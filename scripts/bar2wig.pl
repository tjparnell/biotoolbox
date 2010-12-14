#!/usr/bin/perl

# This script will convert David Nix's bar file to a wig file

use strict;
use Getopt::Long;
use Pod::Usage;
use DirHandle;
use File::Spec;
use Archive::Zip qw( :ERROR_CODES );
use Statistics::Lite qw(mean median sum max);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_db_helper qw(
	$TIM_CONFIG
	open_db_connection
);
use tim_file_helper qw(
	open_to_read_fh
	open_to_write_fh
);
eval {
	# check for bigwig file conversion support
	require tim_db_helper::bigwig;
	tim_db_helper::bigwig->import;
};

print "\n This program will convert bar files to a wig file\n";

### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Get command line options and initialize values
my (
	$infile,
	$outfile,
	$bar_app_path,
	$method,
	$use_track,
	$interbase,
	$bigwig,
	$database,
	$chromo_file,
	$bw_app_path,
	$gz,
	$help
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the directory name of bar files
	'out=s'     => \$outfile, # name of output file 
	'barapp=s'  => \$bar_app_path, # path to David's Bar2Gr file
	'method=s'  => \$method, # method for dealing with duplicate positions
	'track!'    => \$use_track, # boolean to include a track line
	'inter!'    => \$interbase, # shift from interbase to 1-base numbering scheme
	'bw!'       => \$bigwig, # boolean for bigwig output
	'db=s'      => \$database, # name of database to get chromo info
	'chromof=s' => \$chromo_file, # name of a chromosome file
	'bwapp=s'   => \$bw_app_path, # path to wigToBigWig utility
	'gz!'       => \$gz, # boolean to compress output text file
	'help'      => \$help # request help
);

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}





### Check for requirements
unless ($infile) {
	$infile = shift @ARGV or
		die "  No source data file specified! see help\n";
}

# confirm directory
unless (-e $infile and -r $infile) {
	die "  input does not exist or is not readable!\n";
}
my $input_is_dir = 0;
if (-d $infile) {
	# we're working with a directory of files, # otherwise assume a file
	$input_is_dir = 1;
}

# default values
my $method_sub;
if ($method eq 'mean') {
	$method_sub = \&mean();
}
elsif ($method eq 'median') {
	$method_sub = \&median();
}
elsif ($method eq 'sum') {
	$method_sub = \&sum();
}
elsif ($method eq 'max') {
	$method_sub = \&max();
}
else {
	$method = 'sum';
	$method_sub = \&sum;
}
unless (defined $gz) {
	$gz = 0;
}
unless (defined $bigwig) {
	$bigwig = 0;
}
if ($bigwig) {
	unless ($database or $chromo_file) {
		die "  Must define a database name or provide a chromosome file \n" . 
			"to use the bigwig option! see help\n";
	}
}
unless (defined $use_track) {
	if ($bigwig) {
		# if we're generating bigwig file, no track is needed
		$use_track = 0;
	}
	else {
		# default is to to generate a track
		$use_track = 1;
	}
}
unless (defined $interbase) {
	# the USeq and T2 packages typically use interbase (0-base) format
	# default is to shift to 1-base format
	$interbase = 1;
}

# identify application paths
my $java = $TIM_CONFIG->param('applications.java') || `which java` || undef;
unless (defined $java) {
	die " unable to identify java executable!\n";
}
unless ($bar_app_path) {
	$bar_app_path = $TIM_CONFIG->param('applications.Bar2Gr') || undef;
	unless ($bar_app_path) {
		die "  Must define the path to the USeq/T2 java application Bar2Gr! see help\n";
	}
}




### Check input file(s)
print "  checking input file(s)....\n";

# add them to the list
my @to_delete_files;

if ($input_is_dir) {
	# infile is a directory
	# this is typical
	
	# open and read directory
	my $dir = new DirHandle $infile;
	foreach my $file ($dir->read) {
		
		# process only bar files by checking for the bar extension
		if ($file =~ /\.bar/) {
			# check for zip status
			my $newfile = check_bar_file($infile, $file);
				# returns the name of the unzipped file
			
			# record the unzipped file for later cleanup
			if ($newfile) {
				push @to_delete_files, $newfile;
			}
		}
		
	}
	$dir->close;
}

else {
	# a single bar file!? ok, we'll do it
	
	# process bar file
	unless ($infile =~ /\.bar/) {
		die " only .bar files may be processed!\n";
	}
	my $newfile = check_bar_file('.', $infile);
		# returns the name of the unzipped file
	
	# record the unzipped file for later cleanup
	if ($newfile) {
		push @to_delete_files, $newfile;
	}
}






### Convert the bar files

# run David's Bar2Gr utility to convert bar file(s) to text file(s)
print "  converting to gr files...\n";
system $java, '-Xmx1500M', '-jar', $bar_app_path, '-f', $infile;
	# Bar2Gr can accept either a file or a directory of files
	# it will ignore any .bar.zip files present
	# the .gr files it writes use the bar file name as base, appended with 
	# the genome version






### Process the gr files
my @wigfiles;

# A directory of files
if ($input_is_dir) {
	# working with directory of files
	
	# arrays for stranded files
	my @f_files;
	my @r_files;
	my @n_files;
	
	my $dir = new DirHandle $infile;
	foreach my $file ($dir->read) {
		
		# find the text gr files
		if ($file =~ /\.gr$/) {
			# we have a gr file
			
			# generate a path
			my $path = File::Spec->catfile($infile, $file);
			
			# identify the strand, if present
			# strand is usually denoted as _+_ or _-_
			if ($file =~ /_\+_/) {
				# forward strand
				push @f_files, $path;
			}
			elsif ($file =~ /_\-_/) {
				# reverse strand
				push @r_files, $path;
			}
			else {
				# no apparent strand
				push @n_files, $path;
			}
			
			# remember for cleanup
			push @to_delete_files, $path;
		}
	}
	$dir->close;
	
	# process
	if (@f_files) {
		my $wigfile = process_gr_files('forward', $infile, @f_files);
		if (defined $wigfile) {
			push @wigfiles, $wigfile;
			print " generated wig file '$wigfile'\n";
		}
	}
	if (@r_files) {
		my $wigfile = process_gr_files('reverse', $infile, @r_files);
		if (defined $wigfile) {
			push @wigfiles, $wigfile;
			print " generated wig file '$wigfile'\n";
		}
	}
	if (@n_files) {
		my $wigfile = process_gr_files('none', $infile, @n_files);
		if (defined $wigfile) {
			push @wigfiles, $wigfile;
			print " generated wig file '$wigfile'\n";
		}
	}

}

else {
	# a single file
	
	# generate the gr file name
	my $grfile = $infile;
	$grfile =~ s/\.zip$//;
	$grfile =~ s/\.bar$/.gr/; # exchange extension
	
	# confirm gr file existence
	unless (-e $grfile) {
		die " unable to find converted file '$grfile'!\n";
	}
	
	# process according to strand
	if ($infile =~ /_\+_/) {
		# forward strand
		my $wigfile = process_gr_files('forward', $infile, $grfile);
		if (defined $wigfile) {
			push @wigfiles, $wigfile;
			print " generated wig file '$wigfile'\n";
		}
	}
	elsif ($infile =~ /_\-_/) {
		# forward strand
		my $wigfile = process_gr_files('reverse', $infile, $grfile);
		if (defined $wigfile) {
			push @wigfiles, $wigfile;
			print " generated wig file '$wigfile'\n";
		}
	}
	else {
		# forward strand
		my $wigfile = process_gr_files('none', $infile, $grfile);
		if (defined $wigfile) {
			push @wigfiles, $wigfile;
			print " generated wig file '$wigfile'\n";
		}
	}
	
	# remember file names
	push @to_delete_files, $grfile;
}



### Convert to bigWig as requested
if ($bigwig and exists &wig_to_bigwig_conversion) {
	
	
	# open database connection if necessary
	my $db;
	if ($database) {
		$db = open_db_connection($database);
	}
	
	# find wigToBigWig utility
	unless ($bw_app_path) {
		# check for an entry in the configuration file
		$bw_app_path = $TIM_CONFIG->param('applications.wigToBigWig') || 
			undef;
	}
	unless ($bw_app_path) {
		# next check the system path
		$bw_app_path = `which wigToBigWig` || undef;
	}
			
	# determine reference sequence type
	my $ref_seq_type = 
		$TIM_CONFIG->param("$database\.reference_sequence_type") ||
		$TIM_CONFIG->param('default_db.reference_sequence_type');
	
	# perform the conversion
	foreach my $wigfile (@wigfiles) {
		
		# conversion
		my $bw_file = wig_to_bigwig_conversion( {
				'wig'       => $wigfile,
				'db'        => $db,
				'seq_type'  => $ref_seq_type,
				'chromo'    => $chromo_file,
				'bwapppath' => $bw_app_path,
		} );
		
		# check whether successful
		# since we're using an external utility, there's no easy way to 
		# capture errors, see Bio::DB::BigFile for more info
		# therefore, we'll simply check for it's existence
		if ($bw_file) {
			print " Converted '$wigfile' to bigWig '$bw_file'\n";
			push @to_delete_files, $wigfile; # we no longer need the wig file
		}
		else {
			print " Failure to convert '$wigfile' to bigWig. See STDERR for errors\n";
			# we won't delete the wigfile
		}
	}
}
elsif ($bigwig and !exists &wig_to_bigwig_conversion) {
	# conversion to bigwig is requested but not supported
	warn " Support for converting to bigwig format is not available\n" . 
		" Please convert manually. See documentation for more info\n";
}



### Finish up
# clean up files
print " cleaning up temp files...";
foreach my $file (@to_delete_files) {
	unlink $file;
}
print " done\n";




########################   Subroutines   ###################################


### Process zipped bar files
sub check_bar_file {
	
	# assemble the path
	my ($directory, $file) = @_;
	my $path = File::Spec->catfile($directory, $file);
	
	if ($path =~ /\.bar.zip$/) {
		# file is zipped
		# usually these zip archives only have one zip member
		# Bar2Gr doesn't appear to read zipped files, 
		# why David? your other programs do....
		# we'll need to extract it
		
		# open zip archive
		my $zip = Archive::Zip->new();
		unless ($zip->read($path) == AZ_OK) {
			die " unable to read zip archive '$path'!\n";
		}
		
		# identify the barfile
		my $barfile = $file;
		$barfile =~ s/\.zip$//; # strip the zip extension
		my $barpath = File::Spec->catfile($directory, $barfile);
		
		# extract the file
		print "   extracting $barfile...\n";
		unless ($zip->extractMember($barfile, $barpath) == AZ_OK) {
			die " unable to extract member '$barfile' from zip archive '$path'!\n";
		}
		
		# return
		return $barpath;
	}
	else {
		# nothing need be done
		return;
	}
}



### process the gr files
sub process_gr_files {
	
	# collect information
	my $strand = shift @_;
	my $path = shift @_;
	my @grfiles = @_;
	
	# determine base name for output file
	my $name;
	if ($outfile) {
		
		# user specified base name
		$name = $outfile;
	}
	else {
		
		# automatically determine base name from input
		if ($input_is_dir) {
			# specified path is a directory
			# need to identify the last element
			my @directories = File::Spec->splitdir($path);
			while (@directories) {
				my $n = pop @directories;
				if ($n ne q() ) {
					# sometimes splitdir gives empty elements, for whatever 
					# reason, so need to avoid those
					$name = $n;
					last;
				}
			}
		}
		else {
			# specified path is a file
			# need to extract the filename
			my ($volume,$directory,$file) = File::Spec->splitpath($path);
			$file =~ s/\.bar(?:\.zip)?$//; # strip extensions
			$name = $file;
		}
	}
	
	# adjust file name for strand if necessary
	my $wigfile = $name;
	if ($strand eq 'forward') {
		$wigfile .= '_f.wig';
	}
	elsif ($strand eq 'reverse') {
		$wigfile .= '_r.wig';
	}
	else {
		$wigfile .= '.wig';
	}
	if ($gz) {
		$wigfile .= '.gz';
	}
	
	# open the wig file for writing
	my $wig = open_to_write_fh($wigfile, $gz) or 
		warn " unable to open wig file for writing!\n";
	return unless $wig;
	
	# write track line
	if ($use_track) {
		print {$wig} "track type=wiggle_0 name=$name\n";
	}
	
	
	# process the gr files
	foreach my $file (@grfiles) {
		
		# first identify the chromosome
		# this is encoded in the gr file as the first element of the file name
		my $chr;
		if ($file =~ /(chr[\da-zA-Z]+)/) {
			$chr = $1;
		}
		else {
			warn " unable to identify chromosome by regex in file name for '$file'!\n" . 
				" nothing written!\n";
			next;
		}
		
		# open the file
		my $gr_fh = open_to_read_fh($file);
		unless ($gr_fh) {
			warn " unable to open input gr file! nothing done!\n";
			next;
		}
		
		# write the definition line
		print {$wig} "variableStep chrom=$chr\n";
		
		# set reusuable variables
		my @data;
		my $position;
		
		# collect first line, load into variables
		my $firstline = $gr_fh->getline;
		chomp $firstline;
		($position, $data[0]) = split /\t/, $firstline;
		
		# collect subsequent lines
		while (my $line = $gr_fh->getline) {
			chomp $line;
			my ($pos, $score) = split /\t/, $line;
			if ($position == $pos) {
				# another score value for the current position
				push @data, $score;
			}
			else {
				# we've moved to another position
				# time to write the data
				
				# combine score values
				my $newscore = &{$method_sub}(@data);
				
				# shift position
				$position += 1 if $interbase;
				
				# print
				print {$wig} "$position\t$newscore\n";
				
				# reset
				$position = $pos;
				@data = q();
				push @data, $score;
			}
		}
		
		# final write if necessary
		if (@data) {
			# combine score values
			my $newscore = &{$method_sub}(@data);
			
			# print
			print {$wig} "$position\t$newscore\n";
		}
		
		# close the gr file
		$gr_fh->close;
	}
	
	# finished with this wig file
	$wig->close;
	return $wigfile; # return the actual wig file name
}







__END__

=head1 NAME

bar2wig.pl

=head1 SYNOPSIS

bar2wig.pl [--options...] <filename>
  
  Options:
  --in <filename> or <directory>
  --out <filename> 
  --barapp </path/to/Bar2Gr>
  --method [mean | median | sum | max]
  --(no)track
  --(no)inter
  --bw
  --db <database>
  --chromof <filename>
  --bwapp </path/to/wigToBigWig>
  --(no)gz
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input bar file or (more likely) a directory of bar files. 
The bar files may be zipped (*.bar.zip). 

=item --out <filename>

Specify the output filename. By default it uses the base name of the 
input file or the input directory name. The output file will either 
have a .wig or .bw file extension.

=item --barapp </path/to/Bar2Gr>

Specify the full path to David Nix's USeq or T2 application Bar2Gr (it 
is included in both software packages). By default it uses the path 
defined in the biotoolbox configuration file, tim_db_helper.cfg.

=item --method [mean | median | sum | max]

Define the method used to combine multiple data values at a single 
position. Bar files can frequently have one or more values at a single 
position, usually due to either multiple identical oligo probes (microarray 
data) or multiple sequence tags aligning to the same position (next 
generation sequencing data). Typically, with microarray data, the values 
should be averaged (mean), while with sequencing data, values should be summed.
The default method is 'sum'.

=item --(no)track

Wig files typically include a track line at the beginning of the file which 
defines the appearance. However, conversion of a wig file to bigWig 
requires that the track line is absent. Do not include the track line when 
generating bigWig files manually. The default is to include for wig files 
(true), and exclude for bigWig files (false).

=item --(no)inter

Bar files are typically in interbase (0-base) format, whereas wig files are 
typically 1-base format. Shift (or not) the position by 1 bp. Default is 
true.

=item --bw

Flag to indicate that a binary bigWig file should be generated rather than 
a text wig file.

=item --db <database>

Specify the name of a Bio::DB database from which to extract chromosome 
names and sizes. This information is required when generating a bigWig 
file.

=item --chromof <filename>

Alternative to the --db argument, a pre-generated chromosome sizes text 
file may be specified. This text file should consist of two columns, 
delimited by whitespace, consisting of the chromosome name and size.

=item --bwapp </path/to/wigToBigWig>

Specify the full path to Jim Kent's wigToBigWig conversion utility. By 
default it uses the path defined in the biotoolbox configuration file, 
tim_db_helper.cfg. If it is not defined here or in the config file, then 
the system path is searched for the executable. Finally, failing that, 
it will attempt to use Lincoln Stein's Bio::DB::BigFile module for 
conversion, if available.

=item --(no)gz

Specify whether (or not) the output wig file should be compressed with gzip.
This option does not affect bigWig files, but will affect intermediate 
wig files.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will convert binary .bar files to a conventional 
text-based wig file, with the option to further convert to a binary 
bigWig format. 

The bar file is a Java serialized, floating-point, binary file used
natively by David Nix's TiMAT2 (T2) and USeq software packages. It allows
for effective visualization by the Integrated Genome Browser (IGB), but
very few (any?) other genome browsers.

The wig file is a commonly used text-based file format for exchanging and 
representing dense genomic data (microarray and sequencing) that is read 
by nearly all genome browsers. In most genome browsers, the data values are 
downsampled to 8-bit precision (0-255), sufficient for visualization but 
limited for data analysis. 

The bigWig format is an indexed, compressed, binary extension of the wig
format. It maintains numeric precision, allows for rapid statistical
summaries, and can be rapidly and randomly accessed either locally or
remotely. It is supported by the UCSC and GBrowse genome browsers as well 
as biotoolbox scripts.

This program first uses David Nix's Bar2Gr application to convert the bar 
file to a very simplified text format, a .gr file. This is then further 
processed into one or more wig (or bigWig) files. Stranded data (denoted 
by _+_ and _-_ in the bar file names) is written to two stranded output 
files, appended with either '_f' or '_r' to the basename. The wig files 
are in variableStep format.

Conversion from wig to bigWig requires Jim Kent's wigToBigWig utility or 
Lincoln Stein's Bio::DB::BigFile support.


=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  


