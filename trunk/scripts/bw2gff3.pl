#!/usr/bin/perl

# a script to generate GFF3 files for bigwig and bam files

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Copy;
use File::Spec;
use File::Path 'make_path';
use File::Basename qw(fileparse);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper;
eval {use Bio::DB::BigWig};
eval {use Bio::DB::Sam};

print "\n This script will generate a GFF3 file for a bigwig file\n";


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
	$path,
	$source,
	$type,
	$rename,
	$write_metadata,
	$set_name,
	$write_conf,
	$help
);
my @infiles;
my @names;
my @strands;

# Command line options
GetOptions( 
	'in=s'      => \@infiles, # name of input files
	'path=s'    => \$path, # path to move the bigwig file
	'source=s'  => \$source, # the gff source
	'type=s'    => \$type, # the gff type
	'name=s'    => \@names, # the gff name 
	'strand=s'  => \@strands, # indicate the strand for the feature
	'rename'    => \$rename, # rename the file
	'set!'      => \$write_metadata, # write a metadata index file for BigWigSet
	'setname=s' => \$set_name, # name for the bigwigset
	'conf!'     => \$write_conf, # write GBrowse conf stanzas
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



### Check for general required values
# files
if (@infiles) {
	if (scalar @infiles == 1 and $infiles[0] =~ /,/) {
		# a comma delimited list is provided
		my $file = shift @infiles;
		@infiles = split /,/, $file;
	}
}
else {
	# file list was provided on the command line
	@infiles = @ARGV or
		die "  OOPS! No source data files specified! \n use --help\n";
}

# names
if (scalar @names == 1 and $names[0] =~ /,/) {
	# a comma delimited list is provided
	my $name = shift @names;
	@names = split /,/, $name;
	if (scalar @names != scalar @infiles) {
		die " unequal number of names (" . scalar(@names) . ") and files (" . 
			scalar(@infiles) . ") provided!\n";
	}
}

# strands
if (scalar @strands == 1 and $strands[0] =~ /,/) {
	# a separate strand value given for each file as comma delimited list
	my $strand = shift @strands;
	@strands = split /,/, $strand;
	if (scalar @strands != scalar @infiles) {
		die " unequal number of strands (" . scalar(@strands) . ") and files (" . 
			scalar(@infiles) . ") provided!\n";
	}
}
elsif (scalar @strands == 1) {
	# a single strand value
	# assign to each input file name
	my $strand = shift @strands;
	foreach (@infiles) {
		push @strands, $strand;
	}
}
if (@strands) {
	# check the strand values
	foreach my $strand (@strands) {
		unless ($strand =~ /^[frwc\.\-\+]$/) {
			die " unknown strand value '$strand'!\n use --help\n";
		}
	}
}


# target directory
if (defined $path) {
	$path = File::Spec->canonpath($path);
	unless (-e $path) {
		make_path($path) or die "unable to generate target directory: '$path'";
	}
	unless (-w $path) {
		die " target '$path' doesn't seem to be writeable!\n";
	}
}
else {
	# default is to use the current working directory
	$path = File::Spec->curdir();
}

# my bigwigset name
unless ($set_name) {
	# we'll assume that the last directory in the defined path is the 
	# name we'll use for the set
	my @dirs = File::Spec->splitdir($path);
	$set_name = $dirs[-1];
}

# default source
unless ($source) {
	$source = 'data';
}


# prepare metadata and conf output arrays
	# we'll be dumping the metadata in here for writing later after going 
	# through all the input files
	# it will be a simple array of lines to be written (or appended) to a 
	# metadata index file written to the output file path
my @metadata;
	# for the conf array, we'll either be writing individual conf stanzas
	# for each bigwig file or one stanza for the bigwigset
my @confdata;



### Processing each of the input files
while (@infiles) {
	
	# get infile name
	my $infile = shift @infiles;
	print "\n processing '$infile'...\n";
	
	
	### Determine file specific GFF variables
	# determine input file components
	my ($infile_basename, $infile_path, $infile_ext) = 
		fileparse($infile, '.bw', '.bigwig', '.bam');
	
	# determine gff name
	my $name;
	if (@names) {
		# hopefully the same number of names was provided as files!
		$name = shift @names;
	}
	else {
		# default is to use the base name
		$name = $infile_basename;
	}
	$name =~ s/-|\./_/g; # substitute any dashes or periods with underscores
	
	# determine gff type
	my $gfftype;
	if ($type) {
		$gfftype = $type;
	}
	else {
		# default to use the name
		$gfftype = $name;
	}
	$gfftype =~ s/-|\./_/g; # substitute any dashes or periods with underscores
	
	# determine strand
	my $strand;
	if (@strands) {
		$strand =  shift @strands;
	}
	else {
		# default is no strand
		$strand = '.';
	}
	
	### Generate the GFF data structure
	# Initialize the data structure
	my $main_data_ref = initialize_data_structure();
	
	# Determine the target file name
	my $target_basename;
	if ($rename) {
		$target_basename = "$source.$name";
	}
	else {
		$target_basename = $infile_basename;
	}
	my $target_file = 
		File::Spec->catfile($path, "$target_basename$infile_ext");
	
	# Load the chromosome data
	if ($infile_ext eq '.bw' or $infile_ext eq '.bigwig') {
		# source data is a bigwig file
		collect_chromosomes_from_bigwig(
			$infile, $target_file, $main_data_ref, $strand);
	}
	elsif ($infile_ext eq '.bam') {
		# source data is a bam file
		collect_chromosomes_from_bam(
			$infile, $target_file, $main_data_ref, $strand);
	}
	else {
		die " unrecognized input file format!\n";
	}
	
	
	### Move the Input files
	if (-e $target_file) {
		print "  target file already present\n";
	}
	else {
		# copy the file if isn't there already
		copy($infile, $target_file);
		
		# check
		if (-e $target_file) {
			print "  Copied target file '$target_file'\n";
		}
		else {
			warn "  attempted to copy target file but failed!?\n";
		}
	}
	
	
	### Write the GFF3 file
	# if we're creating a bigwigset, then we'll be writing a single GFF file
	# otherwise write individual GFF files for each bigwig file
	if ($write_metadata) {
		# writing a single gff file
		
		# first convert the data structure to GFF for writing
		convert_genome_data_2_gff_data( {
			'data'       => $main_data_ref,
			'version'    => 3,
			'source'     => $source,
			'type'       => $gfftype,
			'name'       => $name,
			'strand'     => 3,
			'tags'       => [4],
		} ) or die " unable to convert data to GFF format!\n";
		
		# write new or append existing GFF file
		my $gff_file = $set_name . '.gff3';
		if (-e $gff_file) {
			# file exists, append to it
			my $gff_fh = open_to_write_fh($gff_file, 0, 1) or 
				die " can't append to GFF file!\n";
			
			# append the table contents to the file
			for my $row (1 .. $main_data_ref->{'last_row'}) {
				print {$gff_fh} join("\t", 
					@{ $main_data_ref->{'data_table'}->[$row] }), "\n";
			}
			$gff_fh->close;
			print "  appended to GFF3 file '$gff_file'\n";
		}
		else {
			# file doesn't exist, write new one
			my $success = write_tim_data_file( {
				'data'       => $main_data_ref,
				'filename'   => $gff_file,
			} );
			if ($success) {
				print "  wrote GFF3 file '$success'\n";
			}
			else {
				print "  unable to write output file!\n";
			}
		}
	}
	else {
		# writing individual GFF files
		
		# convert and write the GFF file
		my $success = convert_and_write_to_gff_file( {
			'data'       => $main_data_ref,
			'filename'   => $name,
			'version'    => 3,
			'source'     => $source,
			'type'       => $gfftype,
			'name'       => $name,
			'strand'     => 3,
			'tags'       => [4],
		} );
		if ($success) {
			print "  wrote GFF3 file '$success'\n";
		}
		else {
			print "  unable to write output file!\n";
		}
	}
	
	
	### Add metadata to index file
	if ($write_metadata) {
		# we'll be adding all of the input bigwig files to the metadata index
		
		# cannot support bam data
		if ($infile_ext eq '.bam') {next}
		
		# add metadata header block
			# this is the file name in square brackets
		push @metadata, "[$target_basename$infile_ext]\n";
		
		# add metadata
		push @metadata, "primary_tag  = $gfftype\n";
		push @metadata, "source       = $source\n";
		push @metadata, "display_name = $name\n";
		if ($strand =~ /^f|w|\+|1/) {
			push @metadata, "strand       = 1\n";
		}
		elsif ($strand =~ /^r|c|\-/) {
			push @metadata, "strand       = -1\n";
		}
		push @metadata, "\n"; # empty line to make things purdy
	}
	
	### Write conf stanza data
	if ($write_conf and !$write_metadata) {
		# we are using the write metadata boolean variable as an indicator 
		# of whether to write an individual conf stanza for each bigwig 
		# file or write a single stanza for the bigwig set
		
		# here we write individual stanzas for each bigwig file
		
		# add the database stanza
		push @confdata, "[$target_basename\_db:database]\n";
		push @confdata, "db_adaptor   = Bio::DB::BigWig\n";
		push @confdata, "db_args      = -bigwig $target_file\n";
		
		# add the basic track stanza
		push @confdata, "\n[$name]\n";
		push @confdata, "database     = $target_basename\_db\n";
		push @confdata, "feature      = summary\n";
		push @confdata, "glyph        = wiggle_whiskers\n";
		push @confdata, "graph_type   = boxes\n";
		push @confdata, "\n";
	}
}


### Write metadata file
if ($write_metadata) {
	my $md_file = File::Spec->catfile($path, 'metadata.index');
	my $md_fh;
	if (-e $md_file) {
		# file already exists!?
		# append to it
		$md_fh = open_to_write_fh($md_file, 0, 1) or 
			die " can't append to metadata index file!\n";
	}
	else {
		# write a new file
		$md_fh = open_to_write_fh($md_file) or 
			die " can't write metadata index file!\n";
	}
	foreach (@metadata) {
		print {$md_fh} $_;
	}
	$md_fh->close;
	print " wrote metadata index file '$md_file'\n";
}



### write the starter conf data
if ($write_conf) {
	# Generate the stanzas for a bigwigset if necessary
	if ($write_metadata) {
		# write a conf stanza for the bigwig set
		
		# add the database stanza
		push @confdata, "[$set_name\_db:database]\n";
		push @confdata, "db_adaptor   = Bio::DB::BigWigSet\n";
		push @confdata, "db_args      = -dir $path\n";
		
		# add the basic track stanza
		push @confdata, "\n[$set_name]\n";
		push @confdata, "database     = $set_name\_db\n";
		push @confdata, "feature      = $type\n";
		push @confdata, "glyph        = wiggle_whiskers\n";
		push @confdata, "graph_type   = boxes\n";
		push @confdata, "\n";
	}
	
	# write the conf stanza file
	my $conf_file = 'conf_stanzas.txt';
	my $conf_fh;
	if (-e $conf_file) {
		# file already exists!?
		# append to it
		$conf_fh = open_to_write_fh($conf_file, 0, 1) or 
			die " can't append to GBrowse configuration stanza file!\n";
	}
	else {
		# write a new file
		$conf_fh = open_to_write_fh($conf_file) or 
			die " can't write GBrowse configuration stanza file!\n";
	}
	foreach (@confdata) {
		print {$conf_fh} $_;
	}
	$conf_fh->close;
	print " wrote sample GBrowse configuration stanza file '$conf_file'\n";
}
print " Finished\n";






#############################  Subroutines  ################################

sub initialize_data_structure {
	
	# generate the data hash
	my %datahash;
	
	# populate the standard data hash keys
	$datahash{'program'}        = $0;
	$datahash{'feature'}        = 'data_features';
	$datahash{'number_columns'} = 5;
	
	# set column metadata
	$datahash{0} = {
		# the chromosome
		'name'     => 'Chromosome',
		'index'    => 0,
	};
	$datahash{1} = {
		# the start position 
		'name'     => 'Start',
		'index'    => 1,
	};
	$datahash{2} = {
		# the stop position
		'name'     => 'Stop',
		'index'    => 2,
	};
	$datahash{3} = {
		# strand
		'name'     => 'Strand',
		'index'    => 3,
	};
	$datahash{4} = {
		# Name
		'name'     => 'File',
		'index'    => 4,
	};
	
	
	# Set the data table
	my @data_table = ( [ qw(
		Chromosome
		Start
		Stop
		Strand
		File
	) ] );
	$datahash{'data_table'} = \@data_table;
	
	# return the reference to the generated data hash
	return \%datahash;
}



sub collect_chromosomes_from_bigwig {
	my ($infile, $target_file, $data_ref, $strand) = @_;
	
	# adjust the File column name
	$data_ref->{4}{'name'} = 'bigwigfile';
	$data_ref->{'data_table'}->[0][4] = 'bigwigfile';
	
	# open the bigwig file
	my $wig = Bio::DB::BigWig->new(-bigwig => $infile) or 
		die " unable to open bigwig data file!\n";
	
	
	# collect the chromosomes
	my @seq_ids = $wig->seq_ids;
	unless (@seq_ids) {
		die " no chromosomes in the bigwig file!?\n";
	}
	
	# fill out the data table
	foreach my $seq_id (@seq_ids) {
		push @{ $data_ref->{'data_table'} }, [ (
			$seq_id,
			1,
			$wig->length($seq_id),
			$strand,
			$target_file,
		) ];
	}
	
	# update
	$data_ref->{'last_row'} = scalar(@seq_ids);
	
	# done
	undef $wig;
}



sub collect_chromosomes_from_bam {
	my ($infile, $target_file, $data_ref, $strand) = @_;
	
	# adjust the File column name
	$data_ref->{4}{'name'} = 'bamfile';
	$data_ref->{'data_table'}->[0][4] = 'bamfile';
	
	# open the bam file
	my $sam = Bio::DB::Sam->new(-bam => $infile) or 
		die " unable to open bam data file!\n";
	
	
	# collect the chromosomes
	my $seq_num = $sam->n_targets;
	unless ($seq_num) {
		die " no chromosomes in the bam file!?\n";
	}
	
	# fill out the data table
	foreach my $tid (0 .. ($seq_num - 1) ) {
		push @{ $data_ref->{'data_table'} }, [ (
			$sam->target_name($tid),
			1,
			$sam->target_len($tid),
			$strand,
			$target_file,
		) ];
	}
	
	# update
	$data_ref->{'last_row'} = $seq_num;
	
	# done
	undef $sam;
}




__END__

=head1 NAME

bw2gff3.pl

=head1 SYNOPSIS

bw2gff3.pl [--options...] <filename1.bw> <filename2.bw> ...
  
  --in <file> or <file1,file2,...>
  --path </destination/path/for/bigwig/files/>
  --source <text>
  --name <text> or <text1,text2,...>
  --type <text>
  --strand [f|r|w|c|+|-|1|0|-1],...
  --rename
  --set
  --conf
  --help
  

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename.bw | filename.bam>

Provide the name(s) of the input bigwig. The list may be 
specified by reiterating the --in option, providing a single comma-delimited 
list, or simply listing all files after the options (e.g. "data_*.bw"). 
Limited support is provided for .bam files.

=item --path </destination/path/for/bigwig/files/>

Provide the destination directory name for the bigwig files. If the
destination does not exist, then it will created. This directory should be
writeable by the user and readable by all (or at least the Apache and MySQL
users). If the input files are not currently located here, they will be
copied to the directory for you. Note that when generating a BigWigSet, a
unique directory for just the indicated files should be provided. The
default path is the current path for the input file.

=item --source <text>

Provide the name for the GFF feature's source for all of the input value. 
The default value is "data". Unique values for each input file is not 
supported.

=item --name <text>

Provide the name(s) for the GFF feature(s). This will be used as the GFF 
feature's name. A unique name should be provided for each file, and may be 
specified as a single comma-delimited list or by reiterating the --name 
option. The default value to use the input file basename.

=item --type <text>

Provide the GFF type for the GFF features for all input files. By default, 
it re-uses the GFF name. Unique values for each input file is not supported.

=item --strand [f|r|w|c|+|-|1|0|-1],...

Indicate which strand the feature will be located. Acceptable values include 
(f)orward, (r)everse, (w)atson, (c)rick, +, -, 1, or -1. By default no strand 
is used. For mulitple input files with different strands, use the option 
repeatedly or provide a comma-delimited list. If only one value is provided, 
it is used for all input files.

=item --rename

Rename the input source file basenames as "$source.$name" when moving to the 
target directory. This may help in organization and clarity of file listings. 
Default is false.

=item --set

Indicate that all of the input bigwig files should be part of a BigWigSet,
which treats all of the bigwig files in the target directory as a single
database. This is primarily for GBrowse, as biotoolbox scripts (currently)
don't use the BigWigSet adaptor. A text file is written in the target
directory with metadata for each bigwig file (feature, source, strand,
name) as described in Bio::DB::BigWigSet documentation. Additional metadata 
may be manually added as necessary. The default is false.

=item --conf

Write sample GBrowse database and track configuration stanzas. Each bigwig 
file will get individual stanzas, unless the --set option is enabled, where 
a single stanza for the BigWigSet is provided. This is helpful when setting 
up GBrowse database and configurations. Default is false.

=item --help

Display this POD help.

=back

=head1 DESCRIPTION

This program will generate a GFF3 file with features representing a 
bigwig file. The features will encompass the entire length of each 
chromosome represented in the data file. The name of the data file and its 
absolute path is stored as a tag value in each feature. This tag value can 
be used by B<tim_db_helper.pm> to collect data from the file with respect to 
various locations and features in the database. 

The source data file is copied to the destination directory. The file may be 
renamed as "source.name" so as to avoid confusion when lots of files are 
dumped into the same directory. 

The bigwig files may also be designated as a BigWigSet, with unique metadata 
assigned to each file (source, type, name, strand). The BigWigSet may be 
treated as a single database with multiple bigwig data sources by GBrowse. A 
metadata index file is written in the target directory as described in 
Bio::DB::BigWigSet.

Multiple source files may be designated, and each may have its own name. 
This facilitates multiple file processing.

The generated GFF3 file is written to the current directory. One GFF file is 
written for each input file. It uses the provided GFF name as the basename 
for the file.

Optionally, sample database and track GBrowse configuration stanzas may also be 
written to the current directory to facilitate setting up GBrowse.



=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112



