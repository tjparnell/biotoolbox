#!/usr/bin/perl

# a script to generate GFF3 files for bigwig and bam files

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Copy;
use Cwd 'abs_path';
use File::Basename qw(fileparse);
use Bio::DB::BigWig;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper;

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
	$strand,
	$help
);
my @infiles;
my @names;

# Command line options
GetOptions( 
	'in=s'      => \@infiles, # name of input files
	'path=s'    => \$path, # path to move the bigwig file
	'source=s'  => \$source, # the gff source
	'type=s'    => \$type, # the gff type
	'name=s'    => \@names, # the gff name 
	'strand=s'  => \$strand, # indicate the strand for the feature
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
	if (scalar @infiles == 1) {
		# only one file provided, but may be comma delimited list
		my $file = shift @infiles;
		@infiles = split /,/, $file;
	}
}
else {
	# file list was provided on the command line
	@infiles = @ARGV or
		die "  OOPS! No source data files specified! \n use $0 --help\n";
}

# names
if (scalar @names == 1) {
	# only one name provided, but it may be comma delimited list
	my $name = shift @names;
	@names = split /,/, $name;
}

# assign gff defaults
unless ($source) {
	$source = 'data';
}
if ($strand) {
	# check the value
	unless ($strand =~ /^[frwc\.]$/) {
		die " unknown strand value!\n use $0 --help\n";
	}
}
else {
	# default
	$strand = '.';
}


### Processing each of the input files
while (@infiles) {
	
	# get infile name
	my $infile = shift @infiles;
	
	
	### Determine file specific GFF variables
	# determine input file components
	my $infile_abs_path = abs_path($infile);
	my ($infile_basename, $infile_path, $infile_ext) = 
		fileparse($infile_abs_path, '.bw', '.bigwig', '.bam');
	#print " path components are\n" .
	#	"  path: $infile_path\n  basename: $infile_basename\n  extension: $infile_ext\n";
	
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
	
	# check path
	unless ($path) {
		# use the current infile's path
		$path = $infile_path;
	}
	unless ($path =~ /\/$/) {
		# make sure there is a trailing "/" if path was provided
		# otherwise it's not really recognized as a directory !??
		# this will most certainly destroy functionality on non-unix machines
		$path .= '/';
	}
	unless (-e $path and -w $path) {
		die " path '$path' either doesn't exist or is not writeable!\n";
	}
	
	
	### Generate the GFF data structure
	# Initialize the data structure
	my $main_data_ref = initialize_data_structure();
	
	# Load the chromosome data
	my $target_file = $path . "$source.$name" . $infile_ext;
	print "\n processing '$infile'...\n";
	if ($infile_ext eq '.bw' or $infile_ext eq '.bigwig') {
		# source data is a bigwig file
		collect_chromosomes_from_bigwig($infile, $target_file, $main_data_ref);
	}
	elsif ($infile_ext eq '.bam') {
		# source data is a bam file
		collect_chromosomes_from_bam($infile, $target_file, $main_data_ref);
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
	my ($infile, $target_file, $data_ref) = @_;
	
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
	my ($infile, $target_file, $data_ref) = @_;
	
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
  --strand [f r w c]
  --help
  

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename.bw | filename.bam>

Provide the name(s) of the input bigwig or bam file(s). The list may be 
specified by reiterating the --in option, providing a single comma-delimited 
list, or simply listing all files after the options (e.g. "data_*.bw").

=item --path </destination/path/for/bigwig/files/>

Provide the destination directory name for the bigwig files. This directory 
should be writeable by the user and readable by all. If the bigwig/bam file 
is not currently located here, the program will copy the file to this 
directory for you. The default path is the current path for the input file.

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

=item --strand [f r w c]

Indicate which strand the feature will be located. Acceptable values include 
(f)orward, (r)everse, (w)atson, or (c)rick. By default no strand is used. 

=item --help

Display this POD help.

=back

=head1 DESCRIPTION

This program will generate a GFF3 file with features representing a 
bigwig or bam file. The features will encompass the entire length of each 
chromosome represented in the data file. The name of the data file and its 
absolute path is stored as a tag value in each feature. This tag value can 
be used by L<tim_db_helper.pm> to collect data from the file with respect to 
various locations and features in the database. 

The source data file is copied to the destination directory. The file is 
renamed as "source.name" so as to avoid confusion when lots of files are 
dumped into the same directory. 

Multiple source files may be designated, and each may have its own name. 
This facilitates multiple file processing.

The generated GFF3 file is generated in the current directory. It use the 
provided GFF name as the basename for the file.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112



