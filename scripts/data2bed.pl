#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use Bio::ToolBox::Data::Stream;
use Bio::ToolBox::utility;
use Bio::ToolBox::big_helper qw(bed_to_bigbed_conversion);
my $VERSION =  '1.61';

print "\n This program will write a BED file\n";

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
	$chr_index,
	$start_index,
	$stop_index,
	$name,
	$score_index,
	$strand_index,
	$zero_based,
	$no_header,
	$ask,
	$bigbed,
	$bb_app_path,
	$database,
	$chromo_file,
	$gz,
	$bgz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'i|in=s'           => \$infile, # the solexa data file
	'o|out=s'          => \$outfile, # name of output file 
	'H|noheader'       => \$no_header, # source has no header line
	'a|ask'            => \$ask, # request help in assigning indices
	'c|chr=i'          => \$chr_index, # index of the chromosome column
	'b|begin|start=i'  => \$start_index, # index of the start position column
	'e|stop|end=i'     => \$stop_index, # index of the stop position coloumn
	'n|name=s'         => \$name, # index for the name column
	's|score=i'        => \$score_index, # index for the score column
	't|strand=i'       => \$strand_index, # index for the strand column
	'0|zero!'          => \$zero_based, # source is 0-based numbering, convert
	'B|bigbed|bb'      => \$bigbed, # generate a binary bigbed file
	'd|db=s'           => \$database, # database for bigbed file generation
	'chromof=s'        => \$chromo_file, # name of a chromosome file
	'bbapp=s'          => \$bb_app_path, # path to bedToBigBed utility
	'z|gz!'            => \$gz, # compress output
	'Z|bgz!'           => \$bgz, # compress with bgzip
	'h|help'           => \$help, # request help
	'v|version'        => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script data2bed.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}



### Check for requirements
unless ($infile) {
	$infile = shift @ARGV or
		die " no input file! use --help for more information\n";
}
if ($bigbed and ($gz or $bgz)) {
	# do not allow compression even if requested when converting to bigbed
	warn " compression not allowed when converting to BigBed format\n";
	$gz = 0;
	$bgz = 0;
}


# define name base or index
my $name_index;
my $name_base;
if (defined $name) {
	if ($name =~ /^(\d+)$/) {
		# looks like an index was provided
		$name_index = $1;
	}
	elsif ($name =~ /(\w+)/i) {
		# text that will be used as the name base when autogenerating
		$name_base = $1;
	}
}



### Load file
my $Input = Bio::ToolBox::Data::Stream->new(file => $infile, noheader => $no_header) or
	die "Unable to open file '$infile'!\n";

if ($bigbed) {
	# identify database if needed
	unless ($database or $chromo_file) {
		$database = $Input->database or 
			die " No database name or chromosome file provided for generating bigbed file!\n";
	}
}



### Determine indices

# Ask user interactively
if ($ask) {
	# the user has specified that we should ask for specific indices
	print " Press Return to accept the suggested index\n";
	
	# request chromosome index
	unless (defined $chr_index) {
		my $suggestion = $Input->chromo_column;
		$chr_index = ask_user_for_index($Input, 
			" Enter the index for the chromosome column [$suggestion]  ");
		$chr_index = defined $chr_index ? $chr_index : $suggestion;
		unless (defined $chr_index) {
			die " No identifiable chromosome column index!\n";
		}
	}
	
	# request start index
	unless (defined $start_index) {
		my $suggestion = $Input->start_column;
		$start_index = ask_user_for_index($Input, 
			" Enter the index for the start column [$suggestion]  ");
		$start_index = defined $start_index ? $start_index : $suggestion;
		unless (defined $start_index) {
			die " No identifiable start position column index!\n";
		}
	}
	
	# request stop index
	unless (defined $stop_index) {
		my $suggestion = $Input->stop_column;
		$stop_index = ask_user_for_index($Input, 
			" Enter the index for the stop or end column [$suggestion]  ");
		$stop_index = defined $stop_index ? $stop_index : $suggestion;
		unless (defined $stop_index) {
			die " No identifiable stop position column index!\n";
		}
	}
	
	# request name index or text
	unless (defined $name) {
		# this is a special input, can't use the ask_user_for_index sub
		# accepts either index or text string
		my $suggestion = $Input->name_column;
		print " Enter the index for the feature name column or\n" . 
			"   the base text for auto-generated names [$suggestion]  ";
		my $in = <STDIN>;
		if ($in =~ /^(\d+)$/) {
			$name_index = $1;
		}
		elsif ($in =~ /(\w+)/) {
			$name_base = $1;
		}
		elsif (defined $suggestion) {
			$name_index = $suggestion;
		}
	}
	
	# request score index
	unless (defined $score_index) {
		my $suggestion = $Input->find_column('^score$');
		$score_index = ask_user_for_index($Input, 
			" Enter the index for the feature score column [$suggestion]  ");
		$score_index = defined $score_index ? $score_index : $suggestion;
	}
	
	# request strand index
	unless (defined $strand_index) {
		my $suggestion = $Input->strand_column;
		$strand_index = ask_user_for_index($Input, 
			" Enter the index for the feature strand column [$suggestion]  ");
		$strand_index = defined $strand_index ? $strand_index : $suggestion;
	}
}
else {
	# or else the indices need to be automatically identified
	unless (
		$Input->feature_type eq 'coordinate' or 
		(defined $chr_index and defined $start_index and defined $stop_index) or 
		($Input->feature_type eq 'named' and ($database or $Input->database)) 
	) {
		die "Not enough information has been provided to convert to bed file.\n" . 
			"Coordinate column names must be recognizable or specified. Use --help\n";
	}
}


# Open output stream
unless ($outfile) {
	$outfile = $Input->path . $Input->basename;
}
my $Output = Bio::ToolBox::Data::Stream->new(
	out     => $outfile,
	bed     => 6,
	gz      => $bgz ? 2 : $gz ? 1 : 0,
) or die " unable to create output file $outfile!";



# Convert the input stream
printf " Converting using \n  - chromosome index %s\n  - start index %s\n" . 
	"  - stop index %s\n  - name index %s\n  - score index %s\n  - strand index %s\n",
	defined $chr_index ? $chr_index : defined $Input->chromo_column ? $Input->chromo_column : '-',
	defined $start_index ? $start_index : defined $Input->start_column ? $Input->start_column : '-',
	defined $stop_index ? $stop_index : defined $Input->stop_column ? $Input->stop_column : '-',
	defined $name_index ? $name_index : defined $name_base ? $name_base : defined $Input->name_column ? $Input->name_column : '-',
	defined $score_index ? $score_index : '-',
	defined $strand_index ? $strand_index : defined $Input->strand_column ? $Input->strand_column : '-';

my $count = 0; # the number of lines processed
my $notsorted = 0; # flag to set 
my ($prev_chr, $prev_start);
if (defined $start_index and substr($Input->name($start_index), -1) eq '0') {
	$zero_based = 1; # name suggests 0-based indexing
}
my $do_feature = $Input->feature_type eq 'named' ? 1 : 0; # get features from db?
while (my $row = $Input->next_row) {
	
	# get the feature from the db if necessary
	my $f = $row->feature if $do_feature;
	
	# check if it's sorted
	if ($bigbed and not $notsorted) {
		$prev_chr ||= $row->value($chr_index);
		$prev_start ||= $row->value($start_index);
		if ($row->value($chr_index) ne $prev_chr) {
			$prev_chr = $row->value($chr_index);
			$prev_start = 1;
		}
		$notsorted++ if $row->value($start_index) < $prev_start;
	}
	
	# build the arguments
	# retrieve information from row object if indices were provided
	my @args;
	if (defined $chr_index) {
		push @args, 'chromo', $row->value($chr_index);
	}
	if (defined $start_index) {
		push @args, 'start';
		push @args, $zero_based ? $row->value($start_index) + 1 : 
			$row->value($start_index);
	}
	if (defined $stop_index) {
		push @args, 'stop', $row->value($stop_index);
	}
	if (defined $strand_index) {
		push @args, 'strand', $row->value($strand_index);
	}
	if (defined $score_index) {
		push @args, 'score', $row->value($score_index);
	}
	if (defined $name_index) {
		push @args, 'name', $row->value($name_index);
	} elsif (defined $name_base) {
		push @args, 'name', sprintf("%s_%07d", $name_base, $count);
	}
			
	# write
	$Output->add_row( $row->bed_string(@args) );
	$count++;
}
$Output->close_fh;


### Convert to BigBed format
if ($bigbed) {
	# requested to continue and generate a binary bigbed file
	printf " wrote %s lines to temporary bed file '%s'\n", 
		format_with_commas($count), $Output->filename;
	print " converting to bigbed file....\n";
	
	# perform the conversion
	my $bb_file = bed_to_bigbed_conversion(
			'bed'       => $Output->filename,
			'db'        => $database,
			'chromo'    => $chromo_file,
			'bbapppath' => $bb_app_path,
	);

	# confirm
	if ($bb_file) {
		print " BigBed file '$bb_file' generated\n";
		unlink $outfile; # remove the bed file
	}
	else {
		die " BigBed file not generated! see standard error\n";
	}
	
}
else {
	printf " wrote %s lines to BED file '%s'\n", format_with_commas($count), 
		$Output->filename;
}


__END__

=head1 NAME

data2bed.pl

A program to convert a data file to a bed file.

=head1 SYNOPSIS

data2bed.pl [--options...] <filename>
  
  File Options:
  -i --in <filename>                    input file: txt, gff, vcf, etc
  -o --out <filename>                   output file name
  -H --noheader                         input file has no header row
  -0 --zero                             file is in 0-based coordinate system
  
  Column indices:
  -a --ask                              interactive selection of columns
  -c --chr <index>                      chromosome column
  -b --begin --start <index>            start coordinate column
  -e --end --stop <index>               stop coordinate column
  -n --name <text | index>              name column or base name text
  -s --score <index>                    score column
  -t --strand <index>                   strand column
  
  BigBed options:
  -B --bb --bigbed                      generate a bigBed file
  -d --db <database>                    database to collect chromosome lengths
  --chromof <filename>                  specify a chromosome file
  --bwapp </path/to/bedToBigBed>        specify path to bedToBigBed
  
  General Options:
  -z --gz                               compress output file
  -Z --bgz                              bgzip compress output file
  -v --version                          print version and exit
  -h --help                             show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 File Options

=over 4

=item --in E<lt>filenameE<gt>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. Genome 
coordinates are required. The first row should be column headers. Text 
files generated by other B<BioToolBox> scripts are acceptable. Files may 
be gzipped compressed.

=item --out E<lt>filenameE<gt>

Specify the output filename. By default it uses the basename of the input 
file.

=item --noheader

The input file does not have column headers, often found with UCSC 
derived annotation data tables. 

=item --zero

Indicate that the source data is already in interbase (0-based) 
coordinates and do not need to be converted. By convention, all 
BioPerl (and, by extension, all biotoolbox) scripts are base 
(1-based) coordinates. Default behavior is to convert.

=back

=head2 Column indices

=over 4

=item --ask

Indicate that the program should interactively ask for the indices for 
feature data. It will present a list of the column 
names to choose from. Enter nothing for non-relevant columns or to 
accept default values.

=item --chr E<lt>column_indexE<gt>

The index of the dataset in the data table to be used 
as the chromosome or sequence ID column in the BED data.

=item --start E<lt>column_indexE<gt>

=item --begin E<lt>column_indexE<gt>

The index of the dataset in the data table to be used 
as the start position column in the BED data.

=item --start E<lt>column_indexE<gt>

=item --end E<lt>column_indexE<gt>

The index of the dataset in the data table to be used 
as the stop or end position column in the BED data.

=item --name E<lt>column_index | base_textE<gt>

Supply either the index of the column in the data table to 
be used as the name column in the BED data, or the base text 
to be used when auto-generating unique feature names. The 
auto-generated names are in the format 'text_00000001'. 
If the source file is GFF3, it will automatically extract the 
Name attribute.

=item --score E<lt>column_indexE<gt>

The index of the dataset in the data table to be used 
as the score column in the BED data.

=item --strand E<lt>column_indexE<gt>

The index of the dataset in the data table to be used
for strand information. Accepted values might include
any of the following: +, -, 1, -1, 0, .

=back

=head2 BigBed options

=over 4

=item --bigbed

=item --bb

Indicate that a binary BigBed file should be generated instead of 
a text BED file. A .bed file is first generated, then converted to 
a .bb file, and then the .bed file is removed.

=item --db E<lt>databaseE<gt>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
or other indexed data file, e.g. Bam or bigWig file, from which chromosome 
length information may be obtained. For more information about using databases, 
see L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. It 
may be supplied from the input file metadata.

=item --chromf E<lt>filenameE<gt>

When converting to a BigBed file, provide a two-column tab-delimited 
text file containing the chromosome names and their lengths in bp. 
Alternatively, provide a name of a database, below.

=item --bbapp </path/to/bedToBigBed>

Specify the path to the UCSC bedToBigBed conversion utility. The 
default is to first check the BioToolBox  configuration 
file C<biotoolbox.cfg> for the application path. Failing that, it will 
search the default environment path for the utility. If found, it will 
automatically execute the utility to convert the bed file.

=back

=head2 General options

=over 4

=item --gz

Specify whether the output file should be compressed with gzip.

=item --bgz

Specify whether the output file should be compressed with block gzip 
(bgzip) for tabix compatibility.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will convert a tab-delimited data file into a BED file,
according to the specifications here
L<http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED>. A minimum 
of three and a maximum of six columns may be generated. Thin and
thick block data (columns greater than 6) are not written. 

Column identification may be specified on the command line, chosen 
interactively, or automatically determined from the column headers. GFF 
source files should have columns automatically identified. 

All lower-numbered columns must be defined before writing higher-numbered 
columns, as per the specification. Dummy data may be filled in for 
Name and/or Score if a higher column is requested. 

Browser and Track lines are not written. 

Following specification, all coordinates are written in interbase
(0-based) coordinates. Base (1-based) coordinates (the BioPerl standard) 
will be converted. 

Score values should be integers within the range 1..1000. Score values 
are not converted in this script. However, the biotoolbox script 
L<manipulate_datasets.pl> has tools to do this if required.

An option exists to further convert the BED file to an indexed, binary BigBed 
format. Jim Kent's bedToBigBed conversion utility must be available, and 
either a chromosome definition file or access to a Bio::DB database is required.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
