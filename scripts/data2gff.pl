#!/usr/bin/perl

# documentation at end of file

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use IO::Prompt::Tiny qw(prompt);
use Bio::ToolBox::Data;
use Bio::ToolBox::utility qw(parse_list ask_user_for_index format_with_commas);

our $VERSION = '2.03';

print "\n This script will convert a data file to a GFF\n\n";

### Quick help
unless (@ARGV) {

	# when no command line options are present
	# print SYNOPSIS
	pod2usage(
		{
			'-verbose' => 0,
			'-exitval' => 1,
		}
	);
}

### Get command line options and initialize values
my (
	$infile,     $outfile,     $no_header,    $chr_index, $start_index,
	$stop_index, $score_index, $strand_index, $use_name,  $id_index,
	$source,     $type,        $tag,          $ask,       $unique,
	$interbase,  $sort_data,   $gz,           $bgz,       $help,
	$print_version,
);

# Command line options
GetOptions(
	'i|in=s'          => \$infile,           # specify the input data file
	'o|out=s'         => \$outfile,          # name of output gff file
	'H|noheader'      => \$no_header,        # source has no header line
	'c|chr=i'         => \$chr_index,        # index of the chromosome column
	'b|begin|start=i' => \$start_index,      # index of the start position column
	'e|stop|end=i'    => \$stop_index,       # index of the stop position coloumn
	's|score=i'       => \$score_index,      # index for the score column
	't|strand=i'      => \$strand_index,     # index for the strand column
	'n|name=s'        => \$use_name,         # index for the name column or the name text
	'd|id=i'          => \$id_index,         # index for the ID column
	'r|source=s'      => \$source,           # text to put in the source column
	'y|type=s'        => \$type,             # test to put in the type column
	'g|tag|tags=s'    => \$tag,              # comma list of tag column indices
	'a|ask'           => \$ask,              # request help in assigning indices
	'unique!'         => \$unique,           # make the names unique
	'0|zero!'         => \$interbase,        # input file is interbase format
	'sort!'           => \$sort_data,        # sort the output file
	'z|gz!'           => \$gz,               # boolean to compress output file
	'Z|bgz!'          => \$bgz,              # compress with bgzip
	'h|help'          => \$help,             # request help
	'v|version'       => \$print_version,    # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {

	# print entire POD
	pod2usage(
		{
			'-verbose' => 2,
			'-exitval' => 1,
		}
	);
}

# Print version
if ($print_version) {
	print " Biotoolbox script data2gff.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

### Check for required values
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		print STDERR " FATAL: no input file! use --help for more information\n";
		exit 1;
	}
}
if ($bgz) {
	$gz        = 2;
	$sort_data = 1;
}

# define name base or index
my ( $name_index, $name_base );
if ( defined $use_name ) {
	if ( $use_name =~ /^(\d+)$/ ) {

		# looks like an index was provided
		$name_index = $1;
	}
	elsif ( $use_name =~ /(\w+)/i ) {

		# text that will be used as the name base when autogenerating
		$name_base = $1;
	}
}

# define type base or index
my ( $type_index, $type_base );
if ( defined $type ) {
	if ( $type =~ /^(\d+)$/ ) {

		# looks like an index was provided
		$type_index = $1;
	}
	elsif ( $type =~ /(\w+)/i ) {

		# text that will be used as the type base when autogenerating
		$type_base = $1;
	}
}

# define source base or index
my ( $source_index, $source_base );
if ( defined $source ) {
	if ( $source =~ /^(\d+)$/ ) {

		# looks like an index was provided
		$source_index = $1;
	}
	elsif ( $source =~ /(\w+)/i ) {

		# text that will be used as the source base when autogenerating
		$source_base = $1;
	}
}

# gff attribute tag indices
my @tag_indices;
if ($tag) {
	@tag_indices = parse_list($tag);
}

### Load file
my $Input = Bio::ToolBox::Data->new(
	in       => $infile,
	noheader => $no_header,
	stream   => 1,
) or die "Unable to open file '$infile'!\n";

### Determine indices
if ($ask) {

	# the user has specified that we should ask for specific indices
	print " Press Return to accept the suggested index\n";

	# request chromosome index
	unless ( defined $chr_index ) {
		my $suggestion = $Input->chromo_column;
		$chr_index = ask_user_for_index( $Input,
			" Enter the index for the chromosome column [$suggestion]: " );
		$chr_index = defined $chr_index ? $chr_index : $suggestion;
		unless ( defined $chr_index ) {
			print STDERR " FATAL: No identifiable chromosome column index!\n";
			exit 1;
		}
	}

	# request start index
	unless ( defined $start_index ) {
		my $suggestion = $Input->start_column;
		$start_index = ask_user_for_index( $Input,
			" Enter the index for the start column [$suggestion]: " );
		$start_index = defined $start_index ? $start_index : $suggestion;
		unless ( defined $start_index ) {
			print STDERR " FATAL: No identifiable start position column index!\n";
			exit 1;
		}
	}

	# request stop index
	unless ( defined $stop_index ) {
		my $suggestion = $Input->stop_column;
		$stop_index = ask_user_for_index( $Input,
			" Enter the index for the stop or end column [$suggestion]: " );
		$stop_index = defined $stop_index ? $stop_index : $suggestion;
		unless ( defined $stop_index ) {
			print STDERR " FATAL: No identifiable stop position column index!\n";
			exit 1;
		}
	}

	# request source text
	unless ( defined $source ) {

		# this is a special input, can't use the ask_user_for_index sub
		# accepts either index or text string
		my $default = $Input->basename;
		my $p  = " Enter the text string or column index for the GFF source [$default]: ";
		my $in = prompt( $p, $default );
		if ( $in =~ /^(\d+)$/ ) {
			$source_index = $1;
		}
		else {
			$source_base = $in;
		}
	}

	# request type text
	unless ( defined $type ) {

		# this is a special input, can't use the ask_user_for_index sub
		# accepts either index or text string
		my $default = $Input->feature || 'feature';
		my $p  = " Enter the text string or column index for the GFF type [$default]: ";
		my $in = prompt( $p, $default );
		if ( $in =~ /^(\d+)$/ ) {
			$type_index = $1;
		}
		else {
			$type_base = $in;
		}
	}

	# request score index
	unless ( defined $score_index ) {
		my $suggestion = $Input->find_column('^score$');
		$score_index = ask_user_for_index( $Input,
			" Enter the index for the feature score column [$suggestion]: " );
		$score_index = defined $score_index ? $score_index : $suggestion;
	}

	# request strand index
	unless ( defined $strand_index ) {
		my $suggestion = $Input->strand_column;
		$strand_index = ask_user_for_index( $Input,
			" Enter the index for the feature strand column [$suggestion]: " );
		$strand_index = defined $strand_index ? $strand_index : $suggestion;
	}

	# request name index or text
	unless ( defined $use_name ) {

		# this is a special input, can't use the ask_user_for_index sub
		# accepts either index or text string
		my $suggestion = $Input->name_column;
		my $p          = " Enter the index for the feature name column or\n"
			. "   the base text for auto-generated names: ";
		my $in = prompt( $p, $suggestion );
		if ( $in =~ /^(\d+)$/ ) {
			$name_index = $1;
		}
		elsif ( $in =~ /\w+/ ) {
			$name_base = $in;
		}
	}

	# request ID index
	unless ( defined $id_index ) {
		my $suggestion = $Input->id_column;
		$suggestion = q() unless $suggestion;
		$id_index   = ask_user_for_index( $Input,
			" Enter the index for the feature unique ID column [$suggestion]: " );
		$id_index = defined $id_index ? $id_index : $suggestion;
	}

	# request tags
	unless ( defined $tag ) {
		@tag_indices = ask_user_for_index( $Input,
			" Enter zero or more column indices for GFF group tags  " );
	}
}
else {
	# or else the indices need to be automatically identified
	unless ( $Input->feature_type eq 'coordinate'
		or ( defined $chr_index and defined $start_index and defined $stop_index )
		or ( $Input->feature_type eq 'named' and $Input->database ) )
	{
		print STDERR
			" FATAL: Not enough information has been provided to convert to GFF file.\n"
			. "Coordinate column names must be recognizable or specified. Use --help\n";
		exit;
	}
}

### Open output data
my $Output = Bio::ToolBox::Data->new( gff => 3, )
	or die " unable to create output data strucutre!\n";

### Convert the input stream
# check some things first
my $do_feature = $Input->feature_type eq 'named' ? 1 : 0;    # get features from db?
if ( defined $start_index and substr( $Input->name($start_index), -1 ) eq '0' ) {

	# start column name suggests it is 0-based
	$interbase = 1;
}
if ( $unique and not( defined $name_index or $name_base ) ) {
	print STDERR
		" FATAL: must provide a name index or name base to make unique feature names!\n";
	exit 1;
}
my $unique_name_counter = {};    # hash for making unique feature names
my $unique_id_counter   = {};    # same for IDs
my $count               = 0;     # the number of lines processed
while ( my $row = $Input->next_row ) {

	# get the feature from the db if necessary
	if ($do_feature) {
		my $f = $row->feature;
	}

	# build the arguments
	# retrieve information from row object if indices were provided
	my @args;
	if ($chr_index) {
		my $c = $row->value($chr_index);
		next if $c eq '.';
		push @args, 'chromo', $c;
	}
	if ($start_index) {
		my $s = $row->value($start_index);
		next    if $s eq '.';
		$s += 1 if $interbase;
		push @args, 'start', $s;
	}
	if ($stop_index) {
		push @args, 'stop', $row->value($stop_index);
	}
	if ($strand_index) {
		push @args, 'strand', $row->value($strand_index);
	}
	if ($score_index) {
		push @args, 'score', $row->value($score_index);
	}
	if ($name_index) {
		my $name =
			$unique
			? generate_unique_name( $row->value($name_index), $unique_name_counter )
			: $row->value($name_index);
		push @args, 'name', $name;
	}
	elsif ($name_base) {
		push @args, 'name', sprintf( "%s_%d", $name_base, $count );
	}
	if ($id_index) {
		my $id = $row->value($id_index);
		push @args, 'id', generate_unique_name( $id, $unique_id_counter );
	}
	if ($type_index) {
		push @args, 'type', $row->value($type_index);
	}
	elsif ($type_base) {
		push @args, 'type', $type_base;
	}
	if ($source_index) {
		push @args, 'source', $row->value($source_index);
	}
	elsif ($source_base) {
		push @args, 'source', $source_base;
	}
	if (@tag_indices) {
		push @args, 'attributes', \@tag_indices;
	}

	# add to output
	my $string = $row->gff_string(@args);
	$Output->add_row($string) if length($string);

	# weirdly, this should work, as the add_row will split the columns of
	# the gff string automatically
	$count++;
}

### Finish
$Input->close_fh;
if ($sort_data) {
	print " Sorting data...\n";
	$Output->gsort_data;
}
unless ($outfile) {
	$outfile = sprintf "%s%s.gff", $Input->path, $Input->basename;
}
$outfile = $Output->write_file(
	filename => $outfile,
	gz       => $gz,
);

printf " Converted %s lines of input data to GFF file '%s'\n",
	format_with_commas($count), $outfile;

sub generate_unique_name {
	my ( $name, $counter ) = @_;
	my $new_name;

	# check uniqueness
	if ( exists $counter->{$name} ) {

		# we've encountered this name before
		# generate a unique name by appending the count number
		$counter->{$name} += 1;
		$new_name = $name . '.' . $counter->{$name};
	}
	else {
		# first time for this name
		# record in the hash
		$new_name = $name;
		$counter->{$name} = 0;
	}
	return $new_name;
}

__END__

=head1 NAME

data2gff.pl

A program to convert a generic data file to GFF format.

=head1 SYNOPSIS

data2gff.pl [--options...] <filename>
  
  File options:
  -i --in <filename>                    input file: txt
  -o --out <filename>                   output file name
  -H --noheader                         input file has no header row
  -0 --zero                             file is in 0-based coordinate system
  
  Column indices:
  -a --ask                              interactive selection of columns
  -c --chr <index>                      chromosome column
  -b --begin --start <index>            start coordinate column
  -e --end --stop <index>               stop coordinate column
  -s --score <index>                    score column
  -t --strand <index>                   strand column
  -n --name <text | index>              name column or base name text
  -d --id <index>                       primary ID column
  -g --tags <index,index,...>           zero or more columns for tag attributes
  -r --source <text | index>            source column or text
  -y --type <text | index>              type column or text
  
  General options:
  --unique                              make IDs unique
  --sort                                sort output by genomic coordinates
  -z --gz                               compress output file
  -Z --bgz                              bgzip compress output file
  -v --version                          print version and exit
  -h --help                             show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 File options

=over 4

=item --in E<lt>filenameE<gt>

Specify an input file containing either a list of database features or 
genomic coordinates for which to convert to GFF format. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. Files may 
be gzipped compressed.

=item --out E<lt>filenameE<gt>

Optionally specify the name of of the output file. The default is to use 
the input file base name. The '.gff' extension is automatically
added if required.

=item --noheader

The input file does not have column headers, often found with UCSC 
derived annotation data tables. 

=item --zero

Input file is in interbase or 0-based coordinates. This should be 
automatically detected for most known file formats, e.g. BED.

=back

=head2 Column indices

=over 4

=item --ask

Indicate that the program should interactively ask for column indices or
text strings for the GFF attributes, including coordinates, source, type, 
etc. It will present a list of the column names to choose from. Enter 
nothing for non-relevant columns or to accept default values.

=item --chr E<lt>column_indexE<gt>

The index of the dataset in the data table to be used 
as the chromosome or sequence ID column in the gff data.

=item --start E<lt>column_indexE<gt>

=item --begin E<lt>column_indexE<gt>

The index of the dataset in the data table to be used 
as the start position column in the gff data.

=item --stop E<lt>column_indexE<gt>

=item --end E<lt>column_indexE<gt>

The index of the dataset in the data table to be used 
as the stop or end position column in the gff data.

=item --score E<lt>column_indexE<gt>

The index of the dataset in the data table to be used 
as the score column in the gff data.

=item --strand E<lt>column_indexE<gt>

The index of the dataset in the data table to be used
for strand information. Accepted values might include
any of the following "+, -, 1, -1, 0, .".

=item --name E<lt>text | column_indexE<gt>

Enter either the text that will be shared name among 
all the features, or the index of the dataset in the data 
table to be used as the name of each gff feature. This 
information will be used in the 'group' column.

=item --id E<lt>column_indexE<gt>

The index of the dataset in the data table to be used
as the unique ID of each gff feature. This information
will be used in the 'group' column of GFF v3 files 
only. The default is to automatically generate a 
unique identifier.

=item --tags E<lt>column_indicesE<gt>

Provide a comma delimited list of column indices that contain 
values to be included as group tags in the GFF features. The 
key will be the column name.

=item --source E<lt>text | column_indexE<gt>

Enter either a text string or a column index representing the 
GFF source that should be used for the features. The default is 
'data'.

=item --type E<lt>text | column_indexE<gt>

Enter either a text string or a column index representing the 
GFF 'type' or 'method' that should be used for the features. If 
not defined, it will use the column name for either 
the 'score' or 'name' column, if defined. As a last resort, it 
will use the most creative method of 'Experiment'.

=back

=head2 General options

=over 4

=item --unique

Indicate whether the feature names should be made unique. A count 
number is appended to the name of subsequent features to make them 
unique. Using a base text string for the name will automatically 
generate unique names.

=item --sort

Sort the output file by genomic coordinates. Automatically enabled 
when compressing with bgzip. 

=item --gz

Indicate whether the output file should be compressed with gzip.

=item --bgz

Specify whether the output file should be compressed with block gzip 
(bgzip) for tabix compatibility.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will convert a data file into a GFF version 3 formatted text file. 
Only simple conversions are performed, where each data line is converted 
to a single feature. Complex features with parent-child relationships (such 
as genes) should be converted with something more advanced.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
