#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::Data::Stream;
use Bio::ToolBox::utility;
my $VERSION =  '1.40';

print "\n This script will convert a data file to a GFF\n\n";


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
	$score_index,
	$strand_index,
	$name_index,
	$name,
	$id_index,
	$source,
	$type,
	$tag,
	$ask,
	$unique,
	$interbase,
	$gz,
	$help,
	$print_version,
);


# Command line options
GetOptions( 
	'in=s'      => \$infile, # specify the input data file
	'out=s'     => \$outfile, # name of output gff file 
	'chr=i'     => \$chr_index, # index of the chromosome column
	'start=i'   => \$start_index, # index of the start position column
	'stop|end=i'=> \$stop_index, # index of the stop position coloumn
	'score=i'   => \$score_index, # index for the score column
	'strand=i'  => \$strand_index, # index for the strand column
	'name=s'    => \$name, # index for the name column or the name text
	'id=i'      => \$id_index, # index for the ID column
	'source=s'  => \$source, # text to put in the source column
	'type=s'    => \$type, # test to put in the type column
	'tag|tags=s'=> \$tag, # comma list of tag column indices
	'ask'       => \$ask, # request help in assigning indices
	'unique!'   => \$unique, # make the names unique
	'zero!'     => \$interbase, # input file is interbase format
	'gz!'       => \$gz, # boolean to compress output file
	'help'      => \$help, # request help
	'version'   => \$print_version, # print the version
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
	print " Biotoolbox script data2gff.pl, version $VERSION\n\n";
	exit;
}



### Check for required values
unless ($infile) {
	$infile = shift @ARGV or
		die "  OOPS! No source data file specified! \n use $0 --help\n";
}
unless (defined $gz) {
	$gz = 0;
}

# define name base or index
my ($name_index, $name_base);
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

# define type base or index
my ($type_index, $type_base);
if (defined $type) {
	if ($type =~ /^(\d+)$/) {
		# looks like an index was provided
		$type_index = $1;
	}
	elsif ($type =~ /(\w+)/i) {
		# text that will be used as the type base when autogenerating
		$type_base = $1;
	}
}

# define source base or index
my ($source_index, $source_base);
if (defined $source) {
	if ($source =~ /^(\d+)$/) {
		# looks like an index was provided
		$source_index = $1;
	}
	elsif ($source =~ /(\w+)/i) {
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
my $Input = Bio::ToolBox::Data::Stream->new(in => $infile) or
	die "Unable to open file '$infile'!\n";

### Determine indices
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
	
	# request source text
	unless (defined $source) {
		# this is a special input, can't use the ask_user_for_index sub
		# accepts either index or text string
		print " Enter the text string or column index for the GFF source (Suggested)  ";
		my $in = <STDIN>;
		if ($in =~ /^(\d+)$/) {
			$source_index = $1;
		}
		elsif ($in =~ /(\w+)/) {
			$source_base = $1;
		}
	}
	
	# request type text
	unless (defined $type) {
		# this is a special input, can't use the ask_user_for_index sub
		# accepts either index or text string
		print " Enter the text string or column index for the GFF type (Suggested)  ";
		my $in = <STDIN>;
		if ($in =~ /^(\d+)$/) {
			$type_index = $1;
		}
		elsif ($in =~ /(\w+)/) {
			$type_base = $1;
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
	
	# request ID index
	unless (defined $id_index) {
		my $suggestion = $Input->id_column;
		$id_index = ask_user_for_index($Input, 
			" Enter the index for the feature unique ID column [$suggestion]  ");
		$id_index = defined $id_index ? $id_index : $suggestion;
	}
	
	# request tags
	unless (defined $tag) {
		@tag_indices = ask_user_for_index($Input, 
			" Enter zero or more column indices for GFF group tags  ");
	}
}
else {
	# or else the indices need to be automatically identified
	unless (
		$Input->feature_type eq 'coordinate' or 
		(defined $chr_index and defined $start_index and defined $stop_index) or 
		($Input->feature_type eq 'named' and $Input->database) 
	) {
		die "Not enough information has been provided to convert to GFF file.\n" . 
			"Coordinate column names must be recognizable or specified. Use --help\n";
	}
}



### Open output stream
unless ($outfile) {
	$outfile = $Input->path . $Input->basename;
}
my $Output = Bio::ToolBox::Data::Stream->new(
	out     => $outfile,
	gff     => 3,        # default, only one
	gz      => $gz
) or die " unable to create output file $outfile!";



### Convert the input stream
# check some things first
my $do_feature = $Input->feature_type eq 'named' ? 1 : 0; # get features from db?
if (defined $start_index and substr($Input->name($start_index), -1) eq '0') {
	# start column name suggests it is 0-based
	$interbase = 1;
}
if ($unique and not (defined $name_index or $name_base)) {
	die " must provide a name index or name base to make unique feature names!\n";
}
my $unique_name_counter = {}; # hash for making unique feature names
my $unique_id_counter   = {}; # same for IDs
my $count = 0; # the number of lines processed
while (my $row = $Input->next_row) {
	
	# get the feature from the db if necessary
	my $f = $row->feature if $do_feature;
	
	# build the arguments
	# retrieve information from row object if indices were provided
	my @args;
	if (defined $chr_index) {
		push @args, 'chromo', $row->value($chr_index);
	}
	if (defined $start_index) {
		my $s = $row->value($start_index);
		$s += 1 if $interbase;
		push @args, 'start', $s;
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
		my $name = $unique ? 
			generate_unique_name($row->value($name_index), $unique_name_counter) : 
			$row->value($name_index);
		push @args, 'name', $name;
	} elsif (defined $name_base) {
		push @args, 'name', sprintf("%s_%07d", $name_base, $count);
	}
	if (defined $id_index) {
		my $id = $row->value($id_index);
		push @args, 'id', generate_unique_name($id, $unique_id_counter);
	}
	if (defined $type_index) {
		push @args, 'type', $row->value($type_index);
	} elsif (defined $type_base) {
		push @args, 'type', $type_base;
	}
	if (defined $source_index) {
		push @args, 'source', $row->value($source_index);
	} elsif (defined $source_base) {
		push @args, 'source', $source_base;
	}
	if (@tag_indices) {
		push @args, 'attributes', \@tag_indices;
	}
			
	# write
	$Output->add_row( $row->gff_string(@args) );
	$count++;
}




### Finish
$Input->close_fh;
$Output->close_fh;
printf " Converted %s lines of input data to GFF file '%s'\n", 
	format_with_commas($count), $Output->filename;


sub generate_unique_name {
	my ($name, $counter) = @_;
	my $new_name;
			
	# check uniqueness
	if (exists $counter->{$name} ) {
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

A script to convert a generic data file to GFF format.

=head1 SYNOPSIS

data2gff.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --ask
  --chr <column_index>
  --start <column_index>
  --stop | --end <column_index>
  --score <column_index>
  --strand <column_index>
  --name <text | column_index>
  --id <column_index>
  --tags <column_index,column_index,...>
  --source <text | column_index>
  --type <text | column_index>
  --unique
  --zero
  --out <filename> 
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify an input file containing either a list of database features or 
genomic coordinates for which to convert to GFF format. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. Files may 
be gzipped compressed.

=item --ask

Indicate that the program should interactively ask for column indices or
text strings for the GFF attributes, including coordinates, source, type, 
etc. It will present a list of the column names to choose from. Enter 
nothing for non-relevant columns or to accept default values.

=item --chr <column_index>

The index of the dataset in the data table to be used 
as the chromosome or sequence ID column in the gff data.

=item --start <column_index>

The index of the dataset in the data table to be used 
as the start position column in the gff data.

=item --stop <column_index>
=item --end <column_index>

The index of the dataset in the data table to be used 
as the stop or end position column in the gff data.

=item --score <column_index>

The index of the dataset in the data table to be used 
as the score column in the gff data.

=item --name <text | column_index>

Enter either the text that will be shared name among 
all the features, or the index of the dataset in the data 
table to be used as the name of each gff feature. This 
information will be used in the 'group' column.

=item --id <column_index>

The index of the dataset in the data table to be used
as the unique ID of each gff feature. This information
will be used in the 'group' column of GFF v3 files 
only. The default is to automatically generate a 
unique identifier.

=item --strand <column_index>

The index of the dataset in the data table to be used
for strand information. Accepted values might include
any of the following "+, -, 1, -1, 0, .".

=item --tags <column_indices>

Provide a comma delimited list of column indices that contain 
values to be included as group tags in the GFF features. The 
key will be the column name.

=item --source <text | column_index>

Enter either a text string or a column index representing the 
GFF source that should be used for the features. The default is 
'data'.

=item --type <text | column_index>

Enter either a text string or a column index representing the 
GFF 'type' or 'method' that should be used for the features. If 
not defined, it will use the column name for either 
the 'score' or 'name' column, if defined. As a last resort, it 
will use the most creative method of 'Experiment'.

=item --unique

Indicate whether the feature names should be made unique. A count 
number is appended to the name of subsequent features to make them 
unique. Using a base text string for the name will automatically 
generate unique names.

=item --zero

Input file is in interbase or 0-based coordinates. This should be 
automatically detected for most known file formats, e.g. BED.

=item --out <filename>

Optionally specify the name of of the output file. The default is to use 
the input file base name. The '.gff' extension is automatically
added if required.

=item --gz

Indicate whether the output file should (not) be compressed with gzip.

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
