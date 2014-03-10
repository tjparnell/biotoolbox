#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::data_helper qw(
	find_column_index
	format_with_commas
);
use Bio::ToolBox::file_helper qw(
	open_tim_data_file
	write_tim_data_file
	open_to_write_fh
	convert_genome_data_2_gff_data
);
my $VERSION = '1.15';

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
	$midpoint,
	$format,
	$zero_based,
	$unique,
	$ask,
	$version,
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
	'midpoint!' => \$midpoint, # boolean to use the midpoint
	'format=i'  => \$format, # format output to indicated number of places
	'zero!'     => \$zero_based, # source is 0-based numbering, convert
	'unique!'   => \$unique, # make the names unique
	'ask'       => \$ask, # request help in assigning indices
	'version=i' => \$version, # the gff version
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
unless ($version) {
	$version = 3;
}
unless (defined $gz) {
	$gz = 0;
}


### Load file
print " Opening data file '$infile'...\n";
my ($in_fh, $metadata_ref) = open_tim_data_file($infile);
unless ($in_fh) {
	die "Unable to open data table!\n";
}


### Determine indices

# Ask user interactively
if ($ask) {
	# the user has specified that we should ask for specific indices
	
	# print the column names
	print " These are the column names in the datafile\n";
	for (my $i = 0; $i < $metadata_ref->{'number_columns'}; $i++) {
		print "   $i\t", $metadata_ref->{$i}{'name'}, "\n";
	}
	print " Note that not all options are not required\n";
	
	# request chromosome index
	unless (defined $chr_index) {
		my $suggestion = find_column_index($metadata_ref, '^chr|seq|refseq');
		print " Enter the index for the chromosome column (Required) [$suggestion]  ";
		my $in = <STDIN>;
		chomp $in;
		if ($in =~ /^\d+$/ and exists $metadata_ref->{$in}) {
			$chr_index = $in;
		}
		elsif (not $in and defined $suggestion) {
			$chr_index = $suggestion;
		}
		else {
			die " No identifiable chromosome column index!\n";
		}
	}
	
	# request start index
	unless (defined $start_index) {
		my $suggestion = find_column_index($metadata_ref, 'start');
		print " Enter the index for the start column (Required) [$suggestion]  ";
		my $in = <STDIN>;
		chomp $in;
		if ($in =~ /^\d+$/ and exists $metadata_ref->{$in}) {
			$start_index = $in;
		}
		elsif (not $in and defined $suggestion) {
			$start_index = $suggestion;
		}
		else {
			die " No identifiable start position column index!\n";
		}
	}
	
	# request stop index
	unless (defined $stop_index) {
		my $suggestion = find_column_index($metadata_ref, 'stop|end');
		print " Enter the index for the stop or end column (Required) [$suggestion]  ";
		my $in = <STDIN>;
		chomp $in;
		if ($in =~ /^\d+$/ and exists $metadata_ref->{$in}) {
			$stop_index = $in;
		}
		elsif (not $in and defined $suggestion) {
			$stop_index = $suggestion;
		}
		else {
			die " No identifiable stop position column index!\n";
		}
	}
	
	# request source text
	unless (defined $source) {
		print " Enter the text string or column index for the GFF source (Suggested)  ";
		$source = <STDIN>;
		chomp $source;
	}
	
	# request type text
	unless (defined $type) {
		print " Enter the text string or column index for the GFF type (Suggested)  ";
		$type = <STDIN>;
		chomp $type;
	}
	
	# request score index
	unless (defined $score_index) {
		my $suggestion = find_column_index($metadata_ref, '^score$');
		print " Enter the index for the feature score column [$suggestion]  ";
		my $in = <STDIN>;
		chomp $in;
		if ($in =~ /^\d+$/ and exists $metadata_ref->{$in}) {
			$score_index = $in;
		}
		elsif (not $in and defined $suggestion) {
			$score_index = $suggestion;
		}
		elsif ($in) {
			print " unrecognized index, skipping\n";
		}
	}
	
	# request strand index
	unless (defined $strand_index) {
		my $suggestion = find_column_index($metadata_ref, 'strand');
		print " Enter the index for the feature strand column [$suggestion]  ";
		my $in = <STDIN>;
		chomp $in;
		if ($in =~ /^\d+$/ and exists $metadata_ref->{exists}) {
			$strand_index = $in;
		}
		elsif (not $in and defined $suggestion) {
			$strand_index = $suggestion;
		}
		elsif ($in) {
			print " unrecognized index, skipping\n";
		}
	}
	
	# request name index or text
	unless (defined $name) {
		my $suggestion = find_column_index($metadata_ref, '^name|id');
		print " Enter the index for the feature name column or the text (Suggested) [$suggestion]  ";
		$name = <STDIN>;
		chomp $name;
		if (not $name and defined $suggestion) {
			$name = $suggestion;
		}
	}
	
	# request ID index
	unless (defined $id_index) {
		print " Enter the index for the feature unique ID column  ";
		my $in = <STDIN>;
		chomp $in;
		if ($in =~ /^\d+$/ and exists $metadata_ref->{$in}) {
			$id_index = $in;
		}
		elsif ($in) {
			print " unrecognized index, skipping\n";
		}
	}
	
	# request tags
	unless (defined $tag) {
		print " Enter a comma-delimited list of column indices for GFF group tags  ";
		$tag = <STDIN>;
		chomp $tag;
		$tag =~ s/\s+//g; # remove whitespace
	}
}

# otherwise attempt to identify indices automatically
else {
	unless (defined $chr_index) {
		$chr_index = find_column_index($metadata_ref, '^chr|seq|refseq');
	}
	unless (defined $start_index) {
		$start_index = find_column_index($metadata_ref, '^start');
	}
	unless (defined $stop_index) {
		$stop_index = find_column_index($metadata_ref, '^stop|end');
	}
	unless (defined $strand_index) {
		$strand_index = find_column_index($metadata_ref, '^strand');
	}
	unless (defined $name or defined $name_index) {
		$name_index = find_column_index($metadata_ref, '^name|ID');
		unless (defined $name_index) {
			$name_index = find_column_index($metadata_ref, 'name|ID$');
		}
	}
	unless (defined $score_index) {
		$score_index = find_column_index($metadata_ref, '^score$');
	}
	unless (defined $type) {
		$type = $metadata_ref->{'feature'} || 'region';
	}
	unless (defined $source) {
		$source = $metadata_ref->{'basename'};
	}
}

if ($metadata_ref->{'bed'}) {
	# bedfiles are 0-based
	$zero_based = 1;
}


# Convert tag text to list of indices
my @tag_indices;
if ($tag) {
	foreach (split /,\s*/, $tag) {
		if (/(\d+)/) {
			push @tag_indices, $1;
		}
	}
}

# Determine the name index or text
if (
	defined $name and 
	$name =~/^\d+$/ and 
	exists $metadata_ref->{$name}
) {
	# looks like a number, so assume it is the index
	$name_index = $name;
}
if ($unique and (not defined $name and not defined $name_index)) {
	die " unable to assign unique feature names without a name index or text!\n";
}
if ($unique and not defined $name_index) {
	# a uniqe name is requested but we don't have an index, yet
	# so make one
	$name_index = $metadata_ref->{'number_columns'};
	$metadata_ref->{$name_index} = {
		'name'   => 'Name',
		'index'  => $name_index,
	};
	push @{ $metadata_ref->{'column_names'} }, 'Name';
	$metadata_ref->{'number_columns'} += 1;
}

# Generate the output file name
unless ($outfile) {
	# re-use the input file name
	$outfile = $metadata_ref->{'path'} . $metadata_ref->{'basename'};
}



### Print the user or auto-generated options
print " converting to GFF using\n";
print "  - '", $metadata_ref->{$chr_index}{name}, "' column for chromosome\n" 
	if defined $chr_index;
print "  - '", $metadata_ref->{$start_index}{name}, "' column for start\n" 
	if defined $start_index;
print "  - '", $metadata_ref->{$stop_index}{name}, "' column for stop\n" 
	if defined $stop_index;
print "  - '", $metadata_ref->{$strand_index}{name}, "' column for strand\n" 
	if defined $strand_index;
print "  - '", $metadata_ref->{$score_index}{name}, "' column for score\n" 
	if defined $score_index;
if (defined $name_index) {
	print "  - '", $metadata_ref->{$name_index}{name}, "' column for name\n" 
}
elsif (defined $name) {
	print "  - '$name' for name\n" 
}
print "  - '", $metadata_ref->{$id_index}{name}, "' column for ID\n" 
	if defined $id_index;
if ($type =~ /^\d+$/ and exists $metadata_ref->{$type}) {
	# type looks like a column index
	print "  - '", $metadata_ref->{$type}{name}, "' column for GFF type\n";
}
else {
	# type must be a text string
	print "  - '$type' for type\n" if defined $type;
}
if ($source =~ /^\d+$/ and exists $metadata_ref->{$source}) {
	# type looks like a column index
	print "  - '", $metadata_ref->{$source}{name}, "' column for GFF source\n";
}
else {
	# type must be a text string
	print "  - '$source' for GFF source\n" if defined $source;
}
print "  - '", join(", ", map { $metadata_ref->{$_}{name} } @tag_indices ), 
	"' columns for group tags\n" if $tag;




### Generate a temporary gff data structure based on the input metadata
# doing this manually rather than calling Bio::ToolBox::data_helper to make a new 
# structure, because we're basing this off the input file metadata, 
# and we populate, empty, and regenerate numerous times as we walk 
# through the file. Ugh, it's complicated, and very old un-optimized, crazy code
my %output_data;
if (exists $metadata_ref->{'program'}) {
	# use originating data file program
	$output_data{'program'} = $metadata_ref->{'program'};
}
else {
	# use current program name
	$output_data{'program'} = $0;
}
if (exists $metadata_ref->{'feature'}) {
	$output_data{'feature'} = $metadata_ref->{'feature'};
}
if (exists $metadata_ref->{'db'}) {
	$output_data{'db'} = $metadata_ref->{'db'};
}
$output_data{'other'} = [ @{ $metadata_ref->{'other'} } ];
$output_data{'data_table'} = [];
regenerate_output_data_hash();



### Global variables
# set output control variables
my $out_fh; # the output file handle
my $count = 0; # the number of lines processed before writing output
my $total_count = 0;

# set unique name counter
my %unique_name_counter;




### Parse input file into GFF
# To avoid exorbitant memory requirements for ginormous files, 
# we will only convert 20000 lines at time. 
while (my $line = $in_fh->getline) {
	
	# add line to the data table
	chomp $line;
	push @{ $output_data{'data_table'} }, [ split /\t/, $line ];
	$output_data{'last_row'} += 1;
	
	# increment counter
	$count++;
	
	# temporarily write output
	if ($count == 20000) {
				
		# generate unique names if requested
		if ($unique) {
			generate_unique_names();
		}
		
		# progressively write out the converted gff data
		write_gff_data();
		
		# regenerate the output data hash
		regenerate_output_data_hash();
		
		# reset count
		$total_count += $count;
		$count = 0;
	}
}





### Write final output
if ($unique) {
	generate_unique_names();
}
write_gff_data();
$total_count += $count;



### Finish
$in_fh->close;
$out_fh->close;
printf " Converted %s lines of input data to GFF file '$outfile'\n", 
	format_with_commas($total_count);
if ($unique) {
	printf " There were %s original unique names\n", 
		format_with_commas(scalar keys %unique_name_counter);
}
# That's it!




sub regenerate_output_data_hash {
	# a subroutine to prepare a temporary output data hash for gff conversion
	
	# reset data table specifics
	$output_data{'number_columns'} = $metadata_ref->{'number_columns'};
	$output_data{'gff'} = 0;
	$output_data{'last_row'} = 0;
	
	# delete the column metadata hashes
	for (my $i = 0; $i < 9; $i++) {
		delete $output_data{$i};
	}
	
	# restore the original column metadata
	for (my $i = 0; $i < $metadata_ref->{'number_columns'}; $i++) {
		$output_data{$i} = { %{ $metadata_ref->{$i} } };
	}
	
	# empty the table
	$output_data{'data_table'} = [ q() ]; # clear the array
	
	# add the original column names
	$output_data{'data_table'}->[0] = [ 
		@{ $metadata_ref->{'column_names'} } 
	];
}


sub write_gff_data {
	# a subroutine to progressively write out the converted gff data
	
	# convert to gff
	my @arguments = (
		'data'     => \%output_data,
		'chromo'   => $chr_index,
		'start'    => $start_index,
		'stop'     => $stop_index,
		'score'    => $score_index,
		'strand'   => $strand_index,
		'source'   => $source, # will interpret as either text or index
		'type'     => $type, # will interpret as either text or index
		'midpoint' => $midpoint,
		'version'  => $version,
		'tags'     => [ @tag_indices ],
		'id'       => $id_index,
		'zero'     => $zero_based,
	);
	if ($unique) {
		# we've generated a new name index
		push @arguments, 'name' => $name_index;
	}
	elsif (defined $name_index) {
		# supplied name index
		push @arguments, 'name' => $name_index;
	}
	else {
		# text name
		push @arguments, 'name' => $name;
	}
	convert_genome_data_2_gff_data( @arguments ) or 
		die " Unable to convert to GFF format!\n";
	
	# format the numbers
	if (defined $format) {
		for (my $row = 1; $row <= $output_data{'last_row'}; $row++) {
			# walk through each row in the data table
			# format the score value to the indicated number of spaces
			if ($format == 0) {
				# no decimal places
				$output_data{'data_table'}->[$row][5] = 
					sprintf( "%.0f", $output_data{'data_table'}->[$row][5] );
			}
			elsif ($format == 1) {
				# 1 decimal place
				$output_data{'data_table'}->[$row][5] = 
					sprintf( "%.1f", $output_data{'data_table'}->[$row][5] );
			}
			elsif ($format == 2) {
				# 2 decimal places
				$output_data{'data_table'}->[$row][5] = 
					sprintf( "%.2f", $output_data{'data_table'}->[$row][5] );
			}
			elsif ($format == 3) {
				# 3 decimal places
				$output_data{'data_table'}->[$row][5] = 
					sprintf( "%.3f", $output_data{'data_table'}->[$row][5] );
			}
		}
		# update metadata
		$output_data{5}{'formatted'} = $format;
	}
	
	# check for file handle
	if ($out_fh) {
		# the output file has been opened and partially written
		# we now only need to write the data table portion and not the 
		# metadata
		
		for my $row (1..$output_data{'last_row'}) {
			print {$out_fh} join(
					"\t", @{ $output_data{'data_table'}->[$row] }
				), "\n";
				 
		}
	}
	else {
		# we will need to open the output file to write if it's not 
		# opened yet
		
		# rather than generating new code for writing the gff file,
		# we will simply use the write_tim_data_file sub
		$outfile = write_tim_data_file(
			'data'      => \%output_data,
			'filename'  => $outfile,
			'gz'        => $gz,
		);
		
		# but now we will have to reopen the file for appended writing
		$out_fh = open_to_write_fh($outfile, $gz, 1);
	}

}


sub generate_unique_names {
	
	for (my $row = 1; $row <= $output_data{'last_row'}; $row++) {
		
		# name scalars
		my $new_name;
		my $current; # current name
		
		# determine current name
		if (defined $output_data{'data_table'}->[$row][$name_index]) {
			$current = $output_data{'data_table'}->[$row][$name_index];
		}
		else {
			$current = $name;
		}
			
		# check uniqueness
		if (exists $unique_name_counter{$current} ) {
			# we've encountered this name before
			# generate a unique name by appending the count number
			$unique_name_counter{ $current } += 1;
			$new_name = $current . '.' . 
				$unique_name_counter{ $current };
		}
		else {
			# first time for this name
			# record in the hash
			$new_name = $current;
			$unique_name_counter{$current} = 0;
		}
		
		# assign the new name
		$output_data{'data_table'}->[$row][$name_index] = $new_name;
	}

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
  --source <text>
  --type <text | column_index>
  --zero
  --format [0,1,2,3]
  --midpoint
  --unique
  --out <filename> 
  --version [2,3]
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. Genome 
coordinates are required. The first row should be column headers. Text 
files generated by other B<BioToolBox> scripts are acceptable. Files may 
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

=item --name <column_index>

Enter either the text that will be shared name among 
all the features, or the index of the dataset in the data 
table to be used as the name of each gff feature. This 
information will be used in the 'group' column.

=item --id <column_index>

The index of the dataset in the data table to be used
as the unique ID of each gff feature. This information
will be used in the 'group' column of GFF v.3 files 
only. The default is to automatically generate a 
unique identifier.

=item --strand <column_index>

The index of the dataset in the data table to be used
for strand information. Accepted values might include
any of the following "f(orward), r(everse), w(atson),
c(rick), +, -, 1, -1, 0, .".

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

=item --zero

Indicate whether the source data is in interbase or 0-based 
coordinates, as is used with UCSC source data or USeq data 
packages. The coordinates will then be converted to 1-based 
coordinates, consistent with the rest of bioperl conventions.
The default is false (will not convert).

=item --format [0,1,2,3]

Indicate the number of decimal places the score value should
be formatted. Acceptable values include 0, 1, 2, or 3 places.
Anything else is ignored.

=item --midpoint

A boolean (1 or 0) value to indicate whether the 
midpoint between the actual 'start' and 'stop' values
should be used instead of the actual values. Default 
is false.

=item --unique

Indicate whether the feature names should be made unique. A count 
number is appended to the name of subsequent features to make them 
unique. This should only be applied to genomic features, and not to 
genomic data values (microarray data, sequencing data, etc). The 
default behavior is false (not unique).

=item --out <filename>

Optionally specify the name of of the output file. The default is to use 
the assigned type value. The '.gff' extension is automatically
added if required.

=item --version [2,3]

Specify the GFF version. The default is version 3.

=item --gz

Indicate whether the output file should (not) be compressed with gzip.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will convert a data file into a GFF formatted text file. 
Only simple conversions are performed, where each data line is converted 
to a single feature. Complex features with parent-child relationships (such 
as genes) should be converted with something more advanced.

The input file should have chromosomal coordinates, i.e. chromosome, 
start, and (optionally) stop or end coordinates. They may be specified 
upon execution or identified automatically. If they are not found, the 
GFF conversion will fail. 

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
