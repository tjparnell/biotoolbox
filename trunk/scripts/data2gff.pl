#!/usr/bin/perl

# A script to convert a generic data file into a gff file
# this presumes it has chromosomal coordinates to convert

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	find_column_index
);
use tim_file_helper qw(
	open_tim_data_file
	write_tim_data_file
	open_to_write_fh
	convert_genome_data_2_gff_data
);

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
	$help
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

# automatically identify sgr format
if ($metadata_ref->{'extension'} =~ /sgr/) {
	$chr_index = 0 unless defined $chr_index;
	$start_index = 1 unless defined $start_index;
	$stop_index = 1 unless defined $stop_index; # same as start
	$score_index = 2 unless defined $score_index;
}

# automatically identify bed format
elsif ($metadata_ref->{'extension'} =~ /bed/) {
	$chr_index = 0 unless defined $chr_index;
	$start_index = 1 unless defined $start_index;
	$stop_index = 2 unless defined $stop_index; # same as start
	
	if (!defined $name_index and $metadata_ref->{'number_columns'} >= 4) {
		$name_index = 3;
	}
	if (!defined $score_index and $metadata_ref->{'number_columns'} >= 5) {
		$score_index = 4;
	}
	if (!defined $strand_index and $metadata_ref->{'number_columns'} >= 6) {
		$strand_index = 5;
	}
}

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
		print " Enter the index for the chromosome column [$suggestion]  ";
		my $in = <STDIN>;
		if ($in =~ /(\d+)/) {
			$chr_index = $1;
		}
		elsif (defined $suggestion) {
			$chr_index = $suggestion;
		}
		else {
			die " No identifiable chromosome column index!\n";
		}
	}
	
	# request start index
	unless (defined $start_index) {
		my $suggestion = find_column_index($metadata_ref, 'start');
		print " Enter the index for the start column [$suggestion]  ";
		my $in = <STDIN>;
		if ($in =~ /(\d+)/) {
			$start_index = $1;
		}
		elsif (defined $suggestion) {
			$start_index = $suggestion;
		}
		else {
			die " No identifiable start position column index!\n";
		}
	}
	
	# request stop index
	unless (defined $stop_index) {
		my $suggestion = find_column_index($metadata_ref, 'stop|end');
		print " Enter the index for the stop or end column [$suggestion]  ";
		my $in = <STDIN>;
		if ($in =~ /(\d+)/) {
			$stop_index = $1;
		}
		elsif (defined $suggestion) {
			$stop_index = $suggestion;
		}
		else {
			die " No identifiable stop position column index!\n";
		}
	}
	
	# request score index
	unless (defined $score_index) {
		print " Enter the index for the feature score column  ";
		my $in = <STDIN>;
		if ($in =~ /(\d+)/) {
			$score_index = $1;
		}
	}
	
	# request name index or text
	unless (defined $name) {
		print " Enter the index for the feature name column or the text  ";
		$name = <STDIN>;
		chomp $name;
	}
	
	# request name index
	unless (defined $id_index) {
		print " Enter the index for the feature ID column  ";
		my $in = <STDIN>;
		if ($in =~ /(\d+)/) {
			$id_index = $1;
		}
	}
	
	# request strand index
	unless (defined $strand_index) {
		print " Enter the index for the feature strand column  ";
		my $in = <STDIN>;
		if ($in =~ /(\d+)/) {
			$strand_index = $1;
		}
	}
	
	# request type text
	unless (defined $type) {
		print " Enter the text string or column index for the GFF type  ";
		$type = <STDIN>;
		chomp $type;
	}
	
	# request source text
	unless (defined $source) {
		print " Enter the text string for the GFF source  ";
		$source = <STDIN>;
		chomp $source;
	}
	
	# request tags
	unless (defined $tag) {
		print " Enter a comma-delimited list of column indices for GFF group tags  ";
		$tag = <STDIN>;
		chomp $tag;
		$tag =~ s/\s+//g; # remove whitespace
	}
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
if (defined $name and $name =~/^\d+$/) {
	$name_index = $name;
}
if ($unique and !defined $name) {
	die " unable to assign unique feature names without a name index or text!\n";
}
if ($unique and !defined $name_index) {
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


### Convert to GFF progressively
# To avoid exorbitant memory requirements for ginormous files, 
# we will only convert 20000 lines at time. 
print " converting to GFF using\n";
print "  - '", $metadata_ref->{$chr_index}{name}, "' for chromosome\n" 
	if defined $chr_index;
print "  - '", $metadata_ref->{$start_index}{name}, "' for start\n" 
	if defined $start_index;
print "  - '", $metadata_ref->{$stop_index}{name}, "' for stop\n" 
	if defined $stop_index;
print "  - '", $metadata_ref->{$strand_index}{name}, "' for strand\n" 
	if defined $strand_index;
print "  - '", $metadata_ref->{$score_index}{name}, "' for score\n" 
	if defined $score_index;
if (defined $name_index) {
	print "  - '", $metadata_ref->{$name_index}{name}, "' for name\n" 
}
elsif (defined $name) {
	print "  - '$name' for name\n" 
}
print "  - '", $metadata_ref->{$id_index}{name}, "' for ID\n" 
	if defined $id_index;
if ($type =~ /^\d+$/ and $type <= $metadata_ref->{'number_columns'}) {
	# type looks like a column index
	print "  - '", $metadata_ref->{$type}{name}, "' for GFF type\n";
}
else {
	# type must be a text string
	print "  - '$type' for type\n" if defined $type;
}
print "  - '$source' for source\n" if defined $source;
print "  - '", join(", ", map { $metadata_ref->{$_}{name} } @tag_indices ), 
	"' for group tags\n" if $tag;


# first, generate a temporary gff data structure based on the input metadata
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

# set output control variables
my $out_fh; # the output file handle
my $count = 0; # the number of lines processed before writing output
my $total_count = 0;

# set unique name counter
my %unique_name_counter;

# parse through the data lines in the input data file
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
print " Converted $total_count lines of input data to GFF file '$outfile'\n";
if ($unique) {
	print " There were ", scalar keys %unique_name_counter, 
		" original unique names\n";
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
	my %arguments = (
		'data'     => \%output_data,
		'chromo'   => $chr_index,
		'start'    => $start_index,
		'stop'     => $stop_index,
		'score'    => $score_index,
		'strand'   => $strand_index,
		'source'   => $source,
		'type'     => $type,
		'midpoint' => $midpoint,
		'version'  => $version,
		'tags'     => [ @tag_indices ],
		'id'       => $id_index,
		'zero'     => $zero_based,
	);
	if ($unique) {
		# we've generated a new name index
		$arguments{'name'} = $name_index;
	}
	elsif (defined $name_index) {
		# supplied name index
		$arguments{'name'} = $name_index;
	}
	else {
		# text name
		$arguments{'name'} = $name;
	}
	convert_genome_data_2_gff_data( \%arguments ) or 
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
		
		# check for filename
		unless ($outfile) {
			# by default, we typically use the method or type name
			# we'll grab it from the first feature, type column
			# we're grabbing it from the data table because the user
			# may not have specified it via $type, whereupon the the
			# gff conversion subroutine will automatically derive one
			$outfile = $output_data{'data_table'}->[1][2];
			# the extension will be changed automatically
		}
		
		# rather than generating new code for writing the gff file,
		# we will simply use the write_tim_data_file sub
		$outfile = write_tim_data_file( {
			'data'      => \%output_data,
			'filename'  => $outfile,
			'gz'        => $gz,
		} );
		
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
  --(no)zero
  --format [0,1,2,3]
  --(no)midpoint
  --(no)unique
  --out <filename> 
  --version [2,3]
  --(no)gz
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a data file. It must be a tab-delimited text file,
preferably in the tim data format as described in 'tim_file_helper.pm', 
although any format should work. The file may be compressed with gzip.

=item --ask

Indicate that the program should interactively ask for the indices for 
feature score, name, and strand. It will present a list of the column 
names to choose from. It will also ask for the source, and type
or method text strings. Enter nothing for non-relevant columns or to 
accept default values.

=item --chr <column_index>

The index of the dataset in the data table to be used 
as the chromosome or sequence ID column in the gff data.

=item --start <column_index>

The index of the dataset in the data table to be used 
as the start position column in the gff data.

=item --start <column_index>
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

=item --source <text>

A scalar value to be used as the text in the 'source' 
column. Default is 'data'.

=item --type <text | column_index>

Enter either a text string to be used as the GFF 'method'
or 'type' column, or a column index representing the 
feature type when one simply won't due. If not defined, 
it will use the name of the dataset used for either the 'score' 
or 'name' column, if defined. As a last resort, it 
will use the most creative method of 'Experiment'.

=item --(no)zero

Indicate whether the source data is in interbase or 0-based 
coordinates, as is used with UCSC source data or USeq data 
packages. The coordinates will then be converted to 1-based 
coordinates, consistent with the rest of bioperl conventions.
The default is false (will not convert).

=item --format [0,1,2,3]

Indicate the number of decimal places the score value should
be formatted. Acceptable values include 0, 1, 2, or 3 places.
Anything else is ignored.

=item --(no)midpoint

A boolean (1 or 0) value to indicate whether the 
midpoint between the actual 'start' and 'stop' values
should be used instead of the actual values. Default 
is false.

=item --(no)unique

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

=item --(no)gz

Indicate whether the output file should (not) be compressed with gzip.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will convert a data file into a GFF formatted text file. 
Only simple conversions are performed, where each data line is converted 
to a single feature. Complex features with parent-child relationships (such 
as genes) should be converted with something else.

The source file should have chromosomal coordinates, i.e. chromosome, 
start, and (optionally) stop or end coordinates. They may be specified 
upon execution or identified automatically. If they are not found, the 
GFF conversion will fail. 

Additional feature information may be specified using appropriate command-line 
arguments, including score values, strand, name, ID, source, type, etc. See 
the L<OPTIONS> for more details. 


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





