#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::data_helper qw(find_column_index);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	open_tim_data_file
	write_tim_data_file
);
my $VERSION = '1.14';


print "\n This script will map oligo data to the genome and generate a GFF file\n";


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
	$oligo_file,
	$data_file,
	$column,
	$type,
	$source,
	$name,
	$strand,
	$midpoint,
	$places,
	$outfile,
	$gz,
	$help,
	$print_version,
);


# Command line options
GetOptions( 
	'oligo=s'     => \$oligo_file, # the oligo feature file name
	'data=s'      => \$data_file, # the oligo data file
	'index|col=i' => \$column, # the index of the data column
	'type=s'      => \$type, # the new gff type
	'source=s'    => \$source, # the new source
	'name=s'      => \$name, # the name for the feature
	'strand!'     => \$strand, # respect or ignore strand information
	'mid!'        => \$midpoint, # convert positions to midpoint
	'places=i'    => \$places, # number of digits to format the data number
	'out=s'       => \$outfile, # the output file name
	'gz!'         => \$gz, # compress output file with gzip
	'help'        => \$help, # request help
	'version'     => \$print_version, # print the version
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
	print " Biotoolbox script map_oligo_data2gff.pl, version $VERSION\n\n";
	exit;
}


### Check for required values
unless ($oligo_file and $data_file) {
	die " Both source oligo data and oligo feature files must be specified!\n";
}


### Load Oligo data file
my ($oligo_data_ref, $oligo_data_metadata) = load_microarray_data();



### Define default data
unless ($name) {
	$name = $oligo_data_metadata->{$column}{'name'};
}
unless ($type) {
	$type = $name;
}
unless ($source) {
	$source = 'lab';
}
unless (defined $gz) {
	$gz = 0;
}
my $formatter;
if (defined $places) {
	# determine how many places to format the number
	if ($places == 0) {
		$formatter = "%.0f";
	}
	elsif ($places == 1) {
		$formatter = "%.1f";
	}
	elsif ($places == 2) {
		$formatter = "%.2f";
	}
	elsif ($places == 3) {
		$formatter = "%.3f";
	}
}
unless ($outfile) {
	$outfile = $name;
}
unless (defined $strand) {
	# default is to discard strand information
	$strand = 0;
}
unless (defined $midpoint) {
	# default is to convert position information to a single midpoint
	$midpoint = 1;
}



### Load the oligo GFF data
my $oligo_feature_ref = load_tim_data_file($oligo_file) or 
	die " unable to load oligo feature file!\n";

unless ($oligo_feature_ref->{gff}) {
	die " oligo feature file does not appear to be GFF format!\n";
}



### Start the conversions

	# we will actually be modifying the loaded oligo GFF data with the 
	# microarray data, placing it in the score column and updating the 
	# group column

# check metadata
if (scalar keys %{ $oligo_data_metadata->{$column} } > 2) {
	# there appears to be more than the basic data in here
	
	# replace the metadata
	$oligo_feature_ref->{5} = $oligo_data_metadata->{$column}; 
		# replace the hash data
		# and update the info
	$oligo_feature_ref->{5}{name} = 'Score';
	$oligo_feature_ref->{5}{index} = 5;
}
$oligo_feature_ref->{5}{'source_data_file'} = $data_file;
$oligo_feature_ref->{5}{'formatted_places'} = $places if defined $places;



# Walk through the data
print " Adding microarray data values to oligo features....\n";

my $table = $oligo_feature_ref->{data_table}; # shortcut ref
my @no_data; # an array of the oligos with no data, to be deleted
for (my $row = 1; $row <= $oligo_feature_ref->{last_row}; $row++) {
	
	# first need to identify the oligo name
	my $oligo_name;
	if ($oligo_feature_ref->{gff} == 3) {
		# gff version 3 file
		my %stuff = split /\s?[;=]\s?/, $table->[$row][8];
			# attempt to split both key=value pairs and multiple pairs
			# simultaneously
		#print "  row $row has " . join(", ", keys %stuff) . "\n";
		$oligo_name = $stuff{'Name'}; 
			# we will use the name instead of the ID since the ID may have 
			# unique identifier numbers appended to it for multicopy oligos
	}
	else {
		# assume gff version 2????
		my @stuff = split /\s*;\s*/, $table->[$row][8];
		my ($class, $id) = split / /, $stuff[0];
		$id =~ s/"//g; # strip any quotation marks
		$oligo_name = $id;
	}
	
	if (exists $oligo_data_ref->{$oligo_name}) {
		# oligo does exist!
		
		# update the score
		if ($formatter) {
			# format the score value and put into the gff table
			$table->[$row][5] = sprintf $formatter, 
				$oligo_data_ref->{$oligo_name};
		}
		else {
			# no formatting, just put into gff table
			$table->[$row][5] = $oligo_data_ref->{$oligo_name};
		}
		
		# update type
		$table->[$row][2] = $type;
		
		# update source
		$table->[$row][1] = $source;
		
		# update position
		if ($midpoint) {
			# convert position to midpoint
			my $position = int( ($table->[$row][3] + $table->[$row][4]) / 2);
			$table->[$row][3] = $position;
			$table->[$row][4] = $position;
		}
		
		# update strand
		unless ($strand) {
			# if true, we'll keep the strand information
			# otherwise we toss strand information
			$table->[$row][6] = '.';
		}
		
		# update group
		$table->[$row][8] = "ID=$name." . $row . ";Name=$name";
			# make the ID unique by appending the row number to the name
		
	}
	else {
		# oligo data does not exist!
		# remember this for future deletion
		push @no_data, $row;
	}
		
}



### Delete the features with no data
if (@no_data) {
	print " ", scalar @no_data, " oligo features did not match a data value!\n";
	
	# delete the rows
	while (@no_data) {
		# proceed from the bottom up
		my $row = pop @no_data;
		splice @{ $table}, $row, 1; # delete the row
	}
	
	$oligo_feature_ref->{'last_row'} = scalar @{ $table } - 1;
}



### Write out the new gff

# update metadata
$oligo_feature_ref->{'gff'} = "3";
$oligo_feature_ref->{'basename'} = undef;
$oligo_feature_ref->{'filename'} = undef;
$oligo_feature_ref->{'extension'} = undef;

my $success = write_tim_data_file(
	'data'       => $oligo_feature_ref,
	'filename'   => $outfile,
	'gz'         => $gz
);
if ($success) {
	print " Wrote file '$success'\n";
}
else {
	print " Unable to write file!\n";
}



#######  Subroutines ##############

sub load_microarray_data {
	
	# open the data file
	my ($data_fh, $data_ref) = open_tim_data_file($data_file) or 
		die " Unable to open oligo data file '$data_file'!\n";
	
	# Determine the column of microarray data
	unless (defined $column) {
		# print the headers
		print "\n These are the columns in the data file.\n";
		for (my $i = 0; $i < $data_ref->{'number_columns'}; $i++) {
			print "   $i\t$data_ref->{$i}{name}\n";
		}
		
		# process the answer
		print " Enter the number of the column with the data   ";
		my $answer = <STDIN>;
		chomp $answer;
		
		# check answer and return
		if (exists $data_ref->{$answer}) {
			# answer appears to be a column index
			$column = $answer;
		}
		else {
			die " Invalid response!\n";
		}
	}
	
	# identify the probe ID column index
	my $probe_i = find_column_index($data_ref, "probe|oligo|id");
	
	# Load the microarray values data into a hash
	my %hash;
	while (my $line = $data_fh->getline) {
		
		# process line
		chomp $line;
		my @data = split /\t/, $line;
		
		# check that the probe is unique
		if (exists $hash{ $data[$probe_i] } ) {
			warn " Probe '" . $data[$probe_i] . 
				"' exists more than once! Using first value only\n";
		}
		else {
			$hash{ $data[$probe_i] } = $data[$column];
		}
	}
	$data_fh->close;
	print " loaded " . scalar(keys %hash) . " microarray probe values from '$data_file'\n";
	
	# return
	return (\%hash, $data_ref);
}







__END__

=head1 NAME

map_oligo_data2gff.pl

A script to assign processed microarray data to genomic coordinates.

=head1 SYNOPSIS

map_oligo_data2gff.pl --oligo <oligo_file.gff> --data <oligo_data.txt> [--options]
  
  Options:
  --oligo <oligo_file.gff>
  --data <oligo_data.txt>
  --index <column_index>
  --name <text>
  --type <text>
  --source <text>
  --strand
  --(no)mid
  --places [0,1,2,3]
  --out <filename> 
  --gz
  --version
  --help
  
=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --oligo <oligo_file.gff>

Specify the file name of the oligo probe GFF file. This file should identify 
all positions that the probes align to in the genome. It may be either v.2 or 
v.3 GFF file, but the oligo probe name or ID must be present in the group 
field.

=item --data <oligo_data.tx>

Specify the file name of a data file. It must be a tab-delimited text file,
preferably in the tim data format as described in Bio::ToolBox::file_helper, 
although any format should work. The file may be compressed with gzip. The 
first column MUST be the oligo or probe unique name or ID.

=item --index <column_index>

Specify the data column index in the data file that will be used in the 
final GFF score column. By default, it lists the column names for the user 
to interactively select.

=item --col <column_index>

Alias to --index.

=item --name <text>

Specify the name of the dataset to be used in the output GFF file. By default 
it uses the column name from the data file.

=item --type <text>

Specify the output GFF type or method field value. By default, it uses the name.

=item --source <text>

Specify the output GFF source field value. By default it is 'lab'.

=item --strand

Indicate whether the original strand information should be kept from the 
original oligo GFF file. The default is false, as most ChIP data is inherently 
not stranded.

=item --(no)mid

Indicate whether the original position information should (not) be converted 
to the midpoint position. The GBrowse xyplot data works best when data is 
present at single points (start = end) rather than regions. Default is true.

=item --places [0,1,2,3]

Indicate the number of decimal places the score value should be formatted. The 
default is no formatting.

=item --out <filename> 

Specify the output filename. The default is to use the name value.

=item --gz

Indicate whether the output file should (not) be compressed with gzip.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This script assigns data values for microarray oligo probes to positions in the 
genome. It essentially merges the information from a data file, consisting of 
unique oligo probe names and data values, with a GFF file referencing the 
genomic positions of each microarray oligo probe. 

The score value may be formatted, and the GFF type, source, name, and strand 
may be set to new values. It will write a GFF v.3 file, regardless of the 
source oligo GFF file version.

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
