#!/usr/bin/perl

# A script to map data around a feature start position

use strict;
use Pod::Usage;
use Getopt::Long;
use Statistics::Lite qw(mean median sum stddevp min max);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	find_column_index
);
use tim_db_helper qw(
	validate_dataset_list
	get_dataset_list
	open_db_connection
	get_new_feature_list
	get_region_dataset_hash
);
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
	write_summary_data
);
#use Data::Dumper;

print "\n This script will map data points relative to a genomic feature\n\n";

### Quick help
unless (@ARGV) { # when no command line options are present
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Get command line options and initialize values

## Initialize values
my (
	$infile, 
	$outfile, 
	$database, 
	$dataset, 
	$feature, 
	$value_type,
	$method,
	$win, 
	$number,
	$position, 
	$strand,
	$set_strand,
	$smooth,
	$sum,
	$log,
	$gz,
	$raw,
	$help
); # command line variables

## Command line options
GetOptions( 
	'out=s'      => \$outfile, # output file name
	'in=s'       => \$infile, # input file name
	'db=s'       => \$database, # database name
	'data=s'     => \$dataset, # dataset name
	'feature=s'  => \$feature, # type of feature
	'value=s'    => \$value_type, # the type of data to collect
	'method=s'   => \$method, # method to combine data
	'window=i'   => \$win, # window size
	'number=i'   => \$number, # number of windows
	'position=s' => \$position, # indicate relative location of the feature
	'strand=s'   => \$strand, # collected stranded data
	'set_strand' => \$set_strand, # enforce an artificial strand
	'smooth!'    => \$smooth, # smooth by interpolation
	'sum!'       => \$sum, # generate average profile
	'raw'        => \$raw, # write raw data
	'log!'       => \$log, # data is in log2 space
	'gz!'        => \$gz, # compress the output file
	'help'       => \$help, # print help
);


# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}






## Check for required values
unless ($database or $infile) {
	die " You must define a database!\n Use --help for more information\n";
}

unless ($outfile) {
	if ($infile) {
		$outfile = $infile;
	}
	else {
		die " You must define an output filename !\n Use --help for more information\n";
	}
}
$outfile =~ s/\.txt$//; # strip extension, we'll add it later

unless ($feature or $infile) {
	die " You must define a feature or use an input file!\n Use --help for more information\n";
}

unless ($win) {
	print " Using default window size of 50 bp\n";
	$win = 50;
}

unless ($number) {
	print " Using default window number of 20 per side\n";
	$number = 20;
}

unless (defined $log) {
	# default is to assume not log2
	$log = 0;
}

if (defined $value_type) {
	# check the region method or type of data value to collect
	unless (
			$value_type eq 'score' or
			$value_type eq 'length' or
			$value_type eq 'count' or
			$value_type eq 'enumerate' # legacy value
	) {
		die " Unknown data value '$value_type'!\n " . 
			"Use --help for more information\n";
	}
	if ($value_type eq 'enumerate') {
		$value_type = 'count'; # make it consistent with tim_db_helper
	}
}
else {
	# default is to take the score
	print " Collecting default data 'score' values\n";
	$value_type = 'score';
}

if (defined $method) {
	# check the requested method
	unless (
			$method eq 'mean' or
			$method eq 'median' or
			$method eq 'sum' or
			$method eq 'min' or
			$method eq 'max' or
			$method eq 'stddev'
	) {
		die " Unknown method '$method'!\n Use --help for more information\n";
	}
}
else {
	# set default method
	if ($value_type eq 'count') {
		print " Using default method of 'sum'\n";
		$method = 'sum';
	}
	else {
		print " Using default method of 'mean'\n";
		$method = 'mean';
	}
}

if (defined $position) {
	# check the position value
	unless (
		$position == 5 or
		$position == 3 or
		$position == 1 or
		$position eq 'm'
	) {
		die " Unknown relative position '$position'!\n";
	}
	if ($position eq 'm') {$position = 1} # change to match internal usage
}
else {
	# default position to use the 5' end
	$position = 5;
}

if (defined $strand) {
	unless (
		$strand eq 'sense' or
		$strand eq 'antisense' or
		$strand eq 'all'
	) {
		die " Unknown strand value '$strand'!\n";
	}
}
else {
	# default
	$strand = 'all';
}

unless (defined $sum) {
	# assume to write a summary file, nearly always want this, at least I do
	$sum = 1;
}




## Generate or load the input dataset
my $main_data_ref; # a reference to the main data structure
if ($infile) {
	# load the gene dataset from existing file
	print " Loading feature set from file $infile....\n";
	$main_data_ref = load_tim_data_file($infile);
	
	# update program name
	unless ($main_data_ref->{'program'} eq $0) {
		$main_data_ref->{'program'} = $0;
	}
} 
else {
	# we will start a new file with a new dataset
	generate_a_new_feature_dataset();
}
unless ($main_data_ref) {
	# check for data
	die " No data loaded! Nothing to do!\n";
}
unless ($database) {
	# define the database from the input file if necessary
	if ( $main_data_ref->{'db'} ) {
		$database = $main_data_ref->{'db'};
	}
	else {
		die " A database must be defined!\n";
	}
}

# simple reference to the data table
my $data_table_ref = $main_data_ref->{'data_table'};

# the number of columns already in the data array
my $startcolumn = $main_data_ref->{'number_columns'}; 





## Collect the data

# Open database connection
my $db = open_db_connection($database);
unless ($db) {
	die " no database connection opened!\n";
}

# Check the dataset
set_and_verify_dataset();

my $start_time = time;

# open a file for writing out raw data files for debugging purposes
if ($raw) { 
	# this is only done if specifically requested
	print " Preparing raw output file....\n";
	open RAWFILE, ">$outfile\_raw.txt";
}

# Collect the nucleosome occupancy
map_relative_data();

# Close the raw file if it was opened
if ($raw) {close RAWFILE}



## Interpolate values
if ($smooth) {
	print " Interpolating missing values....\n";
	go_interpolate_values();
}



## Generate summed data - 
# an average across all features at each position suitable for plotting
if ($sum) {
	print " Generating final summed data....\n";
	my $sumfile = write_summary_data( {
		'data'        => $main_data_ref,
		'filename'    => $outfile,
		'startcolumn' => $startcolumn,
		'dataset'     => $dataset,
		'log'         => $log,
	} );
	if ($sumfile) {
		print " Wrote summary file '$sumfile'\n";
	}
	else {
		print " Unable to write summary file!\n";
	}
}



## Output the data
# reset the program metadata
$main_data_ref->{'program'} = $0;

# we will write a standard tim data file
# appropriate extensions and compression should be taken care of
my $written_file = write_tim_data_file( {
	'data'     => $main_data_ref,
	'filename' => $outfile,
	'gz'       => $gz,
} );
if ($written_file) {
	# success!
	print " wrote file $written_file\n";
}
else {
	# failure! the subroutine will have printed error messages
	print " unable to write output file!\n";
}


## Conclusion
my $timediff = sprintf "%.1f", (time - $start_time)/60;
print " Completed in $timediff minutes\n";



#### Subroutines #######


## validate the given dataset, or ask for one
sub set_and_verify_dataset {
	if ($dataset) {
		
		# check for a remote file
		if ($dataset =~ /^http|ftp/) {
			# a remote file
			# should be good, no verification here though
			return;
		}
		elsif ($dataset =~ /\.(?:bw|bb|bam)$/i) {
			# looks like we have a file 
			if (-e $dataset) {
				# file exists
				$dataset = "file:$dataset";
				return;
			}
			else {
				# maybe it's a funny named dataset?
				if (validate_dataset_list($db, $dataset) ) {
					# returned true, the name of the bad dataset
					die " The requested file or dataset '$dataset' " . 
						"neither exists or is valid!\n";
				}
			}
		}
		else {
			# must be a database feature type
		
			# validate the given dataset
			my $bad_dataset = validate_dataset_list($db, $dataset);
			if ($bad_dataset) {
				die " The requested dataset $bad_dataset is not valid!\n";
			}
			return;
		}
	} 
	
	else {
		# dataset not specified
		# present the dataset list to the user and get an answer
		my %datasethash = get_dataset_list($db);
		print "\n These are the available data sets in the database $database:\n";
		foreach (sort {$a <=> $b} keys %datasethash) {
			# print out the list of microarray data sets
			print "  $_\t$datasethash{$_}\n"; 
		}
		print " Enter the number of the data set you would like to analyze  ";
		my $answer = <STDIN>;
		chomp $answer;
		if (exists $datasethash{$answer}) {
			$dataset = $datasethash{$answer};
		} 
		else {
			die " That number doesn't correspond to a data set!\n";
		}
	}
}



## generate a feature dataset if one was not loaded
sub generate_a_new_feature_dataset {
	# a subroutine to generate a new feature dataset
	
	$main_data_ref = get_new_feature_list( {
			'db'       => $database,
			'features' => $feature,
	} );
	
	# set the current program
	$main_data_ref->{'program'} = $0;
}



## Collect the nucleosome occupancy data
sub map_relative_data {
	
	### Identify columns for feature identification
	# name
	my $name_index = find_column_index($main_data_ref, '^name');
	unless (defined $name_index) {
		die 'unable to identify Name column in data table!';
	}
	# type
	my $type_index = find_column_index($main_data_ref, 'type');
	unless (defined $type_index) {
		die 'unable to identify Type column in data table!';
	}
	# strand if requested
	my $strand_index;
	if ($set_strand) {
		$strand_index = find_column_index($main_data_ref, 'strand');
		unless (defined $strand_index) {
			die 'unable to strand column in data table!';
		}
	}
	
	### Prepare window values
	my $startingpoint = 0 - ($win * $number); 
		# default values will give startingpoint of -1000
	my $endingpoint = $win * $number; 
		# likewise default values will give endingpoint of 1000
	print " Collecting $dataset ";
	if ($log) { 
		print "log2 ";
	}
	print "data from \n   $startingpoint to $endingpoint ";
	if ($position == 3) {
		print "relative to the 3' end in $win bp increments...\n";
	}
	if ($position eq 'm') {
		print "relative to the midpoint in $win bp increments...\n";
	}
	else {
		print "in $win bp increments...\n";
	}
	

	### Prepare and annotate the header names and metadata
	for (my $start = $startingpoint; $start < $endingpoint; $start += $win) {
		# we will be progressing from the starting to ending point
		# in increments of window size
		
		# set the stop position
		my $stop = $start + $win - 1;
		
		# deal with the pesky 0 value
		# since we're working with 1-base coordinates, we don't really have 
		# a real 0 coordinate, so need to skip it as it doesn't really exist
		my $zero_check;
		for my $i ($start .. $stop) {
			$zero_check = 1 if $i == 0;
		}
		if ($zero_check) {
			# adjust the coordinates accordingly
			if ($start == 0) {
				$start += 1;
				$stop += 1;
			}
			elsif ($stop == 0) {
				$stop += 1;
			}
			else {
				# some number in between
				$stop += 1;
			}
		}
		
		# set the new index value, which is equivalent to the number of columns
		my $new_index = $main_data_ref->{'number_columns'};
		
		# the new name
		my $new_name = $start . '..' . $stop;
		
		# set the metadata
		
		my %metadata = (
			'name'        => $new_name,
			'index'       => $new_index,
			'start'       => $start,
			'stop'        => $stop,
			'win'         => $win,
			'log2'        => $log,
			'dataset'     => $dataset,
			'method'      => $method,
			'value'       => $value_type,
		);
		if ($position == 5) {
			$metadata{'relative_position'} = '5prime_end';
		}
		elsif ($position == 3) {
			$metadata{'relative_position'} = '3prime_end';
		}
		else { # midpoint
			$metadata{'relative_position'} = 'center';
		}
		if ($set_strand) {
			$metadata{'strand_implied'} = 1;
		}
		$main_data_ref->{$new_index} = \%metadata;
		
		# set the column header
		$data_table_ref->[0][$new_index] = $new_name;
		
		# update number of columns
		$main_data_ref->{'number_columns'} += 1;
	}
	
	### Collect the data
	for my $row (1..$main_data_ref->{'last_row'}) {
		
		# determine the region
		my $name = $data_table_ref->[$row][$name_index]; # name
		my $type = $data_table_ref->[$row][$type_index]; # type
		if ($raw) { 
			print RAWFILE join "\t", @{ $data_table_ref->[$row] };
		}
		
		# collect the region scores
		my %regionscores = get_region_dataset_hash( {
				'db'          => $db,
				'dataset'     => $dataset,
				'name'        => $name,
				'type'        => $type,
				'start'       => $startingpoint,
				'stop'        => $endingpoint,
				'position'    => $position,
				'value'       => $value_type,
				'strand'      => $strand,
				'set_strand'  => $set_strand ? 
								$data_table_ref->[$row][$strand_index] : undef, 
		} );
		
		# debugging
		if ($raw) {
			print RAWFILE " found ", scalar keys %regionscores, 
				" datapoints for '$name'\n";
			print RAWFILE Dumper(\%regionscores);
		}
		
		# assign the scores to the windows in the region
		for (
			# we will process each window one at a time
			# proceed by the column index for each window
			my $column = $startcolumn; 
			$column < $main_data_ref->{'number_columns'}; 
			$column++
		) {
			# get start and stop
			my $start = $main_data_ref->{$column}{'start'};
			my $stop = $main_data_ref->{$column}{'stop'};
			
			# collect a score at each position in the window
			my @scores;
			for (my $n = $start; $n <= $stop; $n++) {
				# we will walk through the window one bp at a time
				# look for a score associated with the position
				push @scores, $regionscores{$n} if exists $regionscores{$n};
			}
			
			# deal with log scores if necessary
			if ($log) {
				@scores = map { 2 ** $_ } @scores;
			}
			
			# calculate the combined score for the window
			my $winscore;
			if (@scores) {
				# we have scores, so calculate a value based on the method
				
				if ($method eq 'mean') {
					$winscore = mean(@scores);
				}
				elsif ($method eq 'median') {
					$winscore = median(@scores);
				}
				elsif ($method eq 'stddev') {
					$winscore = stddevp(@scores);
				}
				elsif ($method eq 'sum') {
					$winscore = sum(@scores);
				}
				elsif ($method eq 'min') {
					$winscore = sum(@scores);
				}
				elsif ($method eq 'max') {
					$winscore = sum(@scores);
				}
				
				# deal with log2 scores
				if ($log) { 
					# put back in log2 space if necessary
					$winscore = log($winscore) / log(2);
				}
			}
			else {
				# no scores
				# assign a "null" value
				$winscore = '.';
			}
			
			# put the value into the data table
			# we're using a push function instead of explicitly assigning 
			# a position, since the loop is based on relative genomic 
			# position rather than 
			$data_table_ref->[$row][$column] = $winscore;
			if ($raw) { print RAWFILE "\t$winscore" }
		}
		if ($raw) { 
			# finish the raw data line
			print RAWFILE "\n";
		}
	}
	
}



## Interpolate the '.' values with the mean of the neighbors
sub go_interpolate_values {
	
	# determine counts
	my $lastwindow = $main_data_ref->{'number_columns'} - 2; 
		# lastwindow is the index of the second to last column
	
	# walk through each data line and then each window
	for my $row (1..$main_data_ref->{'last_row'}) {
		
		for my $col ($startcolumn..$lastwindow) {
			# walk through the windows of a data row
			# skipping the very first and last windows (columns)
			# we will look for null values
			# if one is found, interpolate from neighbors
			if ($data_table_ref->[$row][$col] eq '.') {
				# we will interpolate the value from the neighbors
				# first, need to check that the neighbors have values
				
				# obtain the neighboring values
				my $before = $data_table_ref->[$row][$col - 1];
				my $after = $data_table_ref->[$row][$col + 1];
				if (
					$before ne '.' and 
					$after ne '.'
				) {
					# neighboring values are not nulls
					# then calculate the interpolated value
					if ($log) { 
						# working with log2 values
						$before = 2 ** $before;
						$after = 2 ** $after;
						my $mean = mean($before, $after);
						$data_table_ref->[$row][$col] = log($mean) / log(2);
					} 
					else { 
						# otherwise calculate straight value
						$data_table_ref->[$row][$col] = mean($before, $after);
					}
				}
			}
			
		}
		
	}
	
}




__END__

=head1 NAME

map_data.pl

A script to map data relative to and flanking a genomic feature

=head1 SYNOPSIS
 
 map_data.pl --db <database> --feature <name> --out <file> [--options]
 map_data.pl --in <in_filename> --out <out_filename> [--options]
  
  Options:
  --db <name|file.gff3>
  --feature [type, type:source]
  --in <filename> 
  --out <filename>
  --data <dataset_name | filename>
  --method [mean|median|min|max|stddev|sum]
  --value [score|count|length]
  --win <integer>
  --num <integer>
  --pos [5|3|m]
  --strand [sense|antisense|all]
  --set_strand
  --(no)sum
  --(no)smooth
  --(no)log
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4


=item --db <name|file.gff3>

Specify the name of the BioPerl SeqFeature::Store database to use as
source. Alternatively, a single GFF3 file may be loaded into a in-memory
database. Specifying the database is required for new feature data files.
For pre-existing input data files, this value may be obtained from the
input file metadata. However, if provided, it overrides the database listed
in the file; this is useful for collecting data from multiple databases.

=item --out <filename>

Specify the output file name. Required for new files; otherwise, 
input files will be overwritten unless specified.

=item --in <filename>

Specify the filename of a data table containing the list of 
features to map nucleosomes around. The file must be in the
'tim_data' format and specify a feature to use. If an input 
file is not specified, then a new list of features will be 
generated from the database.

=item --feature [type, type:source]

Specify the type of feature to map data around. The feature may be 
listed either as GFF type or GFF type:source. The list 
of features will be automatically generated from the database. 

This is optional if the features are defined in the input file.

=item --data <dataset_name | filename>

Specify the name of the data set from which you wish to 
collect data. If not specified, the data set may be chosen
interactively from a presented list. Other
features may be collected, and should be specified using the type 
(GFF type:source), especially when collecting alternative data values. 
Alternatively, the name of a data file may be provided. Supported 
file types include BigWig (.bw), BigBed (.bb), or single-end Bam 
(.bam). The file may be local or remote.

=item --method [mean|median|min|max|stddev|sum]

Specify the method of combining multiple values within each window. The mean, 
median, minimum, maximum, standard deviation, or sum of the values may be 
reported. The default value is mean for score and length values, or sum for 
count values.

=item --value [score|count|length]

Optionally specify the type of data value to collect from the dataset or 
data file. Three values are accepted: score, count, or length. The default 
value type is score. Note that some data sources only support certain 
types of data values. Wig and BigWig files only support score and count; 
BigBed and database features support count and length and optionally 
score; Bam files support basepair coverage (score), count (number of 
alignments), and length.

=item --win <integer>

Specify the window size. Default 50

=item --num <integer>

Specify the number of windows on either side of the feature position 
(total number will be 2 x [num]). Default 20

=item --pos [5|3|m]

Indicate the relative position of the feature around which the 
data is mapped. Three values are accepted: "5" indicates the 
5' prime end is used, "3" indicates the 3' end is used, and "m" or "1"
indicates the middle of the feature is used. The default is to 
use the 5' end, or the start position of unstranded features. 
If the feature "tts" is selected above, the 3' end is 
automatically selected.

=item --strand [sense|antisense|all]

Specify whether stranded data should be collected for each of the 
datasets. Either sense or antisense (relative to the feature) data 
may be collected. The default value is 'all', indicating all 
data will be collected.

=item --set_strand

For features that are not inherently stranded (strand value of 0), 
impose an artificial strand for each feature (1 or -1). This will 
have the effect of enforcing a relative orientation for each feature, 
or to collected stranded data. This requires the presence a 
column in the input data file with a name of "strand". Hence, it 
will not work with newly generated datasets, but only with input 
data files. Default is false.

=item --(no)sum

Indicate that the data should be averaged across all features at
each position, suitable for graphing. A separate text file will 
be written with the suffix '_summed' with the averaged data. 
Default is true.

=item --(no)smooth

Indicate that windows without values should (not) be interpolated
from neighboring values.

=item --(no)log

Dataset values are (not) in log2 space and should be treated 
accordingly. Output values will be in the same space.

=item --help

Display this help


=back

=head1 DESCRIPTION

This program will map data around a genomic feature. Features
may include the transcription start site (TSS) of a transcript, tRNA gene,
or some other genomic feature. The data is collected in a series of windows
flanking the feature start, end, or midpoint position. The number and size of
windows are specified on the command line, or the program will default to
20 windows (on either side of the feature, 40 total) of 50 bp size. These
default window settings corresponds to 1 kb on either side of the feature.
Windows without a value may be interpolated (smoothed) from neigboring 
values, if available.

The default value that is collected is a dataset score (e.g. microarray 
values). However, other values may be collected, including 'count' or 
'length'. Use the --method argument to collect alternative values.

The program writes out a tim data formatted text file.


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






