#!/usr/bin/perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	calculate_score
);
use Bio::ToolBox::utility;
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};
use constant DATASET_HASH_LIMIT => 20001;
		# This constant determines the maximum size of the dataset hash to be 
		# returned from the get_region_dataset_hash(). To increase performance, 
		# the program normally queries the database once for each feature or 
		# region, and a hash returned with potentially a score for each basepair. 
		# This may become unwieldy for very large regions, which may be better 
		# served by separate database queries for each window.
my $VERSION = '1.54';

print "\n A script to collect windowed data flanking a relative position of a feature\n\n";
  

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
	$main_database, 
	$data_database,
	$dataset, 
	$feature, 
	$method,
	$win, 
	$number,
	$position, 
	$strand_sense,
	$set_strand,
	$avoid,
	$avoidtype,
	$long_data,
	$smooth,
	$sum,
	$gz,
	$cpu,
	$help,
	$print_version,
); # command line variables

## Command line options
GetOptions( 
	'out=s'      => \$outfile, # output file name
	'in=s'       => \$infile, # input file name
	'db=s'       => \$main_database, # main or annotation database name
	'ddb=s'      => \$data_database, # data database
	'data=s'     => \$dataset, # dataset name
	'feature=s'  => \$feature, # type of feature
	'method=s'   => \$method, # method to combine data
	'window=i'   => \$win, # window size
	'number=i'   => \$number, # number of windows
	'position=s' => \$position, # indicate relative location of the feature
	'strand=s'   => \$strand_sense, # collected stranded data
	'force_strand|set_strand' => \$set_strand, # enforce an artificial strand
				# force_strand is preferred option, but respect the old option
	'avoid!'     => \$avoid, # avoid conflicting features
	'avtype=s'   => \$avoidtype, # type of feature to avoid
	'long!'      => \$long_data, # collecting long data features
	'smooth!'    => \$smooth, # smooth by interpolation
	'sum!'       => \$sum, # generate average profile
	'gz!'        => \$gz, # compress the output file
	'cpu=i'      => \$cpu, # number of execution threads
	'help'       => \$help, # print help
	'version'    => \$print_version, # print the version
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
	print " Biotoolbox script map_data.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}






## Check for required values
check_defaults();
my $start_time = time;




## Generate or load the input dataset
my $Data;
if ($infile) {
	$Data = Bio::ToolBox::Data->new(
		file       => $infile, 
		parse      => 1,
		feature    => $feature || 'gene',
	) or die " unable to load input file '$infile'\n";
	if ($Data->last_row) {
		printf " Loaded %s features from $infile.\n", format_with_commas( $Data->last_row );
	}
	else {
		die " No features loaded!\n";
	}
	
	# update main database as necessary
	if ($main_database) {
		if (defined $Data->database and $Data->database ne $main_database) {
			# update with new database
			printf " updating main database name from '%s' to '%s'\n", 
				$Data->database, $main_database;
# 			print "   Re-run without --db option if you do not want this to happen\n";
			$Data->database($main_database);
		}
	}
	else {
		$main_database = $Data->database;
	}
	
	# update feature type as necessary
	if (not defined $Data->feature and not defined $Data->type_column and defined $feature) {
		$Data->feature($feature);
	}
}
else {
	# generate a new file
	print " Generating a new feature list from database '$main_database'...\n";
	$Data = Bio::ToolBox::Data->new(
		db      => $main_database,
		feature => $feature,
	) or die " unable to generate new feature list\n";
}
$Data->program("$0, v $VERSION");

# the number of columns already in the data array
my $startcolumn = $Data->number_columns; 

# make sure data table supports avoid option
if ($avoid) {
	unless ($main_database) {
		warn " avoid option not supported without an annotation database! Disabling\n";
		undef $avoid; # 0 complicates things, leave undefined
	}
	if ($avoidtype) {
		if ($avoidtype =~ /,/) {
			# a comma delimited list
			$avoid = [ split(',', $avoidtype) ];
		}
		else {
			# presume a single feature type, take it as is
			$avoid = [ $avoidtype ];
		}
	}
}


## Prepare to collect data
# Open data database
my $ddb;
if (defined $data_database) {
	# specifically defined a data database
	$ddb = open_db_connection($data_database) or 
		die "unable to establish data database connection to $data_database!\n";
}

# Check the dataset
$dataset = verify_or_request_feature_types(
	'db'      => $ddb || $Data->database,
	'feature' => $dataset,
	'single'  => 1,
	'prompt'  => " Enter the number of the feature or dataset from which to" . 
					" collect data   ",
);
unless ($dataset) {
	die " No verifiable dataset provided. Check your file path, database, or dataset.\n";
} 



## Collect the relative data

# check whether it is worth doing parallel execution
if ($cpu > 1) {
	while ($cpu > 1 and $Data->last_row / $cpu < 100) {
		# I figure we need at least 100 lines in each fork split to make 
		# it worthwhile to do the split, otherwise, reduce the number of 
		# splits to something more worthwhile
		$cpu--;
	}
}

if ($cpu > 1) {
	# parallel execution
	print " Forking into $cpu children for parallel data collection\n";
	parallel_execution();
}
else {
	# single process execution
	single_execution();
}


## Generate summed data 
# an average across all features at each position suitable for plotting
if ($sum) {
	print " Generating final summed data....\n";
	my $sumfile = $Data->summary_file(
		'filename'    => $outfile,
		'startcolumn' => $startcolumn,
		'dataset'     => $dataset,
	);
	if ($sumfile) {
		print " Wrote summary file '$sumfile'\n";
	}
	else {
		print " Unable to write summary file!\n";
	}
}


## Output the data
unless ($outfile) {
	$outfile = $Data->path . $Data->basename;
}
my $written_file = $Data->save(
	'filename' => $outfile,
	'gz'       => $gz,
);
if ($written_file) {
	# success!
	print " wrote file $written_file\n";
}
else {
	# failure! the subroutine will have printed error messages
	print " unable to write output file!\n";
}
printf " Completed in %.1f minutes\n", (time - $start_time)/60;





#### Subroutines #######

## check required variables and assign default values
sub check_defaults {
	unless ($main_database or $infile) {
		die " You must define a database or input file!\n Use --help for more information\n";
	}

	unless ($outfile or $infile) {
		die " You must define an output filename !\n Use --help for more information\n";
	}

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

	if (defined $method) {
		# check the requested method
		unless (
				$method eq 'mean' or
				$method eq 'median' or
				$method eq 'sum' or
				$method eq 'min' or
				$method eq 'max' or
				$method eq 'stddev' or
				$method =~ /^\w?count$/
		) {
			die " Unknown method '$method'!\n Use --help for more information\n";
		}
	}
	else {
		$method = 'mean';
	}

	if (defined $position) {
		# check the position value
		unless (
			$position == 5 or
			$position == 3 or
			$position == 4 or
			$position eq 'm'
		) {
			die " Unknown relative position '$position'!\n";
		}
		if ($position eq 'm') {$position = 4} # change to match internal usage
	}
	else {
		# default position to use the 5' end
		$position = 5;
	}

	if (defined $strand_sense) {
		unless (
			$strand_sense eq 'sense' or
			$strand_sense eq 'antisense' or
			$strand_sense eq 'all'
		) {
			die " Unknown strand value '$strand_sense'!\n";
		}
	}
	else {
		# default
		$strand_sense = 'all';
	}

	unless (defined $sum) {
		# assume to write a summary file, nearly always want this, at least I do
		$sum = 1;
	}

	if ($parallel) {
		$cpu ||= 4;
	}
	else {
		# disable cores
		print " disabling parallel CPU execution, no support present\n" if $cpu;
		$cpu = 1;
	}
}


## Run in parallel
sub parallel_execution {
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_start( sub { sleep 1; }); 
		# give a chance for child to start up and open databases, files, etc 
		# without creating race conditions
	
	# generate base name for child processes
	my $child_base_name = $outfile . ".$$"; 

	# Split the input data into parts and execute in parallel in separate forks
	for my $i (1 .. $cpu) {
		$pm->start and next;
		
		#### In child ####
	
		# splice the data structure
		$Data->splice_data($i, $cpu);
		
		# re-open database objects to make them clone safe
		# pass second true to avoid cached database objects
		my $db = $Data->open_database(1);
		if ($data_database) {
			$ddb = open_db_connection($data_database, 1);
		}
		
		# Add the columns for each window 
		# and calculate the relative starting and ending points
		my ($starting_point, $ending_point) = prepare_window_datasets();
	
		# determine long data collection for very large regions
		if ($ending_point - $starting_point > DATASET_HASH_LIMIT) {
			$long_data = 1;
		}
	
		# Select the appropriate method for data collection
		if ($long_data) {
			map_relative_long_data($starting_point, $ending_point);
		}
		else {
			map_relative_data($starting_point, $ending_point);
		}

		# Interpolate values
		if ($smooth) {
			print " Interpolating missing values....\n";
			go_interpolate_values();
		}
		
		# write out result
		my $success = $Data->save(
			'filename' => sprintf("$child_base_name.%03s",$i),
			'gz'       => 0, # faster to write without compression
		);
		if ($success) {
			printf " wrote child file $success\n";
		}
		else {
			# failure! the subroutine will have printed error messages
			die " unable to write file!\n";
			# no need to continue
		}
		
		# Finished
		$pm->finish;
	}
	$pm->wait_all_children;
		
	# reassemble children files into output file
	my @files = glob "$child_base_name.*";
	unless (@files) {
		die "unable to find children files!\n";
	}
	unless (scalar @files == $cpu) {
		die "only found " . scalar(@files) . " child files when there should be $cpu!\n";
	}
	my $count = $Data->reload_children(@files);
	printf " reloaded %s features from children\n", format_with_commas($count);
}


## Run in single thread
sub single_execution {
	
	# Add the columns for each window 
	# and calculate the relative starting and ending points
	my ($starting_point, $ending_point) = prepare_window_datasets();

	# determine long data collection for very large regions
	if ($ending_point - $starting_point > DATASET_HASH_LIMIT) {
		$long_data = 1;
	}

	# Select the appropriate method for data collection
	if ($long_data) {
		map_relative_long_data($starting_point, $ending_point);
	}
	else {
		map_relative_data($starting_point, $ending_point);
	}

	# Interpolate values
	if ($smooth) {
		print " Interpolating missing values....\n";
		go_interpolate_values();
	}
}


## Prepare columns for each window
sub prepare_window_datasets {
	
	# Determine starting and ending points
	my $starting_point = 0 - ($win * $number); 
		# default values will give startingpoint of -1000
	my $ending_point = $win * $number; 
		# likewise default values will give endingpoint of 1000
	
	# Print collection statement
	print " ";
	printf " Collecting data from %d to %d at the %s in %d bp windows...\n", 
		$starting_point, $ending_point, $position == 3 ? "3' end" : 
		$position == 4 ? "midpoint" : "5' end", $win;
	
	# Prepare and annotate the header names and metadata
	for (my $start = $starting_point; $start < $ending_point; $start += $win) {
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
		
		# the new name
		my $new_name = $start . '..' . $stop;
		
		# add new column
		my $new_index = $Data->add_column($new_name);
		
		# set the metadata
		$Data->metadata($new_index, 'start' , $start);
		$Data->metadata($new_index, 'stop' , $stop);
		$Data->metadata($new_index, 'window' , $win);
		$Data->metadata($new_index, 'dataset' , $dataset);
		$Data->metadata($new_index, 'method' , $method);
		if ($position == 5) {
			$Data->metadata($new_index, 'relative_position', '5prime_end');
		}
		elsif ($position == 3) {
			$Data->metadata($new_index, 'relative_position', '3prime_end');
		}
		else { # midpoint
			$Data->metadata($new_index, 'relative_position', 'center');
		}
		if ($set_strand) {
			$Data->metadata($new_index, 'strand_implied', 1);
		}
		if ($strand_sense =~ /sense/) {
			$Data->metadata($new_index, 'strand',  $strand_sense);
		}
		if ($data_database) {
			$Data->metadata($new_index, 'db', $data_database);
		}
		if ($avoid) {
			$Data->metadata($new_index, 'avoid', 
				ref($avoid) eq 'ARRAY' ? join(',', @$avoid) : $avoid
			);
		}
	}
	
	return ($starting_point, $ending_point);
}


## Collect relative scores using single database call
sub map_relative_data {
	my ($starting_point, $ending_point) = @_;
	
	# Collect the data
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		my $regionscores = $row->get_relative_point_position_scores(
			'ddb'         => $ddb,
			'dataset'     => $dataset,
			'position'    => $position,
			'extend'      => $ending_point, # equivalent to the extension
			'stranded'    => $strand_sense,
			'strand'      => $set_strand ? $row->strand : undef, 
			'avoid'       => $avoid,
		);
		
		# assign the scores to the windows in the region
		for (
			# we will process each window one at a time
			# proceed by the column index for each window
			my $column = $startcolumn; 
			$column < $Data->number_columns; 
			$column++
		) {
		
			# record nulls if no data returned
			unless (scalar keys %$regionscores) {
				$row->value($column, '.');
				next;
			}
		
			# get start and stop
			my $start = $Data->metadata($column, 'start');
			my $stop = $Data->metadata($column, 'stop');
		
			# collect a score at each position in the window
			my @scores = 	map { $regionscores->{$_} } 
							grep { $_ >= $start and $_ <= $stop}
							keys %$regionscores;
			
			# put the value into the data table
			$row->value($column, calculate_score($method, \@scores) );
		}
	}
}


## Collect relative data using individual window database calls
sub map_relative_long_data {
	my ($starting_point, $ending_point) = @_;
	
	# Collect the data
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		
		# get feature from the database if necessary
		my $feature = $row->seqfeature || $row;
		
		# calculate a reference point
		my $reference;
		if ($feature->strand >= 0 and $position == 5) {
			# 5' end of forward strand
			$reference = $feature->start;
		}
		elsif ($feature->strand == -1 and $position == 5) {
			# 5' end of reverse strand
			$reference = $feature->stop;
		}
		elsif ($feature->strand >= 0 and $position == 3) {
			# 3' end of forward strand
			$reference = $feature->stop;
		}
		elsif ($feature->strand == -1 and $position == 3) {
			# 3' end of reverse strand
			$reference = $feature->start;
		}
		elsif ($position == 4) {
			# midpoint regardless of strand
			$reference = int( ( ($feature->stop + $feature->start) / 2) + 0.5);
		}
		else {
			# something happened
			die " programming error!? feature " . 
				" at data row $row\n";
		}
		
		# collect the data for every window 
		for (
			my $column = $startcolumn; 
			$column < $Data->number_columns; 
			$column++
		) {
			# we must modify the start and stop position with the adjustments
			# recorded in the current column metadata
			my $score = $row->get_score(
				'db'          => $ddb,
				'dataset'     => $dataset,
				'chromo'      => $feature->seq_id,
				'start'       => $feature->strand >= 0 ? 
									$reference + $Data->metadata($column, 'start') :
									$reference - $Data->metadata($column, 'start'),
				'stop'        => $feature->strand >= 0 ? 
									$reference + $Data->metadata($column, 'stop') : 
									$reference - $Data->metadata($column, 'stop'),
				'strand'      => $feature->strand,
				'method'      => $method,
				'stranded'    => $strand_sense,
			);
			$row->value($column, $score);
		}
	}
}



## Interpolate the '.' values with the mean of the neighbors
sub go_interpolate_values {
	
	# determine counts
	my $lastwindow = $Data->last_column; 
		# lastwindow is the index of the last column
	# walk through each data line and then each window
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		my $col = $startcolumn + 1;
		while ($col < $lastwindow) {
			# walk through the windows of a data row
			# skipping the very first and last windows (columns)
			# we will look for null values
			# if one is found, interpolate from neighbors
			if ($row->value($col) eq '.' and $row->value($col - 1) ne '.') {
				# we will interpolate the value from the neighbors
				# first, need to check that the neighbors have values
				
				# find the next real value
				my $next_i;
				for (my $i = $col + 1; $i <= $lastwindow; $i++) {
					if ($row->value($i) ne '.') {
						$next_i = $i;
						last;
					}
				}
				unless (defined $next_i) {
					$col++;
					next;
				}
				
				# determine fractional value
				my $initial = $row->value($col - 1);
				my $fraction = ($row->value($next_i) - $initial) / 
					($next_i - $col + 1);
				
				# apply fractional values
				for (my $n = $col; $n < $next_i; $n++) {
					$row->value($n, $initial + ($fraction * ($n - $col + 1)) );
				}
				
				# jump ahead
				$col = $next_i;
			}
			$col++;
		}
	}
}



__END__

=head1 NAME

get_relative_data.pl

A script to collect data in bins around a relative position.

=head1 SYNOPSIS
 
get_relative_data.pl --in <in_filename> --out <out_filename> [--options]
  
  Options for existing files:
  --in <filename>                  (txt bed gff gtf refFlat ucsc)
  
  Options for new files:
  --db <name|file>
  --feature <type | type:source | alias>, ...
  
  Options for data collection:
  --ddb <name|file>
  --data <dataset_name | filename>
  --method [mean|median|stddev|min|max|range|sum|          (mean)
            count|pcount|ncount]
  --strand [all|sense|antisense]                            (all)
  --force_strand
  --avoid
  --avtype [type,type,...]
  --long
  
  Bin specification:
  --win <integer>                                           (50)
  --num <integer>                                           (20)
  --pos [5|m|3]                                             (5)
  
  Post-processing:
  --(no)sum                                                 (true)
  --smooth
  
  General Options:
  --out <filename>
  --gz
  --cpu <integer>                                           (2)
  --version
  --help                              show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. Any tab-delimited text 
file with recognizable headers is supported. Gene annotation file 
formats are also supported, including bed, gtf, gff3, refFlat, and 
UCSC native formats such as gene prediction tables are all supported. 
Gene annotation files will be parsed as sequence features. 
Files may be gzipped compressed.

=item --out <filename>

Specify the output file name. Required for new files; otherwise, 
input files will be overwritten unless specified.

=item --db <name | filename>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A database is 
required for generating new data files with features. This option may 
skipped when using coordinate information from an input file (e.g. BED 
file), or when using an existing input file with the database indicated 
in the metadata. For more information about using annotation databases, 
see L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. 

=item --feature [type, type:source]

Specify the type of feature to map data around. The feature may be 
listed either as GFF type or GFF type:source. The list 
of features will be automatically generated from the database. 
This is only required when an input file is not specified. 

=item --ddb <name | filename>

If the data to be collected is from a second database that is separate 
from the annotation database, provide the name of the data database here. 
Typically, a second C<Bio::DB::SeqFeature::Store> or BigWigSet database 
is provided here. 

=item --data <dataset_name | filename>

Specify the name of the data set from which you wish to 
collect data. If not specified, the data set may be chosen
interactively from a presented list. Other
features may be collected, and should be specified using the type 
(GFF type:source), especially when collecting alternative data values. 
Alternatively, the name of a data file may be provided. Supported 
file types include BigWig (.bw), BigBed (.bb), or single-end Bam 
(.bam). The file may be local or remote.

=item --method <text>

Specify the method for combining all of the dataset values within the 
genomic region of the feature. Accepted values include:

=over 4

=item * mean (default)

=item * median

=item * sum

=item * stddev  Standard deviation of the population (within the region)

=item * min

=item * max

=item * range   Returns difference of max and min

=item * count

Counts the number of overlapping items.

=item * pcount (precise count)

Counts the number of items that precisely fall within the query 
region. Partially overlapping are not counted.

=item * ncount (name count)

Counts unique names. Useful when spliced alignments overlap more 
than one exon and you want to avoid double-counting.

=back

=item --strand [sense|antisense|all]

Specify whether stranded data should be collected for each of the 
datasets. Either sense or antisense (relative to the feature) data 
may be collected. The default value is 'all', indicating all 
data will be collected.

=item --force_strand

For features that are not inherently stranded (strand value of 0)
or that you want to impose a different strand, set this option when
collecting stranded data. This will reassign the specified strand for
each feature regardless of its original orientation. This requires the
presence of a "strand" column in the input data file. This option only
works with input file lists of database features, not defined genomic
regions (e.g. BED files). Default is false.

=item --avoid

Indicate whether neighboring features should be avoided when calculating 
values in a window. After collecting the data, each window is checked for 
overlapping features; if the window overlaps another feature, no value 
is reported for that window. This option requires using an annotation 
database (--db option). This is useful to avoid scoring windows that 
overlap a neighboring gene, for example. The default is false (return 
all values regardless of overlap).

=item --avtype [type,type,...]

Provide a feature type (primary_tag or primary_tag:source) or a 
comma-delimited list of types to be used when avoiding neighboring 
features. The default is to avoid features of the same type as that 
of the query, i.e. collecting data around a feature of type 'gene' 
will avoid other 'gene' features. This option allows you to avoid 
other features too, such as 'tRNA' or 'repeat'.

=item --long

Indicate that the dataset from which scores are collected are 
long features (counting genomic annotation for example) and not point 
data (microarray data or sequence coverage). Normally long features are 
only recorded at their midpoint, leading to inaccurate representation at 
some windows. This option forces the program to collect data separately 
at each window, rather than once for each file feature or region and 
subsequently assigning scores to windows. This may result in counting 
features more than once if it overlaps more than one window, a result 
that may or may not be desired. Execution time will likely increase. 
Default is false.

=item --win <integer>

Specify the window size. The default is 50 bp.

=item --num <integer>

Specify the number of windows on either side of the feature position 
(total number will be 2 x [num]). The default is 20, or 1 kb on either 
side of the reference position if the default window size is used.

=item --pos [5|m|3]

Indicate the relative position of the feature around which the 
data is mapped. Three values are accepted: "5" indicates the 
5' prime end is used, "3" indicates the 3' end is used, and "m" 
indicates the middle of the feature is used. The default is to 
use the 5' end, or the start position of unstranded features. 

=item --(no)sum

Indicate that the data should be averaged across all features at
each position, suitable for graphing. A separate text file will 
be written with the suffix '_summed' with the averaged data. 
Default is true (sum).

=item --smooth

Indicate that windows without values should (not) be interpolated
from neighboring values. The default is false (nosmooth).

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --cpu <integer>

Specify the number of CPU cores to execute in parallel. This requires 
the installation of Parallel::ForkManager. With support enabled, the 
default is 2. Disable multi-threaded execution by setting to 1. 

=item --version

Print the version number.

=item --help

Display this help

=back

=head1 DESCRIPTION

This program will collect data around a relative coordinate of a genomic 
feature or region. The data is collected in a series of windows flanking the 
feature start (5' position for stranded features), end (3' position), or 
the midpoint position. The number and size of windows are specified via 
command line arguments, or the program will default to 20 windows on both 
sides of the relative position (40 total) of 50 bp size, corresponding to 
2 kb total (+/- 1 kb). Windows without a value may be interpolated 
(smoothed) from neigboring values, if available.

Stranded data may be collected. If the feature does not have an inherent 
strand, one may be specified to enforce stranded collection or a particular 
orientation. 

When features overlap, or the collection windows of one feature overlaps 
with another feature, then data may be ignored and not collected (--avoid).

=head1 EXAMPLES

These are some examples of some common scenarios for collecting data.

=over 4

=item Collect scores in intervals around start

You want to collect the mean score from a bigWig file in twenty 50 bp 
intervals flanking the start position of each feature in Bed file.

  get_relative_data.pl --data scores.bw --in input.bed

=item Collect scores in intervals around middle

You want to collect median scores in 20 bp intervals extending 500 bp 
from the midpoint of each feature.

  get_relative_data.pl --win 20 --num 25 --pos m --data scores.bw --in \
  input.txt

=item Collect scores in intervals from annotation database

You want to collect scores in intervals around the transcription start 
site of genes in an annotation database, but also avoid intervals that 
may overlap neighboring genes. You want to collect alignment counts 
from a Bam file in a stranded fashion. 

  get_relative_data.pl --db annotation --feature gene --avoid --strand \
  sense --method count --data alignments.bam --out gene_tss
    
=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
