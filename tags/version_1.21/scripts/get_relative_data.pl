#!/usr/bin/perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($Bin);
use Statistics::Lite qw(mean median sum stddevp min max);
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	check_dataset_for_rpm_support
);
use Bio::ToolBox::utility;
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};
use constant LOG2 => log(2);
use constant DATASET_HASH_LIMIT => 5001;
		# This constant determines the maximum size of the dataset hash to be 
		# returned from the get_region_dataset_hash(). To increase performance, 
		# the program normally queries the database once for each feature or 
		# region, and a hash returned with potentially a score for each basepair. 
		# This may become unwieldy for very large regions, which may be better 
		# served by separate database queries for each window.
my $VERSION = 1.21;

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
	$value_type,
	$method,
	$win, 
	$number,
	$position, 
	$strand_sense,
	$set_strand,
	$avoid,
	$long_data,
	$smooth,
	$sum,
	$log,
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
	'value=s'    => \$value_type, # the type of data to collect
	'method=s'   => \$method, # method to combine data
	'window=i'   => \$win, # window size
	'number=i'   => \$number, # number of windows
	'position=s' => \$position, # indicate relative location of the feature
	'strand=s'   => \$strand_sense, # collected stranded data
	'force_strand|set_strand' => \$set_strand, # enforce an artificial strand
				# force_strand is preferred option, but respect the old option
	'avoid!'     => \$avoid, # avoid conflicting features
	'long!'      => \$long_data, # collecting long data features
	'smooth!'    => \$smooth, # smooth by interpolation
	'sum!'       => \$sum, # generate average profile
	'log!'       => \$log, # data is in log2 space
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
	print " Biotoolbox script map_data.pl, version $VERSION\n\n";
	exit;
}






## Check for required values
check_defaults();
my $start_time = time;




## Generate or load the input dataset
my $Data;
if ($infile) {
	$Data = Bio::ToolBox::Data->new(file => $infile) or 
		die " unable to load input file '$infile'\n";
	printf "  Loaded %s '%s' features.\n", 
		format_with_commas( $Data->last_row ), $Data->feature;
	
	# update main database as necessary
	if ($main_database) {
		if ($main_database ne $Data->database) {
			# update with new database
			printf " updating main database name from '%s' to '%s'\n", 
				$Data->database, $main_database;
			print "   Re-run without --db option if you do not want this to happen\n";
			$Data->database($main_database);
		}
	}
	else {
		$main_database = $Data->database;
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

# the number of columns already in the data array
my $startcolumn = $Data->number_columns; 

# make sure data table supports avoid option
if ($avoid) {
	unless ($Data->feature_type eq 'named' and defined $Data->type_column) {
		warn " avoid option not supported with current Data table. Disabling\n";
		$avoid = 0;
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

# Check the RPM method if necessary
my $rpm_read_sum;
if ($method eq 'rpm') {
	print " Checking RPM support for dataset '$dataset'...\n";
	$rpm_read_sum = check_dataset_for_rpm_support($dataset, $ddb, $cpu);
		# this step can be multi-threaded
	if ($rpm_read_sum) {
		printf "   %s total features\n", format_with_commas($rpm_read_sum);
	}
	else {
		die " RPM method not supported! Try something else\n";
	}
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

## Conclusion
printf " Completed in %.1f minutes\n", (time - $start_time)/60;





#### Subroutines #######

## check required variables and assign default values
sub check_defaults {
	unless ($main_database or $infile) {
		die " You must define a database or input file!\n Use --help for more information\n";
	}

	unless ($outfile) {
		if ($infile) {
			$outfile = $infile;
		}
		else {
			die " You must define an output filename !\n Use --help for more information\n";
		}
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

	unless (defined $log) {
		# default is to assume not log2
		$log = 0;
	}

	if (defined $value_type) {
		# check the region method or type of data value to collect
		unless (
				$value_type eq 'score' or
				$value_type eq 'length' or
				$value_type eq 'count'
		) {
			die " Unknown data value '$value_type'!\n " . 
				"Use --help for more information\n";
		}
	}
	else {
		# default is to take the score
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
				$method eq 'stddev' or
				$method eq 'rpm'
		) {
			die " Unknown method '$method'!\n Use --help for more information\n";
		}
	
		if ($method eq 'rpm') {
			# make sure we collect the right values
			$value_type = 'count';
		}
	}
	else {
		# set default method
		if ($value_type eq 'count') {
			$method = 'sum';
		}
		else {
			$method = 'mean';
		}
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
		# conservatively enable 2 cores
		$cpu ||= 2;
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
		if ($data_database) {
			$ddb = open_db_connection($data_database, 1);
		}
		
		# Collect the data
		map_relative_data();

	
		# Interpolate values
		if ($smooth) {
			print " Interpolating missing values....\n";
			go_interpolate_values();
		}
		# convert null values to zero if necessary
		if ($method eq 'sum' or $method eq 'rpm') {
			null_to_zeros();
		}
		
		# write out result
		my $success = $Data->save(
			'filename' => "$child_base_name.$i",
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
	my @args = ("$Bin/join_data_file.pl", "--out", $outfile);
	push @args, '--gz' if $gz;
	push @args, @files;
	system(@args) == 0 or die " unable to execute join_data_file.pl! $?\n";
	unlink @files;
	
	# generate summary file
	if ($sum) {
		# reopen the combined file
		my $Data2 = Bio::ToolBox::Data->new(file => $outfile);
		unless ($Data2) {
			warn " cannot re-open $outfile to generate summary file!\n";
			return;
		}
		print " Generating final summed data....\n";
		my $sumfile = $Data2->summary_file(
			# it will automatically define a new output name
			'startcolumn' => $startcolumn,
			'dataset'     => $dataset,
			'log'         => $log,
		);
		if ($sumfile) {
			print " Wrote summary file '$sumfile'\n";
		}
		else {
			print " Unable to write summary file!\n";
		}
	}
	# done
}


## Run in single thread
sub single_execution {
	
	# Collect the data
	map_relative_data();

	
	# Interpolate values
	if ($smooth) {
		print " Interpolating missing values....\n";
		go_interpolate_values();
	}
	# convert null values to zero if necessary
	if ($method eq 'sum' or $method eq 'rpm') {
		null_to_zeros();
	}


	# Generate summed data - 
	# an average across all features at each position suitable for plotting
	if ($sum) {
		print " Generating final summed data....\n";
		my $sumfile = $Data->summary_file(
			'filename'    => $outfile,
			'startcolumn' => $startcolumn,
			'dataset'     => $dataset,
			'log'         => $log,
		);
		if ($sumfile) {
			print " Wrote summary file '$sumfile'\n";
		}
		else {
			print " Unable to write summary file!\n";
		}
	}


	## Output the data
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
	# done
}


## Prepare columns for each window
sub prepare_window_datasets {
	
	# Determine starting and ending points
	my $startingpoint = 0 - ($win * $number); 
		# default values will give startingpoint of -1000
	my $endingpoint = $win * $number; 
		# likewise default values will give endingpoint of 1000
	
	# Print collection statement
	print " Collecting ";
	if ($log) { 
		print "log2 ";
	}
	print "data from $startingpoint to $endingpoint at the ";
	if ($position == 3) {
		print "3' end"; 
	}
	elsif ($position == 4) {
		print "midpoint";
	}
	else {
		print "5' end";
	}
	print " in $win bp windows...\n";
	
	# Prepare and annotate the header names and metadata
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
		
		# the new name
		my $new_name = $start . '..' . $stop;
		
		# add new column
		my $new_index = $Data->add_column($new_name);
		
		# set the metadata
		$Data->metadata($new_index, 'start' , $start);
		$Data->metadata($new_index, 'stop' , $stop);
		$Data->metadata($new_index, 'window' , $win);
		$Data->metadata($new_index, 'log2' , $log);
		$Data->metadata($new_index, 'dataset' , $dataset);
		$Data->metadata($new_index, 'method' , $method);
		$Data->metadata($new_index, 'value' , $value_type);
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
			$Data->metadata($new_index, 'avoid', 1);
		}
	}
	
	return ($startingpoint, $endingpoint);
}



## Collect the nucleosome occupancy data
sub map_relative_data {
	
	# Add the columns for each window 
	# and calculate the relative starting and ending points
	my ($starting_point, $ending_point) = prepare_window_datasets();
	
	
	# determine long data collection for very large regions
	if ($ending_point - $starting_point > DATASET_HASH_LIMIT) {
		# This could potentially create performance issues where returned hashes   
		# for each feature or interval are too big for efficient data collection.
		# Better to collect data for individual windows using the long data method.
		$long_data = 1;
	}
	
	# Select the appropriate method for data collection
	if ($Data->feature_type eq 'coordinate' and not $long_data) {
		# mapping point data features using genome segments
		map_relative_data_for_regions($starting_point, $ending_point);
	}
	elsif ($Data->feature_type eq 'coordinate' and $long_data) {
		# mapping long data features using genome segments
		map_relative_long_data_for_regions();
	}
	elsif ($Data->feature_type eq 'named' and not $long_data) {
		# mapping point data features using named features
		map_relative_data_for_features($starting_point, $ending_point);
	}
	elsif ($Data->feature_type eq 'named' and $long_data) {
		# mapping long data features using named features
		map_relative_long_data_for_features();
	}
	else {
		die " Unable to identify columns with feature identifiers!\n" .
			" File must have Primary_ID or Name and Type, or Chromo, Start, Stop columns\n";
	}
}


sub map_relative_data_for_features {
	my ($starting_point, $ending_point) = @_;
	
	# Collect the data
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		my %regionscores = $row->get_position_scores(
			'ddb'         => $ddb,
			'dataset'     => $dataset,
			'start'       => $starting_point,
			'stop'        => $ending_point,
			'position'    => $position,
			'value'       => $value_type,
			'stranded'    => $strand_sense,
			'strand'      => $set_strand ? $row->strand : undef, 
			'avoid'       => $avoid,
		);
		
		# record the scores
		record_scores($row, \%regionscores);
	}
}



sub map_relative_long_data_for_features {
	
	# Collect the data
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		
		# get feature from the database
		my $feature = $row->feature;
		unless ($feature) {
			# record a null values
			for (my $c = $startcolumn; $c < $Data->number_columns; $c++) {
				$row->value($c, '.');
			}
			next;
		}
		
		# Collect the scores for each window
		collect_long_data_window_scores(
			$row,
			$feature->seq_id,
			$feature->start,
			$feature->end,
			$set_strand ? $row->strand : $feature->strand,
		);
		
		if ($avoid) {
			warn " avoid option is currently not supported with long data collection! Disabling!\n";
			$avoid = 0;
		}
	}
}

sub map_relative_data_for_regions {
	my ($starting_point, $ending_point) = @_;
	
	# Collect the data
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		
		# calculate new coordinates based on relative adjustments
			# this is a little tricky, because we're working with absolute 
			# coordinates but we want relative coordinates, so we must do 
			# the appropriate conversions
		my ($start, $stop, $region_start);
		
		if ($row->strand >= 0 and $position == 5) {
			# 5' end of forward strand
			$region_start = $row->start;
			$start = $row->start + $starting_point;
			$stop  = $row->start + $ending_point;
		}
		
		elsif ($row->strand == -1 and $position == 5) {
			# 5' end of reverse strand
			$region_start = $row->stop;
			$start = $row->stop - $ending_point;
			$stop  = $row->stop - $starting_point;
		}
		
		elsif ($row->strand >= 0 and $position == 3) {
			# 3' end of forward strand
			$region_start = $row->stop;
			$start = $row->stop + $starting_point;
			$stop  = $row->stop + $ending_point;
		}
		
		elsif ($row->strand == -1 and $position == 3) {
			# 3' end of reverse strand
			$region_start = $row->start;
			$start = $row->start - $ending_point;
			$stop  = $row->start - $starting_point;
		}
		
		elsif ($position == 4) {
			# midpoint regardless of strand
			$region_start = int( ( ($row->stop + $row->start) / 2) + 0.5);
			$start = $region_start + $starting_point;
			$stop  = $region_start + $ending_point;
		}
		
		else {
			# something happened
			die " programming error!? feature " . 
				" at data row $row->row_index\n";
		}
		
		# collect the region scores
		my %regionscores = $row->get_position_scores(
				'db'          => $ddb,
				'dataset'     => $dataset,
				'start'       => $start,
				'stop'        => $stop,
				'value'       => $value_type,
				'stranded'    => $strand_sense,
		);
		
		# convert the regions scores back into relative scores
		my %relative_scores;
		foreach my $position (keys %regionscores) {
			$relative_scores{ $position + $starting_point } = 
				$regionscores{$position};
		}
		
		# record the scores
		record_scores($row, \%relative_scores);
	}
}



sub map_relative_long_data_for_regions {
	
	# Collect the data
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		# Collect the scores for each window
		collect_long_data_window_scores(
			$row, $row->seq_id, $row->start, $row->end, $row->strand);
	}
}



sub record_scores {
	
	# row object and raw scores
	my ($row, $regionscores) = @_;
	
	# assign the scores to the windows in the region
	for (
		# we will process each window one at a time
		# proceed by the column index for each window
		my $column = $startcolumn; 
		$column < $Data->number_columns; 
		$column++
	) {
		# get start and stop
		my $start = $Data->metadata($column, 'start');
		my $stop = $Data->metadata($column, 'stop');
		
		# collect a score at each position in the window
		my @scores;
		for (my $n = $start; $n <= $stop; $n++) {
			# we will walk through the window one bp at a time
			# look for a score associated with the position
			push @scores, $regionscores->{$n} if exists $regionscores->{$n};
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
			elsif ($method eq 'rpm') {
				$winscore = ( sum(@scores) * 1000000 ) / $rpm_read_sum;
			}
			
			# deal with log2 scores
			if ($log) { 
				# put back in log2 space if necessary
				$winscore = log($winscore) / LOG2;
			}
		}
		else {
			# no scores
			# assign a "null" value
			$winscore = $method eq 'sum' ? 0 : '.';
		}
		
		# put the value into the data table
		$row->value($column, $winscore);
	}
}



## Collecting long data in windows
sub collect_long_data_window_scores {
	
	# passed row object and coordinates
	my (
		$row,
		$fchromo,
		$fstart,
		$fstop,
		$fstrand
	) = @_;

	# Translate the actual reference start position based on requested 
	# reference position and region strand
	my $reference;
	if ($fstrand >= 0 and $position == 5) {
		# 5' end of forward strand
		$reference = $fstart;
	}
	elsif ($fstrand == -1 and $position == 5) {
		# 5' end of reverse strand
		$reference = $fstop;
	}
	elsif ($fstrand >= 0 and $position == 3) {
		# 3' end of forward strand
		$reference = $fstop;
	}
	elsif ($fstrand == -1 and $position == 3) {
		# 3' end of reverse strand
		$reference = $fstart;
	}
	elsif ($position == 4) {
		# midpoint regardless of strand
		$reference = int( ( ($fstop + $fstart) / 2) + 0.5);
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
			'chromo'      => $fchromo,
			'start'       => $fstrand >= 0 ? 
								$reference + $Data->metadata($column, 'start') :
								$reference - $Data->metadata($column, 'start'),
			'stop'        => $fstrand >= 0 ? 
								$reference + $Data->metadata($column, 'stop') : 
								$reference - $Data->metdata($column, 'stop'),
			'strand'      => $fstrand,
			'method'      => $method,
			'value'       => $value_type,
			'stranded'    => $strand_sense,
			'log'         => $log,
		);
		$row->value($column, $score);
	}
}



## Interpolate the '.' values with the mean of the neighbors
sub go_interpolate_values {
	
	# determine counts
	my $lastwindow = $Data->number_columns - 1; 
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
				for (my $i = $col + 1; $col <= $lastwindow; $i++) {
					if ($row->value($i) ne '.') {
						$next_i = $i;
						last;
					}
				}
				next unless defined $next_i;
				
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



## Convert null values to proper zero values
sub null_to_zeros {
	# for those methods where we expect a true zero, sum and rpm
	# convert '.' null scores to zero
	
	# this wasn't done before because the null value makes the 
	# interpolation a lot easier without having to worry about zero
	
	# walk through each data line and then each window
	$Data->iterate( sub {
		my $row = shift;
		for (my $c = $startcolumn; $c < $Data->number_columns; $c++) {
			if ($row->value($c) eq '.') {
				$row->value($c, 0);
			}
		}
	} );
}



__END__

=head1 NAME

get_relative_data.pl

A script to collect data in bins around a relative position.

=head1 SYNOPSIS
 
get_relative_data.pl --in <in_filename> --out <out_filename> [--options]
  
  Options for existing files:
  --in <filename> 
  
  Options for new files:
  --db <name|file>
  --feature <type | type:source | alias>, ...
  
  Options for data collection:
  --ddb <name|file>
  --data <dataset_name | filename>
  --method [mean|median|min|max|stddev|sum|rpm]             (mean)
  --value [score|count|length]                              (score)
  --strand [all|sense|antisense]                            (all)
  --force_strand
  --avoid
  --long
  --log
  
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
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. The 
first row should be column headers. Bed files are acceptable, as are 
text files generated by other B<BioToolBox> scripts. Files may be 
gzipped compressed.

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

=item --method [mean|median|min|max|stddev|sum|rpm]

Specify the method of combining multiple values within each window. The mean, 
median, minimum, maximum, standard deviation, or sum of the values may be 
reported. The default value is mean for score and length values, or sum for 
count values.

=item --value [score|count|length]

Optionally specify the type of data value to collect from the dataset or 
data file. Three values are accepted: score, count, or length. Note that 
some data sources only support certain types of data values. Wig and 
BigWig files only support score and count; BigBed database features 
support count and length and optionally score; Bam files support basepair 
coverage (score), count (number of alignments), and length. The default 
value type is score. 

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

Indicate whether features of the same type should be avoided when 
calculating values in a window. Each window is checked for 
overlapping features of the same type; if the window does overlap 
another feature of the same type, no value is reported for the 
window. This option requires using named database features and must 
include a feature GFF type column. The default is false (return all 
values regardless of overlap).

=item --long

Indicate that the dataset from which scores are collected are 
long features (counting genomic annotation for example) and not point 
data (microarray data or sequence coverage). Normally long features are 
only recorded at their midpoint, leading to inaccurate representation at 
some windows. This option forces the program to collect data separately 
at each window, rather than once for each file feature or region and 
subsequently assigning scores to windows. Execution times may be 
longer than otherwise. Default is false.

=item --log

Dataset values are (not) in log2 space and should be treated 
accordingly. Output values will be in the same space. The default is 
false (nolog).

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

The default value that is collected is a dataset score (e.g. microarray 
values). However, other values may be collected, including 'count' or 
'length'. Use the --method argument to collect alternative values.

Stranded data may be collected. If the feature does not have an inherent 
strand, one may be specified to enforce stranded collection or a particular 
orientation. 

When features overlap, or the collection windows of one feature overlaps 
with another feature, then data may be ignored and not collected (--avoid).

The program writes out a tim data formatted text file. It will also 
generate a '*_summed.txt' file, in which each the mean value of all 
features for each window is generated and written as a data row. This 
summed data may be graphed using the biotoolbox script L<graph_profile.pl> 
or merged with other summed data sets for comparison.

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
from a Bam file in a stranded fashion. You also want to plot the profile.

  get_relative_data.pl --db annotation --feature gene --avoid --strand \
  sense --value count --method sum --data alignments.bam --out gene_tss
  
  graph_profile.pl --in gene_tss_summed.txt --min 0 --max 100
  
=back

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
