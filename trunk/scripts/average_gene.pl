#!/usr/bin/perl

# A script to collect values across the body of gene and assign them to an 'average'
# gene body. It collects values into percentage-based bins across the body.

use strict;
use Pod::Usage;
use Getopt::Long;
use Statistics::Lite qw(mean median sum min max stddev);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	find_column_index
);
use tim_db_helper qw(
	open_db_connection
	process_and_verify_dataset
	get_new_feature_list
	get_region_dataset_hash
	check_dataset_for_rpm_support
);
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
	write_summary_data
);
my $VERSION = '1.7.0';

print "\n This script will collect binned values across genes to create an average gene\n\n";

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
	$value_type,
	$stranded,
	$bins,
	$extension,
	$extension_size,
	$min_length,
	$smooth,
	$sum,
	$log,
	$set_strand,
	$raw,
	$help,
	$print_version,
); # command line variables

## Command line options
GetOptions( 
	'in=s'        => \$infile, # input file
	'out=s'       => \$outfile, # name of outfile
	'db=s'       => \$main_database, # main or annotation database name
	'ddb=s'      => \$data_database, # data database
	'data=s'      => \$dataset, # dataset name
	'feature=s'   => \$feature, # what type of feature to work with
	'method=s'    => \$method, # method for collecting the data
	'value=s'     => \$value_type, # type of data to collect
	'strand=s'    => \$stranded, # indicate whether stranded data should be taken
	'bins=i'      => \$bins, # number of bins
	'ext=i'       => \$extension, # number of bins to extend beyond the feature
	'extsize=i'   => \$extension_size, # explicit size of extended bins
	'min=i'       => \$min_length, # minimum feature size
	'smooth!'     => \$smooth, # do not interpolate over missing values
	'sum'         => \$sum, # determine a final average for all the features
	'log!'        => \$log, # dataset is in log2 space
	'set_strand'  => \$set_strand, # enforce an artificial strand
	'raw'         => \$raw, # output raw data
	'help'        => \$help, # print the help
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
	print " Biotoolbox script average_gene.pl, version $VERSION\n\n";
	exit;
}



### Check for required values

unless ($main_database or $infile) {
	# a database must be defined
	# a tim data file used as an input file would also define one
	die " You must define a database or input file!\n";
}

unless ($feature or $infile) {
	# the feature must be defined
	# a tim data file used as an input file would also define one
	die " You must define a feature!\n Use --help for more information\n";
}

unless ($outfile) {
	die " You must define an output filename !\n Use --help for more information\n";
}

if ($stranded) {
	unless (
		$stranded eq 'sense' or 
		$stranded eq 'antisense' or 
		$stranded eq 'all'
	) {
		die " '$stranded' is not recognized for strand\n Use --help for more information\n";
	}
} 
else {
	$stranded = 'all'; # default is no strand information
}

if (defined $value_type) {
	# validate the requested value type
	unless (
		$value_type eq 'score' or
		$value_type eq 'count' or
		$value_type eq 'length'
	) {
		die " unknown value type '$value_type'!\n";
	}
}
else {
	# default value
	print " Collecting default data 'score' values\n";
	$value_type = 'score';
}

if (defined $method) {
	unless (
		$method eq 'mean' or 
		$method eq 'median' or 
		$method eq 'sum' or 
		$method eq 'min' or 
		$method eq 'max' or 
		$method eq 'stddev' or
		$method eq 'rpm'
	) {
		die " '$method' is not recognized for method\n Use --help for more information\n";
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

# assign default values
unless (defined $bins) {
	# dividing gene into 10 (10%) bins seems reasonable to me
	$bins = 10;
} 

unless (defined $extension) {
	# default is no extension
	$extension = 0;
}

unless (defined $log) {
	# default is that data is not in log 2
	$log = 0;
}

unless (defined $smooth) {
	# default is to not include smoothing
	$smooth = 0;
}




### Generate or load the dataset
my $main_data_ref; # a reference to the main data structure
if ($infile) {
	# load the gene dataset from existing file
	print " Loading feature set from file $infile....\n";
	$main_data_ref = load_tim_data_file($infile);
	
	# update program name
	unless ($main_data_ref->{'program'} eq $0) {
		$main_data_ref->{'program'} = $0;
	}
	
	# define database if necessary
	unless ($main_database) {
		$main_database = $main_data_ref->{'db'};
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

# the number of columns already in the data array
my $startcolumn = $main_data_ref->{'number_columns'}; 
# set the reference to the data table
my $data_table_ref = $main_data_ref->{'data_table'};





### Open database connection
unless ($main_database) {
	# define the database 
	if ( $main_data_ref->{'db'} ) {
		# from the input file if possible
		$main_database = $main_data_ref->{'db'};
	}
	elsif ($dataset) {
		# or use the dataset if possible
		$main_database = $dataset;
	}
	else {
		die " You must define a database or dataset!\n" . 
			" Use --help for more information\n";
	}
}
my $mdb = open_db_connection($main_database) or
	die " unable to open database '$main_database'!\n";

# Open data database
my $ddb;
if (defined $data_database) {
	# specifically defined a data database
	$ddb = open_db_connection($data_database) or 
		die "unable to establish data database connection to $data_database!\n";
}
else {
	# reuse the main database connection
	$ddb = $mdb;
}






## Check for the dataset
$dataset = process_and_verify_dataset( {
	'db'      => $ddb,
	'dataset' => $dataset,
	'single'  => 1,
} );

# Check the RPM method if necessary
my $rpm_read_sum;
if ($method eq 'rpm') {
	print " Checking RPM support for dataset '$dataset'...\n";
	$rpm_read_sum = check_dataset_for_rpm_support($dataset, $ddb);
	if ($rpm_read_sum) {
		printf " %s total features\n", format_with_commas($rpm_read_sum);
	}
	else {
		die " RPM method not supported! Try something else\n";
	}
}






## Collect the binned data
my $start_time = time;
if ($raw) {
	open RAWFILE, ">$outfile\_raw.txt";
}
print " Collecting $method data from $dataset in " . 
	($bins + 2 * $extension) . " bins....\n"; 
collect_binned_data();
close RAWFILE;



## Interpolate values
if ($smooth) {
	print " Smoothing data by interpolation....\n";
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

# write main output
my $written_file = write_tim_data_file( {
	# we will write a tim data file
	# appropriate extensions and compression should be taken care of
	'data'     => $main_data_ref,
	'filename' => $outfile,
} );
if ($written_file) {
	print " Wrote data file '$written_file'\n";
}
else {
	print " unable to write data file!\n";
}

# print completion
printf " Completed in %.1f minutes\n", (time - $start_time)/60;



#### Subroutines #######

## generate a feature dataset if one was not loaded
sub generate_a_new_feature_dataset {
	# a subroutine to generate a new feature dataset
	
	$main_data_ref = get_new_feature_list( {
			'db'       => $main_database,
			'features' => $feature,
	} );
	
	# set the current program
	$main_data_ref->{'program'} = $0;
}



## Collect the binned data across the gene
sub collect_binned_data {
	
	## Prepare the metadata and header names
	my $binsize = (100/$bins); 
	prepare_bins($binsize);
	
	
	## Identify columns for feature identification
	my $name = find_column_index($main_data_ref, '^name|id');
	
	my $type = find_column_index($main_data_ref, '^type|class');
	
	my $chromo = find_column_index($main_data_ref, '^chr|seq|ref|ref.?seq');
	
	my $start = find_column_index($main_data_ref, '^start|position');
	
	my $stop = find_column_index($main_data_ref, '^stop|end');
	
	my $strand = find_column_index($main_data_ref, '^strand');
	
	
	
	## Select the appropriate method for data collection
	if (
		defined $name and
		defined $type
	) {
		# using named features
		collect_binned_data_for_features(
			$binsize, $name, $type, $strand);
	}
	
	elsif (
		defined $start and
		defined $stop  and
		defined $chromo
	) {
		# using genome segments
		collect_binned_data_for_regions(
			$binsize, $chromo, $start, $stop, $strand);
	}
	
	else {
		die " Unable to identify columns with feature identifiers!\n" .
			" File must have Name and Type, or Chromo, Start, Stop columns\n";
	}
}	



sub collect_binned_data_for_features {	
	
	## Passed data
	my ($binsize, $name_index, $type_index, $strand_index) = @_;
	
	
	## Collect the data
	for my $row (1..$main_data_ref->{'last_row'}) {
		# walk through each feature
		
		# identify the feature first
		my $name = $data_table_ref->[$row][$name_index]; # name
		my $type = $data_table_ref->[$row][$type_index]; # class
		
		# start printing raw data
		if ($raw) { print RAWFILE join "\t", @{ $data_table_ref->[$row] } } 
		
		# pull gene from database
		my @genes = $mdb->features( 
				-name  => $name,
				-type => $type,
		);
		if (scalar @genes > 1) {
			# there should only be one feature found
			# if more, there's redundant or duplicated data in the db
			# warn the user, this should be fixed
			warn " Found more than one feature of '$type => $name' in " .
				"the database!\n Using the first feature only!\n";
		}
		my $gene = shift @genes; 
		unless ($gene) {
			die " unable to establish db region for $type feature $name!\n";
		}
		
		# define the starting and ending points based on gene length
		my $length = $gene->length;
		
		# check the length
		if (defined $min_length and $length < $min_length) {
			# this feature is too short to divided into bins
			# we will skip this feature
			
			# but first, put in null values
			for my $column ($startcolumn..($main_data_ref->{'number_columns'} - 1) ) {
				$data_table_ref->[$row][$column] = '.';
			}
			
			# move on to next feature
			next;
		}
		
		# calculate the actual extension in bp
		# this is determined from the extension_size or the feature length
		# and multiplied by the number of extensions requested
		my $extra;
		if ($extension) {
			if ($extension_size) {
				# extension is specific bp in size
				$extra = $extension_size * $extension;
			}
			else {
				# extension is dependent on feature length
				$extra = int( ($extension * $binsize * 0.01 * $length) + 0.5);
			}
		}
		
		
		# collect the region scores
		my %regionscores = get_region_dataset_hash( {
					'db'        => $mdb,
					'ddb'       => $ddb,
					'dataset'   => $dataset,
					'value'     => $value_type,
					'name'      => $name,
					'type'      => $type,
					'extend'    => $extra,
					'stranded'  => $stranded,
					'strand'    => $set_strand ? 
						$data_table_ref->[$row][$strand_index] : undef,
		} );
		if ($raw) {
			foreach (sort {$a <=> $b} keys %regionscores) {
				print RAWFILE "\tdb: $regionscores{$_} \@ $_";
			}
		}
		
		# record the scores for each bin
		record_the_bin_values($row, $length, \%regionscores);
	}	
}



sub collect_binned_data_for_regions {
	
	## Passed data
	my ($binsize, $chromo_index, $start_index, $stop_index, $strand_index) = @_;
	
	
	## Collect the data
	for my $row (1..$main_data_ref->{'last_row'}) {
		# walk through each feature
		
		# identify the coordinates
		my $chromo = $data_table_ref->[$row][$chromo_index]; 
		my $start  = $data_table_ref->[$row][$start_index]; 
		my $stop   = $data_table_ref->[$row][$stop_index]; 
		my $strand = $data_table_ref->[$row][$strand_index] || 1; # default +1
		
		# start printing raw data
		if ($raw) { print RAWFILE join "\t", @{ $data_table_ref->[$row] } } 
		
		# determine the segment length
		my $length = $stop - $start + 1;
		
		# check the length
		if (defined $min_length and $length < $min_length) {
			# this feature is too short to divided into bins
			# we will skip this feature
			
			# but first, put in null values
			for my $column ($startcolumn..($main_data_ref->{'number_columns'} - 1) ) {
				$data_table_ref->[$row][$column] = '.';
			}
			
			# move on to next feature
			next;
		}
		
		# the starting and ending points will be calculated from the number of
		# extensions, the binsize (multiply by 0.01 to get fraction), and the gene
		# length. No extensions should give just the length of the gene.
		my ($startingpoint, $endingpoint);
		if ($extension_size) {
			# extension is specific bp in size
			my $extra = $extension_size * $binsize;
			$startingpoint = int( $start - $extra + 0.5 );
			$endingpoint   = int( $stop + $extra + 0.5 );
		}
		else {
			# extension is dependent on feature length
			my $extra = $extension * $binsize * 0.01 * $length;
			$startingpoint = int( $start - $extra + 0.5);
			$endingpoint   = int( $stop + $extra + 0.5);
		}
		
		# collect the region scores
		my %regionscores = get_region_dataset_hash( {
					'db'       => $ddb,
					'dataset'  => $dataset,
					'value'    => $value_type,
					'chromo'   => $chromo,
					'start'    => $startingpoint,
					'stop'     => $endingpoint,
					'absolute' => 1,
					'stranded' => $stranded,
					'strand'   => $strand,
		} );
		if ($raw) {
			foreach (sort {$a <=> $b} keys %regionscores) {
				print RAWFILE "\tdb: $regionscores{$_} \@ $_";
			}
		}
		
		# convert absolute scores to relative scores
		# must be done here because we have the original coordinates
		my %relative_scores;
		if ($strand >= 0) {
			# forward strand
			foreach my $pos (keys %regionscores) {
				$relative_scores{ $pos - $start } = $regionscores{$pos};
			}
		}
		else {
			# reverse strand
			foreach my $pos (keys %regionscores) {
				$relative_scores{ $stop - $pos } = $regionscores{$pos};
			}
		}
		
		# record the scores for each bin
		record_the_bin_values($row, $length, \%relative_scores);
		
	}
}



sub record_the_bin_values {
	
	# get the passed values
	my ($row, $length, $regionscores) = @_;
	
	
	# assign the scores to the bins in the region
	for my $column ($startcolumn..($main_data_ref->{'number_columns'} - 1) ) {
		# we will step through each data column, representing each window (bin)
		# across the feature's region
		# any scores within this window will be collected and the mean 
		# value reported
		
		# convert the window start and stop coordinates (as percentages) to
		# actual bp
		# this depends on whether the binsize is explicitly defined in bp or
		# is a fraction of the feature length
		my ($start, $stop);
		if ($main_data_ref->{$column}{'bin_size'} =~ /bp$/) {
			# the bin size is explicitly defined
			
			# the start and stop points are relative to either the feature
			# start (always 0) or the end (the feature length), depending
			# upon whether the 5' or 3' end of the feature
			
			# determine this by the sign of the start position
			if ($main_data_ref->{$column}{'start'} < 0) {
				# the start position is less than 0, implying the 5' end
				# the reference position will be the feature start, or 0
				$start = $main_data_ref->{$column}{'start'};
				$stop = $main_data_ref->{$column}{'stop'};
			}
			else {
				# the start position is greather than 0, implying the 3' end
				# the reference position will be the feature end, or length
				$start = $main_data_ref->{$column}{'start'} + $length;
				$stop = $main_data_ref->{$column}{'stop'} + $length;
			}
		}
		else {
			# otherwise the bin size is based on feature length
			$start = sprintf "%.0f", ( 
				$main_data_ref->{$column}{'start'} * 0.01 * $length) + 1;
			$stop = sprintf "%.0f", ( 
				$main_data_ref->{$column}{'stop'} * 0.01 * $length);
		}
		
		# collect the scores for this window
		my @scores;
		for (my $n = $start; $n <= $stop; $n++) {
			# we will walk through each bp in the window looking for a score
			push @scores, $regionscores->{$n} if exists $regionscores->{$n};
		}
		
		
		# calculate the value
		my $window_score;
		if (@scores) {
			# we have values in the window
			
			# combine the scores according to the specified method
			if ($method eq 'sum') {
				# either the count or the sum methods require that the 
				# scores be summed
				$window_score = sum(@scores);
				# raw output
				if ($raw) { print RAWFILE "\tfound: " . join(",", @scores) }
				
			}
			
			else {
				# method of mean or median to combine the scores
			
				# convert from log2 if necessary
				if ($log) {
					@scores = map { 2 ** $_ } @scores;
				}
				
				# calculate the score appropriately
				if ($method eq 'mean') {
					$window_score = mean(@scores); 
				}
				elsif ($method eq 'median') {
					$window_score = median(@scores); 
				}
				elsif ($method eq 'min') {
					$window_score = min(@scores); 
				}
				elsif ($method eq 'max') {
					$window_score = max(@scores); 
				}
				elsif ($method eq 'sum') {
					$window_score = sum(@scores);
				}
				elsif ($method eq 'stddev') {
					$window_score = stddev(@scores);
				}
				elsif ($method eq 'rpm') {
					$window_score = ( sum(@scores) * 1000000 ) / $rpm_read_sum;
				}
				
				# raw output
				if ($raw) { print RAWFILE "\tfound: " . join(",", @scores) }
				
				# convert back to log if necessary
				if ($log) {
					if ($window_score != 0) {
						$window_score = log($window_score) / log(2);
					}
					else {
						$window_score = '.';
					}
				}
			}
		} 
		else {
			# no values in this window
			if ($method eq 'sum' or $method eq 'rpm') {
				# score gets 0
				$window_score = 0;
			}
			else {
				# no score gets a null symbol
				$window_score = '.'; 
			}
		}
		
		
		# record the value
		$data_table_ref->[$row][$column] = $window_score;
		if ($raw) { print RAWFILE "\tscore for col $column start $start: $window_score" }
	}
	
	# finish the raw line
	if ($raw) { print RAWFILE "\n" }

}





## Interpolate the '.' values with the mean of the neighbors
sub go_interpolate_values {
	
	# determine the index of the second to last position in the data table array
	my $lastcolumn = $main_data_ref->{'number_columns'} - 2; 
	
	for my $row (1..$main_data_ref->{'last_row'}) {
		# walk through each feature row
		
		for my $column (($startcolumn + 1)..$lastcolumn) {
			# walk through the datasets for a single feature
			# skipping the very first and last bins
			
			# we will check each data element for the non-value '.'
			if ($data_table_ref->[$row][$column] eq '.') {
				# we have a non-value
				# we will interpolate this value from the neighbors
				# check the neighbors for real values
				# and then calculate from the values
				
				# we can tolerate several consecutive non-values
				
				# only one non-value ($column) flanked by values
				if (
					$data_table_ref->[$row][$column - 1] ne '.' and
					$data_table_ref->[$row][$column + 1] ne '.'
				) {
					# take simple average of the two neighbors
					my $new_value;
					if ($log) {
						# log values
						$new_value = mean(
							( 2 ** $data_table_ref->[$row][$column - 1] ), 
							( 2 ** $data_table_ref->[$row][$column + 1] )
						);
						$new_value = log($new_value) / log(2);
					}
					else {
						# non-log values
						$new_value = mean(
							$data_table_ref->[$row][$column - 1], 
							$data_table_ref->[$row][$column + 1]
						);
					}
					$data_table_ref->[$row][$column] = $new_value;
				}
				
				# two non-values ($column and $column+1) flanked by values
				elsif (
					$data_table_ref->[$row][$column - 1] ne '.' and 
					$data_table_ref->[$row][$column + 1] eq '.' and 
					( 
						$data_table_ref->[$row][$column + 2] ne '.' or 
						$data_table_ref->[$row][$column + 2] ne '' # true null
					)
				) {
					# interpolation of two intervening non-values from neighbors
					
					# identify the begin and end values
					my $begin_value = 
						$data_table_ref->[$row][$column - 1];
					my $last_value = 
						$data_table_ref->[$row][$column + 2];
					
					if ($log) {
						# log values
						
						# adjust the begin and end values
						my $last_value = 2 ** $last_value;
						my $begin_value = 2 ** $begin_value;
						
						# calculate the one-third value of the delta
						my $third = ( $last_value - $begin_value ) / 3;
						
						# calculate the intervening non-values
						$data_table_ref->[$row][$column] = 
							log( $begin_value + $third ) / log(2);
						$data_table_ref->[$row][$column+1] = 
							log( $begin_value + (2 * $third) ) / log(2);
					}
					else {
						# non-log values
						
						# calculate the one-third value of the delta
						my $third = ( $last_value - $begin_value ) / 3;
						
						# calculate the intervening non-values
						$data_table_ref->[$row][$column] = 
							$begin_value + $third;
						$data_table_ref->[$row][$column+1] = 
							$begin_value + (2 * $third);
					}
				}
				
				# three non-values ($column, $column+1, and $column+2) flanked 
				# by values
				elsif ($data_table_ref->[$row][$column - 1] ne '.' and 
					$data_table_ref->[$row][$column + 1] eq '.' and 
					$data_table_ref->[$row][$column + 2] eq '.' and 
					(
						$data_table_ref->[$row][$column + 3] ne '.' or 
						$data_table_ref->[$row][$column + 3] ne '' # true null
					)
				) {
					# interpolation of three intervening non-values from neighbors
					
					# identify the begin and end values
					my $begin_value = 
						$data_table_ref->[$row][$column - 1];
					my $last_value = 
						$data_table_ref->[$row][$column + 3];
					
					if ($log) {
						# log values
						
						# adjust the begin and end values
						my $last_value = 2 ** $last_value;
						my $begin_value = 2 ** $begin_value;
						
						# calculate the one-third value of the delta
						my $fourth = ( $last_value - $begin_value ) / 4;
						
						# calculate the intervening non-values
						$data_table_ref->[$row][$column] = 
							log( $begin_value + $fourth ) / log(2);
						$data_table_ref->[$row][$column+1] = 
							log( $begin_value + (2 * $fourth) ) / log(2);
						$data_table_ref->[$row][$column+2] = 
							log( $begin_value + (3 * $fourth) ) / log(2);
					}
					else {
						# non-log values
						
						# calculate the one-third value of the delta
						my $fourth = ( $last_value - $begin_value ) / 4;
						
						# calculate the intervening non-values
						$data_table_ref->[$row][$column] = 
							$begin_value + $fourth;
						$data_table_ref->[$row][$column+1] = 
							$begin_value + (2 * $fourth);
						$data_table_ref->[$row][$column+2] = 
							$begin_value + (3 * $fourth);
					}
				}
				
				# how much longer can I keep this up?
			} # end check
		} # end column
	} # end row
}


### Prepare all of the bin columns and their metadata
sub prepare_bins {
	
	my $binsize = shift;
	
	# the size of the bin in percentage units, default would be 10%
	# each bin will be titled the starting and ending point for that bin in 
	# percentage units
	# for example, -20..-10,-10..0,0..10,10..20
	
	# if $extension is defined, then it will add the appropriate flanking bins,
	# otherwise it should skip them 
	
	# bin(s) on 5' flank
	if ($extension) {
		# 5' bins are requested
		if ($extension_size) {
			# extended bins should be of specific bp size
			for (my $i = $extension; $i > 0; $i--) { 
				my $start = 0 - ($extension_size * $i);
				my $stop = 0 - ($extension_size * ($i - 1));
				_set_metadata($start, $stop, $extension_size, 'bp');
			}
		}
		else {
			# extended bin size will be based on feature length
			for (my $i = $extension; $i > 0; $i--) { 
				my $start = 0 - ($binsize * $i);
				my $stop = 0 - ($binsize * ($i - 1));
				_set_metadata($start, $stop, $binsize, '%');
			}
		}
	}
	
	# bins over the gene body
	for (my $i = 0; $i < $bins; $i++) { 
		my $start = ($i * $binsize );
		my $stop = ($i + 1) * $binsize;
		_set_metadata($start, $stop, $binsize, '%');
	}
	
	# bin(s) on 3' flank
	if ($extension) {
		# 5' bins are requested
		if ($extension_size) {
			# extended bins should be of specific bp size
			for (my $i = 0; $i < $extension; $i++) { 
				my $start = ($extension_size * $i);
				my $stop = ($extension_size * ($i + 1));
				_set_metadata($start, $stop, $extension_size, 'bp');
			}
		}
		else {
			# extended bin size will be based on feature length
			for (my $i = 0; $i < $extension; $i++) { 
				my $start = 100 + ($binsize * $i);
				my $stop = 100 + ($binsize * ($i + 1));
				_set_metadata($start, $stop, $binsize, '%');
			}
		}
	}
}




### Set the metadata for a new data table column (dataset)
sub _set_metadata {
	# the start and stop positions are passed
	my ($start, $stop, $binsize, $unit) = @_;
	
	# set new index
	my $new_index = $main_data_ref->{'number_columns'};
	
	# set new name
	my $name = $start . '..' . $stop . $unit;
	
	# set the metadata using passed and global variables
	$main_data_ref->{$new_index} = {
		'name'      => $name,
		'index'     => $new_index,
		'start'     => $start,
		'stop'      => $stop,
		'log2'      => $log,
		'dataset'   => $dataset,
		'strand'    => $stranded,
		'bin_size'  => $binsize . $unit,
		'method'    => $method,
		'value'     => $value_type,
	};
	if ($set_strand) {
		$main_data_ref->{$new_index}{'strand_implied'} = 1;
	}
	if ($data_database) {
		$main_data_ref->{$new_index}{'db'} = $data_database;
	}
	
	# set the column header
	$data_table_ref->[0][$new_index] = $name;
	
	# increase the column numbers count
	$main_data_ref->{'number_columns'} += 1;
}





__END__

=head1 NAME

average_gene.pl

=head1 SYNOPSIS
 
 average_gene.pl --db <database> --feature <text> --out <file> [--options]
 average_gene.pl --in <file> --out <file> [--options]
  
  Options:
  --in <filename> 
  --out <filename>
  --db <name | filename>
  --ddb <name | filename>
  --feature <type | type:source | alias>,...
  --data <dataset_name | filename>
  --method [mean|median|stddev|min|max|range|sum|rpm|rpkm]
  --value [score|count|length]
  --bins <integer>
  --ext <integer>
  --extsize <integer>
  --min <integer>
  --strand [all|sense|antisense]
  --sum
  --smooth
  --set_strand
  --(no)log
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
text files generated with this program.

=item --out <filename>

Specify the output file name. 

=item --db <name | filename>

Specify the name of a BioPerl database from which to obtain the 
annotation, chromosomal information, and/or data. Typically a 
Bio::DB::SeqFeature::Store database schema is used, either from a 
relational database, SQLite file, or a single GFF3 file to be loaded 
into memory. Alternatively, a BigWigSet directory, or a single BigWig, 
BigBed, or Bam file may be specified. 

A database is required for generating new files. When generating a new 
genome interval file, a bigFile or Bam file listed as a data source 
will be adopted as the database. 

For input files, the database name may be obtained from the file 
metadata. A different database may be specified from that listed in 
the metadata when a different source is desired. 

=item --ddb <name | filename>

Optionally specify the name of an alternate data database from which 
the data should be collected, separate from the primary annotation 
database. The same options apply as to the --db option. 

=item --feature <type | type:source | alias>,...

Specify the type of feature from which to collect values. This is required 
only for new feature tables. Three types of values may be passed: the 
feature type, feature type and source expressed as 'type:source', or an 
alias to one or more feature types. Aliases are specified in the 
C<biotoolbox.cfg> file and provide a shortcut to a list of one or more 
database features. More than one feature may be included as a 
comma-delimited list (no spaces). 
  
=item --data <dataset_name | filename>

Provide the name of the dataset to collect the values. If no 
dataset is specified on the command line, then the program will 
interactively present a list of datasets from the database to select. 

The dataset may be a feature type in a BioPerl Bio::DB::SeqFeature::Store 
or Bio::DB::BigWigSet database. Provide either the feature type or 
type:source. The feature may point to another data file whose path is 
stored in the feature's attribute tag (for example a binary 
Bio::Graphics::Wiggle .wib file, a bigWig file, or Bam file), or the 
features' scores may be used in data collection.

Alternatively, the dataset may be a database file, including bigWig (.bw), 
bigBed (.bb), or Bam alignment (.bam) files. The files may be local or 
remote (specified with a http: or ftp: prefix).

=item --method [mean|median|stddev|min|max|range|sum|rpm|rpkm]

Specify the method for combining all of the dataset values within the 
genomic region of the feature. Accepted values include:
  
  - mean        (default)
  - median
  - sum
  - stddev      Standard deviation of the population (within the region)
  - min
  - max
  - rpm         Reads Per Million mapped, for Bam and BigBed only

=item --value [score|count|length]

Optionally specify the type of data value to collect from the dataset or 
data file. Three values are accepted: score, count, or length. The default 
value type is score. Note that some data sources only support certain 
types of data values. Wig and BigWig files only support score and count; 
BigBed and database features support count and length and optionally 
score; Bam files support basepair coverage (score), count (number of 
alignments), and length.

=item --bins <integer>

Specify the number of bins that will be generated over the length 
of the feature. The size of the feature is a percentage of the 
feature length. The default number is 10, which results in bins of 
size equal to 10% of the feature length. 

=item --ext <integer>

Specify the number of extended bins on either side of the feature. 
The bins are of the same size as determined by the feature 
length and the --bins value. The default is 0. 

=item --extsize <integer>

Specify the exact bin size in bp of the extended bins rather than
using a percentage of feature of length.

=item --min <integer>

Specify the minimum feature size to be averaged. Features with a
length below this value will not be skipped (all bins will have
null values). This is to avoid having bin sizes below the average 
microarray tiling distance. The default is undefined (no limit).

=item --strand [all|sense|antisense]

Specify whether stranded data should be collected. Three values are 
allowed: all datasets should be collected (default), only sense 
datasets, or only antisense datasets should be collected.

=item --sum

Indicate that the data should be averaged across all features at
each position, suitable for graphing. A separate text file will be
written with the suffix '_summed' with the averaged data. The default 
is false.

=item --smooth

Indicate that windows without values should (not) be interpolated
from neighboring values. The default is false.

=item --set_strand

For features that are not inherently stranded (strand value of 0), 
impose an artificial strand for each feature (1 or -1). This will 
have the effect of enforcing a relative orientation for each feature, 
or to collected stranded data. This requires the presence a 
column in the input data file with a name of "strand". Hence, it 
will not work with newly generated datasets, but only with input 
data files. Default is false.

=item --(no)log

Dataset values are (not) in log2 space and should be treated 
accordingly. Output values will be in the same space.

=item --version

Print the version number.

=item --help

This help text.

=back

=head1 DESCRIPTION

This program will collect data across a gene or feature body into numerous 
percentile bins to. It is used to determine if there is a spatial 
distribution preference for the dataset over gene bodies. The number 
of bins may be specified as a command argument (default 10). Additionally, 
extra bins may be extended on either side of the gene (default 0 on either 
side). The bin size is determined as a percentage of gene length.

The program writes out a tim data formatted text file. It will also 
optionally generate a "summed" or average profile for all of the features. 

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



