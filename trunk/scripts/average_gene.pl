#!/usr/bin/perl

# A script to collect values across the body of gene and assign them to an 'average'
# gene body. It collects values into percentage-based bins across the body.

use strict;
use Pod::Usage;
use Getopt::Long;
use Statistics::Lite qw(mean median sum min max);
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
);
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
	write_summary_data
);

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
	$database,
	$dataset,
	$feature,
	$method,
	$value_type,
	$strand,
	$bins,
	$extension,
	$extension_size,
	$min_length,
	$smooth,
	$sum,
	$log,
	$set_strand,
	$raw,
	$help
); # command line variables

## Command line options
GetOptions( 
	'in=s'        => \$infile, # input file
	'out=s'       => \$outfile, # name of outfile
	'db=s'        => \$database, # database name
	'data=s'      => \$dataset, # dataset name
	'feature=s'   => \$feature, # what type of feature to work with
	'method=s'    => \$method, # method for collecting the data
	'value=s'     => \$value_type, # type of data to collect
	'strand=s'    => \$strand, # indicate whether stranded data should be taken
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
unless ($outfile) {
	die " You must define an output filename !\n Use --help for more information\n";
}
unless ($feature or $infile) {
	# the feature must be defined
	# a tim data file used as an input file would also define one
	die " You must define a feature!\n Use --help for more information\n";
}
unless ($database or $infile) {
	# a database must be defined
	# a tim data file used as an input file would also define one
	die " You must define a database or input file!\n";
}
if ($strand) {
	unless (
		$strand eq 'sense' or 
		$strand eq 'antisense' or 
		$strand eq 'all'
	) {
		die " '$strand' is not recognized for strand\n Use --help for more information\n";
	}
} 
else {
	$strand = 'all'; # default is no strand information
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
		$method eq 'count'
	) {
		die " '$method' is not recognized for method\n Use --help for more information\n";
	}
	
	# convenience method, for backwards compatibility
	if ($method eq 'count') {
		$method = 'sum';
		$value_type = 'count';
	}
} 
else {
	# default is mean
	print " Using default method of 'mean'\n";
	$method = 'mean'; 
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



## Generate or load the dataset
my $main_data_ref; # reference to the main data structure hash
if ($infile) {
	# load a pre-existing data file
	print " Loading file $infile....\n";
	$main_data_ref = load_tim_data_file($infile);
	
	# check
	unless ($main_data_ref) {
		die " No file data loaded! Nothing done!\n";
	}
	
	# define database if necessary
	unless ($database) {
		$database = $main_data_ref->{'db'};
	}
	# define feature if necessary
	unless ($feature) {
		$feature = $main_data_ref->{'feature'};
	}
	# update program name
	unless ($main_data_ref->{'program'} eq $0) {
		$main_data_ref->{'program'} = $0;
	}
} 
else {
	# generate a new data file
	print " Generating new feature list....\n";
	$main_data_ref = get_new_feature_list( {
		'db'       => $database,
		'features' => $feature,
		'dubious'  => 1, # skip dubious genes
	} );
	
	# check
	unless ($main_data_ref) {
		die " No feature data generated! Nothing done!\n";
	}
	
	# define the program
	$main_data_ref->{'program'} = $0;
}
# the number of columns already in the data array
my $startcolumn = $main_data_ref->{'number_columns'}; 
# set the reference to the data table
my $data_table_ref = $main_data_ref->{'data_table'};



## Open database connectin
my $db = open_db_connection($database);



## Check for the dataset
$dataset = process_and_verify_dataset( {
	'db'      => $db,
	'dataset' => $dataset,
	'single'  => 1,
} );




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
my $timediff = sprintf "%.1f", (time - $start_time)/60;
print " Completed in $timediff minutes\n";



#### Subroutines #######



## Collect the binned data across the gene
sub collect_binned_data {
	
	## Prepare the metadata and header names
	my $binsize = (100/$bins); 
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
				_set_metadata($start . '%', $stop . '%', $binsize, '%');
			}
		}
	}
	
	# bins over the gene body
	for (my $i = 0; $i < $bins; $i++) { 
		my $start = ($i * $binsize );
		my $stop = ($i + 1) * $binsize;
		_set_metadata($start . '%', $stop . '%', $binsize, '%');
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
				_set_metadata($start . '%', $stop . '%', $binsize, '%');
			}
		}
	}
	
	## Identify columns for feature identification
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
	
	
	
	## Collect the data
	for my $row (1..$main_data_ref->{'last_row'}) {
		# walk through each feature
		
		# identify the feature first
		my $name = $data_table_ref->[$row][$name_index]; # name
		my $type = $data_table_ref->[$row][$type_index]; # class
		
		# start printing raw data
		if ($raw) { print RAWFILE join "\t", @{ $data_table_ref->[$row] } } 
		
		# pull gene from database
		my @genes = $db->features( 
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
		
		# the starting and ending points will be calculated from the number of
		# extensions, the binsize (multiply by 0.01 to get fraction), and the gene
		# length. No extensions should give just the length of the gene.
		# Remember that these end points will be relative to the feature, not
		# the chromosome!
		my $startingpoint = sprintf "%.0f", 
			( 0 - ($extension * $binsize * 0.01 * $length) +1 );
		my $endingpoint = sprintf "%.0f", 
			( $length + ($extension * $binsize * 0.01 * $length) );
		
		# collect the region scores
		my %regionscores = get_region_dataset_hash( {
					'db'       => $db,
					'dataset'  => $dataset,
					'name'     => $name,
					'type'     => $type,
					'start'    => $startingpoint,
					'stop'     => $endingpoint,
					'strand'   => $strand,
					'value'    => $value_type,
					'set_strand' => $set_strand ? 
						$data_table_ref->[$row][$strand_index] : undef,
		} );
		if ($raw) {
			foreach (sort {$a <=> $b} keys %regionscores) {
				print RAWFILE "\tdb: $regionscores{$_} \@ $_";
			}
		}
		
		# assign the scores to the bins in the region
		my $stepsize = sprintf "%.0f", $binsize * 0.01 * $length;
			# this converts the binsize from a percentage to actual bp length
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
				push @scores, $regionscores{$n} if exists $regionscores{$n};
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
				if ($method eq 'count' or $method eq 'sum') {
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
			if ($raw) { print RAWFILE "\tscore for $start: $window_score" }
		}
		
		# finish the raw line
		if ($raw) { print RAWFILE "\n" }
	}	
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
		'strand'    => $strand,
		'bin_size'  => $binsize . $unit,
		'method'    => $method,
		'value'     => $value_type,
	};
	if ($set_strand) {
		$main_data_ref->{$new_index}{'strand_implied'} = 1;
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
  --db <name|file.gff3>
  --feature [type, type:source, alias]
  --out <filename>
  --in <filename> 
  --data <dataset_name | filename>
  --method [mean|median|min|max|sum]
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

Specify the output file name. 

=item --in <filename>

Specify the filename of a data table containing the list of 
features to average. The file must be in the 'tim_data' format and specify 
a feature to use. If not specified, then a new feature list will be 
generated from the database.

=item --feature [type, type:source, alias]

Specify the type of feature from which to collect values. This is required 
for new feature tables. Two types of values may be passed: either a specific 
feature type present in the database, or an alias to one or more features. 
The feature may be specified as either type or type:source. Aliases are 
specified in the C<biotoolbox.cfg> file, and provide a shortcut to a 
list of one or more features. More than feature may be included as a 
comma-delimited list (no commas).
  
=item --data <dataset_name | filename>

Specify the name of the microarray data set for which you wish to 
collect data. If not specified, the data set may be 
chosen interactively from a presented list. When enumerating features, 
the features' type or type:source values should be indicated.
Alternatively, the name of a data file may be provided. Supported 
file types include BigWig (.bw), BigBed (.bb), or single-end Bam 
(.bam). The file may be local or remote.

=item --method [mean|median|min|max|sum]

Specify the method of collecting and combining the data into each bin. 
Three statistical methods for combining score values are allowed: mean, 
median, minimum, maximum, and sum. The defined method does not affect 
the interpolation or summary functions of the program, only initial data 
collection. The default method is 'mean'.

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
microarray tiling distance.

=item --strand [all|sense|antisense]

Specify whether stranded data should be collected. Three values are 
allowed: all datasets should be collected (default), only sense 
datasets, or only antisense datasets should be collected.

=item --sum

Indicate that the data should be averaged across all features at
each position, suitable for graphing. A separate text file will be
written with the suffix '_summed' with the averaged data

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

=item --help

This help text.

=back

=head1 DESCRIPTION

This program will collect data across a gene body into bins to generate an
'average' gene summary. It is used to determine if there is a spatial 
distribution preference for the dataset over gene bodies. The number of bins
may be specified as a command argument (default 10). Additionally, extra bins 
may be extended on either side of the gene (default 0 on either side). The bin 
size is determined as a percentage of gene length.

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



