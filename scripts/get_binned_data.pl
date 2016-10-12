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
use Bio::ToolBox::GeneTools qw(:transcript);
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};
my $VERSION = '1.50';

print "\n This script will collect binned values across features\n\n";

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
	$subfeature,
	$method,
	$stranded,
	$bins,
	$extension,
	$extension_size,
	$min_length,
	$long_data,
	$smooth,
	$sum,
	$set_strand,
	$gz,
	$cpu,
	$help,
	$print_version,
); # command line variables

## Command line options
GetOptions( 
	'in=s'        => \$infile, # input file
	'out=s'       => \$outfile, # name of outfile
	'db=s'        => \$main_database, # main or annotation database name
	'ddb=s'       => \$data_database, # data database
	'data=s'      => \$dataset, # dataset name
	'feature=s'   => \$feature, # what type of feature to work with
	'exons!'      => \$subfeature, # indicate to restrict to subfeatures
	'method=s'    => \$method, # method for collecting the data
	'strand=s'    => \$stranded, # indicate whether stranded data should be taken
	'bins=i'      => \$bins, # number of bins
	'ext=i'       => \$extension, # number of bins to extend beyond the feature
	'extsize=i'   => \$extension_size, # explicit size of extended bins
	'min=i'       => \$min_length, # minimum feature size
	'long!'       => \$long_data, # collecting long data features
	'smooth!'     => \$smooth, # do not interpolate over missing values
	'sum'         => \$sum, # determine a final average for all the features
	'force_strand|set_strand'  => \$set_strand, # enforce an artificial strand
				# force_strand is preferred option, but respect the old option
	'gz!'         => \$gz, # compress the output file
	'cpu=i'       => \$cpu, # number of execution threads
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
check_defaults();
my $start_time = time;
my $collapsed = 0; # global value to indicate transcripts are collapsed



### Generate or load the dataset
my $Data;
if ($infile) {
	$Data = Bio::ToolBox::Data->new(file => $infile, parse => 1) or 
		die " unable to load input file '$infile'\n";
	printf " Loaded %s features from $infile.\n", format_with_commas( $Data->last_row );
	
	# update main database as necessary
	if ($main_database) {
		if ($main_database ne $Data->database) {
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



# Open data database
my $ddb;
if (defined $data_database) {
	# specifically defined a data database
	$ddb = open_db_connection($data_database) or 
		die "unable to establish data database connection to $data_database!\n";
}

# Check for the dataset
$dataset = verify_or_request_feature_types(
	'db'      => $ddb || $Data->database,
	'feature' => $dataset,
	'prompt'  => " Enter the number of the feature or dataset from which to" . 
					" collect data   ",
	'single'  => 1,
);
unless ($dataset) {
	die " No verifiable dataset provided. Check your file path, database, or dataset.\n";
} 




## Collect the relative data
printf " Collecting $method data from $dataset in %s bins....\n",
	($bins + 2 * $extension); 

# check whether it is worth doing parallel execution
if ($cpu > 1) {
	while ($cpu > 1 and $Data->last_row/$cpu < 100) {
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



## Generate summed data - 
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



## Write main output
unless ($outfile) {
	$outfile = $Data->path . $Data->basename;
}
my $written_file = $Data->save(
	'filename' => $outfile,
	'gz'       => $gz,
);
if ($written_file) {
	print " Wrote data file '$written_file'\n";
}
else {
	print " unable to write data file!\n";
}
printf " Completed in %.1f minutes\n", (time - $start_time)/60;
# done



#### Subroutines #######

sub check_defaults {
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

	unless ($outfile or $infile) {
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

	if (defined $method) {
		unless (
			$method eq 'mean' or 
			$method eq 'median' or 
			$method eq 'sum' or 
			$method eq 'min' or 
			$method eq 'max' or 
			$method eq 'stddev' or
			$method =~ /^\w?count/
		) {
			die " '$method' is not recognized for method\n Use --help for more information\n";
		}
	
	} 
	else {
		$method = 'mean';
	}
	
	if ($subfeature and $long_data) {
		die " Long data collection is incompatible with exon subfeatures\n Use --help for more information\n";
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

	unless (defined $smooth) {
		# default is to not include smoothing
		$smooth = 0;
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



## Parallel execution for efficiency
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
		
		# Prepare the metadata and header names
		my $binsize = (100/$bins); 
		prepare_bins($binsize);
		
		# Collect the data
		if ($long_data) {
			collect_binned_long_data($binsize);
		}
		else {
			collect_binned_data($binsize);
		}
		
		# Interpolate values
		if ($smooth) {
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
	
	# Prepare the metadata and header names
	my $binsize = (100/$bins); 
	prepare_bins($binsize);
	
	# Collect the data
	if ($long_data) {
		collect_binned_long_data($binsize);
	}
	else {
		collect_binned_data($binsize);
	}
	
	# Interpolate values
	if ($smooth) {
		print " Smoothing data by interpolation....\n";
		go_interpolate_values();
	}
}



sub collect_binned_data {	
	my $binsize = shift;
	
	# collapse transcripts if needed
	if ($subfeature and not $collapsed) {
		$collapsed = collapse_all_features();
	}
	my $length_i = $Data->find_column('Merged_Transcript_Length');
	
	
	## Collect the data
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		
		# identify the feature first
		my $feature = $row->feature || $row;
		
		# define the starting and ending points based on gene length
		my $length = $length_i ? $row->value($length_i) : 
			$feature->length;
		
		# check the length
		if (defined $min_length and $length < $min_length) {
			# this feature is too short to divided into bins
			# we will skip this feature
			
			# but first, put in null values
			for my $c ($startcolumn..($Data->last_column) ) {
				$row->value($c, '.');
			}
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
		
		my $regionscores = $row->get_region_position_scores(
			'ddb'       => $ddb,
			'dataset'   => $dataset,
			'method'    => $method,
			'extend'    => $extra,
			'stranded'  => $stranded,
			'strand'    => $set_strand ? $row->strand : $feature->strand,
			'exon'      => $subfeature,
			'length'    => $length,
		);
		# record the scores for each bin
		record_the_bin_values($row, $length, $regionscores);
	}	
}


sub collect_binned_long_data {	
	my $binsize = shift;
	
	## Collect the data
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		
		# identify the feature or use the row
		my $feature = $row->feature || $row;
		my $fstart = $feature->start;
		my $fstop = $feature->end;
		my $strand = $set_strand ? $row->strand : $feature->strand;
		my $length = $feature->length;
		
		# collect the scores to the bins in the region
		for my $column ($startcolumn..($Data->last_column) ) {
			# we will step through each data column, representing each window (bin)
			# across the feature's region
			# any scores within this window will be collected and the mean 
			# value reported
		
			# convert the window start and stop coordinates (as percentages) to
			# actual bp
			# this depends on whether the binsize is explicitly defined in bp or
			# is a fraction of the feature length
			my ($start, $stop);
			if ($Data->metadata($column, 'bin_size') =~ /bp$/) {
				# the bin size is explicitly defined
			
				# the start and stop points are relative to either the feature
				# start (always 0) or the end (the feature length), depending
				# upon whether the 5' or 3' end of the feature
			
				# determine this by the sign of the start position
				if ($Data->metadata($column, 'start') < 0 and $strand >= 0) {
					# the start position is less than 0, implying the 5' end
					# the reference position will be the feature start on plus strand
					$start = $fstart + $Data->metadata($column, 'start');
					$stop  = $fstart + $Data->metadata($column, 'stop');
				}
				elsif ($Data->metadata($column, 'start') < 0 and $strand < 0) {
					# the start position is less than 0, implying the 5' end
					# the reference position will be the feature end on minus strand
					$start = $fstop - $Data->metadata($column, 'start');
					$stop  = $fstop - $Data->metadata($column, 'stop');
				}
				elsif ($Data->metadata($column, 'start') >= 0 and $strand >= 0) {
					# the start position is greather than 0, implying the 3' end
					# the reference position will be the feature start on plus strand
					$start = $fstop + $Data->metadata($column, 'start');
					$stop  = $fstop + $Data->metadata($column, 'stop');
				}
				elsif ($Data->metadata($column, 'start') >= 0 and $strand < 0) {
					# the start position is greather than 0, implying the 3' end
					# the reference position will be the feature end on minus strand
					$start = $fstart - $Data->metadata($column, 'start');
					$stop  = $fstart - $Data->metadata($column, 'stop');
				}
				else {
					warn " unable to unable to identify region orientation: start " . 
						$Data->metadata($column, 'start') . ", strand $strand\n";
					return;
				}
			}
			else {
				# otherwise the bin size is based on feature length
				if ($strand >= 0) {
					# forward plus strand
					$start = int( $fstart + 
						($Data->metadata($column, 'start') * 0.01 * $length) + 0.5);
					$stop  = int( $fstart + 
						($Data->metadata($column, 'stop') * 0.01 * $length) - 1 + 0.5);
				}
				else {
					# reverse minus strand
					$start = int( $fstop - 
						($Data->metadata($column, 'start') * 0.01 * $length) + 0.5);
					$stop  = int( $fstop - 
						($Data->metadata($column, 'stop') * 0.01 * $length) + 1 + 0.5);
				}
			}
		
			# collect the data for this bin
			my $score = $row->get_score(
				'ddb'         => $ddb,
				'dataset'     => $dataset,
				'chromo'      => $feature->seq_id,
				'start'       => $start,
				'stop'        => $stop,
				'strand'      => $strand,
				'method'      => $method,
				'stranded'    => $stranded,
			);
			$row->value($column, $score);
		}
	}	
}


sub collapse_all_features {
	# collect all the transcripts and collapse them
	
	# check for collapsed transcript length index
	my $length_i = $Data->find_column('Merged_Transcript_Length');
	unless (defined $length_i) {
		$length_i = $Data->add_column('Merged_Transcript_Length');
	}
	$Data->iterate( sub {
		my $row = shift;
		my $feature = $row->seqfeature;
		my $collSeqFeat = collapse_transcripts($feature);
		$Data->store_seqfeature($row->row_index, $collSeqFeat);
		if ($length_i) {
			$row->value($length_i, get_transcript_length($feature));
		}
	} );
	
	# reorder columns
	my @neworder = (0 .. $startcolumn - 1, $length_i, $startcolumn .. $Data->last_column - 1);
	$Data->reorder_column(@neworder);
	$startcolumn++;
	return 1;
}



sub record_the_bin_values {
	
	# get the passed values
	my ($row, $length, $regionscores) = @_;
	
	
	# assign the scores to the bins in the region
	for my $column ($startcolumn..($Data->last_column) ) {
		# we will step through each data column, representing each window (bin)
		# across the feature's region
		# any scores within this window will be collected and the mean 
		# value reported
		
		# record nulls if no data returned
		unless (scalar keys %$regionscores) {
			$row->value($column, '.');
			next;
		}
		
		# convert the window start and stop coordinates (as percentages) to
		# actual bp
		# this depends on whether the binsize is explicitly defined in bp or
		# is a fraction of the feature length
		my ($start, $stop);
		if ($Data->metadata($column, 'bin_size') =~ /bp$/) {
			# the bin size is explicitly defined
			
			# the start and stop points are relative to either the feature
			# start (always 0) or the end (the feature length), depending
			# upon whether the 5' or 3' end of the feature
			
			# determine this by the sign of the start position
			if ($Data->metadata($column, 'start') < 0) {
				# the start position is less than 0, implying the 5' end
				# the reference position will be the feature start, or 0
				$start = $Data->metadata($column, 'start');
				$stop  = $Data->metadata($column, 'stop');
			}
			else {
				# the start position is greather than 0, implying the 3' end
				# the reference position will be the feature end, or length
				$start = $Data->metadata($column, 'start') + $length;
				$stop  = $Data->metadata($column, 'stop') + $length;
			}
		}
		else {
			# otherwise the bin size is based on feature length
			$start = sprintf "%.0f", ( 
				$Data->metadata($column, 'start') * 0.01 * $length) + 1;
			$stop = sprintf "%.0f", ( 
				$Data->metadata($column, 'stop') * 0.01 * $length);
		}
		
		# collect the scores for this window
		my @scores;
		for (my $n = $start; $n <= $stop; $n++) {
			# we will walk through each bp in the window looking for a score
			push @scores, $regionscores->{$n} if exists $regionscores->{$n};
		}
		
		# calculate the value
		my $window_score = calculate_score($method, \@scores);
		
		# record the value
		$row->value($column, $window_score);
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
	
	# set new name
	my $name = sprintf "%s..%s%s", $start, $stop, $unit;
	
	# set new index
	my $new_index = $Data->add_column($name);
	
	# set the metadata using passed and global variables
	# set the metadata
	$Data->metadata($new_index, 'start' , $start);
	$Data->metadata($new_index, 'stop' , $stop);
	$Data->metadata($new_index, 'dataset' , $dataset);
	$Data->metadata($new_index, 'method' , $method);
	$Data->metadata($new_index, 'bin_size' , $binsize . $unit);
	$Data->metadata($new_index, 'strand' , $stranded);
	if ($set_strand) {
		$Data->metadata($new_index, 'strand_implied', 1);
	}
	if ($data_database) {
		$Data->metadata($new_index, 'db', $data_database);
	}
}





__END__

=head1 NAME

get_binned_data.pl

A program to collect data in bins across a list of features.

=head1 SYNOPSIS
 
 get_binned_data.pl [--options] <filename>
  
  Options for existing files:
  --in <filename>                  (txt bed gff gtf refFlat ucsc)
  
  Options for new files:
  --db <name | filename>
  --feature <type | type:source | alias>, ...
  
  Options for data collection:
  --ddb <name | filename>
  --data <dataset_name | filename>
  --method [mean|median|stddev|min|max|range|sum|          (mean)
            count|pcount|ncount]
  --strand [all|sense|antisense]                            (all)
  --force_strand
  --long
  --(no)log
  
  Bin specification:
  --bins <integer>                                          (10)
  --ext <integer>                                           (0)
  --extsize <integer>
  --min <integer>
  
  Post-processing:
  --sum
  --smooth
  
  General options:
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

=item --db <name | filename>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A database is 
required for generating new data files with features. This option may 
skipped when using coordinate information from an input file (e.g. BED 
file), or when using an existing input file with the database indicated 
in the metadata. For more information about using annotation databases, 
see L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. 

=item --feature <type | type:source | alias>,...

Specify the type of feature from which to collect values. This is required 
only for new feature tables. Three types of values may be passed: the 
feature type, feature type and source expressed as 'type:source', or an 
alias to one or more feature types. Aliases are specified in the 
C<biotoolbox.cfg> file and provide a shortcut to a list of one or more 
database features. More than one feature may be included as a 
comma-delimited list (no spaces). 
  
=item --ddb <name | filename>

If the data to be collected is from a second database that is separate 
from the annotation database, provide the name of the data database here. 
Typically, a second C<Bio::DB::SeqFeature::Store> or BigWigSet database 
is provided here. 

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

=item --strand [all|sense|antisense]

Specify whether stranded data should be collected. Three values are 
allowed: all datasets should be collected (default), only sense 
datasets, or only antisense datasets should be collected.

=item --force_strand

For features that are not inherently stranded (strand value of 0)
or that you want to impose a different strand, set this option when
collecting stranded data. This will reassign the specified strand for
each feature regardless of its original orientation. This requires the
presence of a "strand" column in the input data file. This option only
works with input file lists of database features, not defined genomic
regions (e.g. BED files). Default is false.

=item --long

Indicate that the dataset from which scores are collected are 
long features (counting genomic annotation for example) and not point 
data (microarray data or sequence coverage). Normally long features are 
only recorded at their midpoint, leading to inaccurate representation at 
some windows. This option forces the program to collect data separately 
at each window, rather than once for each file feature or region and 
subsequently assigning scores to windows. Execution times may be 
longer than otherwise. Default is false.

=item --(no)log

Dataset values are (not) in log2 space and should be treated 
accordingly. Output values will be in the same space.

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

=item --sum

Indicate that the data should be averaged across all features at
each position, suitable for graphing. A separate text file will be
written with the suffix '_summed' with the averaged data. The default 
is false.

=item --smooth

Indicate that windows without values should (not) be interpolated
from neighboring values. The default is false.

=item --out <filename>

Specify the output file name. 

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --cpu <integer>

Specify the number of CPU cores to execute in parallel. This requires 
the installation of Parallel::ForkManager. With support enabled, the 
default is 2. Disable multi-threaded execution by setting to 1. 

=item --version

Print the version number.

=item --help

This help text.

=back

=head1 DESCRIPTION

This program will collect data across a gene or feature body into numerous 
percentile bins. It is used to determine if there is a spatial 
distribution preference for the dataset over gene bodies. The number 
of bins may be specified as a command argument (default 10). Additionally, 
extra bins may be extended on either side of the gene (default 0 on either 
side). The bin size is determined as a percentage of gene length.

The program writes out a tim data formatted text file. It will also 
optionally generate a summary or average profile for all of the features. 

=head1 EXAMPLES

These are some examples of some common scenarios for collecting data.

=over 4

=item Collect scores in intervals

You want to collect the mean score from a bigWig file in 10% intervals 
across each feature in a Bed file.

  get_binned_data.pl --data scores.bw --in input.bed

=item Collect scores in intervals plus extended regions

You want to collect the maximum score in 5% intervals across each each 
feature as well as five 100 bp intervals outside of each interval.

  get_binned_data.pl --bins 20 --method max --ext 5 --extsize 100 --data \
  scores.bw --in input.txt

=item Collect scores in intervals for genes

You want to collect stranded alignment counts from a Bam file for genes 
in an annotation database, then generate a profile graph.

  get_binned_data.pl --db annotation --feature gene --strand sense --value \
  count --method sum --data alignments.bam --out gene_profile --sum
  
  graph_profile.pl --in gene_profile_summed.txt --min 0 --max 100

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
