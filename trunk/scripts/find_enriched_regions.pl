#!/usr/bin/perl

# A script to look for enriched regions for a specific microarray data set

use strict;
use Getopt::Long;
use Statistics::Lite qw(mean median stddevp);
use Pod::Usage;
use Data::Dumper;
use File::Basename qw(fileparse);
use Bio::Range;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	generate_tim_data_structure
);
use tim_db_helper qw(
	open_db_connection
	process_and_verify_dataset
	get_region_dataset_hash
	get_chromo_region_score
);
use tim_file_helper qw(
	write_tim_data_file
	convert_and_write_to_gff_file
);
use tim_db_helper::config;
my $VERSION = '1.5.8';

print "\n This script will find enriched regions for a specific data set\n\n";

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

# Initialize values
my (
	$dataset,
	$database,
	$outfile,
	$win,
	$step,
	$sdlimit,
	$threshold,
	$method,
	$deplete,
	$tolerance,
	$feat,
	$genes,
	$trim,
	$sort,
	$log,
	$html,
	$gff,
	$help,
	$print_version,
	$debug,
); # command line variables

# Command line options
GetOptions( 
	'data=s'    => \$dataset, # the dataset to look for enriched regions
	'db=s'      => \$database, # database name
	'out=s'     => \$outfile, # output file name
	'win=i'     => \$win, # size of the window to scan the genome
	'step=i'    => \$step, # step size to move the window along the genome
	'sd=f'      => \$sdlimit, # the number of standard deviations above mean to set as the threshold
	'thresh=s'  => \$threshold, # the explicitly given threshold value
	'method=s'  => \$method, # method of combining values
	'deplete'   => \$deplete, # look for depleted regions instead of enriched
	'tol=i'     => \$tolerance, # tolerance for merging windows
	'feat!'     => \$feat, # collect feature information
	'genes'     => \$genes, # indicate a text file of overlapping genes shoudl be written
	'trim!'     => \$trim, # do trim the windows
	'sort!'     => \$sort, # sort the windows by score
	'log!'      => \$log, # dataset is in log2 space
	'gff'       => \$gff, # write out a gff file
#	'html'      => \$html, # write out a html file with hyperlinks to gbrowse
	'debug'     => \$debug, # limit to chromosome 1 for debugging purposes
	'help'      => \$help, # print help
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
	print " Biotoolbox script find_enriched_regions.pl, version $VERSION\n\n";
	exit;
}


# Check for required flags and assign undefined variables default values
unless ($database) {
	die " You must define a database!\n Use --help for more information\n";
}

$outfile =~ s/\.txt$//; # strip extension, it'll be added later

# window defaults
unless ($win) {
	$win = 250;
}
unless ($step) {
	# default is to use the window size
	$step = $win;
}
unless (defined $tolerance) {
	# default is 1/2 of the window size
	$tolerance = int($win / 2);
}

# threshold default
unless ($threshold) {
	unless ($sdlimit) {
		$sdlimit = 1.5;
		print " Using default limit of 1.5 standard deviation above the mean\n";
	}
}

# set the method of combining scores
if (defined $method) {
	unless ($method eq 'mean' or $method eq 'median') {
		die " unknown method '$method'!\n";
	}
} 
else {
	#  default is average
	$method = 'mean';
}

# set log2 default
unless (defined $log) {
	# default is false
	$log = 0;
}

# set trimming default
unless (defined $trim) {
	$trim = 0;
}





#### Main #####

## Preparing global variables
	# This program predates my development of the tim data file and memory
	# data structures described in 'tim_file_helper.pm'. These structures
	# were bolted on afterwards. As such, the program still uses lots of
	# arrays described immediately below, and only at the end prior to 
	# output is a tim data structure generated.
my @windows; # a temporary of the found enriched windows
	# this is an array of arrays
	# the first array is an array of the found windows, and consists of
	# the second array which is comprised of the following elements
	# $chr, $start, $end, $window_score
my @genelist; # an array of gene names overlapping the regions
my %chrom2length; # a hash to store the chromosome lengths
my $db = open_db_connection($database);


## Begin the search for the enriched windows

# First need to get the data set name if not already provided
# Next, determine the threshold
# Finally, walk through the genome looking for enriched windows. These will be
# stored in the @windows array

# Check or request the dataset
$dataset = process_and_verify_dataset( {
	'db'      => $db,
	'dataset' => $dataset,
	'single'  => 1,
} );

# get a simplified dataset name
my $dataset_name;
if ($dataset =~ /^ (?: http | ftp | file ) : \/* (.+) \. (?: bb | bw | bam ) $/xi) {
	# use the file's basename as the dataset name
	$dataset_name = $1;
}
else {
	# a database feature
	$dataset_name = $dataset;
}

# Generate output file name if necessary
unless ($outfile) {
	$outfile = "$dataset_name\_w$win\_s$step\_t$threshold";
}


## Determine the cutoff value
# the actual value used to determine if a region is enriched or not
unless (defined $threshold) { 
	# otherwise determine cutoff value from the dataset distribution
	print "  Determining threshold....\n";
	$threshold = go_determine_cutoff();
}
my $cutoff = $threshold; 
if ($log) {
	# we assume that the cutoff value is also provided as a log2 number
	# but the calculations require the value to be de-logged
	$cutoff = 2 ** $cutoff;
}

## Find the enriched regions
go_find_enriched_regions();
unless (@windows) { # exit the program if nothing found
	warn " No windows found!\n";
	exit;
}

# DEBUGGING: printing out the intermediate @windows array
if ($debug) {
	open FILE, ">$outfile.debug.post_windows.txt";
	print FILE Dumper(\@windows);
	close FILE;
}



## Merge the windows into larger regions
# this will merge the overlapping windows in @windows and put them into back
go_merge_windows(\@windows);
print "  Merged windows into " . scalar @windows . " windows\n";

# DEBUGGING: printing out the intermediate @windows array
if ($debug) {
	open FILE, ">$outfile.debug.post_merge1.txt";
	print FILE Dumper(\@windows);
	close FILE;
}



## Trim the merged windows of datapoints that are below the threshold
if ($trim) {
	print "  Trimming windows....\n";
	go_trim_windows();
	
	# DEBUGGING: printing out the intermediate @windows array
	if ($debug) {
		open FILE, ">$outfile.debug.post_trim.txt";
		print FILE Dumper(\@windows);
		close FILE;
	}
}


## Double check the merging
# Go back quickly through double-checking that we don't have two neighboring windows
# I still seem to have some slip through....
go_merge_windows(\@windows);
print "  Merged trimmed windows into " . scalar @windows . " windows\n";

# DEBUGGING: printing out the intermediate @windows array
if ($debug) {
	open FILE, ">$outfile.debug.post_merge2.txt";
	print FILE Dumper(\@windows);
	close FILE;
}


## Get score for final window
print "  Calculating final score of merged, trimmed windows....\n";
get_final_window_score();


## Sort the array by the final score of the windows
if ($sort) {
	print "  Sorting windows by score....\n";
	sort_data_by_final_score();
}

## Name the windows
name_the_windows();


## Identify features for merged windows
if ($feat) {
	print "  Identifying associated genomic features....\n";
	get_overlapping_features();
}


## Generate the final primary data hash
# this data hash is compatible with the tim data text format described in
# tim_data_helper.pm and tim_file_helper.pm
# converting to this structure makes it easier for writing files 
# via tim_file_helper
# can you tell that this was bolted on long after writing the original script?
my $main_data_ref = generate_main_data_hash();
unless ($main_data_ref) {
	die " unable to generate main data hash!\n";
}



## Print the output
# write standard output data file
my $write_success = write_tim_data_file( {
	'data'     => $main_data_ref,
	'filename' => $outfile,
} );
if ($write_success) {
	print " Wrote data file '$write_success'\n";
}
else {
	print " unable to write data file!\n";
}

# write html output file
if ($html) { 
	write_html_file();
}

# write gff file
if ($gff) { 
	my $method;
	if ($deplete) {
		$method = 'depleted_region';
	}
	else {
		$method = 'enriched_region';
	}
	my $gff_file = convert_and_write_to_gff_file( {
		'data'     => $main_data_ref,
		'score'    => 5,
		'name'     => 0,
		'source'   => 'find_enriched_regions.pl',
		'method'   => $method,
		'version'  => 3,
		'filename' => $outfile,
	} );
	if ($gff_file) {
		print " Wrote GFF file '$gff_file'\n";
	}
	else {
		print " unable to write GFF file!\n";
	}
}

print "All done!\n\n";




############# Subroutines ###################



### Determine the cutoff values for the dataset
sub go_determine_cutoff {
	
	# collect sample of values from the dataset
	my @chromosomes = $db->seq_ids; 
	unless (@chromosomes) {
		die " unable to identify chromosome sequences in the database!\n";
	}
	
	# select chromosome randomly
	my $n = rand (scalar @chromosomes);
	while ($chromosomes[$n] =~ /chrm|chrmt|chrmt|NA/i) {
		# avoid that mitochrondrial chromosome like the plague!
		$n = rand (scalar @chromosomes);
	}
	
	# collect statistics on the chromosome
	print " Sampling '$dataset_name' values across chromosome '" . 
		$chromosomes[$n] . "'...\n";
	my $chromosome = $db->segment($chromosomes[$n]);
	
	my $mean = get_chromo_region_score( {
		'db'           => $db,
		'dataset'      => $dataset,
		'method'       => 'mean',
		'chromo'       => $chromosomes[$n],
		'start'        => 1,
		'stop'         => $chromosome->length,
		'log'          => $log,
	} );
	unless ($mean) { 
		die " unable to determine mean value for '$dataset'!\n";
	}
	my $stdev = get_chromo_region_score( {
		'db'           => $db,
		'dataset'      => $dataset,
		'method'       => 'stddev',
		'chromo'       => $chromosomes[$n],
		'start'        => 1,
		'stop'         => $chromosome->length,
		'log'          => $log,
	} );
	unless ($stdev) { 
		die " unable to determine stdev value for '$dataset'!\n";
	}
	print "   the mean value is $mean and standard deviation $stdev\n";
	
	# calculate the actual cuttoff value, depending on enriched or depleted
	my $value; 
	if ($deplete) { 
		# look for depleted regions
		# cutoff is the defined multiples of std dev above the mean
		$value = $mean - ($stdev * $sdlimit); 
	} 
	else { 
		# default to look for enriched regions
		# cutoff is the defined multiples of std dev above the mean
		$value = $mean + ($stdev * $sdlimit); 
	}
	
	# conclusion
	print "  Using a threshold of $value ($sdlimit std devs)\n";
	return $value;
}



### Walk through each chromosome sequentially looking for windows of enrichment
sub go_find_enriched_regions {
	
	# print messages
	if ($deplete) {
		print "  Looking for depleted regions ";
	} 
	else {
		print "  Looking for enriched regions ";
	}
	print "using window $method values\n";
	
	
	## collect chromosomes and data
	# get list of chromosomes
	my @chromosomes = $db->seq_ids; 
	unless (@chromosomes) {
		die " unable to retrieve chromosomes from the database!\n";
	}
	# Get the names of chromosomes to avoid
	my @excluded_chromosomes = 
			$TIM_CONFIG->param("$database\.chromosome_exclude") ||
			$TIM_CONFIG->param('default_db.chromosome_exclude');
	
	
	# walk through each chromosome
	foreach my $chr (@chromosomes) {
		
		# check for excluded chromosomes
		my $skip_chr = 0;
		foreach (@excluded_chromosomes) {
			if ($chr eq $_) {
				$skip_chr = 1;
				last;
			}
		}
		next if $skip_chr;
		
		# START DEBUGGING # 
		if ($debug) {
			# LIMIT TO ONE CHROMOSOME
			if ($chr eq 'chr2') {last} 
		}
		# END DEBUGGING #
		
		# collect the dataset values for the current chromosome
		# store in a hash the position (key) and values
		print "  Searching $chr....\n";
		
		# generate a segment representing the chromosome
		# due to fuzzy name matching, we may get more than one back
		my @segments = $db->segment($chr);
		# need to find the right one
		my $chrobj;
		while (@segments) {
			$chrobj = shift @segments;
			last if $chrobj->seq_id eq $chr;
		}
		
		# walk windows along the chromosome and find enriched windows
		my $length = $chrobj->length; # length of the chromosome
		$chrom2length{$chr} = $length; # remember this length for later
		for (my $start = 1; $start < $length; $start += $step) {
			# define the window to look in
			my $end = $start + $win -1;
			# ensure don't go over chromosome length
			if ($end > $length) {
				$end = $length;
			} 
						
			# determine window value
			my $window_score = get_chromo_region_score( {
				'db'         => $db,
				'dataset'    => $dataset,
				'method'     => $method,
				'chromo'     => $chr,
				'start'      => $start,
				'stop'       => $end,
				'log'        => $log,
			} );
			unless ($window_score) {
				#print "no values at $chr:$start..$end!\n"; 
				next;
			}
			if ($log) {
				$window_score = 2 ** $window_score;
			}
			
			# calculate if window passes threshold
			if ($deplete) { 
				# depleted regions
				if ($window_score <= $cutoff) { 
					# score passes our threshold
					push @windows, [$chr, $start, $end, $window_score];
				}
			} 
			else { 
				# enriched regions
				if ($window_score >= $cutoff) { 
					# score passes our threshold
					push @windows, [$chr, $start, $end, $window_score];
				}
			}
		}
		
	}
	print "  Found " . scalar @windows . " windows for $dataset.\n";
}


### Condense the list of overlapping windows
sub go_merge_windows {
 	my $array_ref = shift;
 	
 	# set up new target array and move first item over
 	my @merged; 
 	push @merged, shift @{$array_ref};
 	
 	while ( @$array_ref ) {
 		my $window = shift @{$array_ref};
 		
 		# first check whether chromosomes are equal
 		if ( $merged[-1][0] eq $window->[0] ) {
 			# same chromosome
 			
 			# generate Range objects
 			my $range1 = Bio::Range->new(
 				-start  => $merged[-1][1],
 				-end    => $merged[-1][2] + $tolerance
 				# we add tolerance only on side where merging might occur
 			);
 			my $range2 = Bio::Range->new(
 				-start  => $window->[1],
 				-end    => $window->[2]
 				# we add tolerance only on side where merging might occur
 			);
			
			# check for overlap
			if ( $range1->overlaps($range2) ) {
				# we have overlap
				# merge second range into first
				my ($mstart, $mstop, $mstrand) = $range1->union($range2);
				
				# update the merged window
				$merged[-1][1] = $mstart;
				$merged[-1][2] = $mstop;
				
				# score is no longer relevent
				$merged[-1][3] = '.';
			}
			else {
				# no overlap
				push @merged, $window;
			}
 		}
 		
 		
 		else {
 			# not on same chromosome
 			# move onto old array
 			push @merged, shift @{$array_ref};
 		}
 	}
 	
 	# Put the merged windows back
 	@{$array_ref} = @merged;
 	
}

	
	
	
### Fine trim the merged windows
sub go_trim_windows {
	# since the merged window represents a combined score from a region of 
	# multiple datapoints, the region may include data points on the edges of 
	# the region that actually do not cross the threshold
	# this subroutine will find the true endpoints of the enriched region
	# note that this method won't work well if the dataset is noisy
	
	
	# Walk through the list of merged windows
	foreach my $window (@windows) {
		
		# calculate extended window size
		my $start = $window->[1] - $tolerance;
		my $stop = $window->[2] + $tolerance;
		
		# check sizes so we don't go over limit
		if ($start < 1) {
			$start = 1;
		}
		if ($stop > $chrom2length{ $window->[0] }) {
			$stop = $chrom2length{ $window->[0] };
		}
		
		# get values across the extended window
		my %pos2score = get_region_dataset_hash( {
			'db'       => $db,
			'dataset'  => $dataset,
			'chromo'   => $window->[0],
			'start'    => $start,
			'stop'     => $stop,
			'absolute' => 1,
		} );
		unless (%pos2score) {
			# we should be able to! this region has to have scores!
			die " unable to generate value hash for window $window->[0]:$start..$stop!\n";
		}
		if ($debug) {
			print " window $start..$stop scores at ", 
				join(",", sort {$a <=> $b} keys %pos2score), "\n";
		}
		
		# de-log if necessary
		if ($log) {
			foreach (keys %pos2score) {
				$pos2score{$_} = 2 ** $pos2score{$_};
			}
		}
		
		# look for first position whose value crosses the threshold
		foreach my $pos (sort {$a <=> $b} keys %pos2score) {
			if ($deplete) {
				# looking for depleted regions
				if ( $pos2score{$pos} <= $cutoff) {
					# we found one!
					# assign it to the window start
					$window->[1] = $pos;
					
					# go no further
					last;
				}
			}
			else {
				# lookig for enriched regions
				if ($pos2score{$pos} >= $cutoff) {
					# we found one!
					# assign it to the window start
					$window->[1] = $pos;
					
					# go no further
					last;
				}
			}
		}
		
		# look for last position whose value crosses the threshold
		foreach my $pos (sort {$b <=> $a} keys %pos2score) {
			# reverse sort order
			if ($deplete) {
				# looking for depleted regions
				if ($pos2score{$pos} <= $threshold) {
					# we found one!
					# assign it to the window start
					$window->[2] = $pos;
					
					# go no further
					last;
				}
			}
			else {
				# lookig for enriched regions
				if ($pos2score{$pos} >= $threshold) {
					# we found one!
					# assign it to the window start
					$window->[2] = $pos;
					
					# go no further
					last;
				}
			}
		}
		
		# what to do with single points?
		if ($window->[1] == $window->[2]) { 
			# a single datapoint
			# make it at least a small window half the size of $win
			# this is just to make the data look pretty and avoid a region
			# of 1 bp, which could interfere with the final score later on
			$window->[1] -= int($step/4);
			$window->[2] += int($step/4);
			
			# check sizes so we don't go over limit
			if ($window->[1] < 1) {
				$window->[1] = 1;
			}
			if ($window->[2] > $chrom2length{ $window->[0] }) {
				$window->[2] = $chrom2length{ $window->[0] };
			}
		}
		
	}

}
	

	
### Get the final score for the merged window and other things
sub get_final_window_score {
	for my $i (0..$#windows) {
		# arrays currently have $chr, $start, $end, $region_score
		# we will calculate a new region score, as well as the window size
		
		# determine window size
		my $size = $windows[$i][2] - $windows[$i][1] + 1;
		
		# replace the current score with the window size
		$windows[$i][3] = $size;
		# arrays now have $chr, $start, $end, $size
		
		# re-calculate window score
		$windows[$i][4] = get_chromo_region_score( {
				'db'       => $db,
				'dataset'  => $dataset, 
				'chromo'   => $windows[$i][0],
				'start'    => $windows[$i][1],
				'stop'     => $windows[$i][2],
				'method'   => $method,
				'log'      => $log,
		} );
		
		# arrays now have $chr, $start, $end, $size, $finalscore
	}
}


### Collect the overlapping features of the enriched windows
sub get_overlapping_features {
	
	# walk through the list of windows
	for my $i (0..$#windows) {
		
		# collect the genomic features in the region
		my @features = $db->get_features_by_location(
			-seq_id    => $windows[$i][1],
			-start     => $windows[$i][2],
			-end       => $windows[$i][3],
		);
		
		my (@orf_list, @rna_list, @non_gene_list);
		foreach my $f (@features) {
			
			# collect info
			my $type = $f->primary_tag;
			my $name = $f->display_name;
			
			# determine which category to put the feature in
			# this is pretty generic, how useful is this, really?
			# keeping it this way to avoid breakage and rewrites....
			if ($type =~ /rna/i) {
				push @rna_list, "$type $name";
			}
			elsif ($type =~ /gene|orf/i) {
				push @orf_list, "$type $name";
			}
			else {
				push @non_gene_list, "$type $name";
			}
		}
		
		
		# tack on the list of features
		if (@orf_list) {
			push @{ $windows[$i] }, join(", ", @orf_list);
		}
		else {
			push @{ $windows[$i] }, '.';
		}
		if (@rna_list) {
			push @{ $windows[$i] }, join(", ", @rna_list);
		}
		else {
			push @{ $windows[$i] }, '.';
		}
		if (@non_gene_list) {
			push @{ $windows[$i] }, join(", ", @non_gene_list);
		}
		else {
			push @{ $windows[$i] }, '.';
		}
		# arrays now have $chr, $start, $end, $size, $finalscore,
		# plus, if requested, $orf_list, $rna_list, $non_gene_list
		
		# Record the gene names if requested
		if ($genes) { 
			foreach (@orf_list) {
				# each item is 'type name'
				# only keep the name
				push @genelist, (split / /)[1]; 
			}
			foreach (@rna_list) {
				# each item is 'type name'
				# only keep the name
				push @genelist, (split / /)[1];
			}
		}
	}
}


### Sort the array be final score
sub sort_data_by_final_score {
	
	# sort by the score value and place into a temporary array
	my @temp;
	if ($deplete) {
		# increasing order
		@temp = sort { $a->[4] <=> $b->[4] } @windows;
	}
	else {
		# decreasing order
		@temp = sort { $b->[4] <=> $a->[4] } @windows;
	}
	
	# put back
	@windows = @temp;
}


### Name the windows
sub name_the_windows {
	my $count = 1;
	foreach (@windows) {
		my $name = $dataset_name . '_win' . $count;
		unshift @{ $_ }, $name;
		$count++;
	}
	# arrays now have $name, $chr, $start, $end, $size, $region_score
	# plus, if requested, $orf_list, $rna_list, $non_gene_list
	
}



### Generate the main data hash compatible with tim_file_helper.pm
# this is bolted on after writing the main script for compatibility with 
# tim_file_helper modules
sub generate_main_data_hash {
	
	# determine the data feature
	my $feature;
	if ($deplete) {
		# depleted regions
		$feature = 'depleted_regions';
	}
	else {
		# enriched regions
		$feature = 'enriched_regions';
	}
	
	# generate a new data structure
	my $data = generate_tim_data_structure(
		$feature,
		'WindowID',
		'Chromosome',
		'Start',
		'Stop',
		'Size',
		'Score',
	) or die " unable to generate tim data structure!\n";
	
	# Add metadata
	$data->{'db'} = $database;
	
	# window metadata
		# traditionally with the genome feature datasets, extra pertinant
		# information regarding the window generation goes under the start 
		# metadata
	$data->{2}{'win'} = $win;
	$data->{2}{'step'} = $step;
	if ($trim) {
		$data->{2}{'trimmed'} = 1;
	} else {
		$data->{2}{'trimmed'} = 0;
	}
	
	# score metadata
	$data->{5}{'log2'} = $log;
	$data->{5}{'method'} = $method;
	$data->{5}{'dataset'} = $dataset;
	if ($threshold) {
		$data->{5}{'threshold'} = $threshold;
	} else {
		$data->{5}{'standard_deviation_limit'} = $sdlimit;
		$data->{5}{'threshold'} = $cutoff;
	}
	
	# add feature metadata if it was requested
	if ($feat) {
		# metadata keys
		$data->{6} = {
			# the orf list
			'name'  => 'ORF_Features',
			'index' => 6,
		};
		$data->{7} = {
			# the orf list
			'name'  => 'RNA_Features',
			'index' => 7,
		};
		$data->{8} = {
			# the orf list
			'name'  => 'Non-gene_Features',
			'index' => 8,
		};
		
		# update the number of columns
		$data->{'number_columns'} = 9;
		
	}
	
	# add the column headers
	# the windows array never had column headers, so we'll add them here
	unshift @windows, []; # put in space for the header
	# add the names to the row 0 array
	for (my $i = 0; $i < $data->{'number_columns'}; $i++) {
		$windows[0][$i] = $data->{$i}{'name'};
	}
	
	# associate the @windows array with the data table in the structure
	$data->{'data_table'} = \@windows;
	
	# update last row
	$data->{'last_row'} = scalar(@windows) - 1;
	
	# return the reference to the generated data hash
	return $data;
}




### Write the html output file
# this is old code for writing an html table file of the hits
# it worked with GBrowse 1.x, but has not been vetted with GBrowse 2.x and 
# is unlikely to work with the new browser
# I'm keeping the code here in case I ever want to revisit it in the future
# in the meantime, I've deleted references to it in the POD
sub write_html_file {
	open HTML, ">$outfile.html";
	# print the head
	print HTML "<HTML>\n\n<HEAD>\n\t<TITLE>\n\t$outfile\n\t</TITLE>\n</HEAD>\n";
	print HTML "<BODY BGCOLOR=\"#FFFFFF\" text=\"#000000\">\n\n\n";
	print HTML "<H2>\n$outfile\n</H2>\n";
	# print the parameters
	print HTML "Program $0<p>\n"; # header information marked with #
	print HTML "Database $database<p>\n";
	print HTML "Dataset $dataset<p>\n";
	print HTML "Window $win<p>\n";
	print HTML "Step $step<p>\n";
	if ($threshold) {
		print HTML "Threshold $threshold<p>\n";
	} else {
		print HTML "Standard deviation limit $sdlimit<p>\n";
		print HTML "Threshold $cutoff<p>\n";
	}
	if ($deplete) {
		print HTML "Searching for depleted regions<p>\n";
	} else {
		print HTML "Searching for enriched regions<p>\n";
	}
	print HTML "Method $method<p>\n";
	if ($trim) {
		print HTML "Windows trimmed<p>\n";
	} else {
		print HTML "Windows not trimmed<p>\n";
	}
	if ($feat) {
		print HTML "Features collected<p>\n";
	}
	# print the table headers
	print HTML "\n\n<table border cellspacing=0 cellpadding=3>\n";
	print HTML "<tr bgcolor=\"#ccccff\">\n";
	print HTML "<th align=\"left\">WindowID</th>\n";
	print HTML "<th align=\"left\">Position</th>\n";
	print HTML "<th align=\"left\">Size</th>\n";
	print HTML "<th align=\"left\">$method score</th>\n";
	print HTML "<th align=\"left\">ORF features</th>\n";
	print HTML "<th align=\"left\">RNA features</th>\n";
	print HTML "<th align=\"left\">Non-gene features</th>\n";
	print HTML "</tr>\n\n";
	
	# print the data
	my $gbrowse = "http://localhost/cgi-bin/gbrowse/$database/?"; # hyperlink to machine & gbrowse
	for my $i (0..$#windows) {
		# http://localhost/cgi-bin/gbrowse/cerevisiae/?name=chr1:137066..145763;enable=RSC_ChIP_ypd_244k;h_region=chr1:142066..143368@lightcyan
		my $position = "$windows[$i][1]:$windows[$i][2]..$windows[$i][3]"; # chromosome:start-stop
		
		# set size of selected region for gbrowse
		# data should now have $window_name, $chr, $start, $end, $size, $region_score
		my ($start, $stop);
		if ($windows[$i][4] < 500) { 
			# size is under 500 bp, set 2 kb region
			$start = $windows[$i][2] - 1000;
			$stop = $windows[$i][2] + 999;
		} 
		elsif ($windows[$i][4]  < 1000) { 
			# size is under 1 kb, set 5 kb region
			$start = $windows[$i][2] - 2500;
			$stop = $windows[$i][2] + 2499;
		} 
		elsif ($windows[$i][4]  < 5000) { 
			# size is under 5 kb, set 10 kb region
			$start = $windows[$i][2] - 5000;
			$stop = $windows[$i][2] + 4999;
		} 
		else { 
			# size is really big, set 20 kb region
			$start = $windows[$i][2] - 10000;
			$stop = $windows[$i][2] + 9999;
		}
		if ($start < 0) {$start = 1} # in case we're at left end of chromosome
		
		# generate hypertext reference
		my $region = "name=$windows[$i][1]:$start\..$stop";
		my $track = "enable=$dataset";
		my $highlight = "h_region=$position\@bisque";
		my $link = $gbrowse . "$region;$track;$highlight";
		
		# generate table data
		print HTML "<tr>\n";
		print HTML "<td><a href=\"$link\">$windows[$i][0]</a></td>\n"; 
			# windowID with hyperlink text
		print HTML "<td>$position</td>\n"; # position
		print HTML "<td>$windows[$i][4]</td>\n"; # size
		print HTML "<td>$windows[$i][5]</td>\n"; # score
		print HTML "<td>$windows[$i][6]</td>\n"; # ORF features
		print HTML "<td>$windows[$i][7]</td>\n"; # RNA features
		print HTML "<td>$windows[$i][8]</td>\n"; # non-gene features
		print HTML "</tr>\n\n";
	}
	# close up
	print HTML "</table>\n";
	print HTML "</BODY>\n</HTML>\n";
	close HTML;
	print " Wrote html file $outfile.html\n";
}





__END__

=head1 NAME

find_enriched_regions.pl

=head1 SYNOPSIS
 
 find_enriched_regions.pl --db <db_name> [--options]
 
  Options:
  --db <name|file.gff3>
  --data <dataset | filename>
  --out <filename>
  --win <integer>
  --step <integer>
  --tol <integer>
  --thresh <number>
  --sd <number>
  --method [mean|median]
  --deplete
  --(no)trim
  --(no)sort
  --(no)feat
  --(no)log
  --genes
  --gff
  --gz
  --version
  --help

 
=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <database_name>

=item --db <name|file.gff3>

Specify the name of the BioPerl SeqFeature::Store database to use as
source. Alternatively, a single GFF3 file may be loaded into a in-memory
database. 

=item --data <dataset | filename>

Specify the name of the dataset from which to collect the scores. 
If not specified, then the data set may be chosen interactively 
from a presented list.
Alternatively, the name of a data file may be provided. Supported 
file types include BigWig (.bw), BigBed (.bb), or single-end Bam 
(.bam). The file may be local or remote.

=item --out <filename>

Specify the output file name. If not specified, then it will be 
automatically generated from dataset, window, step, and threshold 
values.

=item --win <integer>

Specify the genomic bin window size in bp. Default value is 250 bp.

=item --step <integer>

Specify the step size for moving the window. Default value is 
equal to the window size.

=item --tol <integer>

Specify the tolerance distance when merging windows together. 
Windows not actually overlapping but within this tolerance 
distance will actually be merged. Default value is 1/2 the 
window size.

=item --thresh <number>

Specify the window score threshold explicitly rather than calculating
a threshold based on standard deviations from the mean.

=item --sd <number>

Specify the multiple of standard deviations above the mean as a
window score threshold. Default is 1.5 standard deviations. Be 
quite careful with this method as it attempts to pull all of the 
datapoints out of the database to calculate the mean and 
standard deviation - which may be acceptable for limited tiling 
microarrays but not acceptable for next generation sequencing 
data with single bp resolution.

=item --method [mean|median]

Specify the method for combining score values within each window 
when determining whether the window exceeds the threshold.
Default method is mean.

=item --deplete

Specify whether depleted regions should be reported instead. 
For example, windows whose scores are 1.5 standard deviations 
below the mean, rather than above.

=item --(no)trim

Indicate that the merged windows should (not) be trimmed of below 
threshold scores on the ends of the window. Normally when windows 
are merged, there may be some data points on the ends of the 
windows whose scores don't actually pass the threshold, but were 
included because the entire window mean (or median) exceeded 
the threshold. This step removes those data points. The default 
behavior is false (notrim).

=item --(no)feat

Indicate that features overlapping the windows should be 
identified. The default behavior is false.

=item --genes

Write out a text file containing a list of the found overlapping genes.

=item --gff

Indicate that a GFF version 3 file should be written out in
addition to the data text file.

=item --(no)log

Flag to indicate that source data is (not) log2, and to calculate 
accordingly and report as log2 data.

=item --gz

Compress the output file through gzip.

=item --version

Print the version number.

=item --help

Display the POD documentation.

=back

=head1 DESCRIPTION

This program will search for regions in the genome that are enriched for a 
particular data set. It walks through each chromosome using a 
window of specified size (default 500 bp) and specified step size (default 100
bp). Data scores within a window that exceed a determined threshold
will be noted. Adjacent windows are merged and then trimmed on the ends to the
minimum thresholded window.

The threshold scores for identifying an enriched region may either be 
explicitly set or automatically determined from the mean and standard 
deviation (SD) of the entire collection of datapoints across the genome. 
This, of course, assumes a normal distribution of datapoint scores, which may 
or may not be suitable for the particular dataset. Note that the 
automatic method may not appropriate for very extremely large datasets 
(e.g. next generation sequencing) as it attempts to calculate the mean and SD 
on all of the datapoints in the database. 

The program writes out a tim data formatted text file consisting of chromosome, 
start, stop, score, and overlapping gene or non-gene genomic features. It 
will optionally write a GFF file.


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





