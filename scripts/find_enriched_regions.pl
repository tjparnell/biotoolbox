#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Statistics::Lite qw(mean median stddevp);
use Pod::Usage;
use File::Basename qw(fileparse);
use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure
	format_with_commas
);
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	get_chromosome_list
	get_region_dataset_hash
	get_chromo_region_score
	get_chromosome_list
);
use Bio::ToolBox::file_helper qw(
	write_tim_data_file
	convert_and_write_to_gff_file
);
use Bio::ToolBox::db_helper::config;
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};
# use Data::Dumper;
my $VERSION = '1.14';

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
	$main_database,
	$data_database,
	$outfile,
	$win,
	$step,
	$strand,
	$sdlimit,
	$threshold,
	$method,
	$value,
	$deplete,
	$tolerance,
	$feat,
	$genes,
	$trim,
	$min_size,
	$sort,
	$log,
	$html,
	$gff,
	$gz,
	$help,
	$cpu,
	$print_version,
	$debug,
); # command line variables

# Command line options
GetOptions( 
	'data=s'    => \$dataset, # the dataset to look for enriched regions
	'db=s'      => \$main_database, # main or feature database name
	'ddb=s'     => \$data_database, # data database name
	'out=s'     => \$outfile, # output file name
	'win=i'     => \$win, # size of the window to scan the genome
	'step=i'    => \$step, # step size to move the window along the genome
	'strand=s'  => \$strand, # specific strand to search
	'sd=f'      => \$sdlimit, # the number of standard deviations above mean to set as the threshold
	'thresh=s'  => \$threshold, # the explicitly given threshold value
	'method=s'  => \$method, # method of combining values
	'value=s'   => \$value, # type of value to collect
	'deplete'   => \$deplete, # look for depleted regions instead of enriched
	'tol=i'     => \$tolerance, # tolerance for merging windows
	'feat!'     => \$feat, # collect feature information
	'genes'     => \$genes, # indicate a text file of overlapping genes shoudl be written
	'trim!'     => \$trim, # do trim the windows
	'min=i'     => \$min_size, # do not trim below this size
	'sort!'     => \$sort, # sort the windows by score
	'log!'      => \$log, # dataset is in log2 space
	'gff'       => \$gff, # write out a gff file
	'gz!'       => \$gz, # write compressed file
	'debug'     => \$debug, # limit to chromosome 1 for debugging purposes
	'cpu=i'     => \$cpu, # number of execution threads
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

$outfile =~ s/\.txt$//; # strip extension, it'll be added later

# window defaults
unless ($win) {
	# collect the window size from biotoolbox.cfg
	$win = $BTB_CONFIG->param("$main_database\.window") || 
			$BTB_CONFIG->param("$data_database\.window") || 
			$BTB_CONFIG->param("default_db.window") || 500;
}
unless ($step) {
	# default is to use the window size
	$step = $win;
}
unless (defined $tolerance) {
	# default is 1/2 of the window size
	$tolerance = int($win / 2);
}
unless ($min_size) {
	# don't go below the window size
	$min_size = $win;
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
	unless ($method =~ /^mean|median|sum$/) {
		die " unknown method '$method'!\n";
	}
} 
else {
	#  default is average
	$method = 'mean';
}

# set the type of score value to collect
if (defined $value) {
	unless ($value =~ /^score|count|length$/) {
		die " unknown method '$value'!\n";
	}
} 
else {
	#  default is score
	$value = 'score';
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

# set strand
if ($strand) {
	if ($strand =~ /^[f\+]/) {
		# forward strand
		$strand = 1;
	}
	elsif ($strand =~ /^[r\-]/i) {
		# reverse strand
		$strand = -1;
	}
	else {
		die " unrecognized strand character!\n";
	}
}
else {
	# default is no strand
	$strand = 0;
}
# set strandedness of data collection
my $strandedness;
if ($strand == 0) {
	# collect data from both strands since strand is not requested
	$strandedness = 'all';
}
else {
	# stranded data, we assume the user wants sense orientation
	# this could become a user defined option in the future....
	$strandedness = 'sense';
}

# check parallel support
if ($parallel) {
	# conservatively enable 2 cores
	$cpu ||= 2;
}
else {
	# disable cores
	print " disabling parallel CPU execution, no support present\n" if $cpu;
	$cpu = 0;
}



#### Main #####

## Preparing global variables
	# This program predates my development of the tim data file and memory
	# data structures described in Bio::ToolBox::file_helper. These structures
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


## Open databases
my ($fdb, $ddb); # feature and data databases
if ($main_database and $data_database) {
	# two separate databases defined
	$fdb = open_db_connection($main_database) or 
		die " unable to establish connection to database '$main_database'!\n";
	$ddb = open_db_connection($data_database) or 
		die " unable to establish connection to database '$data_database'!\n";
} 
elsif ($main_database and !$data_database) {
	# only main database defined
	$fdb = open_db_connection($main_database) or 
		die " unable to establish connection to database '$main_database'!\n";
	$ddb = $fdb; # reuse
} 
elsif (!$main_database and $data_database) {
	# only data database defined
	$ddb = open_db_connection($data_database) or 
		die " unable to establish connection to database '$data_database'!\n";
	if ($feat) {
		warn " no main or feature database defined! disabling search for features\n\n";
		$feat = 0;
	}
} 
elsif (!$main_database and !$data_database and $dataset =~ /\.(?:bw|bb|bam)$/) {
	# dataset is a bigwig, bigbed, or bam file
	# use this as the data database
	$ddb = open_db_connection($dataset) or 
		die " unable to establish connection to database '$dataset'!\n";
	if ($feat) {
		warn " no main or feature database defined! disabling search for features\n\n";
		$feat = 0;
	}
}
else {
	die " no databases defined! see help for more information\n";
}


## Begin the search for the enriched windows

# First need to get the data set name if not already provided
# Next, determine the threshold
# Finally, walk through the genome looking for enriched windows. These will be
# stored in the @windows array

# Check or request the dataset
$dataset = verify_or_request_feature_types(
	'db'      => $ddb,
	'feature' => $dataset,
	'single'  => 1,
	'prompt'  => " Enter the number of the feature or dataset to scan for" . 
					" enrichment   ",
);

# get a simplified dataset name
my $dataset_name;
if ($dataset =~ /^ (?: http | ftp | file ) \: \/* (.+) \. (?: bb | bw | bam ) $/xi) {
	# use the file's basename as the dataset name
	($dataset_name, undef, undef) = fileparse($1);
}
else {
	# a database feature
	$dataset_name = $dataset;
}
my $start_time = time;



## Determine the cutoff value
# the actual value used to determine if a region is enriched or not
unless (defined $threshold) { 
	# otherwise determine cutoff value from the dataset distribution
	print " Determining threshold....\n";
	$threshold = go_determine_cutoff();
}



## Generate output file name if necessary
# we're doing this here instead of before or later because we need 
# the threshold value and the debug option needs an outfile name
unless ($outfile) {
	$outfile = "$dataset_name\_w$win\_s$step\_t$threshold";
}





## Find the enriched regions
# print messages
if ($deplete) {
	print " Looking for depleted regions ";
} 
else {
	print " Looking for enriched regions ";
}
print "using window $method $value values\n";
	
# execute based on the number of CPUs available
if ($cpu > 1) {
	# parallel execution
	print " Forking into $cpu children for parallel execution\n";
	parallel_execution();
}

else {
	# single threaded execution
	single_execution();
}



## Sort the array by the final score of the windows
if ($sort) {
	print " Sorting windows by score....\n";
	sort_data_by_final_score();
}
else {
	@windows = sort { $a->[0] cmp $b->[0] or $a->[1] <=> $b->[1] } @windows;
}

## Name the windows
name_the_windows();


## Identify features for merged windows
if ($feat) {
	print " Identifying associated genomic features....\n";
	get_overlapping_features();
}


## Generate the final primary data hash
# this data hash is compatible with the tim data text format described in
# Bio::ToolBox::data_helper and Bio::ToolBox::file_helper
# converting to this structure makes it easier for writing files 
# via Bio::ToolBox::file_helper
# can you tell that this was bolted on long after writing the original script?
my $main_data_ref = generate_main_data_hash();
unless ($main_data_ref) {
	die " unable to generate main data hash!\n";
}



## Print the output
# write standard output data file
my $write_success = write_tim_data_file(
	'data'     => $main_data_ref,
	'gz'       => $gz,
	'filename' => $outfile,
);
if ($write_success) {
	print " Wrote data file '$write_success'\n";
}
else {
	print " unable to write data file!\n";
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
	my $gff_file = convert_and_write_to_gff_file(
		'data'     => $main_data_ref,
		'score'    => 6,
		'name'     => 0,
		'source'   => 'find_enriched_regions.pl',
		'method'   => $method,
		'version'  => 3,
		'filename' => $outfile,
		'gz'       => $gz,
	);
	if ($gff_file) {
		print " Wrote GFF file '$gff_file'\n";
	}
	else {
		print " unable to write GFF file!\n";
	}
}

printf " Finished in %.1f minutes\n", (time - $start_time)/60;




############# Subroutines ###################



### Determine the cutoff values for the dataset
sub go_determine_cutoff {
	
	# select the largest chromosome
	my $length = 1;
	my $chromosome;
	foreach ( get_chromosome_list($ddb, 1) ) {
		
		# this is an array of chromosome name and length
		my ($name, $size) = @$_;
		
		# check its size
		if ($size > $length) {
			$chromosome = $name;
			$length = $size;
		}
	}
	
	# collect statistics on the chromosome
	print " Sampling '$dataset_name' values across largest chromosome " . 
		"$chromosome...\n";
	my $mean = get_chromo_region_score(
		'db'           => $ddb,
		'dataset'      => $dataset,
		'method'       => 'mean',
		'value'        => $value,
		'chromo'       => $chromosome,
		'start'        => 1,
		'stop'         => $length,
		'log'          => $log,
		'strand'       => $strand,
		'stranded'     => $strandedness,
	);
	unless ($mean) { 
		die " unable to determine mean value for '$dataset'!\n";
	}
	my $stdev = get_chromo_region_score(
		'db'           => $ddb,
		'dataset'      => $dataset,
		'method'       => 'stddev',
		'value'        => $value,
		'chromo'       => $chromosome,
		'start'        => 1,
		'stop'         => $length,
		'log'          => $log,
		'strand'       => $strand,
		'stranded'     => $strandedness,
	);
	unless ($stdev) { 
		die " unable to determine stdev value for '$dataset'!\n";
	}
	print "   the mean value is $mean and standard deviation $stdev\n";
	
	# calculate the actual cuttoff value, depending on enriched or depleted
	my $cutoff; 
	if ($deplete) { 
		# look for depleted regions
		# cutoff is the defined multiples of std dev above the mean
		$cutoff = $mean - ($stdev * $sdlimit); 
	} 
	else { 
		# default to look for enriched regions
		# cutoff is the defined multiples of std dev above the mean
		$cutoff = $mean + ($stdev * $sdlimit); 
	}
	
	# conclusion
	print " Using a threshold of $cutoff ($sdlimit std devs)\n";
	return $cutoff;
}


### Single threaded execution
sub single_execution {
	
	## collect chromosomes and data
	# walk through each chromosome
	foreach ( get_chromosome_list($ddb, 1) ) {
		
		# this is an array of chromosome name and length
		my ($chr, $length) = @$_;
		go_find_enriched_regions($chr, $length);
	}
	print " Found " . format_with_commas(scalar @windows) . " total windows for $dataset.\n";
	
	# check for windows
	unless (@windows) { # exit the program if nothing found
		exit;
	}
}


### Parallel forked execution
sub parallel_execution {
	
	# Initiate the fork manager
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_finish( sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $child_windows) = @_;
		push @windows, @{$child_windows};
	});
	
	# the chromosome list
	my @chromosomes = sort { $b->[1] <=> $a->[1] } get_chromosome_list($ddb, 1);
		# get_chromosome_list() returns arrays of chromosome [name, length]
		# sort decreasing by the length
	# generate an execution list based on size 
		# the idea is to equalize the amount of work to do for each child process
	my @list; # array of arrays, [process][name,length]
	my $i = 1;
	while (@chromosomes) {
		push @{ $list[$i] }, shift @chromosomes;
		$i++;
		$i = 1 if $i > $cpu;
	}
	
	# Count the chromosomes in parallel processess
	foreach my $n (1 .. $cpu) {
		$pm->start and next;

		### In child
		foreach ( @{$list[$n]} ) {
			# process each chromosome in this child list
			go_find_enriched_regions( @{$_} );
		}
		$pm->finish(0, \@windows); 
	}
	$pm->wait_all_children;
	
	# check for windows
	print " Found " . format_with_commas(scalar @windows) . " total windows for $dataset.\n";
	unless (@windows) { # exit the program if nothing found
		exit;
	}
}



### Walk through each chromosome sequentially looking for windows of enrichment
sub go_find_enriched_regions {
	
	# passed chromosome name and length
	my ($chr, $length) = @_;
	
	# collect the dataset values for the current chromosome
	my @chr_windows; # an array of windows for the current chromosome
	
	# walk windows along the chromosome and find enriched windows
	$chrom2length{$chr} = $length; # remember this length for later
	for (my $start = 1; $start < $length; $start += $step) {
		# define the window to look in
		my $end = $start + $win -1;
		# ensure don't go over chromosome length
		if ($end > $length) {
			$end = $length;
		} 
					
		# determine window value
		my $window_score = get_chromo_region_score(
			'db'         => $ddb,
			'dataset'    => $dataset,
			'method'     => $method,
			'value'      => $value,
			'chromo'     => $chr,
			'start'      => $start,
			'stop'       => $end,
			'log'        => $log,
			'strand'     => $strand,
			'stranded'   => $strandedness,
		);
		next unless (defined $window_score);
		
		# calculate if window passes threshold
		if ($deplete) { 
			# depleted regions
			if ($window_score <= $threshold) { 
				# score passes our threshold
				push @chr_windows, [$chr, $start, $end, $window_score];
			}
		} 
		else { 
			# enriched regions
			if ($window_score >= $threshold) { 
				# score passes our threshold
				push @chr_windows, [$chr, $start, $end, $window_score];
			}
		}
	}
	
	# Check for windows
	return unless @chr_windows;
	
	# Merge the windows into larger regions
	go_merge_windows(\@chr_windows);
	
	# Trim the merged windows of datapoints that are below the threshold
	if ($trim) {
		go_trim_windows(\@chr_windows);
	
		# Double check the merging
		# I still seem to have some slip through....
		go_merge_windows(\@chr_windows);
	}
	
	# Get score for final window
	get_final_window_score(\@chr_windows);
	
	# finished
	push @windows, @chr_windows;
	print " Found ", format_with_commas(scalar @windows), " merged", 
		 $trim ? ", trimmed " : " ", "windows passing threshold on chromosome $chr\n"; 
}


### Condense the list of overlapping windows
sub go_merge_windows {
 	
 	my $window_ref = shift;
 	
 	# set up new target array and move first item over
 	my @merged; 
 	push @merged, shift @{$window_ref};
 	
 	while (@{$window_ref}) {
 		my $window = shift @{$window_ref};
 		
		# check for overlap
		if ( $merged[-1]->[2] + $tolerance >= $window->[1] ) {
			# we have overlap with current end and new window start
			# merge them
			$merged[-1]->[2] = $window->[2];
			
			# score no longer relavent
			$merged[-1]->[3] = '.';
			
			if ( $merged[-1]->[1] > $window->[1] ) {
				# in the very rare (unlikely) situation where next window has lower 
				# start than the current window, adjust accordingly
				# merge second range into first
				$merged[-1]->[1] = $window->[1];
			}
		}
		elsif ( $merged[-1]->[1] > $window->[1] ) {
			# in the very rare (unlikely) situation where next window has lower 
			# start than the current window, adjust accordingly
			# merge second range into first
			$merged[-1]->[1] = $window->[1];
		}
		else {
			# no overlap
			push @merged, $window;
		}
 	}
 	
 	# Put the merged windows back
 	push @{$window_ref}, @merged;
}

	
	
	
### Fine trim the merged windows
sub go_trim_windows {
	# since the merged window represents a combined score from a region of 
	# multiple datapoints, the region may include data points on the edges of 
	# the region that actually do not cross the threshold
	# this subroutine will find the true endpoints of the enriched region
	# note that this method won't work well if the dataset is noisy
	
	
 	my $window_ref = shift;
 	
	# Walk through the list of merged windows
	foreach my $window ( @{$window_ref} ) {
		
		# calculate extended window size
		my $start = $window->[1] - $tolerance;
		my $stop  = $window->[2] + $tolerance;
		
		# check sizes so we don't go over limit
		if ($start < 1) {
			$start = 1;
		}
		if ($stop > $chrom2length{ $window->[0] }) {
			$stop = $chrom2length{ $window->[0] };
		}
		
		# get values across the extended window
		my %pos2score = get_region_dataset_hash(
			'db'       => $ddb,
			'dataset'  => $dataset,
			'chromo'   => $window->[0],
			'value'    => $value,
			'start'    => $start,
			'stop'     => $stop,
			'absolute' => 1,
			'strand'   => $strand,
			'stranded' => $strandedness,
		);
		unless (%pos2score) {
			# we should be able to! this region has to have scores!
			warn " unable to generate value hash for window $window->[0]:$start..$stop!\n";
			next;
		}
			
		
		# look for first position whose value crosses the threshold
		foreach my $pos (sort {$a <=> $b} keys %pos2score) {
			if ($deplete) {
				# looking for depleted regions
				if ( $pos2score{$pos} <= $threshold) {
					# we found one!
					# assign it to the window start
					$window->[1] = $pos;
					
					# go no further
					last;
				}
			}
			else {
				# lookig for enriched regions
				if ($pos2score{$pos} >= $threshold) {
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
		
		# check the minimum size of the window
		my $length = $window->[2] - $window->[1] + 1;
		if ($length < $min_size) {
			# find how much to add to either side to make it a minimum size
			my $add = int( ( ($min_size - $length) / 2) + 0.5);
			
			# adjust the window
			$window->[1] -= $add;
			$window->[2] += $add;
			
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
	my $window_ref = shift; 
	
	for my $i (0 .. scalar(@$window_ref)-1) {
		# arrays currently have $chr, $start, $end, $region_score
		# we will calculate a new region score, as well as the window size
		
		# replace the current score with the window size
		$window_ref->[$i][3] = $window_ref->[$i][2] - $window_ref->[$i][1] + 1;
		# arrays now have $chr, $start, $end, $size
		
		# add the strand column
		$window_ref->[$i][4] = $strand;
		
		# re-calculate window score
		$window_ref->[$i][5] = get_chromo_region_score(
				'db'       => $ddb,
				'dataset'  => $dataset, 
				'chromo'   => $window_ref->[$i][0],
				'start'    => $window_ref->[$i][1],
				'stop'     => $window_ref->[$i][2],
				'method'   => $method,
				'value'    => $value,
				'log'      => $log,
				'strand'   => $strand,
				'stranded' => $strandedness,
		);
		
		# arrays now have $chr, $start, $end, $size, $strand, $finalscore
	}
}


### Collect the overlapping features of the enriched windows
sub get_overlapping_features {
	
	# first check the database
	my $fdb_ref = ref $fdb;
	unless ($fdb_ref =~ m/^Bio::DB::SeqFeature::Store/i) {
		warn " feature database is not Bio::DB::SeqFeature::Store!\n" .
			" unable to search for features\n";
		return;
	}
	
	# walk through the list of windows
	for my $i (0..$#windows) {
		
		# collect the genomic features in the region
		my @features = $fdb->get_features_by_location(
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
		# arrays now have $chr, $start, $end, $size, $strand, $finalscore,
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
		@temp = sort { $a->[5] <=> $b->[5] } @windows;
	}
	else {
		# decreasing order
		@temp = sort { $b->[5] <=> $a->[5] } @windows;
	}
	
	# put back
	@windows = @temp;
}


### Name the windows
sub name_the_windows {
	my $count = 1;
	foreach (@windows) {
		unshift @{ $_ }, $dataset_name . '_win' . $count;
		$count++;
	}
	# arrays now have $name, $chr, $start, $end, $size, $strand, $region_score
	# plus, if requested, $orf_list, $rna_list, $non_gene_list
}



### Generate the main data hash compatible with Bio::ToolBox::file_helper
# this is bolted on after writing the main script for compatibility with 
# Bio::ToolBox::file_helper modules
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
		'Strand',
		'Score',
	) or die " unable to generate tim data structure!\n";
	
	# Add metadata
	$data->{'db'} = $main_database || $data_database;
	
	# window metadata
		# traditionally with the genome feature datasets, extra pertinant
		# information regarding the window generation goes under the start 
		# metadata
	$data->{2}{'win'} = $win;
	$data->{2}{'step'} = $step;
	if ($trim) {
		$data->{4}{'trimmed'} = 1;
	} else {
		$data->{4}{'trimmed'} = 0;
	}
	$data->{4}{'min_size'} = $min_size;
	
	# score metadata
	$data->{6}{'log2'} = $log;
	$data->{6}{'method'} = $method;
	$data->{6}{'value'} = $value;
	$data->{6}{'dataset'} = $dataset;
	$data->{6}{'threshold'} = $threshold;
	if ($sdlimit) {
		$data->{6}{'standard_deviation_limit'} = $sdlimit;
	}
	if ($main_database and $data_database) {
		$data->{6}{'data_database'} = $data_database;
	}
	
	# add feature metadata if it was requested
	if ($feat) {
		# metadata keys
		$data->{6} = {
			# the orf list
			'name'  => 'ORF_Features',
			'index' => 6,
		};
		$data->{8} = {
			# the orf list
			'name'  => 'RNA_Features',
			'index' => 7,
		};
		$data->{9} = {
			# the orf list
			'name'  => 'Non-gene_Features',
			'index' => 8,
		};
		
		# update the number of columns
		$data->{'number_columns'} = 10;
		
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




__END__

=head1 NAME

find_enriched_regions.pl

A script to find enriched regions in a dataset using a simple threshold.

=head1 SYNOPSIS
 
 find_enriched_regions.pl --db <db_name> [--options]
 
 find_enriched_regions.pl --data <file> [--options]
 
  Options:
  --db <name | filename>
  --ddb <name | filename>
  --data <dataset | filename>
  --out <filename>
  --win <integer>
  --step <integer>
  --tol <integer>
  --thresh <number>
  --sd <number>
  --method [mean|median|sum]
  --value [score|count|length]
  --strand [f|r]
  --deplete
  --trim
  --min <integer>
  --sort
  --feat
  --log
  --genes
  --gff
  --gz
  --cpu <integer>
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <name>

Specify the name or file of a Bio::DB::SeqFeature::Store database 
or BigWigSet database from which to collect chromosomes, data scores, 
and/or overlapping feature annotations. Features may only be 
collected from a SeqFeature::Store databases. This may be skipped 
if a big file (bigWig, bigBed, or Bam) file is provided as the dataset.

=item --ddb <name>

When data scores are present in a separate database from annotation,
then specify the second data-specific database. The same options 
apply as --db.

=item --data <dataset | filename>

Specify the feature type or primary_tag of the dataset within the 
database from which to collect the scores. If not specified, then 
the data set may be chosen interactively from a presented list.

Alternatively, the name of a single data file may be provided. 
Supported file types include BigWig (.bw), BigBed (.bb), or 
single-end Bam (.bam). The file may be local or remote.

=item --out <filename>

Specify the output file name. If not specified, then it will be 
automatically generated from dataset, window, step, and threshold 
values.

=item --win <integer>

Specify the genomic bin window size in bp. If not specified, 
then the default window size is retrieved from the biotoolbox.cfg 
configuration file. Default is 500 bp.

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

=item --method [mean|median|sum]

Specify the method for combining score values within each window 
when determining whether the window exceeds the threshold.
Default method is mean.

=item --value [score|count|length]

Specify the type of value to collect from the dataset. The 
default value type is score.

--strand [f|r]

Optionally specify a specific strand from which to restrict the 
data collection. This requires that the dataset supports 
stranded data collection (GFF3, Bam, bigBed, BigWigSet). 
Default is to collect from both strands.

=item --deplete

Specify whether depleted regions should be reported instead. 
For example, windows whose scores are 1.5 standard deviations 
below the mean, rather than above.

=item --trim

Indicate that the merged windows should be trimmed of below 
threshold scores on the ends of the window. Normally when windows 
are merged, there may be some data points on the ends of the 
windows whose scores don't actually pass the threshold, but were 
included because the entire window mean (or median) exceeded 
the threshold. This step removes those data points. The default 
behavior is false.

=item --min <integer>

Set the minimum window size in bp when trimming the merged windows. 
The default value is equal to the search window size.

=item --sort

Indicate that the regions should be sorted by their score. 
Sort order is decreasing for enriched regions and increasing 
for depleted regions. Default is false (they should be 
ordered mostly by coordinate).

=item --feat

Indicate that features overlapping the windows should be 
identified. The default behavior is false.

=item --genes

Write out a text file containing a list of the found overlapping genes.

=item --gff

Indicate that a GFF version 3 file should be written out in
addition to the data text file.

=item --log

Flag to indicate that source data is log2, and to calculate 
accordingly and report as log2 data.

=item --gz

Compress the output file through gzip.

=item --cpu <integer>

Specify the number of CPU cores to execute in parallel. This requires 
the installation of Parallel::ForkManager. With support enabled, the 
default is 2. Disable multi-threaded execution by setting to 1. 

=item --version

Print the version number.

=item --help

Display the POD documentation.

=back

=head1 DESCRIPTION

This program will search for regions in the genome that are enriched for a 
particular data set. It walks through each chromosome using a 
window of specified size (default 500 bp) and specified step size (default 
same as window). Data scores within a window that exceed a determined threshold
will be noted. Adjoining windows (within a specific tolerance, default is 
1/2 of window size) are merged. The windows may be optionally trimmed 
of flanking below-threshold positions.

The method for identifying enrichment is based on a very simple criteria: 
a window is kept if the mean (or median) value for the window is greater 
than (or less than for depletion) the threshold. No statistics or False 
Discovery Rate is calculated.

The threshold may be automatically determined based on a calculated 
mean and standard deviation from a sampling of the dataset. For sampling 
purposes, the largest chromosome, scaffold, or sequence defined in the 
database is used. A multiple (default 1.5X) of the standard deviation 
is used to set the threshold.

If an annotation database is provided, gene, ORF, non-coding RNA, or 
other features may optionally be identified overlapping the enriched 
regions.

The program writes out a tab-delimited text file consisting of chromosome, 
start, stop, strand, score, and overlapping gene or non-gene genomic 
features. It will optionally write a GFF3 file.

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
