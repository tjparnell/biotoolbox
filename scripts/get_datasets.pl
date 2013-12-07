#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use File::Spec;
use Statistics::Lite qw(
	sum
	mean
	median
	min
	max
	range
	stddevp
);
use Bio::ToolBox::data_helper qw(
	find_column_index
	format_with_commas
	splice_data_structure
);
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	check_dataset_for_rpm_support
	get_new_feature_list 
	get_new_genome_list 
	get_feature
	get_chromo_region_score
	get_region_dataset_hash
);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	write_tim_data_file
);
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};
my $VERSION = '1.14';


print "\n A program to collect data for a list of features\n\n";


### Display quick help
unless (@ARGV) {
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}

### Get command line options and initialize values

# Initialize values
my (  
	$infile,
	$new,
	$outfile,
	$main_database,
	$data_database,
	$feature,
	$method,
	$value_type,
	$log,
	$stranded,
	$subfeature,
	$extend,
	$start_adj,
	$stop_adj,
	$fstart,
	$fstop,
	$limit,
	$position,
	$win,
	$step,
	$set_strand,
	$gz,
	$cpu,
	$help,
	$print_version,
); 
my @datasets; # an array of names of dataset values to be retrieved

# Command line options
GetOptions( 
	'in=s'       => \$infile, # load a pre-existing file
	'new'        => \$new, # generate a new file
	'out=s'      => \$outfile, # name of new output file 
	'db=s'       => \$main_database, # main or annotation database name
	'ddb=s'      => \$data_database, # data database
	'feature=s'  => \$feature, # name of genomic feature to analyze
	'method=s'   => \$method, # method of collecting & reporting data
	'value=s'    => \$value_type, # type of data to collect
	'data=s'     => \@datasets, # the list of datasets to collect data from
	'log!'       => \$log, # dataset is in log2 space
	'strand=s'   => \$stranded, # indicate strandedness of data
	'exons!'     => \$subfeature, # indicate to restrict to subfeatures
	'extend=i'   => \$extend, # extend the size of the genomic feature
	'start=i'    => \$start_adj, # adjustment to relative position
	'stop=i'     => \$stop_adj, # adjustment relative position
	'fstart=f'   => \$fstart, # fractional start position
	'fstop=f'    => \$fstop, # fractional stop position
	'limit=i'    => \$limit, # size limit to fractionate a feature
	'pos=s'      => \$position, # set the relative feature position
	'win=i'      => \$win, # indicate the size of genomic intervals
	'step=i'     => \$step, # step size for genomic intervals
	'force_strand|set_strand' => \$set_strand, # enforce a specific strand
				# force_strand is preferred option, but respect the old option
	'gz!'        => \$gz, # compress output file
	'cpu=i'      => \$cpu, # number of execution threads
	'help'       => \$help, # request help
	'version'    => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# print help if requested
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script get_datasets.pl, version $VERSION\n\n";
	exit;
}


# Check for required values
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		# we will assume the user wants to make a new file
		$new = 1;
	}
}
if ($new) {
	unless ($outfile) {
		die " You must define an output filename!";
	}
	unless ($feature) {
		die  " You must define an input file or new feature type! see help\n";
	}
}
if (defined $start_adj or defined $stop_adj) {
	unless (defined $start_adj and defined $stop_adj) {
		die " You must define both start and stop coordinate adjustments!\n";
	}
}
if (defined $fstart or defined $fstop) {
	unless (defined $fstart and defined $fstop) {
		die " You must define both fstart and fstop coordinate adjustments!\n";
	}
}


# Assign default values
set_defaults();
my $start_time = time;

# Assign database for new feature lists
if ($new and not defined $main_database) {
	# creating a new feature list requires a main database 
	# otherwise we will postpone this till after loading the input file
	
	if (defined $data_database) {
		# reuse the data database
		$main_database = $data_database;
	}
	elsif (@datasets and $feature eq 'genome') {
		# we could use a dataset file only if we're collecting genome windows
		# take the first element
		$main_database = $datasets[0];
	}
	else {
		die " You must define a database or an appropriate dataset file! see help\n";
	}
}





### Initialize main data, database, and datasets

# Generate or open main feature list and identify columns
my (
	$main_data_ref, 
	$name_i, 
	$type_i, 
	$id_i,
	$chromo_i, 
	$start_i, 
	$stop_i, 
	$strand_i
) = get_main_data_ref();


# Open main database connection
unless ($main_database) {
	# this could've been obtained either from the input file, 
	# command line, or a source data file
	# lacking that, we'll attempt to use the first dataset provided
	# and hope for the best 
	if ($data_database) {
		# use data database if defined instead of main database
		$main_database = $data_database;
	}
	elsif (@datasets) {
		# we hope this some sort of indexed data file like bigWig or Bam
		$main_database = $datasets[0];
		print " no database defined, using $main_database\n";
	}
	else {
		die " no database defined! see help\n";
	}
}
my $mdb = open_db_connection($main_database) or 
	die " unable to establish database connection to '$main_database'!\n";

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

# Check the datasets
unless ($datasets[0] eq 'none') {
	@datasets = verify_or_request_feature_types(
		'db'      => $ddb,
		'feature' => [ @datasets ],
		'prompt'  => " Enter the dataset(s) or feature type(s) from which \n" . 
					" to collect data. Comma delimited or range is acceptable\n",
	);
}


# Working with RPM and RPKM value datasets
# global values
my %dataset2sum; # for tot
# total reads in Bam file when using rpkm method
if ($method eq 'rpm' or $method eq 'rpkm') {
	foreach my $d (@datasets) {
		print " Checking RPM support for dataset '$d'...\n";
		my $sum = check_dataset_for_rpm_support($d, $ddb, $cpu);
		if ($sum) {
			$dataset2sum{$d} = $sum;
			printf "   %s total features\n", format_with_commas($sum);
		}
		else {
			warn " $method method requested but not supported for " .
				"dataset '$d'\n using summed count instead\n";
			$method = 'sum'; 
				# this could negatively impact any subsequent datasets
			last;
		}
	}
}



### Collect the data from each datasets

# check that we have a dataset
if ($datasets[0] eq 'none') {
	print " Nothing to collect!\n";
	if ($new) {
		my $success = write_tim_data_file(
			'data'     => $main_data_ref,
			'filename' => $outfile,
			'gz'       => $gz,
		);
		if ($success) {
			printf " wrote file $success\n";
		}
		else {
			# failure! the subroutine will have printed error messages
			print " unable to write file!\n";
		}
	}
	exit;
}

# check whether it is worth doing parallel execution
if ($cpu > 1) {
	while ($cpu > 1 and $main_data_ref->{'last_row'}/$cpu < 100) {
		# I figure we need at least 100 lines in each fork split to make 
		# it worthwhile to do the split, otherwise, reduce the number of 
		# splits to something more worthwhile
		$cpu--;
	}
}

# execute data collection in 1 or more processes
if ($cpu > 1) {
	# parallel execution
	print " Collecting $method $value_type from datasets @datasets...\n";
	print " Forking into $cpu children for parallel data collection\n";
	parallel_execution();
}

else {
	# single threaded execution
	single_execution();
}


### Finished
printf " Finished in %.1f minutes\n", (time - $start_time)/60;



############# Subroutines ######################################################



### Set default parameters if undefined
sub set_defaults {
	# assign default values
	# these are all global values that could've been assigned on the 
	# command line
	
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
	
	# check datasets
	if ($datasets[0] =~ /,/) {
		# seems to be a comma delimited list, possibly more than one?????
		my @list;
		foreach my $d (@datasets) {
			push @list, (split /,/, $d);
		}
		@datasets = @list;
	}
	
	# check method
	if ($method) {
		# check the method that was defined on the command line
		unless ($method =~ 
			m/^(?:median|mean|stddev|min|max|range|sum|count|enumerate|rpm|rpkm)$/
		) {
			die " unknown method '$method'!";
		}
		
		# set appropriate options for specific methods
		if ($method eq 'count') {
			# convenience method
			$method = 'sum';
			$value_type = 'count';
		}
		elsif ($method eq 'enumerate') {
			# convenience method
			$method = 'sum';
			$value_type = 'count';
		}
		elsif ($method eq 'rpkm') {
			$value_type = 'count';
		}
		elsif ($method eq 'rpm') {
			$value_type = 'count';
		}
	}
	else {
		# set the default to use the mean
		unless ($datasets[0] =~ m/none/) {
			$method = 'mean';
		}
	}
	
	# check the type of value to collect
	if (defined $value_type) {
		# validate the requested value type
		unless ($value_type =~ m/^score|count|length/) {
			die " unknown value type '$value_type'!\n";
		}
	}
	else {
		# default value
		unless ($datasets[0] =~ /none/) {
			$value_type = 'score';
		}
	}
	
	# check strandedness of data to collect
	if (defined $stranded) { 
		# check the strand request that was defined on the command line
		unless ($stranded =~ m/^all|antisense|sense$/i) {
			die " unknown strand '$stranded'!";
		}
	} 
	else {
		# default value
		$stranded = 'all'; 
	}
	
	# check the limit when using fractional start and stop
	if ($fstart and $fstop) {
		# fractional start and stop requested
		unless ($limit) {
			# set a minimum size limit on sub fractionating a feature
			# 1000 bp seems like a reasonable cut off, no?
			$limit = 1000;
		}
	}
	
	# check the relative position
	if (defined $position) {
		# check the position value
		unless ($position =~ m/^5|4|3|m$/i) {
			die " Unknown relative position '$position'!\n";
		}
		$position =~ s/m/4/i # change to match internal usage
	}
	else {
		# default position to use the 5' end
		$position = 5;
	}
	
	# check the output file
	unless ($outfile) {
		# overwrite the input file
		$outfile = $infile;
	}
}



### Collect the feature list and populate the 
sub get_main_data_ref {
	
	my ($data_ref, $name, $type, $id, $chromo, $start, $stop, $strand);
	
	# Generate new input file from the database
	if ($new) { 
		print " Generating a new feature list from database '$main_database'...\n";
		
		# generate feature list based on type of feature
		if ($feature eq 'genome') { 
			# generate a list of genomic intervals
			
			# if window and step size were not set at execution from the 
			# command line, then we'll be using the defaults hard encoded
			# in the module subroutine, which should be window of 500 bp 
			# and step size = window size
			
			# generate the list
			$data_ref = get_new_genome_list(
				'db'       => $main_database, 
				'win'      => $win,
				'step'     => $step,
			);
			
			# assign the column indices
			$chromo = 0;
			$start  = 1;
			$stop   = 2;
		} 
		
		else { 
			# everything else works off a feature list
			
			# generate the gene list
			$data_ref = get_new_feature_list(
				'db'        => $main_database,
				'features'  => $feature,
			);
			
			# assign the column indices
			$id     = 0;
			$name   = 1;
			$type   = 2;
		}
		
		# check
		unless ($data_ref) {
			die "no new feature list generated!";
		}
		
		# assign the current program
		$data_ref->{'program'} = $0;
	} 
	
	
	# Open a file containing the data
	else { 
		print " Loading feature list from file '$infile'...\n";
		$data_ref = load_tim_data_file($infile) or
			die "no file data loaded!";
		
		# identify relevant feature data columns 
		$name   = find_column_index($data_ref, '^name');
		$type   = find_column_index($data_ref, '^type|class');
		$id     = find_column_index($data_ref, '^primary_id');
		$chromo = find_column_index($data_ref, '^chr|seq|ref|ref.?seq');
		$start  = find_column_index($data_ref, '^start|position');
		$stop   = find_column_index($data_ref, '^stop|end');
		$strand = find_column_index($data_ref, '^strand');
		
		# now set missing arguments with values from loaded file
		# database
		unless ($main_database) {
			$main_database = $data_ref->{'db'};
		}
		# feature
		unless ($feature) {
			if ( exists $data_ref->{'feature'} ) {
				# load the feature defined in the file's metadata
				$feature = $data_ref->{'feature'};
			}
		}
		
		# validate the feature
		if ($feature =~ /genome|region|segment|interval/i) {
			# genomic regions?
			
			# must have coordinates
			unless (defined $chromo and defined $start) {
				# no basic coordinates!!!!
				
				# does it at least have name and type?
				if (defined $name and defined $type) {
					# whew! reassign feature to named_feature
					$feature = "Named feature";
				}
				else {
					# nothing recognized!
					die " File suggests genomic regions but coordinate" . 
						" column headers can not be found!\n";
				}
			}
			
			# make sure we have stop index defined
			unless (defined $stop) {
				# we'll use the start index for the stop index in case 
				# it isn't defined
				# this will define a region of 1 bp
				$stop = $start;
			}
		}
		
		elsif ($feature =~ /\w+/) {
			# some named feature
			
			# must have name and type columns
			unless (defined $id or (defined $name and defined $type)) {
				# no name or type defined!?
				
				# do we at least have coordinates?
				if (defined $chromo and defined $start and defined $stop) {
					# whew! reassign to genome feature
					$feature = 'region';
				}
				else {
					# nothing recognized
					die " File suggests '$feature' features but name and" . 
						" type column headers cannot be found!\n";
				}
			}
		}
		
		else {
			# feature was not explicitly defined
			
			# determine empirically from column names
			if (defined $id or (defined $name and defined $type)) {
				# named features!
				$feature = "Named feature";
			}
			elsif (defined $chromo and defined $start) {
				# genomic coordinates
				$feature = 'region';
				
				# make sure we have stop index defined
				unless (defined $stop) {
					# we'll use the start index for the stop index in case 
					# it isn't defined
					# this will define a region of 1 bp
					$stop = $start;
				}
			}
			else {
				# no useable columns defined
				die " File does not have recognizable name, type, and/or" .
					" coordinate column headers!\n Unable to proceed\n";
			}
		}
		
		# print statement on feature
		printf "  Loaded %s '%s' features.\n", 
			format_with_commas( $data_ref->{'last_row'} ),
			$feature;
	}
	
	# done
	return ($data_ref, $name, $type, $id, $chromo, $start, $stop, $strand);
}



sub parallel_execution {
	my $pm = Parallel::ForkManager->new($cpu);
	
	# generate base name for child processes
	my $child_base_name = $outfile . ".$$"; 

	# Split the input data into parts and execute in parallel in separate forks
	for my $i (1 .. $cpu) {
		$pm->start and next;
	
		#### In child ####
	
		# splice the data structure
		splice_data_structure($main_data_ref, $i, $cpu);
		
		# re-open database objects to make them clone safe
		$mdb = open_db_connection($main_database);
		if ($data_database) {
			$ddb = open_db_connection($data_database);
		}
		else {
			$ddb = $mdb;
		}
		
		# collect the dataset
		foreach my $dataset (@datasets) {
			unless ($dataset eq 'none') {
				collect_dataset($dataset);
			}
		}
		
		# write out result
		my $success = write_tim_data_file(
			'data'     => $main_data_ref,
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
	# done
}



sub single_execution {
	foreach my $dataset (@datasets) {
	
		# collect the dataset
		unless ($dataset eq 'none') {
			print " Collecting $method $value_type from dataset '$dataset'...\n";
			collect_dataset($dataset);
		}
		last if $dataset eq 'none';
	}
	
	# write the output file
	# we will rewrite the file after each collection
	# appropriate extensions and compression should be taken care of
	my $success = write_tim_data_file(
		'data'     => $main_data_ref,
		'filename' => $outfile,
		'gz'       => $gz,
	);
	if ($success) {
		printf " wrote file $success\n";
	}
	else {
		# failure! the subroutine will have printed error messages
		print " unable to write file!\n";
	}
}



# Dataset collection
sub collect_dataset {
	
	my $dataset = shift;
	
	# set the new metadata for this new dataset
	my $index = record_metadata($dataset);
	
	# check that we have strand data if necessary
	if ($set_strand) {
		unless (defined $strand_i) {
			die " requested to set strand but a strand column was not found!\n";
		}
	}
	
	# we need to determine how we will collect the data from the features
	# using genomic windows or regions, or named features?
	# are we modifying or extending the coordinates of the region or feature?
	
	# Genomic regions
	if ($feature =~ /genome|region|segment|interval/i) {
		
		# collect the genome dataset based on whether we need to modify 
		# the genomic coordinates or not
		
		if (defined $extend) {
			# extend the region on both sides
			get_extended_genome_dataset($dataset, $index);
		}
		
		elsif (defined $start_adj and defined $stop_adj) {
			# specifically defined relative start and stop positions
			get_adjusted_genome_dataset($dataset, $index);
		}
		
		elsif (defined $fstart and defined $fstop) {
			# use a subfraction of the region
			get_fractionated_genome_dataset($dataset, $index);
		}
		
		else {
			# no modifications to the position
			get_genome_dataset($dataset, $index);
		}
	}
	
	# Named features
	else {
		
		# check that we have the appropriate database opened
		my $mdb_ref = ref $mdb;
		unless ($mdb_ref =~ /^Bio::DB::SeqFeature::Store/) {
			die "\n\n Collecting data for named features is currently only " .
				"supported with\n Bio::DB::SeqFeature::Store databases! " .
				"Please either open an appropriate\n" . 
				" database or include coordinate information in your file\n";
		}
		
		
		# collect the named feature dataset based on whether we need to modify 
		# the genomic coordinates or not
		
		if ($subfeature) {
			# collect feature subfeatures
			
			get_subfeature_dataset($dataset, $index);
		}
		
		elsif (defined $extend) {
			# extend the region on both sides
			get_extended_feature_dataset($dataset, $index);
		}
		
		elsif (defined $start_adj and defined $stop_adj) {
			# specifically defined relative start and stop positions
			get_adjusted_feature_dataset($dataset, $index);
		}
		
		elsif (defined $fstart and defined $fstop) {
			# use a subfraction of the region
			get_fractionated_feature_dataset($dataset, $index);
		}
		
		else {
			# no modifications to the position
			get_feature_dataset($dataset, $index);
		}
	}
	
}



sub get_genome_dataset {
	
	# get passed arguments
	my ($dataset, $index) = @_;
	
	
	# Loop through the genomic regions
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		
		# get region score
		my $score = get_chromo_region_score(
				'db'        => $ddb,
				'dataset'   => $dataset,
				'chromo'    => $main_data_ref->{'data_table'}->[$row][$chromo_i],
				'start'     => $main_data_ref->{'data_table'}->[$row][$start_i],
				'stop'      => $main_data_ref->{'data_table'}->[$row][$stop_i],
				'strand'    => defined $strand_i ?
						$main_data_ref->{'data_table'}->[$row][$strand_i] : 0,
				'value'     => $value_type,
				'method'    => $method,
				'log'       => $main_data_ref->{$index}{'log2'},
				'stranded'  => $stranded,
		);
		unless (defined $score) {
			# this should return a value, otherwise we record an internal null
			$score = '.';
		}
		$main_data_ref->{'data_table'}->[$row][$index] = $score;
	}
}


sub get_extended_genome_dataset {
	
	# get passed arguments
	my ($dataset, $index) = @_;
	
	
	# Loop through the genomic regions
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		
		# get extended coordinates
		my $start = $main_data_ref->{'data_table'}->[$row][$start_i] - $extend;
		my $stop  = $main_data_ref->{'data_table'}->[$row][$stop_i]  + $extend;
		
		# get region score
		my $score = get_chromo_region_score(
				'db'        => $ddb,
				'dataset'   => $dataset,
				'chromo'    => $main_data_ref->{'data_table'}->[$row][$chromo_i],
				'start'     => $start,
				'stop'      => $stop,
				'strand'    => defined $strand_i ?
						$main_data_ref->{'data_table'}->[$row][$strand_i] : 0,
				'value'     => $value_type,
				'method'    => $method,
				'log'       => $main_data_ref->{$index}{'log2'},
				'stranded'  => $stranded,
		);
		unless (defined $score) {
			# this should return a value, otherwise we record an internal null
			$score = '.';
		}
		$main_data_ref->{'data_table'}->[$row][$index] = $score;
	}
}


sub get_adjusted_genome_dataset {
	
	# get passed arguments
	my ($dataset, $index) = @_;
	
	
	# Loop through the genomic regions
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		
		# get strand of the region if defined
		my $strand = defined $strand_i ? 
			$main_data_ref->{'data_table'}->[$row][$strand_i] : 0;
		
		# calculate new relative start and stop positions
		# this depends on both feature orientation and the 
		# relative position requested
		my ($start, $stop);
		
		# adjust relative to the 5' position 
		if ($position == 5) {
		
			if ($strand >= 0) {
				# feature is on the forward, watson strand
				$start = $main_data_ref->{'data_table'}->[$row][$start_i] +
					$start_adj;
				$stop  = $main_data_ref->{'data_table'}->[$row][$start_i] + 
					$stop_adj;
			}
			elsif ($strand < 0) {
				# feature is on the reverse, crick strand
				$start = $main_data_ref->{'data_table'}->[$row][$stop_i] -
					$stop_adj;
				$stop  = $main_data_ref->{'data_table'}->[$row][$stop_i] - 
					$start_adj;
			}
		}
		
		# adjust relative to the 3' position
		elsif ($position == 3) {
			
			if ($strand >= 0) {
				# feature is on the forward, watson strand
				$start = $main_data_ref->{'data_table'}->[$row][$stop_i] +
					$start_adj;
				$stop  = $main_data_ref->{'data_table'}->[$row][$stop_i] + 
					$stop_adj;
			}
			elsif ($strand < 0) {
				# feature is on the reverse, crick strand
				$start = $main_data_ref->{'data_table'}->[$row][$start_i] -
					$stop_adj;
				$stop  = $main_data_ref->{'data_table'}->[$row][$start_i] - 
					$start_adj;
			}
		}
		
		# adjust relative to the midpoint position
		elsif ($position == 4) {
			# determine the midpoint position
			my $middle = int(
				# calculate same way regardless of strand
				(
					$main_data_ref->{'data_table'}->[$row][$stop_i] +
					$main_data_ref->{'data_table'}->[$row][$start_i] 
				) / 2
			);
			if ($strand >= 0) {
				# feature is on the forward, watson strand
				$start = $middle + $start_adj;
				$stop  = $middle + $stop_adj;
			}
			elsif ($strand < 0) {
				# feature is on the reverse, crick strand
				$start = $middle - $stop_adj;
				$stop  = $middle - $start_adj;
			}
		}
		
		# get region score
		my $score = get_chromo_region_score(
				'db'        => $ddb,
				'dataset'   => $dataset,
				'chromo'    => $main_data_ref->{'data_table'}->[$row][$chromo_i],
				'start'     => $start,
				'stop'      => $stop,
				'strand'    => $strand,
				'value'     => $value_type,
				'method'    => $method,
				'log'       => $main_data_ref->{$index}{'log2'},
				'stranded'  => $stranded,
		);
		unless (defined $score) {
			# this should return a value, otherwise we record an internal null
			$score = '.';
		}
		$main_data_ref->{'data_table'}->[$row][$index] = $score;
	}
}



sub get_fractionated_genome_dataset {
	
	# get passed arguments
	my ($dataset, $index) = @_;
	
	
	# Loop through the genomic regions
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		
		# get strand of the region if defined
		my $strand = defined $strand_i ? 
			$main_data_ref->{'data_table'}->[$row][$strand_i] : 0;
		
		# get region length
		my $length = 
			$main_data_ref->{'data_table'}->[$row][$stop_i] -
			$main_data_ref->{'data_table'}->[$row][$start_i] + 1;
			
		
		# calculate new fractional start and stop positions
		# the fraction depends on the length
		# this depends on both feature orientation and the 
		# relative position requested
		my $relative_start = int( ($length * $fstart) + 0.5);
		my $relative_stop  = int( ($length * $fstop) + 0.5);
		my ($start, $stop);
		
		# check whether length exceeds the limit
		if ($length >= $limit) {
			# length exceeds our minimum limit
			# we can take a fractional length
			
			# adjust relative to the 5' position 
			if ($position == 5) {
			
				if ($strand >= 0) {
					# feature is on the forward, watson strand
					$start = $main_data_ref->{'data_table'}->[$row][$start_i] + 
						$relative_start;
					$stop  = $main_data_ref->{'data_table'}->[$row][$start_i] + 
						$relative_stop;
				}
				elsif ($strand < 0) {
					# feature is on the reverse, crick strand
					$start = $main_data_ref->{'data_table'}->[$row][$stop_i] - 
						$relative_stop;
					$stop  = $main_data_ref->{'data_table'}->[$row][$stop_i] - 
						$relative_start;
				}
			}
			
			# adjust relative to the 3' position
			elsif ($position == 3) {
				
				if ($strand >= 0) {
					# feature is on the forward, watson strand
					$start = $main_data_ref->{'data_table'}->[$row][$stop_i] + 
						$relative_start;
					$stop  = $main_data_ref->{'data_table'}->[$row][$stop_i] + 
						$relative_stop;
				}
				elsif ($strand < 0) {
					# feature is on the reverse, crick strand
					$start = $main_data_ref->{'data_table'}->[$row][$start_i] - 
						$relative_stop;
					$stop  = $main_data_ref->{'data_table'}->[$row][$start_i] - 
						$relative_start;
				}
			}
			
			# adjust relative to the midpoint position
			elsif ($position == 4) {
				# determine the midpoint position
				my $middle = int(
					# calculate same way regardless of strand
					(
						$main_data_ref->{'data_table'}->[$row][$stop_i] +
						$main_data_ref->{'data_table'}->[$row][$start_i] 
					) / 2
				);
				if ($strand >= 0) {
					# feature is on the forward, watson strand
					$start = $middle + $relative_start;
					$stop  = $middle + $relative_stop;
				}
				elsif ($strand < 0) {
					# feature is on the reverse, crick strand
					$start = $middle - $relative_stop;
					$stop  = $middle - $relative_start;
				}
			}
		}
		
		# length doesn't exceed minimum limit
		else {
			# simply take the whole fragment
			$start = $main_data_ref->{'data_table'}->[$row][$start_i];
			$stop  = $main_data_ref->{'data_table'}->[$row][$stop_i];
		}
		
		# get region score
		my $score = get_chromo_region_score(
				'db'        => $ddb,
				'dataset'   => $dataset,
				'chromo'    => $main_data_ref->{'data_table'}->[$row][$chromo_i],
				'start'     => $start,
				'stop'      => $stop,
				'strand'    => $strand,
				'value'     => $value_type,
				'method'    => $method,
				'log'       => $main_data_ref->{$index}{'log2'},
				'stranded'  => $stranded,
		);
		unless (defined $score) {
			# this should return a value, otherwise we record an internal null
			$score = '.';
		}
		$main_data_ref->{'data_table'}->[$row][$index] = $score;
	}
}



sub get_feature_dataset {
	
	# get passed arguments
	my ($dataset, $index) = @_;
	
	
	# Loop through the list of features
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		
		# get feature from the database
		my $feature = get_feature(
			'db'    => $mdb,
			'id'    => defined $id_i ? 
				$main_data_ref->{'data_table'}->[$row][$id_i] : undef,
			'name'  => defined $name_i ? 
				$main_data_ref->{'data_table'}->[$row][$name_i] : undef,
			'type'  => defined $type_i ? 
				$main_data_ref->{'data_table'}->[$row][$type_i] : undef,
		);
		unless ($feature) {
			# record a null value and move on
			$main_data_ref->{'data_table'}->[$row][$index] = '.';
			next;
		}
		
		# reassign strand value if requested
		if ($set_strand) {
			$feature->strand($main_data_ref->{'data_table'}->[$row][$strand_i]);
		}
		
		
		# get region score
		my $score = get_chromo_region_score(
				'db'        => $ddb,
				'dataset'   => $dataset,
				'chromo'    => $feature->seq_id,
				'start'     => $feature->start,
				'stop'      => $feature->end,
				'strand'    => $feature->strand,
				'value'     => $value_type,
				'method'    => $method,
				'log'       => $main_data_ref->{$index}{'log2'},
				'stranded'  => $stranded,
		);
		unless (defined $score) {
			# this should return a value, otherwise we record an internal null
			$score = '.';
		}
		$main_data_ref->{'data_table'}->[$row][$index] = $score;
	}
}



sub get_extended_feature_dataset {
	
	# get passed arguments
	my ($dataset, $index) = @_;
	
	
	# Loop through the list of features
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		
		# get feature from the database
		my $feature = get_feature(
			'db'    => $mdb,
			'id'    => defined $id_i ? 
				$main_data_ref->{'data_table'}->[$row][$id_i] : undef,
			'name'  => defined $name_i ? 
				$main_data_ref->{'data_table'}->[$row][$name_i] : undef,
			'type'  => defined $type_i ? 
				$main_data_ref->{'data_table'}->[$row][$type_i] : undef,
		);
		unless ($feature) {
			# record a null value and move on
			$main_data_ref->{'data_table'}->[$row][$index] = '.';
			next;
		}
		
		# reassign strand value if requested
		if ($set_strand) {
			$feature->strand($main_data_ref->{'data_table'}->[$row][$strand_i]);
		}
		
		# get region score
		my $score = get_chromo_region_score(
				'db'        => $ddb,
				'dataset'   => $dataset,
				'chromo'    => $feature->seq_id,
				'start'     => $feature->start - $extend,
				'stop'      => $feature->end + $extend,
				'strand'    => $feature->strand,
				'value'     => $value_type,
				'method'    => $method,
				'log'       => $main_data_ref->{$index}{'log2'},
				'stranded'  => $stranded,
		);
		unless (defined $score) {
			# this should return a value, otherwise we record an internal null
			$score = '.';
		}
		$main_data_ref->{'data_table'}->[$row][$index] = $score;
	}
}



sub get_adjusted_feature_dataset {
	
	# get passed arguments
	my ($dataset, $index) = @_;
	
	
	# Loop through the list of features
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		
		# get feature from the database
		my $feature = get_feature(
			'db'    => $mdb,
			'id'    => defined $id_i ? 
				$main_data_ref->{'data_table'}->[$row][$id_i] : undef,
			'name'  => defined $name_i ? 
				$main_data_ref->{'data_table'}->[$row][$name_i] : undef,
			'type'  => defined $type_i ? 
				$main_data_ref->{'data_table'}->[$row][$type_i] : undef,
		);
		unless ($feature) {
			# record a null value and move on
			$main_data_ref->{'data_table'}->[$row][$index] = '.';
			next;
		}
		
		# reassign strand value if requested
		if ($set_strand) {
			$feature->strand($main_data_ref->{'data_table'}->[$row][$strand_i]);
		}
		
		
		
		# calculate new relative start and stop positions
		# this depends on both feature orientation and the 
		# relative position requested
		my ($start, $stop);
		
		# adjust relative to the 5' position 
		if ($position == 5) {
		
			if ($feature->strand >= 0) {
				# feature is on the forward, watson strand
				$start = $feature->start + $start_adj;
				$stop  = $feature->start + $stop_adj;
			}
			elsif ($feature->strand < 0) {
				# feature is on the reverse, crick strand
				$start = $feature->end - $stop_adj;
				$stop  = $feature->end - $start_adj;
			}
		}
		
		# adjust relative to the 3' position
		elsif ($position == 3) {
			
			if ($feature->strand >= 0) {
				# feature is on the forward, watson strand
				$start = $feature->end + $start_adj;
				$stop  = $feature->end + $stop_adj;
			}
			elsif ($feature->strand < 0) {
				# feature is on the reverse, crick strand
				$start = $feature->start - $stop_adj;
				$stop  = $feature->start - $start_adj;
			}
		}
		
		# adjust relative to the midpoint position
		elsif ($position == 4) {
			# determine the midpoint position
			my $middle = ($feature->start + int( ($feature->length / 2) + 0.5));
			if ($feature->strand >= 0) {
				# feature is on the forward, watson strand
				$start = $middle + $start_adj;
				$stop  = $middle + $stop_adj;
			}
			elsif ($feature->strand < 0) {
				# feature is on the reverse, crick strand
				$start = $middle - $stop_adj;
				$stop  = $middle - $start_adj;
			}
		}
		
		
		# get region score
		my $score = get_chromo_region_score(
				'db'        => $ddb,
				'dataset'   => $dataset,
				'chromo'    => $feature->seq_id,
				'start'     => $start,
				'stop'      => $stop,
				'strand'    => $feature->strand,
				'value'     => $value_type,
				'method'    => $method,
				'log'       => $main_data_ref->{$index}{'log2'},
				'stranded'  => $stranded,
		);
		unless (defined $score) {
			# this should return a value, otherwise we record an internal null
			$score = '.';
		}
		$main_data_ref->{'data_table'}->[$row][$index] = $score;
	}
}




sub get_fractionated_feature_dataset {
	
	# get passed arguments
	my ($dataset, $index) = @_;
	
	
	# Loop through the list of features
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		
		
		# get feature from the database
		my $feature = get_feature(
			'db'    => $mdb,
			'id'    => defined $id_i ? 
				$main_data_ref->{'data_table'}->[$row][$id_i] : undef,
			'name'  => defined $name_i ? 
				$main_data_ref->{'data_table'}->[$row][$name_i] : undef,
			'type'  => defined $type_i ? 
				$main_data_ref->{'data_table'}->[$row][$type_i] : undef,
		);
		unless ($feature) {
			# record a null value and move on
			$main_data_ref->{'data_table'}->[$row][$index] = '.';
			next;
		}
		
		# reassign strand value if requested
		if ($set_strand) {
			$feature->strand($main_data_ref->{'data_table'}->[$row][$strand_i]);
		}
		
		
		
		# calculate new fractional start and stop positions
		# the fraction depends on the length
		# this depends on both feature orientation and the 
		# relative position requested
		my $relative_start = int( ($feature->length * $fstart) + 0.5);
		my $relative_stop  = int( ($feature->length * $fstop) + 0.5);
		my ($start, $stop);
		
		# use the entire feature length if it doesn't match the minimum limit
		if ($feature->length < $limit) {
			$start = $feature->start;
			$stop  = $feature->end;
		}
		
		# adjust relative to the 5' position 
		elsif ($position == 5) {
		
			if ($feature->strand >= 0) {
				# feature is on the forward, watson strand
				$start = $feature->start + $relative_start;
				$stop  = $feature->start + $relative_stop;
			}
			elsif ($feature->strand < 0) {
				# feature is on the reverse, crick strand
				$start = $feature->end - $relative_stop;
				$stop  = $feature->end - $relative_start;
			}
		}
		
		# adjust relative to the 3' position
		elsif ($position == 3) {
			
			if ($feature->strand >= 0) {
				# feature is on the forward, watson strand
				$start = $feature->end + $relative_start;
				$stop  = $feature->end + $relative_stop;
			}
			elsif ($feature->strand < 0) {
				# feature is on the reverse, crick strand
				$start = $feature->start - $relative_stop;
				$stop  = $feature->start - $relative_start;
			}
		}
		
		# adjust relative to the midpoint position
		elsif ($position == 4) {
			# determine the midpoint position
			my $middle = ($feature->start + int( ($feature->length / 2) + 0.5));
			if ($feature->strand >= 0) {
				# feature is on the forward, watson strand
				$start = $middle + $relative_start;
				$stop  = $middle + $relative_stop;
			}
			elsif ($feature->strand < 0) {
				# feature is on the reverse, crick strand
				$start = $middle - $relative_stop;
				$stop  = $middle - $relative_start;
			}
		}
		
		
		# get region score
		my $score = get_chromo_region_score(
				'db'        => $ddb,
				'dataset'   => $dataset,
				'chromo'    => $feature->seq_id,
				'start'     => $start,
				'stop'      => $stop,
				'strand'    => $feature->strand,
				'value'     => $value_type,
				'method'    => $method,
				'log'       => $main_data_ref->{$index}{'log2'},
				'stranded'  => $stranded,
		);
		unless (defined $score) {
			# this should return a value, otherwise we record an internal null
			$score = '.';
		}
		$main_data_ref->{'data_table'}->[$row][$index] = $score;
	}
}



sub get_subfeature_dataset {
	
	# get passed arguments
	my ($dataset, $index) = @_;
	
	
	# Loop through the list of features
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		
		# get feature from the database
		my $feature = get_feature(
			'db'    => $mdb,
			'id'    => defined $id_i ? 
				$main_data_ref->{'data_table'}->[$row][$id_i] : undef,
			'name'  => defined $name_i ? 
				$main_data_ref->{'data_table'}->[$row][$name_i] : undef,
			'type'  => defined $type_i ? 
				$main_data_ref->{'data_table'}->[$row][$type_i] : undef,
		);
		unless ($feature) {
			# record a null value and move on
			$main_data_ref->{'data_table'}->[$row][$index] = '.';
			next;
		}
		
		# reassign strand value if requested
		if ($set_strand) {
			$feature->strand($main_data_ref->{'data_table'}->[$row][$strand_i]);
		}
		
		
		# Collect and identify the subfeatures
		# we could use either exons or CDS subfeatures
		# but we need to separate them from other types of subfeatures
		my @exons;
		my @cdss;
		foreach my $subfeat ($feature->get_SeqFeatures) {
			if ($subfeat->primary_tag =~ /exon/i) {
				push @exons, $subfeat;
			}
			elsif ($subfeat->primary_tag =~ /utr|untranslated/i) {
				push @cdss, $subfeat;
			}
			elsif ($subfeat->primary_tag =~ /cds/i) {
				push @cdss, $subfeat;
			}
			elsif ($subfeat->primary_tag =~ /rna/i) {
				# an RNA subfeature, keep going down another level
				foreach my $f ($subfeat->get_SeqFeatures) {
					if ($f->primary_tag =~ /exon/i) {
						push @exons, $f;
					}
					elsif ($f->primary_tag =~ /utr|untranslated/i) {
						push @cdss, $f;
					}
					elsif ($f->primary_tag =~ /cds/i) {
						push @cdss, $f;
					}
				}
			}
		}
		
		
		# Determine which subfeatures to collect
		# we prefer to use exons because they are easier
		# if exons are not defined then we'll infer them from CDSs and UTRs
		my @features_to_check;
		if (@exons) {
			@features_to_check = @exons;
		}
		elsif (@cdss) {
			@features_to_check = @cdss;
		}
		else {
			# found neither exons, CDSs, or UTRs
			# possibly because there were no subfeatures
			# in this case we just take the whole thing
			push @features_to_check, $feature;
		}
		
		
		# Collect the subfeature values and summed length
		my %feature_values;
		my %duplicate_check;
		my $gene_length = 0;
		foreach my $subfeat (@features_to_check) {
			# we don't want a single value for each subfeature
			# rather we will collect the raw scores and then combine them
			# ourselves later
			
			# since there may be duplicate subfeatures, for example 
			# shared exons form alternate mRNA transcripts, we need to 
			# avoid these
			next if exists $duplicate_check{ $subfeat->start }{ $subfeat->end };
			# then record this one
			$duplicate_check{ $subfeat->start }{ $subfeat->end } = 1;
			
			# we will use the indexed score function, which returns a 
			# hash of postions and scores, but we'll just take the scores
			# we will do this using genomic coordinates rather than 
			# name and feature type, since we can't guarantee that the 
			# subfeatures are indexed separately from the parent in the 
			# database
			my %pos2scores = get_region_dataset_hash(
				'db'        => $ddb,
				'dataset'   => $dataset,
				'chromo'    => $subfeat->seq_id,
				'start'     => $subfeat->start,
				'stop'      => $subfeat->end,
				'strand'    => $subfeat->strand,
				'value'     => $value_type,
				'stranded'  => $stranded,
				'absolute'  => 1, # do not convert to relative positions
			);
			
			# record the values
			# we need to make sure we don't take duplicate scores, thus inappropriately 
			# inflating the score values for duplicate subfeature, for example with 
			# shared exons between alternate transcripts
			# this hash might get really big
			foreach my $p (keys %pos2scores) {
				# we might overwrite identical positions, but they should be same
				$feature_values{$p} = $pos2scores{$p};
			}
			
			# record the subfeature length
				# potential bug if we have overlapping sub features <sigh>
			$gene_length += $subfeat->length;
		}
	
		
		
		# Calculate the final score
		
		# no data collected!? record null or zero value
		unless (%feature_values) {
			if ($method =~ /sum|count|rpm|rpkm/) {
				$main_data_ref->{'data_table'}->[$row][$index] = 0;
			}
			else {
				$main_data_ref->{'data_table'}->[$row][$index] = '.';
			}
			next;
		}
		
		# convert log2 values if necessary
		if ($main_data_ref->{$index}{'log2'}) {
			foreach (keys %feature_values) {
				$feature_values{$_} = 2 ** $feature_values{$_};
			}
		}
		
		# final score is calculated according to the requested method
		my $parent_score;
		if ($method eq 'median') {
			# take the median value
			$parent_score = median(values %feature_values);
		}
		elsif ($method eq 'mean') {
			# or take the mean value
			$parent_score = mean(values %feature_values);
		} 
		elsif ($method eq 'range') {
			# or take the range value
			# this is 'min-max'
			$parent_score = range(values %feature_values);
		}
		elsif ($method eq 'stddev') {
			# or take the standard deviation value
			# we are using the standard deviation of the population, 
			# since these are the only scores we are considering
			$parent_score = stddevp(values %feature_values);
		}
		elsif ($method eq 'min') {
			# or take the minimum value
			$parent_score = min(values %feature_values);
		}
		elsif ($method eq 'max') {
			# or take the maximum value
			$parent_score = max(values %feature_values);
		}
		elsif ($method eq 'count') {
			# count the number of values
			$parent_score = sum(values %feature_values);
		}
		elsif ($method eq 'sum') {
			# sum the number of values
			$parent_score = sum(values %feature_values);
		}
		elsif ($method eq 'rpm') {
			# calculate reads per per million
			$parent_score = ( sum(values %feature_values) * 10^6 ) / $dataset2sum{$dataset};
		}
		elsif ($method eq 'rpkm') {
			# calculate reads per kb per million
			$parent_score = ( sum(values %feature_values) * 10^9 ) / 
								( $gene_length * $dataset2sum{$dataset});
		}
	
		# convert back to log2 if necessary
		if ($main_data_ref->{$index}{'log2'}) { 
			$parent_score = log($parent_score) / log(2);
		}
			
		# record the final score for the parent feature
		$main_data_ref->{'data_table'}->[$row][$index] = $parent_score;
	}
}


# subroutine to record the metadata for a dataset
sub record_metadata {
	my $dataset = shift;
	
	# generate column name
	my $column_name;
	if ($dataset =~ /^file|http|ftp/) {
		# a specified file
		# we just want the file name, split it from the path
		foreach (split /&/, $dataset) {
			my (undef, undef, $file_name) = File::Spec->splitpath($_);
			if ($column_name) {
				$column_name .= '&' . $file_name;
			}
			else {
				$column_name = $file_name;
			}
		}
	}
	else {
		# a feature type, take as is
		$column_name = $dataset;
	}
	
	
	# determine the log2 status of the dataset or not
	# this will adversely affect the math if not set correctly!
	my $logstatus;
	if (defined $log) {
		# defined globally for all datasets by command line argument
		$logstatus = $log;
	} else {
		if ($dataset =~ /log2/i) {
			# if we're smart we'll encode the log2 status in the dataset name
			$logstatus = 1;
		} else {
			# otherwise assume it is non-log
			# unsafe, but what else to do? we'll put the onus on the user
			$logstatus = 0;
		}
	}
	
	
	# determine new index
	my $new_index = $main_data_ref->{'number_columns'};
	# remember that the index counting is 0-based, so the new index is 
	# essentially number_columns - 1 + 1
	$main_data_ref->{'number_columns'} += 1; # update
	
	# generate new metadata hash for this column
	my %metadata = (
		'name'     => $column_name,
		'dataset'  => $dataset,
		'index'    => $new_index,
		'log2'     => $logstatus,
		'method'   => $method, # global argument for dataset combining
		'strand'   => $stranded, # global argument for strand specificity
		'value'    => $value_type, # global argument for value type
	);
	
	# populate metadata hash as necessary from global arguments
	if (defined $extend) {
		$metadata{'extend'} = $extend;
	}
	if (defined $start_adj) {
		$metadata{'start'} = $start_adj;
	}
	if (defined $stop_adj) {
		$metadata{'stop'} = $stop_adj;
	}
	if (defined $fstart) {
		$metadata{'fstart'} = $fstart;
	}
	if (defined $fstop) {
		$metadata{'fstop'} = $fstop;
	}
	if ($position == 3) {
		$metadata{'relative_position'} = "3'end";
	}
	elsif ($position == 1) {
		$metadata{'relative_position'} = "middle";
	}
	if ($limit) {
		$metadata{'limit'} = $limit;
	}
	if ($subfeature) {
		$metadata{'exons'} = 'yes';
	}
	if ($set_strand) {
		$metadata{'forced_strand'} = 'yes';
	}
	
	# add database name if different
	if (defined $data_database) {
		$metadata{'db'} = $data_database;
	}
	elsif ($main_database ne $main_data_ref->{'db'}) {
		$metadata{'db'} = $main_database;
	}
	
	# place metadata hash into main data structure
	$main_data_ref->{$new_index} = \%metadata;
	
	# set column header name
	$main_data_ref->{'data_table'}->[0][$new_index] = $column_name;
	
	# return the index number
	return $new_index;
}






__END__

=head1 NAME

get_datasets.pl

A program to collect data for a list of features

=head1 SYNOPSIS

get_datasets.pl [--options...] [<filename>]
  
  Options:
  --new
  --in <filename>
  --out filename
  --db <name | filename>
  --ddb <name | filename>
  --feature <type | type:source | alias>, ...
  --data <none | file | type>, ...
  --method [mean | median | stddev | min | max | range | sum | rpm | rpkm]
  --value [score | count | length]
  --log
  --strand [all | sense | antisense]
  --exons
  --extend <integer>
  --start=<integer>
  --stop=<integer>
  --fstart=<decimal>
  --fstop=<decimal>
  --pos [5 | m | 3]
  --limit <integer>
  --win <integer>
  --step <integer>
  --force_strand
  --gz
  --cpu <integer>
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --new

Generate a new table of features. Overrides any specified input file. 
Requires a database and feature to be defined.

=item --in <filename>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. The 
first row should be column headers. Bed files are acceptable, as are 
text files generated with this program or get_features.pl.

=item --out <filename>

Specify the output file name. Required for new feature tables; optional for 
current files. If this is argument is not specified then the input file is 
overwritten.

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

=item --feature genome

Specify the type of feature from which to collect values. This is required 
only for new feature tables. Three types of values may be passed: the 
feature type, feature type and source expressed as 'type:source', or an 
alias to one or more feature types. Aliases are specified in the 
C<biotoolbox.cfg> file and provide a shortcut to a list of one or more 
database features. More than one feature may be included as a 
comma-delimited list (no spaces). 

To collect genomic intervals (or regions) simply specify 'genome' as 
the feature type.

=item --data <type1,type2,type3&type4,...>

=item --data <file1,...>

=item --data none

Provide the name of the dataset to collect the values. Use this argument 
repeatedly for each dataset to be collected. Two or more datasets may be
merged into one by delimiting with an ampersand "&" (no spaces!). If no 
dataset is specified on the command line, then the program will 
interactively present a list of datasets from the database to select. 

The dataset may be a feature type in a BioPerl Bio::DB::SeqFeature::Store 
or Bio::DB::BigWigSet database. Provide either the feature type or 
type:source. The feature may point to another data file whose path is 
stored in the feature's attribute tag (for example a binary 
Bio::Graphics::Wiggle .wib file, a bigWig file, or Bam file), or the 
features' scores may be used in data collection.

Alternatively, the dataset may be a database file, including bigWig (.bw), 
bigBed (.bb), useq (.useq), or Bam alignment (.bam) files. The files may 
be local or remote (specified with a http: or ftp: prefix).

To force the program to simply write out the list of collected features 
without collecting data, provide the dataset name of "none".

=item --method [mean | median | stddev | min | max | range | sum | rpm | rpkm]

Specify the method for combining all of the dataset values within the 
genomic region of the feature. Accepted values include:
  
  - mean        (default)
  - median
  - sum
  - stddev      Standard deviation of the population (within the region)
  - min
  - max
  - range       Returns difference of max and min
  - rpm         Reads Per Million mapped, Bam/BigBed only
  - rpkm        Reads Per Kilobase per Million Mapped, Bam/BigBed only

When collecting data using rpkm, the normalized sum of the reads is 
divided by the length of the feature requested (the Kilobase part in rpkm). 
Note that for mRNA or gene features, this will be the sum of the exon 
lengths, not the gene or mRNA.
  
=item --value [score | count | length]

Optionally specify the type of data value to collect from the dataset or 
data file. Three values are accepted: score, count, or length. The default 
value type is score. Note that some data sources only support certain 
types of data values. Wig and BigWig files only support score and count; 
BigBed and database features support count and length and optionally 
score; Bam files support basepair coverage (score), count (number of 
alignments), and length.

=item --log

Indicate the dataset is (not) in log2 space. The log2 status of the dataset is 
critical for accurately mathematically combining the dataset values in the 
feature's genomic region. It may be determined automatically if the dataset 
name includes the phrase "log2".

=item --strand [all | sense | antisense]

Specify whether stranded data should be collected for each of the 
datasets. Either sense or antisense (relative to the feature) data 
may be collected. Note that strand is not supported with some 
data files, including bigWig files (unless specified through a GFF3 feature 
attribute or Bio::DB::BigWigSet database) and Bam files (score coverage
is not but count is). The default value is 'all', indicating all data 
will be collected.  

=item --exons

Optionally indicate that data should be collected only over the exon 
subfeatures of a gene or transcript, rather than the entire gene. 
Subfeatures with a primary_tag of exon are preferentially taken. If exons 
are not defined, then CDS and UTR subfeatures are used, or the entire 
gene or transcript if no appropriate subfeatures are found. Note that 
the options extend, start, stop, fstart, and fstop are ignored. 
Default is false. 

=item --extend <integer>

Optionally specify the bp extension that will be added to both sides of the 
feature's region.

=item --start=<integer>

=item --stop=<integer>

Optionally specify adjustment values to adjust the region to collect values 
relative to the feature position defined by the --pos option (default is 
the 5' position). A negative value is shifted upstream (5' direction), 
and a positive value is shifted downstream. Adjustments are always made 
relative to the feature's strand. Both options must be applied; one is 
not allowed.

=item --fstart=<number>

=item --fstop=<number>

Optionally specify the fractional start and stop position of the region to 
collect values as a function of the feature's length and relative to the 
specified feature position defined by the --pos option (default is 5'). The 
fraction should be presented as a decimal number, e.g. 0.25. Prefix a 
negative sign to specify an upstream position. Both options must be 
applied; one is not allowed. 

=item --pos [5 | m | 3]

Indicate the relative position of the feature with which the 
data is collected when combined with the "start" and "stop" or "fstart" 
and "fstop" options. Three values are accepted: "5" indicates the 
5' prime end is used, "3" indicates the 3' end is used, and "m" 
indicates the middle of the feature is used. The default is to 
use the 5' end, or the start position of unstranded features. 

=item --limit <integer>

Optionally specify the minimum size limit for subfractionating a feature's 
region. Used in combination with fstart and fstop to prevent taking a 
subregion from a region too small to support it. The default is 1000 bp.

=item --win <integer>

When generating a new genome interval list (feature type 'genome'), 
optionally specify the window size. The default size is defined in the 
configuration file, biotoolbox.cfg. 

=item --step <integer>

Optionally indicate the step size when generating a new list of intervals 
across the genome. The default is equal to the window size.

=item --force_strand

For features that are not inherently stranded (strand value of 0)
or that you want to impose a different strand, set this option when
collecting stranded data. This will reassign the specified strand for
each feature regardless of its original orientation. This requires the
presence of a "strand" column in the input data file. This option only
works with input file lists of database features, not defined genomic
regions (e.g. BED files). Default is false.

=item --gz

Indicate whether the output file should (not) be compressed by gzip. 
If compressed, the extension '.gz' is appended to the filename. If a compressed 
file is opened, the compression status is preserved unless specified otherwise.

=item --cpu <integer>

Specify the number of CPU cores to execute in parallel. This requires 
the installation of Parallel::ForkManager. With support enabled, the 
default is 2. Disable multi-threaded execution by setting to 1. 

=item --version

Print the version number.

=item --help

Display the POD documentation for this program.

=back

=head1 DESCRIPTION

This program will collect dataset values from a variety of sources, including 
features in a BioPerl Bio::DB::SeqFeature::Store database, binary wig files 
(.wib) loaded in a database using Bio::Graphics::Wiggle, bigWig files, 
bigBed files, Bam alignment files, or a Bio::DB::BigWigSet database. 

The values are collected for a list of known database features (genes, 
transcripts, etc.) or genomic regions (defined by chromosome, start, and 
stop). The list may be provided as an input file or generated as a new 
list from a database. Output data files may be reloaded for additional 
data collection.

At each feature or interval, multiple data points within the genomic segment 
are combined statistically and reported as a single value for the feature. 
The method for combining datapoints may be specified; the default method is 
the mean of all datapoints.

The coordinates of the features may be adjusted in numerous ways, including 
specifying a specific relative start and stop, a fractional start and stop, 
an extension to both start and stop, and specifying the relative position 
(5' or 3' or midpoint).

Stranded data may be collected, if the dataset supports stranded information. 
Also, two or more datasets may be combined and treated as one. Note that 
collecting stranded data may significantly slow down data collection.

The output file is a standard tim data formatted file, a tab delimited 
file format with each row a genomic feature and each column a dataset. 
Metadata regarding the datasets are stored in comment lines at the beginning 
of the file. The file may be gzipped.

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

