#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use File::Spec;
use Statistics::Lite qw(sum mean median min max range stddevp);
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
my $VERSION = 1.22;


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

# Assign default values
set_defaults();
my $start_time = time;






### Initialize main data, database, and datasets

# Generate or open Data table
my $Data;
if ($infile) {
	$Data = Bio::ToolBox::Data->new(file => $infile) or 
		die " unable to load input file '$infile'\n";
	printf " Loaded %s features from $infile.\n", format_with_commas( $Data->last_row );
	
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
elsif ($new) {
	# generate a new file
	print " Generating a new feature list from database '$main_database'...\n";
	$Data = Bio::ToolBox::Data->new(
		db      => $main_database,
		feature => $feature,
		win     => $win,
		step    => $step,
	) or die " unable to generate new feature list\n";
}

# update program name
unless ($Data->program eq $0) {
	$Data->program($0);
}



# Open data database
my $ddb;
if (defined $data_database) {
	# specifically defined a data database
	$ddb = open_db_connection($data_database) or 
		die "unable to establish data database connection to $data_database!\n";
}

# Check the datasets
unless ($datasets[0] eq 'none') {
	@datasets = verify_or_request_feature_types(
		'db'      => $ddb || $Data->database,
		'feature' => [ @datasets ],
		'prompt'  => " Enter the dataset(s) or feature type(s) from which \n" . 
					" to collect data. Comma delimited or range is acceptable\n",
	);
}
unless (@datasets) {
	die " No verifiable datasets provided. Check your file path, database, or dataset.\n";
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
		my $success = $Data->save(
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
	while ($cpu > 1 and $Data->last_row / $cpu < 100) {
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
		unless ($value_type =~ m/^(?:score|count|length)$/) {
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
		unless ($stranded =~ m/^(?:all|antisense|sense)$/i) {
			die " unknown strand '$stranded'!";
		}
	} 
	else {
		# default value
		$stranded = 'all'; 
	}
	
	# check the relative position
	if (defined $position) {
		# check the position value
		unless ($position =~ m/^(?:5|4|3|m)$/i) {
			die " Unknown relative position '$position'!\n";
		}
		$position =~ s/m/4/i # change to match internal usage
	}
	else {
		# default position to use the 5' end
		$position = 5;
	}
	
	# check the limit when using fractional start and stop
	if ($fstart and $fstop) {
		# fractional start and stop requested
		unless ($limit) {
			# set a minimum size limit on sub fractionating a feature
			# 1000 bp seems like a reasonable cut off, no?
			$limit = 1000;
		}
		if ($position == 4) {
			die " set position to 5 or 3 only when using fractional start and stop\n";
		}
	}
	
	# check the output file
	unless ($outfile) {
		# overwrite the input file
		$outfile = $infile;
	}
	
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
		$Data->splice_data($i, $cpu);
		
		# re-open database objects to make them clone safe
		# pass second true to avoid cached database objects
		if ($data_database) {
			$ddb = open_db_connection($data_database, 1);
		}
		
		# collect the dataset
		foreach my $dataset (@datasets) {
			unless ($dataset eq 'none') {
				collect_dataset($dataset);
			}
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
	# done
}



sub single_execution {
	
	# collect the datasets
	foreach my $dataset (@datasets) {
		unless ($dataset eq 'none') {
			print " Collecting $method $value_type from dataset '$dataset'...\n";
			collect_dataset($dataset);
		}
		last if $dataset eq 'none';
	}
	
	# write the output file
	# we will rewrite the file after each collection
	# appropriate extensions and compression should be taken care of
	my $success = $Data->save(
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
	my $index = add_new_dataset($dataset);
	
	# check that we have strand data if necessary
	if ($set_strand) {
		unless (defined $Data->strand_column) {
			die " requested to set strand but a strand column was not found!\n";
		}
	}
	
	# we need to determine how we will collect the data from the features
	# using genomic windows or regions, or named features?
	# are we modifying or extending the coordinates of the region or feature?
	
	# Genomic regions
	if ($Data->feature =~ /genome|region|segment|interval/i) {
		
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
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		my $score = $row->get_score(
			'db'        => $ddb,
			'dataset'   => $dataset,
			'value'     => $value_type,
			'method'    => $method,
			'log'       => $Data->metadata($index, 'log2'),
			'stranded'  => $stranded,
		);
		$row->value($index, $score);
	} );
}


sub get_extended_genome_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		my $score = $row->get_score(
			# we need to extend the coordinates
			'start'     => $row->start - $extend,
			'stop'      => $row->end - $extend,
			'db'        => $ddb,
			'dataset'   => $dataset,
			'value'     => $value_type,
			'method'    => $method,
			'log'       => $Data->metadata($index, 'log2'),
			'stranded'  => $stranded,
		);
		$row->value($index, $score);
	} );
}


sub get_adjusted_genome_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		
		# adjust coordinates as requested
		# depends on feature strand and relative position
		my ($start, $stop);
		if ($position == 5 and $row->strand >= 0) { 
			# 5' end of forward strand
			$start = $row->start + $start_adj;
			$stop  = $row->start + $stop_adj;
		}
		elsif ($position == 5 and $row->strand < 0) { 
			# 5' end of reverse strand
			$start = $row->end - $stop_adj;
			$stop  = $row->end - $start_adj;
		}
		elsif ($position == 3 and $row->strand >= 0) { 
			# 3' end of forward strand
			$start = $row->end + $start_adj;
			$stop  = $row->end + $stop_adj;
		}
		elsif ($position == 3 and $row->strand < 0) {
			# 3' end of reverse strand
			$start = $row->start - $stop_adj;
			$stop  = $row->start - $start_adj;
		}
		elsif ($position == 4) {
			# middle position
			my $middle = int( ($row->end - $row->start) / 2);
			if ($row->strand >= 0) {
				$start = $middle + $start_adj;
				$stop  = $middle + $stop_adj;
			}
			else {
				$start = $middle - $stop_adj;
				$stop  = $middle - $start_adj;
			}
		}
		
		# now collect score
		my $score = $row->get_score(
			'start'     => $start,
			'stop'      => $stop,
			'db'        => $ddb,
			'dataset'   => $dataset,
			'value'     => $value_type,
			'method'    => $method,
			'log'       => $Data->metadata($index, 'log2'),
			'stranded'  => $stranded,
		);
		$row->value($index, $score);
	} );
}



sub get_fractionated_genome_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		
		# calculate length
		my $length = $row->end - $row->start + 1;
		
		# calculate new fractional start and stop positions
		# the fraction depends on the length
		# this depends on both feature orientation and the 
		# relative position requested
		my $relative_start = int( ($length * $fstart) + 0.5);
		my $relative_stop  = int( ($length * $fstop) + 0.5);
		my ($start, $stop);
		if ($length >= $limit) {
			# length exceeds our minimum limit
			# we can take a fractional length
			
			if ($position == 5 and $row->strand >= 0) { 
				# 5' end of forward strand
				$start = $row->start + $relative_start;
				$stop  = $row->start + $relative_stop;
			}
			elsif ($position == 5 and $row->strand < 0) { 
				# 5' end of reverse strand
				$start = $row->end - $relative_stop;
				$stop  = $row->end - $relative_start;
			}
			elsif ($position == 3 and $row->strand >= 0) { 
				# 3' end of forward strand
				$start = $row->end + $relative_start;
				$stop  = $row->end + $relative_stop;
			}
			elsif ($position == 3 and $row->strand < 0) {
				# 3' end of reverse strand
				$start = $row->start - $relative_stop;
				$stop  = $row->start - $relative_start;
			}
			# midpoint is not accepted
		}
		else {
			# length doesn't meet minimum limit
			# simply take the whole fragment
			$start = $row->start;
			$stop  = $row->end;
		}
		
		# now collect score
		my $score = $row->get_score(
			'start'     => $start,
			'stop'      => $stop,
			'db'        => $ddb,
			'dataset'   => $dataset,
			'value'     => $value_type,
			'method'    => $method,
			'log'       => $Data->metadata($index, 'log2'),
			'stranded'  => $stranded,
		);
		$row->value($index, $score);
	} );
}



sub get_feature_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		my $score = $row->get_score(
			'db'        => $ddb,
			'dataset'   => $dataset,
			'value'     => $value_type,
			'method'    => $method,
			'log'       => $Data->metadata($index, 'log2'),
			'stranded'  => $stranded,
			'strand'    => $set_strand ? $row->strand : undef,
		);
		$row->value($index, $score);
	} );
}



sub get_extended_feature_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		
		# collect score
		my $score = $row->get_score(
			'start'     => $row->start - $extend,
			'stop'      => $row->end + $extend,
			'db'        => $ddb,
			'dataset'   => $dataset,
			'value'     => $value_type,
			'method'    => $method,
			'log'       => $Data->metadata($index, 'log2'),
			'stranded'  => $stranded,
			'strand'    => $set_strand ? $row->strand : undef,
		);
		$row->value($index, $score);
	} );
}



sub get_adjusted_feature_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		
		# we need to collect the feature to adjust coordinates
		my $feature = $row->feature;
		
		# adjust coordinates as requested
		# depends on feature strand and relative position
		my ($start, $stop);
		if ($position == 5 and $feature->strand >= 0) { 
			# 5' end of forward strand
			$start = $feature->start + $start_adj;
			$stop  = $feature->start + $stop_adj;
		}
		elsif ($position == 5 and $feature->strand < 0) { 
			# 5' end of reverse strand
			$start = $feature->end - $stop_adj;
			$stop  = $feature->end - $start_adj;
		}
		elsif ($position == 3 and $feature->strand >= 0) { 
			# 3' end of forward strand
			$start = $feature->end + $start_adj;
			$stop  = $feature->end + $stop_adj;
		}
		elsif ($position == 3 and $feature->strand < 0) {
			# 3' end of reverse strand
			$start = $feature->start - $stop_adj;
			$stop  = $feature->start - $start_adj;
		}
		elsif ($position == 4) {
			# middle position
			my $middle = int( ($feature->end - $feature->start) / 2);
			if ($feature->strand >= 0) {
				$start = $middle + $start_adj;
				$stop  = $middle + $stop_adj;
			}
			else {
				$start = $middle - $stop_adj;
				$stop  = $middle - $start_adj;
			}
		}
		
		# now collect score
		my $score = $row->get_score(
			'start'     => $start,
			'stop'      => $stop,
			'db'        => $ddb,
			'dataset'   => $dataset,
			'value'     => $value_type,
			'method'    => $method,
			'log'       => $Data->metadata($index, 'log2'),
			'stranded'  => $stranded,
			'strand'    => $set_strand ? $row->strand : undef,
		);
		$row->value($index, $score);
	} );
}




sub get_fractionated_feature_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		
		# we need to collect the feature to adjust coordinates
		my $feature = $row->feature;
		
		# calculate length
		my $length = $feature->length;
		
		# calculate new fractional start and stop positions
		# the fraction depends on the length
		# this depends on both feature orientation and the 
		# relative position requested
		my $relative_start = int( ($length * $fstart) + 0.5);
		my $relative_stop  = int( ($length * $fstop) + 0.5);
		my ($start, $stop);
		if ($length >= $limit) {
			# length exceeds our minimum limit
			# we can take a fractional length
			
			if ($position == 5 and $feature->strand >= 0) { 
				# 5' end of forward strand
				$start = $feature->start + $relative_start;
				$stop  = $feature->start + $relative_stop;
			}
			elsif ($position == 5 and $feature->strand < 0) { 
				# 5' end of reverse strand
				$start = $feature->end - $relative_stop;
				$stop  = $feature->end - $relative_start;
			}
			elsif ($position == 3 and $feature->strand >= 0) { 
				# 3' end of forward strand
				$start = $feature->end + $relative_start;
				$stop  = $feature->end + $relative_stop;
			}
			elsif ($position == 3 and $feature->strand < 0) {
				# 3' end of reverse strand
				$start = $feature->start - $relative_stop;
				$stop  = $feature->start - $relative_start;
			}
			# midpoint is not accepted
		}
		else {
			# length doesn't meet minimum limit
			# simply take the whole fragment
			$start = $feature->start;
			$stop  = $feature->end;
		}
		
		# now collect score
		my $score = $row->get_score(
			'start'     => $start,
			'stop'      => $stop,
			'db'        => $ddb,
			'dataset'   => $dataset,
			'value'     => $value_type,
			'method'    => $method,
			'log'       => $Data->metadata($index, 'log2'),
			'stranded'  => $stranded,
			'strand'    => $set_strand ? $row->strand : undef,
		);
		$row->value($index, $score);
	} );
}



sub get_subfeature_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		
		# we need to collect the feature to adjust coordinates
		my $feature = $row->feature;
		
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
		
		
		# Order and collapse the subfeatures
		# we don't want overlapping features
		my @sorted_features =   map {$_->[1]} 
								sort {$a->[0] <=> $b->[0]} 
								map { [$_->start, $_] } 
								@features_to_check;
		my @subfeatures;
		my $current = shift @sorted_features;
		while ($current) {
			if (@sorted_features) {
				my $next = shift @sorted_features;
				if ($current->overlaps($next) ) {
					# overlapping features, so reassign the end point to that of the next
					$current->end($next->end);
				}
				else {
					# no overlap, keep current, next becomes current
					push @subfeatures, $current;
					$current = $next;
				}
			}
			else {
				# no more features
				push @subfeatures, $current;
				undef $current;
			}
		}
		
		
		# Collect the subfeature values and summed length
		my $gene_length = 0;
		my @scores;
		foreach my $subfeat (@subfeatures) {
			# we will repeatedly call get_score on this $row feature but using
			# the subfeature coordinates. as long as we provide everything, nothing 
			# will be taken from the parent feature
			my $score_ref = $row->get_score(
				'db'        => $ddb,
				'dataset'   => $dataset,
				'chromo'    => $subfeat->seq_id,
				'start'     => $subfeat->start,
				'stop'      => $subfeat->end,
				'strand'    => $set_strand ? $row->strand : $subfeat->strand,
				'stranded'  => $stranded,
				'value'     => $value_type,
				'method'    => 'scores', # special method to return array ref of scores
				'log'       => $Data->metadata($index, 'log2'),
			);
			# we don't want a single value for each subfeature
			# rather we will collect the raw scores and then combine them
			# ourselves later
			
			# record the values
			push @scores, @$score_ref;
			
			# record the subfeature length
			$gene_length += $subfeat->length;
		}

		# Calculate the final score
		unless (@scores) {
			# no data collected!? record null or zero value
			if ($method =~ /sum|count|rpm|rpkm/) {
				$row->value($index, 0);
			}
			else {
				$row->value($index, '.');
			}
			next;
		}
		
		# convert log2 values if necessary
		if ($Data->metadata($index, 'log2')) {
			@scores = map {2 ** $_} @scores;
		}
		
		# final score is calculated according to the requested method
		my $parent_score;
		if ($method eq 'median') {
			# take the median value
			$parent_score = median(@scores);
		}
		elsif ($method eq 'mean') {
			# or take the mean value
			$parent_score = mean(@scores);
		} 
		elsif ($method eq 'range') {
			# or take the range value
			# this is 'min-max'
			$parent_score = range(@scores);
		}
		elsif ($method eq 'stddev') {
			# or take the standard deviation value
			# we are using the standard deviation of the population, 
			# since these are the only scores we are considering
			$parent_score = stddevp(@scores);
		}
		elsif ($method eq 'min') {
			# or take the minimum value
			$parent_score = min(@scores);
		}
		elsif ($method eq 'max') {
			# or take the maximum value
			$parent_score = max(@scores);
		}
		elsif ($method eq 'count') {
			# count the number of values
			$parent_score = sum(@scores);
		}
		elsif ($method eq 'sum') {
			# sum the number of values
			$parent_score = sum(@scores);
		}
		elsif ($method eq 'rpm') {
			# calculate reads per per million
			$parent_score = ( sum(@scores) * 1000000 ) / $dataset2sum{$dataset};
		}
		elsif ($method eq 'rpkm') {
			# calculate reads per kb per million
			$parent_score = ( sum(@scores) * 1000000000 ) / 
								( $gene_length * $dataset2sum{$dataset});
		}
	
		# convert back to log2 if necessary
		if ($Data->metadata($index, 'log2')) { 
			$parent_score = log($parent_score) / LOG2;
		}
		
		# record
		$row->value($index, $parent_score);
	} );
}


# subroutine to record the metadata for a dataset
sub add_new_dataset {
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
	
	# add new column
	my $index = $Data->add_column($column_name);
	
	# update metadata 
	$Data->metadata($index, 'dataset', $dataset);
	$Data->metadata($index, 'method', $method);
	$Data->metadata($index, 'value', $value_type);
	$Data->metadata($index, 'log2', $logstatus);
	$Data->metadata($index, 'strand', $stranded);
	$Data->metadata($index, 'extend', $extend)   if defined $extend;
	$Data->metadata($index, 'start', $start_adj) if defined $start_adj;	
	$Data->metadata($index, 'stop', $stop_adj)   if defined $stop_adj;	
	$Data->metadata($index, 'fstart', $fstart)   if defined $fstart;	
	$Data->metadata($index, 'fstop', $fstop)     if defined $fstop;	
	$Data->metadata($index, 'limit', $limit)     if defined $limit;	
	$Data->metadata($index, 'exons', 'yes')      if $subfeature;	
	$Data->metadata($index, 'forced_strand', 'yes') if $set_strand;	
	if ($position == 3) {
		$Data->metadata($index, 'relative_position', "3'end");	
	}
	elsif ($position == 4) {
		$Data->metadata($index, 'relative_position', 'middle');
	}

	# add database name if different
	if (defined $data_database) {
		$Data->metadata($index, 'db', $data_database);
	}
	
	# return the index number
	return $index;
}






__END__

=head1 NAME

get_datasets.pl

A program to collect data for a list of features

=head1 SYNOPSIS

get_datasets.pl [--options...] [<filename>]
  
  Options for existing files:
  --in <filename>
  
  Options for new files:
  --db <name | filename>
  --feature <type | type:source | alias>, ...
  --win <integer>                                           (500)
  --step <integer>                                          (win)
  
  Options for data collection:
  --ddb <name | filename>
  --data <none | file | type>, ...
  --method [mean|median|stddev|min|max|range|sum|rpm|rpkm]  (mean)
  --value [score|count|length]                              (score)
  --strand [all|sense|antisense]                            (all)
  --force_strand
  --exons
  --log
  
  Adjustments to features:
  --extend <integer>
  --start=<integer>
  --stop=<integer>
  --fstart=<decimal>
  --fstop=<decimal>
  --pos [5|m|3]                                             (5)
  --limit <integer>
  
  General options:
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

Specify the output file name. Required for new feature tables; optional for 
current files. If this is argument is not specified then the input file is 
overwritten.

=item --db <name | filename>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A database is 
required for generating new data files with features. This option may 
skipped when using coordinate information from an input file (e.g. BED 
file), or when using an existing input file with the database indicated 
in the metadata. For more information about using annotation databases, 
see L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. 

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

=item --win <integer>

When generating a new genome interval list (feature type 'genome'), 
optionally specify the window size. The default size is defined in the 
configuration file, biotoolbox.cfg. 

=item --step <integer>

Optionally indicate the step size when generating a new list of intervals 
across the genome. The default is equal to the window size.

=item --ddb <name | filename>

If the data to be collected is from a second database that is separate 
from the annotation database, provide the name of the data database here. 
Typically, a second C<Bio::DB::SeqFeature::Store> or BigWigSet database 
is provided here. 

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

=item --method [mean|median|stddev|min|max|range|sum|rpm|rpkm]

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
  
=item --value [score|count|length]

Optionally specify the type of data value to collect from the dataset or 
data file. Three values are accepted: score, count, or length. The default 
value type is score. Note that some data sources only support certain 
types of data values. Wig and BigWig files only support score and count; 
BigBed and database features support count and length and optionally 
score; Bam files support basepair coverage (score), count (number of 
alignments), and length.

=item --strand [all | sense | antisense]

Specify whether stranded data should be collected for each of the 
datasets. Either sense or antisense (relative to the feature) data 
may be collected. Note that strand is not supported with some 
data files, including bigWig files (unless specified through a GFF3 feature 
attribute or Bio::DB::BigWigSet database) and Bam files (score coverage
is not but count is). The default value is 'all', indicating all data 
will be collected.  

=item --force_strand

For features that are not inherently stranded (strand value of 0)
or that you want to impose a different strand, set this option when
collecting stranded data. This will reassign the specified strand for
each feature regardless of its original orientation. This requires the
presence of a "strand" column in the input data file. This option only
works with input file lists of database features, not defined genomic
regions (e.g. BED files). Default is false.

=item --exons

Optionally indicate that data should be collected only over the exon 
subfeatures of a gene or transcript, rather than the entire gene. 
Subfeatures with a primary_tag of exon are preferentially taken. If exons 
are not defined, then CDS and UTR subfeatures are used, or the entire 
gene or transcript if no appropriate subfeatures are found. Note that 
the options extend, start, stop, fstart, and fstop are ignored. 
Default is false. 

=item --log

Indicate the dataset is (not) in log2 space. The log2 status of the dataset is 
critical for accurately mathematically combining the dataset values in the 
feature's genomic region. It may be determined automatically if the dataset 
name includes the phrase "log2".

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

=item --pos [5|m|3]

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

=head1 EXAMPLES

These are some examples of some common scenarios for collecting data.

=over 4

=item Simple mean scores

You want to collect the mean score from a bigWig file for each feature 
in a BED file of intervals.

  get_datasets.pl --data scores.bw --in input.bed

=item Collect normalized counts

You want to collect normalized read counts from a Bam file of alignments 
for each feature in a BED file.

  get_datasets.pl --data alignments.bam --method rpm --in input.bed

=item Collect stranded RNASeq data

You have stranded RNASeq data, and you would like to determine the 
expression level for all genes in an annotation database.
  
  get_datasets.pl --db annotation --feature gene --data rnaseq.bam \
  --strand sense --exons --method rpkm --out expression.txt

=item Restrict to specific region

You have ChIPSeq enrichment scores in a bigWig file and you now want 
to score just the transcription start site of known transcripts in an 
annotation database. Here you will restrict to 500 bp flanking the TSS.
  
  get_datasets.pl --db annotation --feature mRNA --start=-500 \
  --stop=500 --pos 5 --data scores.bw --out tss_scores.txt

=item Count intervals

You have identified all possible transcription factor binding sites in 
the genome and put them in a bigBed file. Now you want to count how 
many exist in each upstream region of each gene.
  
  get_datasets.pl --db annotation --feature gene --start=-5000 \
  --stop=0 --data tfbs.bb --method sum --value count --out tfbs_sums.txt

=item Many datasets at once

You can place multiple bigWig files in a single directory and treat it 
as a data database, known as a BigWigSet. Each file becomes a database 
feature, and you can interactively choose one or more from which to 
collect. Each dataset is appended to the input file as new column.
  
  get_datasets.pl --ddb /path/to/bigwigset --in input.txt

=item Stranded BigWig data

You can generate stranded RNASeq coverage from a Bam file using the 
BioToolBox script bam2wig.pl, which yields rnaseq_f.bw and rnaseq_r.bw 
files. These are automatically interpreted as stranded datasets in a 
BigWigSet context.
  
  get_datasets.pl --ddb /path/to/rnaseq/bigwigset --strand sense \
  --in input.txt

=item Binned coverage across the genome

You are interested in sequencing depth across the genome to look for 
depleted regions. You count reads in 1 kb intervals across the genome.
  
  get_datasets.pl --db genome.fasta --feature genome --win 1000 \
  --data alignments.bam --value count --method sum --out coverage.txt

=item Middle of feature

You are interested in the maximum score in the central 50% of each 
feature.
  
  get_datasets.pl --fstart=0.25 --fstop=0.75 --data scores.bw --in \
  input.txt

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

