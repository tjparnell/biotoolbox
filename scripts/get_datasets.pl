#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use File::Spec;
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	check_dataset_for_rpm_support
	$BAM_ADAPTER
	$BIG_ADAPTER
);
use Bio::ToolBox::utility;

my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};

my $VERSION = '1.60';


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
	$stranded,
	$subfeature,
	$exon_subfeature,
	$extend,
	$start_adj,
	$stop_adj,
	$fstart,
	$fstop,
	$limit,
	$position,
	$fpkm_method,
	$win,
	$step,
	$set_strand,
	$gz,
	$cpu,
	$help,
	$verbose,
	$print_version,
); 
my @datasets; # an array of names of dataset values to be retrieved

# Command line options
GetOptions( 
	'i|in=s'           => \$infile, # load a pre-existing file
	'new'              => \$new, # generate a new file
	'o|out=s'          => \$outfile, # name of new output file 
	'd|db=s'           => \$main_database, # main or annotation database name
	'D|ddb=s'          => \$data_database, # data database
	'f|feature=s'      => \$feature, # name of genomic feature to analyze
	'm|method=s'       => \$method, # method of collecting & reporting data
	'a|data=s'         => \@datasets, # the list of datasets to collect data from
	't|strand=s'       => \$stranded, # indicate strandedness of data
	'u|subfeature=s'   => \$subfeature, # indicate to restrict to subfeatures
	'exons!'           => \$exon_subfeature, # old parameter
	'x|extend=i'       => \$extend, # extend the size of the genomic feature
	'b|begin|start=i'  => \$start_adj, # adjustment to relative position
	'e|end|stop=i'     => \$stop_adj, # adjustment relative position
	'fstart=f'         => \$fstart, # fractional start position
	'fstop=f'          => \$fstop, # fractional stop position
	'limit=i'          => \$limit, # size limit to fractionate a feature
	'p|pos=s'          => \$position, # set the relative feature position
	'fpkm|rpkm=s'      => \$fpkm_method, # set the fpkm method  
	'win=i'            => \$win, # indicate the size of genomic intervals
	'step=i'           => \$step, # step size for genomic intervals
	'force_strand|set_strand' => \$set_strand, # enforce a specific strand
				# force_strand is preferred option, but respect the old option
	'z|gz!'            => \$gz, # compress output file
	'c|cpu=i'          => \$cpu, # number of execution threads
	'h|help'           => \$help, # request help
	'v|version'        => \$print_version, # print the version
	'bam=s'            => \$BAM_ADAPTER, # explicitly set the bam adapter
	'big=s'            => \$BIG_ADAPTER, # explicitly set the big adapter
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
	print " Biotoolbox script get_datasets.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

# Assign default values
set_defaults();
my $start_time = time;
my $collapsed = 0; # global value to indicate transcripts are collapsed





### Initialize main data, database, and datasets

# Generate or open Data table
my $Data;
if ($infile) {
	$Data = Bio::ToolBox::Data->new(
		file       => $infile, 
		parse      => 1,
		feature    => $feature,
		subfeature => $subfeature,
	) or die " unable to load input file '$infile'\n";
	if ($Data->last_row) {
		printf " Loaded %s features from $infile.\n", format_with_commas( $Data->last_row );
	}
	else {
		die " No features loaded!\n";
	}
	
	# update main database as necessary
	if ($main_database) {
		if ($main_database ne $Data->database) {
			# update with new database
			printf " updating main database name from '%s' to '%s'\n", 
				$Data->database, $main_database;
			$Data->database($main_database);
		}
	}
	else {
		$main_database = $Data->database;
	}
	
	# update feature type as necessary
	if (not defined $Data->feature and not defined $Data->type_column and 
		defined $feature
	) {
		$Data->feature($feature);
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
$Data->program("$0, v $VERSION");



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

# Working with RPKM value datasets
my %dataset2sum; # for genomic totals of reads for rpkm determination
if ($fpkm_method eq 'genome') {
	foreach my $d (@datasets) {
		print " Summing total reads for dataset '$d'...\n";
		my $sum = check_dataset_for_rpm_support($d, $cpu);
		if ($sum) {
			$dataset2sum{$d} = $sum;
			printf "   %s total features\n", format_with_commas($sum);
		}
		else {
			die " $method method requested but not supported for " .
				"dataset $d!\n FPKM is only supported with Bam and BigBed.\n";
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
	while ($cpu > 1 and $Data->last_row / $cpu < 500) {
		# We need at least 500 lines in each fork split to make 
		# it worthwhile to do the split, otherwise, reduce the number of 
		# splits to something more worthwhile
		# dang it! it all depends on what we're collecting from
		# bw mean goes super fast, but bam rpkms are really slow
		$cpu--;
	}
}

# execute data collection in 1 or more processes
if ($cpu > 1) {
	# parallel execution
	# print statements here before we fork, less we have duplicate statements!
	print " Collecting $method scores from datasets @datasets...\n";
	print " Forking into $cpu children for parallel data collection\n";
	parallel_execution();
}
else {
	# single threaded execution
	single_execution();
}


### Finished
# write the output file
# we will rewrite the file after each collection
# appropriate extensions and compression should be taken care of
unless ($outfile) {
	$outfile = $Data->path . $Data->basename;
}
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
printf " Finished in %.1f minutes\n", (time - $start_time)/60;



############# Subroutines ######################################################



### Set default parameters if undefined
sub set_defaults {
	# assign default values
	# these are all global values that could've been assigned on the 
	# command line
	
	# Check for required values
	unless ($main_database) {
		$feature ||= 'gene';
	}
	unless ($infile) {
		if (@ARGV and not $feature) {
			# making an assumption that first unnamed variable is input file
			# but only if no new file feature defined
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
		# database for new files checked below
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
		$cpu ||= 4;
	}
	else {
		# disable cores
		print " disabling parallel CPU execution, no support present\n" if $cpu;
		$cpu = 0;
	}
	
	# check datasets
	if (not @datasets and @ARGV) {
		@datasets = @ARGV;
	}
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
			m/^(?:median|mean|stddev|min|max|range|sum|count|pcount|ncount)$/
		) {
			die " unknown method '$method'!";
		}
	}
	else {
		# set the default to use the mean
		$method = 'mean';
	}
	
	# check fpkm method
	if ($fpkm_method) {
		unless ($fpkm_method eq 'region' or $fpkm_method eq 'genome') {
			die " fpkm option must be one of 'region' or 'genome'! see help\n";
		}
		unless ($method =~ /count/) {
			die " method must be a count if you use the fpkm option!\n";
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
			$limit = 10;
		}
		if ($position == 4) {
			die " set position to 5 or 3 only when using fractional start and stop\n";
		}
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
		elsif (defined $infile) {
			# input could be a gene table in gff or ucsc format to be parsed
			# hope for the best here....
		}
		else {
			die " You must define a database or an appropriate dataset file! see help\n";
		}
	}
	
	# Check subfeatures
	if ($exon_subfeature) {
		# legacy option
		$subfeature = 'exon';
	}
	if ($subfeature and $subfeature !~ /^(?:exon|cds|5p_utr|3p_utr)$/) {
		die " unrecognized subfeature option '$subfeature'! Use exon, cds, 5p_utr or 3p_utr\n";
	} 
}



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
		
		# collapse transcripts if needed
		if ($feature =~ /^gene/i and $subfeature eq 'exon') {
			my $c = $Data->collapse_gene_transcripts;
			if ($c != $Data->last_row) {
				printf " Not all row SeqFeatures could be collapsed, %d failed\n", 
					$Data->last_row - $c;
			}
		}
		
		# collect the dataset
		foreach my $dataset (@datasets) {
			unless ($dataset eq 'none') {
				collect_dataset($dataset);
			}
		}
		
		# write out result
		my $success = $Data->save(
			'filename' => sprintf("$child_base_name.%03s",$i),
			'gz'       => 0, # faster to write without compression
		);
		if ($success) {
			printf " wrote child file $success\n";
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
	
	calculate_fpkm_values() if $fpkm_method;
}



sub single_execution {
	
	# collapse transcripts if needed
	if ($feature =~ /^gene/i and $subfeature eq 'exon') {
		my $c = $Data->collapse_gene_transcripts;
		if ($c != $Data->last_row) {
			printf " Not all row SeqFeatures could be collapsed, %d failed\n", 
				$Data->last_row - $c;
		}
	}
		
	# collect the datasets
	foreach my $dataset (@datasets) {
		unless ($dataset eq 'none') {
			print " Collecting $method scores from dataset '$dataset'...\n";
			collect_dataset($dataset);
		}
		last if $dataset eq 'none';
	}
	calculate_fpkm_values() if $fpkm_method;
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
	
	# collect 
	if (defined $start_adj and defined $stop_adj) {
		# specifically defined relative start and stop positions
		get_adjusted_dataset($dataset, $index);
	}
	elsif (defined $fstart and defined $fstop) {
		# use a subfraction of the region
		get_fractionated_genome_dataset($dataset, $index);
	}
	else {
		# everything else
		get_dataset($dataset, $index);
	}
}


sub get_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		my $score = $row->get_score(
			'db'        => $ddb,
			'dataset'   => $dataset,
			'method'    => $method,
			'stranded'  => $stranded,
			'extend'    => $extend,
			'subfeature' => $subfeature,
		);
		$row->value($index, $score);
	} );
}


sub get_adjusted_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		
		# make sure we work with the represented seqfeature if present
		my $feature = $row->seqfeature || $row;
		
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
			'method'    => $method,
			'stranded'  => $stranded,
		);
		$row->value($index, $score);
	} );
}



sub get_fractionated_dataset {
	my ($dataset, $index) = @_;
	
	# collect the scores from the dataset for this index
	$Data->iterate( sub {
		my $row = shift;
		
		# make sure we work with the represented seqfeature if present
		my $feature = $row->seqfeature || $row;
		
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
			'method'    => $method,
			'stranded'  => $stranded,
		);
		$row->value($index, $score);
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
	
	
	# add new column
	my $index = $Data->add_column($column_name);
	
	# update metadata 
	$Data->metadata($index, 'dataset', $dataset);
	$Data->metadata($index, 'method', $method);
	$Data->metadata($index, 'strand', $stranded);
	$Data->metadata($index, 'extend', $extend)   if defined $extend;
	$Data->metadata($index, 'start', $start_adj) if defined $start_adj;	
	$Data->metadata($index, 'stop', $stop_adj)   if defined $stop_adj;	
	$Data->metadata($index, 'fstart', $fstart)   if defined $fstart;	
	$Data->metadata($index, 'fstop', $fstop)     if defined $fstop;	
	$Data->metadata($index, 'limit', $limit)     if defined $limit;	
	$Data->metadata($index, 'subfeature', $subfeature) if $subfeature;	
	$Data->metadata($index, 'forced_strand', 'yes') if $set_strand;	
	$Data->metadata($index, 'total_reads', $dataset2sum{$dataset}) if 
		exists $dataset2sum{$dataset};
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



sub calculate_fpkm_values {
	print " Calculating FPKM values....\n";
	
	# calculate which indices;
	my @indices = ($Data->number_columns - scalar @datasets) .. $Data->last_column;
	
	# identify the length column as necessary
	my $length_i;
	if ($subfeature) {
		if ($subfeature eq 'exon') {
			$length_i = $Data->find_column('Merged_Transcript_Length');
		}
		elsif ($subfeature eq 'cds') {
			$length_i = $Data->find_column('Transcript_CDS_Length');
		}
		elsif ($subfeature eq '5p_utr') {
			$length_i = $Data->find_column('Transcript_5pUTR_Length');
		}
		elsif ($subfeature eq '3p_utr') {
			$length_i = $Data->find_column('Transcript_3pUTR_Length');
		}
		unless (defined $length_i) {
			$length_i = $Data->add_transcript_length($subfeature);
		}
	}
	
	# work through each collected dataset
	foreach my $index (@indices) {
	
		# determine total number of reads for this dataset
		my $total;
		if ($fpkm_method eq 'genome') {
			# this was calculated earlier
			$total = $dataset2sum{ $Data->metadata($index, 'dataset') };
		}
		else {
			# use the total from the counted regions or genes
			$total = 0;
			$Data->iterate( sub {
				my $row = shift;
				$total += $row->value($index);
			} );
		}
		next unless $total;
	
		# add new column
		my $fpkm_index = $Data->add_column($Data->name($index) . '_FPKM');
		$Data->metadata($fpkm_index, 'dataset', $Data->metadata($index, 'dataset'));
		$Data->metadata($fpkm_index, 'fpkm_method', $fpkm_method);
		$Data->metadata($fpkm_index, 'total_count', $total); 
	
		# calculate and store the fpkms
		$Data->iterate( sub {
			my $row = shift;
			my $length = defined $length_i ? $row->value($length_i) : $row->length;
			$length ||= 1; # why would the length be zero?
			my $fpkm = ($row->value($index) * 1_000_000_000) / ($length * $total);
			$row->value($fpkm_index, $fpkm);
		} );
	}
}


__END__

=head1 NAME

get_datasets.pl

A program to collect data for a list of features

=head1 SYNOPSIS

get_datasets.pl [--options...] <filename>

get_datasets.pl [--options...] --in <filename> <data1> <data2...>
  
  Options for data files:
  -i --in <filename>                  input file: txt bed gff gtf refFlat ucsc
  -o --out <filename>                 optional output file, default overwrite 
  
  Options for new files:
  -d --db <name>                      annotation database: mysql sqlite
  -f --feature <type>                 one or more feature types from db or gff
  
  Options for feature "genome":
  --win <integer>                     size of windows across genome (500 bp)
  --step <integer>                    step size of windows across genome 
  
  Options for data collection:
  -D --ddb <name|file>                data or BigWigSet database
  -a --data <dataset|filename>        data from which to collect: bw bam etc
  -m --method [mean|median|stddev|    statistical method for collecting data
            min|max|range|sum|count|   default mean
            pcount|ncount]
  -t --strand [all|sense|antisense]   strand of data relative to feature (all)
  -u --subfeature [exon|cds|          collect over gene subfeatures 
        5p_utr|3p_utr] 
  --force_strand                      use the specified strand in input file
  --fpkm [region|genome]              convert count data to depth normalized
  
  Adjustments to features:
  -s --extend <integer>               extend the feature in both directions
  -b --begin --start <integer>        adjust relative start coordinate
  -e --end --stop <integer>           adjust relative stop coordinate
  -p --pos [5|m|3]                    define the relative position to adjust
  --fstart=<decimal>                  adjust fractional start
  --fstop=<decimal>                   adjust fractional stop
  --limit <integer>                   minimum size to take fractional window
  
  General options:
  -z --gz                             compress output file
  -c --cpu <integer>                  number of threads, default 4
  -v --version                        print version and exit
  -h --help                           show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 Options for data files

=over 4

=item --in E<lt>filenameE<gt>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. Any tab-delimited text 
file with recognizable headers is supported. Gene annotation file 
formats are also supported, including bed, gtf, gff3, refFlat, and 
UCSC native formats such as gene prediction tables are all supported. 
Gene annotation files will be parsed as sequence features. 
Files may be gzipped compressed.

=item --out E<lt>filenameE<gt>

Specify the output file name. Required for new feature tables; optional for 
current files. If this is argument is not specified then the input file is 
overwritten.

=back

=head2 Options for new files

=over 4

=item --db E<lt>name | filenameE<gt>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A database is 
required for generating new data files with features. This option may 
skipped when using coordinate information from an input file (e.g. BED 
file), or when using an existing input file with the database indicated 
in the metadata.  

=item --feature <type | type:source | alias>,...

Specify the type of feature from which to collect values. This is required 
only for new feature tables. Three types of values may be passed: the 
feature type, feature type and source expressed as 'type:source', or an 
alias to one or more feature types. More than one feature may be included 
as a comma-delimited list (no spaces). 

=back

=head2 Options for feature "genome"

=over 4

=item --feature genome

To collect genomic intervals (or regions) simply specify 'genome' as 
the feature type.

=item --win E<lt>integerE<gt>

When generating a new genome interval list (feature type 'genome'), 
optionally specify the window size.  

=item --step E<lt>integerE<gt>

Optionally indicate the step size when generating a new list of intervals 
across the genome. The default is equal to the window size.

=back

=head2 Options for data collection

=over 4

=item --ddb E<lt>nameE<gt>

If the data to be collected is from a second database that is separate 
from the annotation database, provide the name of the data database here. 
Typically, a second L<Bio::DB::SeqFeature::Store> or BigWigSet database 
is provided here. 

=item --data <type1,type2,type3&type4,...>

=item --data <file1,...>

=item --data none

Provide the name of the dataset to collect the values. Use this argument 
repeatedly for each dataset to be collected. Two or more datasets may be
merged into one by delimiting with an ampersand "&" (no spaces!). If no 
dataset is specified on the command line, then the program will 
interactively present a list of datasets from the database to select. 

The dataset may be a feature type in a BioPerl L<Bio::DB::SeqFeature::Store> 
or L<Bio::DB::BigWigSet> database. Provide either the feature type or 
type:source. The feature may point to another data file whose path is 
stored in the feature's attribute tag (for example a binary 
Bio::Graphics::Wiggle F<.wib> file, a bigWig file, or Bam file), or the 
features' scores may be used in data collection.

Alternatively, the dataset may be a database file, including bigWig (.bw), 
bigBed (.bb), useq (.useq), or Bam alignment (.bam) files. The files may 
be local or remote (specified with a http: or ftp: prefix).

To force the program to simply write out the list of collected features 
without collecting data, provide the dataset name of "none".

=item --method E<lt>textE<gt>

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

=item --strand [all | sense | antisense]

Specify whether stranded data should be collected for each of the 
datasets. Either sense or antisense (relative to the feature) data 
may be collected. Note that strand is not supported with some 
data files, including bigWig files (unless specified through a GFF3 feature 
attribute or BigWigSet database) and Bam files (score coverage
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

=item --subfeature [ exon | cds | 5p_utr | 3p_utr ]

Optionally specify the type of subfeature to collect from, rather than 
the entire gene. If the parent feature is gene and the subfeature is exon, 
then all transcripts of the gene will be collapsed. The other subfeatures 
(cds, 5p_utr, and 3p_utr) will not work with gene features but only with 
coding mRNA transcripts. Note that the options extend, start, stop, fstart, 
and fstop are ignored. Default is null. 

=item --exons

Legacy option for specifying --subfeature exon.

=item --fpkm [region|genome]

Optionally indicate that counts should be converted to Fragments Per 
Kilobase per Million mapped (FPKM). This is a method for normalizing 
sequence read depth and is used with Bam (or optionally bigBed) files. 
Two methods exist for normalizing: 

=over 4

=item region 

Uses the sum of counts over all input regions examined and ignores 
non-counted reads

=item genome

Uses the sum of all reads across the genome, regardless of whether 
it was counted in an input region or not

=back

The region method is best used with RNASeq data and a complete gene 
annotation table. The genome method is best used with partial annotation 
tables or other Seq types, such as ChIPSeq. This option can only be used 
with one of the count methods (count, ncount, pcount). The FPKM values 
are appended as additional columns in the output table.

=back

=head2 Adjustments to features

=over 4

=item --extend E<lt>integerE<gt>

Optionally specify the bp extension that will be added to both sides of the 
feature's region.

=item --start E<lt>integerE<gt>

=item --stop E<lt>integerE<gt>

=item --begin E<lt>integerE<gt>

=item --end E<lt>integerE<gt>

Optionally specify adjustment values to adjust the region to collect values 
relative to the feature position defined by the C<--pos> option (default is 
the 5' position). A negative value is shifted upstream (5' direction), 
and a positive value is shifted downstream. Adjustments are always made 
relative to the feature's strand. Both options must be applied; one is 
not allowed.

=item --pos [5|m|3]

Indicate the relative position of the feature with which the 
data is collected when combined with the "start" and "stop" or "fstart" 
and "fstop" options. Three values are accepted: "5" indicates the 
5' prime end is used, "3" indicates the 3' end is used, and "m" 
indicates the middle of the feature is used. The default is to 
use the 5' end, or the start position of unstranded features. 

=item --fstart=<number>

=item --fstop=<number>

Optionally specify the fractional start and stop position of the region to 
collect values as a function of the feature's length and relative to the 
specified feature position defined by the C<--pos> option (default is 5'). The 
fraction should be presented as a decimal number, e.g. 0.25. Prefix a 
negative sign to specify an upstream position. Both options must be 
applied; one is not allowed. 

=item --limit E<lt>integerE<gt>

Optionally specify the minimum size limit for subfractionating a feature's 
region. Used in combination with fstart and fstop to prevent taking a 
subregion from a region too small to support it. The default is 10 bp.

=back

=head2 General options

=over 4

=item --gz

Indicate whether the output file should (not) be compressed by gzip. 
If compressed, the extension F<.gz> is appended to the filename. If a compressed 
file is opened, the compression status is preserved unless specified otherwise.

=item --cpu E<lt>integerE<gt>

Specify the number of CPU cores to execute in parallel. This requires 
the installation of L<Parallel::ForkManager>. With support enabled, the 
default is 4. Disable multi-threaded execution by setting to 1. 

=item --version

Print the version number.

=item --help

Display the POD documentation for this program.

=back

=head1 DESCRIPTION

This program will collect dataset values from a variety of sources, including 
features in a BioPerl Bio::DB::SeqFeature::Store database, binary wig files 
F<.wib> loaded in a database using Bio::Graphics::Wiggle, bigWig files, 
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
  --stop=0 --data tfbs.bb --method count --out tfbs_sums.txt

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
  --data alignments.bam --method count --out coverage.txt

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
it under the terms of the Artistic License 2.0.  

