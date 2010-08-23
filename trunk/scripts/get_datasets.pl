#!/usr/bin/perl -w

# A script to grab the microarray data from the gbrowse db for promoters

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_db_helper qw(
	open_db_connection
	get_dataset_list 
	validate_dataset_list 
	get_new_feature_list 
	get_new_genome_list 
	get_feature_dataset 
	get_genome_dataset 
);
use tim_file_helper;
#use Data::Dumper;


print "\n A program to collect feature data from the database\n\n";


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
my (  # command line variables
	$infile,
	$new,
	$outfile,
	$database,
	$feature,
	$method,
	$log,
	$strand,
	$extend,
	$start,
	$stop,
	$fstart,
	$fstop,
	$limit,
	$position,
	$win,
	$step,
	$dubious,
	$gz,
	$help,
	$doc,
); 
my @datasets; # an array of names of dataset values to be retrieved

# Command line options
GetOptions( 
	'in=s'       => \$infile, # load a pre-existing file
	'new'        => \$new, # generate a new file
	'out=s'      => \$outfile, # name of new output file 
	'db=s'       => \$database, # database name
	'feature=s'  => \$feature, # name of genomic feature to analyze
	'method=s'   => \$method, # method of collecting & reporting data
	'dataset=s'  => \@datasets, # the list of datasets to collect data from
	'log!'       => \$log, # dataset is in log2 space
	'strand=s'   => \$strand, # indicate strandedness of data
	'extend=i'   => \$extend, # extend the size of the genomic feature
	'start=s'    => \$start, # adjustment to relative start position
	'stop=s'     => \$stop, # adjustment relative stop position
	'fstart=f'   => \$fstart, # fractional start position
	'fstop=f'    => \$fstop, # fractional stop position
	'limit=i'    => \$limit, # size limit to fractionate a feature
	'pos=i'      => \$position, # set the relative feature position
	'win=i'      => \$win, # indicate the size of genomic intervals
	'step=i'     => \$step, # step size for genomic intervals
	'dubious!'   => \$dubious, # collect dubious genes
	'gz!'        => \$gz, # compress output file
	'help'       => \$help, # request help
	'doc'        => \$doc, # print POD documentation
) 
# if errors in the getoptions
or pod2usage( {
	'-message' => 'unknown argument(s) given!',
	'-verbose' => 0, 
	'-exitval' => 1,
} );

# print help if requested
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}


# Check for required values
my @errors; # collect the error messages here
unless ($infile or $new) {
	push @errors, " You must specify a file or request a new feature table!";
}
if ($new) {
	unless ($outfile) {
		push @errors, " You must define an output filename";
	}
	unless ($database) {
		push @errors,  " You must define a database!";
	}
	unless ($feature) {
		push @errors,  " You must define a genomic feature!";
	}
}
if (@errors) {
	pod2usage( {
		'-message' => join("\n", @errors) . "\n",
		'-exitval' => 2,
		'-verbose' => 0,
	} );
}




# Set default parameters if undefined
if ($method) {
	# check the method that was defined on the command line
	unless ( _verify_method($method) ) {
		die " unknown method '$method'!";
	}
}
else {
	# set the default to use the mean
	$method = 'mean';
}

if (defined $strand) { 
	# check the strand request that was defined on the command line
	unless (
		$strand eq 'sense' or 
		$strand eq 'antisense' or 
		$strand eq 'none'
	) {
		die " unknown strand '$strand'!";
	}
} 
else {
	# default value
	$strand = 'none'; 
}

if ($fstart and $fstop) {
	# fractional start and stop requested
	unless ($limit) {
		# set a minimum size limit on sub fractionating a feature
		# 1000 bp seems like a reasonable cut off, no?
		$limit = 1000;
	}
}

if (defined $position) {
	# check the position value
	unless (
		$position == 5 or
		$position == 3 or
		$position == 1 or
		$position eq 'm'
	) {
		die " Unknown relative position '$position'!\n";
	}
	if ($position eq 'm') {$position = 1} # change to match internal usage
}
else {
	# default position to use the 5' end
	$position = 5;
}

unless ($outfile) {
	# overwrite the input file
	$outfile = $infile;
}










########## Main #############

### set main variables
my $start_time = time;
my $main_data_ref; # this will be a reference to the primary data hash



### Collect the feature list and populate the 
if ($new) { 
	# generating new input file from the database
	print " Generating a new feature list from database '$database'...\n";
	$main_data_ref = generate_new_list();
	unless ($main_data_ref) {
		die "no new feature list generated!";
	}
	
	# set the current program
	$main_data_ref->{'program'} = $0;
	
} else { 
	# open a file containing the data
	print " Loading feature list from file '$infile'...\n";
	$main_data_ref = load_tim_data_file($infile);
	unless ($main_data_ref) {
		die "no file data loaded!";
	}
	
	# set missing arguments with values from loaded file
	unless ($database) {
		$database = $main_data_ref->{'db'};
	}
	unless ($feature) {
		$feature = $main_data_ref->{'feature'};
	}
}
# 
# # Data Dump new data
# open DATADUMP, ">new_data_dump.txt";
# print DATADUMP Dumper($main_data_ref);
# close DATADUMP;


### Establish db connection
my $db = open_db_connection($database);
unless ($db) {
	die " unable to establish database connection!\n";
}



### Get the data sets

# set start, stop for tRNA
if ($feature eq 'trna') {
	# since tRNAs are so small, we'll automatically set the region to include
	# 150 bp of flanking DNA to ensure we'll get microarray coverage
	unless ($start) {$start = -150}
	unless ($stop) {$stop = 250}
}

# collect the datasets
if (@datasets) {
	# Datasets were requested on the command line, automatic execution
	unless ($datasets[0] eq 'none') {
		# we can skip all data collection if the first dataset is 'none'
		auto_validate_and_collect_datasets();
	}

} else {
	# Interactively request datasets
	interactive_collect_datasets();
}

# # Data Dump new data
# open DATADUMP, ">post_data_dump.txt";
# print DATADUMP Dumper($main_data_ref);
# close DATADUMP;




### Output the data
# we will write a tim data file
# appropriate extensions and compression should be taken care of
my $success = write_tim_data_file( {
	'data'     => $main_data_ref,
	'filename' => $outfile,
	'gz'       => $gz,
} );
if ($success) {
	print " wrote file '$success'\n";
}
else {
	# failure! the subroutine will have printed error messages
	print " unable to write file!\n";
}

my $stop_time = sprintf "%.1f", (time - $start_time)/60;
print " Completed in $stop_time minutes\n";
	


############# Subroutines ###################

## Generate a new gene/feature list
sub generate_new_list {
	
	# generate feature list based on type of feature
	if ($feature eq 'genome') { 
		# generate a list of genomic intervals
		
		# if window and step size were not set at execution from the 
		# command line, then we'll be using the defaults hard encoded
		# in the module subroutine, which should be window of 500 bp 
		# and step size = window size
		
		# generate the list
		return get_new_genome_list( {
			'db'       => $database, 
			'win'      => $win,
			'step'     => $step,
		} );
	} 
	
	else { 
		# everything else works off a feature list
		
		# generate the gene list
		return get_new_feature_list( {
			'db'        => $database,
			'features'  => $feature,
			'dubious'   => $dubious,
		} );
	}
}



# Automatic dataset collection
sub auto_validate_and_collect_datasets {
	# first validate the proffered datasets
	my $bad_dataset = validate_dataset_list($db, @datasets);
	
	if ($bad_dataset) {
		# there is at least one bad dataset
		# just refuse to even look up good data data sets and skip the bad
		# it's easier that way, but not very fair
		die " The following datasets are not valid: $bad_dataset\n";
	} 
	else {
		# all the datasets appear to be good (nothing returned from validator)
		# perform the data lookups
		foreach (@datasets) {
			submit_dataset_request($_); # process the request
		}
	}
}



# Interactively request datasets and collect
sub interactive_collect_datasets {
	
	# first retreive the list of microarray data sets from the database
	my %number_of_dataset = &get_dataset_list($db);
	
	# then present the list to the user 
	print 
		"\n These are the microarray data sets in the database '$database':\n";
	foreach (sort {$a <=> $b} keys %number_of_dataset) {
		# print out the list of microarray data sets
		print "  $_\t$number_of_dataset{$_}\n"; 
	}
	# print instructions
	print 
		" Enter the number(s) of the data set(s) you would like to retrieve\n" . 
		" as a comma-delimited list and/or range (x-y).\n Merge multiple datasets into one using an '&'.\n"
	;
	
	# get an answer
	my $answer = <STDIN>;
	chomp $answer;
	
	# process the user's request for the data set
	$answer =~ s/\s+//g;
	@datasets = _parse_list($answer); # split the user requests
	
	EACH_REQUEST:
	foreach my $data_request (@datasets) {
		# process each request depending on whether datasets need merging
		if ($data_request =~ m/&/) {
			# multiple datasets are to be merged
			
			# convert numbers to names
			my @merged_requests;
			foreach (split /&/, $data_request) {
				# check each number for validity
				# and convert the number to a dataset name
				if (exists $number_of_dataset{$_} ) {
					push @merged_requests, $number_of_dataset{$_};
				}
				# otherwise complain to user
				else {
					warn "unknown number '$_' in request '$data_request'!";
					next EACH_REQUEST;
				}
			}
			
			# perform the dataset request if we reach this point
			submit_dataset_request( join('&', @merged_requests) );
		}
		
		else {
			# only one dataset requested this time
			
			# check the number for validity
			# and convert the number to a dataset name
			if (exists $number_of_dataset{$data_request} ) {
				submit_dataset_request( $number_of_dataset{$data_request} );
			}
			# otherwise complain to user
			else {
				warn "unknown number in request '$data_request'!";
				next EACH_REQUEST;
			}
		}
	}
}



# Subroutine to perform the dataset request
sub submit_dataset_request {
	my $dataset = shift;
	
	# we need to know whether we are working with a log2 dataset or not
	# as it will adversely affect the math!
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
	
	print " Collecting $method data for dataset $dataset...\n";
	
	if ($feature eq 'genome') {
		# collecting a dataset for genomic intervals
		if ( get_genome_dataset( {
				'data'     => $main_data_ref->{'data_table'},
				'db'       => $db,
				'dataset'  => $dataset,
				'method'   => $method,
				'strand'   => $strand,
				'log'      => $logstatus,
			} )
		) {
			print " in ";
			printf "%.1f", ( (time - $start_time)/60 ); 
			print " minutes\n";
			
			# record metadata
			_record_metadata($dataset, $logstatus);
		}
		else {
			# error messages should have printed from the subroutine
			print " failed!\n";
		}
		
		
	}
	
	else {
		# collecting a dataset for list of features
		if ( get_feature_dataset( {
				'data'      => $main_data_ref->{'data_table'},
				'db'        => $db,
				'dataset'   => $dataset,
				'method'    => $method,
				'strand'    => $strand,
				'log'       => $logstatus,
				'extend'    => $extend,
				'start'     => $start,
				'stop'      => $stop,
				'position'  => $position,
				'fstart'    => $fstart,
				'fstop'     => $stop,
				'limit'     => $limit,
			} )
		) {
			print " in ";
			printf "%.1f", ( (time - $start_time)/60 ); 
			print " minutes\n";
			
			# record metadata
			_record_metadata($dataset, $logstatus);
		}
		else {
			# error messages should have printed from the subroutine
			print " failed!\n";
		}
	}
	
}


# subroutine to record the metadata for a dataset
sub _record_metadata {
	my ($dataset, $logstatus) = @_;
	
	# determine new index
	my $new_index = $main_data_ref->{'number_columns'};
	# remember that the index counting is 0-based, so the new index is 
	# essentially number_columns - 1 + 1
	$main_data_ref->{'number_columns'} += 1; # update
	
	# generate new metadata hash for this column
	my %metadata = (
		'name'     => $dataset,
		'dataset'  => $dataset,
		'index'    => $new_index,
		'log2'     => $logstatus,
		'method'   => $method, # global argument for dataset combining
		'strand'   => $strand, # global argument for strand specificity
	);
	
	# populate metadata hash as necessary from global arguments
	if ($extend) {
		$metadata{'extend'} = $extend;
	}
	if ($start) {
		$metadata{'start'} = $start;
	}
	if ($stop) {
		$metadata{'stop'} = $stop;
	}
	if ($fstart) {
		$metadata{'fstart'} = $fstart;
	}
	if ($fstop) {
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
	# add database name if different
	if ($database ne $main_data_ref->{'db'}) {
		$metadata{'db'} = $database;
	}
	
	# place metadata hash into main data structure
	$main_data_ref->{$new_index} = \%metadata;
}



# subroutine to verify methods
sub _verify_method {
	my $method = shift;
	# generate a simple hash of acceptable methods
	# additional methods may be added here
	my %acceptable = (
		'median'     => 1,
		'mean'       => 1,
		'stddev'     => 1,
		'min'        => 1,
		'max'        => 1, 
		'range'      => 1,
		'enumerate'  => 1,
	);
	# return true if this method exists
	return $acceptable{$method};
}


# subroutine to parse a list
sub _parse_list {
	# this subroutine will parse a string into an array
	# it is designed for a string of numbers delimited by commas
	# a range of numbers may be specified using a dash
	# hence 1,2,5-7 would become an array of 1,2,5,6,7
	
	my $string = shift;
	$string =~ s/\s+//g; 
	my @list;
	foreach (split /,/, $string) {
		# check for a range
		if (/\-/) { 
			my ($start, $stop) = split /\-/;
			# add each item in the range to the list
			for (my $i = $start; $i <= $stop; $i++) {
				push @list, $i;
			}
			next;
		} 
		else {
			# ordinary number
			push @list, $_;
		}
	}
	return @list;
}






__END__

=head1 NAME

get_datasets.pl

A script to collect genomic datasets from a Bioperl SeqFeature::Store db.

=head1 SYNOPSIS

get_datasets.pl [--new | --in <filename>] [--options...]
  
  --out filename
  --db name
  --feature [gene, orf, rna, trna, genome, cen, ars]
  --method [mean, median, stddev, min, max, range, enumerate]
  --dataset name
  --(no)log
  --strand [none, sense, antisense]
  --extend integer
  --start integer
  --stop integer
  --fstart decimal
  --fstop decimal
  --limit integer
  --pos [5 | 3 | m]
  --win integer
  --step integer
  --(no)dubious
  --(no)gz
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --new

Generate a new table of features. Required unless an input file is specified.

=item --in

Specify the file name of a previously generated feature dataset.
Required unless "--new" is specified. It should be in the tim data
format that is generated by this program and others, although other
tab-delimited text data formats may be usable. See the file
description below and in L<tim_db_helper.pm>.

=item --out

Specify the output file name. Required for new feature tables; optional for 
current files. If this is argument is not specified then the input file is 
overwritten.

=item --db

Specify the name of the BioPerl gff database to use as source. This is required 
for new feature data files. For pre-existing input data files, this argument 
is optional, but if given it overrides the database listed in the file; this 
is useful for collecting data from multiple databases.

=item --feature

Specify the type of feature from which to collect values. This is required 
for new feature tables. Accepted values include:
  
  - gene       All genes, including ORFs, snRNAs, snoRNAs, ncRNAs
  - orf        Only ORF genes
  - rna        Only ncRNAs, snRNAs, and snoRNAs
  - trna       Only tRNAs. Note that, unless specified otherwise,
               this will also automatically set start and stop 
               options to -150 and +250 to encompass the gene
  - genome     Take the entire genome in windows (default is 500 bp,
               or as defined by the --win option)
  - cen        Take all centromeres
  - ars        Take all Autonomously Replicating Sequences
  - tim_transcript          Transcripts determined by tim using the 
                            script 'map_transcripts.pl' using 
  - perrochi_transcript     Transcripts 
  - miura_transcript        Transcripts 
  - nagalakshmi_transcript  Transcripts 
  - method:source  For custom features, the GFF type of the features 
               may be provided as a comma delimited list, where the 
               type corresponds to the combined GFF method:source.

Finally, one or more actual Bio::DB::SeqFeature::Store type strings
comprised of the GFF's method:source elements may be passed, for
example 'gene:SGD', delimited by commas (no spaces). These will not be
interpreted but used directly.

=item --method

Specify the method for combining all of the dataset values within the 
genomic region of the feature. Accepted values include:
  
  - mean        (default)
  - median
  - stddev      Standard deviation of the population (within the region)
  - min
  - max
  - range       Returns 'min-max'
  - enumerate   Counts the number of features or dataset points
  

=item --dataset

Provide the name of the dataset to collect the values. Use this argument 
repeatedly for each dataset to be collected. Two or more datasets may be
merged into one by delimiting with an "&" (no spaces!). If the dataset is not 
specified on the command line, then the program will interactively present a 
list of datasets from the database to select. To force the program to 
simply write out the list of collected features without collecting data, 
provide the dataset name of "none".

=item --(no)log

Indicate the dataset is (not) in log2 space. The log2 status of the dataset is 
critical for accurately mathematically combining the dataset values in the 
feature's genomic region. It may be determined automatically if the dataset 
name includes the phrase "log2".

=item --strand

Indicate whether stranded microarray data should be collected. Accepted values 
include:
  
  - none (default)
  - sense
  - antisense
  

=item --extend

Optionally specify the bp extension that will be added to both sides of the 
feature's region.

=item --start <integer>

Optionally specify the start position of the region to collect values relative 
to the feature start. Prefix a negative sign to specify 
an upstream position. Specify a negative value on the command line with an 
equals sign, e.g. "--start=-300'. Must be combined with "--stop".

=item --stop <integer>

Optionally specify the stop position of the region to collect values relative 
to the feature start. Must be combined with "--start".

=item --fstart <number>

Optionally specify the fractional start position of the region to collect 
values relative to the feature start (or end if specified). The fraction is 
based on the feature's region length. The fraction should be presented as a 
decimal number, e.g. 0.25. Prefix a negative sign to specify an upstream 
position. Must be combined with "--fstop".

=item --fstop <number>

Optionally specify the fractional stop position of the region to collect 
values relative to the feature start (or end if specified). The fraction is 
based on the feature's region length. The fraction should be presented as a 
decimal number, e.g. 0.25. A value > 1 would include the region downstream 
of the feature. Must be combined with "--fstart".

=item --limit <integer>

Optionally specify the minimum size limit for subfractionating a feature's 
region. Used in combination with fstart and fstop to prevent taking a 
subregion from a region too small to support it. The default is 1000 bp.

=item --pos [5 | 3 | m]

Indicate the relative position of the feature with which the 
data is collected when combined with the "start" and "stop" or "fstart" 
and "fstop" options. Three values are accepted: "5" indicates the 
5' prime end is used, "3" indicates the 3' end is used, and "m" or "1" 
indicates the middle of the feature is used. The default is to 
use the 5' end, or the start position of unstranded features. 

=item --win <integer>

When generating a new genome interval list, optionally specify the window 
size. The default size is 500 bp.

=item --step <integer>

Optionally indicate the step size when generating a new list of intervals 
across the genome. The default is equal to the window size.

=item --(no)dubious

When generating a new feature list of ORFs, indicate that ORFs flagged 
'dubious' should (not) be included. The GFF database must include the attribute 
'Qualifier'. The default behavior is to not include dubious genes.

=item --(no)gz

Indicate whether the output file should (not) be compressed by gzip. 
If compressed, the extension '.gz' is appended to the filename. If a compressed 
file is opened, the compression status is preserved unless specified otherwise.

=item --help

Display the complete POD documentation for this program.

=back

=head1 DESCRIPTION

This program will collect dataset values from a Bioperl SeqFeature::Store 
formatted database.
It will generate a list of genomic features or intervals, collect values from 
one or more datasets for each feature, and write out a data file. Previous 
data files may be re-loaded and additional datasets added to the list of 
features.

At each feature or interval, multiple data points within the genomic segment 
are combined mathematically and reported as a single value for the feature. 
The method for combining datapoints may be specified; the default method is 
the mean of all datapoints.

Datasets in the database to be collected may be chosen interactively from a list 
for convenience. Alternatively, the datasets may be specified upon execution 
at the command line. The program is designed to be run completely without 
interaction, allowing for convenient usage within a shell script.

The data file that is written is a simple tab-delimited text file, with each 
row representing the feature or interval, and each column representing 
either feature identification information or the dataset values. Features 
without a dataset value have a 'null' value, represented by a period (.). 
Metadata about each dataset are stored in comment lines, prefixed by an 
octothorpe (#) at the beginning of the file. These metadata lines are read 
and used by the program when the file is opened for adding subsequent  
datasets. The data file may be compressed using gzip for compact storage.


=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112


=head1 TODO




