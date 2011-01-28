#!/usr/bin/perl

# A script to collect data from the gbrowse db for genomic features

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	parse_list
	find_column_index
);
use tim_db_helper qw(
	open_db_connection
	get_dataset_list 
	validate_dataset_list 
	get_new_feature_list 
	get_new_genome_list 
	get_feature_dataset 
	get_genome_dataset 
);
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
);
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
my (  
	$infile,
	$new,
	$outfile,
	$database,
	$feature,
	$method,
	$value_type,
	$log,
	$strand,
	$subfeature,
	$extend,
	$start,
	$stop,
	$fstart,
	$fstop,
	$limit,
	$position,
	$win,
	$step,
	$set_strand,
	$gz,
	$help,
	$doc,
); 
my @datasets; # an array of names of dataset values to be retrieved
my @datafiles; # a

# Command line options
GetOptions( 
	'in=s'       => \$infile, # load a pre-existing file
	'new'        => \$new, # generate a new file
	'out=s'      => \$outfile, # name of new output file 
	'db=s'       => \$database, # database name
	'feature=s'  => \$feature, # name of genomic feature to analyze
	'method=s'   => \$method, # method of collecting & reporting data
	'value=s'    => \$value_type, # type of data to collect
	'dataset=s'  => \@datasets, # the list of datasets to collect data from
	'dataf=s'    => \@datafiles, # list of data files to collect data from
	'log!'       => \$log, # dataset is in log2 space
	'strand=s'   => \$strand, # indicate strandedness of data
	'exons!'     => \$subfeature, # indicate to restrict to subfeatures
	'extend=i'   => \$extend, # extend the size of the genomic feature
	'start=s'    => \$start, # adjustment to relative start position
	'stop=s'     => \$stop, # adjustment relative stop position
	'fstart=f'   => \$fstart, # fractional start position
	'fstop=f'    => \$fstop, # fractional stop position
	'limit=i'    => \$limit, # size limit to fractionate a feature
	'pos=i'      => \$position, # set the relative feature position
	'win=i'      => \$win, # indicate the size of genomic intervals
	'step=i'     => \$step, # step size for genomic intervals
	'set_strand' => \$set_strand, # enforce a specific strand
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
	
	# set convenience method, for backwards compatibility
	if ($method eq 'count' or $method eq 'enumerate') {
		$method = 'sum';
		$value_type = 'count';
	}
}
else {
	# set the default to use the mean
	print " Using default method of 'mean'\n";
	$method = 'mean';
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

if (defined $strand) { 
	# check the strand request that was defined on the command line
	unless (
		$strand eq 'sense' or 
		$strand eq 'antisense' or 
		$strand eq 'all'
	) {
		die " unknown strand '$strand'!";
	}
} 
else {
	# default value
	$strand = 'all'; 
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
		if ( $main_data_ref->{'feature'} ) {
			# load the feature defined in the file's metadata
			$feature = $main_data_ref->{'feature'};
		}
		else {
			# the file's metadata doesn't have the feature defined
			# most likely because the file doesn't have metadata
			# let's try looking for known column names
			my $name_i = find_column_index($main_data_ref, "name");
			my $type_i = find_column_index($main_data_ref, "type");
			my $chr_i = find_column_index($main_data_ref, "^chr|seq");
			my $start_i = find_column_index($main_data_ref, "^start");
			my $stop_i = find_column_index($main_data_ref, "^stop|end");
			
			# guess 
			if (defined $name_i and defined $type_i) {
				$feature = 'Named feature';
			}
			elsif (defined $chr_i and defined $start_i and defined $stop_i) {
				$feature = 'genome';
			}
			else {
				$feature = 'unknown feature';
			}
		}
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

# check for comma-delimited datasets
if (scalar @datasets == 1 and $datasets[0] =~ /,/) {
	@datasets = split /,/, shift @datasets;
}
if (scalar @datafiles == 1 and $datafiles[0] =~ /,/) {
	@datafiles = split /,/, shift @datafiles;
}

# collect the data files if provided
if (@datafiles) {
	foreach my $file (@datafiles) {
		# look for transfer protocol
		if ($file =~ /^(?:ftp|http):/) {
			# BigWig, BigBed, and Bam files are supposed to be readable
			# across network transfer protocols
			# this is handled by the appropriate modules
			# I have not tested this feature
			push @datasets, $file;
		}
		else {
			# otherwise a local file
			# check first for combined files with &
			if ($file =~ /&/) {
 				# need to prefix each file name
 				push @datasets, join("&", map {"file:$_"} split(/&/, $file) );
 			}
			else {
				push @datasets, "file:$file";
			}
		}
	}
}


# Prcoess the provided datasets
if (@datasets) {
	# Datasets were requested on the command line, automatic execution
	unless (scalar(@datasets) == 1 and $datasets[0] eq 'none') {
		# we can skip all data collection if the only dataset is 'none'
		auto_validate_and_collect_datasets();
	}
} 

else {
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
		} );
	}
}



# Automatic dataset collection
sub auto_validate_and_collect_datasets {
	# first validate the proffered datasets, and then process
	
	foreach my $dataset (@datasets) {
		# do this one at a time
		# we need to determine whether the dataset is database feature
		# or a data file, either local or remote
		if ($dataset =~ /^(?:http|ftp|file):/) {
			# data file, either local or remote
			# no database validation necessary
			# process the request
			submit_dataset_request($dataset);
		}
		else {
			# must be a database feature type
			# validate
			my $bad = validate_dataset_list($db, $dataset);
			if ($bad) {
				warn " The dataset '$dataset' is not valid. Skipping.\n";
				next;
			}
			else {
				# the dataset is valid
				submit_dataset_request($dataset);
			}
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
	@datasets = parse_list($answer); # split the user requests
	
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
		my $column_name = get_genome_dataset( {
				'data'     => $main_data_ref->{'data_table'},
				'db'       => $db,
				'dataset'  => $dataset,
				'method'   => $method,
				'strand'   => $strand,
				'log'      => $logstatus,
				'value'    => $value_type,
		} );
		
		if (defined $column_name) {
			print " in ";
			printf "%.1f", ( (time - $start_time)/60 ); 
			print " minutes\n";
			
			# record metadata
			_record_metadata($dataset, $column_name, $logstatus);
		}
		else {
			# error messages should have printed from the subroutine
			print " failed!\n";
		}
		
		
	}
	
	else {
		# collecting a dataset for list of features
		my $column_name = get_feature_dataset( {
				'data'        => $main_data_ref->{'data_table'},
				'db'          => $db,
				'dataset'     => $dataset,
				'method'      => $method,
				'strand'      => $strand,
				'log'         => $logstatus,
				'extend'      => $extend,
				'start'       => $start,
				'stop'        => $stop,
				'position'    => $position,
				'fstart'      => $fstart,
				'fstop'       => $stop,
				'limit'       => $limit,
				'set_strand'  => $set_strand,
				'value'       => $value_type,
				'subfeature'  => $subfeature,
		} );
		if (defined $column_name) {
			print " in ";
			printf "%.1f", ( (time - $start_time)/60 ); 
			print " minutes\n";
			
			# record metadata
			_record_metadata($dataset, $column_name, $logstatus);
		}
		else {
			# error messages should have printed from the subroutine
			print " failed!\n";
		}
	}
	
}


# subroutine to record the metadata for a dataset
sub _record_metadata {
	my ($dataset, $column_name, $logstatus) = @_;
	
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
		'strand'   => $strand, # global argument for strand specificity
		'value'    => $value_type, # global argument for value type
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
	if ($subfeature) {
		$metadata{'exons'} = 'yes';
	}
	# add database name if different
	if ($database ne $main_data_ref->{'db'}) {
		$metadata{'db'} = $database;
	}
	if ($set_strand) {
		$metadata{'strand_implied'} = 'yes';
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
		'sum'        => 1,
		'count'      => 1,
		'enumerate'  => 1,
	);
	
	# return true if this method exists
	return exists $acceptable{$method};
}







__END__

=head1 NAME

get_datasets.pl

A script to collect genomic datasets from a Bioperl SeqFeature::Store db.

=head1 SYNOPSIS

get_datasets.pl [--options...] [<filename>]
  
  Options:
  --new
  --in <filename>
  --out filename
  --db <name|file.gff3>
  --feature <type | type:source | alias, ...>
  --dataset <none | name, ...>
  --dataf <file1,file2,...>
  --method [mean|median|stddev|min|max|range|sum]
  --value [score|count|length]
  --(no)log
  --strand [all|sense|antisense]
  --exons
  --extend <integer>
  --start <integer>
  --stop <integer>
  --fstart <decimal>
  --fstop <decimal>
  --limit <integer>
  --pos [5|3|m]
  --win <integer>
  --step <integer>
  --set_strand
  --(no)gz
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --new

Generate a new table of features. Overrides any specified input file.

=item --in <filename>

Specify the file name of a previously generated feature dataset.
Required unless "--new" is specified. It should be in the tim data
format that is generated by this program and others, although other
tab-delimited text data formats may be usable. See the file
description below and in C<tim_db_helper.pm>.

=item --out <filename>

Specify the output file name. Required for new feature tables; optional for 
current files. If this is argument is not specified then the input file is 
overwritten.

=item --db <name|file.gff3>

Specify the name of the BioPerl SeqFeature::Store database to use as
source. Alternatively, a single GFF3 file may be loaded into a in-memory
database. Specifying the database is required for new feature data files.
For pre-existing input data files, this value may be obtained from the
input file metadata. However, if provided, it overrides the database listed
in the file; this is useful for collecting data from multiple databases.

=item --feature <type | type:source | alias, ...>

Specify the type of feature from which to collect values. This is required 
for new feature tables. Two types of values may be passed: either a specific 
feature type present in the database, or an alias to one or more features. 
The feature may be specified as either type or type:source. Aliases are 
specified in the C<tim_db_helper.cfg> file, and provide a shortcut to a 
list of one or more features. More than feature may be included as a 
comma-delimited list (no spaces).

=item --dataset <none | name, ...>

Provide the name of the dataset to collect the values. Use this argument 
repeatedly for each dataset to be collected. Two or more datasets may be
merged into one by delimiting with an "&" (no spaces!). If the dataset is not 
specified on the command line, then the program will interactively present a 
list of datasets from the database to select. To force the program to 
simply write out the list of collected features without collecting data, 
provide the dataset name of "none".

=item --dataf <file1,file2,...>

Alternative to using datasets stored in the database, a data file may be
directly referenced and used in the collection of score values. Only data
formats which directly support chromosome:start..stop searches are supported;
these currently include BigWig (.bw), BigBed (.bb), and single-end Bam
(.bam) files. Multiple data files may be specified as a comma delimited
list. Merging files as a single source is also supported using a "&".
Remote files should also be accessible by prefixing with the appropriate
transfer protocol (http: or ftp:).

=item --method [sum|mean|median|stddev|min|max|range]

Specify the method for combining all of the dataset values within the 
genomic region of the feature. Accepted values include:
  
  - sum
  - mean        (default)
  - median
  - stddev      Standard deviation of the population (within the region)
  - min
  - max
  - range       Returns 'min-max'
  
=item --value [score|count|length]

Optionally specify the type of data value to collect from the dataset or 
data file. Three values are accepted: score, count, or length. The default 
value type is score. Note that some data sources only support certain 
types of data values. Wig and BigWig files only support score and count; 
BigBed and database features support count and length and optionally 
score; Bam files support basepair coverage (score), count (number of 
alignments), and length.

=item --(no)log

Indicate the dataset is (not) in log2 space. The log2 status of the dataset is 
critical for accurately mathematically combining the dataset values in the 
feature's genomic region. It may be determined automatically if the dataset 
name includes the phrase "log2".

=item --strand [all | sense | antisense]

Specify whether stranded data should be collected for each of the 
datasets. Either sense or antisense (relative to the feature) data 
may be collected. Note that strand is not (currently) supported with 
coverage (score) from a BAM file, but count is. The default value is 
'all', indicating all data will be collected.  

=item --exons

Optionally indicate that data should be collected only over the exon 
subfeatures of a gene or transcript, rather than the entire gene. 
Subfeatures with a primary_tag of exon are preferentially taken. If exons 
are not defined, then CDS and UTR subfeatures are used, or the entire 
gene or transcript if no appropriate subfeatures are found. Note that the 
data collection method is applied twice, once for each subfeature, and then 
again on all of the subfeature combined values. Also note that the options 
extend, start, stop, fstart, and fstop are ignored. Default is false.

=item --extend <integer>

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
size. The default size is defined in the configuration file, 
tim_db_helper.cfg.

=item --step <integer>

Optionally indicate the step size when generating a new list of intervals 
across the genome. The default is equal to the window size.

=item --set_strand

For features that are not inherently stranded (strand value of 0), 
impose an artificial strand for each feature (1 or -1). This will 
have the effect of enforcing a relative orientation for each feature, 
or to collected stranded data. This requires the presence of a 
column in the input data file with a name of "strand". Hence, it 
will not work with newly generated datasets, but only with input 
data files. Default is false.

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


