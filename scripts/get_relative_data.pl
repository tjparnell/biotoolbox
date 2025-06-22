#!/usr/bin/perl

# documentation at end of file

use warnings;
use strict;
use English      qw(-no_match_vars);
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	calculate_score
	$BAM_ADAPTER
	$BIG_ADAPTER
);
use Bio::ToolBox::utility qw(format_with_commas simplify_dataset_name);

my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};

# This constant determines the maximum size of the dataset hash to be
# returned from the get_region_dataset_hash(). To increase performance,
# the program normally queries the database once for each feature or
# region, and a hash returned with potentially a score for each basepair.
# This may become unwieldy for very large regions, which may be better
# served by separate database queries for each window.
use constant DATASET_HASH_LIMIT => 4999;

our $VERSION = '2.02';

print
	"\n A script to collect windowed data flanking a relative position of a feature\n\n";

### Quick help
unless (@ARGV) {    # when no command line options are present
					# when no command line options are present
					# print SYNOPSIS
	pod2usage(
		{
			'-verbose' => 0,
			'-exitval' => 1,
		}
	);
}

### Get command line options and initialize values

## Initialize values
my (
	$infile,      $parse,     $outfile,      $main_database, $data_database,
	$feature,     $method,    $win,          $number,        $up_number,
	$down_number, $position,  $strand_sense, $set_strand,    $avoid,
	$avoidtype,   $long_data, $smooth,       $sum,           $format,
	$groupcol,    $gz,        $cpu,          $help,          $print_version,
);    # command line variables
my @datasets;

## Command line options
GetOptions(
	'o|out=s'      => \$outfile,          # output file name
	'i|in=s'       => \$infile,           # input file name
	'parse!'       => \$parse,            # parse input file
	'd|db=s'       => \$main_database,    # main or annotation database name
	'D|ddb=s'      => \$data_database,    # data database
	'a|data=s'     => \@datasets,         # dataset name
	'f|feature=s'  => \$feature,          # type of feature
	'm|method=s'   => \$method,           # method to combine data
	'w|window=i'   => \$win,              # window size
	'n|number=i'   => \$number,           # number of windows
	'up=i'         => \$up_number,        # number of windows upstream
	'down=i'       => \$down_number,      # number of windows downstream
	'p|position=s' => \$position,         # indicate relative location of the feature
	't|strand=s'   => \$strand_sense,     # collected stranded data
	'force_strand|set_strand' => \$set_strand,    # enforce an artificial strand
		# force_strand is preferred option, but respect the old option
	'avoid!'     => \$avoid,            # avoid conflicting features
	'avtype=s'   => \$avoidtype,        # type of feature to avoid
	'long!'      => \$long_data,        # collecting long data features
	'smooth!'    => \$smooth,           # smooth by interpolation
	'U|sum!'     => \$sum,              # generate average profile
	'r|format=i' => \$format,           # decimal formatting
	'g|groups'   => \$groupcol,         # write group column file
	'z|gz!'      => \$gz,               # compress the output file
	'c|cpu=i'    => \$cpu,              # number of execution threads
	'h|help'     => \$help,             # print help
	'v|version'  => \$print_version,    # print the version
	'bam=s'      => \$BAM_ADAPTER,      # explicitly set the bam adapter
	'big=s'      => \$BIG_ADAPTER,      # explicitly set the big adapter
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {

	# print entire POD
	pod2usage(
		{
			'-verbose' => 2,
			'-exitval' => 1,
		}
	);
}

# Print version
if ($print_version) {
	print " Biotoolbox script get_relative_data.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

## Check for required values
my $formatter;
check_defaults();
my $start_time = time;

## Generate or load the input dataset
my $Data;
if ($infile) {
	$Data = Bio::ToolBox::Data->new(
		file    => $infile,
		parse   => $parse,
		feature => $feature,
	) or die " unable to load input file '$infile'\n";
	if ( $Data->last_row ) {
		printf " Loaded %s features from $infile.\n",
			format_with_commas( $Data->last_row );
	}
	else {
		print STDERR " FATAL: No features loaded from file '$infile'!\n";
		exit 1;
	}

	# update main database as necessary
	if ($main_database) {
		if ( $main_database ne $Data->database ) {

			# update with new database
			printf " updating main database name from '%s' to '%s'\n",
				$Data->database, $main_database;
			$Data->database($main_database);
		}
	}
	elsif ( $Data->database ) {
		$main_database = $Data->database if $Data->database !~ /^Parsed/;
	}

	# update feature type as necessary
	if (    not defined $Data->feature
		and not defined $Data->type_column
		and defined $feature )
	{
		$Data->feature($feature);
	}

	# set headers to true if file was explicitly not parsed
	# since warnings will occur when we add columns and write to a file
	if ($parse == 0) {
		$Data->headers(1);
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
$Data->program("$PROGRAM_NAME, v $VERSION");

# check peak summit
if ( $position == 9 and $Data->format ne 'narrowPeak' ) {
	print
		" Peak summit indicated as reference point, but input file is not narrowPeak!\n";
	print "  using interval midpoint instead\n";
	$position = 4;
}

# the number of columns already in the data array
my $beginningcolumn = $Data->last_column;
my $startcolumn;    # this is now calculated separately for each dataset

# Check output file name
unless ($outfile) {
	if ( $Data->filename ) {

		# reuse the input filename and overwrite
		$outfile = $Data->filename;
	}
	elsif ( $Data->basename ) {

		# no filename but a basename is indicative of a parsed file
		# so stick a txt extension on it, and write in current directory
		$outfile = sprintf "%s.txt", $Data->basename;
	}
	else {
		print STDERR " FATAL: No output file provided!\n";
		exit 1;
	}
}

# make sure data table supports avoid option
if ($avoid) {
	unless ($main_database) {
		print
" WARNING: avoid option not supported without an annotation database! Disabling\n";
		undef $avoid;    # 0 complicates things, leave undefined
	}
	if ($avoidtype) {
		if ( $avoidtype =~ /,/ ) {

			# a comma delimited list
			$avoid = [ split( /,/, $avoidtype ) ];
		}
		else {
			# presume a single feature type, take it as is
			$avoid = [$avoidtype];
		}
	}
}

## Prepare to collect data
# Open data database
my $ddb;
if ( defined $data_database ) {

	# specifically defined a data database
	$ddb = open_db_connection($data_database)
		or die "unable to establish data database connection to $data_database!\n";
}

# Check the dataset
@datasets = verify_or_request_feature_types(
	'db'      => $ddb || $Data->database,
	'feature' => \@datasets,
	'prompt'  => " Enter the dataset(s) or feature type(s) from which \n"
		. " to collect data. Comma delimited or range is acceptable\n",
);
unless (@datasets) {
	print STDERR
" FATAL: No verifiable dataset(s) provided. Check your file path, database, or dataset.\n";
	exit 1;
}

## Collect the relative data

# Determine starting and ending points
my $starting_point = 0 - ( $win * $up_number );

# default values will give startingpoint of -1000
my $ending_point = $win * $down_number;

# likewise default values will give endingpoint of 1000

# Print collection statement
printf
	" Collecting %s $method data between %d..%d from the %s in %d bp windows from %s\n",
	$long_data ? 'long' : 'hashed',
	$starting_point, $ending_point,
	$position == 5    ? "5' end"
	: $position == 4  ? "midpoint"
	: $position == 10 ? "peak summit"
	: "3' end",
	$win, join( ', ', @datasets ),;

# check whether it is worth doing parallel execution
if ( $cpu > 1 and ( $Data->number_rows / $cpu ) < 100 ) {
	while ( $cpu > 1 and ( $Data->number_rows / $cpu ) < 100 ) {

		# I figure we need at least 100 lines in each fork split to make
		# it worthwhile to do the split, otherwise, reduce the number of
		# splits to something more worthwhile
		$cpu--;
	}
	print " Reducing number of process forks to $cpu due to number of input features\n";
}

if ( $cpu > 1 ) {

	# parallel execution
	print " Forking into $cpu children for parallel data collection...\n";
	parallel_execution();
}
else {
	# single process execution
	single_execution();
}

## Generate summed data
# an average across all features at each position suitable for plotting
if ($sum) {
	print " Generating trimmed-mean summary file....\n";
	my $sumfile = $Data->summary_file(
		'filename' => $outfile,
		'dataset'  => \@datasets,
		'method'   => 'trimmean'
	);
	if ($sumfile) {
		print " Wrote summary file '$sumfile'\n";
	}
	else {
		print " Unable to write summary file!\n";
	}
}

## Output the data
my $written_file = $Data->save(
	'filename' => $outfile,
	'gz'       => $gz,
);
if ($written_file) {

	# success!
	print " wrote file $written_file\n";
}
else {
	# failure! the subroutine will have printed error messages
	print " unable to write output file!\n";
}

# write the column group file
if ($groupcol) {
	my $groupfile = $written_file;
	$groupfile =~ s/\.txt (?:\.gz)? $/.col_groups.txt/x;
	my $fh = Bio::ToolBox::Data->open_to_write_fh($groupfile);
	$fh->print("Name\tDataset\n");
	for ( my $i = $beginningcolumn + 1; $i <= $Data->last_column; $i++ ) {
		my $name    = $Data->name($i);
		my $dataset = $name;
		$dataset =~ s/:[\-\d]+$//;
		$fh->print("$name\t$dataset\n");
	}
	$fh->close;
	print " wrote group column file $groupfile\n";
}

printf " Completed in %.1f minutes\n", ( time - $start_time ) / 60;

#### Subroutines #######

## check required variables and assign default values
sub check_defaults {
	unless ( $main_database or $infile ) {
		print STDERR
" FATAL: You must define a database or input file!\n Use --help for more information\n";
		exit 1;
	}
	$parse = 1 if ( $infile and not defined $parse );

	unless ( $outfile or $infile ) {
		print STDERR
" FATAL: You must define an output filename !\n Use --help for more information\n";
		exit 1;
	}

	# check datasets
	if ( not @datasets and @ARGV ) {
		@datasets = @ARGV;
	}
	if ( @datasets and $datasets[0] =~ /,/ ) {

		# seems to be a comma delimited list, possibly more than one?????
		my @list;
		foreach my $d (@datasets) {
			push @list, ( split /,/, $d );
		}
		@datasets = @list;
	}

	# window size
	unless ($win) {
		print " Using default window size of 50 bp\n";
		$win = 50;
	}

	# number of windows
	if ( $up_number or $down_number ) {

		# either one could be set, assume the other is zero
		$up_number   ||= 0;
		$down_number ||= 0;
	}
	elsif ($number) {
		$up_number   = $number;
		$down_number = $number;
	}
	else {
		print " Using default window number of 20 per side\n";
		$up_number   = 20;
		$down_number = 20;
	}

	if ( defined $method ) {

		# check the requested method
		unless ( $method eq 'mean'
			or $method eq 'median'
			or $method eq 'sum'
			or $method eq 'min'
			or $method eq 'max'
			or $method eq 'stddev'
			or $method =~ /^\w?count$/ )
		{
			print STDERR
				" FATAL: Unknown method '$method'!\n Use --help for more information\n";
			exit 1;
		}
	}
	else {
		$method = 'mean';
	}

	if ( defined $position ) {

		# check the position value
		if ( $position !~ /^ (?: 5 | 3 | 53 | m | p | 4 ) $/x ) {
			print STDERR " FATAL: Unknown relative position '$position'!\n";
			exit 1;
		}

		# change to match internal usage
		if ( $position eq 'm' ) {
			$position = 4;    # midpoint
		}
		elsif ( $position eq 'p' ) {
			$position = 10;    # narrowPeak summit
		}
	}
	else {
		# default position to use the 5' end
		$position = 5;
	}

	if ( defined $strand_sense ) {
		unless ( $strand_sense eq 'sense'
			or $strand_sense eq 'antisense'
			or $strand_sense eq 'all' )
		{
			print STDERR " FATAL: Unknown strand value '$strand_sense'!\n";
			exit 1;
		}
	}
	else {
		# default
		$strand_sense = 'all';
	}

	unless ( defined $sum ) {

		# assume to write a summary file, nearly always want this, at least I do
		$sum = 1;
	}

	# generate formatter
	if ( defined $format ) {
		$formatter = '%.' . $format . 'f';
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

## Run in parallel
sub parallel_execution {
	my $pm = Parallel::ForkManager->new($cpu);

	# generate base name for child processes
	my $child_base_name = $outfile . ".$PID";

	# Split the input data into parts and execute in parallel in separate forks
	for my $i ( 1 .. $cpu ) {
		$pm->start and next;

		#### In child ####

		# splice the data structure
		$Data->splice_data( $i, $cpu );

		# re-open database objects to make them clone safe
		# pass second true to avoid cached database objects
		my $db = $Data->open_database(1);
		if ($data_database) {
			$ddb = open_db_connection( $data_database, 1 );
		}

		# work through each dataset
		foreach my $dataset (@datasets) {

			# new start column for this dataset
			$startcolumn = $Data->number_columns + 1;

			# Add the columns for each window
			# and calculate the relative starting and ending points
			prepare_window_datasets($dataset);

			# determine long data collection for very large regions
			if ( $ending_point - $starting_point > DATASET_HASH_LIMIT ) {
				$long_data = 1;
			}

			# Select the appropriate method for data collection
			if ($long_data) {
				map_relative_long_data($dataset);
			}
			else {
				map_relative_data($dataset);
			}

			# Interpolate values
			if ($smooth) {
				go_interpolate_values();
			}
		}

		# write out result
		my $success = $Data->save(
			'filename' => sprintf( "$child_base_name.%03s", $i ),
			'gz'       => 0,    # faster to write without compression
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
	unless ( scalar @files == $cpu ) {
		die "only found " . scalar(@files) . " child files when there should be $cpu!\n";
	}
	my $count = $Data->reload_children(@files);
	printf " reloaded %s features from children\n", format_with_commas($count);
}

## Run in single thread
sub single_execution {

	# work through each dataset
	foreach my $dataset (@datasets) {

		# new start column for this dataset
		$startcolumn = $Data->number_columns + 1;

		# Add the columns for each window
		# and calculate the relative starting and ending points
		prepare_window_datasets($dataset);

		# determine long data collection for very large regions
		if ( $ending_point - $starting_point > DATASET_HASH_LIMIT ) {
			$long_data = 1;
		}

		# Select the appropriate method for data collection
		if ($long_data) {
			map_relative_long_data($dataset);
		}
		else {
			map_relative_data($dataset);
		}

		# Interpolate values
		if ($smooth) {
			print " Interpolating missing values....\n";
			go_interpolate_values();
		}
	}
}

## Prepare columns for each window
sub prepare_window_datasets {
	my $dataset = shift;

	# generate a simplified new name
	my $new_name = simplify_dataset_name($dataset);

	# Prepare and annotate the header names and metadata
	for ( my $start = $starting_point; $start < $ending_point; $start += $win ) {

		# we will be progressing from the starting to ending point
		# in increments of window size

		# set the stop position
		my $stop = $start + $win - 1;

		# deal with the pesky 0 value
		# since we're working with 1-base coordinates, we don't really have
		# a real 0 coordinate, so need to skip it as it doesn't really exist
		my $zero_check;
		for my $i ( $start .. $stop ) {
			$zero_check = 1 if $i == 0;
		}
		if ($zero_check) {

			# adjust the coordinates accordingly
			if ( $start == 0 ) {
				$start += 1;
				$stop  += 1;
			}
			elsif ( $stop == 0 ) {
				$stop += 1;
			}
			else {
				# some number in between
				$stop += 1;
			}
		}

		# add new column
		my $new_index = $Data->add_column( sprintf( "%s:%d", $new_name, $start ) );

		# set the metadata
		$Data->metadata( $new_index, 'start',          $start );
		$Data->metadata( $new_index, 'stop',           $stop );
		$Data->metadata( $new_index, 'window',         $win );
		$Data->metadata( $new_index, 'dataset',        $dataset );
		$Data->metadata( $new_index, 'method',         $method );
		$Data->metadata( $new_index, 'decimal_format', $format ) if defined $format;
		if ( $position == 5 ) {
			$Data->metadata( $new_index, 'relative_position', '5prime_end' );
		}
		elsif ( $position == 4 ) {    # midpoint
			$Data->metadata( $new_index, 'relative_position', 'center' );
		}
		elsif ( $position == 9 ) {
			$Data->metadata( $new_index, 'relative_position', 'peak_summit' );
		}
		else {
			$Data->metadata( $new_index, 'relative_position', '3prime_end' );
		}
		if ($set_strand) {
			$Data->metadata( $new_index, 'strand_implied', 1 );
		}
		if ( $strand_sense =~ /sense/ ) {
			$Data->metadata( $new_index, 'strand', $strand_sense );
		}
		if ($data_database) {
			$Data->metadata( $new_index, 'db', $data_database );
		}
		if ($avoid) {
			$Data->metadata( $new_index, 'avoid',
				ref($avoid) eq 'ARRAY' ? join( ',', @{$avoid} ) : $avoid );
		}
	}
}

## Collect relative scores using single database call
sub map_relative_data {
	my $dataset = shift;

	# Collect the data
	my $stream = $Data->row_stream;
	while ( my $row = $stream->next_row ) {
		my $regionscores = $row->get_relative_point_position_scores(
			'ddb'      => $ddb,
			'dataset'  => $dataset,
			'position' => $position,
			'extend'   => $ending_point,                     # equivalent to the extension
			'stranded' => $strand_sense,
			'strand'   => $set_strand ? $row->strand : undef,
			'avoid'    => $avoid,
			'method'   => $method    # only really matters for count
		);

		# assign the scores to the windows in the region
		for my $column ( $startcolumn .. $Data->number_columns ) {

			# record nulls if no data returned
			unless ( scalar keys %{$regionscores} ) {
				$row->value( $column, '.' );
				next;
			}

			# get start and stop
			my $start = $Data->metadata( $column, 'start' );
			my $stop  = $Data->metadata( $column, 'stop' );

			# collect a score at each position in the window
			my @scores = map { $regionscores->{$_} }
				grep { $_ >= $start and $_ <= $stop }
				keys %{$regionscores};

			# put the value into the data table
			my $score = calculate_score( $method, \@scores );
			if ( $formatter and looks_like_number($score) ) {
				$score = sprintf( $formatter, $score );
			}
			$row->value( $column, $score );
		}
	}
}

## Collect relative data using individual window database calls
sub map_relative_long_data {
	my $dataset = shift;

	# Collect the data
	my $stream = $Data->row_stream;
	while ( my $row = $stream->next_row ) {

		# calculate a reference point using an internal feature method
		my $reference = $row->calculate_reference($position);

		# collect the data for every window
		for my $column ( $startcolumn .. $Data->number_columns ) {

			# we must modify the start and stop position with the adjustments
			# recorded in the current column metadata
			my $score = $row->get_score(
				'db'      => $ddb,
				'dataset' => $dataset,
				'chromo'  => $row->seq_id,
				'start'   => $row->strand >= 0
				? $reference + $Data->metadata( $column, 'start' )
				: $reference - $Data->metadata( $column, 'start' ),
				'stop' => $row->strand >= 0
				? $reference + $Data->metadata( $column, 'stop' )
				: $reference - $Data->metadata( $column, 'stop' ),
				'strand'   => $row->strand,
				'method'   => $method,
				'stranded' => $strand_sense,
			);
			if ( $formatter and looks_like_number($score) ) {
				$score = sprintf( $formatter, $score );
			}
			$row->value( $column, $score );
		}
	}
}

## Interpolate the '.' values with the mean of the neighbors
sub go_interpolate_values {

	# determine counts
	my $lastwindow = $Data->last_column;

	# lastwindow is the index of the last column
	# walk through each data line and then each window
	my $stream = $Data->row_stream;
	while ( my $row = $stream->next_row ) {
		my $col = $startcolumn + 1;
		while ( $col < $lastwindow ) {

			# walk through the windows of a data row
			# skipping the very first and last windows (columns)
			# we will look for null values
			# if one is found, interpolate from neighbors
			if ( $row->value($col) eq '.' and $row->value( $col - 1 ) ne '.' ) {

				# we will interpolate the value from the neighbors
				# first, need to check that the neighbors have values

				# find the next real value
				my $next_i;
				for my $i ( ( $col + 1 ) .. $lastwindow ) {
					if ( $row->value($i) ne '.' ) {
						$next_i = $i;
						last;
					}
				}
				unless ( defined $next_i ) {
					$col++;
					next;
				}

				# determine fractional value
				my $initial = $row->value( $col - 1 );
				my $fraction =
					( $row->value($next_i) - $initial ) / ( $next_i - $col + 1 );

				# apply fractional values
				for my $n ( $col .. ( $next_i - 1 ) ) {
					my $score = $initial + ( $fraction * ( $n - $col + 1 ) );
					$score = sprintf( $formatter, $score )
						if ( $formatter and $score ne '.' );
					$row->value( $n, $score );
				}

				# jump ahead
				$col = $next_i;
			}
			$col++;
		}
	}
}

__END__

=head1 NAME

get_relative_data.pl

A program to collect data in bins around a relative position.

=head1 SYNOPSIS
 
get_relative_data.pl [--options] --in <filename> --out <filename>
  
get_relative_data.pl [--options] -i <filename> <data1> <data2...>
  
  Options for data files:
  -i --in <filename>                  input file: txt bed gff gtf refFlat ucsc
  -o --out <filename>                 optional output file, default overwrite 
  
  Options for new files
  -d --db <name>                      annotation database: mysql sqlite
  -f --feature <type>                 one or more feature types from db or gff
  
  Options for data collection:
  -D --ddb <name|file>                data or BigWigSet database
  -a --data <dataset|filename>        data from which to collect: bw bam etc
  -m --method [mean|median|stddev|    statistical method for collecting data
            min|max|range|sum|count|   default mean
            pcount|ncount]
  -t --strand [all|sense|antisense]   strand of data relative to feature (all)
  --force_strand                      use the specified strand in input file
  --avoid                             avoid neighboring features
  --avtype [type,type,...]            alternative types of feature to avoid
  --long                              collect each window independently
  -r --format <integer>               number of decimal places for numbers
  
  Bin specification:
  -w --win <integer>                  size of windows, default 50 bp
  -n --num <integer>                  number of windows flanking reference, 20
  --up <integer>                        or number of windows upstream
  --down <integer>                      and number of windows downstream
  -p --pos [5|m|3|p]                  reference position, default 5'
  
  Post-processing:
  -U --sum                            generate summary file
  --smooth                            smoothen sparse data
  
  General Options:
  -g --groups                         write columns group index file for plotting
  -z --gz                             compress output file
  -c --cpu <integer>                  number of threads, default 4
  --noparse                           do not parse input file into SeqFeatures
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

Specify the output file name. Required for new files; otherwise, 
input files will be overwritten unless specified.

=back

=head2 Options for new files

=over 4

=item --db E<lt>name | filenameE<gt>

Specify the name of a L<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A database is 
required for generating new data files with features. This option may 
skipped when using coordinate information from an input file (e.g. BED 
file), or when using an existing input file with the database indicated 
in the metadata.  

=item --feature [type, type:source]

Specify the type of feature to map data around. The feature may be 
listed either as GFF type or GFF type:source. The list 
of features will be automatically generated from the database. 
This is only required when an input file is not specified. 

=back

=head2 Options for data collection

=over 4

=item --ddb E<lt>name | filenameE<gt>

If the data to be collected is from a second database that is separate 
from the annotation database, provide the name of the data database here. 
Typically, a second L<Bio::DB::SeqFeature::Store> or BigWigSet database 
is provided here. 

=item --data E<lt>dataset_name | filenameE<gt>

Provide the name of the dataset to collect the values. If no 
dataset is specified on the command line, then the program will 
interactively present a list of datasets from the data database to select. 

The dataset may be a database file, including bigWig (.bw), 
bigBed (.bb), or Bam alignment (.bam) files. The files may be local or 
remote (specified with a http: or ftp: prefix).

Alternatively, the dataset may be a feature type in a BioPerl 
L<Bio::DB::SeqFeature::Store> or L<Bio::DB::BigWigSet> database. Provide 
either the feature type or C<type:source>. 

More than one datasource may be provided; use multiple data options or list 
the datasets at the end of the command.

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

=item --strand [sense|antisense|all]

Specify whether stranded data should be collected for each of the 
datasets. Either sense or antisense (relative to the feature) data 
may be collected. The default value is 'all', indicating all 
data will be collected.

=item --force_strand

For features that are not inherently stranded (strand value of 0)
or that you want to impose a different strand, set this option when
collecting stranded data. This will reassign the specified strand for
each feature regardless of its original orientation. This requires the
presence of a "strand" column in the input data file. This option only
works with input file lists of database features, not defined genomic
regions (e.g. BED files). Default is false.

=item --avoid

Indicate whether neighboring features should be avoided when calculating 
values in a window. After collecting the data, each window is checked for 
overlapping features; if the window overlaps another feature, no value 
is reported for that window. This option requires using an annotation 
database (--db option). This is useful to avoid scoring windows that 
overlap a neighboring gene, for example. The default is false (return 
all values regardless of overlap).

=item --avtype [type,type,...]

Provide a feature type (primary_tag or primary_tag:source) or a 
comma-delimited list of types to be used when avoiding neighboring 
features. The default is to avoid features of the same type as that 
of the query, i.e. collecting data around a feature of type 'gene' 
will avoid other 'gene' features. This option allows you to avoid 
other features too, such as 'tRNA' or 'repeat'.

=item --long

Indicate that data should be collected independently for each long 
window. This may be enabled automatically if the sum of the entire 
window length passes a predefined threshold. The default for 'short' 
windows is to collect all of the point data from the dataset first, and 
then divide the results into the different windows. Datasets consisting 
of "long" features, for example long alignments, may be counted more 
than once in long mode when they span multiple windows.

=item --format E<lt>integerE<gt>

Specify the number of decimal positions to format the collected scores. 
Default is not to format, often leading to more than the intended 
significant digits.

=back

=head2 Bin specification

=over 4

=item --win E<lt>integerE<gt>

Specify the window size. The default is 50 bp.

=item --num E<lt>integerE<gt>

Specify the number of windows on either side of the feature position 
(total number will be 2 x [num]). The default is 20, or 1 kb on either 
side of the reference position if the default window size is used.

=item --up E<lt>integerE<gt>

=item --down E<lt>integerE<gt>

Alternatively specify the exact number of windows upstream and 
downstream of the reference position. If only one option is set, 
then the other option is assumed to be zero. 

=item --pos [5|m|3]

Indicate the relative position of the feature around which the 
data is mapped. Three values are accepted: "5" indicates the 
5' prime end is used, "3" indicates the 3' end is used, and "m" 
indicates the middle of the feature is used. The default is to 
use the 5' end, or the start position of unstranded features. 

=back

=head2 Post-processing

=over 4

=item --(no)sum

Indicate that the data should be averaged across all features at
each position, suitable for graphing. A separate text file will 
be written with the suffix '_summed' with the averaged data. 
Default is true (sum).

=item --smooth

Indicate that windows without values should (not) be interpolated
from neighboring values. The default is false (nosmooth).

=back

=head2 General options

=over 4

=item --groups

Optionally write a secondary file with the list of column group names and 
their corresponding dataset group. This can be used to assist in designating 
the metadata when plotting files, for example in R with pheatmap. The 
file is named the output basename appended with F<.col_groups.txt>.

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --cpu E<lt>integerE<gt>

Specify the number of CPU cores to execute in parallel. This requires 
the installation of Parallel::ForkManager. With support enabled, the 
default is 2. Disable multi-threaded execution by setting to 1. 

=item --noparse

Prevent input annotation files from being automatically parsed into sequence 
features. Coordinates will be used as is and new data columns will be appended 
to the input file. 

=item --version

Print the version number.

=item --help

Display this help

=back

=head1 DESCRIPTION

This program will collect data around a relative coordinate of a genomic 
feature or region. The data is collected in a series of windows flanking the 
feature start (5' position for stranded features), end (3' position), or 
the midpoint position. The number and size of windows are specified via 
command line arguments, or the program will default to 20 windows on both 
sides of the relative position (40 total) of 50 bp size, corresponding to 
2 kb total (+/- 1 kb). Windows without a value may be interpolated 
(smoothed) from neigboring values, if available.

Stranded data may be collected. If the feature does not have an inherent 
strand, one may be specified to enforce stranded collection or a particular 
orientation. 

When features overlap, or the collection windows of one feature overlaps 
with another feature, then data may be ignored and not collected (--avoid).

=head1 EXAMPLES

These are some examples of some common scenarios for collecting data.

=over 4

=item Collect scores in intervals around start

You want to collect the mean score from a bigWig file in twenty 50 bp 
intervals flanking the start position of each feature in Bed file.

  get_relative_data.pl --data scores.bw --in input.bed

=item Collect scores in intervals around middle

You want to collect median scores in 20 bp intervals extending 500 bp 
from the midpoint of each feature.

  get_relative_data.pl --win 20 --num 25 --pos m --data scores.bw --in \
  input.txt

=item Collect scores in intervals from annotation database

You want to collect scores in intervals around the transcription start 
site of genes in an annotation database, but also avoid intervals that 
may overlap neighboring genes. You want to collect alignment counts 
from a Bam file in a stranded fashion. 

  get_relative_data.pl --db annotation --feature gene --avoid --strand \
  sense --method count --data alignments.bam --out gene_tss
    
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
