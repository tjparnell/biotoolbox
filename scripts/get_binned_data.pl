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
	use_minimum_mapq
	$BAM_ADAPTER
	$BIG_ADAPTER
);
use Bio::ToolBox::utility qw( format_with_commas simplify_dataset_name );
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};

our $VERSION = '2.03';

print "\n This script will collect binned values across features\n\n";

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
	$infile,         $parse,      $outfile,       $main_database,
	$data_database,  $feature,    $subfeature,    $exon_subfeature,
	$method,         $stranded,   $bins,          $extension,
	$extension_size, $min_length, $long_data,     $smooth,
	$sum,            $format,     $groupcol,      $gz,
	$cpu,            $help,       $print_version, $mapq
);    # command line variables
my @datasets;

## Command line options
GetOptions(
	'i|in=s'         => \$infile,             # input file
	'parse!'         => \$parse,              # parse input file
	'o|out=s'        => \$outfile,            # name of outfile
	'd|db=s'         => \$main_database,      # main or annotation database name
	'D|ddb=s'        => \$data_database,      # data database
	'a|data=s'       => \@datasets,           # dataset name
	'f|feature=s'    => \$feature,            # what type of feature to work with
	'u|subfeature=s' => \$subfeature,         # indicate to restrict to subfeatures
	'exons!'         => \$exon_subfeature,    # old parameter
	'm|method=s'     => \$method,             # method for collecting the data
	't|strand=s'     => \$stranded,       # indicate whether stranded data should be taken
	'b|bins=i'       => \$bins,           # number of bins
	'x|ext=i'        => \$extension,      # number of bins to extend beyond the feature
	'X|extsize=i'    => \$extension_size, # explicit size of extended bins
	'min=i'          => \$min_length,     # minimum feature size
	'long!'          => \$long_data,      # collecting long data features
	'smooth!'        => \$smooth,         # do not interpolate over missing values
	'U|sum!'         => \$sum,            # determine a final average for all the features
	'r|format=i'     => \$format,         # decimal formatting
	'mapq=i'         => \$mapq,           # minimum map quality
	'g|groups'       => \$groupcol,       # write group column file
	'z|gz!'          => \$gz,             # compress the output file
	'c|cpu=i'        => \$cpu,            # number of execution threads
	'h|help'         => \$help,           # print the help
	'v|version'      => \$print_version,  # print the version
	'bam=s'          => \$BAM_ADAPTER,    # explicitly set the bam adapter
	'big=s'          => \$BIG_ADAPTER,    # explicitly set the big adapter
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
	print " Biotoolbox script get_binned_data.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

### Check for required values
my $formatter;
check_defaults();
my $start_time = time;
my $length_i;    # global value for the merged transcript length

### Generate or load the dataset
my $Data;
if ($infile) {
	$Data = Bio::ToolBox::Data->new(
		file       => $infile,
		parse      => $parse,
		feature    => $feature,
		subfeature => $subfeature,
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
	if ( $parse == 0 ) {
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

# the number of columns already in the data array
my $beginningcolumn = $Data->last_column + 1;
my $startcolumn;    # this is now calculated separately for each datasets

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

# Open data database
my $ddb;
if ( defined $data_database ) {

	# specifically defined a data database
	$ddb = open_db_connection($data_database)
		or die "unable to establish data database connection to $data_database!\n";
}

# Check for the dataset
@datasets = verify_or_request_feature_types(
	'db'      => $ddb || $Data->database,
	'feature' => \@datasets,
	'prompt'  => " Enter the dataset(s) or feature type(s) from which \n"
		. " to collect data. Comma delimited or range is acceptable\n",
);
unless (@datasets) {
	die
" No verifiable dataset(s) provided. Check your file path, database, or dataset.\n";
}

## Collect the relative data
printf " Collecting $method data in %d bins from %s\n",
	( $bins + 2 * $extension ), join( ", ", @datasets );

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

## Generate summed data -
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

## Write main output
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

# write the column group file
if ($groupcol) {
	my $groupfile = $written_file;
	$groupfile =~ s/\.txt (?:\.gz)? $/.col_groups.txt/x;
	my $fh = Bio::ToolBox::Data->open_to_write_fh($groupfile);
	$fh->print("Name\tDataset\n");
	for my $i ( $beginningcolumn .. $Data->last_column ) {
		my $name    = $Data->name($i);
		my $dataset = $name;
		$dataset =~ s/: [\-\d%bp]+ $//x;
		$fh->print("$name\t$dataset\n");
	}
	$fh->close;
	print " wrote group column file $groupfile\n";
}

printf " Completed in %.1f minutes\n", ( time - $start_time ) / 60;

# done

#### Subroutines #######

sub check_defaults {
	unless ( $main_database or $infile ) {

		# a database must be defined
		# a tim data file used as an input file would also define one
		print STDERR " FATAL: You must define a database or input file!\n";
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

	# strand
	if ($stranded) {
		unless ( $stranded eq 'sense'
			or $stranded eq 'antisense'
			or $stranded eq 'all' )
		{
			print STDERR
" FATAL: '$stranded' is not recognized for strand\n Use --help for more information\n";
			exit 1;
		}
	}
	else {
		$stranded = 'all';    # default is no strand information
	}

	# method
	if ( defined $method ) {
		unless ( $method eq 'mean'
			or $method eq 'median'
			or $method eq 'sum'
			or $method eq 'min'
			or $method eq 'max'
			or $method eq 'stddev'
			or $method =~ /^\w?count/ )
		{
			print STDERR
" FATAL: '$method' is not recognized for method\n Use --help for more information\n";
			exit 1;
		}

	}
	else {
		$method = 'mean';
	}

	# subfeatures
	if ( $subfeature and $long_data ) {
		print STDERR
" FATAL: Long data collection is incompatible with subfeatures\n Use --help for more information\n";
		exit 1;
	}
	if ($exon_subfeature) {

		# legacy option
		$subfeature = 'exon';
	}
	if ( $subfeature and $subfeature !~ /^ (?: exon | cds | 5p_utr | 3p_utr )$/x ) {
		print STDERR
" FATAL: unrecognized subfeature option '$subfeature'! Use exon, cds, 5p_utr or 3p_utr\n";
		exit 1;
	}

	# assign default window bin values
	unless ( defined $bins ) {

		# dividing gene into 10 (10%) bins seems reasonable to me
		$bins = 10;
	}
	unless ( defined $extension ) {

		# default is no extension
		$extension = 0;
	}
	unless ( defined $smooth ) {

		# default is to not include smoothing
		$smooth = 0;
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

	# set minimum alignment mapping quality
	if ( $method =~ /count/ and $mapq ) {
		use_minimum_mapq($mapq);
	}
}

## Parallel execution for efficiency
sub parallel_execution {
	my $pm = Parallel::ForkManager->new($cpu);

	# generate base name for child processes
	my $child_base_name = $outfile . ".$PID";

	# Split the input data into parts and execute in parallel in separate forks
	for my $i ( 1 .. $cpu ) {
		$pm->start and next;

		#### In child ####

		# splice the data structure
		$Data->split_data( $i, $cpu );

		# re-open database objects to make them clone safe
		# pass second true to avoid cached database objects
		my $db = $Data->open_database(1);
		if ($data_database) {
			$ddb = open_db_connection( $data_database, 1 );
		}

		# collapse transcripts if needed
		if (    $feature
			and $feature =~ /gene/i
			and $subfeature
			and $subfeature =~ m/exon | intron/xi )
		{
			$Data->collapse_gene_transcripts;
		}

		# calculate length if necessary
		if ($subfeature) {
			$length_i = $Data->add_transcript_length($subfeature);
		}

		# work through each dataset
		foreach my $dataset (@datasets) {

			# new start column for this dataset
			$startcolumn = $Data->number_columns + 1;

			# Prepare the metadata and header names
			my $binsize = ( 100 / $bins );
			prepare_bins( $binsize, $dataset );

			# Collect the data
			if ($long_data) {
				collect_binned_long_data( $binsize, $dataset );
			}
			else {
				collect_binned_data( $binsize, $dataset );
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

	# collapse transcripts if needed
	if (    $feature
		and $feature =~ /gene/i
		and $subfeature
		and $subfeature =~ m/exon | intron/xi )
	{
		$Data->collapse_gene_transcripts;
	}

	# calculate length if necessary
	if ($subfeature) {
		$length_i = $Data->add_transcript_length($subfeature);
	}

	# work through each dataset
	foreach my $dataset (@datasets) {

		# new start column for this dataset
		$startcolumn = $Data->number_columns + 1;

		# Prepare the metadata and header names
		my $binsize = ( 100 / $bins );
		prepare_bins( $binsize, $dataset );

		# Collect the data
		if ($long_data) {
			collect_binned_long_data( $binsize, $dataset );
		}
		else {
			collect_binned_data( $binsize, $dataset );
		}

		# Interpolate values
		if ($smooth) {
			print " Smoothing data by interpolation....\n";
			go_interpolate_values();
		}
	}
}

sub collect_binned_data {
	my ( $binsize, $dataset ) = @_;

	## Collect the data
	my $stream = $Data->row_stream;
	while ( my $row = $stream->next_row ) {

		# define the starting and ending points based on gene length
		my $length = defined $length_i ? $row->value($length_i) : $row->length;

		# check the length
		if ( defined $min_length and $length < $min_length ) {

			# this feature is too short to divided into bins
			# we will skip this feature

			# but first, put in null values
			for my $c ( $startcolumn .. ( $Data->last_column ) ) {
				$row->value( $c, '.' );
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
				$extra = int( ( $extension * $binsize * 0.01 * $length ) + 0.5 );
			}
		}

		my $regionscores = $row->get_region_position_scores(
			'ddb'        => $ddb,
			'dataset'    => $dataset,
			'method'     => $method,
			'extend'     => $extra,
			'stranded'   => $stranded,
			'subfeature' => $subfeature,
			'length'     => $length,
		);

		# record the scores for each bin
		record_the_bin_values( $row, $length, $regionscores );
	}
}

sub collect_binned_long_data {
	my ( $binsize, $dataset ) = @_;

	## Collect the data
	my $stream = $Data->row_stream;
	while ( my $row = $stream->next_row ) {

		# identify the feature or use the row
		my $fstart = $row->start;
		my $fstop  = $row->end;
		my $strand = $row->strand;
		my $length = $row->length;   # subfeatures not allowed here, so use feature length

		# collect the scores to the bins in the region
		for my $column ( $startcolumn .. ( $Data->last_column ) ) {

			# we will step through each data column, representing each window (bin)
			# across the feature's region
			# any scores within this window will be collected and the mean
			# value reported

			# convert the window start and stop coordinates (as percentages) to
			# actual bp
			# this depends on whether the binsize is explicitly defined in bp or
			# is a fraction of the feature length
			my ( $start, $stop );
			if ( $Data->metadata( $column, 'bin_size' ) =~ /bp$/ ) {

				# the bin size is explicitly defined

				# the start and stop points are relative to either the feature
				# start (always 0) or the end (the feature length), depending
				# upon whether the 5' or 3' end of the feature

				# determine this by the sign of the start position
				if ( $Data->metadata( $column, 'start' ) < 0 and $strand >= 0 ) {

					# the start position is less than 0, implying the 5' end
					# the reference position will be the feature start on plus strand
					$start = $fstart + $Data->metadata( $column, 'start' );
					$stop  = $fstart + $Data->metadata( $column, 'stop' );
				}
				elsif ( $Data->metadata( $column, 'start' ) < 0 and $strand < 0 ) {

					# the start position is less than 0, implying the 5' end
					# the reference position will be the feature end on minus strand
					$start = $fstop - $Data->metadata( $column, 'start' );
					$stop  = $fstop - $Data->metadata( $column, 'stop' );
				}
				elsif ( $Data->metadata( $column, 'start' ) >= 0 and $strand >= 0 ) {

					# the start position is greather than 0, implying the 3' end
					# the reference position will be the feature start on plus strand
					$start = $fstop + $Data->metadata( $column, 'start' );
					$stop  = $fstop + $Data->metadata( $column, 'stop' );
				}
				elsif ( $Data->metadata( $column, 'start' ) >= 0 and $strand < 0 ) {

					# the start position is greather than 0, implying the 3' end
					# the reference position will be the feature end on minus strand
					$start = $fstart - $Data->metadata( $column, 'start' );
					$stop  = $fstart - $Data->metadata( $column, 'stop' );
				}
				else {
					warn " unable to unable to identify region orientation: start "
						. $Data->metadata( $column, 'start' )
						. ", strand $strand\n";
					return;
				}
			}
			else {
				# otherwise the bin size is based on feature length
				if ( $strand >= 0 ) {

					# forward plus strand
					$start =
						int( $fstart +
							( $Data->metadata( $column, 'start' ) * 0.01 * $length ) +
							0.5 );
					$stop =
						int( $fstart +
							( $Data->metadata( $column, 'stop' ) * 0.01 * $length ) - 1 +
							0.5 );
				}
				else {
					# reverse minus strand
					$start =
						int( $fstop -
							( $Data->metadata( $column, 'start' ) * 0.01 * $length ) +
							0.5 );
					$stop =
						int( $fstop -
							( $Data->metadata( $column, 'stop' ) * 0.01 * $length ) + 1 +
							0.5 );
				}
			}

			# collect the data for this bin
			my $score = $row->get_score(
				'ddb'      => $ddb,
				'dataset'  => $dataset,
				'start'    => $start,
				'stop'     => $stop,
				'method'   => $method,
				'stranded' => $stranded,
			);
			if ( $formatter and looks_like_number($score) ) {
				$score = sprintf( $formatter, $score );
			}
			$row->value( $column, $score );
		}
	}
}

sub record_the_bin_values {

	# get the passed values
	my ( $row, $length, $regionscores ) = @_;

	# assign the scores to the bins in the region
	for my $column ( $startcolumn .. ( $Data->last_column ) ) {

		# we will step through each data column, representing each window (bin)
		# across the feature's region
		# any scores within this window will be collected and the mean
		# value reported

		# record nulls if no data returned
		unless ( scalar keys %{$regionscores} ) {
			$row->value( $column, calculate_score( $method, undef ) );
			next;
		}

		# convert the window start and stop coordinates (as percentages) to
		# actual bp
		# this depends on whether the binsize is explicitly defined in bp or
		# is a fraction of the feature length
		my ( $start, $stop );
		if ( $Data->metadata( $column, 'bin_size' ) =~ /bp$/ ) {

			# the bin size is explicitly defined

			# the start and stop points are relative to either the feature
			# start (always 0) or the end (the feature length), depending
			# upon whether the 5' or 3' end of the feature

			# determine this by the sign of the start position
			if ( $Data->metadata( $column, 'start' ) < 0 ) {

				# the start position is less than 0, implying the 5' end
				# the reference position will be the feature start, or 0
				$start = $Data->metadata( $column, 'start' );
				$stop  = $Data->metadata( $column, 'stop' );
			}
			else {
				# the start position is greather than 0, implying the 3' end
				# the reference position will be the feature end, or length
				$start = $Data->metadata( $column, 'start' ) + $length;
				$stop  = $Data->metadata( $column, 'stop' ) + $length;
			}
		}
		else {
			# otherwise the bin size is based on feature length
			$start = sprintf "%.0f",
				( $Data->metadata( $column, 'start' ) * 0.01 * $length ) + 1;
			$stop = sprintf "%.0f",
				( $Data->metadata( $column, 'stop' ) * 0.01 * $length );
		}

		# collect the scores for this window
		my @scores = map { $regionscores->{$_} }
			grep { $_ >= $start and $_ <= $stop }
			keys %{$regionscores};

		# calculate the value
		my $window_score = calculate_score( $method, \@scores );
		if ( $formatter and looks_like_number($window_score) ) {
			$window_score = sprintf( $formatter, $window_score );
		}

		# record the value
		$row->value( $column, $window_score );
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
				for ( my $i = $col + 1; $col <= $lastwindow; $i++ ) {
					if ( $row->value($i) ne '.' ) {
						$next_i = $i;
						last;
					}
				}
				next unless defined $next_i;

				# determine fractional value
				my $initial = $row->value( $col - 1 );
				my $fraction =
					( $row->value($next_i) - $initial ) / ( $next_i - $col + 1 );

				# apply fractional values
				for ( my $n = $col; $n < $next_i; $n++ ) {
					my $score = $initial + ( $fraction * ( $n - $col + 1 ) );
					if ( $formatter and looks_like_number($score) ) {
						$score = sprintf( $formatter, $score );
					}
					$row->value( $n, $score );
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

	my ( $binsize, $dataset ) = @_;

	# the size of the bin in percentage units, default would be 10%
	# each bin will be titled the starting and ending point for that bin in
	# percentage units
	# for example, -20..-10,-10..0,0..10,10..20

	# if $extension is defined, then it will add the appropriate flanking bins,
	# otherwise it should skip them

	# bin(s) on 5' flank
	if ($extension) {

		if ($extension_size) {

			# extended bins should be of specific bp size
			for ( my $i = $extension; $i > 0; $i-- ) {
				my $start = 0 - ( $extension_size * $i );
				my $stop  = 0 - ( $extension_size * ( $i - 1 ) );
				_set_metadata( $start, $stop, $extension_size, 'bp', $dataset );
			}
		}
		else {
			# extended bin size will be based on feature length
			for ( my $i = $extension; $i > 0; $i-- ) {
				my $start = 0 - ( $binsize * $i );
				my $stop  = 0 - ( $binsize * ( $i - 1 ) );
				_set_metadata( $start, $stop, $binsize, '%', $dataset );
			}
		}
	}

	# bins over the gene body
	for ( my $i = 0; $i < $bins; $i++ ) {
		my $start = ( $i * $binsize );
		my $stop  = ( $i + 1 ) * $binsize;
		_set_metadata( $start, $stop, $binsize, '%', $dataset );
	}

	# bin(s) on 3' flank
	if ($extension) {

		if ($extension_size) {

			# extended bins should be of specific bp size
			for ( my $i = 0; $i < $extension; $i++ ) {
				my $start = ( $extension_size * $i );
				my $stop  = ( $extension_size * ( $i + 1 ) );
				_set_metadata( $start, $stop, $extension_size, 'bp', $dataset );
			}
		}
		else {
			# extended bin size will be based on feature length
			for ( my $i = 0; $i < $extension; $i++ ) {
				my $start = 100 + ( $binsize * $i );
				my $stop  = 100 + ( $binsize * ( $i + 1 ) );
				_set_metadata( $start, $stop, $binsize, '%', $dataset );
			}
		}
	}
}

### Set the metadata for a new data table column (dataset)
sub _set_metadata {

	# the start and stop positions are passed
	my ( $start, $stop, $binsize, $unit, $dataset ) = @_;

	# generate a simplified new name
	my $new_name = simplify_dataset_name($dataset);

	# set new index
	my $new_index = $Data->add_column( sprintf( "%s:%d%s", $new_name, $start, $unit ) );

	# set the metadata using passed and global variables
	# set the metadata
	$Data->metadata( $new_index, 'start',          $start );
	$Data->metadata( $new_index, 'stop',           $stop );
	$Data->metadata( $new_index, 'dataset',        $dataset );
	$Data->metadata( $new_index, 'method',         $method );
	$Data->metadata( $new_index, 'bin_size',       $binsize . $unit );
	$Data->metadata( $new_index, 'strand',         $stranded );
	$Data->metadata( $new_index, 'decimal_format', $format ) if defined $format;

	if ($data_database) {
		$Data->metadata( $new_index, 'db', $data_database );
	}

	if ( $dataset =~ / \. (b|cr) am $/xni and $method =~ /count/ ) {
		$Data->metadata( $index, 'mapq', $mapq );
	}
}

__END__

=head1 NAME

get_binned_data.pl

A program to collect data in bins across a list of features.

=head1 SYNOPSIS
 
 get_binned_data.pl [--options] --in <filename> --out <filename>
  
 get_binned_data.pl [--options] -i <filename> <data1> <data2...>
  
  Options for data files:
  -i --in <filename>                  input file: txt bed gff gtf refFlat ucsc
  -o --out <filename>                 optional output file, default overwrite 
  
  Options for new files:
  -d --db <name>                      annotation database: mysql sqlite
  -f --feature <type>                 one or more feature types from db or gff
  
  Options for data collection:
  -D --ddb <name|file>                data or BigWigSet database
  -a --data <dataset|filename>        data from which to collect: bw bam etc
  -m --method [mean|median|stddev|    statistical method for collecting data
        min|max|range|sum|count|      default mean
        pcount|ncount]
  -t --strand [all|sense|antisense]   strand of data relative to feature (all)
  -u --subfeature [exon|cds|          collect over gene subfeatures 
        5p_utr|3p_utr] 
  --long                              collect each window independently
  -r --format <integer>               number of decimal places for numbers
  --mapq <integer>                    minimum map quality of counted alignments
  
  Bin specification:
  -b --bins <integer>                 number of bins feature is divided (10)
  -x --ext <integer>                  number of extended bind outside feature
  -X --extsize <integer>              size of extended bins
  --min <integer>                     minimum size of feature to divide
  
  Post-processing:
  -U --sum                            generate summary file
  --smooth                            smoothen sparse data
  
  General options:
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

Specify the output file name. Default is to overwrite the input text 
file. Required if generating a new file from a database.

=back

=head2 Options for new files

=over 4

=item --db E<lt>nameE<gt>

Specify the name of a L<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be collected rather than providing 
an input file. The name may be that of a MySQL database or a SQLite file. 

=item --feature <type | type:source | alias>,...

Specify the type of feature from which to collect values. This is required 
only for new feature tables. Three types of values may be passed: the 
feature type, feature type and source expressed as 'type:source'. 
More than one feature may be included as a comma-delimited list (no spaces). 
  
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

=item --strand [all|sense|antisense]

Specify whether stranded data should be collected. Three values are 
allowed: all datasets should be collected (default), only sense 
datasets, or only antisense datasets should be collected.

=item --subfeature [ exon | cds | 5p_utr | 3p_utr ]

Optionally specify the type of subfeature to collect from, rather than 
the entire gene. If the parent feature is gene and the subfeature is exon, 
then all transcripts of the gene will be collapsed. The other subfeatures 
(cds, 5p_utr, and 3p_utr) will not work with gene features but only with 
coding mRNA transcripts. Note that the long option is incompatible. 
Default is null. 

=item --exons

Legacy option for specifying --subfeature exon.

=item --long

Indicate that data should be collected independently for each long 
window. This may be enabled automatically if the sum of the entire 
window length passes a predefined threshold. The default for 'short' 
windows is to collect all of the point data from the dataset first, and 
then divide the results into the different windows. Datasets consisting 
of "long" features, for example long alignments, may be counted more 
than once in long mode when they span multiple windows. Not compatible 
when subfeatures are enabled.

=item --format E<lt>integerE<gt>

Specify the number of decimal positions to format the collected scores. 
Default is not to format, often leading to more than the intended 
significant digits.

=item --mapq E<lt>integer<gt>

Specify the minimum mapping quality of alignments to be considered when
counting from a Bam file. Default is 0, which will include all alignments,
including multi-mapping (typically MAPQ of 0). Set to an integer in range
of 0..255. Only affects count methods, including C<count>, C<ncount>, and
C<pcount>. Other methods involving coverage, e.g. C<mean>, do not filter
alignments.

=back

=head2 Bin specification

=over 4

=item --bins E<lt>integerE<gt>

Specify the number of bins that will be generated over the length 
of the feature. The size of the feature is a percentage of the 
feature length. The default number is 10, which results in bins of 
size equal to 10% of the feature length. 

=item --ext E<lt>integerE<gt>

Specify the number of extended bins on either side of the feature. 
The bins are of the same size as determined by the feature 
length and the --bins value. The default is 0. 

=item --extsize E<lt>integerE<gt>

Specify the exact bin size in bp of the extended bins rather than
using a percentage of feature of length.

=item --min E<lt>integerE<gt>

Specify the minimum feature size to be averaged. Features with a
length below this value will not be skipped (all bins will have
null values). This is to avoid having bin sizes below the average 
microarray tiling distance. The default is undefined (no limit).

=back

=head2 Post-processing

=over 4

=item --sum

Indicate that the data should be averaged across all features at
each position, suitable for graphing. A separate text file will be
written with the suffix '_summed' with the averaged data. The default 
is false.

=item --smooth

Indicate that windows without values should (not) be interpolated
from neighboring values. The default is false.

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
default is 4. Disable multi-threaded execution by setting to 1. 

=item --noparse

Prevent input annotation files from being automatically parsed into sequence 
features. Coordinates will be used as is and new data columns will be appended 
to the input file. 

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
in an annotation database.

  get_binned_data.pl --db annotation --feature gene --strand sense \
  --method count --data alignments.bam --out gene_profile --sum
  
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
