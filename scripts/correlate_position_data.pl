#!/usr/bin/perl

# documentation at end of file

use warnings;
use strict;
use English qw(-no_match_vars);
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use List::Util qw(sum0 min max);
use Statistics::Lite qw(median stddevp);
use Statistics::Descriptive;
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
);
use Bio::ToolBox::utility qw(format_with_commas);

my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};

# ANOVA Parameters for P-value Statistics
my $anova;
eval {
	# check for ANOVA support
	require Statistics::ANOVA;
	$anova = 1;
};
my $INDEPENDENCE = 0;
my $PARAMETRIC   = 1;
my $ORDINAL      = 0;

use constant LOG10 => log(10);

our $VERSION = '1.70';

print "\n This program will correlate positions of occupancy between two datasets\n\n";

### Quick help
unless (@ARGV) {

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
my (
	$infile,     $outfile,     $database,    $data_database,
	$refDataSet, $testDataSet, $find_pvalue, $find_shift,
	$radius,     $position,    $set_strand,  $norm_method,
	$gz,         $cpu,         $help,        $print_version,
);

# Command line options
GetOptions(
	'i|in=s'                  => \$infile,           # the input data file
	'o|out=s'                 => \$outfile,          # name of output file
	'd|db=s'                  => \$database,         # name of database
	'D|ddb=s'                 => \$data_database,    # database containing datasets
	'r|ref=s'                 => \$refDataSet,       # reference dataset
	't|test=s'                => \$testDataSet,      # test dataset
	'pval!'                   => \$find_pvalue,      # calculate student t-test
	'shift!'                  => \$find_shift,       # calculate optimum shift
	'radius=i'                => \$radius,           # for collecting data when shifting
	'p|pos=s'                 => \$position,         # set the relative feature position
	'force_strand|set_strand' => \$set_strand,       # enforce an artificial strand
		# force_strand is preferred option, but respect the old option
	'norm=s'    => \$norm_method,      # method of normalization
	'z|gz!'     => \$gz,               # compress output
	'c|cpu=i'   => \$cpu,              # number of execution threads
	'h|help'    => \$help,             # request help
	'v|version' => \$print_version,    # print the version
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
	print " Biotoolbox script correlate_position_data.pl, version $VERSION\n\n";
	exit;
}

### Check for requirements
check_defaults();

### Open input file
my $Data = Bio::ToolBox::Data->new( file => $infile )
	or die " unable to open input file '$infile'!\n";
printf " Loaded %s features from $infile.\n", format_with_commas( $Data->last_row );
$Data->program("$PROGRAM_NAME, v $VERSION");

### Open database connection
if ($database) {

	# update or add database as required
	if ( $Data->database and $database ne $Data->database ) {

		# update with new database
		printf " updating main database name from '%s' to '%s'\n",
			$Data->database, $database;
		print "   Re-run without --db option if you do not want this to happen\n";
		$Data->database($database);
	}
	elsif ( not $Data->database ) {
		$Data->database($database);
	}
}
my $ddb;
if ( defined $data_database ) {

	# specifically defined a data database
	$ddb = open_db_connection($data_database)
		or die "unable to establish data database connection to $data_database!\n";
}

### Validate input datasets
validate_or_request_dataset();

### Collect the data correlations
my $start_time = time;

# check whether it is worth doing parallel execution
if ( $cpu > 1 ) {
	while ( $cpu > 1 and $Data->last_row / $cpu < 100 ) {

		# We need at least 1000 lines in each fork split to make
		# it worthwhile to do the split, otherwise, reduce the number of
		# splits to something more worthwhile
		$cpu--;
	}
}

if ( $cpu > 1 ) {

	# parallel execution
	parallel_execution();
}
else {
	# single process execution
	single_execution();
}

### Write output file
{
	my $success = $Data->write_file(
		'filename' => $outfile,
		'gz'       => $gz,
	);
	if ($success) {
		print " wrote output file $success\n";
	}
	else {
		print " unable to write output file!\n";
	}
}
printf " Finished in %.2f minutes\n", ( time - $start_time ) / 60;

########################   Subroutines   ###################################

sub check_defaults {

	unless ($infile) {
		if (@ARGV) {
			$infile = shift @ARGV;
		}
		else {
			print STDERR " FATAL: no input file! use --help for more information\n";
			exit 1;
		}
	}
	unless ( defined $gz ) {
		$gz = 0;
	}

	if ( defined $position ) {

		# check the position value
		unless ( $position eq 'm'
			or $position == 5
			or $position == 3
			or $position == 4 )
		{
			print STDERR " FATAL: Unknown relative position '$position'!\n";
			exit 1;
		}
		if ( $position eq 'm' ) { $position = 4 }    # change to match internal usage
	}
	else {
		# default position to use the midpoint
		$position = 4;
	}

	if ($norm_method) {
		unless ( $norm_method =~ /^rank|sum$/i ) {
			print STDERR " FATAL: Unrecognized normalization method! see help\n";
			exit 1;
		}
	}

	if ($find_pvalue) {
		unless ($anova) {
			print STDERR
" FATAL: Cannot calculate requested P-value! Pleast install Statistics::ANOVA\n";
			exit 1;
		}
	}

	if ($find_shift) {
		unless ($radius) {
			print STDERR
" FATAL: Must define a radius and reference point to calculate optimum shift\n";
			exit 1;
		}
	}

	if ($parallel) {
		$cpu ||= 4;
	}
	else {
		# disable cores
		print " WARNING: disabling parallel CPU execution, no support present\n" if $cpu;
		$cpu = 0;
	}

	unless ($outfile) {
		$outfile = $infile;
	}

}

sub validate_or_request_dataset {

	# Process the reference dataset
	$refDataSet = verify_or_request_feature_types(
		'db'      => $ddb || $Data->database,
		'feature' => $refDataSet,
		'prompt'  => "Enter the reference data set  ",
		'single'  => 1,
	);
	unless ($refDataSet) {
		die " A valid reference data set must be provided!\n";
	}

	# Process the test dataset
	$testDataSet = verify_or_request_feature_types(
		'db'      => $ddb || $Data->database,
		'feature' => $testDataSet,
		'prompt'  => "Enter the test data set  ",
		'single'  => 1,
	);
	unless ($testDataSet) {
		die " A valid test data set must be provided!\n";
	}
}

sub parallel_execution {
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_start( sub { sleep 1; } );

	# give a chance for child to start up and open databases, files, etc
	# without creating race conditions

	# Prepare new columns
	my ( $r_i, $p_i, $shift_i, $shiftr_i ) = add_new_columns();

	# generate base name for child processes
	my $child_base_name = $outfile . ".$PID";

	# variables for reporting the summary of results
	my @correlations;
	my @optimal_shifts;
	my @optimal_correlations;
	my $count           = 0;
	my $not_enough_data = 0;
	my $no_variance     = 0;
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $result ) = @_;
			push @correlations,         @{ $result->[0] };
			push @optimal_shifts,       @{ $result->[1] };
			push @optimal_correlations, @{ $result->[2] };
			$count           += $result->[3];
			$not_enough_data += $result->[4];
			$no_variance     += $result->[5];
		}
	);

	# Split the input data into parts and execute in parallel in separate forks
	print " Forking into $cpu children for parallel correlation collections\n";
	for my $i ( 1 .. $cpu ) {
		$pm->start and next;

		#### In child ####

		# splice the data structure
		$Data->splice_data( $i, $cpu );

		# re-open database objects to make them clone safe
		my $db = $Data->open_database(1);
		if ($data_database) {
			$ddb = open_db_connection( $data_database, 1 );
		}

		# collect
		my @results = collect_correlations( $r_i, $p_i, $shift_i, $shiftr_i );

		# write output file
		my $success = $Data->write_file(
			'filename' => sprintf( "$child_base_name.%03s", $i ),
			'gz'       => 0,    # faster to write without compression
		);
		if ($success) {
			print " wrote output file $success\n";
		}
		else {
			# failure! the subroutine will have printed error messages
			die " unable to write file!\n";

			# no need to continue
		}

		# Finished
		$pm->finish( 0, \@results );
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
	my $new_count = $Data->reload_children(@files);
	printf " reloaded %s features from children\n", format_with_commas($new_count);

	# summarize the results
	summarize_results( \@correlations, \@optimal_shifts, \@optimal_correlations,
		$count, $not_enough_data, $no_variance );
}

sub single_execution {

	# run in a single thread

	# Prepare new columns
	my ( $r_i, $p_i, $shift_i, $shiftr_i ) = add_new_columns();

	# collect
	print " Collecting correlations....\n";
	my @results = collect_correlations( $r_i, $p_i, $shift_i, $shiftr_i );

	# summarize
	summarize_results(@results);
}

sub collect_correlations {
	my ( $r_i, $p_i, $shift_i, $shiftr_i ) = @_;

	# Set variables for summary analysis
	my @correlations;
	my @optimal_shifts;
	my @optimal_correlations;
	my $count           = 0;
	my $not_enough_data = 0;
	my $no_variance     = 0;

	# check that we can collect information
	unless ( $Data->feature_type eq 'named' or $Data->feature_type eq 'coordinate' ) {
		die " Unable to identify the type of features in the input file."
			. " File must be a recognizable format and/or have column header names\n"
			. " for chromosome, start, stop; or named database features\n";
	}

	# Walk through features
	my $stream = $Data->row_stream;
	while ( my $row = $stream->next_row ) {

		# Collect positioned data
		my %ref_pos2data;
		my %test_pos2data;
		if ($radius) {

			# collecting data in a radius around a reference point

			# collect data
			%ref_pos2data = $row->get_relative_point_position_scores(
				ddb      => $ddb,
				dataset  => $refDataSet,
				position => $position,
				extend   => ( 2 * $radius ),
			);
			%test_pos2data = $row->get_relative_point_position_scores(
				ddb      => $ddb,
				dataset  => $testDataSet,
				position => $position,
				extend   => ( 2 * $radius ),
			);
		}
		else {
			# just collect data over the region
			# collect data
			%ref_pos2data = $row->get_region_position_scores(
				ddb      => $ddb,
				dataset  => $refDataSet,
				position => $position,
			);
			%test_pos2data = $row->get_region_position_scores(
				ddb      => $ddb,
				dataset  => $testDataSet,
				position => $position,
			);
		}

		# Verify minimum data count
		if (   sum0( map {abs} values %ref_pos2data ) == 0
			or sum0( map {abs} values %test_pos2data ) == 0 )
		{
			# not enough data points to work with
			$row->value( $r_i, '.' );
			$row->value( $p_i, '.' ) if ($find_pvalue);
			if ($find_shift) {
				$row->value( $shift_i,  '.' );
				$row->value( $shiftr_i, '.' );
			}
			$not_enough_data++;
			next;
		}

		# Calculate Paired T-Test P value
		# going to try this before normalizing
		if ($find_pvalue) {
			$row->value( $p_i, calculate_ttest( \%ref_pos2data, \%test_pos2data ) );
		}

		# Normalize the data
		# there are couple ways to do this
		if ( $norm_method and $norm_method =~ /rank/i ) {

			# convert to rank values
			# essentially a Spearman's rank correlation
			normalize_values_by_rank( \%ref_pos2data );
			normalize_values_by_rank( \%test_pos2data );
		}
		elsif ( $norm_method and $norm_method =~ /sum/i ) {

			# normalize the sums of both
			# good for read counts
			normalize_values_by_sum( \%ref_pos2data );
			normalize_values_by_sum( \%test_pos2data );
		}

		# Collect the (normalized) data from the position data
		# this may be all of the data or a portion of it, depending
		my @ref_data;
		my @test_data;
		if ($radius) {
			for my $i ( ( 0 - $radius ) .. $radius ) {

				# get values for each position, default 0
				push @ref_data,  $ref_pos2data{$i}  || 0;
				push @test_data, $test_pos2data{$i} || 0;
			}
		}
		else {
			for ( my $i = 1; $i < $row->length; $i++ ) {

				# get values for each position, default 0
				push @ref_data,  $ref_pos2data{$i}  || 0;
				push @test_data, $test_pos2data{$i} || 0;
			}
		}

		# verify the reference data
		if ( min(@ref_data) == max(@ref_data) ) {

			# no variance in the data! cannot calculate Pearson
			$row->value( $r_i, '.' );
			$row->value( $p_i, '.' ) if ($find_pvalue);
			if ($find_shift) {
				$row->value( $shift_i,  '.' );
				$row->value( $shiftr_i, '.' );
			}
			$no_variance++;
			next;
		}

		# Calculate the Pearson correlation between test and reference
		my $r;
		if ( min(@test_data) == max(@test_data) ) {

			# test has no variance
			# report a Pearson of 0
			$r = 0;

			# this may not be a total loss, as there may be data points
			# outside of the window that would generate an optimal Pearson
			# so we will continue rather than bailing
		}
		else {
			# test has variance, calculate Pearson
			$r = calculate_pearson( \@ref_data, \@test_data ) || 0;
		}
		$row->value( $r_i, $r );
		push @correlations, $r if $r;

		# Determine optimal shift if requested
		if ($find_shift) {

			# find optimal shift
			my ( $best_shift, $best_r ) =
				calculate_optimum_shift( \%test_pos2data, \@ref_data, $r, );

			# change strand
			if ($set_strand) {

				# user is imposing a strand
				if ( $row->strand < 0 ) {

					# only change shift direction if reverse strand
					$best_shift *= -1;
				}
			}

			# record in data table
			$row->value( $shift_i,  $best_shift );
			$row->value( $shiftr_i, $best_r > $r ? $best_r : '.' );

			# record for summary analyses
			push @optimal_shifts,       $best_shift if ( $best_r > $r );
			push @optimal_correlations, $best_r     if ( $best_r > $r );
		}

		# Finished with this row
		$count++;
	}

	# Return correlation results for summary
	return ( \@correlations, \@optimal_shifts, \@optimal_correlations,
		$count, $not_enough_data, $no_variance );
}

sub add_new_columns {

	# Required columns
	# the Pearson column
	my $r_i = $Data->add_column('Pearson_correlation');
	$Data->metadata( $r_i, 'reference', $refDataSet );
	$Data->metadata( $r_i, 'test',      $testDataSet );

	# extra information
	$Data->metadata( $r_i, 'data_database', $data_database ) if $data_database;
	$Data->metadata( $r_i, 'normalization', $norm_method )   if $norm_method;
	if ($radius) {
		$Data->metadata( $r_i, 'radius', $radius );
		$Data->metadata( $r_i, 'reference_point',
			$position == 5 ? '5_prime' : $position == 4 ? 'midpoint' : '3_prime' );
	}

	# Optional columns if we're calculating an optimal shift
	my ( $shift_i, $shiftr_i );
	if ($find_shift) {

		# optimal test-reference shift value
		$shift_i = $Data->add_column('optimal_shift');
		$Data->metadata( $shift_i, 'reference',     $refDataSet );
		$Data->metadata( $shift_i, 'test',          $testDataSet );
		$Data->metadata( $shift_i, 'radius',        $radius );
		$Data->metadata( $shift_i, 'data_database', $data_database ) if $data_database;
		$Data->metadata( $shift_i, 'reference_point',
			$position == 5 ? '5_prime' : $position == 4 ? 'midpoint' : '3_prime' );

		# optimal test-reference shift value
		$shiftr_i = $Data->add_column('shift_correlation');
		$Data->metadata( $shiftr_i, 'reference',     $refDataSet );
		$Data->metadata( $shiftr_i, 'test',          $testDataSet );
		$Data->metadata( $shiftr_i, 'radius',        $radius );
		$Data->metadata( $shiftr_i, 'data_database', $data_database ) if $data_database;
		$Data->metadata( $shiftr_i, 'reference_point',
			$position == 5 ? '5_prime' : $position == 4 ? 'midpoint' : '3_prime' );
		$Data->metadata( $shiftr_i, 'normalization', $norm_method ) if $norm_method;
	}

	# Paired Student T-Test P-value
	my $p_i;
	if ($find_pvalue) {
		$p_i = $Data->add_column('ANOVA_Pval');
		$Data->metadata( $p_i, 'reference', $refDataSet );
		$Data->metadata( $p_i, 'test',      $testDataSet );
		$Data->metadata( $p_i, 'transform', '-Log10' );
	}

	# record the indices
	return ( $r_i, $p_i, $shift_i, $shiftr_i );
}

sub calculate_optimum_shift {

	# Calculate an optimal shift that maximizes the Pearson correlation

	my ( $test_pos2data, $ref_data, $r ) = @_;

	# we will calculate a Pearson correlation for each shift
	# calculating window of 2*radius, shifting 2*radius across region
	my %shift2pearson;
	for my $shift ( 1 .. $radius ) {

		# we will test shifts incrementally from 0 to radius in both directions

		## forward shift
		my @test_data;
		for my $i ( ( 0 - $radius + $shift ) .. ( $radius + $shift ) ) {
			push @test_data, $test_pos2data->{$i} || 0;
		}

		# check the test array
		if ( min(@test_data) == max(@test_data) ) {

			# no variance!
			next;
		}

		# calculate Pearson correlation
		my $pearson = calculate_pearson( $ref_data, \@test_data );
		if ( defined $pearson and $pearson > 0.5 ) {
			$shift2pearson{$shift} = $pearson;
		}

		## Reverse shift
		@test_data = ();
		for my $i ( ( 0 - $radius - $shift ) .. ( $radius - $shift ) ) {
			push @test_data, $test_pos2data->{$i} || 0;
		}

		# check the test array
		if ( min(@test_data) == max(@test_data) ) {

			# no variance!
			next;
		}

		# calculate Pearson correlation for reverse shift
		my $pearson1 = calculate_pearson( $ref_data, \@test_data );
		if ( defined $pearson1 and $pearson1 > 0.5 ) {
			$shift2pearson{ $shift * -1 } = $pearson1;
		}
	}

	# identify the best shift
	my $best_r     = $r;
	my $best_shift = 1000000;    # an impossibly high number
	foreach my $shift ( keys %shift2pearson ) {
		if ( $shift2pearson{$shift} == $best_r ) {

			# same as before
			if ( $best_shift != 1000000 and abs($shift) < abs($best_shift) ) {

				# if it's a lower shift, then keep it,
				# but not if it's original
				$best_r     = $shift2pearson{$shift};
				$best_shift = $shift;
			}
		}
		elsif ( $shift2pearson{$shift} > $best_r ) {

			# a good looking r value
			# but we only want those with the smallest shift
			if ( abs($shift) < abs($best_shift) ) {
				$best_r     = $shift2pearson{$shift};
				$best_shift = $shift;
			}
		}
	}
	if ( $best_shift == 1000000 ) {

		# no good shift found!? then that means no shift
		$best_shift = 0;
	}

	# done
	return ( $best_shift, $best_r );
}

sub normalize_values_by_rank {

	# convert values to rank value

	my $data = shift;
	my %ranks;    # a hash for the quantile lookup

	# put values into quantile lookup hash
	my $i = 1;
	foreach ( sort { $a <=> $b } values %{$data} ) {

		# for the quantile hash, the key will be the original data value
		# and the hash value will be the rank in the ordered data array
		# the rank is stored in anonymous array since there may be more
		# than one identical data value and we want the ranks for all of
		# them
		push @{ $ranks{$_} }, $i;
		$i++;
	}

	# replace the current value with the scaled quantile value
	foreach my $k ( keys %{$data} ) {

		# usually there is only one rank value for each data value
		# sometimes, there are two or more identical data values within the
		# the dataset, in which case we will simply take the mean rank

		# we need to look up the rank(s) for the current data value
		# then replace it with the mean rank
		# many times there will only be a single value
		$data->{$k} = mean( @{ $ranks{ $data->{$k} } } );
	}
}

sub normalize_values_by_sum {

	# Scale the data sums to 100

	my $data             = shift;
	my $ref_scale_factor = 100 / sum0( map {abs} values %{$data} );
	map { $data->{$_} = $ref_scale_factor * $data->{$_} } keys %{$data};
}

sub normalize_values_by_median {

	# Scale the data median to 100

	my $data             = shift;
	my $ref_scale_factor = 100 / median( map {abs} values %{$data} );
	map { $data->{$_} = $ref_scale_factor * $data->{$_} } keys %{$data};
}

sub calculate_pearson {

	# positioned data hashes, start, and stop coordinates for evaluating
	my ( $ref, $test ) = @_;

	# calculate correlation
	# 	my $stat = Statistics::LineFit->new();
	# 	$stat->setData($ref, $test) or warn " bad data!\n";
	# 	return $stat->rSquared();
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data($ref);
	my ( $q, $m, $r, $rms ) = $stat->least_squares_fit( @{$test} );
	return $r;
}

sub calculate_ttest {

	# generate the data arrays
	my ( $ref, $test ) = @_;
	my @ref_data;
	my @test_data;

	# load the data arrays
	my $min_value = min( min( keys %{$ref} ), min( keys %{$test} ) );
	my $max_value = max( max( keys %{$ref} ), max( keys %{$test} ) );
	for my $i ( $min_value .. $max_value ) {
		push @ref_data,  $ref->{$i}  || 0;
		push @test_data, $test->{$i} || 0;
	}

	# calculate ANOVA
	my $p_value;
	my $aov = Statistics::ANOVA->new();
	eval {
		$aov->load(
			{
				'test' => \@test_data,
				'ref'  => \@ref_data,
			}
		);
		$aov->anova(
			'independent' => $INDEPENDENCE,
			'parametric'  => $PARAMETRIC,
			'ordinal'     => $ORDINAL,
		);
		$p_value = $aov->{_stat}{p_value};
	};

	$p_value = -1 * ( log($p_value) / LOG10 ) if $p_value;
	return $p_value || '.';
}

sub summarize_results {

	# the results from the correlation analysis
	my ( $correlations, $optimal_shifts, $optimal_correlations,
		$count, $not_enough_data, $no_variance )
		= @_;

	# Summary analyses
	printf " Correlated %s features\n", format_with_commas($count);
	printf " Mean Pearson correlation was %.3f  +/- %.3f\n",
		mean( @{$correlations} ), stddevp( @{$correlations} );

	if ($find_shift) {
		printf " Mean absolute optimal shift was %.0f +/- %.0f bp\n",
			mean( map {abs} @{$optimal_shifts} ),
			stddevp( map {abs} @{$optimal_shifts} );
		printf " Mean optimal Pearson correlation was %.3f  +/- %.3f\n",
			mean( @{$optimal_correlations} ), stddevp( @{$optimal_correlations} );
	}
	if ($not_enough_data) {
		printf " %s features did not have enough data points\n",
			format_with_commas($not_enough_data);
	}
	if ($no_variance) {
		printf " %s features had no variance in the reference data points\n",
			format_with_commas($no_variance);
	}
}

sub mean {
	return sum0(@_) / ( scalar @_ || 1 );
}

__END__

=head1 NAME

correlate_position_data.pl

A script to calculate correlations between two datasets along the length of a feature.

=head1 SYNOPSIS

correlate_position_data.pl [--options] <filename>
  
  Options for data files:
  -i --in <filename>               input file: txt bed etc
  -o --out <filename>              optional output file, default overwrite 
  -d --db <name>                   alternate annotation database
  
  Options for data sources
  -D --ddb <name|file>             data or BigWigSet database
  -r --ref <dataset|filename>      reference data: bw, name, etc
  -t --test <dataset|filename>     test data: bw, name, etc
  
  Options for correlating data
  --pval                           calculate P-value by ANOVA
  --shift                          determine optimal shift to match datasets
  --radius <integer>               radius in bp around reference point to calculate
  -p --pos [5|m|3]                 reference point to measure correlation (m)
  --norm [rank|sum]                normalization method between datasets
  --force_strand                   force an alternate strand
  
  General options:
  -c --cpu <interger>              number of threads (4)
  -z --gz                          compress output with gz
  -v --version                     print version and exit
  -h --help                        show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 Options for data files

=over 4

=item --in <filename>

Specify the input file of features. The features may be comprised of 
name and type, or chromosome, start, and stop. Strand information is 
optional, but is respected if included. A feature list may be 
generated with the program L<get_features.pl>. The file may be 
compressed with gzip.

=item --out <filename>

Specify the output filename. By default it rewrites the input file.

=item --db <name | filename>

Specify the name of a L<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A database is 
required for generating new data files with features. This option may 
skipped when using coordinate information from an input file (e.g. BED 
file), or when using an existing input file with the database indicated 
in the metadata. 

=back

=head2 Options for data sources

=over 4

=item --ddb <name | filename>

If the data to be collected is from a second database that is separate 
from the annotation database, provide the name of the data database here. 
Typically, a second L<Bio::DB::SeqFeature::Store> or BigWigSet database 
is provided here. 

=item --ref <type | filename>

=item --test <type | filename>

Define both the reference and test datasets with which to compare and 
correlate. These may be GFF type or name in a database or BigWigSet, or 
they may be a BigWig or even Bam file. Both options are required. If 
not provided, they may be interactively chosen from the database.

=back

=head2 Options for correlating data

=over 4

=item --pval

Perform an ANOVA analysis between the test and reference data sets and 
report a P-value. By default, this performs a dependent, parametric 
ANOVA. This requires the L<Statistic::ANOVA> module to be installed. 
Please refer to the module documentation for details. If your needs 
require a change to the test, you may edit the parameters at the top 
of this script. For convenience, the P-values are reported as -Log10(P) 
transformed values. The default is false.

=item --shift

Optionally specify whether an optimal shift should be calculated that 
would result in a better Pearson correlation value. The default is 
false.

=item --radius <integer>

Define the radius in basepairs around a reference point to determine 
the window size for the correlation analysis. This value is required 
when calculating an optimal shift (--shift option). The default is to 
take the length of the feature as the window for calculating the 
correlation.  

=item --pos [5|m|3]

Indicate the relative position of the feature to be used as the 
reference point around which the window (determined by the radius 
value) for collecting data will be centered. Three values are 
accepted: "5" indicates the 5' prime end is used, "3" indicates the 
3' end is used, and "m" indicates the middle of the feature is used. 
The default is to use the midpoint. 

=item --norm [rank|sum]

Optionally define a method of normalizing the scores between the 
reference and test data sets prior to calculating the correlation. 
Two methods are currently supported: "rank" converts all values 
to rank values (the mean rank is reported for identical values) 
and essentially calculating a Spearman's rank correlation, while 
"sum" scales all values so that the absolute sums are identical. 
Normalization occurs after missing or zero values are interpolated. 
The default is no normalization.

=item --force_strand

If enabled, a strand orientation will be enforced when determining the 
optimal shift. This does not affect the correlation calculation, only 
the direction of the reported shift. This requires the presence of a 
data column in the input file with strand information. The default is 
no enforcement of strand.

=back

=head2 General options

=over 4

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will calculate statistics between the positioned scores of
two different datasets over a window from an annotated feature or
chromosomal segment. These statistics will help determine whether the
positions or distribution of scores across the window vary or underwent
a positional shift between a test and a reference dataset. For example,
if the enrichment of nucleosome signal from a ChIP experiment shifts in
genomic position, indicating a change in nucleosome position. 

Two statistics may be calculated. First, it will calculate a a Pearson
linear correlation coefficient (r value) between the datasets (default). 
Additionally, an ANOVA analysis may be performed between the datasets and 
generate a P-value. 

By default, the correlation is determined between the data points 
collected over the entire length of the feature. Alternatively, a 
radius and reference point (default is midpoint) may be provided 
that sets the window for collecting scores and calculating a correlation.

In general, to ensure a more reliable Pearson value, fragment ChIP or 
nucleosome coverage should be used rather than point (start or midpoint) 
data, as it will give more reliable results. Fragment coverage is more 
akin to smoothened data and gives better results than interpolated point 
data. 

Normalized read-depth data should be used when possible. If necessary, 
Values can be normalized using one of two methods. The values may be 
converted to rank positions (compare to Kendall's tau), or scaled such 
that the absolute sum values are equal (for example, when working with 
sequence tag read counts).

In addition to calculating a correlation coefficient, an optimal shift 
may also be calculated. This essentially shifts the data, 1 bp at a time, 
in order to identify a shift that would produce a higher correlation. In 
other words, what amount of movement to the left or right would make the 
test data look like the reference data? The window is shifted from -2 
radius to +2 radius relative to the reference point, and the highest 
correlation is reported along with the shift value that generated it. 

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
