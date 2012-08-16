#!/usr/bin/perl

# This script will correlate the positions of occupancy from two datasets
# for a region (nucleosome)

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(sum min max mean stddev);
use Statistics::LineFit;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	find_column_index
	format_with_commas
);
use tim_db_helper qw(
	open_db_connection
	process_and_verify_dataset
	get_region_dataset_hash
);
use tim_file_helper qw(
	load_tim_data_file 
	write_tim_data_file 
);
my $VERSION = '1.9.0';

print "\n This program will correlate positions of occupancy between two datasets\n\n";

### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Get command line options and initialize values
my (
	$infile,
	$outfile,
	$database,
	$data_database,
	$refDataSet,
	$testDataSet,
	$refer,
	$radius,
	$position,
	$set_strand,
	$interpolate,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'        => \$infile, # the input data file
	'out=s'       => \$outfile, # name of output file 
	'db=s'        => \$database, # name of database
	'ddb=s'       => \$data_database, # database containing datasets
	'ref=s'       => \$refDataSet, # reference dataset
	'test=s'      => \$testDataSet, # test dataset
	'radius=i'    => \$radius, # for collecting data
	'pos=s'       => \$position, # set the relative feature position
	'set_strand'  => \$set_strand, # enforce a specific orientation
	'interpolate!' => \$interpolate, # interpolate the position data
	'gz!'         => \$gz, # compress output
	'help'        => \$help, # request help
	'version'     => \$print_version, # print the version
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
	print " Biotoolbox script correlate_position_data.pl, version $VERSION\n\n";
	exit;
}



### Check for requirements
unless ($infile) {
	$infile = shift @ARGV or
		die " no input file! use --help for more information\n";
}
unless (defined $gz) {
	$gz = 0;
}

unless ($radius) {
	$radius = 30;
}

if (defined $position) {
	# check the position value
	unless (
		$position == 5 or
		$position == 3 or
		$position == 4 or
		$position eq 'm'
	) {
		die " Unknown relative position '$position'!\n";
	}
	if ($position eq 'm') {$position = 4} # change to match internal usage
}
else {
	# default position to use the midpoint
	$position = 4;
}

unless (defined $interpolate) {
	$interpolate = 1;
}



### Open input file
my $mainData = load_tim_data_file($infile) or 
	die " unable to open input file '$infile'!\n";



### Open database connection
unless ($database) {
	$database = $mainData->{'db'} or 
		die " no database provided! use --help for more information\n";
}
my $db = open_db_connection($database) or 
	die " unable to open database '$database' connection!\n";
my $ddb;
if ($data_database) {
	# separate database for the refDataSet and testDataSet
	$ddb = open_db_connection($data_database) or 
		die " unable to open data database '$data_database' connection!\n";
}
else {
	# use the same database
	$ddb = $db;
}



### Validate input datasets
validate_or_request_dataset();



### Collect correlations
print " Collecting correlations....\n";
collect_correlations();



### Output file
unless ($outfile) {
	$outfile = $infile;
}
my $success = write_tim_data_file( {
	'data'     => $mainData,
	'filename' => $outfile,
	'gz'       => $gz,
} );
if ($success) {
	print " wrote output file $success\n";
}
else {
	print " unable to write output file!\n";
}
# The End





########################   Subroutines   ###################################

sub validate_or_request_dataset {
	
	# Process the reference dataset
	$refDataSet = process_and_verify_dataset( {
		'db'      => $ddb,
		'dataset' => $refDataSet,
		'prompt'  => "Enter the reference data set  ",
		'single'  => 1,
	} );
	unless ($refDataSet) {
		die " A valid reference data set must be provided!\n";
	}
	
	
	# Process the test dataset
	$testDataSet = process_and_verify_dataset( {
		'db'      => $ddb,
		'dataset' => $testDataSet,
		'prompt'  => "Enter the test data set  ",
		'single'  => 1,
	} );
	unless ($testDataSet) {
		die " A valid test data set must be provided!\n";
	}
}



sub collect_correlations {
	
	# Identify columns
	my $name_i   = find_column_index($mainData, "^name|id");
	my $type_i   = find_column_index($mainData, "^type");
	my $chr_i    = find_column_index($mainData, "^chr|seq");
	my $start_i  = find_column_index($mainData, "^start");
	my $stop_i   = find_column_index($mainData, "^stop|end");
	my $strand_i = find_column_index($mainData, "strand");
	unless (defined $name_i or defined $chr_i) {
		die " unable to identify at least name or chromosome column index!\n";
	}
	if ($set_strand and not $strand_i) {
		warn " unable to identify strand column to enforce strand!\n";
		$set_strand = 0;
	}
	
	# Prepare new columns
	# the rSquared Pearson column
	my $r2_i = $mainData->{'number_columns'};
	$mainData->{$r2_i} = {
		'index'     => $r2_i,
		'name'      => 'Pearson_correlation',
		'reference' => $refDataSet,
		'test'      => $testDataSet,
	};
	$mainData->{'data_table'}->[0][$r2_i] = 'Pearson_correlation';
	$mainData->{'number_columns'} += 1;
	
	# optimal test-reference shift value
	my $shift_i = $mainData->{'number_columns'};
	$mainData->{$shift_i} = {
		'index'     => $shift_i,
		'name'      => 'optimal_shift',
		'reference' => $refDataSet,
		'test'      => $testDataSet,
		'radius'    => $radius,
	};
	$mainData->{'data_table'}->[0][$shift_i] = 'optimal_shift';
	$mainData->{'number_columns'} += 1;
	
	# optimal test-reference shift value
	my $shiftr2_i = $mainData->{'number_columns'};
	$mainData->{$shiftr2_i} = {
		'index'     => $shiftr2_i,
		'name'      => 'shift_correlation',
	};
	$mainData->{'data_table'}->[0][$shiftr2_i] = 'shift_correlation';
	$mainData->{'number_columns'} += 1;
	
	# data database metadata
	if ($data_database) {
		$mainData->{$r2_i}{'data_database'} = $data_database;
		$mainData->{$shift_i}{'data_database'} = $data_database;
	}
	
	# add reference position
	$mainData->{$shift_i}{'reference_point'} = $position == 5 ? 
		'5_prime' : $position == 4 ? 'midpoint' : '3_prime';
	
	
	# Set variables for summary analysis
	my @correlations;
	my @optimal_shifts;
	my @optimal_correlations;
	my $count = 0;
	my $not_enough_data = 0;
	
	
	# Walk through features
	for my $row (1 .. $mainData->{'last_row'}) {
		
		# Determine coordinates
		my ($chromo, $start, $stop, $strand);
		if (defined $name_i and defined $type_i) {
			# using named features
			# retrieve the coordinates from the database
			my @features = $db->features( 
				-name  => $mainData->{'data_table'}->[$row][$name_i],
				-type  => $mainData->{'data_table'}->[$row][$type_i],
			);
			if (scalar @features > 1) {
				# there should only be one feature found
				# if more, there's redundant or duplicated data in the db
				# warn the user, this should be fixed
				warn " Found more than one feature of type " . 
					$mainData->{'data_table'}->[$row][$type_i] . ", name " . 
					$mainData->{'data_table'}->[$row][$name_i] . 
					" in the database!\n Using the first feature only!\n";
			}
			my $feature = shift @features; 
			
			# collect coordinates
			$chromo = $feature->seq_id;
			$start  = $feature->start;
			$stop   = $feature->end;
			if ($set_strand) {
				$strand = $mainData->{'data_table'}->[$row][$strand_i] =~ /-/ ?
					-1 : 1;
			}
			else {
				$strand = $feature->strand;
			}
		}
		else {
			# genomic coordinates defined in the table
			$chromo = $mainData->{'data_table'}->[$row][$chr_i];
			$start  = $mainData->{'data_table'}->[$row][$start_i];
			$stop   = $mainData->{'data_table'}->[$row][$stop_i];
			if (defined $strand_i) {
				$strand = $mainData->{'data_table'}->[$row][$strand_i] =~ /-/ ?
					-1 : 1;
			}
			else {
				$strand = 0;
			}
		}	
		
		# determine reference point
		my $ref_point;
		if ($position == 5) {
			$ref_point = $strand >= 0 ? $start : $stop;
		}
		elsif ($position == 4) {
			$ref_point = int( ( ( $start + $stop ) / 2 ) + 0.5);
		}
		elsif ($position == 3) {
			$ref_point = $strand >= 0 ? $stop : $start;
		}
		
		
		# Collect positioned data
		my %ref_pos2data = get_region_dataset_hash( {
			'db'        => $db,
			'ddb'       => $ddb,
			'dataset'   => $refDataSet,
			'chromo'    => $chromo,
			'start'     => $ref_point - (2 * $radius),
			'stop'      => $ref_point + (2 * $radius),
			'position'  => 4, # reference point is now midpoint based on new coordinates
			'value'     => 'score',
		} );
		my %test_pos2data = get_region_dataset_hash( {
			'db'        => $db,
			'ddb'       => $ddb,
			'dataset'   => $testDataSet,
			'chromo'    => $chromo,
			'start'     => $ref_point - (2 * $radius),
			'stop'      => $ref_point + (2 * $radius),
			'position'  => 4, # reference point is now midpoint based on new coordinates
			'value'     => 'score',
		} );
		
		
		# verify
		if (
			scalar(keys %ref_pos2data) < 5 or 
			scalar(keys %test_pos2data) < 5 or
			sum( map {abs $_} values %ref_pos2data ) == 0 or 
			sum( map {abs $_} values %test_pos2data ) == 0
		) {
			# not enough data points to work with
			# warn "  feature at data row $row returned < 10 data points, skipping\n";
			$mainData->{'data_table'}->[$row][$r2_i]      = '.';
			$mainData->{'data_table'}->[$row][$shift_i]   = '.';
			$mainData->{'data_table'}->[$row][$shiftr2_i] = '.';
			$not_enough_data++;
			next;
		}
			
		
		# Scale the data sums to 100
		# use absolute value just in case there are negatives....
		# I suppose squared values would also work
		# scaling reference data
		my $ref_scale_factor = 100 / sum( map {abs $_} values %ref_pos2data );
		map { 
			$ref_pos2data{$_} = $ref_scale_factor * $ref_pos2data{$_} 
		} keys %ref_pos2data;
		
		# scaling test data
		my $test_scale_factor = 100 / 
			sum( map {abs $_} values %test_pos2data );
		map { 
			$test_pos2data{$_} = $test_scale_factor * $test_pos2data{$_} 
		} keys %test_pos2data;
		
		
		# Interpolate data points
		if ($interpolate) {
			interpolate_values(\%ref_pos2data, 2 * $radius);
			interpolate_values(\%test_pos2data, 2 * $radius);
		}
		
		
		# Calculate the Pearson correlation between test and reference
		# collect data +/- 1 radius from midpoint
		my @ref_data;
		my @test_data;
		for my $i ( (0 - $radius) .. $radius ) {
			# get values for each position, default 0
			push @ref_data, $ref_pos2data{$i} || 0;
			push @test_data, $test_pos2data{$i} || 0;
		}
		my $r2 = calculate_pearson(\@ref_data, \@test_data);
		
		
		# Determine optimal shift
		# calculate a Pearson r^2 correlation for each shift
		# calculating window of 2*radius, shifting 2*radius across region
		my %shift2pearson;
		for my $adjustment ( (-1 * $radius) .. $radius ) {
			# adjustment is from -1radius to +1radius
			# starting point is at -1radius to +1radius
			# so adding together we get first window from -2radius to 0
			# and last window is from 0 to +2radius
			# this will be compared to the reference window, which does not move
			
			# generate the array of test values
			my @test_data;
			for my $i ((0 - $radius + $adjustment) .. ($radius + $adjustment)) {
				push @test_data, $test_pos2data{$i} || 0;
			}
			
			# calculate Pearson correlation
			my $pearson = calculate_pearson( \@ref_data, \@test_data); 
			if (defined $pearson) {
				$shift2pearson{$adjustment} = $pearson;
			}
		}
		
		# identify the best shift
		my $best_r2 = 0; 
		my $best_shift = 1000;
		foreach my $shift (keys %shift2pearson) {
			if ($shift2pearson{$shift} >= $best_r2) {
				# a good looking R^2 value
				# but we only want those with the smallest shift
				if (abs($shift) < abs($best_shift)) {
					$best_r2 = $shift2pearson{$shift};
					$best_shift = $shift;
				}
			}
		}
		if ($best_shift == 1000) {
			# no good shift found!? then that means no shift
			$best_shift = 0;
		}
		
		# Record final data
		# change strand
		if ($set_strand) {
			# user is imposing a strand
			if ($strand < 0) {
				# only change shift direction if reverse strand
				$best_shift *= -1;
			}
		}
		
		# record in data table
		$mainData->{'data_table'}->[$row][$r2_i] = defined $r2 ? $r2 : '.';
		$mainData->{'data_table'}->[$row][$shift_i] = $best_shift;
		$mainData->{'data_table'}->[$row][$shiftr2_i] = $best_r2;
		
		# record for summary analyses
		push @correlations, $r2 if defined $r2;
		push @optimal_shifts, $best_shift;
		push @optimal_correlations, $best_r2;
		$count++;
	}
	
	
	# Summary analyses
	printf " Correlated %s features\n", format_with_commas($count);
	printf " Mean Pearson correlation was %.3f  +/- %.3f\n", 
		mean(@correlations), stddev(@correlations);
	printf " Mean absolute optimal shift was %.0f +/- %.0f bp\n", 
		mean( map {abs $_} @optimal_shifts), 
		stddev( map {abs $_} @optimal_shifts);
	printf " Mean optimal Pearson correlation was %.3f  +/- %.3f\n", 
		mean(@optimal_correlations), stddev(@optimal_correlations);
	if ($not_enough_data) {
		printf " %s features did not have enough data points\n", 
			format_with_commas($not_enough_data);
	}
}


sub interpolate_values {
	
	my ($data, $limit) = @_;
	
	# Fill out the data so that we have all elements filled
	for (my $i = 0 - $limit; $i <= $limit; $i++) {
		unless (exists $data->{$i}) {
			# using default value of zero
			$data->{$i} = 0;
		}
	}
	
	# Interpolate data
	my $i = 0 - $limit; # starting point
	while ($i <= $limit) {
		# walk through the data
		# we're not using a for loop here because we may interpolate
		# more than one element at a time
		
		
		# Look for zero values
		if ($data->{$i} == 0) {
			# we have a zero value
			
			# find the next non-zero value to interpolate
			my $next = $i + 1;
			while ($next < $limit) {
				last if ($data->{$next} != 0); # we found it
				$next++;
			}
			
			# check if we reached the limit
			if ($next == $limit) {
				# that's it, nothing left to interpolate
				# return what we have
				return;
			}
			
			# determine fractional number
			# difference between next non-zero number and previous number
			# divided by the number of elements in between
			my $fraction = ($data->{$next} - $data->{$i - 1}) / 
				(abs($next - $i) + 1);
			
			# apply fractions
			for (my $n = $i; $n < $next; $n++) {
				$data->{$n} = $data->{$i - 1} + ($fraction * (abs($n - $i) + 1));
			}
			
			# reset for next round
			$i = $next + 1;
		}
		
		else {
			# a non-zero number, no interpolation necessary, next please
			$i++;
		}
	}
	# done, just return
}




sub calculate_pearson {
	
	# positioned data hashes, start, and stop coordinates for evaluating
	my ($ref, $test) = @_;
	
	# check that we have non-zero values in both arrays
	return unless (
		sum(map {abs $_} @$ref) > 0 and 
		sum(map {abs $_} @$test) > 0
	);
	
	# calculate correlation
	my $stat = Statistics::LineFit->new();
	$stat->setData($ref, $test) or warn " bad data!\n";
	return $stat->rSquared();
}






__END__

=head1 NAME

correlate_position_data.pl

=head1 SYNOPSIS

correlate_position_data.pl [--options] <filename>
  
  Options:
  --in <filename>
  --out <filename> 
  --db <name | filename>
  --ddb <name | filename>
  --ref <type | filename>
  --test <type | filename>
  --radius <integer>
  --pos [5 | m | 3]
  --set_strand
  --(no)interpolate
  --(no)gz
  --version
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input file of features. The features may be comprised of 
name and type, or chromosome, start, and stop. Strand information is 
optional, but is respected if included. A feature list may be 
generated with the program B<get_features.pl>. The file may be 
compressed with gzip.

=item --out <filename>

Specify the output filename. By default it rewrites the input file.

=item --db <name | filename>

Specify the BioPerl database that contains the features listed in the 
input file and potentially the reference and test datasets. A 
Bio::DB::SeqFeature::Store database is typically used for 
storing annotation. If genomic coordinates are provided instead, then 
any database or bigFile may be provided. The default database is 
obtained from the metadata of the input file. 

=item --ddb <name | filename>

Optionally specify an alternative BioPerl database or BigWigSet directory 
from which the reference and test datasets may be obtained. 

=item --ref <type | filename>

=item --test <type | filename>

Define both the reference and test datasets with which to compare and 
correlate. These may be GFF type or name in a database or BigWigSet, or 
they may be a BigWig or even Bam file. Both options are required. If 
not provided, they may be interactively chosen from the database.

=item --radius <integer>

Define the radius in basepairs to determine the window size for the 
correlation analysis and the limit for sliding the window to determine 
the optimal shift. Default is 30 bp. 

=item --pos [5 | m | 3]

Indicate the relative position of the feature to be used as the 
reference point around which the window for collecting data will 
be centered. Three values are accepted: "5" indicates the 
5' prime end is used, "3" indicates the 3' end is used, and "m" 
indicates the middle of the feature is used. The default is to 
use the midpoint. 

=item --set_strand

If enabled, a strand orientation will be enforced when determining the 
optimal shift. This does not affect the correlation calculation, only 
the direction of the reported shift. This requires the presence of a 
data column in the input file with strand information. The default is 
no enforcement of strand.

=item --(no)interpolate

Interpolate missing or zero positioned values in each window for both 
reference and test data. This will improve the Pearson correlation 
values, especially for sparse data. Enabled by default.

=item --(no)gz

Specify whether (or not) the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will calculate a Pearson correlation coefficient (R^2 value) 
between the positioned scores (occupancy) of two datasets over a window of 
an annotated feature or chromosomal segment. This will determine whether 
the positions or distribution of scores across the window vary between 
two different data sets: a test dataset and a reference dataset. 
The original implementation of this program is to compare nucleosome 
occupancy differences between two datasets and identify shifts in position. 

The window is determined by the radius parameter. The window is set as 
+1 to -1 radius from the feature's reference point. The default reference 
point is the feature midpoint.

To ensure a more reliable Pearson value, the absolute sum of the data 
points in the window are scaled to the same value, and missing values or 
values of zero are interpolated from neighboring values, when possible.

To determine whether a significant shift has occurred, the window for the 
test dataset is shifted, 1 bp at a time, from -2 radius to +2 radius from 
the reference point, and new Pearson correlations are determined. The 
highest correlation found is reported along with the shift value that 
generated it.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  

