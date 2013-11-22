#!/usr/bin/perl

# This script will verify nucleosome mapping
# It will identify overlaps and accuracy

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(min max statsinfo);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_db_helper qw(
	open_db_connection
	get_region_dataset_hash
);
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
);

print "\n This program will verify the mapping of nucleosomes\n";

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
	$help
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the solexa data file
	'help'      => \$help # request help
);

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}



### Check for requirements
unless ($infile) {
	$infile = shift @ARGV or
		die "  OOPS! No source data file specified! \n use --help\n";
}



### Load the nucleosome file
my $data_ref = load_tim_data_file($infile) or die 
	" can't load file '$infile'!\n";
unless ($data_ref->{feature} eq 'nucleosome') {
	die " file feature is not 'nucleosome'! Are you loading the right file?\n";
}
print " Loaded file '$infile' with " . $data_ref->{last_row} . " nucleosomes\n";
my $table = $data_ref->{data_table};


### Open db connection
my $db = open_db_connection($data_ref->{db}) or 
	die " can't connect to database ". $data_ref->{db} . "\n";
my $dataset = $data_ref->{5}{dataset}; # Occupancy column



### Add new columns
my $overlap_i = $data_ref->{'number_columns'};
$data_ref->{$overlap_i} = {
	'name'       => 'overlap_length',
	'index'      => '$overlap_i',
};
my $mapping_i = $overlap_i + 1;
$data_ref->{$mapping_i} = {
	'name'       => 'center_peak_mapping',
	'index'      => '$mapping_i',
};
my $offset_i = $overlap_i + 2;
$data_ref->{$offset_i} = {
	'name'       => 'center_peak_offset',
	'index'      => '$offset_i',
};
$data_ref->{'number_columns'} += 3;
$table->[0][$overlap_i] = 'overlap_length';
$table->[0][$mapping_i] = 'center_peak_mapping';
$table->[0][$offset_i]  = 'center_peak_offset';


### Process the nucleosomes
# initialize variables
my $overlap_count = 0;
my @overlaps;
my @peak_distances;
my $on_target_count = 0;

# identify indices
my $chrom_i = 0;
my $start_i = 1;
my $stop_i = 2;
my $mid_i = 3;

# main loop through data
foreach (my $row = 1; $row <= $data_ref->{last_row}; $row++) {
	
	## Check overlaps
	if ($row < $data_ref->{last_row} and 
		$table->[$row][$chrom_i] eq $table->[$row + 1][$chrom_i]) {
		
		if (
			# next start less than current stop
			$table->[$row + 1][$start_i] < $table->[$row][$stop_i] and 
			# current start greater than previous start (should be)
			$table->[$row + 1][$start_i] > $table->[$row][$start_i]
		) {
			# there is overlap
			my $overlap =  # calculate
				$table->[$row][$stop_i] - $table->[$row + 1][$start_i];
			
			# record
			$table->[$row][$overlap_i] = $overlap;
			push @overlaps, $overlap;
			$overlap_count++;
		}
		else {
			# no overlap
			$table->[$row][$overlap_i] = '.';
		}
	}
	
	## Check midpoint
	# collect the raw data for 101 bp around the mapped midpoint
	# this distance is arbitrary
	my %nuc_data = get_region_dataset_hash( {
		'db'         => $db,
		'dataset'    => $dataset,
		'name'       => $table->[$row][$chrom_i],
		'type'       => 'chromosome',          # this shouldn't be hardcoded
		'start'      => $table->[$row][$mid_i] - 50,
		'stop'       => $table->[$row][$mid_i] + 50,
	});
	
	# find the max peak
	my @distances;
	my $max = max(values %nuc_data);
	foreach my $pos (keys %nuc_data) {
		if ($nuc_data{$pos} == $max) {
			# one of the peaks, likely only one but may be more
			# store the distance between this peak and recorded midpoint
			# these distances are usually positive because of the 
			# way they are mapped - a premature call of a nucleosome
			push @distances, $pos - $table->[$row][$mid_i];
		}
	}
	# record the minimum distance from the recorded midpoint and occupancy peak
	if (scalar @distances == 1) {
		# only one peak found
		if (abs($distances[0]) <= 10) {
			# our tolerance is an arbitrary 10 bp to consider on target
			$on_target_count++;
			$table->[$row][$mapping_i] = 'centered';
		}
		else {
			push @peak_distances, abs($distances[0]);
			$table->[$row][$mapping_i] = 'offset';
		}
		
		# record the distance
		$table->[$row][$offset_i] = $distances[0];
	}
	elsif (scalar @distances > 1) {
		# more than one peak found
		
		# we'll use the closest one
		my $distance = min(map { abs($_) } @distances);
		if ($distance <= 10) {
			# our tolerance is an arbitrary 10 bp to consider on target
			$on_target_count++;
			$table->[$row][$mapping_i] = 'centered';
		}
		else {
			push @peak_distances, $distance;
			$table->[$row][$mapping_i] = 'offset';
		}
		
		# record the distance
		$table->[$row][$offset_i] = $distances[0];
	}
	
}

### digest results

my $percent = sprintf "%.0f%%", ($overlap_count / $data_ref->{last_row}) * 100;

print "  There are $overlap_count overlapping fragments out of a total of " . 
	"$data_ref->{last_row} ($percent)\n";

print "  overlaps:\n" . statsinfo(@overlaps) . "\n";


$percent = sprintf "%.0f%%", ($on_target_count / $data_ref->{last_row}) * 100;
print "\n There were $on_target_count ($percent) perfectly mapped nucleosomes\n";
print "  of the remaining " . ($data_ref->{last_row} - $on_target_count) .
	" nucleosomes, the offset peak distances stats are:\n" . 
	statsinfo(@peak_distances) . "\n";


### Write results
if (
	write_tim_data_file( {
		'data' => $data_ref,
	}) 
) {
	print " wrote data file\n";
}
else {
	print " unable to write file!\n";
}




__END__

=head1 NAME

<name>.pl

=head1 SYNOPSIS

<name>.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --out <filename> 
  --(no)gz
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input file. The program expects the output data file from 
the biotoolbox script 'map_nucleosomes.pl'. The file may be compressed 
with gzip.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will verify the mapping of nucleosomes. It expects as input 
the output data file from 'map_nucleosomes.pl'. Nucleosomes are verified 
by comparing the peak of nucleosome reads collected from the original 
dataset with the recorded midpoint of the mapped nucleosome. If the peak 
is <= 10 bp from the recorded midpoint, then the nucleosome is considered 
centered on the peak and it is properly mapped. If the peak is > 10 bp 
from the recorded midpoint, then it is considered offset, or improperly 
mapped. This is most often due to nucleosomes being called prematurely 
when the dataset is being scanned in windows from left to right. 
Adjusting the parameters for window and buffer in map_nucleosomes.pl can 
limit the number of overlapping nucleosomes.

The same data file is re-written with three additional columns appended, 
overlap_length, center_peak_mapping, and center_peak_offset. 

Overlap_length records the amount of overlap between mapped nucleosomes. 
Ideally this should be 0 as nucleosomes should not overlap; overlapping 
nucleosomes indicate an error in mapping or multiple phasing of nucleosomes. 

Center_peak_mapping records whether the nucleosome was properly mapped or 
not; one of two values is recorded: centered or offset.

Center_peak_offset records the distance in bp between the nucleosome peak 
and the recorded midpoint.



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

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
