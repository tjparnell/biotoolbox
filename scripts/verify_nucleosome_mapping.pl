#!/usr/bin/perl

# This script will verify nucleosome mapping
# It will identify overlaps and accuracy

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(min max mean stddev);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	find_column_index
);
use tim_db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	get_region_dataset_hash
);
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
);
my $VERSION = '1.10';

print "\n This program will verify the mapping of nucleosomes\n\n";

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
	$database,
	$dataset,
	$filter,
	$max_overlap,
	$outfile,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the input data file
	'db=s'      => \$database, # the database
	'data=s'    => \$dataset, # the dataset to verify
	'filter!'   => \$filter, # remover overlapping nucleosomes
	'max=i'     => \$max_overlap, # maximum overlap allowed
	'out=s'     => \$outfile, # the output file name
	'gz!'       => \$gz, # compress output
	'help'      => \$help, # request help
	'version'   => \$print_version, # print the version
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
	print " Biotoolbox script verify_nucleosome_mapping.pl, version $VERSION\n\n";
	exit;
}




### Check for requirements
unless ($infile) {
	$infile = shift @ARGV or
		die "  OOPS! No source data file specified! \n use --help\n";
}
unless (defined $max_overlap) {
	$max_overlap = 30;
}
if ($filter and !$outfile) {
	die " You should specify an output file name when filtering nucleosomes\n";
}

### Load the nucleosome file
my $data_ref = load_tim_data_file($infile) or die 
	" can't load file '$infile'!\n";
unless ($data_ref->{feature} eq 'nucleosome') {
	warn " file feature is not 'nucleosome'! proceeding anyway\n";
}
print " Loaded file '$infile' with " . $data_ref->{last_row} . " nucleosomes\n";
my $table = $data_ref->{data_table};



### Open db connection
my $db;
if ($database) {
	# a specific database was defined on the command line
	$db = open_db_connection($database) or 
		die " can't connect to database $database\n";
}
elsif ($data_ref->{'db'}) {
	# a database was defined in the metadata
	$db = open_db_connection($data_ref->{db}) or 
		die " can't connect to database ". $data_ref->{db} . "\n";
}
else {
	# no database was defined
	# use one of the datasets
	if ($dataset) {
		$db = open_db_connection($dataset) or 
			die " can't connect to a database!\n";
	}
	else {
		# look for a dataset defined in the metadata
		my $index = find_column_index($data_ref, 'Occupancy') || 
					find_column_index($data_ref, 'NucleosomeID') ||
					find_column_index($data_ref, 'score') ||
					undef;
		unless (defined $index) {
			die " Unable to identify Occupancy or Score dataset column!\n";
		}
		
		# pull out the dataset name
		$dataset = $data_ref->{$index}{'dataset'} || 
				   $data_ref->{$index}{'scan_dataset'} || 
				   undef; 
		$dataset =~ s/^\w+\/*://; # strip any prefix
		$db = open_db_connection($dataset) or 
			die " can't connect to a database!\n";
	}
}



### Define the dataset
# defined in the Occupancy 
if ($dataset) {
	# dataset defined on the command line
	$dataset = verify_or_request_feature_types( {
		'db'        => $db,
		'feature'   => $dataset,
		'prompt'    => "Enter the dataset to use for verifying nucleosome maps",
		'single'    => 1,
	} );
}	
else {	
	# use a dataset defined in the input metadata
	# usually defined in the Occupancy metadata, or maybe the NucleosomeID
	# or maybe just the score column if nothing else is present
	my $index = find_column_index($data_ref, 'Occupancy') || 
				find_column_index($data_ref, 'NucleosomeID') ||
				find_column_index($data_ref, 'score') ||
				undef;
	unless (defined $index) {
		die " Unable to identify Occupancy or Score dataset column!\n";
	}
	
	# pull out the dataset name
	$dataset = $data_ref->{$index}{'dataset'} || 
			   $data_ref->{$index}{'scan_dataset'} || 
			   undef; 
	$dataset = verify_or_request_feature_types( {
		'db'        => $db,
		'feature'   => $dataset,
		'prompt'    => "Enter the dataset to use for verifying nucleosome maps",
		'single'    => 1,
	} );
}
unless ($dataset) {
	die " unable to identify the dataset in the input file metadata!\n";
}



### Add new columns
my $overlap_i = $data_ref->{'number_columns'};
$data_ref->{$overlap_i} = {
	'name'       => 'overlap_length',
	'index'      => '$overlap_i',
};
if ($filter) {
	$data_ref->{$overlap_i}{'max_overlap'} = $max_overlap;
}
my $mapping_i = $overlap_i + 1;
$data_ref->{$mapping_i} = {
	'name'       => 'center_peak_mapping',
	'index'      => '$mapping_i',
	'dataset'    => $dataset,
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
# initialize variables for statistical analysis at the end
my $overlap_count = 0;
my @overlaps;
my @peak_distances;
my $on_target_count = 0;

# Identify indices
my $chrom_i = find_column_index($data_ref, '^chrom|chr|seq|ref');
my $start_i = find_column_index($data_ref, '^start');
my $stop_i  = find_column_index($data_ref, '^stop|end');
my $mid_i   = find_column_index($data_ref, '^midpoint|mid');
my $score_i = find_column_index($data_ref, '^score|occupancy');
unless (defined $chrom_i and defined $start_i and defined $stop_i) {
	die " Unable to identify Chromosome, Start, and Stop columns!\n";
}

# main loop through data
foreach (my $row = 1; $row <= $data_ref->{last_row}; $row++) {
	
	## Check overlaps
	if ($row < $data_ref->{last_row} and 
		$table->[$row][$chrom_i] eq $table->[$row + 1][$chrom_i]
	) {
		
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
	# determine midpoint
	my $midpoint = defined $mid_i ? $table->[$row][$mid_i] : 
		int( ( ($table->[$row][$start_i] + $table->[$row][$stop_i]) / 2) + 0.5 );
	
	# collect the raw data for 71 bp around the mapped midpoint
	# this distance is arbitrary but about 1/2 a nucleosome
	my %nuc_data = get_region_dataset_hash( {
		'db'         => $db,
		'dataset'    => $dataset,
		'chromo'     => $table->[$row][$chrom_i],
		'start'      => $midpoint - 35,
		'stop'       => $midpoint + 35,
		'value'      => 'score',
		'absolute'   => 1,
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
			push @distances, $pos - $midpoint;
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
	}
	elsif (scalar @distances > 1) {
		# more than one peak found
		
		# we'll use the closest one
		# sort the distances by their absolute value
		@distances = sort { abs($a) <=> abs($b) } @distances;
		if (abs($distances[0]) <= 10) {
			# our tolerance is an arbitrary 10 bp to consider on target
			$on_target_count++;
			$table->[$row][$mapping_i] = 'centered';
		}
		else {
			push @peak_distances, $distances[0];
			$table->[$row][$mapping_i] = 'offset';
		}
	}
	
	# record the distance
	$table->[$row][$offset_i] = $distances[0];
}


### digest results
printf " There were $overlap_count overlapping nucleosomes out of %s (%.0f%%)\n" . 
	"   the overlap mean was %.0f +/- %.0f bp\n", 
	$data_ref->{last_row}, ($overlap_count / $data_ref->{last_row}) * 100, 
	mean(@overlaps), stddev(@overlaps);


printf "\n There were $on_target_count (%.0f%%) accurately mapped nucleosomes\n" . 
	" and %s offcenter nucleosomes\n" .
	"   the peak offset distance mean was %.0f +/- %.0f bp\n", 
	($on_target_count / $data_ref->{last_row}) * 100, 
	$data_ref->{last_row} - $on_target_count, 
	mean(@peak_distances), stddev(@peak_distances);



### Filter out the overlapping nucleosomes if requested
if ($filter and !defined $score_i) {
	warn " unable to identify score or Occupancy column, cannot filter\n";
}
elsif ($filter) {
	
	# we will identify those that need to be deleted and put them in this array
	my @to_delete;
	my $offset_count = 0;
	my $occupancy_count = 0;
	
	# walk through the list again
	foreach (my $row = 1; $row <= $data_ref->{last_row}; $row++) {
		
		# check for overlap
		next if $table->[$row][$overlap_i] eq '.';
		if ($table->[$row][$overlap_i] > $max_overlap) {
			# overlap with the next nucleosome is too much
			
			# decide which one to delete
			if (
				$table->[$row][$mapping_i] eq 'offset' and 
				$table->[$row + 1][$mapping_i] eq 'offset'
			) {
				# both nucleosomes are offset
				# delete whichever one is offset more
				if ( $table->[$row][$offset_i] > $table->[$row+1][$offset_i] ) {
					# current one is more offset, delete this one
					push @to_delete, $row;
				}
				else {
					# next one is more offset or at least equal
					# delete that one
					push @to_delete, $row+1;
					
					# advance ahead
					$row++;
				}
				$offset_count++;
			}
			elsif ($table->[$row][$mapping_i] eq 'offset') {
				# current nucleosome is offset, so delete that one
				push @to_delete, $row;
				$offset_count++;
			}
			elsif ($table->[$row + 1][$mapping_i] eq 'offset') {
				# next nucleosome is offset, so delete that one
				push @to_delete, $row + 1;
				
				# advance ahead
				$row++;
				$offset_count++;
			}
			else {
				# both nucleosomes are centered, so take whichever one is more occupied
				if ( $table->[$row][$score_i] < $table->[$row+1][$score_i] ) {
					# current one is less occupied, delete this one
					push @to_delete, $row;
				}
				else {
					# next one is less occupied or at least equal
					# delete that one
					push @to_delete, $row+1;
					
					# advance ahead
					$row++;
				}
				$occupancy_count++;
			}
		}
	}
	
	
	# proceed to delete nucleosomes
	my $deleted = 0;
	while (@to_delete) {
		# take from the end, highest rows first
		my $row = pop @to_delete;
		splice( @{$table}, $row, 1);
		$deleted++;
		$data_ref->{last_row}--;
	}
	
	# report
	print " $deleted nucleosomes were filtered out due to extensive overlap\n";
	print "   $offset_count were tossed because they were offset\n";
	print "   $occupancy_count were tossed because of low occupancy\n";
}


### Write results
# check filename
unless ($outfile) {
	$outfile = $infile;
}
my $success;
if (
	$success = write_tim_data_file( {
		'data'     => $data_ref,
		'filename' => $outfile,
		'gz'       => $gz,
	}) 
) {
	print " wrote data file $success\n";
}
else {
	print " unable to write file!\n";
}




__END__

=head1 NAME

verify_nucleosome_mapping.pl

=head1 SYNOPSIS

verify_nucleosome_mapping.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --db <text>
  --data <text | filename>
  --filter
  --max <integer>
  --out <filename> 
  --(no)gz
  --version
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input file. The program expects a text data file with 
at least chromosome, start, and stop coordinates representing 
nucleosome annotation. The output file from the biotoolbox script 
'map_nucleosomes.pl' is ideal, although a Bed, GFF, or other file 
is allowed. The file may be compressed with gzip.

=item --db <text>

Specify the name of a BioPerl database to pull the source data. A 
SeqFeature::Store database or BigWigSet directory may be supplied. 
The default is to use the database defined in the input file metadata 
or the dataset file. 

=item --data <text | filename>

Provide the name of the dataset or data file (bigWig format)
containing the nucleosome midpoint occupancy data with which to
verify nucleosomal positions. If data is obtained from a database,
the type or primary_tag should be provided. Default is to use the 
dataset defined in the input file metadata.

=item --filter

Optionally filter out nucleosomes that exceed a maximum allowed 
overlap. Filtered nucleosomes are deleted from the file. Default is 
no filtering.

=item --max <integer>

Specify the maximum allowed overlap in bp when filtering out 
overlapping nucleosomes. The default is 30 bp.

=item --out <filename>

Optionally specify the output file name. By default it will overwrite 
the input file.

=item --gz

Optionally indicate the file should be compressed when written.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will verify the mapping of nucleosomes. It expects as
input a text file containing the genomic coordinates of annotated
nucleosomes. The output data file from the BioToolBox script
'map_nucleosomes.pl' is ideal, although other text files, including
Bed and GFF, are supported. Nucleosomes are verified by comparing the
peak of nucleosome reads collected from the original dataset with the
recorded midpoint of the mapped nucleosome. If the peak is <= 10 bp
from the recorded midpoint, then the nucleosome is considered centered
on the peak and it is properly mapped. If the peak is > 10 bp from the
recorded midpoint, then it is considered offset, or improperly mapped.
This is most often due to nucleosomes being called prematurely when
the dataset is being scanned in windows from left to right. Adjusting
the parameters for window and buffer in map_nucleosomes.pl can limit
the number of overlapping nucleosomes.

The same data file is re-written or a new file written with three 
additional columns appended, overlap_length, center_peak_mapping, and 
center_peak_offset. 

Overlap_length records the amount of overlap between mapped nucleosomes. 
Ideally this should be 0 as nucleosomes should not overlap; overlapping 
nucleosomes indicates either an error in mapping or multiple phasing of 
nucleosomes. 

Center_peak_mapping records whether the nucleosome was properly mapped or 
not. One of two values is recorded: centered or offset.

Center_peak_offset records the distance in bp between the nucleosome peak 
and the recorded midpoint.

Basic statistics are reported to Standard Output for the overlap and 
offset lengths.

The program can also optionally filter out nucleosomes which exceed a 
set limit of overlap. The overlapping nucleosome to be deleted is 
chosen based on a set of rules: offset nucleosomes are deleted, or if 
both nucleosomes are offset, then the one with the greatest offset is 
deleted. If neither nucleosome is offset, then the nucleosome with the 
lowest occupancy is deleted, or if both occupancies are equal, then 
the rightmost is deleted.

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
