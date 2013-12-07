#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(min max mean stddev);
use Bio::ToolBox::data_helper qw(
	find_column_index
	format_with_commas
);
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	get_region_dataset_hash
);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	write_tim_data_file
);
my $VERSION = '1.14';

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
	$recenter,
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
	'recenter!' => \$recenter, # correct off center nucleosomes
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
	$dataset = verify_or_request_feature_types(
		'db'        => $db,
		'feature'   => $dataset,
		'prompt'    => "Enter the dataset to use for verifying nucleosome maps",
		'single'    => 1,
	);
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
	$dataset = verify_or_request_feature_types(
		'db'        => $db,
		'feature'   => $dataset,
		'prompt'    => "Enter the dataset to use for verifying nucleosome maps",
		'single'    => 1,
	);
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
my $recenter_count  = 0;

# Identify indices
my $chrom_i = find_column_index($data_ref, '^chrom|chr|seq|ref');
my $start_i = find_column_index($data_ref, '^start');
my $stop_i  = find_column_index($data_ref, '^stop|end');
my $mid_i   = find_column_index($data_ref, '^midpoint|mid');
my $score_i = find_column_index($data_ref, '^score|occupancy');
my $name_i  = find_column_index($data_ref, '^NucleosomeID|name|ID');
unless (defined $chrom_i and defined $start_i and defined $stop_i) {
	die " Unable to identify Chromosome, Start, and Stop columns!\n";
}

# main loop through data
foreach (my $row = 1; $row <= $data_ref->{last_row}; $row++) {
	
	## Check midpoint
	# determine midpoint
	my $midpoint = defined $mid_i ? $table->[$row][$mid_i] : 
		int( ( ($table->[$row][$start_i] + $table->[$row][$stop_i]) / 2) + 0.5 );
	
	# collect the raw data for 71 bp around the mapped midpoint
	# this distance is arbitrary but about 1/2 a nucleosome
	my %nuc_data = get_region_dataset_hash(
		'db'         => $db,
		'dataset'    => $dataset,
		'chromo'     => $table->[$row][$chrom_i],
		'start'      => $midpoint - 35,
		'stop'       => $midpoint + 35,
		'value'      => 'score',
		'absolute'   => 1,
	);
	
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
		if (abs($distances[0]) <= 1) {
			# our tolerance is an arbitrary 1 bp to consider on target
			$on_target_count++;
			$table->[$row][$mapping_i] = 'centered';
		}
		elsif ($recenter) {
			# an off-center nucleosome and user requested to correct it
			$table->[$row][$start_i] += $distances[0];
			$table->[$row][$stop_i]  += $distances[0];
			$midpoint += $distances[0];
			if (defined $mid_i) {
				$table->[$row][$mid_i] = $midpoint;
			}
			$table->[$row][$name_i] =~ s/\:\d+$/:$midpoint/;
			$table->[$row][$mapping_i] = 'recentered';
			$distances[0] = 0; # reset for recording below
			$recenter_count++;
		}
		else {
			# record off-center nucleosome
			push @peak_distances, abs($distances[0]);
			$table->[$row][$mapping_i] = 'offset';
		}
	}
	elsif (scalar @distances > 1) {
		# more than one peak found
		
		# we'll use the closest one
		# sort the distances by their absolute value
		@distances = sort { abs($a) <=> abs($b) } @distances;
		if (abs($distances[0]) <= 1) {
			# our tolerance is an arbitrary 1 bp to consider on target
			$on_target_count++;
			$table->[$row][$mapping_i] = 'centered';
		}
		elsif ($recenter) {
			# an off-center nucleosome and user requested to correct it
			$table->[$row][$start_i] += $distances[0];
			$table->[$row][$stop_i]  += $distances[0];
			$midpoint += $distances[0];
			if (defined $mid_i) {
				$table->[$row][$mid_i] = $midpoint;
			}
			$table->[$row][$name_i] =~ s/\:\d+$/:$midpoint/;
			$table->[$row][$mapping_i] = 'recentered';
			$distances[0] = 0; # reset for recording below
			$recenter_count++;
		}
		else {
			# record off-center nucleosome
			push @peak_distances, $distances[0];
			$table->[$row][$mapping_i] = 'offset';
		}
	}
	# record the distance
	$table->[$row][$offset_i] = $distances[0];
	
	# recalculate overlap of previous nucleosome if recentered
	if ($recenter and $table->[$row][$mapping_i] eq 'recentered') {
		# the previous nucleosome would've calculated an overlap value 
		# based on the old position
		# now that we have re-centered the current nucleosome, we need 
		# to go back and re-calculate the overlap
		
		# re-calculate
		if ($table->[$row - 1][$chrom_i] eq $table->[$row][$chrom_i]
		) {
		
			# first adjust the count if necessary
			if ($table->[$row - 1][$overlap_i] =~ /\d+/) {
				$overlap_count--;
				pop @overlaps;
			}
		
			if ( $table->[$row - 1][$stop_i] > $table->[$row][$start_i] ) {
				# previous stop greater than current start
				# there is overlap
				my $overlap =  # calculate
					$table->[$row - 1][$stop_i] - $table->[$row][$start_i];
			
				# record
				$table->[$row - 1][$overlap_i] = $overlap;
				push @overlaps, $overlap;
				$overlap_count++;
			}
			else {
				# no overlap
				$table->[$row - 1][$overlap_i] = '.';
			}
		}
	}
	
	
	## Check overlaps
	if ($row < $data_ref->{last_row} and 
		$table->[$row][$chrom_i] eq $table->[$row + 1][$chrom_i]
	) {
		
		if ( $table->[$row][$stop_i] > $table->[$row + 1][$start_i] ) {
			# current stop greater than next start
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
	
}


### digest results
printf " %s (%.0f%%) nucleosomes were overlapping\n" . 
	"   the overlap mean was %.0f +/- %.0f bp\n", 
	format_with_commas($overlap_count), ($overlap_count / $data_ref->{last_row}) * 100, 
	mean(@overlaps), stddev(@overlaps);

if ($recenter) {
	printf "\n %s (%.0f%%) nucleosomes were accurately mapped\n" . 
		" %s nucleosomes were re-centered on the data peak\n",
		format_with_commas($on_target_count),
		($on_target_count / $data_ref->{last_row}) * 100,
		format_with_commas($recenter_count);
}
else { 
	printf "\n %s (%.0f%%) nucleosomes were accurately mapped\n" . 
		" %s (%.0f%%) nucleosomes were offcenter\n" .
		"   the offset distance mean was %.0f +/- %.0f bp\n", 
		format_with_commas($on_target_count),
		($on_target_count / $data_ref->{last_row}) * 100, 
		format_with_commas($data_ref->{last_row} - $on_target_count), 
		( ($data_ref->{last_row} - $on_target_count) / $data_ref->{last_row}) * 100,
		mean(@peak_distances), stddev(@peak_distances);
}



### Filter out the overlapping nucleosomes if requested
if ($filter and !defined $score_i) {
	warn " unable to identify score or Occupancy column, cannot filter\n";
}
elsif ($filter) {
	
	# we will identify those that need to be deleted and put them in this array
	my @to_delete;
	my $offset_count    = 0;
	my $occupancy_count = 0;
	my $dual_count      = 0;
	
	# walk through the list again
	foreach (my $row = 1; $row <= $data_ref->{last_row}; $row++) {
		
		# check for overlap
		next if $table->[$row][$overlap_i] eq '.';
		if ($table->[$row][$overlap_i] > $max_overlap) {
			# overlap with the next nucleosome is too much
			
			# decide which one to delete
			if ($table->[$row + 1][$overlap_i] > $max_overlap) {
				# next nucleosome is also overlapping
				# we'll presume the next nucleosome is the oddball since it overlaps
				# two nucleosomes, so delete that one and we should be ok
				push @to_delete, $row + 1;
				
				# advance ahead
				$row++;
				$dual_count++;
			}
			elsif (
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
	
	# recalculate overlaps
	undef @overlaps;
	$overlap_count = 0;
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
	}
	
	
	# report
	printf " %s were filtered out due to extensive overlap\n", 
		format_with_commas($deleted);
	printf "   %s were tossed because they were offset\n", 
		format_with_commas($offset_count) if $offset_count;
	printf "   %s were tossed because of low occupancy\n", 
		format_with_commas($occupancy_count) if $occupancy_count;
	printf "   %s were tossed because of overlap with two nucleosomes\n",
		format_with_commas($dual_count) if $dual_count;
	printf " %s nucleosomes are remaining\n", format_with_commas($data_ref->{last_row});
	printf " %s nucleosomes (%.0f%%) are overlapping \n" . 
		"   the overlap mean is %.0f +/- %.0f bp\n", 
		format_with_commas($overlap_count),
		($overlap_count / $data_ref->{last_row}) * 100, 
		mean(@overlaps), stddev(@overlaps);
}


### Write results
# check filename
unless ($outfile) {
	$outfile = $infile;
}
my $success;
if (
	$success = write_tim_data_file(
		'data'     => $data_ref,
		'filename' => $outfile,
		'gz'       => $gz,
	) 
) {
	print " wrote data file $success\n";
}
else {
	print " unable to write file!\n";
}




__END__

=head1 NAME

verify_nucleosome_mapping.pl

A script to verify nucleosome mapping and identify overlaps.

=head1 SYNOPSIS

verify_nucleosome_mapping.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --db <text>
  --data <text | filename>
  --filter
  --max <integer>
  --recenter
  --out <filename> 
  --gz
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

=item --recenter

Optionally re-center those nucleosomes determined to be offset 
from the actual peak in the data. 

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

The program can optionally re-center those nucleosomes considered to 
be offset from their actual data peak. The start, end, midpoint, and 
name is changed to reflect the new position.

The program can also optionally filter out nucleosomes which exceed a 
set limit of overlap. The overlapping nucleosome to be deleted is 
chosen based on a set of rules: offset nucleosomes are deleted, or if 
both nucleosomes are offset, then the one with the greatest offset is 
deleted. If neither nucleosome is offset, then the nucleosome with the 
lowest occupancy is deleted, or if both occupancies are equal, then 
the rightmost is deleted. Overlap statistics are then recalculated after 
filtering. 

The same data file is re-written or a new file written with three 
additional columns appended, overlap_length, center_peak_mapping, and 
center_peak_offset. 

Overlap_length records the amount of overlap between mapped nucleosomes. 
Ideally this should be 0 as nucleosomes should not overlap; overlapping 
nucleosomes indicates either an error in mapping or multiple phasing of 
nucleosomes. 

Center_peak_mapping records whether the nucleosome was properly mapped or 
not. One of three values is recorded: centered, recentered, or offset.

Center_peak_offset records the distance in bp between the observed data 
peak and the recorded midpoint.

Basic statistics are reported to Standard Output for the overlap and 
offset lengths.

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
