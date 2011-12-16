#!/usr/bin/perl

# A script to look for nucleosome movement

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw();
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	generate_tim_data_structure
	find_column_index
);
use tim_db_helper qw(
	open_db_connection
	get_dataset_list
	validate_dataset_list
	get_region_dataset_hash
);
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
	convert_and_write_to_gff_file
);
my $VERSION = '1.5.4';

print "\n This script will identify nucleosome movement\n\n";

### Quick help
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
	$database,
	$dataset,
	$source_data_file,
	$type,
	$source,
	$outfile,
	$gain,
	$loss,
	$delta,
	$win,
	$help,
	$print_version,
); # command line variables

# Command line options
GetOptions( 
	'db=s'     => \$database, # database name
	'data=s'   => \$dataset, # the dataset to look for movement
	'dataf=s'  => \$source_data_file, # a file with the source data
	'type=s'   => \$type, # the GFF type for the movement
	'source=s' => \$source, # the GFF source for the movement
	'out=s'    => \$outfile, # output file name
	'gain=f'   => \$gain, # the amount of nucleosome occupancy gain
	'loss=f'   => \$loss, # the amount of nucleosome occupancy loss
	'delta=f'  => \$delta, # the minimum difference between the gain and loss values
	'win=i'    => \$win, # size of the window to scan the genome
	'help'     => \$help, # print help
	'version'  => \$print_version, # print the version
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
	print " Biotoolbox script find_nucleosome_movement.pl, version $VERSION\n\n";
	exit;
}



# Check for Requirements
unless ($database or $source_data_file) {
	die " You must define a data source\n Use --help for more information\n";
}
unless ($gain and $loss) {
	die " You must define gain and loss values!\n Use --help for more information\n";
}


# Assign default variables
if ($loss > 0) { 
	# convert loss value to negative
	$loss = 0 - $loss;
} 
unless ($delta) {
	# calculate delta value
	$delta = $gain + abs($loss);
}
unless ($win) {
	# set default value of 100 bp
	$win = 100;
} 
unless ($type) {
	if ($dataset) {
		# use the dataset name
		$type = $dataset . '_movement';
	}
	else {
		# very generic, user should've really provided the type
		$type = 'nucleosome_movement';
	}
}
unless ($source) {
	# set default source, the name of this program
	$source = 'find_nucleosome_movement.pl';
}
unless ($outfile) {
	$outfile = $type;
}



# Identify the dataset to use
if ($database) {
	validate_or_request_dataset();
}


# Initialize data structures
my $source_data_ref = initialize_source_data_hash();
my $movement_data_ref = initialize_movement_data_hash();



# Map nucleosome movement
my $movement_number = map_nucleosome_movement();
print " Identified $movement_number movement events\n";



# Output
# a tim data file
my $success = write_tim_data_file( {
	'data'     => $movement_data_ref,
	'filename' => $outfile,
} );
if ($success) {
	print " Wrote data file '$success'\n";
}
else {
	print " Unable to write data file!\n";
}
# a gff file
$success = convert_and_write_to_gff_file( {
	'data'     => $movement_data_ref,
	'filename' => $outfile,
	'version'  => 3,
	'score'    => 4,
	'strand'   => 3,
	'name'     => 5,
	'type'     => $type,
	'source'   => $source,
} );
if ($success) {
	print " Wrote GFF file '$success'\n";
}
else {
	print " Unable to write GFF file!\n";
}

# The End







############################# Subroutines #####################################




###### Validate an offered dataset name or interactively ask for a new one

sub validate_or_request_dataset {
	
	if ($dataset) {
		# first validate the dataset name
		my $bad_dataset = validate_dataset_list($database, $dataset);
		
		# this will return the name(s) of the bad datasets
		if ($bad_dataset) {
			die " The requested dataset '$bad_dataset' is not valid!\n";
		} 
		else {
			# returning nothing from the subroutine is good
			print " Using requested data set $dataset....\n";
		}
	}	
	
	# Otherwise ask for the data set
	else {
		
		# Present the data set list to the user and get an answer
		my %datasethash = get_dataset_list($database); # list of data sets
		print "\n These are the microarray data sets in the database:\n";
		foreach (sort {$a <=> $b} keys %datasethash) {
			# print out the list of microarray data sets
			print "  $_\t$datasethash{$_}\n"; 
		}
		
		# get answer 
		print " Enter the number of the data set you would like to analyze  ";
		my $answer = <STDIN>;
		chomp $answer;
		
		# check answer
		if (exists $datasethash{$answer}) {
			$dataset = $datasethash{$answer};
			print " Using data set '$dataset'....\n";
		} 
		else {
			die " Unrecognized dataset request!\n";
		}
	}
}





###### Initialize the primary output data hash

sub initialize_movement_data_hash {
	# the output data structure for the identified nucleosome movement events
	
	# generate a new data structure
	my $data = generate_tim_data_structure(
		'nucleosome_movement',
		'Chromosome',
		'Start',
		'Stop',
		'Direction',
		'Score',
		'EventID',
	) or die " unable to generate tim data structure!\n";
	
	# add metadata
	$data->{'db'}           = $database || $source_data_file;
	$data->{4}{'dataset'}   = $dataset || $source_data_file;
	$data->{4}{'gain'}      = $gain;
	$data->{4}{'loss'}      = $loss;
	$data->{4}{'delta'}     = $delta;
	$data->{4}{'log2'}      = 1;
	
	# return the reference to the generated data hash
	return $data;
}



###### Initialize the source data hash

sub initialize_source_data_hash {
	# we'll be collecting the source data from either a database or file
	# and putting it into a standard tim data memory structure
	# for ease of use and manipulation
	
	### Generate a new data structure
	my $data = generate_tim_data_structure(
		'microarray_data',
		'Chromosome',
		'Start',
		'Score',
	) or die " unable to generate tim data structure!\n";
	
	# add metadata
	$data->{'db'}           = $database || $source_data_file;
	$data->{2}{'dataset'}   = $dataset || $source_data_file;
	$data->{2}{'log2'}      = 1; # assumed
	
	
	
	### Collect Source data
	if ($database) {
		# Source data from the database
		
		print " Collecting source data from database $database...\n";
		
		# connect
		my $db = open_db_connection($database) or 
			die " unable to connect to database!\n";
		
		# get list of chromosomes
		my @chromosomes = $db->seq_ids; 
		
		# walk through each chromosome
		foreach my $chr (@chromosomes) {
			
			# skip mitochondrial chromosome
			if ($chr =~ /chrm|chrmt/i) {next}
			
			# get the length
			my $segment = $db->segment($chr);
			my $length  = $segment->length;
			
			# collect the dataset values for the current chromosome
			# store in a hash the position (key) and values
			my %chromodata = get_region_dataset_hash( {
				'db'       => $db,
				'dataset'  => $dataset,
				'chromo'   => $chr,
				'start'    => 1,
				'stop'     => $length,
			} ) or die " no data collected for chromosome $chr!";
			# chromdata is organized as position => score
			
			# collect the data from the chromosome
			foreach my $position (sort {$a <=> $b} keys %chromodata) {
				push @{$data->{'data_table'}}, 
					[ ($chr, $position, $chromodata{$position} ) ];
			}
			print "  ...collected ", scalar keys %chromodata, " values from ", 
				"chromosome $chr\n";
		}
	}
	elsif ($source_data_file) {
		# User specified data file
		
		# load the file 
		my $input_ref = load_tim_data_file($source_data_file) or die 
			" unable to load source data from file '$source_data_file'!\n";
		
		# identify the appropriate columns
		my ($chr_index, $start_index, $stop_index, $score_index);
		if ($input_ref->{'gff'}) {
			# a gff file
			# this makes column identification easy
			$chr_index = 0;
			$start_index = 3;
			$stop_index = 4;
			$score_index = 5;
		}
		else {
			# a tim data file
			# let's hope that the columns are identifiable
			$chr_index = find_column_index($input_ref, '^chr|refseq') || undef;
			$start_index = find_column_index($input_ref, '^start|pos') || undef;
			$stop_index = find_column_index($input_ref, '^stop|end') || undef;
			$score_index = find_column_index($input_ref, "score|$dataset") || undef;
		}
		
		# check that we have identified columns
		unless (
			defined $chr_index and 
			defined $start_index 
		) {
			die " unable to identify genomic coordinate data columns in data file!\n";
		}
		unless (defined $score_index) {
			die " unable to identify source dataset in data file!\n";
		}
		
		# Load the data 
		my $indata_ref = $input_ref->{'data_table'}; # quick reference
		# this process depends on whether we have a stop index and 
		# need to calculate midpoint
		if (
			defined $stop_index and
			$indata_ref->[1][$start_index] != $indata_ref->[1][$stop_index]
		) {
			# we must calculate a midpoint between different start and stop
			for (my $row = 1; $row <= $input_ref->{'last_row'}; $row++) {
				# walk through the input data and copy the relevant data
				
				# skip mitochondrial chromosome
				if ($indata_ref->[$row][$chr_index] =~ /chrm|chrmt/i) {next}
				
				# calculate position and store
				my $position = sprintf "%.0f", 
					( $indata_ref->[$row][$start_index] / 
					$indata_ref->[$row][$stop_index] );
				push @{$data->{'data_table'}}, [ (
					$indata_ref->[$row][$chr_index],
					$position,
					$indata_ref->[$row][$score_index]
				) ];
			}
		}
		if (
			defined $stop_index and
			$indata_ref->[1][$start_index] == $indata_ref->[1][$stop_index]
		) {
			# start and stop are equivalent, so no worry
			for (my $row = 1; $row <= $input_ref->{'last_row'}; $row++) {
				# walk through the input data and copy the relevant data
				
				# skip mitochondrial chromosome
				if ($indata_ref->[$row][$chr_index] =~ /chrm|chrmt/i) {next}
				
				# store
				push @{$data->{'data_table'}}, [ (
					$indata_ref->[$row][$chr_index],
					$indata_ref->[$row][$start_index],
					$indata_ref->[$row][$score_index]
				) ];
			}
		}
		else {
			# only have start index
			for (my $row = 1; $row <= $input_ref->{'last_row'}; $row++) {
				# walk through the input data and copy the relevant data
				
				# skip mitochondrial chromosome
				if ($indata_ref->[$row][$chr_index] =~ /chrm|chrmt/i) {next}
				
				# store
				push @{$data->{'data_table'}}, [ (
					$indata_ref->[$row][$chr_index],
					$indata_ref->[$row][$start_index],
					$indata_ref->[$row][$score_index]
				) ];
			}
		}
	}
	
	
	else {
		# How did we get here?
		die " why don't we have any source data!?\n";
	}
	
	# Finish the data structure
	$data->{'last_row'} = scalar(@{$data->{'data_table'}}) - 1;
	
	# Return
	return $data;
}



###### Initialize the source data hash

sub map_nucleosome_movement {
	# The source data structure should ordered by genomic position
	# We will simply march through looking for movement
	
	print " Scanning for movements...\n";
	
	# Initialize count and references
	my $count = 0; # count for the ID
	my $move_ref = $movement_data_ref->{'data_table'};
	
	# Walk through data
	for (my $i = 1; $i < $source_data_ref->{'last_row'}; $i++) {
		
		# Reference elements
		my $row = $source_data_ref->{'data_table'}->[$i];
		my $next = $source_data_ref->{'data_table'}->[$i + 1];
		
		# Look for movement by series of logical checks
		if (
			# same chromo
			$row->[0] eq $next->[0] 
			and
			# within window
			( $next->[1] - $row->[1] ) <= $win
		) {
			# now check values
			# a movement event is defined by three values
				# first value >= gain
				# second value <= loss
				# absolute difference between values >= delta
			
			my $score = abs($row->[2]) + abs($next->[2]);
			
			# direction is defined by order of gain and loss
			# Gain, Loss
			if (
				$row->[2]  >= $gain and
				$next->[2] <= $loss and
				$score >= $delta
			) {
				# a movement leftward
				
				# calculate unique ID
				$count++;
				my $id = $type . sprintf("%07d", $count);
				
				# record the movement
				push @{ $move_ref }, [
					$row->[0], # chromosome
					$row->[1], # start
					$next->[1], # stop
					'-', # left direction
					$score,
					$id
				];
			}
			
			# Loss, Gain
			elsif (
				$row->[2]  <= $loss and
				$next->[2] >= $gain and
				$score >= $delta
			) {
				# a movement rightward
				
				# calculate unique ID
				$count++;
				my $id = $type . sprintf("%07d", $count);
				
				# record the movement
				push @{ $move_ref }, [
					$row->[0], # chromosome
					$row->[1], # start
					$next->[1], # stop
					'+', # right direction
					$score,
					$id
				];
			}
		}
	}
	
	return $count;
}





__END__

=head1 NAME

find_nucleosome_movement.pl

A script to map nucleosome movement based on mononucleosome DNA hybridzed 
to Agilent cerevisiae 244K arrays.

=head1 SYNOPSIS

find_nucleosome_movement.pl [--options...]
  
  --db <database_name>
  --dataf <filename>
  --data <dataset_name>
  --gain <number>
  --loss <number>
  --type <gff_type>
  --source <gff_source>
  --out <filename>
  --delta <number>
  --win <integer>
  --help



=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <database_name>

Specify the name of the BioPerl gff database to pull the source data.

=item --dataf <filename>

Provide the name of a data file containing the dataset. The file 
should provide genomic coordinates (e.g. GFF file). If more than 
one dataset is within the file, then the dataset name should be 
provided.

=item --data <dataset_name>

Provide the name of the dataset containing the nucleosome occupancy data
from which to identify the movements. The dataset should be a ratio 
or normalized difference between experimental and control nucleosome 
occupancies. It is assumed to be in log2 data space.

=item --gain <number>

Provide the minimum score value required to identify a gain event. 
Adjacent concordant gain and loss identifies nucleosome movement.

=item --loss <number>

Provide the minimum score value required to identify a loss event. 
Adjacent concordant gain and loss identifies nucleosome movement. 
The loss number is automatically converted to a negative if 
necessary (assuming the data is in log2 space).

=item --type <gff_type>

Provide the text to be used as the GFF type (or method) used in 
writing the GFF file. It is also used as the base for the event 
name. The default value is the dataset name appended with 
'_movement'.

=item --source <gff_source>

Provide the text to be used as the GFF source used in writing the 
GFF file. The default value is the name of this program.

=item --out <filename>

Specify the output file name. The default is to use the type provided.

=item --delta <number>

Provide the minimum absolute difference value to call a nucleosome 
movement event. The default is the sum of absolute gain and loss 
values.

=item --win <integer>

Provide the maximum window size in bp to look for adjacent probes when 
calling nucleosome movement events. The number should be greater than 
the average microarray probe spacing but less than the length of a 
nucleosome. The default value is 100 bp.

=item --help

Display the POD documentation of the script. 

=back

=head1 DESCRIPTION

This program will look for nucleosome movement. It was designed explicitly
to look for shifts in nucleosome occupancy in medium to high resolution
microarrays, such as the Agilent 244K array for cerevisiae. A nucleosome
movement event is defined as a concordant gain and loss of nucleosome
signal (occupancy) at two neighboring probes. The neighboring probes must
be within a specified window size (default 100 bp). The dataset to be
scanned should be a log2 ratio or normalized difference between
experimental and control nucleosome occupancies. Both minimum gain and loss
values must be met, and optionally the minimum delta value (absolute
difference between the adjacent probe values).

A tim data file will be written and a GFF file suitable for 
loading in GBrowse are written.




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


=head1 TODO






