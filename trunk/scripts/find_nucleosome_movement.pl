#!/usr/bin/perl

# A script to look for nucleosome movement

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw();
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_db_helper qw(
	open_db_connection
	get_dataset_list
	validate_dataset_list
	get_region_dataset_hash
);
use tim_file_helper;

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
	$help
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
) or die " unknown arguments! use --help\n";


# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
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
validate_or_request_dataset();


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
	
	# generate the data hash
	my %datahash;
	
	# populate the standard data hash keys
	$datahash{'program'}        = $0;
	$datahash{'db'}             = $database;
	$datahash{'feature'}        = 'nucleosome_movement';
	$datahash{'gff'}            = 0;
	$datahash{'number_columns'} = 6;
	
	# set column metadata
	$datahash{0} = {
		# the chromosome
		'name'     => 'Chromosome',
		'index'    => 0,
	};
	$datahash{1} = {
		# the start position 
		'name'     => 'Start',
		'index'    => 1,
	};
	$datahash{2} = {
		# the stop position
		'name'     => 'Stop',
		'index'    => 2,
	};
	$datahash{3} = {
		# the direction of movement
		# this will essentially be recorded as strand in the GFF
		'name'     => 'Direction',
		'index'    => 3,
	};
	$datahash{4} = {
		# the score position
		# this will be absolute value of the difference between the two scores
		'name'     => 'Score',
		'index'    => 4,
		'dataset'  => $dataset,
		'gain'     => $gain,
		'loss'     => $loss,
		'delta'    => $delta,
		'log2'     => 1, # assumed
	};
	$datahash{5} = {
		# a unique ID for the movement
		'name'     => 'EventID',
		'index'    => 5,
	};
	
	# Set the data table
	my @data_table = ( [ qw(
		Chromosome
		Start
		Stop
		Direction
		Score
		EventID
	) ] );
	$datahash{'data_table'} = \@data_table;
	
	# return the reference to the generated data hash
	return \%datahash;
}



###### Initialize the source data hash

sub initialize_source_data_hash {
	# we'll be collecting the source data from either a database or file
	# and putting it into a standard tim data memory structure
	# for ease of use and manipulation
	
	### Initialize the data structure
	my %datahash;
	# populate the standard data hash keys
	$datahash{'program'}        = $0;
	$datahash{'db'}             = $database;
	$datahash{'feature'}        = 'microarray_data';
	$datahash{'gff'}            = 0;
	$datahash{'number_columns'} = 3;
	
	# set column metadata
	$datahash{0} = {
		# the chromosome
		'name'     => 'Chromosome',
		'index'    => 0,
	};
	$datahash{1} = {
		# the start position 
		'name'     => 'Start',
		'index'    => 1,
	};
	$datahash{2} = {
		# the score
		'name'     => 'Score',
		'index'    => 2,
		'log2'     => 1, # assumed
		'dataset'  => $dataset,
	};
	
	# Set the data table
	my @data_table = ( [ qw(
		Chromosome
		Start
		Score
	) ] );
	
	
	
	### Collect Source data
	if ($database) {
		# Source data from the database
		
		print " Collecting source data from database $database...\n";
		
		# connect
		my $db = open_db_connection($database) or 
			die " unable to connect to database!\n";
		
		# get list of chromosomes
		my @chromosomes = $db->features(-type => 'chromosome'); 
		
		# walk through each chromosome
		foreach my $chrobj (
			# trying a Schwartzian transformation here
			map $_->[0],
			sort { $a->[1] <=> $b->[1] }
			map [$_, ($_->name =~ /(\d+)/)[0] ], 
			@chromosomes
		) {
			# sort chromosomes by increasing number
			# we're using RE to pull out the digit number in the chromosome name
			# and sorting increasingly by it
			
			# chromosome name
			my $chr = $chrobj->name; # this is actually returning an object, why????
			$chr = "$chr"; # force as string
			
			# skip mitochondrial chromosome
			if ($chr =~ /chrm|chrmt/i) {next}
			
			# collect the dataset values for the current chromosome
			# store in a hash the position (key) and values
			my %chromodata = get_region_dataset_hash( {
				'db'       => $db,
				'dataset'  => $dataset,
				'name'     => $chr,
				'type'    => 'chromosome',
			} ) or die " no data collected for chromosome $chr!";
			# chromdata is organized as position => score
			
			# collect the data from the chromosome
			foreach my $position (sort {$a <=> $b} keys %chromodata) {
				push @data_table, [ ($chr, $position, $chromodata{$position} ) ];
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
			for (my $i = 0; $i < $input_ref->{'number_columns'}; $i++) {
				# check the names of each column
				# looking for chromo, start, stop
				if ($input_ref->{$i}{'name'} =~ /^chrom|refseq/i) {
					$chr_index = $i;
				}
				elsif ($input_ref->{$i}{'name'} =~ /start|position/i) {
					$start_index = $i;
				}
				elsif ($input_ref->{$i}{'name'} =~ /stop|end/i) {
					$stop_index = $i;
				}
				elsif ($input_ref->{$i}{'name'} =~ /score/i) {
					$score_index = $i;
				}
				elsif ($input_ref->{$i}{'name'} eq $dataset) {
					$score_index = $i;
				}
			}
		}
		
		# check that we have the scores
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
				push @data_table, [ (
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
				push @data_table, [ (
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
				push @data_table, [ (
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
	$datahash{'data_table'} = \@data_table;
	$datahash{'last_row'} = scalar(@data_table) - 1;
	
	# Return
	return \%datahash;
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


=head1 TODO






