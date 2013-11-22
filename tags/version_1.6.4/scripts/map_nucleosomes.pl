#!/usr/bin/perl

# A script to map nucleosomes


use strict;
use Pod::Usage;
use Getopt::Long;
use Statistics::Lite qw(count min max range sum mean stddevp);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	generate_tim_data_structure
);
use tim_db_helper qw(
	open_db_connection
	get_dataset_list
	validate_dataset_list
	get_region_dataset_hash
);
use tim_file_helper qw(
	write_tim_data_file
	convert_and_write_to_gff_file
);
#use Data::Dumper;
my $VERSION = '1.5.7';

print "\n This script will map nucleosomes\n\n";

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
	$scan_dataset,
	$tag_dataset,
	$outfile,
	$thresh,
	$window,
	$buffer,
	$bin,
	$gff,
	$type,
	$source,
	$help,
	$print_version,
	$debug
); # command line variables

# Command line options
GetOptions( 
	'db=s'     => \$database, # database name
	'data=s'   => \$dataset, # the dataset to look for movement
	'sdata=s'  => \$scan_dataset, # the dataset for scanning for nucs
	'tdata=s'  => \$tag_dataset, # the dataset containing nuc tag counts
	'out=s'    => \$outfile, # output file name
	'thresh=f' => \$thresh, # the nucleosome signal to call a nuc
	'win=i'    => \$window, # size of the window to scan the genome
	'buf=i'    => \$buffer, # the buffer between previous nuc and next window
	'bin!'     => \$bin, # allow for binning data to find nuc peaks
	'gff!'     => \$gff, # boolean to write gff file
	'type=s'   => \$type, # the GFF type for the movement
	'source=s' => \$source, # the GFF source for the movement
	'help'     => \$help, # print help
	'version'  => \$print_version, # print the version
	'debug'    => \$debug, # a limiter to help debug my program
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
	print " Biotoolbox script map_nucleosomes.pl, version $VERSION\n\n";
	exit;
}


# Check for Requirements
unless ($database) {
	die " You must define a database\n Use --help for more information\n";
}
unless ($thresh) {
	die " You must define the threshold!\n Use --help for more information\n";
}


# Assign default variables
unless ($window) {
	# set default value of 150 bp
	$window = 150;
} 
unless ($buffer) {
	$buffer = 0;
}
unless (defined $bin) {
	$bin = 0;
}



# Establish db connection
my $db = open_db_connection($database) or 
	die " unable to connect to database 'database'!\n";



# Identify the dataset to use
validate_or_request_dataset();



# Initialize data structures
my $nucleosomes_ref = generate_tim_data_structure(
	'nucleosome',
	qw(
		Chromosome
		Start
		Stop
		Midpoint
		NucleosomeID
		Occupancy
		Fuzziness
	)
) or die " unable to generate tim data structure!\n";
$nucleosomes_ref->{'db'}               = $database;
$nucleosomes_ref->{4}{'scan_dataset'}  = $scan_dataset,
$nucleosomes_ref->{4}{'window'}        = $window;
$nucleosomes_ref->{4}{'buffer'}        = $buffer;
$nucleosomes_ref->{4}{'threshold'}     = $thresh;
$nucleosomes_ref->{5}{'dataset'}       = $tag_dataset;
$nucleosomes_ref->{5}{'log2'}          = 0;
$nucleosomes_ref->{6}{'dataset'}       = $tag_dataset;
$nucleosomes_ref->{6}{'log2'}          = 0;



# Map nucleosomes
if ($debug) {
	open DEBUG_FH, ">$outfile.debug.txt";
}
my $found_nucs = map_nucleosomes();
print " Identified $found_nucs nucleosomes\n";
if ($debug) {
	close DEBUG_FH;
}


# Output
# a tim data file
exit unless $found_nucs;
unless ($outfile) {
	$outfile = $scan_dataset . '_nucleosome';
}
my $success = write_tim_data_file( {
	'data'     => $nucleosomes_ref,
	'filename' => $outfile,
} );
if ($success) {
	print " Wrote data file '$success'\n";
}
else {
	print " Unable to write data file!\n";
}
if ($gff) {
	# write a gff file if requested
	unless ($type) {
		if ($dataset) {
			# use the dataset name
			$type = $scan_dataset . '_nucleosome';
		}
		else {
			# default GFF3 compatible SO term
			$type = 'histone_binding_site';
		}
	}
	unless ($source) {
		# set default source, the name of this program
		$source = 'map_nucleosomes.pl';
	}
	$success = convert_and_write_to_gff_file( {
		'data'     => $nucleosomes_ref,
		'filename' => $outfile,
		'version'  => 3,
		'score'    => 5,
		'name'     => 4,
		'type'     => $type,
		'source'   => $source,
		'tag'      => [6],
	} );
	if ($success) {
		print " Wrote GFF file '$success'\n";
	}
	else {
		print " Unable to write GFF file!\n";
	}
}
# The End







############################# Subroutines #####################################




###### Validate an offered dataset name or interactively ask for a new one

sub validate_or_request_dataset {
	
	if ($dataset or $scan_dataset or $tag_dataset) {
		
		# set the scan and tag datasets
		if ($dataset) {
			# one general dataset for both
			$scan_dataset = $dataset;
			$tag_dataset = $dataset;
		}
		elsif (!$scan_dataset or !$tag_dataset) {
			# need to define both
			die " Both scan and tag datasets must be set! See help\n";
		}
		
		# first validate the dataset name
		my $bad_dataset = validate_dataset_list(
			$database, $scan_dataset, $tag_dataset
		);
		
		# this will return the name(s) of the bad datasets
		if ($bad_dataset) {
			die " The requested dataset '$bad_dataset' is not valid!\n";
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
		
		# get scan dataset answer 
		print " Enter the dataset to use for scanning for nucleosome peaks  ";
		my $answer = <STDIN>;
		chomp $answer;
		if (exists $datasethash{$answer}) {
			$scan_dataset = $datasethash{$answer};
		} 
		else {
			die " Unrecognized dataset request!\n";
		}
		
		# get tag dataset answer 
		print " Enter the dataset of tag count data for calculating nucleosome stats [$answer] ";
		my $answer = <STDIN>;
		chomp $answer;
		if ($answer) {
			# user provided a dataset number
			if (exists $datasethash{$answer}) {
				$tag_dataset = $datasethash{$answer};
			} 
			else {
				die " Unrecognized dataset request!\n";
			}
		}
		else {
			# user accepted the default answer, which is the same dataset
			$tag_dataset = $scan_dataset;
		}
	}
}





###### Initialize the primary output data hash


###### Now map those nucleosomes!

sub map_nucleosomes {
	# The principle of this method is essentially to find peaks of  
	# enrichment that are relatively periodic along the length of the 
	# chromosome. We will look for the maximum signal within a window 
	# of a given size that passes a set threshold signal. If it finds 
	# one, then it calls that peak the midpoint of a nucleosome. It 
	# then advances the window, using the endpoint of the just-found 
	# nucleosome as the new starting position. If it does not find a 
	# nucleosome peak that passes the threshold, then it advances the 
	# window to look for another peak. It will use the end point of the
	# previous window as the start point of the next window.
	
	# We will search one chromosome at a time, starting at the beginning 
	# and continuing until the end. Well, how else would we do it? Randomly?
	
	my $found_nucleosomes = 0;
	
	# Walk through each chromosome
	my @chromosomes = $db->seq_ids or 
		die " unable to find sequences in database!\n";
	foreach my $chromo (@chromosomes) {
		
		# skip mitochrondrial chromosome
		if ($chromo =~ /chrm|chrmt/i) {next}
		
		# get chromosome length
		my $chr_segment = $db->segment($chromo);
		my $chr_length = $chr_segment->length;
		
		# Limit to one chromosome if debug
		if ($debug) {
			if ($chromo eq 'chr2') {
				# normally chromosomes are returned in order
				# this should force only chr1 to be done
				last;
			}
			
			# artificially set chr_length to small number
			# $chr_length = 5000;
		}
		
		# Progress report
		print " Scanning chromosome $chromo....\n";
		
		# Move along the chromosome
		my $position = 1;
		while ($position < $chr_length) {
			
			# determine the window stop position
			my $win_stop = $position + $window - 1;
			
			# collect the window scanning scores
			my %window_pos2score = get_region_dataset_hash( {
					'db'       => $db,
					'dataset'  => $scan_dataset,
					'chromo'   => $chromo,
					'start'    => $position,
					'stop'     => $win_stop,
					'value'    => 'score',
			} );
			if ($debug) {
				print DEBUG_FH "### Window $chromo:$position..$win_stop\n";
				print DEBUG_FH "  Window scores:\n  ";
				foreach (sort {$a <=> $b} keys %window_pos2score) {
					print DEBUG_FH "$_:$window_pos2score{$_}, ";
				}
				print DEBUG_FH "\n";
			}
			
			# check whether we actually have a valid nucleosome peak
			# if not, then try again after binning the data
			my $found_peak = 0; 
			my $window_scan_max = max(values %window_pos2score);
			
			if ( $window_scan_max >= $thresh) {
				# we have a peak that passes the threshold
				# yeah, a nucleosome
				$found_peak = 1;
			}
			elsif ($window_scan_max < $thresh and $bin) {
				# no obvious peak, so we will try binning to find a peak
				
				my %binned_pos2score;
				my @positions = sort {$a <=> $b} keys %window_pos2score;
				for (my $i = 0; $i < scalar(@positions); $i += 3) {
					# we will bin three positions into one
					# the sum of three positions' score indexed at the 
					# middle position
					$binned_pos2score{ $positions[ $i + 1 ] } = 
						sum(
							$window_pos2score{ $positions[$i] },
							$window_pos2score{ $positions[$i + 1] },
							$window_pos2score{ $positions[$i + 2] },
						);
				}
				
				# look again for a peak
				$window_scan_max = max(values %binned_pos2score);
				if ($window_scan_max >= $thresh) {
					# we found one
					$found_peak = 1;
					
					# reassign position2score hash
					%window_pos2score = %binned_pos2score;
					if ($debug) {
						print DEBUG_FH "  Binned window scores:\n  ";
						foreach (sort {$a <=> $b} keys %window_pos2score) {
							print DEBUG_FH "$_:$window_pos2score{$_}, ";
						}
					}
				}
			}
			else {
				# no peak
				$position += $window;
				next;
			}
			
			
			# identify the peak position
			my $peak_position;
			if ($found_peak) {
				# we know there is a peak within this window, now must find it
				
				# collect the tag scores for the window
				# the scan data scores may not be reliable, as the statistical 
				# method used to identify peaks may not directly correspond to
				# the actual peak of tags, which is what we really want
				# this is true for False Discovery Rate scores
				my %window_pos2tags = get_region_dataset_hash( {
						'db'       => $db,
						'dataset'  => $tag_dataset,
						'chromo'   => $chromo,
						'start'    => $position,
						'stop'     => $win_stop,
						'value'    => 'score',
				} );
				
				# get all of the positions that correspond to that peak value
				my $tag_peak_value = max(values %window_pos2tags);
				my @peak_positions;
				foreach (sort {$a <=> $b} keys %window_pos2tags) {
					if ($window_pos2tags{$_} == $tag_peak_value) {
						push @peak_positions, $_;
					}
				}
				
				# now identify the best peak position, ideally this is just one
				if (scalar @peak_positions == 1) {
					# one peak makes this easy
					$peak_position = $peak_positions[0];
				}
				
				else {
					# more than one equal peaks
					# if the positions are really close together, take the 
					# mean position
					if (range(@peak_positions) <= 10) {
						# the range between the peaks is 10 bp or less
						# 10 bp is totally arbitrary but reasonable
						$peak_position = int( mean(@peak_positions) + 0.5);
					}
					
					else {
						# they are too far apart, separate nucleosomes?
						# we'll take the one closest to the middle of the range
						my %best;
						my $window_midpoint =  
									int( ( ($position + $win_stop) / 2) + 0.5);
						foreach (reverse @peak_positions) {
							# we're reversing the order to take the rightmost
							# position first
							# that way if the positions are equidistant and the 
							# key gets overwritten, we'll take the leftmost 
							# position
							$best{ abs( $_ -  $window_midpoint ) } = $_;
						}
						
						# take the closest to the middle position
						$peak_position = $best{ min( keys %best ) }
					}
				}
				
				if ($debug) {
					print DEBUG_FH " Peak found at position $peak_position\n";
				}
			}
			
			# process the peak if found
			if (defined $peak_position) {
				# we have a defined nucleosome position!
				
				# Determine the fuzziness of this nucleosome
					# This is essentially the standard deviation of the counts
					# around the midpoint of the nuclesome. Well positioned 
					# nucleosomes will have all their midpoints centered in 
					# position, whereas fuzzy nucleosomes will have their 
					# midpoints spread out across the region. 
					
					# We will need to collect new scores surrounding the 
					# new found nucleosome peak.
					# arbitrarily collecting 50 bp worth of scores from each
					# side of the peak
					# we can't use the previously collected tag scores because
					# the nucleosome peak may have been at the edge of the window
				my %nucleosome_pos2score = get_region_dataset_hash( {
					'db'       => $db,
					'dataset'  => $tag_dataset,
					'chromo'   => $chromo,
					'start'    => $peak_position - 50,
					'stop'     => $peak_position + 50,
					'value'    => 'score',
				} );
				my @nucleosome_positions;
				for (my $i = -37; $i <= 37; $i++) {
					# we will advance along the relative position around the 
					# nucleosome peak and record the number of nucleosome 
					# midpoints found at each relative position
					
					# The size of the relative window to look for fuzziness
					# is essentially 1/2 of a nucleosome length (~74 bp). This
					# is totally arbitrary, but reasonable.
					
					# convert relative to absolute
					my $pos = $peak_position + $i; 
					
					# collect the counts
					if ($nucleosome_pos2score{$pos} > 0) {
						foreach (1 .. $nucleosome_pos2score{$pos}) {
							push @nucleosome_positions, abs($i);
						}
					}
				}
				unless (@nucleosome_positions) {
					# nothing found!? not acceptable
					$position += $window;
					next;
				}
				
				# Calculate standard deviation for fuzziness
				my $fuzziness = sprintf "%.0f", stddevp(@nucleosome_positions);
				
				# Calculate the number of nucleosome midpoints used in calling 
				# this nucleosome, this will be the score
				my $score = count(@nucleosome_positions);
				
				# Determine the nucleosome coordinates
					# we are assuming a standard 147 bp sized nucleosome
					# this may not be accurate for those partial and/or
					# fragile nucleosomes
				my $nuc_start = $peak_position - 73;
				my $nuc_stop = $peak_position + 73;
				
				# Generate name
					# not quite following Pugh's example, we will name the 
					# nucleosome based on chromosome number and midpoint 
					# position (instead of start position)
				$chromo =~ /^(?:chr)?(.+)$/i;
				my $nuc_name = 'N'. $1 . ':' . $peak_position;
				
				# Record the nucleosome information
				push @{ $nucleosomes_ref->{'data_table'} }, [
					$chromo,
					$nuc_start,
					$nuc_stop,
					$peak_position,
					$nuc_name,
					$score,
					$fuzziness,
				];
				
				# Reset the position to the current nucleosome endpoint
				# this may (likely) generate in some overlap if the threshold
				# is set too low
				$position = $nuc_stop + $buffer;
				
				$found_nucleosomes++;
				#print "      nucleosome $found_nucleosomes is $nuc_name with score $score and fuzziness $fuzziness\n";
				#print "  new position is $position\n";
			}
			
			else {
				# no nucleosome found, move on to next window
				$position += $window;
			}
		}
	}
	
	# record the number of nucleosomes
	$nucleosomes_ref->{'last_row'} = $found_nucleosomes;
	return $found_nucleosomes;
}





__END__

=head1 NAME

map_nucleosomes.pl

A script to map nucleosomes.

=head1 SYNOPSIS

map_nucleosomes.pl --db <database> --thresh <number> [--options...]
  
  --db <database_name>
  --data <dataset_name>
  --sdata <dataset_name>
  --tdata <dataset_name>
  --thesh <number>
  --win <integer>
  --buf <integer>
  --(no)bin
  --gff
  --type <gff_type>
  --source <gff_source>
  --out <filename>
  --version
  --help



=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <database_name>

Specify the name of the BioPerl gff database to pull the source data. 

=item --data, --sdata, --tdata <dataset_name>

Provide the name of the dataset(s) containing the nucleosome midpoint
occupancy data from which to identify nucleosomal positions. Two 
datasets are required, the scan dataset (--sdata) and the tag dataset 
(--tdata). The same dataset could be used for both; in which case 
set both to the same using the --data option. See the DESCRIPTION for 
details regarding the datasets. If the datasets are not specified on 
the commandline, then they may be interactively chosen from a list 
from the database.

=item --thresh <number>

Provide the minimum score value required to call a nucleosome position. 
This is only used when scanning the scan dataset. 

=item --win <integer>

Provide the window size for which to scan the chromosome. Setting this  
value too large and overlapping nucleosomes may result, while setting 
this value too low may miss some nucleosomes. The default value is 150 bp. 

=item --buf <integer>

Provide the buffer size in bp which will be added between the end of the 
previous found nucleosome and the beginning of the window to scan for the 
next nucleosome. Setting this value may limit the number of overlapping 
nucleosomes. Default is 0 bp.

=item --bin

Indicate whether the scan data should be binned in an attempt to find 
and map a nucleosome midpoint peak. Binning only occurs when a single 
position fails to pass the threshold value. Bins are generated from 
three adjacent positions and their values summed. This may help  
those genomic regions with very few nucleosome reads. The default is 
to not bin the scan data.

=item --gff

Indicate whether a GFF file should be written in addition to the standard 
text data file. The GFF file version is 3. The default value is false.

=item --type <gff_type>

Provide the text to be used as the GFF type (or method) used in 
writing the GFF file. The default value is the name of the scan 
dataset appended with "_nucleosome".

=item --source <gff_source>

Provide the text to be used as the GFF source used in writing the 
GFF file. The default value is the name of this program.

=item --out <filename>

Specify the output file name. The default is the name of the scan 
dataset appended with "_nucleosome".

=item --version

Print the version number.

=item --help

Display the POD documentation of the script. 

=back

=head1 DESCRIPTION

This program will identify and map the positions of nucleosomes given a 
dataset of nucleosome occupancy data. The dataset should ideally be 
enumerated counts of sequenced nucleosomal fragment midpoints, although 
very high resolution microarray data could also be used in principle. 

Nucleosome calls are made by scanning the chromosomes in windows of 
specified size (default 150 bp, set with --win option) looking for 
the maximum peak that exceeds the minimum threshold value. The position 
of the maximum peak is called as the new nucleosome midpoint. The 
window is then advanced starting at the previous just-identified nucleosome 
endpoint, or at the end of the previous window if no nucleosome was 
identified. This position may be further advanced by setting the 
buffer value, which inserts space between the previous nucleosome end and 
the next window. By advancing the window relative to the previously identified 
nucleosome, the program can adapt to variable nucleosome spacing. 

Two datasets must be provided for the mapping, although the same could be  
used for both functions. The datasets should represent sequenced nucleosome 
fragment midpoints. For the most accurate mapping, the midpoints 
should be mapped at 1 bp resolution, but lower resolutions (5 or 10 bp) 
could work. 

The scan dataset (--sdata) is used to scan for a potential
nucleosome. The dataset may be a simple difference, normalized difference, 
P-value, False Discovery Rate, or simply raw enumerated counts. Probability 
scores should be converted using -10Log10(P) for proper interpretation. 

The tag dataset (--tdata) is used for calculating nucleosome statistics, 
including the precise nucleosome midpoint peak, occupancy, and fuzziness.
The tag dataset must be enumerated counts of nucleosome midpoints and
be whole integers, whether raw counts or a simple difference. Negative counts
from a difference dataset are not counted. 

For genomic regions with few nucleosome reads and no positions that pass the 
threshold value, the window positions may optionally be binned together to 
increase sensitivity (see the --bin option).

Two attributes of each identified nucleosome are calculated, Occupancy and 
Fuzziness. Occupancy represents the sum of the tag counts that support the 
nucleosome position; higher scores indicate a more highly occupied 
nucleosome. Fuzziness indicates how well all of the scores are aligned in 
register. The Fuzziness value is the standard deviation of the population 
of scores from the peak; a high value indicates the nucleosome has a 
fuzzy or variable position, whereas a low value indicates the nucleosome 
is highly positioned. For both attributes, the scores are counted in a 
window representing 1/2 of a nucleosome length (74 bp) centered on the 
defined nucleosome midpoint. The size of this window is arbitrary but 
reasonable.

The identified nucleosomes are set to 147 bp in length, centered at the 
maximum peak. A second script, C<get_actual_nuc_sizes.pl>, will determine 
the actual nucleosome sizes based on paired-end reads in a BAM file.




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






