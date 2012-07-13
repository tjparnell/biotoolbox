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
	process_and_verify_dataset
	get_region_dataset_hash
	get_chromo_region_score
);
use tim_file_helper qw(
	write_tim_data_file
	convert_and_write_to_gff_file
);
#use Data::Dumper;
my $VERSION = '1.8.3';

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


# Check for Requirements and assign default variables
unless ($thresh) {
	die " You must define the threshold!\n Use --help for more information\n";
}

unless ($window) {
	# set default value of 150 bp
	print " Using default window size of 175 bp\n";
	$window = 175;
} 
unless ($buffer) {
	$buffer = 20;
}
unless (defined $bin) {
	$bin = 0;
}



# Establish db connection
my $db;
if ($database) {
	# a database was provided by the user
	$db = open_db_connection($database) or 
		die " unable to connect to database 'database'!\n";
}
else {
	if ($dataset or $scan_dataset) {
		# we might be able to use an input bigwig file
		my $alternate = $dataset || $scan_dataset;
		
		unless ($alternate =~ m/\.bw$/i) {
			die " You must define a database\n Use --help for more information\n";
		}
		$db = open_db_connection($alternate) or 
			die " unable to connect to database 'database'!\n";
	}
	else {
		die " You must define a database\n Use --help for more information\n";
	}
}



# Identify the dataset to use
validate_or_request_dataset();



# Initialize data structures
my $nucleosomes_ref = initialize_nucleosome_data();



### Map nucleosomes
if ($debug) {
	open DEBUG_FH, ">$outfile.debug.txt";
}

my $found_nucs = map_nucleosomes();
print " Identified $found_nucs nucleosomes\n";

if ($debug) {
	close DEBUG_FH;
}


### Output
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

# write a gff file if requested
if ($gff) {
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


sub initialize_nucleosome_data {
	
	# generate the data structure based on what datasets were provided
	my $data;
	if ($tag_dataset) {
		$data = generate_tim_data_structure(
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
	}
	else {
		# no nucleosome statistics will be collected
		$data = generate_tim_data_structure(
			'nucleosome',
			qw(
				Chromosome
				Start
				Stop
				Midpoint
				NucleosomeID
			)
		) or die " unable to generate tim data structure!\n";
	}
	
	# add additional metadata
	$data->{'db'}               = $database;
	$data->{4}{'scan_dataset'}  = $scan_dataset,
	$data->{4}{'window'}        = $window;
	$data->{4}{'buffer'}        = $buffer;
	$data->{4}{'threshold'}     = $thresh;
	if ($tag_dataset) {
		$data->{5}{'dataset'}   = $tag_dataset;
		$data->{6}{'dataset'}   = $tag_dataset;
	}
	
	return $data;
}



###### Validate an offered dataset name or interactively ask for a new one

sub validate_or_request_dataset {
	
	# Process the scan dataset
	if ($dataset and !$scan_dataset) {
		# only one dataset provided, then assign it to scan
		$scan_dataset = $dataset;
	}
	$scan_dataset = process_and_verify_dataset( {
		'db'      => $database,
		'dataset' => $scan_dataset,
		'prompt'  => "Enter the dataset to use for scanning for nucleosome peaks",
		'single'  => 1,
	} );
	unless ($scan_dataset) {
		die " A valid scan dataset must be provided!\n";
	}
	
	
	# Process the tag dataset
	if ($dataset and !$tag_dataset) {
		# only one dataset provided, then assign it to scan
		$tag_dataset = $dataset;
	}
	if ($tag_dataset or $database) {
		$tag_dataset = process_and_verify_dataset( {
			'db'      => $database,
			'dataset' => $tag_dataset,
			'prompt'  => "Enter the dataset of tag count data for calculating nucleosome stats",
			'single'  => 1,
		} );
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
		if ($chromo =~ /^chrm|chrmt|mt|mito/i) {next}
		
		# get chromosome length
		my $chr_segment = $db->segment($chromo);
		my $chr_length = $chr_segment->length;
		
		# Limit to one chromosome if debug
		if ($debug) {
			if ($chromo eq $chromosomes[1]) {
				# limit to only the first chromosome
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
			
			
			# scan the window for a potential nucleosome peak
			my $found_peak = scan_window($chromo, $position, $win_stop);

			
			# identify the peak position
			if ($found_peak) {
				# we found a peak, acurately map it with tag data
				my $peak_position = identify_peak_position(
					$chromo, $position, $win_stop);
				
				# now process the peak and return a new position
				$position = process_nucleosome($chromo, $peak_position);
				
				# increment found counter
				$found_nucleosomes++;
			}
			else {
				# nothing found, move on
				if ($debug) {
					print DEBUG_FH " Did not find a peak\n";
				}
				$position += $window;
				next;
			}
		}
	}
	
	# record the number of nucleosomes
	$nucleosomes_ref->{'last_row'} = $found_nucleosomes;
	return $found_nucleosomes;
}



sub scan_window {
	
	# the window we're looking at
	my ($chromo, $position, $win_stop) = @_;
	
	# collect the window scanning scores
	my %window_pos2score = get_region_dataset_hash( {
			'db'       => $db,
			'dataset'  => $scan_dataset,
			'chromo'   => $chromo,
			'start'    => $position,
			'stop'     => $win_stop,
			'value'    => 'score',
			'absolute' => 1,
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
	my $window_scan_max = max(values %window_pos2score);
	
	if ( $window_scan_max >= $thresh) {
		# we have a peak that passes the threshold
		# yeah, a nucleosome
		return 1;
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
			
			if ($debug) {
				print DEBUG_FH "  Binned window scores:\n  ";
				foreach (sort {$a <=> $b} keys %binned_pos2score) {
					print DEBUG_FH "$_:$window_pos2score{$_}, ";
				}
			}
			
			# return success
			return 1;
		}
		else {
			# still no nucleosome
			return 0;
		}
	}
	else {
		# no peak
		return 0;
	}
}




sub identify_peak_position {
	# we know there is at least one peak within this window, now must find it
	
	# the window we're looking at
	my ($chromo, $position, $win_stop) = @_;
	
	# only process if we have a tag dataset
	if (!$tag_dataset) {
		# just use the position that was found from the scan dataset
		return $position;
	}
	
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
			'absolute' => 1,
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
	my $peak_position;
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
	
	# done
	return $peak_position;
}




## Process the nucleosome and record it
sub process_nucleosome {
	
	my ($chromo, $peak_position) = @_;
	
	# nucleosome statistics
	my $score;
	my $fuzziness;
	
	# we can only determine these statistics if we have a tag dataset
	if ($tag_dataset) {
		# Determine the fuzziness of this nucleosome
			# This is essentially the standard deviation of the counts
			# around the midpoint of the nuclesome. Well positioned 
			# nucleosomes will have all their midpoints centered in 
			# position, whereas fuzzy nucleosomes will have their 
			# midpoints spread out across the region. 
			
			# We will need to collect new scores surrounding the 
			# new found nucleosome peak.
			# arbitrarily collecting 40 bp worth of scores from each
			# side of the peak
			# we can't use the previously collected tag scores because
			# the nucleosome peak may have been at the edge of the window
		my %nucleosome_pos2score = get_region_dataset_hash( {
			'db'       => $db,
			'dataset'  => $tag_dataset,
			'chromo'   => $chromo,
			'start'    => $peak_position - 40,
			'stop'     => $peak_position + 40,
			'value'    => 'score',
			'absolute' => 1,
			# collecting values based on absolute, not relative, coord
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
				foreach (1 .. int(100 * $nucleosome_pos2score{$pos} )) {
					# for each position, we are going to add the 
					# relative position to the array of nucleosome_positions
					# for each tag count
					# to account for rpm counts (decimals) we will 
					# multiply by 100 - this shouldn't affect 
					# real number counts
					push @nucleosome_positions, $i;
				}
			}
		}
		
		# calculate fuzziness
		if (@nucleosome_positions) {
			$fuzziness = sprintf "%.0f", stddevp(@nucleosome_positions);
		}
		else {
			# nothing found!? not acceptable
			$fuzziness = '.';
		}
		
		# Calculate the number of nucleosome midpoints used in calling 
		# this nucleosome, this will be the score
		$score = get_chromo_region_score( {
			'db'       => $db,
			'dataset'  => $tag_dataset,
			'chromo'   => $chromo,
			'start'    => $peak_position - 37,
			'stop'     => $peak_position + 37,
			'value'    => 'score',
			'method'   => 'sum',
		} );
	}
	
	# no tag dataset
	else {
		$score = '.';
		$fuzziness = '.';
	}
	
	# Determine the nucleosome coordinates
		# we are assuming a standard 147 bp sized nucleosome
		# this may not be accurate for those partial and/or
		# fragile nucleosomes
	my $nuc_start = $peak_position - 73;
	my $nuc_stop  = $peak_position + 73;
	
	# Generate name
		# not quite following Pugh's example, we will name the 
		# nucleosome based on chromosome number and midpoint 
		# position (instead of start position)
	$chromo =~ /^(?:chr)?(.+)$/i;
	my $nuc_name = 'Nuc'. $1 . ':' . $peak_position;
	
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
	
	# calculate new position
	my $new_position = $nuc_stop + $buffer;
	
	if ($debug) {
		print DEBUG_FH " Nucleosome $nuc_name found at $nuc_start\..$nuc_stop\n";
		print DEBUG_FH " Advancing position to $new_position\n";
	}
	
	# Return the new position which is the current nucleosome endpoint plus buffer
	# the buffer helps to avoid too many overlapping nucleosomes
	return $new_position;
}



__END__

=head1 NAME

map_nucleosomes.pl

A script to map nucleosomes.

=head1 SYNOPSIS

map_nucleosomes.pl --db <text> --thresh <number> [--options...]
map_nucleosomes.pl --sdata <text|file> --thresh <number> [--options...]
  
  --db <text>
  --sdata <text|file>
  --tdata <text|file>
  --data <text|file>
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

Specify the name of the BioPerl database to pull the source data 
and/or chromsomes. A SeqFeature::Store database may be supplied, 
or a BigWigSet directory. Required unless data is pulled from 
a bigWig file (.bw).

=item --sdata <text|file>

=item --tdata <text|file>

=item --data <text|file>

Provide the name of the dataset(s) and/or data files (bigWig format)
containing the nucleosome midpoint occupancy data from which to
identify nucleosomal positions. If data is obtained from a database,
the type or primary_tag should be provided.

The scan dataset (--sdata) is required, while the tag dataset
(--tdata) is only used to collect nucleosome statistics. The same
dataset could be used for both, in which case both may be set simultaneously 
using the --data option. See the L<DESCRIPTION> for details regarding
the datasets. If the datasets are not specified on the commandline,
then they may be interactively chosen from a list from the database.

=item --thresh <number>

Provide the minimum score value required to call a nucleosome position. 
This is only used when scanning the scan dataset. 

=item --win <integer>

Provide the window size for which to scan the chromosome. Setting this  
value too large and overlapping nucleosomes may result, while setting 
this value too low may miss some nucleosomes. The default value is 175 bp. 

=item --buf <integer>

Provide the buffer size in bp which will be added between the end of the 
previous found nucleosome and the beginning of the window to scan for the 
next nucleosome. Setting this value may limit the number of overlapping 
nucleosomes. Default is 20 bp.

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
specified size (default 175 bp, set with --win option) looking for the
maximum peak that exceeds the minimum threshold value. The position of
the maximum peak is called as the new nucleosome midpoint. The window
is then advanced starting at the previous just-identified nucleosome
endpoint, or at the end of the previous window if no nucleosome was
identified. This position may be further advanced by setting the
buffer value (default 20 bp), which inserts space between the previous
nucleosome end and the next window. By advancing the window relative
to the previously identified nucleosome, the program can adapt to
variable nucleosome spacing and (hopefully) avoid overlapping
nucleosome calls.

This approach works reasonably well if the data shows an inherent, 
periodic pattern of peaks. Noisy datasets derived from partially 
fragmented nucleosomes or low sequencing depth may require some 
statistical smoothing.

The default values (window and buffer) were determined empirically
using yeast nucleosomes and paired-end next generation sequencing.
Values should be optimized and adjusted empirically for new datasets
and confirmed in a genome browser.

Two datasets must be provided for the mapping, although the same could
be used for both functions. The datasets should represent sequenced
nucleosome fragment midpoints. For the most accurate mapping, the
midpoints should be mapped at 1 bp resolution, but lower resolutions
(5 or 10 bp) could also work.

The scan dataset (--sdata) is used to scan for a potential nucleosome.
The dataset may be a simple difference, normalized difference, ratio, 
P-value, Q-Value, False Discovery Rate, or raw or normalized tag 
counts. Probability scores should be converted using -10Log10(P) for
proper interpretation.

The optional tag dataset (--tdata) is used for calculating nucleosome
statistics, including the precise nucleosome midpoint peak, occupancy,
and fuzziness. The tag dataset must be enumerated counts of nucleosome
midpoints. Normalized counts (for example, Reads Per Million mapped or
RPM) or difference counts (experiment minus control) may be used, but
P-value scores should not.

When working with next generation sequencing tags from nucleosomes, tags 
may be enumerated using the BioToolBox script L<bam2wig.pl>.

For genomic regions with few nucleosome reads and no positions that
pass the threshold value, the window positions may optionally be
binned together to increase sensitivity (see the --bin option).

Two attributes of each identified nucleosome are calculated, Occupancy
and Fuzziness. Occupancy represents the sum of the tag counts that
support the nucleosome position; higher scores indicate a more highly
occupied nucleosome. Fuzziness indicates how well all of the scores
are aligned in register. The Fuzziness value is the standard deviation
of the population of scores from the peak; a high value indicates the
nucleosome has a fuzzy or variable position, whereas a low value
indicates the nucleosome is highly positioned. For both attributes,
the scores are counted in a window representing 1/2 of a nucleosome
length (74 bp) centered on the defined nucleosome midpoint. The size
of this window is arbitrary but reasonable.

The identified nucleosomes are artificially set to 147 bp in length,
centered at the maximum peak. A second script,
C<get_actual_nuc_sizes.pl>, can determine the actual nucleosome sizes
based on paired-end reads in a BAM file.

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
