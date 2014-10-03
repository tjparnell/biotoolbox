#!/usr/bin/perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long;
use Statistics::Lite qw(count min max range sum mean stddevp);
use Bio::ToolBox::data_helper qw(generate_tim_data_structure);
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	get_chromosome_list
	get_region_dataset_hash
	get_chromo_region_score
);
use Bio::ToolBox::file_helper qw(
	write_tim_data_file
	convert_and_write_to_gff_file
);
#use Data::Dumper;
my $VERSION = '1.17';

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
	$outfile,
	$threshold,
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
	'out=s'    => \$outfile, # output file name
	'thresh=f' => \$threshold, # the nucleosome signal to call a nuc
	'win=i'    => \$window, # size of the window to scan the genome
	'buf=i'    => \$buffer, # the buffer between previous nuc and next window
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
unless (defined $threshold) {
	die " You must define the threshold!\n Use --help for more information\n";
}

unless ($window) {
	# set default value of 145 bp
	print " Using default window size of 145 bp\n";
	$window = 145;
} 
unless (defined $buffer) {
	print " Using default buffer size of 5 bp\n";
	$buffer = 5;
}



# Establish db connection
my $db;
if ($database) {
	# a database was provided by the user
	$db = open_db_connection($database) or 
		die " unable to connect to database 'database'!\n";
}
else {
	if ($dataset) {
		# we might be able to use an input bigwig file
		
		unless ($dataset =~ m/\.bw$/i) {
			die " You must define a database\n Use --help for more information\n";
		}
		$db = open_db_connection($dataset) or 
			die " unable to connect to database 'database'!\n";
	}
	else {
		die " You must define a database\n Use --help for more information\n";
	}
}



# Identify the dataset to use
$dataset = verify_or_request_feature_types(
		'db'      => $database,
		'feature' => $dataset,
		'prompt'  => "Enter the dataset to use for scanning for nucleosome peaks ",
		'single'  => 1,
);
unless ($dataset) {
	die " A valid scan dataset must be provided!\n";
}




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
	$outfile = $dataset . '_nucleosome';
}
my $success = write_tim_data_file(
	'data'     => $nucleosomes_ref,
	'filename' => $outfile,
);
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
			$type = $dataset . '_nucleosome';
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
	$success = convert_and_write_to_gff_file(
		'data'     => $nucleosomes_ref,
		'filename' => $outfile,
		'version'  => 3,
		'score'    => 5,
		'name'     => 4,
		'type'     => $type,
		'source'   => $source,
		'tag'      => [6],
	);
	if ($success) {
		print " Wrote GFF file '$success'\n";
	}
	else {
		print " Unable to write GFF file!\n";
	}
}
# The End







############################# Subroutines #####################################



###### Initialize the primary output data hash
sub initialize_nucleosome_data {
	
	# generate the data structure based on what datasets were provided
	my $data = generate_tim_data_structure(
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
	
	# add additional metadata
	$data->{'db'}               = $database;
	$data->{4}{'window'}        = $window;
	$data->{4}{'buffer'}        = $buffer;
	$data->{4}{'threshold'}     = $threshold;
	$data->{5}{'dataset'}       = $dataset;
	$data->{6}{'dataset'}       = $dataset;
	
	return $data;
}





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
	my @chromosomes = get_chromosome_list($db, 1) or 
		die " unable to find sequences in database!\n";
		# force the database to discard unwanted chromosomes
	foreach (@chromosomes) {
		
		# we have chromosome name and length
		my ($chromo, $chr_length) = @{$_};
		
		# skip mitochrondrial chromosome
			# just in case it wasn't excluded by default above
		if ($chromo =~ /^chrm|chrmt|mt|mito/i) {next}
		
		# Limit to one chromosome if debug
		if ($debug) {
			if ($chromo eq $chromosomes[2]) {
				# limit to only the first chromosome
				# this should force only chr1 to be done
				last;
			}
			
			# artificially set chr_length to small number for debug purposes
			# $chr_length = 5000;
		}
		
		# Progress report
		print " Scanning chromosome $chromo....\n";
		
		# Move along the chromosome
		my $position = 1;
		while ($position < $chr_length) {
			
			# determine the window stop position
			my $win_stop = $position + $window - 1;
			$win_stop = $chr_length if $win_stop > $chr_length;
			
			# collect the scores
			my %pos2score = get_region_dataset_hash(
					'db'       => $db,
					'dataset'  => $dataset,
					'chromo'   => $chromo,
					'start'    => $position,
					'stop'     => $win_stop,
					'value'    => 'score',
					'absolute' => 1,
			);
			
			if ($debug) {
				print DEBUG_FH "### Window $chromo:$position..$win_stop\n";
				print DEBUG_FH "  Window scores:\n  ";
				foreach (sort {$a <=> $b} keys %pos2score) {
					print DEBUG_FH "$_:$pos2score{$_}, ";
				}
				print DEBUG_FH "\n";
			}
			
			# scan the window for a potential nucleosome peak
			if ( max(values %pos2score) >= $threshold ) {
				# we have a peak, acurately map it with tag data
				# identify the peak position
				my $peak_position = identify_peak_position(
					$chromo, $position, $win_stop, \%pos2score);
				
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




sub identify_peak_position {
	# we know there is at least one peak within this window, now must find it
	
	# the window we're looking at
	my ($chromo, $position, $win_stop, $pos2score) = @_;
	
	# get all of the positions that correspond to that peak value
	my $tag_peak_value = max(values %{$pos2score});
	my @peak_positions;
	foreach (keys %{$pos2score}) {
		if ($pos2score->{$_} == $tag_peak_value) {
			push @peak_positions, $_;
		}
	}
	if ($debug) {
		print DEBUG_FH "  Max peak of $tag_peak_value at @peak_positions\n";
	}
	# now identify the best peak position, ideally this is just one
	my $peak_position;
	if (scalar @peak_positions == 1) {
		# one peak makes this easy
		$peak_position = shift @peak_positions;
	}
	
	else {
		# more than one equal peaks
		# if the positions are really close together, take the 
		# closest position
		if (range(@peak_positions) <= 10) {
			# the range between the peaks is 10 bp or less
			# 10 bp is totally arbitrary but reasonable
			# so take the closest one
			$peak_position = min(@peak_positions);
		}
		
		else {
			# they are too far apart, separate nucleosomes?
			# we'll take the one closest to the middle of the range
			my %best;
			my $window_midpoint =  
						int( ( ($position + $win_stop) / 2) + 0.5);
			foreach (sort {$b <=> $a} @peak_positions) {
				# we're reversing the order to take the rightmost or 
				# highest position first
				# that way if the positions are equidistant and the 
				# key gets overwritten, we'll take the leftmost 
				# lowest position
				$best{ abs( $_ -  $window_midpoint ) } = $_;
			}
			
			# take the position closest to the middle
			$peak_position = $best{ min( keys %best ) }
		}
	}
	
	if ($debug) {
		print DEBUG_FH " Peak found at position $peak_position\n";
	}
	
	# run a sanity check to verify the peak position
	$peak_position = verify_peak_position($chromo, $position, $peak_position);
	
	# done
	return $peak_position;
}



## Verify the peak position is accurate
sub verify_peak_position {
	
	# This sub verifies that the peak position identified in the window scan
	# is accurate.
	# We will re-scan within the vicinity of the supposed peak position.
	# If it is real, then no other peaks should be in the vicinity.
	# Otherwise, we may have prematurely called a nucleosome that is not 
	# accurate, and this sub should fix that.
	# This should greatly reduce the incidence of offset or improperly 
	# mapped nucleosomes found with the script verify_nucleosome_mapping.pl.
	
	# coordinates
	my ($chromo, $scan_start, $peak_position) = @_;
	
	# collect the raw data for +/- 50 bp around the supposed peak position
	# this distance is arbitrarily about 2/3rd nucleosome size
	# tried 1/2 and still getting offset nucleosomes
	# limit this by not going further backwards than the original window scan 
	# start - do not want to overlap previous nucleosome
	my %pos2score = get_region_dataset_hash(
			'db'       => $db,
			'dataset'  => $dataset,
			'chromo'   => $chromo,
			'start'    => $peak_position - 50 < $scan_start ? $scan_start : 
							$peak_position - 50,
			'stop'     => $peak_position + 50,
			'value'    => 'score',
			'absolute' => 1,
	);

	
	# find the max peak
	my @peaks;
	my $max = max(values %pos2score);
	foreach my $pos (keys %pos2score) {
		if ($pos2score{$pos} == $max) {
			# one of the peaks, 
			# very likely there will only be one, but there may be more
			# record all the positions with this peak value
			push @peaks, $pos;
		}
	}
	
	# as in the identify_peak_position, for multiple identical observed peaks
	# we will take the one closest to the original peak, assuming that is best
	if (scalar @peaks == 1) {
		# only one, that's great
		if ($debug and $peaks[0] != $peak_position) {
			print DEBUG_FH "  Reset the peak position from $peak_position" .
				" to $peaks[0] via sanity check\n";
		}
		return shift @peaks;
	}
	else {
		# more than one
		# identify one closest to the original peak
		my %best;
		foreach (sort {$b <=> $a} @peaks) {
			# we're reversing the order to take the rightmost or 
			# highest position first
			# that way if the positions are equidistant and the 
			# key gets overwritten, we'll take the leftmost 
			# position
			$best{ abs( $_ -  $peak_position ) } = $_;
		}
		
		# take the position closest to the original
		my $new_peak = $best{ min( keys %best ) };
		
		if ($debug and $new_peak != $peak_position) {
			print DEBUG_FH "  Reset the peak position from $peak_position" .
				" to $new_peak via sanity check\n";
		}
		return $new_peak;
	}
}


## Process the nucleosome and record it
sub process_nucleosome {
	
	my ($chromo, $peak_position) = @_;
	
	# nucleosome statistics
	my $score;
	my $fuzziness;
	
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
		
		# we will skip moving average here if it was enabled
		# we want real data, not averaged
	my %pos2score = get_region_dataset_hash(
		'db'       => $db,
		'dataset'  => $dataset,
		'chromo'   => $chromo,
		'start'    => $peak_position - 40,
		'stop'     => $peak_position + 40,
		'value'    => 'score',
		'absolute' => 1,
	);
	
	my @nucleosome_positions;
	for (my $i = -37; $i <= 37; $i++) {
		# we will advance along the relative position around the 
		# nucleosome peak and record the number of nucleosome 
		# midpoints found at each relative position
		
		# The size of the relative window to look for fuzziness
		# is a little more than 1/2 of a nucleosome length (81 bp). This
		# is totally arbitrary, but reasonable.
		
		# convert relative to absolute
		my $pos = $peak_position + $i; 
		
		# collect the counts
		if (exists $pos2score{$pos} and $pos2score{$pos} > 0) {
			foreach (1 .. int(100 * $pos2score{$pos} )) {
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
	$score = get_chromo_region_score(
		'db'       => $db,
		'dataset'  => $dataset,
		'chromo'   => $chromo,
		'start'    => $peak_position - 37,
		'stop'     => $peak_position + 37,
		'value'    => 'score',
		'method'   => 'sum',
	);
	
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

map_nucleosomes.pl --data <text|file> --thresh <number> [--options...]
  
  Options:
  --db <text>
  --data <text|file>
  --thesh <number>
  --win <integer>       (145)
  --buf <integer>       (5)
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

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
or other indexed data file, e.g. Bam or bigWig file, from which chromosome 
length information may be obtained. For more information about using databases, 
see L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. 
Required unless data is pulled from a bigWig file (.bw).

=item --data <text|file>

Provide the name of the dataset and/or data files (bigWig format)
containing the nucleosome midpoint occupancy data from which to
identify nucleosomal positions. If data is obtained from a database,
the name or type should be provided.

=item --thresh <number>

Provide the minimum score value required to call a nucleosome position. 
This is only used when scanning the scan dataset. 

=item --win <integer>

Provide the window size for which to scan the chromosome. Setting this  
value too large and overlapping nucleosomes may result, while setting 
this value too low may miss some nucleosomes. The default value is 145 bp. 

=item --buf <integer>

Provide the buffer size in bp which will be added between the end of the 
previous found nucleosome and the beginning of the window to scan for the 
next nucleosome. Setting this value may limit the number of overlapping 
nucleosomes. Default is 5 bp.

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
specified size (default 145 bp, set with --win option) looking for the
maximum peak that exceeds the minimum threshold value. The position of
the maximum peak is called as the new nucleosome midpoint. The window
is then advanced starting at the previous just-identified nucleosome
endpoint, or at the end of the previous window if no nucleosome was
identified. This position may be further advanced by setting the
buffer value (default 5 bp), which inserts space between the previous
nucleosome end and the next window. By advancing the window relative
to the previously identified nucleosome, the program can adapt to
variable nucleosome spacing and (hopefully) avoid overlapping
nucleosome calls.

This approach works reasonably well if the data shows an inherent, 
periodic pattern of peaks. Noisy datasets derived from partially 
fragmented nucleosomes or low sequencing depth may require some 
statistical smoothing.

=head1 DATASETS

This programs expects to work with enumerated midpoint counts of 
nucleosomal sequence fragments at a single bp resolution. Typically, 
nucleosome fragments are sequenced using massively parallel sequencing 
technologies, and the mapped alignments are converted to predicted 
midpoint occupancy counts. The midpoints may be precisely mapped with 
paired-end sequencing, or estimated by 3' end shifting of single-end 
sequence alignments. See the BioToolBox script bam2wig.pl for one 
approach.

=head1 REPORTING

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
based on paired-end reads in a BAM file. Because of stochastic 
positioning of nucleosomes in vivo, it is common to see 10-20% of 
the mapped nucleosomes exhibit predicted overlap with its neighbors. 

=head1 PARAMETERS

The default values (window and buffer) were determined empirically
using yeast nucleosomes and paired-end next generation sequencing.  
Parameters tested included windows of 140 to 180 bp in 5 bp increments, 
and buffers of 0 to 20 bp in 5 bp increments. In general, in order of 
importance, lower threshold, lower buffer sizes (with the exception of 
a buffer of 0 bp), and lower window sizes, lead to higher numbers of 
nucleosomes identified. The percentages of predicted nucleosome overlap 
and offcenter nucleosomes (where the observed tag dataset peak does not 
correspond to the reported nucleosome midpoint) generally increase as 
the number of nucleosomes are identified.

Values should be optimized and adjusted empirically for new datasets
and confirmed in a genome browser. 

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
