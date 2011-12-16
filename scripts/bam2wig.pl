#!/usr/bin/perl

# This script will convert alignments from a Bam file into enumerated 
# point data in a wig format

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper qw(
	open_to_write_fh
);
use tim_data_helper qw(
	format_with_commas
);
use tim_db_helper::config;
eval {
	# check for bam support
	require tim_db_helper::bam;
	tim_db_helper::bam->import;
};
my $VERSION = '1.5.6';
	

print "\n This program will convert bam alignments to enumerated wig data\n";

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
	$position,
	$use_coverage,
	$splice,
	$paired,
	$shift,
	$strand,
	$min_mapq,
	$interpolate,
	$rpm,
	$log,
	$track,
	$bedgraph,
	$bigwig,
	$bwapp,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the solexa data file
	'out=s'     => \$outfile, # name of output file 
	'position=s'=> \$position, # define position
	'coverage!' => \$use_coverage, # calculate coverage
	'splice|split!'   => \$splice, # split splices
	'pe!'       => \$paired, # paired-end alignments
	'shift=i'   => \$shift, # shift coordinates 3'
	'strand=s'  => \$strand, # select specific strands
	'qual=i'    => \$min_mapq, # minimum mapping quality
	'inter|fix!'=> \$interpolate, # positions with no count
	'rpm!'      => \$rpm, # calculate reads per million
	'log=i'     => \$log, # transform count to log scale
	'track!'    => \$track, # write a track line in the wig file
	'bed!'      => \$bedgraph, # write a bedgraph rather than wig file
	'bw!'       => \$bigwig, # generate bigwig file
	'bwapp=s'   => \$bwapp, # utility to generate a bigwig file
	'gz!'       => \$gz, # compress text output
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
	print " Biotoolbox script bam2wig.pl, version $VERSION\n\n";
	exit;
}



### Check for requirements and set defaults
unless ($infile) {
	$infile = shift @ARGV or
		die " no input file! use --help for more information\n";
}

my ($use_start, $use_mid, $use_span);
if ($position) {
	if ($position eq 'start') {
		$use_start = 1;
	}
	elsif ($position eq 'mid') {
		$use_mid = 1;
	}
	elsif ($position eq 'span') {
		$use_span = 1;
	}
	elsif ($position eq 'coverage') { 
		# for backwards compatibility
		$use_coverage = 1;
	}
	else {
		die " unrecognized position value '$position'! see help\n";
	}
}
else {
	# default values
	if ($paired) {
		$use_mid = 1;
	}
	else {
		$use_start = 1;
	}
}

my ($forward, $reverse);
if ($strand) {
	if ($strand eq 'f') {
		$forward = 1;
	}
	elsif ($strand eq 'r') {
		$reverse = 1;
	}
	else {
		warn " using default both strands\n";
	}
}

$shift = 0 unless defined $shift;

if (defined $min_mapq) {
	die " quality score must be < 255!\n" if $min_mapq > 255;
}
else {
	$min_mapq = 0;
}

unless (defined $splice) {
	$splice = 0;
}

if ($paired and $splice) {
	# kind of redundant to have spliced reads with paired-end alignments
	# plus this poses problems with the bam adaptor - getting substr errors
	# in Bio::DB::Bam::AlignWrapper
	warn " disabling splices with paired-end reads\n";
	$splice = 0;
}

unless ($outfile) {
	$outfile = $infile;
	$outfile =~ s/\.bam$//;
}

if (defined $gz) {
	# overide to false if bigwig is true
	$gz = 0 if $bigwig;
} 
else {
	# default is to use compression unless a bigwig file is requested
	# then the file is only temporary anyway
	$gz = $bigwig ? 0 : 1;
}

if ($bigwig) {
	# we need to set some options prior to writing the wig file if 
	# we're going to be writing a bigWig later
	
	# check for the app
	unless ($bwapp) {
		# check the environment path for Kent's conversion utilities
		# we prefer over wig over bedgraph, even though it has higher 
		# memory requirements, the file size is smaller but this can 
		# easily overcome using the bedgraph option, 
		if ($bedgraph) {
			$bwapp = 
				$TIM_CONFIG->param('applications.bedGraphToBigWig') || undef;
			
			# try looking in the environment path using external which 
			unless ($bwapp) {
				$bwapp = `which bedGraphToBigWig`;
				chomp $bwapp;
			}
			
			# can't find anything
			unless ($bwapp) {
				die " unable to find bedGraphToBigWig utility! see help\n";
			}
		}
		else {
			$bwapp = 
				$TIM_CONFIG->param('applications.wigToBigWig') || undef;
			
			# try looking in the environment path using external which 
			unless ($bwapp) {
				$bwapp = `which wigToBigWig`;
				chomp $bwapp;
			}
			
			# can't find anything
			unless ($bwapp) {
				die " unable to find wigToBigWig utility! see help\n";
			}
		}	
		
		# check executable
		if ($bwapp =~ /ToBigWig$/) {
			# looks good
			# remove newline just in case we used the which command
			chomp $bwapp;
			# make sure we don't write a track
			$track = 0;
		}
		else {
			warn " Unable to find bigWig conversion utility!\n" .
				" Generating wig file only\n";
			$bigwig = 0;
		}
	}
}

if ($log) {
	unless ($log == 2 or $log == 10) {
		die " requested log base '$log' not supported! Use 2 or 10\n";
	}
}



### Open files
# Bam file
unless (exists &open_bam_db) {
	die " unable to load Bam file support! Is Bio::DB::Sam installed?\n"; 
}
my $sam = open_bam_db($infile) or die " unable to open bam file '$infile'!\n";
$sam->split_splices($splice);

# low level bam and index
my $bam = $sam->bam;
my $index = $sam->bam_index;

# output file
my $outfh = open_wig_file();



### Calculate Total read numbers
my $total_read_number = 0;
if ($rpm) {
	# this is only required when calculating reads per million
	print " Calculating total number of aligned fragments....\n";
	$total_read_number = sum_total_bam_alignments($sam, $min_mapq, $paired);
	print "   ", format_with_commas($total_read_number), " mapped fragments\n";
}



### Process bam files
# global hash for storing current chromosome variables
my %data; 

# process according to type of data collected and alignment type
if ($use_coverage) {
	# special, speedy, low-level, single-bp coverage 
	process_bam_coverage();
}
elsif ($splice) {
	# single end alignments with splices require special callback
	process_alignments( \&single_end_spliced_callback );
}
elsif ($paired) {
	# paired end alignments require special callback
	process_alignments( \&paired_end_callback );
}
else {
	# single end alignments
	process_alignments( \&single_end_callback );
}

# finish
$outfh->close;
print " wrote wig file '$outfile'\n";

# Convert to BigWig if requested
convert_to_bigwig() if $bigwig;

print " Finished\n";





########################   Subroutines   ###################################

### Open the output file handle 
sub open_wig_file {
	# check extensions
	unless ($outfile =~ /\.wig$/) {
		$outfile .= '.wig';
	}
	if ($gz and $outfile !~ /\.gz$/) {
		$outfile .= '.gz';
	}
	
	# open
	my $fh = open_to_write_fh($outfile, $gz) or 
		die " unable to open output wig file '$outfile'!\n";
	
	# write track line
	if ($track) {
		if ($bedgraph) {
			$fh->print("track type=bedGraph\n");
		}
		else {
			$fh->print("track type=wiggle_0\n");
		}
	}
	
	return $fh;
}




### Collect alignment coverage
sub process_bam_coverage {
	# using the low level bam coverage method, not strand specific
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as 
		
		# get sequence info
		my $seq_length = $sam->target_len($tid);
		my $seq_id = $sam->target_name($tid);
		print " Converting reads on $seq_id...";
		
		# prepare definition line for fixedStep
		unless ($bedgraph) {
			$outfh->print("fixedStep chrom=$seq_id start=1 step=1 span=1\n");
		}
		
		# walk through the chromosome
		my $count = 0;
		my $dump = 50000;
		for (my $start = 0; $start < $seq_length; $start += $dump) {
			
			# the low level interface works with 0-base indexing
			my $end = $start + $dump -1;
			$end = $seq_length if $end > $seq_length;
			
			# using the low level interface for a little more performance
			# we're not binning, asking for 1 bp resolution
			my $coverage = $index->coverage(
				$bam,
				$tid,
				$start,
				$end
			);
			
			# now dump the coverage out to file
			# we're skipping the %data hash and write_wig() functions
			# to eke out more performance
			if ($bedgraph) {
				# we're writing a bedgraph file
				for (my $i = 0; $i < scalar(@{ $coverage }); $i++) {
					$outfh->print( 
						join("\t",
							$seq_id,
							$start + $i, 
							$start + $i + 1, 
							$coverage->[$i]
						) . "\n"
					);
					$count++;
				}
			}
			else {
				# we're writing a fixed step file
				for (my $i = 0; $i < scalar(@{ $coverage }); $i++) {
					$outfh->print($coverage->[$i] . "\n");
					$count++;
				}
			}
		}
		print "  ", format_with_commas($count), " positions were recorded\n";
	}
}


### Walk through the alignments on each chromosome
sub process_alignments {
	my $callback = shift;
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as 
		# a numeric target identifier
		# we can easily convert this to an actual sequence name
		# we will force the conversion to go one chromosome at a time
		
		# sequence name
		my $seq_id = $sam->target_name($tid);
		print " Converting reads on $seq_id...";
		
		# process the reads
		$sam->fetch($seq_id, $callback);
		
		# convert to reads per million if requested
		convert_to_rpm() if $rpm;
		
		# convert to log if requested
		convert_to_log() if $log;
		
		# write current chromo data to wig
		write_wig($seq_id, $sam->target_len($tid));
	}
}


### Callback for processing single-end alignments
sub single_end_callback {
	my $a = shift;
	my $checked = shift; # if true then no need to check alignment
						# only relevent when called from 
						# single_end_spliced_callback()
	
	# check alignment
	unless ($checked) {
		# subfeatures from split splices are not full AlignWrapper objects
		# so they can't be checked for alignment or mapping quality scores
		# skip this test in that case, as the parent was already checked
		
		# check if mapped and mapping quality
		return if $a->unmapped;
		return if $a->qual < $min_mapq;
	}
	
	# collect alignment data
	my $start  = $a->start;
	my $end    = $a->end;
	my $strand = $a->strand;
	return unless $start; # for some reason the alignment doesn't have a start?
	
	# check strand
	if ($forward or $reverse) {
		# stranded data is wanted
		
		# do nothing if strand is not what we want
		return if ($forward and $strand == -1);
		return if ($reverse and $strand == 1);
	}
	
	# shift 3' if requested
	if ($shift) {
		if ($strand == 1) {
			$start += $shift;
			$end += $shift;
		}
		else {
			$start -= $shift;
			$end -= $shift;
		}
	}
	
	# count this tag using the appropriate requested position
	if ($use_start) {
		# record at the 5' position
		if ($strand == 1) {	
			$data{$start} += 1;
		}
		else {
			$data{$end} += 1;
		}
	}
	elsif ($use_mid) {
		# calculate the midpoint position
		my $mid = int( ($start + $end) / 2);
		$data{$mid} += 1;
	}
	elsif ($use_span) {
		# we'll count every position along the alignment
		for (my $i = $start; $i <= $end; $i++) {
			$data{$i} += 1;
		}
	}
	return;
}


### Callback for processing single-end split alignments
sub single_end_spliced_callback {
	my $a = shift;
	
	# check alignment
	return if $a->unmapped;
	
	# check mapping quality
	return if $a->qual < $min_mapq;
	
	# check for subfeatures
	my @subfeatures = $a->get_SeqFeatures;
	if (@subfeatures) {
		# process each subfeature
		foreach my $subf (@subfeatures) {
			single_end_callback($subf, 1)
		}
	}
	else {
		# no subfeatures found
		# treat this as a single read
		single_end_callback($a, 1);
	}
	
	return;
}


### Callback for working with paired-end alignments
sub paired_end_callback {
	my $a = shift;
	
	# check alignment
	return if $a->unmapped;
	return unless $a->proper_pair;
	
	# we only need to process one of the two pairs, 
	# so only take the left (forward strand) read
	return unless $a->strand == 1;
	
	# check mapping quality
		# yes, this only checks the forward strand alignment, but using the 
		# fetch method (for performance reasons) doesn't allow checking both 
		# alignments simultaneously
		# it's a sacrifice
	return if $a->qual < $min_mapq;
	
	# check strand
	if ($forward or $reverse) {
		# stranded data is wanted
		# normally proper paired-end alignments are inherently not stranded
		# but in the case of paired-end RNA-Seq, we do want strand
		# The TopHat aligner records this information as a BAM record 
		# attribute under the flag XS
		my $strand = $a->get_tag_values('XS');
		
		# do nothing if strand is not what we want
		return if ($forward and $strand eq '-');
		return if ($reverse and $strand eq '+');
	}
	
	# collect alignment data
	my $start  = $a->start;
	my $isize  = $a->isize; # insert size
	
	# calculate end
		# I occasionally get errors if I call mate_end method
		# rather trust the reported insert size listed in the original bam file
	my $end = $start + $isize - 1;
	
	# count this tag using the appropriate requested position
	if ($use_start) {
		$data{$start} += 1;
	}
	elsif ($use_mid) {
		my $mid = int( ($start + $end) / 2);
		$data{$mid} += 1;
	}
	elsif ($use_span) {
		# we'll count every position along the alignment
		for (my $i = $start; $i <= $end; $i++) {
			$data{$i} += 1;
		}
	}
	return;
}


### Write the wig data for the current chromosome
sub write_wig {
	my $seq_id = shift;
	my $seq_length = shift;
	my $count = 0;
	
	# begin writing out the data
	if ($interpolate) {
		# we are interpolating the positions that don't have coverage and 
		# writing 0s
		
		if ($bedgraph) {
			# we're writing a bedgraph file
			for (my $i = 1; $i <= $seq_length; $i++) {
				$outfh->print( 
					join("\t",
						$seq_id,
						$i - 1, # bedgraphs are 0-based
						$i,
						$data{$i} ||= 0
					) . "\n"
				);
				$count++;
				delete $data{$i} if exists $data{$i};
			}
		}
		else {
			# we're writing a fixed step file
			$outfh->print("fixedStep chrom=$seq_id start=1 step=1 span=1\n");
			for (my $i = 1; $i <= $seq_length; $i++) {
				my $value = $data{$i} ||= 0;			
				$outfh->print("$value\n");
				$count++;
				delete $data{$i} if exists $data{$i};
			}
		}
	}
	
	else {
		# we are only writing the positions that have a tag 
		
		if ($bedgraph) {
			# we're writing a bedgraph file
			foreach my $i (sort {$a <=> $b} keys %data) {
				# only process properly positioned data
				if ($i > 0 and $i < $seq_length) {
					$outfh->print( 
						join("\t",
							$seq_id,
							$i - 1, # bedgraphs are 0-based
							$i,
							$data{$i}
						) . "\n"
					);
					$count++;
				}
				delete $data{$i};
			}
		}
		else {
			# we're writing a variable step file
			$outfh->print("variableStep chrom=$seq_id span=1\n");
			foreach my $i (sort {$a <=> $b} keys %data) {
				# only process properly positioned data
				if ($i > 0 and $i < $seq_length) {
					$outfh->print("$i\t$data{$i}\n");
					$count++;
				}
				delete $data{$i};
			}
		}
	}
	print "  ", format_with_commas($count), " positions were recorded\n";
	
	# empty the data hash for the next chromosome
	#%data = ();
}


### Run the BigWig conversion utility
sub convert_to_bigwig {
	print " Converting to bigWig...\n";
	
	# make new bw file name
	my $bw_file = $outfile;
	$bw_file =~ s/\.wig$/.bw/;
	
	# generate chromosome information file
	# we'll use the bam sequence header info for this
	my $chr_fh = open_to_write_fh('chromo.info');
	for my $tid (0 .. $sam->n_targets - 1) {
		$chr_fh->print( 
			join("\t", 
				$sam->target_name($tid),
				$sam->target_len($tid)
			) . "\n"
		);
	}
	$chr_fh->close;
	
	# run the utility, trapping errors in a file
	print " Running $bwapp...\n";
	system($bwapp, $outfile, 'chromo.info', $bw_file);
	
	# confirm
	if (-s $bw_file) {
		print " bigwig file '$bw_file' generated\n";
		unlink $outfile; # remove the wig file
	}
	else {
		warn " bigwig file not generated! see standard error\n";
	}
	unlink 'chromo.info';
}


### Convert to log scaling
sub convert_to_log {
	
	# converting all the values to log score
	if ($log == 2) {
		# log2 scale
		foreach my $i (keys %data) {
			next if $data{$i} == 0;
			$data{$i} = log($data{$i}) / log(2);
		}
	}
	
	elsif ($log == 10) {
		# log10 scale
		foreach my $i (keys %data) {
			next if $data{$i} == 0;
			$data{$i} = log($data{$i}) / log(10);
		}
	}
}


### Convert to reads per million
sub convert_to_rpm {
	foreach my $i (keys %data) {
		$data{$i} = ($data{$i} * 1000000) / $total_read_number;
	}
}



__END__

=head1 NAME

bam2wig.pl

=head1 SYNOPSIS

bam2wig.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --out <filename> 
  --position [start|mid|span]
  --coverage
  --splice|split
  --pe
  --shift <integer>
  --strand [f|r]
  --qual <integer>
  --inter|fix
  --rpm
  --log [2|10]
  --(no)track
  --bed
  --bw
  --bwapp </path/to/wigToBigWig or /path/to/bedGraphToBigWig>
  --(no)gz
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input Bam alignment file. The file should be sorted and 
indexed, although it may be indexed automatically

=item --out <filename>

Specify the output filename. By default it uses the base name of the 
input file.

=item --position [start|mid|span]

Specify the position of the alignment coordinate which should be 
recorded. Several positions are accepted: the start (5') position of 
the alignment, the midpoint of the alignment, or at all positions 
along the length of the alignment (span). Note that the  
span option gives coverage but not a true count of the number of 
alignments, unlike start or mid. With paired-end alignments, the 
positions are relative to the entire insert fragment defined by two 
alignments. The default value is start for single-end and mid for 
paired-end alignments.

=item --coverage

Calculate the coverage of the alignments over the genome at single 
base pair resolution. This ignores the position, quality, strand, shift, 
and log options. It uses faster low level interfaces to the Bam file to 
eke out performance. It is equivalent to specifying --position=span, 
--inter, minimum quality of 0, no strand, no rpm, and no log.

=item --splice

=item --split

The Bam file alignments may contain splices, where the 
read is split between two separate alignments. This is most common 
with splice junctions from RNA-Seq data. In this case, treat each 
alignment as a separate tag. This only works with single-end alignments. 
Paired-end spliced alignments are currently treated as single-end 
spliced alignments.

=item --pe

The Bam file consists of paired-end alignments, and only properly 
mapped pairs of alignments will be considered. The default is to 
treat all alignments as single-end.

=item --shift <integer>

Shift the positions of all single-end alignments towards the 3' end by 
the indicated number of basepairs. The value should be 1/2 the average 
length of the insert library sequenced. Useful for ChIP-Seq applications. 
Positions outside the chromosome length are not recorded.

=item --strand [f|r]

Only process those single-end alignments which map to the indicated 
strand. For paired-end RNA-Seq alignments that were generated with 
TopHat, the XS attribute is honored for strand information. 
The default is to take all alignments regardless of strand.

=item --qual <integer>

Set a minimum mapping quality score of alignments to count. The mapping
quality score is a posterior probability that the alignment was mapped
incorrectly, and reported as a -10Log10(P) value, rounded to the nearest
integer (range 0..255). Higher numbers are more stringent. For performance
reasons, when counting paired-end reads, only the left alignment is
checked. The default value is 0 (accept everything).

=item --inter

=item --fix

Specify whether or not to record interpolating positions with count of 0. 
If true, a fixedStep wig file (step=1 span=1) is written, otherwise a 
variableStep wig file is written that only records the positions 
where a tag is found. This will also work with bedGraph output. 
The default behavior is to not record empty positions.

=item --rpm

Convert the data to Reads (or Fragments) Per Million mapped. This is useful 
for comparing read coverage between different datasets. This conversion 
is applied before converting to log, if requested. The default is no 
conversion.

=item --log [2|10]

Transform the count to a log scale. Specify the base number, 2 or 
10. Note that positions with a count of 0 are not converted and 
remain 0, and positions with a count of 1 are converted to a 
log value of 0. Therefore, only really useful with Bam alignment 
files with high count numbers. Default is to not transform the count.

=item --(no)track

Specify whether or not to include a track line in the wig file. In 
general, track lines are not required when further converting to a 
BigWig file.

=item --bed

Specify whether or not to write a bedGraph (chromosome start stop value) 
file or a traditional fixedStep or variableStep wiggle file. The 
default is false.

=item --bw

Specify whether or not the wig file should be further converted into 
an indexed, compressed, binary BigWig file. The default is false.

=item --bwapp </path/to/bedGraphToBigWig or /path/to/wigToBigWig>

Specify the full path to Jim Kent's BigWig conversion utility. Two 
different utilities may be used, bedGraphToBigWig or wigToBigWig, 
depending on the format of the wig file generated. The wigToBigWig 
is preferred only because smaller file sizes are produced, but with 
slightly higher memory requirements. The application paths may be 
set in the biotoolbox.cfg file.

=item --(no)gz

Specify whether (or not) the output file should be compressed with 
gzip. The default is compress the output unless a BigWig file is 
requested.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will enumerate aligned sequence tags and generate a wig, 
or optionally BigWig, file. Each position in the wig file records the 
number of alignments which map to that position. Alignments may be 
counted at the start (5') or midpoint positions, or optionally 
enumerated at every position across the alignment (resulting in a 
coverage map rather than an alignment enumeration). Further, alignments 
may be selected according to strand, and the position may be shifted 
towards the 3' direction (for ChIP-Seq applications).

Note that the memory consumed by the program is roughly proportional to 
the size of the chromosome, particularly for dense read coverage. 
The total number of alignments should not matter. 

Counting and processing the alignments is a lengthy process, and may 
require a long time for large numbers, especially for splices. A 
simple coverage map should be much faster to calculate, but loses the 
options of selecting strand or shifting coordinates.

Conversion to a bigWig file requires the installation of Jim Kent's 
bedGraphToBigWig or wigToBigWig utilities. Conversion from a 
bedGraph file is slightly more memory efficient, but at the expense of 
larger file sizes.

More information about wiggle files can be found at 
http://genome.ucsc.edu/goldenPath/help/wiggle.html, bedGraph at 
http://genome.ucsc.edu/goldenPath/help/bedgraph.html, and bigWig at 
http://genome.ucsc.edu/goldenPath/help/bigWig.html.

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

