#!/usr/bin/perl
$| = 1;

# This script will convert alignments from a Bam file into enumerated 
# point data in a wig format

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(sum min max mean stddev);
use Statistics::LineFit;
use Statistics::Descriptive;
use Data::Dumper;
use FindBin qw($Bin);
#use lib "$Bin/../lib";
use lib '/usr/local/biotoolbox/lib';
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
my $VERSION = '1.7.2';
	

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
	$shift_value,
	$sample_number,
	$strand,
	$bin_size,
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
	'shift!'    => \$shift, # shift coordinates 3'
	'shiftval=i' => \$shift_value, # value to shift coordinates
	'sample=i'  => \$sample_number, # number of samples to test for shift
	'strand!'   => \$strand, # separate strands
	'bin=i'     => \$bin_size, # size of bin to make
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
# global variables
my ($use_start, $use_mid, $use_span, $bin, $record_count);

check_defaults();

# global hash for storing current chromosome counts
my %data1; # for either forward or combined strand
my %data2; # for reverse strand only

# record start time
my $start_time = time;





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

# output files
my ($filenames, $fh1, $fh2) = open_wig_files();




### Calculate Total read numbers
my $total_read_number = 0;
if ($rpm) {
	# this is only required when calculating reads per million
	print " Calculating total number of aligned fragments....\n";
	$total_read_number = sum_total_bam_alignments($sam, $min_mapq, $paired);
	print "   ", format_with_commas($total_read_number), " total mapped fragments\n";
	printf " counted in %.1f minutes\n", (time - $start_time)/60;
}



### Calculate shift value
if ($shift and !$shift_value) {
	print " Calculating 3' shift value...\n";
	$shift_value = determine_shift_value();
	print " Reads will be shifted by $shift_value bp\n";
}



### Process bam file
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




### Finish
$fh1->close;
if (defined $fh2) {
	$fh2->close;
}
foreach (@{ $filenames }) {
	print " wrote wig file '$_'\n";
}

# Convert to BigWig if requested
convert_to_bigwig() if $bigwig;

printf " Finished in %.1f min\n", (time - $start_time)/60;




########################   Subroutines   ###################################

### check required command line options and assign default values
sub check_defaults {
	# checking default and required values from the command line options
	# moved here to make it cleaner
	
	# check input file
	unless ($infile) {
		$infile = shift @ARGV or
			die " no input file! use --help for more information\n";
	}
	unless ($infile =~ /\.bam$/i) {
		die " must provide a .bam file as input!\n";
	}
	
	# check log number
	if ($log) {
		unless ($log == 2 or $log == 10) {
			die " requested log base '$log' not supported! Use 2 or 10\n";
		}
	}
	
	# check mapping quality
	if (defined $min_mapq) {
		die " quality score must be < 255!\n" if $min_mapq > 255;
	}
	else {
		$min_mapq = 0;
	}
	
	# check position
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
	if ($paired and $use_start) {
		warn " using midpoint with paired-end reads\n";
		$use_start = 0;
		$use_mid   = 1;
	}
	
	# check splices
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
	
	# check to shift position or not
	unless (defined $shift) {
		$shift = 0;
	}
	if ($shift and ($paired or $splice or $use_span) ) {
		warn " disabling shift with paired reads\n" if $paired;
		warn " disabling shift with splices enabled\n" if $splice;
		warn " disabling shift with span enabled\n" if $use_span;
		$shift = 0;
	}
	unless ($sample_number) {
		$sample_number = 100;
	}
	
	# check bin size
	if ($bin_size) {
		# set the boolean variable as to whether we're binning above 
		# the default 1 bp resolution 
		$bin = $bin_size > 1 ? 1 : 0;
	}
	else {
		# use the default 1 bp bin
		$bin_size = 1;
		$bin = 0;
	}
	if ($bin) {
		# binned data automatically assumes a fixed step wig
		$interpolate = 1;
	}
	
	# check output file
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
	
	# look for bigwig app
	if ($bigwig) {
		# we need to set some options prior to writing the wig file if 
		# we're going to be writing a bigWig later
		find_bigwig_app();
	}
	
	# determine the method used to record the read count
	# based on the position used, strandedness, and/or shift
	if ($use_coverage) {
		print " recording coverage spanning alignments\n";
	}
	elsif ($use_start and $strand and $shift) {
		$record_count = \&record_stranded_shifted_start;
		print " recording stranded, shifted-start positions\n";
	}
	elsif ($use_start and $strand and !$shift) {
		$record_count = \&record_stranded_start;
		print " recording stranded, start positions\n";
	}
	elsif ($use_start and !$strand and $shift) {
		$record_count = \&record_shifted_start;
		print " recording shifted-start positions\n";
	}
	elsif ($use_start and !$strand and !$shift) {
		$record_count = \&record_start;
		print " recording start positions\n";
	}
	elsif ($use_mid and $strand and $shift and $paired) {
		# this should not happen, shift is disabled with paired
		die " programming error!\n";
	}
	elsif ($use_mid and $strand and $shift and !$paired) {
		$record_count = \&record_stranded_shifted_mid;
		print " recording stranded, shifted-mid positions\n";
	}
	elsif ($use_mid and $strand and !$shift and $paired) {
		$record_count = \&record_stranded_paired_mid;
		print " recording stranded, mid positions of pairs\n";
	}
	elsif ($use_mid and $strand and !$shift and !$paired) {
		$record_count = \&record_stranded_mid;
		print " recording stranded, mid positions\n";
	}
	elsif ($use_mid and !$strand and $shift and $paired) {
		# this should not happen, shift is disabled with paired
		die " programming error!\n";
	}
	elsif ($use_mid and !$strand and $shift and !$paired) {
		$record_count = \&record_shifted_mid;
		print " recording shifted-mid positions\n";
	}
	elsif ($use_mid and !$strand and !$shift and $paired) {
		$record_count = \&record_paired_mid;
		print " recording mid positions of pairs\n";
	}
	elsif ($use_mid and !$strand and !$shift and !$paired) {
		$record_count = \&record_mid;
		print " recording mid position\n";
	}
	elsif ($use_span and $strand and $paired) {
		$record_count = \&record_stranded_paired_span;
		print " recording stranded positions spanning paired alignments\n";
	}
	elsif ($use_span and $strand and !$paired) {
		$record_count = \&record_stranded_span;
		print " recording stranded positions spanning alignments\n";
	}
	elsif ($use_span and !$strand and $paired) {
		$record_count = \&record_paired_span;
		print " recording positions spanning paired alignments\n";
	}
	elsif ($use_span and !$strand and !$paired) {
		$record_count = \&record_span;
		print " recording positions spanning alignments\n";
	}
	else {
		# what else is left!?
		die " programming error!\n";
	}
}


### Find the appropriate executable for generating a bigWig
sub find_bigwig_app {
	# check for the app
	unless ($bwapp) {
		# first check the biotoolbox configuration file for the UCSC 
		# conversion utilities
		# if it's not listed, then check the environment path 
		# we prefer wig over bedgraph, even though it has higher 
		# memory requirements, the file size is smaller
		if ($bedgraph) {
			$bwapp = 
				$TIM_CONFIG->param('applications.bedGraphToBigWig') || undef;
			
			# try looking in the environment path using external which 
			unless ($bwapp) {
				$bwapp = `which bedGraphToBigWig`;
				chomp $bwapp;
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
		}	
	}	
	
	# check the executable
	if ($bwapp =~ /ToBigWig$/) {
		# looks good
		# make sure we don't write a track
		$track = 0;
	}
	else {
		warn " Unable to find bigWig conversion utility!\n" .
			" Generating wig file only\n";
		$bigwig = 0;
		$track  = 0; # be nice, don't write a track
	}
}


### Open the output file handle 
sub open_wig_files {
	
	# prepare output file names
	my @names;
	$outfile =~ s/\.wig(?:\.gz)?$//i; # strip extension if present
	$outfile =~ s/\.bedgraph(?:\.gz)?$//i; # strip extension if present
	my $ext = $bedgraph ? '.bedgraph' : '.wig';
	if ($strand) {
		# we need to two names
		push @names, $outfile . '_f' . $ext;
		push @names, $outfile . '_r' . $ext;
	}
	else {
		# just one name
		push @names, $outfile . $ext;
		
	}
	
	# check compression
	if ($gz) {
		@names = map {$_ .= '.gz'} @names;
	}
	
	# open
	my @fhs;
	foreach my $name (@names) {
		my $fh = open_to_write_fh($name, $gz) or 
			die " unable to open output wig file '$name'!\n";
		
		# write track line
		if ($track) {
			if ($bedgraph) {
				$fh->print("track type=bedGraph\n");
			}
			else {
				$fh->print("track type=wiggle_0\n");
			}
		}
		
		push @fhs, $fh;
	}
	
	# finished, return ref to names array, and 1 or 2 filehandles
	return (\@names, @fhs);
}


### Determine the shift value
sub determine_shift_value {
	
	# remember original values
	my $original_record_count = $record_count;
	$record_count = \&record_stranded_start;
	my $original_bin = $bin;
	$bin = 1;
	my $original_bin_size = $bin_size;
	$bin_size = 10;
	
	# find the biggest chromosome to sample
	my %chrom2size;
	for my $tid (0 .. $sam->n_targets - 1) {
		# key is chromosome name, value is its size
		$chrom2size{$tid} = $sam->target_len($tid);
	}
	my $biggest_id;
	my $biggest_size = 1;
	foreach (keys %chrom2size) {
		if ($chrom2size{$_} > $biggest_size) {
			$biggest_id = $_;
			$biggest_size = $chrom2size{$_};
		}
	}
	my $biggest_chrom = $sam->target_name($biggest_id);
	
	# identify top regions to score
	# we will walk through the largest chromosome looking for the top  
	# 1 kb regions containing the highest unstranded coverage to use
	my %coverage2region;
	
	# look for test regions
	print "  searching for top $sample_number coverage regions on $biggest_chrom to sample... ";
	for (my $start = 1; $start < $biggest_size; $start += 500) {
		
		# the low level interface works with 0-base indexing
		my $end = $start + 499;
		$end = $biggest_size if $end > $biggest_size;
		
		# using the low level interface for a little more performance
		my $coverage = $index->coverage(
			$bam,
			$biggest_id,
			$start - 1,
			$end,
		);
		my $sum_coverage = sum( @{$coverage} );
		next if $sum_coverage == 0;
		
		# check if our coverage exceeds the lowest region
		if (scalar keys %coverage2region < $sample_number) {
			# less than requested regions found so far, so keep it
			# record the start and stop position for this region
			$coverage2region{$sum_coverage} = [$start, $end];
		}
		else {
			# we already have the maximum number
			
			# find the lowest one
			my $current_min = min( keys %coverage2region );
			if ($sum_coverage > $current_min) {
				# it's a new high over the lowest minimum
				
				# remove the previous lowest region
				delete $coverage2region{ $current_min };
				
				# add the current region
				# record the start and end position for this region
				$coverage2region{$sum_coverage} = [$start, $end];
			}
		}
	}
	printf " in %.1f min\n", (time - $start_time)/60;
	
		
	# now determine the optimal shift for each of the test regions
	my @shift_values;
	foreach my $i (sort {$a <=> $b} keys %coverage2region) {
		
		# get the start and end positions
		# we're adjusting them by 500 bp in both directions to actually 
		# sample a 1.5 kb region centered over the original region
		my $start = $coverage2region{$i}->[0] - 500;
		my $end   = $coverage2region{$i}->[1] + 500;
		# just in case we go over
		$start = 1 if $start < 1;
		$end = $biggest_size if $end > $biggest_size;
		
		print "  sampling $biggest_chrom:$start..$end  ";
		
		# collect stranded data from our sample window
		$sam->fetch("$biggest_chrom:$start..$end", \&single_end_callback);
		unless (%data1 and %data2) {
			die " no stranded data collected for calculating shift value!\n";
		}
		
		# generate data arrays
		my @f;
		my @r;
		for (my $i = $start; $i <= $end; $i += 10) {
			if (exists $data1{$i}) {
				push @f, $data1{$i};
				delete $data1{$i};
			}
			else {
				push @f, 0;
			}
			if (exists $data2{$i}) {
				push @r, $data2{$i};
				delete $data2{$i};
			}
			else {
				push @r, 0;
			}
		}
		
		# calculate correlations
		my %r2shift;
		for (my $i = 1; $i <= 50; $i++) {
			# check shift from 10 to 500 bp
			
			# adjust the arrays, mimicking shifting arrays towards the 3'
			unshift @f, 0;
			pop @f;
			shift @r;
			push @r, 0;
			
			# skip the first 30 bp
			next if $i < 4;
			
			# calculate correlation
			my $stat = Statistics::LineFit->new();
			$stat->setData(\@r, \@f) or warn " bad data!\n";
			my $r2 = $stat->rSquared();
			
			# store rsquared
			if ($r2 > 0.2) {
				$r2shift{$r2} = $i * 10;
			}
		}
		
		# determine best shift
		if (%r2shift) {
			my $max_r = max(keys %r2shift);
			push @shift_values, $r2shift{$max_r};
			printf " shift $r2shift{$max_r} bp (r^2 %.3f)\n", $max_r;
		}
		else {
			print "\n";
		}
	}
	
	# determine the optimal shift value
	my $best_value = mean(@shift_values);
	printf "  The mean shift value is %.0f +/- %.0f bp\n", 
		$best_value, stddev(@shift_values);
	
	# restore original values
	$record_count = $original_record_count;
	$bin = $original_bin;
	$bin_size = $original_bin_size;
	
	# make sure data hashes are empty
	%data1 = ();
	%data2 = ();
	
	# done
	return sprintf("%.0f", $best_value);
}

### Collect alignment coverage
sub process_bam_coverage {
	# using the low level bam coverage method, not strand specific
	
	# determine the dump size
	# the dump size indicates how much of the genome we take before we 
	# dump to file
	# this is based on the requested bin_size, default should be 1 bp 
	# but we want to keep the coverage array reasonable size, 10000 elements
	# is plenty. we'll do some math to make it a multiple of bin_size to fit
	my $multiplier = int( 10000 / $bin_size);
	my $dump = $bin_size * $multiplier; 
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as an integer
		
		# get sequence info
		my $seq_length = $sam->target_len($tid);
		my $seq_id = $sam->target_name($tid);
		print " Converting reads on $seq_id...";
		
		# prepare definition line for fixedStep
		unless ($bedgraph) {
			$fh1->print(
				"fixedStep chrom=$seq_id start=1 step=$bin_size span=$bin_size\n"
			);
		}
		
		# walk through the chromosome
		my $count = 0;
		for (my $start = 0; $start < $seq_length; $start += $dump) {
			
			# the low level interface works with 0-base indexing
			my $end = $start + $dump -1;
			$end = $seq_length if $end > $seq_length;
			
			# using the low level interface for a little more performance
			my $coverage = $index->coverage(
				$bam,
				$tid,
				$start,
				$end,
			);
			
			# now dump the coverage out to file
			# we're skipping the %data hash and write_wig() functions
			# to eke out a little more performance
			if ($bedgraph) {
				# we're writing a bedgraph file
				
				for (my $i = 0; $i < scalar(@{ $coverage }); $i += $bin_size) {
					
					# sum the reads within our bin
					my $sum = 0;
					for (my $n = $i; $n < $i + $bin_size; $n++) {
						$sum += $coverage->[$n];
					}
					
					# print the bedgraph line
					$fh1->print( 
						join("\t",
							$seq_id,
							$start + $i, 
							$start + $i + $bin_size, 
							$sum
						) . "\n"
					);
					$count++;
				}
			}
			else {
				# we're writing a fixed step file
				
				for (my $i = 0; $i < scalar(@{ $coverage }); $i += $bin_size) {
					
					# sum the reads within our bin
					my $sum = 0;
					for (my $n = $i; $n < $i + $bin_size; $n++) {
						$sum += $coverage->[$n];
					}
					
					# print the wig line
					$fh1->print("$sum\n");
					$count++;
				}
			}
		}
		print "  ", format_with_commas($count), " positions were recorded";
		printf " in %.1f minutes\n", (time - $start_time)/60;
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
		
		# sequence info
		my $seq_id = $sam->target_name($tid);
		my $seq_length = $sam->target_len($tid);
		print " Converting $seq_id...";
		
		# process the reads in batches across the chromosome
		my $count = 0;
		my $start = 1;
		while ($start < $seq_length) {
			# we'll be processing in windows of 2 Mb, plus a little extra
			# to buffer for paired ends 
			# after we go throug the window, we'll write out the current data
			
			# determine end
			my $end = $start + 2001000; 
			if ($end > $seq_length) {
				$end = $seq_length;
			}
			
			$sam->fetch("$seq_id:$start..$end", $callback);
			
			# convert to reads per million if requested
			convert_to_rpm() if $rpm;
			
			# convert to log if requested
			convert_to_log() if $log;
			
			# write current window data to wig
			$count += write_wig(
				$seq_id, 
				$seq_length,
				$start,
				\%data1,
				$fh1
			);
			if (%data2) {
				$count += write_wig(
					$seq_id, 
					$seq_length,
					$start,
					\%data2,
					$fh2
				);
			}
			
			# ready for next segment
			$start += 2000000;
		}
		
		# print result
		$count = format_with_commas($count);
		if ($strand) {
			print " $count stranded positions recorded";
		}
		else {
			print " $count positions recorded";
		}
		printf " in %.1f minutes\n", (time - $start_time)/60;
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
	
	# record the alignment
	&{ $record_count }($a);
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
	
	# record the alignment
	&{ $record_count }($a);
}


### Write the wig data for the current chromosome
sub write_wig {
	my ($seq_id, $seq_length, $start, $data, $fh) = @_;
	my $count = 0;
	
	# check that we have data to write
	if (!%{ $data } and !$interpolate) {
		# nothing to write for a varStep file!
		# doesn't count for fixedStep
		return $count;
	}
	
	# set the window endpoint
	# add the window size 2 Mb I'm using, minus 1 kb buffer, minus 1 bp
	my $end = $start + 1999999; 
	if ($end > $seq_length) {
		if ($bin) {
			# need to compensate for the step and span parameter of the 
			# wig file and make it a multiple of bin_size
			# it's possible we may lose some small amount of information 
			# at the very end
			$end = $seq_length - ($seq_length % $bin_size) - 2 * $bin_size;
		}
		else {
			$end = $seq_length;
		}
	}
	
	# begin writing out the data
	if ($interpolate) {
		# we are interpolating the positions that don't have coverage and 
		# writing 0s
		
		if ($bedgraph) {
			# we're writing a bedgraph file
			for (my $i = $start; $i <= $end; $i += $bin_size) {
				$fh->print( 
					join("\t",
						$seq_id,
						$i - 1, # bedgraphs are 0-based
						$i + $bin_size - 1, # for the span
						$data->{$i} ||= 0
					) . "\n"
				);
				$count++;
				# delete $data->{$i} if exists $data->{$i};
			}
		}
		else {
			# we're writing a fixed step file
			
			# write a header line if we're at the beginning of a chromosome
			if ($start == 1) {
				$fh->print(
					"fixedStep chrom=$seq_id start=1 step=$bin_size span=$bin_size\n"
				);
			}
			for (my $i = $start; $i <= $end; $i += $bin_size) {
				my $value = $data->{$i} ||= 0;			
				$fh->print("$value\n");
				$count++;
				# delete $data->{$i} if exists $data->{$i};
			}
		}
	}
	
	else {
		# we are only writing the positions that have a tag 
		# these are 1 bp bins, any bins > 1 bp should automatically 
		# be interpolated and written above
		
		if ($bedgraph) {
			# we're writing a bedgraph file
			
			# walk through data
			foreach my $i (sort {$a <=> $b} keys %{ $data }) {
				# only process properly positioned data
				if ($i >= $start and $i <= $end) {
					$fh->print( 
						join("\t",
							$seq_id,
							$i - 1, # bedgraphs are 0-based
							$i,
							$data->{$i}
						) . "\n"
					);
					$count++;
				}
				# delete $data->{$i};
			}
		}
		else {
			# we're writing a variable step file
			
			# write a header line if we're at the beginning of a chromosome
			if ($start == 1) {
				$fh->print("variableStep chrom=$seq_id span=1\n");
			}
			
			# walk through data
			foreach my $i (sort {$a <=> $b} keys %{ $data }) {
				# only process properly positioned data
				if ($i >= $start and $i <= $end) {
					$fh->print("$i\t$data->{$i}\n");
					$count++;
				}
				# delete $data->{$i};
			}
		}
	}
	
	# done
	%{ $data } = ();
	return $count;
}


### Run the BigWig conversion utility
sub convert_to_bigwig {
	print " Converting to bigWig...\n";
	
	# generate chromosome information file
	# we'll use the bam sequence header info for this
	my $chr_file = "$outfile\.chromo.info";
	my $chr_fh = open_to_write_fh($chr_file);
	for my $tid (0 .. $sam->n_targets - 1) {
		$chr_fh->print( 
			join("\t", 
				$sam->target_name($tid),
				$sam->target_len($tid)
			) . "\n"
		);
	}
	$chr_fh->close;
	
	# process each wig file
	my $error = 0;
	foreach my $name (@{ $filenames }) {
	
		# make new bw file name
		my $bw_file = $name;
		$bw_file =~ s/\.wig$/.bw/;
		
		# run the utility, trapping errors in a file
		print " Running $bwapp...\n";
		system($bwapp, $name, $chr_file, $bw_file);
		
		# confirm
		if (-s $bw_file) {
			print " bigwig file '$bw_file' generated\n";
			unlink $name; # remove the wig file
		}
		else {
			warn " bigwig file not generated! see standard error\n";
			$error = 1;
		}
	}
	
	# remove chromosome file
	if ($error) {
		warn " leaving chromosome file $chr_file\n";
	}
	else {
		unlink $chr_file;
	}
}


### Convert to log scaling
sub convert_to_log {
	
	# converting all the values to log score
	if ($log == 2) {
		# log2 scale
		
		# process main data hash
		foreach my $i (keys %data1) {
			$data1{$i} = log( $data1{$i} + 1 ) / log(2);
		}
		
		# check second data hash
		foreach my $i (keys %data2) {
			$data2{$i} = log( $data2{$i} + 1 ) / log(2);
		}
	}
	
	elsif ($log == 10) {
		# log10 scale
		
		# process main data hash
		foreach my $i (keys %data1) {
			$data1{$i} = log( $data1{$i} + 1 ) / log(10);
		}
		
		# check second data hash
		foreach my $i (keys %data2) {
			$data2{$i} = log( $data2{$i} + 1 ) / log(10);
		}
	}
}


### Convert to reads per million
sub convert_to_rpm {
	
	# process main data hash
	foreach my $i (keys %data1) {
		$data1{$i} = ($data1{$i} * 1000000) / $total_read_number;
	}
	
	# process second data hash
	foreach my $i (keys %data2) {
		$data2{$i} = ($data2{$i} * 1000000) / $total_read_number;
	}
}


### Record stranded at shifted start position
sub record_stranded_shifted_start {
	my $a = shift;
	
	# record based on the strand
	if ($a->strand == 1) {
		# forward strand
		
		# calculate record position
		my $position = $a->start + $shift_value;
		
		# determine bin
		if ($bin) {
			$position = $position - ($position % $bin_size) + 1;
		}
		
		# record position in forward strand data hash
		$data1{$position} += 1;
	}
	else {
		# reverse strand
		
		# calculate record position
		my $position = $a->end - $shift_value;
		
		# determine bin
		if ($bin) {
			$position = $position - ($position % $bin_size) + 1;
		}
		
		# record position in reverse strand data hash
		$data2{$position} += 1;
	}
}


### Record stranded at start position
sub record_stranded_start {
	my $a = shift;
	
	# record based on the strand
	if ($a->strand == 1) {
		# forward strand
		
		# record position in forward strand data hash
		if ($bin) {
			# calculate bin
			my $position = $a->start - ($a->start % $bin_size) + 1;
			$data1{$position} += 1;
		}
		else {
			# no bin
			$data1{ $a->start } += 1;
		}
	}
	else {
		# reverse strand
		
		# record position in reverse strand data hash
		if ($bin) {
			# calculate bin
			my $position = $a->end - ($a->end % $bin_size) + 1;
			$data2{$position} += 1;
		}
		else {
			# no bin
			$data1{ $a->end } += 1;
		}
	}
}


### Record at shifted start position
sub record_shifted_start {
	my $a = shift;
	
	# record based on the strand
	my $position;
	if ($a->strand == 1) {
		# forward strand
		$position = $a->start + $shift_value;
	}
	else {
		# reverse strand
		$position = $a->end - $shift_value;
	}
		
	# determine bin
	if ($bin) {
		$position = $position - ($position % $bin_size) + 1;
	}
	
	# record position in main data hash
	$data1{$position} += 1;
}


### Record at start position
sub record_start {
	my $a = shift;
	
	# record based on the strand
	my $position;
	if ($a->strand == 1) {
		# forward strand
		$position = $a->start;
	}
	else {
		# reverse strand
		$position = $a->end;
	}
		
	# determine bin
	if ($bin) {
		$position = $position - ($position % $bin_size) + 1;
	}
	
	# record position in main data hash
	$data1{$position} += 1;
}


### Record stranded at shifted mid position
sub record_stranded_shifted_mid {
	my $a = shift;
	
	# record based on the strand
	if ($a->strand == 1) {
		# forward strand
		
		# calculate record position
		my $position = int( ($a->start + $a->end) / 2) + $shift_value;
		
		# determine bin
		if ($bin) {
			$position = $position - ($position % $bin_size) + 1;
		}
		
		# record position in forward strand data hash
		$data1{$position} += 1;
	}
	else {
		# reverse strand
		
		# calculate record position
		my $position = int( ($a->start + $a->end) / 2) - $shift_value;
		
		# determine bin
		if ($bin) {
			$position = $position - ($position % $bin_size) + 1;
		}
		
		# record position in reverse strand data hash
		$data2{$position} += 1;
	}
}


### Record stranded at mid position for paired end
sub record_stranded_paired_mid {
	my $a = shift;
	
	# calculate record position
	my $position = int( ($a->start + ($a->start + $a->isize -1) ) / 2);
	
	# determine bin
	if ($bin) {
		$position = $position - ($position % $bin_size) + 1;
	}
	
	# determine strand
	# normally proper paired-end alignments are inherently not stranded
	# but in the case of paired-end RNA-Seq, we do want strand
	# The TopHat aligner records this information as a BAM record 
	# attribute under the flag XS
	my $strand = $a->aux_get('XS') || '+';
		# default is going to be the forward strand
	
	# record based on the strand
	if ($strand eq '+') {
		# record position in forward strand data hash
		$data1{$position} += 1;
	}
	elsif ($strand eq '-') {
		# record position in reverse strand data hash
		$data2{$position} += 1;
	}
	else {
		die " unrecognized strand value '$strand' for XS attribute in BAM record\n";
	}
}


### Record stranded at mid position
sub record_stranded_mid {
	my $a = shift;
	
	# calculate record position
	my $position = int( ($a->start + $a->end) / 2);
	
	# determine bin
	if ($bin) {
		$position = $position - ($position % $bin_size) + 1;
	}
	
	# record based on the strand
	if ($a->strand == 1) {
		# record position in forward strand data hash
		$data1{$position} += 1;
	}
	else {
		# record position in reverse strand data hash
		$data2{$position} += 1;
	}
}


### Record at shifted mid position
sub record_shifted_mid {
	my $a = shift;
	
	# calculate shifted position based on strand
	my $position;
	if ($a->strand == 1) {
		# forward strand
		$position = int( ($a->start + $a->end) / 2) + $shift_value;
		
	}
	else {
		# reverse strand
		$position = int( ($a->start + $a->end) / 2) - $shift_value;
	}
	
	# determine bin
	if ($bin) {
		$position = $position - ($position % $bin_size) + 1;
	}
	
	# record position in main data hash
	$data1{$position} += 1;
}


### Record at mid position for paired end
sub record_paired_mid {
	my $a = shift;
	
	# calculate record position
	my $position = int( ($a->start + ($a->start + $a->isize -1) ) / 2);
	
	# determine bin
	if ($bin) {
		$position = $position - ($position % $bin_size) + 1;
	}
	
	# record 
	$data1{$position} += 1;
}


### Record at mid position
sub record_mid {
	my $a = shift;
	
	# calculate record position
	my $position = int( ($a->start + $a->end) / 2);
	
	# determine bin
	if ($bin) {
		$position = $position - ($position % $bin_size) + 1;
	}
	
	# record 
	$data1{$position} += 1;
}


### Record stranded across alignment for paired end
sub record_stranded_paired_span {
	my $a = shift;
	
	# determine the end 
	my $end = $a->start + $a->isize - 1;
	
	# determine strand
	# normally proper paired-end alignments are inherently not stranded
	# but in the case of paired-end RNA-Seq, we do want strand
	# The TopHat aligner records this information as a BAM record 
	# attribute under the flag XS
	my $strand = $a->aux_get('XS') || '+';
		# default is going to be the forward strand
	
	# we'll count every position along the alignment
	if ($strand eq '+') {
		# forward strand
		
		if ($bin) {
			# bin the counts
			
			# walk along every position of the alignment
			for (my $i = $a->start; $i <= $end; $i++) {
				
				# determin the bin
				my $position = $i - ($i % $bin_size) + 1;
				
				# record 
				$data1{$position} += 1;
			}
		}
		else {
			# no need to bin the positions
			
			# walk along every position of the alignment
			for (my $i = $a->start; $i <= $end; $i++) {
				$data1{$i} += 1;
			}
		}
	}
	elsif ($strand eq '-') {
		# reverse strand
		
		if ($bin) {
			# bin the counts
			
			# walk along every position of the alignment
			for (my $i = $a->start; $i <= $a->end; $i++) {
				
				# determin the bin
				my $position = $i - ($i % $bin_size) + 1;
				
				# record 
				$data2{$position} += 1;
			}
		}
		else {
			# no need to bin the positions
			
			# walk along every position of the alignment
			for (my $i = $a->start; $i <= $a->end; $i++) {
				$data2{$i} += 1;
			}
		}
	}
	else {
		die " unrecognized strand value '$strand' for XS attribute in BAM record\n";
	}
}


### Record stranded across alignment
sub record_stranded_span {
	my $a = shift;
	
	# we'll count every position along the alignment
	if ($a->strand == 1) {
		# forward strand
		
		if ($bin) {
			# bin the counts
			
			# walk along every position of the alignment
			for (my $i = $a->start; $i <= $a->end; $i++) {
				
				# determin the bin
				my $position = $i - ($i % $bin_size) + 1;
				
				# record 
				$data1{$position} += 1;
			}
		}
		else {
			# no need to bin the positions
			
			# walk along every position of the alignment
			for (my $i = $a->start; $i <= $a->end; $i++) {
				$data1{$i} += 1;
			}
		}
	}
	else {
		# reverse strand
		
		if ($bin) {
			# bin the counts
			
			# walk along every position of the alignment
			for (my $i = $a->start; $i <= $a->end; $i++) {
				
				# determin the bin
				my $position = $i - ($i % $bin_size) + 1;
				
				# record 
				$data2{$position} += 1;
			}
		}
		else {
			# no need to bin the positions
			
			# walk along every position of the alignment
			for (my $i = $a->start; $i <= $a->end; $i++) {
				$data2{$i} += 1;
			}
		}
	}
}


### Record across alignment for paired end
sub record_paired_span {
	my $a = shift;
	
	# we'll count every position between the alignments
	my $end = $a->start + $a->isize - 1;
	if ($bin) {
		# bin the counts
		
		# walk along every position of the alignment
		for (my $i = $a->start; $i <= $end; $i++) {
			
			# determine the bin
			$position = $i - ($i % $bin_size) + 1;
			
			# record 
			$data1{$position} += 1;
		}
	}
	else {
		# no need to bin the positions
		
		# walk along every position of the alignment
		for (my $i = $a->start; $i <= $end; $i++) {
			$data1{$i} += 1;
		}
	}
}


### Record across alignment
sub record_span {
	my $a = shift;
	
	# we'll count every position along the alignment
	if ($bin) {
		# bin the counts
		
		# walk along every position of the alignment
		for (my $i = $a->start; $i <= $a->end; $i++) {
			
			# determin the bin
			my $position = $i - ($i % $bin_size) + 1;
			
			# record 
			$data1{$position} += 1;
		}
	}
	else {
		# no need to bin the positions
		
		# walk along every position of the alignment
		for (my $i = $a->start; $i <= $a->end; $i++) {
			$data1{$i} += 1;
		}
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
  --bin <integer>
  --shift
  --shiftval <integer>
  --sample <integer>
  --strand
  --qual <integer>
  --inter|fix
  --bed
  --rpm
  --log [2|10]
  --(no)track
  --bw
  --bwapp </path/to/wigToBigWig or /path/to/bedGraphToBigWig>
  --(no)gz
  --version
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
alignments. The default value is the start position for single-end and 
midpoint position for paired-end alignments.

=item --coverage

Quickly calculate the coverage of the alignments over the genome, 
either at single bp resolution (default) or in bins. This method ignores 
the position, quality, strand, shift, and log options. It uses a fast 
low-level interface to the Bam file to eke out performance. It is 
equivalent to specifying --position=span, --fix, --split, 
--nope, --noshift, --nostrand, --qual=0, --norpm, and no log. 

=item --splice

=item --split

The Bam file alignments may contain splices, where the 
read is split between two separate alignments. This is most common 
with splice junctions from RNA-Seq data. In this case, treat each 
alignment as a separate tag. This only works with single-end alignments. 
Splices are disabled for paired-end reads.

=item --pe

The Bam file consists of paired-end alignments, and only properly 
mapped pairs of alignments will be counted. The default is to 
treat all alignments as single-end.

=item --bin <integer>

Specify the window or bin size in which alignments are counted. This 
may make visualization in browsers easier when zoomed out to large 
views. This option is compatible with all modes of counting, including 
coverage. The default is to count at single basepair resolution. 

=item --shift

Specify that the positions of the alignment should be shifted towards 
the 3' end. Useful for ChIP-Seq applications, where only the ends of 
the fragments are counted and often seen as discrete peaks on separate 
strands flanking the true target site. This option is disabled with 
paired-end and spliced reads, and when --position=span is enabled. 

=item --shiftval <integer>

Provide the value in bp that the record position should be shifted. 
The value should be 1/2 the average length of the insert library 
that was sequenced. The default is to automatically and empirically 
determine the appropriate shift value. See below for the approach.

=item --sample <integer>

Indicate the number of top regions to sample when empirically 
determining the shift value. The default is 100.

=item --strand

Indicate that separate wig files should be written for each strand. 
The output file basename is appended with either '_f' or '_r' for 
both files. For paired-end RNA-Seq alignments that were generated with 
TopHat, the XS attribute is honored for strand information. 
The default is to take all alignments regardless of strand.

=item --qual <integer>

Set a minimum mapping quality score of alignments to count. The mapping
quality score is a posterior probability that the alignment was mapped
incorrectly, and reported as a -10Log10(P) value, rounded to the nearest
integer (range 0..255). Higher numbers are more stringent. For performance
reasons, when counting paired-end reads, only the left alignment is
checked. The default value is 0 (accept everything).

=item --fix

=item --inter

Specify whether or not to record interpolating positions with a count 
of 0. If true, a fixedStep wig file is written, otherwise a 
variableStep wig file is written that only records the positions 
where a tag is found. This will also work with bedGraph output. 
This option is automatically enabled when the --bin option is enabled. 
The default behavior is to not record empty positions.

=item --bed

Specify whether or not to write a bedGraph (chromosome start stop value) 
file or a traditional fixedStep or variableStep wiggle file. The 
default is false. 

=item --rpm

Convert the data to Reads (or Fragments) Per Million mapped. This is useful 
for comparing read coverage between different datasets. This conversion 
is applied before converting to log, if requested. The default is no 
conversion.

=item --log [2|10]

Transform the count to a log scale. Specify the base number, 2 or 
10. The counts are increased by 1 before taking a log transformation, 
thus avoiding taking a log of 0. Only really useful with Bam alignment 
files with high count numbers. Default is to not transform the count.

=item --(no)track

Specify whether or not to include a track line in the wig file. This 
option is automatically disabled when further converting to a bigWig 
file.

=item --bw

Specify whether or not the wig file should be further converted into 
an indexed, compressed, binary BigWig file. The default is false.

=item --bwapp < /path/to/wigToBigWig or /path/to/bedGraphToBigWig >

Specify the full path to Jim Kent's bigWig conversion utility. Two 
different utilities may be used, bedGraphToBigWig or wigToBigWig, 
depending on the format of the wig file generated. The wigToBigWig 
is preferred only because smaller file sizes are produced, but with 
slightly higher memory requirements. The application paths may be 
set in the biotoolbox.cfg file.

=item --(no)gz

Specify whether (or not) the output file should be compressed with 
gzip. The default is compress the output unless a BigWig file is 
requested.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will enumerate aligned sequence tags and generate a wig, 
or optionally BigWig, file. Each position in the wig file records the 
number of alignments which map to that position. Alignments may be 
counted at the start (5') or midpoint positions, or optionally 
enumerated at every position across the alignment (resulting in a 
coverage map rather than an alignment enumeration). Additionally, 
alignments may be summed in bins (10, 100, 237, etc bp), which may 
facilitate viewing in genome browsers zoomed out.

Both paired-end and single-end alignments may be counted. Alignments 
with splices (e.g. RNA-Seq) may be counted singly or separately. 
Alignment counts may be separated by strand, facilitating analysis of 
RNA-Seq experiments. 

For ChIP-Seq experiments, the alignment position may be shifted 
in the 3' direction. This effectively merges the separate peaks 
(representing the ends of the enriched fragments) on each strand 
into a single peak centered over the target locus. The shift value 
may be empirically determined from the sequencing data (see below). 

The output wig file may be either a variableStep, fixedStep, or 
bedGraph test format. The wig file may be further converted into a 
compressed, indexed, binary bigWig format, dependent on the availability 
of the appropriate conversion utilities. 

More information about wiggle files can be found at 
http://genome.ucsc.edu/goldenPath/help/wiggle.html, bedGraph at 
http://genome.ucsc.edu/goldenPath/help/bedgraph.html, and bigWig at 
http://genome.ucsc.edu/goldenPath/help/bigWig.html.

=head1 RECOMMENDED SETTINGS

The type of wig file to generate for your Bam sequencing file can vary 
depending on your particular experimental application. Here are a few 
common sequencing applications and my recommended settings for generating 
the wig or bigWig file.

=over

=item Straight coverage

To generate a straight-forward coverage map, similar to what most genome 
browsers display when using a Bam file as source, use the following 
settings:
 
 bam2wig.pl --coverage --in <bamfile>

=item Single-end ChIP-Seq

When sequencing Chromatin Immuno-Precipitation products, one generally 
is more interested in the number of tag counts, rather than coverage. 
Hence, we can simply count the start position of the sequence tags. 

To adjust the positions of tag count peaks to center over the presumed 
site of interest, let the program empirically determine the shift 
value from the sequence data (recommended). Otherwise, if you know 
the mean size of your ChIP eluate fragments, you can use the --shiftval 
option. 

Finally, to compare ChIP-Seq alignments from multiple experiments, 
convert your reads to Reads Per Million Mapped, which will help to 
normalize read counts.
 
 bam2wig.pl --pos start --shift --rpm --in <bamfile>

=item Paired-end ChIP-Seq

If both ends of the ChIP eluate fragments are sequenced, then we do not 
need to calculate a shift value. Instead, we will simply count at the 
midpoint of each properly-mapped sequence pair.
 
 bam2wig.pl --pos mid --pe --rpm --in <bamfile>

=item Unstranded RNA-Seq

With RNA-Sequencing, we may be interested in either coverage (generating 
a transcriptome map) or simple tag counts (differential gene expression), 
so we can count in one of two ways. 

To compare RNA-Seq data from different experiments, convert the read 
counts to Reads Per Million Mapped, which will help to normalize read 
counts.
 
 bam2wig --pos span --rpm --in <bamfile>
 
 bam2wig --pos mid --rpm --in <bamfile>

=item Stranded, single-end RNA-Seq

If the library was generated in such a way as to preserve strand, then 
we can separate the counts based on the strand of the alignment. Note 
that the reported strand may be accurate or flipped, depending upon 
whether first-strand or second-strand synthesized cDNA was sequenced, 
and whether your aligner took this into account. Please check the 
output wig files in a genome browser to verify which one is which, and 
rename appropriately.
 
 bam2wig --pos mid --strand --rpm --in <bamfile>
 
 bam2wig --pos span --strand --rpm --in <bamfile>

=item Stranded, paired-end RNA-Seq

Strand presents a complication when sequencing both ends of the cDNA 
product from a library that preserves orientation. Currently, 
the TopHat aligner can handle stranded, paired-end RNA-Seq alignments. 
Because each pair will align to both strands, the aligner must record 
separately which strand the original fragment should align. The TopHat 
program records an 'XS' attribute for each alignment, and, if present, 
bam2wig.pl will use this to set the strand.
 
 bam2wig --pe --pos mid --strand --rpm --in <bamfile>
 
 bam2wig --pe --pos span --strand --rpm --in <bamfile>

=back

=head1 SHIFT VALUE DETERMINATION

To determine the shift value, the top enriched 500 bp regions (the 
default number is 100) from the largest chromosome are identified by 
their read coverage. Stranded read counts are then collected in 10 bp 
bins over a 1.5 kb region encompassing the identified region. A 
Pearson product-moment correlation coefficient is then reiteratively 
determined between the stranded data as the bins are shifted from 
30 to 500 bp. The shift corresponding to the highest R squared value 
is retained for each sampled region, and the mean best shift for all 
sampled regions is used as the final shift value. This approach works 
best with clean, distinct peaks. Not all sampled regions may return 
a significant R squared value. The peak shift may be evaluated by 
viewing separate, stranded wig files together with the shifted wig 
file in a genome browser.

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

