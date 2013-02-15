#!/usr/bin/env perl

=head1 NAME

bam2wig.pl

A script to enumerate Bam alignments or coverage into a wig file

=head1 SYNOPSIS

bam2wig.pl [--options...] <filename.bam>
  
  Options:
  --in <filename.bam>
  --out <filename> 
  --position [start|mid|span|extend]
  --coverage
  --splice|split
  --pe
  --bin <integer>
  --shift
  --shiftval <integer>
  --sample <integer>
  --chrom <integer>
  --minr <float>
  --strand
  --qual <integer>
  --max <integer>
  --rpm
  --log [2|10]
  --bw
  --bwapp </path/to/wigToBigWig or /path/to/bedGraphToBigWig>
  --(no)gz
  --version
  --help

=cut

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(sum min max mean stddev);
use Statistics::LineFit;
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

# Declare constants for this program
use constant {
	VERSION         => '1.10',
	LOG2            => log(2),
	LOG10           => log(10),
	ALIGN_COUNT_MAX => 200_000, # Maximum number of alignments processed before writing 
	                            # to file. Increasing this number may improve 
	                            # performance at the cost of memory usage.
	BUFFER_MIN      => 1_200, # Leave this much in bp in the coverage buffer when
	                          # writing a bedGraph file to account for additional 
	                          # future coverage. Increase this if alignment length, 
	                          # paired-end span, or 2 x read shift is greater than 
	                          # this value.
};
	
	

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
	$chr_number,
	$correlation_min,
	$strand,
	$bin_size,
	$min_mapq,
	$max_dup,
	$rpm,
	$log,
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
	'chrom=i'   => \$chr_number, # number of chromosomes to sample
	'minr=f'    => \$correlation_min, # R^2 minimum value for shift
	'strand!'   => \$strand, # separate strands
	'bin=i'     => \$bin_size, # size of bin to make
	'qual=i'    => \$min_mapq, # minimum mapping quality
	'max=i'     => \$max_dup, # maximum duplicate positions
	'rpm!'      => \$rpm, # calculate reads per million
	'log=i'     => \$log, # transform count to log scale
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
	print " Biotoolbox script bam2wig.pl, version " . VERSION . "\n\n";
	exit;
}



### Check for requirements and set defaults
# global variables
my ($use_start, $use_mid, $use_span, $bin, $bedgraph, $callback, 
	$write_wig, $convertor, $data_ref);
check_defaults();

# record start time
my $start_time = time;





### Open files
# Bam file
unless (exists &open_bam_db) {
	die " unable to load Bam file support! Is Bio::DB::Sam installed?\n"; 
}
my $sam = open_bam_db($infile) or die " unable to open bam file '$infile'!\n";
my $bam = $sam->bam;
my $index = $sam->bam_index;





### Calculate Total read numbers
my $total_read_number = 0;
if ($rpm) {
	# this is only required when calculating reads per million
	print " Calculating total number of aligned fragments... this may take a while...\n";
	$total_read_number = sum_total_bam_alignments($sam, $min_mapq, $paired);
	print "   ", format_with_commas($total_read_number), " total mapped fragments\n";
	printf " counted in %.1f minutes\n", (time - $start_time)/60;
}



### Calculate shift value
if ($shift and !$shift_value) {
	print " Calculating 3' shift value...\n";
	$shift_value = determine_shift_value();
	
	# precalculate double shift when recording extended position
	if ($position eq 'extend') {
		$shift_value = $shift_value * 2;
		print " Reads will be extended by $shift_value bp\n";
	}
	else {
		print " Reads will be shifted by $shift_value bp\n";
	}
}



### Process bam file
# process according to type of data collected and alignment type
if ($use_coverage) {
	# special, speedy, low-level, single-bp coverage 
	process_bam_coverage();
}
elsif ($splice) {
	# single end alignments with splices require special callbacks
	# we must do this with the high level API
	$sam->split_splices($splice);
	warn " WARNING: enabling splices will increase processing times\n";
	process_split_alignments();
}
else {
	# all other alignments, single, paired, mid, start, span
	process_alignments();
}




### Finish

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
	
	# maximum duplicates
	unless (defined $max_dup) {
		$max_dup = 1000;
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
		elsif ($position eq 'extend') {
			if ($paired) {
				$use_span = 1;
			}
			else {
				$use_span = 1;
				$shift = 1;
			}
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
			if (!$shift and !$strand and !$log and !$rpm and !$splice) {
				$use_coverage = 1;
			}
			else {
				$use_start = 1;
			}
		}
	}
	if ($paired and $use_start) {
		warn " using midpoint with paired-end reads\n";
		undef $use_start;
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
		undef $splice;
	}
	
	# check to shift position or not
	if ($shift and ($paired or $splice) ) {
		warn " disabling shift with paired reads\n" if $paired;
		warn " disabling shift with splices enabled\n" if $splice;
		undef $shift;
	}
	if ($shift) {
		unless ($sample_number) {
			$sample_number = 200;
		}
		unless ($chr_number) {
			$chr_number = 2;
		}
		if (defined $correlation_min) {
			if ($correlation_min <= 0 or $correlation_min >= 1) {
				die " cannot use minimum correlation value of $correlation_min!\n" .
					" use --help for more information\n";
			}
		}
		else {
			$correlation_min = 0.25;
		}
	}
	
	# check bin size
	if ($bin_size) {
		# set the boolean variable as to whether we're binning above 
		# the default 1 bp resolution 
		
		# cannot use bin with spanned features
		if ($use_span) {
			$bin_size = 1;
			undef $bin;
			warn " disabling bin when recording span or extended positions\n". 
				"   a bedGraph file will automatically be written\n";
		}
		else {
			$bin = $bin_size > 1 ? 1 : 0;
		}
	}
	else {
		# default bin size, required to be set
		$bin_size = 1;
	}
	
	# determine output format
	$bedgraph = $use_span ? 1 : 0; 
	
	# look for bigwig app
	if ($bigwig) {
		# we need to set some options prior to writing the wig file if 
		# we're going to be writing a bigWig later
		find_bigwig_app();
	}
	
	# check output file
	unless ($outfile) {
		$outfile = $infile;
		$outfile =~ s/\.bam$//;
	}
	if (defined $gz) {
		# overide to false if bigwig is true
		undef $gz if $bigwig;
	} 
	else {
		# default is to use compression unless a bigwig file is requested
		# then the file is only temporary anyway
		$gz = $bigwig ? 0 : 1;
	}
	$outfile =~ s/\.(?:wig|bdg|bedgraph)(?:\.gz)?$//i; # strip extension if present
	
	# determine the alignment callback method
	# and the wig writing method
	# based on the position used, strandedness, and/or shift
	if ($use_coverage) {
		# no callbacks necessary, special mode
		print " recording coverage spanning alignments\n";
	}
	elsif ($use_start and $strand and $shift) {
		$callback  = \&record_stranded_shifted_start;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording stranded, shifted-start positions\n";
	}
	elsif ($use_start and $strand and !$shift) {
		$callback  = \&record_stranded_start;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording stranded, start positions\n";
	}
	elsif ($use_start and !$strand and $shift) {
		$callback  = \&record_shifted_start;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording shifted-start positions\n";
	}
	elsif ($use_start and !$strand and !$shift) {
		$callback  = \&record_start;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording start positions\n";
	}
	elsif ($use_mid and $strand and $shift and $paired) {
		# this should not happen, shift is disabled with paired
		die " programming error!\n";
	}
	elsif ($use_mid and $strand and $shift and !$paired) {
		$callback  = \&record_stranded_shifted_mid;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording stranded, shifted-mid positions\n";
	}
	elsif ($use_mid and $strand and !$shift and $paired) {
		$callback = \&record_stranded_paired_mid;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording stranded, mid positions of pairs\n";
	}
	elsif ($use_mid and $strand and !$shift and !$paired) {
		$callback = \&record_stranded_mid;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording stranded, mid positions\n";
	}
	elsif ($use_mid and !$strand and $shift and $paired) {
		# this should not happen, shift is disabled with paired
		die " programming error!\n";
	}
	elsif ($use_mid and !$strand and $shift and !$paired) {
		$callback  = \&record_shifted_mid;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording shifted-mid positions\n";
	}
	elsif ($use_mid and !$strand and !$shift and $paired) {
		$callback  = \&record_paired_mid;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording mid positions of pairs\n";
	}
	elsif ($use_mid and !$strand and !$shift and !$paired) {
		$callback  = \&record_mid;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording mid position\n";
	}
	elsif ($use_span and $strand and $paired) {
		$callback  = \&record_stranded_paired_span;
		$write_wig = \&write_bedgraph;
		print " recording stranded positions spanning paired alignments\n";
	}
	elsif ($use_span and $strand and !$shift and !$paired) {
		$callback  = \&record_stranded_span;
		$write_wig = \&write_bedgraph;
		print " recording stranded positions spanning alignments\n";
	}
	elsif ($use_span and !$strand and $paired) {
		$callback  = \&record_paired_span;
		$write_wig = \&write_bedgraph;
		print " recording positions spanning paired alignments\n";
	}
	elsif ($use_span and !$strand and !$shift and !$paired) {
		if (!$rpm and !$log) {
			# this is actually coverage
			undef $use_span;
			undef $bedgraph; # no bedgraph file, coverage uses fixedStep
			$use_coverage = 1;
			print " recording coverage spanning alignments\n";
		}
		else {
			$callback  = \&record_span;
			$write_wig = \&write_bedgraph;
			print " recording positions spanning alignments\n";
		}
	}
	elsif ($use_span and !$strand and $shift and !$paired) {
		$callback  = \&record_extended;
		$write_wig = \&write_bedgraph;
		print " recording extended alignments\n";
	}
	else {
		# what else is left!?
		die " programming error!\n";
	}
	
	# determine the convertor callback
	if ($rpm and !$log) {
		$convertor = sub {
			return ($_[0] * 1_000_000) / $total_read_number;
		};
	}
	elsif ($rpm and $log == 2) {
		# calculate rpm first before log
		$convertor = sub {
			return log( ( ($_[0] * 1_000_000) / $total_read_number) ) + 1 / LOG2;
		};
	}
	elsif ($rpm and $log == 10) {
		# calculate rpm first before log
		$convertor = sub {
			return log( ( ($_[0] * 1_000_000) / $total_read_number) ) + 1 / LOG10;
		};
	}
	elsif (!$rpm and $log == 2) {
		$convertor = sub {
			return log($_[0] + 1) / LOG2;
		};
	}
	elsif (!$rpm and $log == 10) {
		$convertor = sub {
			return log($_[0] + 1) / LOG10;
		};
	}
	elsif (!$rpm and !$log) {
		# this one is easy!
		$convertor = sub {
			return $_[0];
		};
	}
}


### Find the appropriate executable for generating a bigWig
sub find_bigwig_app {
	# check for the app
	unless ($bwapp) {
		# first check the biotoolbox configuration file for the UCSC 
		# conversion utilities
		# if it's not listed, then check the environment path 
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
	unless ($bwapp =~ /ToBigWig$/) {
		warn " Unable to find bigWig conversion utility!\n" .
			" Generating text file only\n";
		$bigwig = 0;
	}
	if ($bedgraph) {
		if ($bwapp =~ /wigToBigWig$/) {
			warn " Wrong utility! Need bedGraphToBigWig.\n";
			$bigwig = 0;
		}
	}
	else {
		if ($bwapp =~ /bedGraphToBigWig$/) {
			warn " Wrong utility! Need wigToBigWig.\n";
			$bigwig = 0;
		}
	}
}


### Open the output file handle 
sub open_wig_file {
	
	# generate name
	my $name = shift;
	$name .= $bedgraph ? '.bdg' : '.wig';
	$name .= '.gz', if $gz;
		
	
	# open
	my $fh = open_to_write_fh($name, $gz) or 
		die " unable to open output wig file '$name'!\n";
		
	# write track line
	if ($bedgraph) {
		$fh->print("track type=bedGraph\n");
	}
	else {
		$fh->print("track type=wiggle_0\n");
	}
	
	# finished, return ref to names array, and 1 or 2 filehandles
	return ($name, $fh);
}


### Determine the shift value
sub determine_shift_value {
	
	# find the biggest chromosome to sample
		# this is assuming all the chromosomes have different sizes ;-)
	my %size2chrom;
	for my $tid (0 .. $sam->n_targets - 1) {
		# key is chromosome size, value is its name
		$size2chrom{ $sam->target_len($tid) } = $tid;
	}
	
	# identify top regions to score
	# we will walk through the largest chromosome(s) looking for the top  
	# 1 kb regions containing the highest unstranded coverage to use
	my %coverage2region;
	
	# look for test regions
	# sampling the number of requested regions from number of requested chromosomes
	print "  searching for the top $sample_number coverage regions\n" .
		"   on the largest $chr_number chromosomes to sample...\n";
	my $chrom_count = 0;
	for my $size (sort {$b <=> $a} keys %size2chrom) {
		# sort largest to smallest
		last if $chrom_count == $chr_number; 
		my $tid = $size2chrom{$size};
		for (my $start = 0; $start < $size; $start += 500) {
		
			my $end = $start + 500;
			$end = $size if $end > $size;
		
			# using the low level interface for a little more performance
			my $coverage = $index->coverage(
				$bam,
				$tid, # need the low level target ID
				$start, 
				$end,
			);
			my $sum_coverage = sum( @{$coverage} );
			next if $sum_coverage == 0;
		
			# check if our coverage exceeds the lowest region
			if (scalar keys %coverage2region < $sample_number) {
				# less than requested regions found so far, so keep it
				# record the coordinates for this region
				$coverage2region{$sum_coverage} = [$tid, $start, $end];
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
					# record the coordinates for this region
					$coverage2region{$sum_coverage} = [$tid, $start, $end];
				}
			}
		}
		$chrom_count++;
	}
	printf "    done in %.1f minutes\n", (time - $start_time)/60;
	
		
	# now determine the optimal shift for each of the test regions
	my @shift_values;
	foreach my $i (sort {$a <=> $b} keys %coverage2region) {
		
		# get the start and end positions
		# we're adjusting them by 500 bp in both directions to actually 
		# sample a 1.5 kb region centered over the original region
		my $tid   = $coverage2region{$i}->[0];
		my $chrom = $sam->target_name($tid);
		my $start = $coverage2region{$i}->[1] - 500;
		my $end   = $coverage2region{$i}->[2] + 500;
		# just in case we go over
		$start = 0 if $start < 1;
		$end = $sam->target_len($tid) if $end > $sam->target_len($tid);
		
		
		# collect stranded data from our sample window
		my %data = {
			'f'       => {},
			'r'       => {},
			'count'   => 0,
			'print'   => 0, # we will not print this data to file
		};
		$index->fetch($bam, $tid, $start, $end, \&record_stranded_start, \%data);
		unless ($data{'count'}) {
			die " no stranded data collected for calculating shift value!\n";
		}
		
		# generate data arrays
		my @f;
		my @r;
		for (my $i = $start; $i <= $end; $i += 10) {
			# walk in the region by increments in 10
			# sum the counts in each 10 bp interval, and push it to the array
			push @f, sum( map { $data{f}{$_} ||= 0 } ($i .. $i+9) );
			push @r, sum( map { $data{r}{$_} ||= 0 } ($i .. $i+9) );
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
			if ($r2 >= $correlation_min) {
				$r2shift{$r2} = $i * 10;
			}
		}
		
		# determine best shift
		my $string = "  sampling $chrom:$start..$end  ";
		if (%r2shift) {
			my $max_r = max(keys %r2shift);
			push @shift_values, $r2shift{$max_r};
			printf "$string shift %s bp (r^2 %.3f)\n", $r2shift{$max_r}, $max_r;
		}
		else {
			print "$string\n";
		}
	}
	
	# determine the optimal shift value
	# we will be using a trimmed mean value to avoid outliers
	@shift_values = sort {$a <=> $b} @shift_values;
	my $cut = int( ( scalar(@shift_values) / 10 ) + 0.5); # take 10%
	$cut++ if ($cut % 2); # make an even number
	splice(@shift_values, 0, $cut);
	splice(@shift_values, scalar(@shift_values) - $cut);
	my $best_value = mean(@shift_values);
	printf "  The mean shift value is %.0f +/- %.0f bp\n", 
		$best_value, stddev(@shift_values);
	
	# done
	return sprintf("%.0f", $best_value);
}

### Collect alignment coverage
sub process_bam_coverage {
	# using the low level bam coverage method, not strand specific
	
	# open wig file
	my ($filename, $fh) = open_wig_file($outfile);
	
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
		
		# prepare definition line for fixedStep
		$fh->print(
			"fixedStep chrom=$seq_id start=1 step=$bin_size span=$bin_size\n"
		);
		
		# walk through the chromosome
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
			if ($bin) {
				for (my $i = 0; $i < scalar(@{ $coverage }); $i += $bin_size) {
				
					# sum the reads within our bin
					my $sum = sum( map {$coverage->[$_]} ($i .. $i + $bin_size - 1));
				
					# print the wig line
					$fh->print("$sum\n");
				}
			}
			else {
				for my $i (0 .. scalar(@$coverage)-1) {
					# print the wig line
					$fh->print("$coverage->[$i]\n");
				}
			}
		}
		printf " Converted reads on $seq_id in %.3f minutes\n", (time - $start_time)/60;
	}
	
	$fh->close;
	print " Wrote file $filename\n";
	
	# convert to bigwig if requested
	if ($bigwig) {
		convert_to_bigwig($filename);
	}
}



### Process alignments
sub process_alignments {
	
	# open wig files
	my ($filename1, $filename2, $fh1, $fh2);
	if ($strand) {
		($filename1, $fh1) = open_wig_file("$outfile\_f");
		($filename2, $fh2) = open_wig_file("$outfile\_r");
	}
	else {
		($filename1, $fh1) = open_wig_file($outfile);
	}
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as 
		# a numeric target identifier
		# we can easily convert this to an actual sequence name
		# we will force the conversion to go one chromosome at a time
		
		# sequence info
		my $seq_id = $sam->target_name($tid);
		my $seq_length = $sam->target_len($tid);
		
		# print chromosome line
		if ($use_start or $use_mid) {
			foreach ($fh1, $fh2) {
				next unless defined $_;
				if ($bin) {
					$_->print(
						"fixedStep chrom=$seq_id start=1 step=$bin_size span=$bin_size\n"
					);
				}
				else {
					$_->print("variableStep chrom=$seq_id\n");
				}
			}
		}
		
		# process the chromosome alignments
		my %data = (
			'f'       => {},
			'r'       => {},
			'fhf'     => $fh1,
			'fhr'     => $fh2,
			'bufferf' => [],
			'bufferr' => [],
			'offsetf' => $bedgraph ? 0 : 1, # bedgraph used 0-base indexing
			'offsetr' => $bedgraph ? 0 : 1, # wig used 1-base indexing
			'count'   => 0,
			'seq_id'  => $seq_id,
			'seq_length' => $seq_length,
			'print'   => 1,
		);
		$index->fetch($bam, $tid, 0, $seq_length, $callback, \%data);
		
		# finish up this chromosome
		&$write_wig(\%data, 1); # final write
		printf " Converted %s alignments on $seq_id in %.3f minutes\n", 
			format_with_commas( $data{'count'}), (time - $start_time)/60;
	}
	
	# finished
	$fh1->close;
	$fh2->close if $fh2;
	print " Wrote file $filename1\n";
	print " Wrote file $filename2\n" if $filename2;
	
	# convert to bigwig if requested
	if ($bigwig) {
		convert_to_bigwig($filename1, $filename2);
	}
}




### Walk through the alignments on each chromosome
sub process_split_alignments {
	
	# open wig files
	my ($filename1, $filename2, $fh1, $fh2);
	if ($strand) {
		($filename1, $fh1) = open_wig_file("$outfile\_f");
		($filename2, $fh2) = open_wig_file("$outfile\_r");
	}
	else {
		($filename1, $fh1) = open_wig_file($outfile);
	}
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as 
		# a numeric target identifier
		# we can easily convert this to an actual sequence name
		# we will force the conversion to go one chromosome at a time
		
		# sequence info
		my $seq_id = $sam->target_name($tid);
		my $seq_length = $sam->target_len($tid);
		
		# print chromosome line
		if ($use_start or $use_mid) {
			foreach ($fh1, $fh2) {
				next unless defined $_;
				if ($bin) {
					$_->print(
						"fixedStep chrom=$seq_id start=1 step=$bin_size span=$bin_size\n"
					);
				}
				else {
					$_->print("variableStep chrom=$seq_id\n");
				}
			}
		}
		
		# process the reads across the chromosome
		# Due to requiring split splices, we need to use the high-level API
		# it's easier to use the high-level than try and re-invent my own code 
		# for dealing with splices <sigh>
		# this will increase processing time as each alignment must be incorporated 
		# into AlignWrapper objects
		
		# since the high level API does not allow for a data structure to be passed 
		# to the callback, we'll have to link it to a global variable
		my %data = (
			'f'       => {},
			'r'       => {},
			'fhf'     => $fh1,
			'fhr'     => $fh2,
			'bufferf' => [],
			'bufferr' => [],
			'offsetf' => $bedgraph ? 0 : 1, # bedgraph used 0-base indexing
			'offsetr' => $bedgraph ? 0 : 1, # wig used 1-base indexing
			'count'   => 0,
			'seq_id'  => $seq_id,
			'seq_length' => $seq_length,
			'print'   => 1,
		);
		$data_ref = \%data;
		
		# fetch using the high level API
		$sam->fetch("$seq_id:1..$seq_length", \&single_end_spliced_callback);
		
		# finish up this chromosome
		&$write_wig(\%data, 1); # final write
		printf " Converted %s alignments on $seq_id in %.3f minutes\n", 
			format_with_commas( $data{'count'}), (time - $start_time)/60;
	}
	
	# finished
	$fh1->close;
	$fh2->close if $fh2;
	print " Wrote file $filename1\n";
	print " Wrote file $filename2\n" if $filename2;
	
	# convert to bigwig if requested
	if ($bigwig) {
		convert_to_bigwig($filename1, $filename2);
	}
}


### Callback for processing single-end split alignments
sub single_end_spliced_callback {
	my $a = shift;
	
	# check alignment
	return if $a->unmapped;
	
	# check for subfeatures
	my @subfeatures = $a->get_SeqFeatures;
	if (@subfeatures) {
		# process each subfeature
		foreach my $subf (@subfeatures) {
			&$callback($subf, $data_ref);
		}
	}
	else {
		# no subfeatures found
		# treat this as a single read
		&$callback($a, $data_ref);
	}
}


### Record stranded at shifted start position
sub record_stranded_shifted_start {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	return if $a->calend <= $shift_value; # cannot have negative positions
	
	# record based on the strand
	if ($a->reversed) {
		# reverse strand
		$data->{r}{ $a->calend - $shift_value }++;
	}
	else {
		# forward strand
		$data->{f}{ $a->pos + 1 + $shift_value }++;
	}
	check_data($data);
}


### Record stranded at start position
sub record_stranded_start {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
	# record based on the strand
	if ($a->reversed) {
		# reverse strand
		$data->{r}{ $a->calend }++;
	}
	else {
		# forward strand
		$data->{f}{ $a->pos + 1 }++;
	}
	check_data($data);
}


### Record at shifted start position
sub record_shifted_start {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	return if $a->calend <= $shift_value; # cannot have negative positions
	
	# shift based on strand, record on forward
	if ($a->reversed) {
		# reverse strand
		$data->{f}{ $a->calend - $shift_value }++;
	}
	else {
		# forward strand
		$data->{f}{ $a->pos + 1 + $shift_value }++;
	}
	check_data($data);
}


### Record at start position
sub record_start {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
	# start based on strand, record on forward
	if ($a->reversed) {
		# reverse strand
		$data->{f}{ $a->calend }++;
	}
	else {
		# forward strand
		$data->{f}{ $a->pos + 1 }++;
	}
	check_data($data);
}


### Record stranded at shifted mid position
sub record_stranded_shifted_mid {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
	# calculate mid position
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	return if $pos <= $shift_value; # cannot have negatives
	
	# record based on the strand
	if ($a->reversed) {
		# reverse strand
		$data->{r}{ $pos - $shift_value }++;
	}
	else {
		# forward strand
		$data->{f}{ $pos + $shift_value }++;
	}
	check_data($data);
}


### Record stranded at mid position
sub record_stranded_mid {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
	# calculate mid position
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	
	# record based on the strand
	if ($a->reversed) {
		# reverse strand
		$data->{r}{ $pos }++;
	}
	else {
		# forward strand
		$data->{f}{ $pos }++;
	}
	check_data($data);
}


### Record at shifted mid position
sub record_shifted_mid {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
	# calculate mid position
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	return if $pos <= $shift_value; # cannot have negatives
	
	# shift based on strand, record on forward
	if ($a->reversed) {
		# reverse strand
		$data->{f}{ $pos - $shift_value }++;
	}
	else {
		# forward strand
		$data->{f}{ $pos + $shift_value }++;
	}
	check_data($data);
}


### Record at mid position for paired end
sub record_paired_mid {
	my ($a, $data) = @_;
	
	# check
	return if $a->reversed; # only take left (forward) alignments, not right
	return unless $a->proper_pair;
	return unless $a->mreversed;
	return unless $a->tid == $a->mtid; # on the same chromosome
	return if $a->qual < $min_mapq;
	
	# calculate mid position of the pair
	my $position = $a->pos + 1 + int( $a->isize / 2 );
	
	# record position in forward strand data hash
	$data->{'f'}{$position}++;
	check_data($data);
}


### Record stranded at mid position for paired end
sub record_stranded_paired_mid {
	my ($a, $data) = @_;
	
	# check
	return if $a->reversed; # only take left (forward) alignments, not right
	return unless $a->proper_pair;
	return unless $a->mreversed;
	return unless $a->tid == $a->mtid; # on the same chromosome
	return if $a->qual < $min_mapq;
	
	# calculate mid position of the pair
	my $position = $a->pos + 1 + int( $a->isize / 2 );
	
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
		$data->{'f'}{$position}++;
	}
	elsif ($strand eq '-') {
		# record position in reverse strand data hash
		$data->{'r'}{$position}++;
	}
	else {
		my $name = $a->qname;
		die " unrecognized strand value '$strand' for XS attribute in BAM record $name\n";
	}
	check_data($data);
}


### Record at mid position
sub record_mid {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
	# record mid position
	$data->{f}{ int( ($a->pos + 1 + $a->calend) / 2) }++;
	
	check_data($data);
}


### Record stranded across alignment for paired end
sub record_stranded_paired_span {
	my ($a, $data) = @_;
	
	# check
	return if $a->reversed; # only take left (forward) alignments, not right
	return unless $a->proper_pair;
	return unless $a->mreversed;
	return unless $a->tid == $a->mtid; # on the same chromosome
	return if $a->qual < $min_mapq;
	
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
		$data->{'f'}{ $a->pos } .= $a->isize . ",";
	}
	elsif ($strand eq '-') {
		# record position in reverse strand data hash
		$data->{'r'}{ $a->pos } .= $a->isize . ",";
	}
	else {
		my $name = $a->qname;
		die " unrecognized strand value '$strand' for XS attribute in BAM record $name\n";
	}
	check_data($data);
}


### Record stranded across alignment
sub record_stranded_span {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
	# record length at the start position based on strand
	if ($a->reversed) {
		# reverse strand
		$data->{r}{ $a->pos } .= $a->calend - $a->pos . ',';
	}
	else {
		# forward strand
		$data->{f}{ $a->pos } .= $a->calend - $a->pos . ',';
	}
	check_data($data);
}


### Record across alignment for paired end
sub record_paired_span {
	my ($a, $data) = @_;
	
	# check
	return if $a->reversed; # only take left (forward) alignments, not right
	return unless $a->proper_pair;
	return unless $a->mreversed;
	return unless $a->tid == $a->mtid; # on the same chromosome
	return if $a->qual < $min_mapq;
	
	# record isize at start
	$data->{'f'}{ $a->pos } .= $a->isize . ",";
	check_data($data);
}


### Record across alignment
sub record_span {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
	# record the length
	$data->{f}{ $a->pos } .= $a->calend - $a->pos . ',';
	check_data($data);
}



### Record extended alignment
sub record_extended {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
	# start based on strand, record on forward
	$data->{f}{ $a->pos } .= "$shift_value,";
	check_data($data);
}


### Increment count and check the data size for writing
sub check_data {
	my $data = $_[0];
	$data->{'count'}++;
	
	# write when we reach buffer maximum number of alignments read
	if ($data->{'print'} and $data->{'count'} % ALIGN_COUNT_MAX == 0) {
		&$write_wig($data, 0);
	}
}


### Write a varStep wig file
sub write_varstep {
	my ($data, $final) = @_;
	
	# do each strand one at a time
	foreach my $s (qw(f r)) {
		next unless keys( %{ $data->{$s} } );
		my $fh     = "fh$s";
		
		# check the maximum position that we cannot go beyond
		# defined either by the minimum buffer value or the end of the chromosome
		my $maximum = $final ?  $data->{'seq_length'} : 
			max(keys %{$data->{$s}}) - BUFFER_MIN; 
		
		# write the data
		foreach my $pos (sort {$a <=> $b} keys %{$data->{$s}}) {
			
			# first check the position and bail when we've gone far enough
			last if ($pos > $maximum);
			
			# write line
			my $score = $data->{$s}{$pos} > $max_dup ? $max_dup : $data->{$s}{$pos};
			$data->{$fh}->print( join("\t", 
				$pos, 
				&$convertor($score) # convert RPM or log before writing
			) . "\n");
			
			# clean up
			delete $data->{$s}{$pos};
		}
		
		# warn about tossed values at the chromosome end
		if ($final and keys(%{ $data->{$s} }) ) {
			warn "  Warning: " . scalar keys(%{ $data->{$s} }) . " data points trimmed" .
				" from " . $data->{'seq_id'} . " end\n";
		}
	}
}


### Write a fixedStep wig file
sub write_fixstep {
	my ($data, $final) = @_;
	
	# do each strand one at a time
	foreach my $s (qw(f r)) {
		next unless keys( %{ $data->{$s} } );
		my $offset = "offset$s";
		my $fh     = "fh$s";
		
		# check the maximum position that we cannot go beyond
		# defined either by the minimum buffer value or the end of the chromosome
		my $maximum = $final ?  $data->{'seq_length'} : 
			max(keys %{$data->{$s}}) - BUFFER_MIN; 
		
		# write the data
		for (
			my $pos = $data->{$offset}; 
			$pos < $maximum - $bin_size; 
			$pos += $bin_size
		) {
			
			# sum the counts in the bin interval
			my $score = sum( 
				map { $data->{$s}{$_} ||= 0 } ($pos .. $pos + $bin_size -1) );
			$score = $score > $max_dup ? $max_dup : $score;
			
			# write line
			$data->{$fh}->print( &$convertor($score) . "\n" );
			
			# clean up
			delete $data->{$s}{$pos};
		}
		
		# warn about tossed values at the chromosome end
		if ($final and keys(%{ $data->{$s} }) ) {
			warn "  Warning: " . scalar keys(%{ $data->{$s} }) . " data points trimmed" .
				" from " . $data->{'seq_id'} . " end\n";
		}
	}
}


### Write a bedGraph file
sub write_bedgraph {
	my ($data, $final) = @_;
	
	# do each strand one at a time
	foreach my $s (qw(f r)) {
		next unless keys( %{ $data->{$s} } );
		my $buffer = "buffer$s";
		my $offset = "offset$s";
		my $fh     = "fh$s";
		
		# check the maximum position that we cannot go beyond
		# defined either by the minimum buffer value or the end of the chromosome
		my $maximum = $final ?  $data->{'seq_length'} : 
			max(keys %{$data->{$s}}) - BUFFER_MIN; 
		
		# convert read lengths to coverage in the buffer array
		my $pos_start_count = scalar(keys %{$data->{$s}});
		foreach my $pos (sort {$a <=> $b} keys %{ $data->{$s} }) {
			
			# first check the position and bail when we've gone far enough
			last if ($pos > $maximum);
			
			# split the lengths, limit to max, and generate coverage
			my @lengths = split(',', $data->{$s}{$pos});
			if (scalar @lengths > $max_dup) {
				splice(@lengths, $max_dup + 1); # delete the extra ones
			}
			foreach my $len (@lengths) {
				# generate coverage
				# we're relying on autovivification of the buffer array here
				# this is what makes Perl both great and terrible at the same time
				# this could also balloon memory usage - oh dear
				for (0 .. $len -1) { $data->{$buffer}->[$pos - $data->{$offset} + $_] += 1 }
			}
			delete $data->{$s}{$pos};
		}
		
		# write the array into bedgraph
		my $current_pos = $data->{$offset};
		my $current_value = shift @{ $data->{$buffer} } || 0;
		my $current_offset = 0;
		while (
			scalar( @{ $data->{$buffer} } ) > BUFFER_MIN or 
				# keep at least minimum length in the buffer
			( $final and scalar @{ $data->{$buffer} } )
				# or we're at the end of the chromosome and need to write everything
		) {
			my $value = shift @{ $data->{$buffer} } || 0;
			if ($value == $current_value) {
				# same value, extend the interval
				$current_offset++;
			}
			else {
				# write out bedgraph line for the current interval of identical values
				
				# generate the end point
				my $end = $current_pos + $current_offset + 1;
				if ($final and $end >= $maximum) {
					# we've reached the end of the chromosome
					$end = $data->{'seq_length'};
					
					# dump anything left in the array to finish up
					warn "  Warning: " . scalar(@{ $data->{$buffer} }) . " bp trimmed" .
						" from " . $data->{'seq_id'} . " end\n";
					$data->{$buffer} = [];
				}
				
				# write line
				$data->{$fh}->print( join("\t", 
					$data->{'seq_id'}, 
					$current_pos,
					$end,
					&$convertor($current_value) # convert RPM or log before writing
				) . "\n");
				
				# reset for the next interval
				$current_pos += $current_offset + 1;
				$current_value = $value;
				$current_offset = 0;
			}
		}
		
		# make sure we don't leave a value hanging
		if ($current_offset) {
			$data->{$fh}->print( join("\t", 
					$data->{'seq_id'}, 
					$current_pos,
					$current_pos + $current_offset + 1,
					&$convertor($current_value) # convert RPM or log before writing
			) . "\n");
		}
		
		# remember the current position for next writing
		$data->{$offset} = $current_pos + $current_offset + 1;
	}
}


### Run the BigWig conversion utility
sub convert_to_bigwig {
	my @filenames = @_;
	
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
	foreach my $name (@filenames) {
		next unless defined $name;
		
		# make new bw file name
		my $bw_file = $name;
		$bw_file =~ s/(?:wig|bdg)$/bw/;
		
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



__END__

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input Bam alignment file. The file should be sorted by 
genomic position and indexed, although it may be indexed automatically.

=item --out <filename>

Specify the output base filename. By default it uses the base name of 
the input file.

=item --position [start|mid|span|extend]

Specify the position of the alignment coordinate which should be 
recorded. Several positions are accepted: 
     
    start     the 5' position of the alignment
    mid       the midpoint of the alignment
    span      along the length of the alignment (coverage)
    extend    along the length of the predicted fragment
              equal to 2 x the shift value (enables --shift) 
    coverage  another way to specify --coverage
     
With paired-end alignments, the positions are relative to the entire 
insert fragment defined by two proper alignments. For single-end 
alignments, the default value is coverage if no other options are 
specified, otherwise start. For paired-end alignments, the midpoint 
position is the default.

=item --coverage

Quickly calculates the coverage of the alignments over the genome, 
either at single bp resolution (default) or in bins. This method ignores 
the position, quality, strand, shift, and log options. It is 
equivalent to specifying --position=span, --fix, --split, 
--nope, --noshift, --nostrand, --qual=0, --max=8000, --norpm, and no log. 

=item --splice

=item --split

The Bam file alignments may contain splices, where the 
read is split between two separate alignments. This is most common 
with splice junctions from RNA-Seq data. In this case, treat each 
alignment as a separate tag. This only works with single-end alignments. 
Splices are disabled for paired-end reads. Note that this will 
increase processing time.

=item --pe

The Bam file consists of paired-end alignments, and only properly 
mapped pairs of alignments will be counted. The default is to 
treat all alignments as single-end.

=item --bin <integer>

Specify the window or bin size in which alignment counts are summed. 
This option is compatible with start, mid, and coverage recording 
options, but is automatically disabled with span and extend recording 
options. The default is to count at single basepair resolution. 

=item --shift

Specify that the positions of the alignment should be shifted towards 
the 3' end. Useful for ChIP-Seq applications, where only the ends of 
the fragments are counted and often seen as separated discrete peaks 
on opposite strands flanking the true target site. This option is 
disabled with paired-end and spliced reads (where it is not needed). 

The extend position option is a special mode where the entire predicted 
ChIP fragment is recorded across the span. The length is 2 x the shift 
value. 

=item --shiftval <integer>

Provide the value in bp that the record position should be shifted. 
The value should be 1/2 the average length of the insert library 
that was sequenced. The default is to automatically and empirically 
determine the appropriate shift value. See below for the approach.

=item --sample <integer>

Indicate the number of top regions to sample when empirically 
determining the shift value. The default is 200.

=item --chrom <integer>

Indicate the number of sequences or chromosomes to sample when 
empirically determining the shift value. The sequences listed 
in the Bam file are taken in order of decreasing length, and 
one or more are taken as a representative sample of the genome. 
The default value is 2. 

=item --minr <float>

Provide the minimum R^2 value to accept a shift value when 
empirically determining the shift value. Enter a decimal value 
between 0 and 1. Higher values are more stringent. The default 
is 0.25.

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

=item --max <integer>

Set a maximum number of duplicate alignments tolerated at a single position. 
The default is 1000. 

=item --rpm

Convert the data to Reads (or Fragments) Per Million mapped. This is useful 
for comparing read coverage between different datasets. This conversion 
is applied before converting to log, if requested. This will increase 
processing time, as the alignments must first be counted. Only mapped reads 
with a minimum mapping quality are counted. All duplicates are counted. 
The default is no RPM conversion. 

=item --log [2|10]

Transform the count to a log scale. Specify the base number, 2 or 
10. The counts are increased by 1 before taking a log transformation, 
thus avoiding taking a log of 0. Only really useful with Bam alignment 
files with high count numbers. Default is to not transform the count.

=item --bw

Specify whether or not the wig file should be further converted into 
an indexed, compressed, binary BigWig file. The default is false.

=item --bwapp < /path/to/wigToBigWig or /path/to/bedGraphToBigWig >

Specify the full path to Jim Kent's bigWig conversion utility. Two 
different utilities may be used, bedGraphToBigWig or wigToBigWig, 
depending on the format of the wig file generated. The application 
paths may be set in the biotoolbox.cfg file.

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
or optionally BigWig, file. Alignments may be counted and recorded 
in several different ways. Strict enumeration may be performed and 
recorded at either the alignment's start or midpoint position. 
Alternatively, either the alignment or fragment may be recorded 
across its span. Finally, a basic unstranded, unshifted, and 
non-transformed alignment coverage may be generated. 

Both paired-end and single-end alignments may be counted. Alignments 
with splices (e.g. RNA-Seq) may be counted singly or separately. 
Alignment counts may be separated by strand, facilitating analysis of 
RNA-Seq experiments. 

For ChIP-Seq experiments, the alignment position may be shifted 
in the 3' direction. This effectively merges the separate peaks 
(representing the ends of the enriched fragments) on each strand 
into a single peak centered over the target locus. Alternatively, 
the entire predicted fragment may be recorded across its span. 
The shift value may be empirically determined from the sequencing 
data (see below). 

The output wig file may be either a variableStep, fixedStep, or 
bedGraph format. The file format is dictated by where the alignment 
position is recorded. Recording start and midpoint at 
single base-pair resolution writes a variableStep wig file. Binned start 
or midpoint counts and coverage are written as a fixedStep wig file. 
Span and extended positions are written as a bedGraph file. 

The wig file may be further converted into a compressed, indexed, binary 
bigWig format, dependent on the availability of the appropriate 
conversion utilities. 

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
performs a 3' shift adjustment to center the fragment's end reads 
over the predicted center and putative target. To adjust the positions 
of tag count peaks, let the program empirically determine the shift 
value from the sequence data (recommended). Otherwise, if you know 
the mean size of your ChIP eluate fragments, you can use the --shiftval 
option. 

Depending on your downstream applications and/or preferences, you 
can record strict enumeration (start positions) or coverage (extend 
position).

Finally, to compare ChIP-Seq alignments from multiple experiments, 
convert your reads to Reads Per Million Mapped, which will help to 
normalize read counts.
 
 bam2wig.pl --pos start --shift --rpm --in <bamfile>
 
 bam2wig.pl --pos extend --rpm --in <bamfile>

=item Paired-end ChIP-Seq

If both ends of the ChIP eluate fragments are sequenced, then we do not 
need to calculate a shift value. Instead, we will simply count at the 
midpoint of each properly-mapped sequence pair, or record the defined 
fragment span.
 
 bam2wig.pl --pos mid --pe --rpm --in <bamfile>
 
 bam2wig.pl --pos span --pe --rpm --in <bamfile>

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
 
 bam2wig --pos span --strand --rpm --in <bamfile>

 bam2wig --pos mid --strand --rpm --in <bamfile>
 
=item Stranded, paired-end RNA-Seq

Strand presents a complication when sequencing both ends of the cDNA 
product from a library that preserves orientation. Currently, 
the TopHat aligner can handle stranded, paired-end RNA-Seq alignments. 
Because each pair will align to both strands, the aligner must record 
separately which strand the original fragment should align. The TopHat 
program records an 'XS' attribute for each alignment, and, if present, 
bam2wig.pl will use this to set the strand.
 
 bam2wig --pe --pos span --strand --rpm --in <bamfile>

 bam2wig --pe --pos mid --strand --rpm --in <bamfile>
 
=back

=head1 SHIFT VALUE DETERMINATION

To determine the shift value, the top 500 bp regions (the default
number is 200 regions) with the highest coverage are collected from
the largest chromosome or sequence represented in the Bam file. The
largest chromosome is used merely as a representative fraction of the
genome for performance reasons rather than random sampling. Stranded
read counts are collected in 10 bp bins over the region and flanking
500 bp regions (1.5 kb total). A Pearson product-moment correlation
coefficient is then reiteratively determined between the stranded data
as the bins are shifted from 30 to 500 bp. The shift corresponding to
the highest R squared value is retained for each sampled region. The
default minimum R^2 value to keep is 0.25, and not all sampled regions
may return a significant R^2 value. A trimmed mean is then calculated
from all of the calculated shift values to be used used as the final
shift value. This approach works best with clean, distinct peaks. The
peak shift may be evaluated by viewing separate, stranded wig files
together with the shifted wig file in a genome browser.

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

