#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(sum min max mean stddev);
use Statistics::LineFit;
use File::Spec;
use Bio::ToolBox::file_helper qw(
	open_to_read_fh 
	open_to_write_fh 
	write_tim_data_file
);
use Bio::ToolBox::data_helper qw(
	format_with_commas
	generate_tim_data_structure
);
use Bio::ToolBox::big_helper qw(wig_to_bigwig_conversion);
eval {
	# check for bam support
	require Bio::ToolBox::db_helper::bam;
	Bio::ToolBox::db_helper::bam->import;
};
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};

# Declare constants for this program
use constant {
	LOG2            => log(2),
	LOG10           => log(10),
};
my $VERSION = '1.14';
	
	

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
	$model,
	$strand,
	$bin_size,
	$min_mapq,
	$max_dup,
	$max_count,
	$rpm,
	$log,
	$bigwig,
	$bwapp,
	$gz,
	$cpu,
	$buffer_min,
	$alignment_count,
	$verbose,
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
	'model!'    => \$model, # write the strand shift model data
	'strand!'   => \$strand, # separate strands
	'bin=i'     => \$bin_size, # size of bin to make
	'qual=i'    => \$min_mapq, # minimum mapping quality
	'max=i'     => \$max_dup, # maximum duplicate positions
	'max_cnt=i' => \$max_count, # maximum pileup count in coverage mode only
	'rpm!'      => \$rpm, # calculate reads per million
	'log=i'     => \$log, # transform count to log scale
	'bw!'       => \$bigwig, # generate bigwig file
	'bwapp=s'   => \$bwapp, # utility to generate a bigwig file
	'gz!'       => \$gz, # compress text output
	'cpu=i'     => \$cpu, # number of cpu cores to use
	'buffer=i'  => \$buffer_min, # minimum buffer length
	'count=i'   => \$alignment_count, # number of alignments before writing to disk
	'verbose!'  => \$verbose, # print sample correlations
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
my ($use_start, $use_mid, $use_span, $bin, $bedgraph, $callback, 
	$split_callback, $write_wig, $convertor, $data_ref, $outbase);
check_defaults();

# record start time
my $start_time = time;





### Open files
# Bam file
unless (exists &open_bam_db) {
	die " unable to load Bam file support! Is Bio::DB::Sam installed?\n"; 
}
my $sam = open_bam_db($infile) or die " unable to open bam file '$infile'!\n";




### Calculate Total read numbers
my $total_read_number = 0;
if ($rpm) {
	# this is only required when calculating reads per million
	print " Pre-counting total number of aligned fragments... \n";
	$total_read_number = sum_total_bam_alignments($sam, $min_mapq, $paired, $cpu);
		# this is multi-threaded as well so pass the cpu number available
	
	# print result
	print "   ", format_with_commas($total_read_number), " total mapped fragments\n";
	printf " counted in %.1f minutes\n", (time - $start_time)/60;
}



### Calculate shift value
if ($shift) {
	if (not $shift_value) {
		print " Calculating 3' shift value...\n";
		$shift_value = determine_shift_value();
	}
	
	# precalculate double shift when recording extended position
	if ($position eq 'extend') {
		$shift_value = $shift_value * 2;
		print " Reads will be extended by $shift_value bp\n";
	}
	else {
		print " Read counts will be shifted by $shift_value bp\n";
	}
}



### Process bam file
# process according to type of data collected and alignment type
if ($cpu > 1) {
	# multiple cpu core execution
	
	if ($use_coverage) {
		# special, speedy, low-level, single-bp coverage 
		parallel_process_bam_coverage();
	}
	else {
		# process alignments individually
		if ($splice) {
			# enable splices, which will use special callbacks
			$sam->split_splices($splice);
			print " WARNING: enabling splices will increase processing times\n";
		}
		parallel_process_alignments();
	}
	
}
else {
	# single-thread execution
	if ($use_coverage) {
		# special, speedy, low-level, single-bp coverage 
		process_bam_coverage();
	}
	else {
		# process alignments individually
		if ($splice) {
			# enable splices, which will use special callbacks
			$sam->split_splices($splice);
			print " WARNING: enabling splices will increase processing times\n";
		}
		process_alignments();
	}
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
	
	# check parallel support
	if ($parallel) {
		# conservatively enable 2 cores
		$cpu ||= 2;
	}
	else {
		# disable cores
		print " disabling parallel CPU execution, no support present\n" if $cpu;
		$cpu = 0;
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
	if ($rpm and $max_dup) {
		print " WARNING: The read count for rpm normalization includes all duplicates." . 
			"\n  Please filter your bam file if you wish to limit duplicates and have\n" . 
			"  an accurate rpm normalization.\n";
	}
	
	# check minimum buffer
	unless (defined $buffer_min) {
		$buffer_min = 1200;
	}
	
	# check alignment count
	unless (defined $alignment_count) {
		$alignment_count = 200000;
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
			if ($shift) {
				# essentially the same thing as extend
				$position eq 'extend';
			}
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
		if ($strand) {
			warn " disabling strand with shift enabled\n";
			$strand = 0;
		}
	}
	if ($shift_value and !$shift) {
		undef $shift_value;
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
	(undef, undef, $outbase) = File::Spec->splitpath($outfile);
	$outfile =~ s/\.(?:wig|bdg|bedgraph)(?:\.gz)?$//i; # strip extension if present
	$outbase =~ s/\.(?:wig|bdg|bedgraph)(?:\.gz)?$//i; # splitpath doesn't do extensions
	
	
	### Determine the alignment callback method
	# and the wig writing method
	# based on the position used, strandedness, and/or shift
	if ($use_coverage) {
		# no callbacks necessary, special mode
		print " recording coverage spanning alignments\n";
	}
	elsif ($strand and $shift) {
		# this should not happen, strand is disabled with shift
		die " programming error!\n";
	}
	elsif ($use_start and $paired) {
		# this should not happen, paired start is changed to paired mid
		die " programming error!\n";
	}
	elsif ($shift and $paired) {
		# this should not happen, shift is disabled with paired
		die " programming error!\n";
	}
	elsif ($use_start and !$strand and !$shift and !$paired) {
		$callback  = \&record_start;
		$split_callback = \&record_split_start;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording start positions\n";
	}
	elsif ($use_start and $strand and !$shift and !$paired) {
		$callback  = \&record_stranded_start;
		$split_callback = \&record_split_stranded_start;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording stranded, start positions\n";
	}
	elsif ($use_start and !$strand and $shift and !$paired) {
		$callback  = \&record_shifted_start;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording shifted-start positions\n";
	}
	elsif ($use_mid and $strand and !$shift and $paired) {
		$callback = \&record_stranded_paired_mid;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording stranded, mid positions of pairs\n";
	}
	elsif ($use_mid and $strand and !$shift and !$paired) {
		$callback = \&record_stranded_mid;
		$split_callback = \&record_split_stranded_mid;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording stranded, mid positions\n";
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
		$split_callback = \&record_split_mid;
		$write_wig = $bin ? \&write_fixstep : \&write_varstep;
		print " recording mid position\n";
	}
	elsif ($use_span and $strand and !$shift and $paired) {
		$callback  = \&record_stranded_paired_span;
		$write_wig = \&write_bedgraph;
		print " recording stranded positions spanning paired alignments\n";
	}
	elsif ($use_span and $strand and !$shift and !$paired) {
		$callback  = \&record_stranded_span;
		$split_callback = \&record_split_stranded_span;
		$write_wig = \&write_bedgraph;
		print " recording stranded positions spanning alignments\n";
	}
	elsif ($use_span and !$strand and !$shift and $paired) {
		$callback  = \&record_paired_span;
		$write_wig = \&write_bedgraph;
		print " recording positions spanning paired alignments\n";
	}
	elsif ($use_span and !$strand and !$shift and !$paired) {
		if (!$rpm and !$log and !$splice) {
			# this is actually coverage
			undef $use_span;
			undef $bedgraph; # no bedgraph file, coverage uses fixedStep
			$use_coverage = 1;
			print " recording coverage spanning alignments\n";
		}
		else {
			$callback  = \&record_span;
			$split_callback = \&record_split_span;
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
	
	
	### Determine the convertor callback
	if ($rpm and !$log) {
		$convertor = sub {
			return ($_[0] * 1_000_000) / $total_read_number;
		};
	}
	elsif ($rpm and $log == 2) {
		# calculate rpm first before log
		$convertor = sub {
			return 0 if $_[0] == 0;
			return log( ($_[0] * 1_000_000) / $total_read_number) / LOG2;
		};
	}
	elsif ($rpm and $log == 10) {
		# calculate rpm first before log
		$convertor = sub {
			return 0 if $_[0] == 0;
			return log( ($_[0] * 1_000_000) / $total_read_number) / LOG10;
		};
	}
	elsif (!$rpm and $log == 2) {
		$convertor = sub {
			return 0 if $_[0] == 0;
			return log($_[0]) / LOG2;
		};
	}
	elsif (!$rpm and $log == 10) {
		$convertor = sub {
			return 0 if $_[0] == 0;
			return log($_[0]) / LOG10;
		};
	}
	elsif (!$rpm and !$log) {
		# this one is easy!
		$convertor = sub {
			return $_[0];
		};
	}
}


### Open the output file handle 
sub open_wig_file {
	my ($name, $track) = @_;
	
	# add extension to filename
	$name .= $bedgraph ? '.bdg' : '.wig';
	$name .= '.gz', if $gz;
		
	
	# open
	my $fh = open_to_write_fh($name, $gz) or 
		die " unable to open output wig file '$name'!\n";
		
	# write track line
	if ($bedgraph) {
		$fh->print("track type=bedGraph\n") if $track;
	}
	else {
		$fh->print("track type=wiggle_0\n") if $track;
	}
	
	# finished
	return ($name, $fh);
}


### Determine the shift value
sub determine_shift_value {
	
	# identify top regions to score
	# we will walk through the largest chromosome(s) looking for the top  
	# 500 bp regions containing the highest unstranded coverage to use
	print "  sampling the top $sample_number coverage regions " .
		"on the largest $chr_number chromosomes...\n";
	
	# first sort the chromosomes by size
		# this is assuming all the chromosomes have different sizes ;-)
		# use a Schwartzian transform
	my @chromosomes = 
		map { $_->[0] }
		sort { $b->[1] <=> $a->[1] }
		map { [$_, $sam->target_len($_)] }
		(0 .. $sam->n_targets - 1);
	
	# result arrays
	my @shift_values;
	my @f_profile;
	my @r_profile;
	my @shifted_profile;
	my @r2_values;
	my @regions; 
	
	# look for high coverage regions to sample
	# do this in multi-threaded fashion if possible
	
	if ($cpu > 1) {
		# do each chromosome in parallel
		print "   Forking into children for parallel scanning\n";
		
		# set up 
		my $pm = Parallel::ForkManager->new($cpu);
		$pm->run_on_finish( sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result) = @_;
			
			# record the chromosome results
			push @shift_values, @{ $result->[0] }; # push the actual values
			push @f_profile, @{ $result->[1] }; 
			push @r_profile, @{ $result->[2] };
			push @shifted_profile, @{ $result->[3] };
			push @r2_values, @{ $result->[4] };
			push @regions, @{ $result->[5] };
		} );
		
		# scan the chromosomes in parallel
		for my $chromosome_index (1 .. $chr_number) {
			$pm->start and next;
			
			### In child
			$sam->clone; # to make it fork safe
			
			# scan
			my $tid = $chromosomes[$chromosome_index - 1]; # convert to 0-based indexing
			print "   Scanning ", $sam->target_name($tid), "\n";
			my %coverage2region; 
			scan_high_coverage($tid, \%coverage2region);
			
			# calculate the correlation
			my $result = calculate_strand_correlation(\%coverage2region);
			
			$pm->finish(0, $result); 
		}
		$pm->wait_all_children;
	}
	else {
		# one chromosome at a time
		for my $chromosome_index (1 .. $chr_number) {
			
			# scan
			my $tid = $chromosomes[$chromosome_index - 1]; # convert to 0-based indexing
			print "   Scanning ", $sam->target_name($tid), "\n";
			my %coverage2region; 
			scan_high_coverage($tid, \%coverage2region);
			
			# calculate the correlation
			my $result = calculate_strand_correlation(\%coverage2region);
			
			# record the results for this chromosome
			push @shift_values, @{ $result->[0] }; # push the actual values
			push @f_profile, @{ $result->[1] }; 
			push @r_profile, @{ $result->[2] };
			push @shifted_profile, @{ $result->[3] };
			push @r2_values, @{ $result->[4] };
			push @regions, @{ $result->[5] };
		}
	}
	printf "  %s regions found with a correlative shift in %.1f minutes\n", 
		format_with_commas(scalar @shift_values), (time - $start_time)/60;
	
	# determine the optimal shift value
	# we will be using a trimmed mean value to avoid outliers
	my $raw_mean = mean(@shift_values);
	my $raw_sd = stddev(@shift_values);
	printf "  The collected mean shift value is %.0f +/- %.0f bp\n", 
		$raw_mean, $raw_sd;
	my $raw_min = $raw_mean - (1.5 * $raw_sd);
	$raw_min = 0 if $raw_min < 0;
	my $raw_max = $raw_mean + (1.5 * $raw_sd);
	my @trimmed_shift_values;
	my @trimmed_f_profile;
	my @trimmed_r_profile;
	my @trimmed_shifted_profile;
	my @trimmed_r2_values;
	my @trimmed_regions;
	foreach my $i (0 .. $#shift_values) {
		if ($shift_values[$i] >= $raw_min and $shift_values[$i] <= $raw_max) {
			push @trimmed_shift_values, $shift_values[$i];
			push @trimmed_f_profile, $f_profile[$i];
			push @trimmed_r_profile, $r_profile[$i];
			push @trimmed_shifted_profile, $shifted_profile[$i];
			push @trimmed_r2_values, $r2_values[$i];
			push @trimmed_regions, $regions[$i];
		}
	}
	my $best_value = sprintf("%.0f", mean(@trimmed_shift_values) );
	printf "  The trimmed mean shift value is %s +/- %.0f bp from %s regions\n", 
		$best_value, stddev(@trimmed_shift_values), 
		format_with_commas(scalar @trimmed_shift_values);
	
	# write out the shift model data file
	if ($model) {
		write_model_file($best_value, \@trimmed_f_profile, \@trimmed_r_profile, 
			\@trimmed_shifted_profile, \@trimmed_r2_values, \@trimmed_regions);
	}
	
	# done
	return $best_value;
}



### Scan chromosome for high coverage regions
sub scan_high_coverage {
	my ($tid, $coverage2region) = @_;
	
	my $size = $sam->target_len($tid);
	
	my $current_min = 0;
	for (my $start = 0; $start < $size; $start += 500) {
	
		my $end = $start + 500;
		$end = $size if $end > $size;
	
		# using the low level interface for a little more performance
		my $coverage = $sam->bam_index->coverage(
			$sam->bam,
			$tid, # need the low level target ID
			$start, 
			$end,
			1, # return mean coverage in a single bin
		);
		my $sum_coverage = $coverage->[0];
		next if $coverage->[0] == 0;
	
		# check if our coverage exceeds the lowest region
		if (scalar keys %{$coverage2region} < $sample_number) {
			# less than requested regions found so far, so keep it
			# record the coordinates for this region
			$coverage2region->{$coverage->[0]} = [$tid, $start, $end];
		}
		else {
			# we already have the maximum number
			
			# check that we have a current minimum value
			unless ($current_min) {
				$current_min = min( keys %{$coverage2region} );
			}
			
			# find the lowest one
			if ($coverage->[0] > $current_min) {
				# it's a new high over the lowest minimum
			
				# remove the previous lowest region
				delete $coverage2region->{ $current_min };
			
				# add the current region
				# record the coordinates for this region
				$coverage2region->{$coverage->[0]} = [$tid, $start, $end];
				$current_min = min( keys %{$coverage2region} );
			}
		}
	}

}



### Calculate the cross strand correlation for the sample regions
sub calculate_strand_correlation {
	my $coverage2region = shift;
	
	my @shift_values;
	my @all_r2_values;
	my @f_profile;
	my @r_profile;
	my @shifted_profile;
	my @regions;
	
	# determine the optimal shift for each of the test regions
	foreach my $i (keys %{$coverage2region}) {
		
		# get the start and end positions
		# we're adjusting them by 400 bp in both directions to actually 
		# sample a 1.3 kb region centered over the original region
		my $tid   = $coverage2region->{$i}->[0];
		my $start = $coverage2region->{$i}->[1] - 400;
		my $end   = $coverage2region->{$i}->[2] + 400;
		my $region = $sam->target_name($tid) . ":$start..$end";
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
		$sam->bam_index->fetch(
			$sam->bam, $tid, $start, $end, \&record_stranded_start, \%data);
		
		# generate data arrays
		my @f;
		my @r;
		for (my $i = $start; $i <= $end; $i += 10) {
			# walk in the region by increments in 10
			# sum the counts in each 10 bp interval, and push it to the array
			push @f, sum( map { $data{f}{$_} ||= 0 } ($i .. $i+9) );
			push @r, sum( map { $data{r}{$_} ||= 0 } ($i .. $i+9) );
		}
		my @original_f = @f;
		my @original_r = @r;
		my $best_r = 0;
		my $best_shift = 0;
		my @best_profile;
		my @r2_values;
		
		# calculate correlations
		for (my $i = 0; $i <= 40; $i++) {
			# check shift from 0 to 400 bp
			
			# adjust the arrays, mimicking shifting arrays towards the 3'
			if ($i) {
				unshift @f, 0;
				pop @f;
				shift @r;
				push @r, 0;
			}
			
			# calculate correlation
			my $stat = Statistics::LineFit->new();
			$stat->setData(\@r, \@f) or warn " bad data!\n";
			my $r2 = $stat->rSquared();
			push @r2_values, $r2;
			
			# check correlation
			if ($r2 >= $correlation_min and $r2 > $best_r) {
				# record new values
				$best_shift = $i * 10;
				$best_r = $r2;
				
				# record the best profile, average of f and r values
				if ($model) {
					# this is only required when reporting the model
					for my $i (0 .. 129) {
						$best_profile[$i] = mean( $f[$i], $r[$i] );
					}
				}
			}
		}
		
		# print result
		if ($verbose) {
			my $string = "  sampling $region  ";
			if ($best_r) {
				printf "$string shift $best_shift bp (r^2 %.3f)\n", $best_r;
			}
			else {
				print "$string\n";
			}
		}
		
		# record result
		if ($best_r) {
			push @shift_values, $best_shift;
			if ($model) {
				push @f_profile, \@original_f;
				push @r_profile, \@original_r;
				push @shifted_profile, \@best_profile;
				push @all_r2_values, \@r2_values;
				push @regions, $region;
			}
		}
	}
	
	return [ \@shift_values, \@f_profile, \@r_profile, 
		\@shifted_profile, \@all_r2_values, \@regions];
}


### Write a text data file with the shift model data
sub write_model_file {
	my ($value, $f_profile, $r_profile, $shifted_profile, $r2_values, $regions) = @_;
	
	### Profile model
	# Prepare the centered profiles from the raw profiles
	# these will be -450 to +450, with 0 being the shifted profile peak
	my @centered_f_profile;
	my @centered_r_profile;
	my @centered_shifted_profile;
	for my $r (0 .. $#{$regions}) {
		# first identify the peak in the shifted profile
		my $peak = 0;
		my $peak_i;
		for my $i (0 .. 129) {
			if ($shifted_profile->[$r][$i] > $peak) {
				$peak = $shifted_profile->[$r][$i];
				$peak_i = $i;
			}
		}
		
		# collect the centered profiles
		for my $i (0 .. 90) {
			my $current_i = $peak_i - 45 + $i;
			if ($current_i >= 0 and $current_i <= 129) {
				# check that the current index is within the raw index range
				$centered_f_profile[$r][$i] = $f_profile->[$r][$current_i];
				$centered_r_profile[$r][$i] = $r_profile->[$r][$current_i];
				$centered_shifted_profile[$r][$i] = $shifted_profile->[$r][$current_i];
			}
			else {
				# otherwise record zeros
				$centered_f_profile[$r][$i] = 0;
				$centered_r_profile[$r][$i] = 0;
				$centered_shifted_profile[$r][$i] = 0;
			}
		}
	}
	
	
	# Prepare the data structure
	my $profile = generate_tim_data_structure(
		'shift_model_profile',
		'Start',
		"$outbase\_F",
		"$outbase\_R",
		"$outbase\_shift"
	);
	$profile->{'db'} = $infile;
	push @{ $profile->{'other'} }, "# Average profile of read start point sums\n";
	push @{ $profile->{'other'} }, 
		"# Only profiles of trimmed shift value samples included\n";
	push @{ $profile->{'other'} }, "# Final shift value calculated as $value bp\n";
	$profile->{2}{'minimum_r2'} = $correlation_min;
	$profile->{2}{'number_of_chromosomes_sampled'} = $chr_number;
	$profile->{2}{'regions_sampled'} = scalar(@$f_profile);
	
	
	# Load data table
	# first we will put the mean value for all the regions
	for my $i (0 .. 90) {
		# generate the start position
		$profile->{'data_table'}->[$i+1][0] = ($i - 45) * 10;
		
		# generate the mean value for each chromosome tested at each position
		$profile->{'data_table'}->[$i+1][1] = mean( map { $centered_f_profile[$_][$i] } 
			(0 .. $#centered_f_profile ) );
		$profile->{'data_table'}->[$i+1][2] = mean( map { $centered_r_profile[$_][$i] } 
			(0 .. $#centered_r_profile ) );
		$profile->{'data_table'}->[$i+1][3] = 
			mean( map { $centered_shifted_profile[$_][$i] } 
			(0 .. $#centered_shifted_profile ) );
	}
	$profile->{'last_row'} = 91;
	
	
	# Add the region specific profiles if verbose if requested
	if ($verbose) {
		for my $r (0 .. $#{$regions}) {
		
			# add column specific metadata
			my $column = $profile->{'number_columns'};
			$profile->{$column} = {
				'name'  => $regions->[$r] . '_F',
				'index' => $column,
			};
			$profile->{$column + 1} = {
				'name'  => $regions->[$r] . '_R',
				'index' => $column + 1,
			};
			$profile->{$column + 2} = {
				'name'  => $regions->[$r] . '_Shift',
				'index' => $column + 2,
			};
			$profile->{'data_table'}->[0][$column]   = $profile->{$column}{'name'};
			$profile->{'data_table'}->[0][$column+1] = $profile->{$column+1}{'name'};
			$profile->{'data_table'}->[0][$column+2] = $profile->{$column+2}{'name'};
		
			# fill in the columns
			for my $i (0 .. 90) {
				$profile->{'data_table'}->[$i+1][$column]   = $centered_f_profile[$r][$i];
				$profile->{'data_table'}->[$i+1][$column+1] = $centered_r_profile[$r][$i];
				$profile->{'data_table'}->[$i+1][$column+2] = 
					$centered_shifted_profile[$r][$i];
			}
		
			$profile->{'number_columns'} += 3;
		}
	}
	
	# Write the model file
	my $profile_file = write_tim_data_file( 
		'data'     => $profile,
		'filename' => "$outfile\_model.txt",
		'gz'       => 0,
	);
	print "  Wrote shift profile model data file $profile_file\n" if $profile_file;
	
	
	### R squared data
	# prepare the data structure
	my $r2_data = generate_tim_data_structure(
		'Shift_correlations',
		'Shift',
		"$outbase\_R2"
	);
	$r2_data->{'db'} = $infile;
	push @{ $r2_data->{'other'} }, "# Average R Squared values for each shift\n";
	push @{ $r2_data->{'other'} }, "# Final shift value calculated as $value bp\n";
	$r2_data->{1}{'minimum_r2'} = $correlation_min;
	$r2_data->{1}{'number_of_chromosomes_sampled'} = $chr_number;
	$r2_data->{1}{'regions_sampled'} = scalar(@$f_profile);
	
	# load data table
	# first we will put the mean value for all the r squared values
	for my $i (0 .. 40) {
		# generate the start position
		$r2_data->{'data_table'}->[$i+1][0] = $i * 10;
		
		# generate the mean value for each chromosome
		$r2_data->{'data_table'}->[$i+1][1] = mean( map { $r2_values->[$_][$i] } 
			(0 ..  $#{$regions}) );
	}
	$r2_data->{'last_row'} = 41;
	
	# add the chromosome specific profiles
	if ($verbose) {
		for my $r (0 .. $#{$regions}) {
		
			# add column specific metadata
			my $column = $r2_data->{'number_columns'};
			$r2_data->{$column} = {
				'name'  => $regions->[$r] . '_R2',
				'index' => $column,
			};
			$r2_data->{'data_table'}->[0][$column]   = $r2_data->{$column}{'name'};
		
			# fill in the columns
			for my $i (0 .. 40) {
				$r2_data->{'data_table'}->[$i+1][$column]   = $r2_values->[$r][$i];
			}
		
			$r2_data->{'number_columns'}++;
		}
	}
	
	# write the r squared file
	my $rSquared_file = write_tim_data_file( 
		'data'     => $r2_data,
		'filename' => "$outfile\_correlations.txt",
		'gz'       => 0,
	);
	print "  Wrote shift correlation data file $rSquared_file\n" if $rSquared_file;
}



### Collect alignment coverage
sub process_bam_coverage {
	# using the low level bam coverage method, not strand specific
	
	# open wig file
	my ($filename, $fh) = open_wig_file($outfile, 1);
	
	# determine the dump size
	# the dump size indicates how much of the genome we take before we 
	# dump to file
	# this is based on the requested bin_size, default should be 1 bp 
	# but we want to keep the coverage array reasonable size, 10000 elements
	# is plenty. we'll do some math to make it a multiple of bin_size to fit
	my $dump = $bin_size * int( 10000 / $bin_size); 
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as an integer
		process_bam_coverage_on_chromosome($fh, $tid, $dump);
	}
	
	$fh->close;
	print " Wrote file $filename\n";
	
	# convert to bigwig if requested
	if ($bigwig) {
		convert_to_bigwig($filename);
	}
}



### Parallel process coverage
sub parallel_process_bam_coverage {
	# using the low level bam coverage method, not strand specific
	# parallel multiple core execution
	
	# determine the dump size
	# the dump size indicates how much of the genome we take before we 
	# dump to file
	# this is based on the requested bin_size, default should be 1 bp 
	# but we want to keep the coverage array reasonable size, 10000 elements
	# is plenty. we'll do some math to make it a multiple of bin_size to fit
	my $dump = $bin_size * int( 10000 / $bin_size); 
	
	# we don't want child files to be compressed, takes too much CPU time
	my $original_gz = $gz;
	$gz = 0;
	
	# prepare ForkManager
	print " Forking into $cpu children for parallel conversion\n";
	my $pm = Parallel::ForkManager->new($cpu);
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as an integer
		
		# run each chromosome in a separate fork
		$pm->start and next;
		
		### in child ###
		# first clone the Bam file object to make it safe for forking
		$sam->clone;
		
		# then write the chromosome coverage in separate chromosome file
		my ($filename, $fh) = open_wig_file($outfile . '#' . sprintf("%05d", $tid));
		process_bam_coverage_on_chromosome($fh, $tid, $dump);
		
		# finished with this chromosome
		$fh->close;
		$pm->finish;
	}
	$pm->wait_all_children;
	
	# merge the children files back into one output file
	print " merging separate chromosome files\n";
	my @files = glob "$outfile#*";
	die "can't find children files!\n" unless (@files);
	$gz = $original_gz;
	my ($filename, $fh) = open_wig_file($outfile, 1);
	while (@files) {
		my $file = shift @files;
		my $in = open_to_read_fh($file);
		while (<$in>) {$fh->print($_)}
		$in->close;
		unlink $file;
	}
	print " Wrote file $filename\n";
	
	# convert to bigwig if requested
	if ($bigwig) {
		convert_to_bigwig($filename);
	}
}



### Process bam coverage for a specific chromosome
sub process_bam_coverage_on_chromosome {
	my ($fh, $tid, $dump) = @_;
	
	# get sequence info
	my $seq_length = $sam->target_len($tid);
	my $seq_id = $sam->target_name($tid);
	
	# prepare definition line for fixedStep
	$fh->print(
		"fixedStep chrom=$seq_id start=1 step=$bin_size span=$bin_size\n"
	);
	
	# set maximum pileup count if requested, default is 8000
	$sam->max_pileup_cnt($max_count) if $max_count;
	
	# walk through the chromosome
	for (my $start = 0; $start < $seq_length; $start += $dump) {
		
		# the low level interface works with 0-base indexing
		my $end = $start + $dump -1;
		$end = $seq_length if $end > $seq_length;
		
		# using the low level interface for a little more performance
		my $coverage = $sam->bam_index->coverage(
			$sam->bam,
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


### Process alignments
sub process_alignments {
	
	# open wig files
	my ($filename1, $filename2, $fh1, $fh2);
	if ($strand) {
		($filename1, $fh1) = open_wig_file("$outfile\_f", 1);
		($filename2, $fh2) = open_wig_file("$outfile\_r", 1);
	}
	else {
		($filename1, $fh1) = open_wig_file($outfile, 1);
	}
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as 
		# a numeric target identifier
		# we can easily convert this to an actual sequence name
		# we will force the conversion to go one chromosome at a time
		
		# processing depends on whether we need to split splices or not
		if ($splice) {
			process_split_alignments_on_chromosome($fh1, $fh2, $tid);
		}
		else {
			process_alignments_on_chromosome($fh1, $fh2, $tid);
		}
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



### Process alignments in parallel
sub parallel_process_alignments {
	
	# we don't want child files to be compressed, takes too much CPU time
	my $original_gz = $gz;
	$gz = 0;
	
	# prepare ForkManager
	print " Forking into $cpu children for parallel conversion\n";
	my $pm = Parallel::ForkManager->new($cpu);
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as 
		# a numeric target identifier
		# run each chromosome in a separate fork
		$pm->start and next;
		
		### in child ###
		# first clone the Bam file object to make it safe for forking
		$sam->clone;
		
		# open chromosome wig files
		my ($filename1, $filename2, $fh1, $fh2);
		if ($strand) {
			($filename1, $fh1) = open_wig_file($outfile . '_f#' . sprintf("%05d", $tid));
			($filename2, $fh2) = open_wig_file($outfile . '_r#' . sprintf("%05d", $tid));
		}
		else {
			($filename1, $fh1) = open_wig_file($outfile . '#' . sprintf("%05d", $tid));
		}
	
		# processing depends on whether we need to split splices or not
		if ($splice) {
			process_split_alignments_on_chromosome($fh1, $fh2, $tid);
		}
		else {
			process_alignments_on_chromosome($fh1, $fh2, $tid);
		}
		
		# finished with this chromosome
		$fh1->close;
		$fh2->close if $fh2;
		$pm->finish;
	}
	$pm->wait_all_children;
	
	# merge the children files back into one output file
	print " merging separate chromosome files\n";
	my ($filename1, $filename2);
	if ($strand) {
		# separate stranded files
		
		my @ffiles = glob "$outfile\_f#*";
		my @rfiles = glob "$outfile\_r#*";
		die "can't find children files!\n" unless (@ffiles and @rfiles);
		
		$gz = $original_gz;
		my $fh;
		
		# combine forward files
		($filename1, $fh) = open_wig_file("$outfile\_f", 1);
		while (@ffiles) {
			my $file = shift @ffiles;
			my $in = open_to_read_fh($file);
			while (<$in>) {$fh->print($_)}
			$in->close;
			unlink $file;
		}
		$fh->close;
		print " Wrote file $filename1\n";
		
		# combine reverse files
		($filename2, $fh) = open_wig_file("$outfile\_r", 1);
		while (@rfiles) {
			my $file = shift @rfiles;
			my $in = open_to_read_fh($file);
			while (<$in>) {$fh->print($_)}
			$in->close;
			unlink $file;
		}
		$fh->close;
		print " Wrote file $filename2\n";
	}
	else {
		# single wig file
		
		my @files = glob "$outfile#*";
		die "can't find children files!\n" unless (@files);
		$gz = $original_gz;
		my $fh;
		($filename1, $fh) = open_wig_file($outfile, 1);
		while (@files) {
			my $file = shift @files;
			my $in = open_to_read_fh($file);
			while (<$in>) {$fh->print($_)}
			$in->close;
			unlink $file;
		}
		$fh->close;
		print " Wrote file $filename1\n";
	}
	
	# convert to bigwig if requested
	if ($bigwig) {
		convert_to_bigwig($filename1, $filename2);
	}
}



### Process alignments for a specific chromosome
sub process_alignments_on_chromosome {
	my ($fh1, $fh2, $tid) = @_;
	
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
	$sam->bam_index->fetch($sam->bam, $tid, 0, $seq_length, $callback, \%data);
	
	# finish up this chromosome
	&$write_wig(\%data, 1); # final write
	printf " Converted %s alignments on $seq_id in %.3f minutes\n", 
		format_with_commas( $data{'count'}), (time - $start_time)/60;
}



### Process split alignments for a specific chromosome
sub process_split_alignments_on_chromosome {
	my ($fh1, $fh2, $tid) = @_;
	
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
	
	# to make things even more complicated, the split alignment subfeatures
	# are no longer bam alignments, but essentially SeqFeature objects
	# which means we must use special callbacks to record them (sigh....)
	
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


### Callback for processing single-end split alignments
sub single_end_spliced_callback {
	my $a = shift;
	
	# check alignment
	return if $a->unmapped;
	return if $a->qual < $min_mapq;
		# we need to check the quality here instead of later
		# because the subfeatures lose this quality score
	
	# check for subfeatures
	my @subfeatures = $a->get_SeqFeatures;
	if (@subfeatures) {
		# process each subfeature
		foreach my $subf (@subfeatures) {
			# use special split callbacks since these are no longer bam alignments
			&$split_callback($subf, $data_ref);
		}
	}
	else {
		# no subfeatures found
		# treat this as a single read
		&$callback($a, $data_ref);
	}
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


### Record split stranded at start position
sub record_split_stranded_start {
	my ($part, $data) = @_;
	
	# record based on the strand
	if ($part->strand == 1) {
		# forward strand
		$data->{f}{ $part->start }++;
	}
	else {
		# reverse strand
		$data->{r}{ $part->end }++;
	}
	check_data($data);
}


### Record at shifted start position
sub record_shifted_start {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
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


### Record split at start position
sub record_split_start {
	my ($part, $data) = @_;
	
	# start based on strand, record on forward
	# these are high level objects
	if ($part->strand == 1) {
		# forward strand
		$data->{f}{ $part->start }++;
	}
	else {
		# reverse strand
		$data->{f}{ $part->end }++;
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


### Record split stranded at mid position
sub record_split_stranded_mid {
	my ($part, $data) = @_;
	
	# calculate mid position
	my $pos = int( ($part->start + $part->end) / 2);
	
	# record based on the strand
	if ($part->strand == 1) {
		# forward strand
		$data->{f}{ $pos }++;
	}
	else {
		# reverse strand
		$data->{r}{ $pos }++;
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


### Record at mid position
sub record_split_mid {
	my ($part, $data) = @_;
	
	# record mid position
	$data->{f}{ int( ($part->start + $part->end) / 2) }++;
	
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


### Record stranded across alignment
sub record_split_stranded_span {
	my ($part, $data) = @_;
	
	# record length at the start position based on strand
	if ($part->strand == 1) {
		# forward strand
		$data->{f}{ $part->start } .= $part->length . ',';
	}
	else {
		# reverse strand
		$data->{r}{ $part->start } .= $part->length . ',';
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


### Record across alignment
sub record_split_span {
	my ($part, $data) = @_;
	
	# record the length
	$data->{f}{ $part->start } .= $part->length . ',';
	check_data($data);
}


### Record extended alignment
sub record_extended {
	my ($a, $data) = @_;
	
	# check
	return if $a->qual < $min_mapq;
	
	# start based on strand, record on forward
	if ($a->reversed) {
		# reverse strand
		# must calculate the start position of the 3 prime extended fragment
		$data->{f}{ $a->calend - $shift_value } .= "$shift_value,";
	}
	else {
		# forward strand
		$data->{f}{ $a->pos } .= "$shift_value,";
	}
	check_data($data);
}


### Increment count and check the data size for writing
sub check_data {
	my $data = $_[0];
	$data->{'count'}++;
	
	# write when we reach buffer maximum number of alignments read
	if ($data->{'print'} and $data->{'count'} % $alignment_count == 0) {
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
			max(keys %{$data->{$s}}) - $buffer_min; 
		
		# deal with negative positions
		if (min( keys %{ $data->{$s} }) <= 0) {
			print "  Warning: " . (abs( min( keys %{ $data->{$s} } )) + 1) . 
				" bp trimmed from the beginning of chromosome " . $data->{'seq_id'} . "\n" 
				if $verbose ;
			
			# delete the negative positions
			foreach ( keys %{ $data->{$s} } ) {
				delete $data->{$s}{$_} if $_ <= 0;
			}
		}
		
		# write the data
		foreach my $pos (sort {$a <=> $b} keys %{$data->{$s}}) {
			
			# first check the position and bail when we've gone off the chromosome end
			last if ($pos > $maximum);
			
			# write line
			my $score = $data->{$s}{$pos};
			if ($max_dup) {
				$score = $score > $max_dup ? $max_dup : $score;
			}
			$data->{$fh}->print( join("\t", 
				$pos, 
				&$convertor($score) # convert RPM or log before writing
			) . "\n");
			
			# clean up
			delete $data->{$s}{$pos};
		}
		
		# warn about tossed values at the chromosome end
		if ($verbose and $final and keys(%{ $data->{$s} }) ) {
			print "  Warning: " . (max( keys(%{ $data->{$s} }) ) - $maximum) . 
				" bp trimmed" . " from the end of chromosome " . $data->{'seq_id'} .
				 "\n";
		}
	}
}


### Write a fixedStep wig file
sub write_fixstep {
	# only binned data is ever written as a fixedStep wig file
	my ($data, $final) = @_;
	
	# do each strand one at a time
	foreach my $s (qw(f r)) {
		next unless keys( %{ $data->{$s} } );
		my $offset = "offset$s";
		my $fh     = "fh$s";
		
		# check the maximum position that we cannot go beyond
		# defined either by the minimum buffer value or the end of the chromosome
		my $maximum = $final ?  $data->{'seq_length'} : 
			max(keys %{$data->{$s}}) - $buffer_min; 
		
		# deal with negative positions
		if (min( keys %{ $data->{$s} }) <= 0) {
			print "  Warning: " . (abs( min( keys %{ $data->{$s} } )) + 1) . 
				" bp trimmed from the beginning of chromosome " . $data->{'seq_id'} . "\n" 
				if $verbose ;
			
			# delete the negative positions
			foreach ( keys %{ $data->{$s} } ) {
				delete $data->{$s}{$_} if $_ <= 0;
			}
		}
		
		# write binned data
		for (
			my $pos = $data->{$offset}; 
			$pos < $maximum - $bin_size; 
			$pos += $bin_size
		) {
		
			# sum the counts in the bin interval
			my $score;
			if ($max_dup) {
				$score = sum(
					map { 
						my $a = $data->{$s}{$_} ||= 0;
						return $a > $max_dup ? $max_dup : $a;
					} ($pos .. $pos + $bin_size -1 )
				);
			}
			else {
				$score = sum(
					map { $data->{$s}{$_} ||= 0 } ($pos .. $pos + $bin_size -1 ) );
			}
		
			# write line
			$data->{$fh}->print( &$convertor($score) . "\n" );
		}
		
		# clean up
		for (
			my $pos = $data->{$offset}; 
			$pos < $maximum; 
			$pos++
		) {
			next unless exists $data->{$s}{$_};
			delete $data->{$s}{$pos};
		}		
		
		# warn about tossed values at the chromosome end
		if ($verbose and $final and keys(%{ $data->{$s} }) ) {
			print "  Warning: " . (max( keys(%{ $data->{$s} }) ) - $maximum) . 
				" bp trimmed" . " from the end of chromosome " . $data->{'seq_id'} .
				 "\n";
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
		
		# reset final to true as necessary
		$final = 1 if ( max(keys %{$data->{$s}}) >= $data->{'seq_length'} );
		
		# check the maximum position that we cannot go beyond
		# defined either by the minimum buffer value or the end of the chromosome
		my $maximum = $final ?  $data->{'seq_length'} : 
			max(keys %{$data->{$s}}) - $buffer_min; 
		
		# print warning about negative positions
		if ($verbose and min( keys %{ $data->{$s} }) < 0) {
			print "  Warning: " . abs( min( keys %{ $data->{$s} } )) . " bp trimmed" .
				" from the beginning of chromosome " . $data->{'seq_id'} . "\n";
		}
		
		# convert read lengths to coverage in the buffer array
		foreach my $pos (sort {$a <=> $b} keys %{ $data->{$s} }) {
			
			# first check the position and bail when we've gone far enough
			last if ($pos > $maximum);
			
			# split the lengths, limit to max, and generate coverage
			my @lengths = split(',', $data->{$s}{$pos});
			if ($max_dup and scalar @lengths > $max_dup) {
				splice(@lengths, $max_dup + 1); # delete the extra ones
			}
			
			# we will take a shortcut if all lengths are the same
			if ( $shift_value or all_equal(\@lengths) ) {
				# all the lengths are equal
				# each position will get the same pileup number
				for (0 .. $lengths[0] - 1) {
					# generate coverage
					# we're relying on autovivification of the buffer array here
					# this is what makes Perl both great and terrible at the same time
					# this could also balloon memory usage - oh dear
					my $p = $pos - $data->{$offset} + $_;
					$data->{$buffer}->[$p] += scalar(@lengths) if $p >= 0;
					# avoid modifying negative positions
				}
			}
			else {
				# not all the lengths are equal, must slog through all of them
				foreach my $len (@lengths) {
					# generate coverage same as above
					for (0 .. $len -1) { 
						my $p = $pos - $data->{$offset} + $_;
						$data->{$buffer}->[$p] += 1 if $p >= 0;
						# avoid modifying modifying negative positions
					}
				}
			}
			delete $data->{$s}{$pos};
		}
		
		# write the array into bedgraph
		my $current_pos = $data->{$offset};
		my $current_value = shift @{ $data->{$buffer} } || 0;
		my $current_offset = 0;
		while (
			scalar( @{ $data->{$buffer} } ) > $buffer_min or 
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
				# remember we're working with 0-base positions here
				my $end = $current_pos + $current_offset + 1;
				if ($final and $end >= $maximum) {
					# we've reached the end of the chromosome
					$end = $data->{'seq_length'};
					
					# dump anything left in the array to finish up
					print "  Warning: " . scalar(@{ $data->{$buffer} }) . " bp trimmed" .
						" from the end of chromosome " . $data->{'seq_id'} . "\n" 
						if $verbose;
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
			my $end = $current_pos + $current_offset + 1;
			if ($final) {
				# interval must not go off the chromosome end
				if ($current_pos < $data->{'seq_length'}) {
					
					# reset end if necessary
					$end = $end > $data->{'seq_length'} ? $data->{'seq_length'} : $end;
					
					# print the hanging value
					$data->{$fh}->print( join("\t", 
							$data->{'seq_id'}, 
							$current_pos,
							$end,
							&$convertor($current_value) # convert RPM or log before writing
					) . "\n");
				}
			}
			else {
				$data->{$fh}->print( join("\t", 
						$data->{'seq_id'}, 
						$current_pos,
						$end,
						&$convertor($current_value) # convert RPM or log before writing
				) . "\n");
			}
		}
		
		# remember the current position for next writing
		$data->{$offset} = $current_pos + $current_offset + 1;
	}
}


### A simple test to confirm if all elements in an array ref are equal
sub all_equal {
	my $array = shift;
	my $value = $array->[0];
	foreach (@$array) {
		return 0 if $_ != $value;
	}
	return 1;
}


### Convert to BigWig format
sub convert_to_bigwig {
	my @files = @_;
	
	foreach my $file (@files) {
		next unless defined $file;
		my $bw_file = wig_to_bigwig_conversion(
			'wig'   => $file,
			'db'    => $sam,
		);
		if ($bw_file) {
			print " Converted to $bw_file\n";
			unlink $file;
		}
		else {
			print " BigWig conversion failed! see standard error for details\n";
		}
	}
}


__END__

=head1 NAME

bam2wig.pl

A script to enumerate Bam alignments or coverage into a wig file.

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
  --model
  --strand
  --qual <integer>
  --max <integer>
  --max_cnt <integer>
  --rpm
  --log [2|10]
  --bw
  --bwapp </path/to/wigToBigWig or /path/to/bedGraphToBigWig>
  --gz
  --cpu <integer>
  --buffer <integer>
  --count <integer>
  --verbose
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input Bam alignment file. The file should be sorted by 
genomic position and indexed, although it may be indexed automatically.

=item --out <filename>

Specify the output base filename. An appropriate extension will be 
added automatically. By default it uses the base name of the 
input file.

=item --position [start|mid|span|extend]

Specify the position of the alignment coordinate which should be 
recorded. Several positions are accepted: 
     
    start     the 5 prime position of the alignment
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
mapped pairs of alignments will be counted. Properly mapped pairs 
include FR reads on the same chromosome, and not FF, RR, RF, or 
pairs aligning to separate chromosomes. The default is to 
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

Provide the value in bp that the recorded position should be shifted. 
The value should be 1/2 the average length of the insert library 
that was sequenced. The default is to empirically determine the 
appropriate shift value. See below for the approach.

=item --sample <integer>

Indicate the number of top regions to sample when empirically 
determining the shift value. The default is 200.

=item --chrom <integer>

Indicate the number of sequences or chromosomes to sample when 
empirically determining the shift value. The reference sequences 
listed in the Bam file header are taken in order of decreasing 
length, and one or more are taken as a representative sample of 
the genome. The default value is 2. 

=item --minr <float>

Provide the minimum R^2 value to accept a shift value when 
empirically determining the shift value. Enter a decimal value 
between 0 and 1. Higher values are more stringent. The default 
is 0.25.

=item --model

Indicate that the shift model profile data should be written to 
file for examination. The average profile, including for each 
sampled chromosome, are reported for the forward and reverse strands, 
as  well as the shifted profile. A standard text file is generated 
using the output base name. The default is to not write the model 
shift data.

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
This uses the alignment start (or midpoint if recording midpoint) position 
to determine duplicity. Note that this has no effect in coverage mode. 
You may want to set a limit when working with random fragments (sonication) 
to avoid PCR bias. Note that setting this value in conjunction with the --rpm 
option may result in lower coverage than anticipated, since the pre-count 
does not account for duplicity. The default is undefined (no limit). 

=item --max_cnt <integer>

In special coverage mode only, this option sets the maximum coverage count 
at any given base. The default is 8000 (set by the bam adaptor).

=item --rpm

Convert the data to Reads (or Fragments) Per Million mapped. This is useful 
for comparing read coverage between different datasets. Only alignments 
that match the minimum mapping quality are counted. Only proper paired-end 
alignments are counted, they are counted as one fragment. The conversion is 
applied before converting to log, if requested. This will increase processing 
time, as the alignments must first be counted. Note that all duplicate reads 
are counted during the pre-count. The default is no RPM conversion. 

=item --log [2|10]

Transform the count to a log scale. Specify the base number, 2 or 
10. Only really useful with Bam alignment files with high count numbers. 
Default is to not transform the count.

=item --bw

Specify whether or not the wig file should be further converted into 
an indexed, compressed, binary BigWig file. The default is false.

=item --bwapp < /path/to/wigToBigWig or /path/to/bedGraphToBigWig >

Specify the full path to Jim Kent's bigWig conversion utility. Two 
different utilities may be used, bedGraphToBigWig or wigToBigWig, 
depending on the format of the wig file generated. The application 
paths may be set in the biotoolbox.cfg file.

=item --gz

Specify whether (or not) the output file should be compressed with 
gzip. The default is compress the output unless a BigWig file is 
requested. Disable with --nogz.

=item --cpu <integer>

Specify the number of CPU cores to execute in parallel. This requires 
the installation of Parallel::ForkManager. With support enabled, the 
default is 2. Disable multi-threaded execution by setting to 1. 

=item --buffer <integer>

Specify the length in bp to reserve as buffer when writing a bedGraph 
file to account for future read coverage. This value must be greater 
than the expected alignment length (including split alignments), 
paired-end span length (especially RNA-Seq), or extended coverage 
(2 x alignment shift). Increasing this value may result in increased 
memory usage, but will avoid errors with duplicate positions 
written to the wig file. The default is 1200 bp. 

=item --count <integer>

Specify the number of alignments processed before writing to file. 
Increasing this count will reduce the number of disk writes and increase 
performance at the cost of increased memory usage. Lowering will 
decrease memory usage. The default is 200,000 alignments.

=item --verbose

Print warnings when read counts go off the end of the chromosomes, 
particularly with shifted read counts. Also print the correlations 
for each sampled region as they are calculated when determining the 
shift value. When writing the model file, data from each region is 
also written. The default is false.

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
in the 3 prime direction. This effectively merges the separate peaks 
(representing the ends of the enriched fragments) on each strand 
into a single peak centered over the target locus. Alternatively, 
the entire predicted fragment may be recorded across its span. 
This extended method of recording is analogous to the approach 
used by the MACS program. The shift value may be empirically 
determined from the sequencing data (see below). If requested, the 
shift model profile may be written to file. Use the BioToolBox 
script C<graph_profile.pl> to graph the data.

The output wig file may be either a variableStep, fixedStep, or 
bedGraph format. The file format is dictated by where the alignment 
position is recorded. Recording start and midpoint at single 
base-pair resolution writes a variableStep wig file. Binned start 
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
performs a 3 prime shift adjustment to center the fragment's end reads 
over the predicted center and putative target. To adjust the positions 
of tag count peaks, let the program empirically determine the shift 
value from the sequence data (recommended). Otherwise, if you know 
the mean size of your ChIP eluate fragments, you can use the --shiftval 
option. 

To evaluate the empirically determined shift value, be sure to include 
the --model option to examine the profiles of stranded and shifted read 
counts and the distribution of cross-strand correlations.

Depending on your downstream applications and/or preferences, you 
can record strict enumeration (start positions) or coverage (extend 
position).

Finally, to compare ChIP-Seq alignments from multiple experiments, 
convert your reads to Reads Per Million Mapped, which will help to 
normalize read counts.
 
 bam2wig.pl --pos start --shift --model --rpm --in <bamfile>
 
 bam2wig.pl --pos extend --model --rpm --in <bamfile>

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
 
You may also wish to increase the buffer size (see --buffer) to account 
for reads that may span very large introns and avoid recording 
duplicate positions.

=back

=head1 SHIFT VALUE DETERMINATION

To empirically determine the shift value, a cross-strand correlation 
method is employed. Regions with the highest read coverage 
are sampled from one or more chromosomes listed in the Bam file. 
The default number of regions is 200 sampled from each of the two 
largest chromosomes. The largest chromosomes are used merely as a 
representative fraction of the genome for performance reasons. Stranded
read counts are collected in 10 bp bins over a 1300 bp region (the 
initial 500 bp high coverage region plus flanking 400 bp). A Pearson 
product-moment correlation coefficient is then reiteratively determined 
between the stranded data sets as the bins are shifted from 0 to 400 bp. 
The shift corresponding to the highest R squared value is recorded for 
each sampled region. The default minimum R squared value to record an 
optimal shift is 0.25, and not all sampled regions may return a 
significant R squared value. After collection, outlier shift values 
E<gt> 1.5 standard deviations from the mean are removed, and the trimmed 
mean is used as the final shift value.

This approach works best with clean, distinct peaks, although even 
noisy data can generate a reasonably good shift model. If requested, 
a text file containing the average read count profiles for the forward 
strand, reverse strand, and shifted data are written so that a model 
graph may be generated. You can generate a visual graph of the shift 
model profiles using the following command:
 
 graph_profile.pl --skip 4 --offset 1 --in <shift_model.txt>
 
The peak shift may also be evaluated by viewing separate, stranded wig 
files together with the shifted wig file in a genome browser.

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
