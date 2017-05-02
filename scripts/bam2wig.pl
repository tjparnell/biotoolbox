#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(sum);
use Bio::ToolBox::db_helper qw(
	open_db_connection
	low_level_bam_coverage
	low_level_bam_fetch
);
use Bio::ToolBox::utility qw(
	format_with_commas
	open_to_read_fh 
	open_to_write_fh 
); 
use Bio::ToolBox::big_helper qw(
	open_wig_to_bigwig_fh 
	generate_chromosome_file
);
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};

my $VERSION = '1.50';
	
	

print "\n This program will convert bam alignments to wig data\n";

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
	$outfile,
	$use_start, 
	$use_mid, 
	$use_span, 
	$use_cspan,
	$use_extend,
	$position,
	$use_coverage,
	$splice,
	$paired,
	$shift,
	$shift_value,
	$extend_value,
	$chr_number,
	$correlation_min,
	$zmin,
	$zmax,
	$model,
	$do_strand,
	$flip,
	$min_mapq,
	$secondary,
	$duplicate,
	$supplementary,
	$max_isize,
	$min_isize,
	$multi_hit_scale,
	$rpm,
	$do_mean,
	$chr_exclude,
	$black_list,
	$dec_precison,
	$bigwig,
	$do_fixstep,
	$do_varstep,
	$do_bedgraph,
	$bwapp,
	$gz,
	$cpu,
	$max_intron,
	$window,
	$verbose,
	$help,
	$print_version,
);
my @bamfiles;
my @scale_values;

# Command line options
GetOptions( 
	'in=s'      => \@bamfiles, # one or more bam files
	'out=s'     => \$outfile, # name of output file 
	'start!'    => \$use_start, # record start point
	'mid!'      => \$use_mid, # record mid point
	'span!'     => \$use_span, # record span
	'cspan!'    => \$use_cspan, # record center span
	'extend!'   => \$use_extend, # extend read
	'coverage!' => \$use_coverage, # calculate coverage
	'position=s' => \$position, # legacy option
	'splice|split!'   => \$splice, # split splices
	'pe!'       => \$paired, # paired-end alignments
	'shift!'    => \$shift, # shift coordinates 3'
	'shiftval=i' => \$shift_value, # value to shift coordinates
	'extval=i'  => \$extend_value, # value to extend reads
	'chrom=i'   => \$chr_number, # number of chromosomes to sample
	'minr=f'    => \$correlation_min, # R minimum value for shift
	'zmin=f'    => \$zmin, # minimum z-score interval for calculating shift
	'zmax=f'    => \$zmax, # maximum z-score interval for calculating shift
	'model!'    => \$model, # write the strand shift model data
	'strand!'   => \$do_strand, # separate strands
	'flip!'     => \$flip, # flip the strands
	'qual=i'    => \$min_mapq, # minimum mapping quality
	'secondary!' => \$secondary, # take secondary alignments
	'duplicate!' => \$duplicate, # include duplicate alignments
	'supplementary!' => \$supplementary, # include supplementary alignments
	'maxsize=i' => \$max_isize, # maximum paired insert size to accept
	'minsize=i' => \$min_isize, # minimum paired insert size to accept
	'fraction!'  => \$multi_hit_scale, # scale by number of hits
	'rpm!'      => \$rpm, # calculate reads per million
	'separate!' => \$do_mean, # rpm scale separately
	'scale=f'   => \@scale_values, # user specified scale value
	'chrskip=s' => \$chr_exclude, # regex for skipping chromosomes
	'blacklist=s' => \$black_list, # file for skipping regions
	'format=i'  => \$dec_precison, # format to number of decimal positions
	'bw!'       => \$bigwig, # generate bigwig file
	'bwapp=s'   => \$bwapp, # utility to generate a bigwig file
	'bdg!'      => \$do_bedgraph, # write a bedgraph output
	'fix!'      => \$do_fixstep, # write a fixedStep output
	'var!'      => \$do_varstep, # write a varStep output
	'gz!'       => \$gz, # compress text output
	'cpu=i'     => \$cpu, # number of cpu cores to use
	'intron=i'  => \$max_intron, # maximum intron size to allow
	'window=i'  => \$window, # window size to control memory usage
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
	print " Biotoolbox script bam2wig.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}



### Check for requirements and set defaults
# more global variables
my ($main_callback, $callback, $wig_writer, $outbase, $chromo_file,
	$binpack, $buflength);
check_defaults();

# record start time
my $start_time = time;





### Open files
my @sams;
foreach (@bamfiles) {
	# this will open each bam file using the high level API
	# with the appropriate installed adapter
	my $sam = open_db_connection($_) or die " unable to open bam file '$_'!\n";
	$sam->split_splices(1) if $splice; 
	push @sams, $sam;
}
# generate the chromosome name list
	# we generate this from the first bam file only, on the assumption that they 
	# all have the same sequences
	# record the chromosome name, rather than internal tid, just in the off chance 
	# that they are not in the same order (!!!???)
my @seq_list;
my %seq_name2length;
for my $tid (0 .. $sams[0]->n_targets - 1) {
	my $chr = $sams[0]->target_name($tid);
	next if $chr_exclude and $chr =~ $chr_exclude;
	push @seq_list, $chr;
	$seq_name2length{$chr} = $sams[0]->target_len($tid);
}
# set the wrapper reference
# this depends on which adapter was opened
my $wrapper_ref;
if ($splice) {
	$wrapper_ref = ref($sams[0]) eq 'Bio::DB::Sam' ? 'Bio::DB::Bam::AlignWrapper' : 
		'Bio::DB::HTS::AlignWrapper';
	eval { require $wrapper_ref; 1 };
}


### Process user provided black lists
my $black_list_hash = process_black_list();


### Calculate shift value
if ($shift or $use_extend) {
	unless ((shift and $shift_value) or ($use_extend and $extend_value)) {
		print " Calculating 3' shift value...\n";
		$shift_value = determine_shift_value();
	}
	
	# precalculate double shift when recording extended position
	if ($use_extend) {
		unless ($extend_value) {
			$extend_value = $shift_value * 2;
		}
		print " Alignments will be extended by $extend_value bp\n";
	}
	else {
		print " Alignments will be shifted by $shift_value bp\n";
	}
}
# set another global value
my $half_extend = int($extend_value / 2);


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
		parallel_process_alignments();
	}
	
}
else {
	# single-thread execution
	print " Install Parallel::ForkManager to significantly increase performance\n" 
		if not $cpu;
	if ($use_coverage) {
		# special, speedy, low-level, single-bp coverage 
		process_bam_coverage();
	}
	else {
		# process alignments individually
		process_alignments();
	}
}



### Finish
unlink $chromo_file if $chromo_file;
printf " Finished in %.3f min\n", (time - $start_time)/60;




########################   Subroutines   ###################################

### check required command line options and assign default values
sub check_defaults {
	# checking default and required values from the command line options
	# moved here to make it cleaner
	
	# check input file
	unless (@bamfiles) {
		die " no input files! use --help for more information\n" unless @ARGV;
		@bamfiles = @ARGV;
	}
	if ($bamfiles[0] =~ /,/) {
		@bamfiles = split /,/, shift @bamfiles;
	}
	
	# check scale values
	if (@scale_values) {
		if ($scale_values[0] =~ /,/) {
			@scale_values = split /,/, shift @scale_values;
		}
		if (scalar(@bamfiles) != scalar(@scale_values)) {
			if (scalar(@scale_values) == 1) {
				# only one, that's ok, use it for all
				my $v = shift @scale_values;
				push @scale_values, $v foreach (@bamfiles);
			}
			else {
				die " number of scale values does not equal number of bam files!\n";
			}
		}
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

	# set missing options if the value was set
	if ($shift_value and !$shift) {
		$shift = 1;
	}
	if ($extend_value and !$use_extend) {
		$use_extend = 1 unless ($use_cspan or $position);
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
		elsif ($position eq 'cspan') {
			$use_cspan = 1;
		}
		elsif ($position eq 'extend') {
			$use_extend = 1;
		}
		elsif ($position eq 'coverage') { 
			# for backwards compatibility
			$use_coverage = 1;
		}
		else {
			die " unrecognized position value '$position'! see help\n";
		}
	}
	my $position_check = $use_start + $use_mid + $use_span + $use_cspan + $use_extend + 
		$use_coverage;
	if ( $position_check > 1) {
		die " Modes are mutually exclusive! Please select only one of\n" . 
			" --start, --mid, --span, --cpsan, --extend, or --coverage\n";
	}
	elsif (not $position_check) {
		die " Please select one of the following modes:\n" . 
			" --start, --mid, --span, --cspan, --extend, or --coverage\n";
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
	$max_intron ||= 50000;
	
	# incompatible options
	if ($paired and $use_start) {
		warn " using start with paired-end reads is not recommended, recording midpoint\n";
		undef $use_start;
		$use_mid   = 1;
	}
	if ($splice and ($use_extend or $extend_value)) {
		warn " disabling splices when extend is defined\n";
		$splice = 0;
	}
	if ($shift and $splice) {
		die " enabling both shift and splices is currently not allowed. Pick one.\n";
	}
	if ($shift and $paired) {
		warn " disabling shift with paired reads\n";
		undef $shift;
	}
	
	# check to shift position or not
	if ($shift or $use_extend) {
		# set parameters for calculating shift value 
		unless ($shift_value or $extend_value) {
			# not provided by user, empirical calculation required
			eval {
				# required for calculating shift
				require Statistics::Descriptive;
			};
			die " Provide a shift value or install the Perl module Statistics::Descriptive\n"
				. " to empirically determine the shift value.\n"
				if $@;
		}
		
		if (defined $correlation_min) {
			if ($correlation_min <= 0 or $correlation_min >= 1) {
				die " cannot use minimum correlation value of $correlation_min!\n" .
					" use --help for more information\n";
			}
		}
		else {
			$correlation_min = 0.5;
		}
		$chr_number ||= 4;
		$zmin ||= 3;
		$zmax ||= 10;
	}
	
	# check mapping quality
	if (defined $min_mapq) {
		die " quality score must be 0-255!\n" if $min_mapq > 255;
	}
	else {
		$min_mapq = 0;
	}
	
	# check paired-end insert size
	unless (defined $max_isize) {
		$max_isize = 600;
	}
	unless (defined $max_isize) {
		$max_isize = 30;
	}
	
	# check flag parameters
	unless (defined $secondary) {
		$secondary = 1;
	}
	unless (defined $supplementary) {
		$supplementary = 1;
	}
	unless (defined $duplicate) {
		$duplicate = 1;
	}
	
	# set decimal formatting
	unless (defined $dec_precison) {
		$dec_precison = 4 if ($rpm or $multi_hit_scale or @scale_values);
	}
	
	# determine binary file packing and length
	if ($multi_hit_scale) {
		# pack as floating point values
		$binpack = 'f';
	}
	else {
		# dealing only with integers here
		# let's hope we never have depth greater than 65,536
		$binpack = 'S';
	}
	$buflength = length(pack($binpack, 1));
	
	# set window length for processing through packed binary chromosome strings
	$window ||= 10000;
		# empirical tests show a window of 1000 to 10000 is best, bigger or smaller 
		# result in longer execution times
	
	
	# check output file
	unless ($outfile) {
		if (scalar @bamfiles == 1) {
			$outfile = $bamfiles[0];
			$outfile =~ s/\.bam$//;
		}
		else {
			die " Please define an output filename when providing multiple bam files!\n";
		}
	}
	$outbase = $outfile;
	$outbase =~ s/\.(?:wig|bdg|bedgraph|bw|bigwig)(?:\.gz)?$//i; # strip extension if present
	
	
	# determine output format
	unless ($do_bedgraph or $do_varstep or $do_fixstep) {
		# pick an appropriate format for the user
		if ($use_span or $use_extend or $use_cspan) {
			$do_bedgraph = 1;
		}
		elsif ($use_start or $use_mid) {
			$do_varstep = 1;
		}
		elsif ($use_coverage) {
			$do_fixstep = 1;
		}
		else {
			die " No recording mode defined to set wig type!\n";
		}
	}
	if ( ($do_bedgraph + $do_varstep + $do_fixstep) > 1) {
		die " Please select only one of --bedgraph, --fixstep, or --varstep\n";
	}
	
	# set wig writer method
	if ($do_bedgraph) {
		$wig_writer = \&write_bedgraph;
	}
	elsif ($do_fixstep) {
		$wig_writer = \&write_fixstep;
	}
	elsif ($do_varstep) {
		$wig_writer = \&write_varstep;
	}
	
	# set the initial main callback for processing alignments
	$main_callback = $paired ? \&pe_callback : \&se_callback;
	
	
	### Determine the alignment recording callback method
	# coverage
	if ($use_coverage) {
		# does not use a callback subroutine
		undef $callback;
	}
	# start
	elsif (not $paired and not $do_strand and not $shift and $use_start) {
		$callback = \&se_start;
	}
	elsif (not $paired and not $do_strand and $shift and $use_start) {
		$callback = \&se_shift_start;
	}
	elsif (not $paired and $do_strand and not $shift and $use_start) {
		$callback = \&se_strand_start;
	}
	elsif (not $paired and $do_strand and $shift and $use_start) {
		$callback = \&se_shift_strand_start;
	}
	# midpoint
	elsif (not $paired and not $do_strand and not $shift and $use_mid) {
		$callback = \&se_mid;
	}
	elsif (not $paired and not $do_strand and $shift and $use_mid) {
		$callback = \&se_shift_mid;
	}
	elsif (not $paired and $do_strand and not $shift and $use_mid) {
		$callback = \&se_strand_mid;
	}
	elsif (not $paired and $do_strand and $shift and $use_mid) {
		$callback = \&se_shift_strand_mid;
	}
	# span
	elsif (not $paired and not $do_strand and not $shift and $use_span) {
		$callback = \&se_span;
	}
	elsif (not $paired and not $do_strand and $shift and $use_span) {
		$callback = \&se_shift_span;
	}
	elsif (not $paired and $do_strand and not $shift and $use_span) {
		$callback = \&se_strand_span;
	}
	elsif (not $paired and $do_strand and $shift and $use_span) {
		$callback = \&se_shift_strand_span;
	}
	# center span
	elsif (not $paired and not $do_strand and not $shift and $use_cspan) {
		$callback = \&se_center_span;
	}
	elsif (not $paired and not $do_strand and $shift and $use_cspan) {
		$callback = \&se_shift_center_span;
	}
	elsif (not $paired and $do_strand and not $shift and $use_cspan) {
		$callback = \&se_strand_center_span;
	}
	elsif (not $paired and $do_strand and $shift and $use_cspan) {
		$callback = \&se_shift_strand_center_span;
	}
	# extend
	elsif (not $paired and not $do_strand and not $shift and $use_extend) {
		$callback = \&se_extend;
	}
	elsif (not $paired and not $do_strand and $shift and $use_extend) {
		$callback = \&se_shift_extend;
	}
	elsif (not $paired and $do_strand and not $shift and $use_extend) {
		$callback = \&se_strand_extend;
	}
	elsif (not $paired and $do_strand and $shift and $use_extend) {
		$callback = \&se_shift_strand_extend;
	}
	# paired-end midpoint
	elsif ($paired and not $do_strand and $use_mid) {
		$callback = \&pe_mid;
	}
	elsif ($paired and $do_strand and $use_mid) {
		$callback = \&pe_strand_mid;
	}
	# paired-end span
	elsif ($paired and not $do_strand and $use_span) {
		$callback = \&pe_span;
	}
	elsif ($paired and $do_strand and $use_span) {
		$callback = \&pe_strand_span;
	}
	# paired-end span
	elsif ($paired and not $do_strand and $use_cspan) {
		$callback = \&pe_center_span;
	}
	elsif ($paired and $do_strand and $use_cspan) {
		$callback = \&pe_strand_center_span;
	}
	else {
		die "programmer error!\n";
	}
}

sub process_black_list {
	if ($black_list) {
		eval {require 'Bio::ToolBox::Data';1;} or 
			die "unable to load Bio::ToolBox::Data!!!\n";
		eval {require 'Set::IntervalTree';1;} or 
			do {
				warn " PROBLEM! Please install Set::IntervalTree to use black lists\n";
				undef $black_list;
			};
		my %black_list_hash = map { $_ => [] } @seq_list;
		my $Data = Bio::ToolBox::Data->new(file => $black_list) or 
			die "unable to read black list file '$black_list'\n";
		$Data->iterate( sub {
			my $row = shift;
			push @{ $black_list_hash{ $row->seq_id } }, 
				[ $row->start - 1, $row->end ]
				if exists $black_list_hash{ $row->seq_id };
		} );
		return \%black_list_hash;
	}
	return;
}

### Determine the shift value
sub determine_shift_value {
	
	# identify top regions to score
	# we will walk through the largest chromosome(s) looking for the top  
	# 500 bp regions containing the highest unstranded coverage to use
	print "  sampling the top coverage regions " .
		"on the largest $chr_number chromosomes...\n";
	
	# first sort the chromosomes by size
		# this is assuming all the chromosomes have different sizes ;-)
		# use a Schwartzian transform
	my $sam = $sams[0];
	my @chromosomes = 
		map { $_->[0] }
		sort { $b->[1] <=> $a->[1] }
		map { [$_, $sam->target_len($_)] }
		(0 .. $sam->n_targets - 1);
	@chromosomes = splice(@chromosomes, 0, $chr_number);
	
	# result arrays
	my @shift_values;
	my @f_profile;
	my @r_profile;
	my @shifted_profile;
	my @r_values;
	my @regions; 
	
	# look for high coverage regions to sample
	# do this in multi-threaded fashion if possible
	
	if ($cpu > 1) {
		# do each chromosome in parallel
		print "   Forking into children for parallel scanning\n";
		printf "   Scanning %s\n", join(", ", map { $sam->target_name($_) } @chromosomes);
		print "SeqID\t#Intervals\tMean\tMax\tStdDev\tMinCutoff\tMaxCutoff\t#Collected\n"
			if $verbose;
		
		# set up 
		my $pm = Parallel::ForkManager->new($cpu);
		$pm->run_on_finish( sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result) = @_;
			
			# record the chromosome results
			push @shift_values, @{ $result->[0] }; # push the actual values
			push @f_profile, @{ $result->[1] }; 
			push @r_profile, @{ $result->[2] };
			push @shifted_profile, @{ $result->[3] };
			push @r_values, @{ $result->[4] };
			push @regions, @{ $result->[5] };
		} );
		
		# scan the chromosomes in parallel
		for my $tid (@chromosomes) {
			$pm->start and next;
			
			### In child
			$sam->clone; # to make it fork safe
			
			# calculate the correlation
			my $result = calculate_strand_correlation($sam, $tid);
			
			$pm->finish(0, $result); 
		}
		$pm->wait_all_children;
	}
	else {
		# one chromosome at a time
		printf "   Scanning %s\n", join(", ", map { $sam->target_name($_) } @chromosomes);
		print "SeqID\t#Intervals\tMean\tMax\tStdDev\tMinCutoff\tMaxCutoff\t#Collected\n"
			if $verbose;
		for my $tid (@chromosomes) {
			
			# calculate the correlation
			my $result = calculate_strand_correlation($tid);
			
			# record the results for this chromosome
			push @shift_values, @{ $result->[0] }; # push the actual values
			push @f_profile, @{ $result->[1] }; 
			push @r_profile, @{ $result->[2] };
			push @shifted_profile, @{ $result->[3] };
			push @r_values, @{ $result->[4] };
			push @regions, @{ $result->[5] };
		}
	}
	printf "  %s regions found with a correlative shift in %.1f minutes\n", 
		format_with_commas(scalar @shift_values), (time - $start_time)/60;
	
	# determine the optimal shift value
	# we will be using a trimmed mean value to avoid outliers
	my $stat = Statistics::Descriptive::Sparse->new;
	$stat->add_data(@shift_values);
	printf "  The collected mean shift value is %.0f +/- %.0f bp\n", 
		$stat->mean, $stat->standard_deviation;
	my $raw_min = $stat->mean - (1.5 * $stat->standard_deviation);
	$raw_min = 0 if $raw_min < 0;
	my $raw_max = $stat->mean + (1.5 * $stat->standard_deviation);
	my @trimmed_shift_values;
	my @trimmed_f_profile;
	my @trimmed_r_profile;
	my @trimmed_shifted_profile;
	my @trimmed_r_values;
	my @trimmed_regions;
	foreach my $i (0 .. $#shift_values) {
		if ($shift_values[$i] >= $raw_min and $shift_values[$i] <= $raw_max) {
			push @trimmed_shift_values, $shift_values[$i];
			push @trimmed_f_profile, $f_profile[$i];
			push @trimmed_r_profile, $r_profile[$i];
			push @trimmed_shifted_profile, $shifted_profile[$i];
			push @trimmed_r_values, $r_values[$i];
			push @trimmed_regions, $regions[$i];
		}
	}
	$stat->clear;
	$stat->add_data(@trimmed_shift_values);
	my $best_value = sprintf("%.0f", $stat->mean);
	printf "  The trimmed mean shift value is %s +/- %.0f bp from %s regions\n", 
		$best_value, $stat->standard_deviation, 
		format_with_commas($stat->count);
	
	# write out the shift model data file
	if ($model) {
		write_model_file($best_value, \@trimmed_f_profile, \@trimmed_r_profile, 
			\@trimmed_shifted_profile, \@trimmed_r_values, \@trimmed_regions);
	}
	
	# done
	return $best_value;
}


sub calculate_strand_correlation {
	my ($sam, $tid) = @_;
	
	my ($collected, $data) = scan_high_coverage($sam, $tid);
	
	my @shift_values;
	my @all_r_values;
	my @f_profile;
	my @r_profile;
	my @shifted_profile;
	my @regions;
	
	# determine the optimal shift for each of the test regions
	foreach my $pos (@$collected) {
		
		# generate region string
		my $region = sprintf "%s:%d..%d", $sam->target_name($tid), $pos * 10, 
			($pos + 50) * 10;
		
		# grab the stranded scores from the collected binned point data 
		my $start = $pos - 50;
		$start = 0 if $start < 0;
		my $stop = $pos + 99;
		$stop = scalar( @{$data->{f}} ) if $stop > scalar( @{$data->{f}} );
		my @f = map { $data->{f}->[$_] || 0 } ($start .. $stop);
		my @r = map { $data->{r}->[$_] || 0 } ($start .. $stop);
		my @original_f = @f;
		my @original_r = @r;
		my $best_r = 0;
		my $best_shift = 0;
		my @best_profile;
		my @r_values;
		
		# calculate correlations
		for (my $i = 0; $i <= 50; $i++) {
			# check shift from 0 to 500 bp
			
			# adjust the arrays, mimicking shifting arrays towards the 3'
			if ($i) {
				unshift @f, 0;
				pop @f;
				shift @r;
				push @r, 0;
			}
			
			# calculate correlation
			my $stat = Statistics::Descriptive::Full->new();
			$stat->add_data(@f);
			my ($q, $m, $r, $rms) = $stat->least_squares_fit(@r);
				# this may produce errors when all values are equal
				# as might happen with high duplicate coverage 
			
			# check correlation
			push @r_values, $r;
			if ($r >= $correlation_min and $r > $best_r) {
				# record new values
				$best_shift = $i * 10;
				$best_r = $r;
				
				# record the best profile, average of f and r values
				if ($model) {
					# this is only required when reporting the model
					for my $i (0 .. 129) {
						$best_profile[$i] = ($f[$i] + $r[$i]) / 2;
					}
				}
			}
		}
		
		# record result
		if ($best_r > $correlation_min) {
			push @shift_values, $best_shift;
			if ($model) {
				push @f_profile, \@original_f;
				push @r_profile, \@original_r;
				push @shifted_profile, \@best_profile;
				push @all_r_values, \@r_values;
				push @regions, $region;
			}
		}
	}
	
	return [ \@shift_values, \@f_profile, \@r_profile, 
		\@shifted_profile, \@all_r_values, \@regions];
}

sub scan_high_coverage {
	my ($sam, $tid) = @_;
	
	# walk entire length of chromosome collecting alignment start data in bins of 10 bp
	my $chr_length = $sam->target_len($tid);
	my %data = (
		f => [],
		r => [],
	);
	low_level_bam_fetch($sam, $tid, 0, $chr_length, \&shift_value_callback, \%data);

	# score 500 bp intervals for coverage
	my %pos2depth;
	for (my $start = 0; $start < int($chr_length/10); $start += 50) {
		my $readsum = sum( map { $data{f}->[$_] || 0 } ($start .. $start + 49) );
		my $readsum += sum( map { $data{r}->[$_] || 0 } ($start .. $start + 49) );
		next if $readsum == 0;
		$pos2depth{$start} = $readsum;
	}
	
	# use Statistics::Descriptive to determine the mean, stddev, and then 
	# loop through the pos2depth hash and select out those +2 z-scores
	# maybe exclude the highest ones? +4 or +5?
	my $stat = Statistics::Descriptive::Sparse->new();
	$stat->add_data(values %pos2depth);
	my $count = $stat->count;
	my $mean = $stat->mean;
	my $max = $stat->max;
	my $sd = $stat->standard_deviation;
	my $mincutoff = $mean + ($zmin * $sd);
	my $maxcutoff = $mean + ($zmax * $sd);
	
	# filter out the best ones
	my @collected;
	foreach my $p (keys %pos2depth) {
		push @collected, $p if $pos2depth{$p} > $mincutoff and $pos2depth{$p} < $maxcutoff;
	}
	if ($verbose) {
		printf "%s\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%d\n", $sam->target_name($tid), 
			$count, $mean, $max, $sd, $mincutoff, $maxcutoff, scalar(@collected); 
	}
	
	return (\@collected, \%data);
}


sub shift_value_callback {
	my ($a, $data) = @_;
	
	# check alignment quality and flags
	return if ($min_mapq and $a->qual < $min_mapq); # mapping quality
	my $flag = $a->flag;
	return if ($flag & 0x0100); # not primary hit
	return if ($flag & 0x0200); # QC failed but still aligned? is this necessary?
	return if ($flag & 0x0400); # marked duplicate
	return if ($flag & 0x0800); # supplementary hit
	
	# record stranded start in 10 bp bins
	if ($a->reversed) {
		my $p = int( $a->calend / 10);
		$data->{r}->[$p] += 1;
	}
	else {
		my $p = int( $a->pos / 10);
		$data->{f}->[$p] += 1;
	}
}


### Write a text data file with the shift model data
sub write_model_file {
	my ($value, $f_profile, $r_profile, $shifted_profile, $r_valuess, $regions) = @_;
	
	# check Data
	my $data_good;
	eval {use Bio::ToolBox::Data; $data_good = 1;};
	unless ($data_good) {
		warn "unable to write model files! Cannot load Data library!\n";
		return;
	}
	
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
	my $profile = Bio::ToolBox::Data->new(
		feature     => 'shift_model_profile',
		datasets    => ['Start', "$outfile\_F", "$outfile\_R", "$outfile\_shift"],
	);
	return unless $profile;
	$profile->add_comment("Average profile of read start point sums");
	$profile->add_comment("Only profiles of trimmed shift value samples included");
	$profile->add_comment("Final shift value calculated as $value bp");
	$profile->metadata(2, 'minimum_r', $correlation_min);
	$profile->metadata(2, 'number_of_chromosomes_sampled', $chr_number);
	$profile->metadata(2, 'regions_sampled', scalar(@$f_profile));
	
	
	# Load data table
	# first we will put the mean value for all the regions
	for my $i (0 .. 90) {
		my $s = ($i - 45) * 10;
		my $f = mean( map { $centered_f_profile[$_][$i] } (0 .. $#centered_f_profile ) );
		my $r = mean( map { $centered_r_profile[$_][$i] } (0 .. $#centered_r_profile ) );
		my $m = mean( map { $centered_shifted_profile[$_][$i] } 
				(0 .. $#centered_shifted_profile ) );
		$profile->add_row( [$s, $f, $r, $m] );
	}
	
	# Write the model file
	my $profile_file = $profile->write_file( 
		'filename' => "$outfile\_model.txt",
		'gz'       => 0,
	);
	print "  Wrote shift profile model data file $profile_file\n" if $profile_file;
	
	
	### R squared data
	# prepare the data structure
	my $Data = Bio::ToolBox::Data->new(
		feature   => 'Shift_correlations',
		datasets  => ['Shift', "$outfile\_R"],
	);
	$Data->add_comment("Average correlation values for each shift");
	$Data->add_comment("Final shift value calculated as $value bp");
	$Data->metadata(1, 'minimum_r', $correlation_min);
	$Data->metadata(1, 'number_of_chromosomes_sampled', $chr_number);
	$Data->metadata(1, 'regions_sampled', scalar(@$f_profile));
	
	# load data table
	# first we will put the mean value for all the r squared values
	for my $i (0 .. 50) {
		# generate the start position
		my $s = $i * 10;
		# generate the mean value for each chromosome
		my $m = mean( map { $r_valuess->[$_][$i] } (0 ..  $#{$regions}) );
		$Data->add_row( [$s, $m] );
	}
	
	# write the r squared file
	my $success = $Data->write_file( 
		'filename' => "$outfile\_correlations.txt",
		'gz'       => 0,
	);
	print "  Wrote shift correlation data file $success\n" if $success;
}


sub mean {
	return sum(@_) / scalar(@_);
}

sub open_wig_file {
	my ($name, $do_bw) = @_;
	
	# open a bigWig file handle if requested
	if ($bigwig and $do_bw) {
		print " Writing directly to bigWig converter\n";
		$name .= '.bw' unless $name =~ /\.bw$/;
		$chromo_file = generate_chromosome_file($sams[0]);
		my $fh = open_wig_to_bigwig_fh(
			file      => $name,
			chromo    => $chromo_file,
			bwapppath => $bwapp,
		);
		if ($fh) {
			return ($name, $fh);
		}
		else {
			# we couldn't open a wigToBigWig filehandle for some reason
			# so default to standard wig file
			# typically a failure of wigToBigWig will bring the entire process
			# down, so this may not be all that necessary
			print " unable to open a bigWig file, writing standard wig file\n";
			$name =~ s/\.bw$//;
		}
	}
	
	# otherwise we open a text wig file
	unless ($name =~ /\.(?:bdg|wig)(?:\.gz)?$/i) {
		$name .= $do_bedgraph ? '.bdg' : '.wig';
	}
	$name .= '.gz' if ($gz and $name !~ /\.gz$/i);
	my $fh = open_to_write_fh($name, $gz) or 
		die " unable to open output wig file '$name'!\n";
		
	# finished
	return ($name, $fh);
}



### Collect alignment coverage
sub process_bam_coverage {
	
	# flag to write temporary binary files or go straight to wig files
	my $do_temp_bin = scalar(@sams) > 1 ? 1 : 0;
	
	# walk through each bam file and chromosome one a time 
	for my $samid (0 .. $#sams) {
		foreach my $seq_id (@seq_list) {
			process_bam_coverage_on_chromosome($samid, $sams[$samid], $seq_id, 
				$do_temp_bin);
		}
	}
	
	# find and merge binary files
	if ($do_temp_bin) {
		print " Merging temporary files from each bam\n";
		
		# find the children
		my %files;
		my @filelist = glob "$outbase.*.temp.bin";
		die " unable to find children files!\n" unless @filelist;
		foreach my $file (@filelist) {
			# each file name is basename.samid.seqid.count.strand.bin.gz
			if ($file =~ /$outbase\.(\d+)\.(.+)\.0\.f\.temp\.bin\Z/) {
				my $samid = $1;
				my $seqid = $2;
				$files{$seqid}{$samid} = $file;
			}
		}
		
		# check to see how we combine the multiple source sam file coverages
		my @norms;
		if ($do_mean) {
			# just need to take an average
			foreach (@sams) { push @norms, (1/scalar(@sams)) }
		}
		# otherwise we just add the samples 
		
		# merge the samples
		foreach my $seq_id (@seq_list) {
			merge_bin_files($seq_id, 'f', 0, $files{$seq_id}, \@norms);
		}
	}
	
	# write final wig file
	return write_final_wig_file();
}



### Parallel process coverage
sub parallel_process_bam_coverage {
	
	# flag to write temporary binary files or go straight to wig files
	my $do_temp_bin = scalar(@sams) > 1 ? 1 : 0;
	
	# generate pool of to do items 
	my @pool;
	for my $s (0 .. $#sams) {
		foreach my $seq_id (@seq_list) {
			push @pool, [$s, $sams[$s], $seq_id];
		}
	}
		
	# prepare ForkManager for working on all chromosomes and bam files in parallel
	print " Forking into $cpu children for parallel conversion\n";
	my $pm = Parallel::ForkManager->new($cpu);
	foreach my $stuff (@pool) {
		$pm->start and next;
		my ($samid, $sam, $seq_id) = @$stuff;
		
		# clone the sam object for safe forking
		$sam->clone;
		process_bam_coverage_on_chromosome($samid, $sam, $seq_id, $do_temp_bin);
		$pm->finish;
	} 
	$pm->wait_all_children;
	
	
	# find and merge binary files
	# we can do this in parallel too!
	if ($do_temp_bin) {
		print " Merging temporary files from each bam\n";
		
		# find the children
		my %files;
		my @filelist = glob "$outbase.*.temp.bin";
		die " unable to find children files!\n" unless @filelist;
		foreach my $file (@filelist) {
			# each file name is basename.samid.seqid.count.strand.bin.gz
			if ($file =~ /$outbase\.(\d+)\.(.+)\.0\.f\.temp\.bin\Z/) {
				my $samid = $1;
				my $seqid = $2;
				$files{$seqid}{$samid} = $file;
			}
		}
		
		# check to see how we combine the multiple source sam file coverages
		my @norms;
		if ($do_mean) {
			# just need to take an average
			foreach (@sams) { push @norms, (1/scalar(@sams)) }
		}
		# otherwise we just add the samples 
		
		# merge the samples
		foreach my $seq_id (@seq_list) {
			$pm->start and next;
			merge_bin_files($seq_id, 'f', 0, $files{$seq_id}, \@norms);
			$pm->finish;
		}
		$pm->wait_all_children;
	}
	
	# merge and write final wig file
	return write_final_wig_file();
}

sub process_bam_coverage_on_chromosome {
	my ($samid, $sam, $seq_id, $do_temp_bin) = @_;
	my $chr_start_time = time; 
	
	# identify chromosome target id
	my ($tid,undef,undef) = $sam->header->parse_region($seq_id);
	my $seq_length = $seq_name2length{$seq_id};
	
	# walk through the chromosome in 1 kb increments
	my $chrom_data;
	for (my $start = 0; $start < $seq_length; $start += 1000) {
		# set endpoint
		my $end = $start + 1000;
		$end = $seq_length if $end > $seq_length;
		
		# using the low level interface for a little more performance
		my $coverage = low_level_bam_coverage($sam, $tid, $start, $end);
		
		# record the coverage
		$chrom_data .= pack("$binpack*", @$coverage);
	}
	
	# write out chromosome binary file, set count to arbitrary 0
	if ($do_temp_bin) {
		write_bin_file($chrom_data, 
			join('.', $outbase, $samid, $seq_id, 0, 'f', 'temp.bin') );
	}
	else {
		&$wig_writer($chrom_data, $binpack, $seq_id, $seq_length, 
			join('.', $outbase, $samid, $seq_id, 0, 'f', 'temp.wig') );
	}
	
	# verbose status line 
	if ($verbose) {
		printf " Generated read coverage on $seq_id in %d seconds\n", 
			time - $chr_start_time;
	}
}

sub process_alignments {
	
	# flag to write temporary binary files or go straight to wig files
	my $do_temp_bin = scalar(@sams) > 1 ? 1 : $rpm ? 1 : @scale_values ? 1 : 0;
	
	# walk through each bam file and chromosome one a time 
	for my $samid (0 .. $#sams) {
		foreach my $seq_id (@seq_list) {
			process_alignments_on_chromosome($samid, $sams[$samid], $seq_id, $do_temp_bin);
		}
	}
	
	if ($verbose) {
		printf " Finished converting alignments in %.3f minutes\n", 
			(time - $start_time) / 60;
	}
	
	# find and merge binary files
	# we can do this in parallel too!
	if ($do_temp_bin) {
		
		# find the children
		my @totals;
		my %seq_totals;
		my %files;
		my @filelist = glob "$outbase.*.temp.bin";
		die " unable to find children files!\n" unless @filelist;
		foreach my $file (@filelist) {
			# each file name is basename.samid.seqid.count.strand.bin.gz
			if ($file =~ /$outbase\.(\d+)\.(.+)\.(\d+)\.([fr])\.temp\.bin\Z/) {
				my $samid = $1;
				my $seq_id = $2;
				my $count = $3;
				my $strand = $4;
				$totals[$samid] += $count;
				$seq_totals{$seq_id} += $count;
				$files{$seq_id}{$strand}{$samid} = $file;
			}
		}
		
		# merging multiple files
		if (scalar(@sams) > 1) {
			print " Merging temporary files from each bam\n" if (scalar(@sams) > 1);
			if ($rpm) {
				print " Normalizing depth\n" if $rpm;
				for my $i (0 .. $#totals) {
					printf "  %s had %s total counted alignments\n", $bamfiles[$i], 
						format_with_commas($totals[$i]);
				}
			}
			# check to see how we combine the multiple source sam file coverages
			my @norms;
			if ($rpm and $do_mean) {
				# we need to scale each sam source individually and take an average
				# calculate normalization factor for each sam file
				@norms = map { 1_000_000 / ($_ * scalar(@sams)) } @totals;
			}
			elsif ($rpm and not $do_mean) {
				# scale them all the same
				my $factor = 1_000_000 / ( sum(@totals) * scalar(@sams) );
				foreach (@sams) {
					push @norms, $factor;
				}
			}
			elsif (not $rpm and $do_mean) {
				# just need to take an average
				foreach (@sams) { push @norms, (1/scalar(@sams)) }
			}
			# otherwise we just add the samples 
		
			# add scale values on top as necessary
			if (@scale_values) {
				for my $i (0 .. $#norms) {
					$norms[$i] *= $scale_values[$i];
				}
			}
			
			# merge the samples
			foreach my $seq_id (@seq_list) {
				foreach my $strand (qw(f r)) {
					next unless defined $files{$seq_id}{$strand};
					merge_bin_files($seq_id, $strand, $seq_totals{$seq_id}, 
						$files{$seq_id}{$strand}, \@norms);
				}
			}
			
			if ($verbose) {
				printf " Finished merging%s in %.3f minutes\n", 
					defined $norms[0] ? " and normalizing" : "",
					(time - $start_time) / 60;
			}
		}
		
		# otherwise we just normalize
		else {
			# determine scaling factor
			my $scale_factor;
			if ($rpm) {
				printf " Normalizing depth based on %s total counted alignments\n",
					format_with_commas($totals[0]); 
				$scale_factor = 1_000_000 / $totals[0];
			}
			if (scalar @scale_values) {
				# user supplied scaling factor
				print " Scaling depth with user-supplied factor\n";
				if ($scale_factor) {
					$scale_factor *= $scale_values[0];
				}
				else {
					$scale_factor = $scale_values[0];
				}
			}
			
			# normalize the wig files
			foreach my $seq_id (@seq_list) {
				foreach my $strand (qw(f r)) {
					next unless defined $files{$seq_id}{$strand};
					normalize_wig_file($files{$seq_id}{$strand}{0}, $scale_factor, $seq_id);
				}
			}
		}
	}
	
	
	# write final wig file
	print " Merging temporary files\n";
	return write_final_wig_file();
}

sub parallel_process_alignments {
	# generate pool of to do items 
	my @pool;
	for my $s (0 .. $#sams) {
		foreach my $seq_id (@seq_list) {
			push @pool, [$s, $sams[$s], $seq_id];
		}
	}
	
	# flag to write temporary binary files or go straight to wig files
	my $do_temp_bin = scalar(@sams) > 1 ? 1 : $rpm ? 1 : @scale_values ? 1 : 0;
	
	# prepare ForkManager for working on all chromosomes and bam files in parallel
	print " Forking into $cpu children for parallel conversion\n";
	my $pm = Parallel::ForkManager->new($cpu) or 
		die "unable to initialize ForkManager object!\n";
	foreach my $stuff (@pool) {
		$pm->start and next;
		my ($samid, $sam, $seq_id) = @$stuff;
		
		# clone the sam object for safe forking
		$sam->clone;
		process_alignments_on_chromosome($samid, $sam, $seq_id, $do_temp_bin);
		$pm->finish;
	} 
	$pm->wait_all_children;
	
	if ($verbose) {
		printf " Finished converting alignments in %.3f minutes\n", 
			(time - $start_time) / 60;
	}
	
	# find and merge binary files
	# we can do this in parallel too!
	if ($do_temp_bin) {
		
		# find the children
		my @totals;
		my %seq_totals;
		my %files;
		my @filelist = glob "$outbase.*.temp.bin";
		die " unable to find children files!\n" unless @filelist;
		foreach my $file (@filelist) {
			# each file name is basename.samid.seqid.count.strand.bin.gz
			if ($file =~ /$outbase\.(\d+)\.(.+)\.(\d+)\.([fr])\.temp\.bin\Z/) {
				my $samid = $1;
				my $seq_id = $2;
				my $count = $3;
				my $strand = $4;
				$totals[$samid] += $count;
				$seq_totals{$seq_id} += $count;
				$files{$seq_id}{$strand}{$samid} = $file;
			}
		}
		
		# merging multiple files
		if (scalar(@sams) > 1) {
			# print processing statements
			print " Merging temporary files from each bam\n" if (scalar(@sams) > 1);
			if ($rpm) {
				print " Normalizing depth\n" if $rpm;
				for my $i (0 .. $#totals) {
					printf "  %s had %s total counted alignments\n", $bamfiles[$i], 
						format_with_commas($totals[$i]);
				}
			}
			
			# check to see how we combine the multiple source sam file coverages
			my @norms;
			if ($rpm and $do_mean) {
				# we need to scale each sam source individually and take an average
				# calculate normalization factor for each sam file
				@norms = map { 1_000_000 / ($_ * scalar(@sams)) } @totals;
				
			}
			elsif ($rpm and not $do_mean) {
				# scale them all the same
				my $factor = 1_000_000 / ( sum(@totals) * scalar(@sams) );
				foreach (@sams) {
					push @norms, $factor;
				}
			}
			elsif (not $rpm and $do_mean) {
				# just need to take an average
				foreach (@sams) { push @norms, (1/scalar(@sams)) }
			}
			# otherwise we just add the samples 
		
			# add scale values on top as necessary
			if (@scale_values) {
				for my $i (0 .. $#norms) {
					$norms[$i] *= $scale_values[$i];
				}
			}
			
			# merge the samples
			foreach my $seq_id (@seq_list) {
				foreach my $strand (qw(f r)) {
					next unless defined $files{$seq_id}{$strand};
					$pm->start and next;
					merge_bin_files($seq_id, $strand, $seq_totals{$seq_id}, 
						$files{$seq_id}{$strand}, \@norms);
					$pm->finish;
				}
			}
			$pm->wait_all_children;
			
			if ($verbose) {
				printf " Finished merging%s in %.3f minutes\n", 
					defined $norms[0] ? " and normalizing" : "",
					(time - $start_time) / 60;
			}
		}
		
		# otherwise we just normalize
		else {
			# determine scaling factor
			my $scale_factor;
			if ($rpm) {
				printf " Normalizing depth based on %s total counted alignments\n",
					format_with_commas($totals[0]); 
				$scale_factor = 1_000_000 / $totals[0];
			}
			if (scalar @scale_values) {
				# user supplied scaling factor
				print " Scaling depth with user-supplied factor\n";
				if ($scale_factor) {
					$scale_factor *= $scale_values[0];
				}
				else {
					$scale_factor = $scale_values[0];
				}
			}
			
			# normalize the wig files
			foreach my $seq_id (@seq_list) {
				foreach my $strand (qw(f r)) {
					next unless defined $files{$seq_id}{$strand};
					$pm->start and next;
					normalize_wig_file($files{$seq_id}{$strand}{0}, $scale_factor, $seq_id);
					$pm->finish;
				}
			}
			$pm->wait_all_children;
		}
	}
	
	
	# write final wig file
	print " Merging temporary files\n";
	return write_final_wig_file();
}

sub process_alignments_on_chromosome {
	my ($samid, $sam, $seq_id, $do_temp_bin) = @_;
	my $chr_start_time = time; 
	
	# identify chromosome target id
	my ($tid,undef,undef) = $sam->header->parse_region($seq_id);
	my $seq_length = $seq_name2length{$seq_id};
	
	# generate data structure for callback
	my $data = {
		f               => [],
		r               => [],
		fpack           => undef,
		rpack           => undef,
		f_offset        => 0,
		r_offset        => 0,
		pair            => {},
		black_list      => undef,
		count           => 0,
		sam             => $sam,
	};
	
	# process black lists for this chromosome
	# since we're using external interval tree module that is not fork-safe, must 
	# recreate interval tree each time
	if ($black_list_hash and scalar @{ $black_list_hash->{$seq_id} }) {
		my $tree = Set::IntervalTree->new;
		foreach (@{ $black_list_hash->{$seq_id} }) {
			# don't need to insert any particular value, just want the interval
			$tree->insert(1, $_->[0], $_->[1]);
		}
		$data->{black_list} = $tree;
	}
	
	# Process alignments on this chromosome
	low_level_bam_fetch($sam, $tid, 0, $seq_length, $main_callback, $data);
	
	# pack remaining data
	my $fdiff = $seq_length - $data->{f_offset};
	$data->{fpack} .= pack("$binpack$fdiff", @{$data->{f}});
	if ($do_strand) {
		my $rdiff = $seq_length - $data->{r_offset};
		$data->{rpack} .= pack("$binpack*", @{$data->{r}});
	}
	
	# write out file
	# we always write the forward strand, and reverse strand if stranded data
	if ($do_temp_bin) {
		# write a temporary binary file for merging later
		write_bin_file($data->{fpack}, 
			join('.', $outbase, $samid, $seq_id, $data->{count}, 'f', 'temp.bin') );
		write_bin_file($data->{rpack}, 
			join('.', $outbase, $samid, $seq_id, $data->{count}, 'r', 'temp.bin') ) 
			if $do_strand;
	}
	else {
		# write a chromosome specific wig file
		&$wig_writer($data->{fpack}, $binpack, $seq_id, $seq_length, 
			join('.', $outbase, $samid, $seq_id, $data->{count}, 'f', 'temp.wig') );
		&$wig_writer($data->{rpack}, $binpack, $seq_id, $seq_length, 
			join('.', $outbase, $samid, $seq_id, $data->{count}, 'r', 'temp.wig') ) 
			if $do_strand;
	}
	
	# verbose status line 
	if ($verbose) {
		printf "  Converted %s alignments on $seq_id in %d seconds\n", 
			format_with_commas( $data->{count}), time - $chr_start_time;
		if ($paired and keys %{$data->{pair}}) {
			printf "   %d orphan paired alignments were left behind and not counted!\n", 
				scalar keys %{$data->{pair}};
		}
	}
}

sub write_bin_file {
	my ($data, $filename) = @_;
	my $fh = open_to_write_fh($filename) or 
		die " unable to write temporary file '$filename'!\n";
	$fh->binmode;
	$fh->print($data);
	$fh->close;
}

sub merge_bin_files {
	my ($seq_id, $strand, $total, $files, $norm_factors) = @_;
	my $merge_start_time = time;
	my $long_window = 100 * $window;
	
	# open filehandles to each binary file
	my %fhs;
	foreach my $samid (keys %$files) {
		my $fh = open_to_read_fh($files->{$samid}) or 
			die sprintf " unable to read temporary file %s!\n", $files->{$samid};
		$fh->binmode;
		$fhs{$samid} = $fh;
	}
	my $first_norm = $norm_factors->[0] || undef;
	
	# march along chromosome in defined windows to keep memory usage down
	# apply normalization as data is loaded into the combined chrom_data array
	my $chrom_data; 
	for (my $pos = 0; $pos < $seq_name2length{$seq_id}; $pos += $long_window) {
		# check length
		my $len = ($pos + $window) > $seq_name2length{$seq_id} ? 
					($seq_name2length{$seq_id} - $pos) : $long_window;
		my $binary_len = $len * $buflength;
		
		# collect the data from the first file handle
		my @win_data;
		if (defined $first_norm) {
			$fhs{0}->read(my $string, $binary_len);
			@win_data = map {$_ * $first_norm} unpack("$binpack*", $string);
		}
		else {
			$fhs{0}->read(my $string, $binary_len);
			@win_data = unpack("$binpack*", $string);
		}
		
		# add to the current window the remaining filehandle data
		foreach my $samid (1 .. $#sams) {
			my $norm = $norm_factors->[$samid] || undef;
			
			# read and unpack current section from binary file
			$fhs{$samid}->read(my $string, $binary_len);
			if (defined $norm) {
				my @data = unpack("$binpack*", $string);
				for my $i (0 .. $#data) {
					$win_data[$i] += ($data[$i] * $norm) if $data[$i]; 
				}
			}
			else {
				my @data = unpack("$binpack*", $string);
				for my $i (0 .. $#data) {
					$win_data[$i] += $data[$i] if $data[$i]; 
				}
			}
		}
		
		# final pack
		$chrom_data .= pack("f$len", @win_data);
	}
	
	# clean up
	foreach my $samid (keys %fhs) {
		$fhs{$samid}->close;
		unlink $files->{$samid};
	}
	
	# now rewrite the merged bin file
	&$wig_writer($chrom_data, 'f', $seq_id, $seq_name2length{$seq_id}, 
		join('.', $outbase, '0', $seq_id, $total, $strand, 'temp.wig') );
	if ($verbose) {
		printf "  Merged%s $seq_id temp files in %d seconds\n", 
			defined $norm_factors->[0] ? " and normalized" : "", time - $merge_start_time;
	}
}

sub normalize_wig_file {
	my ($file, $scale_factor, $seq_id) = @_;
	my $norm_start_time = time;
	my $long_window = 10 * $window;
	
	# open file
	my $fh = open_to_read_fh($file) or 
		die " unable to read temporary $file!\n";
	$fh->binmode;
	
	# march along chromosome in defined windows to keep memory usage down
	# apply normalization as data is loaded into the combined chrom_data array
	my $chrom_data; 
	for (my $pos = 0; $pos < $seq_name2length{$seq_id}; $pos += $long_window) {
		# check length
		my $len = ($pos + $long_window) > $seq_name2length{$seq_id} ? 
					($seq_name2length{$seq_id} - $pos) : $long_window;
		
		# read, unpack, normalize, and re-pack current window from binary file
		my $string;
		$fh->read($string, $len * $buflength);
		$chrom_data .= pack("f*", map {$_ * $scale_factor} unpack("$binpack*", $string));
	}
	
	unlink $file;
	# now write the wig file
	$file =~ s/\.bin$/.wig/;
	&$wig_writer($chrom_data, 'f', $seq_id, $seq_name2length{$seq_id}, $file);
	if ($verbose) {
		printf "  Scaled $seq_id temp file in %d seconds\n", time - $norm_start_time;
	}
}

sub write_final_wig_file {
	
	# find children files
	my @filelist = glob "$outbase.*.temp.wig";
	unless (scalar @filelist) {
		die " can't find children file!\n";
	}
	
	# assemble into a hash
	my %files;
	foreach my $file (@filelist) {
		# each file name is basename.samid.seqid.count.strand.temp.wig.gz
		if ($file =~ /$outbase\.\d+\.(.+)\.\d+\.([fr])\.temp\.wig\Z/) {
			my $seq_id = $1;
			my $strand = $2;
			$files{$seq_id}{$strand} = $file;
		}
	}
	my @f_filelist = map { $files{$_}{f} } @seq_list;
	my @r_filelist = map { $files{$_}{r} } @seq_list;
	
	# write wig files with the appropriate wig writer
	if ($do_strand and !$flip) {
		if ($cpu > 1) {
			# we can fork this!!!!
			my $pm = Parallel::ForkManager->new(2);
			for my $i (1 .. 2) {
				$pm->start and next;
				merge_wig_files("$outbase\_f", @f_filelist) if $i == 1;
				merge_wig_files("$outbase\_r", @r_filelist) if $i == 2;
				$pm->finish;
			}
			$pm->wait_all_children;
		}
		else {
			merge_wig_files("$outbase\_f", @f_filelist);
			merge_wig_files("$outbase\_r", @r_filelist);
		}
	}
	elsif ($do_strand and $flip) {
		if ($cpu > 1) {
			# we can fork this!!!!
			my $pm = Parallel::ForkManager->new(2);
			for my $i (1..2) {
				$pm->start and next;
				merge_wig_files("$outbase\_r", @f_filelist) if $i == 1;
				merge_wig_files("$outbase\_f", @r_filelist) if $i == 2;
				$pm->finish;
			}
			$pm->wait_all_children;
		}
		else {
			merge_wig_files("$outbase\_r", @f_filelist);
			merge_wig_files("$outbase\_f", @r_filelist);
		}
	}
	else {
		merge_wig_files($outbase, @f_filelist);
	}
}

sub write_bedgraph {
	my ($data, $packer, $seq_id, $seq_length, $filename) = @_;	
	
	# set the printf formatter for decimal or integer
	my $formatter = $dec_precison ? 
		"$seq_id\t%d\t%d\t%." . $dec_precison. "f\n" : "$seq_id\t%d\t%d\t%s\n";
	
	# work though chromosome
	my $buflength = length(pack($packer, 1));
	my $out_string;
	my $cpos = 0; # current position
	my $lpos = 0; # last position
	my $cval = 0; # current value
	for (my $pos = 0; $pos < $seq_length; $pos += $window) {
		# check length
		my $len = ($pos + $window) > $seq_length ? ($seq_length - $pos) : $window;
		
		# unpack current window from the passed binary string
		my @win_data = unpack("$packer*", 
			substr($data, $pos * $buflength, $len * $buflength));
		
		# work through current window
		foreach my $value (@win_data) {
			if ($value == $cval) {
				$cpos++;
			}
			else {
				$out_string .= sprintf($formatter, $lpos, $cpos, $cval);
				$lpos = $cpos;
				$cval = $value;
				$cpos++;
			}
		}
	}
	
	# final write
	if ($cpos > $lpos) {
		$out_string .= sprintf($formatter, $lpos, $cpos, $cval);
	}
	
	# write wig file
	my ($filename, $outfh) = open_wig_file($filename, 0);
	$outfh->print($out_string);
	$outfh->close;
}

sub write_fixstep {
	my ($data, $packer, $seq_id, $seq_length, $filename) = @_;	

	# set the printf formatter for decimal or integer
	my $formatter = $dec_precison ? "%." . $dec_precison . "f\n" : "%s\n";
	
	# write fixStep header
	my $out_string = "fixedStep chrom=$seq_id start=1 step=1 span=1\n";
	
	# work though chromosome
	my $buflength = length(pack($packer, 1));
	for (my $pos = 0; $pos < $seq_length; $pos += $window) {
		# check length
		my $len = ($pos + $window) > $seq_length ? ($seq_length - $pos) : $window;
		
		# unpack current window from the passed binary string
		my @win_data = unpack("$packer*", 
			substr($data, $pos * $buflength, $len * $buflength));
		
		# work through current window
		foreach my $value (@win_data) {
			$out_string .= sprintf($formatter, $value);
		}
	}
	
	# write wig file
	my ($filename, $outfh) = open_wig_file($filename, 0);
	$outfh->print($out_string);
	$outfh->close;
}

sub write_varstep {
	my ($data, $packer, $seq_id, $seq_length, $filename) = @_;	

	# set the printf formatter for decimal or integer
	my $formatter = $dec_precison ? "%d\t%." . $dec_precison. "f\n" : "%d\t%s\n";
	
	# write fixStep header
	my $out_string = "variableStep chrom=$seq_id\n";
	
	# work though chromosome
	my $buflength = length(pack($packer, 1));
	for (my $pos = 0; $pos < $seq_length; $pos += $window) {
		# check length
		my $len = ($pos + $window) > $seq_length ? ($seq_length - $pos) : $window;
		
		# unpack current window from the passed binary string
		my @win_data = unpack("$packer*", 
			substr($data, $pos * $buflength, $len * $buflength));
		
		# work through current window
		for (my $i = 0; $i <= $#win_data; $i++) {
			next unless $win_data[$i]; # should be a non-zero value
			$out_string .= sprintf($formatter, $i + $pos + 1, $win_data[$i]);
		}
	}
	
	# write wig file
	my ($filename, $outfh) = open_wig_file($filename, 0);
	$outfh->print($out_string);
	$outfh->close;
}

sub merge_wig_files {
	my ($outfile, @files) = @_;
	
	my ($filename1, $fh) = open_wig_file($outfile, 1);
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

sub se_callback {
	my ($a, $data) = @_;
	
	# check alignment quality and flags
	return if ($min_mapq and $a->qual < $min_mapq); # mapping quality
	my $flag = $a->flag;
	return if (not $secondary and $flag & 0x0100); # secondary alignment
	return if (not $duplicate and $flag & 0x0400); # marked duplicate
	return if ($flag & 0x0200); # QC failed but still aligned? is this necessary?
	return if (not $supplementary and $flag & 0x0800); # supplementary hit
	
	# filter black listed regions
	if (defined $data->{black_list}) {
		my $results = $data->{black_list}->fetch($a->pos, $a->calend);
		return if @$results;
	}
	
	# scale by number of hits
	my $score = 1;
	if ($multi_hit_scale) {
		my $nh = $a->aux_get('NH') || 1;
		$score = $nh > 1 ? 1/$nh : 1;
	}
	$data->{count} += $score;
	
	# pass checks
	if ($splice and $a->cigar_str =~ /N/) {
		# check for splices
		se_spliced_callback($a, $data, $score);
	}
	else {
		&$callback($a, $data, $score);
	}
	
	# check data size
	if (scalar(@{$data->{f}}) > 1_000_000) {
		$data->{fpack} .= pack("$binpack*", splice(@{$data->{f}}, 0, 500_000));
		$data->{f_offset} += 500_000;
	}
	if (scalar(@{$data->{r}}) > 1_000_000) {
		$data->{rpack} .= pack("$binpack*", splice(@{$data->{r}}, 0, 500_000));
		$data->{r_offset} += 500_000;
	}
}

sub pe_callback {
	my ($a, $data) = @_;
	
	# check paired status
	return unless $a->proper_pair;
	return unless $a->tid == $a->mtid; # same chromosome, redundant with proper_pair?
	return if $a->isize > $max_isize;
	return if $a->isize < $min_isize;
	
	# check alignment quality and flags
	return if ($min_mapq and $a->qual < $min_mapq); # mapping quality
	my $flag = $a->flag;
	return if (not $secondary and $flag & 0x0100); # secondary alignment
	return if (not $duplicate and $flag & 0x0400); # marked duplicate
	return if ($flag & 0x0200); # QC failed but still aligned? is this necessary?
	return if (not $supplementary and $flag & 0x0800); # supplementary hit
	
	# filter black listed regions
	if ($data->{black_list}) {
		my $results;
		if ($a->reversed) {
			$results = $data->{black_list}->fetch($a->calend - $a->isize, $a->calend);
		}
		else {
			$results = $data->{black_list}->fetch($a->pos, $a->pos + $a->isize);
		}
		return if @$results;
	}
	
	# look for pair
	my $f; # the forward alignment
	if ($a->reversed) {
		$f = $data->{pair}->{$a->qname} || undef;
		return unless $f; # no pair, no record
		delete $data->{pair}->{$a->qname};
		
		# scale by number of hits
		my $score = 1;
		if ($multi_hit_scale) {
			my $r_nh = $a->aux_get('NH') || 1;
			my $f_nh = $f->aux_get('NH') || 1;
			if ($f_nh < $r_nh) {
				# take the lowest number of hits recorded
				$score = $f_nh > 1 ? 1/$f_nh : 1;
			}
			else {
				$score = $r_nh > 1 ? 1/$r_nh : 1;
			}
		}
		
		# record based only on the forward read
		$data->{count}++;
		&$callback($f, $data, $score);
	}
	else {
		# store until we find it's mate
		$data->{pair}->{$a->qname} = $a;
	}
	
	# check data size
	if (scalar(@{$data->{f}}) > 600_000) {
		$data->{fpack} .= pack("$binpack*", splice(@{$data->{f}}, 0, 500_000));
		$data->{f_offset} += 500_000;
	}
	if (scalar(@{$data->{r}}) > 600_000) {
		$data->{rpack} .= pack("$binpack*", splice(@{$data->{r}}, 0, 500_000));
		$data->{r_offset} += 500_000;
	}
}

sub se_spliced_callback {
	my ($a, $data, $score) = @_;
	my $aw = $wrapper_ref->new($a, $data->{sam});
	my @segments = $aw->get_SeqFeatures;
	my $size = 1;
	my $cigars = $aw->cigar_array;
	foreach my $c (@$cigars) {
		$size = $c->[1] if ($c->[0] eq 'N' and $c->[1] > $size);
	}
	return if $size > $max_intron; # exceed maximum intron size
	if ($use_start) {
		foreach my $segment (@segments) {
			if ($do_strand and $a->reversed) {
				# reverse strand
				$data->{r}->[$segment->start - 1 - $data->{r_offset}] += $score;
			}
			else {
				# otherwise forward strand
				$data->{f}->[$segment->start - 1 - $data->{f_offset}] += $score;
			}
		}
	}
	elsif ($use_span) {
		foreach my $segment (@segments) {
			if ($do_strand and $a->reversed) {
				# reverse strand
				for ($segment->start - 1 - $data->{r_offset} .. 
					$segment->end - $data->{r_offset}
				) {
					$data->{r}->[$_] += $score;
				}
			}
			else {
				# otherwise forward strand
				for ($segment->start - 1  - $data->{f_offset} .. 
					$segment->end - $data->{f_offset}
				) {
					$data->{f}->[$_] += $score;
				}
			}
		}
	}
}

sub se_start {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		$data->{f}->[$a->calend - 1 - $data->{f_offset}] += $score;
	}
	else {
		$data->{f}->[$a->pos - $data->{f_offset}] += $score;
	}
}

sub se_shift_start {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $pos = $a->calend - 1 - $shift_value - $data->{f_offset};
		$data->{f}->[$pos] += $score if $pos >= 0;
	}
	else {
		$data->{f}->[$a->pos + $shift_value - $data->{f_offset}] += $score;
	}
}

sub se_strand_start {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		$data->{r}->[$a->calend - 1 - $data->{r_offset}] += $score;
	}
	else {
		$data->{f}->[$a->pos - $data->{f_offset}] += $score;
	}
}

sub se_shift_strand_start {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $pos = $a->calend -1 - $shift_value - $data->{r_offset};
		$data->{r}->[$pos] += $score if $pos >= 0;
	}
	else {
		$data->{f}->[$a->pos + $shift_value - $data->{f_offset}] += $score;
	}
}

sub se_mid {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	$data->{f}->[$mid - $data->{f_offset}] += $score;
}

sub se_shift_mid {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	if ($a->reversed) {
		my $pos = $mid - $shift_value;
		$data->{f}->[$pos - $data->{f_offset}] += $score if $pos >= 0;
	}
	else {
		$data->{f}->[$mid + $shift_value - $data->{f_offset}] += $score;
	}
}

sub se_strand_mid {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	if ($a->reversed) {
		$data->{r}->[$mid - $data->{r_offset}] += $score;
	}
	else {
		$data->{f}->[$mid - $data->{f_offset}] += $score;
	}
}

sub se_shift_strand_mid {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	if ($a->reversed) {
		my $pos = $mid - $shift_value;
		$data->{r}->[$pos - $data->{r_offset}] += $score if $pos >= 0;
	}
	else {
		$data->{f}->[$mid + $shift_value - $data->{f_offset}] += $score;
	}
}

sub se_span {
	my ($a, $data, $score) = @_;
	foreach ($a->pos - $data->{f_offset} .. $a->calend - 1 - $data->{f_offset}) {
		$data->{f}->[$_] += $score;
	}
}

sub se_strand_span {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		foreach ($a->pos - $data->{r_offset} .. $a->calend - 1 - $data->{r_offset}) {
			$data->{r}->[$_] += $score;
		}
	}
	else {
		foreach ($a->pos - $data->{f_offset} .. $a->calend - 1 - $data->{f_offset}) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_shift_span {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		foreach ($a->pos - $shift_value - $data->{r_offset} .. 
			$a->calend - 1 - $shift_value - $data->{r_offset} 
		) {
			$data->{f}->[$_] += $score if $_ >= 0;
		}
	}
	else {
		foreach ($a->pos + $shift_value - $data->{f_offset} .. 
			$a->calend - 1 + $shift_value - $data->{f_offset}
		) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_shift_strand_span {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		foreach ($a->pos - $shift_value - $data->{r_offset} .. 
			$a->calend - 1 - $shift_value - $data->{r_offset}
		) {
			$data->{r}->[$_] += $score if $_ >= 0;
		}
	}
	else {
		foreach ($a->pos + $shift_value - $data->{f_offset} .. 
			$a->calend - 1 + $shift_value - $data->{f_offset}
		) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_center_span {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	my $start = $mid - $half_extend + 1;
	$start = 0 if $start < 0;
	foreach ($start - $data->{f_offset} .. $mid + $half_extend - $data->{f_offset}) {
		$data->{f}->[$_] += $score;
	}
}

sub se_strand_center_span {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	my $start = $mid - $half_extend + 1;
	$start = 0 if $start < 0;
	if ($a->reversed) {
		foreach ($start - $data->{r_offset} .. 
			$mid + $half_extend - $data->{r_offset}
		) {
			$data->{r}->[$_] += $score;
		}
	}
	else {
		foreach ($mid - $half_extend + 1 - $data->{f_offset} .. 
			$mid + $half_extend - $data->{f_offset}
		) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_shift_center_span{
}

sub se_shift_strand_center_span{
}

sub se_extend {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $start = $a->calend - $extend_value;
		$start = 0 if $start < 0; # avoid negative positions in array
		foreach ($start - $data->{f_offset} .. $a->calend - 1 - $data->{f_offset}) {
			$data->{f}->[$_] += $score;
		}
	}
	else {
		foreach ($a->pos - $data->{f_offset} ..
			 $a->pos + $extend_value - 1 - $data->{f_offset}
		) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_strand_extend {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $start = $a->calend - $extend_value;
		$start = 0 if $start < 0; # avoid negative positions in array
		foreach ($start - $data->{r_offset} .. $a->calend - 1 - $data->{r_offset}) {
			$data->{r}->[$_] += $score;
		}
	}
	else {
		foreach ($a->pos - $data->{f_offset} .. 
			$a->pos + $extend_value - 1 - $data->{f_offset}
		) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_shift_extend {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $start = $a->calend - $shift_value - $extend_value;
		$start = 0 if $start < 0; # avoid negative positions in array
		foreach ($start - $data->{f_offset} .. 
			$a->calend - 1 - $shift_value - $data->{f_offset}
		) {
			$data->{f}->[$_] += $score;
		}
	}
	else {
		foreach ($a->pos + $shift_value - $data->{f_offset} .. 
			$a->pos + $shift_value + $extend_value - 1 - $data->{f_offset}
		) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_shift_strand_extend {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $start = $a->calend - $shift_value - $extend_value;
		$start = 0 if $start < 0; # avoid negative positions in array
		foreach ($start - $data->{r_offset} .. 
			$a->calend - 1 - $shift_value - $data->{r_offset}
		) {
			$data->{r}->[$_] += $score;
		}
	}
	else {
		foreach ($a->pos + $shift_value - $data->{f_offset} .. 
			$a->pos + $shift_value + $extend_value - 1 - $data->{f_offset}
		) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub pe_start {
}

sub pe_strand_start {
}

sub pe_mid {
	my ($a, $data, $score) = @_;
	my $mid = $a->pos + int( $a->isize / 2 );
	$data->{f}->[$mid - $data->{f_offset}] += $score;
}

sub pe_strand_mid {
	my ($a, $data, $score) = @_;
	my $mid = $a->pos + int( $a->isize / 2 );
	my $flag = $a->flag;
	if ($flag & 0x0040) {
		# first read
		if ($a->reversed) {
			$data->{r}->[$mid - $data->{r_offset}] += $score;
		}
		else {
			$data->{f}->[$mid - $data->{f_offset}] += $score;
		}
	}
	elsif ($flag & 0x0080) {
		# second read
		if ($a->reversed) {
			$data->{f}->[$mid - $data->{f_offset}] += $score;
		}
		else {
			$data->{r}->[$mid - $data->{r_offset}] += $score;
		}
	} 
}

sub pe_span {
	my ($a, $data, $score) = @_;
	foreach ($a->pos  - $data->{f_offset} .. $a->pos  - $data->{f_offset} + $a->isize) {
		$data->{f}->[$_] += $score;
	}
}

sub pe_strand_span {
	my ($a, $data, $score) = @_;
	my $flag = $a->flag;
	if ($flag & 0x0040) {
		# first read
		if ($a->reversed) {
			foreach ($a->pos - $data->{r_offset} .. 
				$a->pos - $data->{r_offset} + $a->isize
			) {
				$data->{r}->[$_] += $score;
			}
		}
		else {
			foreach ($a->pos - $data->{f_offset} .. 
				$a->pos - $data->{f_offset} + $a->isize
			) {
				$data->{f}->[$_] += $score;
			}
		}
	}
	elsif ($flag & 0x0080) {
		# second read
		if ($a->reversed) {
			foreach ($a->pos - $data->{f_offset} .. 
				$a->pos - $data->{f_offset} + $a->isize
			) {
				$data->{f}->[$_] += $score;
			}
		}
		else {
			foreach ($a->pos - $data->{r_offset} .. 
				$a->pos - $data->{r_offset} + $a->isize
			) {
				$data->{r}->[$_] += $score;
			}
		}
	} 
}

sub pe_center_span {
	my ($a, $data, $score) = @_;
	my $position = $a->pos + int( $a->isize / 2 );
	foreach ($position - $half_extend - $data->{f_offset} .. 
		$position - $data->{f_offset} + $half_extend
	) {
		$data->{f}->[$_] += $score;
	}
}

sub pe_strand_center_span {
	my ($a, $data, $score) = @_;
	my $position = $a->pos + int( $a->isize / 2 );
	my $flag = $a->flag;
	if ($flag & 0x0040) {
		# first read
		if ($a->reversed) {
			foreach ($position - $data->{r_offset} - $half_extend .. 
				$position - $data->{r_offset} + $half_extend
			) {
				$data->{r}->[$_] += $score;
			}
		}
		else {
			foreach ($position - $data->{f_offset} - $half_extend .. 
				$position - $data->{f_offset} + $half_extend
			) {
				$data->{f}->[$_] += $score;
			}
		}
	}
	elsif ($flag & 0x0080) {
		# second read
		if ($a->reversed) {
			foreach ($position - $data->{f_offset} - $half_extend .. 
				$position - $data->{f_offset} + $half_extend
			) {
				$data->{f}->[$_] += $score;
			}
		}
		else {
			foreach ($position - $data->{r_offset} - $half_extend .. 
				$position - $data->{r_offset} + $half_extend
			) {
				$data->{r}->[$_] += $score;
			}
		}
	} 
}


__END__

=head1 NAME

bam2wig.pl

A script to convert Bam alignments into a wig representation file.

=head1 SYNOPSIS

bam2wig.pl [--options...] <file.bam>

bam2wig.pl --extend --rpm --separate --out file --bw file1.bam file2.bam
  
 Required options:
  --in <filename.bam>           repeat if multiple bams, or comma-delimited list
 
 Reporting options (pick one):
  --start                       record at 5' position
  --mid                         record at midpoint of alignment
  --span                        record across entire alignment or pair
  --extend                      extend alignment (record predicted fragment)
  --cspan                       record a span centered on midpoint
  --coverage                    raw alignment coverage
 
 Alignment reporting options:
  --splice                      split alignment at N splices
  --strand                      record separate strands
  --flip                        flip the strands for convenience
  
 Paired-end alignments:
  --pe                          treat as paired-end alignments
  --minsize <integer>           minimum allowed insertion size (30)
  --maxsize <integer>           maximum allowed insertion size (600)
  
 Alignment filtering options:
  --qual <integer>              minimum mapping quality (0)          
  --nosecondary                 skip secondary alignments (false)
  --noduplicate                 skip marked duplicate alignments (false)
  --nosupplementary             skip supplementary alignments (false)
  --chrskip <regex>             regular expression to skip chromosomes
  --blacklist <file>            interval file of regions to skip (bed, gff, txt)
  
  Shift options:
  --shift                       shift reads in the 3' direction
  --shiftval <integer>          explicit shift value in bp (default is to calculate) 
  --extval <integer>            explicit extension size in bp (default is to calculate)
  --chrom <integer>             number of chromosomes to sample (4)
  --minr <float>                minimum pearson correlation to calculate shift (0.5)
  --zmin <float>                minimum z-score from average to test peak for shift (3)
  --zmax <float>                maximum z-score from average to test peak for shift (10)
  --model                       write peak shift model file for graphing
  
 Score options:
  --rpm                         scale depth to Reads Per Million mapped
  --scale <float>               explicit scaling factor, repeat for each bam file
  --separate                    scale multiple bam files independently before merging
  --fraction                    assign fractional counts to all multi-mapped alignments                    
  --format <integer>            number of decimal positions (4)
 
 Output options:
  --out <filename>              output file name, default is bam file basename
  --bw                          convert to bigWig format
  --bwapp /path/to/wigToBigWig  path to external converter
  --gz                          gzip compress output
  
 Wig format:
  --bdg                         bedGraph, default for span, cspan, extend
  --fix                         fixedStep, default for coverage
  --var                         varStep, default for start, mid
  
 General options:
  --cpu <integer>               number of parallel processes (2)
  --verbose                     report additional information
  --version                     print version information
  --help                        show full documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 Input

=over 4

=item --in <filename>

Specify the input Bam alignment file. More than one file may be 
specified, either with repeated options, a comma-delimited list, 
or simply appended to the command. Bam files will be automatically 
indexed if necessary.

=back

=head2 Reporting Options

=over 4

=item --start

Specify that the 5' position should be recorded in the wig file.

=item --mid

Specify that the midpoint of the alignment (single-end) or fragment 
(paired-end) will be recorded in the wig file.

=item --span

Specify that the entire span of the alignment (single-end) or 
fragment (paired-end) will be recorded in the wig file. 

=item --extend

Specify that the alignment should be extended in the 3' direction 
and that the entire length of the extension be recorded in the wig 
file. The extension may be defined by the user or empirically 
determined.

=item --cspan

Specify that a defined span centered at the alignment (single-end) 
or fragment (paired-end) midpoint will be recorded in the wig file.
The span is defined by the extension value.

=item --coverage

Specify that the raw alignment coverage be calculated and reported 
in the wig file. This utilizes a special low-level operation and 
precludes any alignment filtering or post-normalization methods.

=item --position [start|mid|span|extend|cspan|coverage]

Legacy option for supporting previous versions of bam2wig. 

=back

=head2 Alignment reporting options

=over 4

=item --splice

Indicate that the bam file contains alignments with splices, such as 
from RNASeq experiments. Alignments will be split on cigar N operations 
and each sub fragment will be recorded. This only works with single-end 
alignments, and is disabled for paired-end reads. Only start and span 
recording options are supported.

=item --strand

Indicate that separate wig files should be written for each strand. 
The output file basename is appended with either '_f' or '_r' for 
both files. Strand for paired-end alignments are determined by the 
strand of the first read.

=item --flip

Flip the strand of the output files when generating stranded wig files. 
Do this when RNA-Seq alignments map to the opposite strand of the 
coding sequence, depending on the library preparation method. 

=back

=head2 Paired-end alignments

=over 4

=item --pe

The Bam file consists of paired-end alignments, and only properly 
mapped pairs of alignments will be counted. Properly mapped pairs 
include FR reads on the same chromosome, and not FF, RR, RF, or 
pairs aligning to separate chromosomes. The default is to 
treat all alignments as single-end.

=item --minsize <integer>

Specify the minimum paired-end fragment size in bp to accept for recording. 
Default is 30 bp.

=item --maxsize <integer>

Specify the maximum paired-end fragment size in bp to accept for recording. 
Default is 600 bp.

=back

=head2 Alignment filtering options:

=over 4

=item --qual <integer>

Set a minimum mapping quality score of alignments to count. The mapping 
quality is a range from 0-255, with higher numbers indicating lower 
probability of a mapping error. Multi-mapping alignments often have a 
map quality of 0. The default is 0 (accept everything).

=item --secondary | --nosecondary

Boolean flag to accept or skip secondary alignments, indicated by the 
alignment bit flag 0x100. Secondary alignments typically represent 
alternative mapping locations, or multi-mapping events. By default, all 
alignments are accepted. 

=item --duplicate | --noduplicate

Boolean flag to accept or skip duplicate alignments, indicated by the 
alignment bit flag 0x400. Duplicates alignments may represent a PCR or 
optical duplication. By default, all alignments are accepted. 

=item --supplementary | --nosupplementary

Boolean flag to accept or skip supplementary alignments, indicated by 
the alignment bit flag 0x800. Supplementary alignments are typically 
associated with chimeric fragments. By default, all alignments are 
accepted.

=item --chrskip <regex>

Provide a regular expression to skip certain chromosomes. Perl-based 
regular expressions are employed. Expressions should be quoted or 
properly escaped on the command line. Examples might be 
    
    'chrM'
    'scaffold.+'
    'chr.+alt|chrUn.+|chr.+_random'

=item --blacklist <file>

Provide a file of genomic intervals from which to exclude alignments. 
Examples might include repeats, ribosomal RNA, or heterochromatic regions.
The file should be any text file interpretable by L<Bio::ToolBox::Data> 
with chromosome, start, and stop coordinates, including BED and GFF formats.
Note that this only excludes overlapping alignments, and does not include 
extended alignments.

=back

=head2 Shift options

=item --shift

Specify that the positions of the alignment should be shifted towards 
the 3' end. Useful for ChIP-Seq applications, where only the ends of 
the fragments are counted and often seen as separated discrete peaks 
on opposite strands flanking the true target site. This option is 
disabled with paired-end and spliced reads (where it is not needed). 

=item --shiftval <integer>

Provide the value in bp that the recorded position should be shifted. 
The value should be 1/2 the average length of the library insert size.
The default is to automatically and empirically determine the 
appropriate shift value using cross-strand correlation (recommended). 

=item --extval <integer>

Manually set the length for reads to be extended. By default, the shift 
value is determined empirically and extension is set to 2X the shift 
value. This is also used for the cspan mode.

=item --chrom <integer>

Indicate the number of sequences or chromosomes to sample when 
empirically determining the shift value. The reference sequences 
listed in the Bam file header are taken in order of decreasing 
length, and one or more are taken as a representative sample of 
the genome. The default value is 4. 

=item --minr <float>

Provide the minimum Pearson correlation value to accept a shift 
value when empirically determining the shift value. Enter a decimal value 
between 0 and 1. Higher values are more stringent. The default 
is 0.5.

=item --zmin <float>

Specify the minimum z-score (or number of standard deviations) from 
the chromosomal mean depth to test for a peak shift. Increase this 
number to test for strong robust peaks, which give a better estimations 
of the shift value. Default is 3.

=item --zmax <float> 

Specify the maximum z-score (or number of standard deviations) from 
the chromosomal mean depth to test for a peak shift. This excludes 
erroneous peaks due to repetitive sequence alignments with high coverage. 
Increase this number to include more robust peaks that can give a 
better estimation of the shift value. Default is 10.

=item --model

Indicate that the shift model profile data should be written to 
file for examination. The average profile, including for each 
sampled chromosome, are reported for the forward and reverse strands, 
as  well as the shifted profile. A standard text file is generated 
using the output base name. The default is to not write the model 
shift data.

=back

=head2 Score Options

=over 4

=item --rpm

Convert the data to Reads (or Fragments) Per Million mapped. This is useful 
for comparing read coverage between different datasets. The default is 
no RPM conversion. 

=item --scale <float>

Optionally provide your own scaling factor. This will be multiplied with 
every position when generating the wig file. This may be combined with the 
rpm factor as well. When combining multiple bam files, either a single scale 
factor may be supplied for all files, or individual scale factors may be 
supplied for each bam file. If supplying multiple, use the option multiple 
times or give a comma-delimited list. The values should be in the same order 
as the bam files. 

--separate

When combining multiple bam files with the rpm option, scale each bam file 
separately before adding together. This allows for a more accurate averaging 
of replicates without biasing towards samples with more reads. Otherwise, 
the default behavior is to simply add the files together prior to scaling. 

--fraction

Indicate that multi-mapping alignments should be given fractional counts 
instead of full counts. The number of alignments is determined using the 
NH alignment tag. If a read has 10 alignments, then each alignment is 
given a count of 0.1. 

--format <integer>

Indicate the number of decimal postions reported in the wig file. This 
is only applicable when rpm, scale, or fraction options are provided. 
The default value is 4 decimal positions.

=back

=head2 Output Options

=over 4

=item --out <filename>

Specify the output base filename. An appropriate extension will be 
added automatically. By default it uses the base name of the 
input file.

=item --bw

Specify whether or not the wig file should be further converted into 
an indexed, compressed, binary BigWig file. The default is false.

=item --bwapp /path/to/wigToBigWig

Optionally specify the full path to the UCSC I<wigToBigWig> conversion 
utility. The application path may be set in the .biotoolbox.cfg file 
or found in the default executable path, which makes this option 
unnecessary. 

=item --gz

Specify whether (or not) the output file should be compressed with 
gzip. Disable with --nogz.

=back

=head2 Wig format

=over 4

=item --bdg

Specify that the output wig format is a bedGraph-style wig format. This is 
the default format for extend, span, and cspan modes of operation.

=item --fix

Specify that the output wig format is in fixedStep wig format. This is the 
default format for coverage mode of operation.

=item --var

Specify that the output wig format is in variableStep wig format. This is 
the default format for start and midpoint modes of operation.

=back

=head2 General options

=over 4

=item --cpu <integer>

Specify the number of parallel instances to run simultaneously. This requires 
the installation of L<Parallel::ForkManager>. With support enabled, the 
default is 2. Disable multi-threaded execution by setting to 1. 

=item --verbose

Print extra informational statements during processing. The default is false.

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
This extended method of recording infers the mean size of the 
library fragments, thereby emulating the coverage of paired-end 
sequencing using single-end sequence data. The shift value is 
empirically determined from the sequencing data or 
provided by the user. If requested, the shift model profile may be 
written to file. 

The output wig file may be either a variableStep, fixedStep, or 
bedGraph format. The wig file may be further converted into a 
compressed, indexed, binary bigWig format, dependent on the 
availability of the appropriate conversion utilities. 

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
 
 bam2wig.pl --start --shift --model --rpm --in <bamfile>
 
 bam2wig.pl --extend --model --rpm --in <bamfile>

=item Paired-end ChIP-Seq

If both ends of the ChIP eluate fragments are sequenced, then we do not 
need to calculate a shift value. Instead, we will simply count at the 
midpoint of each properly-mapped sequence pair, or record the defined 
fragment span.
 
 bam2wig.pl --mid --pe --rpm --in <bamfile>
 
 bam2wig.pl --span --pe --rpm --in <bamfile>

=item Unstranded RNA-Seq

With RNA-Sequencing, we may be interested in either coverage (generating 
a transcriptome map) or simple tag counts (differential gene expression), 
so we can count in one of two ways. 

To compare RNA-Seq data from different experiments, convert the read 
counts to Reads Per Million Mapped, which will help to normalize read 
counts.
 
 bam2wig --span --splice --rpm --in <bamfile>
 
 bam2wig --mid --rpm --in <bamfile>

=item Stranded, single-end RNA-Seq

If the library was generated in such a way as to preserve strand, then 
we can separate the counts based on the strand of the alignment. Note 
that the reported strand may be accurate or flipped, depending upon 
whether first-strand or second-strand synthesized cDNA was sequenced, 
and whether your aligner took this into account. Check the Bam 
alignments in a genome browser to confirm the orientation relative to 
coding sequences. If alignments are opposite to the direction of 
transcription, you can include the --flip option to switch the output.
 
 bam2wig --span ---splice -strand --rpm --in <bamfile>

 bam2wig --pos mid --strand --rpm --in <bamfile>
 
=item Paired-end RNA-Seq

Since paired-end mode of bam2wig interprets pairs as a single fragment, 
and splices are disable with paired-end alignments, paired-end RNAseq 
may be best treated as single-end RNAseq.
 
=back

=head1 TEXT REPRESENTATION OF RECORDING ALIGNMENTS

To help users visualize how this program records alignments in a wig 
file, drawn below are 10 alignments, five forward and five reverse. 
They may be interpreted as either single-end or paired-end. Drawn 
below are the numbers that would be recorded in a wig file for various 
parameter settings. Note that alignments are not drawn to scale and 
are drawn for visualization purposes only. Values of X represent 10.

=over 4

=item Alignments

  ....>>>>>>.....................................<<<<<<.............
  .....>>>>>>..................................<<<<<<...............
  ........>>>>>>.......................................<<<<<<.......
  ........>>>>>>.........................................<<<<<<.....
  ..........>>>>>>............................................<<<<<<

=item Starts

  ....11..2.1.......................................1.1.....1.1....1

=item Midpoints

  ......11..2.1..................................1.1.....1.1....1...

=item Stranded Starts

  F...11..2.1.......................................................
  R.................................................1.1.....1.1....1

=item Span (Coverage)

  ....122244433311.............................112222111122221211111

=item Mid Span (extend value 2)

  ......121.2211.................................1111....1111...11..

=item Stranded Span

  F...122244433311..................................................
  R............................................112222111122221211111

=item Shifted Starts (shift value 26)

  ........................1.1...11121.1..1..........................

=item Shifted Span (shift value 26)

  ...................11222211112344365544411........................

=item Extend (extend value 52)

  12223445789999XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX999887665321111

=item Paired-End Midpoints

  ............................1...111..1............................

=item Paired-End Mid span (extend value 6)

  ..........................111123333432111.........................

=item Paired-End Span

  ....12224455555555555555555555555555555555555555555443333332211111

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
