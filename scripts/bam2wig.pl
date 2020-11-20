#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use File::Spec;
use File::Temp;
use List::Util qw(sum0);
use List::MoreUtils qw(natatime);
use Bio::ToolBox::db_helper qw(
	open_db_connection
	low_level_bam_coverage
	low_level_bam_fetch
	$BAM_ADAPTER
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

my $VERSION = '1.68';
	
	

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
	$use_smartpe,
	$use_ends,
	$position,
	$use_coverage,
	$splice,
	$paired,
	$fastpaired,
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
	$nosecondary,
	$noduplicate,
	$nosupplementary,
	$max_isize,
	$min_isize,
	$first,
	$second,
	$multi_hit_scale,
	$splice_scale,
	$rpm,
	$do_mean,
	$chrnorm,
	$chrapply,
	$chr_exclude,
	$black_list,
	$bin_size,
	$dec_precison,
	$bigwig,
	$do_fixstep,
	$do_varstep,
	$do_bedgraph,
	$no_zero,
	$bwapp,
	$gz,
	$cpu,
	$max_intron,
	$window,
	$verbose,
	$tempdir,
	$help,
	$print_version,
);
my @bamfiles;
my @scale_values;

# Command line options
GetOptions( 
	'i|in=s'       => \@bamfiles, # one or more bam files
	'o|out=s'      => \$outfile, # name of output file 
	's|start!'     => \$use_start, # record start point
	'd|mid!'         => \$use_mid, # record mid point
	'a|span!'      => \$use_span, # record span
	'cspan!'       => \$use_cspan, # record center span
	'e|extend!'    => \$use_extend, # extend read
	'smartcov!'    => \$use_smartpe, # smart paired coverage
	'ends!'        => \$use_ends, # record paired-end endpoints 
	'coverage!'    => \$use_coverage, # calculate coverage
	'position=s'   => \$position, # legacy option
	'l|splice|split!'   => \$splice, # split splices
	'p|pe!'        => \$paired, # paired-end alignments
	'P|fastpe!'    => \$fastpaired, # fast paired-end alignments
	'I|shift!'     => \$shift, # shift coordinates 3'
	'H|shiftval=i' => \$shift_value, # value to shift coordinates
	'x|extval=i'   => \$extend_value, # value to extend reads
	'chrom=i'      => \$chr_number, # number of chromosomes to sample
	'minr=f'       => \$correlation_min, # R minimum value for shift
	'zmin=f'       => \$zmin, # minimum z-score interval for calculating shift
	'zmax=f'       => \$zmax, # maximum z-score interval for calculating shift
	'M|model!'       => \$model, # write the strand shift model data
	't|strand!'    => \$do_strand, # separate strands
	'flip!'        => \$flip, # flip the strands
	'q|qual=i'     => \$min_mapq, # minimum mapping quality
	'S|nosecondary'  => \$nosecondary, # skip secondary alignments
	'D|noduplicate' => \$noduplicate, # skip duplicate alignments
	'U|nosupplementary' => \$nosupplementary, # skip supplementary alignments
	'maxsize=i'    => \$max_isize, # maximum paired insert size to accept
	'minsize=i'    => \$min_isize, # minimum paired insert size to accept
	'first!'       => \$first, # only take first read
	'second!'      => \$second, # only take second read
	'fraction!'    => \$multi_hit_scale, # scale by number of hits
	'splfrac!'     => \$splice_scale, # divide counts by number of spliced segments
	'r|rpm!'       => \$rpm, # calculate reads per million
	'm|separate|mean!' => \$do_mean, # rpm scale separately
	'scale=s'      => \@scale_values, # user specified scale value
	'chrnorm=f'    => \$chrnorm, # chromosome-specific normalization
	'chrapply=s'   => \$chrapply, # chromosome-specific normalization regex
	'K|chrskip=s'  => \$chr_exclude, # regex for skipping chromosomes
	'B|blacklist=s' => \$black_list, # file for skipping regions
	'bin=i'        => \$bin_size, # size for binning the data
	'format=i'     => \$dec_precison, # format to number of decimal positions
	'b|bw!'        => \$bigwig, # generate bigwig file
	'bwapp=s'      => \$bwapp, # utility to generate a bigwig file
	'bdg!'         => \$do_bedgraph, # write a bedgraph output
	'fix!'         => \$do_fixstep, # write a fixedStep output
	'var!'         => \$do_varstep, # write a varStep output
	'nozero!'       => \$no_zero, # do not write zero coverage
	'z|gz!'        => \$gz, # compress text output
	'c|cpu=i'      => \$cpu, # number of cpu cores to use
	'intron=i'     => \$max_intron, # maximum intron size to allow
	'window=i'     => \$window, # window size to control memory usage
	'V|verbose!'   => \$verbose, # print sample correlations
	'temp=s'       => \$tempdir, # directory to write temp files
	'adapter=s'    => \$BAM_ADAPTER, # explicitly set the adapter version
	'h|help'       => \$help, # request help
	'v|version'    => \$print_version, # print the version
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
	$binpack, $buflength, $coverage_dump, $coverage_sub);
check_defaults();
print " Writing temp files to $tempdir\n" if $verbose;
my $items = $paired ? 'fragments' : 'alignments';

# record start time
my $start_time = time;





### Open files
my @sams;
printf " Processing files %s...\n", join(", ", @bamfiles);
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
	if ($chr_exclude and $chr =~ /$chr_exclude/i) {
		print "  skipping sequence $chr\n" if $verbose;
		next;
	}
	push @seq_list, $chr;
	$seq_name2length{$chr} = $sams[0]->target_len($tid);
}
# set the wrapper reference
# this depends on which adapter was opened
my $wrapper_ref;
if ($splice or $use_smartpe) {
	$wrapper_ref = ref($sams[0]) eq 'Bio::DB::Sam' ? 'Bio::DB::Bam::AlignWrapper' : 
		ref($sams[0]) eq 'Bio::DB::HTS' ? 'Bio::DB::HTS::AlignWrapper' : 'none';
	eval { require $wrapper_ref; 1 };
}
printf " Using the %s Bam adapter and align wrapper $wrapper_ref\n", ref($sams[0]) if $verbose;

### Process user provided black lists
my $black_list_hash = process_black_list();


### Calculate shift value
if ($shift or $use_extend) {
	unless (($shift and $shift_value) or ($use_extend and $extend_value)) {
		print " Calculating 3' shift value...\n";
		$shift_value = determine_shift_value();
	}
	
	# quietly exit here after determining shift value if no wig is to be generated
	exit unless ($use_start or $use_extend or $use_span or $use_mid or $use_cspan 
				or $use_coverage);
	
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
# set center extension global value
my $half_extend = int($extend_value / 2);
if ($use_cspan) {
	print " $items will be center extended by $half_extend bp both directions\n";
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
		$cpu ||= 4;
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
		elsif ($position eq 'smartcov') { 
			# for backwards compatibility
			$use_smartpe = 1;
			$paired = 1;
			$splice = 1;
		}
		elsif ($position eq 'ends') {
			$use_ends = 1;
			$paired = 1;
		}
		else {
			die " unrecognized position value '$position'! see help\n";
		}
	}
	my $position_check = $use_start + $use_mid + $use_span + $use_cspan + $use_extend + 
		$use_coverage + $use_smartpe + $use_ends;
	if ( $position_check > 1) {
		die " Modes are mutually exclusive! Please select only one of\n" . 
			" --start, --mid, --span, --cpsan, --extend, --smartcov, --ends, or --coverage\n";
	}
	elsif (not $position_check) {
		# we allow no position if user has selected shift so that we can calculate
		# the shift value without running through entire wig conversion
		unless ($shift) {
			die " Please select one of the following modes:\n" . 
				" --start, --mid, --span, --cspan, --extend, --smartcov, --ends, or --coverage\n";
		}
	}
	
	# check splices
	unless (defined $splice) {
		$splice = 0;
	}
	if ($paired and $splice) {
		if ($use_span) {
			print " switching mode to smart paired coverage with paired and splice enabled\n";
			$use_span = 0;
			$use_smartpe = 1;
		}
		elsif ($use_smartpe) {
			# perfect
		}
		else {
			die " incompatible mode with paired and splices enabled! Try --smartcov\n";
		}
	}
	$max_intron ||= 0;
	
	# check paired-end requirements
	$paired = 1 if $fastpaired; # for purposes here, fastpair is paired
	if (($paired or $use_smartpe) and ($first or $second)) {
		$paired = 0; # not necessary
		$use_smartpe = 0;
	}
	if ($use_smartpe) {
		$paired = 1;
		$splice = 1;
		
		# warnings
		if ($fastpaired) {
			warn " disabling fast-paired mode with smart-paired coverage\n";
			$fastpaired = 0;
		}
		if ($splice_scale) {
			warn " disabling splice scaling with smart-paired coverage\n";
			$splice_scale = 0;
		}
		
		# required modules
		eval {
			# required for conveniently assembling overlapping and spliced coverage
			require Set::IntSpan::Fast;
		};
		die " Module Set::IntSpan::Fast is required for smart-paired coverage\n" if $@;
	}
	if ($use_ends) {
		$paired = 1;
	}
	
	# incompatible options
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
	if ($use_ends and $shift) {
		warn " disabling shift when recording paired endpoints\n";
		undef $shift;
	}
	if ($use_ends and $splice) {
		warn " disabling splices when recording paired endpoints\n";
		undef $splice;
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
	
	# center span should have extend value
	if ($use_cspan and !$extend_value) {
		die " please use the --extval option to define an extend value when using center span\n";
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
		if ($use_smartpe) {
			# Smart paired-end coverage doesn't care at all about the size of 
			# the insertion, but it is nevertheless tested in the pe_callback to 
			# accomodate other functions. So set this to a reasonably really high number.
			$max_isize = 100000; # 100 kb should be sufficiently high
		}
		else {
			# vast majority of paired-end fragments are less than 600 bp
			# really, really big fragments are most likely mapping errors
			$max_isize = 600;
		}
	}
	unless (defined $min_isize) {
		$min_isize = 30;
	}
	
	# chromosome-specific normalization
	if ($chrnorm or $chrapply) {
		# must have both of these parameters
		die "missing --chrnorm value a for specific-chromosome normalization!\n" unless $chrnorm;
		die "missing --chrapply regex for specific-chromosome normalization!\n" unless $chrapply;
	}
	
	# check flag parameters
	if ($verbose) {
		printf "  %s secondary 0x100 reads\n", $nosecondary ? 'Skipping' : 'Including';
		printf "  %s duplicate 0x400 reads\n", $noduplicate ? 'Skipping' : 'Including';
		printf "  %s supplementary 0x800 reads\n", $nosupplementary ? 'Skipping' : 'Including';
	}
	
	# set bin size
	if ($bin_size) {
		die "bin size cannot be negative!\n" if $bin_size < 0;
	}
	else {
		# set default to 10 bp for any span or coverage, or 1 bp for point data
		$bin_size = ($use_start or $use_mid or $use_ends) ? 1 : 10;
	}
	
	# determine binary file packing and length
	if ($multi_hit_scale or $splice_scale or $chrnorm or 
		($use_coverage and $bin_size > 1)
	) {
		# pack as floating point values, this is 32 bit
		$binpack = 'f';
	}
	else {
		# dealing only with integers here
		# yes we do occasionally have depth greater than 65,536
		# originally short (16 bits), now long (32 bits)
		$binpack = 'L';
	}
	$buflength = length(pack($binpack, 1));
	
	# set coverage dump size and subroutine code global values
	# this is for processing the coverage dump array
	if ($use_coverage) {
		print " ignoring duplicate read filters with coverage\n" if $noduplicate;
		print " ignoring supplementary read filters with coverage\n" if $nosupplementary;
		print " ignoring secondary read filters with coverage\n" if $nosecondary;
		print " ignoring map quality filter with coverage\n" if $min_mapq;
		print " ignoring paired-end option with coverage\n" if $paired;
		print " ignoring RPM option with coverage\n" if $rpm;
		print " ignoring custom scale option with coverage\n" if @scale_values;
		print " ignoring multi-hit fractional counting with coverage\n" if $multi_hit_scale;
		print " ignoring splice segment fractional counting with coverage\n" if $splice_scale;
		print " ignoring chromosome-specific normalization with coverage\n" if $chrnorm;
		
		if ($bin_size > 1) {
			$coverage_dump = int(1000 / $bin_size) * $bin_size;
			$coverage_dump = $bin_size if $coverage_dump == 0;
			$coverage_sub = sub {
				my ($coverage, $chrom_data) = @_;
				my @binned_scores;
				my $iterator = natatime $bin_size, @$coverage;
				while (my @values = $iterator->()) {
					push @binned_scores, mean(@values);
				}
				$$chrom_data .= pack("$binpack*", @binned_scores);
			};
		}
		else {
			$coverage_dump = 1000;
			$coverage_sub = sub {
				my ($coverage, $chrom_data) = @_;
				$$chrom_data .= pack("$binpack*", @$coverage);
			};
		}
	}
	
	# set window length for processing through packed binary chromosome strings
	$window ||= 10000;
		# empirical tests show a window of 1000 to 10000 is best, bigger or smaller 
		# result in longer execution times
	
	
	# set decimal formatting
	unless (defined $dec_precison) {
		$dec_precison = 4 if ($rpm or $multi_hit_scale or @scale_values);
		$dec_precison = 1 if ($use_coverage and $bin_size > 1);
	}
	
	
	# check output file
	unless ($outfile) {
		if (scalar @bamfiles == 1) {
			$outfile = $bamfiles[0];
			$outfile =~ s/\.(?:b|cr)am$//;
		}
		else {
			die " Please define an output filename when providing multiple bam files!\n";
		}
	}
	(undef, my $outdir, $outbase) = File::Spec->splitpath($outfile);
	$outbase =~ s/\.(?:wig|bdg|bedgraph|bw|bigwig)(?:\.gz)?$//i; # strip extension if present
	$outfile =~ s/\.(?:wig|bdg|bedgraph|bw|bigwig)(?:\.gz)?$//i; # strip extension if present
	
	
	# set output temporary directory
	if ($tempdir) {
		unless (-x $tempdir and -d _ and -w _ ) {
			die " Specified temp directory '$tempdir' either does not exist or is not writeable!\n";
		}
		$tempdir = File::Temp->newdir("bam2wigTEMP_XXXX", DIR=>$tempdir, CLEANUP=>1);
	}
	else {
		# use the base path of the output wig file identified above
		unless ($outdir) {
			$outdir = File::Spec->curdir;
		}
		$tempdir = File::Temp->newdir("bam2wigTEMP_XXXX", DIR=>$outdir, CLEANUP=>1);
	}
	# reassemble outbase file
	$outbase = File::Spec->catfile($tempdir, $outbase);
	
	
	# determine output format
	unless ($do_bedgraph or $do_varstep or $do_fixstep) {
		# pick an appropriate format for the user
		if ($use_span or $use_extend or $use_cspan or $use_coverage or $use_smartpe) {
			if ($bin_size > 1) {
				$do_fixstep = 1;
			}
			else {
				$do_bedgraph = 1;
			}
		}
		elsif ($use_start or $use_mid or $use_ends) {
			if ($bin_size > 1) {
				$do_fixstep = 1;
			}
			else {
				$do_varstep = 1;
			}
		}
	}
	if ( ($do_bedgraph + $do_varstep + $do_fixstep) > 1) {
		die " Please select only one of --bedgraph, --fixstep, or --varstep\n";
	}
	if ($do_varstep and $bin_size > 1) {
		warn " Writing variableStep wig files with bin size of $bin_size bp is not supported!" . 
			"\n Writing fixedStep wig file instead\n";
		$do_varstep = 0;
		$do_fixstep = 1;
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
	$main_callback = $fastpaired ? \&fast_pe_callback : $paired ? \&pe_callback : 
		\&se_callback;
		
	### Determine the alignment recording callback method
	# coverage
	if ($use_coverage) {
		# does not use a callback subroutine
		undef $callback;
		print " Recording raw alignment coverage\n";
	}
	# start
	elsif (not $paired and not $do_strand and not $shift and $use_start) {
		$callback = \&se_start;
		print " Recording single-end alignment starts\n";
	}
	elsif (not $paired and not $do_strand and $shift and $use_start) {
		$callback = \&se_shift_start;
		print " Recording single-end shifted alignment starts\n";
	}
	elsif (not $paired and $do_strand and not $shift and $use_start) {
		$callback = \&se_strand_start;
		print " Recording single-end stranded alignment starts\n";
	}
	elsif (not $paired and $do_strand and $shift and $use_start) {
		$callback = \&se_shift_strand_start;
		print " Recording single-end shifted, stranded alignment starts\n";
	}
	# midpoint
	elsif (not $paired and not $do_strand and not $shift and $use_mid) {
		$callback = \&se_mid;
		print " Recording single-end alignment midpoints\n";
	}
	elsif (not $paired and not $do_strand and $shift and $use_mid) {
		$callback = \&se_shift_mid;
		print " Recording single-end shifted alignment midpoints\n";
	}
	elsif (not $paired and $do_strand and not $shift and $use_mid) {
		$callback = \&se_strand_mid;
		print " Recording single-end stranded alignment midpoints\n";
	}
	elsif (not $paired and $do_strand and $shift and $use_mid) {
		$callback = \&se_shift_strand_mid;
		print " Recording single-end shifted, stranded alignment midpoints\n";
	}
	# span
	elsif (not $paired and not $do_strand and not $shift and $use_span) {
		$callback = \&se_span;
		print " Recording single-end alignment span\n";
	}
	elsif (not $paired and not $do_strand and $shift and $use_span) {
		$callback = \&se_shift_span;
		print " Recording single-end shifted alignment span\n";
	}
	elsif (not $paired and $do_strand and not $shift and $use_span) {
		$callback = \&se_strand_span;
		print " Recording single-end stranded alignment span\n";
	}
	elsif (not $paired and $do_strand and $shift and $use_span) {
		$callback = \&se_shift_strand_span;
		print " Recording single-end shifted, stranded alignment span\n";
	}
	# center span
	elsif (not $paired and not $do_strand and not $shift and $use_cspan) {
		$callback = \&se_center_span;
		print " Recording single-end alignment center-span\n";
	}
	elsif (not $paired and not $do_strand and $shift and $use_cspan) {
		$callback = \&se_shift_center_span;
		print " Recording single-end shifted alignment center-span\n";
	}
	elsif (not $paired and $do_strand and not $shift and $use_cspan) {
		$callback = \&se_strand_center_span;
		print " Recording single-end stranded alignment center-span\n";
	}
	elsif (not $paired and $do_strand and $shift and $use_cspan) {
		$callback = \&se_shift_strand_center_span;
		print " Recording single-end shifted, stranded alignment center-span\n";
	}
	# extend
	elsif (not $paired and not $do_strand and not $shift and $use_extend) {
		$callback = \&se_extend;
		print " Recording single-end extended alignment span\n";
	}
	elsif (not $paired and not $do_strand and $shift and $use_extend) {
		$callback = \&se_shift_extend;
		print " Recording single-end shifted, extended alignment span\n";
	}
	elsif (not $paired and $do_strand and not $shift and $use_extend) {
		$callback = \&se_strand_extend;
		print " Recording single-end stranded, extended alignment span\n";
	}
	elsif (not $paired and $do_strand and $shift and $use_extend) {
		$callback = \&se_shift_strand_extend;
		print " Recording single-end shifted, stranded, extended alignment span\n";
	}
	# paired-end start
	elsif ($paired and not $do_strand and $use_start) {
		$callback = \&pe_start;
		print " Recording paired-end fragment start\n";
	}
	elsif ($paired and $do_strand and $use_start) {
		$callback = \&pe_strand_start;
		print " Recording paired-end stranded fragment start\n";
	}
	# paired-end midpoint
	elsif ($paired and not $do_strand and $use_mid) {
		$callback = \&pe_mid;
		print " Recording paired-end fragment midpoint\n";
	}
	elsif ($paired and $do_strand and $use_mid) {
		$callback = \&pe_strand_mid;
		print " Recording paired-end stranded fragment midpoint\n";
	}
	# paired-end span
	elsif ($paired and not $do_strand and $use_span) {
		$callback = \&pe_span;
		print " Recording paired-end fragment span\n";
	}
	elsif ($paired and $do_strand and $use_span) {
		$callback = \&pe_strand_span;
		print " Recording paired-end stranded, fragment span\n";
	}
	# paired-end center span
	elsif ($paired and not $do_strand and $use_cspan) {
		$callback = \&pe_center_span;
		print " Recording paired-end fragment center-span\n";
	}
	elsif ($paired and $do_strand and $use_cspan) {
		$callback = \&pe_strand_center_span;
		print " Recording paired-end stranded, fragment center-span\n";
	}
	# paired-end smart coverage
	elsif ($paired and not $do_strand and $use_smartpe) {
		$callback = \&smart_pe;
		print " Recording smart paired-end coverage\n";
	}
	elsif ($paired and $do_strand and $use_smartpe) {
		$callback = \&smart_stranded_pe;
		print " Recording stranded, smart paired-end coverage\n";
	}
	elsif ($paired and not $do_strand and $use_ends) {
		$callback = \&pe_ends;
		print " Recording paired-end fragment endpoints\n";
	}
	elsif ($paired and $do_strand and $use_ends) {
		$callback = \&pe_strand_ends;
		print " Recording stranded, paired-end fragment endpoints\n";
	}
	else {
		die "programmer error!\n" unless $shift; # special exception
	}
	
	# summary of wig file being written
	if ($do_bedgraph) {
		print " Writing bedGraph format in $bin_size bp increments\n";
	}
	elsif ($do_varstep) {
		printf " Writing variableStep format in $bin_size bp bins\n";
	}
	elsif ($do_fixstep) {
		printf " Writing fixedStep format in $bin_size bp bins\n";
	}
}

sub process_black_list {
	if ($black_list and -e $black_list) {
		eval {require 'Bio::ToolBox::Data'};
		my $i = 0;
		eval {require Set::IntervalTree; $i = 1;};
		unless ($i) {
			warn " PROBLEM! Please install Set::IntervalTree to use black lists\n";
			undef $black_list;
			return;
		}
		my %black_list_hash = map { $_ => [] } @seq_list;
		my $Data = Bio::ToolBox::Data->new(file => $black_list) or 
			die "unable to read black list file '$black_list'\n";
		$Data->iterate( sub {
			my $row = shift;
			push @{ $black_list_hash{ $row->seq_id } }, 
				[ $row->start - 1, $row->end ]
				if exists $black_list_hash{ $row->seq_id };
		} );
		printf " Loaded %s blacklist regions\n", 
			format_with_commas($Data->last_row);
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
	my @all_r_values;
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
			push @all_r_values, @{ $result->[4] };
			push @regions, @{ $result->[5] };
			push @r_values, @{ $result->[6] };
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
			my $result = calculate_strand_correlation($sam, $tid);
			
			# record the results for this chromosome
			push @shift_values, @{ $result->[0] }; # push the actual values
			push @f_profile, @{ $result->[1] }; 
			push @r_profile, @{ $result->[2] };
			push @shifted_profile, @{ $result->[3] };
			push @all_r_values, @{ $result->[4] };
			push @regions, @{ $result->[5] };
			push @r_values, @{ $result->[6] };
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
	my @trimmed_bestr;
	foreach my $i (0 .. $#shift_values) {
		if ($shift_values[$i] >= $raw_min and $shift_values[$i] <= $raw_max) {
			push @trimmed_shift_values, $shift_values[$i];
			push @trimmed_f_profile, $f_profile[$i];
			push @trimmed_r_profile, $r_profile[$i];
			push @trimmed_shifted_profile, $shifted_profile[$i];
			push @trimmed_r_values, $all_r_values[$i];
			push @trimmed_regions, $regions[$i];
			push @trimmed_bestr, $r_values[$i];
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
			\@trimmed_shifted_profile, \@trimmed_r_values, \@trimmed_regions, 
			\@trimmed_shift_values, \@trimmed_bestr);
	}
	
	# done
	return $best_value;
}


sub calculate_strand_correlation {
	my ($sam, $tid) = @_;
	
	my ($collected, $data) = scan_high_coverage($sam, $tid);
	
	my @shift_values;
	my @bestr_values;
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
				push @bestr_values, $best_r;
			}
		}
	}
	
	return [ \@shift_values, \@f_profile, \@r_profile, 
		\@shifted_profile, \@all_r_values, \@regions, \@bestr_values ];
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
		my $readsum = sum0( map { $data{f}->[$_] || 0 } ($start .. $start + 49) );
		my $readsum += sum0( map { $data{r}->[$_] || 0 } ($start .. $start + 49) );
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
	my ($value, $f_profile, $r_profile, $shifted_profile, $r_valuess, 
		$regions, $region_shifts, $region_bestr) = @_;
	
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
	
	
	### Write regions
	undef $Data;
	$Data = Bio::ToolBox::Data->new(
		feature  => 'Correlated regions',
		datasets => ['Region', 'Shift', 'BestCorrelation']
	);
	for my $i (0 ..  $#{$regions}) {
		$Data->add_row( [ $regions->[$i], $region_shifts->[$i], $region_bestr->[$i] ] );
	}
	$success = $Data->write_file(
		'filename' => "$outfile\_correlated_regions.txt",
		'gz'       => 0
	);
	print "  Wrote correlated regions data file $success\n" if $success;
}


sub mean {
	return sum0(@_) / scalar(@_);
}

sub open_wig_file {
	my ($name, $do_bw) = @_;
	
	# open a bigWig file handle if requested
	if ($bigwig and $do_bw) {
		print " Writing directly to bigWig converter\n";
		$name .= '.bw' unless $name =~ /\.bw$/;
		$chromo_file = generate_chromosome_file($sams[0], $chr_exclude) 
			unless $chromo_file;
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
	my $do_gz = ($gz and $do_bw) ? 1 : 0; 
		# using the do_bw value because that tells us if it's a temp file or not
	$name .= '.gz' if ($do_gz and $name !~ /\.gz$/i);
	my $fh = open_to_write_fh($name, $do_gz) or 
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
		my @filelist = glob("$outbase.*.temp.bin");
		my $expected = scalar(@sams) * scalar(@seq_list) * ($do_strand ? 2 : 1);
		if (scalar(@filelist) != $expected) {
			die sprintf(" unable to find all the temporary binary children files! Only found %d, expected %d\n",
				scalar(@filelist), $expected);
		}
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
		my @filelist = glob("$outbase.*.temp.bin");
		my $expected = scalar(@sams) * scalar(@seq_list) * ($do_strand ? 2 : 1);
		if (scalar(@filelist) != $expected) {
			die sprintf(" unable to find all the temporary binary children files! Only found %d, expected %d\n",
				scalar(@filelist), $expected);
		}
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
	for (my $start = 0; $start < $seq_length; $start += $coverage_dump) {
		# set endpoint
		my $end = $start + $coverage_dump;
		$end = $seq_length if $end > $seq_length;
		
		# using the low level interface for a little more performance
		my $coverage = low_level_bam_coverage($sam, $tid, $start, $end);
		
		# record the coverage
		&$coverage_sub($coverage, \$chrom_data);
	}
	
	# write out chromosome binary file, set count to arbitrary 0
	if ($do_temp_bin) {
		write_bin_file(\$chrom_data, 
			join('.', $outbase, $samid, $seq_id, 0, 'f', 'temp.bin'));
	}
	else {
		&$wig_writer(\$chrom_data, $binpack, $seq_id, $seq_length, 
			join('.', $outbase, $samid, $seq_id, 0, 'f', 'temp.wig'));
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
		printf " Finished converting $items in %.3f minutes\n", 
			(time - $start_time) / 60;
	}
	
	# find and merge binary files
	# we can do this in parallel too!
	if ($do_temp_bin) {
		
		# find the children
		my @totals;
		my %seq_totals;
		my %files;
		my @filelist = glob("$outbase.*.temp.bin");
		my $expected = scalar(@sams) * scalar(@seq_list) * ($do_strand ? 2 : 1);
		if (scalar(@filelist) != $expected) {
			die sprintf(" unable to find all the temporary binary children files! Only found %d, expected %d\n",
				scalar(@filelist), $expected);
		}
		foreach my $file (@filelist) {
			# each file name is basename.samid.seqid.count.strand.bin.gz
			if ($file =~ /$outbase\.(\d+)\.(.+)\.(\d+)\.([fr])\.temp\.bin\Z/) {
				my $samid = $1;
				my $seq_id = $2;
				my $count = $3;
				my $strand = $4;
				$totals[$samid] += $count;
				$seq_totals{$seq_id}{$strand} += $count;
				$files{$seq_id}{$strand}{$samid} = $file;
			}
		}
		
		# merging multiple files
		if (scalar(@sams) > 1) {
			print " Merging temporary files from each bam\n" if (scalar(@sams) > 1);
			if ($rpm) {
				print " Normalizing depth\n" if $rpm;
				for my $i (0 .. $#totals) {
					printf "  %s had %s total counted $items\n", $bamfiles[$i], 
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
				my $factor = 1_000_000 / ( sum0(@totals) * scalar(@sams) );
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
				print " Scaling depth with user-supplied factor\n";
				if (@norms) {
					for my $i (0 .. $#norms) {
						$norms[$i] *= $scale_values[$i];
					}
				}
				else {
					@norms = @scale_values;
				}
			}
			printf "  Normalization factors: %s\n", join(' ', @norms) if $verbose;
			
			# merge the samples
			foreach my $seq_id (@seq_list) {
				foreach my $strand (qw(f r)) {
					next unless defined $files{$seq_id}{$strand};
					merge_bin_files($seq_id, $strand, $seq_totals{$seq_id}{$strand}, 
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
			my $scale_factor = 1;
			if ($rpm) {
				printf " Normalizing depth based on %s total counted $items\n",
					format_with_commas($totals[0]); 
				$scale_factor = 1_000_000 / $totals[0];
			}
			if (scalar @scale_values) {
				# user supplied scaling factor
				print " Scaling depth with user-supplied factor\n";
				$scale_factor *= $scale_values[0];
			}
			
			# normalize the wig files
			foreach my $seq_id (@seq_list) {
				foreach my $strand (qw(f r)) {
					next unless defined $files{$seq_id}{$strand};
					normalize_bin_file($files{$seq_id}{$strand}{0}, $scale_factor, $seq_id);
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
		printf " Finished converting $items in %.3f minutes\n", 
			(time - $start_time) / 60;
	}
	
	# find and merge binary files
	# we can do this in parallel too!
	if ($do_temp_bin) {
		
		# find the children
		my @totals;
		my %seq_totals;
		my %files;
		my @filelist = glob("$outbase.*.temp.bin");
		my $expected = scalar(@sams) * scalar(@seq_list) * ($do_strand ? 2 : 1);
		if (scalar(@filelist) != $expected) {
			die sprintf(" unable to find all the temporary binary children files! Only found %d, expected %d\n",
				scalar(@filelist), $expected);
		}
		foreach my $file (@filelist) {
			# each file name is basename.samid.seqid.count.strand.bin.gz
			if ($file =~ /$outbase\.(\d+)\.(.+)\.(\d+)\.([fr])\.temp\.bin\Z/) {
				my $samid = $1;
				my $seq_id = $2;
				my $count = $3;
				my $strand = $4;
				$totals[$samid] += $count;
				$seq_totals{$seq_id}{$strand} += $count;
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
					printf "  %s had %s total counted $items\n", $bamfiles[$i], 
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
				my $factor = 1_000_000 / ( sum0(@totals) * scalar(@sams) );
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
				print " Scaling depth with user-supplied factor\n";
				if (@norms) {
					for my $i (0 .. $#norms) {
						$norms[$i] *= $scale_values[$i];
					}
				}
				else {
					@norms = @scale_values;
				}
			}
			printf "  Normalization factors: %s\n", join(' ', @norms) if $verbose;
			
			# merge the samples
			foreach my $seq_id (@seq_list) {
				foreach my $strand (qw(f r)) {
					next unless defined $files{$seq_id}{$strand};
					$pm->start and next;
					merge_bin_files($seq_id, $strand, $seq_totals{$seq_id}{$strand}, 
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
			my $scale_factor = 1;
			if ($rpm) {
				printf " Normalizing depth based on %s total counted $items\n",
					format_with_commas($totals[0]); 
				$scale_factor = 1_000_000 / $totals[0];
			}
			if (scalar @scale_values) {
				# user supplied scaling factor
				print " Scaling depth with user-supplied factor\n";
				$scale_factor *= $scale_values[0];
			}
			
			# normalize the wig files
			foreach my $seq_id (@seq_list) {
				foreach my $strand (qw(f r)) {
					next unless defined $files{$seq_id}{$strand};
					$pm->start and next;
					normalize_bin_file($files{$seq_id}{$strand}{0}, $scale_factor, $seq_id);
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
		f_count         => 0,
		r_count         => 0,
		sam             => $sam,
		score           => 1,
	};
	
	# chromosome specific normalization
	if ($chrapply and $seq_id =~ /$chrapply/i) {
		$data->{score} = $chrnorm;
	}
	
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
	my $fdiff = int($seq_length/$bin_size) - $data->{f_offset};
	$data->{fpack} .= pack("$binpack$fdiff", @{$data->{f}});
	if ($do_strand) {
		my $rdiff = int($seq_length/$bin_size) - $data->{r_offset};
		$data->{rpack} .= pack("$binpack$rdiff", @{$data->{r}});
	}
	
	# round the final count to a solid integer as necessary
	if ($multi_hit_scale) {
		$data->{f_count} = sprintf("%.0f", $data->{f_count});
		$data->{r_count} = sprintf("%.0f", $data->{r_count}) if $do_strand;
	}
	
	# write out file
	# we always write the forward strand, and reverse strand if stranded data
	if ($do_temp_bin) {
		# write a temporary binary file for merging later
		write_bin_file(\$data->{fpack}, 
			join('.', $outbase, $samid, $seq_id, $data->{f_count}, 'f', 'temp.bin') );
		write_bin_file(\$data->{rpack}, 
			join('.', $outbase, $samid, $seq_id, $data->{r_count}, 'r', 'temp.bin')
		) if $do_strand;
	}
	else {
		# write a chromosome specific wig file
		&$wig_writer(\$data->{fpack}, $binpack, $seq_id, $seq_length, 
			join('.', $outbase, $samid, $seq_id, $data->{f_count}, 'f', 'temp.wig') );
		&$wig_writer(\$data->{rpack}, $binpack, $seq_id, $seq_length, 
			join('.', $outbase, $samid, $seq_id, $data->{r_count}, 'r', 'temp.wig') 
		) if $do_strand;
	}
	
	# verbose status line 
	if ($verbose) {
		printf "  Converted %s $items on $seq_id in %d seconds\n", 
			format_with_commas( $data->{f_count} + $data->{r_count} ), time - $chr_start_time;
		if ($chrapply and $seq_id =~ /$chrapply/i) {
			printf "   Scaled counts by $chrnorm\n";
		}
		if ($paired and keys %{$data->{pair}}) {
			printf "   %d orphan paired alignments were left behind and not counted!\n", 
				scalar keys %{$data->{pair}};
		}
	}
}

sub write_bin_file {
	my ($data, $filename) = @_;
		# note that $data is a reference
	my $fh = open_to_write_fh($filename) or 
		die " unable to write temporary file '$filename'!\n";
	$fh->binmode;
	$fh->print($$data);
	$fh->close;
}

sub merge_bin_files {
	my ($seq_id, $strand, $total, $files, $norm_factors) = @_;
	my $merge_start_time = time;
	my $long_window = 100 * $window;
	my $seq_bin_length = int($seq_name2length{$seq_id} / $bin_size);
	
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
	for (my $pos = 0; $pos < $seq_bin_length; $pos += $long_window) {
		# check length
		my $len = ($pos + $window) > $seq_bin_length ? 
					($seq_bin_length - $pos) : $long_window;
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
	&$wig_writer(\$chrom_data, 'f', $seq_id, $seq_name2length{$seq_id}, 
		join('.', $outbase, '0', $seq_id, $total, $strand, 'temp.wig') );
	if ($verbose) {
		printf "  Merged%s $seq_id temp files in %d seconds\n", 
			defined $norm_factors->[0] ? " and normalized" : "", time - $merge_start_time;
	}
}

sub normalize_bin_file {
	my ($file, $scale_factor, $seq_id) = @_;
	my $norm_start_time = time;
	my $seq_bin_length = int($seq_name2length{$seq_id} / $bin_size);
	my $long_window = 10 * $window;
	
	# open file
	my $fh = open_to_read_fh($file) or 
		die " unable to read temporary $file!\n";
	$fh->binmode;
	
	# march along chromosome in defined windows to keep memory usage down
	# apply normalization as data is loaded into the combined chrom_data array
	my $chrom_data; 
	for (my $pos = 0; $pos < $seq_bin_length; $pos += $long_window) {
		# check length
		my $len = ($pos + $long_window) > $seq_bin_length ? 
					($seq_bin_length - $pos) : $long_window;
		
		# read, unpack, normalize, and re-pack current window from binary file
		my $string;
		$fh->read($string, $len * $buflength);
		$chrom_data .= pack("f*", map {$_ * $scale_factor} unpack("$binpack*", $string));
	}
	
	# finish
	$fh->close;
	unlink $file;
	
	# now write the wig file
	$file =~ s/\.bin$/.wig/;
	&$wig_writer(\$chrom_data, 'f', $seq_id, $seq_name2length{$seq_id}, $file);
	if ($verbose) {
		printf "  Scaled $seq_id temp file in %d seconds\n", time - $norm_start_time;
	}
}

sub write_final_wig_file {
	
	# find children files
	my @filelist = glob("$outbase.*.temp.wig");
	my $expected = scalar(@seq_list) * ($do_strand ? 2 : 1);
	if (scalar(@filelist) != $expected) {
		die sprintf(" unable to find all the children wig files! Only found %d, expected %d\n",
			scalar(@filelist), $expected);
	}
	
	# assemble into a hash
	my %files;
	my $f_total = 0;
	my $r_total = 0;
	foreach my $file (@filelist) {
		# each file name is basename.samid.seqid.count.strand.temp.wig
		if ($file =~ /$outbase\.\d+\.(.+)\.(\d+)\.([fr])\.temp\.wig\Z/) {
			my $seq_id = $1;
			my $total  = $2;
			my $strand = $3;
			$files{$seq_id}{$strand} = $file;
			$f_total += $total if $strand eq 'f';
			$r_total += $total if $strand eq 'r';
		}
	}
	my @f_filelist = map { $files{$_}{f} } @seq_list;
	my @r_filelist = map { $files{$_}{r} } @seq_list;
	
	# print total alignment summaries
	if ($r_total and $flip) {
		printf " %s total forward $items\n", format_with_commas($r_total);
		printf " %s total reverse $items\n", format_with_commas($f_total);
	}
	elsif ($r_total) {
		printf " %s total forward $items\n", format_with_commas($f_total);
		printf " %s total reverse $items\n", format_with_commas($r_total);
	}
	else {
		printf " %s total $items\n", format_with_commas($f_total);
	}
	
	# write wig files with the appropriate wig writer
	if ($do_strand and !$flip) {
		
		if ($cpu > 1) {
			# we can fork this!!!!
			my $pm = Parallel::ForkManager->new(2);
			for my $i (1 .. 2) {
				$pm->start and next;
				merge_wig_files("$outfile\_f", @f_filelist) if $i == 1;
				merge_wig_files("$outfile\_r", @r_filelist) if $i == 2;
				$pm->finish;
			}
			$pm->wait_all_children;
		}
		else {
			merge_wig_files("$outfile\_f", @f_filelist);
			merge_wig_files("$outfile\_r", @r_filelist);
		}
	}
	elsif ($do_strand and $flip) {
		if ($cpu > 1) {
			# we can fork this!!!!
			my $pm = Parallel::ForkManager->new(2);
			for my $i (1..2) {
				$pm->start and next;
				merge_wig_files("$outfile\_r", @f_filelist) if $i == 1;
				merge_wig_files("$outfile\_f", @r_filelist) if $i == 2;
				$pm->finish;
			}
			$pm->wait_all_children;
		}
		else {
			merge_wig_files("$outfile\_r", @f_filelist);
			merge_wig_files("$outfile\_f", @r_filelist);
		}
	}
	else {
		merge_wig_files($outfile, @f_filelist);
	}
	
	# clean up 
	unlink $chromo_file if $chromo_file;
}

sub write_bedgraph {
	my ($data, $packer, $seq_id, $seq_length, $filename) = @_;	
		# note that $data is a reference
	
	# set the printf formatter for decimal or integer
	my $formatter = length($dec_precison) ? 
		"$seq_id\t%d\t%d\t%." . $dec_precison. "f\n" : "$seq_id\t%d\t%d\t%s\n";
	
	# work though chromosome
	my $seq_bin_length = int($seq_length / $bin_size);
	my $buflength = length(pack($packer, 1));
	my $out_string;
	my $cpos = 0; # current position
	my $lpos = 0; # last position
	my $cval = 0; # current value
	for (my $pos = 0; $pos < $seq_bin_length; $pos += $window) {
		# check length
		my $len = ($pos + $window) > $seq_bin_length ? ($seq_bin_length - $pos) : $window;
		
		# unpack current window from the passed binary string
		my @win_data = unpack("$packer*", 
			substr($$data, $pos * $buflength, $len * $buflength));
		
		# work through current window
		foreach my $value (@win_data) {
			if ($value == $cval) {
				$cpos++;
			}
			else {
				my $end = $cpos * $bin_size;
				$end = $seq_length if $end > $seq_length;
				if ($no_zero and $cval == 0) {
					# we will skip the nonzero intervals 
					# do nothing here
				}
				else {
					$out_string .= sprintf($formatter, $lpos * $bin_size, $end, $cval) 
						if $end; # this avoids writing nonexistent start 0 end 0 lines
				}
				$lpos = $cpos;
				$cval = $value;
				$cpos++;
			}
		}
	}
	
	# final write
	if ($cpos > $lpos) {
		my $end = $cpos * $bin_size;
		$end = $seq_length if $end != $seq_length;
		$out_string .= sprintf($formatter, $lpos * $bin_size, $end, $cval);
	}
	
	
	# write wig file
	my ($filename, $outfh) = open_wig_file($filename, 0);
	$outfh->print($out_string);
	$outfh->close;
}

sub write_fixstep {
	my ($data, $packer, $seq_id, $seq_length, $filename) = @_;	
		# note that $data is a reference

	# set the printf formatter for decimal or integer
	my $formatter = length($dec_precison) ? "%." . $dec_precison . "f\n" : "%s\n";
	
	# write fixStep header
	my $out_string = "fixedStep chrom=$seq_id start=1 step=$bin_size span=$bin_size\n";
	
	# work though chromosome
	my $seq_bin_length = int($seq_length / $bin_size);
	my $buflength = length(pack($packer, 1));
	for (my $pos = 0; $pos < $seq_bin_length; $pos += $window) {
		# check length
		my $len = ($pos + $window) > $seq_bin_length ? ($seq_bin_length - $pos) : $window;
		
		# unpack current window from the passed binary string
		my @win_data = unpack("$packer*", 
			substr($$data, $pos * $buflength, $len * $buflength));
		
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
		# note that $data is a reference

	# set the printf formatter for decimal or integer
	my $formatter = length($dec_precison) ? "%d\t%." . $dec_precison. "f\n" : "%d\t%s\n";
	
	# write fixStep header
	# we are only supporting bin_size of 1 bp for varStep
	my $out_string = "variableStep chrom=$seq_id\n";
	
	# work though chromosome
	my $buflength = length(pack($packer, 1));
	for (my $pos = 0; $pos < $seq_length; $pos += $window) {
		# check length
		my $len = ($pos + $window) > $seq_length ? ($seq_length - $pos) : $window;
		
		# unpack current window from the passed binary string
		my @win_data = unpack("$packer*", 
			substr($$data, $pos * $buflength, $len * $buflength));
		
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
	my ($ofile, @files) = @_;
	
	my ($filename1, $fh) = open_wig_file($ofile, 1);
	while (@files) {
		my $file = shift @files;
		my $in = open_to_read_fh($file); 
		unless ($in) {
			warn "one or more sub-process forks failed! Check your parameters and input file.\nAttempting to clean up\n";
			foreach (@files) {unlink $_;}
			unlink $chromo_file if $chromo_file;
			exit 1;
		}
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
	return if ($nosecondary and $flag & 0x100); # secondary alignment
	return if ($noduplicate and $flag & 0x400); # marked duplicate
	return if ($flag & 0x200); # QC failed but still aligned? is this necessary?
	return if ($nosupplementary and $flag & 0x800); # supplementary hit
	return if ($first and not $flag & 0x40); # first read in pair
	return if ($second and not $flag & 0x80); # second read in pair
	
	# filter black listed regions
	if (defined $data->{black_list}) {
		my $results = $data->{black_list}->fetch($a->pos, $a->calend);
		return if @$results;
	}
	
	# scale by number of hits
	my $score; 
	if ($multi_hit_scale) {
		# preferentially use the number of included hits, then number of hits
		my $nh = $a->aux_get('IH') || $a->aux_get('NH') || 1;
		$score = $nh > 1 ? 1/$nh : 1;
		# record fractional alignment counts
		if ($do_strand and $a->reversed) {
			$data->{r_count} += $score; 
		}
		else {
			$data->{f_count} += $score;
		}
		$score *= $data->{score}; # multiply by chromosome scaling factor
	}
	else {
		 # always record one alignment
		if ($do_strand and $a->reversed) {
			$data->{r_count} += 1; 
		}
		else {
			$data->{f_count} += 1;
		}
		$score = $data->{score}; # probably 1, but may be chromosome scaled
	}
	
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

sub fast_pe_callback {
	my ($a, $data) = @_;
	
	# check paired status
	return unless $a->proper_pair; # both alignments are mapped
	return if $a->reversed; # only look at forward alignments, that's why it's fast
	return unless $a->tid == $a->mtid; # same chromosome?
	my $isize = abs($a->isize); 
	return if $isize > $max_isize;
	return if $isize < $min_isize; 
	
	# check alignment quality and flags
	return if ($min_mapq and $a->qual < $min_mapq); # mapping quality
	my $flag = $a->flag;
	return if ($nosecondary and $flag & 0x0100); # secondary alignment
	return if ($noduplicate and $flag & 0x0400); # marked duplicate
	return if ($flag & 0x0200); # QC failed but still aligned? is this necessary?
	return if ($nosupplementary and $flag & 0x0800); # supplementary hit
	
	# filter black listed regions
	if ($data->{black_list}) {
		my $results = $data->{black_list}->fetch($a->pos, $a->pos + $isize);
		return if @$results;
	}
	
	# scale by number of hits
	my $score;
	if ($multi_hit_scale) {
		my $nh = $a->aux_get('IH') || $a->aux_get('NH') || 1;
		$score = $nh > 1 ? 1/$nh : 1;
		# record fractional alignment counts
		if ($do_strand and $flag & 0x80) {
			# this alignment is forward and second mate, so must be reverse strand
			$data->{r_count} += $score; 
		}
		else {
			# otherwise this forward alignment must be first mate
			# or unstranded analysis defaults to forward strand
			$data->{f_count} += $score;
		}
		$score *= $data->{score}; # multiply by chromosome scaling factor
	}
	else {
		 # always record one alignment
		if ($do_strand and $flag & 0x80) {
			# this alignment is forward and second mate, so must be reverse strand
			$data->{r_count} += 1; 
		}
		else {
			# otherwise this forward alignment is second mate
			# or unstranded analysis defaults to forward strand
			$data->{f_count} += 1;
		}
		$score = $data->{score}; # probably 1, but may be chromosome scaled
	}
	
	# record based on the forward read
	&$callback($a, $data, $score);
	
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

sub pe_callback {
	my ($a, $data) = @_;
	
	# check paired status
	return unless $a->proper_pair; # both alignments are mapped
	return unless $a->tid == $a->mtid; # same chromosome?
	my $isize = $a->isize; 
	if ($a->reversed) {
		# in proper FR pairs reverse alignments are negative
		return if $isize > 0; # pair is RF orientation
		$isize = abs($isize);
	}
	return if $isize > $max_isize;
	return if $isize < $min_isize; 
	
	# check alignment quality and flags
	return if ($min_mapq and $a->qual < $min_mapq); # mapping quality
	my $flag = $a->flag;
	return if ($nosecondary and $flag & 0x0100); # secondary alignment
	return if ($noduplicate and $flag & 0x0400); # marked duplicate
	return if ($flag & 0x0200); # QC failed but still aligned? is this necessary?
	return if ($nosupplementary and $flag & 0x0800); # supplementary hit
	
	# filter black listed regions
	if ($data->{black_list}) {
		my $results;
		if ($a->reversed) {
			$results = $data->{black_list}->fetch($a->calend - $isize, $a->calend);
		}
		else {
			$results = $data->{black_list}->fetch($a->pos, $a->pos + $isize);
		}
		return if @$results;
	}
	
	# look for pair
	if ($a->reversed) {
		# since we look for FR pairs, the F read should already be found
		my $f = $data->{pair}->{$a->qname} or return;
		delete $data->{pair}->{$a->qname};
		
		# scale by number of hits
		my $score;
		if ($multi_hit_scale) {
			my $r_nh = $a->aux_get('IH') || $a->aux_get('NH') || 1;
			my $f_nh = $f->aux_get('IH') || $f->aux_get('NH') || 1;
			if ($f_nh == $r_nh) {
				$score = 1/$f_nh;
			}
			elsif ($f_nh < $r_nh) {
				# take the lowest number of hits recorded
				$score = $f_nh > 1 ? 1/$f_nh : 1;
			}
			else {
				$score = $r_nh > 1 ? 1/$r_nh : 1;
			}
			# record fractional alignment counts
			if ($do_strand and $flag & 0x40) {
				# this alignment is reverse and first mate, so must be reverse strand
				$data->{r_count} += $score; 
			}
			else {
				# otherwise this reverse alignment is second mate
				# or unstranded analysis defaults to forward strand
				$data->{f_count} += $score;
			}
			$score *= $data->{score}; # multiply by chromosome scaling factor
		}
		else {
			# always record one alignment
			if ($do_strand and $flag & 0x40) {
				# this alignment is reverse and first mate, so must be reverse strand
				$data->{r_count} += 1; 
			}
			else {
				# otherwise this reverse alignment is second mate
				# or unstranded analysis defaults to forward strand
				$data->{f_count} += 1;
			}
			$score = $data->{score}; # probably 1, but may be chromosome scaled
		}
		
		# record based primarily on the forward read, but pass reverse read at end
		# for use in smart pairing
		&$callback($f, $data, $score, $a);
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
	
	# check intron size from the cigar string
	if ($max_intron) {
		my $size = 1;
		my $cigars = $aw->cigar_array;
		foreach my $c (@$cigars) {
			# each element is [operation, size]
			$size = $c->[1] if ($c->[0] eq 'N' and $c->[1] > $size);
		}
		return if $size > $max_intron; # exceed maximum intron size
	}
	
	# get exon segments
	my @segments = $aw->get_SeqFeatures;
	
	# adjust score
	if ($splice_scale) {
		$score /= scalar(@segments);
	}
	
	# record
	if ($use_start) {
		foreach my $segment (@segments) {
			if ($do_strand and $a->reversed) {
				# reverse strand
				my $pos = int( $segment->end / $bin_size );
				$data->{r}->[$pos - $data->{r_offset}] += $score;
			}
			else {
				# otherwise forward strand
				my $pos = int( ($segment->start - 1) / $bin_size );
				$data->{f}->[$pos - $data->{f_offset}] += $score;
			}
		}
	}
	elsif ($use_span) {
		foreach my $segment (@segments) {
			my $start = int( ($segment->start - 1) / $bin_size);
			my $end = int( ($segment->end - 1) / $bin_size);
			if ($do_strand and $a->reversed) {
				# reverse strand
				foreach ($start - $data->{r_offset} .. $end - $data->{r_offset}) {
					$data->{r}->[$_] += $score;
				}
			}
			else {
				# otherwise forward strand
				foreach ($start - $data->{f_offset} .. $end - $data->{f_offset}) {
					$data->{f}->[$_] += $score;
				}
			}
		}
	}
}

sub se_start {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		$data->{f}->[int( ($a->calend - 1) / $bin_size) - $data->{f_offset}] += $score;
	}
	else {
		$data->{f}->[int($a->pos / $bin_size) - $data->{f_offset}] += $score;
	}
}

sub se_shift_start {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $pos = int( ($a->calend - 1 - $shift_value) / $bin_size);
		$data->{f}->[$pos - $data->{f_offset}] += $score if $pos >= 0;
	}
	else {
		my $pos = int( ($a->pos + $shift_value) / $bin_size);
		$data->{f}->[$pos - $data->{f_offset}] += $score if $pos >= 0;
	}
}

sub se_strand_start {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $pos = int( ($a->calend - 1) / $bin_size);
		$data->{r}->[$pos - $data->{r_offset}] += $score;
	}
	else {
		my $pos = int( $a->pos / $bin_size );
		$data->{f}->[$a->pos - $data->{f_offset}] += $score;
	}
}

sub se_shift_strand_start {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $pos = int( ($a->calend - 1 - $shift_value) / $bin_size);
		$data->{r}->[$pos - $data->{r_offset}] += $score if $pos >= 0;
	}
	else {
		my $pos = int( ($a->pos + $shift_value) / $bin_size);
		$data->{f}->[$pos - $data->{f_offset}] += $score if $pos >= 0;
	}
}

sub se_mid {
	my ($a, $data, $score) = @_;
	my $pos = int( int( ($a->pos + $a->calend - 1) / 2) / $bin_size);
	$data->{f}->[$pos - $data->{f_offset}] += $score;
}

sub se_shift_mid {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	if ($a->reversed) {
		my $pos = int( ($mid - $shift_value) / $bin_size);
		$data->{f}->[$pos - $data->{f_offset}] += $score if $pos >= 0;
	}
	else {
		my $pos = int( ($mid + $shift_value) / $bin_size);
		$data->{f}->[$pos - $data->{f_offset}] += $score if $pos >= 0;
	}
}

sub se_strand_mid {
	my ($a, $data, $score) = @_;
	my $pos = int( int( ($a->pos + $a->calend -1) / 2) / $bin_size);
	if ($a->reversed) {
		$data->{r}->[$pos - $data->{r_offset}] += $score;
	}
	else {
		$data->{f}->[$pos - $data->{f_offset}] += $score;
	}
}

sub se_shift_strand_mid {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	if ($a->reversed) {
		my $pos = int( ($mid - $shift_value) / $bin_size);
		$data->{r}->[$pos - $data->{r_offset}] += $score if $pos >= 0;
	}
	else {
		my $pos = int( ($mid + $shift_value) / $bin_size);
		$data->{f}->[$pos - $data->{f_offset}] += $score if $pos >= 0;
	}
}

sub se_span {
	my ($a, $data, $score) = @_;
	my $s = int($a->pos / $bin_size) - $data->{f_offset};
	my $e = int( ($a->calend - 1) / $bin_size) - $data->{f_offset};
	foreach ($s .. $e) {
		$data->{f}->[$_] += $score;
	}
}

sub se_strand_span {
	my ($a, $data, $score) = @_;
	my $s = int($a->pos / $bin_size);
	my $e = int( ($a->calend - 1) / $bin_size);
	if ($a->reversed) {
		foreach ($s - $data->{r_offset} .. $e - $data->{r_offset}) {
			$data->{r}->[$_] += $score;
		}
	}
	else {
		foreach ($s - $data->{f_offset} .. $e - $data->{f_offset}) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_shift_span {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $start = int( ($a->pos - $shift_value) / $bin_size);
		my $end = int( ($a->calend - 1 - $shift_value) / $bin_size);
		$start = 0 if $start < 0;
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset} ) {
			$data->{f}->[$_] += $score if $_ >= 0;
		}
	}
	else {
		my $start = int( ($a->pos + $shift_value) / $bin_size);
		my $end = int( ($a->calend - 1 + $shift_value) / $bin_size);
		$start = 0 if $start < 0;
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset} ) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_shift_strand_span {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $start = int( ($a->pos - $shift_value) / $bin_size);
		$start = 0 if $start < 0;
		my $end = int( ($a->calend - 1 - $shift_value) / $bin_size);
		foreach ($start - $data->{r_offset} .. $end - $data->{r_offset} ) {
			$data->{r}->[$_] += $score if $_ >= 0;
		}
	}
	else {
		my $start = int( ($a->pos + $shift_value) / $bin_size);
		my $end = int( ($a->calend - 1 + $shift_value) / $bin_size);
		$start = 0 if $start < 0;
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset} ) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_center_span {
	my ($a, $data, $score) = @_;	
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	my $start = $mid - $half_extend + 1;
	$start = 0 if $start < 0;
	$start = int($start/$bin_size);
	my $end = int( ($mid + $half_extend) / $bin_size);
	foreach ($start - $data->{f_offset} .. $end - $data->{f_offset}) {
		$data->{f}->[$_] += $score;
	}
}

sub se_strand_center_span {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	my $start = $mid - $half_extend + 1;
	$start = 0 if $start < 0;
	$start = int($start/$bin_size);
	my $end = int( ($mid + $half_extend) / $bin_size);
	if ($a->reversed) {
		foreach ($start - $data->{r_offset} .. $end - $data->{r_offset}) {
			$data->{r}->[$_] += $score;
		}
	}
	else {
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset}) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_shift_center_span{
	my ($a, $data, $score) = @_;	
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	if ($a->reversed) {
		$mid -= $shift_value;
	}
	else {
		$mid += $shift_value;
	}
	my $start = $mid - $half_extend + 1;
	$start = 0 if $start < 0;
	$start = int($start/$bin_size);
	my $end = int( ($mid + $half_extend) / $bin_size);
	foreach ($start - $data->{f_offset} .. $end - $data->{f_offset}) {
		$data->{f}->[$_] += $score;
	}
}

sub se_shift_strand_center_span{
	my ($a, $data, $score) = @_;	
	my $mid = int( ($a->pos + $a->calend -1) / 2);
	if ($a->reversed) {
		$mid -= $shift_value;
	}
	else {
		$mid += $shift_value;
	}
	my $start = $mid - $half_extend + 1;
	$start = 0 if $start < 0;
	$start = int($start/$bin_size);
	my $end = int( ($mid + $half_extend) / $bin_size);
	if ($a->reversed) {
		foreach ($start - $data->{r_offset} .. $end - $data->{r_offset}) {
			$data->{r}->[$_] += $score;
		}
	}
	else {
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset}) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_extend {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $start = $a->calend - $extend_value;
		$start = 0 if $start < 0;
		$start = int($start / $bin_size);
		my $end = int( ($a->calend - 1) / $bin_size);
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset}) {
			$data->{f}->[$_] += $score;
		}
	}
	else {
		my $start = int($a->pos / $bin_size);
		my $end = int( ($a->pos + $extend_value - 1) / $bin_size);
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset}) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_strand_extend {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $start = $a->calend - $extend_value;
		$start = 0 if $start < 0;
		$start = int($start / $bin_size);
		my $end = int( ($a->calend - 1) / $bin_size);
		foreach ($start - $data->{r_offset} .. $end - $data->{r_offset}) {
			$data->{r}->[$_] += $score;
		}
	}
	else {
		my $start = int($a->pos / $bin_size);
		$start = 0 if $start < 0;
		my $end = int( ($a->pos + $extend_value - 1) / $bin_size);
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset}) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_shift_extend {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $start = $a->calend - $shift_value - $extend_value;
		$start = 0 if $start < 0; 
		$start = int($start / $bin_size);
		my $end = int( ($a->calend - 1 - $shift_value) / $bin_size);
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset} ) {
			$data->{f}->[$_] += $score;
		}
	}
	else {
		my $start = int( ($a->pos + $shift_value) / $bin_size);
		$start = 0 if $start < 0; 
		my $end = int( ($a->pos + $shift_value + $extend_value - 1) / $bin_size);
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset} ) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub se_shift_strand_extend {
	my ($a, $data, $score) = @_;
	if ($a->reversed) {
		my $start = $a->calend - $shift_value - $extend_value;
		$start = 0 if $start < 0; 
		$start = int($start / $bin_size);
		my $end = int( ($a->calend - 1 - $shift_value) / $bin_size);
		foreach ($start - $data->{r_offset} .. $end - $data->{r_offset} ) {
			$data->{r}->[$_] += $score;
		}
	}
	else {
		my $start = int( ($a->pos + $shift_value) / $bin_size);
		$start = 0 if $start < 0; 
		my $end = int( ($a->pos + $shift_value + $extend_value - 1) / $bin_size);
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset} ) {
			$data->{f}->[$_] += $score;
		}
	}
}

sub pe_start {
	my ($a, $data, $score) = @_;
	my $flag = $a->flag;
	# we always receive the forward read, never reverse, from pe_callback
	# therefore the first and second read flag indicates orientation
	if ($flag & 0x0040) {
		# first read implies forward orientation
		$data->{f}->[int($a->pos / $bin_size) - $data->{f_offset}] += $score;
	}
	elsif ($flag & 0x0080) {
		# second read implies reverse orientation
		my $pos = int(($a->pos + $a->isize) / $bin_size);
		$data->{f}->[$pos - $data->{f_offset}] += $score;
	}
	else {
		die " Paired-end flags are set incorrectly; neither 0x040 or 0x080 are set for paired-end.";
	}
}

sub pe_strand_start {
	my ($a, $data, $score) = @_;
	my $flag = $a->flag;
	# we always receive the forward read, never reverse, from pe_callback
	# therefore the first and second read flag indicates orientation
	if ($flag & 0x0040) {
		# first read implies forward orientation
		$data->{f}->[int($a->pos / $bin_size) - $data->{f_offset}] += $score;
	}
	elsif ($flag & 0x0080) {
		# second read implies reverse orientation
		my $pos = int(($a->pos + $a->isize) / $bin_size);
		$data->{r}->[$pos - $data->{r_offset}] += $score;
	}
	else {
		die " Paired-end flags are set incorrectly; neither 0x040 or 0x080 are set for paired-end.";
	}
}

sub pe_mid {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + int( $a->isize / 2 ) ) / $bin_size);
	$data->{f}->[$mid - $data->{f_offset}] += $score;
}

sub pe_strand_mid {
	my ($a, $data, $score) = @_;
	my $mid = int( ($a->pos + int( $a->isize / 2 ) ) / $bin_size);
	my $flag = $a->flag;
	# we always receive the forward read, never reverse, from pe_callback
	# therefore the first and second read flag indicates orientation
	if ($flag & 0x0040) {
		# first read, forward
		$data->{f}->[$mid - $data->{f_offset}] += $score;
	}
	elsif ($flag & 0x0080) {
		# second read, reverse
		$data->{r}->[$mid - $data->{r_offset}] += $score;
	} 
	else {
		die " Paired-end flags are set incorrectly; neither 0x040 or 0x080 are set for paired-end.";
	}
}

sub pe_span {
	my ($a, $data, $score) = @_;
	my $start = int($a->pos / $bin_size);
	my $end = int( ($a->pos + $a->isize) / $bin_size);
	foreach ($start - $data->{f_offset} .. $end - $data->{f_offset}) {
		$data->{f}->[$_] += $score;
	}
}

sub pe_strand_span {
	my ($a, $data, $score) = @_;
	my $start = int($a->pos / $bin_size);
	my $end = int( ($a->pos + $a->isize) / $bin_size);
	my $flag = $a->flag;
	# we always receive the forward read, never reverse, from pe_callback
	# therefore the first and second read flag indicates orientation
	if ($flag & 0x0040) {
		# first read, forward
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset} ) {
			$data->{f}->[$_] += $score;
		}
	}
	elsif ($flag & 0x0080) {
		# second read, reverse
		foreach ($start - $data->{r_offset} .. $end - $data->{r_offset} ) {
			$data->{r}->[$_] += $score;
		}
	} 
	else {
		die " Paired-end flags are set incorrectly; neither 0x040 or 0x080 are set for paired-end.";
	}
}

sub pe_center_span {
	my ($a, $data, $score) = @_;
	my $position = $a->pos + int( $a->isize / 2 );
	my $start = int( ($position - $half_extend) / $bin_size);
	$start = 0 if $start < 0;
	my $end = int( ($position + $half_extend) / $bin_size);
	foreach ($start - $data->{f_offset} .. $end - $data->{f_offset} ) {
		$data->{f}->[$_] += $score;
	}
}

sub pe_strand_center_span {
	my ($a, $data, $score) = @_;
	my $position = $a->pos + int( $a->isize / 2 );
	my $start = int( ($position - $half_extend) / $bin_size);
	$start = 0 if $start < 0;
	my $end = int( ($position + $half_extend) / $bin_size);
	my $flag = $a->flag;
	# we always receive the forward read, never reverse, from pe_callback
	# therefore the first and second read flag indicates orientation
	if ($flag & 0x0040) {
		# first read, forward
		foreach ($start - $data->{f_offset} .. $end - $data->{f_offset} ) {
			$data->{f}->[$_] += $score;
		}
	}
	elsif ($flag & 0x0080) {
		# second read, reverse
		foreach ($start - $data->{r_offset} .. $end - $data->{r_offset} ) {
			$data->{r}->[$_] += $score;
		}
	} 
	else {
		die " Paired-end flags are set incorrectly; neither 0x040 or 0x080 are set for paired-end.";
	}
}

sub smart_pe {
	my ($f, $data, $score, $r) = @_;
	my $set = Set::IntSpan::Fast->new;
	
	# process both reads, adding to the integer set
	foreach my $a ($f, $r) {
		if ($a->cigar_str =~ /N/) {
			my $aw = $wrapper_ref->new($a, $data->{sam});
			
			# check intron size from the cigar string
			if ($max_intron) {
				my $size = 1;
				my $cigars = $aw->cigar_array;
				foreach my $c (@$cigars) {
					# each element is [operation, size]
					$size = $c->[1] if ($c->[0] eq 'N' and $c->[1] > $size);
				}
				return if $size > $max_intron; # exceed maximum intron size
			}
			
			# record
			foreach my $segment ($aw->get_SeqFeatures) {
				$set->add_range(
					int( ($segment->start - 1) / $bin_size), 
					int( ($segment->end - 1) / $bin_size)
				);
			}
		}
		else {
			$set->add_range( int($a->pos / $bin_size), int(($a->calend - 1) / $bin_size) );
		}
	}
	
	# record
	foreach my $pos ($set->as_array) {
		$data->{f}->[$pos - $data->{f_offset}] += $score;
	}
}

sub smart_stranded_pe {
	my ($f, $data, $score, $r) = @_;
	my $set = Set::IntSpan::Fast->new;
	
	# determine strand based on the forward alignment
	my ($strand, $offset);
	my $flag = $f->flag;
	if ($flag & 0x0040) {
		# it's the first read, therefore fragment is forward
		$strand = 'f';
		$offset = 'f_offset';
	}
	elsif ($flag & 0x0080) {
		# it's the second read, therefore fragment is reverse
		$strand = 'r';
		$offset = 'r_offset';
	}
	
	# process both reads, adding to the integer set
	foreach my $a ($f, $r) {
		if ($a->cigar_str =~ /N/) {
			my $aw = $wrapper_ref->new($a, $data->{sam});
			
			# check intron size from the cigar string
			if ($max_intron) {
				my $size = 1;
				my $cigars = $aw->cigar_array;
				foreach my $c (@$cigars) {
					# each element is [operation, size]
					$size = $c->[1] if ($c->[0] eq 'N' and $c->[1] > $size);
				}
				return if $size > $max_intron; # exceed maximum intron size
			}
			
			# record
			foreach my $segment ($aw->get_SeqFeatures) {
				$set->add_range(
					int( ($segment->start - 1) / $bin_size), 
					int( ($segment->end - 1) / $bin_size)
				);
			}
		}
		else {
			$set->add_range( int($a->pos / $bin_size), int(($a->calend - 1) / $bin_size) );
		}
	}
	
	# record
	foreach my $pos ($set->as_array) {
		$data->{$strand}->[$pos - $data->{$offset}] += $score;
	}
}

sub pe_ends {
	my ($a, $data, $score) = @_;
	# we always receive the forward read, never reverse, from pe_callback
	
	# first position
	$data->{f}->[int($a->pos / $bin_size) - $data->{f_offset}] += $score;
	
	# second position
	my $pos = int(($a->pos + $a->isize) / $bin_size);
	$data->{f}->[$pos - $data->{f_offset}] += $score;
}

sub pe_strand_ends {
	my ($a, $data, $score) = @_;
	# we always receive the forward read, never reverse, from pe_callback
	# therefore the first and second read flag indicates orientation
	
	# calculate both positions
	my $flag = $a->flag;
	if ($flag & 0x0040) {
		# first read implies forward orientation
		# first position is forward
		$data->{f}->[int($a->pos / $bin_size) - $data->{f_offset}] += $score;
		# second position is reverse
		my $pos = int(($a->pos + $a->isize) / $bin_size);
		$data->{r}->[$pos - $data->{r_offset}] += $score;
	}
	elsif ($flag & 0x0080) {
		# second read implies reverse orientation
		# first position is reverse
		$data->{r}->[int($a->pos / $bin_size) - $data->{r_offset}] += $score;
		# second position is forward
		my $pos = int(($a->pos + $a->isize) / $bin_size);
		$data->{f}->[$pos - $data->{f_offset}] += $score;
	}
	else {
		die " Paired-end flags are set incorrectly; neither 0x040 or 0x080 are set for paired-end.";
	}
}

__END__

=head1 NAME

bam2wig.pl

A program to convert Bam alignments into a wig representation file.

=head1 SYNOPSIS

bam2wig.pl [--options...] E<lt>file.bamE<gt>

bam2wig.pl --extend --rpm --mean --out file --bw file1.bam file2.bam
  
 Required options:
  -i --in <filename.bam>           repeat if multiple bams, or comma-delimited list
 
 Reporting options (pick one):
  -s --start                    record at 5' position
  -d --mid                      record at midpoint of alignment or pair
  -a --span                     record across entire alignment or pair
  -e --extend                   extend alignment (record predicted fragment)
  --cspan                       record a span centered on midpoint
  --smartcov                    record paired coverage without overlaps, splices
  --ends                        record paired endpoints
  --coverage                    raw alignment coverage
 
 Alignment reporting options:
  -l --splice                   split alignment at N splices
  -t --strand                   record separate strands as two wig files
  --flip                        flip the strands for convenience
  
 Paired-end alignments:
  -p --pe                       process paired-end alignments, both are checked
  -P --fastpe                   process paired-end alignments, only F are checked
  --minsize <integer>           minimum allowed insertion size (30)
  --maxsize <integer>           maximum allowed insertion size (600)
  --first                       only process paired first read (0x40) as single-end
  --second                      only process paired second read (0x80) as single-end
  
 Alignment filtering options:
  -K --chrskip <regex>          regular expression to skip chromosomes
  -B --blacklist <file>         interval file of regions to skip (bed, gff, txt)
  -q --qual <integer>           minimum mapping quality (0)          
  -S --nosecondary              ignore secondary alignments (false)
  -D --noduplicate              ignore marked duplicate alignments (false)
  -U --nosupplementary          ignore supplementary alignments (false)
  --intron <integer>            maximum allowed intron size in bp (none)
  
  Shift options:
  -I --shift                    shift reads in the 3' direction
  -x --extval <integer>         explicit extension size in bp (default is to calculate)
  -H --shiftval <integer>       explicit shift value in bp (default is to calculate) 
  --chrom <integer>             number of chromosomes to sample (4)
  --minr <float>                minimum pearson correlation to calculate shift (0.5)
  --zmin <float>                minimum z-score from average to test peak for shift (3)
  --zmax <float>                maximum z-score from average to test peak for shift (10)
  -M --model                    write peak shift model file for graphing
  
 Score options:
  -r --rpm                      scale depth to Reads Per Million mapped
  -m --mean                     average multiple bams (default is addition)
  --scale <float>               explicit scaling factor, repeat for each bam file
  --fraction                    assign fractional counts to multi-mapped alignments
  --splfrac                     assign fractional count to each spliced segment
  --format <integer>            number of decimal positions (4)
  --chrnorm <float>             use chromosome-specific normalization factor
  --chrapply <regex>            regular expression to apply chromosome-specific factor
 
 Output options:
  -o --out <filename>           output file name, default is bam file basename
  -b --bw                       convert to bigWig format
  --bin <integer>               bin size for span or extend mode (10)
  --bwapp /path/to/wigToBigWig  path to external converter
  -z --gz                       gzip compress output
  
 Wig format:
  --bdg                         bedGraph, default for span and extend at bin 1
  --fix                         fixedStep, default for bin > 1
  --var                         varStep, default for start, mid
  --nozero                      do not write zero score intervals in bedGraph
  
 General options:
  -c --cpu <integer>            number of parallel processes (4)
  --temp <directory>            directory to write temporary files (output path)
  -V --verbose                  report additional information
  -v --version                  print version information
  -h --help                     show full documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 Input

=over 4

=item --in E<lt>filenameE<gt>

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

=item --smartcov

Smart alignment coverage of paired-end alignments without 
double-counting overlaps or recording gaps (intron splices). 

=item --ends

Record both endpoints of paired-end fragments, i.e. the outermost 
or 5' ends of properly paired fragments. This may be useful with 
ATAC-Seq, Cut&Run-Seq, or other cleavage experiments where you want 
to record the locations of cutting yet retain the ability to filter 
paired-end fragment sizes.

=item --coverage

Specify that the raw alignment coverage be calculated and reported 
in the wig file. This utilizes a special low-level operation and 
precludes any alignment filtering or post-normalization methods. 
Counting overlapping bases in paired-end alignments are dependent on 
the bam adapter (older versions would double-count).

=item --position [start|mid|span|extend|cspan|coverage]

Legacy option for supporting previous versions of bam2wig. 

=back

=head2 Alignment reporting options

=over 4

=item --splice

Indicate that the bam file contains alignments with splices, such as 
from RNASeq experiments. Alignments will be split on cigar N operations 
and each sub fragment will be recorded. This only works with single-end 
alignments, and is disabled for paired-end reads (just treat as single-end). 
Only start and span recording options are supported.

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
pairs aligning to separate chromosomes. Both alignments are required 
to be present before the pair is counted. The default is to treat 
all alignments as single-end.

=item --fastpe

The Bam file consists of paired-end alignments, but to increase processing 
time and be more tolerant of weird pairings, only the forward alignment is 
required and considered; all reverse alignments are ignored. The default is 
to treat all alignments as single-end.

=item --minsize E<lt>integerE<gt>

Specify the minimum paired-end fragment size in bp to accept for recording. 
Default is 30 bp.

=item --maxsize E<lt>integerE<gt>

Specify the maximum paired-end fragment size in bp to accept for recording. 
Default is 600 bp.

=item --first

Take only the first read of a pair, indicated by flag 0x40, and record as 
a single-end alignment. No test of insert size or proper pair status is 
made.

=item --second

Take only the second read of a pair, indicated by flag 0x80, and record as 
a single-end alignment. No test of insert size or proper pair status is 
made.

=back

=head2 Alignment filtering options:

=over 4

=item --qual E<lt>integerE<gt>

Set a minimum mapping quality score of alignments to count. The mapping 
quality is a range from 0-255, with higher numbers indicating lower 
probability of a mapping error. Multi-mapping alignments often have a 
map quality of 0. The default is 0 (accept everything).

=item --nosecondary

Boolean flag to skip secondary alignments, indicated by the 
alignment bit flag 0x100. Secondary alignments typically represent 
alternative mapping locations, or multi-mapping events. By default,  
secondary alignments are included. 

=item ---noduplicate

Boolean flag to skip duplicate alignments, indicated by the 
alignment bit flag 0x400. Duplicates alignments may represent a PCR or 
optical duplication. By default, duplicate alignments are included. 

=item --nosupplementary

Boolean flag to skip supplementary alignments, indicated by 
the alignment bit flag 0x800. Supplementary alignments are typically 
associated with chimeric fragments. By default, supplementary alignments 
are included.

=item --chrskip E<lt>regexE<gt>

Provide a regular expression to skip certain chromosomes. Perl-based 
regular expressions are employed. Expressions should be quoted or 
properly escaped on the command line. Examples might be 
    
    'chrM'
    'scaffold.+'
    'chr.+alt|chrUn.+|chr.+_random'

=item --blacklist E<lt>fileE<gt>

Provide a file of genomic intervals from which to exclude alignments. 
Examples might include repeats, ribosomal RNA, or heterochromatic regions.
The file should be any text file interpretable by L<Bio::ToolBox::Data> 
with chromosome, start, and stop coordinates, including BED and GFF formats.
Note that this only excludes overlapping alignments, and does not include 
extended alignments.

=item --intron E<lt>integerE<gt>

Provide a positive integer as the maximum intron size allowed in an alignment 
when splitting on splices. If an N operation in the CIGAR string exceeds this 
limit, the alignment is skipped. Default is 0 (no filtering).

=back

=head2 Shift options

=over 4

=item --shift

Specify that the positions of the alignment should be shifted towards 
the 3' end. Useful for ChIP-Seq applications, where only the ends of 
the fragments are counted and often seen as separated discrete peaks 
on opposite strands flanking the true target site. This option is 
disabled with paired-end and spliced reads (where it is not needed). 

=item --shiftval E<lt>integerE<gt>

Provide the value in bp that the recorded position should be shifted. 
The value should be 1/2 the average length of the library insert size.
The default is to automatically and empirically determine the 
appropriate shift value using cross-strand correlation (recommended). 

=item --extval E<lt>integerE<gt>

Manually set the length for reads to be extended. By default, the shift 
value is determined empirically and extension is set to 2X the shift 
value. This is also used for the cspan mode.

=item --chrom E<lt>integerE<gt>

Indicate the number of sequences or chromosomes to sample when 
empirically determining the shift value. The reference sequences 
listed in the Bam file header are taken in order of decreasing 
length, and one or more are taken as a representative sample of 
the genome. The default value is 4. 

=item --minr E<lt>floatE<gt>

Provide the minimum Pearson correlation value to accept a shift 
value when empirically determining the shift value. Enter a decimal value 
between 0 and 1. Higher values are more stringent. The default 
is 0.5.

=item --zmin E<lt>floatE<gt>

Specify the minimum z-score (or number of standard deviations) from 
the chromosomal mean depth to test for a peak shift. Increase this 
number to test for strong robust peaks, which give a better estimations 
of the shift value. Default is 3.

=item --zmax E<lt>floatE<gt> 

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

=item --scale E<lt>floatE<gt>

Optionally provide your own scaling factor. This will be multiplied with 
every position when generating the wig file. This may be combined with the 
rpm factor as well. When combining multiple bam files, either a single scale 
factor may be supplied for all files, or individual scale factors may be 
supplied for each bam file. If supplying multiple, use the option multiple 
times or give a comma-delimited list. The values should be in the same order 
as the bam files. 

=item --mean

=item --separate

When processing multiple bam files, this option will take the mean or average 
across all bam files. Without this option, the bam files are simply added. 
When combined with the rpm option, each bam file will be scaled separately 
before taking the average.  

=item --fraction

Indicate that multi-mapping alignments should be given fractional counts 
instead of full counts. The number of alignments is determined using the 
NH alignment tag. If a read has 10 alignments, then each alignment is 
given a count of 0.1. 

=item --splfrac

Indicate that spliced segments should be given a fractional count. This 
allows a count to be assigned to each spliced segment while avoiding 
double-counting. Best used with RNASeq spliced point data (--start or 
--mid); not recommended for --span.

=item --format E<lt>integerE<gt>

Indicate the number of decimal postions reported in the wig file. This 
is only applicable when rpm, scale, or fraction options are provided. 
The default value is 4 decimal positions.

=item --chrnorm E<lt>floatE<gt>

Apply a normalization factor to the counts on specific chromosomes only. 
Usually this is to normalize, for example, variable copy-number chromosomes, 
such as transfected vector sequences, or haploid sex chromosomes. 

=item --chrapply E<lt>regexE<gt>

Specify the Perl-based regular expression to match the chromosome(s) to 
apply the specific normalization factor. For example, 'chrX$' to specify 
the X chromosome only.

=back

=head2 Output Options

=over 4

=item --out E<lt>filenameE<gt>

Specify the output base filename. An appropriate extension will be 
added automatically. By default it uses the base name of the 
input file.

=item --bin E<lt>integerE<gt>

Specify the bin size in bp for the output wig file. In general, specifying 
a larger bin size will decrease the run time and memory requirements in 
exchange for loss of resolution. The default for span, center span, or 
extend mode is 10 bp; all other modes is 1 bp. 

=item --bw

Specify whether or not the wig file should be further converted into 
an indexed, compressed, binary BigWig file. The default is false.

=item --bwapp /path/to/wigToBigWig

Optionally specify the full path to the UCSC I<wigToBigWig> conversion 
utility. The application path may be set in the F<.biotoolbox.cfg> file 
or found in the default executable path, which makes this option 
unnecessary. 

=item --gz

Specify whether (or not) the output file should be compressed with 
gzip. Disable with --nogz.
=item --nozero

When writing bedGraph format, skip (do not write) intervals with a value of 
zero. Does not apply to fixedStep or variableStep formats.

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

=item --cpu E<lt>integerE<gt>

Specify the number of parallel instances to run simultaneously. This requires 
the installation of L<Parallel::ForkManager>. With support enabled, the 
default is 4. Disable multi-threaded execution by setting to 1. 

=item --temp E<lt>directoryE<gt>

Optionally specify an alternate temporary directory path where the temporary 
files will be written. The default is the specified output file path, or the 
current directory. Temporary files will always be written in a subdirectory of 
the path specified with the template "bam2wigTEMP_XXXX".

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
browsers display when using a Bam file as source. B<NOTE> that this mode 
is pure raw coverage, and does not include any filtering methods. The other 
modes allow alignment filtering.
 
 bam2wig.pl --coverage --in <bamfile>

=item Smart paired-end coverage

When you have paired-end alignments and need explicit alignment coverage
without double-counting overlaps (as would occur if you counted as
single-end span) or uncovered insertion (as would occur if you counted as 
paired-end span) and not counting gaps (e.g. intron splices, as would occur
with span mode), use the smart paired-end coverage mode. This properly
assembles coverage from paired-end alignments taking into account overlaps
and gaps.

 bam2wig --smartcov --in <bamfile>

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
 
 bam2wig --span ---splice --strand --rpm --in <bamfile>

 bam2wig --pos mid --strand --rpm --in <bamfile>
 
=item Paired-end RNA-Seq

Use the smart paired-end coverage mode to properly record paired-end 
alignments with splice junctions. 

 bam2wig --smartcov --strand --rpm --in <bamfile>
 
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
