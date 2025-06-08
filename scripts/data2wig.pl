#!/usr/bin/perl

# documentation at end of file

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use List::Util       qw(sum0 max);
use Scalar::Util     qw(looks_like_number);
use Statistics::Lite qw(median);
use Bio::ToolBox::Data::Stream;
use Bio::ToolBox::utility    qw(parse_list ask_user_for_index);
use Bio::ToolBox::big_helper qw(
	get_wig_to_bigwig_app
	check_wigToBigWig_version
	open_wig_to_bigwig_fh
	generate_chromosome_file
	wig_to_bigwig_conversion
);

our $VERSION = '2.02';

print "\n This script will export a data file to a wig file\n\n";

### Quick help
unless (@ARGV) {

	# when no command line options are present
	# print SYNOPSIS
	pod2usage(
		{
			'-verbose' => 0,
			'-exitval' => 1,
		}
	);
}

### Get command line options and initialize values
my (
	$infile,        $outfile,     $fast,           $step,
	$bedgraph,      $step_size,   $span,           $ask,
	$chr_index,     $start_index, $stop_index,     $score_index,
	@score_indices, $no_header,   $attribute_name, $track_name,
	$use_track,     $midpoint,    $interbase,      $format,
	$method,        $bigwig,      $bw_app_path,    $database,
	$chromo_file,   $gz,          $help,           $print_version,
);

# Command line options
GetOptions(
	'i|in=s'              => \$infile,            # name of input file
	'o|out=s'             => \$outfile,           # name of output gff file
	'f|fast!'             => \$fast,              # fast mode without checks and stuff
	'p|step=s'            => \$step,              # wig step method
	'bed|bdg!'            => \$bedgraph,          # write a bedgraph file
	'size=i'              => \$step_size,         # wig step size
	'span=i'              => \$span,              # the wig span size
	'a|ask'               => \$ask,               # request help in assigning indices
	'c|chr=i'             => \$chr_index,         # index for the chromosome column
	'b|begin|start|pos=i' => \$start_index,       # index for the start column
	'e|stop|end=i'        => \$stop_index,        # index for the stop column
	's|index|score=s'     => \$score_index,       # index for the score column
	'H|noheader'          => \$no_header,         # source has no header line
	'attrib=s'            => \$attribute_name,    # gff or vcf attribute name to use
	'name=s'              => \$track_name,        # name string for the track
	'track!'              => \$use_track,         # boolean to include a track line
	'mid!'                => \$midpoint,          # boolean to use the midpoint
	'0|zero|inter!'       => \$interbase,         # shift from interbase
	'format=i'    => \$format,           # format output to indicated number of places
	'm|method=s'  => \$method,           # method for combining duplicate values
	'B|bigwig|bw' => \$bigwig,           # generate a binary bigwig file
	'd|db=s'      => \$database,         # database for bigwig file generation
	'chromof=s'   => \$chromo_file,      # name of a chromosome file
	'bwapp=s'     => \$bw_app_path,      # path to wigToBigWig utility
	'z|gz!'       => \$gz,               # boolean to compress output file
	'h|help'      => \$help,             # request help
	'v|version'   => \$print_version,    # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {

	# print entire POD
	pod2usage(
		{
			'-verbose' => 2,
			'-exitval' => 1,
		}
	);
}

# Print version
if ($print_version) {
	print " Biotoolbox script data2wig.pl, version $VERSION\n\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

### Check for required or default values
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		print " FATAL: No source data file specified! \n use data2wig.pl --help\n";
		exit 1;
	}
}
unless ( defined $use_track ) {

	# default is to write a track
	$use_track = $bigwig ? 0 : 1;
}
$method ||= 'mean';

### Load input file
my $Input = Bio::ToolBox::Data::Stream->new( file => $infile, noheader => $no_header )
	or die "Unable to open file '$infile'!\n";

### Check and/or ask for specific options

check_indices();

check_track_name();

check_step();

my $post_bw_convert;
if ($bigwig) {
	set_bigwig_options();
}

my $method_sub = set_method_sub();

my $score_sub = set_score_sub();

my $printer = set_print_string();

if ($fast) {
	if ($midpoint) {
		print
" WARNING: cannot use midpoint position in fast mode!\nrunning in slow mode...\n";
		$fast = 0;
	}
	if ($attribute_name) {
		print
" WARNING: cannot use GFF or VCF attribute score in fast mode!\nrunning in slow mode...\n";
		$fast = 0;
	}
}

my $start_time = time;

### Open output file
unless ($outfile) {

	# automatically generate output file name based on track name
	$outfile = $track_name;
	$outfile =~ s/\(\) /_/g;    # strip parentheses and spaces from column name
}
my $out_fh;
if ( $bigwig and not $post_bw_convert ) {

	# we will write directly to a bigWig file
	unless ($chromo_file) {
		$chromo_file = generate_chromosome_file($database)
			or die "unable to generate chromosome file needed for bigWig conversion!\n";
	}
	$outfile .= '.bw' unless $outfile =~ /\.bw$/;
	$out_fh = open_wig_to_bigwig_fh(
		file      => $outfile,
		chromo    => $chromo_file,
		bwapppath => $bw_app_path
	) or die "unable to open a wigToBigWig file handle!\n";
}
else {
	# we will write to a wig file
	unless ( $outfile =~ /\.(?: wig | bdg | bedgraph ) (?:\.gz)? $/xi ) {
		$outfile .= $bedgraph ? '.bdg' : '.wig';
	}
	$out_fh = Bio::ToolBox::Data::Stream->open_to_write_fh( $outfile, $gz )
		or die " unable to open output file '$outfile' for writing!\n";

	# write track line
	if ($use_track) {
		print {$out_fh} "track type=wiggle_0 name=$track_name\n";
	}
}

### Start the conversion
printf " converting '%s' as %s wig file...\n",
	scalar(@score_indices)
	? join( ", ", map { $Input->name($_) } @score_indices )
	: $Input->name($score_index),
	$bedgraph ? 'bedGraph' : $step . 'Step';

if ( $fast and $bedgraph ) {
	fast_convert_to_bedgraph();
}
elsif ( $fast and $step eq 'fixed' ) {
	fast_convert_to_fixedStep();
}
elsif ( $fast and $step eq 'variable' ) {
	fast_convert_to_variableStep();
}
elsif ($bedgraph) {
	convert_to_bedgraph();
}
elsif ( $step eq 'fixed' ) {
	convert_to_fixedStep();
}
elsif ( $step eq 'variable' ) {
	convert_to_variableStep();
}

# close
$out_fh->close;

# convert as necessary
if ( $bigwig and $post_bw_convert ) {

	unless ($chromo_file) {
		$chromo_file = generate_chromosome_file($database)
			or die "unable to generate chromosome file needed for bigWig conversion!\n";
	}
	my $final_bw = wig_to_bigwig_conversion(
		wig       => $outfile,
		chromo    => $chromo_file,
		bwapppath => $bw_app_path,
	);
	if ($final_bw) {
		print " Wrote file $final_bw\n";
		unlink $outfile;
	}

	# error message already printed if conversion had failed
}
else {
	printf " Wrote file %s\n", $outfile;
}

# clean up chromosome file
if ( $bigwig and $database and $chromo_file =~ /^chr_sizes_\w{5}$/x ) {
	unlink $chromo_file;
}

# finished
printf " Finished in %.0f seconds\n", ( time - $start_time );

############ Subroutines ###############

sub check_indices {

	# check coordinates
	if ( not $chr_index ) {
		$chr_index = $Input->chromo_column;
	}
	if ( $ask or not defined $chr_index ) {
		$chr_index =
			ask_user_for_index( $Input, " Enter the index for the chromosome column  " );
		unless ( defined $chr_index ) {
			print " FATAL: No identifiable chromosome column index!\n";
			exit 1;
		}
	}
	if ( not $start_index ) {
		$start_index = $Input->start_column;
	}
	if ( $ask or not $start_index ) {
		$start_index = ask_user_for_index( $Input,
			" Enter the index for the start or position column  " );
		unless ($start_index) {
			print " FATAL: No identifiable start column index!\n";
			exit 1;
		}
	}
	if ( substr( $Input->name($start_index), -1 ) eq '0' ) {

		# name suggests it is 0-based indexed
		$interbase = 1;
	}

	# stop column is optional

	# score
	if ( not $score_index ) {
		$score_index = $Input->gff ? 6 : $Input->bed >= 5 ? 5 : undef;
	}
	if ( $ask or not $score_index ) {

		# first look for a generic score index
		$score_index =
			ask_user_for_index( $Input, " Enter the index for the score column  " );
		unless ($score_index) {
			print " FATAL: No identifiable score column index!\n";
			exit 1;
		}
	}
	if ( $score_index =~ /[,\-]/ ) {
		@score_indices = parse_list($score_index);
	}
	elsif ( ref($score_index) eq 'ARRAY' ) {

		# in case returned from interactive
		@score_indices = @{$score_index};
	}
	if (@score_indices) {
		foreach my $i (@score_indices) {
			unless ( $Input->name($i) ) {
				print " FATAL: Score column $i does not exist!\n";
				exit 1;
			}
		}
	}
	else {
		unless ( $Input->name($score_index) ) {
			print " FATAL: Score column $score_index does not exist!\n";
			exit 1;
		}
	}
}

sub check_track_name {

	# determine what the track name will be
	return if $track_name;
	if (@score_indices) {
		$track_name = join( '_', map { $Input->name($_) } @score_indices );
	}
	elsif ($score_index) {

		$track_name = $Input->name($score_index);
		if ( $track_name =~ /^score$/i ) {

			# some sort of bed or gff file standard score column
			$track_name = $Input->basename;
		}
	}
}

sub check_step {

	# In my biotoolbox scripts that generate genomic bins or windows
	# the win(dow) and step values are recorded as metadata under
	# the start column
	# Therefore, we will check for this metadata to confirm and/or
	# determine the step and span parameters

	# check step
	if ($bedgraph) {

		# write a bedgraph file
		$step      = 'bed';
		$use_track = 0;       # no track data
		return;
	}
	elsif ( $step and $step eq 'bed' ) {

		# write a bedgraph file
		$bedgraph  = 1;
		$use_track = 0;       # no track data
		return;
	}
	elsif ( $step and $step eq 'variable' ) {

		# this is ok, we can work with it
	}
	elsif ( $step and $step eq 'fixed' ) {

		# double check that the data file supports this
		# assign the step size as necessary
		if ( defined $Input->metadata( $start_index, 'step' ) ) {
			if ($step_size) {
				if ( $step_size != $Input->metadata( $start_index, 'step' ) ) {
					print
" FATAL: Requested step size $step_size does not match metadata step size!!!\n";
					exit 1;
				}
			}
			else {
				# define it from the metadata
				$step_size = $Input->metadata( $start_index, 'step' );
			}
		}
		elsif ( not $step_size ) {
			print
" WARNING: Fixed step size not defined by user or metadata! Using 'variableStep'\n";
			$step = 'variable';
		}
	}
	else {
		# attempt to determine automatically
		if ( defined $Input->metadata( $start_index, 'step' ) ) {

			# set step size
			$step      = 'fixed';
			$step_size = $Input->metadata( $start_index, 'step' );
			print
				" Automatically generating 'fixedStep' wig with step of $step_size bp\n";
		}
		else {
			print " Automatically generating 'variableStep' wig\n";
			$step = 'variable';
		}
	}

	# check span
	if ($span) {

		# user set it, confirm with metadata if possible
		if ( defined $Input->metadata( $start_index, 'win' )
			and $Input->metadata( $start_index, 'win' ) != $span )
		{
			# the requested span and metadata window size do not match
			print
" FATAL: Requested span size $span does not match metadata window size!!!\n";
			exit 1;
		}
	}
	else {
		# attempt to determine automatically
		if ( defined $Input->metadata( $start_index, 'win' ) and not $midpoint ) {

			# set the span equal to the window size
			$span = $Input->metadata( $start_index, 'win' );
			print " Automatically setting span to $span bp\n";
		}
		else {
			# default span is 1 bp
			$span = 1;
		}
	}

	# confirm span
	if ( $midpoint and $span > 1 ) {

		# cannot use span parameter if the midpoint is requested
		$span = 1;
		print " Reverting span to 1 bp when midpoint is requested\n";
	}
	if ( $step eq 'fixed' and $span and $span > $step_size ) {

		# span greater than the step will confuse wig parsers,
		# including conversion to bigWig files
		# the wig format definition cannot have > 1 value per position
		print " Span of $span bp is greater than step of $step_size bp!\n"
			. "   Wig cannot have overlapping positions. Switching to midpoint\n";
		$span     = 1;
		$midpoint = 1;
	}
}

sub set_bigwig_options {

	# find external utility
	unless ($bw_app_path) {
		$bw_app_path = get_wig_to_bigwig_app();
	}

	# check the utility and options
	if ($bw_app_path) {
		if ( check_wigToBigWig_version($bw_app_path) ) {

			# we can write directly to utility
			# print " $bw_app_path supports stdin, can write directly\n";
			$post_bw_convert = 0;
		}
		else {
			# we cannot write directly to utility
			# print " $bw_app_path does not support stdin, will write temp wig file\n";
			$post_bw_convert = 1;
			$gz              = 0;
		}

		# if we're generating bigwig file, no track is needed
		$use_track = 0;

		# check that we have a source for chromosome info
		unless ( $database or $chromo_file ) {
			$database = $Input->database
				or die
" No database name or chromosome file provided for generating bigwig file!\n";
		}
	}
	else {
		print " No external bigWig utility available, writing wig format\n";
		$bigwig = 0;
	}
}

sub set_method_sub {

	# for combining values from duplicate positions we need a method
	if ( $method eq 'mean' ) {
		return sub { return sum0(@_) / ( scalar(@_) || 1 ); };
	}
	elsif ( $method eq 'median' ) {
		return \&median;
	}
	elsif ( $method eq 'sum' ) {
		return \&sum0;
	}
	elsif ( $method eq 'max' ) {
		return \&max;
	}
	else {
		print STDERR " FATAL: unrecognized method 'method'!\n";
		exit 1;
	}
}

sub set_score_sub {
	if ( $attribute_name and $Input->gff ) {

		# a GFF attribute
		return sub {
			my $row     = shift;
			my $attribs = $row->gff_attributes;
			my $score   = $attribs->{$attribute_name} || 0;
			return if $score eq '.';

			# format as necessary
			$score =~ s/\%$//;    # remove stupid percents if present
			return $score;
		};
	}
	elsif ( $attribute_name and $Input->vcf and defined $score_index ) {

		# a VCF attribute from one sample
		return sub {
			my $row     = shift;
			my $attribs = $row->vcf_attributes;
			my $score   = $attribs->{$score_index}{$attribute_name} || 0;
			return 0 if $score eq '.';

			# format as necessary
			$score =~ s/\%$//;    # remove stupid percents if present
			return $score;
		};
	}
	elsif ( $attribute_name and $Input->vcf and @score_indices ) {

		# a VCF attribute from many samples
		return sub {
			my $row     = shift;
			my $attribs = $row->vcf_attributes;
			my @scores;
			foreach (@score_indices) {
				my $s = $attribs->{$_}{$attribute_name} || 0;
				$s =~ s/\%$//;    # remove stupid percents if present
				if ( looks_like_number($s) ) {
					push @scores, $s;
				}
			}
			return &{$method_sub}(@scores);
		};
	}
	elsif ( @score_indices and $fast ) {

		# collect over multiple score columns from array reference
		return sub {
			my $data = shift;
			my @v    = grep { looks_like_number($_) } map { $data->[$_] } @score_indices;
			return &{$method_sub}(@v);
		};
	}
	elsif ( @score_indices and not $fast ) {

		# collect over multiple score columns from Feature row object
		return sub {
			my $row = shift;
			my @v =
				grep { looks_like_number($_) } map { $row->value($_) } @score_indices;
			return &{$method_sub}(@v);
		};
	}
	elsif ( defined $score_index and $fast ) {

		# collect from a single score column
		return sub {
			return shift->[$score_index] || 0;
		};
	}
	elsif ( defined $score_index and not $fast ) {

		# collect from a single score column
		return sub {
			my $row = shift;
			return $row->value($score_index) || 0;
		};
	}
	else {
		die "programmer error! how did we get here?\n";
	}
}

sub set_print_string {

	# set the printf format string depending on type of wig file
	if ( $step eq 'fixed' ) {
		return defined $format ? '%.' . $format . "f\n" : "%s\n";
	}
	elsif ( $step eq 'variable' ) {
		return defined $format ? '%d %.' . $format . "f\n" : "%s %s\n";
	}
	elsif ( $step eq 'bed' ) {
		return defined $format ? "%s\t%d\t%d\t%." . $format . "f\n" : "%s\t%d\t%d\t%s\n";
	}
}

sub convert_to_fixedStep {

	# keep track of current chromosome name and length
	my $current_chr  = q();    # current chromosome
	my $previous_pos = 0;

	# walk through the data file
	while ( my $row = $Input->next_row ) {

		# coordinates
		my $chromosome = defined $chr_index ? $row->value($chr_index) : $row->seq_id;
		my $start      = calculate_position($row);

		# check coordinate
		next if $start <= 0;    # skip negative or zero coordinates

		# write definition line if necessary
		if ( $chromosome ne $current_chr ) {

			# new chromosome, new definition line
			$out_fh->printf( "fixedStep chrom=%s start=%d step=%d span=%d\n",
				$chromosome, $start, $step_size, $span );

			# reset the current chromosome
			$current_chr  = $chromosome;
			$previous_pos = $start;        # temporary artificial
		}
		elsif ( $start > ( $previous_pos + $span ) ) {

			# skipped a chunk here
			$out_fh->printf( "fixedStep chrom=%s start=%d step=%d span=%d\n",
				$chromosome, $start, $step_size, $span );
		}
		elsif ( $start < $previous_pos ) {
			printf STDERR
" FATAL: input file is not genomically sorted! %d comes after %d at line %d\n",
				$start, $previous_pos, $row->line_number;
			exit 1;
		}
		$previous_pos = $start;

		# print the score
		$out_fh->printf( $printer, &{$score_sub}($row) );
	}
}

sub fast_convert_to_fixedStep {
	my $current_chr = q();    # current chromosome
	unless ($chr_index) {
		$chr_index = $Input->chromo_column;
	}
	unless ( $chr_index and $start_index ) {
		print STDERR " FATAL: coordinate columns not defined!\n";
		exit 1;
	}
	unless ($score_index) {
		print STDERR " FATAL: no score column defined!\n";
		exit 1;
	}

	# use direct file handle and skip the Stream and row objects
	my $fh = $Input->fh;
	while ( my $line = $fh->getline ) {
		chomp $line;
		my @data = split( /\t/, $line );

		# working with raw array, need to shift values to compensate
		# for 1-base indexing of column indexes
		unshift @data, 0;

		my $chromosome = $data[$chr_index];
		my $start      = $data[$start_index];
		$start++ if $interbase;

		# write definition line if necessary
		if ( $chromosome ne $current_chr ) {
			$out_fh->printf( "fixedStep chrom=%s start=%d step=%d span=%d\n",
				$chromosome, $start, $step_size, $span );
			$current_chr = $chromosome;
		}

		# print the score
		$out_fh->printf( $printer, &{$score_sub}( \@data ) );
	}
}

sub convert_to_variableStep {
	my $current_chr  = q();    # current chromosome
	my $previous_pos = 0;      # previous position to avoid duplicates in wig file
	my @scores;                # reusable array for putting multiple data points in
	while ( my $row = $Input->next_row ) {

		# coordinates
		my $chromosome = $chr_index ? $row->value($chr_index) : $row->seq_id;
		my $start      = calculate_position($row);
		next if $start <= 0;    # skip negative or zero coordinates

		# write definition line if necessary
		if ( $chromosome ne $current_chr ) {

			# first check and write orphan scores
			# this might happen if there was only one score on the entire chr
			if (@scores) {
				if ( scalar @scores == 1 ) {

					# print the one score
					$out_fh->printf( $printer, $previous_pos, $scores[0] );
				}
				else {
					# more than one score, combine them
					my $score = &{$method_sub}(@scores);
					$out_fh->printf( $printer, $previous_pos, $score );
				}
				@scores = ();
			}

			# print new definition line and reset for next
			$out_fh->print("variableStep chrom=$chromosome span=$span\n");
			$current_chr  = $chromosome;
			$previous_pos = $start;
		}

		# collect the score
		my $score = &{$score_sub}($row);
		next unless defined $score;

		# check for duplicate positions and write appropriately
		if ( $start < $previous_pos ) {
			printf STDERR
" FATAL: input file is not genomically sorted! %d comes after %d at line %d\n",
				$start, $previous_pos, $row->line_number;
			exit 1;
		}
		elsif ( $start == $previous_pos ) {

			# same position, add to the score list
			push @scores, $score;
		}
		else {
			# we have moved on to the next position
			# now print the previous scores
			if ( scalar @scores == 1 ) {

				# print the one score
				$out_fh->printf( $printer, $previous_pos, $scores[0] );
			}
			elsif ( scalar @scores > 1 ) {

				# more than one score
				my $multi_score = &{$method_sub}(@scores);
				$out_fh->printf( $printer, $previous_pos, $multi_score );
			}

			# reset for next
			$previous_pos = $start;
			@scores       = ($score);
		}
	}
}

sub fast_convert_to_variableStep {
	my $current_chr = q();    # current chromosome
	unless ($chr_index) {
		$chr_index = $Input->chromo_column;
	}
	unless ( $chr_index and $start_index ) {
		print STDERR " FATAL: coordinate columns not defined!\n";
		exit 1;
	}
	unless ($score_index) {
		print STDERR " FATAL: no score column defined!\n";
		exit 1;
	}

	# use direct file handle and skip the Stream and row objects
	my $fh = $Input->fh;
	while ( my $line = $fh->getline ) {
		chomp $line;
		my @data = split( /\t/, $line );

		# working with raw array, need to shift values to compensate
		# for 1-base indexing of column indexes
		unshift @data, 0;

		my $chromosome = $data[$chr_index];
		my $start      = $data[$start_index];
		$start++ if $interbase;

		# write definition line if necessary
		if ( $chromosome ne $current_chr ) {
			$out_fh->print("variableStep chrom=$chromosome span=$span\n");
			$current_chr = $chromosome;
		}

		# collect the score
		$out_fh->printf( $printer, $start, &{$score_sub}( \@data ) );
	}
}

sub convert_to_bedgraph {

	# variables to check for overlap
	my $current_chr  = q();    # current chromosome
	my $previous_pos = 0;      # previous position to avoid overlap
	while ( my $row = $Input->next_row ) {

		# coordinates
		my $chromosome = $chr_index ? $row->value($chr_index) : $row->seq_id;
		my $start      = $row->value($start_index);
		my $stop =
			defined $stop_index
			? $row->value($stop_index)
			: $row->stop || $row->start || $start;

		# adjust start position
		$start-- unless ($interbase);

		# check coordinates
		if ( defined $previous_pos and defined $current_chr ) {

			# check if on the same chromosome
			if ( $current_chr eq $chromosome ) {

				# check for overlap
				if ( $start < $previous_pos ) {
					print
" WARNING: There are overlapping intervals or the file is not sorted by"
						. " coordinates!\n Compare $chromosome:$start"
						. " with previous stop position $previous_pos\n";

					# die if bigwig
					if ($bigwig) {
						print STDERR
" FATAL: bigWig conversion will fail with overlapping coordinates!\n Fix your file!\n";
						exit 1;
					}
				}

				# otherwise it is ok
				$previous_pos = $stop;
			}
			else {
				# new chromosome
				$current_chr  = $chromosome;
				$previous_pos = $stop;
			}
		}
		else {
			$current_chr  = $chromosome;
			$previous_pos = $stop;
		}

		# collect the score
		my $score = &{$score_sub}($row);
		next unless defined $score;

		# write the feature line
		$out_fh->printf( $printer, $chromosome, $start, $stop, $score );
	}
}

sub fast_convert_to_bedgraph {

	unless ($chr_index) {
		$chr_index = $Input->chromo_column;
	}
	unless ($stop_index) {
		$stop_index = $Input->end_column;
	}
	unless ( $chr_index and $start_index and $stop_index ) {
		print " WARNING: coordinate columns not defined!\n";
		exit 1;
	}
	unless ($score_index) {
		print "WARNING: no score column defined!\n";
		exit 1;
	}

	# use direct file handle and skip the Stream and row objects
	my $fh = $Input->fh;
	while ( my $line = $fh->getline ) {
		chomp $line;
		my @data = split( /\t/, $line );

		# working with raw array, need to shift values to compensate
		# for 1-base indexing of column indexes
		unshift @data, 0;

		# adjust start position
		$data[$start_index]-- unless ($interbase);

		# write the feature line
		$out_fh->printf( $printer, $data[$chr_index], $data[$start_index],
			$data[$stop_index], &{$score_sub}( \@data ) );
	}
}

sub calculate_position {
	my $row   = shift;
	my $start = $row->value($start_index);
	$start += 1 if $interbase;
	if ($midpoint) {

		# calculate midpoint and return
		my $end = $stop_index ? $row->value($stop_index) : $row->end;
		$end ||= $start;    # in case no end was defined
		return $start == $end ? $start : int( ( ( $start + $end ) / 2 ) + 0.5 );
	}
	else {
		return $start;
	}
}

__END__

=head1 NAME

data2wig.pl

A program to convert a generic data file into a wig file.

=head1 SYNOPSIS

data2wig.pl [--options...] <filename> 
  
  File options:
  -i --in <filename>                    input file: txt, gff, bed, vcf, etc
  -o --out <filename>                   output file name
  -H --noheader                         input file has no header row
  -0 --zero                             file is in 0-based coordinate system
  
  Column indices:
  -a --ask                              interactive selection of columns
  -s --score <index>                    score column, may be comma list
  -c --chr <index>                      chromosome column
  -b --begin --start <index>            start coordinate column
  -e --end --stop <index>               stop coordinate column
  --attrib <name>                       GFF or VCF attribute name of score
  
  Wig options:
  -p --step [fixed|variable|bed]        type of wig file 
  --bed --bdg                           alternative shortcut for bedGraph
  --size <integer>                      step size for fixedStep
  --span <integer>                      span size for fixed and variable
  
  Conversion options:
  -f --fast                             fast mode, no error checking
  --name <text>                         optional track name
  --(no)track                           generate a track line
  --mid                                 use the midpoint of feature intervals
  --format <integer>                    format decimal points of scores
  -m --method [mean | median | sum | max] combine multiple score columns
  
  BigWig options:
  -B  --bw --bigwig                     generate a bigWig file
  -d --db <database>                    database to collect chromosome lengths
  --chromof <filename>                  specify a chromosome file
  --bwapp </path/to/wigToBigWig>        specify path to wigToBigWig
  
  General options:
  -z --gz                               compress output text files
  -v --version                          print version and exit
  -h --help                             show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 File options

=over 4

=item --in E<lt>filenameE<gt>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. Genome 
coordinates are required. The first row should be column headers. Text 
files generated by other B<BioToolBox> scripts are acceptable. Files may 
be gzipped compressed.

=item --out E<lt>filenameE<gt>

Optionally specify the name of of the output file. The track name is 
used as default. The '.wig' extension is automatically added if required.

=item --noheader

The input file does not have column headers, often found with UCSC 
derived annotation data tables. 

=item --zero

Source data is in interbase coordinate (0-base) system. Shift the 
start position to base coordinate (1-base) system. Wig files are by 
definition 1-based. This is automatically handled for most input  
files. Default is false.

=back

=head2 Column indices

=over 4

=item --ask

Indicate that the program should interactively ask for column indices or
text strings for the GFF attributes, including coordinates, source, type, 
etc. It will present a list of the column names to choose from. Enter 
nothing for non-relevant columns or to accept default values.

=item --score E<lt>column_index or list of column indicesE<gt>

Indicate the column index of the dataset in the data table 
to be used for the score. If a GFF file is used as input, the score column is 
automatically selected. If not defined as an option, then the program will
interactively ask the user for the column index from a list of available
columns. More than one column may be specified, in which case the scores 
are combined using the method specified.

=item --chr E<lt>column_indexE<gt>

Optionally specify the column index of the chromosome or 
sequence identifier. This is required to generate the wig file. It may be 
identified automatically from the column header names.

=item --start E<lt>column_indexE<gt>

=item --begin E<lt>column_indexE<gt>

Optionally specify the column index of the start or chromosome 
position. This is required to generate the wig file. It may be 
identified automatically from the column header names.

=item --stop E<lt>column_indexE<gt>

=item --end E<lt>column_indexE<gt>

Optionally specify the column index of the stop or end 
position. It may be identified automatically from the column header names.

=item --attrib E<lt>attribute_nameE<gt>

Optionally provide the name of the attribute key which represents the score 
value to put into the wig file. Both GFF and VCF attributes are supported. 
GFF attributes are automatically taken from the attribute column (index 9).
For VCF columns, provide the index number of the sample column 
from which to take the value (usually 10 or higher) using the --index option. 
INFO field attributes can also be taken, if desired (use --index 8).

=back

=head2 Wig options

=over 4

=item --step [fixed | variable | bed]

The type of step progression for the wig file. Three wig formats are 
available: 
  - fixedStep: where data points are positioned at equal distances 
        along the chromosome
  - variableStep: where data points are variably positioned along 
        the chromosome. 
  - bed (bedGraph): where scores are associated with intervals 
        defined by start and stop coordinates.
The fixedStep wig file has one column of data (score), the variableStep 
wig file has two columns (position and score), and the bedGraph has 
four columns of data (chromosome, start, stop, score). If the option 
is not defined, then the format is automatically determined from the 
metadata of the file.

=item --bed

=item --bdg

Convenience option to specify a bedGraph file should be written. Same as 
specifying --step=bed.

=item --size E<lt>integerE<gt>

Optionally define the step size in bp for 'fixedStep' wig file. This 
value is automatically determined from the table's metadata, if available. 
If the C<--step> option is explicitly defined as 'fixed', then the step 
size may also be explicitly defined. If this value is not explicitly
defined or automatically determined, the variableStep format is used by
default.

=item --span E<lt>integerE<gt>

Optionally indicate the size of the region in bp to which the data value 
should be assigned. The same size is assigned to all data values in the 
wig file. This is useful, for example, with microarray data where all of 
the oligo probes are the same length and you wish to assign the value 
across the oligo rather than the midpoint. The default is inherently 1 bp. 

=back

=head2 Conversion options

=over 4

=item --fast

Disable checks for overlapping or duplicated intervals, unsorted data, valid
score values, and calculated midpoint positions. Requires setting the
chromosome, start, end (for bedGraph files only), and score column indices. 
B<WARNING:> Use only if you trust your input file format and content.

=item --name E<lt>textE<gt>

The name of the track defined in the wig file. The default is to use 
the name of the chosen score column, or, if the input file is a GFF file, 
the base name of the input file. 

=item --(no)track

Do (not) include the track line at the beginning of the wig file. Wig 
files normally require a track line, but if you will be converting to 
the binary bigwig format, the converter requires no track line. Why it 
can't simply ignore the line is beyond me. This option is automatically 
set to false when the C<--bigwig> option is enabled.

=item --mid

A boolean value to indicate whether the 
midpoint between the actual 'start' and 'stop' values
should be used. The default is to use only the 'start' position. 

=item --format E<lt>integerE<gt>

Indicate the number of decimal places the score value should
be formatted. The default is to not format the score value.

=item --method [mean | median | sum | max]

Define the method used to combine multiple data values at a single 
position. Wig files do not tolerate multiple identical positions. 
Default is mean.

=back

=head2 BigWig options

=over 4

=item --bigwig

=item --bw

Indicate that a binary BigWig file should be generated instead of 
a text wiggle file. 

=item --db E<lt>databaseE<gt>

Specify the name of a L<Bio::DB::SeqFeature::Store> annotation database 
or other indexed data file, e.g. Bam or bigWig file, from which chromosome 
length information may be obtained. It may be supplied from the input 
file metadata.

=item --chromof E<lt>filenameE<gt>

When converting to a BigWig file, provide a two-column tab-delimited 
text file containing the chromosome names and their lengths in bp. 
Alternatively, provide a name of a database, below.

=item --bwapp </path/to/wigToBigWig>

Specify the path to the UCSC wigToBigWig conversion utility. The default 
is to first check the BioToolBox configuration file C<biotoolbox.cfg> for 
the application path. Failing that, it will search the default 
environment path for the utility. If found, it will automatically 
execute the utility to convert the wig file.

=back 

=head2 General options

=over 4

=item --gz

A boolean value to indicate whether the output wiggle 
file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will convert any tab-delimited data text file into a wiggle 
formatted text file. This requires that the file contains not only the 
scores bu also chromosomal coordinates, i.e. chromosome, start, and 
(optionally) stop. The program should automatically detect these 
columns (if appropriately labeled) or they can be specified. An option 
exists to use the midpoint of a region, e.g. microarray probe.

The wig file format is specified by documentation supporting the UCSC 
Genome Browser and detailed here: http://genome.ucsc.edu/goldenPath/help/wiggle.html.
Three formats are supported, 'fixedStep', 'variableStep', and 'bedGraph'. 
The format may be requested or determined empirically from the input file 
metadata. Genomic bin files generated with BioToolBox scripts record 
the window and step values in the metadata, which are used to determine 
the span and step wig values, respectively. The variableStep format 
is otherwise generated by default. The span is, by default, 1 bp.

Wiggle files cannot tolerate multiple datapoints at the same identical 
position, e.g. multiple microarray probes matching a repetitive sequence. 
An option exists to mathematically combine these positions into one value.

Strand is not inherently supported in wig files. If you have stranded data, 
they should be split into separate files. The C<BioToolBox> script 
C<split_data_file.pl> can be used for this purpose.

A binary BigWig file may also be further generated from the  
text wiggle file. The binary format is preferential to the text version 
for a variety of reasons, including fast, random access and no loss in 
data value precision. More information can be found at this location:
http://genome.ucsc.edu/goldenPath/help/bigWig.html. Conversion requires 
BigWig file support, supplied by the external C<wigToBigWig> or 
C<bedGraphToBigWig> utility available from UCSC.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
