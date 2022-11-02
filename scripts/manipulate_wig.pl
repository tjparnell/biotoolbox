#!/usr/bin/perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case bundling);
use IO::Handle;
use Bio::ToolBox;
use Bio::ToolBox::big_helper qw(
	open_wig_to_bigwig_fh
	open_bigwig_to_wig_fh
	generate_chromosome_file
);
my $VERSION = 1.69;

### Quick help
unless (@ARGV) {    # when no command line options are present
					# print SYNOPSIS
	pod2usage(
		{
			'-verbose' => 0,
			'-exitval' => 1,
		}
	);
}

### Options
my $infile;
my $outfile;
my $skip;
my $apply;
my $doNull = 0;
my $deLogValue;
my $doAbsolute = 0;
my $multiplyValue;
my $addValue;
my $logValue;
my $places;
my $minValue;
my $maxValue;
my $noZeroes;
my $doStats;
my $bw2wig_app;
my $wig2bw_app;
my $chromofile;
my $database;
my $help;
my $print_version;

### Command line options
GetOptions(
	'i|input=s'    => \$infile,
	'o|output=s'   => \$outfile,
	'k|skip=s'     => \$skip,
	'y|apply=s'    => \$apply,
	'u|null!'      => \$doNull,
	'd|delog=i'    => \$deLogValue,
	'b|abs!'       => \$doAbsolute,
	'm|multiply=f' => \$multiplyValue,
	'a|add=f'      => \$addValue,
	'l|log=i'      => \$logValue,
	'p|place=i'    => \$places,
	'n|minimum=f'  => \$minValue,         #
	'x|maximum=f'  => \$maxValue,         #
	'z|zero'       => \$noZeroes,
	't|stats!'     => \$doStats,
	'bw2w=s'       => \$bw2wig_app,
	'w2bw=s'       => \$wig2bw_app,
	'chromo=s'     => \$chromofile,
	'db=s'         => \$database,
	'h|help'       => \$help,             # request help
	'v|version'    => \$print_version,    # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

### Print help if requested
if ($help) {

	# print entire POD
	pod2usage(
		{
			'-verbose' => 2,
			'-exitval' => 1,
		}
	);
}

### Print version
if ($print_version) {
	print " Biotoolbox script manipulate_datasets.pl, version $VERSION\n";
	printf " Biotoolbox package version %s\n", Bio::ToolBox->VERSION;
	exit;
}

### Checks
die "no input file provided!\n" unless $infile;
my $doMin = defined $minValue ? 1 : 0;
my $doMax = defined $maxValue ? 1 : 0;
if ($logValue) {
	$logValue = $logValue == 2 ? log(2) : $logValue == 10 ? log(10) : undef;
	die "bad log value!\n" unless defined $logValue;
}
if ( defined $places ) {
	$places = '%.' . $places . 'f';
}

# chromosome skipping regex
my ( $skip_regex, $apply_regex );
if ($skip) {
	$skip_regex = qr($skip);
}
if ($apply) {
	$apply_regex = qr($apply);
}

### Open file handles
# Input
my ( $infh, $outfh );
if ( $infile =~ /^stdin$/i ) {
	$infh = IO::Handle->new;
	$infh->fdopen( fileno(STDIN), 'r' );
}
elsif ( $infile =~ /(?:bw|bigwig)$/i and -e $infile ) {
	$infh = open_bigwig_to_wig_fh(
		bw        => $infile,
		bwapppath => $bw2wig_app,
	) or die "unable to open input bigWig file '$infile'!\n";
}
elsif ( -e $infile ) {
	$infh = Bio::ToolBox->read_file($infile)
		or die "can't open $infile! $!";
}
else {
	die "unrecognized $infile!";
}

# Output
if ( $outfile =~ /^stdout$/i ) {
	$outfh = IO::Handle->new;
	$outfh->fdopen( fileno(STDOUT), 'w' );
}
elsif ( $outfile =~ /(?:bw|bigwig)$/i ) {

	# check for chromosome file
	if ( $chromofile and -e $chromofile ) {

		# user provided a chromosome file
		$outfh = open_wig_to_bigwig_fh(
			bw        => $outfile,
			chromo    => $chromofile,
			bwapppath => $wig2bw_app,
		) or die "unable to open output bigWig file '$outfile'!\n";
	}
	elsif ( $infile =~ /(?:bw|bigwig)$/i or $database ) {

		# we can use the input bigWig as a database source if one isn't provided
		$database ||= $infile;
		$chromofile = generate_chromosome_file( $database, $skip )
			or die "unable to generate chromosome file from '$database'!\n";
		$outfh = open_wig_to_bigwig_fh(
			bw        => $outfile,
			chromo    => $chromofile,
			bwapppath => $wig2bw_app,
		) or die "unable to open output bigWig file '$outfile'!\n";
	}
	else {
		die "unable to open output bigWig file handle without chromosome information!\n";
	}

}
elsif ($outfile) {
	$outfh = Bio::ToolBox->write_file($outfile)
		or die "can't open $outfile! $!";
}

### stats hash
my $stats = {
	count      => 0,
	sumData    => 0,
	sumSquares => 0,
	minVal     => undef,
	maxVal     => undef,
};

### Walk through the file
my $count        = 0;
my $span         = 1;
my $chrom_skip   = 0;
my $chrom_ignore = 0;
my $wig_process_sub;
while ( my $line = $infh->getline ) {

	# look at the first characters to determine the type of line we have
	my $prefix = lc substr( $line, 0, 5 );
	if ( $prefix eq 'track' ) {

		# track line
		$outfh->print($line) if $outfh;
		next;
	}
	elsif ( $prefix eq 'brows' ) {

		# browser line
		$outfh->print($line) if $outfh;
		next;
	}
	elsif ( $prefix eq 'varia' or $prefix eq 'fixed' ) {

		# a step definition line
		if ( $line =~ /chrom=([\w\-\.]+)/ ) {

			# check the chromosome
			my $chrom = $1;
			if ( $skip_regex and $chrom =~ $skip_regex ) {
				$chrom_skip = 1;
			}
			else {
				$chrom_skip = 0;
			}
			if ( $apply and $chrom !~ $apply_regex ) {
				$chrom_ignore = 1;
			}
			else {
				$chrom_ignore = 0;
			}
		}
		if ( $line =~ /span=(\d+)/i ) {

			# capture span size if present
			$span = $1;
		}
		$outfh->print($line) if $outfh;
		next;
	}
	elsif ( substr( $prefix, 0, 1 ) eq '#' ) {

		# comment line
		$outfh->print($line) if $outfh;
		next;
	}

	# skipping current chromosome
	next if $chrom_skip;

	# ignoring current chromosome
	if ($chrom_ignore) {
		$outfh->print($line) if $outfh;
		next;
	}

	# determine format
	unless ( defined $wig_process_sub ) {
		my @data = split /\s+/, $line;
		my $statement;
		if ( scalar @data == 4 ) {
			$statement       = " processing bedGraph...\n";
			$wig_process_sub = \&process_bedGraph;
		}
		elsif ( scalar @data == 2 ) {
			$statement       = " processing variableStep...\n";
			$wig_process_sub = \&process_variableStep;
		}
		elsif ( scalar @data == 1 ) {
			$statement       = " processing fixedStep...\n";
			$wig_process_sub = \&process_fixedStep;
		}
		if ( $outfile =~ /stdout/i ) {
			print STDERR $statement;
		}
		else {
			print STDOUT $statement;
		}
	}

	# process
	chomp $line;
	&$wig_process_sub($line);
}

### close filehandles
$infh->close;
$outfh->close if $outfh;

# remove chromosome file if we generated it
unlink $chromofile
	if (    $outfile =~ /(?:bw|bigwig)$/i
		and $database
		and $chromofile =~ /^chr_sizes_\w{5}$/ );

### Print final messages
my $statMessage;
if ($doStats) {
	my $basecount = $stats->{count};
	my $min       = $stats->{minVal};
	my $max       = $stats->{maxVal};
	my $mean =
		$stats->{count} ? sprintf( "%.05f", $stats->{sumData} / $stats->{count} ) : 0;
	my $stddev = sprintf( "%.05f", sqrt( binVariance() ) );
	$statMessage = <<STATS;
basesCovered: $basecount
mean: $mean
min: $min
max: $max
std: $stddev
STATS
}

if ( $outfile =~ /stdout/i ) {
	print STDERR " converted $count lines, wrote file $outfile\n" if $outfile;
	print STDERR $statMessage                                     if $statMessage;
}
else {
	print STDOUT " converted $count lines, wrote file $outfile\n" if $outfile;
	print STDOUT $statMessage                                     if $statMessage;
}

#### Subroutines

sub process_bedGraph {
	my @data = split "\t", shift;
	return if ( $skip and $data[0] =~ $skip_regex );
	if ( $apply and $data[0] !~ $apply_regex ) {
		$outfh->printf( "%s\t%d\t%d\t%s\n", $data[0], $data[1], $data[2], $data[3] )
			if $outfh;
		return;
	}
	$data[3] = process_score( $data[3] );
	return if not defined $data[3];
	$outfh->printf( "%s\t%d\t%d\t%s\n", $data[0], $data[1], $data[2], $data[3] )
		if ( defined $data[3] and $outfh );
	$count++;
	if ( $doStats and defined $data[3] ) {
		my $length = $data[2] - $data[1];
		$stats->{count}      += $length;
		$stats->{sumData}    += ( $length * $data[3] );
		$stats->{sumSquares} += ( ( $data[3]**2 ) * $length );
		$stats->{minVal} = $data[3] if not defined $stats->{minVal};
		$stats->{minVal} = $data[3] if $data[3] < $stats->{minVal};
		$stats->{maxVal} = $data[3] if not defined $stats->{maxVal};
		$stats->{maxVal} = $data[3] if $data[3] > $stats->{maxVal};
	}
}

sub process_variableStep {
	my @data = split /\s+/, shift;    # could be either tab or space
	$data[1] = process_score( $data[1] );
	$outfh->printf( "%d %s\n", $data[0], $data[1] ) if ( defined $data[1] and $outfh );
	$count++;
	process_step_stats( $data[1] ) if $doStats;
}

sub process_fixedStep {
	my $score = shift;
	$score = process_score($score);
	$outfh->printf( "%s\n", $score ) if ( defined $score and $outfh );
	$count++;
	process_step_stats($score) if $doStats;
}

sub process_score {
	my $v = shift;    # score
	if ( $doNull and $v =~ /^(?:n.?[na])|(?:\-?inf)/i ) { $v = 0 }
	if ($deLogValue)                                    { $v = $deLogValue**$v }
	if ($doAbsolute)                                    { $v = abs($v) }
	if ($multiplyValue)                                 { $v *= $multiplyValue }
	if ($addValue)                                      { $v += $addValue }
	if ($logValue)                   { $v = $v == 0 ? 0 : log($v) / $logValue }
	if ( $doMin and $v < $minValue ) { $v = $minValue }
	if ( $doMax and $v > $maxValue ) { $v = $maxValue }
	if ($places)                     { $v = sprintf( $places, $v ) }
	return undef if ( $noZeroes and $v == 0 );
	return $v;
}

sub process_step_stats {
	return unless defined $_[0];
	for ( 1 .. $span ) {
		$stats->{count}      += 1;
		$stats->{sumData}    += $_[0];
		$stats->{sumSquares} += $_[0]**2;
		$stats->{minVal} = $_[0] if not defined $stats->{minVal};
		$stats->{maxVal} = $_[0] if not defined $stats->{maxVal};
		$stats->{minVal} = $_[0] if $_[0] < $stats->{minVal};
		$stats->{maxVal} = $_[0] if $_[0] > $stats->{maxVal};
	}
}

sub binVariance {
	return 0 unless $stats->{count};
	my $var = $stats->{sumSquares} - $stats->{sumData}**2 / $stats->{count};
	if ( $stats->{count} > 1 ) {
		$var /= $stats->{count} - 1;
	}
	return 0 if $var < 0;
	return $var;
}

__END__

=head1 NAME

manipulate_wig.pl

A progam to manipulate wiggle files.

=head1 SYNOPSIS

manipulate_wig.pl [options] -i <file1.wig> -o <file1.out.wig>

  File Options: 
  -i --in <file>            Input file. Accepts 'stdin'.
  -o --out <file>           Output file. Accepts 'stdout'.
  
  Selection functions:
  -k --skip <regex>         Skip lines where chromosomes match regex 
  -y --apply <regex>        Only apply manipulations to matching chromosomes 
  
  Manipulation functions (in order of execution):
  -u --null                 Convert null, NA, N/A, NaN, inf values to 0
  -d --delog [2|10]         Delog values of given base 
  -b --abs                  Convert to the absolute value 
  -m --mult <float>         Multiply score by the given value
  -a --add <float>          Add the given value to the score
  -l --log [2|10]           Convert to log2 or log10. 
  -n --min <float>          Set the minimum score
  -x --max <float>          Set the maximum score
  -p --place <int>          Format score to decimal positions
  -z --zero                 Discard lines with zero values

  BigWig support:
  --chromo <file>           Chromosome sizes file for writing bigWig
  --db <file>               Indexed file to obtain chromosome info
  --bw2w <path>             Path to UCSC bigWigToWig utility
  --w2bw <path>             Path to UCSC wigToBigWig utility
  
  General functions:
  -t --stats                Calculate statistics 
  -v --version              print version and exit
  -h --help                 show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 File options

=over 4

=item --in E<lt>fileE<gt>

Specify the input wig file. All three formats, variableStep, fixedStep, and 
bedGraph, are supported. Files may be gzipped. BigWig files are supported, 
so long as the UCSC bigWigToWig utility is available. Alternatively, the input 
may be read from standard input by specifying 'stdin' as the file name. 

=item --out E<lt>fileE<gt>

Specify the output wig file. The output format will be the same format as the
input. The file may be gzipped by appending F<.gz> to the name. BigWig files are
supported, so long as the UCSC wigToBigWig utility is available and a chromosome
file is provided. Alternatively, the output may be sent to standard output by
specifying 'stdout' as the file name. 

=back

=head2 Selection functions

=over 4

=item --skip E<lt>regexE<gt>

Selectively skip (discard) lines corresponding to certain chromosomes that 
match the provided regular expression. For example, skip the 
mitochondrial and random contigs, use "chrM|chrUn|random".

=item --apply E<lt>regexE<gt>

Selectively apply manipulation functions to certain chromosomes that match 
provided regular expression, leaving remaining lines untouched. For example, 
to apply a normalization to the X chromosome, use 'chrX'.

=back

=head2 Manipulation functions

=over 4

=item --null

Convert lines with a score of C<null>, C<NA>, C<N/A>, C<NaN>, or C<inf> to 
a value of 0. 

=item --delog [2|10]

Convert lines from log space in the indicated base.

=item --abs

Convert line scores to absolute values.

=item --mult E<lt>floatE<gt>

Multiply line scores by the indicated value.

=item --add E<lt>floatE<gt>

Add the indicated value to each line score.

=item --log [2|10]
 
Convert the line score to a log equivalent in the indicated base space.

=item --min E<lt>floatE<gt>

Set the minimum floor score. Any score below the indicated value 
will be set to the indicated value.

=item --max E<lt>floatE<gt>

Set the maximum ceiling score. Any score above the indicated value 
will be set to the indicated value.

=item --place E<lt>integerE<gt>

Format the score value to the indicated number of decimal positions.

=item --zero 

Discard lines with a score value of zero.

=back

=head2 BigWig support

=over 4

=item chromo E<lt>fileE<gt>

When writing to a bigWig output file, provide a chromosome sizes text 
file for use with the F<wigToBigWig> utility. Alternatively, use a 
database file, below.

=item db E<lt>fileE<gt>

When writing to a bigWig output file, provide an indexed database file, 
such as another bigWig file, Bam, indexed Fasta, etc, for automatically 
generating a chromosome sizes text file to use with the F<wigToBigWig> 
utility. If a bigWig input file was specified, it will be conveniently 
substituted as a database. B<Note> that the C<--skip> option will be 
applied to the generated chromosome file.

=item bw2w E<lt>pathE<gt>

If the UCSC F<bigWigToWig> utility is not in your environment C<PATH>, 
provide the path with this option.

=item w2bw E<lt>pathE<gt>

If the UCSC F<wigToBigWig> utility is not in your environment C<PATH>, 
provide the path with this option.

=back

=head2 General functions

=over 4

=item --stats

Calculate the statistics across the genome at base pair resolution. 
Statistics are calculated after any processing indicated above are 
performed. Only applied chromosomes are calculated. Results are 
printed to STDOUT, or STDERR if standard out is used as an output.

=item --version

Print the version number of the program and exit.

=item --help

Display the POD documentation using perldoc. 

=back

=head1 DESCRIPTION

A program to manipulate the score value of wig files. This will process all 
forms of text based wig files, including fixedStep, variableStep, and bedGraph. 
Files may be gzip compressed. BigWig files are also transparently supported as 
both input and output, provided that the appropriate UCSC utility files are 
available.

B<NOTE:> More than one option may be specified! The options above are the order 
in which the score is manipulated. If they are not in the order you want, you 
may have to pipe to sequential instances. Use 'stdin' and 'stdout' for filenames.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

