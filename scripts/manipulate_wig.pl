#!/usr/bin/perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case bundling);
use IO::Handle;
use Bio::ToolBox::Data::file;
my $VERSION = 1.60;

### Quick help
unless (@ARGV) { # when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
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
my $help;
my $print_version;



### Command line options
GetOptions( 
	'i|input=s'       => \$infile,
	'o|output=s'      => \$outfile,
	'k|skip=s'        => \$skip,
	'y|apply=s'       => \$apply,
	'u|null!'         => \$doNull,
	'd|delog=i'       => \$deLogValue,
	'b|abs!'          => \$doAbsolute,
	'm|multiply=f'    => \$multiplyValue,
	'a|add=f'         => \$addValue,
	'l|log=i'         => \$logValue,
	'p|place=i'       => \$places,
	'n|minimum=f'     => \$minValue, # 
	'x|maximum=f'     => \$maxValue, # 
	'z|zero'          => \$noZeroes,
	't|stats!'        => \$doStats,
	'h|help'          => \$help, # request help
	'v|version'       => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";


### Print help if requested
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

### Print version
if ($print_version) {
	print " Biotoolbox script manipulate_datasets.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
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
if (defined $places) {
	$places = '%.' . $places . 'f';
}

# chromosome skipping regex
my ($skip_regex, $apply_regex);
if ($skip) {
	$skip_regex = qr($skip);
}
if ($apply) {
	$apply_regex = qr($apply);
}





### Open file handles
my ($infh, $outfh);
if ($infile =~ /^stdin$/i) {
	$infh = IO::Handle->new;
	$infh->fdopen(fileno(STDIN), 'r');
}
elsif (-e $infile) {
	$infh = Bio::ToolBox::Data::file->open_to_read_fh($infile) or 
		die "can't open $infile! $!";
}
else {
	die "unrecognized $infile!";
}

if ($outfile =~ /^stdout$/i) {
	$outfh = IO::Handle->new;
	$outfh->fdopen(fileno(STDOUT), 'w');
}
elsif ($outfile) {
	$outfh = Bio::ToolBox::Data::file->open_to_write_fh($outfile) or 
		die "can't open $outfile! $!";
}




### stats hash
my $stats = {
	count         => 0,
	sumData       => 0,
	sumSquares    => 0,
	minVal        => undef,
	maxVal        => undef,
};


### Walk through the file
my $count = 0;
my $span = 1;
my $chrom_skip = 0;
my $chrom_ignore = 0;
my $wig_process_sub;
while (my $line = $infh->getline) {
	# look at the first characters to determine the type of line we have
	my $prefix = lc substr($line,0,5);
	if ($prefix eq 'track') {
		# track line
		$outfh->print($line) if $outfh;
		next;
	}
	elsif ($prefix eq 'brows') {
		# browser line
		$outfh->print($line) if $outfh;
		next;
	}
	elsif ($prefix eq 'varia' or $prefix eq 'fixed') {
		# a step definition line
		if ($line =~ /chrom=([\w\-\.]+)/) {
			# check the chromosome
			my $chrom = $1;
			if ($skip_regex and $chrom =~ $skip_regex) {
				print STDERR "skipping chromosome $chrom\n";
				$chrom_skip = 1;
			}
			else {
				$chrom_skip = 0;
			}
			if ($apply and $chrom !~ $apply_regex) {
				print STDERR "ignoring chromosome $chrom\n";
				$chrom_ignore = 1;
			}
			else {
				$chrom_ignore = 0;
			}
		}
		if ($line =~ /span=(\d+)/i) {
			# capture span size if present
			$span = $1;
		}
		$outfh->print($line) if $outfh;
		next;
	} 
	elsif (substr($prefix,0,1) eq '#') {
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
	unless (defined $wig_process_sub) {
		my @data = split /\s+/, $line;
		if (scalar @data == 4) {
			$wig_process_sub = \&process_bedGraph;
		}
		elsif (scalar @data == 2) {
			$wig_process_sub = \&process_variableStep;
		}
		elsif (scalar @data == 1) {
			$wig_process_sub = \&process_fixedStep;
		}
	}
	
	# process
	chomp $line;
	&$wig_process_sub($line);
}



### close filehandles
$infh->close;
$outfh->close if $outfh;



### Print final messages
my $statMessage;
if ($doStats) {
	my $basecount = $stats->{count};
	my $min   = $stats->{minVal};
	my $max   = $stats->{maxVal};
	my $mean  = $stats->{count} ? sprintf("%.05f", $stats->{sumData} / $stats->{count}) : 0;
	my $stddev = sprintf("%.05f", sqrt(binVariance()) );
	$statMessage = <<STATS; 
basesCovered: $basecount
mean: $mean
min: $min
max: $max
std: $stddev
STATS
}

if ($outfile =~ /stdout/i) {
	print STDERR " converted $count lines, wrote file $outfile\n" if $outfile;
	print STDERR $statMessage if $statMessage;
}
else {
	print STDOUT " converted $count lines, wrote file $outfile\n" if $outfile;
	print STDOUT $statMessage if $statMessage;
}




#### Subroutines

sub process_bedGraph {
	my @data = split "\t", shift;
	return if ($skip and $data[0] =~ $skip_regex);
	if ($apply and $data[0] !~ $apply_regex) {
		$outfh->printf("%s\t%d\t%d\t%s\n", $data[0], $data[1], $data[2], $data[3]) 
			if $outfh;
		return;
	}
	$data[3] = process_score($data[3]);
	return if not defined $data[3];
	$outfh->printf("%s\t%d\t%d\t%s\n", $data[0], $data[1], $data[2], $data[3]) 
		if (defined $data[3] and $outfh);
	$count++;
	if ($doStats and defined $data[3]) {
		my $length = $data[2] - $data[1];
		$stats->{count} += $length;
		$stats->{sumData} += ($length * $data[3]);
		$stats->{sumSquares} += ( ($data[3] ** 2) * $length );
		$stats->{minVal} = $data[3] if not defined $stats->{minVal};
		$stats->{minVal} = $data[3] if $data[3] < $stats->{minVal};
		$stats->{maxVal} = $data[3] if not defined $stats->{maxVal};
		$stats->{maxVal} = $data[3] if $data[3] > $stats->{maxVal};
	}
}

sub process_variableStep {
	my @data = split /\s+/, shift; # could be either tab or space
	$data[1] = process_score($data[1]);
	$outfh->printf("%d %s\n", $data[0], $data[1]) if (defined $data[1] and $outfh);
	$count++;
	process_step_stats($data[1]) if $doStats;
}

sub process_fixedStep {
	my $score = shift;
	$score = process_score($score);
	$outfh->printf("%s\n",$score) if (defined $score and $outfh);
	$count++;
	process_step_stats($score) if $doStats;
}

sub process_score {
	my $v = shift; # score
	if ($doNull and $v =~ /^(?:n.?[na])|(?:\-?inf)/i) {$v = 0}
	if ($deLogValue) {$v = $deLogValue ** $v}
	if ($doAbsolute) {$v = abs($v)}
	if ($multiplyValue) {$v *= $multiplyValue}
	if ($addValue) {$v += $addValue}
	if ($logValue) {$v = $v == 0 ? 0 : log($v) / $logValue}
	if ($doMin and $v < $minValue) {$v = $minValue}
	if ($doMax and $v > $maxValue) {$v = $maxValue}
	if ($places) {$v = sprintf($places, $v)};
	return undef if ($noZeroes and $v == 0);
	return $v;
}

sub process_step_stats {
	return unless defined $_[0];
	for (1 .. $span) {
		$stats->{count} += 1;
		$stats->{sumData} += $_[0];
		$stats->{sumSquares} += $_[0] ** 2;
		$stats->{minVal} = $_[0] if not defined $stats->{minVal};
		$stats->{maxVal} = $_[0] if not defined $stats->{maxVal};
		$stats->{minVal} = $_[0] if $_[0] < $stats->{minVal};
		$stats->{maxVal} = $_[0] if $_[0] > $stats->{maxVal};
	}
}

sub binVariance {
    return 0 unless $stats->{count};
    my $var = $stats->{sumSquares} - $stats->{sumData}**2/$stats->{count};
    if ($stats->{count} > 1) {
	$var /= $stats->{count}-1;
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
  
  Manipulation functions:
  -u --null                 Convert null, NA, N/A, NaN, inf values to 0
  -d --delog [2|10]         Delog values of given base 
  -a --abs                  Convert to the absolute value 
  -m --mult <float>         Multiply score by the given value
  -a --add <float>          Add the given value to the score
  -l --log [2|10]           Convert to log2 or log10. 
  -n --min <float>          Set the minimum score
  -x --max <float>          Set the maximum score
  -p --place <int>          Format score to decimal positions
  -z --zero                 Discard lines with zero values

  General functions:
  -t --stats                Calculate statistics 
  -v --version              print version and exit
  -h --help                 show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <file>

Specify the input wig file. All three formats, variableStep, fixedStep, and 
bedGraph, are supported. Files may be gzipped. Alternatively, the input 
may be read from standard input by specifying 'stdin' as the file name. For 
bigWig files, you may use a pipe as shown below

    bigWigToWig file1.bw stdout | manipulate_wig.pl -i stdin 

=item --out <file>

Specify the output wig file. The output format will be the same format as the 
input. The file may be gzipped by appending F<.gz> to the name. Alternatively, 
the output may be sent to standard output by specifying 'stdout' as the file 
name. For writing to bigWig files, use a pipe as show below

    manipulate_wig.pl -o stdout | wigToBigWig stdin chroms.txt file.bw

=item --skip <regex>

Selectively skip (discard) lines corresponding to certain chromosomes that 
match the provided regular expression. For example, skip the 
mitochondrial and random contigs, use "chrM|chrUn|random".

=item --apply <regex>

Selectively apply manipulation functions to certain chromosomes that match 
provided regular expression, leaving remaining lines untouched. For example, 
to apply a normalization to the X chromosome, use 'chrX'.

=item --null

Convert lines with a score of C<null>, C<NA>, C<N/A>, C<NaN>, or C<inf> to 
a value of 0. 

=item --delog [2|10]

Convert lines from log space in the indicated base.

=item --abs

Convert line scores to absolute values.

=item --mult <float>

Multiply line scores by the indicated value.

=item --add <float>

Add the indicated value to each line score.

=item --log [2|10]
 
Convert the line score to a log equivalent in the indicated base space.

=item --min <float>

Set the minimum floor score. Any score below the indicated value 
will be set to the indicated value.

=item --max <float>

Set the maximum ceiling score. Any score above the indicated value 
will be set to the indicated value.

=item --place <integer>

Format the score value to the indicated number of decimal positions.

=item --zero 

Discard lines with a score value of zero.

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
Files can be gzipped. 

NOTE: More than one option may be specified! The options above are the order 
in which the score is manipulated. If they are not in the order you want, you 
may have to pipe to sequential instances. Use 'stdin' and 'stdout' for filenames.
Use an equal sign to define options with negative values, e.g. --mult=-1

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

