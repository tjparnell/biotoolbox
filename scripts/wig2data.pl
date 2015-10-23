#!/usr/bin/perl

# documentation at end of file 

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::Data::Stream;
my $VERSION = '1.33';

print "\n This program will convert wiggle files to a tabbed text file\n\n";

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
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the solexa data file
	'out=s'     => \$outfile, # name of output file 
	'gz!'       => \$gz, # compress output
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
	print " Biotoolbox script wig2data.pl, version $VERSION\n\n";
	exit;
}




### Check for requirements
unless ($infile) {
	$infile = shift @ARGV or
		die "  OOPS! No source data file specified! \n use $0 --help\n";
}
unless (defined $gz) {
	$gz = 0;
}
unless ($infile =~ /\.wig(?:\.gz)?$/i) {
	die " only .wig files are supported!\n";
}


# Open the file 
print " Converting file '$infile'...\n";
my $in_fh = Bio::ToolBox::Data::Stream->open_to_read_fh($infile) or 
	die " unable to open file!\n";


# Open output Stream
unless ($outfile) {
	$outfile = $infile;
	$outfile =~ s/\.wig/.txt/;
}
my $Output = Bio::ToolBox::Data::Stream->new(
	out     => $outfile,
	columns => [ qw(Chromosome Start Stop Score) ],
	gz      => $gz,
) or die "unable to write $outfile";


# reusuable variables
my (
	$refseq,
	$fixstart,
	$step,
	$span,
); 
my $count = 0;


# main loop
while (my $line = $in_fh->getline) {
	my $start; # specific line variables
	my $stop;
	my $score;
	
	# The wiggle file can have 3 different formats: BED format, variable step, 
	# and fixed step. We need to determine whether each line is a definition
	# line or a data line, based on the line's contents and/or number of 
	# elements. The definition lines will fill the reusable variables above
	# and help in filling out the specific variables.
	
	## check the line's contents
	$line =~ s/[\r\n]+$//; # strip all line endings
	my @data = split /\s+/, $line;
	
	# a track line
	if ($data[0] =~ /track/i) {
		# not much useable information in here for us
		$Output->add_comment($line);
		next;
	}
	
	# a variable step definition line
	elsif ($data[0] =~ /^variablestep$/i) { 
		foreach (@data) {
			if (/chrom=(\w+)/) {$refseq = $1}
			if (/span=(\w+)/) {$span = $1}
		}
		next;
	} 
	
	# a fixed step definition line
	elsif ($data[0] =~ /^fixedstep$/i) { 
		foreach (@data) {
			if (/chrom=(\w+)/) {$refseq = $1}
			if (/span=(\w+)/) {$span = $1}
			if (/start=(\w+)/) {$fixstart = $1}
			if (/step=(\w+)/) {$step = $1}
		}
		next;
	} 
	
	# a BED data line
	elsif (scalar @data == 4) {
		die "input appears to be a bedGraph wig file. Change the extension to .bdg and use as is.\n"
	} 
	
	# a variable step data line
	elsif (scalar @data == 2) { 
		unless ($refseq) { 
			die "Bad formatting! variable step data but chromosome not defined!\n";
		}
		$start = $data[0];
		if ($span) {
			$stop = $start + $span - 1;
		} 
		else {
			$stop = $start;
		}
		$score = $data[1];
	} 
	
	# a fixed step data line
	elsif (scalar @data == 1) { 
		unless ($refseq) { 
			die "Bad formatting! fixed step data but chromosome not defined!\n";
		}
		unless ($fixstart) { 
			die "Bad formatting! fixed step data but start not defined!\n";
		}
		unless ($step) { 
			die "Bad formatting! fixed step data but step not defined!\n";
		}
		$start = $fixstart;
		$fixstart += $step; # prepare for next round
		if ($span) {
			$stop = $start + $span - 1;
		} 
		else {
			$stop = $start;
		}
		$score = $data[0];
	}
	
	# write the line
	$Output->add_row( [
		$refseq,
		$start,
		$stop,
		$score
	] );
	
	$count++;
}

### Finish
$in_fh->close;
$Output->close_fh;
print " wrote $count lines to file '$outfile'\n";


__END__

=head1 NAME

wig2data.pl

A script to convert a text wiggle file to a tab-delimited text file. 

=head1 SYNOPSIS

wig2data.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --out <filename> 
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a wig file. The file may be compressed with gzip.

=item --out <filename>

Specify the output filename. By default it uses the input base name.

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will convert a fixedStep or variableStep wiggle data file into 
a tabbed delimited text data file. The data file will have four columns: 
chromosome, start, stop, and score. 

More information about wig files may be obtained from here:
http://genome.ucsc.edu/goldenPath/help/wiggle.html

The start position will be the coordinate listed in the wig file. If a span
value is indicated in the wiggle file, the GFF stop will equal start plus
span; otherwise stop will equal the start value. The score value may be
formatted to the indicated number of decimal places.
 
=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
