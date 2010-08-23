#!/usr/bin/perl
# This script will convert a bed file to a GFF file for use in the genome browser
# The bed file is described at http://genome.ucsc.edu/FAQ/FAQformat#format1
# the contents include the following in this order:
#   chromosome, start (0-based), stop (exclusive), name, score, strand, etc.


print "\n\t This program will convert *.bed files to a *.gff v2 file\n";

# Input
if ($ARGV[0]) {
	$filename = $ARGV[0];
} else {
	print "\n\t Please type in the bed file name  ";
	chomp($filename = <STDIN>);
}
unless ($filename =~ /\.bed$/) {die "please enter a *.bed file\n"}

# Ask for specific GFF information
print "What is the name for this data?  ";
my $type = <STDIN>;
chomp $type;
if ($type =~ /\s/) {die("Can't have whitespace in $type\n") }
print "Enter new source name [default: data]  ";
my $source = <STDIN>;
chomp $source;
if ($source eq '') {$source = 'data'}

# Do the conversion
open INFILE, $filename;
my @output;
while (my $line = <INFILE>) {
	if ($line =~ /^track/i) {next} # skip the track definition line
	chomp $line;
	my @data = split /\t/, $line;
	my $refseq = $data[0];
	my $start = $data[1] + 1; # need to shift from 0-based indexing
	my $end = $data[2] - 1; # need to shift from exclusive number to an inclusive number
	my $score;
	# score is optional in the bed format
	if ($data[4]) { $score = $data[4] } else { $score = '.' }
	my $strand;
	# strand is optional, but if present is either + or -
	if ($data[5]) { $strand = $data[5] } else { $strand = '.' }
	my $phase = '.'; 
	my ($name, $group);
	# name is optional
	if ($data[3]) { 
		$name = $data[3];
		$group = "$type \"$name\"";
	} else {
		$group = "Experiment \"$type\"";
	}
	push @output, "$refseq\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\t$group\n";
}
close INFILE;

# Output
$filename =~ s/\.bed$//;
open OUTFILE, ">$filename.gff";
print OUTFILE "##gff-version 2\n";
print OUTFILE "# generated using program $0\n";
print OUTFILE "# from source file $filename.bed\n";
print OUTFILE @output;
close OUTFILE;
