#!/usr/bin/perl
# This script will convert a GFF file to a SGR file for use in IGB
# It will use the midpoint between start and end positions (rounded down)
# 


print "\n This program will convert a *.gff file to *.sgr for use in IGB\n";

# Input
if ($ARGV[0]) {
	$filename = $ARGV[0];
} else {
	print " usage: gff2sgr.pl <filename.gff>\n";
	exit 1;
}
$filename =~ s/\.gff$//; # strip the gff extension
open FILEINPUT, "$filename\.gff" or die("Can't open $filename.gff\n");
my @input = <FILEINPUT>;
close FILEINPUT;



# GFF file consists of refseq, source, type, start, end, score, strand, phase, group
my @output;
foreach (@input) {
	if (/^#/) {next}; # skip comment lines
	my @data = split /\t/;
	my $mid = ($data[3] + $data[4]) / 2; # determine midpoint between start and end
	$mid = sprintf "%d", $mid;
	push @output, "$data[0]\t$mid\t$data[5]\n"; # refseq, midpoint, score
}

# Output
open OUTFILE, ">$filename\.sgr";
foreach (@output) {print OUTFILE $_}
close OUTFILE;

