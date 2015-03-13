#!/usr/bin/perl

# documentation at end of file

# this script was copied directly from Bio::Graphics::Glyph::ideogram 

use strict;
my $VERSION = '1.0.0';
my %stains;
my %centros;
my %chrom_ends;

print "\n A script to generate a GBrowse compatible GFF3 cytoband file\n\n";

foreach (@ARGV) {
    if (/^(ftp|http|https):/) {
	$_ = "lynx --dump $_ |gunzip -c|";
    } elsif (/\.gz$/) {
	$_ = "gunzip -c $_ |";
    }
    print STDERR "Processing $_\n";
}

print "##gff-version 3\n";
while(<>)
{
    chomp;
    my($chr,$start,$stop,$band,$stain) = split /\t/;
    $start++;
    $chr = ucfirst($chr);
    if(!(exists($chrom_ends{$chr})) || $chrom_ends{$chr} < $stop)
    {
	$chrom_ends{$chr} = $stop;
    }
    my ($arm) = $band =~ /(p|q)\d+/;
    $stains{$stain} = 1;
    if ($stain eq 'acen')
    {
	$centros{$chr}->{$arm}->{start} = $stop;
	$centros{$chr}->{$arm}->{stop} = $start;
	next;
    }
    $chr =~ s/chr//i;
    print qq/$chr\tUCSC\tcytoband\t$start\t$stop\t.\t.\t.\tParent=$chr;Name=$chr;Alias=$chr$band;stain=$stain;\n/;
}

foreach my $chr(sort keys %chrom_ends)
{
    my $chr_orig = $chr;
    $chr =~ s/chr//i;
    print qq/$chr\tUCSC\tcentromere\t$centros{$chr_orig}->{p}->{stop}\t$centros{$chr_orig}->{q}->{start}\t.\t+\t.\tParent=$chr;Name=$chr\_cent\n/;
}


__END__

=head1 SYNOPSIS

ucsc_cytoband2gff3.pl <cytoBand_file>

This program will convert a UCSC cytoband file into a GFF3 file compatible
with the ideograpm glyph in GBrowse.

This was copied directly from Bio::Graphics::Glyph::ideogram.

=head1 AUTHOR

Gudmundur A. Thorisson E<lt>mummi@cshl.eduE<gt>

Copyright (c) 2001-2006 Cold Spring Harbor Laboratory

=head1 CONTRIBUTORS

Sheldon McKay E<lt>mckays@cshl.edu<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

