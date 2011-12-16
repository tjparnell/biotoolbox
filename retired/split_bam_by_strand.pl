#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Pod::Usage;
my $VERSION = '1.4.2';


print "\n A program to split a bam file into separate files based on strand\n";


# open input file
my $infile = shift @ARGV or 
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );

my $in_bam = Bio::DB::Bam->open($infile, 'r') or 
	die "unable to open bam file '$infile'\n";



# open output files
my $outfile = $infile;
$outfile =~ s/\.bam$//;

my $f_file = $outfile . '.f.bam';
my $f_bam = Bio::DB::Bam->open($f_file, 'w') or 
	die "unable to open output bam file '$f_file' for writing\n";

my $r_file = $outfile . '.r.bam';
my $r_bam = Bio::DB::Bam->open($r_file, 'w') or 
	die "unable to open output bam file '$r_file' for writing\n";



# write headers
my $header = $in_bam->header();
$f_bam->header_write($header);
$r_bam->header_write($header);



# write the reads based on strand
print " sorting alignments...\n";
my $f_count = 0;
my $r_count = 0;
my $unmapped_count = 0;
while (my $alignment = $in_bam->read1() ) {
	
	# skip unmapped reads
	if ($alignment->unmapped) {
		$unmapped_count++;
		next;
	}
	
	# reverse strand
	if ($alignment->reversed) {
		$r_bam->write1($alignment);
		$r_count++;
	}
	
	# forward strand
	else {
		$f_bam->write1($alignment);
		$f_count++;
	}
}


# report results
print "  $f_count forward reads identified\n" if $f_count;
print "  $r_count reverse reads identified\n" if $r_count;
print "  $unmapped_count reads were unmapped and skipped\n" if $unmapped_count;



# close files
undef $in_bam; # is there a close function?????
undef $f_bam;
undef $r_bam;



# make new indices
	# the files should already be sorted, since we were using the low-level API
	# to read and write without affecting order
print " re-indexing...\n";
Bio::DB::Bam->index_build($f_file);
Bio::DB::Bam->index_build($r_file);

print " Wrote files '$f_file' and '$r_file'\n";
print " Done\n";



__END__

=head1 NAME

split_bam_by_strand.pl

A script to split reads in a bam file into two separate bam files based on strand.

=head1 SYNOPSIS

split_bam_by_strand.pl <file.bam>
  
=head1 DESCRIPTION

This program will read a bam file, identify aligned reads, and write the reads 
to one of two output bam files based on the strand to which the alignment 
matches. The output files are named the input base file name appended with 
either '.f' or '.r'. Unmatched files are discarded. The output files are 
then re-indexed for you. Header information, if present, is also retained.

This program is useful, for example, for stranded RNA-Seq, where you want to 
examine each strand separately.

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




