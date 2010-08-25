#!/usr/bin/perl

# a script to combine pileup SNPs from multiple sequencing reads
# and identify unique and common SNPs

use strict;
use IO::File;
use IO::Zlib;

print "\nA script to identify unique and common SNPs between multiple strains\n\n";



# Check input files
unless (@ARGV) {
	print " usage: merge_SNPs.pl <file1> <file2> ...\n";
	exit;
}
if (scalar @ARGV == 1) {
	# only one file provided, this isn't right
	# check for help
	if ($ARGV[0] eq '-h' or $ARGV[0] eq '--help') {
		# help is requested
		exec "perldoc $0";
	}
	else {
		die " More than one SNP file is required\n";
	}
}



# Generate data structures
my @names; # to store the names of the files
my %snps; # structure will be: chr -> pos -> snp -> name = line
my %output; # structure will be: name(s) -> [lines]



# Collect the snp data from all the files
foreach my $file (@ARGV) {
	
	# generate base name by stripping unnecessary stuff
	my $name = $file;
	$name =~ s/\.gz$//;
	$name =~ s/\.txt$//;
	$name =~ s/\.?pileup//;
	$name =~ s/\.?filtered//;
	push @names, $name;
	
	# open file
	my $fh;
	if ($file =~ /\.gz$/) {
		$fh = new IO::Zlib($file, "rb");
	}
	else {
		$fh = new IO::File($file, "r");
	}
	
	# walk through file
	my $count = 0;
	while (my $line = $fh->getline) {
		chomp $line;
		
		# collect chromosome and position, first two elements
		my ($chr, $pos, $snp) = (split /\t/, $line)[0,1,3];
		
		# store the snps
		$snps{$chr}{$pos}{$snp}{$name} = $line;
		
		$count++;
	}
	
	# finish
	$fh->close;
	print " Loaded '$file' with $count SNPs...\n";
}



# Sort through all the SNPs
print " Sorting SNPs....\n";
foreach my $chr (sort { $a cmp $b } keys %snps) {
	# sort chromosomes asciibetically
	
	foreach my $pos (sort {$a <=> $b} keys %{ $snps{$chr} } ) {
		# sort by increasing position
		
		foreach my $snp (sort {$a cmp $b} keys %{ $snps{$chr}{$pos} } ) {
			# sort by the type of snp
			# in all liklihood this will only be one, but just in case
			
			# pull out the sequence names that have this snp
			my @names = sort {$a cmp $b} keys %{ $snps{$chr}{$pos}{$snp} };
			
			# determine the appropriate place to put the line
			my $name = join('_', @names);
			
			# create output key if necessary
			unless (exists $output{$name}) {
				$output{$name} = [];
			}
			
			# place in appropriate key
			# for common SNPs (found in more than one source) we will only be
			# storing a representative line
			# It would be way too messy to try and combine them
			push @{ $output{$name} }, $snps{$chr}{$pos}{$snp}{$names[0]};
		}
	}
}



# Output the separate SNPs
foreach my $name (sort {$a cmp $b} keys %output) {
	# report count
	print " There were ", scalar( @{ $output{$name} } ), " SNPs for '$name'\n";
	
	# print
	my $file = $name . '_SNPs.txt';
	my $fh = new IO::File($file, 'w');
	foreach ( @{ $output{$name} } ) {
		print {$fh} "$_\n";
	}
	$fh->close;
}







__END__

=head1 NAME

merge_SNPs.pl

=head1 SYNOPSIS
 
merge_SNPs.pl <file1> <file2> ...
  
=head1 DESCRIPTION

This program will identify common and unique SNPs between two or more 
sequenced strains. This is useful, for example, in identifying SNPs that 
may be background polymorphisms common to all the strains, versus unique 
SNPs that may be responsible for a mutant phenotype. 

Each strain should have a separate SNP file, generated using the
'varFilter' function of the 'samtools.pl' script, part of the Samtools
distribution L<http://samtools.sourceforge.net>. That script will generate
a list of the sequence variations that differ from the reference genome. 
No verification of the source file format is performed. The files may be 
gzipped.

In this script, the SNPs are sorted into groups based on their occurance in
one or more strains. These groups are then written out to new separate
files, with the file names being a concatenation of the representing
original filenames. The file format is preserved, as each SNP line in the
file is a representative example from one of the original strains.


=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112










