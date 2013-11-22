#!/usr/bin/perl

# a script to split and strip tags from a single fastq file to several files
# 

use strict;
use File::Basename qw(fileparse);
use Bio::SeqIO;
use Bio::Seq::Quality;

print "\n A script to split a fastq file based on barcode tags\n\n";





### Help and arguments

# print help if necessary
unless (@ARGV) {
	die "   usage: split_by_tags.pl <sequence.fastq> <tag1> <tag2> ...\n";
}
if (scalar @ARGV == 1) {
	# only one argument provided
	# check for help
	if ($ARGV[0] eq '-h' or $ARGV[0] eq '--help') {
		# help is requested
		exec "perldoc $0";
	}
	else {
		die " Must provide both sequence file name and list of tags\n";
	}
}

# get arguments
my $file = shift @ARGV;
my @tags = @ARGV;






### Open input sequence file

# parse the file name, we want the basename for the output file names
my ($basename, $path, $extension) = fileparse($file, qw(
		.txt
		.txt.gz
		.fa
		.fa.gz
		.fq
		.fq.gz
		.fastq
		.fastq.gz
) ); # this list hopefully covers most (all?) the possibilities

# check for gzip status
my $gz = 0;
if ($extension =~ /\.gz$/) {
	# filter all the file IO through gzip compression
	# we're doing it this way because BioPerl modules complain about IO::Zlib
	# file objects
	$gz = 1; 
	$file = "gunzip -c $file |";
}


# open the input file
my $seq_stream = Bio::SeqIO->new(
	-file              => $file,
	-format            => 'fastq',
) or die " unable to open input file '$file'\n";





### Prepare output files
my %out;
foreach my $tag (@tags) {
	# determine the tag length, usually 3
	my $len = length $tag;
	
	# open the output filehandle
	my $outfile = $path . $basename . '_' . $tag . $extension;
	my $outseq;
	
	
	if ($gz) {
		$outseq = Bio::SeqIO->new(
			-file             => " | gzip -c > $outfile",
			-format           => 'fastq',
			-quality_header   => 1,
		) or die " unable to open gzipped output file '$outfile'\n";
	} 
	else {
		$outseq = Bio::SeqIO->new(
			-file             => "> $outfile",
			-format           => 'fastq',
			-quality_header   => 1
		) or die " unable to open output file '$outfile'\n";
	}
	
	# store the filehandle
	$out{ $tag } = [ $len, $outfile, $outseq, 0 ]; # length, name, object, count
}





### Process through the input file
my $none = 0; # number of seqs that fail to match a tag
while (my $seq = $seq_stream->next_seq() ) {
	
	# get the raw sequence
	my $sequence = $seq->seq;
	
	# work through the list of tags
	my $check = 0;
	foreach (@tags) {
		if ($sequence =~ /^$_/i) {
			# the sequence matches the tag
			
			# get the subsequence minus the tag
			my $sub_start = $out{$_}->[0] + 1; # tag length + 1
			my $length = $seq->length;
			my $new_seq = $seq->subseq($sub_start, $length);
			my $new_qual = $seq->subqual($sub_start, $length);
			
			# update the sequence
			$seq->seq($new_seq);
			$seq->qual($new_qual);
			
			# print to the output file stream object
			$out{ $_ }->[2]->write_seq($seq) or 
				die " unable to write to file '$out{$_}->[1]'!!!\n";
			
			# update count
			$out{ $_ }->[3] += 1;
			$check = 1;
			last;
		}
	}
	
	# update the no count if necessary
	unless ($check) {
		$none++;
	}
}




### Print the results
foreach (@tags) {
	print "   split out ", $out{$_}->[3], " sequences for tag $_ and wrote " .
		"file '", $out{$_}->[1], "'\n";
}
if ($none) {
	print "   $none sequences did not match a tag\n";
}
print " finished\n";





__END__

=head1 NAME

split_by_tags.pl

=head1 SYNOPSIS
 
split_by_tags.pl <sequence.fastq> <tag1> <tag2> ...
  
=head1 DESCRIPTION

This program will split a single fastq file into multiple fastq files based on 
a barcode tag at the 5' end of the sequence. Usually the barcode is two 
nucleotides, plus a T for linker ligation, for example GGT or AAT. After 
identifying the tag in the sequence, it is stripped from the sequeence. A 
separate file is then written for each tag, with the output file name being the 
basename appended with the tag sequence. Input files may be gzipped, and 
output files will preserve compression status. Input files should have a 
recognizable extension, e.g. .txt, .fa, .fq, .fastq.

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









