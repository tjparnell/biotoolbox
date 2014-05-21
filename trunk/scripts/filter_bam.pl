#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::data_helper qw(parse_list format_with_commas);
my $bam_ok;
eval {
	# check for Bam support
	require Bio::DB::Sam;
	$bam_ok = 1;
};
my $VERSION = '1.19';

print "\n A script to filter a Bam file for specific criteria\n\n";

### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Command line options
my (
	$infile, 
	$outfile, 
	$write_true,
	$write_false,
	$do_align,
	$do_mismatch,
	$do_gap,
	$do_indel,
	$do_mate_proper,
	$do_mate_strand,
	$do_mate_seqid,
	$do_score,
	$do_length,
	$do_index,
	$help,
	$print_version,
);
my @sequences;
my @attributes;
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
	'pass!'      => \$write_true, # write those that pass criteria
	'fail!'      => \$write_false, # write those that fail criteria
	'align!'     => \$do_align, # check alignment
	'mismatch!'  => \$do_mismatch, # check for mismatches
	'gap!'       => \$do_gap, # check for gaps
	'indel!'     => \$do_indel, # check for indels
	'mproper!'   => \$do_mate_proper, # check for proper pair
	'mseqid!'    => \$do_mate_seqid, # check mate seq_id
	'mstrand!'   => \$do_mate_strand, # check mate strand
	'score=i'    => \$do_score, # check alignment score
	'length=s'   => \$do_length, # check length
	'seq=s'      => \@sequences, # check specific sequence
	'attrib=s'   => \@attributes, # check attributes
	'index!'     => \$do_index, # re-index the output files
	'help'       => \$help, # request help
	'version'    => \$print_version, # print the version
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
	print " Biotoolbox script filter_bam.pl, version $VERSION\n\n";
	exit;
}



### Check for required values and set defaults
unless ($bam_ok) {
	die "Module Bio::DB::Sam must be installed to run this script.\n";
}
my ($true_file, $false_file);
my @lengths;
my @filters;
check_defaults();
my $start_time = time;



### Open the Bam files
my ($in_bam, $true_bam, $false_bam) = open_bam_files();



### Filter the Bam file
filter_bam();

finish_bam_files() if $do_index;

printf " Finished in %.1f in minutes\n", (time - $start_time) / 60;







########### Subroutines ################################

sub check_defaults {
	
	# input file
	unless ($infile) {
		if (@ARGV) {
			$infile = shift @ARGV;
		}
		else {
			die " An input BAM file must be specified!\n";
		}
	}
	unless ($infile =~ /\.bam$/i) {
		die " Input file must be a .bam file!\n";
	}

	# assign output names
	unless ($outfile) {
		$outfile = $infile;
		$outfile =~ s/\.bam$//;
		$outfile =~ s/\.sorted//;
		$outfile .= '.filter';
	}
	$write_true = 1 unless defined $write_true;
	if ($write_true and $write_false) {
		$true_file = "$outfile\.pass";
		$true_file .= '.bam' unless $true_file =~ /\.bam/;
		$false_file = "$outfile\.fail";
		$false_file .= '.bam' unless $false_file =~ /\.bam/;
	}
	elsif ($write_true and not $write_false) {
		$true_file = "$outfile";
		$true_file .= '.bam' unless $true_file =~ /\.bam/;
	}
	elsif (not $write_true and $write_false) {
		$false_file = "$outfile";
		$false_file .= '.bam' unless $false_file =~ /\.bam/;
	}
	
	# check alignment, default is true
	$do_align = 1 unless defined $do_align;
	
	# check user provided lists
	if ($do_length) {
		@lengths = parse_list($do_length);
	}
	if (@sequences) {
		for my $i (0 .. $#sequences) {
			my ($pos, $nuc) = split /:/, $sequences[$i];
			$nuc = lc $nuc;
			die "unrecognized nucleotide in requested sequence filter\n" unless 
				$nuc =~ /^[acgt]+$/;
			$sequences[$i] = [$pos, $nuc]
		}
	}
	if (@attributes) {
		for my $i (0 .. $#attributes) {
			my ($key, $value) = split /:/, $attributes[$i];
			my @values;
			@values = split /,/, $value if defined $value;
			$attributes[$i] = [$key, [@values]]
		}
	}
	
	# set the filters
	push @filters, \&filter_for_alignment if $do_align;
	push @filters, \&filter_for_mismatch if $do_mismatch;
	push @filters, \&filter_for_gap if $do_gap;
	push @filters, \&filter_for_indel if $do_indel;
	push @filters, \&filter_for_proper_mate if $do_mate_proper;
	push @filters, \&filter_for_mate_seqid if $do_mate_seqid;
	push @filters, \&filter_for_mate_strand if $do_mate_strand;
	push @filters, \&filter_for_score if $do_score;
	push @filters, \&filter_for_sequence if @sequences;
	push @filters, \&filter_for_attribute if @attributes;
	push @filters, \&filter_for_length if $do_length;
}

sub open_bam_files {
	
	# input bam file
	my $in = Bio::DB::Bam->open($infile) or die " Cannot open input Bam file!\n";
	my $header = $in->header; # must always get before reading alignments

	# output bam files
	my ($true, $false);
	if ($write_true) {
		$true = Bio::DB::Bam->open($true_file, 'w') or 
			die "Cannot open output file $true_file!\n";
		$true->header_write( $header );
	}
	if ($write_false) {
		$false = Bio::DB::Bam->open($false_file, 'w') or 
			die "Cannot open output file $false_file!\n";
		$false->header_write( $header );
	}
	return ($in, $true, $false);
}

sub filter_bam {
	my @counts = (0,0,0); # total, true, false counts
	
	# walk through each alignment in the bam file
	while (my $a = $in_bam->read1) {
		callback($a, \@counts);
	}
	
	
	# report results
	printf " %s (%.1f%%) alignments passed criteria %s", 
		format_with_commas($counts[1]), 
		($counts[1] / $counts[0]) * 100,
		$true_file ? " and were written to $true_file\n" : "\n";
	printf " %s (%.1f%%) alignments failed criteria %s", 
		format_with_commas($counts[2]), 
		($counts[2] / $counts[0]) * 100,
		$false_file ? " and were written to $false_file\n" : "\n";
}

sub callback {
	my ($a, $counts) = @_;
	
	# check the alignment
	my $check = 1; # default is true
	foreach my $filter (@filters) {
		unless ( &{$filter}($a) ) {
			$check = 0;
			last; # no need to go on
		}
	}
	
	# update counts and write as necessary
	$counts->[0]++; # total count
	if ($check) {
		$counts->[1]++;
		$true_bam->write1($a) if ($true_file);
	}
	else {
		$counts->[2]++;
		$false_bam->write1($a) if ($false_file);
	}
}

sub finish_bam_files {
	# close the output bam files and index them
	# they will be sorted as necessary
	if ($write_true) {
		undef $true_bam;
		Bio::DB::Bam->reindex($true_file);
	}
	if ($write_false) {
		undef $false_bam;
		Bio::DB::Bam->reindex($false_file);
	}
}

sub filter_for_alignment {
	# return true if mapped
	return $_[0]->unmapped ? 0 : 1;
}

sub filter_for_mismatch {
	# return true if a mismatch is present
	# this one is a little tricky
	my $a = $_[0];
	my %keys = map {$_ => 1} $a->aux_keys;
	if (exists $keys{'NM'}) {
		# we have an NM tag, use that value
		return $a->aux_get('NM');
	}
	elsif (exists $keys{'MD'}) {
		# we have an MD tag
		# no mismatch means just a number
		return $a->aux_get('MD') !~ /^\d+$/;
	}
	else {
		# the aligner did not include the NM attribute,
		# attempt to derive from CIGAR string
		return $a->cigar_str =~ /X/;
	}
}

sub filter_for_gap {
	# return true if gap is present
	return $_[0]->cigar_str =~ /[N]/;
}

sub filter_for_indel {
	# return true if insertion or deletion is present
	return $_[0]->cigar_str =~ /[ID]/;
}

sub filter_for_proper_mate {
	# return true if part of proper pair
	return $_[0]->proper_pair;
}

sub filter_for_mate_align {
	# return true if mate is mapped
	return $_[0]->munmapped ? 0 : 1;
}

sub filter_for_mate_seqid {
	# return true if mate is on the same reference sequence
	return $_[0]->tid == $_[0]->mtid;
}

sub filter_for_mate_strand {
	# return true if mate is on the opposite strand
	return $_[0]->reversed != $_[0]->mreversed;
}

sub filter_for_score {
	# return true if score is >= score
	return $_[0]->qual >= $do_score;
}

sub filter_for_sequence {
	# return true if query sequence nucleotide equal to that requested at 
	# the specified position
	my $seq = lc $_[0]->qseq;
	my $check = 1; # default is true
	foreach my $s (@sequences) {
		# each $s is an array of ($position, nucleotide)
		my $nuc = $s->[1];
		unless (substr($seq, $s->[0] - 1, 1) =~ /[$nuc]/ ) {
			$check = 0;
			last;
		}
	}
	return $check;
}

sub filter_for_length {
	# return true if sequence length is equal to one that is requested
	my $len = length $_[0]->qseq;
	foreach my $l (@lengths) {
		return 1 if $l == $len;
	}
	return 0;
}

sub filter_for_attribute {
	# return true if value for attribute key equals the requested value
	# there is some tricky logic here to make sure everything fits
	my $a = $_[0];
	my %keys = map {$_ => 1} $a->aux_keys;
	my $check = 0; # default is false
	ATTRIB_LOOP: foreach my $attrib (@attributes) {
		# each $attrib is an array of (key, [@values])
		# first makes sure the key exists
		if (exists $keys{ $attrib->[0] } ) {
			# key exists, so far so good
			
			# then look for each value one at a time
			if (@{ $attrib->[1] }) {
				VALUE_LOOP: foreach my $value ( @{ $attrib->[1] } ) {
					if ($a->aux_get($attrib->[0]) eq $value) {
						# we have a positive match, check
						$check = 1;
						next ATTRIB_LOOP; 
					}
				}
				$check = 0; # if we arrive here, none of the values matched
				last ATTRIB_LOOP;
			}
			else {
				# no values to match, but since the key is present, pass the check
				$check = 1;
			}
		}
		else {
			# key is not present, fail the check
			$check = 0;
			last ATTRIB_LOOP;
		}
	}
	return $check;
}



__END__

=head1 NAME

filter_bam.pl

A script to filter a Bam file for specific criteria.

=head1 SYNOPSIS

filter_bam.pl <file.bam>
  
  Options:
  --in <file.bam>
  --out <filename>
  --(no)pass
  --(no)fail
  --(no)align
  --mismatch
  --gap
  --indel
  --mproper
  --mseqid
  --mstrand
  --score <integer>
  --length <integer>
  --seq <pos:[ATCG]>
  --attrib <key:value>
  --index
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <file.bam>

Specify the file name of a binary Bam file as described for 
Samtools. It does not need to be sorted or indexed.

=item --out <filename>

Optionally specify the base name of the output file. The default is 
to use input base name, appended with '.filter'. If both pass and 
fail files are written, then they are appended with '.pass' and 
'.fail', respectively. 

=item --pass

=item --nopass

Indicate whether (or not) alignments which pass the test criteria 
should be written to an output Bam file. The default is true.

=item --fail

=item --nofail

Indicate whether (or not) alignments which fail the test criteria 
should be written to an output Bam file. The default is false.

=item --align

=item --noalign

Indicate whether (or not) aligned reads should pass. The default is true.

=item --mismatch

Indicate that only alignments with a mismatch should pass. A 
mismatch is indicated by either the C<NM> or C<MD> attributes of 
the alignment, or by the presence of X (mismatch) operations  
in the CIGAR string. Gaps, clipped, or padded sequences are not 
counted.

=item --gap

Indicate that only alignments with a gap should pass. Gaps are 
determined by the presence of N (skipped) operations in the 
CIGAR string.

=item --indel

Indicate that only alignments with either an insertion or deletion 
should pass. Indels are determined by the presence of I (insertion) 
or D (deletion) operations in the CIGAR string.

=item --mproper

Indicate that only alignments that are part of a proper pair 
should pass. Proper pairs are Forward-Reverse alignments on 
the same reference, and do not include Forward-Forward, 
Reverse-Reverse, Reverse-Forward, or separate reference 
sequence alignments.

=item --mseqid

Indicate that only paired alignments that are on the same 
reference sequence should pass. 

=item --mstrand

Indicate that only paired alignments that align to different 
strands should pass, i.e. a Forward-Reverse or Reverse-Forward.

=item --score <integer>

Indicate that only alignments which have a quality score equal 
or greater than that indicated shall pass. The mapping quality 
score is a posterior probability that the alignment was mapped
incorrectly, and reported as a -10Log10(P) value, rounded to the 
nearest integer (range 0..255). Higher numbers are more stringent. 

=item --length <integer>

Indicate that only alignments whose query sequence equals the 
indicated length shall pass. Provide a comma-delimited list and/or 
range of lengths. Note that only the query sequence is checked, 
not the length of the alignment. Multiple lengths are treated 
as a logical OR operation.

=item --seq <pos:[ATCG]>

Indicate that only alignments that have a specific nucleotide at 
a specific position in the query sequence shall pass. Provide a 
position:nucleotide pair, where position is a 1-based integer and 
the nucleotide is one or more of A,C,G, or T. Providing two or 
more nucleotides per position is treated as a logical OR operation. 
Multiple sequence positions may be tested by issuing multiple 
command line options, in which case they are combined in a logical 
AND operation.

=item --attrib <key>

=item --attrib <key:value>

Indicate that only alignments that contain a specific optional 
attribute shall pass. One or more values may also be provided for 
the key, in which case only those alignments which match one of 
the key values shall pass. The values may be provided as a comma 
delimited list separated from the key by a colon. Attribute keys 
are typically two letter codes; see the SAM specification at 
L<http://samtools.sourceforge.net/SAM1.pdf> for a list of standard 
attributes. Two or more key values are combined in a logical OR 
operation. Two or more attribute keys may be tested by specifying 
multiple --attrib command line options; in this case, they are  
combined in a logical AND operation.

=item --index

Optionally re-index the output bam file(s) when finished. If 
necessary, the bam file is sorted by coordinate first. Default is 
false.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will filter the alignments in a Bam file according to a 
series of one or more boolean tests. Alignments which pass all the 
tests are written to an output Bam file. Alignments which do not pass 
one or more filters may be optionally written to a second Bam file.

There are a number of tests that may be applied to the alignments, 
controlled by command line arguments. Please note carefully how the 
test is performed and whether your desired outcome should be the 
pass or fail outcome. When multiple tests are indicated, they are 
combined using a logical AND operation. 

The input and output files are BAM files as described by the Samtools 
project (http://samtools.sourceforge.net). 

=head1 EXAMPLES

Here are a few examples of how to use filters.

=over

=item Alignments that may indicate a SNP

SNPs could be either a mismatch, insertion, or deletion
 
 filter_bam.pl --mismatch --indel --in file.bam

=item RNASeq alignments that could span an intron
 
 filter_bam.pl --gap --in file.bam

=item MNase digested DNA

Chromatin may be digested using MNase, which cuts blunt ends between 
[AT][AT] dinucleotides. To increase the likelihood that sequences were 
derived from MNase digestion, filter for an [AT] nucleotide at the 
first position.
 
 filter_bam.pl --seq 1:AT --in file.bam

=item Alignments indicating chromosomal rearrangement

Paired-end sequencing of genomic DNA where two ends map to separate 
chromosomes or not in a proper forward-reverse arrangement may 
suggest a chromosomal rearrangement. In this case, we want those 
alignments that fail the test.
 
 filter_bam.pl --nopass --fail --mproper --out non_properly_paired --in file.bam
 
 filter_bam.pl --nopass --fail --mseqid --out translocations --in file.bam

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
