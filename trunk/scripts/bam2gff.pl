#!/usr/bin/perl

# a script to convert bam paired_reads to a gff file

use strict;
use Getopt::Long;
use Pod::Usage;
eval {use Bio::DB::Sam};
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	format_with_commas
);
use tim_file_helper qw(
	open_to_write_fh
);

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
	$paired_end,
	$type, 
	$source,
	$minsize, 
	$maxsize, 
	$AT_ends,
	$gz,
	$help, 
);
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
	'pe'         => \$paired_end, # bam is paired end reads
	'type=s'     => \$type, # the GFF type or method
	'source=s'   => \$source, # GFF source field
	'min=i'      => \$minsize, # the minimum size for paired-read segments
	'max=i'      => \$maxsize, # the maximum size for paired-read segments
	'at'         => \$AT_ends, # discard non-AT ends
	'gz!'        => \$gz, # gzip the output file
	'help'       => \$help, # request help
);

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}



### Check for required values and set defaults
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		die "  An input BAM file must be specified!\n";
	}
}

unless ($type) {
	# derive the type from the input file name
	$type = $infile;
	$type =~ s/\.bam$//;
	$type =~ s/\.sorted$//;
	if ($paired_end) {
		$type .= '_paired_reads';
	}
	else {
		$type .= '_reads';
	}
}

unless ($outfile) {
	$outfile = $type;
}

unless ($source) {
	$source = 'Illumina';
}

if ($paired_end) {
	unless (defined $minsize) {
		$minsize = 100; # set default to 100 bp
		print " Using default minimum size of 100 bp\n";
	}
	
	unless (defined $maxsize) {
		$maxsize = 200; # set default to 200 bp
		print " Using default maximum size of 200 bp\n";
	} 
}

unless (defined $gz) {
	# default is no compression
	$gz = 0;
}




### Load the SAM file
print " Opening bam file....\n";
my $sam = Bio::DB::Sam->new( 
	-bam        => $infile,
	-autoindex  => 1,
) or die " unable to open bam file '$infile'!\n";



### Open output gff file to write
unless ($outfile =~ /\.gff3?$/) {
	$outfile .= '.gff3';
}
my $gff_out = open_to_write_fh($outfile, $gz) or 
	die " unable to open output file '$outfile'!\n";
print {$gff_out} "##gff_version 3\n";
print {$gff_out} "# Program $0\n";
print {$gff_out} "# Converted from source file $infile\n";



### Perform the conversion
my $total_count = 0;
if ($paired_end) {
	convert_paired_end_alignments();
}
else {
	convert_single_end_alignments();
}


### Done
print " wrote file '$outfile'\n";



###################### Subroutines ###########################

sub convert_paired_end_alignments {
	
	# write additional header information
	print {$gff_out} "# Mininum size $minsize\n";
	print {$gff_out} "# Maximum size $maxsize\n";
	print {$gff_out} "# Discarded AT ends\n" if $AT_ends;
	
	
	# Initialize counts
	my $undefined_count  = 0;
	my $improper_count   = 0;
	my $diffchromo_count = 0;
	my $non_AT_end_count = 0;
	my $toosmall_count   = 0;
	my $toobig_count     = 0;
	my $just_right_count = 0;
	
	
	# Start the conversion
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as a numeric
		# target identifier
		# we can easily convert this to an actual sequence name
		# we will force the conversion to go one chromosome at a time
		
		# sequence name
		my $seq_id = $sam->target_name($tid);
		
		print " Converting reads on $seq_id...\n";
		my $iterator = $sam->features(
				'-type'     => 'read_pair',
				'-iterator' => 1,
				'-seq_id'   => $seq_id,
		);
		
		CONVERSION_LOOP:
		while (my $pair = $iterator->next_seq() ) {
			$total_count++;
			
			# get individual reads
			my ($left, $right) = $pair->get_SeqFeatures;
			unless (defined $left and defined $right) {
				# one or both alignments are not defined
				# does this mean one of the alignments is unaligned?
				$undefined_count++;
				next;
			}
			unless ($left->proper_pair) {
				# the pair of reads are not properly mapped
				# is this redundant?
				$improper_count++;
				next;
			}
			my $refseq = $left->seq_id;
			if ($refseq ne $right->seq_id) {
				# not on same chromosomes!???? is this part of the improper count?
				$diffchromo_count++;
				next;
			}
			
			# check AT ends
			if ($AT_ends) {
				# both ends must be either A or T to continue, otherwise report
				my $leftseq = $left->query->dna; 
				my $rightseq = $right->query->dna; 
					# this should be reverse-complemented
				unless ($leftseq =~ /^[aAtT]/ and $rightseq =~ /[aAtT]$/) {
					# one of the ends is not correct
					$non_AT_end_count++;
					next;
				}
			}
			
			# check length first
			my $length = $pair->length;
			if ($length < $minsize) {
				$toosmall_count++;
				next CONVERSION_LOOP;
			}
			elsif ($length > $maxsize) {
				$toobig_count++;
				next CONVERSION_LOOP;
			}
			
			# get coordinates
			my $start = $left->start;
			my $end = $right->end;
			
			# generate group
			my $name = $left->query->name;
			$name =~ s/[\:\-\#\(\)\;\,\.]/_/g; # substitute punctuation with _ 
			my $id = $type . '.' . $total_count;
			my $group = "ID=$id; Name=$name";
			
			# print the GFF feature
			print {$gff_out} join("\t", (
				$refseq,
				$source, 
				$type,
				$start,
				$end,
				'.', # score, not used
				'.', # strand, not used
				'.', # phase, not used
				$group
			) ), "\n";
			
			# success
			$just_right_count++;
		}
		
		undef $iterator; # is this really necessary? help with memory or something
	}
	
	
	
	### Finish up 
	# close file
	undef $sam;
	$gff_out->close; 
	
	# print summaries
	print "\n There were " . format_with_commas($total_count) . 
		" total alignment pairs read\n";
	
	print "   " . format_with_commas($undefined_count) . " (". 
		percent_pc($undefined_count) . 
		") pairs had an undefined (unmapped?) alignment\n" 
		if $undefined_count > 0;
	print "   " . format_with_commas($improper_count) . " (". 
		percent_pc($improper_count) . 
		") pairs were improper\n" if $improper_count > 0;
	print "   " . format_with_commas($diffchromo_count) . " (". 
		percent_pc($diffchromo_count) . 
		") pairs were split between different chromosomes\n" 
		if $diffchromo_count > 0;
	
	print "   " . format_with_commas($non_AT_end_count) . " (". 
		percent_pc($non_AT_end_count) . ") pairs had non-AT ends\n" 
		if $non_AT_end_count > 0;
	
	print "   " . format_with_commas($toosmall_count) . " (". 
		percent_pc($toosmall_count) . 
		") pairs were below the lowest size $minsize bp\n";
	print "   " . format_with_commas($toobig_count) . " (". 
		percent_pc($toobig_count) . 
		") pairs were above the highest size $maxsize bp\n";
	print "   " . format_with_commas($just_right_count) . " (". 
		percent_pc($just_right_count) . ") pairs were just right\n";

}



sub convert_single_end_alignments {
	
	# Initialize counts
	my $unmapped_count = 0;
	my $mapped_count   = 0;
	
	
	# Start the conversion
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as a numeric
		# target identifier
		# we can easily convert this to an actual sequence name
		# we will force the conversion to go one chromosome at a time
		
		# sequence name
		my $seq_id = $sam->target_name($tid);
		
		print " Converting reads on $seq_id...\n";
		my $iterator = $sam->features(
				'-type'     => 'match',
				'-iterator' => 1,
				'-seq_id'   => $seq_id,
		);
		
		CONVERSION_LOOP:
		while (my $alignment = $iterator->next_seq() ) {
			$total_count++;
			
			# skip unmapped reads
			if ($alignment->unmapped) {
				$unmapped_count++;
				next;
			}
			
			# get strand
			my $strand;
			if ($alignment->strand > 0) {
				# forward strand
				# I know, this seems counterintuitive at first, but empirical 
				# analysis says otherwise
				$strand = '+';
			}
			else {
				# reverse strand
				$strand = '-';
			}
			
			# generate group
			my $name = $alignment->query->name;
			$name =~ s/[:\-#\(\);,\.]/_/g; # substitute punctuation with _ 
			my $id = $type . '.' . $total_count;
			my $group = "ID=$id; Name=$name";
			
			# print the GFF feature
			print {$gff_out} join("\t", (
				$seq_id,
				$source, 
				$type,
				$alignment->start,
				$alignment->end,
				$alignment->qual, # alignment quaility for the score
				$strand, 
				'.', # phase, not used
				$group
			) ), "\n";
			
			# success
			$mapped_count++;
		}
		
		undef $iterator; # is this really necessary? help with memory or something
	}
	
	
	
	### Finish up 
	# close file
	undef $sam;
	$gff_out->close; 
	
	# print summaries
	print "\n There were " . format_with_commas($total_count) . 
		" total alignment pairs read\n";
	
	print "   " . format_with_commas($unmapped_count) . " (". 
		percent_pc($unmapped_count) . ") reads were unmapped\n" 
		if $unmapped_count > 0;
	print "   " . format_with_commas($mapped_count) . " (". 
		percent_pc($mapped_count) . ") pairs were properly mapped\n" 
		if $mapped_count > 0;

}





sub percent_pc {
	# for calculating the percent of a count out of the total paired read number
	my $count = shift;
	return sprintf "%.1f%%", ($count / $total_count) * 100;
}





__END__

=head1 NAME

bam_paired2gff.pl

A script to convert alignments from a BAM file to a GFF v3 file.

=head1 SYNOPSIS

bam_2gff.pl [--options] <file.bam>
  
  Options:
  --in <filename.bam>
  --out <filename>
  --pe
  --type <gff_type>
  --source <text>
  --min <integer>
  --max <integer>
  --at
  --gz
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a binary BAM file as described for Samtools. 
Use samtools to convert text sam files to binary bam files..

=item --out <filename>

Optionally specify the name of of the output file. The default is to use the
GFF type as the filename. The extension '.gff3' will be added as necessary.

=item --pe

Indicate that the bam file is comprised of paired-end alignments. With 
paired-end alignments, the resulting GFF feature represents the insert 
fragment whose ends were sequenced. This fragment is subject to size 
cutoffs using the min and max options. The default is to treat the 
bam alignments as single-end alignments.

=item --type <gff_type>

Specify a text string to be used as the GFF type. If not specified, it is 
derived from the base file name, minus '.bam' and '.sorted' if present, 
plus '_paired_read' (for paired-end) or '_read' (for single-end).

=item --source <text>

Optionally the text to be used in the GFF source field. The default is 
'Illumina'. 

=item --min <integer>

Specify the minimum size of the paired-end fragment insert that will be 
included in the output GFF file. The default value is 100 bp.

=item --max <integer>

Specify the maximum size of the paired-end fragment insert that will be 
included in the output GFF file. The default value is 200 bp.

=item --at

Boolean option to indicate that only fragments whose paired sequence reads 
end in [AT] nucleotides should be included in the output GFF file. 
Micrococcal nuclease (MNase) cuts (almost?) exclusively at AT dinucleotides; 
this option ensures that the fragment is (more) likely derived from a MNase 
cut. This option is not available with single-end reads.

=item --gz

Compress the output gff file with gzip.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will read a SAM/BAM file containing alignments and
convert them to a GFF version 3 file. By default it generates the 
GFF features from the individual alignments. With paired-end alignments, 
it will optionally generate GFF features representing the insert fragment 
whose ends were sequenced. Only those fragments whose sizes fall within 
the defined minimum and maximum sizes (defaults of 100 and 200 bp, 
respectively) are reported. Additionally, paired-end fragments may be 
checked for the presence of A/T ends, as produced by digestion with 
Micrococcal Nuclease.

The input file is a SAM/BAM file as described by the Samtools project 
(http://samtools.sourceforge.net). The file must be sorted and indexed 
prior to conversion (indexing should be automatic).

The output file is a GFF version 3 file. Each feature represents an 
aligned sequenced pair of acceptable size. 

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112




