#!/usr/bin/perl

# a script to selectively write out paired-end alignments of a given size

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Sam;


print "\n A script to split a paired-end bam file by insert sizes\n\n";

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
	$minsize, 
	$maxsize, 
	$AT_ends,
	$help, 
);
my @size_list;
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
	'min=i'      => \$minsize, # the minimum cutoff size for paired-read segments
	'max=i'      => \$maxsize, # the maximum cutoff size for paired-read segments
	'size=s'     => \@size_list, # a list of sizes to select
	'at'         => \$AT_ends, # discard non-AT ends
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
# input file
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		die "  An input BAM file must be specified!\n";
	}
}

# sizes to select
my @sizes; # an array of arrays of the sizes to include
if (@size_list) {
	foreach (@size_list) {
		my @s = split /-/;
		die " Improperly formatted size range [$_]!\n" if @s != 2;
		push @sizes, [ @s ];
	}
}
else {
	unless (defined $minsize) {
		$minsize = 100; # set default to 100 bp
		print " Using default minimum size of 100 bp\n";
	}
	unless (defined $maxsize) {
		$maxsize = 200; # set default to 200 bp
		print " Using default maximum size of 200 bp\n";
	} 
	push @sizes, [ ($minsize, $maxsize) ];
}

# identify the lowest and highest sizes
my $lowest = $sizes[0][0];
my $highest = $sizes[0][1];
foreach (@sizes) {
	$lowest = $_->[0] if $_->[0] < $lowest;
	$highest = $_->[1] if $_->[1] > $highest;
}

unless ($outfile) {
	$outfile = $infile;
	$outfile =~ s/\.bam$//;
	$outfile =~ s/\.sorted//;
}

	

### Initialize counts
my $pair_count = 0;
my $undefined_count = 0;
my $toosmall_count = 0;
my $toobig_count = 0;
my $improper_count = 0;
my $diffchromo_count = 0;
my $non_AT_end_count = 0;
my $just_right_count = 0;





### Open BAM files
print " Opening bam files....\n";
# we are opening the bam files using the low level bam API

# input file
my $in_sam = Bio::DB::Sam->new( 
	-bam        => $infile,
	-autoindex  => 1,
) or die " unable to open input bam file '$infile'!\n";
	# we are opening the input bam file using the high level sam API
print "   input file '$infile'\n";

# input header
my $header = $in_sam->header();
unless ($header) {
	die "no header in input bam file!!!\n";
}

# output files
foreach (@sizes) {
	# we will open a separate file for each size range
	
	# generate specific file name
	my $bam_file = $outfile . '.' . $_->[0] . '_' . $_->[1] . '.bam';
	
	# open bam file
	my $bam = Bio::DB::Bam->open($bam_file, 'w') 
		or die "unable to open output bam file '$outfile' for writing!\n";
	print "   output file '$bam_file'\n";
	
	# write headers
	$bam->header_write($header);
	
	# store
	$_->[2] = 0; # a count for the number of pairs written to this file
	$_->[3] = $bam_file;
	$_->[4] = $bam;
}
	




### Start conversion

# loop through the chromosomes
for my $tid (0 .. $in_sam->n_targets - 1) {
	# each chromosome is internally represented in the bam file as a numeric
	# target identifier
	# we can easily convert this to an actual sequence name
	# we will force the conversion to go one chromosome at a time
	
	# sequence name
	my $seq_id = $in_sam->target_name($tid);
	
	print " Converting reads on $seq_id...\n";
	my $iterator = $in_sam->features(
			'-type'     => 'read_pair',
			'-iterator' => 1,
			'-seq_id'   => $seq_id,
	);
	
	while (my $pair = $iterator->next_seq() ) {
		$pair_count++;
		
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
		
		# check length against lowest and highest values
		my $length = $pair->length;
		if ($length < $lowest) {
			$toosmall_count++;
			next;
		}
		elsif ($length > $highest) {
			$toobig_count++;
			next;
		}
		
		# fine the right size range and write to the appropriate file
		foreach my $size (@sizes) {
			
			# check sizes
			if (
				$length >= $size->[0] and
				$length <= $size->[1]
			) {
				# within size range so write the alignments
					# these objects are B:D:AlignWrapper objects, not B:D:Alignment 
					# objects, therefore I can't simply just write them out to the 
					# the low level bam object
					# there is no official way to convert from AlignWrapper to 
					# Alignment (even though AlignWrapper wraps around Alignment).
					
					# therefore, I have to manually extract the Alignment object.
					# Using Dumper, identified the Alignment object as the value 
					# under the object's key 'align'
					# and then write that to the output bam file
				foreach ($left, $right) {
					my $align = $_->{align};
					$size->[4]->write1($align); # write to the bam file object
				}
				$size->[2] += 1; # count
			}
		}
		
		# success
		$just_right_count++;
	}	
	
	undef $iterator; # is this really necessary? help with memory or something
}


### Finish up 
# close files
undef $in_sam;
foreach (@sizes) {
	pop @{$_}; # undefine the bam object
}
sleep 3; # does this give the OS a chance to close and write the files?

# resort and index the bam files
foreach my $size (@sizes) {
	
	# sort
	my $new_file = $size->[3];
	$new_file =~ s/\.bam/.sorted/;
	print " re-sorting $size->[3]... ";
	Bio::DB::Bam->sort_core(0, $size->[3], $new_file);
	
	# make new indices
	$new_file .= '.bam'; # sorting would've automatically added the extension
	if (-e $new_file) {
		unlink $size->[3]; # remove the old unsorted output file
		print " re-indexing...\n";
		Bio::DB::Bam->index_build($new_file);
	}
	$size->[3] = $new_file;
}


# print summaries
print "\n There were " . add_commas($pair_count) . 
	" total alignment pairs read\n";

print "   " . add_commas($undefined_count) . " (". percent_pc($undefined_count) . 
	") pairs had an undefined (unmapped?) alignment\n" if $undefined_count > 0;
print "   " . add_commas($improper_count) . " (". percent_pc($improper_count) . 
	") pairs were improper\n" if $improper_count > 0;
print "   " . add_commas($diffchromo_count) . " (". percent_pc($diffchromo_count) . 
	") pairs were split between different chromosomes\n" if $diffchromo_count > 0;

print "   " . add_commas($non_AT_end_count) . " (". percent_pc($non_AT_end_count) . 
	") pairs had non-AT ends\n" if $non_AT_end_count > 0;

print "   " . add_commas($toosmall_count) . " (". percent_pc($toosmall_count) . 
	") pairs were below the lowest size $lowest bp\n";
print "   " . add_commas($toobig_count) . " (". percent_pc($toobig_count) . 
	") pairs were above the highest size $highest bp\n";
print "   " . add_commas($just_right_count) . " (". percent_pc($just_right_count) . 
	") pairs were just right\n";

foreach (@sizes) {
	print " " . add_commas($_->[2]) . " (" . percent_pc($_->[2]) . 
		") pairs were written to file '$_->[3]'\n";
}
print "\n";

### Subroutines

sub percent_pc {
	# for calculating the percent of pair_count (total)
	my $count = shift;
	return sprintf "%.2f%%", ($count / $pair_count) * 100;
}

sub add_commas {
	# for formatting a number with commas
	my $number = shift;
	my @digits = split //, $number;
	my @formatted;
	while (@digits) {
		if (@digits > 3) {
			unshift @formatted, pop @digits;
			unshift @formatted, pop @digits;
			unshift @formatted, pop @digits;
			unshift @formatted, ',';
		}
		else {
			while (@digits) {
				unshift @formatted, pop @digits;
			}
		}
	}
	return join "", @formatted;
}


__END__

=head1 NAME

split_bam_paired_by_size.pl

=head1 SYNOPSIS

split_bam_paired_by_size.pl [--options] <file.bam>
  
  Options:
  --in <file.bam>
  --min <integer>
  --max <integer>
  --size <min-max>
  --out <filename>
  --at
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <file.bam>

Specify the file name of a binary BAM file as described for Samtools. 
Use samtools to convert text sam files to binary bam files. The file 
should be indexed; this program should be able to do it for you 
automatically if it is not (assuming the bam directory is writeable).

=item --min <integer>

Optionally specify the minimum size of fragment to include in the output
GFF file. The default value is 100 bp.

=item --max <integer>

Optionally specify the maximum size of fragment to include in the output
GFF file. The default value is 200 bp.

=item --size <min-max>

When multiple size ranges are desired, they may be specified using the 
size option. Define the minimum and maximum size as a range separated by 
a dash (no spaces). Use this option repeatedly for multiple size ranges. 
The size option takes precedence over the min and max options.

=item --out <filename>

Optionally specify the base name of the output file. The default is to use 
base input name. The output file names are appended with '.$min_$max'. 

=item --at

Boolean option to indicate that only fragments whose paired sequence reads 
end in a [AT] nucleotide should be included in the output GFF file. 
Micrococcal nuclease (MNase) cuts (almost?) exclusively at AT dinucleotides; 
this option ensures that the fragment is more likely derived from a MNase 
cut.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will read a BAM file containing paired-read alignments and
write a new BAM file containing only those successfully aligned paired 
reads whose insert size fall within the set minimum and maximum 
lengths inclusively (defaults of 100 and 200 bp, respectively). 

Multiple size ranges may be specified, and pairs within each range are 
written to separate files. 

Additionally, fragments may be checked for the presence of an A/T nucleotide 
at the 5' end of the sequence, as is usually produced by digestion with 
Micrococcal Nuclease.

The input and output files are BAM files as described by the Samtools 
project (http://samtools.sourceforge.net). 

The input file should be sorted and indexed prior to sizing. The output file 
will also be automatically re-sorted and re-indexed for you.

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







