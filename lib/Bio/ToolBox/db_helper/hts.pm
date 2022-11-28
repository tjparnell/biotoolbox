package Bio::ToolBox::db_helper::hts;

use warnings;
use strict;
use Carp;
use English qw(-no_match_vars);
use Bio::ToolBox::db_helper::constants;
use Bio::ToolBox::db_helper::alignment_callbacks;
use Bio::DB::HTS;
require Exporter;

my $parallel;
eval {
	# check for parallel support, when counting bam alignments
	require Parallel::ForkManager;
	$parallel = 1;
};

our $VERSION = '1.70';

# Exported names
our @ISA    = qw(Exporter);

## no critic
## this is never intended to be used directly by end users
## and exporting everything is required
our @EXPORT = qw(
	open_bam_db
	open_indexed_fasta
	check_bam_index
	write_new_bam_file
	collect_bam_scores
	sum_total_bam_alignments
);
## use critic

# Hash of Bam chromosomes
my %BAM_CHROMOS;

# sometimes user may request a chromosome that's not in the bigfile
# that could lead to an exception
# we will record the chromosomes list in this hash
# $BAM_CHROMOS{bamfile}{chromos}
# we also record the chromosome name variant with or without chr prefix
# to accommodate different naming conventions

# Opened Bam db objects
my %OPENED_BAM;

# a cache for opened Bam files
# caching here is only for local purposes of collecting scores
# db_helper also provides caching of db objects but with option to force open in
# the case of forking processes - we don't have that here

### Open a bam database connection
sub open_bam_db {
	my $bamfile = shift;

	# check the path
	my $path = $bamfile;
	$path =~ s/^file://;    # strip the file prefix if present

	# check for bam index
	check_bam_index($path);

	# open the bam database object
	# we specifically do not cache the bam object or chromosome names here
	my $sam;
	eval { $sam = Bio::DB::HTS->new( -bam => $path ); };
	if ($sam) {
		return $sam;
	}
	elsif ($EVAL_ERROR or $OS_ERROR) {
		carp "ERROR: Unable to open Bam file: $EVAL_ERROR $OS_ERROR";
		return;
	}
	else {
		return;
	}
}

### Open an indexed fasta file
sub open_indexed_fasta {
	my $fasta = shift;
	eval { require Bio::DB::HTS::Faidx };    # this should be available
	my $fai = Bio::DB::HTS::Faidx->new($fasta);

	# this should automatically build the fai index if possible
	return $fai if defined $fai;
}

### Check for a bam index
sub check_bam_index {

	# HTSlib can accept either .bam.bai or .bai
	my $bamfile = shift;

	# check for possible index names
	return if -e "$bamfile.bai";     # .bam.bai
	return if -e "$bamfile.crai";    # .cram.crai
	return if -e "$bamfile.csi";     # .bam.csi
	my $alt_index = $bamfile;
	$alt_index =~ s/m$/i/i;          # change bam to bai or cram to crai
	return if -e $alt_index;

	# existing index can't be found, so make one
	Bio::DB::HTSfile->reindex($bamfile);
}

### Write a new bam file
sub write_new_bam_file {
	my $file = shift;
	$file .= '.bam' unless $file =~ /\.bam$/i;
	my $bam = Bio::DB::HTSfile->open( $file, 'wb' );
	carp "ERROR: unable to open bam file $file!" unless $bam;
	return $bam;
}

### Collect Bam scores
sub collect_bam_scores {

	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;

	# initialize score structures
	# which one is used depends on the return type variable
	my %pos2data;       # either position => count or position => [scores]
	my $scores = [];    # just scores

	# look at each bamfile
	# usually there is only one, but there may be more than one
	for ( my $b = DATA; $b < scalar @{ $param }; $b++ ) {

		## Open the Bam File
		my $bamfile = $param->[$b];
		my $bam     = $OPENED_BAM{$bamfile} || undef;
		unless ($bam) {

			# open and cache the bam file
			$bam = open_bam_db($bamfile)
				or croak "FATAL: Unable to open bam file '$bamfile'!";
			$OPENED_BAM{$bamfile} = $bam;

			# record the chromosomes and possible variants
			$BAM_CHROMOS{$bamfile} = {};
			foreach my $s ( $bam->seq_ids ) {
				$BAM_CHROMOS{$bamfile}{$s} = $s;
				if ( $s =~ /^chr(.+)$/ ) {
					$BAM_CHROMOS{$bamfile}{$1} = $s;
				}
				else {
					$BAM_CHROMOS{$bamfile}{"chr$s"} = $s;
				}
			}
		}

		# first check that the chromosome is present
		my $chromo = $BAM_CHROMOS{$bamfile}{ $param->[CHR] } or next;

		# convert coordinates into low level coordinates
		# consumed by the low level Bam API
		my ( $tid, $zstart, $end ) = $bam->header->parse_region(
			sprintf( "%s:%d-%d", $chromo, $param->[STRT], $param->[STOP] ) );

		## Collect the data according to the requested value type
		# we will either use simple coverage or alignments (count)
		if ( $param->[METH] =~ /count/ ) {

			# Need to collect and count alignments

			## Set the callback and a callback data structure
			my $callback = assign_callback($param);
			my %data     = (
				'scores' => $scores,
				'index'  => \%pos2data,
				'start'  => $param->[STRT],
				'stop'   => $param->[STOP],
			);

			# get the alignments
			# we are using the low level API to eke out performance
			$bam->hts_index->fetch( $bam->hts_file, $tid, $zstart, $end, $callback,
				\%data );
		}
		else {
			## Coverage
			# I am assuming everything else is working with read coverage

			# generate the coverage, this will ignore strand
			my $coverage = $bam->hts_index->coverage(
				$bam->hts_file,
				$tid,
				$zstart,    # 0-based coordinates
				$end,
			);

			# convert the coverage data
			# by default, this should return the coverage at 1 bp resolution
			if ( scalar @{ $coverage } ) {

				# check whether we need to index the scores
				if ( $param->[RETT] == 2 ) {
					for ( my $i = $param->[STRT]; $i <= $param->[STOP]; $i++ ) {

						# move the scores into the position score hash
						$pos2data{$i} += $coverage->[ $i - $param->[STRT] ];
					}
				}
				else {
					$scores = $coverage;
				}
			}
		}
	}

	## Return collected data
	if ( $param->[RETT] == 2 ) {
		return wantarray ? %pos2data : \%pos2data;
	}
	else {
		return wantarray ? @{ $scores } : $scores;
	}
}

### Determine total number of alignments in a bam file
sub sum_total_bam_alignments {

	# Passed arguments;
	my $sam_file = shift;
	my $min_mapq = shift || 0;    # by default we take all alignments
	my $paired   = shift || 0;    # by default we assume all alignments are single-end
	my $cpu      = shift || 2;    # number of forks to execute in parallel
	$cpu = 1 unless ($parallel);
	unless ($sam_file) {
		carp 'ERROR: no Bam file or bam db object passed!';
		return;
	}

	# Open Bam file if necessary
	my $bam;
	my $bam_ref = ref $sam_file;
	if ( $bam_ref =~ /Bio::DB::HTS/x ) {

		# we have an opened sam db object
		$bam = $sam_file;
	}
	else {
		# we have a name of a sam file
		$bam = open_bam_db($sam_file);
		return unless ($bam);
	}

	# prepare the counting subroutine
	my $counter = sub {
		my $tid    = shift;
		my $number = 0;

		# process the reads according to single or paired-end
		# paired end alignments
		if ($paired) {
			$bam->hts_index->fetch(
				$bam->hts_file,
				$tid, 0,
				$bam->target_len($tid),
				sub {
					my ( $a, $n ) = @_;

					# check alignment
					my $flag = $a->flag;
					return if $flag & 0x4;       # unmapped
					return if $flag & 0x0100;    # secondary alignment
					return if $flag & 0x0400;    # marked duplicate
					return if $flag & 0x0800;    # supplementary hit

					# check paired alignment
					return unless $flag & 0x2;      # proper_pair;
					return if $flag & 0x8;          # mate unmapped
					return if $flag & 0x10;         # reversed, only count left alignments
					return if $a->qual < $min_mapq;

					# count this fragment
					${$n}++;
				},
				\$number
			);
		}

		# single end alignments
		else {
			$bam->hts_index->fetch(
				$bam->hts_file,
				$tid, 0,
				$bam->target_len($tid),
				sub {
					my ( $a, $n ) = @_;

					# check alignment
					my $flag = $a->flag;
					return if $flag & 0x4;            # unmapped
					return if $flag & 0x0100;         # secondary alignment
					return if $flag & 0x0400;         # marked duplicate
					return if $flag & 0x0800;         # supplementary hit
					return if $a->qual < $min_mapq;

					# count this fragment
					${$n}++;
				},
				\$number
			);
		}
		return $number;
	};

	# Count the alignments on each chromosome
	my $total_read_number = 0;
	if ( $cpu > 1 ) {

		# count each chromosome in multiple parallel threads to speed things up

		# generate relatively equal lists of chromosome for each process based on length
		my @chromosomes = map { $_->[0] }
			sort { $b->[1] <=> $a->[1] }
			map { [ $_, $bam->target_len($_) ] } ( 0 .. $bam->n_targets - 1 );
		my @list;    # array of arrays, [process][chromosome id]
		my $i = 1;
		while (@chromosomes) {
			push @{ $list[$i] }, shift @chromosomes;
			$i++;
			$i = 1 if $i > $cpu;
		}

		# we will use Parallel ForkManager for convenience
		my $pm = Parallel::ForkManager->new($cpu);
		$pm->run_on_finish(
			sub {
				my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $count ) =
					@_;
				$total_read_number += ${$count};
			}
		);

		# Count the chromosomes in parallel processess
		foreach my $n ( 1 .. $cpu ) {
			$pm->start and next;

			### In child
			$bam->clone;    # to make it fork safe
			my $count = 0;
			foreach ( @{ $list[$n] } ) {

				# count each chromosome in this process list
				$count += &{$counter}($_);
			}
			$pm->finish( 0, \$count );
		}
		$pm->wait_all_children;
	}

	else {
		# loop through all the chromosomes in one execution thread
		for my $tid ( 0 .. $bam->n_targets - 1 ) {

			# each chromosome is internally represented in the bam file as
			# a numeric target identifier
			$total_read_number += &{$counter}($tid);
		}
	}

	# done
	return $total_read_number;
}

1;

__END__

=head1 NAME

Bio::ToolBox::db_helper::hts

=head1 DESCRIPTION

This module provides support for binary bam alignment files to the 
L<Bio::ToolBox> package through the HTSlib C library. 

=head1 USAGE

The module requires L<Bio::DB::HTS> to be installed, which in turn 
requires the HTSlib version 1.x.

In general, this module should not be used directly. Use the methods 
available in L<Bio::ToolBox::db_helper> or <Bio::ToolBox::Data>.  

All subroutines are exported by default.

=head2 Available subroutines

All subroutines are exported by default.

=over

=item open_bam_db

This subroutine will open a Bam database connection. Pass either the 
local path to a Bam file (F<.bam> extension) or the URL of a remote Bam 
file. A remote bam file must be indexed. A local bam file may be 
automatically indexed upon opening if the user has write permissions 
in the parent directory. 

It will return the opened database object.

=item open_indexed_fasta

This will open an indexed fasta file using the L<Bio::DB::HTS::Faidx> 
module. It requires a F<.fai> file to built, and one should be 
automatically built if it is not present. This provides a very fast 
interface to fetch genomic sequences, but no other support is 
provided. Pass the path to an uncompressed genomic fasta file 
(multiple sequences in one file is supported, but separate chromosome 
sequence files are not). The fasta index object is returned.

=item check_bam_index

This subroutine will check whether a bam index file is present and, 
if not, generate one. The L<Bio::DB::HTS> module uses the samtools 
style index extension, F<.bam.bai>, as opposed to the Picard style 
extension, F<.bai>. If a F<.bai> index is present, it will copy the 
file as F<.bam.bai> index. Unfortunately, a F<.bai> index cannot be 
used directly.

This method is called automatically prior to opening a bam database. 

=item write_new_bam_file

This subroutine will open a new empty Bam file. Pass the name of the 
new file as the argument. It will return a L<Bio::DB::HTSfile> object to 
which you can write a header followed by alignments. Be sure you know 
what to do before using this method! 

=item collect_bam_scores

This subroutine will collect only the data values from a binary bam file 
for the specified database region. The positional information of the 
scores is not retained.

Collected data values may be restricted to strand by specifying the desired 
strandedness (sense, antisense, or all), 
depending on the method of data collection. Collecting scores, or basepair 
coverage of alignments over the region of interest, does not currently support 
stranded data collection (as of this writing). However, enumerating 
alignments (count method) and collecting alignment lengths do support 
stranded data collection. Alignments are checked to see whether their midpoint 
is within the search interval before counting or length collected. 

As of version 1.30, paired-end bam files are properly handled with regards 
to strand; Strand is determined by the orientation of the first mate. However, 
pairs are still counted as two alignments, not one. To avoid this, use the 
value_type of 'ncount' and count the number of unique alignment names. 
(Previous versions treated all paired-end alignments as single-end alignments, 
severely limiting usefulness.)

The subroutine is passed a parameter array reference. See 
L</"Data Collection Parameters Reference"> below for details.

The subroutine returns an array or array reference of the requested dataset 
values found within the region of interest. 

=item collect_bam_position_scores

This subroutine will collect the score values from a binary bam file 
for the specified database region keyed by position. 

The subroutine is passed a parameter array reference. See 
L</"Data Collection Parameters Reference"> below for details.

The subroutine returns a hash or hash reference of the defined dataset values 
found within the region of interest keyed by position. The feature midpoint 
is used as the key position. When multiple features are found at the same 
position, a simple mean (for length data methods) or sum 
(for count methods) is returned. The ncount value type is not supported 
with positioned scores.

=item sum_total_bam_alignments

This subroutine will sum the total number of properly mapped alignments 
in a bam file. Pass the subroutine one to four arguments in the following 
order. 

=over 4

=item 1. Bam file path or object 

The name of the Bam file which should be counted. Alternatively, an 
opened L<Bio::DB::HTS> object may also be given. Required.

=item 2. Minimum mapping quality (integer)

Optionally pass the minimum mapping quality of the reads to be 
counted. The default is 0, where all alignments are counted.
Maximum is 255. See the SAM specification for details.

=item 3. Paired-end (boolean)

Optionally pass a boolean value (1 or 0) indicating whether 
the Bam file represents paired-end alignments. Only proper 
alignment pairs are counted. The default is to treat all 
alignments as single-end.

=item 4. Number of forks (integer)

Optionally pass the number of parallel processes to execute 
when counting alignments. Walking through a Bam file is 
time consuming but can be easily parallelized. The module 
L<Parallel::ForkManager> is required, and the default is a 
conservative two processes when it is installed.

=back
       
The subroutine will return the number of alignments.

=back

=head2 Data Collection Parameters Reference

The data collection subroutines are passed an array reference of parameters. 
The recommended  method for data collection is to use the
L<Bio::ToolBox::db_helper/get_segment_score> method. 

The parameters array reference includes these items:

=over 4

=item 1. chromosome

=item 1. start coordinate

=item 3. stop coordinate

=item 4. strand

Should be standard BioPerl representation: -1, 0, or 1.

=item 5. strandedness

A scalar value representing the desired strandedness of the data 
to be collected. Acceptable values include "sense", "antisense", 
or "all". Only those scores which match the indicated 
strandedness are collected.

=item 6. score method

Acceptable values include score, count, pcount, and ncount.

   * score returns the basepair coverage of alignments over the 
   region of interest
   
   * count returns the number of alignments that overlap the 
   search region. 
   
   * pcount, or precise count, returns the count of alignments 
   whose start and end fall within the region. 
   
   * ncount, or named count, returns an array of alignment read  
   names. Use this to avoid double-counting paired-end reads by 
   counting only unique names. Reads are taken if they overlap 
   the search region.

=item 7. Database object.

Not used here.

=item 8. Paths to bam file

Subsequent bam files may also be provided as additional list items.
Opened Bam file objects are cached. 

=back

=head1 SEE ALSO

L<Bio::ToolBox::Data::Feature>, L<Bio::ToolBox::db_helper>, L<Bio::DB::HTS>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

