package Bio::ToolBox::db_helper::bam;

# modules
require Exporter;
use strict;
use Carp;
use File::Copy;
use Bio::ToolBox::db_helper::constants;
use Bio::ToolBox::db_helper::alignment_callbacks;
use Bio::DB::Sam;
our $parallel;
eval {
	# check for parallel support, when counting bam alignments
	require Parallel::ForkManager;
	$parallel = 1;
};
our $VERSION = '1.51';

# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	open_bam_db
	open_indexed_fasta
	check_bam_index
	write_new_bam_file
	collect_bam_scores
	sum_total_bam_alignments
);

# Hash of Bam chromosomes
our %BAM_CHROMOS;
	# sometimes user may request a chromosome that's not in the bigfile
	# that could lead to an exception
	# we will record the chromosomes list in this hash
	# $BAM_CHROMOS{bamfile}{chromos}
	# we also record the chromosome name variant with or without chr prefix
	# to accommodate different naming conventions

# Opened Bam db objects
our %OPENED_BAM;
	# a cache for opened Bam files 
	# caching here is only for local purposes of collecting scores
	# db_helper also provides caching of db objects but with option to force open in
	# the case of forking processes - we don't have that here

# The true statement
1; 


### Open a bam database connection
sub open_bam_db {
	my $bamfile = shift;
	
	# check the path
	my $path = $bamfile;
	$path =~ s/^file://; # strip the file prefix if present
	
	# check for bam index
	check_bam_index($path);
	
	# open the bam database object
	my $sam;
	eval {
		$sam = Bio::DB::Sam->new(-bam => $path);
	};
	return unless $sam;
	# we specifically do not cache the bam object or chromosome names here
	
	return $sam;
}


### Open an indexed fasta file
sub open_indexed_fasta {
	my $fasta = shift;
	if ($fasta =~ /\.gz$/) {
		die " Bio::DB::Sam::Fai doesn't support compressed fasta files! Please decompress\n";
	}
	my $fai;
	eval {$fai = Bio::DB::Sam::Fai->load($fasta)};
		# this should automatically build the fai index if possible
	return $fai if defined $fai;
}


### Check for a bam index 
sub check_bam_index {
	# the old samtools always expects a .bam.bai index file
	# I find that relying on -autoindex yields a flaky Bio::DB::Sam object that 
	# doesn't always work as expected. Best to create the index BEFORE opening
	
	my $bamfile = shift;
	return if ($bamfile =~ /^(?:http|ftp)/i); # I can't do much with remote files
	
	# we will check the modification time to make sure index is newer
	my $bam_mtime = (stat($bamfile))[9];
	
	# optional index names
	my $bam_index = "$bamfile.bai"; # .bam.bai
	my $alt_index = $bamfile;
	$alt_index =~ s/bam$/bai/i; # picard uses .bai instead of .bam.bai as samtools does
	
	# check for existing index
	if (-e $bam_index) {
		if ( (stat($bam_index))[9] < $bam_mtime) {
			# index is older than bam file
			print " index $bam_index is old. Attempting to update time stamp.\n";
			my $now = time;
			utime($now, $now, $bam_index) || Bio::DB::Bam->reindex($bamfile);
		}
	}
	elsif (-e $alt_index) {
		if ( (stat($alt_index))[9] < $bam_mtime) {
			# index is older than bam file
			print " index $alt_index is old.\n";
		}
		# reuse this index
		copy($alt_index, $bam_index);
	}
	else {
		# make a new index
		Bio::DB::Bam->reindex($bamfile);
	}
}



### Write a new bam file
sub write_new_bam_file {
	my $file = shift;
	$file .= '.bam' unless $file =~ /\.bam$/i;
	my $bam = Bio::DB::Bam->open($file, 'w');
	carp "unable to open bam file $file!\n" unless $bam;
	return $bam;
}



### Collect Bam scores
sub collect_bam_scores {
	
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	
	# initialize score structures
	# which one is used depends on the return type variable
	my %pos2data; # either position => count or position => [scores]
	my $scores = []; # just scores
	
	# look at each bamfile
	# usually there is only one, but there may be more than one
	for (my $b = DATA; $b < scalar @$param; $b++) {
	
		## Open the Bam File
		my $bamfile = $param->[$b];
		my $bam = $OPENED_BAM{$bamfile} || undef;
		unless ($bam) {
			# open and cache the bam file
			$bam = open_bam_db($bamfile) or 
				croak " Unable to open bam file '$bamfile'! $!\n";
			$OPENED_BAM{$bamfile} = $bam;
	
			# record the chromosomes and possible variants
			$BAM_CHROMOS{$bamfile} = {};
			foreach my $s ($bam->seq_ids) {
				$BAM_CHROMOS{$bamfile}{$s} = $s;
				if ($s =~ /^chr(.+)$/) {
					$BAM_CHROMOS{$bamfile}{$1} = $s;
				}
				else {
					$BAM_CHROMOS{$bamfile}{"chr$s"} = $s;
				}
			}
		}
			
		# first check that the chromosome is present
		my $chromo = $BAM_CHROMOS{$bamfile}{$param->[CHR]} or next;
		
		# convert coordinates into low level coordinates
		# consumed by the low level Bam API
		my ($tid, $zstart, $end) = $bam->header->parse_region(
			sprintf("%s:%d-%d", $chromo, $param->[STRT], $param->[STOP]) );
	
		
		## Collect the data according to the requested value type
		# we will either use simple coverage or alignments (count)
		if ($param->[METH] =~ /count/) {
			# Need to collect and count alignments
			
			## Set the callback and a callback data structure
			my $callback = assign_callback($param);
			my %data = (
				'scores' => $scores,
				'index'  => \%pos2data,
				'start'  => $param->[STRT],
				'stop'   => $param->[STOP],
			);
			
			# get the alignments
			# we are using the low level API to eke out performance
			$bam->bam_index->fetch($bam->bam, $tid, $zstart, $end, $callback, \%data);
		}
		else {
			## Coverage
			# I am assuming everything else is working with read coverage
			
			# generate the coverage, this will ignore strand
			my $coverage = $bam->bam_index->coverage(
				$bam->bam,
				$tid,
				$zstart, # 0-based coordinates
				$end,
			);
			
			# convert the coverage data
			# by default, this should return the coverage at 1 bp resolution
			if (scalar @$coverage) {
				# check whether we need to index the scores
				if ($param->[RETT] == 2) {
					for (my $i = $param->[STRT]; $i <= $param->[STOP]; $i++) {
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
	
	# process the ncount arrays
	if ($param->[RETT] == 2 and $param->[METH] eq 'ncount') {
		foreach my $position (keys %pos2data) {
			my %name2count;
			foreach (@{$pos2data{$position}}) { $name2count{$_} += 1 }
			$pos2data{$position} = scalar(keys %name2count);
		}
	}
	
	## Return collected data
	if ($param->[RETT] == 2) {
		return wantarray ? %pos2data : \%pos2data;
	}
	else {
		return wantarray ? @$scores : $scores;
	}
}


### Determine total number of alignments in a bam file
sub sum_total_bam_alignments {
	
	# Passed arguments;
	my $sam_file = shift;
	my $min_mapq = shift || 0; # by default we take all alignments
	my $paired   = shift || 0; # by default we assume all alignments are single-end
	my $cpu      = shift || 2; # number of forks to execute in parallel
	$cpu = 1 unless ($parallel);
	unless ($sam_file) {
		carp " no Bam file or bam db object passed!\n";
		return;
	}
	
	
	# Open Bam file if necessary
	my $sam;
	my $sam_ref = ref $sam_file;
	if ($sam_ref =~ /Bio::DB::Sam/) {
		# we have an opened sam db object
		$sam = $sam_file;
	}
	else {
		# we have a name of a sam file
		$sam = open_bam_db($sam_file);
		return unless ($sam);
	}
	
	# prepare the counting subroutine
	my $counter = sub {
		my $tid = shift;
		my $number = 0;
		
		# process the reads according to single or paired-end
		# paired end alignments
		if ($paired) {
			$sam->bam_index->fetch(
				$sam->bam, 
				$tid, 
				0, 
				$sam->target_len($tid), 
				sub {
					my ($a, $number) = @_;
					
					# check paired alignment
					return unless $a->proper_pair;
					return if $a->reversed; # only count left alignments
					return if $a->qual < $min_mapq;
					
					# count this fragment
					$$number++;
				}, 
				\$number
			);
		}
		
		# single end alignments
		else {
			$sam->bam_index->fetch(
				$sam->bam, 
				$tid, 
				0, 
				$sam->target_len($tid), 
				sub {
					my ($a, $number) = @_;
					
					# check alignment
					return if $a->unmapped;
					return if $a->qual < $min_mapq;
					
					# count this fragment
					$$number++;
				}, 
				\$number
			);
		}
		return $number;
	};
	
	# Count the alignments on each chromosome
	my $total_read_number = 0;
	if ($cpu > 1) {
		# count each chromosome in multiple parallel threads to speed things up 
		
		# generate relatively equal lists of chromosome for each process based on length
		my @chromosomes = map { $_->[0] }
			sort { $b->[1] <=> $a->[1] }
			map { [$_, $sam->target_len($_)] } 
			(0 .. $sam->n_targets - 1);
		my @list; # array of arrays, [process][chromosome id]
		my $i = 1;
		while (@chromosomes) {
			push @{ $list[$i] }, shift @chromosomes;
			$i++;
			$i = 1 if $i > $cpu;
		}
		
		# we will use Parallel ForkManager for convenience
		my $pm = Parallel::ForkManager->new($cpu);
		$pm->run_on_finish( sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $count) = @_;
			$total_read_number += $$count;
		});
		
		# Count the chromosomes in parallel processess
		foreach my $n (1 .. $cpu) {
			$pm->start and next;
	
			### In child
			$sam->clone; # to make it fork safe
			my $count = 0;
			foreach ( @{$list[$n]} ) {
				# count each chromosome in this process list
				$count += &{$counter}($_);
			}
			$pm->finish(0, \$count); 
		}
		$pm->wait_all_children;
	}
	
	else {
		# loop through all the chromosomes in one execution thread
		for my $tid (0 .. $sam->n_targets - 1) {
			# each chromosome is internally represented in the bam file as 
			# a numeric target identifier
			$total_read_number += &{$counter}($tid);
		}
	}
	
	# done
	return $total_read_number;
}



__END__

=head1 NAME

Bio::ToolBox::db_helper::bam

=head1 DESCRIPTION

This module provides support for binary bam alignment files to the 
L<Bio::ToolBox> package through the samtools C library. 

=head1 USAGE

The module requires L<Bio::DB::Sam> to be installed, which in turn 
requires the samtools C library version 1.20 or less to be installed.

In general, this module should not be used directly. Use the methods 
available in L<Bio::ToolBox::db_helper> or <Bio::ToolBox::Data>.  

All subroutines are exported by default.

=head2 Available subroutines

All subroutines are exported by default.

=over

=item open_bam_db()

This subroutine will open a Bam database connection. Pass either the 
local path to a Bam file (.bam extension) or the URL of a remote Bam 
file. A remote bam file must be indexed. A local bam file may be 
automatically indexed upon opening if the user has write permissions 
in the parent directory. 

It will return the opened database object.

=item open_indexed_fasta()

This will open an indexed fasta file using the L<Bio::DB::Sam::Fai> 
module. It requires a F<.fa.fai> file to built, and one should be 
automatically built if it is not present. This provides a very fast 
interface to fetch genomic sequences, but no other support is 
provided. Pass the path to an uncompressed genomic fasta file 
(multiple sequences in one file is supported, but separate chromosome 
sequence files are not). The fasta index object is returned.

=item check_bam_index()

This subroutine will check whether a bam index file is present and, 
if not, generate one. The L<Bio::DB::Sam> module uses the samtools 
style index extension, F<.bam.bai>, as opposed to the picard style 
extension, F<.bai>. If a F<.bai> index is present, it will copy the 
file as F<.bam.bai> index. Unfortunately, a F<.bai> index cannot be 
used directly.

This method is called automatically prior to opening a bam database. 

=item write_new_bam_file()

This subroutine will open a new empty Bam file. Pass the name of the 
new file as the argument. It will return a Bio::DB::Bam object to 
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

The subroutine is passed a parameter array reference. See below for details.

The subroutine returns an array or array reference of the requested dataset 
values found within the region of interest. 

=item collect_bam_position_scores

This subroutine will collect the score values from a binary bam file 
for the specified database region keyed by position. 

The subroutine is passed a parameter array reference. See below for details.

The subroutine returns a hash or hash reference of the defined dataset values 
found within the region of interest keyed by position. The feature midpoint 
is used as the key position. When multiple features are found at the same 
position, a simple mean (for length data methods) or sum 
(for count methods) is returned. The ncount value type is not supported 
with positioned scores.

=item sum_total_bam_alignments()

This subroutine will sum the total number of properly mapped alignments 
in a bam file. Pass the subroutine one to four arguments in the following 
order. 

=over 4

=item 1. Bam file path or object 

The name of the Bam file which should be counted. Alternatively, an 
opened Bio::DB::Sam object may also be given. Required.

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
Parallel::ForkManager is required, and the default is a 
conservative two processes when it is installed.

=back
       
The subroutine will return the number of alignments.

=back

=head2 Data Collection Parameters Reference

The data collection subroutines are passed an array reference of parameters. 
The recommended  method for data collection is to use get_segment_score() method from 
L<Bio::ToolBox::db_helper>. 

The parameters array reference includes these items:

=over 4

=item 1. The chromosome or seq_id

=item 1. The start position of the segment to collect 

=item 3. The stop or end position of the segment to collect 

=item 4. The strand of the segment to collect

Should be standard BioPerl representation: -1, 0, or 1.

=item 5. The strandedness of the data to collect 

A scalar value representing the desired strandedness of the data 
to be collected. Acceptable values include "sense", "antisense", 
or "all". Only those scores which match the indicated 
strandedness are collected.

=item 6. The method for combining scores.

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
   
=item 7. A database object.

Not used here.

=item 8 and higher. Paths to one or more Bam files

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

