package Bio::ToolBox::db_helper::useq;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(mean);
use Bio::DB::USeq;
our $VERSION = '1.50';


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_useq_scores
	collect_useq_position_scores
	open_useq_db
);

# Hash of USeq chromosomes
our %USEQ_CHROMOS;
	# sometimes user may request a chromosome that's not in the useq file
	# that could lead to an exception
	# we will record the chromosomes list in this hash
	# $USEQ_CHROMOS{useqfile}{chromos}
	# we also record the chromosome name variant with or without chr prefix
	# to accommodate different naming conventions

# Opened USeq db objects
our %OPENED_USEQ;
	# a cache for opened USeq databases, primarily for collecting scores
	# caching here is only for local purposes of collecting scores
	# db_helper also provides caching of db objects but with option to force open in
	# the case of forking processes - we don't have that here

# The true statement
1; 



sub collect_useq_scores {
	
	# passed options as array ref
	# chromosome, start, stop, strand, strandedness, method, value, db, dataset
	my $options = shift;
	
	# adjust strand method
	my $strand;
	if ($options->[4] eq 'antisense') {
		$strand = $options->[3] * -1;
	}
	elsif ($options->[4] eq 'all') {
		# Bio::DB::USeq will translate this properly, and collect from 
		# both strands as necessary
		$strand = 0;
	}
	else {
		# default
		$strand = $options->[3];
	}
	
	# unlikely there are more than one useq file, but just in case
	my @scores;
	for (my $i = 8; $i < scalar @$options; $i++) {
		
		# open a new db object
		my $useq = _get_useq($options->[$i]);
		
		# check chromosome first
		my $chromo = $USEQ_CHROMOS{$options->[$i]}{$options->[0]} or next;
	
		# need to collect the scores based on the type of score requested
		
		if ($options->[5] eq 'score') {
			# need to collect scores
			push @scores, $useq->scores(
				-seq_id     => $chromo,
				-start      => $options->[1], 
				-end        => $options->[2],
				-strand     => $strand,
			);
		}
		elsif ($options->[5] eq 'count') {
			# need to collect features across the region
			my $iterator = $useq->get_seq_stream(
				-seq_id     => $chromo,
				-start      => $options->[1], 
				-end        => $options->[2],
				-strand     => $strand,
			);
			return unless $iterator;
			
			# collect the lengths of each feature
			while (my $f = $iterator->next_seq) {
				$scores[0] += 1;
			}
		}
		elsif ($options->[5] eq 'pcount') {
			# need to collect features across the region
			my $iterator = $useq->get_seq_stream(
				-seq_id     => $chromo,
				-start      => $options->[1], 
				-end        => $options->[2],
				-strand     => $strand,
			);
			return unless $iterator;
			
			# collect the lengths of each feature
			while (my $f = $iterator->next_seq) {
				$scores[0] += 1 if 
					($f->start >= $start and $f->end <= $stop);
			}
		}
		elsif ($options->[5] eq 'length') {
			# need to collect features across the region
			my $iterator = $useq->get_seq_stream(
				-seq_id     => $chromo,
				-start      => $options->[1], 
				-end        => $options->[2],
				-strand     => $strand,
			);
			return unless $iterator;
			
			# collect the lengths of each feature
			while (my $f = $iterator->next_seq) {
				push @scores, $f->length;
			}
		}
		else {
			confess " unrecognized method $method!";
		}
	}
	
	return @scores;
}



sub collect_useq_position_scores {
	
	# passed options as array ref
	# chromosome, start, stop, strand, strandedness, method, value, db, dataset
	my $options = shift;
	
	# adjust strand method
	my $strand;
	if ($options->[4] eq 'antisense') {
		$strand = $options->[3] * -1;
	}
	elsif ($options->[4] eq 'all') {
		# Bio::DB::USeq will translate this properly, and collect from 
		# both strands as necessary
		$strand = 0;
	}
	else {
		# default
		$strand = $options->[3];
	}
	
	# unlikely there are more than one useq file, but just in case
	my %pos2score;
	for (my $i = 8; $i < scalar @$options; $i++) {
		
		# open a new db object
		my $useq = _get_useq($options->[$i]);
		
		# check chromosome first
		my $chromo = $USEQ_CHROMOS{$options->[$i]}{$options->[0]} or next;
	
		# collect the features overlapping the region
		my $iterator = $useq->get_seq_stream(
			-seq_id     => $chromo,
			-start      => $options->[1], 
			-end        => $options->[2],
			-strand     => $strand,
		);
		return unless $iterator;
		
		# collect the lengths of each feature
		while (my $f = $iterator->next_seq) {
			
			# determine position to record
			my $position;
			if ($f->start == $f->end) {
				# just one position recorded
				$position = $f->start;
			}
			else {
				# calculate the midpoint
				$position = int( 
					( ($f->start + $f->end) / 2) + 0.5
				);
			}
			
			# check the position
			next unless (
				# want to avoid those whose midpoint are not technically 
				# within the region of interest
				$position >= $options->[1] and $position <= $options->[2]
			);
			
			# record the value
			if ($options->[5] eq 'score') {
				push @{ $pos2score{$position} }, $f->score;
			}
			elsif ($options->[5] eq 'count') {
				$pos2score{$position} += 1;
			}
			elsif ($options->[5] eq 'pcount') {
				$pos2score{$position} += 1 if 
					($f->start >= $options->[1] and $f->end <= $options->[2]);
			}
			elsif ($options->[5] eq 'length') {
				push @{ $pos2score{$position} }, $f->length;
			}
		}
	}
	
	# combine multiple datapoints at the same position
	if ($options->[5] eq 'score' or $options->[5] eq 'length') {
		# each value is an array of one or more datapoints
		# we will take the simple mean
		foreach my $position (keys %pos2score) {
			$pos2score{$position} = mean( @{$pos2score{$position}} );
		}
	}
	
	# return collected data
	return wantarray ? %pos2score : \%pos2score;
}



sub open_useq_db {
	
	# path
	my $useqfile = shift;
	my $path = $useqfile;
	$path =~ s/^file://; # clean up file prefix if present
	
	# open
	my $useq;
	eval {
		$useq = Bio::DB::USeq->new($path);
	};
	return unless $useq;
	
	return $useq;
}



### Internal subroutine for getting the cached USeq object
sub _get_useq {
	my $useqfile = shift;
	
	return $OPENED_USEQ{$useqfile} if exists $OPENED_USEQ{$useqfile};
	
	# open and cache the USeq object
	my $useq = open_useq_db($useqfile) or 
		croak " Unable to open USeq file '$useqfile'! $!\n";
	$OPENED_USEQ{$useqfile} = $useq;
	
	# record the chromosomes and possible variants
	$USEQ_CHROMOS{$useqfile} = {};
	foreach my $s ($useq->seq_ids) {
		$USEQ_CHROMOS{$useqfile}{$s} = $s;
		if ($s =~ /^chr(.+)$/) {
			$USEQ_CHROMOS{$useqfile}{$1} = $s;
		}
		else {
			$USEQ_CHROMOS{$useqfile}{"chr$s"} = $s;
		}
	}
	return $useq;
}



__END__

=head1 NAME

Bio::ToolBox::db_helper::useq

=head1 DESCRIPTION

This module provides support for USeq files to the L<Bio::ToolBox> package. 
Useq files are zip archives representing either intervals or scores. They 
may be used similarly to either bigWig or bigBed files. More information 
about useq files may be found at L<http://useq.sourceforge.net/useqArchiveFormat.html>.
USeq files use the extension F<.useq>.

=head1 USAGE

The module requires L<Bio::DB::USeq> to be installed.

In general, this module should not be used directly. Use the methods 
available in L<Bio::ToolBox::db_helper> or <Bio::ToolBox::Data>.  

All subroutines are exported by default.

=over

=item open_useq_db()

This subroutine will open a useq database connection. Pass the local 
path to a useq file (.useq extension). It will return the opened 
Bio::DB::USeq database object.

=item collect_useq_scores()

This subroutine will collect only the data values from a binary useq file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed a parameter array reference. See below for details.

The subroutine returns an array or array reference of the requested dataset 
values found within the region of interest. 

=item collect_useq_position_scores()

This subroutine will collect the score values from a binary useq file 
for the specified database region keyed by position. 

The subroutine is passed a parameter array reference. See below for details.

The subroutine returns a hash or hash reference of the defined dataset values 
found within the region of interest keyed by position. The feature midpoint 
is used as the key position. When multiple features are found at the same 
position, a simple mean (for score or length data methods) or sum 
(for count methods) is returned.

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

Not used here. 

=item 7. The value type of data to collect

Acceptable values include score, count, pcount, ncount, and length.

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
   
   length returns the lengths of all overlapping alignments 

=item 8. A database object.

Not used here.

=item 9 and higher. Paths to one or more USeq files

Opened USeq file objects are cached. 

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
