package Bio::ToolBox::db_helper::useq;

use warnings;
use strict;
use Carp;
use List::Util qw(min max sum);
use Statistics::Lite qw(median);
use Bio::ToolBox::db_helper::constants;
use Bio::DB::USeq;
require Exporter;

our $VERSION = '1.51';

# Exported names
our @ISA    = qw(Exporter);

## no critic
## this is never intended to be used directly by end users
## and exporting everything is required
our @EXPORT = qw(
	collect_useq_scores
	collect_useq_position_scores
	open_useq_db
);
## use critic

# Hash of USeq chromosomes
my %USEQ_CHROMOS;

# sometimes user may request a chromosome that's not in the useq file
# that could lead to an exception
# we will record the chromosomes list in this hash
# $USEQ_CHROMOS{useqfile}{chromos}
# we also record the chromosome name variant with or without chr prefix
# to accommodate different naming conventions

# Opened USeq db objects
my %OPENED_USEQ;

# a cache for opened USeq databases, primarily for collecting scores
# caching here is only for local purposes of collecting scores
# db_helper also provides caching of db objects but with option to force open in
# the case of forking processes - we don't have that here

sub collect_useq_scores {

	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;

	# adjust strand method
	my $strand;
	if ( $param->[STND] eq 'antisense' ) {
		$strand = $param->[STR] * -1;
	}
	elsif ( $param->[STND] eq 'all' ) {

		# Bio::DB::USeq will translate this properly, and collect from
		# both strands as necessary
		$strand = 0;
	}
	else {
		# default
		$strand = $param->[STR];
	}

	# unlikely there are more than one useq file, but just in case
	my @scores;
	for ( my $d = DATA; $d < scalar @{ $param }; $d++ ) {

		# open a new db object
		my $useq = _get_useq( $param->[$d] );

		# check chromosome first
		my $chromo = $USEQ_CHROMOS{ $param->[$d] }{ $param->[CHR] } or next;

		# need to collect the scores based on the type of score requested
		if ( $param->[METH] eq 'count' ) {

			# need to collect features across the region
			my $iterator = $useq->get_seq_stream(
				-seq_id => $chromo,
				-start  => $param->[STRT],
				-end    => $param->[STOP],
				-strand => $strand,
			);
			return unless $iterator;

			# count each feature
			while ( my $f = $iterator->next_seq ) {
				push @scores, 1;
			}
		}
		elsif ( $param->[METH] eq 'ncount' ) {

			# need to collect features across the region
			my $iterator = $useq->get_seq_stream(
				-seq_id => $chromo,
				-start  => $param->[STRT],
				-end    => $param->[STOP],
				-strand => $strand,
			);
			return unless $iterator;

			# store the names
			while ( my $f = $iterator->next_seq ) {
				push @scores, $f->display_name || $f->primary_id;

				# if no display name, a primary_id should automatically be generated
			}
		}
		elsif ( $param->[METH] eq 'pcount' ) {

			# need to collect features across the region
			my $iterator = $useq->get_seq_stream(
				-seq_id => $chromo,
				-start  => $param->[STRT],
				-end    => $param->[STOP],
				-strand => $strand,
			);
			return unless $iterator;

			# precisely count each feature
			while ( my $f = $iterator->next_seq ) {
				push @scores, 1
					if ( $f->start >= $param->[STRT] and $f->end <= $param->[STOP] );
			}
		}
		else {
			# everything else is just scores
			push @scores,
				$useq->scores(
					-seq_id => $chromo,
					-start  => $param->[STRT],
					-end    => $param->[STOP],
					-strand => $strand,
				);
		}
	}

	return wantarray ? @scores : \@scores;
}

sub collect_useq_position_scores {

	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;

	# adjust strand method
	my $strand;
	if ( $param->[STND] eq 'antisense' ) {
		$strand = $param->[STR] * -1;
	}
	elsif ( $param->[STND] eq 'all' ) {

		# Bio::DB::USeq will translate this properly, and collect from
		# both strands as necessary
		$strand = 0;
	}
	else {
		# default
		$strand = $param->[STR];
	}

	# unlikely there are more than one useq file, but just in case
	my %pos2score;
	for ( my $d = DATA; $d < scalar @{ $param }; $d++ ) {

		# open a new db object
		my $useq = _get_useq( $param->[$d] );

		# check chromosome first
		my $chromo = $USEQ_CHROMOS{ $param->[$d] }{ $param->[CHR] } or next;

		# collect the features overlapping the region
		my $iterator = $useq->get_seq_stream(
			-seq_id => $chromo,
			-start  => $param->[STRT],
			-end    => $param->[STOP],
			-strand => $strand,
		);
		return unless $iterator;

		# collect each feature
		while ( my $f = $iterator->next_seq ) {

			# determine position to record
			my $position;
			if ( $f->start == $f->end ) {

				# just one position recorded
				$position = $f->start;
			}
			else {
				# calculate the midpoint
				$position = int( ( ( $f->start + $f->end ) / 2 ) + 0.5 );
			}

			# check the position
			next
				unless (
					# want to avoid those whose midpoint are not technically
					# within the region of interest
					$position >= $param->[STRT] and $position <= $param->[STOP]
				);

			# record the value
			if ( $param->[METH] eq 'count' ) {
				$pos2score{$position} += 1;
			}
			elsif ( $param->[METH] eq 'ncount' ) {
				$pos2score{$position} ||= [];
				push @{ $pos2score{$position} }, $f->display_name
					|| $f->primary_id;
			}
			elsif ( $param->[METH] eq 'pcount' ) {
				$pos2score{$position} += 1
					if ( $f->start >= $param->[STRT] and $f->end <= $param->[STOP] );
			}
			else {
				# everything else we take the score
				push @{ $pos2score{$position} }, $f->score;
			}
		}
	}

	# combine multiple datapoints at the same position
	if ( $param->[METH] eq 'ncount' ) {
		foreach my $position ( keys %pos2score ) {
			my %name2count;
			foreach ( @{ $pos2score{$position} } ) { $name2count{$_} += 1 }
			$pos2score{$position} = scalar( keys %name2count );
		}
	}
	elsif ( $param->[METH] eq 'count' or $param->[METH] eq 'pcount' ) {

		# do nothing, these aren't arrays
	}
	elsif ( $param->[METH] eq 'mean' ) {
		foreach my $position ( keys %pos2score ) {
			$pos2score{$position} =
				sum( @{ $pos2score{$position} } ) / scalar( @{ $pos2score{$position} } );
		}
	}
	elsif ( $param->[METH] eq 'median' ) {
		foreach my $position ( keys %pos2score ) {
			$pos2score{$position} = median( @{ $pos2score{$position} } );
		}
	}
	elsif ( $param->[METH] eq 'min' ) {
		foreach my $position ( keys %pos2score ) {
			$pos2score{$position} = min( @{ $pos2score{$position} } );
		}
	}
	elsif ( $param->[METH] eq 'max' ) {
		foreach my $position ( keys %pos2score ) {
			$pos2score{$position} = max( @{ $pos2score{$position} } );
		}
	}
	elsif ( $param->[METH] eq 'sum' ) {
		foreach my $position ( keys %pos2score ) {
			$pos2score{$position} = sum( @{ $pos2score{$position} } );
		}
	}
	else {
		# just take the mean for everything else
		foreach my $position ( keys %pos2score ) {
			$pos2score{$position} =
				sum( @{ $pos2score{$position} } ) / scalar( @{ $pos2score{$position} } );
		}
	}

	# return collected data
	return wantarray ? %pos2score : \%pos2score;
}

sub open_useq_db {

	# path
	my $useqfile = shift;
	my $path     = $useqfile;
	$path =~ s/^file://;    # clean up file prefix if present

	# open
	my $useq;
	eval { $useq = Bio::DB::USeq->new($path); };
	return unless $useq;

	return $useq;
}

### Internal subroutine for getting the cached USeq object
sub _get_useq {
	my $useqfile = shift;

	return $OPENED_USEQ{$useqfile} if exists $OPENED_USEQ{$useqfile};

	# open and cache the USeq object
	my $useq = open_useq_db($useqfile)
		or croak " Unable to open USeq file '$useqfile'! $!\n";
	$OPENED_USEQ{$useqfile} = $useq;

	# record the chromosomes and possible variants
	$USEQ_CHROMOS{$useqfile} = {};
	foreach my $s ( $useq->seq_ids ) {
		$USEQ_CHROMOS{$useqfile}{$s} = $s;
		if ( $s =~ /^chr(.+)$/ ) {
			$USEQ_CHROMOS{$useqfile}{$1} = $s;
		}
		else {
			$USEQ_CHROMOS{$useqfile}{"chr$s"} = $s;
		}
	}
	return $useq;
}

1;

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
available in L<Bio::ToolBox::db_helper> or L<Bio::ToolBox::Data>.  

All subroutines are exported by default.

=over

=item open_useq_db

This subroutine will open a useq database connection. Pass the local 
path to a useq file (F<.useq> extension). It will return the opened 
L<Bio::DB::USeq> database object.

=item collect_useq_scores

This subroutine will collect only the data values from a binary useq file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed a parameter array reference. See 
L</"Data Collection Parameters Reference"> below for details.

The subroutine returns an array or array reference of the requested dataset 
values found within the region of interest. 

=item collect_useq_position_scores

This subroutine will collect the score values from a binary useq file 
for the specified database region keyed by position. 

The subroutine is passed a parameter array reference. See 
L</"Data Collection Parameters Reference"> below for details.

The subroutine returns a hash or hash reference of the defined dataset values 
found within the region of interest keyed by position. The feature midpoint 
is used as the key position. When multiple features are found at the same 
position, a simple mean (for score methods) or sum 
(for count methods) is returned.

=back

=head2 Data Collection Parameters Reference

The data collection subroutines are passed an array reference of parameters. 
The recommended  method for data collection is to use the 
L<Bio::ToolBox::db_helper/get_segment_score> method. 

The parameters array reference includes these items:

=over 4

=item 1. chromosome

=item 2. start coordinate

=item 3. stop coordinate 

Coordinates are in BioPerl-style 1-base system.

=item 4. strand

Should be standard BioPerl representation: -1, 0, or 1.

=item 5. strandedness

A scalar value representing the desired strandedness of the data 
to be collected. Acceptable values include "sense", "antisense", 
or "all". Only those scores which match the indicated 
strandedness are collected.

=item 6. score method

Acceptable values include score, count, ncount, and pcount.

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

=item 8. Path to USeq files

Additional USeq files may be appended to the list when merging. 
Opened USeq file objects are cached. 

=back

=head1 SEE ALSO

L<Bio::ToolBox::Data::Feature>, L<Bio::ToolBox::db_helper>, L<Bio::DB::USeq>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
