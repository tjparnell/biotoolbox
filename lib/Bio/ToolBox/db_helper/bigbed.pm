package Bio::ToolBox::db_helper::bigbed;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(mean);
use Bio::DB::BigBed;
use constant {
	CHR  => 0,  # chromosome
	STRT => 1,  # start
	STOP => 2,  # stop
	STR  => 3,  # strand
	STND => 4,  # strandedness
	METH => 5,  # method
	RETT => 6,  # return type
	DB   => 7,  # database object
	DATA => 8,  # first dataset, additional may be present
};
our $VERSION = '1.50';


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_bigbed_scores
	collect_bigbed_position_scores
	open_bigbed_db
	sum_total_bigbed_features
);

# Hash of Bigfile chromosomes
our %BIGBED_CHROMOS;
	# sometimes user may request a chromosome that's not in the bigfile
	# that could lead to an exception
	# we will record the chromosomes list in this hash
	# $BIGBED_CHROMOS{bigfile}{chromos}
	# we also record the chromosome name variant with or without chr prefix
	# to accommodate different naming conventions

# Opened bigBed db objects
our %OPENED_BB;
	# a cache for opened BigBed databases, primarily for collecting scores
	# caching here is only for local purposes of collecting scores
	# db_helper also provides caching of db objects but with option to force open in
	# the case of forking processes - we don't have that here

# The true statement
1; 


### Collect BigBed scores only
sub collect_bigbed_scores {
	
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	
	# look at each bedfile
	# usually there is only one, but for stranded data there may be 
	# two bedfiles (+ and -), so we'll check each bed file for strand info
	my @scores;
	for (my $d = DATA; $d < scalar @$param; $d++) {
	
		# open the bedfile
		my $bb = _get_bb($param->[$d]);
			
		# first check that the chromosome is present
		my $chromo = $BIGBED_CHROMOS{$param->[$d]}{$param->[CHR]} or next;
		
		# collect the features overlapping the region
			# we are using the high level API rather than the low-level
			# since we getting the individual scores from each bed element
		my $bb_stream = $bb->features(
			-seq_id   => $chromo, 
			-start    => $param->[STRT], 
			-end      => $param->[STOP],
			-iterator => 1,
		);
		
		# process each feature
		while (my $bed = $bb_stream->next_seq) {
			
			# First check whether the strand is acceptable
			if (
				$param->[STND] eq 'all' # all data is requested
				or $bed->strand == 0 # unstranded data
				or ( 
					# sense data
					$param->[STR] == $bed->strand 
					and $param->[STND] eq 'sense'
				) 
				or (
					# antisense data
					$param->[STR] != $bed->strand  
					and $param->[STND] eq 'antisense'
				)
			) {
				# we have acceptable data to collect
			
				# store the appropriate datapoint
				if ($param->[METH] eq 'count') {
					push @scores, 1;
				}
				elsif ($param->[METH] eq 'pcount') {
					push @scores, 1 if ($bed->start >= $param->[STRT] and 
						$bed->end <= $param->[STOP]);
				}
				elsif ($param->[METH] eq 'ncount') {
					push @scores, $bed->display_name || $bed->primary_id;
				}
				else {
					push @scores, $bed->score;
				}
			}
		}
	}

	# return collected data
	
	return wantarray ? @scores : \@scores;
}




### Collect positioned BigBed scores
sub collect_bigbed_position_scores {
	
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	
	# look at each bedfile
	# usually there is only one, but there may be more
	my %pos2data;
	for (my $i = DATA; $i < scalar @$param; $i++) {
	
		# open the bedfile
		my $bb = _get_bb($param->[$i]);
			
		# first check that the chromosome is present
		my $chromo = $BIGBED_CHROMOS{$param->[$i]}{$param->[CHR]} or next;
		
		# collect the features overlapping the region
		my $bb_stream = $bb->features(
			-seq_id   => $chromo, 
			-start    => $param->[STRT], 
			-end      => $param->[STOP],
			-iterator => 1,
		);
		
		# process each feature
		while (my $bed = $bb_stream->next_seq) {
			
			# First check whether the strand is acceptable
			if (
				$param->[STND] eq 'all' # all data is requested
				or $bed->strand == 0 # unstranded data
				or ( 
					# sense data
					$param->[STR] == $bed->strand 
					and $param->[STND] eq 'sense'
				) 
				or (
					# antisense data
					$param->[STR] != $bed->strand  
					and $param->[STND] eq 'antisense'
				)
			) {
				# we have acceptable data to collect
			
				# determine position to record
				my $position;
				if ($bed->start == $bed->end) {
					# just one position recorded
					$position = $bed->start;
				}
				else {
					# calculate the midpoint
					$position = int( 
						( ($bed->start + $bed->end) / 2) + 0.5
					);
				}
				
				# check the position
				next unless (
					# want to avoid those whose midpoint are not technically 
					# within the region of interest
					$position >= $param->[STRT] and $position <= $param->[STOP]
				);
				
				# store the appropriate datapoint
				# for score and length, we're putting these into an array
				if ($param->[METH] eq 'count') {
					$pos2data{$position} += 1;
				}
				elsif ($param->[METH] eq 'pcount') {
					$pos2data{$position} += 1 if 
						($bed->start <= $param->[STRT] and $bed->end <= $param->[STOP]);
				}
				elsif ($param->[METH] eq 'ncount') {
					$pos2data{$position} ||= [];
					push @{ $pos2data{$position} }, $bed->display_name || 
						$bed->primary_id;
 				}
				else {
					# everything else we just take the score
					push @{ $pos2data{$position} }, $bed->score + 0;
				}
			}
		}
	}

	# combine multiple datapoints at the same position
	if ($param->[METH] eq 'ncount') {
		foreach my $position (keys %pos2data) {
			my %name2count;
			foreach (@{$pos2data{$position}}) { $name2count{$_} += 1 }
			$pos2data{$position} = scalar(keys %name2count);
		}
	}
	elsif ($param->[METH] eq 'count' or $param->[METH] eq 'pcount') {
		# do nothing, these aren't arrays
	}
	elsif ($param->[METH] eq 'mean') {
		foreach my $position (keys %pos2data) {
			$pos2data{$position} = mean( @{$pos2data{$position}} );
		}
	}
	elsif ($param->[METH] eq 'median') {
		foreach my $position (keys %pos2data) {
			$pos2data{$position} = median( @{$pos2data{$position}} );
		}
	}
	elsif ($param->[METH] eq 'min') {
		foreach my $position (keys %pos2data) {
			$pos2data{$position} = min( @{$pos2data{$position}} );
		}
	}
	elsif ($param->[METH] eq 'max') {
		foreach my $position (keys %pos2data) {
			$pos2data{$position} = max( @{$pos2data{$position}} );
		}
	}
	elsif ($param->[METH] eq 'sum') {
		foreach my $position (keys %pos2data) {
			$pos2data{$position} = sum( @{$pos2data{$position}} );
		}
	}
	else {
		# just take the mean for everything else
		foreach my $position (keys %pos2data) {
			$pos2data{$position} = mean( @{$pos2data{$position}} );
		}
	}
	
	# return collected data
	return wantarray ? %pos2data : \%pos2data;
}



### Open a bigBed database connection
sub open_bigbed_db {
	
	# check path
	my $bedfile = shift;
	my $path = $bedfile;
	$path =~ s/^file://; # clean up file prefix if present
	
	# open
	my $bb;
	eval {
		$bb = Bio::DB::BigBed->new($path);
	};
	return unless $bb;
	
	return $bb;
}



### Sum the total number of features in the bigBed file
sub sum_total_bigbed_features {
	
	# Passed arguments;
	my $bb_file = shift;
	unless ($bb_file) {
		carp " no BigBed file or BigBed db object passed!\n";
		return;
	}
	
	
	# Open BigBed file if necessary
	my $bb;
	my $bb_ref = ref $bb_file;
	if ($bb_ref =~ /Bio::DB::BigBed/) {
		# we have an opened bigbed db object
		$bb = $bb_file;
	}
	else {
		# we have a name of a bigbed file
		# open it but do not remember it
		$bb = open_bigbed_db($bb_file);
		return unless ($bb);
	}
	
	# Return the item count for the bigBed
	# wow, this is easy!
	return $bb->bf->bigBedItemCount;
		# this is the itemCount available to the low level bigFile object
		# the validCount method from summary or bin features simple records
		# the number of covered bases, not entirely useful here
}


### Internal subroutine for getting the cached bigbed object
sub _get_bb {
	my $bbfile = shift;
	
	return $OPENED_BB{$bbfile} if exists $OPENED_BB{$bbfile};
	
	# open and cache the bigWig object
	my $bb = open_bigbed_db($bbfile) or 
		croak " Unable to open bigBed file '$bbfile'! $!\n";
	$OPENED_BB{$bbfile} = $bb;
	
	# record the chromosomes and possible variants
	$BIGBED_CHROMOS{$bbfile} = {};
	foreach my $s ($bb->seq_ids) {
		$BIGBED_CHROMOS{$bbfile}{$s} = $s;
		if ($s =~ /^chr(.+)$/) {
			$BIGBED_CHROMOS{$bbfile}{$1} = $s;
		}
		else {
			$BIGBED_CHROMOS{$bbfile}{"chr$s"} = $s;
		}
	}
	return $bb;
}


__END__

=head1 NAME

Bio::ToolBox::db_helper::bigbed

=head1 DESCRIPTION

This module provides support for binary BigBed files to the 
L<Bio::ToolBox> package. 

=head1 USAGE

The module requires L<Bio::DB::BigBed> to be installed, which in turn 
requires the UCSC Kent C library to be installed.

In general, this module should not be used directly. Use the methods 
available in L<Bio::ToolBox::db_helper> or <Bio::ToolBox::Data>.  

All subroutines are exported by default.

=head2 Available subroutines

=over

=item collect_bigbed_scores

This subroutine will collect only the data values from a binary bigbed file 
for the specified database region. The positional information of the 
scores is not retained.

The subroutine is passed a parameter array reference. See below for details.

The subroutine returns an array or array reference of the requested dataset 
values found within the region of interest. 

=item collect_bigbed_position_scores

This subroutine will collect the score values from a binary bigBed file 
for the specified database region keyed by position. 

The subroutine is passed a parameter array reference. See below for details.

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. The feature midpoint is used 
as the key position. When multiple features are found at the same 
position, a simple mean (for score or length data methods) or sum 
(for count methods) is returned.

=item open_bigbed_db()

This subroutine will open a BigBed database connection. Pass either the 
local path to a bigBed file (.bb extension) or the URL of a remote bigBed 
file. It will return the opened database object.

The opened BigBed object is cached for later use. If you do not want this 
(for example, when forking), pass a second true argument.

=item sum_total_bigbed_features()

This subroutine will sum the total number of bed features present in a 
BigBed file. This may be useful, for example, in calculating fragments 
(reads) per million mapped values when the bigbed file represents 
sequence alignments.

Pass either the name of a bigBed file (.bb), either local or remote, or an 
opened BigBed database object. A scalar value of the total number of features 
is returned.

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

Acceptable values include score, count, and pcount.

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

=item 8 and higher. Paths to one or more BigBed files

Opened BigBed file objects are cached. Both local and remote files are 
supported.

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



