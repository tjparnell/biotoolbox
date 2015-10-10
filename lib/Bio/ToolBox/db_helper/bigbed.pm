package Bio::ToolBox::db_helper::bigbed;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(mean);
use Bio::DB::BigBed;
our $VERSION = 1.33;


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
	
	# pass the required information
	unless (scalar @_ >= 7) {
		confess " At least seven arguments must be passed to collect BigBed scores!\n";
	}
	my ($chromo, $start, $stop, $strand, $stranded, $method, @bed_features) = @_;
		# method can be score, count, or length
	
	# initialize the score array
	# this will record score, count, or lengths per the method
	my @scores;
	
	# look at each bedfile
	# usually there is only one, but for stranded data there may be 
	# two bedfiles (+ and -), so we'll check each bed file for strand info
	foreach my $bedfile (@bed_features) {
	
		# open the bedfile
		my $bb = _get_bb($bedfile);
			
		# first check that the chromosome is present
		$chromo = $BIGBED_CHROMOS{$bedfile}{$chromo} or next;
		
		# collect the features overlapping the region
		my $bb_stream = $bb->features(
			-seq_id   => $chromo, 
			-start    => $start, 
			-end      => $stop,
			-iterator => 1,
		);
		
		# process each feature
		while (my $bed = $bb_stream->next_seq) {
			
			# First check whether the strand is acceptable
			if (
				$stranded eq 'all' # all data is requested
				or $bed->strand == 0 # unstranded data
				or ( 
					# sense data
					$strand == $bed->strand 
					and $stranded eq 'sense'
				) 
				or (
					# antisense data
					$strand != $bed->strand  
					and $stranded eq 'antisense'
				)
			) {
				# we have acceptable data to collect
			
				# store the appropriate datapoint
				if ($method eq 'score') {
					push @scores, $bed->score;
				}
				elsif ($method eq 'count') {
					$scores[0] += 1;
				}
				elsif ($method eq 'pcount') {
					$scores[0] += 1 if 
						($bed->start >= $start and $bed->end <= $stop);
				}
				elsif ($method eq 'length') {
					push @scores, $bed->length;
				}
			}
		}
	}

	# return collected data
	return @scores;
}




### Collect positioned BigBed scores
sub collect_bigbed_position_scores {
	
	# pass the required information
	unless (scalar @_ >= 7) {
		confess " At least seven arguments must be passed to collect BigBed position scores!\n";
	}
	my ($chromo, $start, $stop, $strand, $stranded, $method, @bed_features) = @_;
		# method can be score, count, or length
	
	# set up hash, either position => count or position => [scores]
	my %bed_data;
	
	# look at each bedfile
	# usually there is only one, but there may be more
	foreach my $bedfile (@bed_features) {
	
		# open the bedfile
		my $bb = _get_bb($bedfile);
			
		# first check that the chromosome is present
		$chromo = $BIGBED_CHROMOS{$bedfile}{$chromo} or next;
		
		# collect the features overlapping the region
		my $bb_stream = $bb->features(
			-seq_id   => $chromo, 
			-start    => $start, 
			-end      => $stop,
			-iterator => 1,
		);
		
		# process each feature
		while (my $bed = $bb_stream->next_seq) {
			
			# First check whether the strand is acceptable
			if (
				$stranded eq 'all' # all data is requested
				or $bed->strand == 0 # unstranded data
				or ( 
					# sense data
					$strand == $bed->strand 
					and $stranded eq 'sense'
				) 
				or (
					# antisense data
					$strand != $bed->strand  
					and $stranded eq 'antisense'
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
					$position >= $start and $position <= $stop
				);
				
				# store the appropriate datapoint
				# for score and length, we're putting these into an array
				if ($method eq 'score') {
					# perform addition to force the score to be a scalar value
					push @{ $bed_data{$position} }, $bed->score + 0;
				}
				elsif ($method eq 'count') {
					$bed_data{$position} += 1;
				}
				elsif ($method eq 'pcount') {
					$bed_data{$position} += 1 if 
						($bed->start <= $start and $bed->end <= $stop);
				}
				elsif ($method eq 'length') {
					# I hope that length is supported, but not sure
					# may have to calculate myself
					push @{ $bed_data{$position} }, $bed->length;
				}
			}
		}
	}

	# combine multiple datapoints at the same position
	if ($method eq 'score' or $method eq 'length') {
		# each value is an array of one or more datapoints
		# we will take the simple mean
		foreach my $position (keys %bed_data) {
			$bed_data{$position} = mean( @{$bed_data{$position}} );
		}
	}
	
	# return collected data
	return %bed_data;
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

This module supports the use of bigBed file in the biotoolbox scripts.
It is used to collect the dataset scores from a binary 
bigBed file (.bb). The file may be local or remote.

Scores may be restricted to strand by specifying the desired strandedness. 
For example, to collect transcription data over a gene, pass the strandedness 
value 'sense'. If the strand of the region database object (representing the 
gene) matches the strand of the bed feature, then the data for that bed 
feature is collected.  

For loading bigbed files into a Bio::DB database, see the biotoolbox perl 
script 'big_filegff3.pl'.

=head1 USAGE

The module requires Lincoln Stein's Bio::DB::BigBed to be installed. 

Load the module at the beginning of your program.

	use Bio::ToolBox::db_helper::bigbed;

It will automatically export the name of the subroutines. 

=over

=item collect_bigbed_scores

This subroutine will collect only the data values from a binary bigbed file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed seven or more arguments in the following order:

=over 4

=item 1. The chromosome or seq_id

=item 2. The start position of the segment to collect 

=item 3. The stop or end position of the segment to collect 

=item 4. The strand of the feature or segment.

The BioPerl strand values must be used, i.e. -1, 0, or 1.

=item 5. The strandedness of the bed elements to collect.

A scalar value representing the desired strandedness of the data 
to be collected. Acceptable values include "sense", "antisense", 
or "all". Only those scores which match the indicated 
strandedness are collected.

=item 6. The value type of the data to collect.
 
Acceptable values include score, count, pcount, and length.

   score returns the score of each bed element within the 
   region. Make sure the BigBed elements contain a score 
   column.
   
   count returns the number of elements that overlap the 
   search region. 
   
   pcount, or precise count, returns the count of elements 
   that only fall within the region and do not extend beyond
   the search region.
   
   length returns the lengths of all overlapping elements. 

=item 7. The paths to one or more BigBed files

Always provide the BigBed path. Opened BigBed file objects are 
cached. Both local and remote files are supported.

=back

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_bigbed_position_scores

This subroutine will collect the score values from a binary bigBed file 
for the specified database region keyed by position. 

The subroutine is passed the same arguments as collect_bigbed_scores().

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

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  



