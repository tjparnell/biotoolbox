package Bio::ToolBox::db_helper::useq;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(mean);
use Bio::DB::USeq;
our $VERSION = '1.15';


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

# Opened USeq db objects
our %OPENED_USEQ;
	# a cache for opened USeq databases, primarily for collecting scores

# The true statement
1; 



sub collect_useq_scores {
	
	# pass the required information
	unless (scalar @_ >= 7) {
		confess " At least seven arguments must be passed to collect useq scores!\n";
	}
	my ($chromo, $start, $stop, $strand, $stranded, $method, @useqs) = @_;
		# method can be score, count, or length
	
	# initialize the score array
	# this will record score, count, or lengths per the method
	my @scores;
	
	# adjust strand method
	# do not need to adjust if stranded is sense
	if ($stranded eq 'antisense') {
		$strand = $strand * -1;
	}
	elsif ($stranded eq 'all') {
		# Bio::DB::USeq will translate this properly, and collect from 
		# both strands as necessary
		$strand = 0;
	}
	
	# unlikely there are more than one useq file, but just in case
	foreach my $useqfile (@useqs) {
		
		# open a new db object
		my $useq;
		if (exists $OPENED_USEQ{$useqfile}) {
			# use a cached object
			$useq = $OPENED_USEQ{$useqfile};
		}
		else {
			# open and cache the bigWig object
			$useq = open_useq_db($useqfile) or 
				croak " Unable to open USeq file '$useqfile'! $!\n";
			$OPENED_USEQ{$useqfile} = $useq;
			%{ $USEQ_CHROMOS{$useqfile} } = map { $_ => 1 } $useq->seq_ids;
		}
		
		# check chromosome first
		next unless exists $USEQ_CHROMOS{$useqfile}{$chromo};
	
		# need to collect the scores based on the type of score requested
		
		if ($method eq 'score') {
			# need to collect scores
			my @region_scores = $useq->scores(
				-seq_id     => $chromo,
				-start      => $start,
				-end        => $stop,
				-strand     => $strand,
			);
			push @scores, @region_scores;
		}
		elsif ($method eq 'count') {
			# need to collect features across the region
			my $iterator = $useq->get_seq_stream(
				-seq_id     => $chromo,
				-start      => $start,
				-end        => $stop,
				-strand     => $strand,
			);
			return unless $iterator;
			
			# collect the lengths of each feature
			while (my $f = $iterator->next_seq) {
				push @scores, 1;
			}
		}
		elsif ($method eq 'length') {
			# need to collect features across the region
			my $iterator = $useq->get_seq_stream(
				-seq_id     => $chromo,
				-start      => $start,
				-end        => $stop,
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
	
	# pass the required information
	unless (scalar @_ >= 7) {
		confess " At least seven arguments must be passed to collect useq scores!\n";
	}
	my ($chromo, $start, $stop, $strand, $stranded, $method, @useqs) = @_;
		# method can be score, count, or length
	
	# initialize the score array
	# this will record score, count, or lengths per the method
	my %pos2score;
	
	# adjust strand method
	# do not need to adjust if stranded is sense
	if ($stranded eq 'antisense') {
		$strand = $strand * -1;
	}
	elsif ($stranded eq 'all') {
		# Bio::DB::USeq will translate this properly, and collect from 
		# both strands as necessary
		$strand = 0;
	}
	
	# unlikely there are more than one useq file, but just in case
	foreach my $useqfile (@useqs) {
		
		# open a new db object
		my $useq;
		if (exists $OPENED_USEQ{$useqfile}) {
			# use a cached object
			$useq = $OPENED_USEQ{$useqfile};
		}
		else {
			# open and cache the bigWig object
			$useq = open_useq_db($useqfile) or 
				croak " Unable to open USeq file '$useqfile'! $!\n";
			$OPENED_USEQ{$useqfile} = $useq;
			%{ $USEQ_CHROMOS{$useqfile} } = map { $_ => 1 } $useq->seq_ids;
		}
		
		# check chromosome first
		next unless exists $USEQ_CHROMOS{$useqfile}{$chromo};
	
		# collect the features overlapping the region
		my $iterator = $useq->get_seq_stream(
			-seq_id     => $chromo,
			-start      => $start,
			-end        => $stop,
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
				$position >= $start and $position <= $stop
			);
			
			# record the value
			if ($method eq 'score') {
				push @{ $pos2score{$position} }, $f->score;
			}
			elsif ($method eq 'count') {
				$pos2score{$position} += 1;
			}
			elsif ($method eq 'length') {
				push @{ $pos2score{$position} }, $f->length;
			}
		}
	}
	
	# combine multiple datapoints at the same position
	if ($method eq 'score' or $method eq 'length') {
		# each value is an array of one or more datapoints
		# we will take the simple mean
		foreach my $position (keys %pos2score) {
			$pos2score{$position} = mean( @{$pos2score{$position}} );
		}
	}
	
	# return collected data
	return %pos2score;
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



__END__

=head1 NAME

Bio::ToolBox::db_helper::useq

=head1 DESCRIPTION

This module supports the use of useq file in the Bio::ToolBox distribution.
Useq files are zip archives representing either intervals or scores. They 
may be used similarly to either bigWig or bigBed files. More information 
about useq files may be found at L<http://useq.sourceforge.net/useqArchiveFormat.html>.
USeq files use the extension F<.useq>.

Scores from useq files may be collected using this module. Either a single 
score from an interval, or a hash of scores associated with positions across 
an interval. 

Scores may be restricted to strand by specifying the desired strandedness. 
For example, to collect transcription data over a gene, pass the strandedness 
value 'sense'. If the strand of the region database object (representing the 
gene) matches the strand of the bed feature, then the data for that bed 
feature is collected.  

=head1 USAGE

The module requires the Bio::DB::USeq package to be installed. 

Load the module at the beginning of your program.

	use Bio::ToolBox::db_helper::useq;

It will automatically export the name of the subroutines. 

=over

=item collect_useq_scores()

This subroutine will collect only the data values from a binary useq file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed seven or more arguments in the following order:
    
    1) The chromosome or seq_id
    2) The start position of the segment to collect 
    3) The stop or end position of the segment to collect 
    4) The strand of the original feature (or region), -1, 0, or 1.
    5) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       or "all". Only those scores which match the indicated 
       strandedness are collected.
    6) The method or type of data collected. 
       Acceptable values include 'score', 'count' (returns the number 
       of features found), or 'length' (returns the lengths of the 
       features found).  
    7) The paths to one or more USeq files. It's unlikely to collect 
       from more than one, but, hey, if you want to....

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_useq_position_scores()

This subroutine will collect the score values from a binary useq file 
for the specified database region keyed by position. 

The subroutine is passed the same arguments as collect_useq_scores().

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. The feature midpoint is used 
as the key position. When multiple features are found at the same 
position, a simple mean (for score or length data methods) or sum 
(for count methods) is returned.

=item open_useq_db()

This subroutine will open a useq database connection. Pass the local 
path to a useq file (.useq extension). It will return the opened 
Bio::DB::USeq database object.

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
