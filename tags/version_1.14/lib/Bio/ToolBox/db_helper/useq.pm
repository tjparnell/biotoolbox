package Bio::ToolBox::db_helper::useq;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(mean);
use Bio::DB::USeq;
our $VERSION = '1.14';


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_useq_scores
	collect_useq_position_scores
	open_useq_db
);

# Hashes of opened file objects
our %OPENED_USEQFILES; # opened useq file objects

# Hash of Bigfile chromosomes
our %USEQ_CHROMOS;
	# sometimes user may request a chromosome that's not in the useq file
	# that could lead to an exception
	# we will record the chromosomes list in this hash
	# $USEQ_CHROMOS{useqfile}{chromos}

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
		my $useq = open_useq_db($useqfile);
		
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
		my $useq = open_useq_db($useqfile);
		
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
	
	my $useqfile = shift;
	my $useq;
	
	# check whether the file has been opened or not
	if (exists $OPENED_USEQFILES{$useqfile} ) {
		# this file is already opened, use it
		$useq = $OPENED_USEQFILES{$useqfile};
	}
	
	else {
		# this file has not been opened yet, open it
		my $path = $useqfile;
		$path =~ s/^file://; # clean up file prefix if present
		eval {
			$useq = Bio::DB::USeq->new($path);
		};
		return unless $useq;
		
		# store the opened object for later use
		$OPENED_USEQFILES{$useqfile} = $useq;
		
		# collect the chromosomes for this useq
		%{ $USEQ_CHROMOS{$useqfile} } = map { $_ => 1 } $useq->seq_ids;
	}
	
	return $useq;
}