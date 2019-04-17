package Bio::ToolBox::db_helper::big;

# modules
require Exporter;
use strict;
use Carp;
use List::Util qw(min max sum);
use Bio::ToolBox::db_helper::constants;
use Bio::DB::Big;
our $VERSION = '1.66';

# Initialize CURL buffers
BEGIN {
	# not clear if this should be done only once or if it's harmless to re-init
	# for every new file, so I guess best to just do it here at the very beginning
	# initialization is only for remote files 
	Bio::DB::Big->init();
}

# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	open_bigwig_db
	collect_bigwig_score
	collect_bigwig_scores
	collect_bigwig_position_scores
	open_bigbed_db
	collect_bigbed_scores
	collect_bigbed_position_scores
	open_bigwigset_db
	collect_bigwigset_score
	collect_bigwigset_scores
	collect_bigwigset_position_scores
	sum_total_bigbed_features
);

# Hash of Bigfile chromosomes
my %BIG_CHROMOS;
	# sometimes user may request a chromosome that's not in the bigfile
	# that could lead to an exception
	# we will record the chromosomes list in this hash
	# $BIG_CHROMOS{bigfile}{altchromo} = chromo
	# we also record the chromosome name variant with or without chr prefix
	# to accommodate different naming conventions

# Hash of Bigfile chromosome lengths
my %BIG_CHROMOLENGTHS;
	# since libBigWig doesn't internally clip chromosome lengths
	# we will cache chromosome lengths and check them before it leads to an exception
	# $BIG_CHROMOLENGTHS{bigfile}{chromo} = length

# Opened bigFile db objects
my %OPENED_BIG;
	# a cache for opened Bigfile databases, primarily for collecting scores
	# caching here is only for local purposes of collecting scores
	# db_helper also provides caching of db objects but with option to force open in
	# the case of forking processes - we don't have that here

# BigWigSet bigWig IDs
my %BIGWIGSET_WIGS; 
	# cache for the bigwigs from a BigWigSet used in a query
	# we want to use low level bigWig access which isn't normally 
	# available from the high level BigWigSet, so we identify the 
	# bigWigs from the bigWigSet and cache them here 


# The true statement
1; 



#### BigWig Subroutines

sub open_bigwig_db {
	my $path = shift;
	my $bw = _open_big($path) or 
		croak " Unable to open bigWig file '$path'! $@\n";
	unless ($bw->is_big_wig) {
		croak " $path is not a bigWig file!\n";
	}
	return $bw; # we do not cache here
}


sub collect_bigwig_score {
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset(s)
	my $param = shift;
	
	# this will use the built-in summary statistics
	# the following methods are available: mean, min, max, std, cov
	
	# adjust old method name for compatibility 
	$param->[METH] = 'std' if $param->[METH] eq 'stddev';
		
	# check how many features we have
	if (scalar @$param == 9) {
		# only one bw, great!
		my $bw = _get_bigwig($param->[DATA]);
		my $chromo = $BIG_CHROMOS{$param->[DATA]}{$param->[CHR]} or return;
		$param->[STRT] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
			$param->[STRT] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
		$param->[STOP] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
			$param->[STOP] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
		my $s = $bw->get_stats($chromo, $param->[STRT] - 1, $param->[STOP], 1, 
			$param->[METH]);
		return $s->[0];
	}
	else {
		# we have multiple bigwigs
		my @scores;
		for (my $d = DATA; $d < scalar @$param; $d++) {
			my $bw = _get_bigwig($param->[$d]);
			my $chromo = $BIG_CHROMOS{$param->[$d]}{$param->[CHR]} or next;
			$param->[STRT] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
				$param->[STRT] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
			$param->[STOP] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
				$param->[STOP] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
			my $s = $bw->get_stats($chromo, $param->[STRT] - 1, $param->[STOP], 1, 
				$param->[METH]);
			push @scores, $s->[0];
		}
		if ($param->[METH] eq 'min') {
			return min(@scores);
		}
		elsif ($param->[METH] eq 'max') {
			return max(@scores);
		}
		else {
			confess sprintf " how did we get here? unable to calculate %s from multiple bigwigs!\n", 
				$param->[METH];
		}
	}
}


sub collect_bigwig_scores {
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset(s)
	my $param = shift;
	
	# check how many features we have
	if (scalar @$param == 9) {
		# only one, great!
		my $bw = _get_bigwig($param->[DATA]);
		my $chromo = $BIG_CHROMOS{$param->[DATA]}{$param->[CHR]} or return;
		$param->[STRT] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
			$param->[STRT] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
		$param->[STOP] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
			$param->[STOP] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
		my $raw_scores = $bw->get_values($chromo, $param->[STRT] - 1, $param->[STOP]);
		my @scores = grep {defined} @$raw_scores;
		return wantarray ? @scores : \@scores;
	}
	else {
		# we have multiple bigwigs
		my @scores;
		for (my $d = DATA; $d < scalar @$param; $d++) {
			my $bw = _get_bigwig($param->[$d]);
			my $chromo = $BIG_CHROMOS{$param->[$d]}{$param->[CHR]} or next;
			$param->[STRT] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
				$param->[STRT] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
			$param->[STOP] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
				$param->[STOP] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
			my $raw = $bw->get_values($chromo, $param->[STRT] - 1, $param->[STOP]);
			push @scores, grep {defined} @$raw;
		}
		return wantarray ? @scores : \@scores;
	}
}


sub collect_bigwig_position_scores {
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset(s)
	my $param = shift;
	my %pos2score;
	
	# check how many features we have
	if (scalar @$param == 9) {
		# only one, great!
		my $bw = _get_bigwig($param->[DATA]);
		my $chromo = $BIG_CHROMOS{$param->[DATA]}{$param->[CHR]} or return;
		$param->[STRT] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
			$param->[STRT] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
		$param->[STOP] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
			$param->[STOP] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
		my $intervals = $bw->get_intervals($chromo, $param->[STRT] - 1, $param->[STOP]);
		
		# record intervals into hash
		foreach my $i (@$intervals) {
			for (my $p = $i->{start} + 1; $p <= $i->{end}; $p++) {
				$pos2score{$p} = $i->{value};
			}
		}
	}
	else {
		# we have multiple bigwigs
		my %duplicates; # hash of duplicate positions, position => number
		
		# collect from each one
		for (my $d = DATA; $d < scalar @$param; $d++) {
			my $bw = _get_bigwig($param->[$d]);
			my $chromo = $BIG_CHROMOS{$param->[$d]}{$param->[CHR]} or next;
			$param->[STRT] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
				$param->[STRT] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
			$param->[STOP] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
				$param->[STOP] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
			my $intervals = $bw->get_intervals($chromo, $param->[STRT] - 1, $param->[STOP]);
		
			# record intervals into hash
			foreach my $i (@$intervals) {
				for (my $p = $i->{start} + 1; $p <= $i->{end}; $p++) {
					# check every position to see if it's a duplicate
					if (exists $pos2score{$p} ) {
						if (exists $duplicates{$p} ) {
							# append an incrementing number at the end
							$duplicates{$p} += 1; # increment first
							my $new = sprintf("%d.%d", $p, $duplicates{$p});
							$pos2score{$new} = $i->{value};
						}
						else {
							# first time duplicate
							my $new = $p . '.1';
							$pos2score{$new} = $i->{value};
							$duplicates{$p} = 1;
						}
					}
					else {
						$pos2score{$p} = $i->{value};
					}
				}
			}
		}
		
		# check for duplicate positions - we may not have any
		if (%duplicates) {
			_remove_duplicate_positions(\%pos2score, \%duplicates);
		}
	}
	return wantarray ? %pos2score : \%pos2score;
}





#### BigBed Subroutines

sub open_bigbed_db {
	my $path = shift;
	my $bb = _open_big($path) or 
		croak " Unable to open bigBed file '$path'! $@\n";
	unless ($bb->is_big_bed) {
		croak " $path is not a bigBed file!\n";
	}
	return $bb;
}


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
		my $bb = _get_bigbed($param->[$d]);
			
		# first check chromosome is present
		my $chromo = $BIG_CHROMOS{$param->[$d]}{$param->[CHR]} or next;
		$param->[STRT] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
			$param->[STRT] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
		$param->[STOP] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
			$param->[STOP] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
		
		# collect the features overlapping the region
			# we are using the high level API rather than the low-level
			# since we getting the individual scores from each bed element
		my $bb_stream = Bio::ToolBox::db_helper::big::BedIteratorWrapper->new(
			$bb, $chromo, $param->[STRT], $param->[STOP]);
		
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


sub collect_bigbed_position_scores {
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	
	# look at each bedfile
	# usually there is only one, but there may be more
	my %pos2data;
	for (my $i = DATA; $i < scalar @$param; $i++) {
	
		# open the bedfile
		my $bb = _get_bigbed($param->[$i]);
			
		# first check that the chromosome is present
		my $chromo = $BIG_CHROMOS{$param->[$i]}{$param->[CHR]} or next;
		$param->[STRT] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
			$param->[STRT] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
		$param->[STOP] = $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo} if 
			$param->[STOP] > $BIG_CHROMOLENGTHS{$param->[DATA]}{$chromo};
		
		# collect the features overlapping the region
		my $bb_stream = Bio::ToolBox::db_helper::big::BedIteratorWrapper->new(
			$bb, $chromo, $param->[STRT], $param->[STOP]);
		
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
			$pos2data{$position} = sum( @{$pos2data{$position}} ) / 
									scalar( @{$pos2data{$position}} );
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
			$pos2data{$position} = sum( @{$pos2data{$position}} ) / 
									scalar( @{$pos2data{$position}} );
		}
	}
	
	# return collected data
	return wantarray ? %pos2data : \%pos2data;
}

sub sum_total_bigbed_features {
	# there is no easy way to do this with this adapter, except to literally 
	# walk through the entire file.
	# well, we do this with bam files, I guess we could do the same here
	# honestly, who uses this????? it's legacy. skip for now until someone complains
	return undef;
}



#### BigWigSet Subroutines

sub open_bigwigset_db {
	my $path = shift;
	return Bio::ToolBox::db_helper::big::BigWigSet->new($path);
}


sub collect_bigwigset_score {
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	
	# lookup the bigWig files based on the parameters
	my $ids = _lookup_bigwigset_wigs($param);
	return unless scalar(@$ids) > 0;
	croak("multiple selected bigWig files from a BigWigSet is not supported with single score method")
		if scalar(@$ids) > 1;
	push @$param, @$ids;
	
	# use the low level single bigWig API 
	return collect_bigwig_score($param);
}


sub collect_bigwigset_scores {
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	
	# lookup the bigWig files based on the parameters
	my $ids = _lookup_bigwigset_wigs($param);
	return unless scalar(@$ids) > 0;
	push @$param, @$ids;
	
	# use the low level single bigWig API 
	return collect_bigwig_scores($param);
}


sub collect_bigwigset_position_scores {
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	
	# lookup the bigWig files based on the parameters
	my $ids = _lookup_bigwigset_wigs($param);
	return unless scalar(@$ids) > 0;
	push @$param, @$ids;
	
	# use the low level single bigWig API 
	return collect_bigwig_position_scores($param);
}




#### Internal

sub _open_big {
	my $path = shift;
	$path =~ s/^file://; # clean up file prefix if present
	my $big;
	eval {
		$big = Bio::DB::Big->open($path);
	};
	return $big if $big;
	return;	
}

sub _get_bigwig {
	my $file = shift;
	return $OPENED_BIG{$file} if exists $OPENED_BIG{$file};
	
	# open and cache the bigFile object
	my $bw = _open_big($file) or 
		croak " Unable to open big file '$file'! $@\n";
	unless ($bw->is_big_wig) {
		croak " $file is not a bigWig file!\n";
	}
	$OPENED_BIG{$file} = $bw;
	_record_seqids($file, $bw);
	return $bw;
}

sub _get_bigbed {
	my $file = shift;
	return $OPENED_BIG{$file} if exists $OPENED_BIG{$file};
	
	# open and cache the bigFile object
	my $bb = _open_big($file) or 
		croak " Unable to open big file '$file'! $@\n";
	unless ($bb->is_big_bed) {
		croak " $file is not a bigBed file!\n";
	}
	$OPENED_BIG{$file} = $bb;
	_record_seqids($file, $bb);
	return $bb;
}

sub _record_seqids {
	my ($file, $big) = @_;
	$BIG_CHROMOS{$file} = {};
	my $chroms = $big->chroms(); # returns hash seq_id => length
	foreach my $c (keys %$chroms) {
		my $chr = $chroms->{$c}{name};
		my $len = $chroms->{$c}{length};
		# store length
		$BIG_CHROMOLENGTHS{$file}{$chr} = $len;
		$BIG_CHROMOS{$file}{$chr} = $chr; # itself
		# now store alternate names, with or without chr prefix
		if ($chr =~ /^chr(.+)$/i) {
			$BIG_CHROMOS{$file}{$1} = $chr;
		}
		else {
			$BIG_CHROMOS{$file}{"chr$chr"} = $chr;
		}
	}
}

sub _remove_duplicate_positions {
	
	# collect the feature and hashes
	my ($pos2data, $duplicates) = @_;
	
	# remove the duplicates
		# we will combine all of them with a simple mean, what else to do?
	foreach my $pos (keys %{ $duplicates } ) {
		my $num = $duplicates->{$pos};
		
		# collect all the values
		my @values;
		push @values, $pos2data->{$pos}; # initial value
		for my $i (1..$num) {
			push @values, $pos2data->{ "$pos\.$i" };
			delete $pos2data->{ "$pos\.$i" };
		}
		$pos2data->{$pos} = sum(@values) / scalar(@values);
	}
}

sub _lookup_bigwigset_wigs {
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	# the datasets, could be either types or names, unfortunately
	my @types = splice(@$param, DATA);
	
	# we cache the list of looked up bigwigs to avoid doing this over and over
	my $lookup = sprintf("%s_%s_%s", join('_', @types), $param->[STND], 
		$param->[STR]);
	return $BIGWIGSET_WIGS{$lookup} if exists $BIGWIGSET_WIGS{$lookup};
	
	# filter first by the name or type
	my @ids = $param->[DB]->filter_bigwigs(@types);
	
	# then check for strand if necessary
	if ($param->[STND] ne 'all' and $param->[STR] != 0) {
    	# looks like we are collecting stranded data
    	# try to filter again based on strand attribute
    	my $strand_to_keep;
    	if ($param->[STND] eq 'sense') {
			$strand_to_keep = $param->[STR];
		}
		elsif ($param->[STND] eq 'antisense') {
			$strand_to_keep =   $param->[STR] == -1 ? 1 :
								$param->[STR] == 1 ? -1 : 0;
		}
		else {
			confess sprintf "bad strandedness value: %s", $param->[STND];
		}
		@ids = $param->[DB]->filter_bigwigs_by_strand($strand_to_keep, @ids);
    }
	
	# map the names back to full bigwig paths
	my @paths = map { $param->[DB]->get_bigwig_path($_) } @ids;
	
	# cache and return
	$BIGWIGSET_WIGS{$lookup} = \@paths;
	return \@paths;
}



package Bio::ToolBox::db_helper::big::BigWigSet;
# this package borrows concepts from Bio::DB::BigWigSet by Lincoln Stein
# in order to make it compatible with Bio::DB::Big
# it is NOT a drop-in replacement or equivalent
# if you need more functionality, install Bio::DB::BigWigSet
use Carp;
use IO::Dir;
use IO::File;
use File::Spec;

sub new {
	my ($class, $dir) = @_;
	croak "must call method new with a directory path!" unless $dir;
	croak "BigWigSet '$dir' is not a directory path!" unless -d $dir;
	my $self = {
		dir      => $dir,
		bwfiles  => {},
		metadata => {}
	};
	
	# read directory
	my $D = IO::Dir->new($dir) or croak "unable to open $dir! $!";
	while (my $f = $D->read) {
		my $f_path = File::Spec->catfile($dir, $f);
		next unless -f $f_path;
		if ($f =~ /\.(?:bw|bigwig)$/i) {
			# bigwig file
			$self->{bwfiles}{$f} = $f_path;
		}
		elsif ($f =~ /^meta.*\.txt$/i) {
			# a genuine metadata file - excellent!
			# let's open it and parse the contents
			my $fh = IO::File->new($f_path) or 
				croak "unable to read metadata file '$f_path'! $!";
			my $current_bw; # the current bigwig file we're describing
			while (my $line = $fh->getline) {
				next if $line =~ /^#/;
				next if $line !~ /\w+/;
				chomp $line;
				if ($line =~ /\[(.+\.(?:bw|bigwig))\]/i) {
					# start of a new metadata block for the current bw file
					$current_bw = $1;
					next;
				}
				elsif ($line =~ /^([\w\-\.]+)\s*=\s*([\w\-\.]+)/) {
					# a metadata line
					croak "malformed metadata file! no bw file stanza header before metadata line!"
						unless defined $current_bw;
					$self->{metadata}{$current_bw}{$1} = $2;
				} 
				# ignore everything else
			}
			$fh->close;
		}
		# ignore all other files
	}
	undef $D;
	
	# fill out the metadata
	foreach my $f (keys %{$self->{bwfiles}}) {
		# a genuine metadata file isn't absolutely required
		# therefore we will always put in our own metadata gleaned from file name 
		# do a regex to pull out basename and possibly strand information from filename
		my $name = $f;
		$name =~ s/\.(?:bw|bigwig)$//i; # remove extension for the name
		my ($type, $strand);
		if ($f =~ /^(.+)[\._](?:f|for|forward|plus|\+)\.(?:bw|bigwig)$/i) {
			$type = $1;
			$strand = 1;
		}
		elsif ($f =~ /^(.+)[\._](?:r|rev|reverse|minus|\-)\.(?:bw|bigwig)$/i) {	
			$type = $1;
			$strand = -1;
		}
		else {
			$type = $name;
			$strand = 0;
		}
		# be careful and do not overwrite pre_existing key value
		$self->{metadata}{$f}{name}   ||= $name;
		$self->{metadata}{$f}{type}   ||= $type;
		$self->{metadata}{$f}{strand} ||= $strand;
	}
	
	bless $self, $class;
	return $self;
}

sub bigwig_names {
	my $self = shift;
	my @b = keys %{ $self->{bwfiles} };
	return wantarray ? @b : \@b;
}

sub bigwigs {
	# to maintain compatiblity with the old BigWigSet, we will return the full path
	# of all the bigwig files
	my $self = shift;
	my @b = values %{ $self->{bwfiles} };
	return wantarray ? @b : \@b;
}

sub get_bigwig {
	my ($self, $bwfile) = @_;
	return unless $bwfile;
	my $path = $self->get_bigwig_path($bwfile);
	if ($path) {
		return Bio::ToolBox::db_helper::big::open_bigwig_db($path);
	}
	else {
		carp "unrecognized bigWig file name '$bwfile'!\n";
		return;
	}
}

sub get_bigwig_path {
	my ($self, $bwfile) = @_;
	return unless $bwfile;
	return $self->{bwfiles}{$bwfile} || undef;
}

sub metadata {
	my $self = shift;
	return $self->{metadata};
}

sub filter_bigwigs {
	my ($self, @names) = @_;
	
	# we start with all of the bigwigs available, then filter out
	my $start_list = $self->bigwig_names;
	my $md = $self->metadata;
	my @filtered;
	
	foreach my $name (@names) {
		# there may be more than one name of a dataset passed
		# in all likelihood, probably not, but just in case....
		
		# we are filtering on basically any value that matches and are not really 
		# restricting on a specific attribute key
		# so type, name, display_name, or primary_tag are all checked in that order
		# let's hope this is adequate and won't create too many problems
		# I'm betting on the fact that bigwigset databases are so obscure that 
		# most users won't even bother with a genuine metadata file
		# in fact, why am I even bothering with this at all?
		# because they are a cool concept and I occasionally use them. huh.
		foreach my $b (@$start_list) {
			if (exists $md->{$b}{type}) {
				push @filtered, $b if $name eq $md->{$b}{type};
			} 
			elsif (exists $md->{$b}{name}) {
				push @filtered, $b if $name eq $md->{$b}{name};
			} 
			elsif (exists $md->{$b}{display_name}) {
				push @filtered, $b if $name eq $md->{$b}{display_name};
			} 
			elsif (exists $md->{$b}{primary_tag}) {
				push @filtered, $b if $name eq $md->{$b}{primary_tag};
			} 
		}
	}
	return wantarray ? @filtered : \@filtered;
}

sub filter_bigwigs_by_strand {
	my $self = shift;
	my $strand = shift;
	
	# check what files were passed
	my @files;
	if (scalar(@_) == 1 and ref($files[0]) eq 'ARRAY') {
		@files = @$_[0];
	}
	elsif (scalar(@_) == 0) {
		@files = $self->bigwig_names;
	}
	else {
		@files = @_;
	}
	
	# filter the bigWigs based on strand
	# we keep anything that matches the given strand
	# because the old BigWigSet adapter didn't really handle strand properly, 
	# and certainly not standard BioPerl 1,0,-1 values, I had adopted plus, none, minus
	# now I'm doomed as I have to support this too. ugh.
	my @keepers;
	my $md = $self->metadata;
	if ($strand >= 0) {
		foreach my $f (@files) {
			my $str = $md->{$f}{strand} || 0;
			if ($str =~ /[a-z]/) {
				# ugh old school
				$str =  $str eq 'plus' ? 1 :
						$str eq 'minus' ? -1 :
						$str eq 'none' ? 0 : 0;
			} 
			push @keepers, $f if $str >= 0;
		}
	}
	elsif ($strand < 0) {
		foreach my $f (@files) {
			my $str = $md->{$f}{strand} || 0;
			if ($str =~ /[a-z]/) {
				# ugh old school
				$str =  $str eq 'plus' ? 1 :
						$str eq 'minus' ? -1 :
						$str eq 'none' ? 0 : 0;
			} 
			push @keepers, $f if $str < 0;
		}
	}
	return @keepers;
}


package Bio::ToolBox::db_helper::big::BedIteratorWrapper;
use Carp;
use Bio::ToolBox::SeqFeature;

sub new {
	my $class = shift;
	
	# check the bigBed object
	my $bb = shift;
	unless (ref($bb) eq 'Bio::DB::Big::File' and $bb->is_big_bed) {
		confess "passed big object is not a bigBed file!";
	}
	
	# get coordinates
	my ($seqid, $start, $end) = @_;
	confess "no coordinates!" unless (defined $seqid and defined $start and defined $end);
	$start -= 1; # compensate for 0-based coordinates
	
	# create an iterator
	# include all string entries and just one per iteration
	# return 10 at a time, not too many, not too little
	my $iterator = $bb->get_entries_iterator($seqid, $start, $end, 1, 10);
	
	# return object wrapper
	my $self = {
		bb      => $bb,
		iter    => $iterator,
		seqid   => $seqid,
		start   => $start,
		end     => $end,
		entries => [],
	};
	return bless $self, $class;
}

sub next_seq {
	my $self = shift;
	my $entry = shift @{ $self->{entries} } || undef;
	unless ($entry) {
		$self->{entries} = $self->{iter}->next || undef;
		return unless $self->{entries};
		$entry = shift @{ $self->{entries} } || undef;
	}
	return unless $entry;
	my @bits = split('\t', $entry->{string});
	return Bio::ToolBox::SeqFeature->new(
		-seq_id         => $self->{seqid},
		-start          => $entry->{start} + 1, # compensate for 0-based coordinates
		-end            => $entry->{end},
		-display_name   => $bits[0] || undef,
		-score          => $bits[1] || 1,
		-strand         => $bits[2] || 0,
	);
}


__END__

=head1 NAME

Bio::ToolBox::db_helper::big

=head1 DESCRIPTION

This module provides support for binary BigWig and BigBed files to 
the L<Bio::ToolBox> package. It also provides minimal support for a 
directory of one or more bigWig files as a combined database, known as a 
BigWigSet. 

=head1 USAGE

The module requires L<Bio::DB::Big> to be installed, which in turn 
requires the L<libBigWig | https://github.com/dpryan79/libBigWig> C 
library to be installed. This provides a simpler and easier-to-install 
library compared to the UCSC Kent C libraries.

In general, this module should not be used directly. Use the methods 
available in L<Bio::ToolBox::Data> or L<Bio::ToolBox::db_helper>.  

All subroutines are exported by default.

=head2 Available subroutines

=over

=item collect_bigwig_score()

This subroutine will collect a single value from a binary bigWig file. 
It uses the low-level summary method to collect the statistical 
information and is therefore significantly faster than the other 
methods, which rely upon parsing individual data points across the 
region.

The subroutine is passed a parameter array reference. See below for details.

The object will return either a valid score or a null value.

=item collect_bigwigset_score()

Similar to collect_bigwig_score() but using a BigWigSet database of 
BigWig files. Unlike individual BigWig files, BigWigSet features support 
stranded data collection if a strand attribute is defined in the metadata 
file. 

The subroutine is passed a parameter array reference. See below for details.
    
=item collect_bigwig_scores()

This subroutine will collect only the score values from a binary BigWig file 
for the specified database region. The positional information of the 
scores is not retained.

The subroutine is passed a parameter array reference. See below for details.

The subroutine returns an array or array reference of the requested dataset 
values found within the region of interest. 

=item collect_bigwigset_scores()

Similar to collect_bigwig_scores() but using a BigWigSet database of 
BigWig files. Unlike individual BigWig files, BigWigSet features support 
stranded data collection if a strand attribute is defined in the metadata 
file. 

The subroutine is passed a parameter array reference. See below for details.

=item collect_bigwig_position_scores()

This subroutine will collect the score values from a binary BigWig file 
for the specified database region keyed by position. 

The subroutine is passed a parameter array reference. See below for details.

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. Note that only one value is 
returned per position, regardless of the number of dataset features 
passed. Usually this isn't a problem as only one dataset is examined at a 
time.

=item collect_bigwigset_position_scores()

Similar to collect_bigwig_position_scores() but using a BigWigSet database 
of BigWig files. Unlike individual BigWig files, BigWigSet features support 
stranded data collection if a strand attribute is defined in the metadata 
file. 

The subroutine is passed a parameter array reference. See below for details.

=item open_bigwig_db()

This subroutine will open a BigWig database connection. Pass either the 
local path to a bigWig file (.bw extension) or the URL of a remote bigWig 
file. It will return the opened database object.

=item open_bigwigset_db()

This subroutine will open a BigWigSet database connection using a directory 
of BigWig files and one metadata index file, as described in 
Bio::DB::BigWigSet. Essentially, this treats a directory of BigWig files as 
a single database with each BigWig file representing a different feature 
with unique attributes (type, source, strand, etc). 

Pass the subroutine a scalar value representing the local path to the 
directory. It presumes a feature_type of 'region', as expected by the other 
Bio::ToolBox::db_helper subroutines and modules. It will return the opened database 
object.

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

Acceptable values include mean, min, and max when collecting 
single score over a genomic segment. This uses the built-in 
statistic zoom levels of the bigWig.

=item 7. A database object.

Pass the opened L<Bio::DB::BigWigSet> database object when working 
with BigWigSets.

=item 8 and higher. BigWig file names or BigWigSet database types.

Opened BigWig objects are cached. Both local and remote BigWig files 
are supported. 

=back

=head1 SUPPORT MODULES

This includes two additional object-oriented modules for supporting 
BigWigSets and bigBed SeqFeature iteration. 

=head2 Bio::ToolBox::db_helper::big::BigWigSet

This provides support for a BigWigSet, which is not natively supported 
by the L<Bio::DB::Big> adapter, and is based on the concepts from the 
L<Bio::DB::BigWigSet> adapter. However, it is B<NOT> a drop-in replacement, 
only a few methods are provided, and only a few of these are similar to 
the original adapter. 

This adapter will still read a C<INI>-style F<metadata.txt> 
file as described in L<Bio::DB::BigWigSet> for metadata. Briefly, this 
file format is similar to below

    [file1.bw]
    name = mydata
    type = ChIPSeq
    
    [file2.bw]
    name = mydata2
    type = ChIPSeq
    
Each bigWig file in the directory should have a stanza entry with the 
path and file name in the stanza header in square brackets. Metadata 
is included as simple key = value pairs, where keys can be typical 
SeqFeature attributes, including C<display_name> or C<name>, 
C<primary_tag> or C<type>, and C<strand>.

B<NOTE:> Metadata text files are ideal, but not required. If a 
metadata file is not present, appropriate metadata will be determined 
from the bigWig file names, using the basename as the metadata C<name> 
and possibly extracting the strand from the end of the filename, if 
it ends in a F<_f> or F<_r>. 

The following methods are available.

=over 4

=item new

Generate a new BigWigSet object. The path of the directory must be 
passed as an argument. The contents of the directory will be read, 
bigWig files located, metadata files (if any) read and processed. 
Fasta sequences are not supported.

=item bigwig_names

Returns an array or array reference of the bigWig file names. These 
are just the file names, without the path.

=item bigwigs

Returns an array or array reference of the full path for the bigWig 
files. This is identical to the L<Bio::DB::BigWigSet> method.

=item get_bigwig

Given a bigwig name, this will return an opened bigwig database 
L<Bio::DB::Big> object. This is identical to the L<Bio::DB::BigWigSet> 
method.

=item get_bigwig_path

Given a bigwig name, this will return the full path to the corresponding 
bigWig file.

=item metadata

This will return a hash reference pointing to the metadata hash structure, 
the keys of which are bigWig names, and the values are hash references 
for the metadata key = value metadata pairs. This is identical to the 
L<Bio::DB::BigWigSet> method.

=item filter_bigwigs

Provide one or more names, primary tags, or types to filter the bigWig 
files in the Set. For the purposes of this simple method, no distinction is 
made whether the filtering criteria is a C<display_name>, C<primary_tag>, 
or C<type>. The provided text strings will be used to search all the 
metadata values, and the names of files with exact matches are returned.
For purposes of filtering, the following metadata keys are searched in the 
following order: C<type>, C<name>, C<display_name>, and C<primary_tag>, and 
the first match is kept.

=item filter_bigwigs_by_strand

Pass first the strand, and then optionally a list of bigWig names, 
perhaps the results from filter_bigwigs(). If no names were passed, all 
the names in the BigWigSet will be considered. The names of files whose 
strand matches the given strand will be returned.

=back

=head2 Bio::ToolBox::db_helper::big::BedIteratorWrapper

This is an object wrapper around a bigBed database for retrieving items 
and returning them as convenient SeqFeature objects. 
Only the first 3 to 6 standard BED columns are supported: seq_id, 
start, stop, name, score, and strand. 

The following methods are provided.

=over 4

=item new

Pass the new method the following items: opened L<Bio::DB::Big> bigBed object, 
chromosome, start, and end coordinates. 

=item next_seq

This will return the next available feature in the established search interval 
as a L<Bio::ToolBox::SeqFeature> object. The method name is consistent with 
other L<Bio::Perl> compatible objects and iterators. Yeah, it sucks, and not 
very apropos to the actual function. Oh well. 

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

