package Bio::ToolBox::db_helper::bigwig;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(min max mean);
use Bio::DB::BigWig qw(binMean binStdev);
use Bio::DB::BigWigSet;
use constant {
	CHR  => 0,  # chromosome
	STRT => 1,  # start
	STOP => 2,  # stop
	STR  => 3,  # strand
	STND => 4,  # strandedness
	METH => 5,  # method
	VAL  => 6,  # value type
	DB   => 7,  # database object
	DATA => 8,  # first dataset, additional may be present
};
our $VERSION = '1.50';


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_bigwig_score
	collect_bigwig_scores
	collect_bigwig_position_scores
	collect_bigwigset_score
	collect_bigwigset_scores
	collect_bigwigset_position_scores
	open_bigwig_db
	open_bigwigset_db
);


# Hash of Bigfile chromosomes
our %BIGWIG_CHROMOS;
	# sometimes user may request a chromosome that's not in the bigfile
	# that could lead to an exception
	# we will record the chromosomes list in this hash
	# $BIGWIG_CHROMOS{bigfile}{chromos}
	# we also record the chromosome name variant with or without chr prefix
	# to accommodate different naming conventions

# Opened bigWig db objects
our %OPENED_BW;
	# a cache for opened BigWig databases, primarily for collecting scores
	# caching here is only for local purposes of collecting scores
	# db_helper also provides caching of db objects but with option to force open in
	# the case of forking processes - we don't have that here

# Methods Cache lookup
	# pre-define the methods code for processing summary objects to avoid 
	# excessive if-elsif tests every data cycle
our %ONE_SUMMARY_METHOD = (
	'mean'  => sub {return binMean(shift) },
	'sum'   => sub {return shift->{sumData} },
	'min'   => sub {return shift->{minVal} },
	'max'   => sub {return shift->{maxVal} },
	'count' => sub {return shift->{validCount} },
	'stddev'=> sub {return binStdev(shift) },
);
our %MULTI_SUMMARY_METHOD = (
	'mean'  => sub {
				my $summaries = shift;
				my $sum = 0;
				my $count = 0;
				foreach my $s (@$summaries) {
					$sum += $s->{sumData};
					$count += $s->{validCount};
				}
				return $count ? $sum/$count : '.';
			},
	'sum'   => sub {
				my $summaries = shift;
				my $sum = 0;
				foreach my $s (@$summaries) {
					$sum += $s->{sumData};
				}
				return $sum;
			},
	'min'   => sub {
				my $summaries = shift;
				return min( map {$_->{minVal}} @$summaries);
			},
	'max'   => sub {
				my $summaries = shift;
				return max( map {$_->{maxVal}} @$summaries);
			},
	'count' => sub {
				my $summaries = shift;
				my $count = 0;
				foreach my $s (@$summaries) {
					$count += $s->{validCount};
				}
				return $count;
			},
	'stddev'=> sub {croak "cannot calculate stddev value from multiple bigWig Summary objects!"},
);

# The true statement
1; 



### Modules ###

### Collect single BigWig score
sub collect_bigwig_score {
	
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, value, db, dataset
	my $param = shift;
	
	# Collecting summary features
	# we will collect a summary object for each requested wig feature  
	
	# Walk through each requested feature
	# There is likely only one
	my @summaries;
	for (my $d = DATA; $d < scalar @$param; $d++) {
		my $bw = _get_bw($param->[$d]);
		
		# first check that the chromosome is present
		my $chromo = $BIGWIG_CHROMOS{$param->[$d]}{$param->[CHR]} or next;
		
		# use a low level method to get a single summary hash for 1 bin
		my $sumArrays = $bw->bigWigSummaryArrayExtended(
			# chromo, 0-based start, stop, bin
			$chromo, $param->[STRT] - 1, $param->[STOP], 1
		);
		push @summaries, $sumArrays->[0];
	}
	
	# now process the summary features
	return _process_summaries($param->[METH], \@summaries);
}


### Collect multiple BigWig scores
sub collect_bigwig_scores {
	
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, value, db, dataset
	my $param = shift;
	
	# Walk through each requested feature
	# There is likely only one
	my @scores; 
	for (my $d = DATA; $d < scalar @$param; $d++) {
	
		# open db object
		my $bw = _get_bw($param->[$d]);
		
		# first check that the chromosome is present
		my $chromo = $BIGWIG_CHROMOS{$param->[$d]}{$param->[CHR]} or next;
		
		# initialize low level stream for this segment
		my $list = $bw->bigWigIntervalQuery(
			$chromo, $param->[STRT] - 1, $param->[STOP] );
		my $f = $list->head;
		while ($f) {
			# print $f->start,"..",$f->end,": ",$f->value,"\n";
			push @scores, $f->value; 
			$f = $f->next;
		}
	} 
	
	# return collected data
	return wantarray ? @scores : \@scores;
}


### Collect positioned BigWig scores
sub collect_bigwig_position_scores {
	
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, value, db, dataset
	my $param = shift;
	
	# initialize 
	my %pos2data; # hash of position => score
	my %duplicates; # hash of duplicate positions, position => number
	
	# Walk through each requested feature
	# There is likely only one
	for (my $d = DATA; $d < scalar @$param; $d++) {
	
		# Open the BigWig file
		my $bw = _get_bw($param->[$d]);
		
		# first check that the chromosome is present
		my $chromo = $BIGWIG_CHROMOS{$param->[$d]}{$param->[CHR]} or next;
		
		# initialize low level stream for this segment
		my $list = $bw->bigWigIntervalQuery(
			$chromo, $param->[STRT] - 1, $param->[STOP] );
		my $f = $list->head;
		while ($f) {
			# print $f->start,"..",$f->end,": ",$f->value,"\n";
			# process the feature
			for (my $pos = $f->start + 1; $pos <= $f->end; $pos++) {
				# remember low level features are 0-based
				# check for duplicate positions
				# duplicates should not normally exist for a single wig file
				# but we may be working with multiple wigs that are being combined
				if (exists $pos2data{$pos} ) {
					if (exists $duplicates{$pos} ) {
						# append an incrementing number at the end
						$duplicates{$pos} += 1; # increment first
						my $new = sprintf("%d.%d", $pos, $duplicates{$pos});
						$pos2data{$new} = $f->value;
					}
					else {
						# first time duplicate
						my $new = $pos . '.1';
						$pos2data{$new} = $f->value;
						$duplicates{$pos} = 1;
					}
				}
				else {
					$pos2data{$pos} = $f->value;
				}
			}
			$f = $f->next;
		}
	} 
	
	# check for duplicate positions
	if (%duplicates) {
		_remove_duplicate_positions(\%pos2data, \%duplicates);
	}
	
	# return collected data
	return wantarray ? %pos2data : \%pos2data;
}



sub collect_bigwigset_score {
	
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, value, db, dataset
	my $param = shift;
	
	# Confirm the chromosome
	my $chromo = $BIGWIG_CHROMOS{ ($param->[DB]->bigwigs)[0] }{$param->[CHR]} 
		|| undef;
	unless (defined $chromo) {
		# chromosome is not present, at least in the first bigwig in the set
		# return null
		return $param->[METH] =~ /sum|count/ ? 0 : '.';
	}
	
	# Reset which feature_type to collect from the database
	# this is normally set to region when we opened the bigwigset db
	# but it may be changed by a previous method
	# we now want summary feature_type
	$param->[DB]->feature_type('summary');
	
	# Collecting summary features
	# we will collect a summary object for each requested wig feature  
	my @summaries;
	
	# Work through all the requested feature types
	# the BigWigSet feature request doesn't work well with multiple features
	# so we'll do them one at a time
	for (my $d = DATA; $d < scalar @$param; $d++) {
		my $type = $param->[$d];
		
		# Collect the summary features
		# no low-level joy here, use the high level API
		$type =~ s/\:.+$//; # strip the source if present
			# features method only works with primary_tag, not full 
			# primary_tag:source type
		my @features = $param->[DB]->features(
			-seq_id   => $chromo,
			-start    => $param->[STRT],
			-end      => $param->[STOP],
			-type     => $type,
		);
		# if the type doesn't work, then try display_name instead
		unless (@features) {
			@features = $param->[DB]->features(
				-seq_id   => $chromo,
				-start    => $param->[STRT],
				-end      => $param->[STOP],
				-name     => $type,
			);
		}
			# since we're collecting summary features, we will only get one 
			# per bigwig file that matches the request
			# no need for a seqfeature stream
		
		# Determine which features to take based on strandedness
		
		# Stranded features
		if (
			$param->[STR] != 0 and
			($param->[STND] eq 'sense' or $param->[STND] eq 'antisense')
		) {
			# we will have to check the strand for each object
			# feature objects we collect don't have the standard strand set
			# instead, we will have to get the attribute tag named strand
			
			# check each feature
			foreach my $f (@features) {
				
				# get the feature strand
				my $fstrand = 0; # default
				if ($f->has_tag('strand') ) {
					($fstrand) = $f->get_tag_values('strand');
				}
				
				# collect summary if strand is appropriate
				if (
					$fstrand == 0 or
					(
						# sense data
						$param->[STR] == $fstrand 
						and $param->[STND] eq 'sense'
					) 
					or (
						# antisense data
						$param->[STR] != $fstrand  
						and $param->[STND] eq 'antisense'
					)
				) {
					# we have acceptable data to collect
					push @summaries, $f->score;
				}
			}
		}
		
		# Non-stranded features
		else {
			# take all the features found
			
			# keep all the summaries
			foreach my $f (@features) {
				push @summaries, $f->score;
			}
		}
	}
	
	# now process the summary features
	return _process_summaries($param->[METH], \@summaries);
}




sub collect_bigwigset_scores {
	
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, value, db, dataset
	my $param = shift;
	
	# Confirm the chromosome
	my $chromo = $BIGWIG_CHROMOS{ ($param->[DB]->bigwigs)[0] }{$param->[CHR]}; 
	return unless (defined $chromo);
		# chromosome is not present, at least in the first bigwig in the set
		# return nothing
	
	# Reset which feature_type to collect from the database
	# this is normally set to region when we opened the bigwigset db
	# but it may be changed by a previous method
	# we now want region feature_type
	$param->[DB]->feature_type('region');
	
	# initialize collection array
	my @scores; 
	
	# Go through each feature type requested
	for (my $i = 8; $i < scalar @$param; $i++) {
		my $type = $param->[$i];
		
		# Collect the feature scores
		# no low-level joy here, use the high level API
		$type =~ s/\:.+$//; # strip the source if present
			# features method only works with primary_tag, not full 
			# primary_tag:source type
			# since the default feature_type for the bigwigset database is 
			# region, we will get lots of features returned, one for each 
			# datapoint
			# use an iterator to process them
		my $iterator = $param->[DB]->get_seq_stream(
			-seq_id   => $chromo,
			-start    => $param->[STRT],
			-end      => $param->[STOP],
			-type     => $type,
		);
		
		# check that we have a feature stream
		my $feature = $iterator->next_seq;
		unless ($feature) {
			# uh oh, no feature! perhaps we didn't correctly identify the 
			# correct bigwig in the set
			# try again using display_name instead
			$iterator = $param->[DB]->get_seq_stream(
				-seq_id   => $chromo,
				-start    => $param->[STRT],
				-end      => $param->[STOP],
				-name     => $type,
			);
			$feature = $iterator->next_seq;
		}
		
		# Determine which features to take based on strandedness
		
		# Stranded features
		if (
			$param->[STR] != 0 and
			($param->[STND] eq 'sense' or $param->[STND] eq 'antisense')
		) {
			# we will have to check the strand for each object
			# feature objects we collect don't have the standard strand set
			# instead, we will have to get the attribute tag named strand
			
			# check each feature
			while ($feature) {
				
				# get the feature strand
				my $fstrand = 0; # default
				if ($feature->has_tag('strand') ) {
					($fstrand) = $feature->get_tag_values('strand');
				}
				
				# collect score if strand is appropriate
				if (
					$fstrand == 0 or
					(
						# sense data
						$param->[STR] == $fstrand 
						and $param->[STND] eq 'sense'
					) 
					or (
						# antisense data
						$param->[STR] != $fstrand  
						and $param->[STND] eq 'antisense'
					)
				) {
					# we have acceptable data to collect
					push @scores, $feature->score;
				}
				
				# prepare for the next feature
				$feature = $iterator->next_seq;
			}
		}
		
		# Non-stranded features
		else {
			# take all the features found
			while ($feature) {
				push @scores, $feature->score;
				$feature = $iterator->next_seq;
			}
		}
	}
	
	# Finished
	return wantarray ? @scores : \@scores;;
}



sub collect_bigwigset_position_scores {
	
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, value, db, dataset
	my $param = shift;
	
	# Confirm the chromosome
	my $chromo = $BIGWIG_CHROMOS{ ($param->[DB]->bigwigs)[0] }{$param->[CHR]} 
		|| undef;
	return unless (defined $chromo);
		# chromosome is not present, at least in the first bigwig in the set
		# return nothing
	
	# Reset which feature_type to collect from the database
	# this is normally set to region when we opened the bigwigset db
	# but it may be changed by a previous method
	# we now want region feature_type
	$param->[DB]->feature_type('region');
	
	# initialize collection hash, position => score
	my %pos2data; 
	my %duplicates;
	
	# Go through each feature type requested
	for (my $i = 8; $i < scalar @$param; $i++) {
		my $type = $param->[$i];
		
		# Collect the feature scores
		$type =~ s/\:.+$//; # strip the source if present
			# features method only works with primary_tag, not full 
			# primary_tag:source type
			# since the default feature_type for the bigwigset database is 
			# region, we will get lots of features returned, one for each 
			# datapoint
			# use an iterator to process them
		# no low-level joy here, use the high level API
		my $iterator = $param->[DB]->get_seq_stream(
			-seq_id   => $chromo,
			-start    => $param->[STRT],
			-end      => $param->[STOP],
			-type     => $type,
		);
		
		# check that we have a feature stream
		my $feature = $iterator->next_seq if $iterator;
		unless ($feature) {
			# uh oh, no feature! perhaps we didn't correctly identify the 
			# correct bigwig in the set
			# try again using display_name instead
			$iterator = $param->[DB]->get_seq_stream(
				-seq_id   => $chromo,
				-start    => $param->[STRT],
				-end      => $param->[STOP],
				-name     => $type,
			);
			$feature = $iterator->next_seq;
		}
		
		# Determine which features to take based on strandedness
		
		# Stranded features
		if (
			$param->[STR] != 0 and
			($param->[STND] eq 'sense' or $param->[STND] eq 'antisense')
		) {
			# we will have to check the strand for each object
			# feature objects we collect don't have the standard strand set
			# instead, we will have to get the attribute tag named strand
			
			# Check each feature
			
			# Stranded features
			while ($feature) {
				
				# get the feature strand
				my $fstrand = 0; # default
				if ($feature->has_tag('strand') ) {
					($fstrand) = $feature->get_tag_values('strand');
				}
				
				# collect score if strand is appropriate
				if (
					$param->[STR] == 0 or
					(
						# sense data
						$param->[STR] == $fstrand 
						and $param->[STND] eq 'sense'
					) 
					or (
						# antisense data
						$param->[STR] != $fstrand  
						and $param->[STND] eq 'antisense'
					)
				) {
					# acceptable data point
					_process_position_score_feature($feature, \%pos2data, \%duplicates);
				}
				$feature = $iterator->next_seq;
			}
		}
		
		# Non-stranded features
		else {
			# take all the features found
			while ($feature) {
				# process
				_process_position_score_feature($feature, \%pos2data, \%duplicates);
				$feature = $iterator->next_seq;
			}
		}
	}
	
	# Remove duplicates
	if (%duplicates) {
		_remove_duplicate_positions(\%pos2data, \%duplicates);
	}
	
	
	# Finished
	return %pos2data;
}




### Open a bigWig database connection
sub open_bigwig_db {
	
	# path
	my $wigfile = shift;
	my $path = $wigfile;
	$path =~ s/^file://; # clean up file prefix if present
	
	# open
	my $bw;
	eval {
		$bw = Bio::DB::BigWig->new( -bigwig => $path);
	};
	return unless $bw;
	
	return $bw;
}




### Open a bigWigSet database connection
sub open_bigwigset_db {
	
	my $directory = shift;
	
	# check for trailing slash
	# this seems to interfere with generating the list of files, leading 
	# to duplicates: both raw files as well as contents from metadata
	$directory =~ s/\/$//; # strip trailing slash
	
	# open the database connection 
	# we're using the region feature type because that's what the rest of 
	# Bio::ToolBox::db_helper modules expect and work with
	my $bws;
	eval {
		$bws = Bio::DB::BigWigSet->new(
					-dir            => $directory,
					-feature_type   => 'region',
		);
	};
	return unless $bws;
	
	# check that we haven't just opened a new empty bigwigset object
	my @paths = $bws->bigwigs;
	return unless @paths; # no valid bigWig files, not valid
	
	# we have bigwig files, must be a valid bigwigset directory
	
	# collect the chromosomes from the first bigwig
	# we will assume all of the bigwigs have the same chromosomes!
	_record_seqids($paths[0], $bws->get_bigwig($paths[0]) );
	
	# check for potential implied strandedness
	my $md = $bws->metadata;
	foreach my $i (keys %$md) {
		my $f = $md->{$i}{'dbid'}; # the file path
		if ($f =~ /(?:f|for|forward|top|plus|\+)\.bw$/i) {
			# implied forward strand 
			unless (exists $md->{$i}{'strand'}) {
				$bws->set_bigwig_attributes($f, {'strand' => 1});
			}
		}
		elsif ($f =~ /(?:r|rev|reverse|bottom|minus|\-)\.bw$/i) {
			# implied reverse strand 
			unless (exists $md->{$i}{'strand'}) {
				$bws->set_bigwig_attributes($f, {'strand' => -1});
			}
		}
	}
	
	return $bws;
}


### Internal subroutine for getting the cached bigwig object
sub _get_bw {
	my $bwfile = shift;
	
	return $OPENED_BW{$bwfile} if exists $OPENED_BW{$bwfile};
	
	# open and cache the bigWig object
	my $bw = open_bigwig_db($bwfile) or 
		croak " Unable to open bigWig file '$bwfile'! $!\n";
	$OPENED_BW{$bwfile} = $bw;
	_record_seqids($bwfile, $bw);
	return $bw;
}

	# record the chromosomes and possible variants
sub _record_seqids {
	my ($bwfile, $bw) = @_;
	$BIGWIG_CHROMOS{$bwfile} = {};
	foreach my $s ($bw->seq_ids) {
		$BIGWIG_CHROMOS{$bwfile}{$s} = $s;
		if ($s =~ /^chr(.+)$/) {
			$BIGWIG_CHROMOS{$bwfile}{$1} = $s;
		}
		else {
			$BIGWIG_CHROMOS{$bwfile}{"chr$s"} = $s;
		}
	}
}



### Internal subroutine for processing summary features
# for personal use only
sub _process_summaries {
	my ($method, $summaries) = @_;
	
	# calculate a value based on the number of summary features
	if (scalar @$summaries == 1) {
		# great! only one summary returned
		return &{ $ONE_SUMMARY_METHOD{$method} }($summaries->[0]);
	}
	elsif (scalar @$summaries > 1) {
		# more than one summary
		return &{ $MULTI_SUMMARY_METHOD{$method} }($summaries);
	}
	else {
		# no summaries
		return $method =~ /sum|count/ ? 0 : '.';
	}
}



### Internal subroutine for processing position and score from features
# for personal use only
sub _process_position_score_feature {
	
	# collect the feature and hashes
	my ($f, $pos2data, $duplicates) = @_;
	
	# process across the length of this feature
		# for most wig features this is almost certainly a single point
		# but just in case this is spanned data, we will record the 
		# value at every position
	for (my $pos = $f->start; $pos <= $f->end; $pos++) {
		
		# check for duplicate positions
		# duplicates should not normally exist for a single wig file
		# but we may be working with multiple wigs that are being combined
		if (exists $pos2data->{$pos} ) {
			
			if (exists $duplicates->{$pos} ) {
				# we have lots of duplicates at this position!
				
				# append an incrementing number at the end
				$duplicates->{$pos} += 1; # increment first
				my $new = sprintf("%d.%d", $pos, $duplicates->{$pos});
				$pos2data->{$new} = $f->score;
			}
			else {
				# first time duplicate
				my $new = $pos . '.1';
				$pos2data->{$new} = $f->score;
				$duplicates->{$pos} = 1;
			}
		}
		else {
			$pos2data->{$pos} = $f->score;
		}
	}
}


### Internal subroutine for removing duplicate positions from pos2data hash
# for personal use only
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
		$pos2data->{$pos} = mean(@values);
	}
}



__END__

=head1 NAME

Bio::ToolBox::db_helper::bigwig

=head1 DESCRIPTION

This module provides support for binary BigWig files to the 
L<Bio::ToolBox> package. It also supports a directory of one 
or more bigWig files as a combined database, known as a 
BigWigSet. 

=head1 USAGE

The module requires L<Bio::DB::BigWig> to be installed, which in turn 
requires the UCSC Kent C library to be installed.

In general, this module should not be used directly. Use the methods 
available in L<Bio::ToolBox::db_helper> or <Bio::ToolBox::Data>.  

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

=item collect_bigwigset_position_score()

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

Acceptable values include mean, min, max, stddev, sum, and count.
Used when collecting a single value over a genomic segnment.

=item 7. The value type of data to collect

Typically only 'score' is accepted here.

=item 8. A database object.

Pass the opened L<Bio::DB::BigWigSet> database object when working 
with BigWigSets.

=item 9 and higher. BigWig file names or BigWigSet database types.

Opened BigWig objects are cached. Both local and remote BigWig files 
are supported. 

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



