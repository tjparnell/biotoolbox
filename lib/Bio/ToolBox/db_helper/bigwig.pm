package Bio::ToolBox::db_helper::bigwig;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(min max mean);
use Bio::ToolBox::db_helper::constants;
use Bio::DB::BigWig qw(binMean binStdev);
use Bio::DB::BigWigSet;
our $VERSION = '1.51';


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

# BigWigSet bigWig IDs
our %BIGWIGSET_WIGS; 
	# cache for the bigwigs from a BigWigSet used in a query
	# we want to use low level bigWig access which isn't normally 
	# available from the high level BigWigSet, so we identify the 
	# bigWigs from the bigWigSet and cache them here 


# Methods Cache lookup
	# pre-define the methods code for processing summary objects to avoid 
	# excessive if-elsif tests every data cycle
our %ONE_SUMMARY_METHOD = (
	'mean'  => sub {return unless $_[0]; return binMean($_[0]) },
	'sum'   => sub {return 0 unless $_[0]; return $_[0]->{sumData} },
	'min'   => sub {return unless $_[0]; return $_[0]->{minVal} },
	'max'   => sub {return unless $_[0]; return $_[0]->{maxVal} },
	'count' => sub {return 0 unless $_[0]; return $_[0]->{validCount} },
	'stddev'=> sub {return unless $_[0]; return binStdev($_[0]) },
);
our %MULTI_SUMMARY_METHOD = (
	'mean'  => sub {
				my $summaries = shift;
				return unless $summaries->[0];
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
				return 0 unless $summaries->[0];
				my $sum = 0;
				foreach my $s (@$summaries) {
					$sum += $s->{sumData};
				}
				return $sum;
			},
	'min'   => sub {
				my $summaries = shift;
				return unless $summaries->[0];
				return min( map {$_->{minVal}} @$summaries);
			},
	'max'   => sub {
				my $summaries = shift;
				return unless $summaries->[0];
				return max( map {$_->{maxVal}} @$summaries);
			},
	'count' => sub {
				my $summaries = shift;
				return 0 unless $summaries->[0];
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
	# chromosome, start, stop, strand, strandedness, method, db, dataset
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
		my $sumArrays = $bw->bf->bigWigSummaryArrayExtended(
			# chromo, 0-based start, stop, # bins
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
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	
	# we don't have support for pcount and ncount, so make sure we change those
	$param->[METH] = 'count' if $param->[METH] =~ /\wcount/;
	
	# Walk through each requested feature
	# There is likely only one
	my @scores; 
	for (my $d = DATA; $d < scalar @$param; $d++) {
	
		# open db object
		my $bw = _get_bw($param->[$d]);
		
		# first check that the chromosome is present
		my $chromo = $BIGWIG_CHROMOS{$param->[$d]}{$param->[CHR]} or next;
		
		# initialize low level stream for this segment
		my $list = $bw->bf->bigWigIntervalQuery(
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
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	
	# we don't have support for pcount and ncount, so make sure we change those
	$param->[METH] = 'count' if $param->[METH] =~ /\wcount/;
	
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
		my $list = $bw->bf->bigWigIntervalQuery(
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
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;
	
	# lookup the bigWig files based on the parameters
	my $ids = _lookup_bigwigset_wigs($param);
	return unless scalar(@$ids) > 0;
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
	
	# check for potential implied strandedness based on the file name
		# because the database features are not true SeqFeature objects, 
		# the strand method isn't really supported as expected with SeqFeature
		# instead it is treated as an attribute, and typical strand convention 
		# of -1, 0, 1 gets messed up with regex matching, so we'll silently 
		# use minus, none, and plus as strand attribute values
	my $md = $bws->metadata;
	foreach my $i (keys %$md) {
		my $f = $md->{$i}{dbid}; # the file path
		if ($f =~ /(?:f|for|forward|top|plus|\+)\.bw$/i) {
			# implied forward strand 
			if (exists $md->{$i}{strand}) {
				$md->{$i}{strand} = 'plus' if $md->{$i}{strand} eq '1';
			}
			else {
				$bws->set_bigwig_attributes($f, {strand => 'plus'});
			}
		}
		elsif ($f =~ /(?:r|rev|reverse|bottom|minus|\-)\.bw$/i) {
			# implied reverse strand 
			if (exists $md->{$i}{strand}) {
				$md->{$i}{strand} = 'minus' if $md->{$i}{strand} eq '-1';
			}
			else {
				$bws->set_bigwig_attributes($f, {strand => 'minus'});
			}
		}
		else {
			# check for strand anyway
			if (exists $md->{$i}{strand}) {
				$md->{$i}{strand} = 'plus' if $md->{$i}{strand} eq '1';
				$md->{$i}{strand} = 'minus' if $md->{$i}{strand} eq '-1';
				$md->{$i}{strand} = 'none' if $md->{$i}{strand} eq '0';
			}
			else {
				# assign a non-strand just in case
				$md->{$i}{strand} = 'none';
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


### Internal subroutine to identify and cache selected bigWigs from a BigWigSet
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
	
	# I want to access the bigWigs through the low level API for speed
	# but first I need to find out which ones to use
	# use internal methods to filter the bigWigs in the same manner 
	# that get_seq_stream does in BigWigSet
	my @ids = $param->[DB]->bigwigs;
	
	# filter first by type
	@ids = $param->[DB]->_filter_ids_by_type(\@types, \@ids );
	
	# try filtering by name if that doesn't work
	unless (@ids) {
    	my @bwids = $param->[DB]->bigwigs;
    	foreach my $t (@types) {
    		my @found = $param->[DB]->_filter_ids_by_name($t, \@bwids);
			push @ids, @found if (@found);
    	}
    }
	
	# then check for strand
	if ($param->[STND] ne 'all' and $param->[STR] != 0) {
    	# looks like we are collecting stranded data
    	# try to filter again based on strand attribute
    	# remember to accept 0 strand attributes as well
    	# create an array of acceptable strand values to filter on
    		# remember that strand is an ordinary attribute and does not 
    		# behave like the typical SeqFeature strand attribute
    	my @strands;
    	if ($param->[STND] eq 'sense') {
			@strands = 
				$param->[STR] == -1 ? qw(minus none) :
				$param->[STR] == 1  ? qw(plus none) : qw(minus plus none);
		}
		elsif ($param->[STND] eq 'antisense') {
			@strands = 
				$param->[STR] == -1 ? qw(plus none) :
				$param->[STR] == 1  ? qw(minus none) : qw(minus plus none);
		}
		else {
			confess sprintf "bad strandedness value: %s", $param->[STND];
		}
    	# then filter based on attribute
    	@ids = $param->[DB]->_filter_ids_by_attribute(
    		{strand => \@strands}, \@ids);
    }
	
	# cache and return
	$BIGWIGSET_WIGS{$lookup} = \@ids;
	return \@ids;
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
Methods of pcount and ncount are technically supported, but are 
treated the same as count.

=item 7. A database object.

Pass the opened L<Bio::DB::BigWigSet> database object when working 
with BigWigSets.

=item 8 and higher. BigWig file names or BigWigSet database types.

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



