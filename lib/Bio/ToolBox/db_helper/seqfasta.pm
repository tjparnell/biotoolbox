package Bio::ToolBox::db_helper::seqfasta;

use warnings;
use strict;
use Carp;
use English qw(-no_match_vars);
use Module::Load;    # for dynamic loading during runtime
use List::Util qw(min max sum);
use Statistics::Lite qw(median);
use Bio::ToolBox::db_helper::constants;
use Bio::ToolBox::db_helper::config;
use Bio::DB::Fasta;
use Bio::DB::SeqFeature::Store;
require Exporter;

our $VERSION   = '1.70';

# Exported names
our @ISA    = qw(Exporter);

## no critic
## this is never intended to be used directly by end users
## and exporting everything is required
our @EXPORT = qw(
	open_fasta_db
	open_store_db
	collect_store_scores
);
## use critic

my $WIGGLE_OK = 0;
my $CONFIG_OK = 0;

sub open_fasta_db {
	my $database = shift;
	my $db;

	eval {
		# to prevent annoying error messages, we silence warnings
		# we are often trying to open without absolute confirmation this will work
		# hence working inside an eval
		local $SIG{__WARN__} = sub { };
		$db = Bio::DB::Fasta->new($database);
	};
	unless ($db) {
		warn "unable to open fasta file as database! Make sure the file's directory\n"
			. "writable, and if a *.index file exists, try deleting it so it can be rebuilt\n";
	}
	return $db;
}

sub open_store_db {
	my $database = shift;
	my $db;

	# first determine type of database we're dealing with
	# by checking for an extension

	# GFF3 file to be loaded into memory
	if ( $database =~ m/\.gff3? (?:\.gz)? $/xi ) {

		# gff3 file can be gzipped
		# this might take a while, so print a statement
		print " Loading file into memory database...\n";
		eval {
			$db = Bio::DB::SeqFeature::Store->new(
				-adaptor => 'memory',
				-gff     => $database,
			);
		};
	}

	elsif ( $database =~ m/\. (?: db | sqlite )$/xi ) {
		eval {
			$db = Bio::DB::SeqFeature::Store->new(
				-adaptor => 'DBI::SQLite',
				-dsn     => $database,
			);
		};
	}

	else {
		# open the connection using parameters from the configuration file
		# we'll try to use database specific parameters first, else use
		# the db_default parameters
		my $adaptor = $BTB_CONFIG->param( $database . '.adaptor' )
			|| $BTB_CONFIG->param('default_db.adaptor');
		my $user = $BTB_CONFIG->param( $database . '.user' )
			|| $BTB_CONFIG->param('default_db.user');
		my $pass =
			   $BTB_CONFIG->param( $database . '.pass' )
			|| $BTB_CONFIG->param('default_db.pass')
			|| undef;

		# check for empty passwords
		# config::simple passes an empty array when nothing was defined
		if ( ref($pass) eq 'ARRAY' and scalar @{ $pass } == 0 ) { $pass = undef }

		# set up the dsn
		# it can be specifically defined
		my $dsn = $BTB_CONFIG->param( $database . '.dsn' ) || undef;
		unless ( defined $dsn ) {

			# or dsn can be generated with the dsn_prefix
			$dsn = $BTB_CONFIG->param( $database . '.dsn_prefix' )
				|| $BTB_CONFIG->param('default_db.dsn_prefix');
			$dsn .= $database;
		}

		# establish the database connection
		eval {
			# to prevent annoying error messages from B:DB:SF:S
			local $SIG{__WARN__} = sub { };

			# attempt a connection
			$db = Bio::DB::SeqFeature::Store->new(
				-adaptor => $adaptor,
				-dsn     => $dsn,
				-user    => $user,
				-pass    => $pass,
			);
		};
	}

	# return opened database object or nothing if unsuccessful
	return unless $db;
	return $db;
}

sub collect_store_scores {

	# pass the required information
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;

	# database feature types
	my @types = splice( @{ $param }, DATA );

	# set up iterator from database
	my $iterator = $param->[DB]->get_seq_stream(
		-seq_id      => $param->[CHR],
		-start       => $param->[STRT],
		-end         => $param->[STOP],
		-primary_tag => \@types,
	);
	return unless $iterator;    # we should always get an iterator back, even if the
								# iterator returns nothing useful for us

	# collect the first feature
	my $feature = $iterator->next_seq;
	unless ($feature) {

		# that's odd, we should get a feature
		# try rebuilding the iterator with an alternate chromosome name
		my $chromo = $param->[CHR] =~ /^chr(.+)$/i ? $1 : 'chr' . $param->[CHR];
		$iterator = $param->[DB]->get_seq_stream(
			-seq_id      => $chromo,
			-start       => $param->[STRT],
			-end         => $param->[STOP],
			-primary_tag => \@types,
		);
		$feature = $iterator->next_seq or return;
	}

	# places to stick the scores
	my @scores;
	my %pos2data;

	## First check for potential Wig Data
	# the legacy wiggle adaptor uses a GFF seqfeature db to store file paths
	# to the individual chromosome wib files
	# if this is the case, we must pass this on to the wiggle adaptor
	if ( $feature->has_tag('wigfile') ) {

		# data is in wig format, or at least the first datapoint is

		# determine the type of wigfile
		my ($wigfile) = $feature->get_tag_values('wigfile');

		## Bio::Graphics wib file
		if ( $wigfile =~ /\.wib$/ ) {

			# data is in old-style binary wiggle format
			# based on the Bio::Graphics::Wiggle adaptor

			# get the full list of features to pass off to the
			# helper subroutine
			push @{ $param }, $feature;
			while ( my $f = $iterator->next_seq ) {
				push @{ $param }, $f;
			}

			# check that we have wiggle support
			unless ($WIGGLE_OK) {
				eval {
					my $class = 'Bio::ToolBox::db_helper::wiggle';
					autoload $class;
					$WIGGLE_OK = 1;
				};
			}
			if ($WIGGLE_OK) {

				# get the dataset scores using Bio::ToolBox::db_helper::wiggle

				if ( $param->[RETT] == 2 ) {

					# warn " using collect_wig_position_scores() from tag\n";
					return collect_wig_position_scores($param);
				}
				else {
					# warn " using collect_wig_scores() from tag\n";
					return collect_wig_scores($param);
				}
			}
			else {
				croak " Wiggle support is not enabled! $EVAL_ERROR\n";
			}
		}

		## BigWig file
		elsif ( $wigfile =~ /\.bw$/ ) {

			# data is in bigwig format
			# we are abandoning this support
			die ' Supporting bigWig files via seqfeature attribute is no longer '
				. "supported.\n Please use bigWig files directly\n";
		}
		else {
			croak " Unrecognized wigfile attribute '$wigfile'!"
				. " Unable to continue!\n";
		}
	}

	## Database Data
	# Working with data stored directly in the database
	# this is more straight forward in collection

	# Walk through the datapoints
	while ($feature) {

		# Check which data to take based on strand
		if (
			$param->[STND] eq 'all'     # all data is requested
			or $feature->strand == 0    # unstranded data
			or (
				# sense data
				$param->[STR] == $feature->strand and $param->[STND] eq 'sense'
			)
			or (
				# antisense data
				$param->[STR] != $feature->strand and $param->[STND] eq 'antisense'
			)
			)
		{
			# we have acceptable data to collect

			# data is in the database
			# much easier to collect

			# store data in either indexed hash or score array
			if ( $param->[RETT] == 2 ) {

				# determine position to record
				my $position;
				if ( $feature->start == $feature->end ) {

					# just one position recorded
					$position = $feature->start;
				}
				else {
					# calculate the midpoint
					$position = int( ( $feature->start + $feature->end ) / 2 );
				}

				# store the appropriate value
				if ( $param->[METH] eq 'count' ) {
					$pos2data{$position} += 1;
				}
				elsif ( $param->[METH] eq 'pcount' ) {
					$pos2data{$position} += 1
						if (    $feature->start >= $param->[STRT]
							and $feature->end <= $param->[STOP] );
				}
				elsif ( $param->[METH] eq 'ncount' ) {
					push @{ $pos2data{$position} },
						   $feature->display_name
						|| $feature->primary_id
						|| sprintf( "%s:%s-%s:%s",
							$feature->seq_id, $feature->start,
							$feature->end,    $feature->strand );
				}
				else {
					push @{ $pos2data{$position} }, $feature->score;
				}
			}

			else {
				# just store the score in the array
				# store the appropriate value
				if ( $param->[METH] eq 'count' ) {
					push @scores, 1;
				}
				elsif ( $param->[METH] eq 'pcount' ) {
					push @scores, 1
						if (    $feature->start >= $param->[STRT]
							and $feature->end <= $param->[STOP] );
				}
				elsif ( $param->[METH] eq 'ncount' ) {
					push @scores,
						   $feature->display_name
						|| $feature->primary_id
						|| sprintf( "%s:%s-%s:%s",
							$feature->seq_id, $feature->start,
							$feature->end,    $feature->strand );
				}
				else {
					push @scores, $feature->score;
				}
			}
		}

		# prepare for next
		$feature = $iterator->next_seq || undef;
	}

	# post-process the collected position->score values
	# combine multiple values recorded at the same position
	if ( $param->[RETT] == 2 ) {
		if ( $param->[METH] eq 'ncount' ) {
			foreach my $position ( keys %pos2data ) {
				my %name2count;
				foreach ( @{ $pos2data{$position} } ) { $name2count{$_} += 1 }
				$pos2data{$position} = scalar( keys %name2count );
			}
		}
		elsif ( $param->[METH] eq 'count' or $param->[METH] eq 'pcount' ) {

			# do nothing, these aren't arrays
		}
		elsif ( $param->[METH] eq 'mean' ) {
			foreach my $position ( keys %pos2data ) {
				$pos2data{$position} = sum( @{ $pos2data{$position} } ) /
					scalar( @{ $pos2data{$position} } );
			}
		}
		elsif ( $param->[METH] eq 'median' ) {
			foreach my $position ( keys %pos2data ) {
				$pos2data{$position} = median( @{ $pos2data{$position} } );
			}
		}
		elsif ( $param->[METH] eq 'min' ) {
			foreach my $position ( keys %pos2data ) {
				$pos2data{$position} = min( @{ $pos2data{$position} } );
			}
		}
		elsif ( $param->[METH] eq 'max' ) {
			foreach my $position ( keys %pos2data ) {
				$pos2data{$position} = max( @{ $pos2data{$position} } );
			}
		}
		elsif ( $param->[METH] eq 'sum' ) {
			foreach my $position ( keys %pos2data ) {
				$pos2data{$position} = sum( @{ $pos2data{$position} } );
			}
		}
		else {
			# just take the mean for everything else
			foreach my $position ( keys %pos2data ) {
				$pos2data{$position} = sum( @{ $pos2data{$position} } ) /
					scalar( @{ $pos2data{$position} } );
			}
		}
	}

	# return the appropriate data
	if ( $param->[RETT] == 2 ) {
		return wantarray ? %pos2data : \%pos2data;
	}
	else {
		return wantarray ? @scores : \@scores;
	}
}

1;

__END__


=head1 NAME

Bio::ToolBox::db_helper::seqfasta

=head1 DESCRIPTION

This module provides support for L<Bio::DB::SeqFeature::Store> and 
L<Bio::DB::Fasta> BioPerl database adaptors. Other BioPerl-style 
databases may have limited support.

=head1 USAGE

The module requires BioPerl to be installed, which includes database adaptors 
L<Bio::DB::SeqFeature::Store> and L<Bio::DB::Fasta>. The database storage 
adapter, e.g MySQL, SQLite, etc., will have additional requirements.

In general, this module should not be used directly. Use the methods 
available in L<Bio::ToolBox::db_helper> or <Bio::ToolBox::Data>.  

All subroutines are exported by default.

=head2 Opening databases

=over 4

=item open_store_db

Opens a L<Bio::DB::SeqFeature::Store> database. The connection parameters 
are typically stored in a configuration file, C<.biotoolbox.cfg>. Multiple 
database containers are supported, including MySQL, SQLite, and in-memory. 

Pass the name of the database or database file.

=item open_fasta_db

Opens a L<Bio::DB::Fasta> database. Either a single fasta file or a directory 
of fasta files may be provided. Pass the path name to the file or directory.

=back

=head2 Collecting scores

Scores from seqfeature objects stored in the database may be retrieved. The 
scores may be collected as is, or they may be associated with genomic positions 
(indexed scores). Scores may be restricted to strand by specifying the desired 
strandedness. For example, to collect transcription data over a gene, pass the 
strandedness value 'sense'. If the strand of the region database object 
(representing the gene) matches the strand of the wig file data feature, then 
the data is collected.

Legacy wig file support uses GFF SeqFeature databases to store the file paths 
of the binary wiggle (F<.wib>) files. If the seqfeature objects returned from the 
database include the wigfile attribute, then these objects are forwarded on to 
the L<Bio::ToolBox::db_helper::wiggle> adaptor for appropriate score collection.

=over

=item collect_store_scores

This subroutine will collect only the score values from database features
for the specified database region. The positional information of the 
scores is not retained, and the values may be further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed a parameter array reference. See 
L<"Data Collection Parameters Reference"> below for details.

The subroutine returns an array or array reference of the requested dataset 
values found within the region of interest. 

=item collect_wig_position_scores

This subroutine will collect the score values form features in the database 
for the specified region keyed by position. 

The subroutine is passed a parameter array reference. See 
L</"Data Collection Parameters Reference"> below for details.

The subroutine returns a hash or hash reference of the requested dataset 
values found within the region of interest keyed by position. Note that only 
one value is returned per position, regardless of the number of dataset features 
passed.

=back

=head2 Data Collection Parameters Reference

The data collection subroutines are passed an array reference of parameters. 
The recommended  method for data collection is to use the 
L<Bio::ToolBox::db_helper/get_segment_score> method. 

The parameters array reference includes these items:

=over 4

=item 1. chromosome

=item 1. start coordinate

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

This should be a L<Bio::DB::SeqFeature::Store> database.

=item 8. database type

Set the C<type> or C<primary_tag> of the dataset within the database 
to collect from. Additional dataset items may be added to the list 
when merging data.

=back

=head1 SEE ALSO

L<Bio::ToolBox::Data::Feature>, L<Bio::ToolBox::db_helper>, 
L<Bio::DB::SeqFeature::Store>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

