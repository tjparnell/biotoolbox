package Bio::ToolBox::db_helper::seqfasta;

require Exporter;
use Carp;
use strict;
use Module::Load; # for dynamic loading during runtime
use Statistics::Lite qw(mean);
use Bio::ToolBox::db_helper::config;
use Bio::DB::Fasta;
use Bio::DB::SeqFeature::Store;

our $VERSION = '1.25';
our $WIGGLE_OK = 0;

# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	open_fasta_db
	open_store_db
	collect_store_scores
	collect_store_position_scores
);

1;

sub open_fasta_db {
	my $database = shift;
	my $db;
	
	eval {
		# to prevent annoying error messages, we silence warnings
		# we are often trying to open without absolute confirmation this will work
		# hence working inside an eval
		local $SIG{__WARN__} = sub {}; 
		$db = Bio::DB::Fasta->new($database);
	};
	unless ($db) {
		warn "unable to open fasta file as database! Make sure the file's directory\n" .
			"writable, and if a *.index file exists, try deleting it so it can be rebuilt\n";
	}
	return $db;
}


sub open_store_db {
	my $database = shift;
	my $db;
	
	# first determine type of database we're dealing with
	# by checking for an extension
	
	# GFF3 file to be loaded into memory
	if ($database =~ /\.gff3?(?:\.gz)?$/i) {
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
	
	elsif ($database =~ /\.(?:db|sqlite)$/i) {
		eval {
			$db = Bio::DB::SeqFeature::Store->new(
				-adaptor  => 'DBI::SQLite',
				-dsn      => $database,
			);
		};
	}
	
	else {
		# open the connection using parameters from the configuration file
		# we'll try to use database specific parameters first, else use 
		# the db_default parameters
		my $adaptor = $BTB_CONFIG->param($database . '.adaptor') || 
			$BTB_CONFIG->param('default_db.adaptor');
		my $user = $BTB_CONFIG->param($database . '.user') || 
			$BTB_CONFIG->param('default_db.user');
		my $pass = $BTB_CONFIG->param($database . '.pass') ||
			$BTB_CONFIG->param('default_db.pass') || undef;
		
		# check for empty passwords
		# config::simple passes an empty array when nothing was defined
		if (ref $pass eq 'ARRAY' and scalar @$pass == 0) {$pass = undef}
		
		# set up the dsn
		# it can be specifically defined
		my $dsn = $BTB_CONFIG->param($database . '.dsn') || undef;
		unless (defined $dsn) {
			# or dsn can be generated with the dsn_prefix
			$dsn = $BTB_CONFIG->param($database . '.dsn_prefix') || 
				$BTB_CONFIG->param('default_db.dsn_prefix');
			$dsn .= $database;
		}
		
		# establish the database connection
		eval {
			# to prevent annoying error messages from B:DB:SF:S
			local $SIG{__WARN__} = sub {}; 
		
			# attempt a connection
			$db = Bio::DB::SeqFeature::Store->new(
				-adaptor => $adaptor,
				-dsn     => $dsn,
				-user    => $user,
				-pass    => $pass,
			);
		};
	}
	
	unless ($db) {
		warn "Failed to open Bio::DB::SeqFeature::Store database '$database'";
	}
	
	# return opened database object or nothing if unsuccessful 
	return $db;
}



sub collect_store_scores {
	# collect only the scores
	return _collect_store_data(0, @_);
}


sub collect_store_position_scores {
	# collect positioned scores
	return _collect_store_data(1, @_);
}



sub _collect_store_data {
	# pass the required information
	unless (scalar @_ >= 9) {
		confess " At least eight arguments must be passed to collect db scores!\n";
	}
	my ($do_index, $db, $chromo, $start, $stop, $strand, $strandedness, 
		$value_type, @types) = @_;
	
	my $db_ref = ref $db;
	unless ($db_ref =~ /Bio::DB/) {
		croak "An opened database object must be passed!";
	}
	
	# set up iterator
	my $iterator = $db->get_seq_stream(
		-seq_id      => $chromo,
		-start       => $start,
		-end         => $stop,
		-primary_tag => [@types],
	);
	return unless $iterator;
	
	# collect the first feature
	my $feature = $iterator->next_seq;
	return unless $feature;
	
	# deal with features that might not be from the chromosome we want
	# sometimes chromosome matching is sloppy and we get something else
	# this really should not happen....
	while ($feature->seq_id ne $chromo) {
		$feature = $iterator->next_seq || undef;
	}
	return unless $feature;
	
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
		if ($wigfile =~ /\.wib$/) {
			# data is in old-style binary wiggle format
			# based on the Bio::Graphics::Wiggle adaptor
			
			# get the full list of features to pass off to the 
			# helper subroutine
			my @features;
			push @features, $feature;
			while (my $f = $iterator->next_seq) {
				push @features, $f;
			}
			
			# check that we have wiggle support
			unless ($WIGGLE_OK) {
				eval {
					my $class = 'Bio::ToolBox::db_helper::wiggle';
					load $class;
					$class->import;
					$WIGGLE_OK = 1;
				};
			}
			if ($WIGGLE_OK) {
				# get the dataset scores using Bio::ToolBox::db_helper::wiggle
				
				if ($do_index) {
					# warn " using collect_wig_position_scores() from tag\n";
					return collect_wig_position_scores(
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type,
						@features
					);
				}
				else {
					# warn " using collect_wig_scores() from tag\n";
					return collect_wig_scores(
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type,
						@features
					);
				}
			}
			else {
				croak " Wiggle support is not enabled! " . 
					"Is Bio::Graphics::Wiggle installed?\n";
			}
		}
		
		## BigWig file
		elsif ($wigfile =~ /\.bw$/) {
			# data is in bigwig format
			# we are abandoning this support
			die " Supporting bigWig files via seqfeature attribute is no longer " . 
				"supported.\n Please use bigWig files directly\n";
		}	
		else {
			croak " Unrecognized wigfile attribute '$wigfile'!" . 
				" Unable to continue!\n";
		}
	}
	
	
	## Database Data
	# Working with data stored directly in the database
	# this is more straight forward in collection
	
	# Walk through the datapoints
	while ($feature) {
	
		# Check which data to take based on strand
		if (
			$strandedness eq 'all' # all data is requested
			or $strand == 0 # region is unstranded
			or $feature->strand == 0 # unstranded data
			or ( 
				# sense data
				$strand == $feature->strand 
				and $strandedness eq 'sense'
			) 
			or (
				# antisense data
				$strand != $feature->strand  
				and $strandedness eq 'antisense'
			)
		) {
			# we have acceptable data to collect
		
			# data is in the database
			# much easier to collect
			
			# store data in either indexed hash or score array
			if ($do_index) {
			
				# determine position to record
				my $position;
				if ($feature->start == $feature->end) {
					# just one position recorded
					$position = $feature->start;
				}
				else {
					# calculate the midpoint
					$position = int( 
						($feature->start + $feature->end) / 2
					);
				}
				
				# store the appropriate value
				if ($value_type eq 'score') {
					push @{ $pos2data{$position} }, $feature->score;
				}
				elsif ($value_type eq 'count') {
					$pos2data{$position} += 1;
				}
				elsif ($value_type eq 'pcount') {
					$pos2data{$position} += 1 if 
						($feature->start >= $start and $feature->end <= $stop);
				}
				elsif ($value_type eq 'length') {
					push @{ $pos2data{$position} }, $feature->length;
				}
			}
			
			else {
				# just store the score in the array
				
				# store the appropriate value
				if ($value_type eq 'score') {
					push @scores, $feature->score;
				}
				elsif ($value_type eq 'count') {
					$scores[0] += 1;
				}
				elsif ($value_type eq 'pcount') {
					$scores[0] += 1 if 
						($feature->start >= $start and $feature->end <= $stop);
				}
				elsif ($value_type eq 'length') {
					push @scores, $feature->length;
				}
			}
		}
		
		# prepare for next
		$feature = $iterator->next_seq || undef;
	}
	
	# post-process the collected position->score values 
	# combine multiple values recorded at the same position
	if (
		$do_index and 
		($value_type eq 'score' or $value_type eq 'length')
	) {
		# each 'value' is an array of one or more scores or lengths 
		# from the datapoints collected above
		# the mean value is the best we can do right now for 
		# combining the data
		# really would prefer something else
		# we don't have a true method to utilize
		foreach my $position (keys %pos2data) {
			$pos2data{$position} = mean( @{$pos2data{$position}} );
		}
	}
	
	# return the appropriate data	
	return $do_index ? %pos2data : @scores;
}



__END__


=head1 NAME

Bio::ToolBox::db_helper::seqfasta

=head1 DESCRIPTION

This module supports opening L<Bio::DB::SeqFeature::Store> and 
L<Bio::DB::Fasta> BioPerl database adaptors. It also supports collecting 
feature scores from L<Bio::DB::SeqFeature::Store> databases. Unsupported 
BioPerl-style database adaptors that support generic methods may also 
be used, although success may vary.

=head2 Opening databases

For Fasta databases, either a single fasta file or a directory of fasta 
files may be provided.

For SeqFeature Store databases, the connection parameters are stored in 
a configuration file, C<.biotoolbox.cfg>. Multiple database containers 
are supported, including MySQL, SQLite, and in-memory. 

=head2 Collecting scores

Scores from seqfeature objects stored in the database may be retrieved. The 
scores may be collected as is, or they may be associated with genomic positions 
(indexed scores). Scores may be restricted to strand by specifying the desired 
strandedness. For example, to collect transcription data over a gene, pass the 
strandedness value 'sense'. If the strand of the region database object 
(representing the gene) matches the strand of the wig file data feature, then 
the data is collected.

Legacy wig file support uses GFF SeqFeature databases to store the file paths 
of the binary wiggle (.wib) files. If the seqfeature objects returned from the 
database include the wigfile attribute, then these objects are forwarded on to 
the L<Bio::ToolBox::db_helper::wiggle> adaptor for appropriate score collection.

=head1 USAGE

The module requires the BioPerl adaptors L<Bio::DB::SeqFeature::Store> and 
L<Bio::DB::Fasta>. 

Load the module at the beginning of your program.

	use Bio::ToolBox::db_helper::seqfasta;

It will automatically export the name of the subroutines. 

=over

=item collect_store_scores

This subroutine will collect only the score values from database features
for the specified database region. The positional information of the 
scores is not retained, and the values may be further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed eight or more arguments in the following order:
    
    1) The opened database object. A database name or file is not ok.
    2) The chromosome name
    3) The start position of the segment to collect from
    4) The stop or end position of the segment to collect from
    5) The strand of the original feature (or region), -1, 0, or 1.
    6) A scalar value representing the desired strandedness of the data 
       to be collected. Only those scores which match the indicated 
       strandedness are collected. Acceptable values include 
        "sense", 
        "antisense", 
        "none" or "no".
    7) The type of data collected. Acceptable values include 
       'score' (returns the score), 
       'count' (the number of defined positions with scores), or 
       'length' (the wig step is used here).  
    8) One or more feature types or primary_tags to perform the 
       database search. If nothing is provided, then usually everything 
       in the database is returned!

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_wig_position_scores

This subroutine will collect the score values form features in the database 
for the specified region keyed by position. 

The subroutine is passed the same arguments as collect_wig_scores().

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. Note that only one value is 
returned per position, regardless of the number of dataset features 
passed.

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

