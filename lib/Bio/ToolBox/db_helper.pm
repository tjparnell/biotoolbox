package Bio::ToolBox::db_helper;

use warnings;
use strict;
use Carp qw(carp cluck croak confess);
use Module::Load;    # for dynamic loading during runtime
use List::Util qw(min max sum0 uniq);
use Statistics::Lite qw(median range stddevp);
use Bio::ToolBox::db_helper::constants;
use Bio::ToolBox::utility;
require Exporter;

our $VERSION = '1.70';

# check values for dynamically loaded helper modules
# these are loaded only when needed during runtime to avoid wasting resources
our $BAM_OK      = 0;
our $BIGBED_OK   = 0;
our $BIGWIG_OK   = 0;
our $SEQFASTA_OK = 0;
our $USEQ_OK     = 0;
our $BAM_ADAPTER = undef;    # preference for which bam adapter to use
our $BIG_ADAPTER = undef;

# define reusable variables
my %TOTAL_READ_NUMBER;       # for rpm calculations
my $PRIMARY_ID_WARNING;      # for out of date primary IDs
my %OPENED_DB;               # cache for some opened Bio::DB databases
my %DB_METHODS;              # cache for database score methods

# score calculators
my %SCORE_CALCULATOR_SUB = (
	'mean' => sub {
		my $s = shift;
		return sum0(@{ $s }) / ( scalar(@{ $s }) || 1 );
	},
	'sum' => sub {
		my $s = shift;
		return sum0(@{ $s });
	},
	'median' => sub {
		my $s = shift;
		return '.' unless scalar(@{ $s });
		return median(@{ $s });
	},
	'min' => sub {
		my $s = shift;
		return '.' unless scalar(@{ $s });
		return min(@{ $s });
	},
	'max' => sub {
		my $s = shift;
		return '.' unless scalar(@{ $s });
		return max(@{ $s });
	},
	'count' => sub {
		my $s = shift;
		return scalar(@{ $s });
	},
	'pcount' => sub {
		my $s = shift;
		return scalar(@{ $s });
	},
	'ncount' => sub {

		# Convert names into unique counts
		my $s = shift;
		my %name2count;
		foreach my $n (@{ $s }) {
			if ( ref $n eq 'ARRAY' ) {

				# this is likely from a ncount indexed hash
				foreach (@{ $n }) {
					$name2count{$_} += 1;
				}
			}
			else {
				$name2count{$n} += 1;
			}
		}
		return scalar( keys %name2count );
	},
	'range' => sub {

		# the range value is 'min-max'
		my $s = shift;
		return '.' unless scalar(@{ $s });
		return range(@{ $s });
	},
	'stddev' => sub {

		# we are using the standard deviation of the population,
		# since these are the only scores we are considering
		my $s = shift;
		return '.' unless scalar(@{ $s });
		return stddevp(@{ $s });
	},
	'rpm' => sub {
		confess " The rpm methods have been deprecated due to complexity and "
			. "the variable way of calculating the value. Collect counts and "
			. "calculate your preferred way.\n";
	},
	'rpkm' => sub {
		confess " The rpm methods have been deprecated due to complexity and "
			. "the variable way of calculating the value. Collect counts and "
			. "calculate your preferred way.\n";
	}
);

# Exported names
our @ISA       = qw(Exporter);
our @EXPORT_OK = qw(
	$BAM_ADAPTER
	$BIG_ADAPTER
	use_bam_adapter
	use_big_adapter
	open_db_connection
	get_dataset_list
	verify_or_request_feature_types
	check_dataset_for_rpm_support
	get_new_feature_list
	get_new_genome_list
	validate_included_feature
	get_db_feature
	get_segment_score
	calculate_score
	get_chromosome_list
	low_level_bam_coverage
	low_level_bam_fetch
	get_genomic_sequence
);

### Open a connection to the SeqFeature Store MySQL database

sub use_bam_adapter {
	my $a = shift || undef;
	$BAM_ADAPTER = $a if $a;
	return $BAM_ADAPTER;
}

sub use_big_adapter {
	my $a = shift || undef;
	$BIG_ADAPTER = $a if $a;
	return $BIG_ADAPTER;
}

sub open_db_connection {
	my $database = shift;
	my $no_cache = shift || 0;
	unless ($database) {

		# 		cluck 'no database name passed!';
		return;
	}

	# first check if it is a database reference
	my $db_ref = ref $database;
	if ( $db_ref =~ /DB|big::BigWigSet/ ) {

		# the provided database is already an open database object
		# nothing to open, return as is
		return $database;
	}

	# check to see if we have already opened it
	if ( exists $OPENED_DB{$database} and not $no_cache ) {

		# return the cached database connection
		# but NOT if user explicitly requested no cached databases
		# DO NOT reuse database objects if you have forked!!! Bad things happen
		return $OPENED_DB{$database};
	}

	# skip parsed databases
	return if $database =~ /^Parsed:/;    # we don't open parsed annotation files

	# remove file prefix, just in case
	$database =~ s/^file://;

	### Attempt to open the database
	# we go through a series of checks to determine if it is remote, local,
	# an indexed big data file, SQLite file, etc
	# when all else fails, try to open a SQL connection
	my $db;
	my $error;

	# check if it is a remote file
	if ( $database =~ /^(?:https?|ftp)/i ) {

		# a remote Bam database
		if ( $database =~ /\.bam$/i ) {

			# open using Bam adaptor
			_load_bam_helper_module() unless $BAM_OK;
			if ($BAM_OK) {
				$db = open_bam_db($database);
				unless ($db) {
					$error =
						" ERROR: could not open remote Bam file" . " '$database'! $!\n";
				}
			}
			else {
				$error = " Bam database cannot be loaded because\n"
					. " Bio::DB::Sam or Bio::DB::HTS is not installed\n";
			}
		}

		# a remote BigBed database
		elsif ( $database =~ /\.(?:bb|bigbed)$/i ) {

			# open using BigBed adaptor
			$BIGBED_OK = _load_bigbed_helper_module() unless $BIGBED_OK;
			if ($BIGBED_OK) {
				$db = open_bigbed_db($database);
				unless ($db) {
					$error = " ERROR: could not open remote BigBed file"
						. " '$database'! $!\n";
				}
			}
			else {
				$error = " BigBed database cannot be loaded because\n"
					. " Bio::DB::Big or Bio::DB::BigBed is not installed\n";
			}
		}

		# a remote BigWig database
		elsif ( $database =~ /\.(?:bw|bigwig)$/i ) {

			# open using BigWig adaptor
			$BIGWIG_OK = _load_bigwig_helper_module() unless $BIGWIG_OK;
			if ($BIGWIG_OK) {
				$db = open_bigwig_db($database);
				unless ($db) {
					$error = " ERROR: could not open remote BigWig file"
						. " '$database'! $!\n";
				}
			}
			else {
				$error = " BigWig database cannot be loaded because\n"
					. " Bio::DB::Big or Bio::DB::BigWig is not installed\n";
			}
		}

		# a remote useq file
		elsif ( $database =~ /\.useq$/ ) {

			# uh oh! remote useq files are not supported
			$error = " ERROR: remote useq files are not supported!\n";
		}

		# a remote fasta file???
		elsif ( $database =~ /\.fa(?:sta)?$/i ) {

			# uh oh! remote fasta files are not supported
			$error = " ERROR: remote fasta files are not supported!\n";
		}

		# a presumed remote directory, presumably of bigwig files
		else {
			# open using BigWigSet adaptor
			$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::bigwig')
				unless $BIGWIG_OK;
			if ($BIGWIG_OK) {
				$db = open_bigwigset_db($database);
				unless ($db) {
					$error = " ERROR: could not open presumed remote "
						. "BigWigSet directory '$database'! $!\n";
				}
			}
			else {
				$error = " Presumed BigWigSet database cannot be loaded because\n"
					. " Bio::DB::BigWigSet is not installed\n";
			}
		}

	}

	# a local existing file
	elsif ( -f $database ) {

		# a Bam database
		if ( $database =~ /\.bam$/i ) {

			# open using appropriate bam adaptor
			_load_bam_helper_module() unless $BAM_OK;
			if ($BAM_OK) {
				undef $@;
				$db = open_bam_db($database);
				unless ($db) {
					$error =
						" ERROR: could not open local Bam file" . " '$database'! $@\n";
				}
			}
			else {
				$error = " Bam database cannot be loaded because\n"
					. " Bio::DB::Sam or Bio::DB::HTS is not installed\n";
			}
		}

		# a BigBed database
		elsif ( $database =~ /\.(?:bb|bigbed)$/i ) {

			# open using BigBed adaptor
			$BIGBED_OK = _load_bigbed_helper_module() unless $BIGBED_OK;
			if ($BIGBED_OK) {
				undef $@;
				$db = open_bigbed_db($database);
				unless ($db) {
					$error =
						" ERROR: could not open local BigBed file" . " '$database'! $@\n";
				}
			}
			else {
				$error = " BigBed database cannot be loaded because\n"
					. " Big::DB::Big or Bio::DB::BigBed is not installed\n";
			}
		}

		# a BigWig database
		elsif ( $database =~ /\.(?:bw|bigwig)$/i ) {

			# open using BigWig adaptor
			$BIGWIG_OK = _load_bigwig_helper_module() unless $BIGWIG_OK;
			if ($BIGWIG_OK) {
				undef $@;
				$db = open_bigwig_db($database);
				unless ($db) {
					$error =
						" ERROR: could not open local BigWig file" . " '$database'! $@\n";
				}
			}
			else {
				$error = " BigWig database cannot be loaded because\n"
					. " Big::DB::Big or Bio::DB::BigWig is not installed\n";
			}
		}

		# a useq database
		elsif ( $database =~ /\.useq$/i ) {

			# open using USeq adaptor
			$USEQ_OK = _load_helper_module('Bio::ToolBox::db_helper::useq')
				unless $USEQ_OK;
			if ($USEQ_OK) {
				$db = open_useq_db($database);
				unless ($db) {
					$error =
						" ERROR: could not open local useq file" . " '$database'! $!\n";
				}
			}
			else {
				$error = " Useq database cannot be loaded because\n"
					. " Bio::DB::USeq is not installed\n";
			}
		}

		# a Fasta File
		elsif ( $database =~ /\.fa(?:sta)?(?:\.gz)?$/i ) {

			# first try a modern fai indexed adapter
			_load_bam_helper_module() unless $BAM_OK;
			if ($BAM_OK) {
				$db = open_indexed_fasta($database);
				unless ($db) {
					$error .= " ERROR: could not open indexed fasta file '$database'!\n";
				}
			}
			else {
				# open using the old slow BioPerl Fasta adaptor
				$SEQFASTA_OK = _load_helper_module('Bio::ToolBox::db_helper::seqfasta')
					unless $SEQFASTA_OK;
				if ($SEQFASTA_OK) {
					$db = open_fasta_db($database);
					unless ($db) {
						$error .= " ERROR: could not open fasta file '$database'!\n";
						if ( -e "$database\.index" ) {
							$error .= "   Try deleting $database\.index and try again\n";
						}
					}
				}
				else {
					$error .= " Fasta file could not be loaded because neither"
						. " Bio::DB::HTS, Bio::DB::Sam, or Bio::DB::Fasta is installed\n";
				}
			}
		}

		# a gff3 or sqlite database
		elsif ( $database =~ /\.(?:gff3?|gff3?\.gz|db|sqlite)$/ ) {
			$SEQFASTA_OK = _load_helper_module('Bio::ToolBox::db_helper::seqfasta')
				unless $SEQFASTA_OK;
			if ($SEQFASTA_OK) {
				$db = open_store_db($database);
				unless ($db) {
					$error =
						" ERROR: could not load SeqFeature database file '$database'!\n";
				}
			}
			else {
				$error .= " Module Bio::DB::SeqFeature::Store is required to load "
					. "GFF and database files\n";
			}
		}

		# a cram database
		elsif ( $database =~ /\.cram$/i ) {

			# open using HTS bam adaptor only
			$BAM_ADAPTER ||= 'hts';
			if ( $BAM_ADAPTER eq 'sam' ) {
				$error .= " ERROR: Only HTS adapters support Cram files!\n";
			}
			_load_bam_helper_module() unless $BAM_OK;
			if ($BAM_OK) {
				undef $@;
				$db = open_bam_db($database);
				unless ($db) {
					$error =
						" ERROR: could not open local Cram file" . " '$database'! $@\n";
				}
			}
			else {
				$error = " Cram database cannot be loaded because\n"
					. " Bio::DB::HTS is not installed\n";
			}
		}

		# something unrecognized?
		else {
			$error .=
" ERROR! Cannot identify database type based on extension for $database!\n";
		}
	}

	# a directory, presumably of bigwig files
	elsif ( -d $database ) {

		# try opening using the BigWigSet adaptor
		$BIGWIG_OK = _load_bigwig_helper_module() unless $BIGWIG_OK;
		if ($BIGWIG_OK) {
			$db = open_bigwigset_db($database);
			unless ($db) {
				$error = " ERROR: could not open local BigWigSet "
					. "directory '$database'!\n";
				$error .= "   Does directory contain bigWig .bw files?\n";
			}
		}
		else {
			$error = " Presumed BigWigSet database cannot be loaded because\n"
				. " Bio::DB::Big or Bio::DB::BigWigSet is not installed\n";
		}

		# try opening with the Fasta adaptor
		unless ($db) {
			$SEQFASTA_OK = _load_helper_module('Bio::ToolBox::db_helper::seqfasta')
				unless $SEQFASTA_OK;
			if ($SEQFASTA_OK) {
				$db = open_fasta_db($database);
				unless ($db) {
					$error .= " ERROR: could not open fasta directory '$database'!\n";
					$error .= "   Does directory contain fasta files? If it contains a"
						. " directory.index file,\n   try deleting it and try again.\n";
				}
			}
			else {
				$error .= " Module Bio::DB::Fasta is required to open presumed fasta "
					. "directory\n";
			}
		}
	}

	# otherwise assume the name of a SeqFeature::Store database in the configuration
	else {
		# attempt to open using information from the configuration
		# using default connection information as necessary
		$SEQFASTA_OK = _load_helper_module('Bio::ToolBox::db_helper::seqfasta')
			unless $SEQFASTA_OK;
		if ($SEQFASTA_OK) {
			$db = open_store_db($database);
			unless ($db) {
				$error .= " unable to open relational database named '$database'\n";
			}
		}
		else {
			$error .= " Module Bio::DB::SeqFeature::Store is required to connect "
				. "to databases\n";
		}
	}

	# conditional return
	if ($db) {

		# cache the opened connection for use later
		$OPENED_DB{$database} = $db unless $no_cache;

		# return as appropriate either both object and name or just object
		return $db;
	}
	else {
		print STDERR $error;
		return;
	}
}

### Retrieve a database name from an db object
sub get_db_name {
	my $db = shift;
	return unless $db;
	my $db_ref = ref $db;
	return $db unless $db_ref;    # presumption that non-object is just a database name
	my $db_name;
	if ( $db_ref =~ /^Bio::DB::SeqFeature::Store/ ) {

		# a SeqFeature database, using any DBI adapter
		$db_name = $db->{'dbh'}->{'name'};

		# dig through the object internals to identify the original
		# name of the database
		# this should be relatively well documented through DBI
		# but could break in the future since it's not official API
	}
	elsif ( $db_ref eq 'Bio::DB::Sam' ) {

		# a Samtools Bam database
		# old versions <=1.39 don't have the bam_path method, so it's easier
		# to explicitly grab it from the object internals
		$db_name = $db->{bam_path};
	}
	elsif ( $db_ref eq 'Bio::DB::HTS' ) {

		# a HTSlib Bam database
		$db_name = $db->hts_path;
	}

	# determining the database name from other sources is
	# either not possible or not easy, so won't bother unless
	# there is a really really good need
	return $db_name;
}

### Retrieve a list of the microrarray data sets from the db
sub get_dataset_list {

	my $database         = shift;
	my $use_all_features = shift;

	# Open a db connection
	my $db = open_db_connection($database);
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	my $db_ref = ref $db;

	# process the database types, according to the type of database
	my @types;

	# a SeqFeature database
	if ( $db_ref =~ /^Bio::DB::SeqFeature::Store/ ) {

		# collect the database types,
		# sort first by source, then by method
		@types = (
			map $_->[2],
			sort { ( $a->[0] cmp $b->[0] ) or ( $a->[1] cmp $b->[1] ) }
				map [ $_->source, $_->method, $_ ],
			$db->types
		);
	}

	# a BigWigSet database
	elsif ( $db_ref =~ /BigWigSet/i ) {

		# get the metadata
		my $metadata = $db->metadata;

		# collect
		foreach my $file ( keys %{$metadata} ) {

			# get the type for each file
			my ( $primary, $type, $name );

			# get the appropriate tags
			foreach my $attribute ( keys %{ $metadata->{$file} } ) {
				if ( $attribute =~ m/^type/i ) {
					$type = $metadata->{$file}{$attribute};
				}
				elsif ( $attribute =~ m/name/i ) {
					$name = $metadata->{$file}{$attribute};
				}
				elsif ( $attribute =~ m/^primary_tag|method$/i ) {
					$primary = $metadata->{$file}{$attribute};
				}
			}

			# store the type
			push @types, $type || $primary || $name;
		}

		# sort the types in alphabetical order
		# and discard duplicate types - which may occur with stranded entries
		@types = uniq( sort { $a cmp $b } @types );
	}

	# some other database
	else {
		carp " no dataset lists for database type $db_ref!\n";
	}

	# finish
	return @types;
}

### Process and verify a dataset
sub verify_or_request_feature_types {

	# Retrieve passed values
	my %args = @_;    # the passed argument values as a hash reference

	# Check for single option
	$args{'single'} ||= 0;

	# Collect the datasets
	my @datasets;
	$args{'feature'} ||= undef;
	if ( ref $args{'feature'} eq 'ARRAY' ) {

		# an anonymous array of datasets
		@datasets = @{ $args{'feature'} };
	}
	elsif ( defined $args{'feature'} ) {
		push @datasets, $args{'feature'};
	}

	# else it is null and @datasets remains empty, prompting user input

	# feature limits
	# this is a regex to limit the feature types presented to the user
	# otherwise everything in the database is presented to the user
	my $limit = $args{'limit'} ||= undef;

	# Open database object and collect features
	$args{'db'} ||= undef;
	my $db = $args{'db'} ? open_db_connection( $args{'db'} ) : undef;

	# Initialize main output arrays
	my @good_datasets;
	my @bad_datasets;
	my %db_features;    # hash for features found in the database, use when needed

	# Check provided datasets
	if (@datasets) {

		# check for multiple comma-delimited datasets
		my @list_to_check;
		foreach my $item (@datasets) {
			if ( $item =~ /,/ ) {

				# this one has a comma, therefore it has more than dataset
				push @list_to_check, split( /,/, $item );
			}
			else {
				# a singleton
				push @list_to_check, $item;
			}
		}

		# now verify the datasets
		foreach my $dataset (@list_to_check) {

			# check for a remote file
			if ( $dataset =~ /^(?: http | ftp) .+ \. (?: bam | bw | bb) $/xi ) {

				# a remote file
				# assume it is good, no verification here though
				# it will either work or won't work
				push @good_datasets, $dataset;
			}

			# a local file
			elsif ( $dataset =~ /\.(?:bam|bw|bigwig|bb|bigbed|useq)$/i ) {

				# presume we have a local indexed data file

				# user may have requested two or more files to be merged
				# these should be combined with an ampersand
				# check each one
				my @files;
				foreach my $file ( split /\&/, $dataset ) {
					if ( -e $file ) {

						# file exists
						push @files, "file:$file";
					}
					else {
						# file doesn't exist! can't use this set of files
						@files = ();
						last;
					}
				}
				if (@files) {
					push @good_datasets, join( "&", @files );
				}
				else {
					push @bad_datasets, $dataset;
				}
			}

			# a feature type in a database
			else {
				# must be a database feature type

				# check for database features hash
				unless (%db_features) {
					if ($db) {

						# collect the database features as needed
						# I want both method:source and method in the hash
						foreach my $type ( get_dataset_list($db) ) {
							my ( $method, $source ) = split /:/, $type;
							if ($source) {
								$db_features{$type}   = 1;
								$db_features{$method} = 1;
							}
							else {
								$db_features{$type} = 1;
							}
						}

						# verify
						unless (%db_features) {
							carp " provided database has no feature types "
								. "to verify dataset(s) against!\n";
						}
					}
					else {
						# we need a database
						carp " unable to verify dataset without database";
					}
				}

				# validate the given dataset
				# user may have requested two or more features to be merged
				# these should be combined with an ampersand
				# check each one
				my $check = 0;
				foreach my $d ( split /&/, $dataset ) {

					# validate this dataset
					if ( exists $db_features{$d} ) {
						$check++;
					}
					else {
						# at least one of these is not good, fail check
						$check = 0;
						last;
					}
				}

				# report the verification
				if ($check) {
					push @good_datasets, $dataset;
				}
				else {
					push @bad_datasets, $dataset;
				}
			}
		}
	}

	# User must select datasets
	else {
		# dataset not specified
		# present the dataset list to the user and get an answer

		# get the dataset list
		if ($db) {

			# collect the database features as needed
			my $i = 1;
			foreach my $type ( get_dataset_list($db) ) {
				if ($limit) {
					my ( $p, $s ) = split /:/, $type;  # split into primary_tag and source
						# only keep those types that match our limiter
					next unless $p =~ /$limit/i;
				}
				$db_features{$i} = $type;
				$i++;
			}

			# verify
			unless (%db_features) {
				carp " provided database has no feature types " . "to collect!\n";
				return;
			}
		}
		else {
			# we need a database
			carp " no database provided from which to collect features!\n";
			return;
		}

		# present the list
		print "\n These are the available data sets in the database:\n";
		foreach ( sort { $a <=> $b } keys %db_features ) {

			# print out the list of microarray data sets
			print "  $_\t$db_features{$_}\n";
		}

		# prompt the user
		$args{'prompt'} ||= undef;
		if ( $args{'prompt'} ) {

			# provided custom prompt
			print $args{'prompt'};
		}
		else {
			# generic prompt
			if ( $args{'single'} ) {
				print " Enter one number for the data set or feature   ";
			}
			else {
				print " Enter the number(s) or range of the data set(s) or"
					. " feature(s)   ";
			}
		}

		# get answer from the user
		my $answer = <STDIN>;
		chomp $answer;
		my @answer_list = parse_list($answer);

		# take the first one if requested
		if ( $args{'single'} ) {
			unless ( scalar @answer_list == 1 ) {
				splice( @answer_list, 1 );
			}
		}

		# verify the answer list
		foreach my $answer (@answer_list) {

			# check for merged datasets
			if ( $answer =~ /&/ ) {

				# a merged dataset
				my @list  = split /&/, $answer;
				my $check = 1;

				# check all are good
				foreach (@list) {
					unless ( exists $db_features{$_} ) {
						$check = 0;
					}
				}

				# if all are good
				if ($check) {
					push @good_datasets, join( "&", map { $db_features{$_} } @list );
				}
				else {
					push @bad_datasets, $answer;
				}
			}

			else {
				# a single dataset
				# check if it is good

				if ( exists $db_features{$answer} ) {
					push @good_datasets, $db_features{$answer};
				}
				else {
					push @bad_datasets, $db_features{$answer};
				}
			}
		}
	}

	# Print bad results
	if (@bad_datasets) {
		print " The following datasets could not be verified:\n";
		foreach (@bad_datasets) {
			print "      $_\n";
		}
	}

	# Return good results
	if ( $args{'single'} ) {
		return $good_datasets[0];
	}
	else {
		return @good_datasets;
	}
}

### Process and verify a dataset
sub check_dataset_for_rpm_support {

	# get passed dataset and databases
	my $dataset = shift;
	my $cpu     = shift;

	# Calculate the total number of reads
	# this uses the global variable $rpkm_read_sum
	my $rpm_read_sum;

	if ( exists $TOTAL_READ_NUMBER{$dataset} ) {

		# this dataset has already been summed
		# no need to do it again
		$rpm_read_sum = $TOTAL_READ_NUMBER{$dataset};
	}

	elsif ( $dataset =~ /\.bam$/ ) {

		# a bam file dataset

		_load_bam_helper_module() unless $BAM_OK;
		if ($BAM_OK) {

			# Bio::ToolBox::db_helper::bam was loaded ok
			# sum the number of reads in the dataset
			$rpm_read_sum = sum_total_bam_alignments( $dataset, 0, 0, $cpu );
		}
		else {
			carp " Bam support is not available! "
				. "Is Bio::DB::Sam or Bio::DB::HTS installed?\n";
			return;
		}
	}

	elsif ( $dataset =~ /\.bb$/ ) {

		# a bigbed file dataset

		$BIGBED_OK = _load_helper_module('Bio::ToolBox::db_helper::bigbed')
			unless $BIGBED_OK;
		if ($BIGBED_OK) {

			# Bio::ToolBox::db_helper::bigbed was loaded ok
			# sum the number of features in the dataset
			$rpm_read_sum = sum_total_bigbed_features($dataset);
		}
		else {
			carp " BigBed support is not available! " . "Is Bio::DB::BigBed installed?\n";
			return;
		}
	}

	else {
		# some other non-supported dataset
		return;
	}

	# return the sum value if we've made it this far
	$TOTAL_READ_NUMBER{$dataset} = $rpm_read_sum;
	return $rpm_read_sum;
}

### Generate a new list of features
sub get_new_feature_list {

	# Retrieve passed values
	my %args = @_;    # the passed argument values

	# check data object
	my $data = $args{data} || undef;
	unless ($data) {
		confess "must pass a 'data' key and Bio::ToolBox::Data object!";
		return;
	}
	unless ( ref $data eq 'Bio::ToolBox::Data' ) {
		confess 'must pass a Bio::ToolBox::Data object!';
		return;
	}

	# Open a db connection
	$args{'db'} ||= undef;
	my $db = open_db_connection( $args{'db'} );
	unless ($db) {
		carp 'no database connected!';
		return;
	}

	# Verify a SeqFeature::Store database
	my $db_ref = ref $db;
	unless ( $db_ref =~ /^Bio::DB::SeqFeature::Store/ ) {
		carp "Database type $db_ref doesn't support generating feature lists!\n";
		return;
	}

	# Features to search for
	my $searchFeature = $args{'feature'} || $args{'features'} || undef;
	unless ($searchFeature) {
		carp "no search feature types passed!";
		return;
	}
	my @classes = split /,/, $searchFeature;    # it may or may not be a list

	# chromosomes to skip
	my $chr_exclude = $args{'chrskip'} || undef;

	# Add table columns
	my $pid_i  = $data->add_column('Primary_ID');
	my $name_i = $data->add_column('Name');
	my $type_i = $data->add_column('Type');

	# Set up the database iterator
	print "   Searching for $searchFeature features\n";
	my $iterator =
		$db->get_seq_stream( -types => scalar @classes ? \@classes : $searchFeature, );
	unless ($iterator) {

		# there should be some features found in the database
		carp "could not get feature iterator for database";
		return;
	}

	# Walk through the collected features
	my $total_count = 0;    # total found features
	while ( my $feature = $iterator->next_seq ) {
		$total_count++;

		# skip genes from excluded chromosomes
		next if ( defined $chr_exclude and $feature->seq_id =~ $chr_exclude );

		# Record the feature information
		# in the B::DB::SF::S database, the primary_ID is a number unique to the
		# the specific database, and is not portable between databases
		$data->add_row(
			[ $feature->primary_id, $feature->display_name, $feature->type, ] );
	}

	# print result of search
	printf "   Found %s features in the database\n", format_with_commas($total_count);
	printf "   Kept %s features.\n", format_with_commas( $data->last_row );

	# return the new data structure
	return 1;
}

### Generate a new list genomic windows
sub get_new_genome_list {

	# Collect the passed arguments
	my %args = @_;

	# check data object
	my $data = $args{data} || undef;
	unless ($data) {
		confess "must pass a 'data' key and Bio::ToolBox::Data object!";
	}
	unless ( ref $data eq 'Bio::ToolBox::Data' ) {
		confess 'must pass a Bio::ToolBox::Data object!';
	}

	# Open a db connection
	$args{'db'} ||= undef;
	my $db = open_db_connection( $args{'db'} );
	unless ($db) {
		carp 'no database connected!';
		return;
	}

	# Determine win and step sizes
	$args{'win'}  ||= 500;            # hard coded default size
	$args{'step'} ||= $args{'win'};

	# Prepare data structures
	my $chr_i   = $data->add_column('Chromosome');
	my $start_i = $data->add_column('Start');
	my $stop_i  = $data->add_column('Stop');
	$data->metadata( $start_i, 'win',  $args{'win'} );
	$data->metadata( $start_i, 'step', $args{'step'} );

	# Collect the chromosomes
	my $chr_exclude = $args{'skip'} || $args{'chrskip'} || undef;
	my @chromosomes = get_chromosome_list( $db, $chr_exclude );
	unless (@chromosomes) {
		carp " no sequence IDs were found in the database!\n";
		return;
	}

	# Load exclusion intervals
	my %exclusion_tree;
	$args{exclude} ||= $args{blacklist} || undef;
	if ( exists $args{exclude} and defined $args{exclude} ) {
		if ( ref( $args{exclude} ) eq 'Bio::ToolBox::Data' ) {
			if ( _load_helper_module('Set::IntervalTree') ) {

				# iterate through the data object of intervals
				# prepare an Interval Tree for each chromosome
				# and load the exclusion interval into it
				$args{exclude}->iterate(
					sub {
						my $row = shift;
						unless ( exists $exclusion_tree{ $row->seq_id } ) {
							$exclusion_tree{ $row->seq_id } =
								Set::IntervalTree->new();
						}
						$exclusion_tree{ $row->seq_id }
							->insert( 1, $row->start - 1, $row->end );
					}
				);
			}
			else {
				carp " Set::IntervalTree must be installed to use exclusion intervals!";
			}
		}
		else {
			confess " Exclusion data must be a Bio::ToolBox::Data object!";
		}
	}

	# Collect the genomic windows
	print "   Generating $args{win} bp windows in $args{step} bp increments\n";
	foreach (@chromosomes) {

		# name and length as sub-array in each element
		my ( $chr, $length ) = @{$_};
		my $Tree = $exclusion_tree{$chr} || undef;

		# iterate through chromosome
		for ( my $start = 1; $start <= $length; $start += $args{step} ) {
			my $end = $start + $args{win} - 1;
			if ( $end > $length ) {

				# fix end to the length of the chromosome
				$end = $length;
			}

			# skip if excluded
			next if ( $Tree and scalar( @{ $Tree->fetch( $start - 1, $end ) } ) >= 1 );
			$data->add_row( [ $chr, $start, $end ] );
		}
	}
	print "   Kept " . $data->{'last_row'} . " windows.\n";

	# Return the data structure
	return 1;
}

sub get_db_feature {
	my %args = @_;

	# Open a db connection
	$args{db} ||= undef;
	unless ( $args{db} ) {
		croak 'no database provided for getting a feature!';
	}
	my $db     = open_db_connection( $args{db} );
	my $db_ref = ref $db;

	# check for values and internal nulls
	$args{'id'}   ||= undef;
	$args{'name'} ||= undef;
	$args{'type'} ||= undef;
	undef $args{'name'} if $args{'name'} eq '.';
	undef $args{'id'}   if $args{'id'} eq '.';
	undef $args{'type'} if $args{'type'} eq '.';

	# quick method for feature retrieval
	if ( defined $args{'id'} and $db->can('fetch') ) {

		# we have a unique primary_id to identify the feature
		# usually this pulls out the feature directly
		my $feature = $db->fetch( $args{'id'} ) || undef;

		# check that the feature we pulled out is what we want
		my $check = $feature ? 1 : 0;
		if ($check) {
			$check = 0
				if ( defined $args{'name'} and $feature->display_name ne $args{'name'} );
			$check = 0 if ( defined $args{'type'} and $feature->type ne $args{'type'} );
		}

		# return if things match up as best we can tell
		if ($check) {
			return $feature;
		}
		else {
			# the primary_ids are out of date
			unless ($PRIMARY_ID_WARNING) {
				warn "CAUTION: Some primary IDs in Input file appear to be out of date\n";
				$PRIMARY_ID_WARNING++;
			}
		}
	}

	# otherwise use name and type
	return
		unless $args{'name'};   # otherwise db will return all features! Not what we want!
	my @features = $db->features(
		-name    => $args{'name'},    # we assume this name will be unique
		-aliases => 0,
		-type    => $args{'type'},
	);

	# if none are found
	unless (@features) {

		# try again with aliases enabled
		@features = $db->features(
			-name    => $args{'name'},
			-aliases => 1,
			-type    => $args{'type'},
		);
	}
	unless ( @features and $args{'name'} =~ /[;,\|]/ ) {

		# I used to append aliases to the end of the name in old versions of biotoolbox
		# crazy, I know, but just in case this happened, let's try breaking them apart
		my $name2 = ( split( /\s*[;,\|]\s*/, $args{'name'} ) )[0];

		# multiple names present using common delimiters ;,|
		# take the first name only, assume others are aliases that we don't need
		@features = $db->features(
			-name    => $name2,
			-aliases => 1,
			-type    => $args{'type'},
		);
	}

	# check the number of features returned
	if ( scalar @features > 1 ) {

		# there should only be one feature found
		# if more are returned, there's redundant or duplicated data in the db

		# first check whether addition of aliases may improve accuracy
		if ( $args{'name'} =~ /;/ ) {

			# we have presumed aliases in the name value
			my $check = $args{'name'};
			$check =~ s/\s*;\s*/;/g;    # strip spaces just in case

			# check each feature and try to match name and aliases
			my @candidates;
			foreach my $f (@features) {
				my $f_name =
					join( ';', $f->display_name, ( $f->get_tag_values('Alias') ) );
				if ( $check eq $f_name ) {

					# found one
					push @candidates, $f;
				}
			}

			if ( scalar @candidates == 1 ) {

				# we found it!
				return shift @candidates;
			}
			elsif ( scalar @candidates > 1 ) {

				# hopefully we've improved a little bit?
				@features = @candidates;
			}
		}

		# warn the user, this should be fixed
		printf "  Found %s %s features named '%s' in the database! Using first one.\n",
			scalar(@features), $args{'type'}, $args{'name'};
	}
	elsif ( !@features ) {
		printf "  Found no %s features named '%s' in the database!\n", $args{'type'},
			$args{'name'};
		return;
	}

	# done
	return shift @features;
}

### Get a dataset score for a single region
sub get_segment_score {

	# parameters passed as an array
	# we will be passing this array on as a reference to the appropriate
	# imported helper subroutine
	# chromosome, start, stop, strand, strandedness, method, return type, db, dataset
	confess "incorrect number of parameters passed!" unless scalar @_ == 9;

	# check the database
	$_[DB] = open_db_connection( $_[DB] ) if ( $_[DB] and not ref( $_[DB] ) );

	# check for combined datasets
	if ( $_[DATA] =~ /&/ ) {
		push @_, ( split '&', pop @_ );
	}

	# determine method
	# yes, we're only checking the first dataset, but they should all
	# be the same type
	my $db_method = $DB_METHODS{ $_[METH] }{ $_[RETT] }{ $_[DB] }{ $_[DATA] }
		|| _lookup_db_method( \@_ );

	# return type values
	# 0 = calculate score
	# 1 = score array
	# 2 = hash position scores
	# immediately return calculated score if appropriate
	if ( $_[RETT] > 0 ) {

		# immediately return either indexed hash or array of scores
		return &{$db_method}( \@_ );
	}
	else {
		# calculate a score
		my $scores = &{$db_method}( \@_ );

		# this might be an array reference of scores that need to be combined
		# or it could be a single scalar score which just needs to be returned
		if ( ref $scores ) {
			return calculate_score( $_[METH], $scores );
		}
		else {
			if ( not defined $scores ) {
				return 0 if $_[METH] =~ /count|sum/;
				return '.';
			}
			return $scores;
		}
	}
}

sub calculate_score {
	my ( $method, $scores ) = @_;
	$scores ||= [];    # just in case
	if ( exists $SCORE_CALCULATOR_SUB{$method} ) {
		return &{ $SCORE_CALCULATOR_SUB{$method} }($scores);
	}
	else {
		confess " unrecognized method '$method'!";
	}
}

sub get_chromosome_list {

	# options
	my $database    = shift;
	my $chr_exclude = shift || undef;

	# Open a db connection
	my $db = open_db_connection($database);
	unless ($db) {
		carp 'no database connected!';
		return;
	}

	# Check for BigWigSet database
	# these need to be handled a little differently
	if ( ref $db =~ /BigWigSet/ ) {

		# BigWigSet databases are the only databases that don't
		# support the seq_ids method
		# instead we have to look at one of the bigwigs in the set
		my $bw_file = ( $db->bigwigs )[0];
		$db = open_db_connection($bw_file);
	}

	# Collect chromosome lengths
	# we want to maintain the original order of the chromosomes so we
	# won't be putting it into a hash
	# instead an array of arrays
	my @chrom_lengths;

	# Database specific approaches to collecting chromosomes
	# I used to have one generic approach, but idiosyncrasies and potential
	# bugs make me use different approaches for better consistency
	my $type = ref $db;

	# SeqFeature::Store
	if ( $type =~ /^Bio::DB::SeqFeature::Store/ ) {
		for my $chr ( $db->seq_ids ) {

			# check for excluded chromosomes
			next if ( defined $chr_exclude and $chr =~ /$chr_exclude/i );

			# get chromosome size
			my ($seqf) = $db->get_features_by_name($chr);
			my $length = $seqf ? $seqf->length : 0;

			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}

	# libBigWig Bigfile, including both bigWig and bigBed
	elsif ( $type eq 'Bio::DB::Big::File' ) {

		# this doesn't follow the typical BioPerl convention
		# it's a hash, so randomly sorted!!!!! Never will be in same order as file!!!!
		my $chroms = $db->chroms();

		my @list;
		foreach ( values %{ $chroms } ) {

			# check for excluded chromosomes
			next if ( defined $chr_exclude and $_->{name} =~ /$chr_exclude/i );
			push @list, [ $_->{name}, $_->{length} ];
		}

		# sort consistently by a sane system
		@chrom_lengths = sane_chromo_sort(@list);
	}

	# UCSC kent Bigfile
	elsif ( $type eq 'Bio::DB::BigWig' or $type eq 'Bio::DB::BigBed' ) {
		foreach my $chr ( $db->seq_ids ) {

			# check for excluded chromosomes
			next if ( defined $chr_exclude and $chr =~ /$chr_exclude/i );

			# get chromosome size
			my $length = $db->length($chr);

			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}

	# Samtools Bam
	elsif ( $type eq 'Bio::DB::Sam' or $type eq 'Bio::DB::HTS' ) {
		for my $tid ( 0 .. $db->n_targets - 1 ) {

			# each chromosome is internally represented in the bam file as
			# a numeric target identifier
			# we can easily convert this to an actual sequence name
			# we will force the conversion to go one chromosome at a time

			# sequence info
			my $chr    = $db->target_name($tid);
			my $length = $db->target_len($tid);

			# check for excluded chromosomes
			next if ( defined $chr_exclude and $chr =~ /$chr_exclude/i );

			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}

	# Fast through modern HTS fai index
	elsif ( $type eq 'Bio::DB::HTS::Faidx' ) {
		for my $chr ( $db->get_all_sequence_ids ) {

			# check for excluded
			next if ( defined $chr_exclude and $chr =~ /$chr_exclude/i );

			# get length and store it
			my $length = $db->length($chr);
			push @chrom_lengths, [ $chr, $length ];
		}
	}

	# Fasta through old BioPerl adapter
	elsif ( $type eq 'Bio::DB::Fasta' ) {
		for my $chr ( $db->get_all_ids ) {

			# check for excluded chromosomes
			next if ( defined $chr_exclude and $chr =~ /$chr_exclude/i );

			# get chromosome size
			my $seq    = $db->get_Seq_by_id($chr);
			my $length = $seq ? $seq->length : 0;

			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}

	# other Bioperl adapter?????
	elsif ( $db->can('seq_ids') ) {
		foreach my $chr ( $db->seq_ids ) {

			# check for excluded chromosomes
			next if ( defined $chr_exclude and $chr =~ /$chr_exclude/i );

			# generate a segment representing the chromosome
			# due to fuzzy name matching, we may get more than one back
			my @segments = $db->segment($chr);

			# need to find the right one
			my $segment;
			while (@segments) {
				$segment = shift @segments;
				last if $segment->seq_id eq $chr;
			}

			# check segment
			unless ($segment) {
				carp " No genome segment for seq_id $chr!!?? skipping\n";
				next;
			}

			# get the chromosome length
			my $length = $segment->length;

			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}

	else {
		carp
" Unable to identify chromosomes in database type '$type'. Try using a different database file or adapter\n";
		return;
	}

	# Return
	unless (@chrom_lengths) {
		carp " no chromosome sequences identified in database!\n";
		return;
	}
	return @chrom_lengths;
}

sub low_level_bam_fetch {
	confess "incorrect number of parameters passed!" unless scalar @_ == 6;
	my ( $sam, $tid, $start, $stop, $callback, $data ) = @_;

	# run the the low level bam fetch based on which adapter is being used
	unless ($BAM_ADAPTER) {
		$BAM_ADAPTER = ref $sam =~ /hts/i ? 'hts' : 'sam';
	}
	if ( $BAM_ADAPTER eq 'hts' ) {

		# using Bio::DB::HTS
		my $index = $sam->hts_index;
		return unless $index;
		return $index->fetch( $sam->hts_file, $tid, $start, $stop, $callback, $data );
	}
	elsif ( $BAM_ADAPTER eq 'sam' ) {

		# using Bio::DB::Sam
		my $index = $sam->bam_index;
		return unless $index;
		return $index->fetch( $sam->bam, $tid, $start, $stop, $callback, $data );
	}
	else {
		confess "no bam adapter loaded!\n";
	}
}

sub low_level_bam_coverage {
	confess "incorrect number of parameters passed!" unless scalar @_ == 4;
	my ( $sam, $tid, $start, $stop ) = @_;

	# run the the low level bam coverage based on which adapter is being used
	unless ($BAM_ADAPTER) {
		$BAM_ADAPTER = ref $sam =~ /hts/i ? 'hts' : 'sam';
	}
	if ( $BAM_ADAPTER eq 'hts' ) {

		# using Bio::DB::HTS
		my $index = $sam->hts_index;
		return unless $index;
		return $index->coverage( $sam->hts_file, $tid, $start, $stop );
	}
	elsif ( $BAM_ADAPTER eq 'sam' ) {

		# using Bio::DB::Sam
		my $index = $sam->bam_index;
		return unless $index;
		return $index->coverage( $sam->bam, $tid, $start, $stop );
	}
	else {
		confess "no bam adapter loaded!\n";
	}
}

sub get_genomic_sequence {
	confess "incorrect number of parameters passed!" unless scalar @_ == 4;
	my ( $db, $chrom, $start, $stop ) = @_;

	# check database
	my $type = ref $db;
	$db = open_db_connection($db) if not $type;

	# return sequence based on the type of database adapter we're using
	if ( $type eq 'Bio::DB::HTS::Faidx' ) {
		return $db->get_sequence_no_length("$chrom:$start-$stop");
	}
	elsif ( $type eq 'Bio::DB::Sam::Fai' ) {
		return $db->fetch("$chrom:$start-$stop");
	}
	elsif ( $db->can('seq') ) {

		# BioPerl database including Bio::DB::Fasta and Bio::DB::SeqFeature::Store
		return $db->seq( $chrom, $start, $stop );
	}
	else {
		confess "unsupported database $type for collecting genomic sequence!\n";
	}
}

### deprecated methods
# just in case
sub validate_included_feature {
	confess "validate_included_feature() is no longer a valid method. "
		. "Please update your script to the current API.\n";
}

### Internal subroutine to retrieve the scores from an established region object
sub _lookup_db_method {

	# parameters passed as an array reference
	my $param = shift;

	# determine the appropriate score method
	my $score_method;
	if ( $param->[DATA] =~ /^file|http|ftp/ ) {

		# collect the data according to file type

		# BigWig Data file
		if ( $param->[DATA] =~ /\.(?:bw|bigwig)$/i ) {

			# file is in bigwig format
			# get the dataset scores using Bio::ToolBox::db_helper::bigwig
			# this uses either Bio::DB::BigWig or Bio::DB::Big

			# check that we have bigwig support
			$BIGWIG_OK = _load_bigwig_helper_module() unless $BIGWIG_OK;
			if ($BIGWIG_OK) {
				if ( $param->[RETT] == 2 ) {
					$score_method = \&collect_bigwig_position_scores;
				}
				elsif ( $param->[RETT] == 1 ) {
					$score_method = \&collect_bigwig_scores;
				}
				elsif ( $param->[METH] =~ /min|max|mean/ ) {
					$score_method = \&collect_bigwig_score;
				}
				elsif ( $param->[METH] =~ /sum|count/ ) {

					# difference between ucsc and libBigWig libraries
					$score_method =
						$BIG_ADAPTER eq 'ucsc'
						? \&collect_bigwig_score
						: \&collect_bigwig_scores;
				}
				else {
					$score_method = \&collect_bigwig_scores;
				}
			}
			else {
				croak
" BigWig support is not enabled! Is Bio::DB::Big or Bio::DB::BigFile installed?";
			}
		}

		# BigBed Data file
		elsif ( $param->[DATA] =~ /\.(?:bb|bigbed)$/i ) {

			# data is in bigbed format
			# get the dataset scores using Bio::ToolBox::db_helper::bigbed
			# this uses either Bio::DB::BigBed or Bio::DB::Big

			# check that we have bigbed support
			$BIGBED_OK = _load_bigbed_helper_module() unless $BIGBED_OK;
			if ($BIGBED_OK) {
				if ( $param->[RETT] == 2 ) {
					$score_method = \&collect_bigbed_position_scores;
				}
				else {
					$score_method = \&collect_bigbed_scores;
				}
			}
			else {
				croak
" BigBed support is not enabled! Is Bio::DB::Big or Bio::DB::BigFile installed?";
			}
		}

		# BAM data file
		elsif ( $param->[DATA] =~ /\.bam$/i ) {

			# data is in bam format
			# get the dataset scores using Bio::ToolBox::db_helper::bam
			# this uses the Bio::DB::Sam or Bio::DB::HTS adaptor

			# check that we have Bam support
			_load_bam_helper_module() unless $BAM_OK;
			if ($BAM_OK) {
				$score_method = \&collect_bam_scores;
			}
			else {
				croak " Bam support is not enabled! "
					. "Is Bio::DB::HTS or Bio::DB::Sam installed?\n";
			}
		}

		# USeq Data file
		elsif ( $param->[DATA] =~ /\.useq$/i ) {

			# data is in useq format
			# get the dataset scores using Bio::ToolBox::db_helper::useq
			# this uses the Bio::DB::USeq adaptor

			# check that we have bigbed support
			$USEQ_OK = _load_helper_module('Bio::ToolBox::db_helper::useq')
				unless $USEQ_OK;
			if ($USEQ_OK) {
				if ( $param->[RETT] == 2 ) {
					$score_method = \&collect_useq_position_scores;
				}
				else {
					$score_method = \&collect_useq_scores;
				}
			}
			else {
				croak " USeq support is not enabled! " . "Is Bio::DB::USeq installed?\n";
			}
		}

		# Unsupported Data file
		else {
			confess sprintf " Unsupported or unrecognized file type for file '%s'!\n",
				$param->[DATA];
		}
	}

	### BigWigSet database
	elsif ( ref( $param->[DB] ) =~ m/BigWigSet/ ) {

		# calling features from a BigWigSet database object
		# this uses either Bio::DB::BigWig or Bio::DB::Big

		# check that we have bigwig support
		# duh! we should, we probably opened the stupid database!
		$BIGWIG_OK = _load_bigwig_helper_module() unless $BIGWIG_OK;
		croak
" BigWigSet support is not enabled! Is Bio::DB::Big or Bio::DB::BigFile installed?"
			unless $BIGWIG_OK;

		# the data collection depends on the method
		if ( $param->[RETT] == 2 ) {
			$score_method = \&collect_bigwigset_position_scores;
		}
		elsif ( $param->[RETT] == 1 ) {
			$score_method = \&collect_bigwigset_scores;
		}
		elsif ( $param->[METH] =~ /min|max|mean|sum|count/ ) {

			# difference between ucsc and libBigWig libraries
			# the ucsc library works with summaries and we can handle multiple of these
			# but the big adapter doesn't
			$score_method =
				$BIG_ADAPTER eq 'ucsc'
				? \&collect_bigwigset_score
				: \&collect_bigwigset_scores;
		}
		else {
			$score_method = \&collect_bigwigset_scores;
		}
	}

	### BioPerl style database
	elsif ( ref( $param->[DB] ) =~ m/^Bio::DB/ ) {

		# a BioPerl style database, including Bio::DB::SeqFeature::Store
		# most or all Bio::DB databases support generic get_seq_stream() methods
		# that return seqfeature objects, which we can use in a generic sense

		# check that we have support
		# duh! we should, we probably opened the stupid database!
		$SEQFASTA_OK = _load_helper_module('Bio::ToolBox::db_helper::seqfasta')
			unless $SEQFASTA_OK;
		if ($SEQFASTA_OK) {

			# get the dataset scores using Bio::ToolBox::db_helper::seqfasta
			# check that we support methods
			unless ( $param->[DB]->can('get_seq_stream') ) {
				confess sprintf
"unsupported database! cannot use %s as it does not support get_seq_stream method or equivalent",
					ref( $param->[DB] );
			}
			$score_method = \&collect_store_scores;
		}
		else {
			croak " SeqFeature Store support is not enabled! "
				. "Is BioPerl and Bio::DB::SeqFeature::Store properly installed?\n";
		}
	}

	### Some other database?
	else {
		confess "no recognizeable dataset provided!" unless $param->[DATA];
		confess "no database passed!"                unless $param->[DB];
		confess sprintf "something went wrong! parameters: %s", join ', ', @{ $param };
	}

	### Cache and return the results
	$DB_METHODS{ $param->[METH] }{ $param->[RETT] }{ $param->[DB] }{ $param->[DATA] } =
		$score_method;
	return $score_method;
}

sub _load_helper_module {
	my $class   = shift;
	my $success = 0;
	eval {
		# using Module::Load to dynamically load modules and import default functions
		autoload $class;
		$success = 1;
	};
	return $success;
}

sub _load_bam_helper_module {
	if ($BAM_ADAPTER) {
		if ( $BAM_ADAPTER =~ /sam/i ) {
			$BAM_OK      = _load_helper_module('Bio::ToolBox::db_helper::bam');
			$BAM_ADAPTER = 'sam';    # for internal consistency
		}
		elsif ( $BAM_ADAPTER =~ /hts/i ) {
			$BAM_OK      = _load_helper_module('Bio::ToolBox::db_helper::hts');
			$BAM_ADAPTER = 'hts';    # for internal consistency
		}
		elsif ( $BAM_ADAPTER =~ /none/i ) {

			# basically for testing purposes, don't use a module
			return 0;
		}
		else {
			# unrecognized
			$@ = 'unrecognized';
		}
		if ( not $BAM_OK ) {
			print "Requested '$BAM_ADAPTER' adapter could not be loaded: $@\n";
		}
	}
	else {
		# try hts first, then sam
		# be sure to set BAM_ADAPTER upon success
		$BAM_OK = _load_helper_module('Bio::ToolBox::db_helper::hts');
		if ($BAM_OK) {
			$BAM_ADAPTER = 'hts';
		}
		else {
			$BAM_OK      = _load_helper_module('Bio::ToolBox::db_helper::bam');
			$BAM_ADAPTER = 'sam' if $BAM_OK;
		}
	}
	return $BAM_OK;
}

sub _load_bigwig_helper_module {
	if ($BIG_ADAPTER) {
		if ( $BIG_ADAPTER =~ /ucsc|kent/i ) {
			$BIGWIG_OK   = _load_helper_module('Bio::ToolBox::db_helper::bigwig');
			$BIG_ADAPTER = 'ucsc';    # for internal consistency
		}
		elsif ( $BIG_ADAPTER =~ /big/i ) {
			$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::big');
			$BIGBED_OK = $BIGWIG_OK;                                         # bigbed too!
			$BIG_ADAPTER = 'big';    # for internal consistency
		}
		elsif ( $BIG_ADAPTER =~ /none/i ) {

			# basically for testing purposes, don't use a module
			return 0;
		}
		else {
			# unrecognized
			$@ = 'unrecognized';
		}
		if ( not $BIGWIG_OK ) {
			print "Requested '$BIG_ADAPTER' adapter could not be loaded: $@\n";
		}
	}
	else {
		# we have to try each one out
		# try the modern big adapter first
		$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::big');
		if ($BIGWIG_OK) {

			# success!
			$BIGBED_OK   = $BIGWIG_OK;            # bigbed too!
			$BIG_ADAPTER = 'big' if $BIGWIG_OK;
		}
		else {
			$BIGWIG_OK   = _load_helper_module('Bio::ToolBox::db_helper::bigwig');
			$BIG_ADAPTER = 'ucsc' if $BIGWIG_OK;
		}
	}
	return $BIGWIG_OK;
}

sub _load_bigbed_helper_module {
	if ($BIG_ADAPTER) {
		if ( $BIG_ADAPTER =~ /ucsc|kent/i ) {
			$BIGBED_OK   = _load_helper_module('Bio::ToolBox::db_helper::bigbed');
			$BIG_ADAPTER = 'ucsc';    # for internal consistency
		}
		elsif ( $BIG_ADAPTER =~ /big/i ) {
			$BIGBED_OK = _load_helper_module('Bio::ToolBox::db_helper::big');
			$BIGWIG_OK = $BIGBED_OK;                                         # bigwig too!
			$BIG_ADAPTER = 'big';    # for internal consistency
		}
		elsif ( $BIG_ADAPTER =~ /none/i ) {

			# basically for testing purposes, don't use a module
			return 0;
		}
		elsif ( $BAM_ADAPTER =~ /\w+/ ) {

			# unrecognized
			$@ = 'unrecognized';
		}
		if ( not $BIGWIG_OK ) {
			print "Requested '$BIG_ADAPTER' adapter could not be loaded: $@\n";
		}
	}
	else {
		# we have to try each one out
		# try the modern big adapter first
		$BIGBED_OK = _load_helper_module('Bio::ToolBox::db_helper::big');
		if ($BIGBED_OK) {

			# success!
			$BIGWIG_OK   = $BIGBED_OK;            # bigwig too!
			$BIG_ADAPTER = 'big' if $BIGBED_OK;
		}
		else {
			$BIGBED_OK   = _load_helper_module('Bio::ToolBox::db_helper::bigbed');
			$BIG_ADAPTER = 'ucsc' if $BIGBED_OK;
		}
	}
	return $BIGBED_OK;
}

1;

__END__

=head1 NAME

Bio::ToolBox::db_helper - helper interface to various database formats

=head1 DESCRIPTION

In most cases, this module does not need to be used directly. The 
methods available to L<Bio::ToolBox::Data> and L<Bio::ToolBox::Data::Feature> 
provide convenient access to the methods described here.

These are helper subroutines to work with relevant databases that can be 
accessed through BioPerl modules. These include the L<Bio::DB::SeqFeature::Store> 
relational database as well as L<Bio::DB::Fasta>, L<Bio::DB::Sam>, 
L<Bio::DB::HTS>, L<Bio::DB::BigWig>, L<Bio::DB::BigWigSet>, and 
L<Bio::DB::BigBed> databases. The primary functions 
included opening connections to these databases, verifying features found 
within the database, generating lists of features or genomic intervals for 
subsequent analysis, and collecting features and/or scores within regions 
of interest. Collected scores may be summarized using a variety of statistical 
methods.

When collecting scores or features, the data may hosted in a variety of 
formats and locations. These include:

=over 4

=item SeqFeature database

Genomic feature annotation, including genes, transcripts, exons, etc, 
may be stored in a L<Bio::DB::SeqFeature::Store> database. These 
databases are backed by either a relational database, including 
MySQL, PostGreSQL, or SQLite. Small GFF3 files may also be loaded in 
an in-memory database. These are typically loaded from GFF3 
files; see the adapter documentation for more information. Dataset 
scores, such as from microarray, may also be stored as the score 
value in the source GFF file. 

References to local, binary, indexed files may also be included as 
attributes to features stored in the database. This is legacy support 
from the L<GBrowse genome browser|http://gmod.org/gbrowse>, and should 
not be used in preference to other methods. Supported files include the 
legacy binary wig files (F<.wib>, see L<Bio::Graphics::Wiggle>) using the 
C<wigfile> attribute, or bigWig files using the C<wigfile> or C<bigwigfile> 
attribute. The attribute values must be full paths. 

=item BigWig files

BigWig files (F<.bw> and F<.bigwig>) are compressed, binary, indexed 
versions of text wig files and may be accessed either locally or remotely. 
They support extremely  fast score retrieval from any genomic location of 
any size without sacrificing resolution (spatial and numeric). See
L<bigWig description|http://genome.ucsc.edu/goldenPath/help/bigWig.html> for more information. 
BigWig files are supported by either the L<Bio::DB::Big> or the 
L<Bio::DB::BigWig> adapter, based on either L<libBigWig|https://github.com/dpryan79/libBigWig> 
or the L<UCSC|http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads> libraries, 
respectively; see their respective documentation for more information.

=item Directory of BigWig files

A directory containing two or more BigWig files is assembled into a 
BigWigSet, allowing for metadata, such as strand, to be associated with 
BigWig files. Additional metadata beyond the filename may be included in 
text file F<metadata.txt> within the directory. See the 
L<Bio::DB::BigWigSet> adapter documentation for more information. When using  
L<Bio::DB::Big> adapters, BigWigSet support is natively supported by 
the BioToolBox package.

=item BigBed files

BigBed files are compressed, binary, indexed versions of text BED files. See
L<bigBed description|http://genome.ucsc.edu/goldenPath/help/bigBed.html> for more information.
Both local and remote files may be accessed. BigBed files are supported by 
either the L<Bio::DB::Big> or the L<Bio::DB::BigBed> adapter, based on either L<libBigWig|https://github.com/dpryan79/libBigWig> 
or the L<UCSC|http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads> libraries, 
respectively; see their respective documentation for more information.

=item Bam files

Bam files (F<.bam>) are compressed, binary, indexed versions of the text SAM file, 
or sequence alignment map. They are used with next generation sequencing 
technologies. They support individual alignment retrieval as well as 
read depth coverage. Two different Bam file adapters are supported. 

The L<Bio::DB::HTS> adapter is an interface to the Bam (or Cram) file. 
This is based on the modern HTSlib C library (version E<gt>= 1.0), successor 
to the original samtools library. See L<http://samtools.github.io> for 
more information.

The L<Bio::DB::Sam> adapter is an older interface to the Bam alignment 
file. This is based on the samtools C library version E<lt>= 0.1.19. While 
this adapter is still supported, the above should be used for new installs.

=item Cram files

Cram files (F<.cram>) are similar to Bam files, but are much smaller due to 
only storing sequence differences for each alignment. As such, they require 
an indexed, reference fasta to regenerate the complete alignment. Cram files 
are only supported by the L<Bio::DB::HTS> adapter.

=item USeq files

USeq files (F<.useq>) are compressed, binary, indexed files that support 
multiple information types, including region intervals, BED type 
annotations, or wig type scores distributed across the genome. They 
support rapid, random access across the genome and are comparable to 
both BigWig and BigBed files. See L<http://useq.sourceforge.net/useqArchiveFormat.html> 
for more information. Files may only be local. These files are supported 
by the L<Bio::DB::USeq> adapter. 

=item Fasta files

Fasta files (F<.fa> or F<.fasta>) may be opened. Fasta files are indexed to 
rapidly and randomly fetch genomic sequence. Three different adapters are 
available for indexing the fasta file. 

The L<Bio::DB::HTS::Faidx> adapter is preferentially used, and requires a 
F<.fai> index file. 

The L<Bio::DB::Sam::Fai> adapter may alternatively used, and also requires a 
F<.fai> index file.

The older style L<Bio::DB::Fasta> adapter is much slower, but will index 
either a single genomic fasta file or a directory of individual chromosome 
or contig fasta files. 

Additionally, genomic sequence stored in a L<Bio::DB::SeqFeature::Store> 
annotation database may also be used.

=back

While the functions within this library may appear to be simply a rehashing of 
the methods and functions in L<Bio::DB::SeqFeature::Store> and other modules, they 
either provide a simpler function to often used database methodologies or are 
designed to work intimately with the BioToolBox data format file and data structures 
(see L<Bio::ToolBox::file_helper>). One key advantage to these functions is the ability
to work with datasets that are stranded (transcriptome data, for example).

Historically, this module was initially written to use L<Bio::DB::GFF> for 
database usage. It has since been re-written to use L<Bio::DB::SeqFeature::Store> 
database schema. Additional functionality has also been added for other databases 
and file formats. 

Complete usage and examples for the functions are provided below.

=head1 USAGE

Call the module at the beginning of your perl script and include the module 
names to export. None are exported by default.

  use Bio::ToolBox::db_helper qw(
	  get_new_feature_list 
	  get_db_feature
  );

This will export the indicated subroutine names into the current namespace. 
Their usage is detailed below. 

=head1 EXPORTED SUBROUTINES

=over 4

=item C<$BAM_ADAPTER> variable

=item C<$BIG_ADAPTER> variable

These variables control which Bam and bigWig/bigBed adapters are used
during execution in installations where both adapters are present. These
variables are not set initially; adapters are chosen automatically when an
adapter is requested during execution and modules are evaluated for loading. 

=item use_bam_adapter

=item use_big_adapter

These are convenience methods for explicitly setting and retrieving the 
C<$BAM_ADAPTER> and C<$BIG_ADAPTER> variables, respectively. 

Optionally pass the value to set. It will always return the value being used. 
Values include the following:

=over 4

=item Bam adapter

   * hts     Bio::DB::HTS (default if available)
   * sam     Bio::DB::Sam
   * none    no adapters allowed

=item Big adapter

   * big     Bio::DB::Big (default if available)
   * ucsc    Bio::DB::BigWig and Bio::DB::BigBed
   * none    no adapters allowed

=back

=item open_db_connection

	my $db_name = 'cerevisiae';
	my $db = open_db_connection($db_name);
	
	my $file = '/path/to/file.bam';
	my $db = open_db_connection($file);
	
	# within a forked child process
	# reusing the same variable to simplify code
	# pass second true value to avoid cache
	$db = open_db_connection($file, 1); 

This is a simplified, generalized method that will open a connection to a 
database or indexed file using any one of a number of different 
BioPerl-style adapters. For local and remote files, the appropriate 
adapter is determined by the file extension. Relational database 
connection and authentication information is checked in a configuration 
file. Once identified, the appropriate Perl module adapter is loaded 
during run time automatically, saving the user from knowing in advance  
to C<use> the appropriate module in the script. 

Pass the name of a relational database or the path of the database file to 
the subroutine. The opened database object is returned. If it fails, then 
an error message should be generated and nothing is returned.

B<Important!> If you are forking your perl process, B<always> re-open your 
database objects in the child process, and pass a second true value 
to avoid using the cached database object. By default, opened databases are 
cached to improve efficiency, but this will be disastrous when crossing forks. 

=over 4

=item SeqFeature Store database

Provide the name of the relational database to open. Parameters for 
connecting to a relational database are stored in the BioToolBox 
configuration file, F<.biotoolbox.cfg>. These include database adaptors, 
user name, password, etc. Further information regarding the configuration 
file may be found in L<Bio::ToolBox::db_helper::config>.

For SQLite databases, provide the path to the F<.sqlite> or F<.db> file.

For in-memory databases, provide the path to the F<.gff> or F<.gff3> file. 
The larger the file, the more memory and processing time is consumed as 
the file is parsed into memory. This may be fine for small model organisms 
like yeast, but not for vertebrate annotation.  

=item Bam file database 

Provide the path to the F<.bam> file. This may be a local file, or a remote 
URL (generally supported by the adapter). If a local F<.bam.bai> index file 
is not present, it will automatically be built prior to opening; this may 
fail if the bam file is not sorted. Remote files are not indexed.

=item Cram file database

Provide the path to the local F<.cram> file. Currently, supplementary fasta 
files are not supported, so the Cram file header must point to a valid 
reference. Cram files will be indexed automatically if an index is not available.

=item BigWig file database

Provide the path to a local F<.bw> or F<.bigwig> file, or URL to a 
remote file. 

=item BigWigSet database

A local or remote directory of one or more BigWig files that can treated 
collectively as a database. A special text file may be included in the 
directory to assign metadata to each BigWig file, including attributes such 
as type, source, display name, strand, etc. See L<Bio::DB::BigWigSet> for 
more information on the formatting of the metadata file. If a metadata 
file is not included (it's not necessary), then the filenames are used 
as the name and type. Strand may be parsed from the filenames if the 
basename ends in F<.f> and F<.r>.

=item BigBed database

Provide the path to a local F<.bb> or F<.bigbed> file, or URL to a 
remote file. 

=item USeq database

Provide the path to a local F<.bb> or F<.bigbed> file, or URL to a 
remote file. 

=item Fasta database

Provide the path to a local fasta file for rapidly fetching genomic 
sequence. If the fasta file is not indexed, then it will be indexed 
automatically upon opening. 

=back

=item get_db_name

  my $name = get_db_name($db);

This method will attempt to get the name of an opened Database 
object if, for some reason, it's unknown, e.g. the user only provided an 
opened db object. This only works for some databases, at least those I've 
bothered to track down and find a usable API call to use.

=item get_dataset_list

	my $db_name = 'cerevisiae';
	my @types = get_dataset_list($db_name);

This method will retrieve a list of the available features stored in the 
database and returns an array of the features' types.

Pass either the name of the database or an established database object. 
Supported databases include both L<Bio::DB::SeqFeature::Store> and 
L<Bio::DB::BigWigSet> databases. 

For L<Bio::DB::SeqFeature::Store> databases, the type is represented as 
"C<type:source>", corresponding to the third and second GFF columns, 
respectively. The types are sorted alphabetically first by source, 
then by method.

For L<Bio::DB::BigWigSet> databases, the C<type>, C<primary_tag>, C<method>, 
or C<display_name> attribute may be used, in that respective order of 
availability. The list is sorted alphabetically.

=item verify_or_request_feature_types

	# verify a dataset, either file or database type
	my $db_name = 'cerevisiae';
	my $dataset = 'microaray_data';
	my @features = verify_or_request_feature_types(
		db      => $db_name,
		feature => $dataset,
	);
	unless (@features) {
		die " no valid feature provided!\n";
	}
	
	# select a list of features
	my $db_name = 'cerevisiae';
	my @types = (); # user not provided
	@types = verify_or_request_feature_types(
		db      => $db_name,
		feature => $types,
	);
	# user will be promoted to select from database list
	if (@types) {
		print "using types " . join(", ", @types);
	}
	else {
		die " no valid features selected!\n";
	}

This subroutine will process a list of feature types or data sources to be 
used for data or feature collection. There are two modes of action. If a 
list was provided, it will process and verify the list. If no list was 
provided, then the user will be prompted to select from a list of available 
feature types in the provided database.

The provided list may be a list of feature types or a list of single 
file data sources, including Bam, bigWig, or bigBed files. If the list 
includes feature types, they are verified as existing in the provided 
database. If the list includes file names, the local files are checked 
and the names prefixed with "C<file:>" for use downstream.

If no list was provided, then a list of available feature types in the 
provided database will be presented to the user, and the user prompted 
to make a selection. One or more types may be selected, and a single 
item may be enforced if requested. The response is filtered through 
the parse_list method from L<Bio::ToolBox::utility>, so a mix of single 
numbers or a range of numbers may be accepted. The responses are then 
validated.

Two types of databases are supported: L<Bio::DB::SeqFeature::Store> and 
L<Bio::DB::BigWigSet> databases. For SeqFeature Store databases, the 
type is comprised of "C<method:source>", corresponding to the third and 
second GFF columns. For BigWigSet databases, types correspond to either 
the C<type>, C<method>, C<primary_tag>, or C<display_name> attributes, 
in that order of availability. 

For feature types or files that pass validation, they are returned as a 
list. Types of files that do not pass validation are printed to C<STDOUT>. 

To use this subroutine, pass an array with the following keys 
and values. Not every key is required.

=over 4

=item db

The name of the database or a reference to an established BioPerl
database object. A L<Bio::DB::SeqFeature::Store> or
L<Bio::DB::BigWigSet> database are typically used. Only required when
no features are provided.

=item feature

Pass either a single dataset name as a scalar or an anonymous array
reference of a list of dataset names. These may have been provided as
a command line option and need to be verified. If nothing is passed,
then a list of possible datasets will be presented to the user to be
chosen.

=item prompt

Provide a phrase to be prompted to the user to help in selecting
datasets from the list. If none is provided, a generic prompt will be
used.

=item single

A Boolean value (1 or 0) indicating whether only a single dataset is
allowed when selecting datasets from a presented list. If true, only
one dataset choice is accepted. If false, one or more dataset choices
are allowed. Default is false.

=item limit

Optionally provide a word or regular expression that matches the
feature type (C<primary_tag> only; C<source_tag>, if present, is ignored).
For example, provide "geneE<verbar>mrna" to only present gene and mRNA
features to the user. This is only applicable when a user must select
from a database list. The default is to list all available feature
types.

=back

=item check_dataset_for_rpm_support($dataset, [$cpu])

   # count the total number of alignments
   my $dataset = '/path/to/file.bam';
   my $total_count = check_dataset_for_rpm_support($dataset);
   
   # use multithreading
   my $total_count = check_dataset_for_rpm_support($dataset, $cpu);

This subroutine will check a dataset for RPM, or Reads Per Million mapped, 
support. Only two types of database files support this, Bam files and 
BigBed files. If the dataset is either one of these, or the name of a 
database feature which points to one of these files, then it will 
calculate the total number of mapped alignments (Bam file) or features 
(BigBed file). It will return this total number. If the dataset does 
not support RPM (because it is not a Bam or BigBed file, for example), 
then it will return undefined.

Pass this subroutine one or two values. The first is the name of the 
dataset. Ideally it should be validated using verify_or_request_feature_types() 
and have an appropriate prefix (file, http, or ftp). 

For multi-threaded execution pass a second value, the number of CPU cores 
available to count Bam files. This will speed up counting bam files 
considerably. The default is 2 for environments where L<Parallel::ForkManager> 
is installed, or 1 where it is not.

=item get_new_feature_list

	my $Data = Bio::ToolBox::Data->new;
	my $db_name = 'cerevisiae';
	my %data = get_new_feature_list(
		data      => $Data,
		db        => $db_name,
		features  => 'genes',
	);

This subroutine will generate a new feature list collected from an annotation 
database, such as a L<Bio::DB::SeqFeature::Store> database. It will populate 
an empty L<Bio:ToolBox::Data> object data table with a list of the collected 
features. Columns include the Primary ID, Name, and feature Type. 

The subroutine is passed an array containing three required arguments. The keys 
include

=over 4

=item data

An empty L<Bio::ToolBox::Data> object. This is required.

=item db

The name of the database or a reference to an established database object. 
Currently, this must be a L<Bio::DB::SeqFeature::Store> database.
This is required.

=item features

A scalar value containing a name representing the feature type(s) of
feature(s) to collect. The type may be a C<primary_tag> or C<primary_tag:source_tag>.
It may be a single element or a comma-delimited list. This is required.

=item chrskip

Provide a regular-expression compatible string for skipping features from 
specific chromosomes, for example mitochondrial or unmapped contigs.

=back

Status messages, including the number of features found and kept, are always 
printed to C<STDOUT>. 

=item get_new_genome_list

	my $Data = Bio::ToolBox::Data->new();
	my $db_name = 'cerevisiae';
	my $window_size = 500;
	my $step_size = 250;
	my %data = get_new_genome_list(
		data      => $Data,
		db        => $db_name,
		win       => $window_size,
		step      => $step_size,
	);

This subroutine will generate a new list of genomic windows or intervals. 
The genome is split into intervals of a specified size with the specified 
step size.

The subroutine is passed an array containing the required arguments. 

=over 4

=item data

An empty L<Bio::ToolBox::Data> object. This is required.

=item db

The name of the database or a reference to an established database
object. Any BioPerl adapter, including those described in the L</DESCRIPTION>, 
are supported. Required. 

=item win

A scalar value containing an integer representing the size of the
window in basepairs. The default value is 500 bp.

=item step

A scalar value containing an integer representing the step size for
advancing the window across the genome. The default is the window size.

=item chrskip

=item skip

=item exclude

Provide a regular expression compatible string for excluding specific 
chromosomes or classes of chromosomes, such as the mitochondrial or 
unmapped contigs. 

=back

Status messages are printed to C<STDOUT>. 

=item get_db_feature

    my $db = open_db_connection('annotation.db');
    my $seqfeature = get_db_feature(
        db      => $db,
        type    => 'gene',
        name    => 'ABC1',
    );

This subroutine will retrieve a specific feature from a L<Bio::DB::SeqFeature::Store> 
database for subsequent analysis, manipulation, and/or score retrieval using the 
L</get_segment_score> methods. It relies upon unique information to pull out a 
single, unique feature.

Several attributes may be used to pull out the feature, including the feature's 
unique database primary ID, name andE<sol>or aliases, and GFF type (C<primary_tag> 
andE<sol>or C<source>). The L</get_new_feature_list> subroutine will generate a list 
of features with their unique identifiers. 

The C<primary_id> attribute is preferentially used as it provides the best 
performance. However, it is not always portable between databases or even a 
database re-load. In that case, the C<display_name> and C<type> are used to identify 
potential features. Note that the C<display_name> may not always be unique in the 
database. In this case, the addition of aliases may help. If all else fails, a new 
feature list should be generated. 

To get a feature, pass an array of arguments.

=over 4

=item db

The name of the L<Bio::DB::SeqFeature::Store> database or a reference
to an established database object. 

=item id

Provide the C<primary_id> tag. In the L<Bio::DB::SeqFeature::Store>
database schema this is a (usually) non-portable, unique identifier
specific to a database. It provides the fastest lookup.

=item name

A scalar value representing the feature C<display_name>. Aliases may be
appended with semicolon delimiters. 

=item type

Provide the feature type, which is typically expressed as
C<primary_tag:source>. Alternatively, provide just the C<primary_tag> only.

=back

While it is possible to identify features with any two attributes 
(or possibly just name or ID), the best performance is obtained with 
all three together. The first SeqFeature object is returned if multiple 
are found.

=item get_segment_score

    my $score = get_segment_score(
    	$seqfeature->seq_id,    # chromosome
    	$seqfeature->start,     # start
    	$seqfeature->end,       # end
    	$seqfeature->strand,    # strand
    	'sense',                # strandedness
    	'mean',                 # method
    	0,                      # return type
    	$db,                    # database
    	'scores.bw'             # datasets
    );

Generic method for collecting score(s) from a database. This is the 
primary interface to all of the database handling packages, including 
Bam, bigWig, bigBed, and USeq files, as well as SeqFeature databases. 
By passing common parameters here, the appropriate database or file 
handling packages are loaded during run time and the appropriate 
methods are called. Essentially, this is where the magic happens.

Pass the method an array of nine (or more) parameters. The parameters 
are not keyed as a hash, and the order is explicit. However, see the 
L<Bio::ToolBox::db_helper::constants> module for named index positions. 
The order is as follows.

=over 4

=item * 0 Chromosome

String representing chromosome or sequence ID.

=item * 1 Start

Start coordinate (integer, 1-based)

=item * 2 Stop

Stop or end coordinate (integer, 1-based)

=item * 3 Strand

BioPerl style strand value, including -1, 0, and 1.

=item * 4 Strandedness

A string indicating strandedness: sense, antisense, all.

=item * 5 Method

A string representing the method of data collection. 

Current acceptable methods include mean, median, sum, min, max, 
stddev, count, ncount (named count), and pcount (precise count). 
See the individual db_helper sub classes for more details. Not all 
adapters support all methods.

=item * 6 Return type

An integer (0, 1, 2) representing the type of value to return.

A return type of 0 indicates a single score should be returned, 
calculated using the specified method. A return type of 1 is an 
array or array reference of all scores found. A return type of 2 
is a hash or hash reference of coordinate => score. 

=item * 7 Dataset database

This is either a L<Bio::DB::SeqFeature::Store> or a L<Bio::DB::BigWigSet> 
database containing multiple score dataset features. If your dataset 
does not contain multiple different feature types, pass C<undefined>.

Pass either a string containing the name of database (which will be opened) 
or an already opened database object. 

=item * 8 Dataset

This is the verified dataset from which the scores will be passed. The 
dataset should be verified using the L</verify_or_request_feature_types> 
subroutine. For indexed file datasets, e.g. F<.bam> or F<.bw> files, 
this will the verified local path or URL. For a database, this will be 
a database type, sucha as C<primary_tag> or C<primary_tag:source>. 
or the path to a data file, e.g. Bam, bigWig, bigBed, or USeq file. 

If multiple datasets are to be combined, simply append them to the 
array (index 9, 10, etc).

=back

The returned item is dependent on the value of the return type code.

=item calculate_score

    my $scores = get_segment_score($chr,$start,$stop,$strand,
    	'sense','ncount',1,undef,'file:/path/to/dataset.bam');
    my $score = calculate_score('ncount', $scores);

This subroutine will perform the math to calculate a single score 
from an array of values. Current acceptable methods include mean, 
median, sum, min, max, stddev, count, ncount, and pcount.

Pass the subroutine two items: the name of the mathematical method 
enumerated above, and an array reference of the values to work on.
A single score will be returned.

=item get_chromosome_list

	my $db = open_db_connection('cerevisiae');
	# get all chromosomes in the database
	my @chromosomes = get_chromosome_list($db);
	foreach (@chromosomes) {
		my $name = $_->[0];
		my $length = $_->[1];
		print "chromosome $name is $length bp\n";
	}

This subroutine will collect a list of chromosomes or reference sequences 
in a Bio::DB database and return the list along with their sizes in bp. 
Many BioPerl-based databases are supported, including 
L<Bio::DB::SeqFeature::Store>, L<Bio::DB::Fasta>, L<Bio::DB::Sam>, 
L<Bio::DB::HTS>, L<Bio::DB::BigWig>, L<Bio::DB::BigWigSet>, 
L<Bio::DB::BigBed>, or any others that support the C<seq_ids> method. 
L<Bio::DB::Big> adapters for big files are also supported, though see note 
below. See the L</open_db_connection> subroutine for more information.

Pass the subroutine either the name of the database, the path to the 
database file, or an opened database object.

Optionally pass a second value, a regular expression compatible string 
for skipping specific chromosomes or chromosome classes, such as mitochondrial 
or unmapped contigs. The default is to return all chromosomes. 

The subroutine will return an array, with each element representing each 
chromosome or reference sequence in the database. Each element is an anonymous 
array of two elements, the chromosome name and length in bp.

B<NOTE>: L<Bio::DB::Big> databases for bigWig and bigBed files do 
not return chromosomes in the original order as the file, and are returned 
in a sorted manner that may not be the original order.

=item low_level_bam_coverage

	my $coverage = low_level_bam_coverage($sam, $tid, $start, $stop);

This is a convenience method for running the low level bam coverage method. 
Since both L<Bio::DB::Sam> and L<Bio::DB::HTS> bam file adapters are 
supported, and each has slight variation in the API syntax, this method helps 
to abstract the actual method and use the appropriate syntax depending on 
which adapter is loaded. It is best if the C<$sam> object was opened using the 
L</open_db_connection> method, or that C<$BAM_ADAPTER> is set.

B<NOTE> that this is the B<low level> coverage method based on the index object, 
and not the similarly named high level API method. Read the adapter 
documentation for proper usage.

=item low_level_bam_fetch

	my $success = low_level_bam_fetch($sam, $tid, $start, $stop, $callback, $data);

This is a convenience method for running the low level bam fetch method. 
Since both L<Bio::DB::Sam> and L<Bio::DB::HTS> bam file adapters are 
supported, and each has slight variation in the API syntax, this method helps 
to abstract the actual method and use the appropriate syntax depending on 
which adapter is loaded. It is best if the C<$sam> object was opened using the 
L</open_db_connection> method, or that C<$BAM_ADAPTER> is set.

B<NOTE> that this is the B<low level> fetch method based on the index object, 
and not the similarly named high level API method. Read the adapter 
documentation for proper usage.

=item get_genomic_sequence

	my $sequence = get_genomic_sequence($db, $chrom, $start, $stop);

This is a convenience wrapper function for fetching genomic sequence from 
three different supported fasta database adapters, which, of course, all 
have different APIs, including the L<Bio::DB::Fasta>, L<Bio::DB::HTS::Faidx>, 
and L<Bio::DB::Sam::Fai> adapters, as well as L<Bio::DB::SeqFeature::Store> 
and possibly other BioPerl databases. 

Pass the opened database object, chromosome, start, and stop coordinates. This 
assumes BioPerl standard 1-base coordinates. Only the forward sequence is 
retrieved. The sequence is returned as a simple string.

=back

=head1 INTERNAL SUBROUTINES

These are not intended for normal general consumption but are documented here 
so that at least some of us know what is going on.

=over 4

=item _lookup_db_method

Internal method for determining which database adapter and method call 
should be used given a set of parameters. This result is cached for 
future data queries (likely hundreds or thousands of repeated calls). 
The method is stored in the global hash C<%DB_METHODS>;

Pass the parameter array reference from L</get_segment_score>. This 
should load the appropriate db_helper sub module during run time.

=item _load_helper_module

Subroutine to load and import a db_helper module during run time.
Sets the appropriate global variable if successful. 

=item _load_bam_helper_module

Subroutine to determine which Bam adapter db_helper module to load, 
either the older samtools-based adapter or the newer HTSlib-based adapter.
Uses the global and exportable variable C<$BAM_ADAPTER> to set a 
preference for which adapter to use. Use 'sam' or 'hts' or some 
string containing these names. Do B<NOT> change this variable after 
loading the helper module; it will not do what you think it will do.
If both are available and a preference is not set then the hts helper 
is the default.

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

