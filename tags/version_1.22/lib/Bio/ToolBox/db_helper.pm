package Bio::ToolBox::db_helper;

use strict;
require Exporter;
use Carp qw(carp cluck croak confess);
use Module::Load; # for dynamic loading during runtime
use Statistics::Lite qw(
	sum
	mean
	median
	min
	max
	range
	stddevp
);
use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure 
);
use Bio::ToolBox::db_helper::config;
use Bio::ToolBox::utility;
use constant LOG2 => log(2);

our $VERSION = 1.22;


# check values for dynamically loaded helper modules
# these are loaded only when needed during runtime to avoid wasting resources
our $BAM_OK      = 0;
our $BIGBED_OK   = 0;
our $BIGWIG_OK   = 0;
our $FASTA_OK    = 0;
our $SEQSTORE_OK = 0;
our $USEQ_OK     = 0;
our $WIGGLE_OK   = 0;

# define reusable variables
our $TAG_EXCEPTIONS; # for repeated use with validate_included_feature()
our %total_read_number; # for rpm calculations
our $primary_id_warning; # for out of date primary IDs
our %OPENED_DB; # cache for some opened Bio::DB databases

# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw();
our @EXPORT_OK = qw(
	open_db_connection
	get_dataset_list 
	verify_or_request_feature_types 
	check_dataset_for_rpm_support 
	get_new_feature_list 
	get_new_genome_list 
	validate_included_feature 
	get_feature 
	get_chromo_region_score 
	get_region_dataset_hash 
	get_chromosome_list 
);


# The true statement
1; 

=head1 NAME

Bio::ToolBox::db_helper

=head1 DESCRIPTION

These are helper subroutines to work with relevant databases that can be 
accessed through BioPerl modules. These include the Bio::DB::SeqFeature::Store 
relational database as well as Bio::DB::Fasta, Bio::DB::Bam, Bio::DB::BigWig, 
Bio::DB::BigWigSet, and Bio::DB::BigBed databases. The primary functions 
included opening connections to these databases, verifying features found 
within the database, generating lists of features or genomic intervals for 
subsequent analysis, and collecting features and/or scores within regions 
of interest. Collected scores may be summarized using a variety of statistical 
methods.

When collecting scores or features, the data may hosted in a variety of 
formats and locations. These include:

=over 4

=item SeqFeature::Store Database

Full features may be stored, including genes, transcripts, exons, etc. 
Simple datasets, such as from microarray, may also be stored as the 
score value in the source GFF file. 

References to local, binary, indexed files may also be included as 
attributes to features stored in the database. Supported files 
include binary wig files (.wib, see Bio::Graphics::Wiggle) using the 
'wigfile' attribute, or bigWig files using the 'wigfile' or 'bigwigfile' 
attribute. The attribute values must be full paths. 

SeqFeature::Store databases are usually hosted by a relational database 
server (MySQL or PostGreSQL), SQLite file, or an in-memory database 
(for small GFF3 files only).

=item BigWig files

BigWig files are compressed, binary, indexed versions of text wig files 
and may be accessed either locally or remotely. They support extremely  
fast score retrieval from any genomic location of any size without 
sacrificing resolution (spatial and numeric).

=item Directory of BigWig files

A directory containing one or more BigWig files is assembled into a 
BigWigSet, allowing for metadata, such as strand, to be associated with 
BigWig files. 

=item BigBed files

BigBed files are compressed, binary, indexed versions of text BED files 
and may be accessed either locally or remotely. They support extremely  
fast score and feature retrieval from any genomic location.

=item Bam files

Bam files are compressed, binary, indexed versions of the text SAM file, 
or sequence alignment map. They are used with next generation sequencing 
technologies. They support individual alignment retrieval as well as 
read depth coverage. 

=item USeq files

USeq files are compressed, binary, indexed files that support BED type 
annotations or wig type scores distributed across the genome. They 
support rapid, random access across the genome and are comparable to 
both BigWig and BigBed files.

=back

While the functions within this library may appear to be simply a rehashing of 
the methods and functions in Bio::DB::SeqFeature::Store and other modules, they 
either provide a simpler function to often used database methodologies or are 
designed to work intimately with the BioToolBox data format file and data structures 
(see C<Bio::ToolBox::file_helper>). One key advantage to these functions is the ability
to work with datasets that are stranded (transcriptome data, for example).

Historically, this module was initially written to use Bio::DB::GFF for 
database usage. It has since been re-written to use Bio::DB::SeqFeature::Store 
database schema. Additional functionality has also been added for other databases 
and file formats. 

Complete usage and examples for the functions are provided below.

=head1 USAGE

Call the module at the beginning of your perl script and include the module 
names to export. None are exported by default.

  use Bio::ToolBox::db_helper qw(
	  get_new_feature_list 
	  get_feature
  );

This will export the indicated subroutine names into the current namespace. 
Their usage is detailed below. 

=cut






################################################################################
################           General subroutines             #####################
################################################################################


### Open a connection to the SeqFeature Store MySQL database

=head1 EXPORTED SUBROUTINES

=over 4

=item open_db_connection($database, $no_cache)

This module will open a connection to a BioPerl style database.
It returns an object that represents the connection. Several 
different types of databases are supported.

=over 4

=item Bio::DB::SeqFeature::Store database

These may be represented by a relational database (e.g. MySQL database), 
a SQLite database file (file.sqlite or file.db), or a single GFF3 file 
(file.gff) that can be loaded into an in-memory database. In-memory databases 
should only be used with small files as they demand a lot of memory.

Parameters for connecting to a relational database are stored in the BioToolBox 
configuration file, C<.biotoolbox.cfg>. These include database adaptors, 
user name, password, etc. Information regarding the configuration file may 
be found within the file itself. 

=item Bio::DB::Sam database 

A self-contained database represented by a sorted, indexed Bam file 
(file.bam). See http://samtools.sourceforge.net for more details. Files 
may be either local or remote (prefixed with http:// or ftp://).

=item Bio::DB::BigWig database

A self-contained database of scores represented by a BigWig (file.bw). See
L<http://genome.ucsc.edu/goldenPath/help/bigWig.html> for more information.
Files may be either local or remote (prefixed with http:// or ftp://).

=item Bio::DB::BigWigSet database

A local or remote directory of one or more BigWig files that can treated 
collectively as a database. A special text file may be included in the 
directory to assign metadata to each BigWig file, including attributes such 
as type, source, display name, strand, etc. See L<Bio::DB::BigWigSet> for 
more information on the formatting of the metadata file.

=item Bio::DB::BigBed database

A self-contained database of regions represented by a BigBed (file.bb). See
L<http://genome.ucsc.edu/goldenPath/help/bigBed.html> for more information.
Files may be either local or remote (prefixed with http:// or ftp://).

=item Bio::DB::USeq database

A self-contained database file of indexed regions or scores represented by 
a useq archive file (file.useq). See 
L<http://useq.sourceforge.net/useqArchiveFormat.html> for more information. 
Files may only be local.

=item Bio::DB::Fasta Database

A database of fasta sequences. A single multi-fasta or a directory of 
fasta files may be specified. The directory or parent directory must be 
writeable by the user to write a small index file. 

=back

Pass the name of a relational database or the path of the database file to 
the subroutine. The opened database object is returned. If it fails, then 
an error message should be generated and nothing is returned.

B<Important!> If you are forking your perl process, B<always> re-open your 
database objects in the child process, and pass a second true value 
to avoid using the cached database object. By default, opened databases are 
cached to improve efficiency, but this will be disastrous when crossing forks. 

Example:

	my $db_name = 'cerevisiae';
	my $db = open_db_connection($db_name);
	
	my $file = '/path/to/file.bam';
	my $db = open_db_connection($file);
	
	# within a forked child process
	# reusing the same variable to simplify code
	# pass second true value to avoid cache
	$db = open_db_connection($file, 1); 


=cut

sub open_db_connection {
	my $database = shift;
	my $no_cache = shift || 0;
	unless ($database) {
		cluck 'no database name passed!';
		return;
	}
	
	# first check if it is a database reference
	my $db_ref = ref $database;
	if ($db_ref =~ /^Bio::DB/) {
		# the provided database is already an open database object
		# nothing to open, return as is
		
		# determine the name if possible
		my $db_name;
		if ($db_ref =~ /^Bio::DB::SeqFeature::Store/) {
			# a SeqFeature database, using any DBI adapter
			$db_name = $database->{'dbh'}->{'name'}; 
				# dig through the object internals to identify the original 
				# name of the database
				# this should be relatively well documented through DBI
				# but could break in the future since it's not official API
		}
		elsif ($db_ref eq 'Bio::DB::Sam') {
			# a Bam database
			$db_name = $database->{'bam_path'};
		}
			# determining the database name from other sources is
			# either not possible or not easy, so won't bother unless
			# there is a really really good need
		
		# return as appropriate either both object and name or just object
		return wantarray ? ($database, $db_name) : $database;
	}
	
	# check to see if we have already opened it
	if (exists $OPENED_DB{$database} and not $no_cache) {
		# return the cached database connection
		# but only if user explicitly requested no cached databases
		# DO NOT reuse databases if you have forked!!!
		return wantarray ? ($OPENED_DB{$database}, $database) : $OPENED_DB{$database};
	}
	
	
	### Attempt to open the database
	# we go through a series of checks to determine if it is remote, local, 
	# an indexed big data file, SQLite file, etc
	# when all else fails, try to open a SQL connection
	my $db;
	my $error;
	
	# check if it is a remote file
	if ($database =~ /^http|ftp/i) {
		
		# a remote Bam database
		if ($database =~ /\.bam$/i) {
			# open using BigWig adaptor
			$BAM_OK = _load_helper_module('Bio::ToolBox::db_helper::bam') unless $BAM_OK;
			if ($BAM_OK) {
				$db = open_bam_db($database);
				unless ($db) {
					$error = " ERROR: could not open remote Bam file" .
						" '$database'! $!\n";
				}
			}
			else {
				$error = " Bam database cannot be loaded because\n" . 
					" Bio::DB::Sam is not installed\n";
			}
		}
		
		# a remote BigBed database
		elsif ($database =~ /\.bb$/i) {
			# open using BigBed adaptor
			$BIGBED_OK = _load_helper_module('Bio::ToolBox::db_helper::bigbed') 
				unless $BIGBED_OK;
			if ($BIGBED_OK) {
				$db = open_bigbed_db($database);
				unless ($db) {
					$error = " ERROR: could not open remote BigBed file" .
						" '$database'! $!\n";
				}
			}
			else {
				$error = " BigBed database cannot be loaded because\n" . 
					" Bio::DB::BigBed is not installed\n";
			}
		}
		
		# a remote BigWig database
		elsif ($database =~ /\.bw$/i) {
			# open using BigWig adaptor
			$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::bigwig') 
				unless $BIGWIG_OK;
			if ($BIGWIG_OK) {
				$db = open_bigwig_db($database);
				unless ($db) {
					$error = " ERROR: could not open remote BigWig file" .
						" '$database'! $!\n";
				}
			}
			else {
				$error = " BigWig database cannot be loaded because\n" . 
					" Bio::DB::BigWig is not installed\n";
			}
		}
		
		# a remote useq file
		elsif ($database =~ /\.useq$/) {
			# uh oh! remote useq files are not supported
			$error = " ERROR: remote useq files are not supported!\n";
		}
		
		# a presumed remote directory, presumably of bigwig files
		else {
			# open using BigWigSet adaptor
			$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::bigwig') 
				unless $BIGWIG_OK;
			if ($BIGWIG_OK) {
				$db = open_bigwigset_db($database);
				unless ($db) {
					$error = " ERROR: could not open presumed remote " .
						"BigWigSet directory '$database'! $!\n";
				}
			}
			else {
				$error = " Presumed BigWigSet database cannot be loaded because\n" . 
					" Bio::DB::BigWigSet is not installed\n";
			}
		}
	
	}
	
	# a directory, presumably of bigwig files
	elsif (-d $database) {
		# try opening using the BigWigSet adaptor
		$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::bigwig') 
			unless $BIGWIG_OK;
		if ($BIGWIG_OK) {
			$db = open_bigwigset_db($database);
			unless ($db) {
				$error = " ERROR: could not open local BigWigSet " . 
					"directory '$database'!\n";
				$error .= "   Does directory contain bigWig .bw files?\n";
			}
		}
		else {
			$error = " Presumed BigWigSet database cannot be loaded because\n" . 
				" Bio::DB::BigWigSet is not installed\n";
		}
		
		# try opening with the Fasta adaptor
		unless ($db) {
			$FASTA_OK = _load_helper_module('Bio::DB::Fasta') unless $FASTA_OK;
			if ($FASTA_OK) {
				undef $@; # in case it was set by something else
				eval {
					# to prevent annoying error messages
					local $SIG{__WARN__} = sub {}; 
					$db = Bio::DB::Fasta->new($database);
				};
				unless ($db) {
					$error .= " ERROR: could not open fasta directory '$database'!\n";
					$error .= "   Does directory contain fasta files? If it contains a" . 
						" directory.index file,\n   try deleting it and try again.\n";
				}
			}
			else {
				$error .= " Module Bio::DB::Fasta is required to open presumed fasta " . 
					"directory\n";
			}
		}
	}
	
	# check for a known file type
	elsif ($database =~ /gff|bw|bb|bam|useq|db|sqlite|fa|fasta/i) {
		
		# first check that it exists
		if (-e $database and -r _) {
		
			# a single gff3 file that we can load into memory
			if ($database =~ /\.gff3?(?:\.gz)?$/i) {
				$SEQSTORE_OK = _load_helper_module('Bio::DB::SeqFeature::Store') 
					unless $SEQSTORE_OK;
				if ($SEQSTORE_OK) {
					print " Loading file into memory database...\n";
					undef $@;
					eval {
						$db = Bio::DB::SeqFeature::Store->new(
							-adaptor => 'memory',
							-gff     => $database,
						);
					};
					unless ($db) {
						$error = " ERROR: could not load file '$database' into memory! $@\n";
					}
				}
				else {
					$error .= " Module Bio::DB::SeqFeature::Store is required to load " . 
						"GFF3 files into memory\n";
				}
			}
			
			# a SQLite database
			elsif ($database =~ /\.(?:sqlite|db)$/i) {
				# open using SQLite adaptor
				$SEQSTORE_OK = _load_helper_module('Bio::DB::SeqFeature::Store') 
					unless $SEQSTORE_OK;
				if ($SEQSTORE_OK) {
					undef $@;
					eval {
						$db = Bio::DB::SeqFeature::Store->new(
							-adaptor  => 'DBI::SQLite',
							-dsn      => $database,
						);
					};
					unless ($db) {
						$error = " ERROR: could not open SQLite file '$database'! $@\n";
					}
				}
				else {
					$error .= " Module Bio::DB::SeqFeature::Store is required to access" . 
						" database\n";
				}
			}
			
			# a Bam database
			elsif ($database =~ /\.bam$/i) {
				# open using BigWig adaptor
				$BAM_OK = _load_helper_module('Bio::ToolBox::db_helper::bam') 
					unless $BAM_OK;
				if ($BAM_OK) {
					undef $@;
					$db = open_bam_db($database);
					unless ($db) {
						$error = " ERROR: could not open local Bam file" .
							" '$database'! $@\n";
					}
				}
				else {
					$error = " Bam database cannot be loaded because\n" . 
						" Bio::DB::Sam is not installed\n";
				}
			}
			
			# a BigBed database
			elsif ($database =~ /\.bb$/i) {
				# open using BigBed adaptor
				$BIGBED_OK = _load_helper_module('Bio::ToolBox::db_helper::bigbed') 
					unless $BIGBED_OK;
				if ($BIGBED_OK) {
					undef $@;
					$db = open_bigbed_db($database);
					unless ($db) {
						$error = " ERROR: could not open local BigBed file" .
							" '$database'! $@\n";
					}
				}
				else {
					$error = " BigBed database cannot be loaded because\n" . 
						" Bio::DB::BigBed is not installed\n";
				}
			}
			
			# a BigWig database
			elsif ($database =~ /\.bw$/i) {
				# open using BigWig adaptor
				$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::bigwig') 
					unless $BIGWIG_OK;
				if ($BIGWIG_OK) {
					undef $@;
					$db = open_bigwig_db($database);
					unless ($db) {
						$error = " ERROR: could not open local BigWig file" .
							" '$database'! $@\n";
					}
				}
				else {
					$error = " BigWig database cannot be loaded because\n" . 
						" Bio::DB::BigWig is not installed\n";
				}
			}
			
			# a useq database
			elsif ($database =~ /\.useq$/i) {
				# open using USeq adaptor
				$USEQ_OK = _load_helper_module('Bio::ToolBox::db_helper::useq') 
					unless $USEQ_OK;
				if ($USEQ_OK) {
					$db = open_useq_db($database);
					unless ($db) {
						$error = " ERROR: could not open local useq file" .
							" '$database'! $!\n";
					}
				}
				else {
					$error = " Useq database cannot be loaded because\n" . 
						" Bio::DB::USeq is not installed\n";
				}
			}
			
			# a Fasta File
			elsif ($database =~ /\.fa(?:sta)?$/i) {
				# open using the Fasta adaptor
				$FASTA_OK = _load_helper_module('Bio::DB::Fasta') unless $FASTA_OK;
				if ($FASTA_OK) {
					undef $@;
					eval {
						# to prevent annoying error messages
						local $SIG{__WARN__} = sub {}; 
						$db = Bio::DB::Fasta->new($database);
					};
					unless ($db) {
						$error = " ERROR: could not open fasta file '$database'!\n";
						if (-e "$database\.index") {
							$error .= "   Try deleting $database\.index and try again\n";
						}
					}
				}
				else {
					$error .= " Fasta file could not be loaded because Bio::DB::Fasta" . 
						" is not installed\n";
				}
			}
		}
		
		# file does not exist or can be read
		else {
			if (not -e _) {
				# file does not exist
				$error = " ERROR: file '$database' does not exist!\n";
			}
			else {
				# file must not be readable then
				$error = " ERROR: file '$database' can not be read!\n";
			}
		}
	}
	
	# unrecognized real file
	elsif (-e $database) {
		# file exists, I just don't recognize the extension
		$error = " File '$database' type is not recognized\n";
	}
	
	# otherwise assume the name of a SeqFeature::Store database in the configuration
	
	# attempt to open using information from the configuration
	# using default connection information as necessary
	unless ($db) {
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
		$SEQSTORE_OK = _load_helper_module('Bio::DB::SeqFeature::Store') 
			unless $SEQSTORE_OK;
		if ($SEQSTORE_OK) {
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
		
			unless ($db) {
				$error .= " ERROR: unknown $adaptor database '$database'\n";
				if ($dsn =~ /mysql|pg/i) {
					$error .= "   using user '$user' and password\n";
				}
			}
		}
		else {
			$error .= " Module Bio::DB::SeqFeature::Store is required to connect " . 
				"to $adaptor databases\n";
		}
	}
	
	# conditional return
	if ($db) {
		# cache the opened connection for use later
		$OPENED_DB{$database} = $db unless $no_cache;
		
		# return as appropriate either both object and name or just object
		return wantarray ? ($db, $database) : $db;
	} 
	else {
		$error .= " no database could be found or connected!\n";
		warn $error;
		return;
	}
}



### Retrieve a list of the microrarray data sets from the db

=item get_dataset_list

This subroutine will retrieve a list of the available features stored in the 
database and returns an array of the features' types.

For Bio::DB::SeqFeature::Store databases, the type is represented as 
"type:source", corresponding to the third and second GFF columns, 
respectively. The types are sorted alphabetically first by source, 
then by method.

For Bio::DB::BigWigSet databases, the type, primary_tag, method, or 
display_name attribute may be used, in that respective order of 
availability. The list is sorted alphabetically.

Pass either the name of the database or an established database object. 
Supported databases include both Bio::DB::SeqFeature::Store and 
Bio::DB::BigWigSet databases. 

Example:

	my $db_name = 'cerevisiae';
	my @types = get_dataset_list($db_name);

=cut

sub get_dataset_list {
	
	my $database = shift;
	my $use_all_features = shift;
	
	# Open a db connection 
	my ($db, $db_name) = open_db_connection($database);
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	my $db_ref = ref $db;
	
	# process the database types, according to the type of database
	my @types;
	
	# a SeqFeature database
	if ($db_ref =~ m/Bio::DB::SeqFeature::Store/) {
		
		# collect the database types, 
		# sort first by source, then by method
		@types = (
			map $_->[2],
			sort { ($a->[0] cmp $b->[0]) or ($a->[1] cmp $b->[1]) } 
			map [$_->source, $_->method, $_],
			$db->types
		);
	}
	
	# a BigWigSet database
	elsif ($db_ref eq 'Bio::DB::BigWigSet') {
		
		# get the metadata
		my $metadata = $db->metadata;
		
		# collect
		foreach my $file (keys %{ $metadata }) {
			# get the type for each file
			my ($primary, $type, $name);
			
			# get the appropriate tags
			foreach my $attribute (keys %{ $metadata->{$file} } ) {
				if ($attribute =~ m/^primary_tag|method$/i) {
					$primary = $metadata->{$file}{$attribute};
				}
				elsif ($attribute =~ m/^type/i) {
					$type = $metadata->{$file}{$attribute};
				}
				elsif ($attribute =~ m/^display_name/i) {
					$name = $metadata->{$file}{$attribute};
				}
			}
				
			# store the type
			push @types, $type || $primary || $name;
		}
		
		# sort the types in alphabetical order
		@types = sort {$a cmp $b} @types;
	}
	
	# some other database
	else {
		carp " no dataset lists for database type $db_ref!\n";
	}
	
	# finish
	return @types;
}



### Process and verify a dataset

=item verify_or_request_feature_types()

This subroutine will process a list of feature types or data sources to be 
used for data or feature collection. There are two modes of action. If a 
list was provided, it will process and verify the list. If no list was 
provided, then the user will be prompted to select from a list of available 
feature types in the provided database.

The provided list may be a list of feature types or a list of single 
file data sources, including Bam, bigWig, or bigBed files. If the list 
includes feature types, they are verified as existing in the provided 
database. If the list includes file names, the local files are checked 
and the names prefixed with "file:" for use downstream.

If no list was provided, then a list of available feature types in the 
provided database will be presented to the user, and the user prompted 
to make a selection. One or more types may be selected, and a single 
item may be enforced if requested. The response is filtered through 
the parse_list() method from C<Bio::ToolBox::db_helper>, so a mix of single 
numbers or a range of numbers may be accepted. The responses are then 
validated.

Two types of databases are supported: L<Bio::DB::SeqFeature::Store> and 
L<Bio::DB::BigWigSet> databases. For SeqFeature Store databases, the 
type is comprised of "method:source", corresponding to the third and 
second GFF columns. For BigWigSet databases, types correspond to either 
the type, method, primary_tag, or display_name attributes, in that 
order of availability. 

For feature types or files that pass validation, they are returned as a 
list. Types of files that do not pass validation are printed to STDOUT. 

To use this subroutine, pass an array with the following keys 
and values. Not every key is required.

  db       => The name of the database or a reference to an 
              established BioPerl database object. A 
              Bio::DB::SeqFeature::Store or Bio::DB::BigWigSet 
              database are typically used. Only required when 
              no features are provided.
  feature  => Pass either a single dataset name as a scalar or an 
              anonymous array reference of a list of dataset names. 
              These may have been provided as a command line option and 
              need to be verified. If nothing is passed, then a list of 
              possible datasets will be presented to the user to be 
              chosen.
  prompt   => Provide a phrase to be prompted to the user to help in 
              selecting datasets from the list. If none is provided, a 
              generic prompt will be used.
  single   => A Boolean value (1 or 0) indicating whether only a single 
              dataset is allowed when selecting datasets from a 
              presented list. If true, only one dataset choice is 
              accepted. If false, one or more dataset choices are 
              allowed. Default is false.

The subroutine will return a list of the accepted datasets. It will print 
bad dataset names to standard out.

Example:

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
	

=cut

sub verify_or_request_feature_types {
	
	# Retrieve passed values
	my %args = @_; # the passed argument values as a hash reference
	
	# Check for single option
	$args{'single'} ||= 0;
	
	# Collect the datasets
	my @datasets;
	$args{'feature'} ||= undef;
	if (ref $args{'feature'} eq 'ARRAY') {
		# an anonymous array of datasets
		@datasets = @{ $args{'feature'} };
	}
	elsif (defined $args{'feature'}) {
		push @datasets, $args{'feature'};
	}
	# else it is null and @datasets remains empty, prompting user input
	
	
	# Open database object and collect features
	$args{'db'} ||= undef;
	my $db = $args{'db'} ? open_db_connection( $args{'db'} ) : undef;
	
	# Initialize main output arrays
	my @good_datasets;
	my @bad_datasets;
	my %db_features; # hash for features found in the database, use when needed
	
	# Check provided datasets
	if (@datasets) {
		
		# check for multiple comma-delimited datasets
		my @list_to_check;
		foreach my $item (@datasets) {
			if ($item =~ /,/) {
				# this one has a comma, therefore it has more than dataset
				push @list_to_check, split(/,/, $item);
			}
			else {
				# a singleton
				push @list_to_check, $item;
			}
		}
		
		# now verify the datasets
		foreach my $dataset (@list_to_check) {
			
			# check for a remote file
			if ($dataset =~ /^(?: http | ftp) .+ \. (?: bam | bw | bb) $/xi) {
				# a remote file
				# assume it is good, no verification here though
				# it will either work or won't work
				push @good_datasets, $dataset;
			}
			
			# a local file
			elsif ($dataset =~ /\.(?:bw|bb|bam|useq)$/i) {
				# presume we have a local indexed data file
				
				# user may have requested two or more files to be merged
				# these should be combined with an ampersand
				# check each one 
				my @files;
				foreach my $file (split /\&/, $dataset) {
					if (-e $file) {
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
					push @good_datasets, join("&", @files);
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
							my ($method, $source) = split /:/, $type;
							if ($source) {
								$db_features{$type} = 1;
								$db_features{$method} = 1;
							}
							else {
								$db_features{$type} = 1;
							}
						}
						
						# verify
						unless (%db_features) {
							carp " provided database has no feature types " . 
								"to verify dataset(s) against!\n";
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
				foreach my $d (split /&/, $dataset) {
					# validate this dataset
					if (exists $db_features{$d}) {
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
				$db_features{$i} = $type;
				$i++;
			}
			
			# verify
			unless (%db_features) {
				carp " provided database has no feature types " . 
					"to collect!\n";
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
		foreach (sort {$a <=> $b} keys %db_features) {
			# print out the list of microarray data sets
			print "  $_\t$db_features{$_}\n"; 
		}
		
		# prompt the user
		$args{'prompt'} ||= undef;
		if ($args{'prompt'}) {
			# provided custom prompt
			print $args{'prompt'};
		}
		else {
			# generic prompt
			if ($args{'single'}) {
				print " Enter one number for the data set or feature   ";
			}
			else {
				print " Enter the number(s) or range of the data set(s) or" . 
					" feature(s)   ";
			}
		}
		
		# get answer from the user
		my $answer = <STDIN>;
		chomp $answer;
		my @answer_list = parse_list($answer);
		
		# take the first one if requested
		if ($args{'single'}) {
			unless (scalar @answer_list == 1) {
				splice(@answer_list, 1);
			}
		}
		
		# verify the answer list
		foreach my $answer (@answer_list) {
			
			# check for merged datasets
			if ($answer =~ /&/) {
				# a merged dataset
				my @list = split /&/, $answer;
				my $check = 1;
				
				# check all are good
				foreach (@list) {
					unless (exists $db_features{$_}) {
						$check = 0;
					}
				}
				
				# if all are good
				if ($check) {
					push @good_datasets, 
						join( "&", map { $db_features{$_} } @list);
				}
				else {
					push @bad_datasets, $answer;
				}
			}
			
			else {
				# a single dataset
				# check if it is good
				
				if (exists $db_features{$answer}) {
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
	if ($args{'single'}) {
		return $good_datasets[0];
	}
	else {
		return @good_datasets;
	}
}





### Process and verify a dataset

=item check_dataset_for_rpm_support()

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
and have an appropriate prefix (file, http, or ftp). If it does not 
have a prefix, then it is assumed to be a database feature. The second 
passed feature is the name of a BioPerl database or an opened database 
object. This database will be checked for the indicated dataset, and 
the first returned feature checked for an attribute pointing to a 
supported file.

For multi-threaded execution pass a third value, the number of CPU cores 
available to count Bam files. The default is 2 for environments where 
Parallel::ForkManager is installed, or 1 where it is not.

=cut

sub check_dataset_for_rpm_support {
	
	# get passed dataset and databases
	my $dataset  = shift;
	my $database = shift;
	my $cpu      = shift;
	
	# check that we haven't done this already
	
	
	# Check the dataset to see if supports RPM or RPKM method
	# if so, then calculate the total number of reads
	# this uses the global variable $rpkm_read_sum
	my $rpm_read_sum;
	
	if (exists $total_read_number{$dataset}) {
		# this dataset has already been summed
		# no need to do it again
		
		$rpm_read_sum = $total_read_number{$dataset};
	}
	
	elsif ($dataset =~ /\.bam$/) {
		# a bam file dataset
		
		if ($BAM_OK) {
			# Bio::ToolBox::db_helper::bam was loaded ok
			# sum the number of reads in the dataset
			$rpm_read_sum = sum_total_bam_alignments($dataset, 0, 0, $cpu);
		}
		else {
			carp " Bam support is not available! " . 
				"Is Bio::DB::Sam installed?\n";
			return;
		}
	}
	
	elsif ($dataset =~ /\.bb$/) {
		# a bigbed file dataset
		
		if ($BIGBED_OK) {
			# Bio::ToolBox::db_helper::bigbed was loaded ok
			# sum the number of features in the dataset
			$rpm_read_sum = sum_total_bigbed_features($dataset);
		}
		else {
			carp " BigBed support is not available! " . 
				"Is Bio::DB::BigBed installed?\n";
			return;
		}
	}
	
	elsif ($dataset !~ /^file|ftp|http/i) {
		# a database feature
		# this feature might point to a bam or bigbed file
		
		# Check the database
		my $db; # the database object to be used
		if (defined $database) {
			my $db_ref = ref $database;
			if ($db_ref =~ /^Bio::DB/) {
				$db = $database;
			}
			else {
				# the name of a database was passed, create a database connection
				$db = open_db_connection($database);
				unless ($db) {
					carp " No database to check feature '$dataset' for RPM support!\n";
					return;
				}
			}
		}
		else {
			carp " No database to check feature '$dataset' for RPM support!\n";
			return;
		}
		
		# get a sample of the features from the database
		my @features;
		if ($dataset =~ /&/) {
			# in case we have a combined datasets
			@features = $db->features(-type => [ split /&/, $dataset ]);
		}
		else {
			@features = $db->features(-type => $dataset);
		}
		unless (@features) {
			 # no feature found
			 # basically nothing to do
			 return;
		}
		
		# look for the database file in the attributes
		if ($features[0]->has_tag('bamfile')) {
			# specifying a bam file
			my ($bamfile) = $features[0]->get_tag_values('bamfile');
			
			if ($BAM_OK) {
				# Bio::ToolBox::db_helper::bam was loaded ok
				# sum the number of reads in the dataset
				$rpm_read_sum = sum_total_bam_alignments($bamfile, 0, 0, $cpu);
			}
			else {
				carp " Bam support is not available! " . 
					"Is Bio::DB::Sam installed?\n";
				return;
			}
		}
		
		elsif ($features[0]->has_tag('bigbedfile')) {
			# specifying a bigbed file
			my ($bedfile) = $features[0]->get_tag_values('bigbedfile');
			
			if ($BIGBED_OK) {
				# Bio::ToolBox::db_helper::bigbed was loaded ok
				# sum the number of features in the dataset
				$rpm_read_sum = sum_total_bigbed_features($bedfile);
			}
			else {
				carp " BigBed support is not available! " . 
					"Is Bio::DB::BigBed installed?\n";
				return;
			}
		}
		
		else {
			# can't find an appropriate dataset
			# nothing to do
			return;
		}
	}
	
	else {
		# some other non-supported dataset
		return;
	}
	
	# return the sum value if we've made it this far
	$total_read_number{$dataset} = $rpm_read_sum;
	return $rpm_read_sum;
}







################################################################################
################           Feature subroutines             #####################
################################################################################


### Generate a new list of features


=item get_new_feature_list 

This subroutine will generate a new feature list collected from the database. 
Once the list of genomic features is generated, then data may be collected
for each item in the list. 

The subroutine will generate and return a data hash as described in 
C<Bio::ToolBox::file_helper>. The data table will have two or three columns. The 
feature name and type:source are listed in columns one and two, respectively.
If the features have an Alias tag, then a third column is included with 
a comma delimited list of the feature aliases.

The subroutine is passed an array containing the arguments. 
The keys include

  Required:
  db       => The name of the database or a reference to an 
              established database object. 
  features => A scalar value containing a name representing the 
              type(s) of feature(s) to collect. This name will be 
              parsed into an actual list with the internal subroutine 
              _features_to_classes(). Refer to that documentation for 
              a list of appropriate features.

The subroutine will return a reference to the data hash. It will print 
status messages to STDOUT. 

Example

	my $db_name = 'cerevisiae';
	my %data = get_new_feature_list(
		'db'        => $db_name,
		'features'  => 'genes',
	);


=cut


sub get_new_feature_list {

	# Retrieve passed values
	my %args = @_; # the passed argument values
	
	# Open a db connection 
	$args{'db'} ||= undef;
	my ($db, $db_name) = open_db_connection($args{'db'});
	unless ($db) {
		carp 'no database connected!';
		return;
	}

	# Verify a SeqFeature::Store database
	my $db_ref = ref $db;
	unless ($db_ref =~ /^Bio::DB::SeqFeature::Store/) {
		carp "Database type $db_ref doesn't support generating feature lists!\n";
		return;
	}
	
	
	# Translate the features into a list of classes
	my @classes = _features_to_classes($args{'features'});
	unless (@classes) {
		carp "no or unknown features passed!";
		return;
	}
	
	# Generate data structures
	my $new_data = generate_tim_data_structure(
		$args{'features'},
		'Primary_ID',
		'Name',
		'Type'
	);
	unless ($new_data) {
		cluck " cannot generate tim data structure!\n";
		return;
	}
	
	# name of the database
	$new_data->{'db'} = $db_name; 
	
	# List of types
	if (scalar @classes > 1) {
		$new_data->{1}->{'include'} = join(",", @classes);
	}
	
	# Get the names of chromosomes to avoid
	my @excluded_chromosomes = 
		$BTB_CONFIG->param("$db_name\.chromosome_exclude");
	unless (@excluded_chromosomes) {
		@excluded_chromosomes = 
			$BTB_CONFIG->param('default_db.chromosome_exclude');
	}
	my %excluded_chr_lookup = map {$_ => 1} @excluded_chromosomes;
	
	# Set up the database iterator
	print "   Searching for " . join(", ", @classes) . "\n";
	my $iterator = $db->get_seq_stream(
			-types    => \@classes
	); 
	unless ($iterator) {
		# there should be some features found in the database
		carp "could not get feature iterator for database";
		return;
	}
	
	
	# Walk through the collected features
	my $total_count = 0; # total found features
	while (my $feature = $iterator->next_seq) {
		$total_count++;
		
		# skip genes from excluded chromosomes
		next if exists $excluded_chr_lookup{ $feature->seq_id };
		
		# skip anything that matches the tag exceptions
		next unless validate_included_feature($feature);
		
		# Record the feature information
		# in the B::DB::SF::S database, the primary_ID is a number unique to the 
		# the specific database, and is not portable between databases
		push @{ $new_data->{'data_table'} }, [
			$feature->primary_id, 
			$feature->display_name,  
			$feature->type,
		];
		$new_data->{'last_row'} += 1;
	}
	
	# print result of search
	printf "   Found %s features in the database\n", format_with_commas($total_count);
	printf "   Kept %s features.\n", format_with_commas($new_data->{'last_row'});
	
	# return the new data structure
	return $new_data;
}



### Generate a new list genomic windows

=item get_new_genome_list 

This subroutine will generate a new list of genomic windows. The genome
is split into intervals of a specified size that is moved along the 
genome in specified step sizes.

The subroutine will generate and return a data hash as described in 
C<Bio::ToolBox::file_helper>. The data table will have 3 columns, including 
Chromosome, Start, and Stop.

The subroutine is passed an array containing the arguments. 
The keys include

  Required:
  db       => The name of the database or a reference to an 
              established database object. 
  Optional: 
  win      => A scalar value containing an integer representing the
              size of the window in basepairs. The default value 
              is defined in biotoolbox.cfg file.
  step     => A scalar value containing an integer representing the
              step size for advancing the window across the genome. 
              The default is the window size.

The subroutine will return a reference to the data hash. It will print 
status messages to STDOUT. 

Example

	my $db_name = 'cerevisiae';
	my $window_size = 500;
	my $step_size = 250;
	my %data = get_new_genome_list(
		'db'        => $db_name,
		'win'       => $window_size,
		'step'      => $step_size,
	);


=cut

sub get_new_genome_list {

	# Collect the passed arguments
	my %args = @_; 
	
	
	# Open a db connection 
	$args{'db'} ||= undef;
	my ($db, $db_name) = open_db_connection($args{'db'});
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	
	# Determine win and step sizes
	$args{'win'} ||= undef;
	unless ($args{'win'}) {
		$args{'win'} = 
			$BTB_CONFIG->param("$db_name\.window") ||
			$BTB_CONFIG->param('default_db.window');
		print "  Using default window size $args{win} bp\n";
	}
	$args{'step'} ||= $args{'win'};
		
	
	# Generate data structures
	my $new_data = generate_tim_data_structure(
		'genome',
		'Chromosome',
		'Start',
		'Stop'
	);
	unless ($new_data) {
		cluck " cannot generate tim data structure!\n";
		return;
	}
	my $feature_table = $new_data->{'data_table'}; 
	
	# Begin loading basic metadata information
	$new_data->{'db'}      = $db_name; # the db name
	$new_data->{1}{'win'}  = $args{'win'}; # under the Start metadata
	$new_data->{1}{'step'} = $args{'step'};
	
	
	
	# Collect the chromosomes
	# include option to exclude those listed in biotoolbox.cfg that
	# we don't want
	my @chromosomes = get_chromosome_list($db, 1);
	unless (@chromosomes) {
		carp " no sequence IDs were found in the database!\n";
		return;
	}
	
	
	# Collect the genomic windows
	print "   Generating $args{win} bp windows in $args{step} bp increments\n";
	foreach (@chromosomes) {
		
		# name and length as sub-array in each element
		my ($chr, $length) = @{$_};
		
		for (my $start = 1; $start <= $length; $start += $args{step}) {
			# set the end point
			my $end = $start + $args{win} - 1; 
			
			if ($end > $length) {
				# fix end to the length of the chromosome
				$end = $length;
			} 
			
			# add to the output list
			push @{$feature_table}, [ $chr, $start, $end,];
			$new_data->{'last_row'}++;
		}
	}
	print "   Kept " . $new_data->{'last_row'} . " windows.\n";
	
	# Return the data structure
	return $new_data;
}


=item validate_included_feature

This subroutine will validate a database feature to make sure it is 
useable. It will check feature attributes and compare them against 
a list of attributes and values to be avoided. The list of unwanted 
attributes and values is stored in the BioToolBox configuration file 
biotoolbox.cfg. 

Pass the subroutine a Bio::DB::SeqFeature::Store database feature. 
It will return true (1) if the feature passes validation and false 
(undefined) if it contains an excluded attribute and value.

=cut

sub validate_included_feature {
	
	# feature to check
	my $feature = shift;
	
	# get the list of feature exclusion tags
	unless (defined $TAG_EXCEPTIONS) {
		$TAG_EXCEPTIONS = $BTB_CONFIG->get_block('exclude_tags');
	}
	
	# Check the tag exceptions
	# we will check for the presence of the exception tag
	# if the feature tag value matches the exception tag value
	# then the feature should be excluded
	foreach my $key (keys %{ $TAG_EXCEPTIONS }) {
		if ($feature->has_tag($key)) {
			# first check that the feature has this tag
			my @tag_values = $feature->get_tag_values($key);
			if (ref $TAG_EXCEPTIONS->{$key} eq 'ARRAY') {
				# there's more than one exception value!
				# need to check them all
				foreach my $exception (@{ $TAG_EXCEPTIONS->{$key} }) {
					if (grep {$_ eq $exception} @tag_values) {
						# this feature should be excluded
						return;
					}
				}
			}
			else {
				# single tag exception value
				if (grep {$_ eq $TAG_EXCEPTIONS->{$key}} @tag_values) {
					# this feature should be excluded
					return;
				}
			}
		}
	}
	
	# if we've made it thus far, the feature is good
	return 1;
}




=item get_feature

This subroutine will retrieve a specific feature from a Bio::DB::SeqFeature::Store 
database for subsequent analysis, manipulation, and/or score retrieval using the 
get_chromo_region_score() or get_region_dataset_hash() methods. It relies upon 
unique information to pull out a single, unique feature.

Several attributes may be used to pull out the feature, including the feature's 
unique database primary ID, name and/or aliases, and GFF type (primary_tag and/or 
source). The get_new_feature_list() subroutine will generate a list of features 
with their unique identifiers. 

The primary_id attribute is preferentially used as it provides the best 
performance. However, it is not portable between databases or even re-loading. 
In that case, the display_name and type are used to identify potential features. 
Note that the display_name may not be unique in the database. In this case, the 
addition of aliases may help. If all else fails, a new feature list should be 
generated. 

To get a feature, pass an array of arguments.
The keys include

  Required:
  db       => The name of the Bio::DB::SeqFeature::Store database or 
              a reference to an established database object. 
  id       => Provide the primary_id tag. In the 
              Bio::DB::SeqFeature::Store database schema this is a 
              (usually) non-portable, unique identifier specific to a 
              database. It provides the fastest lookup.
  name     => A scalar value representing the feature display_name. 
              Aliases may be appended with semicolon delimiters. 
  type     => Provide the feature type, which is typically expressed 
              as primary_tag:source. Alternatively, provide just the 
              primary_tag only.

While it is possible to identify features with any two attributes 
(or possibly just name or ID), the best performance is obtained with 
all three together.

The first SeqFeature object is returned if found.

=cut

sub get_feature {
	my %args = @_;
	
	# Open a db connection 
	$args{'db'} ||= undef;
	my $db = open_db_connection($args{'db'});
	unless ($db) {
		croak 'no database connection for getting a feature!';
	}
	my $db_ref = ref $db;
	
	# get the name of the feature
	my $name = $args{'name'} || undef; 
	$name = (split(/\s*[;,\|]\s*/, $name))[0] if $name =~ /[;,\|]/;
		 # multiple names present using common delimiters ;,|
		 # take the first name only, assume others are aliases that we don't need
	
	# check for values and internal nulls
	$args{'id'} = exists $args{'id'} ? $args{'id'} : undef;
	$args{'type'} ||= undef;
	undef $name         if $name eq '.';
	undef $args{'id'}   if $args{'id'} eq '.';
	undef $args{'type'} if $args{'type'} eq '.';
	
	# quick method for feature retrieval
	if (defined $args{'id'} and $db_ref =~ /SeqFeature::Store/) {
		# we have a unique primary_id to identify the feature
		# fetch method only works with Bio::DB::SeqFeature::Store databases
		my $feature = $db->fetch($args{'id'}) || undef; 
		
		# check that the feature we pulled out is what we want
		my $check = $feature ? 1 : 0;
		if ($check) {
			$check = 0 if (defined $name and $feature->display_name ne $name);
			$check = 0 if (defined $args{'type'} and $feature->type ne $args{'type'});
		}
		
		# return if things match up as best we can tell
		if ($check) {
			return $feature;
		}
		else {
			# the primary_ids are out of date
			unless ($primary_id_warning) {
				warn "CAUTION: Some primary IDs in Input file appear to be out of date\n";
				$primary_id_warning++;
			}
		}
	}
	
	# otherwise use name and type 
	my @features = $db->features( 
		-name       => $name, # we assume this name will be unique
		-aliases    => 0, 
		-type       => $args{'type'},
	);
	
	# if none are found
	unless (@features) {
		# try again with aliases enabled
		@features = $db->features( 
			-name       => $name,
			-aliases    => 1, 
			-type       => $args{'type'},
		);
	}
	
	# check the number of features returned
	if (scalar @features > 1) {
		# there should only be one feature found
		# if more are returned, there's redundant or duplicated data in the db
		
		# first check whether addition of aliases may improve accuracy
		if ($args{'name'} =~ /;/) {
			# we have presumed aliases in the name value
			my $check = $args{'name'};
			$check =~ s/\s*;\s*/;/g; # strip spaces just in case
			
			# check each feature and try to match name and aliases
			my @candidates;
			foreach my $f (@features) {
				my $f_name = join(';', $f->display_name, ($f->get_tag_values('Alias')));
				if ($check eq $f_name) {
					# found one
					push @candidates, $f;
				}
			}
			
			if (scalar @candidates == 1) {
				# we found it!
				return shift @candidates;
			}
			elsif (scalar @candidates > 1) {
				# hopefully we've improved a little bit?
				@features = @candidates;
			}
		}
		
		# warn the user, this should be fixed
		printf "  Found %s %s features named '$name' in the database! Using first one.\n", 
			scalar(@features), $args{'type'};
	}
	elsif (!@features) {
		printf "  Found no %s features named '$name' in the database!\n", $args{'type'};
		return;
	}
	
	# done 
	return shift @features; 
}





################################################################################
##################           Score subroutines             #####################
################################################################################



### Get a dataset score for a single region

=item get_chromo_region_score 

This subroutine will retrieve a dataset value for a single specified 
region in the genome. The region is specified with chromosomal coordinates:
chromosome name, start, and stop. It will collect all dataset values within the
window, combine them with the specified method, and return the single value.

The subroutine is passed an array containing the arguments. 
The keys include

  Required:
  db       => The name of the database or a reference to an 
              established database object. Optional if an 
              indexed dataset file (Bam, BigWig, etc.) is provided.
  dataset  => The name of the dataset in the database to be 
              collected. The name should correspond to a feature 
              type in the database, either as type or type:source. 
              The name should be verified using the 
              subroutine verify_or_request_feature_types() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces. Alternatively, specify a data file name. 
              A local file should be prefixed with 'file:', while 
              a remote file should be prefixed with the transfer 
              protocol (ftp: or http:).
  method   => The method used to combine the dataset values found
              in the defined region. Acceptable values include 
              sum, mean, median, range, stddev, min, max, rpm, 
              rpkm, and scores. See _get_segment_score() 
              documentation for more info.
  chromo   => The name of the chromosome (reference sequence)
  start    => The start position of the region on the chromosome
  stop     => The stop position of the region on the chromosome
  end      => Alias for stop
  
  Optional:
  strand   => The strand of the region (-1, 0, or 1) on the 
              chromosome. The default is 0, or unstranded.
  value    => Specify the type of value to collect. Acceptable 
              values include score, count, or length. The default 
              value type is score. 
  log      => Boolean value (1 or 0) indicating whether the dataset 
              values are in log2 space or not. If undefined, the 
              dataset name will be checked for the presence of the 
              phrase 'log2' and set accordingly. This argument is
              critical for accurately mathematically combining 
              dataset values in the region.
  stranded => Indicate whether the dataset values from a specific 
              strand relative to the feature should be collected. 
              Acceptable values include sense, antisense, or all.
              Default is 'all'.
  rpm_sum  => When collecting rpm or rpkm values, the total number  
              of alignments may be provided here. Especially 
              useful when collecting via parallel forked processes, 
              otherwise each fork will sum the number of alignments 
              independently, an expensive proposition. 
         	  
The subroutine will return the region score if successful.

Examples

	my $db = open_db_connection('cerevisiae');
	my $score = get_chromo_region_score(
		'db'      => $db,
		'method'  => 'mean',
		'dataset' => $dataset,
		'chr'     => $chromo,
		'start'   => $startposition,
		'stop'    => $stopposition,
		'log'     => 1,
	);
	


=cut

sub get_chromo_region_score {
	
	# retrieve passed values
	my %args = @_; 
	
	# check the data source
	unless ($args{'dataset'}) {
		confess " no dataset requested!";
	}
	
	# Open a db connection 
	$args{'db'} ||= $args{'ddb'} || undef;
	unless ($args{'db'} or $args{'dataset'} =~ /^(?:file|http|ftp)/) {
		# database is only really necessary if we are not using an indexed file dataset
		confess "no database provided!";
	}
	my $db;
	if ($args{'db'}) {
		$db = open_db_connection( $args{'db'} ) or 
			confess "cannot open database!";
	}
	
	# establish coordinates
	$args{'chromo'} ||= $args{'seq'} || $args{'seq_id'} || undef;
	$args{'start'}    = exists $args{'start'} ? $args{'start'} : 1;
	$args{'start'}    = 1 if ($args{'start'} <= 0);
	$args{'stop'}   ||= $args{'end'} || 1;
	$args{'strand'}   = exists $args{'strand'} ? $args{'strand'} : 0;
	
	unless ($args{chromo} and $args{start} and $args{stop}) {
		my $s = sprintf "%s:%s..%s", $args{chromo}, $args{start}, $args{stop};
		cluck "one or more provided genomic coordinates ($s) are invalid!\n";
		return;
	};
	if ($args{'stop'} < $args{'start'}) {
		# coordinates are flipped, reverse strand
		if ($args{'stop'} <= 0) {
			cluck "invalid stop coordinate $args{stop} provided\n";
			return;
		}
		my $stop = $args{'start'};
		$args{'start'} = $args{'stop'};
		$args{'stop'}  = $stop;
		$args{'strand'} = -1;
	}
	
	# define default values as necessary
	$args{'value'}    ||= 'score';
	$args{'stranded'} ||= 'all';
	$args{'log'}      ||= undef;
	unless (defined $args{'log'}) {
		# we need to know whether we are working with a log2 dataset or not
		# as it will adversely affect the math!
		if ($args{'dataset'} =~ /log2/i) {
			# if we're smart we'll encode the log2 status in the dataset name
			# but chances are, we're not that smart
			$args{'log'} = 1;
		} else {
			# otherwise assume it is non-log
			# unsafe, but what else to do? we'll put the onus on the user
			$args{'log'} = 0;
		}
	}
	
	# set RPM sum value if necessary
	if (exists $args{'rpm_sum'} and defined $args{'rpm_sum'}) {
		unless (exists $total_read_number{ $args{'dataset'} }) {
			$total_read_number{ $args{'dataset'} } = $args{'rpm_sum'};
		}
	}
	
	# get the scores for the region
	# pass to internal subroutine to combine dataset values
	return _get_segment_score(
				$db,
				$args{'chromo'},
				$args{'start'},
				$args{'stop'},
				$args{'strand'}, 
				$args{'dataset'},
				$args{'value'},
				$args{'method'}, 
				$args{'stranded'},  
				$args{'log'},
	);
}





### Retrieve hash of dataset values for a region


=item get_region_dataset_hash 

This subroutine will retrieve dataset values or feature attributes from
features located within a defined region and return them as a hash.
The (start) positions will be the hash keys, and the corresponding dataset 
values or attributes will be the hash values. The region is defined based on 
a genomic feature in the database. The region corresponding to the entire 
feature is selected by default. Alternatively, it may be adjusted by passing 
appropriate arguments.

Different dataset values may be collected. The default is to collect 
score values of the dataset features found within the region (e.g. 
microarray values). Alternatively, a count of found dataset features
may be returned, or the lengths of the found dataset features. When
lengths are used, the midpoint position of the feature is used in the
returned hash rather than the start position.

The returned hash is keyed by relative coordinates and their scores. For 
example, requesting a region from -200 to +200 of a feature (using the 
start and stop options, below) will return a hash whose keys are relative 
to the feature start position, i.e. the keys will >= -200 and <= 200. 
Absolute coordinates relative to the reference sequence or chromosome 
may be optionally returned instead.

The subroutine is passed an array containing the arguments. 
The keys include

  Required:
  db       => The name of the annotation database or a reference to 
              an established database object. 
  dataset  => The name of the dataset in the database to be 
              collected. The name should correspond to a feature 
              type in the database, either as type or type:source. 
              The name should be verified using the 
              subroutine verify_or_request_feature_types() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces. Alternatively, specify a data file name. 
              A local file should be prefixed with 'file:', while 
              a remote file should be prefixed with the transfer 
              protocol (ftp: or http:).
  Required for database features:
  id       => The Primary ID of the genomic feature. This is 
              database specific and may be used alone to identify 
              a genomic feature, or in conjunction with name and type.
  name     => The name of the genomic feature. Required if the Primary 
              ID is not provided. 
  type     => The type of the genomic feature. Required if the Primary 
              ID is not provided. 
  Required for coordinate positions:
  chromo   => The chromosome or sequence name (seq_id). This may be 
              used instead of name and type to specify a genomic 
              segment. This must be used with start and stop options, 
              and optionally strand options.
  Optional:
  start    => Indicate an integer value representing the start  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with "stop".
  stop|end => Indicate an integer value representing the stop  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with "start".
  ddb      => The name of the data-specific database or a reference 
              to an established database. Use when the data and 
              annotation are in separate databases.
  extend   => Indicate an integer value representing the number of 
              bp the feature's region should be extended on both
              sides.
  position => Indicate the relative position of the feature from 
              which the "start" and "stop" positions are calculated.
              Three values are accepted: "5", which denotes the 
              5' end of the feature, "3" which denotes the 
              3' end of the feature, or "4" which denotes the 
              middle of the feature. This option is only used in 
              conjunction with "start" and "stop" options. The 
              default value is "5".
  strand   => For those features or regions that are NOT 
              inherently stranded (strand 0), artificially set the 
              strand. Three values are accepted: -1, 0, 1. This 
              will overwrite any pre-existing value (it will not, 
              however, migrate back to the database).
  stranded => Indicate whether the dataset values from a specific 
              strand relative to the feature should be collected. 
              Acceptable values include sense, antisense, or all.
              Default is all.
  value    => Indicate which attribute will be returned. Acceptable 
              values include "score", "count", or "length". The  
              default behavior will be to return the score values.
  avoid    => Boolean value to indicate that other features of the 
              same type should be avoided. This only works if name 
              and type was provided. Any positioned scores which 
              overlap the other feature(s) are not returned. The 
              default is false (return all values).
  absolute => Boolean value to indicate that absolute coordinates 
              should be returned, instead of transforming to 
              relative coordinates, which is the default.
          	  
The subroutine will return the hash if successful.

Example

	my $db = open_db_connection('cerevisiae');
	my %region_scores = get_region_dataset_hash(
		'db'      => $db,
		'dataset' => $dataset,
		'name'    => $name,
		'type'    => $type,
	);
	


=cut

sub get_region_dataset_hash {
	
	# retrieve passed values
	my %args = @_; 
	
	### Initialize parameters
	
	# check the data source
	$args{'dataset'} ||= undef;
	unless ($args{'dataset'}) {
		confess " no dataset requested!";
	}
	
	# Open a db connection 
	$args{'db'} ||= undef;
	my $db;
	if ($args{'db'}) {
		$db = open_db_connection( $args{'db'} ) or 
			confess "cannot open database!";
	}
	
	# Open the data database if provided
	$args{'ddb'} ||= undef;
	my $ddb;
	if ($args{'ddb'}) {
		$ddb = open_db_connection( $args{'ddb'} ) or
			confess "requested data database could not be opened!\n";
	}
	else {
		# reuse something else
		if ($db) {
			$ddb = $db;
		}
		elsif ($args{'dataset'} =~ /^(?:file|http|ftp)/) {
			$ddb = $args{'dataset'};
		}
		else {
			confess "no database or indexed dataset supplied!";
		}
	}
	
	# confirm options and check we have what we need 
	$args{'name'}   ||= undef;
	$args{'type'}   ||= undef;
	$args{'id'}     ||= undef;
	$args{'chromo'} ||= $args{'seq'} || $args{'seq_id'} || undef;
	$args{'start'}    = exists $args{'start'} ? $args{'start'} : 1;
	$args{'stop'}   ||= $args{'end'} || 1;
	$args{'strand'}   = exists $args{'strand'} ? $args{'strand'} : undef;
	unless (
		(defined $args{'name'} and defined $args{'type'}) or 
		(defined $args{'chromo'} and $args{'start'} and $args{'stop'})
	) {
		cluck "the feature name and type or genomic coordinates are missing!";
		return;
	};
	if ($args{'stop'} < $args{'start'}) {
		# coordinates are flipped, reverse strand
		my $stop = $args{'start'};
		$args{'start'} = $args{'stop'};
		$args{'stop'}  = $stop;
		$args{'strand'} = -1;
	}
	
	# assign other defaults
	$args{'stranded'} ||= 'all';
	$args{'value'}    ||= 'score';
	$args{'position'} ||= 5;
	$args{'extend'}   ||= 0;
	$args{'avoid'}    ||= 0;
	$args{'absolute'} ||= 0;
	
	
	# the final coordinates
	my $fref_pos; # to remember the feature reference position
	my $fchromo;
	my $fstart;
	my $fstop;
	my $fstrand;
	
	
	
	### Define the chromosomal region segment
	# we will use the primary database to establish the intitial feature
	# and determine the chromosome, start and stop
	
	
	# Extend a named database feature
	if (
		( $args{'id'} or ( $args{'name'} and $args{'type'} ) ) and 
		$args{'extend'}
	) {
		
		# first define the feature to get endpoints
		confess "database required to use named features" unless $db;
		my $feature = get_feature(
			'db'    => $db,
			'id'    => $args{'id'},
			'name'  => $args{'name'},
			'type'  => $args{'type'},
		) or return; 
		
		# determine the strand
		$fstrand   = defined $args{'strand'} ? $args{'strand'} : $feature->strand;
		
		# record the feature reference position and strand
		if ($args{'position'} == 5 and $fstrand >= 0) {
			$fref_pos = $feature->start;
		}
		elsif ($args{'position'} == 3 and $fstrand >= 0) {
			$fref_pos = $feature->end;
		}
		elsif ($args{'position'} == 5 and $fstrand < 0) {
			$fref_pos = $feature->end;
		}
		elsif ($args{'position'} == 3 and $fstrand < 0) {
			$fref_pos = $feature->start;
		}
		elsif ($args{'position'} == 4) {
			# strand doesn't matter here
			$fref_pos = $feature->start + int(($feature->length / 2) + 0.5);
		}
		
		# record final coordinates
		$fchromo = $feature->seq_id;
		$fstart  = $feature->start - $args{'extend'};
		$fstop   = $feature->end + $args{'extend'};
	} 
		
	# Specific start and stop coordinates of a named database feature
	elsif (
		( $args{'id'} or ( $args{'name'} and $args{'type'} ) ) and 
		$args{'start'} and $args{'stop'}
	) {
		# first define the feature to get endpoints
		confess "database required to use named features" unless $db;
		my $feature = get_feature(
			'db'    => $db,
			'id'    => $args{'id'},
			'name'  => $args{'name'},
			'type'  => $args{'type'},
		) or return; 
		
		# determine the strand
		$fstrand   = defined $args{'strand'} ? $args{'strand'} : $feature->strand;
		
		# determine the cooridnates based on the identified feature
		if ($args{'position'} == 5 and $fstrand >= 0) {
			# feature is on forward, top, watson strand
			# set segment relative to the 5' end
			
			# record final coordinates
			$fref_pos  = $feature->start;
			$fchromo   = $feature->seq_id;
			$fstart    = $feature->start + $args{'start'};
			$fstop     = $feature->start + $args{'stop'};
		}
		
		elsif ($args{'position'} == 5 and $fstrand < 0) {
			# feature is on reverse, bottom, crick strand
			# set segment relative to the 5' end
			
			# record final coordinates
			$fref_pos  = $feature->end;
			$fchromo   = $feature->seq_id;
			$fstart    = $feature->end - $args{'stop'};
			$fstop     = $feature->end - $args{'start'};
		}
		
		elsif ($args{'position'} == 3 and $fstrand >= 0) {
			# feature is on forward, top, watson strand
			# set segment relative to the 3' end
			
			# record final coordinates
			$fref_pos = $feature->end;
			$fchromo   = $feature->seq_id;
			$fstart    = $feature->end + $args{'start'};
			$fstop     = $feature->end + $args{'stop'};
		}
		
		elsif ($args{'position'} == 3 and $fstrand < 0) {
			# feature is on reverse, bottom, crick strand
			# set segment relative to the 3' end
			
			# record final coordinates
			$fref_pos = $feature->start;
			$fchromo   = $feature->seq_id;
			$fstart    = $feature->start - $args{'stop'};
			$fstop     = $feature->start - $args{'start'};
		}
		
		elsif ($args{'position'} == 4) {
			# feature can be on any strand
			# set segment relative to the feature middle
			
			# record final coordinates
			$fref_pos = $feature->start + int(($feature->length / 2) + 0.5);
			$fchromo   = $feature->seq_id;
			$fstart    = $fref_pos + $args{'start'};
			$fstop     = $fref_pos + $args{'stop'};
		}
	}
	
	# an entire named database feature
	elsif ( $args{'id'} or ( $args{'name'} and $args{'type'} ) ) {
		
		# first define the feature to get endpoints
		confess "database required to use named features" unless $db;
		my $feature = get_feature(
			'db'    => $db,
			'id'    => $args{'id'},
			'name'  => $args{'name'},
			'type'  => $args{'type'},
		) or return; 
		
		# determine the strand
		$fstrand   = defined $args{'strand'} ? $args{'strand'} : $feature->strand;
		
		# record the feature reference position and strand
		if ($args{'position'} == 5 and $fstrand >= 0) {
			$fref_pos = $feature->start;
		}
		elsif ($args{'position'} == 3 and $fstrand >= 0) {
			$fref_pos = $feature->end;
		}
		elsif ($args{'position'} == 5 and $fstrand < 0) {
			$fref_pos = $feature->end;
		}
		elsif ($args{'position'} == 3 and $fstrand < 0) {
			$fref_pos = $feature->start;
		}
		elsif ($args{'position'} == 4) {
			# strand doesn't matter here
			$fref_pos = $feature->start + int(($feature->length / 2) + 0.5);
		}
		
		# record final coordinates
		$fchromo   = $feature->seq_id;
		$fstart    = $feature->start;
		$fstop     = $feature->end;
	}
	
	# a genomic region
	elsif ( $args{'chromo'} and defined $args{'start'} and defined $args{'stop'} ) {
		# coordinates are easy
		
		$fchromo   = $args{'chromo'};
		if ($args{'extend'}) {
			# user wants to extend
			$fstart    = $args{'start'} - $args{'extend'};
			$fstop     = $args{'stop'}  + $args{'extend'};
		}
		else {
			$fstart    = $args{'start'};
			$fstop     = $args{'stop'};
		}
		$fstart = 1 if $fstart <= 0;
		
		# determine the strand
		$fstrand   = defined $args{'strand'} ? $args{'strand'} : 0; # default is no strand
		
		# record the feature reference position and strand
		if ($args{'position'} == 5 and $fstrand >= 0) {
			$fref_pos = $args{'start'};
		}
		elsif ($args{'position'} == 3 and $fstrand >= 0) {
			$fref_pos = $args{'stop'};
		}
		elsif ($args{'position'} == 5 and $fstrand < 0) {
			$fref_pos = $args{'stop'};
		}
		elsif ($args{'position'} == 3 and $fstrand < 0) {
			$fref_pos = $args{'start'};
		}
		elsif ($args{'position'} == 4) {
			# strand doesn't matter here
			$fref_pos = $args{'start'} + 
				int( ( ($args{'stop'} - $args{'start'} + 1) / 2) + 0.5);
		}
	}
	
	# or else something is wrong
	else {
		confess " programming error! not enough information provided to" .
			" identify database feature!\n";
	}
	
	# sanity check for $fstart
	$fstart = 1 if $fstart < 1;
	
	### Data collection
	my %datahash = _get_segment_score(
		$ddb, # using the data database here
		$fchromo,
		$fstart,
		$fstop,
		$fstrand, 
		$args{'dataset'}, 
		$args{'value'},
		'indexed', # method
		$args{'stranded'}, 
		0, # log value
	);
	
	
	### Check for conflicting features
	if ($args{'avoid'} and $args{'type'}) {
		# we need to look for any potential overlapping features of the 
		# same type and remove those scores
		
		# get the overlapping features of the same type
		my @overlap_features = $db->features(
			-seq_id  => $fchromo,
			-start   => $fstart,
			-end     => $fstop,
			-type    => $args{'type'}
		);
		if (@overlap_features) {
			# there are one or more feature of the same type in this 
			# region
			# one of them is likely the one we're working with
			# but not necessarily - user may be looking outside original feature
			# the others are not what we want and therefore need to be 
			# avoided
			foreach my $feat (@overlap_features) {
				
				# skip the one we want
				next if ($feat->display_name eq $args{'name'});
				
				# now eliminate those scores which overlap this feature
				foreach my $position (keys %datahash) {
					
					# delete the scored position if it overlaps with 
					# the offending feature
					if (
						$position >= $feat->start and
						$position <= $feat->end
					) {
						delete $datahash{$position};
					}
				}
			}
		}
		
	}
	
	
	
	### Convert the coordinates to relative positions
		# previous versions of this function that used Bio::DB::GFF returned 
		# the coordinates as relative positions, e.g. -200..200
		# to maintain this compatibility we will convert the coordinates to 
		# relative positions
		# most downstream applications of this function expect this, and it's
		# a little easier to work with. Just a little bit, though....
	if ($args{'absolute'}) {
		# do not convert to relative positions
		return %datahash;
	}
	else {
		my %relative_datahash;
		if ($fstrand >= 0) {
			# forward strand
			foreach my $position (keys %datahash) {
				# relative position is real position - reference
				$relative_datahash{ $position - $fref_pos } = $datahash{$position};
			}
		}
		elsif ($fstrand < 0) {
			# reverse strand
			foreach my $position (keys %datahash) {
				# the relative position is -(real position - reference)
				$relative_datahash{ $fref_pos - $position } = $datahash{$position};
			}
		}
		
		# return the collected dataset hash
		return %relative_datahash;
	}
}




=item get_chromosome_list

This subroutine will collect a list of chromosomes or reference sequences 
in a Bio::DB database and return the list along with their sizes in bp. 
Many BioPerl-based databases are supported, including 
Bio::DB::SeqFeature::Store, Bio::DB::Fasta, Bio::DB::Sam, Bio::DB::BigWig, 
Bio::DB::BigWigSet, and Bio::DB::BigBed, or any others that support the 
"seq_ids" method. See the L<open_db_connection> subroutine for more 
information.

Pass the subroutine either the name of the database, the path to the 
database file, or an opened database object.

Optionally pass a second value, a boolean argument to limit and exclude 
unwanted chromosomes as defined by the "chromosome_exclude" option in 
the BioToolBox configuration file, C<biotoolbox.cfg>. A true value limits 
chromosomes, and false includes all chromosomes. The default is to return 
all chromosomes. Sometimes some sequences are simply not wanted in 
analysis, like the mitochondrial chromosome or unmapped contigs.

The subroutine will return an array, with each element representing each 
chromosome or reference sequence in the database. Each element is an anonymous 
array of two elements, the chromosome name and length in bp.

Example
	
	my $db = open_db_connection('cerevisiae');
	# get all chromosomes in the database
	my @chromosomes = get_chromosome_list($db);
	foreach (@chromosomes) {
		my $name = $_->[0];
		my $length = $_->[1];
		print "chromosome $name is $length bp\n";
	}

=cut

sub get_chromosome_list {
	
	# options
	my $database = shift;
	my $limit = shift || 0;
	
	# Open a db connection 
	my ($db, $db_name) = open_db_connection($database);
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	# Check for BigWigSet database
	# these need to be handled a little differently
	if (ref $db eq 'Bio::DB::BigWigSet') {
		# BigWigSet databases are the only databases that don't 
		# support the seq_ids method
		# instead we have to look at one of the bigwigs in the set
		my $bw_file = ($db->bigwigs)[0];
		$db = open_db_connection($bw_file);
	}
	
	# Get chromosome exclusion list
	# these are chromosomes that we do not want to include in the final 
	# list
	# they must be explicitly requested to ignore
	my %excluded_chr_lookup;
	if ($limit) {
		my @excluded_chromosomes = 
			$BTB_CONFIG->param("$db_name\.chromosome_exclude");
		unless (@excluded_chromosomes) {
			@excluded_chromosomes = 
				$BTB_CONFIG->param('default_db.chromosome_exclude');
		}
		%excluded_chr_lookup = map {$_ => 1} @excluded_chromosomes;
	}
		
	
	# Collect chromosome lengths
	# we want to maintain the original order of the chromosomes so we 
	# won't be putting it into a hash
	# instead an array of arrays
	my @chrom_lengths;
	
	# Database specific approaches to collecting chromosomes
	# I used to have one generic approach, but idiosyncrasies and potential 
	# bugs make me use different approaches for better consistency
	
	# Bigfile
	if (ref $db eq 'Bio::DB::BigWig' or ref $db eq 'Bio::DB::BigBed') {
		foreach my $chr ($db->seq_ids) {
			
			# check for excluded chromosomes
			if (exists $excluded_chr_lookup{$chr} ) {
				next;
			}
			
			# get chromosome size
			my $length = $db->length($chr);
			
			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}
	
	# Bam
	elsif (ref $db eq 'Bio::DB::Sam') {
		for my $tid (0 .. $db->n_targets - 1) {
			# each chromosome is internally represented in the bam file as 
			# a numeric target identifier
			# we can easily convert this to an actual sequence name
			# we will force the conversion to go one chromosome at a time
			
			# sequence info
			my $chr    = $db->target_name($tid);
			my $length = $db->target_len($tid);
			
			# check for excluded chromosomes
			if (exists $excluded_chr_lookup{$chr} ) {
				next;
			}
			
			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}
	
	# Fasta
	elsif (ref $db eq 'Bio::DB::Fasta') {
		for my $chr ($db->get_all_ids) {
			
			# check for excluded chromosomes
			if (exists $excluded_chr_lookup{$chr} ) {
				next;
			}
			
			# get chromosome size
			my $seq = $db->get_Seq_by_id($chr);
			my $length = $seq ? $seq->length : 0;
			
			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}
	
	# SeqFeature::Store or other Bioperl
	else {
		foreach my $chr ($db->seq_ids) {
			
			# check for excluded chromosomes
			if (exists $excluded_chr_lookup{$chr} ) {
				next;
			}
			
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
	
	
	# Return
	unless (@chrom_lengths) {
		carp " no chromosome sequences identified in database!\n";
		return;
	}
	return @chrom_lengths;
}




################################################################################
###############           Internal subroutines             #####################
################################################################################


### Internal subroutine to convert a feature category into a list of classes

=back

=head1 INTERNAL SUBROUTINES

These are not intended for normal consumption but are documented here 
so that we know what is going on.

=over 4

=item _features_to_classes

This internal subroutine provides a conveniant look up and conversion of a 
single-word description of a category of features into a list of actual
feature classes in the database. For example, the word 'gene' may include
all ORFs, snRNAs, snoRNAs, and ncRNAs.

Pass the subroutine the feature category name as a scalar value. The 
actual list of feature types will be collected and returned as an array. 
Multiple values may be passed as a comma-delimited string (no spaces).

The aliases and feature lists are specified in the Bio::ToolBox::db_helper 
configuration file, biotoolbox.cfg. Additional lists and aliases 
may be placed there. The lists are database specific, or they can be 
added to the default database.

If the passed category name does not match an alias in the config file, 
then it is assumed to be a feature in the database. No further validation 
will be done (if it's not valid, simply no features would be returned from 
the database). 

Also, feature types may be passed as the GFF's method:source, in which case 
they are assumed to be valid and not checked.

=cut

sub _features_to_classes {
	my $feature = shift;
	my @types;
		
	my $alias2types = $BTB_CONFIG->get_block('features');
	if (exists $alias2types->{$feature} ) {
		# looks like the feature is an alias for a list of features
		# defined in the config file
		if (ref $alias2types->{$feature} eq 'ARRAY') {
			# there's a list of them
			@types = @{ $alias2types->{$feature} };
		}
		else {
			# only one
			$types[0] = $alias2types->{$feature};
		}
	}
	else { 
		# We'll assume this is a specific feature in the database.
		# It may be provided as type:source or type.
		# We'll pass these on directly to the originating subroutine
		# more than one may be passed delimited by commas, but no spaces
		@types = split /,/, $feature;
	} 
	
	return @types;
}




### Internal subroutine to retrieve the scores from an established region object

=item _get_segment_score

This internal subroutine is used to collect the dataset scores for an 
established genomic region. It works with a variety of data sources. 

First, the data may be stored directly in the Bio::DB::SeqFeature::Store 
(using the original GFF score value). Second, the feature may reference 
a data file as an attribute (e.g., wigfile=/path/to/file.wib). Finally, 
the name(s) of a data file may be passed from which to collect the data. 
Supported data files include BigWig (.bw), BigBed (.bb), USeq (.useq), 
and Bam (.bam). A Bio::DB::BigWigSet database is also supported.

The subroutine is passed an array of ten specific values, all of which 
must be defined and presented in this order. These values include
  
  [0] The opened database object that may be used to generate the 
      the segment and contain the data to collect. If the dataset 
      request is from a big file (bigWig, bigBed, Bam), then this 
      database will likely not be used. Otherwise a 
      Bio::DB::SeqFeature::Store or Bio::DB::BigWigSet database 
      object should be passed.
  [1] The chromosome or seq_id of the segment to examine
  [2] The start position of the segment
  [3] The stop or end position of the segment
  [4] The strand of the region (or original feature) -1, 0, or 1.
  [5] The dataset name for filename. Multiple datasets may be included, 
      delimited with an ampersand (&). Multiple datasets are merged into 
      one, unless excluded by strand. Local data source files should be 
      prepended with 'file:', while remote data source files should be 
      prepended with the transfer protocol (http: or ftp:).
  [6] The data type to be collecting. In most cases, the score value 
      is used, but others may be collected. Accepted values include
         
         score
         count
         length
         
  [7] The method of combining all of the dataset values found in the 
      segment into a single value. Accepted values include
         
         sum
         mean
         median
         min
         max
         range (returns difference between max and min)
         stddev (returns the standard deviation of a population)
         indexed (returns hash of postion => score)
         rpm (returns reads per million mapped, only valid with 
              bam and bigbed databases)
         rpkm (same as rpm but normalized for length in kb)
         scores (returns an array reference of all the raw scores)
         
  [8] The strandedness of acceptable data. Genomic segments 
      established from an inherently stranded database feature 
      (e.g. ORF) have a specific strandedness. If the dataset strand 
      does not match the specified criteria for strandedness, then it 
      is ignored. If the dataset does not have a specified strand, 
      then it is used regardless of the specified criteria for 
      strandedness. Currently, only transcription data is stranded. 
      Accepted values include
         
         sense       take only values on the same strand
         antisense   take only values on the opposite strand
         all         take all values
         
  [9] The log status of the dataset. Many microarray datasets are in 
      log2 format, and must be de-logged to perform accurate 
      mathematical calculations, such as mean or median. Supply a 
      boolean (0 or 1) value.
      
  
The subroutine returns a single numeric value if appropriate dataset values
are identified in the genomic segment. If not, then a non-value (.) is 
returned. Use of a period as a non-value avoids un-defined errors in 
some subsequent programmatic manipulations and provides a convenient human 
visual marker.

=cut

sub _get_segment_score {
	
	# get passed arguments
	my (
		$db,
		$chromo,
		$start,
		$stop,
		$strand, 
		$dataset, 
		$value_type,
		$method, 
		$strandedness, 
		$log
	) = @_;
	
	# define
	my @scores; # array of collected scores
	my %pos2data; # hash of position to scores
	my $dataset_type; # remember what type of database the data is from
	my $iterator; # seqfeature stream object for reiterating db features
	my $db_type = $db ? ref $db : undef; # source of the originating db 
	
	my @datasetlist = split /[&,]/, $dataset; 
		# multiple datasets may be combined into a single search, for example
		# transcriptome data that is on f and r strands. These are given as
		# ampersand or comma delimited lists
	
	
	### Determine where we are going to get the data
		# first check whether the provided dataset(s) look like a data file
		# next check whether the database segment object came from a BigWigSet
		# finally assume it is a SeqFeature database object
		# then look for a wigfile, bigwigfile, or bamfile attribute
		# finally then just take the score directly from the database objects
	
	### Data source files provided
	if ($datasetlist[0] =~ /^file|http|ftp/) {
		
		# collect the data according to file type
		
		# BigWig Data file
		if ($datasetlist[0] =~ /\.bw$/i) {
			# file is in bigwig format
			# this uses the Bio::DB::BigWig adaptor
			
			# check that we have bigwig support
			$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::bigwig') 
				unless $BIGWIG_OK;
			if ($BIGWIG_OK) {
				# get the dataset scores using Bio::ToolBox::db_helper::bigwig
				
				# the data collection depends on the method
				if ($value_type eq 'score' and 
					$method =~ /min|max|mean|sum|count/
				) {
					# we can use the low-level, super-speedy, summary method 
					# warn " using collect_bigwig_score() with file\n";
					return collect_bigwig_score(
						$chromo,
						$start,
						$stop,
						$method,
						@datasetlist
					);
				}
				
				elsif ($value_type eq 'count' and $method eq 'sum') {
					# we can use the low-level, super-speedy, summary method 
					# warn " using collect_bigwig_score() with file\n";
					return collect_bigwig_score(
						$chromo,
						$start,
						$stop,
						'count', # special method
						@datasetlist
					);
				}
				
				elsif ($method eq 'indexed') {
					# collect hash of position => scores
					# warn " using collect_bigwig_position_score() with file\n";
					return collect_bigwig_position_scores(
						$chromo,
						$start,
						$stop,
						@datasetlist
					);
				}
				
				else {
					# use the longer region collection method
					# warn " using collect_bigwig_scores() with file\n";
					@scores = collect_bigwig_scores(
						$chromo,
						$start,
						$stop,
						@datasetlist
					);
					$dataset_type = 'bw';
				}
			}
			else {
				croak " BigWig support is not enabled! " . 
					"Is Bio::DB::BigWig installed?\n";
			}
		}		
		
		# BigBed Data file
		elsif ($datasetlist[0] =~ /\.bb$/i) {
			# data is in bigbed format
			# this uses the Bio::DB::BigBed adaptor
			
			# check that we have bigbed support
			$BIGBED_OK = _load_helper_module('Bio::ToolBox::db_helper::bigbed') 
				unless $BIGBED_OK;
			if ($BIGBED_OK) {
				# get the dataset scores using Bio::ToolBox::db_helper::bigbed
				
				if ($method eq 'indexed') {
					# warn " using collect_bigbed_position_scores() with file\n";
					return collect_bigbed_position_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type, 
						@datasetlist
					);
				}
				
				else {
					# warn " using collect_bigbed_scores() with file\n";
					@scores = collect_bigbed_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type, 
						@datasetlist
					);
					$dataset_type = 'bb';
				}
			}
			else {
				croak " BigBed support is not enabled! " . 
					"Is Bio::DB::BigBed installed?\n";
			}
		}
		
		# BAM data file
		elsif ($datasetlist[0] =~ /\.bam$/i) {
			# data is in bam format
			# this uses the Bio::DB::Sam adaptor
			
			# check that we have Bam support
			$BAM_OK = _load_helper_module('Bio::ToolBox::db_helper::bam') unless $BAM_OK;
			if ($BAM_OK) {
				# get the dataset scores using Bio::ToolBox::db_helper::bam
				
				if ($method eq 'indexed') {
					# warn " using collect_bam_position_scores() with file\n";
					return collect_bam_position_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type, 
						@datasetlist
					);
				}
				else {
					# warn " using collect_bam_scores() with file\n";
					@scores = collect_bam_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type, 
						@datasetlist
					);
					$dataset_type = 'bam';
				}
			}
			else {
				croak " Bam support is not enabled! " . 
					"Is Bio::DB::Sam installed?\n";
			}
		}
		
		# USeq Data file
		elsif ($datasetlist[0] =~ /\.useq$/i) {
			# data is in useq format
			# this uses the Bio::DB::USeq adaptor
			
			# check that we have bigbed support
			$USEQ_OK = _load_helper_module('Bio::ToolBox::db_helper::useq') 
				unless $USEQ_OK;
			if ($USEQ_OK) {
				# get the dataset scores using Bio::ToolBox::db_helper::useq
				
				if ($method eq 'indexed') {
					# warn " using collect_useq_position_scores() with file\n";
					return collect_useq_position_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type, 
						@datasetlist
					);
				}
				
				else {
					# warn " using collect_useq_scores() with file\n";
					@scores = collect_useq_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type, 
						@datasetlist
					);
					$dataset_type = 'useq';
				}
			}
			else {
				croak " USeq support is not enabled! " . 
					"Is Bio::DB::USeq installed?\n";
			}
		}
		
		# Unsupported Data file
		else {
			confess " Unsupported file type for file '$datasetlist[0]!\n";
		}
		
	}
	
	
	### BigWigSet database
	elsif ($db_type =~ m/^Bio::DB::BigWigSet/) {
		# calling features from a BigWigSet database object
		
		# we may be able to take advantage of a special low-level 
		# super-speedy interface based on the BigWigSet summary feature
		
		# the data collection depends on the method
		if ($value_type eq 'score' and 
			$method =~ /min|max|mean|sum|count/
		) {
			# we can use the low-level, super-speedy, summary method 
			# warn " using collect_bigwigset_score()\n";
			return collect_bigwigset_score(
				$db,
				$chromo,
				$start,
				$stop,
				$strand, 
				$strandedness, 
				$method,
				@datasetlist
			);
		}
		elsif ($value_type eq 'count' and $method eq 'sum') {
			# we can use the low-level, super-speedy, summary method 
			# warn " using collect_bigwigset_score()\n";
			return collect_bigwigset_score(
				$db,
				$chromo,
				$start,
				$stop,
				$strand, 
				$strandedness, 
				'count', # special method
				@datasetlist
			);
		}
		
		elsif ($value_type eq 'score' and $method eq 'indexed') {
			# want positioned score data
			# warn " using collect_bigwigset_position_score()\n";
			return collect_bigwigset_position_scores(
				$db,
				$chromo,
				$start,
				$stop,
				$strand, 
				$strandedness, 
				@datasetlist
			);
		}
		
		else {
			# simply collect a list of the scores
			# warn " using collect_bigwigset_scores()\n";
			@scores = collect_bigwigset_scores(
				$db,
				$chromo,
				$start,
				$stop,
				$strand, 
				$strandedness, 
				@datasetlist
			);
			$dataset_type = 'bw';
		}
	}
		
	
	### SeqFeature database
	elsif ($db_type =~ m/^Bio::DB::SeqFeature/) {
		# a SeqFeature database
		# normal collection
		$iterator = $db->get_seq_stream(
			-seq_id      => $chromo,
			-start       => $start,
			-end         => $stop,
			-primary_tag => [@datasetlist],
		);
	}
	
	
	### Some other database?
	else {
		# some other Bio::DB database????
		# until I code in every single possibility
		# let's just try a basic features method using whatever the 
		# default type is and hope for the best
		# warn "using other database $db_type!\n";
		confess "no database passed!" unless $db;
		$iterator = $db->get_seq_stream(
			-seq_id      => $chromo,
			-start       => $start,
			-end         => $stop,
		);
	}
		
	
	
	### Process database SeqFeature objects
	if ($iterator) {
		# We have a seqfeature object stream
		# First check whether we're dealing with a datafile pointed to by
		# an attribute tag 
		# Failing that, assume it's the seqfeature objects themselves we want
		
		# collect the first feature
		my $feature = $iterator->next_seq;
		return _return_null($method) unless $feature;
		
		# deal with features that might not be from the chromosome we want
		# sometimes chromosome matching is sloppy and we get something else
		while ($feature->seq_id ne $chromo) {
			$feature = $iterator->next_seq || undef;
		}
		return _return_null($method) unless $feature;
		
		## Wig Data
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
				$WIGGLE_OK = _load_helper_module('Bio::ToolBox::db_helper::wiggle') 
					unless $WIGGLE_OK;
				if ($WIGGLE_OK) {
					# get the dataset scores using Bio::ToolBox::db_helper::wiggle
					
					if ($method eq 'indexed') {
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
						@scores = collect_wig_scores(
							$start,
							$stop,
							$strand, 
							$strandedness, 
							$value_type,
							@features
						);
						$dataset_type = 'wig';
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
				# this uses the Bio::DB::BigWig adaptor
				
				# collect the wigfile paths
				# also check strand while we're at it
				my @wigfiles;
				while ($feature) {
					
					# check if we can take this feature
					if (
						$strandedness eq 'all' # stranded data not requested
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
						# we can take this file, it passes the strand test
						my ($file) = $feature->get_tag_values('wigfile');
						push @wigfiles, "file:$file";
					}
					
					# prepare for next
					$feature = $iterator->next_seq || undef;
				}
				
				# if no wigfiles are found, return empty handed
				# should only happen if the strands don't match
				return _return_null($method) unless (@wigfiles);
				
				# check that we have bigwig support
				$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::bigwig') 
					unless $BIGWIG_OK;
				if ($BIGWIG_OK) {
					
					# the data collection depends on the method
					if ($value_type eq 'score' and 
						$method =~ /min|max|mean|sum|count/
					) {
						# we can use the low-level, super-speedy, summary method 
						# warn " using collect_bigwig_score() from tag\n";
						return collect_bigwig_score(
							$chromo,
							$start,
							$stop,
							$method,
							@wigfiles
						);
					}
					
					elsif ($value_type eq 'count' and $method eq 'sum') {
						# we can use the low-level, super-speedy, summary method 
						# warn " using collect_bigwig_score() from tag\n";
						return collect_bigwig_score(
							$chromo,
							$start,
							$stop,
							'count', # special method
							@wigfiles
						);
					}
					
					elsif ($method eq 'indexed') {
						# warn " using collect_bigwig_position_scores() from tag\n";
						return collect_bigwig_position_scores(
							$chromo,
							$start,
							$stop,
							@wigfiles
						);
					}
					
					else {
						# use the longer region collection method
						# warn " using collect_bigwig_scores() from tag\n";
						@scores = collect_bigwig_scores(
							$chromo,
							$start,
							$stop,
							@wigfiles
						);
						$dataset_type = 'bw';
					}
				}
				else {
					croak " BigWig support is not enabled! " . 
						"Is Bio::DB::BigWig installed?\n";
				}
			}
			else {
				croak " Unrecognized wigfile attribute '$wigfile'!" . 
					" Unable to continue!\n";
			}
		}
		
		
		## Database Data
		else {
			# Working with data stored directly in the database
			# this is more straight forward in collection
			
			# Walk through the datapoints
			# warn " using database\n";
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
					if ($method eq 'indexed') {
					
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
							push @scores, 1;
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
				$method eq 'indexed' and 
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
			
			$dataset_type = 'db';
		} 
	
	} # end database collection
	
	
	
	### Determine region score from collected scores
	# We have collected the positioned scores
	# Now return the appropriate values
	
	# indexed scores
	if ($method eq 'indexed') {
		# requested indexed position scores
		# we will simply return the data hash
		# regardless whether there are scores or not
		return %pos2data;
	}
	
	# scores
	elsif ($method eq 'scores') {
		# just the scores are requested
		# return an array reference
		return \@scores;
	}
	
	# single region score
	else {
		# check that we have scores
		return _return_null($method) unless (@scores);
		
		# requested a single score for this region
		# we need to combine the data
		my $region_score;
		
		# first deal with log2 values if necessary
		if ($log) {
			@scores = map {2 ** $_} @scores;
		}
		
		# determine the region score according to method
		# we are using subroutines from Statistics::Lite
		if ($method eq 'median') {
			# take the median value
			$region_score = median(@scores);
		}
		elsif ($method eq 'mean') {
			# or take the mean value
			$region_score = mean(@scores);
		} 
		elsif ($method eq 'range') {
			# or take the range value
			# this is 'min-max'
			$region_score = range(@scores);
		}
		elsif ($method eq 'stddev') {
			# or take the standard deviation value
			# we are using the standard deviation of the population, 
			# since these are the only scores we are considering
			$region_score = stddevp(@scores);
		}
		elsif ($method eq 'min') {
			# or take the minimum value
			$region_score = min(@scores);
		}
		elsif ($method eq 'max') {
			# or take the maximum value
			$region_score = max(@scores);
		}
		elsif ($method eq 'count') {
			# count the number of values
			$region_score = scalar(@scores);
		}
		elsif ($method eq 'sum') {
			# sum the number of values
			$region_score = sum(@scores);
		}
		elsif ($method =~ /rpk?m/) {
			# convert to reads per million mapped
			# this is only supported by bam and bigbed db, checked above
			
			# total the number of reads if necessary
			unless (exists $total_read_number{$dataset} ) {
				
				# check the type of database
				if ($dataset_type eq 'bam') {
					# a bam database
					
					$total_read_number{$dataset} = 
						sum_total_bam_alignments($dataset);
					print "\n [total alignments: ", 
						format_with_commas( $total_read_number{$dataset} ), 
						"]\n";
				}
				
				elsif ($dataset_type eq 'bb') {
					# bigBed database
					
					$total_read_number{$dataset} = 
						sum_total_bigbed_features($dataset);
					print "\n [total features: ", 
						format_with_commas( $total_read_number{$dataset} ), 
						"]\n";
				}
			}	
			
			# calculate the region score according to the method
			if ($method eq 'rpkm') {
				$region_score = 
					( sum(@scores) * 1000000000 ) / 
					( ($stop - $start + 1) * $total_read_number{$dataset} );
			}
			elsif ($method eq 'rpm') {
				$region_score = 
					( sum(@scores) * 1000000 ) / $total_read_number{$dataset};
			}
			else {
				# this dataset doesn't support rpm methods
				# use the sum method instead
				$region_score = sum(@scores);
			}
		}
		else {
			# somehow bad method snuck past our checks
			confess " unrecognized method '$method'!";
		}
	
		# convert back to log2 if necessary
		if ($log) { 
			$region_score = log($region_score) / LOG2;
		}
		
		# finished
		return $region_score;
	}
}


=item _return_null

Internal method for returning a 0 or internal null '.' character based 
on the method being used.

=cut

sub _return_null {
	my $method = shift;
	
	# return based on the method
	if ($method eq 'sum') { 
		return 0;
	}
	elsif ($method eq 'count') { 
		return 0;
	}
	else {
		# internal null value
		return '.';
	}
}


sub _load_helper_module {
	my $class = shift;
	my $success = 0;
	eval {
		load $class; # using Module::Load to dynamically load modules
		$class->import; # must be done particularly for my modules
		$success = 1;
	};
	return $success;
}



=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  

