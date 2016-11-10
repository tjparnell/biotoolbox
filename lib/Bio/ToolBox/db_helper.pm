package Bio::ToolBox::db_helper;
our $VERSION = '1.50';

use strict;
require Exporter;
use Carp qw(carp cluck croak confess);
use Module::Load; # for dynamic loading during runtime
use List::Util qw(min max sum);
use Statistics::Lite qw(median range stddevp);
use Bio::ToolBox::db_helper::config;
use Bio::ToolBox::utility;
use constant {
	CHR  => 0,  # chromosome
	STRT => 1,  # start
	STOP => 2,  # stop
	STR  => 3,  # strand
	STND => 4,  # strandedness
	METH => 5,  # method
	RETT => 6,  # return type
	DB   => 7,  # database object
	DATA => 8,  # first dataset, additional may be present
};


# check values for dynamically loaded helper modules
# these are loaded only when needed during runtime to avoid wasting resources
our $BAM_OK      = 0;
our $BIGBED_OK   = 0;
our $BIGWIG_OK   = 0;
our $SEQFASTA_OK = 0;
our $USEQ_OK     = 0;

# define reusable variables
our $TAG_EXCEPTIONS; # for repeated use with validate_included_feature()
our %total_read_number; # for rpm calculations
our $primary_id_warning; # for out of date primary IDs
our %OPENED_DB; # cache for some opened Bio::DB databases
our %DB_METHODS; # cache for database score methods

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
	get_db_feature 
	get_segment_score 
	calculate_score
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
	  get_db_feature
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
# 		cluck 'no database name passed!';
		return;
	}
	
	# first check if it is a database reference
	my $db_ref = ref $database;
	if ($db_ref =~ /DB/) {
		# the provided database is already an open database object
		# nothing to open, return as is
		return $database;
	}
	
	# check to see if we have already opened it
	if (exists $OPENED_DB{$database} and not $no_cache) {
		# return the cached database connection
		# but NOT if user explicitly requested no cached databases
		# DO NOT reuse database objects if you have forked!!! Bad things happen
		return $OPENED_DB{$database};
	}
	
	# skip parsed databases
	return if $database =~ /^Parsed:/; # we don't open parsed annotation files
	
	
	### Attempt to open the database
	# we go through a series of checks to determine if it is remote, local, 
	# an indexed big data file, SQLite file, etc
	# when all else fails, try to open a SQL connection
	my $db;
	my $error;
	
	# check if it is a remote file
	if ($database =~ /^(?:https?|ftp)/i) {
		
		# a remote Bam database
		if ($database =~ /\.bam$/i) {
			# open using Bam adaptor
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
		elsif ($database =~ /\.(?:bb|bigbed)$/i) {
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
		elsif ($database =~ /\.(?:bw|bigwig)$/i) {
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
			$SEQFASTA_OK = _load_helper_module('Bio::ToolBox::db_helper::seqfasta') 
				unless $SEQFASTA_OK;
			if ($SEQFASTA_OK) {
				$db = open_fasta_db($database);
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
	elsif ($database =~ /gff|bw|bb|bam|useq|db|sqlite|fa|fasta|bigbed|bigwig/i) {
		
		# first check that it exists
		if (-e $database) {
		
			# a Bam database
			if ($database =~ /\.bam$/i) {
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
			elsif ($database =~ /\.(?:bb|bigbed)$/i) {
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
			elsif ($database =~ /\.(?:bw|bigwig)$/i) {
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
				$SEQFASTA_OK = _load_helper_module('Bio::ToolBox::db_helper::seqfasta') 
					unless $SEQFASTA_OK;
				if ($SEQFASTA_OK) {
					$db = open_fasta_db($database);
					unless ($db) {
						$error .= " ERROR: could not open fasta file '$database'!\n";
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
			
			# a gff3 or sqlite database 
			elsif ($database =~ /\.(?:gff3?|gff3?\.gz|db|sqlite)$/) {
				$SEQFASTA_OK = _load_helper_module('Bio::ToolBox::db_helper::seqfasta') 
					unless $SEQFASTA_OK;
				if ($SEQFASTA_OK) {
					$db = open_store_db($database);
					unless ($db) {
						$error = " ERROR: could not load SeqFeature database file '$database'!\n";
					}
				}
				else {
					$error .= " Module Bio::DB::SeqFeature::Store is required to load " . 
						"GFF and database files\n";
				}
			}
			
			else {
				$error .= " Programmer error! Cannot interpret database type for $database!\n";
			}
		}
		
		# file does not exist or can be read
		else {
			# file does not exist
			$error = " ERROR: file '$database' does not exist!\n";
		}
	}
	
	# unrecognized real file
	elsif (-e $database) {
		# file exists, I just don't recognize the extension
		$error = " File '$database' type and/or extension is not recognized\n";
	}
	
	# otherwise assume the name of a SeqFeature::Store database in the configuration
	
	# attempt to open using information from the configuration
	# using default connection information as necessary
	unless ($db) {
		$SEQFASTA_OK = _load_helper_module('Bio::ToolBox::db_helper::seqfasta') 
			unless $SEQFASTA_OK;
		if ($SEQFASTA_OK) {
			$db = open_store_db($database);
			unless ($db) {
				$error .= " unable to open relational database named '$database'\n";
			}
		}
		else {
			$error .= " Module Bio::DB::SeqFeature::Store is required to connect " . 
				"to databases\n";
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
		$error .= " no database could be found or connected!\n";
		print STDERR $error;
		return;
	}
}


### Retrieve a database name from an db object

=item get_db_name

This subroutine will attempt to get the name of an opened Database 
object if for some reason it's unknown, i.e. user only provided an 
opened db object. Only works for some databases, at least those I've 
bothered to track down and find a usable API call to use.

=cut

sub get_db_name {
	my $db = shift;
	return unless $db;
	my $db_ref = ref($db);
	return $db unless $db_ref; # presumption that non-object is just a database name 
	my $db_name;
	if ($db_ref =~ /^Bio::DB::SeqFeature::Store/) {
		# a SeqFeature database, using any DBI adapter
		$db_name = $db->{'dbh'}->{'name'}; 
			# dig through the object internals to identify the original 
			# name of the database
			# this should be relatively well documented through DBI
			# but could break in the future since it's not official API
	}
	elsif ($db_ref eq 'Bio::DB::Sam') {
		# a Bam database
		$db_name = $db->{'bam_path'};
	}
	# determining the database name from other sources is
	# either not possible or not easy, so won't bother unless
	# there is a really really good need
	return $db_name;
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
	my $db = open_db_connection($database);
	my $db_name = get_db_name($database);
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	my $db_ref = ref $db;
	
	# process the database types, according to the type of database
	my @types;
	
	# a SeqFeature database
	if ($db_ref =~ /^Bio::DB::SeqFeature::Store/) {
		
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
the parse_list() method from L<Bio::ToolBox::utility>, so a mix of single 
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
  limit    => Optionally provide a word or regular expression that matches
              the feature type (primary_tag only; source_tag, if present, 
              is ignored). For example, provide "gene|mrna" to only 
              present gene and mRNA features to the user. This is only 
              applicable when a user must select from a database list. 
              The default is to list all available feature types.

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
			elsif ($dataset =~ /\.(?:bam|bw|bigwig|bb|bigbed|useq)$/i) {
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
				if ($limit) {
					my ($p, $s) = split /:/, $type; # split into primary_tag and source
					# only keep those types that match our limiter
					next unless $p =~ /$limit/i;
				}
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
and have an appropriate prefix (file, http, or ftp). 

For multi-threaded execution pass a second value, the number of CPU cores 
available to count Bam files. This will speed up counting bam files 
considerably. The default is 2 for environments where Parallel::ForkManager 
is installed, or 1 where it is not.

=cut

sub check_dataset_for_rpm_support {
	
	# get passed dataset and databases
	my $dataset  = shift;
	my $cpu      = shift;
	
	# Calculate the total number of reads
	# this uses the global variable $rpkm_read_sum
	my $rpm_read_sum;
	
	if (exists $total_read_number{$dataset}) {
		# this dataset has already been summed
		# no need to do it again
		$rpm_read_sum = $total_read_number{$dataset};
	}
	
	elsif ($dataset =~ /\.bam$/) {
		# a bam file dataset
		
		$BAM_OK = _load_helper_module('Bio::ToolBox::db_helper::bam') unless $BAM_OK;
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
		
		$BIGBED_OK = _load_helper_module('Bio::ToolBox::db_helper::bigbed') unless 
			$BIGBED_OK; 
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
	
	# check data object
	my $data = $args{data} || undef;
	unless ($data) {
		confess "must pass a 'data' key and Bio::ToolBox::Data object!";
		return;
	}
	unless (ref($data) eq 'Bio::ToolBox::Data') {
		confess 'must pass a Bio::ToolBox::Data object!';
		return;
	}
	
	# Open a db connection 
	$args{'db'} ||= undef;
	my $db = open_db_connection($args{'db'});
	my $db_name = get_db_name($args{'db'});
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
	
	# Add table columns
	my $pid_i  = $data->add_column('Primary_ID');
	my $name_i = $data->add_column('Name');
	my $type_i = $data->add_column('Type');
	
	# List of types
	if (scalar @classes > 1) {
		$data->metadata($name_i, 'include', join(",", @classes));
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
		$data->add_row( [
			$feature->primary_id, 
			$feature->display_name,  
			$feature->type,
		] );
	}
	
	# print result of search
	printf "   Found %s features in the database\n", format_with_commas($total_count);
	printf "   Kept %s features.\n", format_with_commas($data->last_row);
	
	# return the new data structure
	return 1;
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
	
	# check data object
	my $data = $args{data} || undef;
	unless ($data) {
		confess "must pass a 'data' key and Bio::ToolBox::Data object!";
		return;
	}
	unless (ref($data) eq 'Bio::ToolBox::Data') {
		confess 'must pass a Bio::ToolBox::Data object!';
		return;
	}
		
	# Open a db connection 
	$args{'db'} ||= undef;
	my $db = open_db_connection($args{'db'});
	my $db_name = get_db_name($args{'db'});
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
		
	
	# Prepare data structures
	my $chr_i   = $data->add_column('Chromosome');
	my $start_i = $data->add_column('Start');
	my $stop_i  = $data->add_column('Stop');
	$data->metadata($start_i, 'win', $args{'win'}); 
	$data->metadata($start_i, 'step', $args{'step'});
	
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
			my $end = $start + $args{win} - 1; 
			if ($end > $length) {
				# fix end to the length of the chromosome
				$end = $length;
			} 
			$data->add_row( [ $chr, $start, $end] );
		}
	}
	print "   Kept " . $data->{'last_row'} . " windows.\n";
	
	# Return the data structure
	return 1;
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




=item get_db_feature

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

sub get_db_feature {
	my %args = @_;
	
	# Open a db connection 
	$args{db} ||= undef;
	unless ($args{db}) {
		croak 'no database provided for getting a feature!';
	}
	my $db = open_db_connection($args{db});
	my $db_ref = ref $db;
	
	# get the name of the feature
	my $name = $args{'name'} || undef; 
	
	# check for values and internal nulls
	$args{'id'} = exists $args{'id'} ? $args{'id'} : undef;
	$args{'type'} ||= undef;
	undef $name         if $name and $name eq '.';
	undef $args{'id'}   if $args{'id'} and $args{'id'} eq '.';
	undef $args{'type'} if $args{'type'} and $args{'type'} eq '.';
	
	# quick method for feature retrieval
	if (defined $args{'id'} and $db->can('fetch')) {
		# we have a unique primary_id to identify the feature
		# usually this pulls out the feature directly
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
	return unless $name; # otherwise db will return all features! Not what we want!
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
	unless (@features and $name =~ /[;,\|]/) {
		# I used to append aliases to the end of the name in old versions of biotoolbox
		# crazy, I know, but just in case this happened, let's try breaking them apart
		my $name2 = (split(/\s*[;,\|]\s*/, $name))[0];
			 # multiple names present using common delimiters ;,|
			 # take the first name only, assume others are aliases that we don't need
			@features = $db->features( 
				-name       => $name2,
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
  avoid    => Provide an array reference of database feature types 
              that should be avoided. Any positioned scores which 
              overlap the other feature(s) are not returned. The 
              default is to return all values. A boolean value of 
              1 can also be passed, in which case the same type 
              of feature as the search feature will be used. This
              was the original implementation (v1.35 and below). 
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



sub get_segment_score {
	# parameters passed as an array
	# we will be passing this array on as a reference to the appropriate 
	# imported helper subroutine
	# chromosome, start, stop, strand, strandedness, method, return type, db, dataset
	confess "incorrect number of parameters passed!" unless scalar @_ >= 9;
	
	# check the database
	$_[DB] = open_db_connection($_[DB]) if ($_[DB] and not ref($_[DB]));
	
	# determine method
	my $db_method = $DB_METHODS{$_[METH]}{$_[RETT]}{$_[DB]}{$_[DATA]} || 
		_lookup_db_method(\@_);
	
	# return type values
		# 0 = calculate score
		# 1 = score array
		# 2 = hash position scores
	# immediately return calculated score if appropriate
	if ($_[RETT] > 0) {
		# immediately return either indexed hash or array of scores
		return &{$db_method}(\@_);
	}
	else {
		# calculate a score 
		my $scores = &{$db_method}(\@_);
		# this might be an array reference of scores that need to be combined
		# or it could be a single scalar score which just needs to be returned
		if (ref($scores)) {
			return calculate_score($_[METH], $scores);
		}
		else {
			return $scores;
		}
	}
}


sub calculate_score {
	my ($method, $scores) = @_;
	
	# empty score
	if (not defined $scores or scalar @$scores == 0) {
		return 0 if $method =~ /count|sum/;
		return '.';	
	}
	
	# calculate a score based on the method
	if ($method eq 'mean') {
		my $n = scalar(@$scores);
		return 0 if $n == 0;
		return $scores->[0] if $n == 1;
		return sum(@$scores)/$n;
	} 
	elsif ($method eq 'sum') {
		return sum(@$scores);
	}
	elsif ($method eq 'median') {
		return median(@$scores);
	}
	elsif ($method eq 'min') {
		return min(@$scores);
	}
	elsif ($method eq 'max') {
		return max(@$scores);
	}
	elsif ($method eq 'count' or $method eq 'pcount') {
		return scalar(@$scores);
	}
	elsif ($method eq 'ncount') {
		# Convert names into unique counts
		my %name2count;
		foreach (@$scores) { $name2count{$_} += 1 }
		return scalar(keys %name2count);
	}
	elsif ($method eq 'range') {
		# the range value is 'min-max'
		return range(@$scores);
	}
	elsif ($method eq 'stddev') {
		# we are using the standard deviation of the population, 
		# since these are the only scores we are considering
		return stddevp(@$scores);
	}
	elsif ($method =~ /rpk?m/) {
		confess " The rpm methods have been deprecated due to complexity and " .
			"the variable way of calculating the value. Collect counts and " . 
			"calculate your preferred way.\n";
	}
	else {
		# somehow bad method snuck past our checks
		confess " unrecognized method '$method'!";
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
	my $db = open_db_connection($database);
	my $db_name = get_db_name($database);
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
	
	# SeqFeature::Store
	if (ref $db =~ /^Bio::DB::SeqFeature::Store/) {
		for my $chr ($db->seq_ids) {
			
			# check for excluded chromosomes
			next if (exists $excluded_chr_lookup{$chr} );
			
			# get chromosome size
			my ($seqf) = $db->get_features_by_name($chr);
			my $length = $seqf ? $seqf->length : 0;
			
			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}
	
	# Bigfile
	elsif (ref $db eq 'Bio::DB::BigWig' or ref $db eq 'Bio::DB::BigBed') {
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
			next if (exists $excluded_chr_lookup{$chr} );
			
			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}
	
	# Fasta
	elsif (ref $db eq 'Bio::DB::Fasta') {
		for my $chr ($db->get_all_ids) {
			
			# check for excluded chromosomes
			next if (exists $excluded_chr_lookup{$chr} );
			
			# get chromosome size
			my $seq = $db->get_Seq_by_id($chr);
			my $length = $seq ? $seq->length : 0;
			
			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}
	
	# other Bioperl
	else {
		foreach my $chr ($db->seq_ids) {
			
			# check for excluded chromosomes
			next if (exists $excluded_chr_lookup{$chr} );
			
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

sub _lookup_db_method {
	# parameters passed as an array reference
	my $param = shift;
	
	# determine the appropriate score method
	my $score_method;
	if ($param->[DATA] =~ /^file|http|ftp/) {
		# collect the data according to file type
		
		# BigWig Data file
		if ($param->[DATA] =~ /\.(?:bw|bigwig)$/i) {
			# file is in bigwig format
			# get the dataset scores using Bio::ToolBox::db_helper::bigwig
			# this uses the Bio::DB::BigWig adaptor
			
			# check that we have bigwig support
			$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::bigwig') 
				unless $BIGWIG_OK;
			if ($BIGWIG_OK) {
				if ($param->[RETT] == 2) {
					$score_method = \&collect_bigwig_position_scores;
				}
				elsif ($param->[METH] =~ /min|max|mean|sum|count/) {
					$score_method = \&collect_bigwig_score;
				}
				else {
					$score_method = \&collect_bigwig_scores;
				}
			}
			else {
				croak " BigWig support is not enabled! " . 
					"Is Bio::DB::BigWig installed?\n";
			}
		}		
		
		# BigBed Data file
		elsif ($param->[DATA] =~ /\.(?:bb|bigbed)$/i) {
			# data is in bigbed format
			# get the dataset scores using Bio::ToolBox::db_helper::bigbed
			# this uses the Bio::DB::BigBed adaptor
			
			# check that we have bigbed support
			$BIGBED_OK = _load_helper_module('Bio::ToolBox::db_helper::bigbed') 
				unless $BIGBED_OK;
			if ($BIGBED_OK) {
				if ($param->[RETT] == 2) {
					$score_method = \&collect_bigbed_position_scores;
				}
				else {
					$score_method = \&collect_bigbed_scores;
				}
			}
			else {
				croak " BigBed support is not enabled! " . 
					"Is Bio::DB::BigBed installed?\n";
			}
		}
		
		# BAM data file
		elsif ($param->[DATA] =~ /\.bam$/i) {
			# data is in bam format
			# get the dataset scores using Bio::ToolBox::db_helper::bam
			# this uses the Bio::DB::Sam adaptor
			
			# check that we have Bam support
			$BAM_OK = _load_helper_module('Bio::ToolBox::db_helper::bam') unless $BAM_OK;
			if ($BAM_OK) {
				$score_method = \&collect_bam_scores;
			}
			else {
				croak " Bam support is not enabled! " . 
					"Is Bio::DB::Sam installed?\n";
			}
		}
		
		# USeq Data file
		elsif ($param->[DATA] =~ /\.useq$/i) {
			# data is in useq format
			# get the dataset scores using Bio::ToolBox::db_helper::useq
			# this uses the Bio::DB::USeq adaptor
			
			# check that we have bigbed support
			$USEQ_OK = _load_helper_module('Bio::ToolBox::db_helper::useq') 
				unless $USEQ_OK;
			if ($USEQ_OK) {
				if ($param->[RETT] == 2) {
					$score_method = \&collect_useq_position_scores;
				}
				else {
					$score_method = \&collect_useq_scores;
				}
			}
			else {
				croak " USeq support is not enabled! " . 
					"Is Bio::DB::USeq installed?\n";
			}
		}
		
		# Unsupported Data file
		else {
			confess sprintf " Unsupported or unrecognized file type for file '%s'!\n", 
				$param->[DATA];
		}
	}
	
	
	### BigWigSet database
	elsif (ref($param->[DB]) =~ m/^Bio::DB::BigWigSet/) {
		# calling features from a BigWigSet database object
		
		# check that we have bigwig support
		# duh! we should, we probably opened the stupid database!
		$BIGWIG_OK = _load_helper_module('Bio::ToolBox::db_helper::bigwig') 
			unless $BIGWIG_OK;
		croak " BigWigSet support is not enabled! Is Bio::DB::BigFile installed?" 
			unless $BIGWIG_OK;
		
		# the data collection depends on the method
		if ($param->[RETT] == 2) {
			$score_method = \&collect_bigwigset_position_scores;
		}
		elsif ($param->[METH] =~ /min|max|mean|sum|count/) {
			$score_method = \&collect_bigwigset_score;
		}
		else {
			$score_method = \&collect_bigwigset_scores;
		}
	}
		

	### BioPerl style database
	elsif (ref($param->[DB]) =~ m/^Bio::DB/) {
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
			unless ($param->[DB]->can('get_seq_stream')) {
				confess sprintf "unsupported database! cannot use %s as it does not support get_seq_stream method or equivalent", 
					ref($param->[DB]);
			}
			$score_method = \&collect_store_scores;
		}
		else {
			croak " SeqFeature Store support is not enabled! " . 
				"Is BioPerl and Bio::DB::SeqFeature::Store properly installed?\n";
		}
	}
	
	
	### Some other database?
	else {
		confess "no recognizeable dataset provided!" unless $param->[DATA];
		confess "no database passed!" unless $param->[DB];
		confess "something went wrong! parameters: @$param";
	}
		
	
	### Cache and return the results
	$DB_METHODS{$param->[METH]}{$param->[RETT]}{$param->[DB]}{$param->[DATA]} = 
		$score_method;
	return $score_method;
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
it under the terms of the Artistic License 2.0.  

