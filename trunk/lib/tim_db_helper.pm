package tim_db_helper;

require Exporter;
use strict;
use Carp;
use FindBin qw($Bin);
use Config::Simple;
use Bio::DB::SeqFeature::Store;
use Statistics::Lite qw(
	mean
	median
	min
	max
	range
	stddevp
);


# check whether these optional modules are available
our $WIGGLE_OK = 0;
our $BIGWIG_OK = 0;
eval { 
	# check for wiggle support
	use Bio::Graphics::Wiggle;
	$WIGGLE_OK = 1;
}; 
eval { 
	# check for BigWig support
	use Bio::DB::BigWig;
	$BIGWIG_OK = 1;
}; 

# Hashes of opened file objects
our %OPENED_WIGFILES; # opened wigfile objects
	# in empirical testing, this doesn't really seem to speed things up
	# like I thought it would
	# oh well, keep it anyway????

# Configuration File import for database info
our $TIM_CONFIG;
if (exists $ENV{'TIM_DB_HELPER'}) {
	 $TIM_CONFIG = Config::Simple->new($ENV{'TIM_DB_HELPER'}) or 
		die Config::Simple->error();
}	
elsif (-e "$ENV{HOME}/tim_db_helper.cfg") {
	 $TIM_CONFIG = Config::Simple->new("$ENV{HOME}/tim_db_helper.cfg") or
	 	die Config::Simple->error();
}
else {
	warn "\n#### Using default configuration file '$Bin/../lib/tim_db_config.cfg'####\n";
	$TIM_CONFIG = Config::Simple->new("$Bin/../lib/tim_db_helper.cfg") or 
	 	die Config::Simple->error();
}
our $TAG_EXCEPTIONS; # for repeated use with validate_included_feature()

# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw();
our @EXPORT_OK = qw(
	$TIM_CONFIG
	open_db_connection
	get_dataset_list 
	validate_dataset_list 
	get_new_feature_list 
	get_new_genome_list 
	validate_included_feature 
	get_feature_dataset 
	get_genome_dataset 
	get_chromo_region_score 
	get_region_dataset_hash 
);


# The true statement
1; 

=head1 NAME

tim_db_helper

=head1 DESCRIPTION

These are helper subroutines to work with microarray and/or next generation
sequencing data stored in a Bioperl SeqFeature Store database. The included
subroutines allow to connect to the database and retrieve data relative to
various genomic features within the database, including genes, transcripts,
transcription start sites, genomic bins, etc. The data may collected and
summarized in a variety of methods, including mean, median, enumeration,
etc.

The data may be stored in the database in a variety of mechanisms. The
simplest is to store the data directly into the database as the score value
of the GFF features stored in the database.   Alternatively, the data may
be stored in a binary wig files and referenced in the attribute tag of the
data's feature. Two types of wig files are supported: a scaled binary file
(.wib file) supported by the module L<Bio::Graphics::Wiggle>, and a binary
BigWig file (.bw file) supported by the module L<Bio::DB::BigWig>. The
BigWig file format is much preferred as it maintains spatial resolution of
the original data and does not lose precision by scaling to 8-bit values,
unlike the .wib file format.

While these functions may appear to be simply a rehashing of the methods
and functions in Bio::DB::SeqFeature::Store, they either provide a simpler
function to often used database methodologies or are designed to work
intimately with the tim data format file and data structures (see
L<tim_file_helper.pm>). One key advantage to these functions is the ability
to work with datasets that are stranded (transcriptome data, for example).

A note on the database. While the Bio::DB::SeqFeature::Store database is
designed to work with GFF3 - formatted data, these functions make some
assumptions that we are working with non-standard GFF3 data. Specifically,
the GFF3 feature's method (third column) is supposed to be a standard
Sequence Ontology term. I'm instead using this column to identify different
datasets. For example, technically one should use a method of
'microarray_oligo' and use the feature name as the identifier of the
dataset, such as 'RSC_ypd_244k'. Instead, these functions assume the method
is being used as the unique identifier, so 'RSC_ypd_244k' is set as both
the method and the feature name in the source GFF3 files. While it would be
possible to use the display_name entirely to select unique datasets, it is
more convenient and accessible to identify through the method.

Historically, this module was initially written to use Bio::DB::GFF for 
database usage. It has since been re-written to use Bio::DB::SeqFeature::Store.

To accomodate multiple different databases and settings, database 
configurations are stored in a separate 'tim_db_helper.cfg' file. This is a 
simple INI style text file that stores various variables, including 
database connection parameters and feature groups. The file may be located 
in your home root directory, or located anywhere and referenced in your 
Environment settings under the key 'TIM_DB_HELPER'. If the file is not 
found then the default file located in the biotoolbox lib directory is 
used. See the internal documentation of the tim_db_helper.cfg file for more 
details. 

Complete usage and examples for the functions are provided below.

=head1 USAGE

Call the module at the beginning of your perl script and include the module 
names to export. 

  
  use tim_db_helper qw(
	  get_new_feature_list 
	  get_feature_dataset 
  );
  

This will export the indicated subroutine names into the current namespace. 
Their usage is detailed below. The configuration object may also be imported 
into the program's namespace as C<$TIM_CONFIG> to allow access to the local 
database configuration.


=over

=cut






################################################################################
################           General subroutines             #####################
################################################################################


### Open a connection to the SeqFeature Store MySQL database

=item open_db_connection

This module will open a connection to the Bioperl SeqFeature Store 
database. It returns an object that represents the connection. The database 
may either be a relational database (e.g. MySQL database) or a GFF3 file 
which can be loaded into a small in-memory database.

Pass the name of a relational database or the name of the GFF3 file to be 
loaded into memory. Other parameters for connecting to the database are 
stored in a configuration file, C<tim_db_helper.cfg>. These include 
database adaptors, user name, password, etc.

Example:

	my $db_name = 'cerevisiae';
	my $db = open_db_connection($db_name);



=cut

sub open_db_connection {
	my $database = shift;
	unless ($database) {
		carp 'no database name passed!';
		return;
	}
	
	# determine whether a database name or a file for memory db
	my $db;
	if (
		$database =~ /\.gff3?(?:\.gz)?$/ and 
		-e $database
	) {
		# it appears database is an actual file
		
		# open using a memory adaptor
		print " Loading file into memory database...\n";
		$db = Bio::DB::SeqFeature::Store->new(
			-adaptor => 'memory',
			-gff     => $database,
		);
	}
	
	else {
		# a name of a relational database
		
		# open the connection using parameters from the configuration file
		# we'll try to use database specific parameters first, else use 
		# the db_default parameters
		my $adaptor = $TIM_CONFIG->param($database . '.adaptor') || 
			$TIM_CONFIG->param('default_db.adaptor');
		my $user = $TIM_CONFIG->param($database . '.user') || 
			$TIM_CONFIG->param('default_db.user');
		my $pass = $TIM_CONFIG->param($database . '.pass') ||
			$TIM_CONFIG->param('default_db.pass');
		my $dsn = $TIM_CONFIG->param($database . '.dsn_prefix') ||
			$TIM_CONFIG->param('default_db.dsn_prefix');
		$dsn .= $database;
		
		# establish the database connection
		$db = Bio::DB::SeqFeature::Store->new(
			-adaptor => $adaptor,
			-dsn     => $dsn,
			-user    => $user,
			-pass    => $pass,
		);
	}
	
	# conditional return
	if ($db) {
		return $db;
	} 
	else {
		carp 'unable to establish a database connection!';
		return;
	}
}



### Retrieve a list of the microrarray data sets from the db

=item get_dataset_list

This subroutine will retrieve a list of the available features stored in the 
database and returns a hash of the feature's GFF types, represented as 
"type:source", corresponding to the third and second GFF columns, respectively.
The hash is keyed with an incrementing number, and the value is the GFF type 
of the dataset. A hash is returned rather than a list to help facilitate 
presenting and having the user select an item from the list. The list of
available features are sorted asciibetically before numbered.

Pass either the name of the database or an established database object.

By default, the list of feature types are filtered by the source. Features 
whose source are listed in the C<source_exclude> array of the 
C<tim_db_helper.cfg> file are excluded from the final hash. These usually 
include sources from official genomic authorities, such as 'SGD', 'GeneDB', 
'UCSC', 'Ensembl', etc. In this way, only special features (e.g. microarray 
datasets) are included in the list.

To include all features without filtering, pass a second true argument 
(1, 'all', etc.).

Example:

	my $db_name = 'cerevisiae';
	my %microarray_dataset = get_dataset_list($db_name);
	foreach (sort {$a <=> $b} keys %microarray_dataset) {
		# print number in ascending order, dataset name
		print $_, $microarray_dataset{$_};
	}
	
	my %all_features = get_dataset_list($db_name, 'all');


=cut

sub get_dataset_list {
	
	my $database = shift;
	my $use_all_features = shift;
	
	# Open a db connection 
	# determine whether we have just a database name or an opened object
	my $db; # the database object to be used
	my $db_name; # the name of the database, for use with config param
	if ($database) {
		my $db_ref = ref $database;
		if ($db_ref =~ /Bio::DB/) {
			# a db object returns the name of the package
			# this appears to be a bioperl db object
			$db = $database;
			$db_name = $db->{'dbh'}->{'name'}; 
				# dig through the object internals to identify the original 
				# name of the database
				# this should be relatively well documented through DBI
				# but could break in the future since it's not official API
		}
		else {
			# the name of a database was passed, create a database connection
			$db_name = $database;
			$db = open_db_connection($db_name);
		}
	}
	else {
		carp 'no database name passed!';
		return;
	}
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	# get sources to skip
		# usually these are features from an official genome authority
	my %source2skip;
	foreach ($TIM_CONFIG->param($db_name . '.source_exclude') ) {
		# database specific source exclusions
		$source2skip{$_} = 1;
	}
	unless (keys %source2skip) {
		# no database specific exclusions, we'll read default then
		foreach ($TIM_CONFIG->param('default_db.source_exclude') ) {
			$source2skip{$_} = 1;
		}
	}
		
	# process the database types
	my %dataset;
	my $i = 1;
	foreach my $type (
		map $_->[1],
		sort {$a->[0] cmp $b->[0]} 
		map [$_->method, $_],
		$db->types
	) {
		# sort the types asciibetically by method
		
		my $source = $type->source;
		
		# add the type to the list
		if ($use_all_features) {
			# all types should be added to the list
			$dataset{$i} = $type;
			$i++;
		}
		else {
			# check and skip unwanted sources
			unless (exists $source2skip{$source}) {
				# keep if it's not on the unwanted list
				$dataset{$i} = $type;
				$i++;
			}
		}
	}
	return %dataset;
}


### Validate a list of microarray data sets against those in a db

=item validate_dataset_list

This subroutine will validate that a list of microarray data set names exist 
within the given database. This is to help with, for example, checking that 
a dataset name written on a commandline is spelled correctly and actually 
exists within the given database. This is why the above subroutine, 
get_dataset_list(), is so helpful as it avoids having to validate existance and
spelling. 

The subroutine is passed an array. The first element of the array must be 
either the name of the database or an established database object reference. 
The subsequent elements are the names of the datasets to be verified.

The subroutine returns a scalar string consisting of the names of the bad
dataset names (to be passed on to the user). The list is separated by a 
comma and space ", ". 

This subroutine could do with a much better interface and needs a re-write and
re-work of the API.

Example:
	
	my @dataset_names = qw(
		H3_ChIP_44k
		H3k4me3_ChIP_44k
		H3k7ac_ChIP_44k   
	); 
	# 3rd dataset should be H3k9ac_ChIP_44k
	my $db_name = 'cerevisiae';
	my $bad_dataset = validate_dataset_list(
		$db_name,
		@dataset_names
	);
	print $bad_dataset; # prints 'H3k7ac_ChIP_44k'

=cut

sub validate_dataset_list {
	my $database = shift;
	unless (scalar @_ > 0) { die "no datasets to validate!\n"}
	
	# Open a db connection 
	# determine whether we have just a database name or an opened object
	my $db; # the database object to be used
	if (defined $database) {
		my $db_ref = ref $database;
		if ($db_ref =~ /Bio::DB/) {
			# a db object returns the name of the package
			# this appears to be a bioperl db object
			$db = $database;
		}
		else {
			# the name of a database
			$db = open_db_connection( $database );
		}
	}
	else {
		carp 'no database name passed!';
		return;
	}
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	# generate a hash of all the data types in the database
	my %datasethash;
	foreach ($db->types) {
		# obtain a list of all data types in the database
		# these are typename objects representing the method and source 
		# of each type of feature in the database
		# we will put these in the hash
		
		# first as a method:source string
		my $type = "$_"; # returns 'method:source'
		$datasethash{$type} = 1; # the value is not important here
		
		# second, put in the just the method, since we often use that alone
		my $primary = $_->method;
		$datasethash{$primary} = 1;
	}
	
	# now go through the list of datasets to check the name
	my @baddatasets; # an array of the names that fail validation
	foreach my $dataset (@_) {
		# we may have combined datasets indicated by a &
		if ($dataset =~ /&/) {
			foreach (split /&/, $dataset) {
				unless (exists $datasethash{$_}) {
					push @baddatasets, $_;
				}
			}
		} else { # only one dataset
			unless (exists $datasethash{$dataset}) {
				push @baddatasets, $dataset;
			}
		}
	}
	
	# return the name of bad datasets
	return join(", ", @baddatasets);
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
tim_file_helper.pm. The data table will have two or three columns. The 
feature name and type:source are listed in columns one and two, respectively.
If the features have an Alias tag, then a third column is included with 
a comma delimited list of the feature aliases.

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  db =>       The name of the database or a reference to an 
              established database object. 
  features => A scalar value containing a name representing the 
              type(s) of feature(s) to collect. This name will be 
              parsed into an actual list with the internal subroutine 
              _features_to_classes(). Refer to that documentation for 
              a list of appropriate features.
  Optional: 
  mito =>     A boolean value (1 or 0) indicating whether features
              from the mitochrondrial genome should be included.
              The default value is false.
  dubious =>  A boolean value (1 or 0) indicating whether genes 
              flagged in the database as 'dubious' should be 
              included. The default is false (not kept).

The subroutine will return a reference to the data hash. It will print 
status messages to STDOUT. 

Example

	my $db_name = 'cerevisiae';
	my %data = get_new_feature_list( {
		'db'        => $db_name,
		'features'  => 'genes',
		'mito'      => 0,
		'dubious'   => 0,
	} );


=cut


sub get_new_feature_list {

	# Retrieve passed values
	my $arg_ref = shift; # the passed argument values as a hash reference
	
	# Open a db connection 
	# determine whether we have just a database name or an opened object
	my $db; # the database object to be used
	my $db_name; # the name of the database, for use with config param
	if (defined $arg_ref->{'db'}) {
		my $db_ref = ref $arg_ref->{'db'};
		if ($db_ref =~ /Bio::DB/) {
			# a db object returns the name of the package
			# this appears to be a bioperl db object
			$db = $arg_ref->{'db'};
			$db_name = $db->{'dbh'}->{'name'}; 
				# dig through the object internals to identify the original 
				# name of the database
				# this should be relatively well documented through DBI
				# but could break in the future since it's not official API
		}
		else {
			# the name of a database was passed, create a database connection
			$db_name = $arg_ref->{'db'};
			$db = open_db_connection($db_name);
		}
	}
	else {
		carp 'no database name passed!';
		return;
	}
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	# Translate the features into a list of classes
	my @classes = _features_to_classes($arg_ref->{'features'});
	unless (@classes) {
		carp "no or unknown features passed!";
		return;
	}
	
	# Check for including mitochondrial genes
	my $mito = $arg_ref->{'mito'} || 0; 
	
	# Generate data structures
	my %data_hash; # the primary data structure that everything will go into
	my @feature_table; # an array of arrays to put the feature list into
	# we're keeping the data table separate for the time being for simplicity
	
	# begin loading basic metadata information
	$data_hash{'db'} = $db_name; # name of the database
	$data_hash{'feature'} = $arg_ref->{'features'};
	$data_hash{'gff'} = 0;
	$data_hash{'number_columns'} = 2;
	
	# Generate the table header
	push @feature_table, [ 
		'Name', 
		'Type', 
	];
	
	# Generate column metadata
	$data_hash{0} = {
			'name'  => 'Name',
			'index' => 0,
	};
	$data_hash{1} = {
			'name'  => 'Type',
			'index' => 1,
	};
	
	# List of types
	if (scalar @classes > 1) {
		$data_hash{1}->{'include'} = join(",", @classes);
	}
	
	# Collect the genes from the database
	print "   Searching for " . join(", ", @classes) . "\n";
	my @featurelist; # an array of found feature objects in the database
	@featurelist = $db->features(
			-types => \@classes
	); 
	unless (@featurelist) {
		# there should be some features found in the database
		carp "no features found in database!";
		return;
	}
	print "   Found " . scalar @featurelist . " features in the database.\n";
	
	
	# Check for aliases
	for (my $i = 0; $i < 50; $i++) {
		# we're checking the first 50 features looking for an Alias tag
		# checking that many because not all features may have the tag
		# we like to have Aliases, because it makes interpreting gene names
		# a little easier
		last unless (defined $featurelist[$i]);
		
		if ($featurelist[$i]->has_tag('Alias')) {
			
			# add an Alias column to the data table
			push @{ $feature_table[0] }, 'Aliases';
			$data_hash{2} = {
					'name'  => 'Aliases',
					'index' => 2,
			};
			$data_hash{'number_columns'} = 3;
			last;
		}
	}
	
	# Process the features
	FEATURE_COLLECTION_LIST:
	foreach my $feature (@featurelist) {
		
		# skip the mitochondrial genes
		unless ($mito) { 
			next FEATURE_COLLECTION_LIST if 
				$feature->seq_id =~ /^chrm|chrmt|mt|mit/i;
		}
		
		# skip anything that matches the tag exceptions
		unless ( validate_included_feature($feature) ) {
			next FEATURE_COLLECTION_LIST;
		}
		
		
		# Record the feature information
		my @data = (
			$feature->display_name, 
			$feature->type 
		);
		
		# Add alias info if available
		if (exists $data_hash{2}) {
			# we only look for Alias info if we have a column for it
			if ($feature->has_tag('Alias')) {
				push @data, join(q( ), $feature->get_tag_values('Alias'));
			}
			else {
				push @data, '.'; # internal null value
			}
		}
		
		# Record information
		push @feature_table, \@data;
	}
	
	
	# print result of search
	print "   Kept " . scalar @feature_table . " features.\n";
	
	# sort the table
	my @feature_table_sorted;
	my $header = shift @feature_table; # move header
	if ($arg_ref->{'features'} =~ m/transcript/i) {
		# sort differently if we're working with transcripts
		
		@feature_table_sorted = sort { 
			# sort first by parent type, then by name
			( $a->[3] cmp $b->[3] ) and ( $a->[0] cmp $b->[0] )
		} @feature_table; 
	}
	else {
		# non-transcript sort
		
		@feature_table_sorted = sort { 
			# sort first by type, then by name
			( $a->[1] cmp $b->[1] ) and ( $a->[0] cmp $b->[0] )
		} @feature_table; 
	}
	undef @feature_table; # empty the original array
	unshift @feature_table_sorted, $header;
	
	# put the feature_table into the data hash
	$data_hash{'data_table'} = \@feature_table_sorted;
	
	# record the last row
	$data_hash{'last_row'} = scalar @feature_table_sorted - 1;
	
	# return the reference to this hash
	return \%data_hash;
}



### Generate a new list genomic windows

=item get_new_genome_list 

This subroutine will generate a new list of genomic windows. The genome
is split into intervals of a specified size that is moved along the 
genome in specified step sizes.

The subroutine will generate and return a data hash as described in 
tim_file_helper.pm. The data table will have 3 columns, including 
Chromosome, Start, and Stop.

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  db =>       The name of the database or a reference to an 
              established database object. 
  Optional: 
  win =>      A scalar value containing an integer representing the
              size of the window in basepairs. Default is 500 bp.
  step =>     A scalar value containing an integer representing the
              step size for advancing the window across the genome. 
              The default is the window size.
  mito =>     A boolean value (1 or 0) indicating whether features
              from the mitochrondrial genome should be included.
              The default value is false.

The subroutine will return a reference to the data hash. It will print 
status messages to STDOUT. 

Example

	my $db_name = 'cerevisiae';
	my $window_size = 500;
	my $step_size = 250;
	my %data = get_new_genome_list( {
		'db'        => $db_name,
		'win'       => $window_size,
		'step'      => $step_size,
		'mito'      => 1,
	} );


=cut

sub get_new_genome_list {

	# Collect the passed arguments
	my $arg_ref = shift; 
	
	# Open a db connection 
	# determine whether we have just a database name or an opened object
	my $db; # the database object to be used
	my $db_name; # the name of the database, for use with config param
	if (defined $arg_ref->{'db'}) {
		my $db_ref = ref $arg_ref->{'db'};
		if ($db_ref =~ /Bio::DB/) {
			# a db object returns the name of the package
			# this appears to be a bioperl db object
			$db = $arg_ref->{'db'};
			$db_name = $db->{'dbh'}->{'name'}; 
				# dig through the object internals to identify the original 
				# name of the database
				# this should be relatively well documented through DBI
				# but could break in the future since it's not official API
		}
		else {
			# the name of a database was passed, create a database connection
			$db_name = $arg_ref->{'db'};
			$db = open_db_connection($db_name);
		}
	}
	else {
		carp 'no database name passed!';
		return;
	}
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	# Determine win and step sizes
	my ($win, $step);
	if ($arg_ref->{'win'}) {
		$win = $arg_ref->{'win'};
	}
	else {
		$win = 500;
		print "  Using default window size 500 bp\n";
	}
	if ($arg_ref->{'step'}) {
		$step = $arg_ref->{'step'};
	}
	else {
		$step = $win;
	}
	
	
	# Generate data structures
	my %data_hash; # the primary data structure that everything will go into
	my @feature_table; # an array of arrays to put the feature list into
	# we're keeping the data table separate for the time being for simplicity
	
	# Begin loading basic metadata information
	$data_hash{'db'} = $db_name; # the db name
	$data_hash{'feature'} = 'genome';
	$data_hash{'number_columns'} = 3;
	$data_hash{'gff'} = 0;
	
	# Prepare column metadata
	# the first three columns are identical regardless of feature type
	# put column headers at the beginning of the array
	push @feature_table, [ qw(
			Chromosome
			Start
			Stop
	) ];
	# column metadata
	$data_hash{0} = {
			'name'  => 'Chromosome',
			'index' => 0,
	};
	$data_hash{1} = {
			'name'  => 'Start',
			'index' => 1,
			'win'   => $win,
			'step'  => $step,
	};
	$data_hash{2} = {
			'name'  => 'Stop',
			'index' => 2,
	};
	
	
	# Collect the chromosomes from the db
		# sort the chromosomes
		# this only works with numbered chromosomes, not Roman numerals
		# additionally skip mitochrondiral chromosome named as chrm or chrmt
		# we're using an elongated pseudo Schwartzian Transform
	my @temp_chromos;
	my $mitochr; # special placeholder for the mitochrondrial chromosome
	my $ref_seq_type = $TIM_CONFIG->param("$db_name\.reference_sequence_type") ||
		$TIM_CONFIG->param('default_db.reference_sequence_type');
			# the annotation gff may have the reference sequences labeled
			# as various types, such as chromosome, sequence, 
			# contig, scaffold, etc
			# this is set in the configuration file
			# this could pose problems if more than one is present
	foreach ( $db->features(-type => $ref_seq_type) ) {
		my $name = $_->display_name;
		my $primary = $_->primary_id;
			# this should be the database primary ID, used in the database 
			# schema to identify features loaded from the initial GFF
			# we're using this for sorting
		if ($name =~ /^chrm|chrmt|mt|mit/i) {
			# skip mitochondrial chromosome, invariably named chrM or chrMT
			$mitochr = $_; # keep it for later in case we want it
			next;
		}
		push @temp_chromos, [$primary, $_];
	}
	my @chromosomes = map $_->[1], sort { $a->[0] <=> $b->[0] } @temp_chromos;
		# the chromosomes should now be in the same order as they occured
		# in the original annotation GFF file
	if (exists $arg_ref->{'mito'} and $arg_ref->{'mito'} and $mitochr) {
		# push mitochrondrial chromosome at the end of the list if requested
		push @chromosomes, $mitochr;
	}
	unless (@chromosomes) {
		carp " no chromosomes found! Please check the reference sequence type\n" .
			" in the configuration file and database\n";
		return;
	}
	
	# Collect the genomic windows
	print "   Generating $win bp windows in $step bp increments\n";
	foreach my $chrobj (@chromosomes) {
		my $chr = $chrobj->name; # chromosome name
		my $length = $chrobj->length;
		for (my $start = 1; $start <= $length; $start += $step) {
			# set the end point
			my $end = $start + $win - 1; 
			
			if ($end > $length) {
				# fix end to the length of the chromosome
				$end = $length;
			} 
			
			# add to the output list
			push @feature_table, [ $chr, $start, $end,];
		}
	}
	print "   Kept " . (scalar @feature_table) . " windows.\n";
	
	
	# Put the feature_table into the data hash
	$data_hash{'data_table'} = \@feature_table;
	
	# Record the last row
	$data_hash{'last_row'} = scalar @feature_table - 1;
	
	# Return the reference to this hash
	return \%data_hash;
}



sub validate_included_feature {
	
	# feature to check
	my $feature = shift;
	
	# get the list of feature exclusion tags
	unless (defined $TAG_EXCEPTIONS) {
		$TAG_EXCEPTIONS = $TIM_CONFIG->get_block('exclude_tags');
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






################################################################################
##################           Score subroutines             #####################
################################################################################


### Retrieve values from a microarray data set for a list of features

=item get_feature_dataset 

This subroutine will retrieve data set values for each item in a feature list. 
The values are numeric, typically from a microarray data set in the database.
The features are genomic features in the database, such as genes or 
transcripts, and should be generated using the get_new_feature_list() 
subroutine. The genomic region corresponding to the entire feature may be used
or adjusted by setting specific parameters. The values within the defined 
region are collected for the specified dataset(s) and combined into a single 
value using the specified method. 

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  db =>       The name of the database or a reference to an 
              established database object. 
  dataset =>  The name(s) of the datasets in the database to be 
              collected. The name(s) should be verified using the 
              subroutine validate_dataset_list() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces.
  data =>     A scalar reference to the data table containing the
              list of features. This should be a reference to the
              key 'data_table' in the data hash described in 
              load_tim_data_file(), not to the entire data hash. 
              Note that the column metadata should be updated 
              separately by the calling program.
  method =>   The method used to combine the dataset values found
              in the defined region. Acceptable values include 
              mean, median, range, stddev, min, max, and enumerate. 
              See _get_segment_score() documentation for more info.
  Optional:
  log =>      Boolean value (1 or 0) indicating whether the dataset 
              values are in log2 space or not. If undefined, the 
              dataset name will be checked for the presence of the 
              phrase 'log2' and set accordingly. This argument is
              critical for accurately mathematically combining 
              dataset values in the region.
  strand =>   Indicate whether the dataset values from a specific 
              strand relative to the feature should be collected. 
              Acceptable values include sense, antisense, no,
              or none. Default is 'none'.
  extend =>   Indicate an integer value representing the number of 
              bp the feature's region should be extended on both
              sides.
  start =>    Indicate an integer value representing the start  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with 'stop'.
  stop =>     Indicate an integer value representing the stop  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with 'start'.
  fstart =>   Indicate the fraction of feature length at which the
              region will start. The value should be a decimal 
              number and not a whole percentage (this is not 
              verified!). Include a negative sign to force upstream 
              of the reference point. Must be combined with 'fstop'.
  fstop =>    Indicate the fraction of feature length at which the
              region will stop. The value should be a decimal 
              number and not a whole percentage (this is not 
              verified!). Include a negative sign to force upstream 
              of the reference point. Must be combined with 'fstart'.
  limit =>    Indicate a minimum size of the feature before "fstart" 
              and "fstop" are utilized, otherwise the entire feature 
              is used. Useful to avoid taking a subfraction of very 
              small features that may be below the resolution of the 
              dataset. Default 1000 bp.
  position => Indicate the relative position of the feature from 
              which the "start" and "stop" positions are calculated.
              Three values are accepted: "5", which denotes the 
              default 5' end of the feature, "3" which denotes the 
              3' end of the feature, or "1" which denotes the 
              middle of the feature. This option is used only in 
              conjunction with "start" and "stop"  or "fstart" and 
              "fstop" options. The default is "5".
          	  
The subroutine will return a true value (1) if successful.

Examples

	my $db_name = 'cerevisiae';
	my $data = $data_hash->{'data_table'};
	
	# entire feature, method average
	my $success = get_feature_dataset( {
		'db'      => $db_name,
		'data'    => $data,
		'method'  => 'mean',
		'dataset' => $dataset,
		'log'     => 1,
	} );
	
	# promoter region, method median
	my $success = get_feature_dataset( {
		'db'      => $db_name,
		'data'    => $data,
		'method'  => 'median',
		'dataset' => $dataset,
		'start'   => -300,
		'stop'    => 100
		'log'     => 1,
	} );
	
	# maximum expression of middle 50% 
	$dataset='Steinmetz_polyA_actDf&Steinmetz_polyA_actDr';
	my $success = get_feature_dataset( {
		'db'      => $db_name,
		'data'    => $data,
		'method'  => 'max',
		'dataset' => $dataset,
		'fstart'  => 0.25,
		'fstop'   => 0.75,
		'strand'  => 'sense',
	} );
	


=cut

sub get_feature_dataset {

	# retrieve passed values
	my $arg_ref = shift; # the passed values as a hash reference
	
	# check required values
	unless ($arg_ref->{'dataset'}) {
		carp " no dataset requested!";
		return;
	}
	unless ($arg_ref->{'method'}) {
		carp " no method defined!";
		return;
	}
	my $data_table_ref = $arg_ref->{'data'};
	unless ($data_table_ref) {
		carp "no feature list defined!";
		return;
	}
	
	
	# define default values as necessary
	my $log = $arg_ref->{'log'};
	unless (defined $log) {
		# we need to know whether we are working with a log2 dataset or not
		# as it will adversely affect the math!
		if ($arg_ref->{'dataset'} =~ /log2/i) {
			# if we're smart we'll encode the log2 status in the dataset name
			# but chances are, we're not that smart
			$log = 1;
		} else {
			# otherwise assume it is non-log
			# unsafe, but what else to do? we'll put the onus on the user
			$log = 0;
		}
	}
	my $stranded = $arg_ref->{'strand'} || 'none';
	my $relative_pos = $arg_ref->{'position'} || 5;
	my $limit = $arg_ref->{'limit'} || 1000;
		# this is an arbitrary default value that seems reasonable for
		# moderate resolution microarrays
		# it's only needed if using the fractional start & stop
	
	# define other variables from the argument hash
	my $extend = $arg_ref->{'extend'} || undef;
	my $start = $arg_ref->{'start'} || undef;
	my $stop = $arg_ref->{'stop'} || undef;
	my $fstart = $arg_ref->{'fstart'} || undef;
	my $fstop = $arg_ref->{'fstop'} || undef;
	
	
	
	# Open a db connection 
	# determine whether we have just a database name or an opened object
	my $db; # the database object to be used
	if (defined $arg_ref->{'db'}) {
		my $db_ref = ref $arg_ref->{'db'};
		if ($db_ref =~ /Bio::DB/) {
			# a db object returns the name of the package
			# this appears to be a bioperl db object
			$db = $arg_ref->{'db'};
		}
		else {
			# the name of a database
			$db = open_db_connection( $arg_ref->{'db'} );
		}
	}
	else {
		carp 'no database name passed!';
		return;
	}
	unless ($db) {
		carp 'no database connected!';
		return;
	}
		
	# Identify the identifying indexes
	my $name_index;
	my $type_index; 
	# we need to identify which columns in the feature table are name and class
	# these should be columns 0 and 1, respectively, but just in case...
	# and to avoid hard coding index values in case someone changes them
	for (my $i = 0; $i < scalar @{ $data_table_ref->[0] }; $i++ ) {
		if ($data_table_ref->[0][$i] =~ /^name$/i) {
			$name_index = $i;
		}
		elsif ($data_table_ref->[0][$i] =~ /^class$/i) {
			$type_index = $i;
		}
		elsif ($data_table_ref->[0][$i] =~ /^type$/i) {
			$type_index = $i;
		}
	}
	unless (defined $name_index and defined $type_index) {
		carp 'unable to identify Name and/or Type columns';
		return;
	}
	
	# generate column name
	push @{ $data_table_ref->[0] }, $arg_ref->{'dataset'};
	
	# number of features in the list
	my $feat_num = scalar @{ $data_table_ref } - 1; 
		# we don't have access to entire data structure, just the data table
		# hence we need to determine the last row directly from the data table
	
	# Begin data collection
	# we will first check for arguments which modify or adjust the size of the
	# feature region. We check this first before stepping through the list to 
	# avoid repeatedly testing at every step. Hence, we have some code 
	# redundancy, since we rewrite the for loop multiple times. I suppose I 
	# could have generalized into another subroutine, but that might take 
	# more effort and overhead
	
	if (defined $extend) {
		# if an extension is specified to the feature region
		
		FEATURE_DATA_COLLECTION:
		for my $i (1 .. $feat_num) {
			# collecting feature identification 
			my $name = $data_table_ref->[$i][$name_index];
			my $type = $data_table_ref->[$i][$type_index];
			
			# first define the the feature
			my @features = $db->features( 
					-name  => $name,
					-type => $type,
			);
			if (scalar @features > 1) {
				# there should only be one feature found
				# if more, there's redundant or duplicated data in the db
				# warn the user, this should be fixed
				warn " Found more than one feature of '$type => $name' in " .
					"the database!\n Using the first feature only!\n";
			}
			my $feature = shift @features; 
			
			# now re-define the region based on the extended coordinates
				# this is regardless of feature orientation
			my $region = $db->segment( 
				-name      => $name,
				-type      => $type,
				-start     => $feature->start - $extend,
				-end       => $feature->end + $extend,
				-strand    => $feature->strand,
				-absolute  => 1,
			);
			
			
			# confirm region
			unless ($region) { 
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$arg_ref->{'dataset'},
						$arg_ref->{'method'}, 
						$stranded, 
						$log,
			);

		}
	
	}
	
	elsif (
			defined $start and 
			defined $stop and
			$relative_pos == 3
	) {
		# if specific start and stop coordinates from the 3'end are requested
		
		FEATURE_DATA_COLLECTION:
		for my $i (1 .. $feat_num) {
			# collecting feature identification 
			my $name = $data_table_ref->[$i][$name_index];
			my $type = $data_table_ref->[$i][$type_index];
			
			# first define the feature to get length
			my @features = $db->features( 
					-name  => $name,
					-type => $type,
			);
			if (scalar @features > 1) {
				# there should only be one feature found
				# if more, there's redundant or duplicated data in the db
				# warn the user, this should be fixed
				warn " Found more than one feature of '$type => $name' in " .
					"the database!\n Using the first feature only!\n";
			}
			my $feature = shift @features; 
			
			# now re-define the region based on the extended coordinates
			# relative to the 3' end
			my $region;
			if ($feature->strand >= 0) {
				# feature is on the forward, watson strand
				$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $feature->end + $start,
					-end       => $feature->end + $stop,
					-strand    => $feature->strand,
					-absolute  => 1,
				);
			}
			elsif ($feature->strand < 0) {
				# feature is on the reverse, crick strand
				$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $feature->start - $stop,
					-end       => $feature->start - $start,
					-strand    => $feature->strand,
					-absolute  => 1,
				);
			}

			
			# confirm region
			unless ($region) { 
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$arg_ref->{'dataset'},
						$arg_ref->{'method'}, 
						$stranded, 
						$log,
			);
		}
	}
	
	elsif (
			defined $start and 
			defined $stop and 
			$relative_pos == 5
	) {
		# if specific start and stop coordinates from the 5'end are requested
		
		FEATURE_DATA_COLLECTION:
		for my $i (1 .. $feat_num) {
			# collecting feature identification 
			my $name = $data_table_ref->[$i][$name_index];
			my $type = $data_table_ref->[$i][$type_index];
			
			# define feature
			my @features = $db->features( 
					-name  => $name,
					-type => $type,
			);
			if (scalar @features > 1) {
				# there should only be one feature found
				# if more, there's redundant or duplicated data in the db
				# warn the user, this should be fixed
				warn " Found more than one feature of '$type => $name' in " .
					"the database!\n Using the first feature only!\n";
			}
			my $feature = shift @features; 
			
			# now re-define the region based on the extended coordinates
			# relative to the 5' end of the feature
			my $region;
			if ($feature->strand >= 0) {
				# feature is on the forward, watson strand
			
				$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $feature->start + $start,
					-end       => $feature->start + $stop,
					-strand    => $feature->strand,
					-absolute  => 1,
				);
			}
			elsif ($feature->strand < 0) {
				# feature is on the reverse, crick strand
			
				$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $feature->end - $stop,
					-end       => $feature->end - $start,
					-strand    => $feature->strand,
					-absolute  => 1,
				);
			}
			
			# confirm region
			unless ($region) { 
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$arg_ref->{'dataset'},
						$arg_ref->{'method'}, 
						$stranded, 
						$log,
			);
		}
	}
	
	
	elsif (
			defined $fstart and 
			defined $fstop and
			$relative_pos == 3
	) {
		# if percentile values for start and stop are specified
		# define region coordinates on 3' end of feature
		# I can't imagine using the percentile - based positions 
		# from the 3' end, but you never know....
			
		FEATURE_DATA_COLLECTION:
		for my $i (1 .. $feat_num) {
			# collecting feature identification 
			my $name = $data_table_ref->[$i][$name_index];
			my $type = $data_table_ref->[$i][$type_index];
			
			# define region coordinates on 3' end of feature
			# first define the region of the feature to get length
			my @features = $db->features( 
					-name  => $name,
					-type => $type,
			);
			if (scalar @features > 1) {
				# there should only be one feature found
				# if more, there's redundant or duplicated data in the db
				# warn the user, this should be fixed
				warn " Found more than one feature of '$type => $name' in " .
					"the database!\n Using the first feature only!\n";
			}
			my $feature = shift @features; 
			my $length = $feature->length;
			
			# Define the region
			my $region;
			# confirm that feature exceeds our minimum size limit
			if ($length >= $limit) {
				# the feature is long enough to fractionate
				my $relative_start = sprintf "%.0f", 
									$length * $fstart;
				my $relative_stop = sprintf "%.0f", 
									$length * $fstop;
				
				# define relative to the 3' end
				if ($feature->strand >= 0) {
					# feature is on the forward, watson strand
					$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $feature->end + $relative_start,
					-end       => $feature->end + $relative_stop,
					-strand    => $feature->strand,
					-absolute  => 1,
					);
				}
				elsif ($feature->strand < 0) {
					# feature is on the reverse, crick strand
					$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $feature->start - $relative_stop,
					-end       => $feature->start - $relative_start,
					-strand    => $feature->strand,
					-absolute  => 1,
					);
				}
			}
			else {
				# feature length is too small to fractionate
				# so we'll simply take the entire region
				$region = $feature->segment;
			}	 

			# confirm region
			unless ($region) { 
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$arg_ref->{'dataset'},
						$arg_ref->{'method'}, 
						$stranded, 
						$log,
			);
		}
	}
					
	elsif (
			defined $fstart and 
			defined $fstop and 
			$relative_pos == 5
	) {
		# if percentile values for start and stop are specified
		# define region coordinates on 5' end of feature
			
		FEATURE_DATA_COLLECTION:
		for my $i (1 .. $feat_num) {
			# collecting feature identification 
			my $name = $data_table_ref->[$i][$name_index];
			my $type = $data_table_ref->[$i][$type_index];
			
			# first define the feature to get length
			my @features = $db->features( 
					-name  => $name,
					-type => $type,
			);
			if (scalar @features > 1) {
				# there should only be one feature found
				# if more, there's redundant or duplicated data in the db
				# warn the user, this should be fixed
				warn " Found more than one feature of '$type => $name' in " .
					"the database!\n Using the first feature only!\n";
			}
			my $feature = shift @features; 
			my $length = $feature->length;
			
			# Define the region
			my $region;
			# confirm that feature exceeds our minimum size limit
			if ($length >= $limit) {
				# the feature is long enough to fractionate
				my $relative_start = sprintf "%.0f", 
									$length * $fstart;
				my $relative_stop = sprintf "%.0f", 
									$length * $fstop;
				
				# define relative to the 3' end
				if ($feature->strand >= 0) {
					# feature is on the forward, watson strand
				
					$region = $db->segment( 
						-name      => $name,
						-type      => $type,
						-start     => $feature->start + $relative_start,
						-end       => $feature->start + $relative_stop,
						-strand    => $feature->strand,
						-absolute  => 1,
					);
				}
				elsif ($feature->strand < 0) {
					# feature is on the reverse, crick strand
				
					$region = $db->segment( 
						-name      => $name,
						-type      => $type,
						-start     => $feature->end - $relative_stop,
						-end       => $feature->end - $relative_start,
						-strand    => $feature->strand,
						-absolute  => 1,
					);
				}
			}
			else {
				# feature length is too small to fractionate
				# so we'll simply take the entire region
				$region = $feature->segment;
			}	 

			# confirm region
			unless ($region) { 
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$arg_ref->{'dataset'},
						$arg_ref->{'method'}, 
						$stranded, 
						$log,
			);
		}
	}
		
	elsif (
			defined $start and 
			defined $stop and 
			$relative_pos == 1
	) {
		# if specific start and stop coordinates from the feature middle
		
		FEATURE_DATA_COLLECTION:
		for my $i (1 .. $feat_num) {
			# collecting feature identification 
			my $name = $data_table_ref->[$i][$name_index];
			my $type = $data_table_ref->[$i][$type_index];
			
			# define feature
			my @features = $db->features( 
					-name  => $name,
					-type => $type,
			);
			if (scalar @features > 1) {
				# there should only be one feature found
				# if more, there's redundant or duplicated data in the db
				# warn the user, this should be fixed
				warn " Found more than one feature of '$type => $name' in " .
					"the database!\n Using the first feature only!\n";
			}
			my $feature = shift @features; 
			
			# now re-define the region 
			# relative to the middle of the feature
			# strand does not matter here
			my $middle = ($feature->start + int( $feature->length / 2) );
			
			my $region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $middle + $start,
					-end       => $middle + $stop,
					-strand    => $feature->strand,
					-absolute  => 1,
			);
			
			# confirm region
			unless ($region) { 
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$arg_ref->{'dataset'},
						$arg_ref->{'method'}, 
						$stranded, 
						$log,
			);
		}
	}
	
	
	elsif (
			defined $fstart and 
			defined $fstop and 
			$relative_pos == 1
	) {
		# if fraction values for start and stop are specified
		# define region coordinates relative to middle of the feature
			
		FEATURE_DATA_COLLECTION:
		for my $i (1 .. $feat_num) {
			# collecting feature identification 
			my $name = $data_table_ref->[$i][$name_index];
			my $type = $data_table_ref->[$i][$type_index];
			
			# first define the feature to get length
			my @features = $db->features( 
					-name  => $name,
					-type => $type,
			);
			if (scalar @features > 1) {
				# there should only be one feature found
				# if more, there's redundant or duplicated data in the db
				# warn the user, this should be fixed
				warn " Found more than one feature of '$type => $name' in " .
					"the database!\n Using the first feature only!\n";
			}
			my $feature = shift @features; 
			my $length = $feature->length;
			my $middle = ($feature->start + int( $length / 2) );
			
			# Define the region
			my $region;
			# confirm that feature exceeds our minimum size limit
			if ($length >= $limit) {
				# the feature is long enough to fractionate
				my $relative_start = sprintf "%.0f", 
									$length * $fstart;
				my $relative_stop = sprintf "%.0f", 
									$length * $fstop;
				
				# define relative to the middle
				# strand does not matter
				$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $middle + $relative_start,
					-end       => $middle + $relative_stop,
					-strand    => $feature->strand,
					-absolute  => 1,
				);
			}
			else {
				# feature length is too small to fractionate
				# so we'll simply take the entire region
				$region = $feature->segment;
			}	 

			# confirm region
			unless ($region) { 
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$arg_ref->{'dataset'},
						$arg_ref->{'method'}, 
						$stranded, 
						$log,
			);
		}
	}
	
	
	else {
		# default is to simply take the whole feature region
		
		FEATURE_DATA_COLLECTION:
		for my $i (1 .. $feat_num) {
			# collecting feature identification 
			my $name = $data_table_ref->[$i][$name_index];
			my $type = $data_table_ref->[$i][$type_index];
			
			# first get the feature
			my @features = $db->features( 
					-name  => $name,
					-type => $type,
			);
			if (scalar @features > 1) {
				# there should only be one feature found
				# if more, there's redundant or duplicated data in the db
				# warn the user, this should be fixed
				warn " Found more than one feature of '$type => $name' in " .
					"the database!\n Using the first feature only!\n";
			}
			my $feature = shift @features; 
			
			# then establish the segment
			my $region = $feature->segment(); 
				# this should automatically take care of strand
				# we could probably call the segment directly, but want to 
				# have mechanism in place in case more than one feature was
				# present
			
			# confirm region
			unless ($region) { 
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$arg_ref->{'dataset'},
						$arg_ref->{'method'}, 
						$stranded, 
						$log,
			);
		}
		
	}
	
	# Finish
	return 1;
}





### Retrieve values from a microarray data set for genomic windows

=item get_genome_dataset 

This subroutine will retrieve data set values across the genome. It uses a 
list of genomic windows or intervals as generated by the subroutine 
get_new_genome_list(). It will collect all dataset values within each 
window, combine them with the specified method, and record the single value
in the data list. 

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  db =>       The name of the database or a reference to an 
              established database object. 
  dataset =>  The name(s) of the datasets in the database to be 
              collected. The name(s) should be verified using the 
              subroutine validate_dataset_list() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces.
  data =>     A scalar reference to the data table containing the
              list of features. This should be a reference to the
              key 'data_table' in the data hash described in 
              load_tim_data_file(), not to the entire data hash. 
              Note that the column metadata should be updated 
              separately by the calling program.
  method =>   The method used to combine the dataset values found
              in the defined region. Acceptable values include 
              mean, median, range, stddev, min, max, and enumerate. 
              See _get_segment_score() documentation for more info.
  Optional:
  log =>      Boolean value (1 or 0) indicating whether the dataset 
              values are in log2 space or not. If undefined, the 
              dataset name will be checked for the presence of the 
              phrase 'log2' and set accordingly. This argument is
              critical for accurately mathematically combining 
              dataset values in the region.
  strand =>   Indicate whether the dataset values from a specific 
              strand relative to the feature should be collected. 
              Acceptable values include sense, antisense, no,
              or none. Default is 'none'.
          	  
The subroutine will return a true value (1) if successful.

Examples

	my $db_name = 'cerevisiae';
	my $data = \$data_hash->{'data_table'};
	my $success = get_genome_dataset( {
		'db'      => $db_name,
		'data'    => $data,
		'method'  => 'mean',
		'dataset' => $dataset,
		'log'     => 1,
	} );
	


=cut

sub get_genome_dataset {

	# retrieve passed values
	my $arg_ref = shift; # the passed values as a hash reference
	
	# check required values
	unless ($arg_ref->{'dataset'}) {
		carp " no dataset requested!";
		return;
	}
	unless ($arg_ref->{'method'}) {
		carp " no method defined!";
		return;
	}
	my $data_table_ref = $arg_ref->{'data'};
	unless ($data_table_ref) {
		carp "no feature list defined!";
		return;
	}
	
	
	# define default values as necessary
	my $log = $arg_ref->{'log'};
	unless (defined $log) {
		# we need to know whether we are working with a log2 dataset or not
		# as it will adversely affect the math!
		if ($arg_ref->{'dataset'} =~ /log2/i) {
			# if we're smart we'll encode the log2 status in the dataset name
			# but chances are, we're not that smart
			$log = 1;
		} else {
			# otherwise assume it is non-log
			# unsafe, but what else to do? we'll put the onus on the user
			$log = 0;
		}
	}
	my $stranded = $arg_ref->{'strand'};
	unless (defined $stranded) {
		$stranded = 'none';
	}
	
	# Open a db connection 
	# determine whether we have just a database name or an opened object
	my $db; # the database object to be used
	if (defined $arg_ref->{'db'}) {
		my $db_ref = ref $arg_ref->{'db'};
		if ($db_ref =~ /Bio::DB/) {
			# a db object returns the name of the package
			# this appears to be a bioperl db object
			$db = $arg_ref->{'db'};
		}
		else {
			# the name of a database
			$db = open_db_connection( $arg_ref->{'db'} );
		}
	}
	else {
		carp 'no database name passed!';
		return;
	}
	unless ($db) {
		carp 'no database connected!';
		return;
	}
		
	# generate column name
	push @{ $data_table_ref->[0] }, $arg_ref->{'dataset'};
	
	
	# Automatically identify column indices
	my ($chr_index, $start_index, $stop_index); 
	# we need to identify which columns in the feature table are for 
	# chromosome, start, and stop. these should be columns 0, 1, and 2,
	# respectively, but just in case...
	# and to avoid hard coding index values in case someone changes them
	for (my $i = 0; $i < scalar @{ $data_table_ref->[0] }; $i++ ) {
		if ($data_table_ref->[0][$i] =~ /^chr/i) {
			# we may only be using
			$chr_index = $i;
		}
		elsif ($data_table_ref->[0][$i] =~ /^start$/i) {
			$start_index = $i;
		}
		elsif ($data_table_ref->[0][$i] =~ /^stop$/i) {
			$stop_index = $i;
		}
	}
	unless (
		defined $chr_index and
		defined $start_index and
		defined $stop_index
	) {
		carp 'unable to identify Chromosome, Start, and/or Stop columns';
		return;
	}
	
	my $n = scalar @{$data_table_ref} - 1;
	for my $i (1..$n) {
		
		# define the region
		my $region = $db->segment( 
					-name  => $data_table_ref->[$i][$chr_index],
					-start => $data_table_ref->[$i][$start_index],
					-end   => $data_table_ref->[$i][$stop_index],
		);
		
		# print warning if not found
		unless ($region) { 
			my $error = "Window at table position $i, chromosome " .
					$data_table_ref->[$i][$chr_index] . " position " .
					$data_table_ref->[$i][$start_index] . " not found!";
			carp "$error";
			next;
		}
 		
 		# get the scores for the region
		# pass to internal subroutine to combine dataset values
 		push @{ $data_table_ref->[$i] }, _get_segment_score(
					$region, 
					$arg_ref->{'dataset'},
					$arg_ref->{'method'}, 
					$stranded, 
					$log,
		);

	}
	
	# Finish
	undef %OPENED_WIGFILES;
	return 1;
}


### Get a dataset score for a single region

=item get_chromo_region_score 

This subroutine will retrieve a dataset value for a single specified 
region in the genome. The region is specified with chromosomal coordinates:
chromosome name, start, and stop. It will collect all dataset values within the
window, combine them with the specified method, and return the single value.

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  db =>       The name of the database or a reference to an 
              established database object. 
  dataset =>  The name(s) of the datasets in the database to be 
              collected. The name(s) should be verified using the 
              subroutine validate_dataset_list() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces.
  method =>   The method used to combine the dataset values found
              in the defined region. Acceptable values include 
              'mean', 'median', 'range', 'stddev', 'min', and 'max'. 
              See _get_segment_score() documentation for more info.
  chr =>      The name of the chromosome (reference sequence)
  start =>    The start position of the region on the chromosome
  stop =>     The stop position of the region on the chromosome
  end =>      Alias for stop
  
  Optional:
  log =>      Boolean value (1 or 0) indicating whether the dataset 
              values are in log2 space or not. If undefined, the 
              dataset name will be checked for the presence of the 
              phrase 'log2' and set accordingly. This argument is
              critical for accurately mathematically combining 
              dataset values in the region.
  strand =>   Indicate whether the dataset values from a specific 
              strand relative to the feature should be collected. 
              Acceptable values include sense, antisense, no,
              or none. Default is 'none'.
          	  
The subroutine will return the region score if successful.

Examples

	my $db = open_db_connection('cerevisiae');
	my $score = get_chromo_region_score( {
		'db'      => $db,
		'method'  => 'mean',
		'dataset' => $dataset,
		'chr'     => $chromo,
		'start'   => $startposition,
		'stop'    => $stopposition,
		'log'     => 1,
	} );
	


=cut

sub get_chromo_region_score {
	
	# retrieve passed values
	my $arg_ref = shift; 
	
	# Open a db connection 
	# determine whether we have just a database name or an opened object
	my $db; # the database object to be used
	if (defined $arg_ref->{'db'}) {
		my $db_ref = ref $arg_ref->{'db'};
		if ($db_ref =~ /Bio::DB/) {
			# a db object returns the name of the package
			# this appears to be a bioperl db object
			$db = $arg_ref->{'db'};
		}
		else {
			# the name of a database
			$db = open_db_connection( $arg_ref->{'db'} );
		}
	}
	else {
		carp 'no database name passed!';
		return;
	}
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	# get coordinates
	my $chromo = $arg_ref->{'chr'} || $arg_ref->{'seq'} || $arg_ref->{'seq_id'}
		|| undef;
	my $start = $arg_ref->{'start'} || undef;
	my $stop = $arg_ref->{'stop'} || $arg_ref->{'end'} || undef;
	unless ($chromo and $start and $stop) {
		carp "one or more genomic region coordinates are missing!";
		return;
	};
	
	# define default values as necessary
	my $log = $arg_ref->{'log'};
	unless (defined $log) {
		# we need to know whether we are working with a log2 dataset or not
		# as it will adversely affect the math!
		if ($arg_ref->{'dataset'} =~ /log2/i) {
			# if we're smart we'll encode the log2 status in the dataset name
			# but chances are, we're not that smart
			$log = 1;
		} else {
			# otherwise assume it is non-log
			# unsafe, but what else to do? we'll put the onus on the user
			$log = 0;
		}
	}
	my $stranded = $arg_ref->{'strand'};
	unless (defined $stranded) {
		$stranded = 'none';
	}
	
	# define the chromosomal region segment
	my $region = $db->segment(
			$chromo,  
			$start, 
			$stop, 
	);
	unless ($region) {
		carp "unable to define genomic region";
		return;
	}
		
	# get the scores for the region
	# pass to internal subroutine to combine dataset values
	return _get_segment_score(
				$region, 
				$arg_ref->{'dataset'},
				$arg_ref->{'method'}, 
				$stranded, 
				$log,
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

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  db =>       The name of the database or a reference to an 
              established database object. 
  dataset =>  The name(s) of the datasets in the database to be 
              collected. The name(s) should be verified using the 
              subroutine validate_dataset_list() prior to passing.
              Multiple datasets may be given, joined by "&", with no
              spaces. If collecting attributes of features (e.g 
              "count" or "length"), then the name of the feature 
              (as either the GFF "type" or 'type:method') should
              be passed.
  name =>     The name of the genomic feature.
  type =>     The type of the genomic feature.
  Optional:
  extend =>   Indicate an integer value representing the number of 
              bp the feature's region should be extended on both
              sides.
  start =>    Indicate an integer value representing the start  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with "stop".
  stop =>     Indicate an integer value representing the stop  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with "start".
  position => Indicate the relative position of the feature from 
              which the "start" and "stop" positions are calculated.
              Three values are accepted: "5", which denotes the 
              5' end of the feature, "3" which denotes the 
              3' end of the feature, or "1" which denotes the 
              middle of the feature. This option is only used in 
              conjunction with "start" and "stop" options. The 
              default value is "5".
  strand =>   Indicate whether the dataset values from a specific 
              strand relative to the feature should be collected. 
              Acceptable values include sense, antisense, no,
              or none. Default is "none".
  method =>   Indicate which attribute will be returned. Acceptable 
              values include "score", "count", or "length". The  
              default behavior will be to return the score values.
          	  
The subroutine will return the hash if successful.

Example

	my $db = open_db_connection('cerevisiae');
	my %region_scores = get_region_dataset_hash( {
		'db'      => $db,
		'dataset' => $dataset,
		'name'    => $name,
		'type'    => $type,
	} );
	


=cut

sub get_region_dataset_hash {
	
	# retrieve passed values
	my $arg_ref = shift; 
	
	# Open a db connection 
	# determine whether we have just a database name or an opened object
	my $db; # the database object to be used
	if (defined $arg_ref->{'db'}) {
		my $db_ref = ref $arg_ref->{'db'};
		if ($db_ref =~ /Bio::DB/) {
			# a db object returns the name of the package
			# this appears to be a bioperl db object
			$db = $arg_ref->{'db'};
		}
		else {
			# the name of a database
			$db = open_db_connection( $arg_ref->{'db'} );
		}
	}
	else {
		carp 'no database name passed!';
		return;
	}
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	
	# confirm options and assign defaults
	my $name = $arg_ref->{'name'} || undef;
	my $type = $arg_ref->{'type'} || $arg_ref->{'class'} || undef;
	unless (
		defined $name and 
		defined $type 
	) {
		carp "the feature name and/or type are missing!";
		return;
	};
	my $stranded = $arg_ref->{'strand'} || 'none';
	unless (
		$stranded eq 'none' 
		or $stranded eq 'sense' 
		or $stranded eq 'antisense'
	) {
		carp "unknown strand value '$stranded'!";
		return;
	}
	my $method = $arg_ref->{'method'} || 'score';
	my $relative_pos = $arg_ref->{'position'} || 5;
	
	
	# Define the chromosomal region segment
	my $region;
	my $fstart; # to remember the feature start
	if ($arg_ref->{'extend'}) {
		# if an extension is specified to the feature region
		# first define the feature
		my @features = $db->features( 
				-name  => $name,
				-type  => $type,
		);
		if (scalar @features > 1) {
			# there should only be one feature found
			# if more, there's redundant or duplicated data in the db
			# warn the user, this should be fixed
			warn " Found more than one feature of '$type => $name' in " .
				"the database!\n Using the first feature only!\n";
		}
		my $feature = shift @features; 
		
		# record the feature start
		if ($feature->strand >= 0) {
			# forward strand
			$fstart = $feature->start;
		}
		else {
			# reverse strand
			$fstart = $feature->end;
		}
		
		# get length
		my $length = $feature->length;
		
		# now re-define the region based on the extended coordinates
		$region = $db->segment( 
				-name      => $name,
				-type      => $type,
				-start     => $feature->start - $arg_ref->{'extend'},
				-end       => $feature->end + $arg_ref->{'extend'},
				-strand    => $feature->strand,
				-absolute  => 1,
		);
	} 
		
	elsif (
			defined $arg_ref->{'start'} and 
			defined $arg_ref->{'stop'}
	) {
		# if specific start and stop coordinates are requested
		
		# first define the feature to get endpoints
		my @features = $db->features( 
				-name  => $name,
				-type => $type,
		);
		if (scalar @features > 1) {
			# there should only be one feature found
			# if more, there's redundant or duplicated data in the db
			# warn the user, this should be fixed
			warn " Found more than one feature of '$type => $name' in " .
				"the database!\n Using the first feature only!\n";
		}
		my $feature = shift @features; 
		
		# next define the region relative to the feature start or end
		if ($feature->strand >= 0 and $relative_pos == 3) {
			# feature is on forward, top, watson strand
			# set segment relative to the 3' end
			
			# record feature start
			$fstart = $feature->end;
			
			$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $feature->end + $arg_ref->{'start'},
					-end       => $feature->end + $arg_ref->{'stop'},
					-strand    => $feature->strand,
					-absolute  => 1,
			);
		}
		elsif ($feature->strand < 0 and $relative_pos == 3) {
			# feature is on reverse, bottom, crick strand
			# set segment relative to the 3' end
			
			# record feature start
			$fstart = $feature->start;
			
			$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $feature->start - $arg_ref->{'stop'},
					-end       => $feature->start - $arg_ref->{'start'},
					-strand    => $feature->strand,
					-absolute  => 1,
			);
		}
		elsif ($feature->strand >= 0 and $relative_pos == 5) {
			# feature is on forward, top, watson strand
			# set segment relative to the 5' end
			
			# record feature start
			$fstart = $feature->start;
			
			$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $feature->start + $arg_ref->{'start'},
					-end       => $feature->start + $arg_ref->{'stop'},
					-strand    => $feature->strand,
					-absolute  => 1,
			);
		}
		elsif ($feature->strand < 0 and $relative_pos == 5) {
			# feature is on reverse, bottom, crick strand
			# set segment relative to the 5' end
			
			# record feature start
			$fstart = $feature->end;
			
			$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $feature->end - $arg_ref->{'stop'},
					-end       => $feature->end - $arg_ref->{'start'},
					-strand    => $feature->strand,
					-absolute  => 1,
			);
		}
		elsif ($relative_pos == 1) {
			# feature is on any strand
			# set segment relative to the feature middle
			
			# determine the feature midpoint and record it
			$fstart = $feature->start + sprintf("%.0f", ($feature->length / 2));
			
			$region = $db->segment( 
					-name      => $name,
					-type      => $type,
					-start     => $fstart + $arg_ref->{'start'},
					-end       => $fstart + $arg_ref->{'stop'},
					-strand    => $feature->strand,
					-absolute  => 1,
			);
		}
	}
	
	else {
		# default is to simply take the whole gene region
		# first get the feature
		my @features = $db->features( 
				-name  => $name,
				-type => $type,
		);
		if (scalar @features > 1) {
			# there should only be one feature found
			# if more, there's redundant or duplicated data in the db
			# warn the user, this should be fixed
			warn " Found more than one feature of '$type => $name' in " .
				"the database!\n Using the first feature only!\n";
		}
		my $feature = shift @features; 
		
		# record the feature start
		if ($feature->strand >= 0) {
			# forward strand
			$fstart = $feature->start;
		}
		else {
			# reverse strand
			$fstart = $feature->end;
		}
		
		# then establish the segment
		$region = $feature->segment(); 
			# this should automatically take care of strand
			# we could probably call the segment directly, but want to 
			# have mechanism in place in case more than one feature was
			# present
	}
	
	# check region
	unless ($region) { # print warning if not found
		carp " Region for '$name' not found!\n";
		return;
	}
	
	# retrieve the list of features/dataset points in the region
	my @datasetlist = split /[&,]/, $arg_ref->{'dataset'}; 
		# multiple datasets may be combined into a single search, for example
		# transcriptome data that is on f and r strands. These are given as
		# ampersand or comma delimited lists
	my @datapoints = $region->features(-types => [@datasetlist]); 
	
	# Collect the score/attributes for each feature
	my %datahash; # the hash to put the collected data into
		
	# Walk through the datapoints
	foreach my $datapoint (@datapoints) {
	
		# Check which data to take based on strand
		if (
			$stranded eq 'none' # stranded data not requested
			or $datapoint->strand == 0 # unstranded data
			or ( 
				# sense data
				$datapoint->strand == $region->strand 
				and $stranded eq 'sense'
			) 
			or (
				# antisense data
				$datapoint->strand != $region->strand  
				and $stranded eq 'antisense'
			)
		) {
			# we have acceptable data to collect
			
			# First check whether the data is in the database or wig file
			
			# Wig file
			if ($datapoint->has_tag('wigfile') ) {
				# data is in a wigfile
				# we'll need to open the file and read the data
				
				# this block is pretty identical, but incompatible, with 
				# the internal subroutine _get_wig_data()
				# the reason I can't used that subroutine is because I need 
				# to correlate position with score, and since position is 
				# not directly encoded in the wigfile, I have to regenerate 
				# it using the wigfile step
				
				# first check the method
				unless ($method eq 'score') {
					carp " Only scores may be retrieved from wigfile data!\n";
					return;
				}
				
				# get the parameters for collecting wig data
				my $start = $region->start;
				my $end = $region->end;
				
				my @wigdata = _get_wig_data($datapoint, $start, $end);
				
				# regenerate position information and put the data in the hash
				my @wigfiles = $datapoint->get_tag_values('wigfile');
				my $wig = $OPENED_WIGFILES{ $wigfiles[0] };
				my $step = $wig->step;
				my $i = 0;
				for (my $pos = $start; $pos <= $end; $pos += $step) {
					if (defined $wigdata[$i]) {
						$datahash{ $pos } = $wigdata[$i];
					}
					$i++;
				}
			}
			
			# BigWig File
			elsif ($datapoint->has_tag('bigwigfile') ) {
				# data is in a BigWig file
				# we'll need to open the file and read the data
				
				# first check the method
				unless ($method eq 'score') {
					carp " Only scores may be retrieved from wigfile data!\n";
					return;
				}
				
				# get the data
				%datahash = _get_bigwig_data(
					$datapoint,			# feature
					$region->seq_id,	# chromo
					$region->start,		# start
					$region->end,		# end
					1 					# request position too 
				);
			}
			
			# Database
			else {
				# data is in the database
				# much easier to collect
				
				# determine position to record
				my $position;
				if ($datapoint->start == $datapoint->end) {
					# just one position recorded
					$position = $datapoint->start;
				}
				else {
					# calculate the midpoint
					$position = int( 
						($datapoint->start + $datapoint->end) / 2
					);
				}
				
				if ($method eq 'score') {
					$datahash{$position} = $datapoint->score + 0;
				}
				elsif ($method eq 'count') {
					$datahash{$position} += 1;
				}
				elsif ($method eq 'length') {
					$datahash{$position} = $datapoint->length;
				}
			}
		}
	}
	
	# Convert the coordinates to relative positions
		# previous versions of this function that used Bio::DB::GFF returned 
		# the coordinates as relative positions, e.g. -200..200
		# to maintain this compatibility we will convert the coordinates to 
		# relative positions
		# most downstream applications of this function expect this, and it's
		# a little easier to work with. Just a little bit, though....
	my %relative_datahash;
	if ($region->strand >= 0) {
		# forward strand
		foreach my $position (keys %datahash) {
			$relative_datahash{ $position - $fstart } = $datahash{$position};
		}
	}
	elsif ($region->strand < 0) {
		# reverse strand
		foreach my $position (keys %datahash) {
			$relative_datahash{ $fstart - $position } = $datahash{$position};
		}
	}

	
	# return the collected dataset hash
	return %relative_datahash;
}






################################################################################
###############           Internal subroutines             #####################
################################################################################


### Internal subroutine to convert a feature category into a list of classes



=item _features_to_classes

This internal subroutine provides a conveniant look up and conversion of a 
single-word description of a category of features into a list of actual
feature classes in the database. For example, the word 'gene' may include
all ORFs, snRNAs, snoRNAs, and ncRNAs.

Pass the subroutine the feature category name as a scalar value. Accepted
values and their associated database classes include the following. Additional
values and lists could easily be added if necessary.
  
  gene (ORF, gene, ncRNA, snRNA, snoRNA)
  orf (ORF, gene)
  rna (ncRNA, snRNA, snoRNA)
  trna (tRNA)
  cen (centromere)
  ars (ARS)
  all (ORF, gene, ncRNA, snRNA, snoRNA, tRNA, rRNA, centromere, ARS, 
       Promoter, Terminator, long_terminal_repeat, retrotransposon)
  
Finally, a Bio::DB::SeqFeature::Store type string consisting of the gff feature's 
method:source as may be passed and be used. No attempt at validation will 
be done (if it's not valid, simply no features would be returned from the 
database). Multiple types may be passed, delimited by commas (no spaces). 

The subroutine will return an array consisting of the classes.


=cut

sub _features_to_classes {
	my $feature = shift;
	my @types;
		
	my $alias2types = $TIM_CONFIG->get_block('features');
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
established genomic region. It should transparently work with either 
data stored directly in the database (using the GFF score value) or 
stored in a binary wig file (referenced in the database using the 
parent GFF feature's attribute 'wigfile').

It is passed an array of five specific values, all of which must be 
defined and presented in this order. These values include
  
  [0] The genomic segment as a gff database object, the establishment 
      of which is the responsibility of the calling subroutine.
  [1] The dataset name. Multiple datasets may be included, delimited 
      with an ampersand (&). Multiple datasets are merged into one, 
      unless excluded by strand.
  [2] The method of combining all of the dataset values found in the 
      segment into a single value. Accepted values include
         
         enumerate (count the number of features within the region)
         count (same as enumerate)
         mean
         median
         min
         max
         range (returns 'min-max')
         stddev (returns the standard deviation of a population)
         
  [3] The strandedness of acceptable data. Genomic segments 
      established from an inherently stranded database feature 
      (e.g. ORF) have a specific strandedness. If the dataset strand 
      does not match the specified criteria for strandedness, then it 
      is ignored. If the dataset does not have a specified strand, 
      then it is used regardless of the specified criteria for 
      strandedness. Currently, only transcription data is stranded. 
      Accepted values include
         
         sense (take only values with the same strand)
         antisense (take only values on the opposite strand)
         none (take all values)
         
  [4] The log status of the dataset. Many microarray datasets are in 
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
	my ($region, $dataset, $method, $strandedness, $log) = @_;
	
	my @datasetlist = split /[&,]/, $dataset; 
		# multiple datasets may be combined into a single search, for example
		# transcriptome data that is on f and r strands. These are given as
		# ampersand or comma delimited lists
	
	# retrieve all datapoints in the region
	my @datapoints = $region->features(
					-type => [@datasetlist],
	); 
	
	# Check for the enumeration method
	if ($method eq 'enumerate' or $method eq 'count') {
		# this is simple, we simply count the number of features
		my $count = 0; 
		
		# check for strandedness
		if ($strandedness eq 'sense') {
			# count only those features on the same strand
			my $region_strand = $region->strand; 
			
			# check each datapoint for a match to the region strand
			foreach (@datapoints) {
				my $strand = $_->strand;
				if ($strand == 0) { 
					# there is no strand - great! - just take it
					$count++;
				} 
				elsif ($strand == $region_strand) { 
					# strand matches the region strand
					$count++;
				}
			}
		}
		elsif ($strandedness eq 'antisense') {
			# count the antisense features
				
			my $region_strand = $region->strand; 
		
			# check each datapoint for a match to the region strand
			foreach (@datapoints) {
				my $strand = $_->strand;
				if ($strand == 0) {
					# there is no strand - great! - just take it
					$count++;
				} elsif ($strand != $region_strand) {
					# strand doesn't match the region strand
					# that means it must be the opposite strand, so take it
					$count++;
				}
			}
		} 
		elsif ($strandedness eq 'none' or $strandedness eq 'no') {
			# Strandedness should be ignored for this dataset
			$count = scalar(@datapoints);
		} 
		
		return $count;
	}
			
	# The rest of the subroutine deals with numeric score data
	unless (@datapoints) {return '.'} # make sure we have datapoints
	
	# define
	my @scores; # an array for the dataset values
	my $region_score; # the score for the genomic segment
	
	# Check whether we're dealing with wig data or database data
	
	# Wig Data
	if ( $datapoints[0]->has_tag('wigfile') ) {
		# data is in wig format, or at least the first datapoint is
		# if not all the data is in wig format, we're screwed, 'cus we can't mix
		# normally this should only return one datapoint (maybe two if stranded) 
		# the datapoint feature should run across the entire chromosome
		# but then it's not really a dataPOINT, is it? ;-)
		# This uses the Bio::Graphics::Wiggle module
		
		# need to check for stranded data
		if ($strandedness eq 'sense') {
			# collect the sense data
			
			# collect region info
			my $region_strand = $region->strand; 
			my $start = $region->start;
			my $end = $region->end;
			
			foreach my $datapoint (@datapoints) {
				my $strand = $datapoint->strand;
				if ($strand == 0) { 
					# there is no strand - great! - just take it
					my @s = _get_wig_data($datapoint, $start, $end);
					foreach (@s) {
						if (defined $_) {push @scores, $_}
					}
				} 
				elsif ($strand == $region_strand) { 
					# strand matches the region strand
					my @s = _get_wig_data($datapoint, $start, $end);
					foreach (@s) {
						if (defined $_) {push @scores, $_}
					}
				}
			}
		}
		elsif ($strandedness eq 'antisense') {
			# collect the antisense data
			
			# collect region info
			my $region_strand = $region->strand; 
			my $start = $region->start;
			my $end = $region->end;
			
			foreach my $datapoint (@datapoints) {
				my $strand = $datapoint->strand;
				if ($strand == 0) { 
					# there is no strand - great! - just take it
					my @s = _get_wig_data($datapoint, $start, $end);
					foreach (@s) {
						if (defined $_) {push @scores, $_}
					}
				} 
				elsif ($strand != $region_strand) { 
					# strand doesn't match the region strand
					# that means it must be the opposite strand, so take it
					my @s = _get_wig_data($datapoint, $start, $end);
					foreach (@s) {
						if (defined $_) {push @scores, $_}
					}
				}
			}
		} 
		elsif ($strandedness eq 'none' or $strandedness eq 'no') {
			# strandedness should be ignored for this dataset
			
			# collect region info
			my $start = $region->start;
			my $end = $region->end;
			
			foreach my $datapoint (@datapoints) {
				my @s = _get_wig_data($datapoint, $start, $end);
				foreach (@s) {
					if (defined $_) {push @scores, $_}
				}
			}
		} 
	}
	
	
	# BigWig Data
	elsif ( $datapoints[0]->has_tag('bigwigfile') ) {
		# data is in bigwig format
		# this uses the Bio::DB::BigWig adaptor
		
		# need to check for stranded data
		if ($strandedness eq 'sense') {
			# collect the sense data
			
			# collect region info
			my $region_strand = $region->strand; 
			my $start = $region->start;
			my $end = $region->end;
			my $chromo = $region->seq_id;
			
			foreach my $datapoint (@datapoints) {
				my $strand = $datapoint->strand;
				if ($strand == 0) { 
					# there is no strand - great! - just take it
					push @scores, _get_bigwig_data(
						$datapoint, $chromo, $start, $end, 0);
				} 
				elsif ($strand == $region_strand) { 
					# strand matches the region strand
					push @scores, _get_bigwig_data(
						$datapoint, $chromo, $start, $end, 0);
				}
			}
		}
		elsif ($strandedness eq 'antisense') {
			# collect the antisense data
			
			# collect region info
			my $region_strand = $region->strand; 
			my $start = $region->start;
			my $end = $region->end;
			my $chromo = $region->seq_id;
			
			foreach my $datapoint (@datapoints) {
				my $strand = $datapoint->strand;
				if ($strand == 0) { 
					# there is no strand - great! - just take it
					push @scores, _get_bigwig_data(
						$datapoint, $chromo, $start, $end, 0);
				} 
				elsif ($strand != $region_strand) { 
					# strand doesn't match the region strand
					# that means it must be the opposite strand, so take it
					push @scores, _get_bigwig_data(
						$datapoint, $chromo, $start, $end, 0);
				}
			}
		} 
		elsif ($strandedness eq 'none' or $strandedness eq 'no') {
			# strandedness should be ignored for this dataset
			
			# collect region info
			my $start = $region->start;
			my $end = $region->end;
			my $chromo = $region->seq_id;
			
			foreach my $datapoint (@datapoints) {
				push @scores, _get_bigwig_data(
					$datapoint, $chromo, $start, $end, 0);
			}
		} 
	}
	
	
	# Database Data
	else {
		# Working with data stored in the database
		# this is more straight forward in collection
		
		# Check for stranded data
		if ($strandedness eq 'sense') {
			# take only values that are on the same strand
			
			my $region_strand = $region->strand; 
			# get the absolute (not relative) strand for the region
		
			# check each datapoint for a match to the region strand
			foreach (@datapoints) {
				my $strand = $_->strand;
				if ($strand == 0) { 
					# there is no strand - great! - just take it
					push @scores, $_->score;
				} 
				elsif ($strand == $region_strand) { 
					# strand matches the region strand
					push @scores, $_->score;
				}
			} 
			
		}
		
		elsif ($strandedness eq 'antisense') {
			# take the antisense value
				
			my $region_strand = $region->strand; 
			# get the absolute (not relative) strand for the region
		
			# check each datapoint for a match to the region strand
			foreach (@datapoints) {
				my $strand = $_->strand;
				if ($strand == 0) {
					# there is no strand - great! - just take it
					push @scores, $_->score;
				} elsif ($strand != $region_strand) {
					# strand doesn't match the region strand
					# that means it must be the opposite strand, so take it
					push @scores, $_->score;
				}
			}
		} 
		
		elsif ($strandedness eq 'none' or $strandedness eq 'no') {
			# Strandedness should be ignored for this dataset
			foreach (@datapoints) {
				push @scores, $_->score;
			}
		} 
		
		else {
			# somehow bad strand value got past our checks
			croak " unrecognized strand value '$strandedness'!";
		}
	}
	
	# Determine final region score if we have scores
	if (@scores) {
		# we have dataset values from the region
		
		# deal with log2 values if necessary
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
		else {
			# somehow bad method snuck past our checks
			croak " unrecognized method '$method'!";
		}
	
		# convert back to log2 if necessary
		if ($log) { 
			$region_score = log($region_score) / log(2);
		}
	}
	
	else {
		# we have no dataset scores from this region
		# return a 'null' or non-value
		$region_score = '.';
	}
	
	return $region_score;
}


### Internal subroutine to retrieve wig file data

=item _get_wig_data

This internal subroutine is used to collect the dataset scores from a binary 
wig file that is referenced in the database. A single feature representing the 
dataset is present across each chromosome. The feature references the location 
of the binary file representing the dataset scores in an attribute. The file is 
read using the Bio::Graphics::Wiggle module, and the values extracted from the 
region of interest. To speed up the program and avoid repetitive opening and 
closing of the files, the opened wig file object is stored in a global 
hash in case it is needed again.

The subroutine is passed three arguments in the following order:
    
    1) the dataset object for the chromosome, which contains 
       the 'wigfile' attribute
    2) the start coordinate of the region of interest
    3) the end coordinate of the region of interest

The subroutine returns an array of the defined dataset values found within 
the region of interest.

=cut


sub _get_wig_data {
	
	# get passed arguments
	my ($datapoint, $start, $end) = @_;
	
	# check that we have wiggle support
	unless ($WIGGLE_OK) {
		carp " Wiggle support is not enabled. Can't use Bio::Graphics::Wiggle\n";
		return;
	}
	
	# collect wig data
	my @scores;
	if ($datapoint->has_tag('wigfile') ) {
		
		# get wigfile name
		my @wigfiles = $datapoint->get_tag_values('wigfile');
		my $wigfile = shift @wigfiles;
		
		# check for opened wigfile
		my $wig;
		if (exists $OPENED_WIGFILES{$wigfile} ) {
			# this file is already opened, use it
			$wig = $OPENED_WIGFILES{$wigfile};
		}
		else {
			# this file has not been opened yet, open it
			$wig = Bio::Graphics::Wiggle->new($wigfile,0);
			unless ($wig) {
				carp " unable to open data wigfile '$wigfile'";
				return;
			}
			
			# store the opened object for later use
			$OPENED_WIGFILES{$wigfile} = $wig;
		}
		
		# adjust as necessary to avoid wig errors
		if ($start < $wig->start) {
			# adjust the start position
			$start = $wig->start;
		}
		elsif ($start > $wig->end) {
			# nothing we can do here, no values
			return;
		}
		if ($end > $wig->end) {
			# adjust the end position
			$end = $wig->end;
		}
		elsif ($end < $wig->start) {
			# nothing we can do here, no values
			return;
		}
		
		# collect the wig values
		push @scores,  @{ $wig->values($start => $end) };
	}
	return @scores;
}

### Internal subroutine to retrieve wig file data

=item _get_bigwig_data

This internal subroutine is used to collect the dataset scores from a binary 
BigWig file that is referenced in the database. A single feature representing 
the dataset is present across each chromosome. The feature references the 
location of the binary file representing the dataset scores in an attribute. 
The file is opend and read using the Bio::DB::BigWig module, and the values 
and their positions are extracted from the region of interest. To speed up 
the program and avoid repetitive opening and closing of the files, the opened 
wig file object is stored in a global hash in case it is needed again.

The subroutine is passed four arguments in the following order:
    
    1) the dataset object for the chromosome, which contains 
       the 'bigwigfile' attribute
    2) the chromosome name of interest
    3) the start coordinate of the region of interest
    4) the end coordinate of the region of interest
    5) boolean value indicating whether position data should be included

The subroutine returns an array. If the position data is requested, the 
array contains interleaving position and score values suitable for loading 
into a hash. Otherwise, only score values are included in the array.


=cut

sub _get_bigwig_data {
	
	# get passed arguments
	my ($datapoint, $chromo, $start, $end, $position_req) = @_;
	
	# check that we have bigwig support
	unless ($BIGWIG_OK) {
		carp " BigWig support is not enabled. Can't use Bio::DB::BigWig\n";
		return;
	}
	
	# collect bigwig data
	my @collected_data;
	if ($datapoint->has_tag('bigwigfile') ) {
		
		# get bigwigfile name
		my @wigfiles = $datapoint->get_tag_values('bigwigfile');
		my $wigfile = shift @wigfiles;
		
		# check for opened bigwigfile
		my $wig;
		if (exists $OPENED_WIGFILES{$wigfile} ) {
			# this file is already opened, use it
			$wig = $OPENED_WIGFILES{$wigfile};
		}
		else {
			# this file has not been opened yet, open it
			$wig = Bio::DB::BigWig->new($wigfile);
			unless ($wig) {
				carp " unable to open data BigWig file '$wigfile'";
				return;
			}
			
			# store the opened object for later use
			$OPENED_WIGFILES{$wigfile} = $wig;
		}
		
		# collect the features
		my @features = $wig->get_features_by_location(
			$chromo, $start, $end);
		if ($position_req) {
			# position data is requested to be included
			foreach (@features) {
				push @collected_data, $_->start;
				push @collected_data, $_->score;
			}
		}
		else {
			# we don't need position data
			foreach (@features) {
				push @collected_data, $_->score;
			}
		}
	}
	
	return @collected_data;
}


=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

