package tim_db_helper;

require Exporter;
use strict;
use Carp;
use File::Spec;
use Config::Simple;
use Bio::DB::SeqFeature::Store;
use Statistics::Lite qw(
	sum
	mean
	median
	min
	max
	range
	stddevp
);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(generate_tim_data_structure);


# check for wiggle support
our $WIGGLE_OK = 0;
eval {
	require tim_db_helper::wiggle;
	tim_db_helper::wiggle->import;
};
unless ($@) {
	$WIGGLE_OK = 1;
}; 
$@ = undef;

# check for BigWig support
our $BIGWIG_OK = 0;
eval { 
	require tim_db_helper::bigwig;
	tim_db_helper::bigwig->import;
};
unless ($@) {
	$BIGWIG_OK = 1;
}; 
$@ = undef;

# check for BigBed support
our $BIGBED_OK = 0;
eval { 
	require tim_db_helper::bigbed;
	tim_db_helper::bigbed->import;
};
unless ($@) {
	$BIGBED_OK = 1;
}; 
$@ = undef;

# check for Bam support
our $BAM_OK = 0;
eval { 
	require tim_db_helper::bam;
	tim_db_helper::bam->import;
};
unless ($@) {
	$BAM_OK = 1;
}; 
$@ = undef;


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
	warn "\n#### Using default configuration file:\n   '$Bin/../lib/tim_db_config.cfg'\n####\n";
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
(.wib file) supported by the module C<Bio::Graphics::Wiggle>, and a binary
BigWig file (.bw file) supported by the module C<Bio::DB::BigWig>. The
BigWig file format is much preferred as it maintains spatial resolution of
the original data and does not lose precision by scaling to 8-bit values,
unlike the .wib file format.

While these functions may appear to be simply a rehashing of the methods
and functions in Bio::DB::SeqFeature::Store, they either provide a simpler
function to often used database methodologies or are designed to work
intimately with the tim data format file and data structures (see
C<tim_file_helper.pm>). One key advantage to these functions is the ability
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
			$TIM_CONFIG->param('default_db.pass') || undef;
		
		# check for empty passwords
		# config::simple passes an empty array when nothing was defined
		if (ref $pass eq 'ARRAY' and scalar @$pass == 0) {$pass = undef}
		
		# set up the dsn
		# it can be specifically defined
		my $dsn = $TIM_CONFIG->param($database . '.dsn') || undef;
		unless (defined $dsn) {
			# or dsn can be generated with the dsn_prefix
			$dsn = $TIM_CONFIG->param($database . '.dsn_prefix') || 
				$TIM_CONFIG->param('default_db.dsn_prefix');
			$dsn .= $database;
		}
		
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
	unless (scalar @_ > 0) { 
		carp "no datasets to validate!\n";
		return "no datasets!";
	}
	
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
		carp "no database name passed!\n";
		return "no database";
	}
	unless ($db) {
		carp "no database connected!\n";
		return "no database connection!";
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
	my $new_data = generate_tim_data_structure(
		$arg_ref->{'features'},
		'Name',
		'Type'
	);
	unless ($new_data) {
		carp " cannot generate tim data structure!\n";
		return;
	}
	my $feature_table = $new_data->{'data_table'}; 
	
	# name of the database
	$new_data->{'db'} = $db_name; 
	
	# List of types
	if (scalar @classes > 1) {
		$new_data->{1}->{'include'} = join(",", @classes);
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
			push @{ $feature_table->[0] }, 'Aliases';
			$new_data->{2} = {
					'name'  => 'Aliases',
					'index' => 2,
			};
			$new_data->{'number_columns'} = 3;
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
		if (exists $new_data->{2}) {
			# we only look for Alias info if we have a column for it
			if ($feature->has_tag('Alias')) {
				push @data, join(q( ), $feature->get_tag_values('Alias'));
			}
			else {
				push @data, '.'; # internal null value
			}
		}
		
		# Record information
		push @{$feature_table}, \@data;
		$new_data->{'last_row'} += 1;
	}
	
	
	# print result of search
	print "   Kept " . $new_data->{'last_row'} . " features.\n";
	
	# sort the table
	my @feature_table_sorted;
	my $header = shift @{$feature_table}; # move header
	@feature_table_sorted = sort { 
		# sort first by type, then by name
		( $a->[1] cmp $b->[1] ) and ( $a->[0] cmp $b->[0] )
	} @{$feature_table}; 
	unshift @feature_table_sorted, $header;
	
	# put the feature_table into the data hash
	$new_data->{'data_table'} = \@feature_table_sorted;
	
	# return the new data structure
	return $new_data;
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
              size of the window in basepairs. The default value 
              is defined in tim_db_helper.cfg file.
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
		$win = 
			$TIM_CONFIG->param("$db_name\.window") ||
			$TIM_CONFIG->param('default_db.window');
		print "  Using default window size $win bp\n";
	}
	if ($arg_ref->{'step'}) {
		$step = $arg_ref->{'step'};
	}
	else {
		$step = $win;
	}
	
	
	# Generate data structures
	my $new_data = generate_tim_data_structure(
		'genome',
		'Chromosome',
		'Start',
		'Stop'
	);
	unless ($new_data) {
		carp " cannot generate tim data structure!\n";
		return;
	}
	my $feature_table = $new_data->{'data_table'}; 
	
	# Begin loading basic metadata information
	$new_data->{'db'}      = $db_name; # the db name
	$new_data->{1}{'win'}  = $win; # under the Start metadata
	$new_data->{1}{'step'} = $step;
	
	
	# Collect the chromosomes from the db
		# sort the chromosomes
		# this only works with numbered chromosomes, not Roman numerals
		# additionally skip mitochrondiral chromosome named as chrm or chrmt
		# we're using an elongated pseudo Schwartzian Transform
	my @temp_chromos;
	my $mitochr; # special placeholder for the mitochrondrial chromosome
	my $ref_seq_type = 
		$TIM_CONFIG->param("$db_name\.reference_sequence_type") ||
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
		carp " no '$ref_seq_type' features found! Please check the reference\n" .
			" sequence type in the configuration file and database\n";
		return;
	}
	
	# Check the potential size of the data array
	# I'm afraid of collecting too much data that overwhelms available memory
	my $total_size = 0;
	foreach my $chrobj (@chromosomes) {
		$total_size += $chrobj->length;
	}
	if ( ($total_size / $win) > 1_000_000) {
		warn " data array is predicted to be $total_size lines long!\n" . 
			" this may require considerable amounts of memory\n" .
			" you may wish to consider using 'generate_genomic_bins.pl'\n";
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
			push @{$feature_table}, [ $chr, $start, $end,];
			$new_data->{'last_row'}++;
		}
	}
	print "   Kept " . $new_data->{'last_row'} . " windows.\n";
	
	
	# Return the data structure
	return $new_data;
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
  db       => The name of the database or a reference to an 
              established database object. 
  dataset  => The name of the dataset in the database to be 
              collected. The name should correspond to a feature 
              type in the database, either as type or type:source. 
              The name should be verified using the 
              subroutine validate_dataset_list() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces. Alternatively, specify a data file name. 
              A local file should be prefixed with 'file:', while 
              a remote file should be prefixed with the transfer 
              protocol (ftp: or http:).
  data     => A scalar reference to the data table containing the
              list of features. This should be a reference to the
              key 'data_table' in the data structure described in 
              tim_data_helper.pm, not to the entire data hash. 
              Note that the column metadata should be updated 
              separately by the calling program.
  method   => The method used to combine the dataset values found
              in the defined region. Acceptable values include 
              sum, mean, median, range, stddev, min, max, and count. 
              See _get_segment_score() documentation for more info.
  Optional:
  value    => Specify the type of value to collect. Acceptable 
              values include score, count, or length. The default 
              value type is score. 
  log      => Boolean value (1 or 0) indicating whether the dataset 
              values are in log2 space or not. If undefined, the 
              dataset name will be checked for the presence of the 
              phrase 'log2' and set accordingly. This argument is
              critical for accurately mathematically combining 
              dataset values in the region.
  strand   => Indicate whether the dataset values from a specific 
              strand relative to the feature should be collected. 
              Acceptable values include sense, antisense, no,
              or none. Default is 'none'.
  extend   => Indicate an integer value representing the number of 
              bp the feature's region should be extended on both
              sides.
  start    => Indicate an integer value representing the start  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with 'stop'.
  stop     => Indicate an integer value representing the stop  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with 'start'.
  fstart   => Indicate the fraction of feature length at which the
              region will start. The value should be a decimal 
              number and not a whole percentage (this is not 
              verified!). Include a negative sign to force upstream 
              of the reference point. Must be combined with 'fstop'.
  fstop    => Indicate the fraction of feature length at which the
              region will stop. The value should be a decimal 
              number and not a whole percentage (this is not 
              verified!). Include a negative sign to force upstream 
              of the reference point. Must be combined with 'fstart'.
  limit    => Indicate a minimum size of the feature before "fstart" 
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
  set_strand => For those features that are NOT inherently stranded 
              (strand 0), artificially set the features' strand 
              values. Pass a true value (1) to assign strand values 
              to each value. There must be a column in the passed 
              data table containing the strand values (-1, 0, 1) 
              and whose header name contains "strand". 
  subfeature => When set to a true value, the subfeatures of the 
              feature are identified and values collected for each 
              one. The subfeature values are then combined and 
              recorded. Note that the method of combining values is 
              applied twice: once for each subfeature, and then with 
              all the subfeatures. Exon subfeatures are 
              preferentially used first, followed by CDS and UTR 
              subfeatures. Options extend, start, stop, fstart, and 
              fstop are ignored.
          	  
The subroutine will return a true value (the new column name) if 
successful. It will return nothing and print an error message to 
STDERR if not successful.

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
	
	# maximum expression of middle 50%, from bam file
	$dataset='Steinmetz_polyA.bam';
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
	my $value_type = $arg_ref->{'value'} || 'score';
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
	my $stranded    = $arg_ref->{'strand'} || 'all';
	my $relative_pos = $arg_ref->{'position'} || 5;
	my $limit       = $arg_ref->{'limit'} || 1000;
		# this is an arbitrary default value that seems reasonable for
		# moderate resolution microarrays
		# it's only needed if using the fractional start & stop
	
	# define other variables from the argument hash
	my $extend      = $arg_ref->{'extend'} || undef;
	my $start       = $arg_ref->{'start'} || undef;
	my $stop        = $arg_ref->{'stop'} || undef;
	my $fstart      = $arg_ref->{'fstart'} || undef;
	my $fstop       = $arg_ref->{'fstop'} || undef;
	my $set_strand  = $arg_ref->{'set_strand'} || 0;
	my $get_subfeat = $arg_ref->{'subfeature'} || 0;
	
	
	
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
	my $strand_index;
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
		if ($set_strand and $data_table_ref->[0][$i] =~ /^strand$/i) {
			$strand_index = $i;
		}
	}
	unless (defined $name_index and defined $type_index) {
		carp 'unable to identify Name and/or Type columns';
		return;
	}
	
	# generate column name
	my $column_name;
	if ($arg_ref->{'dataset'} =~ /^file|http|ftp/) {
		# a specified file
		# we just want the file name, split it from the path
		foreach (split /&/, $arg_ref->{'dataset'}) {
			my (undef, undef, $file_name) = File::Spec->splitpath($_);
			if ($column_name) {
				$column_name .= '&' . $file_name;
			}
			else {
				$column_name = $file_name;
			}
		}
	}
	else {
		# a feature type, take as is
		$column_name = $arg_ref->{'dataset'};
	}
	push @{ $data_table_ref->[0] }, $column_name;
	
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
	
	if ($get_subfeat) {
		# subfeatures are requested
		# we will collect the scores for each subfeature first, then combine
		# the subfeature scores
		
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
			elsif (!@features) {
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			my $feature = shift @features; 
			
			# reassign strand value if requested
			if ($set_strand) {
				$feature->strand($data_table_ref->[$i][$strand_index]);
			}
			
			# get strand
			my $fstrand = $feature->strand;
			
			# separate subfeatures by type
			my @exons;
			my @cdss;
			foreach my $subfeat ($feature->get_SeqFeatures) {
				if ($subfeat->primary_tag =~ /exon/i) {
					push @exons, $subfeat;
				}
				elsif ($subfeat->primary_tag =~ /utr|untranslated/i) {
					push @cdss, $subfeat;
				}
				elsif ($subfeat->primary_tag =~ /cds/i) {
					push @cdss, $subfeat;
				}
				elsif ($subfeat->primary_tag =~ /rna/i) {
					# an RNA subfeature, keep going down another level
					foreach my $f ($subfeat->get_SeqFeatures) {
						if ($f->primary_tag =~ /exon/i) {
							push @exons, $f;
						}
						elsif ($f->primary_tag =~ /utr|untranslated/i) {
							push @cdss, $f;
						}
						elsif ($f->primary_tag =~ /cds/i) {
							push @cdss, $f;
						}
					}
				}
			}
			
			# determine which subfeatures to collect
			my @features_to_check;
			if (@exons) {
				@features_to_check = @exons;
			}
			elsif (@cdss) {
				@features_to_check = @cdss;
			}
			else {
				# found neither exons, CDSs, or UTRs
				# possibly because there were no subfeatures
				# in this case we just take the whole thing
				push @features_to_check, $feature;
			}
			
			# collect the subfeature values
			my @subf_values;
			foreach my $subfeat (@features_to_check) {
			
				# establish the segment
				my $region = $subfeat->segment(); 
				
				# get the score for the subfeature region
				# pass to internal subroutine to combine dataset values
				push @subf_values, _get_segment_score(
							$region, 
							$fstrand,
							$arg_ref->{'dataset'},
							$value_type,
							$arg_ref->{'method'}, 
							$stranded, 
							$log,
				);
			}
			
			
			# deal with log2 values if necessary
			if ($log) {
				@subf_values = map {2 ** $_} @subf_values;
			}
			
			# collect the final score
			# we are using subroutines from Statistics::Lite
			my $parent_score;
			if ($arg_ref->{'method'} eq 'median') {
				# take the median value
				$parent_score = median(@subf_values);
			}
			elsif ($arg_ref->{'method'} eq 'mean') {
				# or take the mean value
				$parent_score = mean(@subf_values);
			} 
			elsif ($arg_ref->{'method'} eq 'range') {
				# or take the range value
				# this is 'min-max'
				$parent_score = range(@subf_values);
			}
			elsif ($arg_ref->{'method'} eq 'stddev') {
				# or take the standard deviation value
				# we are using the standard deviation of the population, 
				# since these are the only scores we are considering
				$parent_score = stddevp(@subf_values);
			}
			elsif ($arg_ref->{'method'} eq 'min') {
				# or take the minimum value
				$parent_score = min(@subf_values);
			}
			elsif ($arg_ref->{'method'} eq 'max') {
				# or take the maximum value
				$parent_score = max(@subf_values);
			}
			elsif ($arg_ref->{'method'} eq 'count') {
				# count the number of values
				$parent_score = sum(@subf_values);
			}
			elsif ($arg_ref->{'method'} eq 'sum') {
				# sum the number of values
				$parent_score = sum(@subf_values);
			}
			else {
				# somehow bad method snuck past our checks
				croak " unrecognized method '$arg_ref->{'method'}'!";
			}
		
			# convert back to log2 if necessary
			if ($log) { 
				$parent_score = log($parent_score) / log(2);
			}
				
			# record the final score for the parent feature
			push @{ $data_table_ref->[$i] }, $parent_score;
		}	
		
	}
	
	elsif (defined $extend) {
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
			elsif (!@features) {
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			my $feature = shift @features; 
			
			# reassign strand value if requested
			if ($set_strand) {
				$feature->strand($data_table_ref->[$i][$strand_index]);
			}
			
			# get strand
			my $fstrand = $feature->strand;
			
			# now re-define the region based on the extended coordinates
				# this is regardless of feature orientation
			my $region = $db->segment( 
				$feature->seq_id,
				$feature->start - $extend,
				$feature->end + $extend,
			);
			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$fstrand,
						$arg_ref->{'dataset'},
						$value_type,
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
			elsif (!@features) {
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			my $feature = shift @features; 
			
			# reassign strand value if requested
			if ($set_strand) {
				$feature->strand($data_table_ref->[$i][$strand_index]);
			}
			
			# get strand
			my $fstrand = $feature->strand;
			
			# now re-define the region based on the extended coordinates
			# relative to the 3' end
			my $region;
			if ($feature->strand >= 0) {
				# feature is on the forward, watson strand
				$region = $db->segment( 
					$feature->seq_id,
					$feature->end + $start,
					$feature->end + $stop,
				);
			}
			elsif ($feature->strand < 0) {
				# feature is on the reverse, crick strand
				$region = $db->segment( 
					$feature->seq_id,
					$feature->start - $stop,
					$feature->start - $start,
				);
			}

			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$fstrand,
						$arg_ref->{'dataset'},
						$value_type,
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
			elsif (!@features) {
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			my $feature = shift @features; 
			
			# reassign strand value if requested
			if ($set_strand) {
				$feature->strand($data_table_ref->[$i][$strand_index]);
			}
			
			# get strand
			my $fstrand = $feature->strand;
			
			# now re-define the region based on the extended coordinates
			# relative to the 5' end of the feature
			my $region;
			if ($feature->strand >= 0) {
				# feature is on the forward, watson strand
			
				$region = $db->segment( 
					$feature->seq_id,
					$feature->start + $start,
					$feature->start + $stop,
				);
			}
			elsif ($feature->strand < 0) {
				# feature is on the reverse, crick strand
			
				$region = $db->segment( 
					$feature->seq_id,
					$feature->end - $stop,
					$feature->end - $start,
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
						$fstrand,
						$arg_ref->{'dataset'},
						$value_type,
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
			elsif (!@features) {
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			my $feature = shift @features; 
			my $length = $feature->length;
			
			# reassign strand value if requested
			if ($set_strand) {
				$feature->strand($data_table_ref->[$i][$strand_index]);
			}
			
			# get strand
			my $fstrand = $feature->strand;
			
			# Define the region
			my $region;
			# confirm that feature exceeds our minimum size limit
			if ($length >= $limit) {
				# the feature is long enough to fractionate
				my $relative_start = int( ($length * $fstart) + 0.5);
									
				my $relative_stop = int( ($length * $fstop) + 0.5);
				
				# define relative to the 3' end
				if ($feature->strand >= 0) {
					# feature is on the forward, watson strand
					$region = $db->segment( 
						$feature->seq_id,
						$feature->end + $relative_start,
						$feature->end + $relative_stop,
					);
				}
				elsif ($feature->strand < 0) {
					# feature is on the reverse, crick strand
					$region = $db->segment( 
						$feature->seq_id,
						$feature->start - $relative_stop,
						$feature->start - $relative_start,
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
						$fstrand,
						$arg_ref->{'dataset'},
						$value_type,
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
			elsif (!@features) {
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			my $feature = shift @features; 
			my $length = $feature->length;
			
			# reassign strand value if requested
			if ($set_strand) {
				$feature->strand($data_table_ref->[$i][$strand_index]);
			}
			
			# get strand
			my $fstrand = $feature->strand;
			
			# Define the region
			my $region;
			# confirm that feature exceeds our minimum size limit
			if ($length >= $limit) {
				# the feature is long enough to fractionate
				my $relative_start = int( ($length * $fstart) + 0.5);
				my $relative_stop = int( ($length * $fstop) + 0.5);
				
				# define relative to the 3' end
				if ($feature->strand >= 0) {
					# feature is on the forward, watson strand
				
					$region = $db->segment( 
						$feature->seq_id,
						$feature->start + $relative_start,
						$feature->start + $relative_stop,
					);
				}
				elsif ($feature->strand < 0) {
					# feature is on the reverse, crick strand
				
					$region = $db->segment( 
						$feature->seq_id,
						-start     => $feature->end - $relative_stop,
						-end       => $feature->end - $relative_start,
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
						$fstrand,
						$arg_ref->{'dataset'},
						$value_type,
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
			elsif (!@features) {
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			my $feature = shift @features; 
			
			# reassign strand value if requested
			if ($set_strand) {
				$feature->strand($data_table_ref->[$i][$strand_index]);
			}
			
			# get strand
			my $fstrand = $feature->strand;
			
			# now re-define the region 
			# relative to the middle of the feature
			# strand does not matter here
			my $middle = ($feature->start + int( ($feature->length / 2) + 0.5));
			
			my $region = $db->segment( 
					$feature->seq_id,
					$middle + $start,
					$middle + $stop,
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
						$fstrand,
						$arg_ref->{'dataset'},
						$value_type,
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
			elsif (!@features) {
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			my $feature = shift @features; 
			my $length = $feature->length;
			my $middle = ($feature->start + int( ($length / 2) + 0.5) );
			
			# reassign strand value if requested
			if ($set_strand) {
				$feature->strand($data_table_ref->[$i][$strand_index]);
			}
			
			# get strand
			my $fstrand = $feature->strand;
			
			# Define the region
			my $region;
			# confirm that feature exceeds our minimum size limit
			if ($length >= $limit) {
				# the feature is long enough to fractionate
				my $relative_start = int( ($length * $fstart) + 0.5);
				my $relative_stop = int( ($length * $fstop) + 0.5);
				
				# define relative to the middle
				# strand does not matter
				$region = $db->segment( 
					$feature->seq_id,
					$middle + $relative_start,
					$middle + $relative_stop,
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
						$fstrand,
						$arg_ref->{'dataset'},
						$value_type,
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
			elsif (!@features) {
				carp "feature $name at table position $i not found!";
				next FEATURE_DATA_COLLECTION;
			}
			my $feature = shift @features; 
			
			# reassign strand value if requested
			if ($set_strand) {
				$feature->strand($data_table_ref->[$i][$strand_index]);
			}
			
			# get strand
			my $fstrand = $feature->strand;
			
			# then establish the segment
			my $region = $feature->segment(); 
				# this should automatically take care of strand
				# we could probably call the segment directly, but want to 
				# have mechanism in place in case more than one feature was
				# present
			
			# get the scores for the region
			# pass to internal subroutine to combine dataset values
			push @{ $data_table_ref->[$i] }, _get_segment_score(
						$region, 
						$fstrand,
						$arg_ref->{'dataset'},
						$value_type,
						$arg_ref->{'method'}, 
						$stranded, 
						$log,
			);
		}
		
	}
	
	# Finish
	return $column_name;
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
  dataset =>  The name of the dataset in the database to be 
              collected. The name should correspond to a feature 
              type in the database, either as type or type:source. 
              The name should be verified using the 
              subroutine validate_dataset_list() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces. Alternatively, specify a data file name. 
              A local file should be prefixed with 'file:', while 
              a remote file should be prefixed with the transfer 
              protocol (ftp: or http:).
  data =>     A scalar reference to the data table containing the
              list of features. This should be a reference to the
              key 'data_table' in the data structure described in 
              tim_data_helper.pm, not to the entire data hash. 
              Note that the column metadata should be updated 
              separately by the calling program.
  method =>   The method used to combine the dataset values found
              in the defined region. Acceptable values include 
              sum, mean, median, range, stddev, min, max, and count. 
              See _get_segment_score() documentation for more info.
  Optional:
  value    => Specify the type of value to collect. Acceptable 
              values include score, count, or length. The default 
              value type is score. 
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
          	  
The subroutine will return a true value (the new column name) if 
successful. It will return nothing and print an error message to 
STDERR if not successful.

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
	my $value_type = $arg_ref->{'value'} || 'score';
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
	my $stranded = $arg_ref->{'strand'} || 'all';
	
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
	my $column_name;
	if ($arg_ref->{'dataset'} =~ /^file|http|ftp/) {
		# a specified file
		# we just want the file name, split it from the path
		foreach (split /&/, $arg_ref->{'dataset'}) {
			my (undef, undef, $file_name) = File::Spec->splitpath($_);
			if ($column_name) {
				$column_name .= '&' . $file_name;
			}
			else {
				$column_name = $file_name;
			}
		}
	}
	else {
		# a feature type, take as is
		$column_name = $arg_ref->{'dataset'};
	}
	push @{ $data_table_ref->[0] }, $column_name;
	
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
					0, # no strand
					$arg_ref->{'dataset'},
					$value_type,
					$arg_ref->{'method'}, 
					$stranded, 
					$log,
		);

	}
	
	# Finish
	return $column_name;
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
  dataset =>  The name of the dataset in the database to be 
              collected. The name should correspond to a feature 
              type in the database, either as type or type:source. 
              The name should be verified using the 
              subroutine validate_dataset_list() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces. Alternatively, specify a data file name. 
              A local file should be prefixed with 'file:', while 
              a remote file should be prefixed with the transfer 
              protocol (ftp: or http:).
  dataf =>    Specify a supported data file from which to collect 
              the data, rather than from the database. 
              Supported file formats include bigwig (.bw), 
              bigbed (.bb), or BAM (.bam). 
  method =>   The method used to combine the dataset values found
              in the defined region. Acceptable values include 
              sum, mean, median, range, stddev, min, and max. 
              See _get_segment_score() documentation for more info.
  chr =>      The name of the chromosome (reference sequence)
  start =>    The start position of the region on the chromosome
  stop =>     The stop position of the region on the chromosome
  end =>      Alias for stop
  
  Optional:
  value    => Specify the type of value to collect. Acceptable 
              values include score, count, or length. The default 
              value type is score. 
  log =>      Boolean value (1 or 0) indicating whether the dataset 
              values are in log2 space or not. If undefined, the 
              dataset name will be checked for the presence of the 
              phrase 'log2' and set accordingly. This argument is
              critical for accurately mathematically combining 
              dataset values in the region.
          	  
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
	
	# check the data source
	unless ($arg_ref->{'dataset'}) {
		carp " no dataset requested!";
		return;
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
	my $value_type = $arg_ref->{'value'} || 'score';
	my $log = $arg_ref->{'log'} || undef;
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
				0, # no strand with genomic region
				$arg_ref->{'dataset'},
				$value_type,
				$arg_ref->{'method'}, 
				'all',  # all stranded datasets
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
  db       => The name of the database or a reference to an 
              established database object. 
  dataset  => The name of the dataset in the database to be 
              collected. The name should correspond to a feature 
              type in the database, either as type or type:source. 
              The name should be verified using the 
              subroutine validate_dataset_list() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces. Alternatively, specify a data file name. 
              A local file should be prefixed with 'file:', while 
              a remote file should be prefixed with the transfer 
              protocol (ftp: or http:).
  name     => The name of the genomic feature.
  type     => The type of the genomic feature.
  Optional:
  extend   => Indicate an integer value representing the number of 
              bp the feature's region should be extended on both
              sides.
  start    => Indicate an integer value representing the start  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with "stop".
  stop     => Indicate an integer value representing the stop  
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
  set_strand => For those features or regions that are NOT 
              inherently stranded (strand 0), artificially set the 
              strand. Three values are accepted: -1, 0, 1. This 
              will overwrite any pre-existing value (it will not, 
              however, migrate back to the database).
  strand   => Indicate whether the dataset values from a specific 
              strand relative to the feature should be collected. 
              Acceptable values include sense, antisense, or all.
              Default is all.
  value    => Indicate which attribute will be returned. Acceptable 
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
	
	### Initialize parameters
	
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
	
	# check the data source
	unless ($arg_ref->{'dataset'}) {
		carp " no dataset requested!";
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
	my $stranded = $arg_ref->{'strand'} || 'all';
	my $value_type = $arg_ref->{'value'} || 'score';
	my $relative_pos = $arg_ref->{'position'} || 5;
	my $set_strand = $arg_ref->{'set_strand'} || undef;
	
	
	
	
	
	### Define the chromosomal region segment
	
	my $region;
	my $fstart; # to remember the feature start
	my $fstrand; # to remember the feature strand
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
		
		# change strand if requested
		if (defined $set_strand) {
			$feature->strand($set_strand);
		}
		
		# record the feature start and strand
		$fstart = $feature->strand >= 0 ? $feature->start : $feature->end;
		$fstrand = $feature->strand;
		
		# get length
		my $length = $feature->length;
		
		# now re-define the region based on the extended coordinates
		$region = $db->segment( 
				$feature->seq_id,
				$feature->start - $arg_ref->{'extend'},
				$feature->end + $arg_ref->{'extend'},
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
		
		# change strand if requested
		if (defined $set_strand) {
			$feature->strand($set_strand);
		}
		
		# next define the region relative to the feature start or end
		if ($feature->strand >= 0 and $relative_pos == 3) {
			# feature is on forward, top, watson strand
			# set segment relative to the 3' end
			
			# record feature start
			$fstart = $feature->end;
			$fstrand = $feature->strand;
			
			$region = $db->segment( 
					$feature->seq_id,
					$feature->end + $arg_ref->{'start'},
					$feature->end + $arg_ref->{'stop'},
			);
		}
		elsif ($feature->strand < 0 and $relative_pos == 3) {
			# feature is on reverse, bottom, crick strand
			# set segment relative to the 3' end
			
			# record feature start
			$fstart = $feature->start;
			$fstrand = $feature->strand;
			
			$region = $db->segment( 
					$feature->seq_id,
					$feature->start - $arg_ref->{'stop'},
					$feature->start - $arg_ref->{'start'},
			);
		}
		elsif ($feature->strand >= 0 and $relative_pos == 5) {
			# feature is on forward, top, watson strand
			# set segment relative to the 5' end
			
			# record feature start
			$fstart = $feature->start;
			$fstrand = $feature->strand;
			
			$region = $db->segment( 
					$feature->seq_id,
					$feature->start + $arg_ref->{'start'},
					$feature->start + $arg_ref->{'stop'},
			);
		}
		elsif ($feature->strand < 0 and $relative_pos == 5) {
			# feature is on reverse, bottom, crick strand
			# set segment relative to the 5' end
			
			# record feature start
			$fstart = $feature->end;
			$fstrand = $feature->strand;
			
			$region = $db->segment( 
					$feature->seq_id,
					$feature->end - $arg_ref->{'stop'},
					$feature->end - $arg_ref->{'start'},
			);
		}
		elsif ($relative_pos == 1) {
			# feature can be on any strand
			# set segment relative to the feature middle
			
			# determine the feature midpoint and record it
			$fstart = $feature->start + int(($feature->length / 2) + 0.5);
			$fstrand = $feature->strand;
			
			$region = $db->segment( 
					$feature->seq_id,
					$fstart + $arg_ref->{'start'},
					$fstart + $arg_ref->{'stop'},
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
		
		# change strand if requested
		if (defined $set_strand) {
			$feature->strand($set_strand);
		}
		
		# record the feature start and strand
		$fstart = $feature->strand >= 0 ? $feature->start : $feature->end;
		$fstrand = $feature->strand;
		
		# then establish the segment
		$region = $feature->segment(); 
			# this should automatically take care of strand
			# we could probably call the segment directly, but want to 
			# have mechanism in place in case more than one feature was
			# present
	}
	
	# check region
	unless ($region) { 
		# print warning if not found
		carp " Region for '$name' not found!\n";
		return;
	}
	
	
	
	
	### Data collection
	
	# Collect the score/attributes for each feature
	my %datahash; # the hash to put the collected data into
	
	# retrieve the list of features/dataset points in the region
	my @datasetlist = split /[&,]/, $arg_ref->{'dataset'}; 
		# multiple datasets may be combined into a single search, for example
		# transcriptome data that is on f and r strands. These are given as
		# ampersand or comma delimited lists
	
	# Data source files were provided
	if ($datasetlist[0] =~ /^file|ftp|http/) {
				
		# BigWig file
		if ($datasetlist[0] =~ /\.bw$/) {
		
			# collect if possible
			if ($BIGWIG_OK) {
				# the system has Bio::DB::BigWig loaded
				# get the dataset scores
				%datahash = collect_bigwig_position_scores(
					$region, $fstrand, $stranded, @datasetlist);
			}
			else {
				croak " Bio::DB::BigWig is not available for retrieving BigWig data!\n";
			}
		}
		
		# BigBed file
		elsif ($datasetlist[0] =~ /\.bb$/) {
			
			# collect if possible
			if ($BIGBED_OK) {
				# the system has Bio::DB::BigBed loaded
				# get the dataset scores
				%datahash = collect_bigbed_position_scores(
					$region, $fstrand, $stranded, $value_type, @datasetlist);
			}
			else {
				croak " Bio::DB::BigBed is not available for retrieving BigBed data!\n";
			}
		}
		
		# Bam file
		elsif ($datasetlist[0] =~ /\.bam$/) {
			
			# collect if possible
			if ($BAM_OK) {
				# the system has Bio::DB::Sam loaded
				# get the dataset scores
				%datahash = collect_bam_position_scores(
					$region, $fstrand, $stranded, $value_type, @datasetlist);
			}
			else {
				croak " Bio::DB::Sam is not available for retrieving Bam data!\n";
			}
		}
		
		# Unsupported file
		else {
			croak " Unsupported data file requested!\n";
		}
	}
	
	# working with database features
	else {
		# retrieve the list of features/dataset points in the region
		my @datapoints = $region->features(-types => [@datasetlist]); 
		
		# Check that we have collected datapoints
		unless (@datapoints) {
			# nothing found, return empty handed
			return %datahash;
		}
		
		# Check for wigfiles
		if ($datapoints[0]->has_tag('wigfile') ) {
			# these datapoints point to a binary wig file
			
			# collect if possible
			if ($WIGGLE_OK) {
				# get the dataset scores using tim_db_helper::wiggle
				%datahash = collect_wig_position_scores(
					$region, $fstrand, $stranded, @datapoints);
			}
			else {
				croak " Bio::Graphics::Wiggle is not available for retrieving wig data!\n";
			}
		}
		
		elsif ($datapoints[0]->has_tag('bigwigfile') ) {
			# these datapoints point to a binary bigwig file
			
			# collect if possible
			if ($BIGWIG_OK) {
				# the system has Bio::DB::BigWig loaded
				# get the dataset scores
				%datahash = collect_bigwig_position_scores(
					$region, $fstrand, $stranded, @datapoints);
			}
			else {
				croak " Bio::DB::BigWig is not available for retrieving BigWig data!\n";
			}
		}
		
		elsif ($datapoints[0]->has_tag('bigbedfile') ) {
			# these datapoints point to a binary bigbed file
			
			# collect if possible
			if ($BIGBED_OK) {
				# the system has Bio::DB::BigBed loaded
				# get the dataset scores
				%datahash = collect_bigbed_position_scores(
					$region, $fstrand, $stranded, $value_type, @datapoints);
			}
			else {
				croak " Bio::DB::BigBed is not available for retrieving BigBed data!\n";
			}
		}
		
		elsif ($datapoints[0]->has_tag('bamfile') ) {
			# these datapoints point to a binary bam file
			
			# collect if possible
			if ($BAM_OK) {
				# the system has Bio::DB::Sam loaded
				# get the dataset scores
				%datahash = collect_bam_position_scores(
					$region, $fstrand, $stranded, $value_type, @datapoints);
			}
			else {
				croak " Bio::DB::Sam is not available for retrieving Bam data!\n";
			}
		}
		
		else {
			# these datapoints must otherwise represent features in the database
			# collect score or attributes directly
			
			# Walk through the datapoints
			foreach my $datapoint (@datapoints) {
			
				# Check which data to take based on strand
				if (
					$stranded eq 'all' # all data is requested
					or $datapoint->strand == 0 # unstranded data
					or ( 
						# sense data
						$fstrand == $datapoint->strand 
						and $stranded eq 'sense'
					) 
					or (
						# antisense data
						$fstrand != $datapoint->strand  
						and $stranded eq 'antisense'
					)
				) {
					# we have acceptable data to collect
				
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
					
					if ($value_type eq 'score') {
						# perform addition to force the score to be a scalar value
						push @{ $datahash{$position} }, $datapoint->score + 0;
					}
					elsif ($value_type eq 'count') {
						$datahash{$position} += 1;
					}
					elsif ($value_type eq 'length') {
						push @{ $datahash{$position} }, $datapoint->length;
					}
				}
			}
			
			# post-process the collected values 
			# combine multiple values recorded at the same position
			if ($value_type eq 'score' or $value_type eq 'length') {
				# each 'value' is an array of one or more scores or lengths 
				# from the datapoints collected above
				# we will take the simple mean
				foreach my $position (keys %datahash) {
					$datahash{$position} = mean( @{$datahash{$position}} );
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
	my %relative_datahash;
	if ($fstrand >= 0) {
		# forward strand
		foreach my $position (keys %datahash) {
			$relative_datahash{ $position - $fstart } = $datahash{$position};
		}
	}
	elsif ($fstrand < 0) {
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

Pass the subroutine the feature category name as a scalar value. The 
actual list of feature types will be collected and returned as an array. 
Multiple values may be passed as a comma-delimited string (no spaces).

The aliases and feature lists are specified in the tim_db_helper 
configuration file, tim_db_helper.cfg. Additional lists and aliases 
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
established genomic region. It works with a variety of data sources. 

First, the data may be stored directly in the Bio::DB::SeqFeature::Store 
(using the original GFF score value). Second, the feature may reference 
a data file as an attribute (e.g., wigfile=/path/to/file.wib). Finally, 
the name(s) of a data file may be passed from which to collect the data.

The subroutine is passed an array of five specific values, all of which 
must be defined and presented in this order. These values include
  
  [0] The genomic segment as a gff database object, the establishment 
      of which is the responsibility of the calling subroutine.
  [1] The strand of the region (or original feature), as the strand is 
      lost when generating the segment. Acceptable values are 
      -1, 0, or 1.
  [2] The dataset name for filename. Multiple datasets may be included, 
      delimited with an ampersand (&). Multiple datasets are merged into 
      one, unless excluded by strand. Local data source files should be 
      prepended with 'file:', while remote data source files should be 
      prepended with the transfer protocol (http: or ftp:).
  [3] The data type to be collecting. In most cases, the score value 
      is used, but others may be collected. Accepted values include
         
         score
         count
         length
         
  [4] The method of combining all of the dataset values found in the 
      segment into a single value. Accepted values include
         
         sum
         mean
         median
         min
         max
         range (returns difference between max and min)
         stddev (returns the standard deviation of a population)
         
  [5] The strandedness of acceptable data. Genomic segments 
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
         
  [6] The log status of the dataset. Many microarray datasets are in 
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
		$region, 
		$region_strand, 
		$dataset, 
		$value_type,
		$method, 
		$strandedness, 
		$log
	) = @_;
	
	# define
	my @scores; # an array for the dataset values
	my $region_score; # the score for the genomic segment
	
	my @datasetlist = split /[&,]/, $dataset; 
		# multiple datasets may be combined into a single search, for example
		# transcriptome data that is on f and r strands. These are given as
		# ampersand or comma delimited lists
	
	# check for data source files
	if ($datasetlist[0] =~ /^(?:file|http|ftp):/) {
		
		# collect the data according to file type
		
		# BigWig Data file
		if ($datasetlist[0] =~ /\.bw$/i) {
			# file is in bigwig format
			# this uses the Bio::DB::BigWig adaptor
			
			# check that we have bigwig support
			unless ($BIGWIG_OK) {
				croak " BigWig support is not enabled. Can't use Bio::DB::BigWig\n";
			}
			
			# collect the scores using tim_db_helper::bigwig
			@scores = collect_bigwig_scores(
				$region, $region_strand, $strandedness, @datasetlist);
		}		
		
		# BigBed Data file
		elsif ($datasetlist[0] =~ /\.bb$/i) {
			# data is in bigbed format
			# this uses the Bio::DB::BigBed adaptor
			
			# check that we have bigbed support
			unless ($BIGBED_OK) {
				croak " BigBed support is not enabled. Can't use Bio::DB::BigBed\n";
			}
			
			# collect the scores using tim_db_helper::bigbed
			@scores = collect_bigbed_scores(
				$region, 
				$region_strand, 
				$strandedness, 
				$value_type, 
				@datasetlist
			);
		}
		
		# BAM data file
		elsif ($datasetlist[0] =~ /\.bam$/i) {
			# data is in bam format
			# this uses the Bio::DB::Sam adaptor
			
			# check that we have bigbed support
			unless ($BAM_OK) {
				croak " Bam support is not enabled. Can't use Bio::DB::Sam\n";
			}
			
			# collect the scores using tim_db_helper::bam
			@scores = collect_bam_scores(
				$region, 
				$region_strand, 
				$strandedness, 
				$value_type, 
				@datasetlist
			);
		}
		
		# Unsupported Data file
		else {
			croak " Unsupported file type for file '$datasetlist[0]!\n";
		}
		
	}
	
	# otherwise all other data is from the database
	else {
	
	
		# retrieve all datapoints in the region
		my @datapoints = $region->features(
						-type => [@datasetlist],
		); 
		
		# make sure we have datapoints
		unless (@datapoints) {return '.'} 
		
		# Check whether we're dealing with wig data or database data
		
		# Wig Data
		if ( $datapoints[0]->has_tag('wigfile') ) {
			# data is in wig format, or at least the first datapoint is
			
			# check that we have wiggle support
			unless ($WIGGLE_OK) {
				croak " Wiggle support is not enabled. Can't use Bio::Graphics::Wiggle\n";
				return;
			}
			
			# collect the scores using tim_db_helper::wiggle
			@scores = collect_wig_scores(
				$region, $region_strand, $strandedness, @datapoints);
		}
		
		
		# BigWig Data
		elsif ( $datapoints[0]->has_tag('bigwigfile') ) {
			# data is in bigwig format
			# this uses the Bio::DB::BigWig adaptor
			
			# check that we have bigwig support
			unless ($BIGWIG_OK) {
				croak " BigWig support is not enabled. Can't use Bio::DB::BigWig\n";
			}
			
			# collect the scores using tim_db_helper::bigwig
			@scores = collect_bigwig_scores(
				$region, $region_strand, $strandedness, @datapoints);
		}
		
		
		# BigBed Data
		elsif ( $datapoints[0]->has_tag('bigbedfile') ) {
			# data is in bigbed format
			# this uses the Bio::DB::BigBed adaptor
			
			# check that we have bigbed support
			unless ($BIGBED_OK) {
				croak " BigBed support is not enabled. Can't use Bio::DB::BigBed\n";
			}
			
			# collect the scores using tim_db_helper::bigbed
			@scores = collect_bigbed_scores(
				$region, 
				$region_strand, 
				$strandedness, 
				$value_type, 
				@datapoints
			);
		}
		
		# Bam Data
		elsif ( $datapoints[0]->has_tag('bamfile') ) {
			# data is in bam format
			# this uses the Bio::DB::Sam adaptor
			
			# check that we have bam support
			unless ($BAM_OK) {
				croak " Bam support is not enabled. Can't use Bio::DB::Sam\n";
			}
			
			# collect the scores using tim_db_helper::bam
			@scores = collect_bam_scores(
				$region, 
				$region_strand, 
				$strandedness, 
				$value_type, 
				@datapoints
			);
		}
		
		
		# Database Data
		else {
			# Working with data stored in the database
			# this is more straight forward in collection
			
			# we'll store the score value for each found feature
			# features that don't have score return 'undef', but that'll still be 
			# added to the array
			
			# Check for stranded data
			if ($strandedness eq 'sense') {
				# take only values that are on the same strand
				
				# check each datapoint for a match to the region strand
				foreach my $f (@datapoints) {
					if (
						$f->strand == 0 or 
						$f->strand == $region_strand
					) { 
						# there is no strand - great! - just take it
						# or strand matches the region strand
						if ($value_type eq 'score') {
							push @scores, $f->score;
						}
						elsif ($value_type eq 'count') {
							push @scores, 1;
						}
						elsif ($value_type eq 'length') {
							push @scores, $f->length;
						}
					} 
				} 
				
			}
			
			elsif ($strandedness eq 'antisense') {
				# take the antisense value
					
				# check each datapoint for a match to the region strand
				foreach my $f (@datapoints) {
					my $strand = $_->strand;
					if (
						$f->strand == 0 or 
						$f->strand != $region_strand
					) { 
						# there is no strand - great! - just take it
						# or strand doesn't match the region strand
						# that means it must be the opposite strand, so take it
						if ($value_type eq 'score') {
							push @scores, $f->score;
						}
						elsif ($value_type eq 'count') {
							push @scores, 1;
						}
						elsif ($value_type eq 'length') {
							push @scores, $f->length;
						}
					} 
				}
			} 
			
			elsif ($strandedness eq 'all') {
				# Strandedness should be ignored for this dataset
				foreach my $f (@datapoints) {
					if ($value_type eq 'score') {
						push @scores, $f->score;
					}
					elsif ($value_type eq 'count') {
						push @scores, 1;
					}
					elsif ($value_type eq 'length') {
						push @scores, $f->length;
					}
				}
			} 
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
		elsif ($method eq 'count') {
			# count the number of values
			$region_score = scalar(@scores);
		}
		elsif ($method eq 'sum') {
			# sum the number of values
			$region_score = sum(@scores);
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
		
		if ($method eq 'count' or $method eq 'sum') {
			# return zero
			$region_score = 0;
		}
		else {
			# return a 'null' or non-value
			$region_score = '.';
		}
	}
	
	return $region_score;
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

