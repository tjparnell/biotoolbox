#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use FindBin qw($Bin);
use Bio::ToolBox::db_helper::config qw(add_database);
use Bio::ToolBox::db_helper qw(open_db_connection);
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

# check for additional requirements 
my $sql;
eval {
	require DBD::SQLite;
	$sql = 1;
};



my $VERSION = '1.15';

print "\n This program will set up an annotation database\n\n";

### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Get command line options and initialize values
my (
	$ucscdb,
	$path,
	$keep,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'db=s'      => \$ucscdb, # the UCSC database shortname
	'path=s'    => \$path, # the optional path for the database
	'keep!'     => \$keep, # keep the annotation files
	'help'      => \$help, # request help
	'version'   => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " BioToolBox script db_setup.pl, version $VERSION\n\n";
	exit;
}



### Check for requirements
unless ($sql) {
	die " Please install Perl module DBD::SQLite to set up a database\n";
}
unless ($ucscdb) {
	$ucscdb = shift @ARGV or
		die " no database name provided! use --help for more information\n";
}
if ($path) {
	$path = File::Spec->canonpath($path);
	unless (-e $path) {
		mkdir $path or die "unable to make database path $path\n$!\n"; 
	}
	unless (-w _) {
		die " $path is not writable!\n";
	}
}
else {
	# determine a path
	if (-e File::Spec->catdir($ENV{HOME}, 'Library')) {
		$path = File::Spec->catdir($ENV{HOME}, 'Library');
	}
	elsif (-e File::Spec->catdir($ENV{HOME}, 'lib')) {
		$path = File::Spec->catdir($ENV{HOME}, 'lib');
	}
	else {
		# make one for us
		$path = File::Spec->catdir($ENV{HOME}, 'lib');
		mkdir $path or die "unable to make database path $path\n$!\n"; 
	}
}



### Get UCSC annotation
print "##### Fetching annotation from UCSC. This may take a while ######\n";
system(File::Spec->catdir($Bin, 'ucsc_table2gff3.pl'), '--db', $ucscdb, '--ftp', 'all', 
	'--gz') == 0 or die "unable to execute ucsc_table2gff3 script!\n";
my @gff = glob("$ucscdb*.gff3.gz");
my @source = glob("$ucscdb*.txt.gz");
unless (@gff) {
	die "unable to find new GFF3 files!\n";
}



### Build database
print "##### Building database. This may take a while ######\n";
my $database = File::Spec->catdir($path, "$ucscdb.sqlite");
my $temp = File::Spec->tmpdir();

# create a database
my $store = Bio::DB::SeqFeature::Store->new(
    -dsn        => $database,
    -adaptor    => 'DBI::SQLite',
    -tmpdir     => $temp,
    -write      => 1,
    -create     => 1,
    -compress   => 0,
) or die " Cannot create a SeqFeature database connection!\n";

# load the database
my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(
    -store              => $store,
    -sf_class           => 'Bio::DB::SeqFeature',
    -verbose            => 1,
    -tmpdir             => $temp,
    -fast               => 1,
    -ignore_seqregion   => 0,
    -index_subfeatures  => 1,
    -noalias_target     => 0,
    -summary_stats      => 0,
)or die " Cannot create a GFF3 loader for the database!\n";

# on signals, give objects a chance to call their DESTROY methods
# borrowed from bp_seqfeature_load.pl
$SIG{TERM} = $SIG{INT} = sub {  undef $loader; undef $store; die "Aborted..."; };
$loader->load(@gff);



### Check database
my $db;
if (-e $database and -s _) {
	$db = open_db_connection($database);
}
if ($db) {
	print "\n##### Created database $database ######\n";
	my $a = add_database(
		'name'      => $ucscdb,
		'dsn'       => $database,
		'adaptor'   => 'DBI::SQLite',
	);
	if ($a) {
		print <<SUCCESS
The database configuration was added to the BioToolBox configuration 
file. You may use the database in any BioToolBox script with the 
option --db $ucscdb.

You can check the database now by running
  print_feature_types.pl $ucscdb
SUCCESS
	;
	}
}
else {
	print "##### Something went wrong! Database could not be opened #####\n";
	unlink $database if -e $database;
}

### Clean up
if ($db and not $keep) {
	unlink @gff;
	unlink @source;
}



__END__

=head1 NAME

db_setup.pl

=head1 SYNOPSIS

db_setup.pl [--options...] <UCSC database>
  
  Options:
  --db <UCSC database>
  --path </path/to/store/database> 
  --keep
  --version
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <UCSC database>

Provide the short UCSC database name for the species and version you want 
to use. See L<http://genome.ucsc.edu/FAQ/FAQreleases.html> for a current 
list of available UCSC genomes. Examples include hg19, mm10, danRer7, 
sacCer3, etc.

=item --path </path/to/store/database>

Specify the optional path to store the SQLite database file. The default 
path is the C<~/lib>.

=item --keep

Keep the downloaded UCSC tables and converted GFF3 files. Default is to 
delete them.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will simplify the task of generating an annotation database. You 
provide the short name of the UCSC database for the species and genome version 
you are interested in, and the script will automatically download gene annotation 
and build a I<Bio::DB::SeqFeature::Store> database for use with BioToolBox 
scripts. 

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  

