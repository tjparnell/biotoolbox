#!/usr/bin/perl

# documentation at end of file

use warnings;
use strict;
use English      qw(-no_match_vars);
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use File::Spec;
use FindBin                         qw($Bin);
use Bio::ToolBox::db_helper::config qw(add_database);
use Bio::ToolBox::db_helper         qw(open_db_connection);
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

# check for additional requirements
my $sql = 0;
eval {
	require DBD::SQLite;
	$sql = 1;
};

our $VERSION = '2.00';

print "\n This program will set up an annotation database\n\n";

### Quick help
unless (@ARGV) {

	# when no command line options are present
	# print SYNOPSIS
	pod2usage(
		{
			'-verbose' => 0,
			'-exitval' => 1,
		}
	);
}

### Get command line options and initialize values
my ( $ucscdb, $path, $keep, $verbose, $help, $print_version, );
my @tables;

# Command line options
GetOptions(
	'd|db=s'     => \$ucscdb,           # the UCSC database shortname
	'p|path=s'   => \$path,             # the optional path for the database
	't|table=s'  => \@tables,           # which tables to collect
	'k|keep!'    => \$keep,             # keep the annotation files
	'V|verbose!' => \$verbose,          # show db loading
	'h|help'     => \$help,             # request help
	'v|version'  => \$print_version,    # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {

	# print entire POD
	pod2usage(
		{
			'-verbose' => 2,
			'-exitval' => 1,
		}
	);
}

# Print version
if ($print_version) {
	print " BioToolBox script db_setup.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

### Check for requirements
unless ($sql) {
	die " Please install Perl module DBD::SQLite to set up a database\n";
}
unless ($ucscdb) {
	$ucscdb = shift @ARGV
		or die " no database name provided! use --help for more information\n";
}
if ($path) {
	$path = File::Spec->rel2abs($path);
	unless ( -e $path ) {
		mkdir $path or die "unable to make database path $path $OS_ERROR\n";
	}
	unless ( -w _ ) {
		die " $path is not writable!\n";
	}
}
else {
	# determine a path
	if ( -e File::Spec->catdir( $ENV{HOME}, 'Library' ) ) {
		$path = File::Spec->catdir( $ENV{HOME}, 'Library' );
	}
	elsif ( -e File::Spec->catdir( $ENV{HOME}, 'lib' ) ) {
		$path = File::Spec->catdir( $ENV{HOME}, 'lib' );
	}
	else {
		# make one for us
		$path = File::Spec->catdir( $ENV{HOME}, 'lib' );
		mkdir $path or die "unable to make database path $path $OS_ERROR\n";
	}
}
if (@tables) {
	if ( $tables[0] =~ /,/ ) {
		my $t = shift @tables;
		@tables = split /,/, $t;
	}
}
else {
	@tables = qw(refgene knowngene);
}
my $start_time = time;

### Get UCSC annotation
print "##### Fetching annotation from UCSC. This may take a while ######\n";
system( File::Spec->catdir( $Bin, 'ucsc_table2gff3.pl' ),
	'--db', $ucscdb, '--ftp', join( ',', @tables ), '--gz' ) == 0
	or die "unable to execute ucsc_table2gff3 script!\n";
my @gff    = glob("$ucscdb*.gff3.gz");
my @source = glob("$ucscdb*.txt.gz");
unless (@gff) {
	die "unable to find new GFF3 files!\n";
}

### Build database
print "##### Building database. This may take a while ######\n";
my $database = File::Spec->catdir( $path, "$ucscdb.sqlite" );
my $temp     = File::Spec->tmpdir();

# create a database
my $store = Bio::DB::SeqFeature::Store->new(
	-dsn      => $database,
	-adaptor  => 'DBI::SQLite',
	-tmpdir   => $temp,
	-write    => 1,
	-create   => 1,
	-compress => 0,               # compression seems to be broken, cannot read db
) or die " Cannot create a SeqFeature database connection!\n";

# load the database
my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(
	-store             => $store,
	-sf_class          => 'Bio::DB::SeqFeature',
	-verbose           => $verbose,
	-tmpdir            => $temp,
	-fast              => 1,
	-ignore_seqregion  => 0,
	-index_subfeatures => 1,
	-noalias_target    => 0,
	-summary_stats     => 0,
) or die " Cannot create a GFF3 loader for the database!\n";

# on signals, give objects a chance to call their DESTROY methods
# borrowed from bp_seqfeature_load.pl
local $SIG{TERM} = local $SIG{INT} =
	sub { undef $loader; undef $store; die "Aborted..."; };
$loader->load(@gff);

### Check database
my $db;
if ( -e $database and -s _ ) {
	$db = open_db_connection($database);
}
if ($db) {
	print "\n##### Created database $database ######\n";
	printf " Finished in %.1f minutes\n\n", ( time - $start_time ) / 60;
	my $a = add_database(
		'name'    => $ucscdb,
		'dsn'     => $database,
		'adaptor' => 'DBI::SQLite',
	);
	if ($a) {
		print <<SUCCESS;
The database configuration was added to the BioToolBox configuration 
file. You may use the database in any BioToolBox script with the 
option --db $ucscdb.

You can check the database now by running
  db_types.pl $ucscdb
SUCCESS
	}
}
else {
	print "##### Something went wrong! Database could not be opened #####\n";
	unlink $database if -e $database;
}

### Clean up
if ( $db and not $keep ) {
	unlink @gff;
	unlink @source;
}

__END__

=head1 NAME

db_setup.pl

A program to setup a SeqFeature::Store SQLite database from UCSC data

=head1 SYNOPSIS

db_setup.pl [--options...] <UCSC database>
  
  Options:
  -d --db <UCSC database>
  -p --path </path/to/store/database> 
  -t --table [refGene|knownGene|all]
  -k --keep
  -V --verbose
  -v --version
  -h --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db E<lt>UCSC databaseE<gt>

Provide the short UCSC database name for the species and version you want 
to use. See L<http://genome.ucsc.edu/FAQ/FAQreleases.html> for a current 
list of available UCSC genomes. Examples include hg19, mm10, danRer7, 
sacCer3, etc.

=item --path E<lt>/path/to/store/databaseE<gt>

Specify the optional path to store the SQLite database file. The default 
path is the C<~/lib>.

=item --table [refGene|knownGene|all]

Provide one or more UCSC tables to load into the database. They may be 
specified as comma-delimited list (no spaces) or as separate, repeated 
arguments. The default is refGene and knownGene (if available).

=item --keep

Keep the downloaded UCSC tables and converted GFF3 files. Default is to 
delete them.

=item --verbose

Show realtime database loading progress.

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
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

