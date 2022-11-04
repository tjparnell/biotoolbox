package Bio::ToolBox::db_helper::config;

use warnings;
use strict;
use Carp;
use File::Spec;
use Config::Simple;
require Exporter;

our $VERSION = '1.70';

# variables
my $default;
my $default_path;
my $config_path;
our $BTB_CONFIG;

# Load the config file before exporting it
BEGIN {

	# the default configuration
	$default = <<DEFAULT;
#### BioToolBox Default Config file ####

[default_db]
user                     = nobody
pass                     = hello
adaptor                  = DBI::mysql
dsn_prefix               = dbi:mysql:

[example_remote]
user                     = me
pass                     = secret
dsn                      = dbi:mysql:host=127.0.0.1;port=3306;database=example

[example_sqlite]
adaptor                  = DBI::SQLite
dsn                      = /path/to/example.sqlite

[applications]

DEFAULT

	# possible paths
	$default_path = File::Spec->catdir( $ENV{HOME}, '.biotoolbox.cfg' );
	my $var = $ENV{'BIOTOOLBOX'} || undef;

	# check for file possible paths
	if ( defined $var and -s $var ) {
		$config_path = $var;
	}
	elsif ( -s $default_path ) {
		$config_path = $default_path;
	}

	# Open the configuration file
	if ($config_path) {
		$BTB_CONFIG = Config::Simple->new($config_path);
	}
	else {
		# no path, open empty object
		# this should still work, it just won't return anything useful
		$BTB_CONFIG = Config::Simple->new( syntax => 'ini' );
	}
}

## no critic
## this is never intended to be used directly by end users
## and exporting config mut always be done
our @ISA       = qw(Exporter);
our @EXPORT    = qw($BTB_CONFIG);
our @EXPORT_OK = qw(add_database add_program);
## use critic

sub add_database {
	my %args = @_;

	# check
	croak 'no name provided for new database configuration' unless ( $args{name} );
	croak 'no dsn provided for new database configuration'
		unless ( $args{dsn} || $args{dsn_prefix} );

	# check that we have config file to update
	unless ($config_path) {
		_write_new_config();
	}

	# set name
	my $name = $args{name};
	delete $args{name};

	# add configuration
	$BTB_CONFIG->set_block( $name, \%args );
	return _rewrite_config();
}

sub add_program {
	my $path = shift;
	unless ( -e $path and -x _ ) {
		carp "$path either does not exist or is not executable\n";
		return;
	}

	# check that we have config file to update
	unless ($config_path) {
		_write_new_config();
	}

	# add parameter
	my ( $vol, $dir, $file ) = File::Spec->splitpath($path);
	$BTB_CONFIG->param( 'applications.' . $file, $path );
	return _rewrite_config();
}

sub _write_new_config {
	open( FH, '>', $default_path )
		or croak "Cannot write biotoolbox configuration file $default_path!\n$!\n";
	print FH $default;
	close FH;
	$config_path = $default_path;
	$BTB_CONFIG->read($default_path);
}

sub _rewrite_config {

	# write new config file
	my $updated;
	if ( -w $config_path ) {
		$updated = $BTB_CONFIG->write;
	}
	unless ($updated) {

		# attempt to write in users own directory, which takes precedence
		my $file = File::Spec->catdir( $ENV{HOME}, '.biotoolbox.cfg' );
		$updated = $BTB_CONFIG->write($file);
		if ($updated) {
			$config_path = $file;
			return 1;
		}
		else {
			carp sprintf "unable to write updated configuration to %s!\n%s\n",
				$file, $BTB_CONFIG->error;
			return;
		}
	}
}

1;

__END__

=head1 NAME

Bio::ToolBox::db_helper::config

=head1 DESCRIPTION

This module accesses the biotoolbox configuration file. This file stores 
multiple database connection settings. It also stores the paths to various 
helper applications. 

The default location for the file is in the user's home directory. 
Alternatively, the location of the file may be referenced through an 
Environment setting under the key C<BIOTOOLBOX>.

Versions prior to 1.54 automatically wrote a config file in every user's 
home directory, whether it was needed or wanted or not. With version 
1.54, a config file is B<only> written when necessary by adding a program 
path or database entry.

The file is intended to be edited by the user for their custom installation. 
The file is a simple INI style text file. The format is detailed below.

=head1 FORMAT

The BioToolBox configuration file F<.biotoolbox.cfg> is a simple INI-style 
text file. Blocks of settings begin with a C<[header]> line, followed by lines 
of key = value. The value may be single text or comma delimited lists. Examples of 
settings include the following.

=over 4

=item default_db

These are default settings that are shared by all databases.

  [default_db]
  user                     = nobody
  pass                     = hello
  adaptor                  = DBI::mysql
  dsn_prefix               = dbi:mysql:

The user and password keys are for authenticating to a relational 
database. B<WARNING!>
For sanity and security sake, please, B<PLEASE>, generate a read-only 
user for relational database access. Do NOT use a privileged account. 
Any password written here is for all to see and is merely a convenience.

The adaptor key specifies the module driver for connecting to the relational 
database containing the L<Bio::DB::SeqFeature::Store> database. Acceptable 
values include C<DBI::mysql>, C<DBI::Pg>, or C<DBI::SQLite>. 

The dsn key defines the string for connecting to the database. For example, 
to connect to a mysql database 'genome' on localhost through a socket

	dbi:mysql:genome

to connect to a remote mysql database

	dbi:mysql:host=127.0.0.1;port=3306;database=genome

The dsn_prefix key simply drops the database name, allowing it to be used 
with any database name.

See the documentation for L<Bio::DB::SeqFeature::Store> for syntax of 
C<adaptor> and C<dsn_prefix> keys. 

Multiple database sections may be included. Simply name the section after the 
name of the database. Database specific keys may be included, and missing 
keys will default to the C<default_db> values. 

=item Applications

Some BioToolBox scripts require helper programs. Enter the name of the 
helper program and the full path of its location. Executable programs 
may also be automatically found in the system path.

  [applications]
  wigToBigWig      = /usr/local/bin/wigToBigWig
  bedToBigBed      = /usr/local/bin/bedToBigBed

=back

=head1 USAGE

The module exports a single L<Config::Simple> object (C<$BTB_CONFIG>) representing 
the biotoolbox settings in the configuration file. Please refer to the 
documentation for L<Config::Simple> for usage.

If an existing configuration file is not present, then it will write a new 
default file in the user's home directory. I make the assumption that the 
user has write privileges in their own home directory. It will fail otherwise.

Two subroutines may be optionally exported for assistance in manipulating the 
configuration file.

=over 4

=item add_database

Easily add a new database configuration block to the file. Pass an array of 
key =E<gt> values to be added as a new configuration block. You must include a 
name =E<gt> short name to label the block. This should be unique in the file, 
if you are adding a new database. Otherwise, it will probably clobber the 
pre-existing configuration block (which may be what you want it to do).
Follow the FORMAT examples for details on what to provide. 

It will attempt to rewrite the configuration file, if the user has write 
privileges. If not, then it will attempt to write a new file in the user's 
home root directory. It will return true upon success. 

   my $success = add_database(
      'name'    => 'hg19',
      'adaptor' => 'DBI::SQLite',
      'dsn'     => '/path/to/hg19.sqlite',
   );

=item add_program

Easily add the path to a binary executable for future reference. This is 
a little bit faster to find in here than searching for it through the 
system PATH.

It will attempt to rewrite the configuration file, if the user has write 
privileges. If not, then it will attempt to write a new file in the user's 
home root directory. It will return true upon success. 

   my $success = add_program('/path/to/wigToBigWig');

=back

=head1 SEE ALSO

L<Bio::ToolBox::Data>, L<Bio::ToolBox::db_helper>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

