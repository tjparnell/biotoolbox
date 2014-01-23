package Bio::ToolBox::db_helper::config;

use strict;
require Exporter;
use File::Spec;
use Carp;
use Config::Simple;
our $VERSION = '1.14';

# variables
our $default;
our $config_path;
our $BTB_CONFIG;



# Load the config file before exporting it
BEGIN {
	
	# the default configuration
	$default = <<DEFAULT
#### BioToolBox Default Config file ####

[default_db]
user                     = nobody
pass                     = hello
adaptor                  = DBI::mysql
dsn_prefix               = dbi:mysql:
chromosome_exclude       = chrMito, chrMT, 2-micron
window                   = 500

[example_remote]
user                     = me
pass                     = secret
dsn                      = dbi:mysql:host=127.0.0.1;port=3306;database=example

[example_sqlite]
adaptor                  = DBI::SQLite
dsn                      = /path/to/example.sqlite

[features]
rna         = ncRNA, snRNA, snoRNA, tRNA
orf         = gene, ORF
repeat      = repeat_region, long_terminal_repeat, transposable_element_gene, LTR_retrotransposon

[exclude_tags]
orf_classification    = Dubious

[applications]

DEFAULT
;
	
	# possible paths
	my $new  = File::Spec->catdir($ENV{HOME}, '.biotoolbox.cfg');
	my $old  = File::Spec->catdir($ENV{HOME}, 'biotoolbox.cfg');
	my $var = $ENV{'BIOTOOLBOX'} || undef;
	
	# check for file in home directory
	if (-s $new) {
		$config_path = $new;
	}
	
	elsif (-s $old) {
		$config_path = $old;

		# upgrade the name
		my $m;
		eval {
			use File::Copy;
			$m = move($old, $new);
		};
		if ($m) {
			$config_path = $new;
			warn "### Updated $old to $new ###\n";
		}
	}
	
	# check for environment variable
	elsif (defined $var and -s $var) {
		$config_path = $var;
	}	

	# finally, when all else fails, write a new file
	else {
		open(FH, '>', $new) or confess 
			"Cannot write biotoolbox configuration file $new!\n$!\n";
		print FH $default;
		close FH;
		$config_path = $new;
	}
	
	# Open the configuration file
	$BTB_CONFIG = Config::Simple->new($config_path);
	confess "could not read biotoolbox configuration file $config_path\n" 
		unless $BTB_CONFIG;
}

# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw($BTB_CONFIG);
our @EXPORT_OK = qw(add_database add_program);


# The true statement
1; 

sub add_database {
	my %args = @_;
	
	# check
	croak "no name provided for new database configuration\n" unless ($args{name});
	croak "no dsn provided for new database configuration\n" unless 
		($args{dsn} || $args{dsn_prefix});
	
	# set name
	my $name = $args{name};
	delete $args{name};
	
	# add configuration
	$BTB_CONFIG->set_block($name, \%args);
	return _rewrite_config();
}

sub add_program {
	my $path = shift;
	unless (-e $path and -x _) {
		carp "$path either does not exist or is not executable\n";
		return;
	}
	
	# add parameter
	my ($vol, $dir, $file) = File::Spec->splitpath($path);
	$BTB_CONFIG->param('applications.' . $file, $path);
	return _rewrite_config();
}

sub _rewrite_config {
	# write new config file
	my $updated;
	if (-w $config_path) {
		$updated = $BTB_CONFIG->write;
	}
	unless ($updated) {
		# attempt to write in users own directory, which takes precedence
		my $file = File::Spec->catdir($ENV{HOME}, '.biotoolbox.cfg');
		$updated = $BTB_CONFIG->write($file);
		if ($updated) {
			$config_path = $file;
			return 1;
		}
		else {
			carp "unable to write updated configuration to $file!\n" . 
			$BTB_CONFIG->error . "\n";
			return;
		}
	}
}

=head1 NAME

Bio::ToolBox::db_helper::config

=head1 DESCRIPTION

This module accesses the biotoolbox configuration file. This file stores 
multiple database connection settings, as well as default behaviors when 
accessing information from the database, such as excluded attribute tags, 
reference sequence GFF type, etc. It also stores the paths to various 
helper applications. 

The default location for the file is in the user's home directory. 
Alternatively, the location of the file may be referenced through an 
Environment setting under the key 'BIOTOOLBOX'.

The file is intended to be edited by the user for their custom installation. 
The file is a simple INI style text file. The format is detailed below.

=head1 FORMAT

The BioToolBox configuration file is a simple INI-style text file. Blocks 
of settings begin with a [header] line, followed by lines of key = value. 
The value may be single text or comma delimited lists. Examples of settings 
include the following.

=over 4

=item default_db

These are default settings that are shared by all databases.
  
  [default_db]
  user                     = nobody
  pass                     = hello
  adaptor                  = DBI::mysql
  dsn_prefix               = dbi:mysql:
  chromosome_exclude       = chrMito, chrMT, 2-micron
  window                   = 500

The user and password keys are for authenticating to a relational 
database. B<WARNING!!!!!!!!!!!!>
For sanity and security sake, please, PLEASE, generate a read-only 
user for relational database access. Do NOT use a privileged account. 
Any password written here is for all to see and is merely a convenience.

The adaptor key specifies the module driver for connecting to the relational 
database containing the Bio::DB::SeqFeature::Store database. Acceptable 
values include DBI::mysql, DBI::Pg, or DBI::SQLite. 

The dsn key defines the string for connecting to the database. For example, 
to connect to a mysql database 'genome' on localhost through a socket
	
	dbi:mysql:genome

to connect to a remote mysql database
	
	dbi:mysql:host=127.0.0.1;port=3306;database=genome

The dsn_prefix key simply drops the database name, allowing it to be used 
with any database name.

See the documentation for Bio::DB::SeqFeature::Store for syntax of 
adaptor and dsn_prefix keys. 

The chromosome_exclude key provides a list of chromosomes to avoid when 
generating either a list of genomic window intervals or genes. For 
example, the mitochondrial chromosome is usually not included when 
performing analyses. 

The window is the size in bp when generating genomic window intervals. It 
is used by the Bio::ToolBox::db_helper::get_new_genome_list() function.

Multiple database sections may be included. Simply name the section after the 
name of the database. Database specific keys may be included, and missing 
keys will default to the default_db values. 

=item Feature Alias Lists

These are aliases for one or more GFF feature types when searching 
for these features in the database.

Specify as either the GFF "type" or "type:source". These represent GFF 
columns 3 and 2:3, respectively.
  
  [features]
  rna         = ncRNA, snRNA, snoRNA, tRNA
  orf         = gene, ORF
  repeat      = repeat_region, long_terminal_repeat, transposable_element_gene, LTR_retrotransposon

=item Exclude Tags

Some features in the database you simply don't want in your list. For 
example, in the cerevisiaie GFF3 annotation, dubious genes are included 
as gene features, but have the GFF "orf_classification" tag value of 
"Dubious". I don't want these features. Ever. These tags are checked 
using the Bio::ToolBox::db_helper::get_new_feature_list() function.

Specify the tag key and the tag value(s) to be excluded
  
  [exclude_tags]
  orf_classification    = Dubious

=item Applications

Some BioToolBox scripts require helper programs. Enter the name of the 
helper program and the full path of its location. Executable programs 
may also be automatically found in the system path.
  
  [applications]
  wigToBigWig      = /usr/local/bin/wigToBigWig
  java             = /usr/bin/java
  Bar2Gr           = /usr/local/USeq/Apps/Bar2Gr

=back

=head1 USAGE

The module exports a single Config::Simple object ($BTB_CONFIG) representing 
the biotoolbox settings in the configuration file. Please refer to the 
documentation for Config::Simple for usage.

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

