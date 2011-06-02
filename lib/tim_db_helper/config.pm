package tim_db_helper::config;

# modules
use strict;
require Exporter;
use Carp;
use Config::Simple;
use FindBin qw($Bin);


### Get the Configuration File 
our $TIM_CONFIG;

# check for environment variable
if (exists $ENV{'BIOTOOLBOX'}) {
	 $TIM_CONFIG = Config::Simple->new($ENV{'BIOTOOLBOX'}) or 
		croak Config::Simple->error();
}	

# check for file in home directory
elsif (-e "$ENV{HOME}/biotoolbox.cfg") {
	 $TIM_CONFIG = Config::Simple->new("$ENV{HOME}/biotoolbox.cfg") or
	 	croak Config::Simple->error();
}

# the old environment variable was TIM_DB_HELPER
# check for that just in case
if (exists $ENV{'TIM_DB_HELPER'}) {
	 $TIM_CONFIG = Config::Simple->new($ENV{'TIM_DB_HELPER'}) or 
		croak Config::Simple->error();
}	

# the old file name that I used to use was tim_db_helper.cfg
# check for that just in case
elsif (-e "$ENV{HOME}/tim_db_helper.cfg") {
	 $TIM_CONFIG = Config::Simple->new("$ENV{HOME}/tim_db_helper.cfg") or
	 	croak Config::Simple->error();
}

# finally, when all else fails, just use the supplied default
else {
	$TIM_CONFIG = Config::Simple->new("$Bin/../biotoolbox.cfg") or 
	 	croak Config::Simple->error();
}

# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw($TIM_CONFIG);


# The true statement
1; 

=head1 NAME

tim_db_helper::config

=head1 DESCRIPTION

This module accesses the biotoolbox configuration file. This file stores 
multiple database connection settings, as well as default behaviors when 
accessing information from the database, such as excluded attribute tags, 
reference sequence GFF type, etc. It also stores the paths to various 
helper applications. 

The default location for the file is in the root biotoolbox directory. 
Alternatively, the file may be stored in the root of the user's
home directory. Or the location of the file may be referenced through an 
Environment setting under the key 'BIOTOOLBOX'.

The file is intended to be edited by the user for their custom installation. 
The file is a simple INI style text file. Documentation on the format may 
be found within the file itself.

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


