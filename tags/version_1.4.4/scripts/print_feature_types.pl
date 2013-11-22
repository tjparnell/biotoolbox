#!/usr/bin/perl

# a quick and dirty program to print out the feature types in the current database

use strict;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_db_helper qw(
	open_db_connection
);

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
	$dbname,
	$help
);

# Command line options
GetOptions( 
	'db=s'      => \$dbname, # the database name
	'help'      => \$help # request help
);

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}




# Check for database
unless ($dbname) {
	$dbname = shift @ARGV or 
		die " Must provide a database name!\n";
}



# Open database
my $db = open_db_connection($dbname) or die " can't open db\n";

# Get the database types
my $count = 0;
my %source2type;
foreach ($db->types) {
	
	# the type is essentially method:source
	# get individual values
	my $source = $_->source;
	my $type = $_->method;
	
	# store the type in an array under the source
	if (exists $source2type{$source}) {
		push @{ $source2type{$source} }, $type;
	}
	else {
		$source2type{$source} = [ ($type) ];
	}
	$count++;
}
print " Found $count feature types in database '$dbname'\n";


# Print the database types by source type
foreach my $source (sort {$a cmp $b} keys %source2type) {
	print "  There are ", scalar @{$source2type{$source}}, " feature types ", 
		"with source '$source'\n";
	foreach (sort {$a cmp $b} @{$source2type{$source}} ) {
		print "     $_\n";
	}
}


print "That's all\n";



__END__

=head1 NAME

print_feature_types.pl

=head1 SYNOPSIS

print_feature_types.pl <database>
  
  --db <database>
  --help
  

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <database>

Specify the name of the Bio::DB::SeqFeature::Store database.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will print a list of all of the known feature types present 
in a Bio::DB::SeqFeature::Store database. The types are organized into 
groups by their source tag.

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

