#!/usr/bin/perl

# a quick and dirty program to print out the feature types in the current database

use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_db_helper qw(
	open_db_connection
);

my $dbname = shift @ARGV or die "enter the name of the database following program name\n";
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
