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

my @sgdf;
my @pombef;
my @dataf;
my @solexaf;
my @otherf;

my $count = 0;
foreach ($db->types) {
	if (/SGD$/) {
		push @sgdf, $_;
	}
	elsif (/geneDB$/) {
		push @pombef, $_;
	}
	elsif (/data$/) {
		push @dataf, $_;
	}
	elsif (/solexa|illumina/i) {
		push @solexaf, $_;
	}
	else {
		push @otherf, $_;
	}
	$count++;
}
print " found $count feature types in database $dbname\n";


if (@sgdf) {
	print "These are " . scalar @sgdf . " SGD feature types in the database:\n";
	foreach (sort {$a cmp $b} @sgdf) { print "  $_\n"}
}
if (@pombef) {
	print "These are " . scalar @pombef . " GeneDB feature types in the database:\n";
	foreach (sort {$a cmp $b} @pombef) { print "  $_\n"}
}
if (@dataf) {
	print "These are " . scalar @dataf . " data feature types in the database:\n";
	foreach (sort {$a cmp $b} @dataf) { print "  $_\n"}
}
if (@solexaf) {
	print "These are " . scalar @solexaf . " Solexa feature types in the database:\n";
	foreach (sort {$a cmp $b} @solexaf) { print "  $_\n"}
}
if (@otherf) {
	print "These are " . scalar @otherf . " remaining feature types in the database:\n";
	foreach (sort {$a cmp $b} @otherf) { print "  $_\n"}
}
print "That's all\n";
