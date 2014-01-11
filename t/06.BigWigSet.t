#!/usr/bin/env perl

# Test script for Bio::ToolBox::Data 
# working with BigWigSet data

use strict;
use Test;
use FindBin '$Bin';

BEGIN {
	if (eval {require Bio::DB::BigWigSet; 1}) {
		plan tests => 16;
	}
	else {
		plan skip_all => 'Optional module Bio::DB::BigWigSet not available';
	}
	$ENV{'BIOTOOLBOX'} = "$Bin/Data/biotoolbox.cfg";
}

use lib "$Bin/../lib";
use Bio::ToolBox::Data;


my $dataset = "$Bin/Data/sample3";

### Open a test file
my $Data = Bio::ToolBox::Data->new(file => "$Bin/Data/sample.bed");
ok($Data);

# add a database
$Data->database($dataset);
ok($Data->database, $dataset);
my $db = $Data->open_database;
ok($db);



### Initialize row stream
my $stream = $Data->row_stream;

# First row is YAL047C
my $row = $stream->next_row;
ok($row->name, 'YAL047C');

# try a segment
my $segment = $row->segment;
ok($segment);
ok($segment->start, 54989);

# score count sum
my $score = $row->get_score(
	'db'       => $dataset,
	'dataset'  => 'sample3',
	'value'    => 'count',
	'method'   => 'sum',
);
# print "count sum for ", $row->name, " is $score\n";
ok($score, 433);

# score mean coverage
$score = $row->get_score(
	'db'       => $db,
	'dataset'  => 'sample3',
	'value'    => 'score',
	'method'   => 'mean',
);
# print "mean coverage for ", $row->name, " is $score\n";
ok(sprintf("%.2f", $score), 1.19);




### Move to the next row
$row = $stream->next_row;
ok($row->start, 57029);
ok($row->strand, -1);

$score = $row->get_score(
	'dataset'  => 'sample3',
	'value'    => 'score',
	'method'   => 'median',
	'stranded' => 'all',
);
# print "both strands score median for ", $row->name, " is $score\n";
ok(sprintf("%.2f", $score), 1.69);

# try stranded data collection
$score = $row->get_score(
	'dataset'  => 'sample3',
	'value'    => 'score',
	'method'   => 'median',
	'stranded' => 'sense',
);
# print "sense score median for ", $row->name, " is $score\n";
ok(sprintf("%.2f", $score), 2.74);

$score = $row->get_score(
	'dataset'  => 'sample3',
	'value'    => 'score',
	'method'   => 'median',
	'stranded' => 'antisense',
);
# print "antisense score median for ", $row->name, " is $score\n";
ok(sprintf("%.2f", $score), 0.38);



### Try positioned score index
my %pos2scores = $row->get_position_scores(
	'dataset'  => 'sample3',
	'value'    => 'score',
	'stranded' => 'sense',
);
ok(scalar keys %pos2scores, 44);
# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
ok(sprintf("%.2f", $pos2scores{55}), 4.16);
ok(sprintf("%.2f", $pos2scores{255}), 2.03);

