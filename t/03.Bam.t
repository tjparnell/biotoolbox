#!/usr/bin/env perl

# Test script for Bio::ToolBox::Data 
# working with Bam data

use strict;
use Test;
use FindBin '$Bin';

BEGIN {
	if (eval {require Bio::DB::Sam; 1}) {
		plan tests => 22;
	}
	else {
		plan skip_all => 'Optional module Bio::DB::Sam not available';
	}
	$ENV{'BIOTOOLBOX'} = "$Bin/Data/biotoolbox.cfg";
}

use lib "$Bin/../lib";
use Bio::ToolBox::Data;


my $dataset = "$Bin/Data/sample1.bam";

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

# read count sum
my $score = $row->get_score(
	'db'       => $dataset,
	'dataset'  => $dataset,
	'value'    => 'count',
	'method'   => 'sum',
);
# print "count sum for ", $row->name, " is $score\n";
ok($score, 453);

# mean coverage
$score = $row->get_score(
	'db'       => $db,
	'dataset'  => $dataset,
	'value'    => 'score',
	'method'   => 'mean',
);
# print "mean coverage for ", $row->name, " is $score\n";
ok(sprintf("%.2f", $score), 16.33);



### Move to the next row
$row = $stream->next_row;
ok($row->start, 57029);
ok($row->strand, -1);

# try stranded data collection
$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'count',
	'method'   => 'sum',
	'stranded' => 'all',
);
# print "all read count sum for ", $row->name, " is $score\n";
ok($score, 183);

$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'count',
	'method'   => 'sum',
	'stranded' => 'sense',
);
# print "sense read count sum for ", $row->name, " is $score\n";
ok($score, 86);

$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'count',
	'method'   => 'sum',
	'stranded' => 'antisense',
);
# print "antisense read count sum for ", $row->name, " is $score\n";
ok($score, 97);

$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'score',
	'method'   => 'mean',
	'stranded' => 'sense',
);
# print "sense mean coverage for ", $row->name, " is $score\n";
ok(sprintf("%.2f", $score), 29.38);

$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'score',
	'method'   => 'mean',
	'stranded' => 'antisense',
);
# print "antisense mean coverage for ", $row->name, " is $score\n";
ok(sprintf("%.2f", $score), 29.38);



### Move to third row
# test row positioned score using bam file
$row = $stream->next_row;
ok($row->name, 'YAL044W-A');

my %pos2scores = $row->get_position_scores(
	'dataset'  => $dataset,
	'value'    => 'count',
);
ok(scalar keys %pos2scores, 110);
# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
ok($pos2scores{1}, 1);
ok($pos2scores{20}, 2);

%pos2scores = $row->get_position_scores(
	'dataset'  => $dataset,
	'value'    => 'count',
	'absolute' => 1,
	'stranded' => 'antisense',
);
# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
ok(scalar keys %pos2scores, 63);
ok($pos2scores{57556}, 2);
ok($pos2scores{57840}, 1);





