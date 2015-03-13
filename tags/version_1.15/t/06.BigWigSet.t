#!/usr/bin/env perl

# Test script for Bio::ToolBox::Data 
# working with BigWigSet data

use strict;
use Test::More;
use FindBin '$Bin';

BEGIN {
	if (eval {require Bio::DB::BigWigSet; 1}) {
		plan tests => 18;
	}
	else {
		plan skip_all => 'Optional module Bio::DB::BigWigSet not available';
	}
	$ENV{'BIOTOOLBOX'} = "$Bin/Data/biotoolbox.cfg";
}

use lib "$Bin/../lib";
require_ok 'Bio::ToolBox::Data' or 
	BAIL_OUT "Cannot load Bio::ToolBox::Data";


my $dataset = "$Bin/Data/sample3";

### Open a test file
my $Data = Bio::ToolBox::Data->new(file => "$Bin/Data/sample.bed");
isa_ok($Data, 'Bio::ToolBox::Data', 'BED Data');

# add a database
$Data->database($dataset);
is($Data->database, $dataset, 'get database');
my $db = $Data->open_database;
isa_ok($db, 'Bio::DB::BigWigSet', 'connected database');



### Initialize row stream
my $stream = $Data->row_stream;
isa_ok($stream, 'Bio::ToolBox::Data::Iterator', 'row stream iterator');

# First row is YAL047C
my $row = $stream->next_row;
is($row->name, 'YAL047C', 'row name');

# try a segment
my $segment = $row->segment;
isa_ok($segment, 'Bio::DB::BigWigSet::Segment', 'row segment');
is($segment->start, 54989, 'segment start');

# score count sum
my $score = $row->get_score(
	'db'       => $dataset,
	'dataset'  => 'sample3',
	'value'    => 'count',
	'method'   => 'sum',
);
# print "count sum for ", $row->name, " is $score\n";
is($score, 433, 'row sum of count');

# score mean coverage
$score = $row->get_score(
	'db'       => $db,
	'dataset'  => 'sample3',
	'value'    => 'score',
	'method'   => 'mean',
);
# print "mean coverage for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), 1.19, 'row mean score');



### Move to the next row
$row = $stream->next_row;
is($row->start, 57029, 'row start position');
is($row->strand, -1, 'row strand');

$score = $row->get_score(
	'dataset'  => 'sample3',
	'value'    => 'score',
	'method'   => 'median',
	'stranded' => 'all',
);
# print "both strands score median for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), 1.69, 'row median score');

# try stranded data collection
$score = $row->get_score(
	'dataset'  => 'sample3',
	'value'    => 'score',
	'method'   => 'median',
	'stranded' => 'sense',
);
# print "sense score median for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), 2.74, 'row sense median score');

$score = $row->get_score(
	'dataset'  => 'sample3',
	'value'    => 'score',
	'method'   => 'median',
	'stranded' => 'antisense',
);
# print "antisense score median for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), 0.38, 'row antisense median score');



### Try positioned score index
my %pos2scores = $row->get_position_scores(
	'dataset'  => 'sample3',
	'value'    => 'score',
	'stranded' => 'sense',
);
is(scalar keys %pos2scores, 44, 'number of positioned scores');
# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is(sprintf("%.2f", $pos2scores{55}), 4.16, 'score at position 55');
is(sprintf("%.2f", $pos2scores{255}), 2.03, 'score at position 255');

