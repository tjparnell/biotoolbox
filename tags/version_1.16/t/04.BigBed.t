#!/usr/bin/env perl

# Test script for Bio::ToolBox::Data 
# working with BigBed data

use strict;
use Test::More;
use FindBin '$Bin';

BEGIN {
	if (eval {require Bio::DB::BigBed; 1}) {
		plan tests => 24;
	}
	else {
		plan skip_all => 'Optional module Bio::DB::BigBed not available';
	}
	$ENV{'BIOTOOLBOX'} = "$Bin/Data/biotoolbox.cfg";
}

use lib "$Bin/../lib";
require_ok 'Bio::ToolBox::Data' or 
	BAIL_OUT "Cannot load Bio::ToolBox::Data";


my $dataset = "$Bin/Data/sample1.bb";

### Open a test file
my $Data = Bio::ToolBox::Data->new(file => "$Bin/Data/sample.bed");
isa_ok($Data, 'Bio::ToolBox::Data', 'BED Data');

# add a database
$Data->database($dataset);
is($Data->database, $dataset, 'get database');
my $db = $Data->open_database;
isa_ok($db, 'Bio::DB::BigBed', 'connected database');



### Initialize row stream
my $stream = $Data->row_stream;
isa_ok($stream, 'Bio::ToolBox::Data::Iterator', 'row stream iterator');

# First row is YAL047C
my $row = $stream->next_row;
is($row->name, 'YAL047C', 'row name');

# try a segment
my $segment = $row->segment;
isa_ok($segment, 'Bio::DB::BigFile::Segment', 'row segment');
is($segment->start, 54989, 'segment start');

# read count sum
my $score = $row->get_score(
	'db'       => $dataset,
	'dataset'  => $dataset,
	'value'    => 'count',
	'method'   => 'sum',
);
# print "count sum for ", $row->name, " is $score\n";
is($score, 453, 'row sum of read count score');

# mean coverage
$score = $row->get_score(
	'db'       => $db,
	'dataset'  => $dataset,
	'value'    => 'score',
	'method'   => 'mean',
);
# print "mean coverage for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), 143.81, 'row mean score');



### Move to the next row
$row = $stream->next_row;
is($row->start, 57029, 'row start position');
is($row->strand, -1, 'row strand');

# try stranded data collection
$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'count',
	'method'   => 'sum',
	'stranded' => 'all',
);
# print "all read count sum for ", $row->name, " is $score\n";
is($score, 183, 'row sum of count score for all strands');

$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'count',
	'method'   => 'sum',
	'stranded' => 'sense',
);
# print "sense read count sum for ", $row->name, " is $score\n";
is($score, 86, 'row sum of count score for sense strand');

$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'count',
	'method'   => 'sum',
	'stranded' => 'antisense',
);
# print "antisense read count sum for ", $row->name, " is $score\n";
is($score, 97, 'row sum of count score for antisense strand');

$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'score',
	'method'   => 'mean',
	'stranded' => 'sense',
);
# print "sense mean coverage for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), 146.88, 'row mean score for sense strand');

$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'score',
	'method'   => 'mean',
	'stranded' => 'antisense',
);
# print "antisense mean coverage for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), 146.53, 'row mean score for sense strand');



### Move to third row
# test row positioned score using bam file
$row = $stream->next_row;
is($row->name, 'YAL044W-A', 'row name');

my %pos2scores = $row->get_position_scores(
	'dataset'  => $dataset,
	'value'    => 'count',
);
is(scalar keys %pos2scores, 111, 'number of positioned scores');
# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is($pos2scores{2}, 1, 'positioned score at 1');
is($pos2scores{20}, 2, 'positioned score at 20');

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
is(scalar keys %pos2scores, 64, 'number of positioned scores');
is($pos2scores{57556}, 2, 'positioned score at 57556');
is($pos2scores{57840}, 1, 'positioned score at 57840');
