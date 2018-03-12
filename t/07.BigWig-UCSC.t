#!/usr/bin/perl -w

# Test script for Bio::ToolBox::Data 
# working with BigWig data

use strict;
use Test::More;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	if (eval {require Bio::DB::BigWig; 1}) {
		plan tests => 24;
	}
	else {
		plan skip_all => 'Optional module Bio::DB::BigWig not available';
	}
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile($Bin, "Data", "biotoolbox.cfg");
}

require_ok 'Bio::ToolBox::Data' or 
	BAIL_OUT "Cannot load Bio::ToolBox::Data";
use_ok( 'Bio::ToolBox::db_helper', 'get_chromosome_list' );


my $dataset = File::Spec->catfile($Bin, "Data", "sample2.bw");

### Open a test file
my $infile = File::Spec->catfile($Bin, "Data", "sample.bed");
my $Data = Bio::ToolBox::Data->new(file => $infile);
isa_ok($Data, 'Bio::ToolBox::Data', 'BED Data');

# add a database
$Data->database($dataset);
is($Data->database, $dataset, 'get database');
my $db = $Data->open_database;
isa_ok($db, 'Bio::DB::BigWig', 'connected database');

# check chromosomes
my @chromos = get_chromosome_list($db);
is(scalar @chromos, 1, 'number of chromosomes');
is($chromos[0][0], 'chrI', 'name of first chromosome');
is($chromos[0][1], 230208, 'length of first chromosome');


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

# score count sum
my $score = $row->get_score(
	'db'       => $dataset,
	'dataset'  => $dataset,
	'method'   => 'count',
);
# print "count sum for ", $row->name, " is $score\n";
is($score, 49, 'row sum of count') or 
	diag("if this test fails, try updating your UCSC kent source library and rebuild");

# score mean coverage
$score = $row->get_score(
	'db'       => $db,
	'dataset'  => $dataset,
	'method'   => 'mean',
);
# print "mean coverage for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), -0.12, 'row mean score') or 
	diag("if this test fails, try updating your UCSC kent source library and rebuild");



### Move to the next row
$row = $stream->next_row;
is($row->start, 57029, 'row start position');
is($row->strand, -1, 'row strand');

$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'count',
);
# print "score count sum for ", $row->name, " is $score\n";
is($score, 7, 'row count sum') or 
	diag("if this test fails, try updating your UCSC kent source library and rebuild");

$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'score',
	'method'   => 'median',
);
# print "score median for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), '0.50', 'row median score') or 
	diag("if this test fails, try updating your UCSC kent source library and rebuild");




### Try positioned score index
my %pos2scores = $row->get_region_position_scores(
	'dataset'  => $dataset,
	'value'    => 'score',
);
is(scalar keys %pos2scores, 7, 'number of positioned scores');
# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is(sprintf("%.2f", $pos2scores{8}), '0.50', 'positioned score at 8');
is(sprintf("%.2f", $pos2scores{142}), 0.58, 'positioned score at 142');
undef %pos2scores;



### Try relative positioned score index
my $pos2scores = $row->get_relative_point_position_scores(
	'dataset'  => $dataset,
	'value'    => 'score',
	'position' => 5,
	'extend'   => 200,
);
is(scalar keys %$pos2scores, 9, 'number of relative positioned scores');
# print "found ", scalar keys %$pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %$pos2scores) {
# 	print "  $_ => $pos2scores->{$_}\n";
# }
is(sprintf("%.2f", $pos2scores->{-114}), 0.41, 'relative positioned score at -114');
is(sprintf("%.2f", $pos2scores->{96}), 0.48, 'relative positioned score at 96');


# BigWig does not support stranded data collection
# save that for BigWigSet


