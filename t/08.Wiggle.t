#!/usr/bin/env perl

# Test script for Bio::ToolBox::Data 
# working with binary wiggle data

use strict;
use Test::More;
use FindBin '$Bin';

BEGIN {
	if (eval {require Bio::Graphics::Wiggle; 1}) {
		plan tests => 17;
	}
	else {
		plan skip_all => 'Optional module Bio::Graphics not available';
	}
	$ENV{'BIOTOOLBOX'} = "$Bin/Data/biotoolbox.cfg";
}

use lib "$Bin/../lib";
require_ok 'Bio::ToolBox::Data' or 
	BAIL_OUT "Cannot load Bio::ToolBox::Data";


### Prepare the GFF database
my $database = "$Bin/Data/sample2.gff3";
unless (-e $database) {
	open(my $fh, ">" ,$database);
	print $fh <<GFF
##gff-version 3
chrI	SGD	chromosome	1	230218	.	.	.	ID=chrI;dbxref=NCBI:NC_001133;Name=chrI
chrI	tim	sample2	1	230218	.	.	.	Name=sample2;wigfile=$Bin/Data/sample2.wib
GFF
;
	close $fh;
}


### Open a test file
my $Data = Bio::ToolBox::Data->new(file => "$Bin/Data/sample.bed");
isa_ok($Data, 'Bio::ToolBox::Data', 'BED Data');

# add a database
$Data->database($database);
is($Data->database, $database, 'get database');
my $db = $Data->open_database;
isa_ok($db, 'Bio::DB::SeqFeature::Store', 'connected database');



### Initialize row stream
my $stream = $Data->row_stream;
isa_ok($stream, 'Bio::ToolBox::Data::Iterator', 'row stream iterator');

# First row is YAL047C
my $row = $stream->next_row;
is($row->name, 'YAL047C', 'row name');

# try a segment
my $segment = $row->segment;
isa_ok($segment, 'Bio::DB::SeqFeature::Segment', 'row segment');
is($segment->start, 54989, 'segment start');

# score count sum
my $score = $row->get_score(
	'db'       => $database,
	'dataset'  => 'sample2',
	'value'    => 'count',
	'method'   => 'sum',
);
# print "count sum for ", $row->name, " is $score\n";
is($score, 49, 'row sum of count');

# score mean coverage
$score = $row->get_score(
	'db'       => $db,
	'dataset'  => 'sample2',
	'value'    => 'score',
	'method'   => 'mean',
);
# print "mean coverage for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), -0.14, 'row mean score');



### Move to the next row
$row = $stream->next_row;
is($row->start, 57029, 'row start position');
is($row->strand, -1, 'row strand');

# try stranded data collection
$score = $row->get_score(
	'dataset'  => 'sample2',
	'value'    => 'count',
	'method'   => 'sum',
);
# print "score count sum for ", $row->name, " is $score\n";
is($score, 7, 'row count sum');

$score = $row->get_score(
	'dataset'  => 'sample2',
	'value'    => 'score',
	'method'   => 'median',
);
# print "score median for ", $row->name, " is $score\n";
is(sprintf("%.2f", $score), '0.49', 'row median score');




### Try positioned score index
my %pos2scores = $row->get_position_scores(
	'dataset'  => 'sample2',
	'value'    => 'score',
);
is(scalar keys %pos2scores, 7, 'number of positioned scores');
# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is(sprintf("%.2f", $pos2scores{8}), '0.49', 'positioned score at 8');
is(sprintf("%.2f", $pos2scores{142}), 0.58, 'positioned score at 142');


END {
	unlink $database;
}

