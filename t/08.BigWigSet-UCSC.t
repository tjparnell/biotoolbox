#!/usr/bin/perl -w

# Test script for Bio::ToolBox::Data
# working with BigWigSet data

use strict;
use Test::More;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	if ( eval { require Bio::DB::BigWigSet; 1 } ) {
		plan tests => 43;
	}
	else {
		plan skip_all => 'Optional module Bio::DB::BigWigSet not available';
	}
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, "Data", "biotoolbox.cfg" );
}

require_ok 'Bio::ToolBox::Data'
	or BAIL_OUT "Cannot load Bio::ToolBox::Data";
use_ok( 'Bio::ToolBox::db_helper', 'get_chromosome_list' );

my $dataset = File::Spec->catfile( $Bin, "Data", "sample3" );

### Open a test file
my $infile = File::Spec->catfile( $Bin, "Data", "sample.bed" );
my $Data   = Bio::ToolBox::Data->new( file => $infile );
isa_ok( $Data, 'Bio::ToolBox::Data', 'BED Data' );

# add a database
is( $Data->big_adapter('ucsc'), 'ucsc', 'set preferred database adapter to ucsc' );
$Data->database($dataset);
is( $Data->database, $dataset, 'get database' );
my $db = $Data->open_database;
isa_ok( $db, 'Bio::DB::BigWigSet', 'connected database' );

# check chromosomes
my @chromos = get_chromosome_list($db);
is( scalar @chromos, 1,      'number of chromosomes' );
is( $chromos[0][0],  'chrI', 'name of first chromosome' );
is( $chromos[0][1],  230208, 'length of first chromosome' );

### Initialize row stream
my $stream = $Data->row_stream;
isa_ok( $stream, 'Bio::ToolBox::Data::Iterator', 'row stream iterator' );

# First row is YAL047C
my $row = $stream->next_row;
is( $row->name, 'YAL047C', 'row name' );

# try a segment
my $segment = $row->segment;
isa_ok( $segment, 'Bio::DB::BigWigSet::Segment', 'row segment' );
is( $segment->start, 54989, 'segment start' );

# score count sum
my $score = $row->get_score(
	'db'      => $dataset,
	'dataset' => 'sample3',
	'method'  => 'count',
);

# print "count sum for ", $row->name, " is $score\n";
is( $score, 435, 'row sum of count' )
	or diag("if this test fails, try updating your UCSC kent source library and rebuild");

# this number is erroneous, the actual value is 434, but the summaries count returns
# 435 for whatever reason - sigh

# score mean coverage
$score = $row->get_score(
	'db'      => $db,
	'dataset' => 'sample3',
	'method'  => 'mean',
);

# print "mean coverage for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), 1.19, 'row mean score' )
	or diag("if this test fails, try updating your UCSC kent source library and rebuild");

# postion scores
my %pos2scores2 = $row->get_region_position_scores( 'dataset' => 'sample3', );

# print "found ", scalar keys %pos2scores2, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores2) {
# 	print "  $_ => $pos2scores2{$_}\n";
# }
is( scalar( keys %pos2scores2 ),           434,  'position scores' );
is( sprintf( "%.2f", $pos2scores2{1671} ), 2.72, 'score at position 1671' );
is( $pos2scores2{1787},                    0,    'score at position 1787' );

# min score
$score = $row->get_score(
	'db'      => $db,
	'dataset' => 'sample3',
	'method'  => 'min',
);
is( $score, 0, 'minimum score' );

# max score
$score = $row->get_score(
	'db'      => $db,
	'dataset' => 'sample3',
	'method'  => 'max',
);
is( sprintf( "%.2f", $score ), 4.57, 'maximum score' );

### Move to the next row
$row = $stream->next_row;
is( $row->start,  57029, 'row start position' );
is( $row->strand, -1,    'row strand' );

$score = $row->get_score(
	'dataset'  => 'sample3',
	'method'   => 'median',
	'stranded' => 'all',
);

# print "both strands score median for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), 1.69, 'row median score' )
	or diag("if this test fails, try updating your UCSC kent source library and rebuild");

# try stranded data collection
$score = $row->get_score(
	'dataset'  => 'sample3',
	'method'   => 'median',
	'stranded' => 'sense',
);

# print "sense score median for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), 2.74, 'row sense median score' )
	or diag("if this test fails, try updating your UCSC kent source library and rebuild");

$score = $row->get_score(
	'dataset'  => 'sample3',
	'method'   => 'median',
	'stranded' => 'antisense',
);

# print "antisense score median for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), 0.38, 'row antisense median score' )
	or diag("if this test fails, try updating your UCSC kent source library and rebuild");

### Try positioned score index
my %pos2scores = $row->get_region_position_scores(
	'dataset'  => 'sample3',
	'stranded' => 'sense',
);
is( scalar keys %pos2scores, 44, 'number of positioned scores' );

# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( sprintf( "%.2f", $pos2scores{55} ),  4.16, 'score at position 55' );
is( sprintf( "%.2f", $pos2scores{255} ), 2.03, 'score at position 255' );

### Try relative positioned score index
%pos2scores = $row->get_relative_point_position_scores(
	'dataset'  => 'sample3',
	'stranded' => 'sense',
	'position' => 5,
	'extend'   => 200,
);
is( scalar keys %pos2scores, 49, 'number of relative positioned scores' );

# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( sprintf( "%.2f", $pos2scores{55} ),  4.16, 'score at relative position 55' );
is( sprintf( "%.2f", $pos2scores{-25} ), 1.73, 'score at relative position -25' );

### Generate new genomic window file
undef $Data;
undef $row;
$Data = Bio::ToolBox::Data->new(
	feature => 'genome',
	db      => $dataset,
	win     => 500
);
isa_ok( $Data, 'Bio::ToolBox::Data', 'new genome window file' );
is( $Data->feature,        'genome',     'Data feature name' );
is( $Data->feature_type,   'coordinate', 'Data feature type is coordinate' );
is( $Data->number_columns, 3,            'Data number of columns' );
is( $Data->number_rows,    461,          'Data number of rows' );
$row = $Data->get_row(1);
isa_ok( $row, 'Bio::ToolBox::Data::Feature', 'First row object' );
is( $row->start,  1,   'First row start coordinate' );
is( $row->stop,   500, 'First row stop coordinate' );
is( $row->length, 500, 'First row length' );
$row = $Data->get_row(461);
isa_ok( $row, 'Bio::ToolBox::Data::Feature', 'Last row object' );
is( $row->start,  230001, 'Last row start coordinate' );
is( $row->length, 208,    'Last row length' );

