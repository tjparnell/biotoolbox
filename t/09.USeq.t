#!/usr/bin/env perl

# Test script for Bio::ToolBox::db_helper::useq
# working with USeq data

use Test2::V0 -no_srand => 1;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	if ( eval { require Bio::DB::USeq; 1 } ) {
		plan(46);
	}
	else {
		skip_all('Optional module Bio::DB::USeq not available');
	}
	## no critic
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, "Data", "biotoolbox.cfg" );
	## use critic
}

use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(get_chromosome_list);

# Data file
my $dataset = File::Spec->catfile( $Bin, "Data", "sample3.useq" );

### Open a test file
my $infile = File::Spec->catfile( $Bin, "Data", "sample.bed" );
my $Data   = Bio::ToolBox::Data->new( file => $infile );
isa_ok( $Data, ['Bio::ToolBox::Data'], 'got a bed Data object' );

# add a database
$Data->database($dataset);
is( $Data->database, $dataset, 'get database' );
my $db = $Data->open_database;
isa_ok( $db, ['Bio::DB::USeq'], 'got a USeq database object' );

# check chromosomes
my @chromos = get_chromosome_list($db);
is( scalar @chromos, 1,      'number of chromosomes' );
is( $chromos[0][0],  'chrI', 'name of first chromosome' );
is( $chromos[0][1],  60997,  'length of first chromosome' );

### Initialize row stream
my $stream = $Data->row_stream;
isa_ok( $stream, ['Bio::ToolBox::Data::Iterator'], 'got a row stream iterator object' );

# First row is YAL047C
my $row = $stream->next_row;
is( $row->name, 'YAL047C', 'row name' );

# try a segment
my $segment = $row->segment;
isa_ok( $segment, ['Bio::DB::USeq::Segment'], 'got a USeq Segment object' );
is( $segment->start, 54989, 'segment start' );

# score count
my $score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'count',
);

# print "count sum for ", $row->name, " is $score\n";
is( $score, 437, 'row score count' );

# score precise count
$score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'pcount',
);

# print "precise count for ", $row->name, " is $score\n";
is( $score, 431, 'row precise count' );

# score name count
$score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'ncount',
);

# print "name count for ", $row->name, " is $score\n";
is( $score, 437, 'row name count' );

# score mean coverage
$score = $row->get_score(
	'db'      => $db,
	'dataset' => $dataset,
	'method'  => 'mean',
);

# print "mean coverage for ", $row->name, " is $score\n";
is( sprintf( "%.1f", $score ), 1.2, 'row mean score' );

### Move to the next row
$row = $stream->next_row;
is( $row->start,  57029, 'row start position' );
is( $row->strand, -1,    'row strand' );

$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'median',
	'stranded' => 'all',
);

# print "both strands score median for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), 1.67, 'row median score' );

# try stranded data collection
$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'median',
	'stranded' => 'sense',
);

# print "sense score median for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), 2.71, 'row sense median score' );

$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'median',
	'stranded' => 'antisense',
);

# print "antisense score median for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), 0.38, 'row antisense median score' );

### Try positioned score index
my %pos2scores = $row->get_region_position_scores(
	'dataset'  => $dataset,
	'stranded' => 'sense',
);

# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 44, 'number of positioned scores' );
is( sprintf( "%.2f", $pos2scores{55} ),  4.16, 'score at position 55' );
is( sprintf( "%.2f", $pos2scores{255} ), 2.03, 'score at position 255' );

# positioned count
%pos2scores = $row->get_region_position_scores(
	'dataset'  => $dataset,
	'stranded' => 'sense',
	'method'   => 'count',
);

# print "found ", scalar keys %pos2scores, " positions with counts\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 44, 'number of positioned counts' );
is( $pos2scores{47},         1,  'count at position 47' );
is( $pos2scores{167},        1,  'count at position 167' );

# positioned pcount
%pos2scores = $row->get_region_position_scores(
	'dataset'  => $dataset,
	'stranded' => 'sense',
	'method'   => 'pcount',
);

# print "found ", scalar keys %pos2scores, " positions with precise counts\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 42, 'number of positioned precise counts' );
is( $pos2scores{15},         1,  'precise count at position 15' );
is( $pos2scores{127},        1,  'precise count at position 127' );

# positioned ncount
%pos2scores = $row->get_region_position_scores(
	'dataset'  => $dataset,
	'stranded' => 'sense',
	'method'   => 'ncount',
);

# print "found ", scalar keys %pos2scores, " positions with named counts\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 44, 'number of positioned named counts' );
is( $pos2scores{7},          1,  'named count at position 7' );
is( $pos2scores{335},        1,  'named count at position 335' );

### Try relative positioned score index
%pos2scores = $row->get_relative_point_position_scores(
	'dataset'  => $dataset,
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

# Generate new genomic window file
undef $Data;
undef $row;
$Data = Bio::ToolBox::Data->new(
	feature => 'genome',
	db      => $dataset,
	win     => 500
);
isa_ok( $Data, ['Bio::ToolBox::Data'], 'got a new genome window Data object' );
is( $Data->feature,        'genome',     'Data feature name' );
is( $Data->feature_type,   'coordinate', 'Data feature type is coordinate' );
is( $Data->number_columns, 3,            'Data number of columns' );
is( $Data->number_rows,    122,          'Data number of rows' );
$row = $Data->get_row(1);
isa_ok( $row, ['Bio::ToolBox::Data::Feature'], 'First row Feature object' );
is( $row->start,  1,   'First row start coordinate' );
is( $row->stop,   500, 'First row stop coordinate' );
is( $row->length, 500, 'First row length' );
$row = $Data->get_row(122);
isa_ok( $row, ['Bio::ToolBox::Data::Feature'], 'Last row Feature object' );
is( $row->start,  60501, 'Last row start coordinate' );
is( $row->length, 497,   'Last row length' );

