#!/usr/bin/env perl

# Test script for Bio::ToolBox::db_helper::big
# working with BigWig data

use Test2::V0 -no_srand => 1;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	if ( eval { require Bio::DB::Big; 1 } ) {
		plan(35);
	}
	else {
		skip_all('Optional module Bio::DB::Big not available');
	}
	## no critic
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, "Data", "biotoolbox.cfg" );
	## use critic
}

use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(get_chromosome_list);

# Data file
my $dataset = File::Spec->catfile( $Bin, "Data", "sample2.bw" );

### Open a test file
my $infile = File::Spec->catfile( $Bin, "Data", "sample.bed" );
my $Data   = Bio::ToolBox::Data->new( file => $infile );
isa_ok( $Data, ['Bio::ToolBox::Data'], 'got a bed Data object' );

# add a database
is( $Data->big_adapter('big'), 'big', 'set preferred database adapter to big' );
$Data->database($dataset);
is( $Data->database, $dataset, 'get database' );
my $db = $Data->open_database;
isa_ok( $db, ['Bio::DB::Big::File'], 'got a BigWig database object' );

# check chromosomes
my @chromos = get_chromosome_list($db);
is( scalar @chromos, 1,      'number of chromosomes' );
is( $chromos[0][0],  'chrI', 'name of first chromosome' );
is( $chromos[0][1],  230208, 'length of first chromosome' );

### Initialize row stream
my $stream = $Data->row_stream;
isa_ok( $stream, ['Bio::ToolBox::Data::Iterator'], 'got a row stream iterator object' );

# First row is YAL047C
my $row = $stream->next_row;
is( $row->name, 'YAL047C', 'row name' );

# segments not used with Bio::DB::Big::File

# score count sum
my $score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'count',
);

# print "count sum for ", $row->name, " is $score\n";
is( $score, 49, 'row sum of count' );

# score mean coverage
$score = $row->get_score(
	'db'      => $db,
	'dataset' => $dataset,
	'method'  => 'mean',
);

# print "mean coverage for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), -0.12, 'row mean score' );

# score min coverage
$score = $row->get_score(
	'db'      => $db,
	'dataset' => $dataset,
	'method'  => 'min',
);
is( sprintf( "%.2f", $score ), -0.62, 'row min score' );

# score max coverage
$score = $row->get_score(
	'db'      => $db,
	'dataset' => $dataset,
	'method'  => 'max',
);
is( sprintf( "%.2f", $score ), 0.52, 'row max score' );

### Move to the next row
$row = $stream->next_row;
is( $row->start,  57029, 'row start position' );
is( $row->strand, -1,    'row strand' );

$score = $row->get_score(
	'dataset' => $dataset,
	'method'  => 'count',
);

# print "score count sum for ", $row->name, " is $score\n";
is( $score, 7, 'row count sum' )
	or diag("if this test fails, try updating your UCSC kent source library and rebuild");

$score = $row->get_score(
	'dataset' => $dataset,
	'value'   => 'score',
	'method'  => 'median',
);

# print "score median for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), '0.50', 'row median score' )
	or diag("if this test fails, try updating your UCSC kent source library and rebuild");

### Try positioned score index
my %pos2scores = $row->get_region_position_scores(
	'dataset' => $dataset,
	'value'   => 'score',
);
is( scalar keys %pos2scores, 7, 'number of positioned scores' );

# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( sprintf( "%.2f", $pos2scores{8} ),   '0.50', 'positioned score at 8' );
is( sprintf( "%.2f", $pos2scores{142} ), 0.58,   'positioned score at 142' );
undef %pos2scores;

### Try relative positioned score index
my $pos2scores = $row->get_relative_point_position_scores(
	'dataset'  => $dataset,
	'value'    => 'score',
	'position' => 5,
	'extend'   => 200,
);
is( scalar keys %{$pos2scores}, 9, 'number of relative positioned scores' );

# print "found ", scalar keys %{$pos2scores}, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %{$pos2scores}) {
# 	print "  $_ => $pos2scores->{$_}\n";
# }
is( sprintf( "%.2f", $pos2scores->{-114} ), 0.41, 'relative positioned score at -114' );
is( sprintf( "%.2f", $pos2scores->{96} ),   0.48, 'relative positioned score at 96' );

# BigWig does not support stranded data collection
# save that for BigWigSet

### Generate new genomic window file
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
is( $Data->number_rows,    461,          'Data number of rows' );
$row = $Data->get_row(1);
isa_ok( $row, ['Bio::ToolBox::Data::Feature'], 'First row Feature object' );
is( $row->start,  1,   'First row start coordinate' );
is( $row->stop,   500, 'First row stop coordinate' );
is( $row->length, 500, 'First row length' );
$row = $Data->get_row(461);
isa_ok( $row, ['Bio::ToolBox::Data::Feature'], 'Last row Feature object' );
is( $row->start,  230001, 'Last row start coordinate' );
is( $row->length, 208,    'Last row length' );

