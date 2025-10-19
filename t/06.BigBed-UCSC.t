#!/usr/bin/perl -w

# Test script for Bio::ToolBox::Data
# working with BigBed data

use strict;
use Test::More;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	if ( eval { require Bio::DB::BigBed; 1 } ) {
		plan tests => 50;
	}
	else {
		plan skip_all => 'Alternate module Bio::DB::BigBed not available';
	}
	## no critic
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, "Data", "biotoolbox.cfg" );
	## use critic
}

require_ok 'Bio::ToolBox::Data'
	or BAIL_OUT "Cannot load Bio::ToolBox::Data";
use_ok( 'Bio::ToolBox::db_helper', 'check_dataset_for_rpm_support',
	'get_chromosome_list' );

my $dataset = File::Spec->catfile( $Bin, "Data", "sample1.bb" );

### Open a test file
my $infile = File::Spec->catfile( $Bin, "Data", "sample.bed" );
my $Data   = Bio::ToolBox::Data->new( file => $infile );
isa_ok( $Data, 'Bio::ToolBox::Data', 'BED Data' );

# add a database
is( $Data->big_adapter('ucsc'), 'ucsc', 'set preferred database adapter to ucsc' );
$Data->database($dataset);
is( $Data->database, $dataset, 'get database' );
my $db = $Data->open_database;
isa_ok( $db, 'Bio::DB::BigBed', 'connected database' );

# check chromosomes
my @chromos = get_chromosome_list($db);
is( scalar @chromos, 1,      'number of chromosomes' );
is( $chromos[0][0],  'chrI', 'name of first chromosome' );
is( $chromos[0][1],  230208, 'length of first chromosome' );

# check total mapped alignments
my $total = check_dataset_for_rpm_support($dataset);
is( $total, 1414, "number of features in BigBed" );

### Initialize row stream
my $stream = $Data->row_stream;
isa_ok( $stream, 'Bio::ToolBox::Data::Iterator', 'row stream iterator' );

# First row is YAL047C
my $row = $stream->next_row;
is( $row->name, 'YAL047C', 'row name' );

# try a segment
my $segment = $row->segment;
isa_ok( $segment, 'Bio::DB::BigFile::Segment', 'row segment' );
is( $segment->start, 54989, 'segment start' );

# read count sum
my $score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'count',
);

# print "count sum for ", $row->name, " is $score\n";
is( $score, 453, 'row read count score' );

# mean coverage
$score = $row->get_score(
	'db'      => $db,
	'dataset' => $dataset,
	'method'  => 'mean',
);

# print "mean coverage for ", $row->name, " is $score\n";
is( sprintf( "%.1f", $score ), 143.8, 'row mean score' );

# read precise count sum
$score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'pcount',
);

# print "count sum for ", $row->name, " is $score\n";
is( $score, 414, 'row precise count' );

# read named count sum
$score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'ncount',
);

# print "count sum for ", $row->name, " is $score\n";
is( $score, 453, 'row named count' );

### Move to the next row
$row = $stream->next_row;
is( $row->start,  57029, 'row start position' );
is( $row->strand, -1,    'row strand' );

# try stranded data collection
$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'count',
	'stranded' => 'all',
);

# print "all read count sum for ", $row->name, " is $score\n";
is( $score, 183, 'row sum of count score for all strands' );

$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'count',
	'stranded' => 'sense',
);

# print "sense read count sum for ", $row->name, " is $score\n";
is( $score, 86, 'row sum of count score for sense strand' );

$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'count',
	'stranded' => 'antisense',
);

# print "antisense read count sum for ", $row->name, " is $score\n";
is( $score, 97, 'row sum of count score for antisense strand' );

$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'mean',
	'stranded' => 'sense',
);

# print "sense mean coverage for ", $row->name, " is $score\n";
is( sprintf( "%.1f", $score ), 146.9, 'row mean score for sense strand' );

$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'mean',
	'stranded' => 'antisense',
);

# print "antisense mean coverage for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), 146.53, 'row mean score for sense strand' );

### Move to third row
# test row positioned score using bam file
$row = $stream->next_row;
is( $row->name, 'YAL044W-A', 'row name' );

my %pos2scores = $row->get_region_position_scores(
	'dataset' => $dataset,
	'method'  => 'count',
);

# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 111, 'number of positioned scores' );
is( $pos2scores{2},          1,   'positioned score at 1' );
is( $pos2scores{20},         2,   'positioned score at 20' );

%pos2scores = $row->get_region_position_scores(
	'dataset' => $dataset,
	'method'  => 'pcount',
);

# print "found ", scalar keys %pos2scores, " positioned precise counts \n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 12, 'number of positioned precise counts' );
is( $pos2scores{2},          1,  'positioned precise count at 1' );
is( $pos2scores{15},         2,  'positioned precise count at 15' );

%pos2scores = $row->get_region_position_scores(
	'dataset' => $dataset,
	'method'  => 'ncount',
);

# print "found ", scalar keys %pos2scores, " positioned named counts\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 111, 'number of positioned named counts' );
is( $pos2scores{7},          2,   'positioned named count at 7' );
is( $pos2scores{36},         1,   'positioned named count at 36' );

%pos2scores = $row->get_region_position_scores(
	'dataset'  => $dataset,
	'method'   => 'count',
	'absolute' => 1,
	'stranded' => 'antisense',
);

# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 64, 'number of positioned scores' );
is( $pos2scores{57556},      2,  'positioned score at 57556' );
is( $pos2scores{57840},      1,  'positioned score at 57840' );

# Generate new genomic window file
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

