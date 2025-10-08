#!/usr/bin/perl -w

# Test script for Bio::ToolBox::Data
# working with Bam data

use strict;
use Test::More;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	if ( eval { require Bio::DB::Sam; 1 } ) {
		plan tests => 71;
	}
	else {
		plan skip_all => 'Optional module Bio::DB::Sam not available';
	}
	## no critic
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, "Data", "biotoolbox.cfg" );
	## use critic
}

require_ok 'Bio::ToolBox::Data'
	or BAIL_OUT "Cannot load Bio::ToolBox::Data";
use_ok(
	'Bio::ToolBox::db_helper', 'check_dataset_for_rpm_support',
	'get_chromosome_list',     'get_genomic_sequence',
	'use_minimum_mapq'
);

my $dataset = File::Spec->catfile( $Bin, "Data", "sample1.bam" );
my $fasta   = File::Spec->catfile( $Bin, "Data", 'sequence.fa' );

### Open a test file
my $infile = File::Spec->catfile( $Bin, "Data", "sample.bed" );
my $Data   = Bio::ToolBox::Data->new( file => $infile );
isa_ok( $Data, 'Bio::ToolBox::Data', 'BED Data' );

# add a database
is( $Data->bam_adapter('sam'), 'sam', 'set preferred database adapter to sam' );
$Data->database($dataset);
is( $Data->database, $dataset, 'get database' );
my $db = $Data->open_database;
isa_ok( $db, 'Bio::DB::Sam', 'connected database' );

# check chromosomes
my @chromos = get_chromosome_list($db);
is( scalar @chromos, 1,      'bam number of chromosomes' );
is( $chromos[0][0],  'chrI', 'bam name of first chromosome' );
is( $chromos[0][1],  230208, 'bam length of first chromosome' );

# check total mapped alignments
my $total = check_dataset_for_rpm_support($dataset);
is( $total, 1414, "number of mapped alignments in bam" );

# check mapping quality
is( use_minimum_mapq(),    0,   'default minimum mapq' );
is( use_minimum_mapq(100), 100, 'update minimum mapq' );
$total = check_dataset_for_rpm_support($dataset);

# this does not change because the value is cached
is( $total,              1414, "number of total alignments with high mapq" );
is( use_minimum_mapq(0), 0,    'reset minimum mapq' );    # reset for below

# check fasta
my $fdb = $Data->open_new_database($fasta);
isa_ok( $fdb, 'Bio::DB::Sam::Fai', 'Sam Fai fasta database' );
my $seq = get_genomic_sequence( $fdb, 'chrI', 257, 275 );
is( $seq, 'ACCCTACCATTACCCTACC', 'fetched fasta sequence' );
unlink "$fasta.fai";

### Initialize row stream
my $stream = $Data->row_stream;
isa_ok( $stream, 'Bio::ToolBox::Data::Iterator', 'row stream iterator' );

# First row is YAL047C
my $row = $stream->next_row;
is( $row->name, 'YAL047C', 'row name' );

# try a segment
my $segment = $row->segment;
isa_ok( $segment, 'Bio::DB::Sam::Segment', 'row segment' );
is( $segment->start, 54989, 'segment start' );

# read count sum with default mapq
my $score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'count',
);
is( $score, 453, 'row sum of read count score with low mapq' );

# check with high mapq
use_minimum_mapq(100);
$score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'count',
);
is( $score, 429, 'read count score with high mapq' );
use_minimum_mapq(0);    # reset again

# mean coverage
$score = $row->get_score(
	'db'      => $db,
	'dataset' => $dataset,
	'method'  => 'mean',
);
is( sprintf( "%.2f", $score ), 16.33, 'row mean coverage with low mapq' );

# check mean coverage with high mapq – should not change
use_minimum_mapq(100);
$score = $row->get_score(
	'db'      => $db,
	'dataset' => $dataset,
	'method'  => 'mean',
);
is( sprintf( "%.2f", $score ), 16.33, 'row mean coverage with high mapq' );
use_minimum_mapq(0);

# read precise count sum
$score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'pcount',
);
is( $score, 414, 'row sum of read precise count score with low mapq' );

# read precise count sum with high mapq
use_minimum_mapq(100);
$score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'pcount',
);
is( $score, 398, 'row sum of read precise count with high mapq' );
use_minimum_mapq(0);

# read ncount sum
$score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'ncount',
);
is( $score, 453, 'row read name count score with low mapq' );

# read ncount sum with high mapq
use_minimum_mapq(100);
$score = $row->get_score(
	'db'      => $dataset,
	'dataset' => $dataset,
	'method'  => 'ncount',
);
is( $score, 429, 'row read name count score with high mapq' );
use_minimum_mapq(0);

### Move to the next row
$row = $stream->next_row;
is( $row->start,  57029, 'row start position' );
is( $row->strand, -1,    'row strand' );

# stranded alignment count
$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'count',
	'stranded' => 'all',
);
is( $score, 183, 'row sum of count score for all strands' );

$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'count',
	'stranded' => 'sense',
);
is( $score, 86, 'row sum of count score for sense strand' );

$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'count',
	'stranded' => 'antisense',
);
is( $score, 97, 'row sum of count score for antisense strand' );

# stranded coverage
# low level coverage does not support stranded collection so values should be identical
$score = $row->get_score(
	'dataset'  => $dataset,
	'method'   => 'mean',
	'stranded' => 'sense',
);
is( sprintf( "%.2f", $score ), 29.38, 'row mean coverage for sense strand' );

$score = $row->get_score(
	'dataset'  => $dataset,
	'value'    => 'score',
	'method'   => 'mean',
	'stranded' => 'antisense',
);
is( sprintf( "%.2f", $score ), 29.38, 'row mean coverage for antisense strand' );

### Move to third row
# test row positioned score using bam file
$row = $stream->next_row;
is( $row->name, 'YAL044W-A', 'row name' );

# positioned scores using standard count
my %pos2scores = $row->get_region_position_scores(
	'dataset' => $dataset,
	'method'  => 'count',
);

# print " > found ", scalar keys %pos2scores, " positions with reads at MAPQ 0\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  > $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 150, 'number of positioned scores' );
is( $pos2scores{101},        2,   'positioned count at 101' );
is( $pos2scores{-21},        2,   'positioned count at -21' );

# repeat with high mapping quality
use_minimum_mapq(100);
%pos2scores = $row->get_region_position_scores(
	'dataset' => $dataset,
	'method'  => 'count',
);

# print " > found ", scalar keys %pos2scores, " positions with reads at MAPQ 30\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  > $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 145, 'number of positioned scores with high mapq' );
is( $pos2scores{101},        1,   'positioned count at 101 with high mapq' );
is( $pos2scores{-21},        2,   'positioned count at -21 with high mapq' );
use_minimum_mapq(0);

# positioned count with precise count
%pos2scores = $row->get_region_position_scores(
	'dataset' => $dataset,
	'method'  => 'pcount',
);

# print " > found ", scalar keys %pos2scores, " positions with precise reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  > $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 89, 'number of precise positioned scores' );
is( $pos2scores{72},         2,  'pcount at position 72' );
is( $pos2scores{118},        1,  'pcount at position 118' );

# repeat pcount positioned data with high mapping quality
use_minimum_mapq(100);
%pos2scores = $row->get_region_position_scores(
	'dataset' => $dataset,
	'method'  => 'pcount',
);

# print " > found ", scalar keys %pos2scores, " positions with precise reads at high mapq\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  > $_ => $pos2scores{$_}\n";
# }
is( scalar keys %pos2scores, 87, 'number of precise positioned scores with high mapq' );
is( $pos2scores{72},         1,  'pcount at pos 72 with high mapq' );
is( ( exists $pos2scores{118} ), q(), 'pcount missing at pos 118 with high mapq' );
use_minimum_mapq(0);

# positioned count data with names
%pos2scores = $row->get_region_position_scores(
	'dataset' => $dataset,
	'method'  => 'ncount',
);

# print " > found ", scalar keys %pos2scores, " positions of named reads at MAPQ 0\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	printf "  > $_ => %s\n", join(',', @{$pos2scores{$_}});
# }
is( scalar keys %pos2scores, 150, 'number of named positioned scores' );
is( $pos2scores{6}->[0], 'HWI-EAS240_0001:7:64:6158:10466#0/1', 'positioned named at 6' );
is( scalar @{ $pos2scores{101} }, 2, 'positioned name count at 101' );

# repeat ncount positioned data with high mapping quality
use_minimum_mapq(100);
%pos2scores = $row->get_region_position_scores(
	'dataset' => $dataset,
	'method'  => 'ncount',
);

# print " > found ", scalar keys %pos2scores, " positions of named reads at MAPQ 100\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	printf "  > $_ => %s\n", join(',', @{$pos2scores{$_}});
# }
is( scalar keys %pos2scores, 145, 'number of named positioned scores with high MAPQ' );
is(
	$pos2scores{6}->[0], 'HWI-EAS240_0001:7:64:6158:10466#0/1',
	'positioned named at 6 with high MAPQ'
);
is( scalar @{ $pos2scores{101} }, 1, 'positioned name count at 101 with high MAPQ' );
use_minimum_mapq(0);

# Fetch alignments
my $alignment_data = { mapq => [] };
my $callback       = sub {
	my ( $a, $data ) = @_;
	push @{ $data->{mapq} }, $a->qual;
};
my $f = $row->fetch_alignments(
	'db'       => $db,
	'data'     => $alignment_data,
	'callback' => $callback,
);
is( $f, 1, 'Raw alignment fetch' );

# printf "found %d alignments\n", scalar(@{$alignment_data->{mapq}});
# for my $i (0..183) {
# 	printf "$i\t%d\n", $alignment_data->{mapq}->[$i];
# }
is( scalar( @{ $alignment_data->{mapq} } ), 184, 'Raw fetch mapq number' );
is( $alignment_data->{mapq}->[0],           150, 'First alignment mapq' );
is( $alignment_data->{mapq}->[66],          0,   '66th alignment mapq' );
is( $alignment_data->{mapq}->[80],          95,  '80th alignment mapq' );

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

