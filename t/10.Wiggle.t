#!/usr/bin/env perl

# Test script for Bio::ToolBox::db_helper::wiggle
# working with binary wiggle data

use Test2::V0 -no_srand => 1;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	if ( eval { require Bio::Graphics::Wiggle; 1 } ) {
		plan(16);
	}
	else {
		skip_all('Optional module Bio::Graphics not available');
	}
	## no critic
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, "Data", "biotoolbox.cfg" );
	## use critic
}

use IO::File;
use Bio::ToolBox::Data;

### Prepare the GFF database
my $database = File::Spec->catfile( $Bin, "Data", "sample2_wib.gff3" );
my $wigfile  = File::Spec->catfile( $Bin, "Data", "sample2.wib" );
unless ( -e $database ) {
	my $fh = IO::File->new( $database, "w" )
		or die "unable to write temporary GFF file!";
	my $gff = <<OUT;
##gff-version 3
chrI	SGD	chromosome	1	230218	.	.	.	ID=chrI;dbxref=NCBI:NC_001133;Name=chrI
chrI	tim	sample2	1	230218	.	.	.	Name=sample2;wigfile=$wigfile
OUT
	$fh->print($gff);
	$fh->close;
}

### Open a test file
my $infile = File::Spec->catfile( $Bin, "Data", "sample.bed" );
my $Data   = Bio::ToolBox::Data->new( file => $infile );
isa_ok( $Data, ['Bio::ToolBox::Data'], 'got a bed Data object' );

# add a database
$Data->database($database);
is( $Data->database, $database, 'get database' );
my $db = $Data->open_database;
isa_ok( $db, ['Bio::DB::SeqFeature::Store'], 'got a SeqFeature Store database object' );

### Initialize row stream
my $stream = $Data->row_stream;
isa_ok( $stream, ['Bio::ToolBox::Data::Iterator'], 'got row stream Iterator ojbect' );

# First row is YAL047C
my $row = $stream->next_row;
is( $row->name, 'YAL047C', 'row name' );

# try a segment
my $segment = $row->segment;
isa_ok( $segment, ['Bio::DB::SeqFeature::Segment'], 'got a row Segment object' );
is( $segment->start, 54989, 'segment start' );

# score count sum
my $score = $row->get_score(
	'db'      => $database,
	'dataset' => 'sample2',
	'method'  => 'count',
);

# print "count sum for ", $row->name, " is $score\n";
is( $score, 49, 'row sum of count' );

# score mean coverage
$score = $row->get_score(
	'db'      => $db,
	'dataset' => 'sample2',
	'method'  => 'mean',
);

# print "mean coverage for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), -0.14, 'row mean score' );

### Move to the next row
$row = $stream->next_row;
is( $row->start,  57029, 'row start position' );
is( $row->strand, -1,    'row strand' );

# try stranded data collection
$score = $row->get_score(
	'dataset' => 'sample2',
	'method'  => 'count',
);

# print "score count sum for ", $row->name, " is $score\n";
is( $score, 7, 'row count sum' );

$score = $row->get_score(
	'dataset' => 'sample2',
	'method'  => 'median',
);

# print "score median for ", $row->name, " is $score\n";
is( sprintf( "%.2f", $score ), '0.49', 'row median score' );

### Try positioned score index
my %pos2scores = $row->get_region_position_scores( 'dataset' => 'sample2', );
is( scalar keys %pos2scores, 7, 'number of positioned scores' );

# print "found ", scalar keys %pos2scores, " positions with reads\n";
# foreach (sort {$a <=> $b} keys %pos2scores) {
# 	print "  $_ => $pos2scores{$_}\n";
# }
is( sprintf( "%.2f", $pos2scores{8} ),   '0.49', 'positioned score at 8' );
is( sprintf( "%.2f", $pos2scores{142} ), 0.58,   'positioned score at 142' );

END {
	unlink $database;
}

