#!/usr/bin/perl -w

# Test script for Bio::ToolBox::Data
# working with a Bio::DB::SeqFeature::Store database
# using in memory database for simplicity here....

use strict;
use Test::More;
use File::Spec;
use FindBin '$Bin';
use Statistics::Lite qw(min max);

BEGIN {
	if ( eval { require Bio::DB::SeqFeature::Store::memory; 1 } ) {
		plan tests => 88;
	}
	else {
		plan skip_all => 'Bio::DB::SeqFeature::Store not available';
	}
	## no critic
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, "Data", "biotoolbox.cfg" );
	## use critic
}

require_ok 'Bio::ToolBox::Data'
	or BAIL_OUT "Cannot load Bio::ToolBox::Data";

### Collect features from a database
my $infile = File::Spec->catfile( $Bin, "Data", "chrI.gff3" );
my $Data   = Bio::ToolBox::Data->new(
	'db'      => 'mem_test_file',
	'feature' => 'gene:SGD',
);
isa_ok( $Data, 'Bio::ToolBox::Data', 'db collected gene table' );

# check metadata
is( $Data->feature,      'gene:SGD', 'feature' );
is( $Data->feature_type, 'named',    'feature type' );
is( $Data->program,      undef,      'program' );
is( $Data->gff,          0,          'gff format' );
is( $Data->bed,          0,          'bed format' );
is( $Data->filename,     q(),        'filename' );

# check data table
is( $Data->number_columns, 3,  'number of columns' );
is( $Data->last_row,       33, 'last row' );  # this now includes dubious genes as of 1.54

# columns
is( $Data->id_column,   1, 'primary ID column' );
is( $Data->name_column, 2, 'identify name column' );
is( $Data->type_column, 3, 'identify type column' );

# test database
is( $Data->database, 'mem_test_file', 'database name' );
my $db = $Data->open_database;
isa_ok( $db, 'Bio::DB::SeqFeature::Store', 'opened database' );

# since we are dealing with a memory database, the features returned
# are not in a predictable order suitable for testing
# so we will sort the data table by increasing name
$Data->sort_data( 2, 'i' );

# test row_stream
my $stream = $Data->row_stream;
isa_ok( $stream, 'Bio::ToolBox::Data::Iterator', 'row stream iterator' );

# look at first row
my $row = $stream->next_row;
isa_ok( $row, 'Bio::ToolBox::Data::Feature', 'row Feature' );
ok( $row->row_index, 'row index' );           # could be anything
is( $row->type, 'gene:SGD', 'Feature type' );
is( $row->name, 'YAL043C',  'Feature name' );
cmp_ok( $row->id, '>=', 1, 'primary id' );    # in memory db, could be anything?
is( $row->start, 58695, 'Feature start' );    # now automatically collected from db
is( $row->stop,  61052, 'Feature stop' );     # now automatically collected from db

# test SeqFeature
my $feature = $row->feature;
isa_ok( $feature, 'Bio::DB::SeqFeature', 'db SeqFeature from row Feature' );
is( $row->start, 58695, 'Feature start, defined' );    # uses db feature for start
is( $row->stop,  61052, 'Feature stop, defined' );     # uses db feature for stop
my $segment = $row->segment;
isa_ok( $segment, 'Bio::DB::SeqFeature::Segment', 'db segment from row Feature' );

### Using a second test database for data collection
# verify data from a data database
my $ddb      = File::Spec->catfile( $Bin, "Data", "sample2.gff3" );
my $verified = $Data->verify_dataset( 'data:sample2', $ddb );
is( $verified, 'data:sample2', 'dataset in data database verified' );

# collect mean score
my $score = $row->get_score(
	ddb      => $ddb,
	dataset  => 'data',
	'method' => 'mean',
);

# print "mean score with data and ddb is $score\n";
is( sprintf( "%.2f", $score ), 0.18, 'mean score across feature' );

# collect max score
$score = $row->get_score(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	'method' => 'max',
);

# print "max score is $score\n";
is( sprintf( "%.2f", $score ), 0.46, 'max score across feature' );

# collect median score
$score = $row->get_score(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	'method' => 'median',
);

# print "median score is $score\n";
is( sprintf( "%.2f", $score ), 0.19, 'median score across feature' );

# collect 5' mean score
$score = $row->get_score(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	'method' => 'mean',
	start    => $row->stop - 1000,
	stop     => $row->stop,
);

# print "5 prime mean score is $score\n";
is( sprintf( "%.2f", $score ), '0.30', 'mean score across 5 prime feature' );

# collect sum count score
$score = $row->get_score(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	'method' => 'count',
);

# print "median score is $score\n";
is( $score, 60, 'count across feature' );

# collect pcount score
$score = $row->get_score(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	'method' => 'pcount',
);

# print "pcount score is $score\n";
is( $score, 60, 'pcount across feature' );

# collect ncount score
$score = $row->get_score(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	'method' => 'ncount',
);

# print "ncount score is $score\n";
is( $score, 60, 'ncount across feature' );

# collect position scores
my %score1 = $row->get_region_position_scores(
	ddb     => $ddb,
	dataset => 'data:sample2',
);

# print "position_score is \n" . print_hash(\%score1);
is( min( keys %score1 ), 58,   'min position in positioned score hash' );
is( max( keys %score1 ), 2342, 'max position in positioned score hash' );
is( $score1{686},        0.26, 'score at position 686 in positioned score hash' );

%score1 = $row->get_region_position_scores(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	'method' => 'ncount',
);

# print "position_score ncount is \n" . print_hash(\%score1);
is( min( keys %score1 ), 58,   'min position in positioned ncount hash' );
is( max( keys %score1 ), 2342, 'max position in positioned ncount hash' );
is( $score1{686},        1,    'score at position 686 in positioned ncount hash' );

# move to next gene to get better relative scores
$row = $stream->next_row;
isa_ok( $row, 'Bio::ToolBox::Data::Feature', 'row Feature' );
is( $row->name, 'YAL044C', 'Feature name' );

# collect 5' relative position scores
%score1 = $row->get_relative_point_position_scores(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	extend   => 500,
	position => 5,
);
is( min( keys %score1 ), -467, 'min 5prime relative position in positioned score hash' );
is( max( keys %score1 ), 467,  'max 5prime relative position in positioned score hash' );
is( $score1{218}, 0.62,
	'score at position 218 in 5prime relative positioned score hash' );

# collect 3' relative position scores
%score1 = $row->get_relative_point_position_scores(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	extend   => 500,
	position => 3,
);

# print "3' position_score is \n" . print_hash(\%score1);
is( min( keys %score1 ), -430, 'min 3prime relative position in positioned score hash' );
is( max( keys %score1 ), 481,  'max 3prime relative position in positioned score hash' );
is( $score1{187}, 0.75,
	'score at position 187 in 3prime relative positioned score hash' );

# collect region extended relative position scores
%score1 = $row->get_region_position_scores(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	extend   => 500,
	absolute => 1,
);

# print "extended position_score is \n" . print_hash(\%score1);
is( min( keys %score1 ),
	57469, 'min extended relative position in positioned score hash' );
is( max( keys %score1 ),
	58929, 'max extended relative position in positioned score hash' );
is( $score1{58532}, 0.34,
	'score at position 58532 in extended relative positioned score hash' );
is( scalar( keys %score1 ), 39,
	'number of keys extended relative positioned score hash' );

# collect avoided extended relative position scores
%score1 = $row->get_region_position_scores(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	extend   => 1000,
	absolute => 1,
	avoid    => 1,
);

# print "avoided extended position_score is \n" . print_hash(\%score1);
is( min( keys %score1 ),
	56969, 'min avoided extended relative position in positioned score hash' );
is( max( keys %score1 ),
	58683, 'max avoided extended relative position in positioned score hash' );
is( $score1{58532}, 0.34,
	'score at position 58532 in avoided extended relative positioned score hash' );
is( scalar( keys %score1 ),
	26, 'number of keys avoided extended relative positioned score hash' );

# collect absolute position scores
%score1 = $row->get_relative_point_position_scores(
	ddb      => $ddb,
	dataset  => 'data:sample2',
	extend   => 500,
	position => 4,
	'method' => 'count',
);

# print "position_score is \n" . print_hash(\%score1);
is( $score1{-6}, 1, 'count at relative position -5 in positioned score hash' );

### Repeat score collection with advanced annotation
undef(%score1);
undef($score);
undef($feature);
undef($row);
undef($stream);
undef($Data);
$Data = Bio::ToolBox::Data->new(
	'db'      => $ddb,
	'feature' => 'mRNA:tim',

	# we need to work with transcripts here, not genes
);
isa_ok( $Data, 'Bio::ToolBox::Data', 'db collected gene table' );

# check metadata
is( $Data->feature,      'mRNA:tim', 'feature' );
is( $Data->feature_type, 'named',    'feature type' );

# check data table
is( $Data->number_columns, 3, 'number of columns' );
is( $Data->last_row,       4, 'last row' );

# columns
is( $Data->id_column,   1, 'primary ID column' );
is( $Data->name_column, 2, 'identify name column' );
is( $Data->type_column, 3, 'identify type column' );

# check first feature
$Data->sort_data( 2, 'i' );
$stream = $Data->row_stream;
isa_ok( $stream, 'Bio::ToolBox::Data::Iterator', 'row stream iterator' );
$row = $stream->next_row;
isa_ok( $row, 'Bio::ToolBox::Data::Feature', 'row Feature' );
is( $row->type, 'mRNA:tim', 'Feature type' );
is( $row->name, 'YAL030W',  'Feature name' );

# test SeqFeature
$feature = $row->feature;
isa_ok( $feature, 'Bio::DB::SeqFeature', 'db SeqFeature from row Feature' );

# test subfeature scores, should use same database
$score = $row->get_score(
	dataset  => 'data',
	'method' => 'count',
);
is( $score, 13, 'count score across mRNA' );
$score = $row->get_score(
	dataset    => 'data',
	'method'   => 'count',
	subfeature => 'exon',
);
is( $score, 9, 'count score across subfeature exons' );
$score = $row->get_score(
	dataset    => 'data',
	'method'   => 'count',
	subfeature => 'cds',
);
is( $score, 8, 'count score across subfeature CDS' );
$score = $row->get_score(
	dataset    => 'data',
	'method'   => 'pcount',
	subfeature => '5p_utr',
);

# we do a precise count here because a quirk of BDBSFS adapter returns EVERYTHING found
# when nothing fits in the coordinate search window. bug or feature? ugh.
is( $score, 0, 'pcount score across subfeature 5p UTR' );
$score = $row->get_score(
	dataset    => 'data:sample2',
	'method'   => 'pcount',
	subfeature => '3p_utr',
);

# another bug? should return 1 with method count, but it returns all????
is( $score, 1, 'pcount score across subfeature 3p UTR' );

%score1 = $row->get_region_position_scores(
	dataset    => 'data:sample2',
	subfeature => 'CDS',
);

# print "CDS position_score is \n" . print_hash(\%score1);
is( $score1{3},   0.482, 'score at position 3 in CDS region position scores' );
is( $score1{161}, 0.392, 'score at position 161 in CDS region position scores' );

# test next feature
$row = $stream->next_row;
is( $row->name, 'YAL043C', 'Feature name' );
$score = $row->get_score(
	dataset    => 'data',
	'method'   => 'count',
	subfeature => 'exon',
);
is( $score, 61, 'count score across exon' );
$score = $row->get_score(
	dataset    => 'data',
	'method'   => 'count',
	subfeature => 'cds',
);
is( $score, 60, 'count score across multi-subfeature CDS' );

%score1 = $row->get_region_position_scores(
	dataset    => 'data:sample2',
	subfeature => 'cds',
);

# print "CDS position_score is \n" . print_hash(\%score1);
is( $score1{322},  0.461, 'score at position 321 in CDS region position scores' );
is( $score1{2314}, 0.262, 'score at position 2313 in CDS region position scores' );

### Extra test for checking non-conforming chromosome name
# using sample.bed file as in other tests
undef $Data;
undef $stream;
$Data =
	Bio::ToolBox::Data->new( file => File::Spec->catfile( $Bin, "Data", "sample.bed" ) );
isa_ok( $Data, 'Bio::ToolBox::Data', 'BED Data' );
$stream = $Data->row_stream;
$row    = $stream->next_row;
$score  = $row->get_score(
	db       => $ddb,
	dataset  => 'data',
	'method' => 'mean',
);
is( sprintf( "%.2f", $score ), -0.14, 'mean score across first bed feature' );
$row   = $stream->next_row;
$score = $row->get_score(
	db       => $ddb,
	dataset  => 'data',
	'method' => 'median',
);
is( sprintf( "%.1f", $score ),
	0.5, 'median score across second bed feature with alt chromo' );

### Test for sequence retrieval using Bio::DB::Fasta
# reusing the Data object immediately above
$Data->bam_adapter('none');    # to ensure we use BioPerl's Bio::DB::Fasta
my $fasta = File::Spec->catfile( $Bin, "Data", "sequence.fa" );
print "opening $fasta...\n";
my $fdb = Bio::ToolBox::Data->open_database($fasta);
isa_ok( $fdb, 'Bio::DB::Fasta', 'BioPerl fasta adapter database' );

# reuse row object from above, but we will specify our own coordinates since our
# fasta file is too small and doesn't cover the bed coordinates
my $sequence = $row->get_sequence(
	db     => $fdb,
	seq_id => 'chrI',
	start  => 257,
	end    => 275,
	strand => 1,
);
is( $sequence, 'ACCCTACCATTACCCTACC', 'fetched fasta sequence' );
undef $fdb;
unlink "$fasta.index";

sub print_hash {
	my $hash = shift;
	foreach ( sort { $a <=> $b } keys %{$hash} ) {
		print "    $_\t=> " . $hash->{$_} . "\n";
	}
}
