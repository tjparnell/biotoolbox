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
	if (eval {require Bio::DB::SeqFeature::Store::memory; 1}) {
		plan tests => 60;
	}
	else {
		plan skip_all => 'Bio::DB::SeqFeature::Store::memory or DB_File not available';
	}
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile($Bin, "Data", "biotoolbox.cfg");
}

require_ok 'Bio::ToolBox::Data' or 
	BAIL_OUT "Cannot load Bio::ToolBox::Data";


### Collect features from a database
my $infile = File::Spec->catfile($Bin, "Data", "chrI.gff3");
my $Data = Bio::ToolBox::Data->new(
	'db'      => $infile,
	'feature' => 'gene:SGD',
);
isa_ok($Data, 'Bio::ToolBox::Data', 'db collected gene table');

# check metadata
is($Data->feature, 'gene:SGD', 'feature');
is($Data->feature_type, 'named', 'feature type');
is($Data->program, undef, 'program');
is($Data->gff, 0, 'gff format');
is($Data->bed, 0, 'bed format');
is($Data->filename, undef, 'filename');

# check data table
is($Data->number_columns, 3, 'number of columns');
is($Data->last_row, 26, 'last row');

# columns
is($Data->id_column, 0, 'primary ID column');
is($Data->name_column, 1, 'identify name column');
is($Data->type_column, 2, 'identify type column');

# test database
is($Data->database, $infile, 'database name');
my $db = $Data->open_database;
isa_ok($db, 'Bio::DB::SeqFeature::Store', 'opened database');

# since we are dealing with a memory database, the features returned 
# are not in a predictable order suitable for testing
# so we will sort the data table by increasing name
$Data->sort_data(1, 'i');

# test row_stream
my $stream = $Data->row_stream;
isa_ok($stream, 'Bio::ToolBox::Data::Iterator', 'row stream iterator');

# look at first row
my $row = $stream->next_row;
isa_ok($row, 'Bio::ToolBox::Data::Feature', 'row Feature');
ok($row->row_index, 'row index'); # could be anything
is($row->type, 'gene:SGD', 'Feature type');
is($row->name, 'YAL043C', 'Feature name');
cmp_ok($row->id, '>=', 1, 'primary id'); # in memory db, could be anything?
is($row->start, undef, 'Feature start, undefined'); # not in table, returns nothing
is($row->stop, undef, 'Feature stop, undefined'); # not in table, returns nothing

# test SeqFeature
my $feature = $row->feature;
isa_ok($feature, 'Bio::DB::SeqFeature', 'db SeqFeature from row Feature');
is($row->start, 58695, 'Feature start, defined'); # uses db feature for start
is($row->stop, 61052, 'Feature stop, defined'); # uses db feature for stop
my $segment = $row->segment;
isa_ok($segment, 'Bio::DB::SeqFeature::Segment', 'db segment from row Feature');

# verify data from a data database
my $ddb = File::Spec->catfile($Bin, "Data", "sample2.gff3");
my $verified = $Data->verify_dataset('data:sample2', $ddb);
is($verified, 'data:sample2', 'dataset in data database verified');

# collect mean score
my $score = $row->get_score(
	ddb     => $ddb,
	dataset => 'data',
	'method'  => 'mean',
);
# print "mean score with data and ddb is $score\n";
is(sprintf("%.2f", $score), 0.18, 'mean score across feature');
 
# collect max score
$score = $row->get_score(
	ddb     => $ddb,
	dataset => 'data:sample2',
	'method'  => 'max',
);
# print "max score is $score\n";
is(sprintf("%.2f", $score), 0.46, 'max score across feature');

# collect median score
$score = $row->get_score(
	ddb     => $ddb,
	dataset => 'data:sample2',
	'method'  => 'median',
);
# print "median score is $score\n";
is(sprintf("%.2f", $score), 0.19, 'median score across feature');

# collect 5' mean score
$score = $row->get_score(
	ddb     => $ddb,
	dataset => 'data:sample2',
	'method'  => 'mean',
	start   => $row->stop - 1000,
	stop    => $row->stop,
);
# print "5 prime mean score is $score\n";
is(sprintf("%.2f", $score), '0.30', 'mean score across 5 prime feature');

# collect sum count score
$score = $row->get_score(
	ddb     => $ddb,
	dataset => 'data:sample2',
	'method'  => 'sum',
	value   => 'count',
);
# print "median score is $score\n";
is($score, 60, 'sum count across feature');

# collect mean length score
$score = $row->get_score(
	ddb     => $ddb,
	dataset => 'data:sample2',
	'method'  => 'mean',
	value   => 'length',
);
# print "median score is $score\n";
is($score, 1, 'mean length of scores across feature');

# collect position scores
my %score1 = $row->get_position_scores(
	ddb     => $ddb,
	dataset => 'data:sample2',
);
# print "position_score is \n" . print_hash(\%score1);
is(min(keys %score1), 58, 'min position in positioned score hash');
is(max(keys %score1), 2342, 'max position in positioned score hash');
is($score1{686}, 0.26, 'score at position 686 in positioned score hash');



# move to next gene to get better relative scores
$row = $stream->next_row;
isa_ok($row, 'Bio::ToolBox::Data::Feature', 'row Feature');
is($row->name, 'YAL044C', 'Feature name');

# collect 5' relative position scores
%score1 = $row->get_position_scores(
	ddb     => $ddb,
	dataset => 'data:sample2',
	start   => -500,
	stop    => 500,
	position => 5,
);
# print "5' position_score is \n" . print_hash(\%score1);
is(min(keys %score1), -467, 'min 5prime relative position in positioned score hash');
is(max(keys %score1), 467, 'max 5prime relative position in positioned score hash');
is($score1{218}, 0.62, 'score at position 218 in 5prime relative positioned score hash');

# collect 3' relative position scores
%score1 = $row->get_position_scores(
	ddb     => $ddb,
	dataset => 'data:sample2',
	start   => -500,
	stop    => 500,
	position => 3,
);
# print "3' position_score is \n" . print_hash(\%score1);
is(min(keys %score1), -430, 'min 3prime relative position in positioned score hash');
is(max(keys %score1), 481, 'max 3prime relative position in positioned score hash');
is($score1{187}, 0.75, 'score at position 187 in 3prime relative positioned score hash');

# collect absolute position scores
%score1 = $row->get_position_scores(
	ddb     => $ddb,
	dataset => 'data:sample2',
	start   => -500,
	stop    => 500,
	position => 5,
	absolute => 1,
);
# print "position_score is \n" . print_hash(\%score1);
is(min(keys %score1), 57995, 'absolute min position in positioned score hash');
is(max(keys %score1), 58929, 'absolute max position in positioned score hash');
is($score1{58212}, 0.52, 'score at absolute position 58212 in positioned score hash');

# collect extended relative position scores
%score1 = $row->get_position_scores(
	ddb     => $ddb,
	dataset => 'data:sample2',
	extend  => 500,
	absolute => 1,
);
# print "extended position_score is \n" . print_hash(\%score1);
is(min(keys %score1), 57469, 'min extended relative position in positioned score hash');
is(max(keys %score1), 58929, 'max extended relative position in positioned score hash');
is($score1{58532}, 0.34, 'score at position 58532 in extended relative positioned score hash');
is(scalar(keys %score1), 39, 'number of keys extended relative positioned score hash');

# collect extended relative position scores
%score1 = $row->get_position_scores(
	ddb     => $ddb,
	dataset => 'data:sample2',
	extend  => 1000,
	absolute => 1,
	avoid   => 1,
);
# print "avoided extended position_score is \n" . print_hash(\%score1);
is(min(keys %score1), 56969, 'min avoided extended relative position in positioned score hash');
is(max(keys %score1), 58683, 'max avoided extended relative position in positioned score hash');
is($score1{58532}, 0.34, 'score at position 58532 in avoided extended relative positioned score hash');
is(scalar(keys %score1), 26, 'number of keys avoided extended relative positioned score hash');

# collect absolute position scores
%score1 = $row->get_position_scores(
	ddb     => $ddb,
	dataset => 'data:sample2',
	start   => -500,
	stop    => 500,
	position => 4,
	value   => 'count',
);
# print "position_score is \n" . print_hash(\%score1);
is($score1{-5}, 1, 'count at relative position -5 in positioned score hash');



# extra test for checking non-conforming chromosome name
# using sample.bed file as in other tests
undef $Data;
undef $stream;
$Data = Bio::ToolBox::Data->new(file => File::Spec->catfile($Bin, "Data", "sample.bed"));
isa_ok($Data, 'Bio::ToolBox::Data', 'BED Data');
$stream = $Data->row_stream;
$row = $stream->next_row;
$score = $row->get_score(
	db      => $ddb,
	dataset => 'data',
	'method'  => 'mean',
);
is(sprintf("%.2f", $score), -0.14, 'mean score across first bed feature');
$row = $stream->next_row;
$score = $row->get_score(
	db      => $ddb,
	dataset => 'data',
	'method'  => 'median',
);
is(sprintf("%.2f", $score), 0.49, 'mean score across second bed feature with alt chromo');





sub print_hash {
	my $hash = shift;
	foreach (sort {$a <=> $b} keys %$hash) {
		print "    $_\t=> " . $hash->{$_} . "\n";
	}
}