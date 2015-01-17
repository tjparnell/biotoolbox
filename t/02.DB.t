#!/usr/bin/perl -w

# Test script for Bio::ToolBox::Data
# working with a Bio::DB::SeqFeature::Store database
# using in memory database for simplicity here....

use strict;
use Test::More;
use FindBin '$Bin';
use Data::Dumper;

BEGIN {
	plan tests => 33;
	$ENV{'BIOTOOLBOX'} = "$Bin/Data/biotoolbox.cfg";
}

use lib "$Bin/../lib";
require_ok 'Bio::ToolBox::Data' or 
	BAIL_OUT "Cannot load Bio::ToolBox::Data";


### Collect features from a database
my $Data = Bio::ToolBox::Data->new(
	'db'      => "$Bin/Data/chrI.gff3.gz",
	'feature' => 'gene:SGD',
);
isa_ok($Data, 'Bio::ToolBox::Data', 'db collected gene table');

# check metadata
is($Data->feature, 'gene:SGD', 'feature');
is($Data->feature_type, 'named', 'feature type');
like($Data->program, qr/02\.DB\.t/, 'program');
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
is($Data->database, "$Bin/Data/chrI.gff3.gz", 'database name');
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
my $ddb = "$Bin/Data/sample2.gff3";
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
 
# collect mean score
$score = $row->get_score(
	db      => $ddb,
	dataset => 'data',
	'method'  => 'mean',
);
# print "mean score is $score\n";
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

# collect position scores
my %score1 = $row->get_position_scores(
	ddb     => $ddb,
	dataset => 'data:sample2',
);
# print "position_score is ", Dumper(\%score1), "\n";



END {
# 	unlink "$Bin/Data/DB_results.txt";
}
