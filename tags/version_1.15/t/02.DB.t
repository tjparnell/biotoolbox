#!/usr/bin/env perl

# Test script for Bio::ToolBox::Data
# working with a Bio::DB::SeqFeature::Store database
# using in memory database for simplicity here....

use strict;
use Test;
use FindBin '$Bin';

BEGIN {
	plan tests => 15;
	$ENV{'BIOTOOLBOX'} = "$Bin/Data/biotoolbox.cfg";
}

use lib "$Bin/../lib";
use Bio::ToolBox::Data;

### Open a test file
my $Data = Bio::ToolBox::Data->new(
	'db'      => "$Bin/Data/chrI.gff3.gz",
	'feature' => 'gene:SGD',
);
ok($Data);

# test number_columns
my $col_number = $Data->number_columns;
ok($col_number, 3);

# test last_row
my $last_row = $Data->last_row;
ok($last_row, 92);

# test column_name
ok($Data->name_column, 1);

# test column_type
ok($Data->type_column, 2);

# test database
ok($Data->database, "$Bin/Data/chrI.gff3.gz");
my $db = $Data->open_database;
ok($db);

# test add_column
my $cds_index = $Data->add_column('CDS_count');
ok($cds_index, 3);

# test row_stream
my $stream = $Data->row_stream;
ok($stream);

# first row feature
my $row = $stream->next_row;
ok($row);

# test row feature
ok($row->type, 'gene:SGD');

# test SeqFeature
my $feature = $row->feature;
ok($feature);
# ok($feature->start, 87286);
# ok($feature->display_name, 'YAL030W');

# test get_score CDS count
my $score = $row->get_score(
	'db'      => $db,
	'dataset' => 'cds',
	'method'  => 'sum',
	'value'   => 'count',
);
ok($score >= 1);
$row->value($cds_index, $score);

# add remaining values
while ($row = $stream->next_row) {
	my $score = $row->get_score(
		'db'      => $db,
		'dataset' => 'cds',
		'method'  => 'sum',
		'value'   => 'count',
	);
	$row->value($cds_index, $score);
}

# test save file
my $file = $Data->save(filename => "$Bin/Data/DB_results");
ok($file, "$Bin/Data/DB_results.txt");
ok(-e "$Bin/Data/DB_results.txt");



END {
	unlink "$Bin/Data/DB_results.txt";
}
