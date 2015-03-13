#!/usr/bin/env perl

# Test script for Bio::ToolBox::Data

use strict;
use Test;
use FindBin '$Bin';

BEGIN {
	plan tests => 26;
	$ENV{'BIOTOOLBOX'} = "$Bin/Data/biotoolbox.cfg";
}

use lib "$Bin/../lib";
use Bio::ToolBox::Data;

### Open a test file
my $Data = Bio::ToolBox::Data->new(
	file => "$Bin/Data/chrI.gff3.gz",
);
ok($Data);

# test gff
ok($Data->gff, 3);

# test number_columns
my $col_number = $Data->number_columns;
ok($col_number, 9);

# test last_row
my $last_row = $Data->last_row;
ok($last_row, 309);

# test find_column
my $group_index = $Data->find_column('Group');
ok($group_index, 8);

# test column chromosome
my $chromo_index = $Data->chromo_column;
ok($chromo_index, 0);

# test row_stream
my $stream = $Data->row_stream;
ok($stream);

# first row feature
my $row = $stream->next_row;
ok($row);
ok($row->value($chromo_index), 'chrI');
ok($row->start, 1);
ok($row->end, 230218);

# second row feature
$row = $stream->next_row;
ok($row->value(2), 'repeat_region');
ok($row->end, 62);

# change value
$row->value(4, 100);

# check the changed value
ok($row->end, 100);
ok($Data->value(2,4), 100);

# test delete row
$Data->delete_row(1);
ok($Data->last_row, 308);
ok($Data->value(1,4), 100);

# test delete_column
$Data->delete_column(8);
ok($Data->number_columns, 8);

# test add_column
$Data->add_column('Name');
ok($Data->number_columns, 9);

# add the names
$stream = $Data->row_stream;
my $i = 1;
while ($row = $stream->next_row) {
	$row->value(8, "Feature$i");
	$i++;
}
ok($Data->value(308,8), 'Feature308');

# test reorder_column
$Data->reorder_column(0,3,4,8);
ok($Data->number_columns, 4);
ok($Data->value(308,3), 'Feature308');

# test splice function
$Data->splice_data(2,4); # second quarter of the data table
ok($Data->last_row, 77);
ok($Data->value(77,3), 'Feature154');

# test save file
my $file = $Data->save(filename => "$Bin/Data/chrI.bed");
ok($file, "$Bin/Data/chrI.bed");
ok(-e "$Bin/Data/chrI.bed");

END {
	unlink "$Bin/Data/chrI.bed";
}
