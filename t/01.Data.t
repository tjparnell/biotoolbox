#!/usr/bin/perl -w

# Test script for Bio::ToolBox::Data

use strict;
use Test;
use FindBin '$Bin';

BEGIN {
	plan tests => 47;
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

# test column names
my @column_names = $Data->list_columns;
ok(scalar @column_names, 9);
ok($column_names[7], 'Phase');
ok($Data->name(7), 'Phase');

# metadata
ok($Data->metadata(1, 'name'), 'Source');
$Data->metadata(2, 'accuracy', 'bogus');
my $md = $Data->metadata(2);
ok($md);
ok(ref $md, 'HASH');
ok($md->{'accuracy'}, 'bogus');

# column values
my $cv = $Data->column_values(3);
ok($cv);
ok(scalar @$cv, 310);
ok($cv->[1], 1);

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
ok($Data->value(1,8), undef);
$Data->delete_column(8);

# add the names as a new column
my @names = qw(Name);
for (my $i = 1; $i <= $Data->last_row; $i++) {
	push @names, "Feature$i";
}
my $index = $Data->add_column(\@names);
ok($index, 8);
ok($Data->value(308,8), 'Feature308');

# copy a column
$index = $Data->copy_column(8);
ok($index, 9);
ok($Data->number_columns, 10);
ok($Data->metadata(9, 'name'), 'Name');

# change and copy metadata
$Data->name(9, 'DuplicateName');
$Data->copy_metadata(2,9);
$md = $Data->metadata(9);
ok($md->{name}, 'DuplicateName');
ok($md->{accuracy}, 'bogus');
ok($md->{'index'}, 9);

# test reorder_column
$Data->reorder_column(0,3,4,8);
ok($Data->number_columns, 4);
ok($Data->value(308,3), 'Feature308');

# test iterate function
my $offset = 1;
my $start_i = $Data->start_column;
ok($start_i, 1);
my $iterate_success = $Data->iterate( sub {
	my $row = shift;
	my $new_start = $row->start - $offset;
	$row->value($start_i, $new_start);
} );
ok($iterate_success);
ok($Data->value(1, $start_i), 0);

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
