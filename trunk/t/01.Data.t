#!/usr/bin/perl -w

# Test script for Bio::ToolBox::Data and Bio::ToolBox::Data::Stream

use strict;
use Test::More;
use FindBin '$Bin';

BEGIN {
	plan tests => 104;
	$ENV{'BIOTOOLBOX'} = "$Bin/Data/biotoolbox.cfg";
}

use lib "$Bin/../lib";
require_ok 'Bio::ToolBox::Data' or 
	BAIL_OUT "Cannot load Bio::ToolBox::Data";
require_ok 'Bio::ToolBox::Data::Stream' or 
	BAIL_OUT "Cannot load Bio::ToolBox::Data::Stream";

### Open a test file
my $Data = Bio::ToolBox::Data->new(
	file => "$Bin/Data/chrI.gff3.gz",
);
isa_ok($Data, 'Bio::ToolBox::Data', 'GFF3 Data');

# test general metadata
is($Data->gff, 3, 'gff version');
is($Data->bed, 0, 'bed version');
is($Data->program, '', 'program name');
is($Data->feature, 'region', 'general feature');
is($Data->feature_type, 'coordinate', 'feature type');
is($Data->database, '', 'database');
is($Data->filename, "$Bin/Data/chrI.gff3.gz", 'filename');
is($Data->basename, 'chrI', 'basename');
is($Data->extension, '.gff3.gz', 'extension');
is($Data->path, "$Bin/Data/", 'path');

# test comments
my @comments = $Data->comments;
is(scalar @comments, 3, 'comments array');
is($comments[0], '# date Tue Feb  8 19:50:12 2011', 'first comment');
$Data->delete_comment(0);
@comments = $Data->comments;
is(scalar @comments, 2, 'comments array after delete');
$Data->add_comment('this is a comment');
@comments = $Data->comments;
is(scalar @comments, 3, 'comments array after adding');
is($comments[2], 'this is a comment', 'added comment');

# test number_columns
is($Data->number_columns, 9, 'number of columns');

# test last_row
is($Data->last_row, 79, 'last row index');

# test columns
is($Data->chromo_column, 0, 'chromosome column');
is($Data->start_column, 3, 'start column');
is($Data->stop_column, 4, 'stop column');
is($Data->strand_column, 6, 'strand column');
is($Data->type_column, 2, 'type column');

# test find_column
is($Data->find_column('Group'), 8, 'find column Group');


# test column names
my @column_names = $Data->list_columns;
is(scalar @column_names, 9, 'number of column names');
is($column_names[7], 'Phase', 'name of column 7');
is($Data->name(7), 'Phase', 'name of column 7 again');

# column metadata
is($Data->metadata(1, 'name'), 'Source', 'column name via metadata value');
$Data->metadata(2, 'accuracy', 'bogus');
my $md = $Data->metadata(2);
ok($md, 'metadata success');
isa_ok($md, 'HASH', 'metadata is a hash');
is($md->{'accuracy'}, 'bogus', 'set metadata value is correct');

# column values
my $cv = $Data->column_values(3);
ok($cv, 'column values');
is(scalar @$cv, 80, 'number of column values');
is($cv->[1], 1, 'check specific column value');

# test row_stream
my $stream = $Data->row_stream;
isa_ok($stream, 'Bio::ToolBox::Data::Iterator', 'Iterator object');

# first row feature
my $row = $stream->next_row;
isa_ok($row, 'Bio::ToolBox::Data::Feature', 'Feature object');
is($row->value(0), 'chrI', 'row object value of chromo index');
is($row->start, 1, 'row object start value');
is($row->end, 230218, 'row object end value');

# second row feature
$row = $stream->next_row;
is($row->value(2), 'repeat_region', 'next row object value');
is($row->end, 62, 'row object end value');

# change value
$row->value(4, 100);

# check the changed value
is($row->end, 100, 'checked changed row object end value');
is($Data->value(2,4), 100, 'checked changed value in data table');

# test delete row
$Data->delete_row(1);
is($Data->last_row, 78, 'last row index after deleting 1 row');
is($Data->value(1,4), 100, 'data table changed value');

# test delete_column
$Data->delete_column(8);
is($Data->number_columns, 8, 'number of columns after deleting column');

# test add_column
my $added = $Data->add_column('Name');
is($added, 8, 'returned index of added column');
is($Data->number_columns, 9, 'number of columns after adding column');
is($Data->value(1,8), undef, 'check value of added column');
$Data->delete_column(8);

# add array of names as a new column
my @new_column = qw(Name);
for (my $i = 1; $i <= $Data->last_row; $i++) {
	push @new_column, "Feature$i";
}
my $index = $Data->add_column(\@new_column);
is($index, 8, 'added column index');
is($Data->value(78,8), 'Feature78', 'checked column value after adding new column values');

# copy a column
$index = $Data->copy_column(8);
is($index, 9, 'index of copied column');
is($Data->number_columns, 10, 'new number of columns');
is($Data->metadata(9, 'name'), 'Name', 'Name of new column');

# change and copy metadata
$Data->name(9, 'DuplicateName');
$Data->copy_metadata(2,9);
$md = $Data->metadata(9);
is($md->{name}, 'DuplicateName', 'metadata of changed column name');
is($md->{accuracy}, 'bogus', 'metadata of copied column');
is($md->{'index'}, 9, 'metadata of copied column');

# sort table
$Data->sort_data(8, 'd');
is($Data->value(1,8), 'Feature9', 'check first name after reverse sort');
is($Data->value(1,3), 538, 'check first start after reverse sort');
is($Data->value(78,8), 'Feature1', 'check last name after reverse sort');
is($Data->value(78,3), 1, 'check last start after reverse sort');

# genomic sort rows
$Data->gsort_data;
is($Data->value(1,8), 'Feature2', 'check first name after genomic sort');
is($Data->value(1,3), 1, 'check first start after genomic sort');
is($Data->value(78,8), 'Feature77', 'check last name after genomic sort');
is($Data->value(78,3), 58695, 'check last start after genomic sort');

# test reorder_column
$Data->reorder_column(0,3,4,8);
is($Data->number_columns, 4, 'number of columns after reordering');
is($Data->value(78,3), 'Feature77', 'value in data table after reordering');

# test iterate function
my $offset = 1;
my $start_i = $Data->start_column;
is($start_i, 1, 'start column index');
my $iterate_success = $Data->iterate( sub {
	my $row = shift;
	my $new_start = $row->start - $offset;
	$row->value($start_i, $new_start);
} );
ok($iterate_success, 'iterate success');
is($Data->value(1, $start_i), 0, 'data table value after iteration');

# test splice function
$Data->splice_data(2,2); # second half of the data table
is($Data->last_row, 39, 'last row index after splicing');
is($Data->value(39,3), 'Feature77', 'data table value after splicing');

# test save file
my $file = $Data->save(filename => "$Bin/Data/chrI.bed");
is($file, "$Bin/Data/chrI.bed", 'output file name success');
ok(-e "$Bin/Data/chrI.bed", 'output file exists');





### Testing Stream object ###

# open the bed file we just wrote
my $Stream = Bio::ToolBox::Data::Stream->new(
	file    => $file,
);
isa_ok($Stream, 'Bio::ToolBox::Data::Stream', 'Stream object');
is($Stream->bed, 4, 'bed value');
is($Stream->gff, 0, 'gff value');
is($Stream->feature, 'region', 'feature');
is($Stream->feature_type, 'coordinate', 'feature_type');
is($Stream->extension, '.bed', 'extension');

# columns
is($Stream->chromo_column, 0, 'chromosome column');
is($Stream->start_column, 1, 'start column');
is($Stream->stop_column, 2, 'stop column');
is($Stream->name_column, 3, 'name column');

# check column names
@column_names = $Stream->list_columns;
is(scalar @column_names, 4, 'number of columns');
is($column_names[3], 'Name', 'name of column');
is($Stream->name(2), 'End', 'name of column again');

# iterate
my $f = $Stream->next_row;
isa_ok($f, 'Bio::ToolBox::Data::Feature', 'next row Feature object');

# check feature
is($f->seq_id, 'chrI', 'feature seq_id');
is($f->start, 35155, 'feature start position transformed');
is($f->stop, 36303, 'feature stop position');
is($f->name, 'Feature41', 'feature name');

$Stream->close_fh;
undef $Stream;

# open again differently
$Stream = Bio::ToolBox::Data->new(
	stream      => 1,
	file        => $file,
);
isa_ok($Stream, 'Bio::ToolBox::Data::Stream', 'Stream object');

# create output file
my $file1 = $file;
$file1 =~ s/\.bed$/_2.bed/;
my $outStream = $Stream->duplicate($file1);
isa_ok($outStream, 'Bio::ToolBox::Data::Stream', 'duplicated Stream object');
is($Stream->basename, 'chrI', 'in Stream basename');
is($outStream->basename, 'chrI_2', 'out Stream basename');

# duplicate file
while (my $row = $Stream->next_row) {
	# just write the same thing, no need to modify
	$outStream->write_row($row);
}
$Stream->close_fh;
$outStream->close_fh;
is(-s $file, -s $file1, "duplicate file sizes");

# duplicate file again
my $file2 = $file1;
$file2 =~ s/_2\.bed/_3.bed/;
$Stream = Bio::ToolBox::Data::Stream->new(file => $file);
$outStream = $Stream->duplicate($file2);
while (my $row = $Stream->next_row) {
	my @a = $row->row_values;
	$outStream->write_row(\@a);
}
$Stream->close_fh;
$outStream->close_fh;
is(-s $file, -s $file2, "duplicate file sizes again");

# reload the duplicate files
my $reloaded = $Data->reload_children($file, $file1, $file2);
is($reloaded, 117, 'reloaded children files');
ok($Data->verify, 'verify data structure');
is($Data->last_row, 117, 'reloaded last row index');
is($Data->number_columns, 4, 'reloaded number of columns');

END {
	# unlink($file, $outfile, $outfile2);
}
