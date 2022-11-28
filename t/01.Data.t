#!/usr/bin/perl -w

# Test script for Bio::ToolBox::Data and Bio::ToolBox::Data::Stream

use strict;
use Test::More;
use Test::Warn;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	plan tests => 228;
	## no critic
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, "Data", "biotoolbox.cfg" );
	## use critic
}

require_ok 'Bio::ToolBox::Data'
	or BAIL_OUT "Cannot load Bio::ToolBox::Data";
require_ok 'Bio::ToolBox::Data::Stream'
	or BAIL_OUT "Cannot load Bio::ToolBox::Data::Stream";

### Open a test file
my $infile = File::Spec->catfile( $Bin, "Data", "chrI.gff3" );
my $Data   = Bio::ToolBox::Data->new( file => $infile, );
isa_ok( $Data, 'Bio::ToolBox::Data', 'GFF3 Data' );

# test general metadata
is( $Data->gff,          3,            'gff version' );
is( $Data->bed,          0,            'bed version' );
is( $Data->format,       'gff3',       'GFF3 format' );
is( $Data->program,      undef,        'program name' );
is( $Data->feature,      'region',     'general feature' );
is( $Data->feature_type, 'coordinate', 'feature type' );
is( $Data->interbase,    0,            'coordinates are not interbase' );
is( $Data->database,     undef,        'database' );
is( $Data->filename,     $infile,      'filename' );
is( $Data->basename,     'chrI',       'basename' );
is( $Data->extension,    '.gff3',      'extension' );

# is($Data->path, File::Spec->catfile($Bin, "Data", ''), 'path');
# 	the path is hard to test, because the returned value includes a trailing slash
# 	but getting File::Spec to include one artificially is hard
# 	this fails on Windows

# test comments
my @comments = $Data->comments;
is( scalar @comments, 3,                                 'comments array' );
is( $comments[0],     '# date Tue Feb  8 19:50:12 2011', 'first comment' );
$Data->delete_comment(0);
@comments = $Data->comments;
is( scalar @comments, 2, 'comments array after delete' );
$Data->add_comment('this is a comment');
@comments = $Data->comments;
is( scalar @comments, 3,                   'comments array after adding' );
is( $comments[2],     'this is a comment', 'added comment' );

# test number_columns
is( $Data->number_columns, 9, 'number of columns' );

# test last_row
is( $Data->number_rows, 79, 'number of rows' );
is( $Data->last_row,    79, 'last row index' );

# test columns
is( $Data->chromo_column, 1, 'chromosome column' );
is( $Data->start_column,  4, 'start column' );
is( $Data->stop_column,   5, 'stop column' );
is( $Data->strand_column, 7, 'strand column' );
is( $Data->type_column,   3, 'type column' );

# change column indexes
is( $Data->chromo_column(2), 2, 'changed chromosome column index' );
is( $Data->start_column(2),  2, 'changed start column index' );
is( $Data->stop_column(2),   2, 'changed stop column index' );
is( $Data->strand_column(2), 2, 'changed strand column index' );
is( $Data->type_column(2),   2, 'changed type column index' );

# change them all back to normal before continuing
$Data->chromo_column(1);
$Data->start_column(4);
$Data->stop_column(5);
$Data->strand_column(7);
$Data->type_column(3);

# test find_column
is( $Data->find_column('Group'), 9, 'find column Group' );

# test column names
my @column_names = $Data->list_columns;
is( scalar @column_names, 9,        'number of column names' );
is( $column_names[7],     'Phase',  'name of column 7' );
is( $Data->name(7),       'Strand', 'name of column 7 again' );

# column metadata
is( $Data->metadata( 1, 'name' ), 'Chromosome', 'column name via metadata value' );
$Data->metadata( 2, 'accuracy', 'bogus' );
my $md = $Data->metadata(2);
ok( $md, 'metadata success' );
isa_ok( $md, 'HASH', 'metadata is a hash' );
is( $md->{'accuracy'}, 'bogus', 'set metadata value is correct' );

# column values
my $cv = $Data->column_values(4);
ok( $cv, 'column values' );
is( scalar @{$cv}, 80, 'number of column values' );
is( $cv->[1],      1,  'check specific column value' );

# test duplicate
my $Dupe = $Data->duplicate;
isa_ok( $Dupe, 'Bio::ToolBox::Data', 'Duplicated object' );
is( $Dupe->gff,            3,            'Dupe gff version' );
is( $Dupe->program,        undef,        'Dupe program name' );
is( $Dupe->feature,        'region',     'Dupe general feature' );
is( $Dupe->feature_type,   'coordinate', 'Dupe feature type' );
is( $Dupe->database,       undef,        'Dupe database' );
is( $Dupe->filename,       q(),          'Dupe filename' );
is( $Dupe->number_columns, 9,            'Dupe number of columns' );
is( $Dupe->last_row,       0,            'Dupe last row index' );

# test row_stream
my $stream = $Data->row_stream;
isa_ok( $stream, 'Bio::ToolBox::Data::Iterator', 'Iterator object' );

# first row feature
my $row = $stream->next_row;
isa_ok( $row, 'Bio::ToolBox::Data::Feature', 'Feature object' );
is( $row->value(1),   'chrI',          'row object value of index 1' );
is( $row->seq_id,     'chrI',          'row object chromosome value' );
is( $row->start,      1,               'row object start value' );
is( $row->end,        230218,          'row object end value' );
is( $row->coordinate, 'chrI:1-230218', 'row object coordinate string' );

# add row feature
my $added_row_i = $Dupe->add_row($row);
is( $added_row_i,         1,      'Dupe add_row Feature object' );
is( $Dupe->last_row,      1,      'Dupe added last row index' );
is( $Dupe->value( 1, 4 ), 1,      'Dupe row value at 4 (start)' );
is( $Dupe->value( 1, 5 ), 230218, 'Dupe row value at 5 (end)' );

# second row feature
$row = $stream->next_row;
is( $row->value(3), 'repeat_region', 'next row object value at 3 (type)' );
is( $row->end,      62,              'row object end value' );

# check gff attribute
my $gff_att = $row->gff_attributes;
isa_ok( $gff_att, 'HASH', 'row GFF attributes hash' );
is( $gff_att->{Name}, 'TEL01L-TR', 'row GFF attribute Name' );
$gff_att->{Note} = 'I hereby claim this telomeric repeat to be mine';
is( $row->rewrite_gff_attributes, 1, 'rewrite row GFF attributes' );
is(
	$row->value(9),
'ID=TEL01L-TR; Name=TEL01L-TR; Note=I%20hereby%20claim%20this%20telomeric%20repeat%20to%20be%20mine',
	'rewritten row GFF attribute'
);

# change end value
$row->value( 5, 100 );

# check the changed value
is( $row->end,            62,  'checked row object cached end value' );
is( $Data->value( 2, 5 ), 100, 'checked changed value in data table' );
is( $row->end(100),       100, 'row end value changed via high level' );

# test delete row
$Data->delete_row(1);
is( $Data->last_row,      78,  'last row index after deleting 1 row' );
is( $Data->value( 1, 5 ), 100, 'data table changed value' );

# test delete_column
$Data->delete_column(8);
is( $Data->number_columns, 8, 'number of columns after deleting column' );

# test add_column
my $added = $Data->add_column('Name');
is( $added,                9,     'returned index of added column' );
is( $Data->number_columns, 9,     'number of columns after adding column' );
is( $Data->value( 1, 9 ),  undef, 'check value of added column' );
$Data->delete_column(9);

# add array of names as a new column
my @new_column = qw(Name);
for my $i ( 1 .. $Data->last_row ) {
	push @new_column, "Feature$i";
}
my $index = $Data->add_column( \@new_column );
is( $index, 9, 'added column index' );
is( $Data->value( 78, 9 ),
	'Feature78', 'checked column value after adding new column values' );

# copy a column
$index = $Data->copy_column(9);
is( $index,                        10,     'index of copied column' );
is( $Data->number_columns,         10,     'new number of columns' );
is( $Data->metadata( 10, 'name' ), 'Name', 'Name of new column' );

# change and copy metadata
$Data->name( 10, 'DuplicateName' );
$Data->copy_metadata( 2, 10 );
$md = $Data->metadata(10);
is( $md->{name},     'DuplicateName', 'metadata of changed column name' );
is( $md->{accuracy}, 'bogus',         'metadata of copied column' );
is( $md->{'index'},  10,              'metadata of copied column' );

# sort table
$Data->sort_data( 9, 'd' );
is( $Data->value( 1,  9 ), 'Feature9', 'check first name after reverse sort' );
is( $Data->value( 1,  4 ), 538,        'check first start after reverse sort' );
is( $Data->value( 78, 9 ), 'Feature1', 'check last name after reverse sort' );
is( $Data->value( 78, 4 ), 1,          'check last start after reverse sort' );

# genomic sort rows
$Data->gsort_data;    # Data still tagged as gff, which influences gsort method
is( $Data->value( 1,  9 ), 'Feature2',  'check first name after genomic sort' );
is( $Data->value( 1,  4 ), 1,           'check first start after genomic sort' );
is( $Data->value( 78, 9 ), 'Feature77', 'check last name after genomic sort' );
is( $Data->value( 78, 4 ), 58695,       'check last start after genomic sort' );

# test reorder_column
$Data->reorder_column( 1, 4, 5, 9 );
is( $Data->number_columns, 4,           'number of columns after reordering' );
is( $Data->value( 78, 4 ), 'Feature77', 'value in data table after reordering' );

# test iterate function
my $offset = 1;
is( $Data->interbase, 0, 'interbase value of 0' );
my $start_i = $Data->start_column;
is( $start_i, 2, 'start column index' );
my $iterate_success = $Data->iterate(
	sub {
		my $row2      = shift;
		my $new_start = $row2->start - $offset;
		$row2->value( $start_i, $new_start );
	}
);
ok( $iterate_success, 'iterate success' );
is( $Data->value( 1, $start_i ), 0, 'data table value after iteration' );
is( $Data->interbase(1),         1, 'interbase value of 1' );

# test splice function
$Data->splice_data( 2, 2 );    # second half of the data table
is( $Data->last_row,       39,          'last row index after splicing' );
is( $Data->value( 39, 4 ), 'Feature77', 'data table value after splicing' );

# test save file
is( $Data->gff(0), 0, 'reset gff value to 0' );
is( $Data->bed(4), 4, 'reset bed value to 4' );
my $outfile = File::Spec->catdir( $Bin, "Data", "chrI.bed" );
my $file    = $Data->save( filename => $outfile );
is( $file, $outfile, 'output file name success' );
ok( -e $outfile, 'output file exists' );

# clean up
undef $row;
undef $stream;
undef $Data;

### reopen saved Bed file
# many of the same types of tests as above, just with bed file
$Data = Bio::ToolBox::Data->new( file => $file );

# metadata tests
isa_ok( $Data, 'Bio::ToolBox::Data', 'Bed Data' );
is( $Data->gff,            0,            'gff version' );
is( $Data->bed,            4,            'bed version' );
is( $Data->format,         'bed4',       'bed format' );
is( $Data->program,        undef,        'program name' );
is( $Data->feature,        'region',     'general feature' );
is( $Data->feature_type,   'coordinate', 'feature type' );
is( $Data->database,       undef,        'database' );
is( $Data->filename,       $file,        'filename' );
is( $Data->basename,       'chrI',       'basename' );
is( $Data->extension,      '.bed',       'extension' );
is( $Data->number_columns, 4,            'number of columns' );
is( $Data->number_rows,    39,           'number of rows' );
is( $Data->last_row,       39,           'last row index' );
is( $Data->interbase,      1,            'coordinates are interbase' );
is( $Data->chromo_column,  1,            'chromosome column' );
is( $Data->start_column,   2,            'start column' );
is( $Data->stop_column,    3,            'stop column' );
is( $Data->strand_column,  undef,        'strand column' );
is( $Data->type_column,    undef,        'type column' );

# stream tests and row feature
$stream = $Data->row_stream;
$row    = $stream->next_row;
isa_ok( $row, 'Bio::ToolBox::Data::Feature', 'Feature object' );
is( $row->value(1),   'chrI',             'row object value of chromo index' );
is( $row->start,      35155,              'row object start value' );
is( $row->end,        36303,              'row object end value' );
is( $row->name,       'Feature41',        'row object feature name' );
is( $row->type,       'region',           'row object type value' );
is( $row->coordinate, 'chrI:35154-36303', 'row object coordinate string' );

# grab row feature directly and change attributes using high API functions - v1.68
$row = $Data->get_row(25);
isa_ok( $row, 'Bio::ToolBox::Data::Feature', 'Direct Feature object' );
is( $row->value(1),        'chrI', 'Feature chromosome actual value' );
is( $row->seq_id,          'chrI', 'Feature chromosome' );
is( $row->seq_id('chrX'),  'chrX', 'Change chromosome via high level' );
is( $row->seq_id,          'chrX', 'Feature changed chromosome' );
is( $Data->value( 25, 1 ), 'chrX', 'Feature chromosome actual value' );

is( $row->value(2),        52800, 'Feature start actual value' );
is( $row->start,           52801, 'Feature start coordinate' );
is( $row->start(52901),    52901, 'Change start coordinate via high level' );
is( $row->start,           52901, 'Feature changed start coordinate' );
is( $Data->value( 25, 2 ), 52900, 'Feature start actual value' );

is( $row->value(3),        54789, 'Feature end actual value' );
is( $row->stop,            54789, 'Feature end coordinate' );
is( $row->stop(54589),     54589, 'Change end coordinate via high level' );
is( $row->end,             54589, 'Feature changed end coordinate' );
is( $Data->value( 25, 3 ), 54589, 'Feature end actual value' );

is( $row->value(4),        'Feature63', 'Feature name actual value' );
is( $row->name,            'Feature63', 'Feature name' );
is( $row->name('bob'),     'bob',       'Change name via high level' );
is( $row->name,            'bob',       'Feature changed name' );
is( $Data->value( 25, 4 ), 'bob',       'Feature name actual value' );

is( $row->value(5), '.', 'Feature actual strand value (nonexistent)' );
is( $row->strand,   0,   'Feature strand (implied)' );
warning_is(
	sub { $row->strand(1) },
	'ERROR: No Strand column to update!',
	'Attempt strand change via high level'
);
is( $row->strand,          0,     'Check attempted strand change' );
is( $Data->value( 25, 5 ), undef, 'Feature actual changed strand value' );

is( $row->type, 'region', 'Feature type (implied)' );
warning_is(
	sub { $row->type('gene') },
	'ERROR: No Type column to update!',
	'Attempt type change via high level'
);
isnt( $row->type, 'gene', 'Check attempted type change' );
is( $row->type, 'region', 'Check actual type value' );

# check calculating reference point
is( $row->calculate_reference(5), 52901, '5\' reference position' );
is( $row->calculate_reference(3), 54589, '3\' reference position' );
is( $row->calculate_reference(4), 53745, 'midpoint reference position' );
my $args = {
	position        => 4,
	practical_start => 1001,
	practical_stop  => 2000
};
is( $row->calculate_reference($args),
	1501, 'midpoint reference position of given positions' );

undef $row;
undef $stream;
undef $Data;

### Testing Stream object ###

# open the bed file we just wrote
my $Stream = Bio::ToolBox::Data::Stream->new( in => $file, );
isa_ok( $Stream, 'Bio::ToolBox::Data::Stream', 'Stream Bed object' );
is( $Stream->bed,          4,            'bed value' );
is( $Stream->gff,          0,            'gff value' );
is( $Stream->feature,      'region',     'feature' );
is( $Stream->feature_type, 'coordinate', 'feature_type' );
is( $Stream->extension,    '.bed',       'extension' );

# columns
is( $Stream->chromo_column, 1, 'chromosome column' );
is( $Stream->start_column,  2, 'start column' );
is( $Stream->stop_column,   3, 'stop column' );
is( $Stream->name_column,   4, 'name column' );

# check column names
@column_names = $Stream->list_columns;
is( scalar @column_names, 4,      'number of columns' );
is( $column_names[3],     'Name', 'name of column' );
is( $Stream->name(3),     'End',  'name of column again' );

# iterate
my $f = $Stream->next_row;
isa_ok( $f, 'Bio::ToolBox::Data::Feature', 'next row Feature object' );

# check feature
is( $f->seq_id,     'chrI',             'feature seq_id' );
is( $f->start,      35155,              'feature start position transformed' );
is( $f->stop,       36303,              'feature stop position' );
is( $f->midpoint,   35729,              'feature midpoint position' );
is( $f->peak,       35729,              'feature peak position, default to midpoint' );
is( $f->name,       'Feature41',        'feature name' );
is( $f->coordinate, 'chrI:35154-36303', 'feature coordinate string' );

$Stream->close_fh;
undef $Stream;

# open again differently
$Stream = Bio::ToolBox::Data->new(
	stream => 1,
	in     => $file,
);
isa_ok( $Stream, 'Bio::ToolBox::Data::Stream', 'Stream object' );

# create output file
my $file1 = $file;
$file1 =~ s/\.bed$/_2.bed/;
my $outStream = $Stream->duplicate($file1);
isa_ok( $outStream, 'Bio::ToolBox::Data::Stream', 'duplicated Stream object' );
is( $Stream->basename,    'chrI',   'in Stream basename' );
is( $outStream->basename, 'chrI_2', 'out Stream basename' );

# duplicate file
while ( my $row2 = $Stream->next_row ) {

	# just write the same thing, no need to modify
	# write as Feature objects
	$outStream->write_row($row2);
}
$Stream->close_fh;
$outStream->close_fh;
is( -s $file, -s $file1, "duplicate file sizes" );

# duplicate file again
# specify out stream as a new empty bed file
my $file2 = $file1;
$file2 =~ s/_2\.bed/_3.bed/;
$Stream    = Bio::ToolBox::Data::Stream->new( in  => $file );
$outStream = Bio::ToolBox::Data::Stream->new( out => $file2, bed => 4 );
while ( my $row2 = $Stream->next_row ) {

	# write as arrays
	my @a = $row2->row_values;
	$outStream->write_row( \@a );
}
$Stream->close_fh;
$outStream->close_fh;
cmp_ok( -s $file2, '<', -s $file, "smaller file size due to lack of comments" );

# reload the duplicate files
# this should effectively delete the child files
$Data = Bio::ToolBox::Data->new();
isa_ok( $Data, 'Bio::ToolBox::Data', 'new empty Data object' );
is( $Data->number_columns, 0, 'number of columns' );
is( $Data->last_row,       0, 'last row index' );

my $reloaded = $Data->reload_children( $file, $file1, $file2 );
is( $reloaded, 117, 'reloaded children files' );
ok( $Data->verify, 'verify data structure' );
is( $Data->last_row,       117, 'reloaded last row index' );
is( $Data->number_columns, 4,   'reloaded number of columns' );

### Open a narrowPeak test file
undef $Data;
undef $row;
$infile = File::Spec->catfile( $Bin, "Data", "H3K4me3.narrowPeak" );
$Data   = Bio::ToolBox::Data->new( file => $infile, );
isa_ok( $Data, 'Bio::ToolBox::Data', 'narrowPeak Data' );

# test general metadata
is( $Data->gff,                   0,             'gff version' );
is( $Data->bed,                   10,            'bed version' );
is( $Data->format,                'narrowPeak',  'narrowPeak format' );
is( $Data->feature,               'region',      'general feature' );
is( $Data->feature_type,          'coordinate',  'feature type' );
is( $Data->extension,             '.narrowPeak', 'narrowPeak extension' );
is( $Data->number_columns,        10,            'number of columns' );
is( $Data->start_column,          2,             'start column' );
is( $Data->stop_column,           3,             'stop column' );
is( $Data->find_column('pValue'), 8,             'find column pValue' );
is( $Data->find_column('peak'),   10,            'find column peak' );
$row = $Data->get_row(1);
isa_ok( $row, 'Bio::ToolBox::Data::Feature', 'first peak interval Feature object' );
is( $row->peak,     11908866, 'peak interval peak coordinate' );
is( $row->midpoint, 11909060, 'peak interval midpoint' );

### Open a gappedPeak test file
undef $Data;
$infile = File::Spec->catfile( $Bin, "Data", "H3K27ac.bed" );
$Data   = Bio::ToolBox::Data->new( file => $infile, );
isa_ok( $Data, 'Bio::ToolBox::Data', 'gappedPeak bed Data' );

# test general metadata
is( $Data->gff,                   0,            'gff version' );
is( $Data->bed,                   15,           'bed version' );
is( $Data->format,                'gappedPeak', 'gappedPeak format' );
is( $Data->feature,               'region',     'general feature' );
is( $Data->feature_type,          'coordinate', 'feature type' );
is( $Data->extension,             '.bed',       'gappedPeak bed extension' );
is( $Data->number_columns,        15,           'number of columns' );
is( $Data->start_column,          2,            'start column' );
is( $Data->stop_column,           3,            'stop column' );
is( $Data->find_column('pValue'), 14,           'find column pValue' );

