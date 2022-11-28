#!/usr/bin/perl -w

# Test script for Bio::ToolBox::utility

use strict;
use Test::More;

BEGIN {
	plan tests => 19;
}

use_ok(
	'Bio::ToolBox::utility', qw(
		parse_list
		format_with_commas
		ask_user_for_index
		simplify_dataset_name
		sane_chromo_sort
	)
) or BAIL_OUT "Cannot load Bio::ToolBox::utility";

### parse list
my $in     = '1-10';
my @expect = ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );
my @got    = parse_list($in);
is_deeply( \@got, \@expect, "Parse list $in" );

$in     = '1,2,6,8-10';
@expect = ( 1, 2, 6, 8, 9, 10 );
@got    = parse_list($in);
is_deeply( \@got, \@expect, "Parse list $in" );

### format_with_commas
$in = 1234567890;
my $expect = '1,234,567,890';
my $got    = format_with_commas($in);
is( $got, $expect, 'Format number 1,234,567,890' );

$in     = $in / 100;
$expect = '12,345,678.9';
$got    = format_with_commas($in);
is( $got, $expect, 'Format number 12,345,678.9' );

$in /= 10;
$in *= -1;
$expect = '-1,234,567.89';
$got    = format_with_commas($in);
is( $got, $expect, 'Format number -1,234,567.89' );

### ask_user_for_index

# this requires user interaction
# I'm sure there are Test modules that fake this
# not worth it at this time

### simplify_dataset_name
$in     = 'file:test1.bw';
$expect = 'test1';
$got    = simplify_dataset_name($in);
is( $got, $expect, "Simplify dataset name $in" );

$in     = 'file:/my/data/path/test1.bw';
$expect = 'test1';
$got    = simplify_dataset_name($in);
is( $got, $expect, "Simplify dataset name $in" );

$in     = 'file:/my/data/path/test1.rpm.bw';
$expect = 'test1';
$got    = simplify_dataset_name($in);
is( $got, $expect, "Simplify dataset name $in" );

$in     = 'file:/my/data/path/test1_log2FE.bw';
$expect = 'test1';
$got    = simplify_dataset_name($in);
is( $got, $expect, "Simplify dataset name $in" );

$in     = 'file:/my/data/path/test1-dup.bw';
$expect = 'test1';
$got    = simplify_dataset_name($in);
is( $got, $expect, "Simplify dataset name $in" );

$in     = 'https://my.server.com/my/data/path/test1.rmdup.bam';
$expect = 'test1';
$got    = simplify_dataset_name($in);
is( $got, $expect, "Simplify dataset name $in" );

$in     = '/my/data/path/test1.dup.extend150.rpm.bw';
$expect = 'test1';
$got    = simplify_dataset_name($in);
is( $got, $expect, "Simplify dataset name $in" );

$in     = '/my/data/path/test1.dup.extend150.rpm.bw';
$expect = 'test1';
$got    = simplify_dataset_name($in);
is( $got, $expect, "Simplify dataset name $in" );

$in     = 'test1&test2';
$expect = 'test1&test2';
$got    = simplify_dataset_name($in);
is( $got, $expect, "Simplify dataset name $in" );

$in     = 'path/test1_f.bw&path/test2_r.bw';
$expect = 'test1_f&test2_r';
$got    = simplify_dataset_name($in);
is( $got, $expect, "Simplify dataset name $in" );

### sane_chromo_sort
my @list     = qw(chr1 chr10 chr11 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 MT);
my @expect_l = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 MT);
my @got_l    = sane_chromo_sort(@list);
is_deeply( \@got_l, \@expect_l, 'Sanely sort chromosome list 1' );

@list = qw(1 10 chr11 chr2 chr3 chr4 5 chr6 chr7 chr8 chr9 contig.124 contig.14 Y X mito);
@expect_l =
	qw(1 chr2 chr3 chr4 5 chr6 chr7 chr8 chr9 10 chr11 X Y mito contig.14 contig.124);
@got_l = sane_chromo_sort(@list);
is_deeply( \@got_l, \@expect_l, 'Sanely sort chromosome list 2' );

@list     = qw(2-micron scaffold-1 M scaffold-10 scaffold-3 I III IX VIII IV XII X Y );
@expect_l = qw(I III IV VIII IX X XII Y M scaffold-1 scaffold-3 scaffold-10 2-micron);
@got_l    = sane_chromo_sort(@list);
is_deeply( \@got_l, \@expect_l, 'Sanely sort chromosome list 3' );

