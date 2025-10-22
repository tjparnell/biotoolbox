#!/usr/bin/env perl

# Test script for Bio::ToolBox scripts

use Test2::V0 -no_srand => 1;
use Test::Script;
plan(22);
use English qw(-no_match_vars);
use File::Spec;
use FindBin '$Bin';

BEGIN {
	## no critic
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, 'Data', 'biotoolbox.cfg' );
	## use critic
}

## check that all scripts compile correctly
script_compiles( 'scripts/bam2wig.pl', 'bam2wig compiles ok' );
script_compiles(
	'scripts/correlate_position_data.pl',
	'correlate_position_data compiles ok'
);
script_compiles( 'scripts/data2bed.pl',         'data2bed compiles ok' );
script_compiles( 'scripts/data2fasta.pl',       'data2fasta compiles ok' );
script_compiles( 'scripts/data2gff.pl',         'data2gff compiles ok' );
script_compiles( 'scripts/data2wig.pl',         'data2wig compiles ok' );
script_compiles( 'scripts/db_types.pl',         'db_types compiles ok' );
script_compiles( 'scripts/get_binned_data.pl',  'get_binned_data compiles ok' );
script_compiles( 'scripts/get_datasets.pl',     'get_datasets compiles ok' );
script_compiles( 'scripts/get_feature_info.pl', 'get_feature_info compiles ok' );
script_compiles( 'scripts/get_features.pl',     'get_features compiles ok' );
script_compiles( 'scripts/get_gene_regions.pl', 'get_gene_regions compiles ok' );
script_compiles(
	'scripts/get_intersecting_features.pl',
	'get_intersecting_features compiles ok'
);
script_compiles( 'scripts/get_relative_data.pl',   'get_relative_data compiles ok' );
script_compiles( 'scripts/join_data_file.pl',      'join_data_file compiles ok' );
script_compiles( 'scripts/manipulate_datasets.pl', 'manipulate_datasets compiles ok' );
script_compiles( 'scripts/manipulate_wig.pl',      'manipulate_wig compiles ok' );
script_compiles( 'scripts/merge_datasets.pl',      'merge_datasets compiles ok' );
script_compiles( 'scripts/pull_features.pl',       'pull_features compiles ok' );
script_compiles( 'scripts/split_data_file.pl',     'split_data_file compiles ok' );
script_compiles( 'scripts/ucsc_table2gff3.pl',     'ucsc_table2gff3 compiles ok' );

SKIP: {
	eval { require Bio::DB::SeqFeature::Store };
	skip( 'Bio::DB::SeqFeature::Store not installed', 1 ) if $EVAL_ERROR;
	script_compiles( 'scripts/db_setup.pl', 'db_setup compiles ok' );
}
