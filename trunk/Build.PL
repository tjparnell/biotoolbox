#!/usr/bin/perl

use strict;
use warnings;
use Module::Build;



### Prepare global variables for installation

# initial build options
my $options = get_base_options();

# Check for additional Modules to install
check_useq();
check_sam();
check_big();



### Build the script
my $build = Module::Build->new(%$options);
$build->create_build_script;
exit 0;


###### subroutines to assist in building ########


sub get_base_options {
	my %options = (
		build_class       => 'Module::Build',
		module_name       => 'Bio::ToolBox',
		license           => 'perl',
		dist_version_from => 'lib/Bio/ToolBox.pm',
		dist_author       => 'Timothy Parnell <parnell.tj@gmail.com>',
		dist_abstract     => 'Tools for querying and analysis of genomic data',
		configure_requires => {
			'Module::Build'       => 0,
		},
		requires          => {
			'Config::Simple'          => 4.58,
			'Statistics::Lite'        => 3.2,
			'Statistics::Descriptive' => 3.0,
			'Bio::Root::Version'      => '1.0069001',
		},
		recommends        => {
			'Algorithm::Cluster'      => 0,
			'DBD::SQLite'             => 0,
			'File::Which'             => 0,
			'GD::Graph'               => 1.44,
			'GD::Graph::smoothlines'  => 1.6,
			'Statistics::LineFit'     => 0.06,
			'Net::FTP'                => 0,
			'Parallel::ForkManager'   => 1.02,
		},
		script_files => [
			'scripts/average_gene.pl',
			'scripts/bam2gff_bed.pl',
			'scripts/bam2wig.pl',
			'scripts/bar2wig.pl',
			'scripts/big_file2gff3.pl',
			'scripts/bin_genomic_data.pl',
			'scripts/change_chr_prefix.pl',
			'scripts/change_cluster_order.pl',
			'scripts/compare_subfeature_scores.pl',
			'scripts/convert_yeast_genome_version.pl',
			'scripts/correlate_position_data.pl',
			'scripts/CpG_calculator.pl',
			'scripts/data2bed.pl',
			'scripts/data2fasta.pl',
			'scripts/data2frequency.pl',
			'scripts/data2gff.pl',
			'scripts/data2wig.pl',
			'scripts/db_setup.pl',
			'scripts/db_types.pl',
			'scripts/filter_bam.pl',
			'scripts/find_enriched_regions.pl',
			'scripts/get_actual_nuc_sizes.pl',
			'scripts/get_bam_seq_stats.pl',
			'scripts/get_datasets.pl',
			'scripts/get_ensembl_annotation.pl',
			'scripts/get_features.pl',
			'scripts/get_feature_info.pl',
			'scripts/get_gene_regions.pl',
			'scripts/get_intersecting_features.pl',
			'scripts/get_relative_data.pl',
			'scripts/gff3_to_ucsc_table.pl',
			'scripts/graph_data.pl',
			'scripts/graph_histogram.pl',
			'scripts/graph_profile.pl',
			'scripts/intersect_nucs.pl',
			'scripts/intersect_SNPs.pl',
			'scripts/join_data_file.pl',
			'scripts/locate_SNPs.pl',
			'scripts/manipulate_datasets.pl',
			'scripts/map_nucleosomes.pl',
			'scripts/map_oligo_data2gff.pl',
			'scripts/map_transcripts.pl',
			'scripts/merge_datasets.pl',
			'scripts/novo_wrapper.pl',
			'scripts/process_microarray.pl',
			'scripts/pull_features.pl',
			'scripts/run_cluster.pl',
			'scripts/split_bam_by_isize.pl',
			'scripts/split_data_file.pl',
			'scripts/ucsc_cytoband2gff3.pl',
			'scripts/ucsc_table2gff3.pl',
			'scripts/verify_nucleosome_mapping.pl',
			'scripts/wig2data.pl',
		],
	);
	return \%options;
}


sub check_useq {
	if ($] >= 5.010000) {
		# Bio::DB::USeq requires perl 5.10
		$options->{'requires'}{'Bio::DB::USeq'} = 0.20;
	}
}


sub check_sam {
	
	# check to see if it is installed
	my $sam_ok;
	eval {require Bio::DB::Sam; $sam_ok = 1;};
	if ($sam_ok) {
		# Bam support is currently installed
		# request a minimum version
		# if they don't meet this minimum, let's hope the user
		# knows how to rectify it.....
		$options->{'requires'}{'Bio::DB::Sam'} = 1.36;
	}
	else {
		$options->{'recommends'}{'Bio::DB::Sam'} = 1.36;
	}
}


sub check_big {
	
	# check to see if it is installed
	my $big_ok;
	eval {require Bio::DB::BigFile; $big_ok = 1;};
	if ($big_ok) {
		# BigFile support is currently installed
		# request a minimum version
		# if they don't meet this minimim, let's hope the user
		# knows how to rectify it.....
		$options->{'requires'}{'Bio::DB::BigFile'} = 1.07;
	}
	else {
		$options->{'recommends'}{'Bio::DB::BigFile'} = 1.07;
	}
}



