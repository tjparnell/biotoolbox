#!/usr/bin/env perl

use strict;
use warnings;
use Config;
use File::Spec;
use Module::Build;



### Prepare global variables for installation

# initial build options
my $options = get_base_options();

# check environment 
my $make_present = check_make();
my @paths = check_paths();



### Check for additional Modules to install
# I could just try to make all of these recommended modules
# but I am trying to make things easier on the user, so 
check_cluster();
check_db();
check_gd();
check_sam();
check_big();



### Build the script
my $build = Module::Build->new(%$options);
$build->create_build_script;
exit 0;


###### subroutines to assist in building ########


sub get_base_options {
	my %options = (
		# build_class       => 'Module::Build',
		module_name       => 'Bio::ToolBox',
		license           => 'perl',
		dist_version_from => 'lib/Bio/ToolBox.pm',
		dist_author       => 'Timothy Parnell <parnell.tj@gmail.com>',
		dist_abstract     => 'Tools for querying and analysis of genomic data',
		requires          => {
			'Config::Simple'          => 0,
			'Statistics::Lite'        => 3.2,
			'Statistics::Descriptive' => 3.0,
			'Statistics::LineFit'     => 0,
			'Parallel::ForkManager'   => 0,
			'Bio::Root::Version'      => '1.0069001',
			'Bio::DB::USeq'           => 0.20, # only because this is pure perl
		},
		recommends        => {
			'Archive::Zip'            => 0,
			'Bio::Graphics'           => 0,
			'File::Which'             => 0,
			'Net::FTP'                => 0,
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
			'scripts/print_feature_types.pl',
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


sub check_make {
	# I am making the assumption this is going to be in the 
	# typical system bin directory
	return (-e '/usr/bin/make' && -x _);
}


sub check_paths {
	# I am going to hard encode every conceivable path that seems likely
	# these are assuming standard unix paths
	# we are assuming that standard lib and inc paths are below these
	# as well as other package directories, such as samtools and kent
	my @paths;
	foreach (
		File::Spec->curdir,          # current directory, just in case
		$ENV{HOME},                  # home
		'/usr',                      # root unix library
		'/usr/local',                # local unix library
		'/opt',                      # optional tree in linux, MacPort on OS X
		'/opt/local',
		'/sw',                       # Fink on Mac OS X
		'/home/BioApps/SamTools',    # HCI Bio servers
	) {
		push @paths, $_ if -e $_;
	}
	return @paths;
}

sub check_cluster {
	return unless $make_present;
	$options->{'requires'}{'Algorithm::Cluster'} = 0;
}


sub check_db {
	return unless $make_present;
	
	# add SQLite as a fully featured and simple database
	$options->{'requires'}{'DBI'} = 0;
	$options->{'requires'}{'DBD::SQLite'} = 0;
	
	# add mysql as an advanced option
	$options->{'recommends'}{'DBD::mysql'} = 0;
}


sub check_gd {
	return unless $make_present;
	
	my $gd_config = `gdlib-config --all`;
	if ($gd_config) {
		# pretty good evidence the libgd libraries are installed
		$options->{'requires'}{'GD'} = 0;
			# this is sometimes problematic on Mac OS X without a force????
			# due to faulty tests??? not actual problems
		$options->{'requires'}{'GD::Graph'} = 0;
		$options->{'requires'}{'GD::Graph::smoothlines'} = 0;
	}
}


sub check_sam {
	return unless $make_present;
	
	# check to see if it is installed
	eval {require Bio::DB::Sam; 1;};
	if ($@) {
		# Bam support is not installed
		# first check if we can find a current samtools directory
		my $samtools;
		PATH: foreach my $path (@paths) {
			my $p = File::Spec->catdir($path, 'samtools*');
			foreach my $candidate (reverse glob($p)) {
				if (-e File::Spec->catdir($candidate, 'libbam.a')) {
					$samtools = $candidate;
					last PATH;
				}
			}
		}
		if ($samtools) {
			# we have libbam.a, I'm going to assume the include files 
			# are present too
			# go ahead and require it, I think we are good
			print "Found samtools library at $samtools\n";
			system('export', "SAMTOOLS=$samtools");
			$options->{'requires'}{'Bio::DB::Sam'} = 1.36;
		}
		else {
# 			warn <<SAM
# 
# A samtools library could not be found for building the Bio::DB::Sam 
# module for working with Bam files. You can install it manually or 
# run the install/install_samtools.sh to install it.
# SAM
# ;
			$options->{'recommends'}{'Bio::DB::Sam'} = 1.36;
		}	
	}
	else {
		# Bam support is currently installed
		# request a minimum version
		# if they don't meet this minimum, let's hope the user
		# knows how to rectify it.....
		$options->{'requires'}{'Bio::DB::Sam'} = 1.36;
	}
}


sub check_big {
	return unless $make_present;
	
	# check to see if it is installed
	eval {require Bio::DB::BigFile; 1;};
	if ($@) {
		# Big file support is not installed
		# first check if we can find a UCSC kent source directory
		my $kent;
		PATH: foreach my $path (@paths) {
			foreach my $mach (qw(x86_64 i686 i386 ppc alpha sparc)) {
				my $p = File::Spec->catdir($path, 'kent', 'src', 'lib', $mach);
				if (-e File::Spec->catdir($p, 'jkweb.a')) {
					$kent = $p;
					last PATH;
				}
			}
		}
		if ($kent) {
			print <<KENT

A UCSC kent source directory was found at $kent
so I will add Bio::DB::BigFile as a requirement. If the build 
fails, then you will need to build Bio::DB::BigFile manually.
KENT
;
			$options->{'requires'}{'Bio::DB::BigFile'} = 1.07;
		}
		else {
# 			print <<KENT
# 
# A UCSC kent source library could not be found for building the 
# Bio::DB::BigFile module for working with BigWig and BigBed files. 
# You can install it manually or run the included script 
# to install it.
# KENT
		}	
	}
	else {
		# BigFile support is currently installed
		# request a minimum version
		# if they don't meet this minimim, let's hope the user
		# knows how to rectify it.....
		$options->{'requires'}{'Bio::DB::BigFile'} = 1.07;
	}
}



