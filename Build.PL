#!/usr/bin/perl

use strict;
use warnings;
use File::Spec;
use Module::Build;



### Prepare global variables for installation

# initial build options
my $options = get_base_options();

# check environment 
my @paths = check_paths();
my $make_present = check_make();



### Check for additional Modules to install
# I am trying to make things easier on the user, so checking for 
# various external libraries and build environment before deciding 
# to require or just recommend
check_useq();
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
			'Statistics::LineFit'     => 0.06,
			'Parallel::ForkManager'   => 1.02,
			'Bio::Root::Version'      => '1.0069001',
			'Net::FTP'                => 0,
			'File::Which'             => 0,
		},
		recommends        => {
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


sub check_paths {
	# I am going to hard encode every conceivable path that seems likely
	# these are assuming standard unix paths
	# we are assuming that standard lib and inc paths are below these
	# as well as other package directories, such as samtools and kent
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
	push @paths, File::Spec->path; # environment paths
	return @paths;
}


sub check_make {
	my $found = 0;
	foreach my $path (@paths) {
		$found++ if (-e File::Spec->catdir($path, 'make') and -x _);
		$found++ if (-e File::Spec->catdir($path, 'bin', 'make') and -x _);
		last if $found;
	}
	unless ($found) {
		warn <<NOMAKE

############################################################################
A compatible C compiler was not found. This is required to install many of 
the required Perl modules, as well many optional modules that require 
external libraries and compiling C extensions.

This installation will likely fail.

On Linux, you can install GCC using your package manager.
On Mac OS X, you can install XCode, plus the optional Command Line Tools.
############################################################################

NOMAKE
;
	}
	return $found;
}


sub check_useq {
	if ($] >= 5.010000) {
		# Bio::DB::USeq requires perl 5.10
		$options->{'requires'}{'Archive::Zip'}  = 0;
		$options->{'requires'}{'Bio::DB::USeq'} = 0.20;
	}
}


sub check_cluster {
	if ($make_present) {
		$options->{'requires'}{'Algorithm::Cluster'} = 0;
	}
	else {
		$options->{'recommends'}{'Algorithm::Cluster'} = 0;
	}
}


sub check_db {
	if ($make_present) {
		# add SQLite as a fully featured and simple database
		$options->{'requires'}{'DBI'} = 0;
		$options->{'requires'}{'DBD::SQLite'} = 0;
	
		# add mysql as an advanced option
		foreach my $path (@paths) {
			if (
				(-e File::Spec->catdir($path, 'mysql_config') and -x _) or
				(-e File::Spec->catdir($path, 'bin', 'mysql_config') and -x _)
			) {
				# pretty good evidence that MySQL is installed
				$options->{'recommends'}{'DBD::mysql'} = 0;
				last;
			}
		}
	}
	else {
		$options->{'recommends'}{'DBI'} = 0;
		$options->{'recommends'}{'DBD::SQLite'} = 0;
	}
}


sub check_gd {
	
	# check whether GD is installed
	my $gd_ok;
	eval {require 'GD'; $gd_ok = 1};
	if ($gd_ok) {
		# GD support is installed
		$options->{'requires'}{'GD'} = 2.44;
		$options->{'requires'}{'GD::Graph'} = 1.44;
		$options->{'requires'}{'GD::Graph::smoothlines'} = 1.6;
		$options->{'recommends'}{'Bio::Graphics'} = 2.29;
	}
	else {
		if ($make_present) {
			# look for whether external libgd libraries are installed
			foreach my $path (@paths) {
				if (
					(-e File::Spec->catdir($path, 'gdlib-config') and -x _) or
					(-e File::Spec->catdir($path, 'bin', 'gdlib-config') and -x _)
				) {
					# pretty good evidence the libgd libraries are installed
					$options->{'requires'}{'GD'} = 2.44;
					$options->{'requires'}{'GD::Graph'} = 1.44;
					$options->{'requires'}{'GD::Graph::smoothlines'} = 1.6;
					$options->{'recommends'}{'Bio::Graphics'} = 2.29;
					last;
				}
			}
		}
		else {
			$options->{'recommends'}{'GD'} = 2.44;
			$options->{'recommends'}{'GD::Graph'} = 1.44;
			$options->{'recommends'}{'GD::Graph::smoothlines'} = 1.6;
			$options->{'recommends'}{'Bio::Graphics'} = 2.29;
		}
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
		if ($samtools and $make_present) {
			# we have libbam.a, I'm going to assume the include files 
			# are present too
			# go ahead and require it, I think we are good
			print "Found samtools library at $samtools\n";
			eval {system('export', "SAMTOOLS=$samtools")};
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
		# Big file support is not installed
		# first check if we can find a UCSC kent source directory
		my $kent;
		my $machtype;
		PATH: foreach my $path (@paths) {
			my $p = File::Spec->catdir($path, 'kent', 'src');
			foreach my $mach (qw(x86_64 i686 i386 ppc alpha sparc)) {
				if (-e File::Spec->catdir($p, 'lib', $mach, 'jkweb.a')) {
					$kent = $p;
					$machtype = $mach;
					last PATH;
				}
			}
		}
		if ($kent and $make_present) {
			# we found the kent source directory
			# so it could be possible to build BigFile
			print <<KENT

A UCSC kent source directory was found at $kent
so I will add Bio::DB::BigFile as a requirement. If the build 
fails, then you will need to build Bio::DB::BigFile manually.
KENT
;
			$options->{'requires'}{'Bio::DB::BigFile'} = 1.07;
			eval {
				# try exporting the paths for the module to find later
				# this assumes a typical bash shell in a unix environment
				# which is a reasonable 
				system('export', "KENT_SRC=$kent");
				system('export', "MACHTYPE=$machtype");
			};
		}
		else {
			$options->{'recommends'}{'Bio::DB::BigFile'} = 1.07;
# 			print <<KENT
# 
# A UCSC kent source library could not be found for building the 
# Bio::DB::BigFile module for working with BigWig and BigBed files. 
# You can install it manually or run the included script 
# to install it.
# KENT
		}	
	}
}



