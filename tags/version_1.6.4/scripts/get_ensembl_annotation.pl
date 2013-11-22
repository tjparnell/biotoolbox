#!/usr/bin/perl

# a script to get ensEMBL annotation and write it out as GFF3


use strict;
use IO::File;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqFeature::Lite;
# use Data::Dumper;

# Check for Bio::EnsEMBL
my $bio_ensembl = 0;
eval {
	require Bio::EnsEMBL::Registry;
	Bio::EnsEMBL::Registry->import;
	$bio_ensembl = 1;
};
my $VERSION = '1.4.1';
	
print "\n A script to fetch genomic annotation from public Ensembl databases\n\n";


### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Get command line options and initialize values
my (
	$species,
	$get_chromo,
	$get_protein_genes,
	$get_rna_genes,
	$get_misc_rna_genes,
	$get_snrna_genes,
	$get_snorna_genes,
	$get_mirna_genes,
	$get_rrna_genes,
	$get_trna_genes,
	$group,
	$host,
	$user,
	$pass,
	$prefix,
	$outfile,
	$printdb,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'species=s' => \$species, # the species to look up
	'chromo!'   => \$get_chromo, # collect chromosome info
	'protein!'  => \$get_protein_genes, 
	'rna!'      => \$get_rna_genes,
	'miscrna!'  => \$get_misc_rna_genes,
	'snrna!'    => \$get_snrna_genes,
	'snorna!'   => \$get_snorna_genes,
	'mirna!'    => \$get_mirna_genes,
	'rrna!'     => \$get_rrna_genes,
	'trna!'     => \$get_trna_genes,
	'group=s'   => \$group, # the database group
	'host=s'    => \$host, # host address
	'user=s'    => \$user, # user name to log in
	'pass=s'    => \$pass, # password to log in with
	'prefix!'   => \$prefix, # prefix chromosome names with chr
	'out=s'     => \$outfile, # name of output file 
	'printdb'   => \$printdb, # print all of the available databases
	'help'      => \$help, # request help
	'version'   => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script get_ensembl_annotation.pl, version $VERSION\n\n";
	exit;
}



# Check EnsEMBL
unless ($bio_ensembl) {
	die "\n Bio::EnsEMBL modules are not installed. Please see help\n";
}


### Check for requirements and set defaults
# required
unless ($species or $printdb) {
	die " Must define a species! See documentation\n";
}
unless ($outfile) {
	$outfile = $species . '_annotation';
}

# connectivity
unless ($host) {
	$host = 'ensembldb.ensembl.org';
}
unless ($user) {
	$user = 'anonymous';
}
unless ($group) {
	$group = 'core';
}

# set default collection types
unless (defined $get_chromo) {
	$get_chromo = 1;
}
unless (defined $get_protein_genes) {
	$get_protein_genes = 1;
}
if ($get_rna_genes) {
	# we will set all the individual rna genes to true unless already defined
	$get_misc_rna_genes = 1 unless defined $get_misc_rna_genes;
	$get_snrna_genes    = 1 unless defined $get_snrna_genes;
	$get_snorna_genes   = 1 unless defined $get_snorna_genes;
	$get_mirna_genes    = 1 unless defined $get_mirna_genes;
	$get_rrna_genes     = 1 unless defined $get_rrna_genes;
	$get_trna_genes     = 1 unless defined $get_trna_genes;
}




### Establish ensembl connection
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
		-host    => $host,
		-user    => $user,
		-pass    => $pass,
) or die " Can't connect to registry!\n";




### Print available database names
if ($printdb) {
	# print all the slice adaptors
	my @db_adaptors = @{ $registry->get_all_DBAdaptors() };
	
	print " These are the available public databases:\n\n";
	foreach my $db_adaptor (@db_adaptors) {
		my $db_connection = $db_adaptor->dbc();
	
		printf(
			"species/group\t%s/%s\n database\t%s\n\n",
			$db_adaptor->species(), 
			$db_adaptor->group(),
			$db_connection->dbname()
		);
	}
	exit;
}



### Connect to database
my $slice_adaptor = $registry->get_adaptor($species, $group, 'slice') or
	die " Can't get slice adaptor!\n";
my $db_connection = $slice_adaptor->dbc;
my $db_name = $db_connection->dbname;
print " Success in connecting to database $db_name\n";





### Open output
# open filehandle
unless ($outfile =~ /\.gff3?$/) {
	$outfile .= '.gff3';
}
my $gff_fh = new IO::File $outfile, 'w';
unless ($gff_fh) {
	die " unable to open file handle for '$outfile'!\n";
}

# write headers
$gff_fh->print("##gff-version 3\n");
$gff_fh->print("# EnsEMBL data for species $species\n");
$gff_fh->print("# Collected from database $db_name\n");








### Begin siphoning data

# Get chromosomes
my $slices_ref = $slice_adaptor->fetch_all('toplevel');
print " Identified " . scalar(@{ $slices_ref }) . " toplevel chromosomes/contigs\n";

# Process chromosomes
my $phase; # global variable for keeping track of CDS phases
foreach my $slice (@{ $slices_ref }) {
	
	print "  Collecting features for chromosome " . 
		$slice->seq_region_name . "\n";
	
	# record how many features we write
	my $chr_feature_count = 0;
	
	# record chromsomes if requested
	if ($get_chromo) {
		
		# set chromosome name
		my $name = $slice->seq_region_name;
		if ($slice->coord_system_name eq 'chromosome' and $prefix) {
			$name = 'chr' . $name;
		}
		
		# generate the chromosome SeqFeature object
		my $chromo_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $name,
				-source        => $db_name,
				-primary_tag   => $slice->coord_system_name, 
				-start         => 1,
				-end           => $slice->seq_region_length,
				-strand        => 0,
				-phase         => '.',
				-display_name  => $name,
				-primary_id    => $name,
		);
		
		# write the object
		$chromo_sf->version(3);
		$gff_fh->print( $chromo_sf->gff_string . "\n");
		$chr_feature_count++;
	}
	
# 	#### DEBUGGING LIMITER ####
# 	unless ($slice->seq_region_name eq '11') {
# 		next;
# 	}
# 	###########################
	
	# collect the protein_coding genes
	if ($get_protein_genes) {
		$chr_feature_count += process_genes($slice, 'protein_coding');
	}
	
	# collect misc RNA genes
	if ($get_misc_rna_genes) {
		$chr_feature_count += process_genes($slice, 'misc_RNA');
	}
	
	# collect snRNA genes
	if ($get_snrna_genes) {
		$chr_feature_count += process_genes($slice, 'snRNA');
	}
	
	# collect snoRNA genes
	if ($get_snorna_genes) {
		$chr_feature_count += process_genes($slice, 'snoRNA');
	}
	
	# collect miRNA genes
	if ($get_mirna_genes) {
		$chr_feature_count += process_genes($slice, 'miRNA');
	}
	
	# collect rRNA genes
	if ($get_rrna_genes) {
		$chr_feature_count += process_genes($slice, 'rRNA');
	}
	
	# collect tRNA genes
	if ($get_trna_genes) {
		$chr_feature_count += process_genes($slice, 'tRNA');
	}
	
	# finished with this chromosome
	if ($chr_feature_count) {
		# print directive to close out all previous genes
		$gff_fh->print("###\n"); 
		print " Wrote $chr_feature_count features for chromosome " . 
			$slice->seq_region_name . "\n";
	}
}





### Finish
$gff_fh->close;
print " Finished! Wrote output file '$outfile'.\n";








########################   Subroutines   ###################################


sub process_genes {
	my ($slice, $biotype) = @_;
	
	# get chromosome name
	my $chr = $slice->seq_region_name;
	if ($slice->coord_system_name eq 'chromosome' and $prefix) {
		$chr = 'chr' . $chr;
	}
	
	# retrieve a list of genes
	my $genes_ref = $slice->get_all_Genes_by_type($biotype, undef, 1);
		# all genes of $biotype, no restriction on source logic name, 
		# retrieve all information
	print "    Pulled " . scalar(@{ $genes_ref }) . " $biotype genes....\n";
	
	# process the genes
	my $count = 0;
	foreach my $gene ( @{ $genes_ref } ) {
		if ($biotype eq 'protein_coding') {
			process_coding_gene($gene, $chr);
		}
		else {
			process_rna_gene($gene, $chr);
		}
		$count++;
	}
	
	return $count;
}




sub process_coding_gene {
	my ($gene, $chr) = @_;
	
	# create the SeqFeature gene object
	my $gene_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $chr,
				-source        => $gene->source || 'ensembl',
				-primary_tag   => 'gene',
				-start         => $gene->start,
				-end           => $gene->end,
				-strand        => $gene->strand,
				-phase         => '.',
				-display_name  => $gene->stable_id,
				-primary_id    => $gene->stable_id,
	);
	$gene_sf->add_tag_value('status', $gene->status);
	
	# get additional information
	my $alias = $gene->external_name;
	if (defined $alias) {
		$alias =~ s/ {3,}//g; # strip excessive spaces
		$gene_sf->add_tag_value('Alias', $alias);
	}
	my $note_text = $gene->description;
	if (defined $note_text) {
		$note_text =~ s/ {3,}//g; # strip excessive spaces
		$gene_sf->add_tag_value('Note', $note_text);
	}
	
	# work through the transcripts
	foreach my $transcript (@{ $gene->get_all_Transcripts }) {
		
		# generate SeqFeature Transcript object
		my $trnscpt_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $chr,
				-source        => $gene->source || 'ensembl',
				-primary_tag   => 'mRNA',
				-start         => $transcript->start,
				-end           => $transcript->end,
				-strand        => $transcript->strand,
				-phase         => '.',
				-display_name  => $transcript->stable_id,
				-primary_id    => $transcript->stable_id,
		);
		
		# add alias
		if (defined $alias) {
			$trnscpt_sf->add_tag_value('Alias', $alias);
		}
		
		# get transcription start/stop
		my $coding_start = $transcript->coding_region_start || undef;
		my $coding_stop  = $transcript->coding_region_end || undef;
		$phase = 0; # reset phase to zero for all new transcripts
		
		# add the exons
		my $exons_ref = $transcript->get_all_Exons;
		add_exon_codon_subfeatures(
			$trnscpt_sf, $exons_ref, $coding_start, $coding_stop);
		
		# process CDS and UTRs
		if ($transcript->strand > 0) {
			# forward strand
			process_forward_exons(
				$trnscpt_sf, $exons_ref, $coding_start, $coding_stop);
		}
		else {
			# reverse strand
			process_reverse_exons(
				$trnscpt_sf, $exons_ref, $coding_start, $coding_stop);
		
		}
		
		# add the transcript
		$gene_sf->add_SeqFeature($trnscpt_sf);
	}
	
	
	# print the gene feature and its subfeatures
	$gene_sf->version(3);
	$gff_fh->print( $gene_sf->gff_string(1) . "\n");
		# the gff_string method is undocumented in the POD, but is a 
		# valid method. Passing 1 should force a recursive action to 
		# print parent and children.
}


sub add_exon_codon_subfeatures {
	my ($trnscpt_sf, $exons_ref, $coding_start, $coding_stop) = @_;
	
	# generate the start and stop codons
	if ($trnscpt_sf->strand == 1) {
		# forward strand
		
		# in some situations, there may be protein_coding gene transcript
		# without a coding start or stop !?
		if (defined $coding_start and defined $coding_stop) {
			# start codon
			$trnscpt_sf->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
					-seq_id        => $trnscpt_sf->seq_id,
					-source        => $trnscpt_sf->source,
					-primary_tag   => 'start_codon',
					-start         => $coding_start,
					-end           => $coding_start + 2,
					-strand        => 1,
					-phase         => 0,
					-display_name  => $trnscpt_sf->display_name . '.start_codon',
					-primary_id    => $trnscpt_sf->display_name . '.start_codon',
				)
			);
			
			# stop codon
			$trnscpt_sf->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
					-seq_id        => $trnscpt_sf->seq_id,
					-source        => $trnscpt_sf->source,
					-primary_tag   => 'stop_codon',
					-start         => $coding_stop - 2,
					-end           => $coding_stop,
					-strand        => 1,
					-phase         => 0,
					-display_name  => $trnscpt_sf->display_name . '.stop_codon',
					-primary_id    => $trnscpt_sf->display_name . '.stop_codon',
				)
			);
		}
	}
	
	else {
		# reverse strand
		
		# in some situations, there may be protein_coding gene transcript
		# without a coding start or stop !?
		if (defined $coding_start and defined $coding_stop) {
		# stop codon
			$trnscpt_sf->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
					-seq_id        => $trnscpt_sf->seq_id,
					-source        => $trnscpt_sf->source,
					-primary_tag   => 'stop_codon',
					-start         => $coding_start,
					-end           => $coding_start + 2,
					-strand        => -1,
					-phase         => 0,
					-display_name  => $trnscpt_sf->display_name . '.stop_codon',
					-primary_id    => $trnscpt_sf->display_name . '.stop_codon',
				)
			);
			
			# start codon
			$trnscpt_sf->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
					-seq_id        => $trnscpt_sf->seq_id,
					-source        => $trnscpt_sf->source,
					-primary_tag   => 'start_codon',
					-start         => $coding_stop - 2,
					-end           => $coding_stop,
					-strand        => -1,
					-phase         => 0,
					-display_name  => $trnscpt_sf->display_name . '.start_codon',
					-primary_id    => $trnscpt_sf->display_name . '.start_codon',
				)
			);
		}
	}
	
	# add exons
	my $exon_count = 1;
	foreach my $exon ( @{ $exons_ref } ) {
		
		# add the exon as a subfeature
		$trnscpt_sf->add_SeqFeature( 
			Bio::SeqFeature::Lite->new(
				-seq_id        => $trnscpt_sf->seq_id,
				-source        => $trnscpt_sf->source,
				-primary_tag   => 'exon',
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => $exon->strand,
				-display_name  => $exon->stable_id,
				-primary_id    => $trnscpt_sf->display_name . ".exon.$exon_count",
			)
		);
		
		# we're giving a unique id based on transcript name appended with 
		# exon and incrementing number
		$exon_count++;
	}
}



sub process_forward_exons {
	my ($transcript, $exons_ref, $cdsStart, $cdsStop) = @_;
	
	# exon counter
	# we need to make the primary ids unique, and since exons may be 
	# shared between multiple transcripts, we can't use the exon's id
	# we'll use the transcript id plus an incrementing number to make unique
	my $ex_count = 1; 
	
	
	foreach my $exon ( @{ $exons_ref } ) {
		
		# we need to determine whether the exon is UTR, CDS, or split both
		
		#### No CDS whatsoever ####
		if (!defined $cdsStart and !defined $cdsStop) {
			# nothing to write
			# there should already be an exon written
		}
		
		#### 5'UTR only ####
		elsif (
			$exon->start < $cdsStart # cdsStart
			and
			$exon->end < $cdsStart
		) {
			# the exon start/end is entirely before the cdsStart
			# we have a 5'UTR
			
			# build the utr SeqFeature object
			my $exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => 1, # we know this already
				-phase         => '.',
				-primary_tag   => 'five_prime_UTR',
				-primary_id    => $transcript->primary_id . ".utr.$ex_count",
				-display_name  => $exon->stable_id,
			);
			
			# associate exon with the parent transcript
			$transcript->add_SeqFeature($exon_sf);
		}
		
		
		#### Split 5'UTR and CDS ####
		elsif (
			$exon->start < $cdsStart 
			and
			$exon->end >= $cdsStart
		) {
			# the start codon is in this exon
			# we need to make two features, 
			# one for the utr half, other for the cds
			
			# build the utr half of the object
			my $utr_exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $cdsStart - 1,
				-strand        => 1,
				-phase         => '.',
				-primary_tag   => 'five_prime_UTR',
				-primary_id    => $transcript->primary_id . ".utr.$ex_count",
				-display_name  => $exon->stable_id . '.utr',
			);
			# since we're actually building two objects here, we have to 
			# complete the process of building the utr object and associate 
			# it with the transcript object
			# the cds half of the exon will be finished below
			
			# associate add the utr half to the parent transcript
			$transcript->add_SeqFeature($utr_exon_sf);
			
			# now build the cds half of the object
			my $exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $cdsStart,
				-end           => $exon->end,
				-strand        => 1,
				-phase         => $phase,
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds.$ex_count",
				-display_name  => $exon->stable_id . '.cds',
			);
			
			# associate exon with the parent transcript
			$transcript->add_SeqFeature($exon_sf);
			
			# reset phase for next CDS
			# phase + (3 - (length % 3)), readjust to 0..2 if necessary
			# adapted from Barry Moore's gtf2gff3.pl script
			$phase = $phase + (3 - ( ($exon->end - $cdsStart + 1) % 3) );
			$phase -=3 if $phase > 2;
		}
		
		#### CDS only ####
		elsif (
			$exon->start >= $cdsStart 
			and
			$exon->end <= $cdsStop
		) {
			# we are in the CDS
			# this will also work with genes that have no UTRs, where the 
			# the cdsStart == txStart, and same with End
			
			my $exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => 1,
				-phase         => $phase,
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds.$ex_count",
				-display_name  => $exon->stable_id,
			);
			
			# associate exon with the parent transcript
			$transcript->add_SeqFeature($exon_sf);
			
			# reset phase for next CDS
			# phase + (3 - (length % 3)), readjust to 0..2 if necessary
			# adapted from Barry Moore's gtf2gff3.pl script
			$phase = $phase + (3 - ( ($exon->end - $exon->start + 1) % 3) );
			$phase -=3 if $phase > 2;
		}
		
		
		#### Split CDS and 3'UTR ####
		elsif (
			$exon->start <= $cdsStop 
			and
			$exon->end > $cdsStop
		) {
			# the stop codon is in this exon
			# we need to make two features, 
			# one for the cds half, other for the utr
			
			# build the cds half of the object
			my $exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $cdsStop,
				-strand        => 1,
				-phase         => $phase,
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds.$ex_count",
				-display_name  => $exon->stable_id . '.cds',
			);
			
			# now build the utr half of the object
			my $utr_exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $cdsStop + 1,
				-end           => $exon->end,
				-strand        => 1,
				-phase         => '.',
				-primary_tag   => 'three_prime_UTR',
				-primary_id    => $transcript->primary_id . ".utr.$ex_count",
				-display_name  => $exon->stable_id . '.utr',
			);
			
			# associate exons with the parent transcript
			$transcript->add_SeqFeature($exon_sf);
			$transcript->add_SeqFeature($utr_exon_sf);
			
			# reset phase for next CDS - huh?
			# phase + (3 - (length % 3)), readjust to 0..2 if necessary
			# adapted from Barry Moore's gtf2gff3.pl script
			$phase = $phase + (3 - ( ($cdsStop - $exon->start + 1) % 3) );
			$phase -=3 if $phase > 2;
		}
		
		#### 3'UTR only ####
		elsif (
			$exon->start > $cdsStop 
			and
			$exon->end > $cdsStop
		) {
			# the exon start/end is entirely after the cdsStop
			# we have a 3'UTR
			
			# build the utr object
			my $exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => 1,
				-phase         => '.',
				-primary_tag   => 'three_prime_UTR',
				-primary_id    => $transcript->primary_id . ".utr.$ex_count",
				-display_name  => $exon->stable_id,
			);
			
			# associate exon with the parent transcript
			$transcript->add_SeqFeature($exon_sf);
		}
		
		#### Something's wrong ####
		else {
			# just in case I goofed something up
			warn " programming error! the exon coordinates don't match up " .
				"with CDS coordinates for transcript '" . 
				$transcript->display_name . "'\n " . 
				" Forward strand exon coordinates " . $exon->start . ".." .
				$exon->end . "\n" .
				" CDS coordinates $cdsStart .. $cdsStop\n";
		}
		
		# increment exon counter
		$ex_count++;
	}
}





sub process_reverse_exons {
	my ($transcript, $exons_ref, $cdsStart, $cdsStop) = @_;
	
	# exon counter
	# we need to make the primary ids unique, and since exons may be 
	# shared between multiple transcripts, we can't use the exon's id
	# we'll use the transcript id plus an incrementing number to make unique
	my $ex_count = 1; 
	
	
	foreach my $exon ( @{ $exons_ref } ) {
		
		# we need to determine whether the exon is UTR, CDS, or split both
		
		#### No CDS whatsoever ####
		if (!defined $cdsStart and !defined $cdsStop) {
			# nothing to write
			# there should already be an exon written
		}
		
		#### 3'UTR only ####
		elsif (
			$exon->start < $cdsStart # cdsStart
			and
			$exon->end < $cdsStart
		) {
			# the exon start/end is entirely before the cdsStart
			# we have a 3'UTR
			
			# build the utr object
			my $exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => -1,
				-phase         => '.',
				-primary_tag   => 'three_prime_UTR',
				-primary_id    => $transcript->primary_id . ".utr.$ex_count",
				-display_name  => $exon->stable_id,
			);
			
			# associate exon with the parent transcript
			$transcript->add_SeqFeature($exon_sf);
			
		}
		
		
		#### Split 3'UTR and CDS ####
		elsif (
			$exon->start < $cdsStart 
			and
			$exon->end >= $cdsStart
		) {
			# the (stop) codon is in this exon
			# we need to make two features, 
			# one for the utr half, other for the cds
			
			# build the cds half of the object
			my $exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $cdsStart,
				-end           => $exon->end,
				-strand        => -1,
				-phase         => $phase,
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds.$ex_count",
				-display_name  => $exon->stable_id . '.cds',
			);
			
			# now build the utr half of the object
			my $utr_exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $cdsStart - 1,
				-strand        => -1,
				-phase         => '.',
				-primary_tag   => 'three_prime_UTR',
				-primary_id    => $transcript->primary_id . ".utr.$ex_count",
				-display_name  => $exon->stable_id . '.utr',
			);
			# since we're actually building two objects here, we have to 
			# complete the process of building the utr object and associate 
			# it with the transcript object
			# the cds half of the exon will be finished below
			
			# associate exons with the parent transcript
			$transcript->add_SeqFeature($exon_sf);
			$transcript->add_SeqFeature($utr_exon_sf);
			
			# reset phase for next CDS - huh?
			# phase + (3 - (length % 3)), readjust to 0..2 if necessary
			# adapted from Barry Moore's gtf2gff3.pl script
			$phase = $phase + (3 - ( ($exon->end - $cdsStart + 1) % 3) );
			$phase -=3 if $phase > 2;
		}
		
		#### CDS only ####
		elsif (
			$exon->start >= $cdsStart 
			and
			$exon->end <= $cdsStop
		) {
			# we are in the CDS
			# this will also work with genes that have no UTRs, where the 
			# the cdsStart == txStart, and same with End
			
			my $exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => -1,
				-phase         => $phase,
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds.$ex_count",
				-display_name  => $exon->stable_id,
			);
			
			# associate exon with the parent transcript
			$transcript->add_SeqFeature($exon_sf);
			
			# reset phase for next CDS
			# phase + (3 - (length % 3)), readjust to 0..2 if necessary
			# adapted from Barry Moore's gtf2gff3.pl script
			$phase = $phase + (3 - ( ($exon->end - $exon->start + 1) % 3) );
			$phase -=3 if $phase > 2;
		}
		
		
		#### Split CDS and 5'UTR ####
		elsif (
			$exon->start <= $cdsStop 
			and
			$exon->end > $cdsStop
		) {
			# the (start) codon is in this exon
			# we need to make two features, 
			# one for the cds half, other for the utr
			
			# build the utr half of the object
			my $utr_exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $cdsStop + 1,
				-end           => $exon->end,
				-strand        => -1,
				-phase         => '.',
				-primary_tag   => 'five_prime_UTR',
				-primary_id    => $transcript->primary_id . ".utr.$ex_count",
				-display_name  => $exon->stable_id . '.utr',
			);
						
			# now build the cds half of the object
			my $exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $cdsStop,
				-strand        => -1,
				-phase         => $phase,
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds.$ex_count",
				-display_name  => $exon->stable_id . '.cds',
			);
			
			# associate exons with the parent transcript
			$transcript->add_SeqFeature($utr_exon_sf);
			$transcript->add_SeqFeature($exon_sf);
			
			# reset phase for next CDS
			# phase + (3 - (length % 3)), readjust to 0..2 if necessary
			# adapted from Barry Moore's gtf2gff3.pl script
			$phase = $phase + (3 - ( ($cdsStop - $exon->start + 1) % 3) );
			$phase -=3 if $phase > 2;
		}
		
		#### 5'UTR only ####
		elsif (
			$exon->start > $cdsStop 
			and
			$exon->end > $cdsStop
		) {
			# the exon start/end is entirely after the cdsStop
			# we have a 5'UTR
			
			# build the utr object
			my $exon_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => -1,
				-phase         => '.',
				-primary_tag   => 'five_prime_UTR',
				-primary_id    => $transcript->primary_id . ".utr.$ex_count",
				-display_name  => $exon->stable_id,
			);
			
			# associate exon with the parent transcript
			$transcript->add_SeqFeature($exon_sf);
		}
			
		#### Something's wrong ####
		else {
			# just in case I goofed something up
			warn " programming error! the exon coordinates don't match up " .
				"with CDS coordinates for transcript '" . 
				$transcript->display_name . "'\n " . 
				" Reverse strand exon coordinates " . $exon->start . ".." .
				$exon->end . "\n" .
				" CDS coordinates $cdsStart .. $cdsStop\n";
		}
		
		# increment exon counter
		$ex_count++;
	}
}




sub process_rna_gene {
	my ($gene, $chr) = @_;
	
	# create the SeqFeature gene object
	my $gene_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $chr,
				-source        => $gene->source || 'ensembl',
				-primary_tag   => 'gene',
				-start         => $gene->start,
				-end           => $gene->end,
				-strand        => $gene->strand,
				-phase         => '.',
				-display_name  => $gene->stable_id,
				-primary_id    => $gene->stable_id,
	);
	$gene_sf->add_tag_value('status', $gene->status);
	
	# add additional information
	my $alias = $gene->external_name;
	if (defined $alias) {
		$alias =~ s/ {3,}//g; # strip excessive spaces
		$gene_sf->add_tag_value('Alias', $alias);
	}
	my $note_text = $gene->description;
	if (defined $note_text) {
		$note_text =~ s/ {3,}//g; # strip excessive spaces
		$gene_sf->add_tag_value('Note', $note_text);
	}
	
	# work through the transcript(s)
	# there is likely only one transcript for RNA genes, but just in case....
	foreach my $transcript (@{ $gene->get_all_Transcripts }) {
		
		# generate SeqFeature Transcript object
		my $trnscpt_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $chr,
				-source        => $gene->source || 'ensembl',
				-primary_tag   => $gene->biotype,
				-start         => $transcript->start,
				-end           => $transcript->end,
				-strand        => $transcript->strand,
				-phase         => '.',
				-display_name  => $transcript->stable_id,
				-primary_id    => $transcript->stable_id,
		);
		
		# add alias
		if ($alias) {
			$trnscpt_sf->add_tag_value('Alias', $alias);
		}
		
		# add the exon(s), again, typically expect only one
		my $exons_ref = $transcript->get_all_Exons;
		my $exon_count = 1;
		foreach my $exon ( @{ $exons_ref } ) {
			
			# add the exon as a subfeature
			$trnscpt_sf->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
					-seq_id        => $trnscpt_sf->seq_id,
					-source        => $trnscpt_sf->source,
					-primary_tag   => 'noncoding_exon',
					-start         => $exon->start,
					-end           => $exon->end,
					-strand        => $exon->strand,
					-display_name  => $exon->stable_id,
					-primary_id    => $trnscpt_sf->display_name . ".exon.$exon_count",
				)
			);
			
			# we're giving a unique id based on transcript name appended with 
			# exon and incrementing number
			$exon_count++;
		}
		
		# associate the transcript with the parent gene feature
		$gene_sf->add_SeqFeature($trnscpt_sf);
	}
	
	# print the gene feature and its subfeatures
	$gene_sf->version(3);
	$gff_fh->print( $gene_sf->gff_string(1) . "\n" );
		# the gff_string method is undocumented in the POD, but is a 
		# valid method. Passing 1 should force a recursive action to 
		# print parent and children.
}





__END__

=head1 NAME

get_ensembl_annotation.pl

=head1 SYNOPSIS

get_ensembl_annotation.pl [--options...] --species <text>
  
  Options:
  --species <text>
  --out <filename>
  --(no)chromo
  --(no)protein
  --(no)rna
  --(no)miscrna
  --(no)snrna
  --(no)snorna
  --(no)mirna
  --(no)rrna
  --(no)trna
  --prefix
  --group <text>
  --host <host.address>
  --user <text>
  --pass <text>
  --printdb
  --version
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --species <text>

Enter the species name for the database to connect to. Common aliases, 
e.g. human, may be acceptable. Check the EnsEMBL website for the species 
databases available or print the list using the --printdb option below.

=item --out <filename>

Specify the output filename. By default, it uses the species name. An 
extension will be added if necessary.

=item --(no)chromo

Boolean flag to indicate whether or not to write (or not) features for 
the toplevel chromosomes/scaffolds/contigs/sequences in the database. 
The GFF type is the generic term 'chromosome'. The default is true.

=item --no(protein)

Boolean flag to indicate whether or not to collect protein-coding genes 
from the database. The default is true.

=item --(no)rna

Boolean flag to indicate whether or not to collect all non-coding RNA 
genes, including misc_RNA, snRNA, snoRNA, rRNA, tRNA, miRNA, and piRNA, 
from the database. This option may be superseded by setting the 
individual RNA options. For example, setting both --rna and --notrna 
will collect all RNA types except tRNAs. The default is false. 

=item --(no)miscrna

Boolean flag to indicate whether or not to collect miscellenaeous 
noncoding RNA genes.

=item --(no)snrna

Boolean flag to indicate whether or not to collect small nuclear 
RNA (snRNA) genes.

=item --(no)snorna

Boolean flag to indicate whether or not to collect small nucleolar 
RNA (snoRNA) genes.

=item --(no)mirna

Boolean flag to indicate whether or not to collect micro 
RNA (miRNA) genes.

=item --(no)rrna

Boolean flag to indicate whether or not to collect ribosomal RNA 
(rRNA) genes.

=item --(no)trna

Boolean flag to indicate whether or not to collect transfer RNA   
(tRNA) genes.

=item --prefix

Boolean flag to prefix chromosome names with 'chr'. Some people and/or 
downstream applications prefer chromosomes to be named as 'chr1' 
instead of '1'. Only coordinate systems of type 'chromosome' are 
changed, not scaffolds, contigs, etc. Default is false.

=item --group <text>

Specify the name of the database group with which to connect. The default 
value is 'core'. See EnsEMBL documentation for more information.

=item --host <host.address>

Specify the Internet address of the EnsEMBL public MySQL database host. 
The default value is 'ensembldb.ensembl.org'.

=item --user <text>

Specify the user name to connect as to the EnsEMBL public database. 
The default value is 'anonymous'.

=item --pass <text>

Specify the password to use when connecting to the database. The default 
value is none.

=item --printdb

Print all of the available database names, species, and groups from the 
connected EnsEMBL database. The program will exit after printing the list.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will connect to the public EnsEMBL MySQL database and 
retrieve the latest genome annotation for a given species. It will 
generate a GFF3 annotation file suitable for loading into a Bio::DB 
database or genome browser. The GFF3 features generated are 
multi-level nested gene->mRNA->CDS features. It will optionally 
generate features for the top-level sequences (chromosomes, contigs, 
scaffolds, etc.) and non-coding RNA genes (snRNA, tRNA, rRNA, etc.).

This program requires EnsEMBL's Perl API modules to connect to their public 
MySQL servers. It is not available through CPAN, unfortunately, but you can 
find installation instructions at http://www.ensembl.org/info/docs/api/api_installation.html.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  


