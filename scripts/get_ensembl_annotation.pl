#!/usr/bin/perl

# documentation at end of file

use strict;
use IO::File;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqFeature::Lite;
# use Data::Dumper;

# Check for Bio::EnsEMBL
my $bio_ensembl = 0;
eval {
	# API version
	require Bio::EnsEMBL::ApiVersion;
	Bio::EnsEMBL::ApiVersion->import;
	$bio_ensembl = software_version();
	
	# main registry
	require Bio::EnsEMBL::Registry;
	Bio::EnsEMBL::Registry->import;
};
my $VERSION = '1.20';
	
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
	$do_cds,
	$do_utr,
	$do_codon,
	$group,
	$host,
	$port,
	$user,
	$pass,
	$ucsc,
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
	'cds!'      => \$do_cds, # include CDS in output
	'utr!'      => \$do_utr, # include UTRs in output
	'codon!'    => \$do_codon, # include start & stop codons in output
	'group=s'   => \$group, # the database group
	'host=s'    => \$host, # host address
	'port=i'    => \$port, # IP port number
	'user=s'    => \$user, # user name to log in
	'pass=s'    => \$pass, # password to log in with
	'ucsc|prefix!' => \$ucsc, # use UCSC style chromosome naming conventions
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
	printf " Bio::EnsEMBL API version %s\n", $bio_ensembl ? $bio_ensembl : 'not installed';
	exit;
}



# Check EnsEMBL
if ($bio_ensembl) {
	print " Using Bio::EnsEMBL API version $bio_ensembl\n";
}
else {
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
unless ($port) {
	$port = 5306;
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
unless (defined $get_rna_genes) {
	$get_rna_genes = 1;
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
unless (defined $do_cds) {
	$do_cds = 1;
	unless (defined $do_codon) {
		$do_codon = 0;
	}
}
unless (defined $do_utr) {
	$do_utr = 1;
}
my $start_time = time;



### Establish ensembl connection
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
		-host    => $host,
		-user    => $user,
		-pass    => $pass,
		-port    => $port,
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
my $csa = Bio::EnsEMBL::Registry->get_adaptor($species, $group, "coordsystem" );
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

my $top_level_coordsys = $csa->fetch_all('toplevel') if $csa;
my $assembly = $top_level_coordsys->[0]->version || undef;
$gff_fh->print("##genome_build $assembly\n") if $assembly; 






### Begin siphoning data

# Get chromosomes
my $slices_ref = $slice_adaptor->fetch_all('toplevel');
print " Identified " . scalar(@{ $slices_ref }) . " toplevel chromosomes/contigs\n";

# reorder chromosomes
$slices_ref = reorder_slices($slices_ref);

# Process chromosomes
my $phase; # global variable for keeping track of CDS phases
my %id2counts; # global hash for keeping track of unique gene names
foreach my $slice (@{ $slices_ref }) {
	
	print "  Collecting features for chromosome " . 
		$slice->seq_region_name . "\n";
	
	# record how many features we write
	my $chr_feature_count = 0;
	
	# record chromsomes if requested
	if ($get_chromo) {
		
		# set chromosome name
		my $name = $slice->seq_region_name;
		if ($slice->coord_system_name eq 'chromosome' and $ucsc) {
			$name =~ s/^MT$/M/; # rename mitochondrial chromosome
			$name = 'chr' . $name;
		}
		
		# generate the chromosome SeqFeature object
		my $chromo_sf = Bio::SeqFeature::Lite->new(
				-seq_id        => $name,
				-source        => $assembly || $db_name,
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
		print " Wrote $chr_feature_count features for chromosome " . 
			$slice->seq_region_name . "\n";
	}
	
	# debugging limiter
# 	last;
}





### Finish
$gff_fh->close;
printf " Finished in %.1f minutes. Wrote file '$outfile' \n", 
	(time - $start_time)/60;








########################   Subroutines   ###################################


sub reorder_slices {
	my $slices = shift;
	
	# hashes by type of name => slice
	my %chr_numeric2slice;
	my %contig_numeric2slice;
	my %chr_string2slice;
	my %contig_string2slice;
	
	# go through the list
	foreach my $slice (@{ $slices }) {
		my $name = $slice->seq_region_name;
		if ($slice->coord_system_name eq 'chromosome') {
			# chromosomes
			if ($name =~ /(\d+)$/) {
				# numeric
				$chr_numeric2slice{$1} = $slice;
			}
			else {
				# string
				$chr_string2slice{$name} = $slice;
			}
		}
		else {
			# contigs etc
			if ($name =~ /(\d+)$/) {
				# numeric
				$contig_numeric2slice{$1} = $slice;
			}
			else {
				# string
				$contig_string2slice{$name} = $slice;
			}
		}
	}
	
	# sort the slices, chromosome first, then contigs
	# taking numeric ones first, then strings
	my @sorted;
	foreach (sort {$a <=> $b} keys %chr_numeric2slice) {
		push @sorted, $chr_numeric2slice{$_};
	}
	foreach (sort {$a cmp $b} keys %chr_string2slice) {
		push @sorted, $chr_string2slice{$_};
	}
	foreach (sort {$a <=> $b} keys %contig_numeric2slice) {
		push @sorted, $contig_numeric2slice{$_};
	}
	foreach (sort {$a cmp $b} keys %contig_string2slice) {
		push @sorted, $contig_string2slice{$_};
	}
	
	# finished
	return \@sorted;
}


sub process_genes {
	my ($slice, $biotype) = @_;
	
	# get chromosome name
	my $chr = $slice->seq_region_name;
	if ($slice->coord_system_name eq 'chromosome' and $ucsc) {
		$chr =~ s/^MT$/M/; # rename mitochondrial chromosome
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
				-display_name  => $gene->external_name || $gene->stable_id,
				-primary_id    => $gene->stable_id,
	);
	
	# add stable ID as an alias
	if ($gene->external_name) {
		$gene_sf->add_tag_value('Alias', $gene->stable_id);
	}
	
	# get additional information
	$gene_sf->add_tag_value('status', $gene->status);
	my $note_text = $gene->description;
	if (defined $note_text) {
		$note_text =~ s/ {2,}/ /g; # strip excessive spaces
		$gene_sf->add_tag_value('Note', $note_text);
	}
	
	# get external db cross reference
	my $xref = $gene->display_xref;
	if ($xref) {
		my $db = $xref->dbname;
		my $id = $xref->primary_id;
		$gene_sf->add_tag_value('Dbxref', "$db:$id");
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
				-display_name  => $transcript->external_name || $transcript->stable_id,
				-primary_id    => $transcript->stable_id,
		);
		
		# add stable ID as an alias
		if ($transcript->external_name) {
			$trnscpt_sf->add_tag_value('Alias', $transcript->stable_id);
		}
		
		# use the parent gene common name as an alias
		$trnscpt_sf->add_tag_value('Alias', $gene_sf->display_name);
		
		# get external db cross reference
		my $xref = $transcript->display_xref;
		if ($xref) {
			my $db = $xref->dbname;
			my $id = $xref->primary_id;
			$trnscpt_sf->add_tag_value('Dbxref', "$db:$id");
		}
	
		# add the exons
		my $exons_ref = $transcript->get_all_Exons;
		add_exons($trnscpt_sf, $exons_ref);
		
		# get transcription start/stop
		my $coding_start = $transcript->coding_region_start || undef;
		my $coding_stop  = $transcript->coding_region_end || undef;
		$phase = 0; # reset phase to zero for all new transcripts
		
		# check whether we have coding potential or not
		if (!defined $coding_start and !defined $coding_stop) {
			# do not need to write CDS, UTRs, or codons
			# change the transcript primary_tag to ncRNA
			$trnscpt_sf->primary_tag('ncRNA');
		}
		else {
			# process CDS and UTRs
			if ($do_cds or $do_utr) {
				add_cds_and_utrs($trnscpt_sf, $exons_ref, $coding_start, $coding_stop);
			}
		
			# add codons
			if ($do_codon) {
				add_codons($trnscpt_sf, $exons_ref, $coding_start, $coding_stop);
			}
		}
		
		# add the transcript
		$gene_sf->add_SeqFeature($trnscpt_sf);
	}
	
	
	# print the gene feature and its subfeatures
	$gene_sf->version(3);
	$gff_fh->print( $gene_sf->gff_string(1), "\n###\n");
		# the gff_string method is undocumented in the POD, but is a 
		# valid method. Passing 1 should force a recursive action to 
		# print parent and children.
		# also print the close directive
}


sub add_exons {
	my ($transcript, $exons_ref) = @_;
	
	# exon counter
	# we need to make the primary ids unique, and since exons may be 
	# shared between multiple transcripts, we can't use the exon's id
	# we'll use the transcript id plus an incrementing number to make unique
	# it will take considerable effort to write unified exons with multiple parentage
	# perhaps if I was writing GFF3 directly without going through SeqFeature objects....
	my $exon_count = 1;
	foreach my $exon ( @{ $exons_ref } ) {
		
		# add the exon as a subfeature
		$transcript->add_SeqFeature( 
			Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'exon',
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => $exon->strand,
				-display_name  => $exon->stable_id,
				-primary_id    => $transcript->primary_id . ".exon$exon_count",
			)
		);
		
		# we're giving a unique id based on transcript name appended with 
		# exon and incrementing number
		$exon_count++;
	}
}


sub add_codons {

	my ($transcript, $exons_ref, $coding_start, $coding_stop) = @_;
	
	# in some situations, there may be protein_coding gene transcript
	# without a coding start or stop !?
	if (defined $coding_start and defined $coding_stop) {
		
		# generate the start and stop codons
		my ($start_codon_start, $start_codon_stop, $stop_codon_start, $stop_codon_stop);
		if ($transcript->strand == 1) {
			# forward strand
			$start_codon_start = $coding_start;
			$start_codon_stop  = $coding_start + 2;
			$stop_codon_start  = $coding_stop -2;
			$stop_codon_stop   = $coding_stop;
		}
		else {
			# reverse strand
			$start_codon_start = $coding_stop -2;
			$start_codon_stop  = $coding_stop;
			$stop_codon_start  = $coding_start;
			$stop_codon_stop   = $coding_start + 2;
		} 
		
		# start codon
		$transcript->add_SeqFeature( 
			Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'start_codon',
				-start         => $start_codon_start,
				-end           => $start_codon_stop,
				-strand        => $transcript->strand,
				-phase         => 0,
				-display_name  => $transcript->primary_id . '.start_codon',
				-primary_id    => $transcript->primary_id . '.start_codon',
			)
		);
		
		# stop codon
		$transcript->add_SeqFeature( 
			Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'stop_codon',
				-start         => $stop_codon_start,
				-end           => $stop_codon_stop,
				-strand        => $transcript->strand,
				-phase         => 0,
				-display_name  => $transcript->primary_id . '.stop_codon',
				-primary_id    => $transcript->primary_id . '.stop_codon',
			)
		);
	}
}


sub add_cds_and_utrs {

	my ($transcript, $exons_ref, $cdsStart, $cdsStop) = @_;
	
	# exon counter
	# we need to make the primary ids unique, and since exons may be 
	# shared between multiple transcripts, we can't use the exon's id
	# we'll use the transcript id plus an incrementing number to make unique
	my $ex_count = 1; 
	foreach my $exon ( @{ $exons_ref } ) {
		
		# we need to determine whether the exon is UTR, CDS, or split both
		
		#### 5'UTR, CDS, and 3'UTR all in one exon
		if (
			$exon->start < $cdsStart 
			and
			$exon->end > $cdsStop
		) {
			# all three subfeatures are in this exon
			# we need to make two UTRs and one CDS
			
			# build the left UTR object
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new( 
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $cdsStart - 1,
				-strand        => $transcript->strand,
				-phase         => '.',
				-primary_tag   => $transcript->strand == 1 ? 
									'five_prime_UTR' :
									'three_prime_UTR',
				-primary_id    => $transcript->strand == 1 ? 
									$transcript->primary_id . ".5utr$ex_count" :
									$transcript->primary_id . ".3utr$ex_count",
				-display_name  => $transcript->strand == 1 ? 
									$exon->stable_id . ".5utr" :
									$exon->stable_id . ".3utr",
				)
			) if $do_utr;
			
			# build the right UTR object
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $cdsStop + 1,
				-end           => $exon->end,
				-strand        => $transcript->strand,
				-phase         => '.',
				-primary_tag   => $transcript->strand == 1 ? 
									'three_prime_UTR':
									'five_prime_UTR' ,
				-primary_id    => $transcript->strand == 1 ? 
									$transcript->primary_id . ".3utr$ex_count" :
									$transcript->primary_id . ".5utr$ex_count",
				-display_name  => $transcript->strand == 1 ? 
									$exon->stable_id . ".3utr" :
									$exon->stable_id . ".5utr",
				)
			) if $do_utr;
			
			# now build the cds object
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $cdsStart,
				-end           => $cdsStop,
				-strand        => $transcript->strand,
				-phase         => $phase,
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds$ex_count",
				-display_name  => $exon->stable_id . ".cds",
				)
			) if $do_cds;
			
			# do not need to reset the phase
		}
		
		#### UTR only ####
		elsif (
			$exon->start < $cdsStart # cdsStart
			and
			$exon->end < $cdsStart
		) {
			# the exon start/end is entirely before the cdsStart
			# we have either a 5'UTR for forward strand or 3'UTR for reverse strand
			
			# build the utr SeqFeature object
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => $transcript->strand, 
				-phase         => '.',
				-primary_tag   => $transcript->strand == 1 ? 
									'five_prime_UTR' :
									'three_prime_UTR',
				-primary_id    => $transcript->primary_id . ".utr$ex_count",
				-display_name  => $exon->stable_id . ".utr",
				)
			) if $do_utr;
		}
		
		
		#### Split UTR and CDS ####
		elsif (
			$exon->start < $cdsStart 
			and
			$exon->end >= $cdsStart
		) {
			# the start or stop codon is in this exon
			# we need to make two features, 
			# 5'UTR for forward strand, or 3'UTR for reverse strand, plus CDS
			
			# build the utr half of the object
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $cdsStart - 1,
				-strand        => $transcript->strand,
				-phase         => '.',
				-primary_tag   => $transcript->strand == 1 ? 
									'five_prime_UTR' :
									'three_prime_UTR',
				-primary_id    => $transcript->primary_id . ".utr$ex_count",
				-display_name  => $exon->stable_id . '.utr',
				)
			) if $do_utr;
			
			# build the cds half of the object
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $cdsStart,
				-end           => $exon->end,
				-strand        => $transcript->strand,
				-phase         => $phase,
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds$ex_count",
				-display_name  => $exon->stable_id . ".cds",
				)
			) if $do_cds;
			
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
			
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => $transcript->strand,
				-phase         => $phase,
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds$ex_count",
				-display_name  => $exon->stable_id . '.cds',
				)
			) if $do_cds;
			
			# reset phase for next CDS
			# phase + (3 - (length % 3)), readjust to 0..2 if necessary
			# adapted from Barry Moore's gtf2gff3.pl script
			$phase = $phase + (3 - ( ($exon->end - $exon->start + 1) % 3) );
			$phase -=3 if $phase > 2;
		}
		
		
		#### Split CDS and UTR ####
		elsif (
			$exon->start <= $cdsStop 
			and
			$exon->end > $cdsStop
		) {
			# the start or stop codon is in this exon
			# we need to make two features, 
			# 3'UTR for forward strand or 5'UTR for reverse strand, plus CDS
			
			# build the cds half of the object
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $cdsStop,
				-strand        => $transcript->strand,
				-phase         => $phase,
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds$ex_count",
				-display_name  => $exon->stable_id . '.cds',
				)
			) if $do_cds;
			
			# now build the utr half of the object
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $cdsStop + 1,
				-end           => $exon->end,
				-strand        => $transcript->strand,
				-phase         => '.',
				-primary_tag   => $transcript->strand == 1 ? 
									'three_prime_UTR':
									'five_prime_UTR' ,
				-primary_id    => $transcript->primary_id . ".utr$ex_count",
				-display_name  => $exon->stable_id . '.utr',
				)
			) if $do_utr;
			
			# reset phase for next CDS - huh?
			# phase + (3 - (length % 3)), readjust to 0..2 if necessary
			# adapted from Barry Moore's gtf2gff3.pl script
			$phase = $phase + (3 - ( ($cdsStop - $exon->start + 1) % 3) );
			$phase -=3 if $phase > 2;
		}
		
		#### UTR only ####
		elsif (
			$exon->start > $cdsStop 
			and
			$exon->end > $cdsStop
		) {
			# the exon start/end is entirely after the cdsStop
			# we have either a 3'UTR for forward strand or 5'UTR for reverse strand
			
			# build the utr object
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source_tag,
				-start         => $exon->start,
				-end           => $exon->end,
				-strand        => $transcript->strand,
				-phase         => '.',
				-primary_tag   => $transcript->strand == 1 ? 
									'three_prime_UTR':
									'five_prime_UTR' ,
				-primary_id    => $transcript->primary_id . ".utr$ex_count",
				-display_name  => $exon->stable_id . '.utr',
				)
			) if $do_utr;
		}
		
		#### Something's wrong ####
		else {
			# just in case I goofed something up
			warn " programming error! the exon coordinates don't match up " .
				"with CDS coordinates for transcript '" . 
				$transcript->display_name . "'\n " . 
				" Exon coordinates " . $exon->start . ".." .
				$exon->end . "\n" .
				" CDS coordinates $cdsStart .. $cdsStop, strand " . $exon->strand . "\n";
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
				-display_name  => $gene->external_name || $gene->stable_id,
				-primary_id    => $gene->stable_id,
	);
	
	# add stable ID as an alias
	if ($gene->external_name) {
		$gene_sf->add_tag_value('Alias', $gene->stable_id);
	}
	
	# get additional information
	$gene_sf->add_tag_value('status', $gene->status);
	my $note_text = $gene->description;
	if (defined $note_text) {
		$note_text =~ s/ {2,}/ /g; # strip excessive spaces
		$gene_sf->add_tag_value('Note', $note_text);
	}
	
	# get external db cross reference
	my $xref = $gene->display_xref;
	if ($xref) {
		my $db = $xref->dbname;
		my $id = $xref->primary_id;
		$gene_sf->add_tag_value('Dbxref', "$db:$id");
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
		
		# add parent display_name as alias
		$trnscpt_sf->add_tag_value('Alias', $gene_sf->display_name);
		
		# get external db cross reference
		my $xref = $transcript->display_xref;
		if ($xref) {
			my $db = $xref->dbname;
			my $id = $xref->primary_id;
			$trnscpt_sf->add_tag_value('Dbxref', "$db:$id");
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
	$gff_fh->print( $gene_sf->gff_string(1), "\n###\n");
		# the gff_string method is undocumented in the POD, but is a 
		# valid method. Passing 1 should force a recursive action to 
		# print parent and children.
		# also print the close directive
}





__END__

=head1 NAME

get_ensembl_annotation.pl

A script to retrieve Ensembl annotation and write it out as GFF3.

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
  --(no)cds
  --(no)utr
  --codon
  --ucsc
  --group <text>
  --host <host.address>
  --port <integer>
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
e.g. human, may be acceptable. This is typically provided in the format 
as "genus_species", e.g. "homo_sapiens". Check the EnsEMBL website for 
the species databases available, or print the list using the --printdb 
option below.

=item --out <filename>

Specify the output filename. By default, it uses the species name. An 
extension will be added if necessary.

=item --(no)chromo

Boolean flag to indicate whether to write (or not) features for 
the toplevel chromosomes/scaffolds/contigs/sequences in the database. 
The default is true.

=item --no(protein)

Boolean flag to indicate whether or not to collect protein-coding genes 
from the database. The default is true.

=item --(no)rna

Boolean flag to indicate whether or not to collect all non-coding RNA 
genes, including misc_RNA, snRNA, snoRNA, rRNA, tRNA, miRNA, and piRNA, 
from the database. This option may be superseded by setting the 
individual RNA options. For example, setting both --rna and --notrna 
will collect all RNA types except tRNAs. The default is true. 

=item --(no)miscrna

Boolean flag to indicate whether or not to collect miscellenaeous 
noncoding RNA genes. The default is true. 

=item --(no)snrna

Boolean flag to indicate whether or not to collect small nuclear 
RNA (snRNA) genes. The default is true. 

=item --(no)snorna

Boolean flag to indicate whether or not to collect small nucleolar 
RNA (snoRNA) genes. The default is true. 

=item --(no)mirna

Boolean flag to indicate whether or not to collect micro 
RNA (miRNA) genes. The default is true. 

=item --(no)rrna

Boolean flag to indicate whether or not to collect ribosomal RNA 
(rRNA) genes. The default is true. 

=item --(no)trna

Boolean flag to indicate whether or not to collect transfer RNA   
(tRNA) genes. The default is true. 

=item --(no)cds

Boolean flag to indicate whether or not to include CDS information for 
mRNA transcripts and genes. Default is true.

=item --(no)utr

Boolean flag to indicate whether or not to include UTR features for 
mRNA transcripts and genes. Default is true.

=item --codon

Boolean flag to indicate whether or not to include start and stop codons 
for mRNA transcripts and genes. Default is false.

=item --ucsc

Boolean flag to prefix chromosome names with 'chr' in the style of 
UCSC genome annotation. Only coordinate systems of type 'chromosome' 
are changed, not scaffolds, contigs, etc. Default is false.

=item --group <text>

Specify the name of the database group with which to connect. The default 
value is 'core'. See EnsEMBL documentation for more information.

=item --host <host.address>

Specify the Internet address of the EnsEMBL public MySQL database host. 
The default value is 'ensembldb.ensembl.org'.

=item --port <integer>

Specify the IP port address for the MySQL server. Default is 5306.

=item --user <text>

Specify the user name with which to connect to the EnsEMBL public database. 
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

Note that EnsEMBL releases new Perl API modules with each database 
release. If you do not see the latest genome version (compared to what 
is available on the web), you should update your EnsEMBL Perl modules. 
The API version should be printed at the beginning of execution.

=head1 REQUIREMENTS

This program requires EnsEMBL's Perl API modules to connect to their public 
MySQL servers. It is not available through CPAN, unfortunately, but you can 
find installation instructions at 
L<http://www.ensembl.org/info/docs/api/api_installation.html>. 

=head1 DATABASE ACCESS

If you having difficulties connecting, check the server and port numbers 
at L<http://www.ensembl.org/info/data/mysql.html>.

To connect to the Ensembl Genomes public mysql server rather than the 
default, please specify the host as "mysql.ebi.ac.uk" and port 4157. 
See L<http://www.ensemblgenomes.org/info/data_access> for up to date 
information.

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
