#!/usr/bin/perl

# a script to get ensEMBL annotation


use strict;
use IO::File;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Registry;
use Bio::SeqFeature::Generic;
use Bio::Annotation::Comment;
use Bio::Tools::GFF;



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
	$host,
	$user,
	$outfile,
	$help
);

# Command line options
GetOptions( 
	'species=s' => \$species, # the species to look up
	'chromo!'   => \$get_chromo, # collect chromosome info
	'protein!'  => \$get_protein_genes, 
	'rna!'      => \$get_rna_genes,
	'host=s'    => \$host, # host address
	'user=s'    => \$user, # user name to log in
	'out=s'     => \$outfile, # name of output file 
	'help'      => \$help # request help
);

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}


### Check for requirements and set defaults
# required
unless ($species) {
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


# set default collection types
unless (defined $get_chromo) {
	$get_chromo = 1;
}
unless (defined $get_protein_genes) {
	$get_protein_genes = 1;
}
unless (defined $get_rna_genes) {
	# leave this false at least until I write the code!!!!
	$get_protein_genes = 0;
}




### Establish connection
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
		-host    => $host,
		-user    => $user,
) or die " Can't connect to registry!\n";

my $slice_adaptor = $registry->get_adaptor($species, 'core', 'slice') or
	die " Can't get slice adaptor!\n";
my $db_connection = $slice_adaptor->dbc;
my $db_name = $db_connection->dbname;

print " Success in connecting to database $db_name\n";



### Open output
# open filehandle
	# we want a handle to write extra comments
	# the GFF output object doesn't handle comments
unless ($outfile =~ /\.gff3?$/) {
	$outfile .= '.gff3';
}
my $gff_fh = new IO::File $outfile, 'w';
unless ($gff_fh) {
	die " unable to open file handle for '$outfile'!\n";
}

# open gff output object
my $gff =Bio::Tools::GFF->new(
	-fh           => $gff_fh,
	-gff_version  => 3,
) or die " Unable to open output file '$outfile'!\n";

# write headers
$gff_fh->print("##gff-version 3\n");
$gff->{'_first'} = 0; # pretend we've already written the first feature
					# this avoids writing two gff declaration lines
$gff_fh->print("# Ensembl data for species $species \n");
$gff_fh->print("# Collected from database $db_name\n");








### Begin siphoning data

# Get chromosomes
my $slices_ref = $slice_adaptor->fetch_all('toplevel');
print " Identified " . scalar(@{ $slices_ref }) . " toplevel chromosomes/contigs\n";

# Process chromosomes
my $phase; # global variable for keeping track of CDS phases
foreach my $slice (@{ $slices_ref }) {
	
	# get basic info
	my $chr = $slice->seq_region_name();
	my $length = $slice->seq_region_length();
	
	# record chromsomes if requested
	if ($get_chromo) {
		
		# generate the chromosome SeqFeature object
		my $chromo_sf = Bio::SeqFeature::Generic->new(
				-seq_id        => $chr,
				-source        => $db_name,
				-primary_tag   => $slice->coord_system_name, 
				-start         => 1,
				-end           => $length,
				-strand        => 0,
				-frame         => '.',
				-display_name  => $chr,
		);
		$chromo_sf->add_tag_value('ID', $chr);
		$chromo_sf->add_tag_value('Name', $chr);
		
		# write the object
		print_my_features_and_subfeatures($chromo_sf);
	}
	
	#### DEBUGGING LIMITER ####
# 	unless ($chr eq '25') {
# 		next;
# 	}
	###########################
	
	# collect the protein_coding genes
	if ($get_protein_genes) {
	
		# retrieve a list of genes
		my $genes_ref = $slice->get_all_Genes_by_type('protein_coding', 'ensembl', 1);
			# all protein coding genes, no analysis name, retrieve all information
		print " Collected " . scalar(@{ $genes_ref }) . 
			" genes from chromosome $chr\n";
		
		# process the genes
		foreach my $gene ( @{ $genes_ref } ) {
			process_coding_gene($gene, $chr);
		}
	}
	
	# collect RNA genes
	if ($get_rna_genes) {
		warn " getting RNA genes not supported yet!\n";
	}
	
	# print directive to close out all previous genes
	$gff_fh->print("###\n"); 
}





### Finish
$gff = undef;
$gff_fh->close;
print " Finished! Wrote output file '$outfile'.\n";








########################   Subroutines   ###################################



sub process_coding_gene {
	my ($gene, $chr) = @_;
	
	# create the SeqFeature gene object
	my $gene_id = $gene->stable_id;
	my $gene_name = $gene->external_name || $gene_id;
	my $gene_sf = Bio::SeqFeature::Generic->new(
				-seq_id        => $chr,
				-source        => $gene->source || 'ensembl',
				-primary_tag   => 'gene',
				-start         => $gene->start,
				-end           => $gene->end,
				-strand        => $gene->strand,
				-frame         => '.',
				-display_name  => $gene_name,
	);
	$gene_sf->add_tag_value('ID', $gene_id);
	$gene_sf->add_tag_value('Name', $gene_name);
	$gene_sf->add_tag_value('status', $gene->status);
	my $gene_note = Bio::Annotation::Comment->new(
		-text  => $gene->description || q(),
	);
	$gene_sf->add_tag_value('Note', $gene_note);
	
	# work through the transcripts
	foreach my $transcript (@{ $gene->get_all_Transcripts }) {
		
		# generate SeqFeature Transcript object
		my $trnscpt_sf = Bio::SeqFeature::Generic->new(
				-seq_id        => $chr,
				-source        => $gene->source || 'ensembl',
				-primary_tag   => 'mRNA',
				-start         => $transcript->start,
				-end           => $transcript->end,
				-strand        => $transcript->strand,
				-frame         => '.',
				-display_name  => $transcript->stable_id,
		);
		$trnscpt_sf->add_tag_value('ID', $transcript->stable_id);
		$trnscpt_sf->add_tag_value('Parent', $gene_id);
		$trnscpt_sf->add_tag_value('Name', $transcript->stable_id);
		
		# get transcription start/stop
		my $coding_start = $transcript->coding_region_start;
		my $coding_stop  = $transcript->coding_region_end;
		$phase = 0; # reset phase to zero for all new transcripts
		
		# add the exons
		my $exons_ref = $transcript->get_all_Exons;
		if ($transcript->strand > 0) {
			# forward strand
			foreach my $exon ( @{ $exons_ref } ) {
				process_forward_exon(
					$trnscpt_sf, $exon, $coding_start, $coding_stop);
			}
		}
		else {
			# reverse strand
			foreach my $exon ( @{ $exons_ref } ) {
				process_reverse_exon(
					$trnscpt_sf, $exon, $coding_start, $coding_stop);
			}
		
		}
		
		# add the transcript
		$gene_sf->add_SeqFeature($trnscpt_sf);
	}
	
	# print the gene feature and its subfeatures
	print_my_features_and_subfeatures($gene_sf);
}


sub process_forward_exon {
	my ($transcript, $exon, $cdsStart, $cdsStop) = @_;
	
	# get the parent id
	my ($parent_tag) = $transcript->get_tag_values('ID');
	
	# we need to determine whether the exon is UTR, CDS, or split both
	my $exon_sf; # the SeqFeature object for the exon
	
	#### 5'UTR only ####
	if (
		$exon->start < $cdsStart # cdsStart
		and
		$exon->end < $cdsStart
	) {
		# the exon start/end is entirely before the cdsStart
		# we have a 5'UTR
		
		# build the utr SeqFeature object
		$exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $exon->start,
			-end           => $exon->end,
			-strand        => 1, # we know this already
			-frame         => '.',
			-primary_tag   => 'five_prime_UTR',
		);
		
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
		my $utr_exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $exon->start,
			-end           => $cdsStart - 1,
			-strand        => 1,
			-frame         => '.',
			-primary_tag   => 'five_prime_UTR',
		);
		# since we're actually building two objects here, we have to 
		# complete the process of building the utr object and associate 
		# it with the transcript object
		# the cds half of the exon will be finished below
		
		# add the utr half id
		$utr_exon_sf->add_tag_value('ID', $exon->stable_id . ":utr");
		
		# add the utr half parent
		$utr_exon_sf->add_tag_value('Parent', $parent_tag);

		# associate add the utr half to the parent transcript
		$transcript->add_SeqFeature($utr_exon_sf);
		
		# now build the cds half of the object
		$exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $cdsStart,
			-end           => $exon->end,
			-strand        => 1,
			-frame         => $phase,
			-primary_tag   => 'CDS',
		);
		
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
		
		$exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $exon->start,
			-end           => $exon->end,
			-strand        => 1,
			-frame         => $phase,
			-primary_tag   => 'CDS',
		);
		
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
		$exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $exon->start,
			-end           => $cdsStop,
			-strand        => 1,
			-frame         => $phase,
			-primary_tag   => 'CDS',
		);
		
		# now build the utr half of the object
		my $utr_exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $cdsStop + 1,
			-end           => $exon->end,
			-strand        => 1,
			-frame         => '.',
			-primary_tag   => 'three_prime_UTR',
		);
		
		# add the utr half id
		$utr_exon_sf->add_tag_value('ID', $exon->stable_id . ":utr");
		
		# add the utr half parent
		$utr_exon_sf->add_tag_value('Parent', $parent_tag);

		# associate add the utr half to the parent transcript
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
		$exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $exon->start,
			-end           => $exon->end,
			-strand        => 1,
			-frame         => '.',
			-primary_tag   => 'three_prime_UTR',
		);
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
	
	
	# Finish the exon object
	if ($exon_sf) {
		# add id
		$exon_sf->add_tag_value('ID', $exon->stable_id);
		
		# add parent
		$exon_sf->add_tag_value('Parent', $parent_tag);
	
		# associate exon with the parent transcript
		$transcript->add_SeqFeature($exon_sf);
	}
	else {
		warn " unable to generate exon seqfeature for " . 
			$transcript->display_name . "!\n";
	}
}





sub process_reverse_exon {
	my ($transcript, $exon, $cdsStart, $cdsStop) = @_;
	
	# get the parent id
	my ($parent_tag) = $transcript->get_tag_values('ID');
	
	# we need to determine whether the exon is UTR, CDS, or split both
	my $exon_sf; # the SeqFeature object for the exon
	
	#### 3'UTR only ####
	if (
		$exon->start < $cdsStart # cdsStart
		and
		$exon->end < $cdsStart
	) {
		# the exon start/end is entirely before the cdsStart
		# we have a 3'UTR
		
		# build the utr object
		$exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $exon->start,
			-end           => $exon->end,
			-strand        => -1,
			-frame         => '.',
			-primary_tag   => 'three_prime_UTR',
		);
		
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
		$exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $cdsStart,
			-end           => $exon->end,
			-strand        => -1,
			-frame         => $phase,
			-primary_tag   => 'CDS',
		);
		
		# now build the utr half of the object
		my $utr_exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $exon->start,
			-end           => $cdsStart - 1,
			-strand        => -1,
			-frame         => '.',
			-primary_tag   => 'three_prime_UTR',
		);
		# since we're actually building two objects here, we have to 
		# complete the process of building the utr object and associate 
		# it with the transcript object
		# the cds half of the exon will be finished below
		
		# add the utr half id
		$utr_exon_sf->add_tag_value('ID', $exon->stable_id . ":utr");
		
		# add the utr half parent
		$utr_exon_sf->add_tag_value('Parent', $parent_tag);

		# associate add the utr half to the parent transcript
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
		
		$exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $exon->start,
			-end           => $exon->end,
			-strand        => -1,
			-frame         => $phase,
			-primary_tag   => 'CDS',
		);
		
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
		my $utr_exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $cdsStop + 1,
			-end           => $exon->end,
			-strand        => -1,
			-frame         => '.',
			-primary_tag   => 'five_prime_UTR',
		);
		
		# add the utr half id
		$utr_exon_sf->add_tag_value('ID', $exon->stable_id . ":utr");
		
		# add the utr half parent
		$utr_exon_sf->add_tag_value('Parent', $parent_tag);

		# associate add the utr half to the parent transcript
		$transcript->add_SeqFeature($utr_exon_sf);
		
		# now build the cds half of the object
		$exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $exon->start,
			-end           => $cdsStop,
			-strand        => -1,
			-frame         => $phase,
			-primary_tag   => 'CDS',
		);
		
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
		$exon_sf = Bio::SeqFeature::Generic->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source_tag,
			-start         => $exon->start,
			-end           => $exon->end,
			-strand        => -1,
			-frame         => '.',
			-primary_tag   => 'five_prime_UTR',
		);
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
	
	# Finish the exon object
	if ($exon_sf) {
		# add id
		$exon_sf->add_tag_value('ID', $exon->stable_id);
		
		# add parent
		$exon_sf->add_tag_value('Parent', $parent_tag);
	
		# associate exon with the parent transcript
		$transcript->add_SeqFeature($exon_sf);
	}
	else {
		warn " unable to generate exon seqfeature for " . 
			$transcript->display_name . "!\n";
	}
}




sub print_my_features_and_subfeatures {
	# Bio::Tools::GFF will not recursively print all subfeatures
	# so we fake it here
	
	my $feat = shift;
	
	# print the passed SeqFeature object
	$gff->write_feature($feat);
	
	# print the subfeatures
	foreach ($feat->get_SeqFeatures) {
		print_my_features_and_subfeatures($_);
	}
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
  --host <host.address>
  --user <text>
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --species <text>

Enter the species name for the database to connect to. Common aliases, 
e.g. human, may be acceptable. Check the Ensembl website for the species 
databases available.

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

Boolean flag to indicate whether or not to collect non-coding RNA genes, 
including misc_RNA, snRNA, snoRNA, rRNA, tRNA, miRNA, and piRNA, from 
the database. The default is false.

=item --host <host.address>

Specify the Internet address of the Ensembl public MySQL database host. 
The default value is 'ensembldb.ensembl.org'.

=item --user <text>

Specify the user name to connect as to the public database. 
The default value is 'anonymous'.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will connect to the public Ensembl MySQL database and 
retrieve the latest genome annotation for a given species. It will 
generate a GFF3 annotation file suitable for loading into a Bio::DB 
database or genome browser. The GFF3 features generated are 
multi-level nested gene->mRNA->CDS features. It will optionally 
generate features for the top-level sequences (chromosomes, contigs, 
scaffolds, etc.) and non-coding RNA genes (snRNA, tRNA, rRNA, etc.).



=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112


