#!/usr/bin/perl -w

# Test script for Bio::ToolBox::GeneTools modules

use strict;
use Test::More;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	plan tests => 85;
	## no critic
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, "Data", "biotoolbox.cfg" );
	## use critic
	use_ok('Bio::ToolBox::Parser::ucsc');
	use_ok( 'Bio::ToolBox::GeneTools', qw(:all) );
}

my $ucscfile = File::Spec->catfile( $Bin, "Data", "ensGene.txt" );
my $enssrc   = File::Spec->catfile( $Bin, "Data", "ensemblSource.txt" );

my $ucsc = Bio::ToolBox::Parser::ucsc->new(
	file    => $ucscfile,
	enssrc  => $enssrc,
	do_gene => 1,
	do_cds  => 1,
	do_exon => 1,
	do_utr  => 1,
);

# parse first feature line
my $parsed = $ucsc->parse_table();
is( $parsed, 1, "parsed table" );
my $gene = $ucsc->next_top_feature;
isa_ok( $gene, 'Bio::ToolBox::SeqFeature', 'first gene object' );
is( $gene->seq_id,       'chr20',           'gene seq_id' );
is( $gene->start,        388142,            'gene start' );
is( $gene->stop,         411610,            'gene stop' );
is( $gene->primary_tag,  'gene',            'gene primary_tag' );
is( $gene->display_name, 'ENSG00000125826', 'gene display_name' );
is( $gene->primary_id,   'ENSG00000125826', 'gene primary_id' );

# gene transcript functions
my @transcripts = get_transcripts($gene);
is( scalar @transcripts, 13, 'get_transcripts method' );

# filter transcript
my $filt_gene1 = filter_transcript_biotype( $gene, 'processed_transcript' );
isa_ok( $filt_gene1, 'Bio::ToolBox::SeqFeature', 'first filtered gene transcript' );
my @filt_gene1_trx = get_transcripts($filt_gene1);
is( scalar @filt_gene1_trx, 2, 'number filtered processed_transcripts' );

# gene exons
my @common_exons = get_common_exons($gene);
is( scalar @common_exons, 0, 'common gene exons' );
my @uncommon_exons = get_uncommon_exons($gene);
is( scalar @uncommon_exons,    14,     'uncommon gene exons' );
is( $uncommon_exons[0]->start, 388697, 'first uncommon gene exon start' );
my @alt_exons = get_alt_exons($gene);
is( scalar @alt_exons,    19,     'alt gene exons' );
is( $alt_exons[0]->start, 388142, 'first alt gene exon start' );

# gene introns
my @common_introns = get_common_introns($gene);
is( scalar @common_introns, 0, 'common gene introns' );
my @uncommon_introns = get_uncommon_introns($gene);
is( scalar @uncommon_introns,    14,     'uncommon gene introns' );
is( $uncommon_introns[0]->start, 389424, 'first uncommon gene intron start' );
my @alt_introns = get_alt_introns($gene);
is( scalar @alt_introns,    8,      'alt gene introns' );
is( $alt_introns[0]->start, 388316, 'first alt gene intron start' );

# collapsed transcript
my $collapsedT = collapse_transcripts(@transcripts);
isa_ok( $collapsedT, 'Bio::ToolBox::SeqFeature', 'collapsed transcript object' );
is( get_transcript_length($collapsedT), 3839, 'collapsed transcript length' );
my @collapsedT_exons = get_exons($collapsedT);
is( scalar @collapsedT_exons, 14, 'collapsed transcript exon number' );
my @collapsedT_cds = get_cds($collapsedT);
is( scalar @collapsedT_cds, 0, 'collapsed transcript cds number' );
my @collapsedT_introns = get_introns($collapsedT);
is( scalar @collapsedT_introns,             13,    'collapsed transcript intron number' );
is( get_cdsStart($collapsedT),              undef, 'collapsed transcript CDS start' );
is( get_cdsEnd($collapsedT),                undef, 'collapsed transcript CDS stop' );
is( get_transcript_cds_length($collapsedT), 0,     'collapsed transcript CDS length' );

# first transcript
my $t1 = shift @transcripts;
isa_ok( $t1, 'Bio::ToolBox::SeqFeature', 'first transcript object' );
is( is_coding($t1),             1,                 'transcript1 is_coding' );
is( $t1->primary_tag,           'mRNA',            'transcript1 primary_tag' );
is( $t1->primary_id,            'ENST00000411647', 'transcript1 primary_id' );
is( get_transcript_length($t1), 593,               'transcript1 get_transcript_length' );
my @t1_exons = get_exons($t1);
is( scalar @t1_exons, 5, 'transcript1 exon number' );
my @t1_cds = get_cds($t1);
is( scalar @t1_cds, 4, 'transcript1 cds number' );
my @t1_introns = get_introns($t1);
is( scalar @t1_introns,             4,      'transcript1 intron number' );
is( get_cdsStart($t1),              389402, 'transcript1 CDS start' );
is( get_cdsEnd($t1),                398466, 'transcript1 CDS stop' );
is( get_transcript_cds_length($t1), 352,    'transcript1 CDS length' );
my @t1_utrs = get_utrs($t1);
is( scalar @t1_utrs,          2,                'transcript UTR number' );
is( $t1_utrs[0]->start,       388142,           'first UTR start' );
is( $t1_utrs[1]->start,       389335,           'second UTR start' );
is( $t1_utrs[1]->end,         389401,           'second UTR end' );
is( $t1_utrs[1]->primary_tag, 'five_prime_UTR', 'second UTR primary tag' );

my @fivep_utrs = get_5p_utrs($t1);
is( scalar @fivep_utrs,    2,      '5 prime UTR number' );
is( $fivep_utrs[0]->start, 388142, '5 prime UTR start' );
my @threep_utrs = get_3p_utrs($t1);
is( scalar @threep_utrs,               0,   '3 prime UTR number' );
is( get_transcript_5p_utr_length($t1), 241, '5 prime UTR length' );
is( get_transcript_3p_utr_length($t1), 0,   '3 prime UTR length' );
is( get_transcript_utr_length($t1),    241, 'combined UTR length' );

# second transcript
my $t2 = shift @transcripts;
isa_ok( $t2, 'Bio::ToolBox::SeqFeature', 'second transcript object' );
is( is_coding($t2),             0,                      'transcript2 is_coding' );
is( $t2->primary_tag,           'processed_transcript', 'transcript2 primary_tag' );
is( $t2->primary_id,            'ENST00000465226',      'transcript2 primary_id' );
is( get_transcript_length($t2), 372, 'transcript2 get_transcript_length' );
my @t2_exons = get_exons($t2);
is( scalar @t2_exons, 2, 'transcript2 exon number' );
my @t2_cds = get_cds($t2);
is( scalar @t2_cds, 0, 'transcript2 cds number' );
my @t2_introns = get_introns($t2);
is( scalar @t2_introns,             1,     'transcript2 intron number' );
is( get_cdsStart($t2),              undef, 'transcript2 CDS start' );
is( get_cdsEnd($t2),                undef, 'transcript2 CDS stop' );
is( get_transcript_cds_length($t2), 0,     'transcript2 CDS length' );
my @t2_utrs = get_utrs($t2);
is( scalar @t2_utrs, 0, 'transcript2 UTR number' );

# third transcript
my $t3 = shift @transcripts;
is( is_coding($t3),             0,                 'transcript3 is_coding' );
is( $t3->primary_tag,           'transcript',      'transcript3 primary_tag' );
is( $t3->primary_id,            'ENST00000382214', 'transcript3 primary_id' );
is( get_transcript_length($t3), 2752,              'transcript3 get_transcript_length' );
my @t3_exons = get_exons($t3);
is( scalar @t3_exons,               12,     'transcript3 exon number' );
is( get_cdsStart($t3),              389402, 'transcript3 CDS start' );
is( get_cdsEnd($t3),                407963, 'transcript3 CDS stop' );
is( get_transcript_cds_length($t3), 1011,   'transcript3 CDS length' );
my @t3_utrs = get_utrs($t3);
is( scalar @t3_utrs, 5, 'transcript3 UTR number' );

@fivep_utrs = get_5p_utrs($t3);
is( scalar @fivep_utrs,    1,      '5 prime UTR number' );
is( $fivep_utrs[0]->start, 388694, '5 prime UTR start' );
@threep_utrs = get_3p_utrs($t3);
is( scalar @threep_utrs,               4,      '3 prime UTR number' );
is( $threep_utrs[0]->start,            407964, '3 prime UTR start' );
is( get_transcript_5p_utr_length($t3), 708,    '5 prime UTR length' );
is( get_transcript_3p_utr_length($t3), 1033,   '3 prime UTR length' );
is( get_transcript_utr_length($t3),    1741,   'combined UTR length' );

# export gene as UCSC refFlat
my $expected_rf = <<END;
ENSG00000125826	ENST00000411647	chr20	+	388141	398466	389401	398466	5	388141,389334,390524,398169,398375,	388315,389423,390669,398263,398466,
ENSG00000125826	ENST00000465226	chr20	+	388183	389395	389395	389395	2	388183,389155,	388315,389395,
ENSG00000125826	ENST00000382214	chr20	+	388693	411610	411610	411610	12	388693,390524,398169,398375,399990,400201,401514,402770,407956,409134,409594,410993,	389423,390669,398263,398574,400112,400375,401650,402882,408136,409233,409738,411610,
ENSG00000125826	ENST00000415942	chr20	+	388696	411610	411610	411610	11	388696,390524,398169,398375,399990,400201,402770,407956,409134,409594,410993,	389423,390669,398263,398574,400112,400375,402882,408136,409233,409738,411610,
ENSG00000125826	ENST00000356286	chr20	+	388696	411610	389401	411074	12	388696,390524,398169,398375,399990,400201,401514,402770,407956,409134,409594,410993,	389423,390669,398263,398574,400112,400375,401675,402882,408136,409233,409738,411610,
ENSG00000125826	ENST00000475269	chr20	+	388721	391408	389401	391206	3	388721,390524,391055,	389423,390669,391408,
ENSG00000125826	ENST00000441733	chr20	+	388790	398467	389401	398467	4	388790,390527,398169,398375,	389423,390669,398263,398467,
ENSG00000125826	ENST00000353660	chr20	+	388814	411610	389382	411074	11	388814,398169,398375,399990,400201,401514,402770,407956,409134,409594,410993,	389423,398263,398574,400112,400375,401675,402882,408136,409233,409738,411610,
ENSG00000125826	ENST00000400245	chr20	+	388948	391408	391408	391408	3	388948,390473,391055,	389423,390669,391408,
ENSG00000125826	ENST00000382181	chr20	+	388956	411610	398463	411074	10	388956,398169,398375,399990,400201,402770,407956,409134,409594,410993,	389423,398263,398574,400112,400375,402882,408136,409233,409738,411610,
ENSG00000125826	ENST00000400247	chr20	+	388986	391408	389382	391206	2	388986,391055,	389423,391408,
ENSG00000125826	ENST00000414880	chr20	+	390527	400286	390527	400286	6	390527,397863,398169,398375,399990,400201,	390669,397986,398263,398574,400112,400286,
ENSG00000125826	ENST00000468272	chr20	+	402797	411610	411610	411610	4	402797,407956,409134,410993,	402882,408136,409738,411610,
END
is( ucsc_string($gene), $expected_rf, 'export gene as refFlat' );

# export gene as GTF
my @gtf = split /\n/, gtf_string($gene);
is( scalar(@gtf), 171, 'number of exported gene gtf lines' );

# just check the first transcript
my $expected_gtf = <<END;
chr20	EnsGene	exon	388142	388315	.	+	.	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
chr20	EnsGene	exon	389335	389423	.	+	.	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
chr20	EnsGene	start_codon	389402	389404	.	+	0	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
chr20	EnsGene	CDS	389402	389423	.	+	0	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
chr20	EnsGene	exon	390525	390669	.	+	.	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
chr20	EnsGene	CDS	390525	390669	.	+	2	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
chr20	EnsGene	exon	398170	398263	.	+	.	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
chr20	EnsGene	CDS	398170	398263	.	+	1	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
chr20	EnsGene	exon	398376	398466	.	+	.	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
chr20	EnsGene	CDS	398376	398466	.	+	0	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
chr20	EnsGene	stop_codon	398464	398466	.	+	0	gene_id "ENSG00000125826"; transcript_id "ENST00000411647"; gene_name "ENSG00000125826"; transcript_name "ENST00000411647"; transcript_biotype "protein_coding";
END
my $check_gtf = join( "\n", splice( @gtf, 0, 11 ) ) . "\n";
is( $check_gtf, $expected_gtf, 'first gtf transcript' );

# export as bed12
my $expected_bed = <<END;
chr20	388141	398466	ENST00000411647	1000	+	389401	398466	0	5	174,89,145,94,91	0,1193,2383,10028,10234
END
my $bed = bed12_string($t1);
is( $bed, $expected_bed, 'export transcript as bed12' );

