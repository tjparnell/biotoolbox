#!/usr/bin/perl -w

# Test script for Bio::ToolBox::GeneTools modules

use strict;
use Test::More;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	plan tests => 79;
}

BEGIN {
	use_ok('Bio::ToolBox::parser::ucsc');
	use_ok('Bio::ToolBox::GeneTools', qw(:all));
}

$ENV{'BIOTOOLBOX'} = File::Spec->catfile($Bin, "Data", "biotoolbox.cfg");
my $ucscfile = File::Spec->catfile($Bin, "Data", "ensGene.txt");

my $ucsc = Bio::ToolBox::parser::ucsc->new(
	file    => $ucscfile,
	do_gene => 1,
	do_cds  => 1,
# 	do_utr  => 1,
);

# parse first feature line
my $parsed = $ucsc->parse_table();
is($parsed, 1, "parsed table");
my $gene = $ucsc->next_top_feature;
isa_ok($gene, 'Bio::ToolBox::SeqFeature', 'first gene object');
is($gene->seq_id, 'chr20', 'gene seq_id');
is($gene->start, 388142, 'gene start');
is($gene->stop, 411610, 'gene stop');
is($gene->primary_tag, 'gene', 'gene primary_tag');
is($gene->display_name, 'ENSG00000125826', 'gene display_name');
is($gene->primary_id, 'ENSG00000125826', 'gene primary_id');

# gene transcript functions
my @transcripts = get_transcripts($gene);
is(scalar @transcripts, 13, 'get_transcripts method');

# gene exons
my @common_exons = get_common_exons($gene);
is(scalar @common_exons, 0, 'common gene exons');
my @uncommon_exons = get_uncommon_exons($gene);
is(scalar @uncommon_exons, 14, 'uncommon gene exons');
is($uncommon_exons[0]->start, 388697, 'first uncommon gene exon start');
my @alt_exons = get_alt_exons($gene);
is(scalar @alt_exons, 19, 'alt gene exons');
is($alt_exons[0]->start, 388142, 'first alt gene exon start');

# gene introns
my @common_introns = get_common_introns($gene);
is(scalar @common_introns, 0, 'common gene introns');
my @uncommon_introns = get_uncommon_introns($gene);
is(scalar @uncommon_introns, 14, 'uncommon gene introns');
is($uncommon_introns[0]->start, 389424, 'first uncommon gene intron start');
my @alt_introns = get_alt_introns($gene);
is(scalar @alt_introns, 8, 'alt gene introns');
is($alt_introns[0]->start, 388316, 'first alt gene intron start');


# collapsed transcript
my $collapsedT = collapse_transcripts(@transcripts);
isa_ok($collapsedT, 'Bio::ToolBox::SeqFeature', 'collapsed transcript object');
is(get_transcript_length($collapsedT), 3839, 'collapsed transcript length');
my @collapsedT_exons = get_exons($collapsedT);
is(scalar @collapsedT_exons, 14, 'collapsed transcript exon number');
my @collapsedT_cds = get_cds($collapsedT);
is(scalar @collapsedT_cds, 0, 'collapsed transcript cds number');
my @collapsedT_introns = get_introns($collapsedT);
is(scalar @collapsedT_introns, 13, 'collapsed transcript intron number');
is(get_cdsStart($collapsedT), undef, 'collapsed transcript CDS start');
is(get_cdsEnd($collapsedT), undef, 'collapsed transcript CDS stop');
is(get_transcript_cds_length($collapsedT), 0, 'collapsed transcript CDS length');


# first transcript
my $t1 = shift @transcripts;
isa_ok($t1, 'Bio::ToolBox::SeqFeature', 'first transcript object');
is(is_coding($t1), 1, 'transcript1 is_coding');
is($t1->primary_tag, 'mRNA', 'transcript1 primary_tag');
is($t1->primary_id, 'ENST00000411647', 'transcript1 primary_id');
is(get_transcript_length($t1), 593, 'transcript1 get_transcript_length');
my @t1_exons = get_exons($t1);
is(scalar @t1_exons, 5, 'transcript1 exon number');
my @t1_cds = get_cds($t1);
is(scalar @t1_cds, 4, 'transcript1 cds number');
my @t1_introns = get_introns($t1);
is(scalar @t1_introns, 4, 'transcript1 intron number');
is(get_cdsStart($t1), 389402, 'transcript1 CDS start');
is(get_cdsEnd($t1), 398466, 'transcript1 CDS stop');
is(get_transcript_cds_length($t1), 352, 'transcript1 CDS length');
my @t1_utrs = get_utrs($t1);
is(scalar @t1_utrs, 2, 'transcript UTR number');
is($t1_utrs[0]->start, 388142, 'first UTR start');
is($t1_utrs[1]->start, 389335, 'second UTR start');
is($t1_utrs[1]->end, 389401, 'second UTR end');
is($t1_utrs[1]->primary_tag, 'five_prime_UTR', 'second UTR primary tag');

my @fivep_utrs = get_5p_utrs($t1);
is(scalar @fivep_utrs, 2, '5 prime UTR number');
is($fivep_utrs[0]->start, 388142, '5 prime UTR start');
my @threep_utrs = get_3p_utrs($t1);
is(scalar @threep_utrs, 0, '3 prime UTR number');
is(get_transcript_5p_utr_length($t1), 241, '5 prime UTR length');
is(get_transcript_3p_utr_length($t1), 0, '3 prime UTR length');
is(get_transcript_utr_length($t1), 241, 'combined UTR length');


# second transcript
my $t2 = shift @transcripts;
isa_ok($t2, 'Bio::ToolBox::SeqFeature', 'second transcript object');
is(is_coding($t2), 0, 'transcript2 is_coding');
is($t2->primary_tag, 'ncRNA', 'transcript2 primary_tag');
is($t2->primary_id, 'ENST00000465226', 'transcript2 primary_id');
is(get_transcript_length($t2), 372, 'transcript2 get_transcript_length');
my @t2_exons = get_exons($t2);
is(scalar @t2_exons, 2, 'transcript2 exon number');
my @t2_cds = get_cds($t2);
is(scalar @t2_cds, 0, 'transcript2 cds number');
my @t2_introns = get_introns($t2);
is(scalar @t2_introns, 1, 'transcript2 intron number');
is(get_cdsStart($t2), undef, 'transcript2 CDS start');
is(get_cdsEnd($t2), undef, 'transcript2 CDS stop');
is(get_transcript_cds_length($t2), 0, 'transcript2 CDS length');
my @t2_utrs = get_utrs($t2);
is(scalar @t2_utrs, 0, 'transcript2 UTR number');



# third transcript
my $t3 = shift @transcripts;
is(is_coding($t3), 1, 'transcript3 is_coding');
is($t3->primary_tag, 'mRNA', 'transcript3 primary_tag');
is($t3->primary_id, 'ENST00000382214', 'transcript3 primary_id');
is(get_transcript_length($t3), 2752, 'transcript3 get_transcript_length');
my @t3_exons = get_exons($t3);
is(scalar @t3_exons, 12, 'transcript3 exon number');
is(get_cdsStart($t3), 389402, 'transcript3 CDS start');
is(get_cdsEnd($t3), 407963, 'transcript3 CDS stop');
is(get_transcript_cds_length($t3), 1011, 'transcript3 CDS length');
my @t3_utrs = get_utrs($t3);
is(scalar @t3_utrs, 5, 'transcript3 UTR number');

@fivep_utrs = get_5p_utrs($t3);
is(scalar @fivep_utrs, 1, '5 prime UTR number');
is($fivep_utrs[0]->start, 388694, '5 prime UTR start');
@threep_utrs = get_3p_utrs($t3);
is(scalar @threep_utrs, 4, '3 prime UTR number');
is($threep_utrs[0]->start, 407964, '3 prime UTR start');
is(get_transcript_5p_utr_length($t3), 708, '5 prime UTR length');
is(get_transcript_3p_utr_length($t3), 1033, '3 prime UTR length');
is(get_transcript_utr_length($t3), 1741, 'combined UTR length');



