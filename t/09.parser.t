#!/usr/bin/perl -w

# Test script for Bio::ToolBox::parser modules

use strict;
use Test::More;
use File::Spec;
use FindBin '$Bin';

BEGIN {
	plan tests => 91;
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile($Bin, "Data", "biotoolbox.cfg");
}

require_ok 'Bio::ToolBox::parser::gff' or 
	BAIL_OUT "Cannot load Bio::ToolBox::parser::gff";
require_ok 'Bio::ToolBox::parser::ucsc' or 
	BAIL_OUT "Cannot load Bio::ToolBox::parser::ucsc";



######## GFF parsing ########

# open and check parser
my $gff = Bio::ToolBox::parser::gff->new(
	file => File::Spec->catfile($Bin, "Data", "chrI.gff3")
);
isa_ok($gff, 'Bio::ToolBox::parser::gff', 'GFF3 Parser');
my $fh = $gff->fh;
isa_ok($fh, 'IO::File', 'IO filehandle');
is($gff->version, 3, 'GFF version');
my @skips = $gff->skip(qw(repeat_region telomere nucleotide_match binding_site 
	ARS long_terminal_repeat region));
is(scalar @skips, 7, 'skipped items');

# parse first feature line
my $f = $gff->next_feature;
isa_ok($f, 'Bio::SeqFeature::Lite', 'first feature object');
is($f->seq_id, 'chrI', 'feature seq_id');
is($f->start, 1, 'feature start');
is($f->stop, 230218, 'feature stop');
is($f->primary_tag, 'chromosome', 'feature primary_tag');
is($f->display_name, 'chrI', 'feature display_name');
is($f->primary_id, 'chrI', 'feature primary_id');
undef $f;

# parse next top feature
$f = $gff->next_top_feature;
isa_ok($f, 'Bio::SeqFeature::Lite', 'next top feature object');
is($f->seq_id, 'chrI', 'feature2 seq_id');
is($f->start, 335, 'feature2 start');
is($f->stop, 649, 'feature2 stop');
is($f->primary_tag, 'gene', 'feature2 primary_tag');
is($f->display_name, 'YAL069W', 'feature2 display_name');
is($f->primary_id, 'YAL069W', 'feature2 primary_id');
my @subf = $f->get_SeqFeatures;
is(scalar @subf, 1, 'feature2 subfeatures');
is($subf[0]->primary_tag, 'CDS', 'feature2 subfeature primary_tag');
is($subf[0]->display_name, 'YAL069W', 'feature2 subfeature display_name');
undef $f;

# top features
my @tops = $gff->top_features;
is(scalar @tops, 32, 'number of top features');

# find gene
$f = $gff->find_gene('YAL055W');
isa_ok($f, 'Bio::SeqFeature::Lite', 'found gene object');
is($f->start, 42177, 'feature3 start');
is($f->stop, 42719, 'feature3 stop');
is($f->primary_tag, 'gene', 'feature3 primary_tag');
is($f->display_name, 'YAL055W', 'feature3 display_name');
is($f->primary_id, 'YAL055W', 'feature3 primary_id');
undef $f;



######## UCSC parsing ########

# open and check parser
my $ucsc = Bio::ToolBox::parser::ucsc->new(
	file => File::Spec->catfile($Bin, "Data", "ensGene.genePred")
);
isa_ok($ucsc, 'Bio::ToolBox::parser::ucsc', 'UCSC Parser');
$fh = $ucsc->fh;
isa_ok($fh, 'IO::File', 'IO filehandle');
is($ucsc->do_gene, 1, 'ucsc do_gene');
is($ucsc->do_cds, 0, 'ucsc do_cds');
is($ucsc->do_utr, 0, 'ucsc do_utr');
is($ucsc->do_codon, 0, 'ucsc do_codon');
is($ucsc->do_name, 0, 'ucsc do_name');
is($ucsc->share, 1, 'ucsc share');
is($ucsc->source, 'EnsGene', 'ucsc source_tag');

# load extra information
my $source_number = $ucsc->load_extra_data(
	File::Spec->catfile($Bin, "Data", "ensemblSource.txt"),
	'ensemblsource'
);
is($source_number, 17, 'ensemblSource extra data');
my $name_number = $ucsc->load_extra_data(
	File::Spec->catfile($Bin, "Data", "ensemblToGeneName.txt"),
	'ensname'
);
is($name_number, 17, 'ensemblToGeneName extra data');

# parse first feature line
$f = $ucsc->next_feature;
isa_ok($f, 'Bio::SeqFeature::Lite', 'first gene object');
is($f->seq_id, 'chr20', 'gene seq_id');
is($f->start, 388142, 'gene start');
is($f->stop, 398466, 'gene stop');
is($f->primary_tag, 'gene', 'gene primary_tag');
is($f->display_name, 'RBCK1', 'gene display_name');
is($f->primary_id, 'ENSG00000125826', 'gene primary_id');

# transcript from first gene
my @transcripts = $f->get_SeqFeatures;
is(scalar @transcripts, 1, 'gene transcripts');
my $t = shift @transcripts;
is($t->start, 388142, 'transcript start');
is($t->stop, 398466, 'transcript stop');
is($t->primary_tag, 'mRNA', 'transcript primary_tag');
is($t->display_name, 'ENST00000411647', 'transcript display_name');

# first transcript exons
my @exons = sort {$a->start <=> $b} $t->get_SeqFeatures; # make sure in order
is(scalar @exons, 5, 'transcript exons');
my $e = shift @exons;
is($e->start, 388142, 'exon start');
is($e->stop, 388315, 'exon stop');
is($e->primary_tag, 'exon', 'exon primary_tag');
undef $f;
undef $t;
undef $e;

# reload the table
my $reload = $ucsc->parse_table( File::Spec->catfile($Bin, "Data", "ensGene.genePred") );
is($reload, 1, "ucsc parse table");

# top features
my @top = $ucsc->top_features;
is(scalar @top, 5, 'ucsc top features');

# first top feature
$f = shift @top;
isa_ok($f, 'Bio::SeqFeature::Lite', 'first top gene object');
is($f->seq_id, 'chr20', 'gene seq_id');
is($f->start, 388142, 'gene start');
is($f->stop, 411610, 'gene stop');
is($f->primary_tag, 'gene', 'gene primary_tag');
is($f->display_name, 'RBCK1', 'gene display_name');
is($f->primary_id, 'ENSG00000125826', 'gene primary_id');

# transcript from first gene
@transcripts = $f->get_SeqFeatures;
is(scalar @transcripts, 13, 'gene transcripts');
$t = shift @transcripts;
is($t->start, 388142, 'first transcript start');
is($t->stop, 398466, 'first transcript stop');
is($t->primary_tag, 'mRNA', 'first transcript primary_tag');
is($t->display_name, 'ENST00000411647', 'first transcript display_name');
undef $t;

# look at the last transcript
$t = pop @transcripts; 
is($t->start, 402798, 'last transcript start');
is($t->stop, 411610, 'last transcript stop');
is($t->primary_tag, 'retained_intron', 'last transcript primary_tag');
is($t->display_name, 'ENST00000468272', 'last transcript display_name');

# last transcript exons
@exons = sort {$a->start <=> $b} $t->get_SeqFeatures; # make sure in order
is(scalar @exons, 4, 'last transcript exons');
$e = shift @exons;
is($e->start, 402798, 'exon start');
is($e->stop, 402882, 'exon stop');
is($e->primary_tag, 'exon', 'exon primary_tag');
undef $f;
undef $t;
undef $e;

# find a gene
$f = $ucsc->find_gene('Y_RNA');
isa_ok($f, 'Bio::SeqFeature::Lite', 'first feature object');
is($f->seq_id, 'chr20', 'feature seq_id');
is($f->start, 431307, 'feature start');
is($f->stop, 431406, 'feature stop');
is($f->primary_tag, 'gene', 'feature primary_tag');
is($f->display_name, 'Y_RNA', 'feature display_name');
is($f->primary_id, 'ENSG00000206797', 'feature primary_id');
undef $f;

# get summary counts
my %counts = $ucsc->counts;
# foreach (keys %counts) {print "$_ => $counts{$_}\n"}
is(scalar keys %counts, 10, "count hash keys");
is($counts{gene}, 5, "count hash gene number");
is($counts{mrna}, 10, "count hash mRNA number");
is($counts{snrna}, 1, "count hash snRNA number");
is($counts{other}, 5, "count hash other number");




