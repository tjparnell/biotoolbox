#!/usr/bin/perl -w

# Test script for Bio::ToolBox::parser modules

use strict;
use Test::More;
use File::Spec;
use FindBin '$Bin';

my $lite = 0;
if (eval {require Bio::SeqFeature::Lite; 1}) {
	$lite = 1;
	plan tests => 224;
}
else {
	plan tests => 135;
}
$ENV{'BIOTOOLBOX'} = File::Spec->catfile($Bin, "Data", "biotoolbox.cfg");

require_ok 'Bio::ToolBox::parser::gff' or 
	BAIL_OUT "Cannot load Bio::ToolBox::parser::gff";
require_ok 'Bio::ToolBox::parser::ucsc' or 
	BAIL_OUT "Cannot load Bio::ToolBox::parser::ucsc";
require_ok 'Bio::ToolBox::Data' or 
	BAIL_OUT "Cannot load Bio::ToolBox::Data";


# with the introduction of Bio::ToolBox::SeqFeature class, we can test the 
# parsers with two different seqfeature classes, or possibly more.
# Thus this test script actually tests two things: the parsers and my seqfeature 
# object class.
my $gfffile = File::Spec->catfile($Bin, "Data", "chrI.gff3");
my $ucscfile = File::Spec->catfile($Bin, "Data", "ensGene.txt");
my $outfile = File::Spec->catfile($Bin, "Data", "tempout.txt");

### Testing with standard BioPerl Bio::SeqFeature::Lite class
if ($lite) {
	test_gff('Bio::SeqFeature::Lite');
	test_ucsc('Bio::SeqFeature::Lite');
}

### Testing with internal Bio::ToolBox::SeqFeature class
test_gff('Bio::ToolBox::SeqFeature');
test_ucsc('Bio::ToolBox::SeqFeature');


### Testing Data with gene table parsing
test_parsed_gff_table();
test_parsed_ucsc_table();


sub test_gff {
	my $sfclass = shift;
	print " >> Testing GFF parser with $sfclass\n";
	
	# open and check parser
	my $gff = Bio::ToolBox::parser::gff->new(
		file => $gfffile,
		class => $sfclass,
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
	isa_ok($f, $sfclass, 'first feature object');
	is($f->seq_id, 'chrI', 'feature seq_id');
	is($f->start, 1, 'feature start');
	is($f->stop, 230218, 'feature stop');
	is($f->primary_tag, 'chromosome', 'feature primary_tag');
	is($f->display_name, 'chrI', 'feature display_name');
	is($f->primary_id, 'chrI', 'feature primary_id');
	undef $f;

	# parse next top feature
	$f = $gff->next_top_feature;
	isa_ok($f, $sfclass, 'next top feature object');
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
	isa_ok($f, $sfclass, 'found gene object');
	is($f->start, 42177, 'feature3 start');
	is($f->stop, 42719, 'feature3 stop');
	is($f->primary_tag, 'gene', 'feature3 primary_tag');
	is($f->display_name, 'YAL055W', 'feature3 display_name');
	is($f->primary_id, 'YAL055W', 'feature3 primary_id');
	undef $f;
}




######## UCSC parsing ########

sub test_ucsc {
	# get the seqfeature class to test
	my $sfclass = shift;
	print " >> Testing UCSC parser with $sfclass\n";
	
	# open and check parser
	my $ucsc = Bio::ToolBox::parser::ucsc->new(
		file => $ucscfile,
		class => $sfclass,
	);
	isa_ok($ucsc, 'Bio::ToolBox::parser::ucsc', 'UCSC Parser');
	my $fh = $ucsc->fh;
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
	my $f = $ucsc->next_feature;
	isa_ok($f, $sfclass, 'first gene object');
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
	my $reload = $ucsc->parse_table($ucscfile);
	is($reload, 1, "ucsc parse table");

	# top features
	my @top = $ucsc->top_features;
	is(scalar @top, 5, 'ucsc top features');

	# first top feature
	$f = shift @top;
	isa_ok($f, $sfclass, 'first top gene object');
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
	isa_ok($f, $sfclass, 'first feature object');
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

}


sub test_parsed_gff_table {
	my $Data = Bio::ToolBox::Data->new();
	isa_ok($Data, 'Bio::ToolBox::Data', 'New Data object');
	my $flavor = $Data->taste_file($gfffile);
	is($flavor, 'gff', 'GFF file flavor');
	my $p = $Data->parse_table($gfffile);
	is($p, 1, 'parsed GFF table');
	is($Data->number_columns, 3, 'number of columns');
	is($Data->last_row, 43, 'number of rows');
	is($Data->database, "Parsed:$gfffile", 'database source');
	is($Data->name(0), 'Primary_ID', 'First column name');
	is($Data->name(1), 'Name', 'Second column name');
	is($Data->name(2), 'Type', 'Third column name');
	is($Data->value(4,0), 'YAL069W', 'Fifth row ID');
	is($Data->value(4,1), 'YAL069W', 'Fifth row Name');
	is($Data->value(4,2), 'gene:SGD', 'Fifth row Type');
	my $f = $Data->get_seqfeature(4);
	isa_ok($f, 'Bio::ToolBox::SeqFeature', 'Fifth row SeqFeature object');
	is($f->display_name, 'YAL069W', 'SeqFeature display name');
	is($f->get_tag_values('orf_classification'), 'Dubious', 'SeqFeature attribute');
	is($f->start, 335, 'SeqFeature start position');
	my @subf = $f->get_SeqFeatures;
	is(scalar @subf, 1, 'Number of SeqFeature sub features');
	my $f2 = shift @subf;
	is($f2->type, 'CDS:SGD', 'Sub feature type');
	is($f2->name, 'YAL069W', 'Sub feature name');
	
	# attempt to reload the file
	my $s = $Data->save($outfile);
	is($s, $outfile, 'Saved temporary file');
	undef $Data;
	$Data = Bio::ToolBox::Data->new(
		file    => $outfile,
		parse   => 1,
	);
	isa_ok($Data, 'Bio::ToolBox::Data', 'Reloaded file');
	is($Data->basename, 'tempout', 'file basename');
	my $f3 = $Data->get_seqfeature(4);
	isa_ok($f3, 'Bio::ToolBox::SeqFeature', 'Reloaded fifth row SeqFeature object');
	is($f3->display_name, 'YAL069W', 'SeqFeature display name');
	is($f3->get_tag_values('orf_classification'), 'Dubious', 'SeqFeature attribute');
	unlink $outfile;
}


sub test_parsed_ucsc_table {
	my $Data = Bio::ToolBox::Data->new();
	isa_ok($Data, 'Bio::ToolBox::Data', 'New Data object');
	my $flavor = $Data->taste_file($ucscfile);
	is($flavor, 'ucsc', 'UCSC file flavor');
	my $p = $Data->parse_table($ucscfile);
	is($p, 1, 'parsed UCSC table');
	is($Data->number_columns, 3, 'number of columns');
	is($Data->last_row, 5, 'number of rows');
	is($Data->database, "Parsed:$ucscfile", 'database source');
	is($Data->name(0), 'Primary_ID', 'First column name');
	is($Data->name(1), 'Name', 'Second column name');
	is($Data->name(2), 'Type', 'Third column name');
	is($Data->value(1,0), 'ENSG00000125826', 'First row ID');
	is($Data->value(1,1), 'ENSG00000125826', 'First row Name');
	is($Data->value(1,2), 'gene:EnsGene', 'First row Type');
	my $f = $Data->get_seqfeature(1);
	isa_ok($f, 'Bio::ToolBox::SeqFeature', 'First row SeqFeature object');
	is($f->display_name, 'ENSG00000125826', 'SeqFeature display name');
	is($f->start, 388142, 'SeqFeature start position');
	my @subf = $f->get_SeqFeatures;
	is(scalar @subf, 13, 'Number of SeqFeature sub features');
	my $f2 = shift @subf;
	is($f2->type, 'mRNA:EnsGene', 'Sub feature type');
	is($f2->name, 'ENST00000411647', 'Sub feature name');
}



