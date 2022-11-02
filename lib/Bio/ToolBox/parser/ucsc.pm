package Bio::ToolBox::parser::ucsc;
our $VERSION = '1.70';

=head1 NAME

Bio::ToolBox::parser::ucsc - Parser for UCSC genePred, refFlat, etc formats

=head1 SYNOPSIS

  use Bio::ToolBox::Parser;
  my $filename = 'file.refFlat';
  
  my $Parser = Bio::ToolBox::Parser->new(
  	file    => $filename,
  	do_gene => 1,
  	do_exon => 1,
  ) or die "unable to open file!\n";
  # the Parser will taste the file and open the appropriate 
  # subclass parser, ucsc in this case
  
  while (my $feature = $Parser->next_top_feature() ) {
	# each $feature is a parent SeqFeature object, usually a gene
  	printf "%s:%d-%d\n", $f->seq_id, $f->start, $f->end;
	
	# subfeatures such as transcripts, exons, etc are nested within
	my @children = $feature->get_SeqFeatures();
  }

=head1 DESCRIPTION

This is the UCSC specific parser subclass to the L<Bio::ToolBox::Parser>
object, and as such inherits generic methods from the parent. 

This is a parser for converting UCSC-style gene prediction flat file formats into 
BioPerl-style L<Bio::SeqFeatureI> compliant objects, complete with nested objects 
representing transcripts, exons, CDS, UTRs, start- and stop-codons. Full control 
is available on what to parse, e.g. exons on, CDS and codons off. Additional gene 
information can be added by supplying additional tables of information, such as 
common gene names and descriptions, available from the UCSC repository. 

=head2 Table formats supported

Supported files are tab-delimited text files obtained from UCSC and described 
at L<http://genome.ucsc.edu/FAQ/FAQformat.html#format9>. Formats are identified 
by the number of columns, rather than specific file extensions, column name 
headers, or other metadata. Therefore, unmodified tables should only be used 
for correct parsing. Some errors are reported for incorrect lines. Unadulterated 
files can safely be downloaded from L<http://hgdownload.soe.ucsc.edu/downloads.html>.
Files obtained from the UCSC Table Browser can also be used with caution. Files 
may be gzip compressed.

File formats supported include the following.

=over 4

=item * Gene Prediction (genePred), 10 columns

=item * Gene Prediction with RefSeq gene Name (refFlat), 11 columns

=item * Extended Gene Prediction (genePredExt), 15 columns

=item * Extended Gene Prediction with bin (genePredExt), 16 columns

=item * knownGene table, 12 columns

=back

=head2 Supplemental information

The UCSC gene prediction tables include essential information, but not detailed 
information, such as common gene names, description, protein accession IDs, etc. 
This additional information can be associated with the genes or transcripts during 
parsing if the appropriate tables are supplied. These tables can be obtained from 
the UCSC download site L<http://hgdownload.soe.ucsc.edu/downloads.html>.

Supported tables include the following.

=over 4 

=item * refSeqStatus, for refGene, knownGene, and xenoRefGene tables

=item * refSeqSummary, for refGene, knownGene, and xenoRefGene tables

=item * ensemblToGeneName, for ensGene tables

=item * ensemblSource, for ensGene tables

=item * kgXref, for knownGene tables

=back

=head2 Implementation

For an implementation of this module to generate GFF3 formatted files from UCSC 
data sources, see the L<Bio::ToolBox> script L<ucsc_table2gff3.pl>.

=head1 METHODS

=head2 Initalize the parser object

In most cases, users should initialize an object using the generic 
L<Bio::ToolBox::Parser> object.

=over 4

=item new

Initiate a UCSC table parser object. Pass a single value (a table file name) 
to open a table and parse its objects. Alternatively, pass an array of key 
value pairs to control how the table is parsed. Options include the following.

=over 4

=item file

=item table

Provide a file name for a UCSC gene prediction table. The file may be gzip 
compressed. 

=item source

Pass a string to be added as the source tag value of the SeqFeature objects. 
The default value is 'UCSC'. If the file name has a recognizable name, 
such as 'refGene' or 'ensGene', it will be used instead.

=item do_gene

Pass a boolean (1 or 0) value to combine multiple transcripts with the same gene 
name under a single gene object. Default is true.

-item do_exon

=item do_cds

=item do_utr

=item do_codon

Pass a boolean (1 or 0) value to parse certain subfeatures, including exon, 
CDS, five_prime_UTR, three_prime_UTR, stop_codon, and start_codon features. 
Default is false.

=item do_name

Pass a boolean (1 or 0) value to assign names to subfeatures, including exons, 
CDSs, UTRs, and start and stop codons. Default is false.

=item share

Pass a boolean (1 or 0) value to recycle shared subfeatures (exons and UTRs) 
between multiple transcripts of the same gene. This results in reduced 
memory usage, and smaller exported GFF3 files. Default is true. 

=item refseqsum

=item refseqstat

=item kgxref

=item ensembltogene

=item ensemblsource

Pass the appropriate file name for additional information.

=item class

Pass the name of a L<Bio::SeqFeatureI> compliant class that will be used to 
create the SeqFeature objects. The default is to use L<Bio::ToolBox::SeqFeature>, 
which is lighter-weight and consumes less memory. A suitable BioPerl alternative
is L<Bio::SeqFeature::Lite>.

=back

=back

=head2 Other methods

See L<Bio::ToolBox::Parser> for generic methods for accessing the 
features. Below are some specific methods to this subclass.

=over 4

=item load_extra_data($file, $type)

	my $file = 'hg19_refSeqSummary.txt.gz';
	my success = $ucsc->load_extra_data($file, 'summary');

Pass two values, the file name of the supplemental file and the type 
of supplemental data. Values can include the following 

=over 4

=item * refseqstatus or status

=item * refseqsummary or summary

=item * kgxref

=item * ensembltogene or ensname

=item * ensemblsource or enssrc

=back

The number of transcripts with information loaded from the supplemental 
data file is returned.

=item counts

This method will return a hash of the number of genes and RNA types that 
have been parsed.

=item typelist

This method will return a comma-delimited list of the feature types or 
C<primary_tag>s found in the parsed file. If a file has not yet been 
parsed, it will return a generic list of expected (typical) feature 
types. Otherwise, it will return the feature types observed in the 
parsed file.

=back

=head1 SEE ALSO

L<Bio::ToolBox::parser::gff>, L<Bio::ToolBox::parser::bed>, 
L<Bio::ToolBox::SeqFeature>

=cut

use strict;
use Carp qw(carp cluck croak confess);
use base 'Bio::ToolBox::Parser';
use Bio::ToolBox::Data;
use Bio::ToolBox::parser::ucsc::builder;

1;

sub new {
	my $class = shift;
	return $class->SUPER::new(@_);
}

sub open_file {
	my $self     = shift;
	my $filename = shift || $self->file || undef;

	# check file
	# Unlike the gff and bed parsers, the ucsc parser can handle opening
	# multiple files, primarily to recycle various UCSC support files.
	# This is historical precedent from earlier versions, unfortunately.
	# This may unnecessarily complicate and/or inflate memory....
	unless ($filename) {
		cluck("no file name passed!");
		return;
	}

	# Check file format type
	my $filetype = $self->filetype || undef;
	unless ($filetype) {
		( my $flavor, $filetype ) = Bio::ToolBox::Data->taste_file($filename);
		unless ( $flavor eq 'ucsc' ) {
			confess "File is not a UCSC-format file!!! How did we get here?";
		}
		$self->{filetype} = $filetype;
	}
	if ( $filetype eq 'genePred' ) {

		# turn off gene processing for simple genePred which has no gene names
		$self->do_gene(0);
	}

	# The ucsc parser does not have convertor subroutines, unlike the gff and bed
	# parsers. Rather, parsing is done by the number of elements in each line
	# and parsed accordingly in the builder object. Slightly inefficient but
	# functional.

	# Open filehandle object
	my $fh = Bio::ToolBox::Data->open_to_read_fh($filename)
		or croak " cannot open file '$filename'!\n";

	# reset source as necessary
	if ( $filename =~ /ensgene/i and $self->source eq 'UCSC' ) {
		$self->source('EnsGene');
	}
	elsif ( $filename =~ /xenorefgene/i and $self->source eq 'UCSC' ) {
		$self->source('xenoRefGene');
	}
	elsif ( $filename =~ /refgene/i and $self->source eq 'UCSC' ) {
		$self->source('refGene');
	}
	elsif ( $filename =~ /refseq/i and $self->source eq 'UCSC' ) {
		$self->source('refSeq');
	}
	elsif ( $filename =~ /knowngene/i and $self->source eq 'UCSC' ) {
		$self->source('knownGene');
	}

	# Check existing data
	# the reason this parser can handle reading a second file without having to make a
	# new parser object (like gff and bed) is so that we can potentially recycle the
	# extra UCSC information
	if ( $self->fh ) {

		# close existing
		$self->{fh}->close;

		# go ahead and clear out existing data
		$self->{gene2seqf}    = {};
		$self->{top_features} = [];
		$self->{loaded}       = {};
		$self->{id2count}     = {};
		$self->{counts}       = {};
		$self->{'eof'}        = 0;
		$self->{line_count}   = 0;
	}

	# Finish
	$self->{fh} = $fh;
	return 1;
}

sub load_extra_data {
	my ( $self, $file, $type ) = @_;
	unless ($file) {
		cluck "no file name passed!";
		return;
	}

	# check the type
	if ( $type =~ /ensembltogene|ensname/i ) {
		$type = 'ensembltogene';
	}
	elsif ( $type =~ /ensemblsource|enssrc/i ) {
		$type = 'ensemblsource';
	}
	elsif ( $type =~ /refseqstat|status/i ) {
		$type = 'refseqstat';
	}
	elsif ( $type =~ /refseqsum|summary/i ) {
		$type = 'refseqsum';
	}
	elsif ( $type =~ /kgxref/i ) {
		$type = 'kgxref';
	}
	else {
		carp "unknown type '$type' to load extra data";
		return;
	}

	my $fh = Bio::ToolBox::Data->open_to_read_fh($file);
	unless ($fh) {
		carp "unable to open file '$file'! $!";
		return;
	}

	# load ensembl data
	my $count = 0;
	if ( $type =~ /ensembl/ ) {

		# we will store gene name in position 0, and source in position 1
		my $index = $type eq 'ensembltogene' ? 0 : 1;
		while ( my $line = $fh->getline ) {

			# process line
			chomp $line;
			next if ( $line =~ /^#/ );
			my @line_data = split /\t/, $line;
			if ( scalar @line_data != 2 ) {
				carp " file $file doesn't seem right!? Line has "
					. scalar @line_data
					. " elements!\n";
				return;
			}

			# store data into hash
			$self->{'ensembldata'}{ $line_data[0] }->[$index] = $line_data[1];
			$count++;
		}
	}

	# load various refSeq data
	else {
		# we just store the line data based on the gene ID
		# each table has different elements
		# here they are for reference

		### refSeqStatus table
# 0	mrnaAcc	RefSeq gene accession name
# 1	status	Status ('Unknown', 'Reviewed', 'Validated', 'Provisional', 'Predicted', 'Inferred')
# 2	molecule type ('DNA', 'RNA', 'ds-RNA', 'ds-mRNA', 'ds-rRNA', 'mRNA', 'ms-DNA', 'ms-RNA', 'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'ss-DNA', 'ss-RNA', 'ss-snoRNA', 'tRNA', 'cRNA', 'ss-cRNA', 'ds-cRNA', 'ms-rRNA')	values	molecule type

		### refSeqSummary table
# 0	RefSeq mRNA accession
# 1	completeness	FullLength ('Unknown', 'Complete5End', 'Complete3End', 'FullLength', 'IncompleteBothEnds', 'Incomplete5End', 'Incomplete3End', 'Partial')
# 1	summary	 	text	values	Summary comments

		### kgXref table
		# 0	kgID	Known Gene ID
		# 1	mRNA	mRNA ID
		# 2	spID	SWISS-PROT protein Accession number
		# 3	spDisplayID	 SWISS-PROT display ID
		# 4	geneSymbol	Gene Symbol
		# 5	refseq	 RefSeq ID
		# 6	protAcc	 NCBI protein Accession number
		# 7	description	Description

		### ensemblToGeneName table
		# 0 Ensembl transcript ID
		# 1 gene name

		# load the table
		while ( my $line = $fh->getline ) {
			chomp $line;
			next if ( $line =~ /^#/ );
			my @line_data = split /\t/, $line;

			# the unique id should be the first element in the array
			# take it off the array, since it doesn't need to be stored there too
			my $id = shift @line_data;

			# check for duplicate lines
			if ( exists $self->{$type}{$id} ) {
				warn "  $type line for identifier $id exists twice!\n";
				next;
			}

			# store data into hash
			$self->{$type}{$id} = [@line_data];
			$count++;
		}
	}

	$fh->close;
	return $count;
}

sub typelist {
	my $self = shift;
	my @items;
	foreach my $k ( keys %{ $self->{counts} } ) {
		push @items, $k if $self->{counts}{$k} > 0;
	}
	if (@items) {
		return join( ',', @items );
	}
	else {
		# return generic list
		return $self->do_gene ? 'gene,mRNA,ncRNA,exon,CDS' : 'mRNA,ncRNA,exon,CDS';
	}
}

sub next_feature {
	my $self = shift;

	# check that we have an open filehandle
	unless ( $self->fh ) {
		croak("no UCSC file loaded to parse!");
	}
	return if $self->{'eof'};

	while ( my $line = $self->fh->getline ) {
		if ( $line !~ /\w+/ ) {
			$self->{line_count}++;
			next;
		}
		if ( $line =~ /^#/ ) {
			push @{ $self->{comments} }, $line;
			$self->{line_count}++;
			next;
		}
		chomp $line;
		my @linedata = split( "\t", $line );
		my $builder  = Bio::ToolBox::parser::ucsc::builder->new( \@linedata, $self );
		$self->{line_count}++;
		unless ($builder) {

			# builder will print its own error message if fails
			warn " unable to parse line number ", $self->{line_count}, "\n";
			next;
		}

		# generate and return the feature from the line
		if ( $self->do_gene ) {
			return $builder->build_gene;
		}
		else {
			return $builder->build_transcript;
		}
	}

	# presumed end of file
	$self->{'eof'} = 1;
	return;
}

*parse_file = \&parse_table;

sub parse_table {
	my $self = shift;
	if (@_) {
		$self->open_file(shift) or return;
	}
	unless ( $self->fh ) {
		carp "must open a file first!";
		return;
	}
	return if ( $self->{'eof'} );

	#### Main Loop
	printf "  Parsing %s format file....\n", $self->filetype;
	while ( my $feature = $self->next_feature ) {

		# record this top feature
		my $id = $feature->primary_id;

		# the builder should have checked and made this ID unique
		$self->{loaded}{$id} = $feature;

		# check chromosome
		my $s = $feature->seq_id;
		unless ( exists $self->{seq_ids}{$s} ) {
			$self->{seq_ids}{$s} = $feature->end;
		}
		$self->{seq_ids}{$s} = $feature->end if $feature->end > $self->{seq_ids}{$s};
	}

	# after parsing the entire file, add all the loaded features to a sorted
	# array of the top features
	# we pseudo-genomic sort by Schwartzian transform because UCSC tables
	# can't be trusted to be sorted
	push @{ $self->{top_features} }, map { $_->[2] }
		sort { $a->[0] cmp $b->[0] or $a->[1] <=> $b->[1] }
		map { [ $_->seq_id, $_->start, $_ ] } values %{ $self->{loaded} };

	# finished parsing
	return 1;
}

sub counts {
	my $self   = shift;
	my %counts = %{ $self->{counts} };
	return wantarray ? %counts : \%counts;
}

__END__

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

