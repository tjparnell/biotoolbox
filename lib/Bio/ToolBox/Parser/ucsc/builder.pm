package Bio::ToolBox::parser::ucsc::builder;

use warnings;
use strict;

our $VERSION = '1.70';

sub new {
	my ( $class, $linedata, $ucsc ) = @_;
	my %self;
	my $format;

	# we're identifying the type of table based on the number of columns
	# may not be the best or accurate method, but it generally works for valid formats

	### Extended Gene Prediction Table ###
	if ( scalar @{ $linedata } == 16 ) {

		# an extended gene prediction table, e.g. refGene, ensGene, xenoRefGene
		# as downloaded from the UCSC Table Browser or FTP site
		# includes the bin number as the first column

		# 0  bin
		# 1  name
		# 2  chrom
		# 3  strand
		# 4  txStart
		# 5  txEnd
		# 6  cdsStart
		# 7  cdsEnd
		# 8  exonCount
		# 9  exonStarts
		# 10 exonEnds
		# 11 score
		# 12 name2
		# 13 cdsStartStat
		# 14 cdsEndStat
		# 15 exonFrames

		$format           = 'genePredExt with bin';
		$self{name}       = $linedata->[1];
		$self{chrom}      = $linedata->[2];
		$self{strand}     = $linedata->[3];
		$self{txStart}    = $linedata->[4] + 1;
		$self{txEnd}      = $linedata->[5];
		$self{cdsStart}   = $linedata->[6] + 1;
		$self{cdsEnd}     = $linedata->[7];
		$self{exonCount}  = $linedata->[8];
		$self{exonStarts} = $linedata->[9];
		$self{exonEnds}   = $linedata->[10];
		$self{name2}      = $linedata->[12];

		if (exists $ucsc->{ensembldata}->{ $linedata->[1] } ) {
			if ( $ucsc->{ensembldata}->{ $linedata->[1] }->[0] ) {
				$self{gene_name} = $ucsc->{ensembldata}->{ $linedata->[1] }->[0];
			}
			else {
				$self{gene_name} = $linedata->[12];
			}
		}
		else {
			$self{gene_name}     = $linedata->[12];
		}

		if ( exists $ucsc->{refseqsum}->{ $linedata->[1] } ) {
			$self{note}          = $ucsc->{refseqsum}->{ $linedata->[1] }->[1] || q();
			$self{completeness}  = $ucsc->{refseqsum}->{ $linedata->[1] }->[0] || q();
		}
		else {
			$self{note}          = q();
			$self{completeness}  = q();
		}

		if ( exists $ucsc->{refseqstat}->{ $linedata->[1] } ) {
			$self{status}        = $ucsc->{refseqstat}->{ $linedata->[1] }->[0] || q();
		}
		else {
			$self{status}        = q();
		}

		if ( $linedata->[1] =~ /^N[MR]_\d+/ ) {
			$self{refseq}        = $linedata->[1];
		}
	}
	### Extended Gene Prediction Table ###
	elsif ( scalar @{ $linedata } == 15 ) {

		# an extended gene prediction table, e.g. refGene, ensGene, xenoRefGene
		# without the bin value

		# 0  name
		# 1  chrom
		# 2  strand
		# 3  txStart
		# 4  txEnd
		# 5  cdsStart
		# 6  cdsEnd
		# 7  exonCount
		# 8  exonStarts
		# 9 exonEnds
		# 10 score
		# 11 name2
		# 12 cdsStartStat
		# 13 cdsEndStat
		# 14 exonFrames

		$format           = 'genePredExt';
		$self{name}       = $linedata->[0];
		$self{chrom}      = $linedata->[1];
		$self{strand}     = $linedata->[2];
		$self{txStart}    = $linedata->[3] + 1;
		$self{txEnd}      = $linedata->[4];
		$self{cdsStart}   = $linedata->[5] + 1;
		$self{cdsEnd}     = $linedata->[6];
		$self{exonCount}  = $linedata->[7];
		$self{exonStarts} = $linedata->[8];
		$self{exonEnds}   = $linedata->[9];
		$self{name2}      = $linedata->[11];

		if (exists $ucsc->{ensembldata}->{ $linedata->[0] } ) {
			if ( $ucsc->{ensembldata}->{ $linedata->[0] }->[0] ) {
				$self{gene_name} = $ucsc->{ensembldata}->{ $linedata->[0] }->[0];
			}
			else {
				$self{gene_name} = $linedata->[11];
			}
		}
		else {
			$self{gene_name}     = $linedata->[11];
		}

		if ( exists $ucsc->{refseqsum}->{ $linedata->[0] } ) {
			$self{note}          = $ucsc->{refseqsum}->{ $linedata->[0] }->[1] || q();
			$self{completeness}  = $ucsc->{refseqsum}->{ $linedata->[0] }->[0] || q();
		}
		else {
			$self{note}          = q();
			$self{completeness}  = q();
		}

		if ( exists $ucsc->{refseqstat}->{ $linedata->[0] } ) {
			$self{status}        = $ucsc->{refseqstat}->{ $linedata->[0] }->[0] || q();
		}
		else {
			$self{status}        = q();
		}

		if ( $linedata->[0] =~ /^N[MR]_\d+/ ) {
			$self{refseq}        = $linedata->[0];
		}
	}
	### Known Gene Table ###
	elsif ( scalar @{ $linedata } == 12 ) {

# 0 name	known gene identifier
# 1 chrom	Reference sequence chromosome or scaffold
# 2 strand	+ or - for strand
# 3 txStart	Transcription start position
# 4 txEnd	Transcription end position
# 5 cdsStart	Coding region start
# 6 cdsEnd	Coding region end
# 7 exonCount	Number of exons
# 8 exonStarts	Exon start positions
# 9 exonEnds	Exon end positions
# 10 proteinID	UniProt display ID for Known Genes, UniProt accession or RefSeq protein ID for UCSC Genes
# 11 alignID	Unique identifier for each (known gene, alignment position) pair

		$format           = 'knownGene';
		$self{name}       = $linedata->[0];
		$self{chrom}      = $linedata->[1];
		$self{strand}     = $linedata->[2];
		$self{txStart}    = $linedata->[3] + 1;
		$self{txEnd}      = $linedata->[4];
		$self{cdsStart}   = $linedata->[5] + 1;
		$self{cdsEnd}     = $linedata->[6];
		$self{exonCount}  = $linedata->[7];
		$self{exonStarts} = $linedata->[8];
		$self{exonEnds}   = $linedata->[9];
		$self{name2}      = $linedata->[0];

		if ( exists $ucsc->{kgxref}->{ $linedata->[0] } ) {
			$self{name}         = $ucsc->{kgxref}->{ $linedata->[0] }->[0];
			$self{gene_name}    = 
				$ucsc->{kgxref}->{ $linedata->[0] }->[3] || # geneSymbol
				$ucsc->{kgxref}->{ $linedata->[0] }->[0] || # mRNA id
				$ucsc->{kgxref}->{ $linedata->[0] }->[4] || # refSeq id
				$linedata->[0];                             # ugly default
			$self{note}         = $ucsc->{kgxref}->{ $linedata->[0] }->[6]    || q();
			$self{refseq}       = $ucsc->{kgxref}->{ $linedata->[0] }->[4]    || q();
			$self{spid}         = $ucsc->{kgxref}->{ $linedata->[0] }->[1]    || q();
			$self{spdid}        = $ucsc->{kgxref}->{ $linedata->[0] }->[2]    || q();
			$self{protacc}      = $ucsc->{kgxref}->{ $linedata->[0] }->[5]    || q();
		}
		else {
			$self{gene_name}    = $linedata->[0];                                             # ugly default
			$self{note}         = q();
			$self{refseq}       = q();
		}
		
		if (
			exists $self{refseq}
			and $self{refseq}
			and exists $ucsc->{refseqstat}->{ $self{refseq} }
		) {
			$self{status}       = $ucsc->{refseqstat}->{ $self{refseq} }->[0] || q();
			$self{completeness} = $ucsc->{refseqsum}->{ $self{refseq} }->[0]  || q();
		}

	}
	### refFlat or Gene Prediction Table ###
	elsif ( scalar @{ $linedata } == 11 ) {

		# 0  name2 or gene name
		# 1  name or transcript name
		# 2  chrom
		# 3  strand
		# 4  txStart
		# 5  txEnd
		# 6  cdsStart
		# 7  cdsEnd
		# 8  exonCount
		# 9  exonStarts
		# 10 exonEnds

		$format = 'refFlat';
		$self{gene_name}    = $linedata->[0];
		$self{name2}        = $linedata->[0];
		$self{name}         = $linedata->[1];
		$self{chrom}        = $linedata->[2];
		$self{strand}       = $linedata->[3];
		$self{txStart}      = $linedata->[4] + 1;
		$self{txEnd}        = $linedata->[5];
		$self{cdsStart}     = $linedata->[6] + 1;
		$self{cdsEnd}       = $linedata->[7];
		$self{exonCount}    = $linedata->[8];
		$self{exonStarts}   = $linedata->[9];
		$self{exonEnds}     = $linedata->[10];

		if ( exists $ucsc->{ensembldata}->{ $linedata->[1] } ) {
			if ( $ucsc->{ensembldata}->{ $linedata->[1] }->[0] ) {
				$self{gene_name} = $ucsc->{ensembldata}->{ $linedata->[1] }->[0];
			}
			else {
				$self{gene_name} = $linedata->[0];
			}
		}

		if ( exists $ucsc->{refseqsum}->{ $linedata->[1] } ) {
			$self{note}          = $ucsc->{refseqsum}->{ $linedata->[1] }->[1]  || q();
			$self{completeness}  = $ucsc->{refseqsum}->{ $linedata->[1] }->[0]  || q();
		}
		else {
			$self{note}          = q();
			$self{completeness}  = q();
		}

		if ( exists $ucsc->{refseqstat}->{ $linedata->[1] } ) {
			$self{status}        = $ucsc->{refseqstat}->{ $linedata->[1] }->[0] || q();
		}
		else {
			$self{status}        = q();
		}

		if ( $linedata->[1] =~ /^N[MR]_\d+/ ) {
			$self{refseq}        = $linedata->[1];
		}
	}
	### Gene Prediction Table ###
	elsif ( scalar @{ $linedata } == 10 ) {

		# a simple gene prediction table, e.g. refGene, ensGene, xenoRefGene

		# 0  name
		# 1  chrom
		# 2  strand
		# 3  txStart
		# 4  txEnd
		# 5  cdsStart
		# 6  cdsEnd
		# 7  exonCount
		# 8  exonStarts
		# 9  exonEnds

		$format             = 'genePred';
		$self{name}         = $linedata->[0];
		$self{chrom}        = $linedata->[1];
		$self{strand}       = $linedata->[2];
		$self{txStart}      = $linedata->[3] + 1;
		$self{txEnd}        = $linedata->[4];
		$self{cdsStart}     = $linedata->[5] + 1;
		$self{cdsEnd}       = $linedata->[6];
		$self{exonCount}    = $linedata->[7];
		$self{exonStarts}   = $linedata->[8];
		$self{exonEnds}     = $linedata->[9];
		$self{name2}        = $linedata->[0];       # re-use transcript name
		$self{gene_name}    = $linedata->[0];       # re-use transcript name
		
		if ( exists $ucsc->{refseqsum}->{ $linedata->[0] } ) {
			$self{note}         = $ucsc->{refseqsum}->{ $linedata->[0] }->[1]  || q();
			$self{completeness} = $ucsc->{refseqsum}->{ $linedata->[0] }->[0]  || q();
		}
		else {
			$self{note}         = q();
			$self{completeness} = q();
		}

		if ( exists $ucsc->{refseqstat}->{ $linedata->[0] } ) {
			$self{status}       = $ucsc->{refseqstat}->{ $linedata->[0] }->[0] || q();
		}
		else {
			$self{status}       = q();
		}

		if ( $linedata->[0] =~ /^N[MR]_\d+/ ) {
			$self{refseq} = $linedata->[0];
		}
	}
	else {
		# unrecognized line format
		printf STDERR "ERROR: unrecognized format, line has %d elements\n",
			scalar @{ $linedata };
		return;
	}

	# verify
	my @errors;
	push @errors, 'strand'     unless $self{strand}     =~ /^[\+\-\.]$/;
	push @errors, 'txStart'    unless $self{txStart}    =~ /^\d+$/;
	push @errors, 'txEnd'      unless $self{txEnd}      =~ /^\d+$/;
	push @errors, 'cdsStart'   unless $self{cdsStart}   =~ /^\d+$/;
	push @errors, 'cdsEnd'     unless $self{cdsEnd}     =~ /^\d+$/;
	push @errors, 'exonCount'  unless $self{exonCount}  =~ /^\d+$/;
	push @errors, 'exonStarts' unless $self{exonStarts} =~ /^[\d,]+$/;
	push @errors, 'exonEnds'   unless $self{exonEnds}   =~ /^[\d,]+$/;

	if (@errors) {
		printf STDERR "ERROR: line format for %s has the following errors: %s\n", 
			$format, join(", ", @errors);
		return;
	}

	# fix values
	$self{strand}     = $self{strand} eq '+' ? 1 : $self{strand} eq '-' ? -1 : 0;
	$self{exonStarts} = [ map { $_ += 1 } ( split /,/, $self{exonStarts} ) ];
	$self{exonEnds}   = [ ( split /,/, $self{exonEnds} ) ];

	# Attempt to identify the transcript type
	if ( $self{cdsStart} - 1 == $self{cdsEnd} ) {

		# there appears to be no coding potential when
		# txEnd = cdsStart = cdsEnd
		# if you'll look, all of the exon phases should also be -1

		# otherwise, we may be able to infer some certain
		# types from the gene name
		if ( $self{name2} =~ /^mir/i ) {

			# a noncoding gene whose name begins with mir is likely a micro RNA
			$self{type} = 'miRNA';
		}
		elsif ( $self{name2} =~ /^snr/i ) {

			# a noncoding gene whose name begins with snr is likely a snRNA
			$self{type} = 'snRNA';
		}
		elsif ( $self{name2} =~ /^sno/i ) {

			# a noncoding gene whose name begins with sno is likely a snoRNA
			$self{type} = 'snoRNA';
		}
		else {
			# a generic ncRNA
			$self{type} = 'ncRNA';
		}
	}
	else {
		# the transcript has an identifiable CDS so likely a mRNA
		$self{type} = 'mRNA';
	}

	# add the ucsc object
	$self{ucsc} = $ucsc;

	return bless \%self, $class;
}

sub name {
	return shift->{name};
}

sub name2 {
	return shift->{name2};
}

sub gene_name {
	return shift->{gene_name};
}

sub chrom {
	return shift->{chrom};
}

sub txStart {
	return shift->{txStart};
}

sub txEnd {
	return shift->{txEnd};
}

sub strand {
	return shift->{strand};
}

sub cdsStart {
	return shift->{cdsStart};
}

sub cdsEnd {
	return shift->{cdsEnd};
}

sub exonCount {
	return shift->{exonCount};
}

sub exonStarts {
	return shift->{exonStarts};
}

sub exonEnds {
	return shift->{exonEnds};
}

sub type {
	return shift->{type};
}

sub refseq {
	my $self = shift;
	if ( exists $self->{refseq} and $self->{refseq} ) {
		return $self->{refseq};
	} else {
		return q();
	};
}

sub note {
	my $self = shift;
	if ( exists $self->{note} and $self->{note} ) {
		return $self->{note};
	} else {
		return q();
	};
}

sub status {
	my $self = shift;
	if ( exists $self->{status} and $self->{status} ) {
		return $self->{status};
	} else {
		return q();
	};
}

sub completeness {
	my $self = shift;
	if ( exists $self->{completeness} and $self->{completeness} ) {
		return $self->{completeness};
	} else {
		return q();
	};
}

sub ucsc {
	return shift->{ucsc};
}

sub build_gene {
	my $self = shift;

	# shortcuts
	my $ucsc     = $self->ucsc;
	my $id2count = $ucsc->{id2count};

	# find a pre-existing gene to update or build a new one
	my $gene;
	my $name = $self->gene_name;
	if ( exists $ucsc->{gene2seqf}->{$name} ) {

		# we found a gene with the same name
		# pull out the gene seqfeature(s) array reference
		# there may actually be more than one gene
		my $genes = $ucsc->{gene2seqf}->{$name};

		# check that the current transcript intersects with the gene
		# sometimes we can have two separate transcripts with the
		# same gene name, but located on opposite ends of the chromosome
		# part of a gene family, but unlikely the same gene 200 Mb in
		# length
		foreach my $g (@{ $genes }) {
			if (
				# overlap method borrowed from Bio::RangeI
				$g->seq_id eq $self->chrom
				and ( $g->strand == $self->strand )
				and not( $g->start > $self->txEnd
					or $g->end < $self->txStart )
				)
			{
				# gene and transcript overlap on the same strand
				# we found the intersecting gene
				$gene = $g;
				last;
			}
		}
	}

	# process the gene
	if ($gene) {

		# update coordinates as necessary
		if ( ( $self->txStart ) < $gene->start ) {
			$gene->start( $self->txStart );
		}
		if ( $self->txEnd > $gene->end ) {
			$gene->end( $self->txEnd );
		}
	}
	else {
		# build a new gene object
		$gene = $ucsc->{sfclass}->new(
			-seq_id       => $self->chrom,
			-source       => $ucsc->source,
			-primary_tag  => 'gene',
			-start        => $self->txStart,
			-end          => $self->txEnd,
			-strand       => $self->strand,
			-phase        => '.',
			-display_name => $self->gene_name,
		);

		# Add a unique primary ID
		my $id = $self->name2;
		if ( exists $id2count->{ lc $id } ) {

			# we've encountered this gene ID before

			# then make name unique by appending the count number
			$id2count->{ lc $id } += 1;
			$id .= '.' . $id2count->{ lc $id };
		}
		else {
			# this is the first transcript with this id
			# set the id counter
			$id2count->{ lc $id } = 0;
		}
		$gene->primary_id($id);

		# Add to the gene2seqf hash - this is based on gene name, not unique id
		if ( exists $ucsc->{gene2seqf}{$name} ) {

			# more than one gene with this name, add it to the list
			push @{ $ucsc->{gene2seqf}{$name} }, $gene;
		}
		else {
			# inaugural member for the gene name
			$ucsc->{gene2seqf}{$name} = [$gene];
		}

		# Add an alias
		if ( $self->name2 ne $self->gene_name ) {
			$gene->add_tag_value( 'Alias', $self->name2 );
		}

		# add to the gene count
		$ucsc->{counts}{'gene'} += 1;
	}

	# now build the transcript for the gene
	my $transcript = $self->build_transcript($gene);
	$gene->add_SeqFeature($transcript);
	$transcript->add_tag_value( 'Parent', $gene->primary_id );

	# update extra attributes as necessary
	$self->update_attributes($gene);

	# finished
	return $gene;
}

sub build_transcript {
	my ( $self, $gene ) = @_;    # gene is not required

	# shortcuts
	my $ucsc        = $self->ucsc;
	my $id2count    = $ucsc->{id2count};
	my $ensembldata = $ucsc->{ensembldata};

	# Uniqueify the transcript ID but only if a gene is not present
	my $id = $self->name;
	if ( exists $id2count->{ lc $id } ) {

		# we've encountered this transcript ID before

		# now need to make ID unique by appending a number
		$id2count->{ lc $id } += 1;
		$id .= '.' . $id2count->{ lc $id };
	}
	else {
		# this is the first transcript with this id
		$id2count->{ lc $id } = 0;
	}


	# Generate the transcript SeqFeature object
	my $transcript = $ucsc->{sfclass}->new(
		-seq_id       => $self->chrom,
		-source       => $ucsc->source,
		-primary_tag  => $self->type,
		-start        => $self->txStart,
		-end          => $self->txEnd,
		-strand       => $self->strand,
		-phase        => '.',
		-display_name => $self->name,
		-primary_id   => $id,
	);

	# add gene name as an alias
	if ( $self->gene_name ne $self->name2 ) {
		$transcript->add_tag_value( 'Alias', $self->gene_name );
	}

	# adjust the primary_tag and biotype values as necessary
	if ( exists $ensembldata->{ $self->name } ) {
		my $t = $ensembldata->{ $self->name }->[1] || q();
		if ( $t =~ /protein.coding/xi ) {
			$transcript->primary_tag('mRNA');
			$transcript->add_tag_value( 'biotype', $t );
		}
		elsif ( $t =~ /(?: rna | transcript )/xi ) {
			$transcript->primary_tag($t); # this is redundant????
			$transcript->add_tag_value( 'biotype', $t );
		}
		elsif ($t) {
			$transcript->primary_tag('transcript');
			$transcript->add_tag_value( 'biotype', $t );
		}
	}

	# update extra attributes as necessary
	$self->update_attributes($transcript);

	# add transcript specific attributes
	if ( defined $self->completeness ) {
		$transcript->add_tag_value( 'completeness', $self->completeness );
	}
	if ( defined $self->status ) {
		$transcript->add_tag_value( 'status', $self->status );
	}
	# add the exons
	if ( $ucsc->do_exon ) {
		$self->add_exons( $transcript, $gene );
	}

	# add CDS, UTRs, and codons if necessary
	if ( $self->cdsStart - 1 != $self->cdsEnd ) {

		if ( $ucsc->do_utr ) {
			$self->add_utrs( $transcript, $gene );
		}

		if ( $ucsc->do_codon ) {
			$self->add_codons( $transcript, $gene );
		}

		if ( $ucsc->do_cds ) {
			$self->add_cds($transcript);
		}
	}

	# record the type of transcript
	$ucsc->{counts}{ $transcript->primary_tag } += 1;

	# transcript is complete
	return $transcript;
}

sub update_attributes {
	my ( $self, $seqf ) = @_;

	# add Note if possible
	if ( $self->note ) {
		$self->add_unique_attribute( $seqf, 'Note', $self->note );
	}

	# add refSeq identifier if possible
	if ( $self->refseq ) {
		$self->add_unique_attribute( $seqf, 'Dbxref', 'RefSeq:' . $self->refseq );
	}

	# add SwissProt identifier if possible
	if ( exists $self->{spid} and $self->{spid} ) {
		$self->add_unique_attribute( $seqf, 'Dbxref', 'Swiss-Prot:' . $self->{spid} );
	}

	# add SwissProt display identifier if possible
	if ( exists $self->{spdid} and $self->{spdid} ) {
		$self->add_unique_attribute( $seqf, 'swiss-prot_display_id', $self->{spdid} );
	}

	# add NCBI protein access identifier if possible
	if ( exists $self->{protacc} and $self->{protacc} ) {
		$self->add_unique_attribute( $seqf, 'Dbxref', 'RefSeq:' . $self->{protacc} );
	}
}

sub add_unique_attribute {
	my ( $self, $seqf, $tag, $value ) = @_;

	# look for a pre-existing identical tag value
	my $check = 1;
	foreach ( $seqf->get_tag_values($tag) ) {
		if ( $_ eq $value ) {
			$check = 0;
			last;
		}
	}

	# add it if our value is unique
	$seqf->add_tag_value( $tag, $value ) if $check;
}

sub add_exons {
	my ( $self, $transcript, $gene ) = @_;
	my $ucsc = $self->ucsc;

	# Add the exons
EXON_LOOP:
	for ( my $i = 0; $i < $self->exonCount; $i++ ) {

		# first look for existing
		if ( $ucsc->share and $gene ) {
			my $exon = $self->find_existing_subfeature(
				$gene, 'exon',
				$self->exonStarts->[$i],
				$self->exonEnds->[$i]
			);
			if ($exon) {

				# we found an existing exon to reuse
				# associate with this transcript
				$transcript->add_SeqFeature($exon);
				next EXON_LOOP;
			}
		}

		# transform index for reverse strands
		# this will allow numbering from 5'->3'
		my $number;
		if ( $transcript->strand == 1 ) {

			# forward strand
			$number = $i;
		}
		else {
			# reverse strand
			$number = abs( $i - $self->exonCount + 1 );
		}

		# build the exon seqfeature
		my $exon = $ucsc->{sfclass}->new(
			-seq_id      => $transcript->seq_id,
			-source      => $transcript->source,
			-primary_tag => 'exon',
			-start       => $self->exonStarts->[$i],
			-end         => $self->exonEnds->[$i],
			-strand      => $transcript->strand,
			-primary_id  => $transcript->primary_id . ".exon$number",
		);

		# add name if requested
		if ( $ucsc->do_name ) {
			$exon->display_name( $transcript->display_name . ".exon$number" );
		}

		# associate with transcript
		$transcript->add_SeqFeature($exon);
	}
}

sub add_utrs {
	my ( $self, $transcript, $gene ) = @_;
	my $ucsc = $self->ucsc;

	# we will scan each exon and look for a potential utr and build it
	my @utrs;
UTR_LOOP:
	for ( my $i = 0; $i < $self->exonCount; $i++ ) {

		# transform index for reverse strands
		# this will allow numbering from 5'->3'
		my $number;
		if ( $transcript->strand == 1 ) {

			# forward strand
			$number = $i;
		}
		else {
			# reverse strand
			$number = abs( $i - $self->exonCount + 1 );
		}

		# identify UTRs
		# we will identify by comparing the cdsStart and cdsStop relative
		# to the exon coordinates
		# the primary tag is determined by the exon strand orientation
		my ( $start, $stop, $tag );

		# in case we need to build two UTRs
		my ( $start2, $stop2, $tag2 );

		# Split 5'UTR, CDS, and 3'UTR all on the same exon
		if (    $self->exonStarts->[$i] < $self->cdsStart
			and $self->exonEnds->[$i] > $self->cdsEnd )
		{
			# the CDS is entirely within the exon, resulting in two UTRs
			# on either side of the exon
			# we must build two UTRs

			# the left UTR
			$start = $self->exonStarts->[$i];
			$stop  = $self->cdsStart - 1;
			$tag   = $transcript->strand == 1 ? 'five_prime_UTR' : 'three_prime_UTR';

			# the right UTR
			$start2 = $self->cdsEnd + 1;
			$stop2  = $self->exonEnds->[$i];
			$tag2   = $transcript->strand == 1 ? 'three_prime_UTR' : 'five_prime_UTR';
		}

		# 5'UTR forward, 3'UTR reverse
		elsif ( $self->exonStarts->[$i] < $self->cdsStart
			and $self->exonEnds->[$i] < $self->cdsStart )
		{
			# the exon start/end is entirely before the cdsStart
			$start = $self->exonStarts->[$i];
			$stop  = $self->exonEnds->[$i];
			$tag   = $transcript->strand == 1 ? 'five_prime_UTR' : 'three_prime_UTR';
		}

		# Split 5'UTR & CDS on forward, 3'UTR & CDS
		elsif ( $self->exonStarts->[$i] < $self->cdsStart
			and $self->exonEnds->[$i] >= $self->cdsStart )
		{
			# the start/stop codon is in this exon
			# we need to make the UTR out of a portion of this exon
			$start = $self->exonStarts->[$i];
			$stop  = $self->cdsStart - 1;
			$tag   = $transcript->strand == 1 ? 'five_prime_UTR' : 'three_prime_UTR';
		}

		# CDS only
		elsif ( $self->exonStarts->[$i] >= $self->cdsStart
			and $self->exonEnds->[$i] <= $self->cdsEnd )
		{
			# CDS only exon
			next UTR_LOOP;
		}

		# Split 3'UTR & CDS on forward, 5'UTR & CDS
		elsif ( $self->exonStarts->[$i] <= $self->cdsEnd
			and $self->exonEnds->[$i] > $self->cdsEnd )
		{
			# the stop/start codon is in this exon
			# we need to make the UTR out of a portion of this exon
			$start = $self->cdsEnd + 1;
			$stop  = $self->exonEnds->[$i];
			$tag   = $transcript->strand == 1 ? 'three_prime_UTR' : 'five_prime_UTR';
		}

		# 3'UTR forward, 5'UTR reverse
		elsif ( $self->exonStarts->[$i] > $self->cdsEnd
			and $self->exonEnds->[$i] > $self->cdsEnd )
		{
			# the exon start/end is entirely after the cdsStop
			# we have a 3'UTR
			$start = $self->exonStarts->[$i];
			$stop  = $self->exonEnds->[$i];
			$tag   = $transcript->strand == 1 ? 'three_prime_UTR' : 'five_prime_UTR';
		}

		# Something else?
		else {
			printf STDERR
"Warning: A malformed UTR that doesn't match known criteria: cdsStart %d, cdsEnd %d, exonStart %d, exonEnd %d\n",
				$self->cdsStart, $self->cdsEnd, $self->exonStarts->[$i],
				$self->exonEnds->[$i];
			next UTR_LOOP;
		}

		## Generate the UTR objects
		my $utr;

		# look for existing utr
		if ( $ucsc->share and $gene ) {
			$utr = $self->find_existing_subfeature( $gene, $tag, $start, $stop );
		}

		# otherwise build the UTR object
		unless ($utr) {
			$utr = $ucsc->{sfclass}->new(
				-seq_id      => $transcript->seq_id,
				-source      => $transcript->source,
				-start       => $start,
				-end         => $stop,
				-strand      => $transcript->strand,
				-phase       => '.',
				-primary_tag => $tag,
				-primary_id  => $transcript->primary_id . ".utr$number",
			);
			$utr->display_name( $transcript->display_name . ".utr$number" )
				if $ucsc->do_name;
		}

		# store this utr seqfeature in a temporary array
		push @utrs, $utr;

		# build a second UTR object as necessary
		if ($start2) {
			my $utr2;

			# look for existing utr
			if ( $ucsc->share and $gene ) {
				$utr2 = $self->find_existing_subfeature( $gene, $tag2, $start2, $stop2 );
			}

			# otherwise build the utr
			unless ($utr2) {
				$utr2 = $ucsc->{sfclass}->new(
					-seq_id      => $transcript->seq_id,
					-source      => $transcript->source,
					-start       => $start2,
					-end         => $stop2,
					-strand      => $transcript->strand,
					-phase       => '.',
					-primary_tag => $tag2,
					-primary_id  => sprintf("%s.utr%da", $transcript->primary_id, $number),
				);
				$utr2->display_name( sprintf("%s.utr%da", $transcript->display_name, $number) )
					if $ucsc->do_name;
			}

			# store this utr seqfeature in a temporary array
			push @utrs, $utr2;
		}
	}

	# associate found UTRs with the transcript
	foreach my $utr (@utrs) {
		$transcript->add_SeqFeature($utr);
	}
}

sub add_cds {
	my ( $self, $transcript ) = @_;
	my $ucsc = $self->ucsc;

	# we will NOT collapse CDS features since we cannot guarantee that a shared
	# CDS will have the same phase, since phase is dependent on the translation
	# start

	# we will scan each exon and look for a potential CDS and build it
	my @cdss;
	my $phase = 0;    # initialize CDS phase and keep track as we process CDSs
CDS_LOOP:
	for ( my $i = 0; $i < $self->exonCount; $i++ ) {

		# transform index for reverse strands
		my $j;
		if ( $transcript->strand == 1 ) {

			# forward strand
			$j = $i;
		}
		else {
			# reverse strand
			# flip the index for exon starts and stops so that we
			# always progress 5' -> 3'
			# this ensures the phase is accurate from the start codon
			$j = abs( $i - $self->exonCount + 1 );
		}

		# identify CDSs
		# we will identify by comparing the cdsStart and cdsStop relative
		# to the exon coordinates
		my ( $start, $stop );

		# Split 5'UTR, CDS, and 3'UTR all on the same exon
		if (    $self->exonStarts->[$j] < $self->cdsStart
			and $self->exonEnds->[$j] > $self->cdsEnd )
		{
			# exon contains the entire CDS
			$start = $self->cdsStart;
			$stop  = $self->cdsEnd;
		}

		# 5'UTR forward, 3'UTR reverse
		elsif ( $self->exonStarts->[$j] < $self->cdsStart
			and $self->exonEnds->[$j] < $self->cdsStart )
		{
			# no CDS in this exon
			next CDS_LOOP;
		}

		# Split 5'UTR & CDS on forward, 3'UTR & CDS
		elsif ( $self->exonStarts->[$j] < $self->cdsStart
			and $self->exonEnds->[$j] >= $self->cdsStart )
		{
			# the start/stop codon is in this exon
			# we need to make the CDS out of a portion of this exon
			$start = $self->cdsStart;
			$stop  = $self->exonEnds->[$j];
		}

		# CDS only
		elsif ( $self->exonStarts->[$j] >= $self->cdsStart
			and $self->exonEnds->[$j] <= $self->cdsEnd )
		{
			# entire exon is CDS
			$start = $self->exonStarts->[$j];
			$stop  = $self->exonEnds->[$j];
		}

		# Split 3'UTR & CDS on forward, 5'UTR & CDS
		elsif ( $self->exonStarts->[$j] <= $self->cdsEnd
			and $self->exonEnds->[$j] > $self->cdsEnd )
		{
			# the stop/start codon is in this exon
			# we need to make the CDS out of a portion of this exon
			$start = $self->exonStarts->[$j];
			$stop  = $self->cdsEnd;
		}

		# 3'UTR forward, 5'UTR reverse
		elsif ( $self->exonStarts->[$j] > $self->cdsEnd
			and $self->exonEnds->[$j] > $self->cdsEnd )
		{
			# the exon start/end is entirely after the cdsStop
			# we have entirely 5' or 3'UTR, no CDS
			next CDS_LOOP;
		}

		# Something else?
		else {
			printf STDERR
"Warning: A malformed CDS that doesn't match known criteria: cdsStart %d, cdsEnd %d, exonStart %d, exonEnd %d\n",
				$self->cdsStart, $self->cdsEnd, $self->exonStarts->[$j], 
				$self->exonEnds->[$j];
			next CDS_LOOP;
		}

		# build the CDS object
		my $cds = $ucsc->{sfclass}->new(
			-seq_id       => $transcript->seq_id,
			-source       => $transcript->source,
			-start        => $start,
			-end          => $stop,
			-strand       => $transcript->strand,
			-phase        => $phase,
			-primary_tag  => 'CDS',
			-primary_id   => sprintf("%s.cds%d", $transcript->primary_id, $i),
			-display_name => sprintf("%s.cds%d", $transcript->display_name, $i),
		);

		# the id and name still use $i for labeling to ensure numbering from 0

		# store this utr seqfeature in a temporary array
		push @cdss, $cds;

		# reset the phase for the next CDS
		# phase + (3 - (length % 3)), readjust to 0..2 if necessary
		# adapted from Barry Moore's gtf2gff3.pl script
		$phase = $phase + ( 3 - ( $cds->length % 3 ) );
		$phase -= 3 if $phase > 2;
	}

	# associate found UTRs with the transcript
	foreach my $cds (@cdss) {
		$transcript->add_SeqFeature($cds);
	}
}

sub add_codons {
	my ( $self, $transcript, $gene ) = @_;
	my $ucsc = $self->ucsc;

	# generate the start and stop codons
	my ( $start_codon, $stop_codon );
	if ( $transcript->strand == 1 ) {

		# forward strand

		# share codons if possible
		if ( $ucsc->share and $gene ) {
			$start_codon = $self->find_existing_subfeature( $gene, 'start_codon',
				$self->cdsStart, $self->cdsStart + 2 );
			$stop_codon =
				$self->find_existing_subfeature( $gene, 'stop_codon', $self->cdsEnd - 2,
					$self->cdsEnd );
		}

		# start codon
		unless ($start_codon) {
			$start_codon = $ucsc->{sfclass}->new(
				-seq_id      => $transcript->seq_id,
				-source      => $transcript->source,
				-primary_tag => 'start_codon',
				-start       => $self->cdsStart,
				-end         => $self->cdsStart + 2,
				-strand      => 1,
				-phase       => 0,
				-primary_id  => $transcript->primary_id . '.start_codon',
			);
			$start_codon->display_name( $transcript->display_name . '.start_codon' )
				if $ucsc->do_name;
		}

		# stop codon
		unless ($stop_codon) {
			$stop_codon = $ucsc->{sfclass}->new(
				-seq_id      => $transcript->seq_id,
				-source      => $transcript->source,
				-primary_tag => 'stop_codon',
				-start       => $self->cdsEnd - 2,
				-end         => $self->cdsEnd,
				-strand      => 1,
				-phase       => 0,
				-primary_id  => $transcript->primary_id . '.stop_codon',
			);
			$stop_codon->display_name( $transcript->display_name . '.stop_codon' )
				if $ucsc->do_name;
		}
	}

	else {
		# reverse strand

		# share codons if possible
		if ( $ucsc->share and $gene ) {
			$stop_codon = $self->find_existing_subfeature( $gene, 'stop_codon',
				$self->cdsStart, $self->cdsStart + 2 );
			$start_codon =
				$self->find_existing_subfeature( $gene, 'start_codon', $self->cdsEnd - 2,
					$self->cdsEnd );
		}

		# stop codon
		unless ($stop_codon) {
			$stop_codon = $ucsc->{sfclass}->new(
				-seq_id      => $transcript->seq_id,
				-source      => $transcript->source,
				-primary_tag => 'stop_codon',
				-start       => $self->cdsStart,
				-end         => $self->cdsStart + 2,
				-strand      => -1,
				-phase       => 0,
				-primary_id  => $transcript->primary_id . '.stop_codon',
			);
			$stop_codon->display_name( $transcript->display_name . '.stop_codon' )
				if $ucsc->do_name;
		}

		# start codon
		unless ($start_codon) {
			$start_codon = $ucsc->{sfclass}->new(
				-seq_id       => $transcript->seq_id,
				-source       => $transcript->source,
				-primary_tag  => 'start_codon',
				-start        => $self->cdsEnd - 2,
				-end          => $self->cdsEnd,
				-strand       => -1,
				-phase        => 0,
				-primary_id   => $transcript->primary_id . '.start_codon',
				-display_name => $transcript->primary_id . '.start_codon',
			);
			$start_codon->display_name( $transcript->display_name . '.start_codon' )
				if $ucsc->do_name;
		}
	}

	# associate with transcript
	$transcript->add_SeqFeature($start_codon);
	$transcript->add_SeqFeature($stop_codon);
}

sub find_existing_subfeature {
	my ( $self, $gene, $type, $start, $stop ) = @_;
	return unless $gene;

	# we will try to find a pre-existing subfeature at identical coordinates
	foreach my $transcript ( $gene->get_SeqFeatures() ) {

		# walk through transcripts
		foreach my $subfeature ( $transcript->get_SeqFeatures() ) {

			# walk through subfeatures of transcripts
			if (    $subfeature->primary_tag eq $type
				and $subfeature->start == $start
				and $subfeature->end == $stop )
			{
				# we found a match
				return $subfeature;
			}
		}
	}
	return;
}

1;

__END__

=head1 NAME

Bio::ToolBox::parser::ucsc::builder - Build gene and transcript SeqFeatures from UCSC formats

=head1 USAGE

Don't use this directly. This is an accessory module for the L<Bio::ToolBox::parser::ucsc> 
module.

It will take a line from a UCSC-formatted gene table (refFlat, genePred, knownGene) 
and assemble a hierarchical gene or transcript SeqFeature object.

=head1 METHODS

=over 4

=item new

Pass an ARRAY reference representing the tab-delimited fields from a line 
read from UCSC-formatted file and the L<Bio::ToolBox::parser::ucsc> object
itself. The method of parsing is determined by the number of elements.

=item build_transcript

=item build_gene

Methods for building the representative transcript and gene from the line data.

=item add_exons

=item add_utrs

=item add_cds

=item add_codons

Add relative subfeatures to a transcript.

=item find_existing_subfeature

Find the same subfeature when sharing is enabled.

=item update_attributes

=item add_unique_attribute

Add attributes.

=item name

=item name2

=item gene_name

=item chrom

=item txStart

=item txEnd

=item strand

=item cdsStart

=item cdsEnd

=item exonCount

=item exonStarts

=item exonEnds

=item type

=item refseq

=item note

=item status

=item completeness

=item ucsc

Read access to parsed values in the object.

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

