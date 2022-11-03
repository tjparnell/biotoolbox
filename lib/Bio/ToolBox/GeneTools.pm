package Bio::ToolBox::GeneTools;

use warnings;
use strict;
use Carp qw(carp cluck croak confess);
require Exporter;

our $VERSION = '1.67';

### Export
our @ISA       = qw(Exporter);
our @EXPORT_OK = qw(
	get_exons
	get_alt_exons
	get_common_exons
	get_uncommon_exons
	get_alt_common_exons
	get_introns
	get_alt_introns
	get_common_introns
	get_uncommon_introns
	get_alt_common_introns
	get_transcripts
	collapse_transcripts
	get_transcript_length
	is_coding
	get_cds
	get_cdsStart
	get_cdsEnd
	get_start_codon
	get_stop_codon
	get_transcript_cds_length
	get_utrs
	get_transcript_utr_length
	get_5p_utrs
	get_3p_utrs
	get_transcript_5p_utr_length
	get_transcript_3p_utr_length
	gff_string
	gtf_string
	ucsc_string
	bed_string
	filter_transcript_support_level
	filter_transcript_gencode_basic
	filter_transcript_biotype
);
our %EXPORT_TAGS = (
	all  => \@EXPORT_OK,
	exon => [
		qw(
			get_exons
			get_alt_exons
			get_common_exons
			get_uncommon_exons
			get_alt_common_exons
		)
	],
	intron => [
		qw(
			get_introns
			get_alt_introns
			get_common_introns
			get_uncommon_introns
			get_alt_common_introns
		)
	],
	transcript => [
		qw(
			get_transcripts
			collapse_transcripts
			get_transcript_length
		)
	],
	cds => [
		qw(
			is_coding
			get_cds
			get_cdsStart
			get_cdsEnd
			get_start_codon
			get_stop_codon
			get_transcript_cds_length
		)
	],
	utr => [
		qw(
			get_utrs
			get_5p_utrs
			get_3p_utrs
			get_transcript_utr_length
			get_transcript_5p_utr_length
			get_transcript_3p_utr_length
		)
	],
	export => [
		qw(
			gff_string
			gtf_string
			ucsc_string
			bed_string
		)
	],
	filter => [
		qw(
			filter_transcript_biotype
			filter_transcript_gencode_basic
			filter_transcript_support_level
		)
	]
);

######## Exon Methods

sub get_exons {

	# initialize
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	my @exons;
	my @cdss;
	my @transcripts;

	# go through the subfeatures
	foreach my $subfeat ( $transcript->get_SeqFeatures ) {
		my $type = $subfeat->primary_tag;
		if ( $type =~ /exon/i ) {
			push @exons, $subfeat;
		}
		elsif ( $type =~ /cds|utr|untranslated/i ) {
			push @cdss, $subfeat;
		}
		elsif ( $type =~ /rna|transcript/i ) {
			push @transcripts, $subfeat;
		}
	}

	# check which array we'll use
	# prefer to use actual exon subfeatures, but those may not be defined
	my @list;
	if (@exons) {
		@list = map { $_->[0] }
			sort { $a->[1] <=> $b->[1] }
			map { [ $_, $_->start ] } @exons;
	}
	elsif (@cdss) {

		# duplicate the CDSs as exons
		foreach (@cdss) {
			my $e = $_->duplicate;
			$e->primary_tag('exon');    # reset tag
			$e->phase('.');             # no phase
			push @list, $e;
		}

		# make sure to merge adjacent exons
		@list = map { $_->[0] }         # must sort first
			sort { $a->[1] <=> $b->[1] }
			map { [ $_, $_->start ] } @list;
		for ( my $i = 0; $i < scalar @list; $i++ ) {
			if ( defined( $list[ $i + 1 ] ) ) {
				if ( $list[ $i + 1 ]->start - $list[$i]->end <= 1 ) {

					# need to merge
					$list[$i]->end( $list[ $i + 1 ]->end );
					splice( @list, $i + 1, 1 );    # remove the merged
					$i--;
				}
			}
		}
	}
	elsif (@transcripts) {
		foreach my $t (@transcripts) {

			# there are possibly duplicates in here if there are alternate transcripts
			# should we remove them?
			my @e = get_exons($t);
			push @list, @e;
		}
		@list = map { $_->[0] }
			sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] }
			map { [ $_, $_->start, $_->end ] } @list;
	}
	else {
		# nothing found!
		return;
	}

	return wantarray ? @list : \@list;
}

sub get_alt_exons {
	my $ac_exons = get_alt_common_exons(@_);
	my @alts;
	foreach my $k ( keys %$ac_exons ) {
		next if $k eq 'common';
		next if $k eq 'uncommon';
		push @alts, @{ $ac_exons->{$k} };
	}
	if (@alts) {

		# re-sort in genomic order
		@alts = map { $_->[1] }
			sort { $a->[0] <=> $b->[0] }
			map { [ $_->start, $_ ] } @alts;
	}
	return wantarray ? @alts : \@alts;
}

sub get_common_exons {
	my $ac_exons = get_alt_common_exons(@_);
	return wantarray ? @{ $ac_exons->{common} } : $ac_exons->{common};
}

sub get_uncommon_exons {
	my $ac_exons = get_alt_common_exons(@_);
	return wantarray ? @{ $ac_exons->{uncommon} } : $ac_exons->{uncommon};
}

sub get_alt_common_exons {
	return _get_alt_common_things( 1, @_ );
}

######## Intron Methods

sub get_introns {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	my @introns;

	# find the exons and/or CDSs
	my @exons = get_exons($transcript);
	return unless @exons;
	return if ( scalar(@exons) == 1 );

	# identify the last exon index position
	my $last = scalar(@exons) - 1;

	# forward strand
	if ( $transcript->strand >= 0 ) {

		# each intron is created based on the previous exon
		for ( my $i = 0; $i < $last; $i++ ) {
			my $e = $exons[$i];
			my $i = $e->new(
				-seq_id       => $e->seq_id,
				-start        => $e->end + 1,
				-end          => $exons[ $i + 1 ]->start - 1,   # up to start of next exon
				-strand       => $transcript->strand,
				-primary_tag  => 'intron',
				-source_tag   => $transcript->source_tag,
				-primary_id   => $transcript->display_name . ".intron$i",
				-display_name => $transcript->display_name . ".intron$i",
			);
			push @introns, $i;
		}
	}

	# reverse strand
	else {
		# each intron is created based on the previous exon
		# ordering from 5' to 3' end direction for convenience in naming
		for ( my $i = $last; $i > 0; $i-- ) {
			my $e = $exons[$i];
			my $i = $e->new(
				-seq_id      => $e->seq_id,
				-start       => $exons[ $i - 1 ]->end + 1,              # end of next exon
				-end         => $e->start - 1,
				-strand      => $transcript->strand,
				-primary_tag => 'intron',
				-source_tag  => $transcript->source_tag,
				-primary_id  => $transcript->display_name . ".intron$i",
				-display_name => $transcript->display_name . ".intron$i",
			);
			push @introns, $i;
		}

		# reorder the introns based on start position
		if (@introns) {
			@introns = map { $_->[0] }
				sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] }
				map { [ $_, $_->start, $_->end ] } @introns;
		}
	}

	# finished
	return wantarray ? @introns : \@introns;
}

sub get_alt_introns {
	my $ac_introns = get_alt_common_introns(@_);
	my @alts;
	foreach my $k ( keys %$ac_introns ) {
		next if $k eq 'common';
		next if $k eq 'uncommon';
		push @alts, @{ $ac_introns->{$k} };
	}
	if (@alts) {

		# re-sort in genomic order
		@alts = map { $_->[1] }
			sort { $a->[0] <=> $b->[0] }
			map { [ $_->start, $_ ] } @alts;
	}
	return wantarray ? @alts : \@alts;
}

sub get_common_introns {
	my $ac_introns = get_alt_common_introns(@_);
	return wantarray ? @{ $ac_introns->{common} } : $ac_introns->{common};
}

sub get_uncommon_introns {
	my $ac_introns = get_alt_common_introns(@_);
	return wantarray ? @{ $ac_introns->{uncommon} } : $ac_introns->{uncommon};
}

sub get_alt_common_introns {
	return _get_alt_common_things( 0, @_ );
}

sub _get_alt_common_things {

	# internal subroutine to get either exons or introns
	my $do_exon = shift;    # true for exon, false for intron
	my @transcripts;
	return unless @_;
	if ( scalar @_ == 1 ) {

		# someone passed a gene, get the transcripts
		confess "not a SeqFeature object!" unless ref( $_[0] ) =~ /seqfeature/i;
		@transcripts = get_transcripts( $_[0] );
	}
	elsif ( scalar @_ > 1 ) {

		# presume these are transcripts?
		@transcripts = @_;
	}

	# hash of transcript to things
	my %tx2things = (
		common   => [],
		uncommon => []
	);

	# no transcripts provided?
	return \%tx2things unless @transcripts;

	# only one transcript provided?
	if ( scalar @transcripts == 1 ) {

		# all exons are common by definition
		# 		my $name = $transcripts[0]->display_name;
		my @things =
			$do_exon ? get_exons( $transcripts[0] ) : get_introns( $transcripts[0] );
		$tx2things{common} = \@things;
		return \%tx2things;
	}

	# get things and put them in has based on coordinates
	my %pos2things;
	foreach my $t (@transcripts) {
		my @things = $do_exon ? get_exons($t) : get_introns($t);
		foreach my $e (@things) {
			my $new_e = $e->duplicate;
			$pos2things{ $e->start }{ $e->end }{ $t->display_name } = $new_e;
		}
		$tx2things{ $t->display_name } = [];
	}

	# put things into categories based on commonality
	# associate things with unique transcripts, common, or uncommon sets
	my $trx_number = scalar @transcripts;
	foreach my $s ( sort { $a <=> $b } keys %pos2things ) {    # sort on start
		foreach my $e ( sort { $a <=> $b } keys %{ $pos2things{$s} } ) {    # sort on stop
			my @names = keys %{ $pos2things{$s}{$e} };
			if ( scalar @names == 1 ) {

				# only 1 thing, must be an alternate
				push @{ $tx2things{ $names[0] } }, $pos2things{$s}{$e}{ $names[0] };
			}
			elsif ( scalar @names == $trx_number ) {

				# common to all transcripts, take the first one as example
				push @{ $tx2things{common} }, $pos2things{$s}{$e}{ $names[0] };
			}
			else {
				# common to some but not all transcripts, so uncommon
				push @{ $tx2things{uncommon} }, $pos2things{$s}{$e}{ $names[0] };
			}
		}
	}
	return \%tx2things;
}

######## Transcript Methods

sub get_transcripts {
	my $gene = shift;
	return                             unless $gene;
	confess "not a SeqFeature object!" unless ref $gene =~ /seqfeature/i;
	return $gene if ( $gene->primary_tag =~ /rna|transcript/i );
	my @transcripts;
	my @exons;
	my @other;
	foreach my $subf ( $gene->get_SeqFeatures ) {
		if ( $subf->primary_tag =~ /rna|transcript|\bprocessed/i ) {
			push @transcripts, $subf;
		}
		elsif ( $subf->primary_tag =~ /^(?:cds|exon|\w+codon)$/i ) {
			push @exons, $subf;
		}
		else {
			# wierdo subfeature types like unprocessed_pseudogene
			push @other, $subf;
		}
	}
	if ( not @transcripts and @exons ) {

		# some weirdly formatted annotation files skip the transcript
		# looking at you SGD
		my $transcript = $gene->new(
			-seq_id      => $gene->seq_id,
			-start       => $gene->start,
			-end         => $gene->end,
			-strand      => $gene->strand,
			-primary_tag => 'transcript',
			-source      => $gene->source_tag,
			-name        => $gene->display_name,
			-segments    => \@exons,
		);
		push @transcripts, $transcript;
	}
	elsif ( not @transcripts and not @exons and @other ) {

		# well, what choice do we have?
		# we can assume these are transcripts, because what else could they be?
		@transcripts = @other;
	}
	@transcripts = map { $_->[0] }
		sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] }
		map { [ $_, $_->start, $_->length ] } @transcripts;
	return wantarray ? @transcripts : \@transcripts;
}

sub collapse_transcripts {
	my @transcripts;
	return unless @_;
	my $example = $_[0];    # parent transcript seqfeature to model new transcript on
	if ( scalar @_ == 1 ) {

		# someone passed a gene, get the transcripts
		@transcripts = get_transcripts( $_[0] );
		return unless @transcripts;
		return $transcripts[0] if scalar @transcripts == 1;
	}
	elsif ( scalar @_ > 1 ) {
		@transcripts = @_;
	}

	# get all useable exons to collapse
	my @exons;
	foreach my $t (@transcripts) {
		my @e = get_exons($t);
		push @exons, @e;
	}

	# check that we have exons - weirdo files may just have CDS!!!????
	unless (@exons) {
		foreach my $t (@transcripts) {
			my @e = get_cds($t);
			push @exons, @e;
		}
	}
	return unless (@exons);

	# sort all the exons
	my @sorted = map { $_->[0] }
		sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] }
		map { [ $_, $_->start, $_->end ] } @exons;

	# build new exons from the original - don't keep cruft
	my $next = shift @sorted;
	my @new;
	$new[0] = $next->new(
		-seq_id      => $next->seq_id,
		-start       => $next->start,
		-end         => $next->end,
		-strand      => $next->strand,
		-primary_tag => 'exon',
	);

	# work through remaining exons, adding and merging as necessary
	while (@sorted) {
		$next = shift @sorted;
		my ( $ns, $ne ) = ( $next->start, $next->end );          # new start & end
		my ( $os, $oe ) = ( $new[-1]->start, $new[-1]->end );    # old start & end
		if ( $ns == $os and $ne > $oe ) {

			# same beginning, further end
			$new[-1]->end($ne);
		}
		elsif ( $ns > $os and $ns < $oe and $ne > $oe ) {

			# overlapping start, further end
			$new[-1]->end($ne);
		}
		elsif ( $ns > $oe ) {

			# completely new exon
			push @new,
				$next->new(
					-seq_id      => $next->seq_id,
					-start       => $ns,
					-end         => $ne,
					-strand      => $next->strand,
					-primary_tag => 'exon',
				);
		}

		# all other possibilities we can skip
	}

	# return the assembled transcript
	return $example->new(
		-seq_id      => $example->seq_id,
		-start       => $new[0]->start,
		-end         => $new[-1]->end,
		-strand      => $example->strand,
		-primary_tag => 'transcript',
		-source      => $example->source_tag,
		-name        => $example->display_name,
		-segments    => \@new,
	);
}

sub get_transcript_length {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	if ( $transcript->primary_tag =~ /gene$/i ) {

		# someone passed a gene object!!!!
		my @lengths;
		foreach my $t ( get_transcripts($transcript) ) {
			push @lengths, get_transcript_length($t);
		}

		# return the longest transcript length
		return ( sort { $b <=> $a } @lengths )[0];
	}
	my $total = 0;
	foreach my $e ( get_exons($transcript) ) {
		$total += $e->length;
	}
	return $total;
}

######## CDS Methods

sub is_coding {
	my $transcript = shift;
	return                             unless $transcript;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	if ( $transcript->primary_tag =~ /gene$/i ) {

		# someone passed a gene, check its subfeatures
		my $code_potential = 0;
		foreach ( $transcript->get_SeqFeatures ) {
			$code_potential += is_coding($_);
		}
		return $code_potential;
	}
	return 1 if $transcript->primary_tag =~ /mrna/i;              # assumption
	return 1 if $transcript->source      =~ /protein.?coding/i;
	if ( $transcript->has_tag('transcript_biotype') ) {

		# ensembl type GTFs
		my ($biotype) = $transcript->get_tag_values('transcript_biotype');
		return $biotype =~ /protein.?coding/i ? 1 : 0;
	}
	elsif ( $transcript->has_tag('biotype') ) {

		# ensembl type GFFs
		my ($biotype) = $transcript->get_tag_values('biotype');
		return $biotype =~ /protein.?coding/i ? 1 : 0;
	}
	elsif ( $transcript->has_tag('gene_biotype') ) {

		# ensembl type GTFs
		# must be careful here, gene_biotype of course pertains to gene,
		# and not necessarily this particular transcript
		my ($biotype) = $transcript->get_tag_values('gene_biotype');
		return 1 if $biotype =~ /protein.?coding/i;
	}
	foreach ( $transcript->get_SeqFeatures ) {

		# old fashioned way
		return 1 if $_->primary_tag eq 'CDS';
	}
	return 0;
}

sub get_cds {
	my $transcript = shift;
	return                             unless $transcript;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	my @cds;
	foreach my $subfeat ( $transcript->get_SeqFeatures ) {
		push @cds, $subfeat if $subfeat->primary_tag eq 'CDS';
	}
	return unless @cds;
	@cds = map { $_->[0] }
		sort { $a->[1] <=> $b->[1] }
		map { [ $_, $_->start ] } @cds;
	return wantarray ? @cds : \@cds;
}

sub get_cdsStart {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	my $cds = get_cds($transcript);
	return unless $cds;
	if ( $transcript->strand >= 0 ) {
		return $cds->[0]->start;
	}
	else {
		# stop codons may or may not be not included
		my $codon = get_stop_codon($transcript);
		if ($codon) {
			return $codon->start < $cds->[0]->start ? $codon->start : $cds->[0]->start;
		}
		else {
			return $cds->[0]->start;
		}
	}
}

sub get_cdsEnd {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	my $cds = get_cds($transcript);
	return unless $cds;
	if ( $transcript->strand >= 0 ) {
		my $codon = get_stop_codon($transcript);
		if ($codon) {
			return $codon->end > $cds->[-1]->end ? $codon->end : $cds->[-1]->end;
		}
		else {
			return $cds->[-1]->end;
		}
	}
	else {
		return $cds->[-1]->end;
	}
}

sub get_transcript_cds_length {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	my $total = 0;
	foreach my $subf ( $transcript->get_SeqFeatures ) {
		next unless $subf->primary_tag eq 'CDS';
		$total += $subf->length;
	}
	return $total;
}

sub get_start_codon {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	my $start_codon;

	# look for existing one
	foreach my $subfeat ( $transcript->get_SeqFeatures ) {
		$start_codon = $subfeat if $subfeat->primary_tag =~ /start.?codon/i;
	}
	return $start_codon if $start_codon;

	# otherwise we have to build one
	my $cdss = get_cds($transcript);
	return unless $cdss;
	if ( $transcript->strand >= 0 ) {
		$start_codon = $transcript->new(
			-seq_id      => $transcript->seq_id,
			-source      => $transcript->source,
			-primary_tag => 'start_codon',
			-start       => $cdss->[0]->start,
			-end         => $cdss->[0]->start + 2,
			-strand      => 1,
			-phase       => 0,
			-primary_id  => $transcript->primary_id . '.start_codon',
		);
	}
	else {
		$start_codon = $transcript->new(
			-seq_id      => $transcript->seq_id,
			-source      => $transcript->source,
			-primary_tag => 'start_codon',
			-start       => $cdss->[-1]->end - 2,
			-end         => $cdss->[-1]->end,
			-strand      => -1,
			-phase       => 0,
			-primary_id  => $transcript->primary_id . '.start_codon',
		);
	}
	return $start_codon;
}

sub get_stop_codon {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	my $stop_codon;

	# look for existing one
	foreach my $subfeat ( $transcript->get_SeqFeatures ) {
		$stop_codon = $subfeat if $subfeat->primary_tag =~ /stop.?codon/i;
	}
	return $stop_codon if $stop_codon;

	# otherwise we have to build one
	# this entirely presumes that the stop codon is inclusive to the last cds
	# this is the case with GFF3 and UCSC tables, but not GTF
	my $cdss = get_cds($transcript);
	return unless $cdss;
	if ( $transcript->strand >= 0 ) {
		$stop_codon = $transcript->new(
			-seq_id      => $transcript->seq_id,
			-source      => $transcript->source,
			-primary_tag => 'stop_codon',
			-start       => $cdss->[-1]->end - 2,
			-end         => $cdss->[-1]->end,
			-strand      => 1,
			-phase       => 0,
			-primary_id  => $transcript->primary_id . '.stop_codon',
		);
	}
	else {
		$stop_codon = $transcript->new(
			-seq_id      => $transcript->seq_id,
			-source      => $transcript->source,
			-primary_tag => 'stop_codon',
			-start       => $cdss->[0]->start,
			-end         => $cdss->[0]->start + 2,
			-strand      => -1,
			-phase       => 0,
			-primary_id  => $transcript->primary_id . '.stop_codon',
		);
	}
	return $stop_codon;
}

sub get_utrs {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;

	# collect the various types of subfeatures
	my @exons;
	my @cdss;
	my @utrs;
	my @transcripts;
	foreach my $subfeat ( $transcript->get_SeqFeatures ) {
		my $type = $subfeat->primary_tag;
		if ( $type =~ /exon/i ) {
			push @exons, $subfeat;
		}
		elsif ( $type =~ /cds/i ) {
			push @cdss, $subfeat;
		}
		elsif ( $type =~ /utr|untranslated/i ) {
			push @utrs, $subfeat;
		}
		elsif ( $type =~ /rna|transcript/i ) {
			push @transcripts, $subfeat;
		}
	}

	# collect the utrs into final list
	my @list;
	if (@utrs) {

		# good, we don't have to do any work
		@list = @utrs;
	}
	elsif (@transcripts) {

		# someone must've passed us a gene
		foreach my $t (@transcripts) {
			my @u = get_utrs($t);
			push @list, @u;
		}
	}
	elsif ( @exons and @cdss ) {

		# calculate the utrs ourselves

		@exons = map { $_->[0] }
			sort { $a->[1] <=> $b->[1] }
			map { [ $_, $_->start ] } @exons;
		@cdss = map { $_->[0] }
			sort { $a->[1] <=> $b->[1] }
			map { [ $_, $_->start ] } @cdss;
		my $firstCDS = $cdss[0];
		my $lastCDS  = $cdss[-1];
		while (@exons) {
			my $exon = shift @exons;
			if ( $exon->end < $firstCDS->start ) {

				# whole exon is UTR
				my $utr = $exon->duplicate;
				$utr->primary_tag(
					$transcript->strand >= 0 ? 'five_prime_UTR' : 'three_prime_UTR' );
				$utr->display_name( $exon->display_name . '.utr' );
				push @list, $utr;
			}
			elsif ( $exon->overlaps($firstCDS) ) {

				# partial UTR on left side
				my $pieces = $exon->subtract($firstCDS);    # array ref of pieces
				next unless $pieces;
				my $utr = $pieces->[0];    # we will want the first one if there are two
				$utr->primary_tag(
					$transcript->strand >= 0 ? 'five_prime_UTR' : 'three_prime_UTR' );
				$utr->display_name( $exon->display_name . '.utr' );
				$utr->strand( $exon->strand );
				$utr->source( $exon->strand );
				push @list, $utr;
			}
			elsif ( $exon->start > $firstCDS->end and $exon->end < $lastCDS->start ) {

				# CDS exon
				next;
			}
			elsif ( $exon->overlaps($lastCDS) ) {

				# partial UTR
				my $pieces = $exon->subtract($lastCDS);    # array ref of pieces
				next unless $pieces;
				my $utr = $pieces->[-1];    # we will want the second one if there are two
				$utr->primary_tag(
					$transcript->strand >= 0 ? 'three_prime_UTR' : 'five_prime_UTR' );
				$utr->display_name( $exon->display_name . '.utr' );
				$utr->strand( $exon->strand );
				$utr->source( $exon->strand );
				push @list, $utr;
			}
			elsif ( $exon->start > $lastCDS->end ) {

				# whole exon is UTR
				my $utr = $exon->duplicate;
				$utr->primary_tag(
					$transcript->strand >= 0 ? 'three_prime_UTR' : 'five_prime_UTR' );
				$utr->display_name( $exon->display_name . '.utr' );
				push @list, $utr;
			}
			else {
				# geometric error
				croak " programmer geometric error!";
			}
		}
	}
	else {
		# nothing usable found to identify UTRs
		return;
	}

	# we have our list
	return wantarray ? @list : \@list;
}

sub get_transcript_utr_length {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;
	my $utrs  = get_utrs($transcript);
	my $total = 0;
	foreach my $utr (@$utrs) {
		$total += $utr->length;
	}
	return $total;
}

sub get_5p_utrs {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;

	# get all UTRs
	my $utrs = get_utrs($transcript);
	return unless scalar(@$utrs);

	my @fivers = grep { $_->primary_tag =~ /5|five/i } @$utrs;
	return wantarray ? @fivers : \@fivers;
}

sub get_3p_utrs {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref $transcript =~ /seqfeature/i;

	# get all UTRs
	my $utrs = get_utrs($transcript);
	return unless scalar(@$utrs);

	my @threes = grep { $_->primary_tag =~ /3|three/i } @$utrs;
	return wantarray ? @threes : \@threes;
}

sub get_transcript_5p_utr_length {
	my $transcript = shift;
	my $utrs       = get_5p_utrs($transcript);
	my $total      = 0;
	foreach my $utr (@$utrs) {
		$total += $utr->length;
	}
	return $total;
}

sub get_transcript_3p_utr_length {
	my $transcript = shift;
	my $utrs       = get_3p_utrs($transcript);
	my $total      = 0;
	foreach my $utr (@$utrs) {
		$total += $utr->length;
	}
	return $total;
}

#### Export methods

sub gff_string {

	# Bio::ToolBox::SeqFeature and Bio::SeqFeature::Lite objects have this method
	# otherwise this will die
	return shift->gff_string(@_);
}

sub gtf_string {
	my $feature = shift;
	my $gene    = shift || undef;    # only present when recursing
	confess "not a SeqFeature object!" unless ref $feature =~ /seqfeature/i;

	# process a gene
	if ( $feature->primary_tag =~ /gene$/i and not defined $gene ) {
		my $string;
		foreach my $t ( get_transcripts($feature) ) {
			$string .= gtf_string( $t, $feature );
		}
		return $string;
	}

	# check that we have transcribed feature with exons
	my @exons = get_exons($feature);
	return unless @exons;    # no exon subfeatures? must not be a transcript....

	# mandatory identifiers
	my ( $gene_id, $gene_name, $gene_biotype );
	if ($gene) {
		$gene_id   = $gene->primary_id || $gene->display_name;
		$gene_name = $gene->display_name;
		($gene_biotype) =
			   $gene->get_tag_values('gene_biotype')
			|| $gene->get_tag_values('biotype')
			|| undef;
	}
	else {
		# these attributes might still be present for transcripts
		($gene_id)      = $feature->get_tag_values('gene_id')      || undef;
		($gene_name)    = $feature->get_tag_values('gene_name')    || undef;
		($gene_biotype) = $feature->get_tag_values('gene_biotype') || undef;
	}
	my $trx_id   = $feature->primary_id || $feature->display_name;
	my $trx_name = $feature->display_name;
	my $group    = sprintf(
		"gene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\"; transcript_name \"%s\";",
		$gene_id, $trx_id, $gene_name, $trx_name );

	# add additional transcript attributes that might be interesting to keep
	if ($gene_biotype) {
		$group .= " gene_biotype \"$gene_biotype\";";
	}
	my ($biotype) = $feature->get_tag_values('transcript_biotype')
		|| $feature->get_tag_values('biotype');
	if ($biotype) {
		$group .= " transcript_biotype \"$biotype\";";
	}
	my ($tsl) = $feature->get_tag_values('transcript_support_level');
	if ($tsl) {
		$group .= " transcript_support_level \"$tsl\";";
	}

	# skip transcript as it is technically not part of the GTF standard....

	# convert exon subfeatures collected above
	my $string;
	my @cds = get_cds($feature);
	push @cds, get_stop_codon($feature);
	push @cds, get_start_codon($feature);
	@exons = map { $_->[0] }
		sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] }
		map { [ $_, $_->start, $_->end ] } ( @exons, @cds );
	foreach my $subf (@exons) {
		$string .= join(
			"\t",
			(
				$subf->seq_id     || '.',
				$subf->source_tag || '.',
				$subf->primary_tag,
				$subf->start,
				$subf->end,
				defined $subf->score ? $subf->score : '.',
				$subf->strand < 0    ? '-'          : '+',
				defined $subf->phase ? $subf->phase : '.',
				"$group\n"
			)
		);
	}

	# finished
	return $string;
}

sub ucsc_string {
	my $feature = shift;
	confess "not a SeqFeature object!" unless ref $feature =~ /seqfeature/i;
	my @ucsc_list;

	# process according to type
	if ( $feature->primary_tag =~ /gene$/i ) {

		# a gene object, we will need to process it's transcript subfeatures
		foreach my $transcript ( get_transcripts($feature) ) {
			my $ucsc = _process_ucsc_transcript( $transcript, $feature );
			push @ucsc_list, $ucsc if $ucsc;
		}
	}
	elsif ( $feature->primary_tag =~ /rna|transcript/i ) {

		# some sort of RNA transcript
		my $ucsc = _process_ucsc_transcript($feature);
		push @ucsc_list, $ucsc if $ucsc;
	}
	else {
		return;
	}

	# return strings
	my $string;
	foreach my $ucsc (@ucsc_list) {
		$string .= sprintf(
			"%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n",
			$ucsc->{'name2'},
			$ucsc->{'name'},
			$ucsc->{'chr'},
			$ucsc->{'strand'},
			$ucsc->{'txStart'},
			$ucsc->{'txEnd'},
			$ucsc->{'cdsStart'},
			$ucsc->{'cdsEnd'},
			$ucsc->{'exonCount'},
			join( ",", @{ $ucsc->{'exonStarts'} } ) . ',',
			join( ",", @{ $ucsc->{'exonEnds'} } ) . ',',
		);
	}
	return $string;
}

sub bed_string {
	my $feature = shift;
	confess "not a SeqFeature object!" unless ref $feature =~ /seqfeature/i;
	my @ucsc_list;

	# process according to type
	if ( $feature->primary_tag =~ /gene$/i ) {

		# a gene object, we will need to process it's transcript subfeatures
		foreach my $transcript ( get_transcripts($feature) ) {
			my $ucsc = _process_ucsc_transcript( $transcript, $feature );
			push @ucsc_list, $ucsc if $ucsc;
		}
	}
	elsif ( $feature->primary_tag =~ /rna|transcript/i ) {

		# some sort of RNA transcript
		my $ucsc = _process_ucsc_transcript($feature);
		push @ucsc_list, $ucsc if $ucsc;
	}
	else {
		return;
	}

	# return strings
	my $string;
	foreach my $ucsc (@ucsc_list) {

		# exon sizes
		my @sizes;
		for ( my $i = 0; $i < scalar( @{ $ucsc->{'exonStarts'} } ); $i++ ) {
			push @sizes, $ucsc->{'exonEnds'}->[$i] - $ucsc->{'exonStarts'}->[$i];
		}
		$string .= sprintf(
			"%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n",
			$ucsc->{'chr'},
			$ucsc->{'txStart'},
			$ucsc->{'txEnd'},
			$ucsc->{'name'},
			$ucsc->{'score'},
			$ucsc->{'strand'},
			$ucsc->{'cdsStart'},
			$ucsc->{'cdsEnd'},
			0,
			$ucsc->{'exonCount'},
			join( ",", @sizes ),
			join( ",", map { $_ - $ucsc->{'txStart'} } @{ $ucsc->{'exonStarts'} } ),
		);
	}
	return $string;
}

#### filter methods

sub filter_transcript_support_level {
	my $gene = shift;
	return unless $gene;
	my $min_tsl = shift || 'best';
	my @list    = qw(1 2 3 4 5 NA Missing);

	# get transcripts
	my @transcripts;
	if ( ref $gene =~ /seqfeature/i ) {
		@transcripts = get_transcripts($gene);
	}
	elsif ( ref $gene eq 'ARRAY' ) {
		@transcripts = @$gene;
	}
	else {
		return;
	}
	return unless @transcripts;

	# categorize transcripts
	my %results = map { $_ => [] } @list;
	foreach my $t (@transcripts) {
		my ($tsl) = $t->get_tag_values('transcript_support_level');
		$tsl ||= 'Missing';    # default in case nothing is present
		push @{ $results{$tsl} }, $t;
	}

	# take appropriate transcripts
	my @keepers;
	if ( $min_tsl eq 'best' ) {

		# progress through all levels and keep only the highest set
		foreach my $tsl (@list) {
			if ( scalar @{ $results{$tsl} } ) {
				@keepers = @{ $results{$tsl} };
				last;
			}
		}
	}
	elsif ( $min_tsl =~ /^best(\d)$/ ) {

		# progress through all levels and take all up to specified max
		my $max = $1;
		foreach my $tsl ( 1 .. 5 ) {
			if ( scalar @{ $results{$tsl} } ) {
				if ( $tsl <= $max ) {
					push @keepers, @{ $results{$tsl} };
				}
				elsif ( not @keepers and $tsl > $max ) {

					# go ahead and take it if we have no other option
					push @keepers, @{ $results{$tsl} };
				}
			}
		}

		# still take the NA and Missing if nothing else
		if ( scalar @keepers == 0 and scalar @{ $results{'NA'} } ) {
			@keepers = @{ $results{'NA'} };
		}
		elsif ( scalar @keepers == 0 and scalar @{ $results{'Missing'} } ) {
			@keepers = @{ $results{'Missing'} };
		}
	}
	elsif ( $min_tsl =~ /^\d$/ ) {

		# take only the specified level
		@keepers = @{ $results{$min_tsl} };
	}
	elsif ( $min_tsl eq 'NA' ) {

		# take only the NA guys
		@keepers = @{ $results{'NA'} };
	}
	else {
		confess "unrecognized minimum TSL value '$min_tsl' Check the documentation!";
	}
	@keepers = @{ $results{'Missing'} } unless @keepers;

	# return
	return _return_filtered_transcripts( $gene, \@keepers );
}

sub filter_transcript_gencode_basic {
	my $gene = shift;
	return unless $gene;

	# get transcripts
	my @transcripts;
	if ( ref $gene =~ /seqfeature/i ) {
		@transcripts = get_transcripts($gene);
	}
	elsif ( ref $gene eq 'ARRAY' ) {
		@transcripts = @$gene;
	}
	else {
		return;
	}
	return unless @transcripts;

	# take appropriate transcripts
	my @keepers;
	foreach my $t (@transcripts) {
		my ($basic) = $t->get_tag_values('tag');
		if ( $basic and $basic eq 'basic' ) {
			push @keepers, $t;
		}
	}

	# return
	return _return_filtered_transcripts( $gene, \@keepers );
}

sub filter_transcript_biotype {
	my $gene = shift;
	return unless $gene;
	my $check = shift;
	confess "no biotype value to check provided!" unless $check;

	# get transcripts
	my @transcripts;
	if ( ref $gene =~ /seqfeature/i ) {
		@transcripts = get_transcripts($gene);
	}
	elsif ( ref $gene eq 'ARRAY' ) {
		@transcripts = @$gene;
	}
	else {
		return;
	}
	return unless @transcripts;

	# take appropriate transcripts
	my @keepers;
	foreach my $t (@transcripts) {
		my ($value) = $t->get_tag_values('transcript_biotype')
			|| $t->get_tag_values('biotype');
		$value ||= $t->primary_tag;
		if ( $value and $value =~ /$check/i ) {
			push @keepers, $t;
		}
	}

	# return
	return _return_filtered_transcripts( $gene, \@keepers );
}

#### internal methods

sub _process_ucsc_transcript {
	my $transcript = shift;
	my $gene       = shift || undef;

	# initialize ucsc hash
	my $ucsc = {
		'name'       => $transcript->display_name || $transcript->primary_id,
		'name2'      => undef,
		'chr'        => $transcript->seq_id,
		'strand'     => $transcript->strand < 0 ? '-' : '+',
		'txStart'    => $transcript->start - 1,
		'txEnd'      => $transcript->end,
		'cdsStart'   => undef,
		'cdsEnd'     => undef,
		'exonCount'  => 0,
		'exonStarts' => [],
		'exonEnds'   => [],
		'score'      => $transcript->score || 1000,
	};

	# determine gene name
	if ($gene) {

		# use provided gene name
		$ucsc->{'name2'} = $gene->display_name || $gene->primary_id;
	}
	else {
		# reuse the name
		$ucsc->{'name2'} = $ucsc->{'name'};
	}

	# record CDS points
	if ( is_coding($transcript) ) {
		my $start = get_cdsStart($transcript);
		my $stop  = get_cdsEnd($transcript);
		if ( $start and $stop ) {
			$ucsc->{cdsStart} = $start - 1;
			$ucsc->{cdsEnd}   = $stop;
		}
		else {
			# if we don't have cds start/stop, then assign to transcript end
			# as if it is a noncoding transcript, regardless of primary_tag
			$ucsc->{cdsStart} = $transcript->end;
			$ucsc->{cdsEnd}   = $transcript->end;
		}
	}
	else {
		# non-coding transcript sets the cds points to the transcript end
		$ucsc->{cdsStart} = $transcript->end;
		$ucsc->{cdsEnd}   = $transcript->end;
	}

	# record the exons
	foreach my $exon ( get_exons($transcript) ) {
		push @{ $ucsc->{'exonStarts'} }, $exon->start - 1;
		push @{ $ucsc->{'exonEnds'} },   $exon->end;
		$ucsc->{'exonCount'} += 1;
	}

	return $ucsc;
}

sub _return_filtered_transcripts {
	my ( $gene, $keepers ) = @_;

	if ( ref $gene =~ /seqfeature/i ) {

		# first check if we were only given a transcript
		if ( $gene->primary_tag =~ /transcript|rna/i ) {

			# we must have been given a single transcript to check, so return it
			$keepers->[0] ||= undef;
			return $keepers->[0];
		}

		# we can't delete subfeatures, so we're forced to create a new
		# parent gene and reattach the filtered transcripts
		my %attributes = $gene->attributes;
		return $gene->new(
			-seq_id      => $gene->seq_id,
			-start       => $gene->start,
			-end         => $gene->end,
			-strand      => $gene->strand,
			-primary_tag => $gene->primary_tag,
			-source      => $gene->source_tag,
			-name        => $gene->display_name,
			-id          => $gene->primary_id,
			-attributes  => \%attributes,
			-segments    => $keepers,
		);
	}
	else {
		return $keepers;
	}
}

1;

__END__

=head1 NAME

Bio::ToolBox::GeneTools - SeqFeature agnostic methods for working with gene models

=head1 SYNOPSIS

    use Bio::ToolBox::GeneTools qw(:all);
    
    my $gene; # a SeqFeatureI compliant gene object obtained elsewhere
              # for example, from Bio::DB::SeqFeature::Store database
              # or parsed from a GFF3, GTF, or UCSC-style gene table using 
              # Bio::ToolBox parsers
    
    if (is_coding($gene)) { # boolean test
    	
    	# collect all exons from all transcripts in gene
    	my @exons = get_exons($gene);
    	
    	# find just the alternate exons used only once
    	my @alternate_exons = get_alt_exons($gene);
    	
    	# collect UTRs, which may not be defined in the original source
    	my @utrs;
    	foreach my $t (get_transcripts($gene)) {
    		my @u = get_utrs($t);
    		push @utrs, @u;
    	}
    }


=head1 DESCRIPTION

This module provides numerous exportable functions for working with gene 
SeqFeature models. This assumes that the gene models follow the BioPerl 
L<Bio::SeqFeatureI> convention with nested SeqFeature objects representing the 
gene, transcript, and exons. For example, 

    gene
      transcript
        exon
        CDS

Depending upon how the SeqFeatures were generated or defined, subfeatures 
may or may not be defined or be obvious. For example, UTRs or introns may 
not be present. Furthermore, the C<primary_tag> or type may not follow 
Sequence Ontology terms. Regular expressions are deployed to handle 
varying naming schemes and exceptions.

These functions should work with most or all L<Bio::SeqFeatureI> compliant 
objects. It has been tested with L<Bio::ToolBox::SeqFeature>, 
L<Bio::SeqFeature::Lite>, and L<Bio::DB::SeqFeature> classes.

New SeqFeature objects that are generated use the same class for 
simplicity and expediency. 

=head1 METHODS

=head2 Function Import

None of the functions are exported by default. Specify which ones you want 
when you import the module. Alternatively, use one of the tags below.

=over 4

=item :all

Import all of the methods.

=item :exon

Import all of the exon methods, including L</get_exons>, L</get_alt_exons>, 
L</get_common_exons>, L</get_uncommon_exons>, and L</get_alt_common_exons>.

=item :intron

Import all of the intron methods, including L</get_introns>, L</get_alt_introns>, 
L</get_common_introns>, L</get_uncommon_introns>, and L</get_alt_common_introns>.

=item :transcript

Import the transcript related methods, including L</get_transcripts>, 
L</get_transcript_length>, and L</collapse_transcripts>.

=item :cds

Import the CDS pertaining methods, including L</is_coding>, L</get_cds>, 
L</get_cdsStart>, L</get_cdsEnd>, L</get_transcript_cds_length>, and L</get_utrs>.

=item :export

Import all of the export methods, including L</gff_string>, L</gtf_string>, 
L</ucsc_string>, and L</bed_string>;

=item :filter

Import all of the transcript filter methods, including L</filter_transcript_biotype>,
L</filter_transcript_gencode_basic>, and L</filter_transcript_support_level>.

=back

=head2 Exon Methods

Functions to get a list of exons from a gene or transcript

=over 4

=item get_exons

	my @exons = get_exons($gene);
	my @exons = get_exons($transcript);

This will return an array or array reference of all the exon subfeatures in 
the SeqFeature object, either gene or transcript. No discrimination whether 
they are used once or more than once. Non-defined exons can be assembled from 
CDS and/or UTR subfeatures. Exons are sorted by start coordinate.

=item get_alt_exons

	my @alternate_exons = get_alt_exons($gene);

This will return an array or array reference of all the exon subfeatures for 
a multi-transcript gene that are used only once in all of the transcripts.

=item get_common_exons

	my @common_exons = get_common_exons($gene);

This will return an array or array reference of all the exon subfeatures for 
a multi-transcript gene that are used in all of the transcripts.

=item get_uncommon_exons

	my @uncommon_exons = get_uncommon_exons($gene);

This will return an array or array reference of all the exon subfeatures for 
a multi-transcript gene that are used in some of the transcripts, i.e. more 
than one but not all.

=item get_alt_common_exons

	my %exon_hash = get_alt_common_exons($gene);

This will return a hash reference with several keys, including "common", 
"uncommon", and each of the transcript IDs. Each key value is an array 
reference with the exons for that category. The "common" will be all 
common exons, "uncommon" will be uncommon exons (used more than once but 
less than all), and each transcript ID will include their specific alternate 
exons (used only once).

For genes with only a single transcript, all exons will be marked as "common" 
for simplicity, although technically they could all be considered "alternate" 
since they're only used once.

=back

=head2 Intron Methods

Functions to get a list of introns from a gene or transcript. Introns are 
not usually defined in gene annotation files, but are inferred from the 
exons and total gene or transcript length. In this case, new SeqFeature 
elements are generated for each intron.

=over 4

=item get_introns

	my @introns = get_introns($gene);
	my @introns = get_introns($transcript);

This will return an array or array reference of all the intron subfeatures in 
the SeqFeature object, either gene or transcript. No discrimination whether 
they are used once or more than once. Non-defined introns can be assembled from 
CDS andE<sol>or UTR subfeatures. Introns are sorted by start coordinate.

=item get_alt_introns

	my @alternate_introns = get_alt_introns($gene);

This will return an array or array reference of all the intron subfeatures for 
a multi-transcript gene that are used only once in all of the transcripts.

=item get_common_introns

	my @common_introns = get_common_introns($gene);

This will return an array or array reference of all the intron subfeatures for 
a multi-transcript gene that are used in all of the transcripts.

=item get_uncommon_introns

	my @uncommon_introns = get_uncommon_introns($gene);

This will return an array or array reference of all the intron subfeatures for 
a multi-transcript gene that are used in some of the transcripts, i.e. more 
than one but not all.

=item get_alt_common_introns

	my %intron_hash = get_alt_common_introns($gene);

This will return a hash reference with several keys, including "common", 
"uncommon", and each of the transcript IDs. Each key value is an array 
reference with the introns for that category. The "common" will be all 
common introns, "uncommon" will be uncommon introns (used more than once but 
less than all), and each transcript ID will include their specific alternate 
introns (used only once).

For genes with only a single transcript, all introns will be marked as "common" 
for simplicity, although technically they could all be considered "alternate" 
since they're only used once.

=back

=head2 Transcript Methods

These methods work on transcripts, typically alternate transcripts from a 
gene SeqFeature.

=over 4

=item get_transcripts

	my @transcripts = get_transcripts($gene);

Returns an array or array reference of the transcripts associated with a 
gene feature.

=item collapse_transcripts

	my $new_transcript = collapse_transcripts($gene);
	my $new_transcript = collapse_transcripts($transcript1, $transcript2);

This method will collapse all of the transcripts associated with a gene 
SeqFeature into a single artificial transcript, merging exons as necessary 
to maximize exon length and minimize introns. This is useful when 
performing, for example, RNASeq analysis on genes. A single SeqFeature 
transcript object is returned containing the merged exon subfeatures. 

Pass either a gene or a list of transcripts to collapse.

=item get_transcript_length

	my $length = get_transcript_length($transcript);

Calculates and returns the transcribed length of a transcript, i.e 
the sum of its exon lengths. B<Warning!> If you pass a gene object, you 
will get the maximum of all transcript exon lengths, which may not be 
what you anticipate!

=back

=head2 CDS methods

These methods calculate values related to the coding sequence of the 
mRNA transcripts. 

=over 4

=item is_coding

	if( is_coding($transcript) ) {
		# do something
	}

This method will return a boolean value (1 or 0) if the passed transcript object 
appears to be a coding transcript. GFF and GTF files are not always immediately 
clear about the type of transcript; there are (unfortunately) multiple ways 
to encode the feature as a protein coding transcript: C<primary_tag>, 
C<source_tag>, GFF attribute, presence of CDS subfeatures, etc. 
This method checks all of these possibilities. B<Note>: If you pass a 
multi-transcript gene, only one transcript need to be coding to pass a true 
value.

=item get_cds

	my @cds = get_cds($transcript);

Returns the CDS subfeatures of the given transcript, if they are 
defined. Returns either an array or array reference.

=item get_cdsStart

	my $start = get_cdsStart($trancript);

Returns the start coordinate of the CDS for the given transcript.
Note that this is the leftmost (smallest) coordinate of the CDS 
and not necessarily the coordinate of the start codon, similar to 
what the UCSC gene tables report. Use the transcript strand to 
determine the 5' end.

=item get_cdsEnd

	my $end = get_cdsEnd($trancript);

Returns the stop coordinate of the CDS for the given transcript.
B<Note> that this is the rightmost (largest) coordinate of the CDS 
and not necessarily the coordinate of the stop codon, similar to 
what the UCSC gene tables report. Use the transcript strand to 
determine which is the C<3'> and C<5'> end.

=item get_start_codon

	my $start_codon = get_start_codon($trancript);

Returns a SeqFeature object representing the start codon. If one is 
not explicitly defined in the hierarchy, then a new object is generated.

=item get_stop_codon

	my $stop_codon = get_stop_codon($transcript);

Returns a SeqFeature object representing the stop codon. If one is 
not defined in the hierarchy, then a new object is created. B<Note> that 
this assumes that the stop codon is inclusive to the defined CDS, which is 
the case with GFF3 and UCSC gene table derived features. On the other hand, 
features derived from GTF is defined with the stop codon exclusive to the CDS. 
This shouldn't matter with GTF, however, since GTF usually explicitly includes 
stop codon features, whereas the other two formats do not.

=item get_transcript_cds_length

	my $length = get_transcript_cds_length($transcript);

Calculates and returns the length of the coding sequence for a 
transcript, i.e. the sum of the CDS lengths.

=item get_utrs

	my @utrs = get_utrs($trancript);

Returns both C<5'> and C<3'> untranslated regions of the transcript. If these are 
not defined in the SeqFeature subfeature hierarchy, then the coordinates will be 
determined from from the exon and CDS subfeatures, if available, and new SeqFeature 
objects generated. Non-coding transcripts will not return anything. 

=item get_5p_utrs

	my @5p_utrs = get_5p_utrs($trancript);

Returns only the C<5'> untranslated regions of the transcript.

=item get_3p_utrs($transcript)

	my @3p_utrs = get_3p_utrs($trancript);

Returns only the C<3'> untranslated regions of the transcript.

=back

=head2 Export methods

These methods are used for exporting a gene andE<sol>or transcript model into 
a text string based on the specified format. 

=over 4

=item gff_string

	my $string .= gff_string($gene, 1);
	my $string .= gff_string($transcript, 1);

This is just a convenience method. SeqFeature objects based on 
L<Bio::SeqFeature::Lite>, L<Bio::DB::SeqFeature>, or L<Bio::ToolBox::SeqFeature>
have a C<gff_string> method, and this will simply call that method. SeqFeature 
objects that do not have this method will, of course, cause the script to 
terminate. 

Pass the seqfeature object and a boolean value to recursively append all 
subfeatures (e.g. exons) to the string. In most cases, this will generate a 
GFF3-style string.

L<Bio::ToolBox::Data::Feature> also provides a simplified gff_string method.

=item gtf_string

	my $string .= gtf_string($gene);
	my $string .= gtf_string($transcript);

This will export a gene or transcript model as a series of GTF formatted 
text lines, following the defined Gene Transfer Format (also known as GFF 
version 2.5). It will ensure that each feature is properly tagged with the 
C<gene_id> and C<transcript_id> attributes. 

This method will automatically recurse through all subfeatures.

=item ucsc_string

	my $string = ucsc_string($gene);

This will export a gene or transcript model as a refFlat formatted Gene 
Prediction line (11 columns). See L<http://genome.ucsc.edu/FAQ/FAQformat.html#format9>
for details. Multiple transcript genes are exported as multiple text lines 
concatenated together.

=item bed_string

	my $string = bed_string($gene);

This will export a gene or transcript model as a UCSC Bed formatted transcript 
line (12 columns). See L<http://genome.ucsc.edu/FAQ/FAQformat.html#format1>
for details. Multiple transcript genes are exported as multiple text lines 
concatenated together. Note that gene information is not preserved with Bed12 
files; only the transcript name is used. The C<RGB> value is set to 0.

=back

=head2 Filter methods

These methods are used to filter genes.

=over 4

=item filter_transcript_support_level

	my $new_gene = filter_transcript_support_level($gene, 'best2');
	my @good_transcripts = filter_transcript_support_level(\@transcripts);

This will filter a gene object for transcripts that match or exceed the 
provided transcript support level. This assumes that the transcripts 
contain the attribute tag 'transcript_support_level', which are present in 
Ensembl provided GFF3 and GTF annotation files. The values are a digit (1-5), 
or 'NA', where 1 is experimentally supported and 5 is entirely predicted 
with no experimental evidence. See 
L<Ensembl TSL glossary entry|http://www.ensembl.org/info/website/glossary.html> 
for details. 

Pass a gene SeqFeature object with one or more transcript subfeatures. 
Alternatively, an array reference of transcripts could be passed as well.

A level may be provided as a second argument. The default is 'best'.

=over 4

=item best

Only the transcripts with the highest existing value will be retained.

=item bestE<lt>digitE<gt>

All transcripts up to the indicated level are retained. For example, 
'best3' would indicate that transcripts with support levels 1, 2, and 3 
would be retained. 

=item E<lt>digitE<gt>

Only transcripts at the given level are retained.

=item NA

Only transcripts with 'NA' as the value are retained. These are typically 
pseudogenes or single-exon transcripts.

=back

If none of the transcripts have the attribute, then all are returned 
(nothing is filtered). 

If a gene object was provided, a new gene object will be returned with 
only the retained transcripts as subfeatures. If an array reference of 
transcripts was provided, then an array reference of the filtered 
transcripts is returned.

=item filter_transcript_gencode_basic

	my $new_gene = filter_transcript_gencode_basic($gene);
	my @good_transcripts = filter_transcript_gencode_basic(\@transcripts);

This will filter a gene object for transcripts for the Ensembl GENCODE 
tag "basic", which indicates that a transcript is tagged as GENCODE Basic 
transcript. 

If a gene object was provided, a new gene object will be returned with 
only the retained transcripts as subfeatures. If an array reference of 
transcripts was provided, then an array reference of the filtered 
transcripts is returned.
 
=item filter_transcript_biotype

	my $new_gene = filter_transcript_gencode_basic($gene, $biotype);
	my @good_transcripts = filter_transcript_gencode_basic(\@transcripts, 'miRNA');

This will filter a gene object for transcripts for specific biotype values 
using the C<transcript_biotype> or C<biotype> attribute tags, commonly found 
in Ensembl annotation.

If a gene object was provided, a new gene object will be returned with 
only the retained transcripts as subfeatures. If an array reference of 
transcripts was provided, then an array reference of the filtered 
transcripts is returned.

=back

=head1 SEE ALSO

L<Bio::ToolBox::SeqFeature>, L<Bio::ToolBox::parser::ucsc>, 
L<Bio::ToolBox::parser::gff>, L<Bio::ToolBox::parser::bed>, L<Bio::Tools::GFF>,
L<Bio::SeqFeature::Lite>, L<Bio::DB::SeqFeature>, L<Bio::SeqFeatureI>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
