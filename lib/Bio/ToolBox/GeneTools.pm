package Bio::ToolBox::GeneTools;
our $VERSION = '1.44';

=head1 NAME

Bio::ToolBox::GeneTools - SeqFeature agnostic methods for working with gene models

=head1 SYNOPSIS

    use Bio::ToolBox::GeneTools qw(:all);
    
    my $gene; # a SeqFeatureI compliant gene object obtained elsewhere
              # for example, from Bio::DB::SeqFeature::Store database
              # or parsed from a GFF3, GTF, or UCSC-style gene table using 
              # Bio::ToolBox::parser::(gff,ucsc) parsers
    
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
not be present. Furthermore, the primary_tag or type may not follow 
Sequence Ontology terms. Regular expressions are deployed to handle 
varying naming schemes and exceptions.

These functions should work with most or all <Bio::SeqFeatureI> compliant 
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

Import all of the exon methods, including get_exons(), get_alt_exons(), 
get_common_exons(), get_uncommon_exons(), and get_alt_common_exons().

=item :intron

Import all of the intron methods, including get_introns(), get_alt_introns(), 
get_common_introns(), get_uncommon_introns(), and get_alt_common_introns().

=item :transcript

Import the transcript related methods, including get_transcripts(), 
get_transcript_length(), and collapse_transcripts().

=item :cds

Import the CDS pertaining methods, including is_coding(), get_cds(), 
get_cdsStart(), get_cdsEnd(), get_transcript_cds_length(), and get_utrs().

=back

=head2 Exon Methods

Functions to get a list of exons from a gene or transcript

=over 4

=item get_exons($gene)

=item get_exons($transcript)

This will return an array or array reference of all the exon subfeatures in 
the SeqFeature object, either gene or transcript. No discrimination whether 
they are used once or more than once. Non-defined exons can be assembled from 
CDS and/or UTR subfeatures. Exons are sorted by start coordinate.

=item get_alt_exons($gene)

This will return an array or array reference of all the exon subfeatures for 
a multi-transcript gene that are used only once in all of the transcripts.

=item get_common_exons($gene)

This will return an array or array reference of all the exon subfeatures for 
a multi-transcript gene that are used in all of the transcripts.

=item get_uncommon_exons($gene)

This will return an array or array reference of all the exon subfeatures for 
a multi-transcript gene that are used in some of the transcripts, i.e. more 
than one but not all.

=item get_alt_common_exons($gene)

This will return a hash reference with several keys, including "common", 
"uncommon", and each of the transcript IDs. Each key value is an array 
reference with the exons for that category. The "common" will be all 
common exons, "uncommon" will be uncommon exons (used more than once but 
less than all), and each transcript ID will include their specific alternate 
exons (used only once).

=back

=head2 Intron Methods

Functions to get a list of introns from a gene or transcript. Introns are 
not usually defined in gene annotation files, but are inferred from the 
exons and total gene or transcript length. In this case, new SeqFeature 
elements are generated for each intron.

=over 4

=item get_introns($gene)

=item get_introns($transcript)

This will return an array or array reference of all the intron subfeatures in 
the SeqFeature object, either gene or transcript. No discrimination whether 
they are used once or more than once. Non-defined introns can be assembled from 
CDS and/or UTR subfeatures. Introns are sorted by start coordinate.

=item get_alt_introns($gene)

This will return an array or array reference of all the intron subfeatures for 
a multi-transcript gene that are used only once in all of the transcripts.

=item get_common_introns($gene)

This will return an array or array reference of all the intron subfeatures for 
a multi-transcript gene that are used in all of the transcripts.

=item get_uncommon_introns($gene)

This will return an array or array reference of all the intron subfeatures for 
a multi-transcript gene that are used in some of the transcripts, i.e. more 
than one but not all.

=item get_alt_common_introns($gene)

This will return a hash reference with several keys, including "common", 
"uncommon", and each of the transcript IDs. Each key value is an array 
reference with the introns for that category. The "common" will be all 
common introns, "uncommon" will be uncommon introns (used more than once but 
less than all), and each transcript ID will include their specific alternate 
introns (used only once).

=back

=head2 Transcript Methods

These methods work on transcripts, typically alternate transcripts from a 
gene SeqFeature.

=over 4

=item get_transcripts($gene)

Returns an array or array reference of the transcripts associated with a 
gene feature.

=item collapse_transcripts($gene)

=item collapse_transcripts($transcript1, $transcript2, ...)

This method will collapse all of the transcripts associated with a gene 
SeqFeature into a single artificial transcript, merging exons as necessary 
to maximize exon length and minimize introns. This is useful when 
performing, for example, RNASeq analysis on genes. A single SeqFeature 
transcript object is returned containing the merged exon subfeatures. 

=item get_transcript_length($transcript)

Calculates and returns the transcribed length of a transcript, i.e 
the sum of its exon lengths.

=back

=head2 CDS methods

These methods calculate values related to the coding sequence of the 
mRNA transcripts. 

=over 4

=item is_coding($gene)

=item is_coding($transcript)

This method will return a boolean value if the passed transcript object 
appears to be a coding transcript. GFF and GTF files are not always immediately 
clear about the type of transcript; there are (unfortunately) multiple ways 
to encode the feature as a protein coding transcript: primary_tag, source_tag, 
attribute, CDS subfeatures, etc. This method tries to determine this.

=item get_cds($transcript)

Returns the CDS subfeatures of the given transcript, if they are 
defined. Returns either an array or array reference.

=item get_cdsStart($transcript)

Returns the start coordinate of the CDS for the given transcript.
Note that this is the leftmost (smallest) coordinate of the CDS 
and not necessarily the coordinate of the start codon, similar to 
what the UCSC gene tables report. Use the transcript strand to 
determine the 5' end.

=item get_cdsEnd($transcript)

Returns the stop coordinate of the CDS for the given transcript.
Note that this is the rightmost (largest) coordinate of the CDS 
and not necessarily the coordinate of the stop codon, similar to 
what the UCSC gene tables report. Use the transcript strand to 
determine the 3' end.

=item get_start_codon($transcript)

Returns a SeqFeature object representing the start codon. If one is 
not defined in the hierarchy, then a new object is created.

=item get_stop_codon($transcript)

Returns a SeqFeature object representing the stop codon. If one is 
not defined in the hierarchy, then a new object is created. Not that 
this assumes that the stop codon is inclusive to the defined CDS.

=item get_transcript_cds_length($transcript)

Calculates and returns the length of the coding sequence for a 
transcript, i.e. the sum of the CDS lengths.

=item get_utrs($transcript)

Returns the 5' and 3' untranslated regions of the transcript. If these are 
not defined in the SeqFeature subfeature hierarchy, then they will be calculated 
from the exon and CDS subfeatures, if available. Non-coding transcripts will not 
return anything. 

=back

=head2 Export methods

These methods are used for exporting a gene and/or transcript model into 
a text string based on the specified format. 

=over 4

=item gff_string($gene)

This is just a convenience method. SeqFeature objects based on 
L<Bio::SeqFeature::Lite>, L<Bio::DB::SeqFeature>, or L<Bio::ToolBox::SeqFeature>
have a gff_string() method, and this will simply call that method. SeqFeature 
objects that do not have this method will, of course, cause the script to 
terminate. 

L<Bio::ToolBox::Feature> also provide a gff_string method.

=item ucsc_string($gene)

This will export a gene or transcript model as a refFlat formatted gene 
Prediction line (11 columns). See L<http://genome.ucsc.edu/FAQ/FAQformat.html#format9>
for details. Multiple transcript genes are exported as multiple text lines 
concatenated together.

=back

=cut

use strict;
use Carp qw(carp cluck croak confess);
require Exporter;

### Export
our @ISA = qw(Exporter);
our @EXPORT = qw();
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
	gff_string
	ucsc_string
);
our %EXPORT_TAGS = (
	all => \@EXPORT_OK,
	exon => [ qw(
		get_exons
		get_alt_exons
		get_common_exons
		get_uncommon_exons
		get_alt_common_exons
	) ],
	intron => [ qw(
		get_introns
		get_alt_introns
		get_common_introns
		get_uncommon_introns
		get_alt_common_introns
	) ],
	transcript => [ qw(
		get_transcripts
		collapse_transcripts
		get_transcript_length
	) ],
	cds => [ qw(
		is_coding
		get_cds
		get_cdsStart
		get_cdsEnd
		get_start_codon
		get_stop_codon
		get_transcript_cds_length
		get_utrs
	) ],
	export => [ qw(
		gff_string
		ucsc_string
	) ],
);


### The True Statement
1; 



######## Exon Methods

sub get_exons {
	# initialize
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref($transcript) =~ /seqfeature/i;
	my @exons;
	my @cdss;
	my @transcripts;
	
	# go through the subfeatures
	foreach my $subfeat ($transcript->get_SeqFeatures) {
		my $type = $subfeat->primary_tag;
		if ($type =~ /exon/i) {
			push @exons, $subfeat;
		}
		elsif ($type =~ /cds|utr|untranslated/i) {
			push @cdss, $subfeat;
		}
		elsif ($type =~ /rna|transcript/i) {
			push @transcripts, $subfeat;
		}
	}
	
	# check which array we'll use
	# prefer to use actual exon subfeatures, but those may not be defined
	my @list;
	if (@exons) {
		@list = map { $_->[0] } 
				sort { $a->[1] <=> $b->[1] }
				map { [$_, $_->start] } 
				@exons;
	}
	elsif (@cdss) {
		# duplicate the CDSs as exons
		@list = map { _duplicate($_) } @cdss;
		foreach (@list) {$_->primary_tag('exon')} # reset tag
		
		# make sure to merge adjacent exons
		@list = map { $_->[0] } # must sort first
				sort { $a->[1] <=> $b->[1] }
				map { [$_, $_->start] } 
				@list;
		for (my $i = 0; $i < scalar @list; $i++) {
			if (defined ($list[ $i + 1 ]) ) {
				if ($list[$i+1]->start - $list[$i]->end <= 1) {
					# need to merge
					$list[$i]->end( $list[$i+1]->end );
					splice(@list, $i + 1, 1); # remove the merged
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
				map { [$_, $_->start, $_->end] } 
				@list;
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
	foreach my $k (keys %$ac_exons) {
		next if $k eq 'common';
		next if $k eq 'uncommon';
		push @alts, @{ $ac_exons->{$k} };
	}
	if (@alts) {
		# re-sort in genomic order
		@alts = map {$_->[1]}
				sort {$a->[0] <=> $b->[0]}
				map { [$_->start, $_] } @alts;
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
	return _get_alt_common_things(1, @_);
}


######## Intron Methods

sub get_introns {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref($transcript) =~ /seqfeature/i;
	my @introns;
	
	# find the exons and/or CDSs
	my @exons = get_exons($transcript);
	return unless @exons;
	return if (scalar(@exons) == 1);
	
	# identify the last exon index position
	my $last = scalar(@exons) - 1;
	
	# forward strand
	if ($transcript->strand >= 0) {
		# each intron is created based on the previous exon
		for (my $i = 0; $i < $last; $i++) {
			my $e = $exons[$i];
			my $i = $e->new(
				-seq_id       => $e->seq_id,
				-start        => $e->end + 1,
				-end          => $exons[$i + 1]->start - 1, # up to start of next exon
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
		for (my $i = $last; $i > 0; $i--) {
			my $e = $exons[$i];
			my $i = $e->new(
				-seq_id       => $e->seq_id,
				-start        => $exons[$i - 1]->end + 1, # end of next exon
				-end          => $e->start - 1,
				-strand       => $transcript->strand,
				-primary_tag  => 'intron',
				-source_tag   => $transcript->source_tag,
				-primary_id   => $transcript->display_name . ".intron$i",
				-display_name => $transcript->display_name . ".intron$i",
			);
			push @introns, $i;
		}
		
		# reorder the introns based on start position
		if (@introns) {
			@introns = map { $_->[0] }
					sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] }
					map { [$_, $_->start, $_->end] } 
					@introns;
		}
	}
	
	# finished
	return wantarray ? @introns : \@introns;
}

sub get_alt_introns {
	my $ac_introns = get_alt_common_introns(@_);
	my @alts;
	foreach my $k (keys %$ac_introns) {
		next if $k eq 'common';
		next if $k eq 'uncommon';
		push @alts, @{ $ac_introns->{$k} };
	}
	if (@alts) {
		# re-sort in genomic order
		@alts = map {$_->[1]}
				sort {$a->[0] <=> $b->[0]}
				map { [$_->start, $_] } @alts;
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
	return _get_alt_common_things(0, @_);
}

sub _get_alt_common_things {
	# internal subroutine to get either exons or introns
	my $do_exon = shift; # true for exon, false for intron
	my @transcripts;
	return unless @_;
	if (scalar @_ == 1) {
		# someone passed a gene, get the transcripts
		confess "not a SeqFeature object!" unless ref($_[0]) =~ /seqfeature/i;
		@transcripts = get_transcripts($_[0]);
	}
	elsif (scalar @_ > 1) {
		# presume these are transcripts?
		@transcripts = @_;
	}
	
	# hash of transcript to things
	my %tx2things = (
		common => [],
		uncommon => []
	);
	
	# no transcripts provided?
	return \%tx2things unless @transcripts;
	
	# only one transcript provided?
	if (scalar @transcripts == 1) {
		# all exons are alternate by definition
		my $name = $transcripts[0]->display_name;
		my @things = $do_exon ? get_exons($transcripts[0]) : get_introns($transcripts[0]);
		$tx2things{$name} = \@things;
		return \%tx2things;
	}
	
	# get things and put them in has based on coordinates
	my %pos2things;
	foreach my $t (@transcripts) {
		my @things = $do_exon ? get_exons($t) : get_introns($t);
		foreach my $e (@things) {
			$pos2things{$e->start}{$e->end}{$t->display_name} = _duplicate($e);
		}
		$tx2things{ $t->display_name } = [];
	}
	
	# put things into categories based on commonality
	# associate things with unique transcripts, common, or uncommon sets
	my $trx_number = scalar @transcripts;
	foreach my $s (sort {$a <=> $b} keys %pos2things) {               # sort on start
		foreach my $e (sort {$a <=> $b} keys %{ $pos2things{$s} }) {  # sort on stop
			my @names = keys %{ $pos2things{$s}{$e} };
			if (scalar @names == 1) {
				# only 1 thing, must be an alternate
				push @{ $tx2things{$names[0]} }, $pos2things{$s}{$e}{$names[0]};
			}
			elsif (scalar @names == $trx_number) {
				# common to all transcripts, take the first one as example
				push @{ $tx2things{common} }, $pos2things{$s}{$e}{$names[0]};
			}
			else {
				# common to some but not all transcripts, so uncommon
				push @{ $tx2things{uncommon} }, $pos2things{$s}{$e}{$names[0]};
			}
		}
	}
	return \%tx2things;
}



######## Transcript Methods

sub get_transcripts {
	my $gene = shift;
	confess "not a SeqFeature object!" unless ref($gene) =~ /seqfeature/i;
	return $gene if ( $gene->primary_tag !~ /gene/i and 
		$gene->primary_tag =~ /rna|transcript/i);
	my @transcripts;
	my @exons;
	foreach my $subf ($gene->get_SeqFeatures) {
		if ($subf->primary_tag =~ /rna|transcript/i) {
			push @transcripts, $subf;
		}
		elsif ($subf->primary_tag =~ /^(?:cds|exon)$/i) {
			push @exons, $subf;
		}
	}
	if (not @transcripts and @exons) {
		# some weirdly formatted annotation files skip the transcript
		# looking at you SGD
		my $transcript = $gene->new(
			-seq_id         => $gene->seq_id,
			-start          => $gene->start,
			-end            => $gene->end,
			-strand         => $gene->strand,
			-primary_tag    => 'transcript',
			-source         => $gene->source_tag,
			-name           => $gene->display_name,
			-segments       => \@exons,
		);
		push @transcripts, $transcript;
	}
	@transcripts = map { $_->[0] }
		sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] }
		map { [$_, $_->start, $_->length] } 
		@transcripts;
	return wantarray ? @transcripts : \@transcripts;
}


sub collapse_transcripts {
	my @transcripts;
	return unless @_;
	my $example = $_[0]; # parent transcript seqfeature to model new transcript on
	if (scalar @_ == 1) {
		# someone passed a gene, get the transcripts
		@transcripts = get_transcripts($_[0]);
		return if scalar @transcripts == 1;
	}
	elsif (scalar @_ > 1) {
		@transcripts = @_;
	}

	# get all useable exons to collapse
	my @exons;
	foreach my $t (@transcripts) {
		my @e = get_exons($t);
		push @exons, @e;
	}
	
	# sort all the exons
	my @sorted = 	map { $_->[0] }
					sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] }
					map { [$_, $_->start, $_->end] } 
					@exons;
	
	# build new exons from the original - don't keep cruft
	my $next = shift @sorted;
	my @new;
	$new[0] = $next->new(
		-seq_id         => $next->seq_id,
		-start 			=> $next->start,
		-end   			=> $next->end,
		-primary_tag  	=> 'exon',
	);
	
	# work through remaining exons, adding and merging as necessary
	while (@sorted) {
		$next = shift @sorted;
		my ($ns, $ne) = ($next->start, $next->end); # new start & end
		my ($os, $oe) = ($new[-1]->start, $new[-1]->end); # old start & end
		if ($ns == $os and $ne > $oe) {
			# same beginning, further end
			$new[-1]->end($ne);
		}
		elsif ($ns > $os and $ns < $oe and $ne > $oe) {
			# overlapping start, further end
			$new[-1]->end($ne);
		}
		elsif ($ns > $oe) {
			# completely new exon
			push @new, $next->new(
				-seq_id         => $next->seq_id,
				-start 			=> $ns,
				-end   			=> $ne,
				-primary_tag  	=> 'exon',
			);
		}
		# all other possibilities we can skip
	}
	
	# return the assembled transcript
	return $example->new(
		-seq_id         => $example->seq_id,
		-start          => $new[0]->start,
		-end            => $new[-1]->end,
		-strand         => $example->strand,
		-primary_tag    => 'transcript',
		-source         => $example->source_tag,
		-name           => $example->display_name,
		-segments       => \@new,
	);
}

sub get_transcript_length {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref($transcript) =~ /seqfeature/i;
	my $total = 0;
	foreach my $e (get_exons($transcript)) {
		$total += $e->length;
	}
	return $total;
}



######## CDS Methods

sub is_coding {
	my $transcript = shift;
	return unless $transcript;
	if ($transcript->primary_tag =~ /gene/i) {
		# someone passed a gene, check its subfeatures
		my $code_potential = 0;
		foreach ($transcript->get_SeqFeatures) {
			$code_potential += is_coding($_);
		}
		return $code_potential;
	}
	return 1 if $transcript->primary_tag =~ /mrna/i; # assumption
	return 1 if $transcript->source =~ /protein.?coding/i;
	if ($transcript->has_tag('biotype')) {
		# ensembl type GFFs
		my ($biotype) = $transcript->get_tag_values('biotype');
		return 1 if $biotype =~ /protein.?coding/i;
	}
	elsif ($transcript->has_tag('gene_biotype')) {
		# ensembl type GTFs
		my ($biotype) = $transcript->get_tag_values('gene_biotype');
		return 1 if $biotype =~ /protein.?coding/i;
	}
	foreach ($transcript->get_SeqFeatures) {
		# old fashioned way
		return 1 if $_->primary_tag eq 'CDS';
	}
	return 0;
}

sub get_cds {
	my $transcript = shift;
	my @cds;
	foreach my $subfeat ($transcript->get_SeqFeatures) {
		push @cds, $subfeat if $subfeat->primary_tag eq 'CDS';
	}
	return unless @cds;
	@cds = map { $_->[0] } 
			sort { $a->[1] <=> $b->[1] }
			map { [$_, $_->start] } 
			@cds;
	return wantarray ? @cds : \@cds;
}

sub get_cdsStart {
	my $transcript = shift;
	my $cds = get_cds($transcript);
	return unless $cds;
	if ($transcript->strand >= 0) {
		return $cds->[0]->start;
	}
	else {
		# stop codons may or may not be not included
		my $codon = get_stop_codon($transcript);
		return $codon->start < $cds->[0]->start ? 
			$codon->start : $cds->[0]->start;
	}
}

sub get_cdsEnd {
	my $transcript = shift;
	my $cds = get_cds($transcript);
	return unless $cds;
	if ($transcript->strand >= 0) {
		my $codon = get_stop_codon($transcript);
		return $codon->end > $cds->[-1]->end ? 
			$codon->end : $cds->[-1]->end;
	}
	else {
		return $cds->[-1]->end;
	}
}

sub get_transcript_cds_length {
	my $transcript = shift;
	my $total = 0;
	foreach my $subf ($transcript->get_SeqFeatures) {
		next unless $subf->primary_tag eq 'CDS';
		$total += $subf->length;
	}
	return $total;
}

sub get_start_codon {
	my $transcript = shift;
	return unless is_coding($transcript);
	my $start_codon;
	
	# look for existing one
	foreach my $subfeat ($transcript->get_SeqFeatures) {
		$start_codon =  $subfeat if $subfeat->primary_tag =~ /start.?codon/i;
	}
	return $start_codon if $start_codon;
	
	# otherwise we have to build one
	my $cdss = get_cds($transcript);
	return unless $cdss;
	if ($transcript->strand >= 0) {
		$start_codon = $transcript->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'start_codon',
				-start         => $cdss->[0]->start,
				-end           => $cdss->[0]->start + 2,
				-strand        => 1,
				-phase         => 0,
				-primary_id    => $transcript->primary_id . '.start_codon',
		);
	}
	else {
		$start_codon = $transcript->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'start_codon',
				-start         => $cdss->[-1]->end - 2,
				-end           => $cdss->[-1]->end,
				-strand        => -1,
				-phase         => 0,
				-primary_id    => $transcript->primary_id . '.start_codon',
		);
	}
	return $start_codon;
}

sub get_stop_codon {
	my $transcript = shift;
	return unless is_coding($transcript);
	my $stop_codon;
	
	# look for existing one
	foreach my $subfeat ($transcript->get_SeqFeatures) {
		$stop_codon =  $subfeat if $subfeat->primary_tag =~ /stop.?codon/i;
	}
	return $stop_codon if $stop_codon;
	
	# otherwise we have to build one
	# this entirely presumes that the stop codon is inclusive to the last cds
	my $cdss = get_cds($transcript);
	return unless $cdss;
	if ($transcript->strand >= 0) {
		$stop_codon = $transcript->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'stop_codon',
				-start         => $cdss->[-1]->end - 2,
				-end           => $cdss->[-1]->end,
				-strand        => 1,
				-phase         => 0,
				-primary_id    => $transcript->primary_id . '.stop_codon',
		);
	}
	else {
		$stop_codon = $transcript->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'stop_codon',
				-start         => $cdss->[0]->start,
				-end           => $cdss->[0]->start + 2,
				-strand        => -1,
				-phase         => 0,
				-primary_id    => $transcript->primary_id . '.stop_codon',
		);
	}
	return $stop_codon;
}

sub get_utrs {
	my $transcript = shift;
	confess "not a SeqFeature object!" unless ref($transcript) =~ /seqfeature/i;
	return unless is_coding($transcript);
	
	# collect the various types of subfeatures
	my @exons;
	my @cdss;
	my @utrs;
	my @transcripts;
	foreach my $subfeat ($transcript->get_SeqFeatures) {
		my $type = $subfeat->primary_tag;
		if ($type =~ /exon/i) {
			push @exons, $subfeat;
		}
		elsif ($type =~ /cds/i) {
			push @cdss, $subfeat;
		}
		elsif ($type =~ /utr|untranslated/i) {
			push @utrs, $subfeat;
		}
		elsif ($type =~ /rna|transcript/i) {
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
	elsif (@exons and @cdss) {
		# calculate the utrs ourselves
		
		@exons = map { $_->[0] }
				sort { $a->[1] <=> $b->[1] }
				map { [$_, $_->start] } 
				@exons;
		@cdss = map { $_->[0] }
				sort { $a->[1] <=> $b->[1] }
				map { [$_, $_->start] } 
				@cdss;
		my $firstCDS = $cdss[0];
		my $lastCDS  = $cdss[-1];
		while (@exons) {
			my $exon = shift @exons;
			if ($exon->end < $firstCDS->start) {
				# whole exon is UTR
				my $utr = _duplicate($exon);
				$utr->primary_tag( 
					$transcript->strand >= 0 ? 'five_prime_UTR' : 'three_prime_UTR' );
				$utr->display_name( $exon->display_name . '.utr' );
				push @list, $utr;
			}
			elsif ($exon->overlaps($firstCDS)) {
				# partial UTR on left side
				my $pieces = $exon->subtract($firstCDS); # array ref of pieces
				next unless $pieces;
				my $utr = $pieces->[0]; # we will want the first one if there are two
				$utr->primary_tag( 
					$transcript->strand >= 0 ? 'five_prime_UTR' : 'three_prime_UTR' );
				$utr->display_name($exon->display_name . '.utr');
				$utr->strand($exon->strand);
				$utr->source($exon->strand);
				push @list, $utr;
			}
			elsif ($exon->start > $firstCDS->end and $exon->end < $lastCDS->start) {
				# CDS exon
				next;
			}
			elsif ($exon->overlaps($lastCDS)) {
				# partial UTR
				my $pieces = $exon->subtract($lastCDS); # array ref of pieces
				next unless $pieces;
				my $utr = $pieces->[-1]; # we will want the second one if there are two
				$utr->primary_tag( 
					$transcript->strand >= 0 ? 'five_prime_UTR' : 'three_prime_UTR' );
				$utr->display_name($exon->display_name . '.utr');
				$utr->strand($exon->strand);
				$utr->source($exon->strand);
				push @list, $utr;
			}
			elsif ($exon->start > $lastCDS->end) {
				# whole exon is UTR
				my $utr = _duplicate($exon);
				$utr->primary_tag( 
					$transcript->strand >= 0 ? 'five_prime_UTR' : 'three_prime_UTR' );
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


#### Export methods

sub gff_string {
	# Bio::ToolBox::SeqFeature and Bio::SeqFeature::Lite objects have this method
	# otherwise this will die
	return shift->gff_string(@_);
}

sub ucsc_string {
	my $feature = shift;
	confess "not a SeqFeature object!" unless ref($feature) =~ /seqfeature/i;
	my @ucsc_list;
	
	# process according to type
	if ($feature->primary_tag =~ /gene/i) {
		# a gene object, we will need to process it's transcript subfeatures
		foreach my $transcript (get_transcripts($feature)) {
			my $ucsc = _process_ucsc_transcript($transcript, $feature);
			push @ucsc_list, $ucsc if $ucsc;
		}
	}
	elsif ($feature->primary_tag =~ /rna|transcript/i) {
		# some sort of RNA transcript
		my $ucsc = _process_ucsc_transcript($feature);
		push @ucsc_list, $ucsc if $ucsc;
	}
	
	# return strings
	my $string;
	foreach my $ucsc (@ucsc_list) {
		$string .= sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
			$ucsc->{'name2'},
			$ucsc->{'name'},
			$ucsc->{'chr'},
			$ucsc->{'strand'},
			$ucsc->{'txStart'},
			$ucsc->{'txEnd'},
			$ucsc->{'cdsStart'},
			$ucsc->{'cdsEnd'},
			$ucsc->{'exonCount'},
			join(",", @{$ucsc->{'exonStarts'}} ) . ',', 
			join(",", @{$ucsc->{'exonEnds'}} ) . ',',
		);
	}
	return $string;
}





#### internal methods

# internal method to duplicate a seqfeature object with just the essential stuff
# sometimes we just don't want to bring along all this cruft....
sub _duplicate {
	my $f = shift;
	return $f->new(
		-seq_id         => $f->seq_id,
		-start          => $f->start,
		-end            => $f->end,
		-strand         => $f->strand,
		-primary_tag    => $f->primary_tag,
		-source         => $f->source_tag,
	);
}

sub _process_ucsc_transcript {
	my $transcript = shift;
	my $gene = shift || undef;
	
	# initialize ucsc hash
	my $ucsc = {
		'name'         => $transcript->display_name || $transcript->primary_id,
		'name2'        => undef,
		'chr'          => $transcript->seq_id,
		'strand'       => $transcript->strand < 0 ? '-' : '+', 
		'txStart'      => $transcript->start - 1,
		'txEnd'        => $transcript->end,
		'cdsStart'     => undef,
		'cdsEnd'       => undef,
		'exonCount'    => 0,
		'exonStarts'   => [],
		'exonEnds'     => [],
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
	if (is_coding($transcript)) {
		$ucsc->{cdsStart} = get_cdsStart($transcript) - 1;
		$ucsc->{cdsEnd}   = get_cdsEnd($transcript);
	}
	else {
		# non-coding transcript sets the cds points to the transcript end
		$ucsc->{cdsStart} = $transcript->end;
		$ucsc->{cdsEnd}   = $transcript->end;
	}
	
	# record the exons
	foreach my $exon (get_exons($transcript)) {
		push @{ $ucsc->{'exonStarts'} }, $exon->start - 1;
		push @{ $ucsc->{'exonEnds'} }, $exon->end;
		$ucsc->{'exonCount'} += 1;
	}
	
	return $ucsc;
}



__END__

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
