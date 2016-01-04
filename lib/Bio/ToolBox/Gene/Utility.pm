package Bio::ToolBox::Gene::Utility;
our $VERSION = '1.35';

=head1 NAME

Bio::ToolBox::Gene::Utility - Methods for working with gene SeqFeatures

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=over 4

=item is_coding($transcript)

This method will return a boolean value if the passed transcript object 
appears to be a coding transcript. GFF and GTF files are not always immediately 
clear about the type of transcript; there are (unfortunately) multiple ways 
to encode the feature as a protein coding transcript: primary_tag, source_tag, 
attribute, CDS subfeatures, etc. This method tries to determine this.

=back

=cut

use strict;
use Carp qw(carp cluck croak confess);
require Exporter;


### Variables
# Export
our @ISA = qw(Exporter);
our @EXPORT = qw(
);
our @EXPORT_OK = qw(
);

### The True Statement
1; 



### Coding Status
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

sub get_exons {
	
	# initialize
	my $transcript = shift;
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
			push @transcripts;
		}
	}
	
	# check which array we'll use
	# prefer to use actual exon subfeatures, but those may not be defined
	my @list;
	if (@exons) {
		@list = @exons;
	}
	elsif (@cdss) {
		@list = @cdss;
	}
	elsif (@transcripts) {
		foreach my $t (@transcripts) {
			# there are possibly duplicates in here if there are alternate transcripts
			# should we remove them?
			push @list, collect_exons($t);
		}
	}
	else {
		# nothing found!
		return;
	}
	
	# return sorted list by start position
	return  map { $_->[0] }
			sort { $a->[1] <=> $b->[1] }
			map { [$_, $_->start] } 
			@list;
}

sub get_alt_exons {
	my $ca_exons = get_alt_common_exons(@_);
	my @alts;
	foreach my $k (keys %$ca_exons) {
		next if $k eq 'common';
		next if $k eq 'uncommon';
		push @alts, @{ $ca_exons->{$k} };
	}
	return wantarray ? @alts : \@alts;
}

sub get_common_exons {
	my $ca_exons = get_alt_common_exons(@_);
	my @com = @{ $ca_exons->{common} };
	return wantarray ? @com : \@com;
}

sub get_uncommon_exons {
	my $ca_exons = get_alt_common_exons(@_);
	my @com = @{ $ca_exons->{uncommon} };
	return wantarray ? @com : \@com;
}

sub get_alt_common_exons {
	return _get_alt_common_things('exon', @_);
}


sub get_introns {
	my $transcript = shift;
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
		@introns = map { $_->[0] }
				sort { $a->[1] <=> $b->[1] }
				map { [$_, $_->start] } 
				@introns;
	}
	
	# finished
	return @introns;
}


sub get_alt_introns {
	my $ca_introns = get_alt_common_introns(@_);
	my @alts;
	foreach my $k (keys %$ca_introns) {
		next if $k eq 'common';
		next if $k eq 'uncommon';
		push @alts, @{ $ca_introns->{$k} };
	}
	return wantarray ? @alts : \@alts;
}

sub get_common_introns {
	my $ca_introns = get_alt_common_introns(@_);
	my @com = @{ $ca_introns->{common} };
	return wantarray ? @com : \@com;
}

sub get_uncommon_introns {
	my $ca_introns = get_alt_common_introns(@_);
	my @com = @{ $ca_introns->{uncommon} };
	return wantarray ? @com : \@com;
}

sub get_alt_common_introns {
	return _get_alt_common_things('intron', @_);
}

sub _get_alt_common_things {
	# internal subroutine to get either exons or introns
	my $type = shift; # exon or intron
	my @transcripts;
	return unless @_;
	if (scalar @_ == 1) {
		# someone passed a gene, get the transcripts
		@transcripts = get_transcripts($_[0]);
		return if scalar @transcripts == 1;
	}
	elsif (scalar @_ > 1) {
		@transcripts = @_;
	}
	
	# type of thing to get
	my $do_exon = $type eq 'exon' ? 1 : 0;
	
	# hash of transcript to things
	my %tx2things = (
		common => [],
		uncommon => []
	);
	
	# get things and put them in has based on coordinates
	my %pos2things;
	foreach my $t (@transcripts) {
		my @things = $do_exon ? get_exons($t) : get_introns($t);
		foreach my $e (@things) {
			push @{ $pos2things{$e->start}{ $e->end} }, [ $t->display_name, $e ];
		}
		$tx2things{ $t->display_name } = [];
	}
	
	# put things into categories based on commonality
	# associate things with unique transcripts, common, or uncommon sets
	my $trx_number = scalar @transcripts;
	foreach my $s (sort {$a <=> $b} keys %pos2things) {               # sort on start
		foreach my $e (sort {$a <=> $b} keys %{ $pos2things{$s} }) {  # sort on stop
			my $n = scalar @{ $pos2things{$s}{$e} };
			if ($n == 1) {
				# only 1 thing, must be an alternate
				push @{ $tx2things{ $pos2things{$s}{$e}[0][0] } }, $pos2things{$s}{$e}[0][1];
			}
			elsif ($n == $trx_number) {
				# common to all transcripts
				push @{ $tx2things{common} }, $pos2things{$s}{$e}[0][1];
			}
			else {
				# common to some but not all transcripts, so uncommon
				push @{ $tx2things{uncommon} }, $pos2things{$s}{$e}[0][1];
			}
		}
	}
	return \%tx2things;
}


sub collapse_transcripts {

}

sub get_transcript_length {

}

sub get_cdsStart {

}

sub get_cdsEnd {

}

sub get_transcripts {
	my $gene = shift;
	my @transcripts;
	foreach my $subf ($gene->get_SeqFeatures) {
		push @transcripts, $subf if $subf->primary_tag =~ /rna|transcript/i;
	}
	return map { $_->[0] }
		sort { $a->[1] <=> $b->[1] }
		map { [$_, $_->start] } 
		@transcripts;
}


__END__

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
