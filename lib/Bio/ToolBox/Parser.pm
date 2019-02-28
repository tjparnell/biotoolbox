package Bio::ToolBox::Parser;
my $VERSION = 1.70;

=head1 NAME

Bio::ToolBox::Parser - generic parsing tool for GFF, UCSC, BED

=head1 SYNOPSIS
  
  # obtain an annotation file
  use Bio::ToolBox::Parser;
  my $filename = shift @ARGV; # could be any annotation format
  
  # open in parser
  my $parser = Bio::ToolBox::Parser->new(
  	file    => $filename,
  ) or die "unable to open $filename!\n";
  # file is tasted and appropriate parser automatically selected
  # returns parser object if recognized
  # could be one of Bio::ToolBox::parser::bed, 
  # Bio::ToolBox::parser::gff, or Bio::ToolBox::parser::ucsc
  
  # do something with parser
  while (my $feature = $parser->next_top_feature() ) {
	# each $feature is a SeqFeature object
	my @children = $feature->get_SeqFeatures();
  }
  
  # alternatively open a specific parser
  use Bio::ToolBox::parser::gff;
  my $parser = Bio::ToolBox::parser::gff->new(
  	file    => $filename,
  	do_gene => 1,
  	do_exon => 1,
  ) or die "unable to open gff file!\n";


=head1 DESCRIPTION

This module is a generic wrapper around the three main annotation file 
parsers. It will taste test the file and choose the appropriate parser and 
open it automatically. These parsers include 

=over 4

=item L<Bio::ToolBox::parser::bed>

Parses most Bed file formats, including 3-12 column Bed formats, and some 
specific Encode Bed6+? formats, including C<narrowPeak>, C<broadPeak>, and 
C<gappedPeak>. 

=item L<Bio::ToolBox::parser::gff>

Parses any GFF flavor, including GTF and GFF3.

=item L<Bio::ToolBox::parser::ucsc>

Parses some of the common UCSC annotation table formats, including 
C<refFlat>, C<genePred>, C<genePredExt>, and C<knownGene>. Support for 
some additional UCSC metadata tables is available.

=back

Files are parsed entirely into memory, assembling gene components (transcripts, 
exon, CDS, UTR, etc) into hierarchical, top-level SeqFeature objects as 
appropriate. These SeqFeature objects can then be iterated through in a loop, 
acting on each one as appropriate. The default SeqFeature class is 
L<Bio::ToolBox::SeqFeature>, an efficient L<Bio::SeqFeatureI> compliant 
object class.

=head1 METHODS






=cut

use strict;
use Module::Load;
use Bio::ToolBox::Data;
use Carp qw(carp cluck croak confess);
use Bio::ToolBox::SeqFeature; # alternative to Bio::SeqFeature::Lite

1;

sub new {
	my $class = shift;
	
	# passed arguments
	my %args;
	if (scalar(@_) == 1) {
		$args{file} = shift;
	}
	else {
		%args = @_;
	}
	$args{file} ||= $args{table} || undef;
	
	# check parser flavor
	my $flavor;
	if (ref($class)) {
		$class = ref($class);
	}
	if ($class =~ /Bio::ToolBox::parser::(\w+)/) {
		# we got the flavor directly from a flavored parser new function
		$flavor = $1;
		# presume this parser is already loaded
	}
	elsif ($args{file}) {
		# find the flavor from the file
		$flavor = Bio::ToolBox::Data->taste_file($args{file}) or return;
		$class = 'Bio::ToolBox::parser::' . $flavor;
		eval {load $class};
		if ($@) {
			carp "unable to load $class! cannot parse $flavor!";
			return;
		}
	}
	else {
		# uh oh! what flavor to use????
		# delay dying here? 
		confess("Unknown flavor! See documentation. ");
	}
	
	# initialize self
	my $self = {
		'fh'            => undef,
		'version'       => 0,
		'do_gene'       => 1, 
		'do_exon'       => 0,
		'do_cds'        => 0, 
		'do_utr'        => 0, 
		'do_codon'      => 0,
		'do_name'       => 0,
		'share'         => 0,
		'simplify'      => 0,
		'source'        => undef,
		'typelist'      => '',
		'seq_ids'       => {},
		'top_features'  => [],
		'eof'           => 0,
		'comments'      => [],
		'convertor_sub' => undef,
		'sfclass'       => 'Bio::ToolBox::SeqFeature',
	};
	bless $self, $class;
	
	# parser specific keys
	if ($flavor eq 'bed') {
		$self->{stream}         = undef;
		$self->{line_count}     = 0;
		$self->{bed}            = undef;
	}
	elsif ($flavor eq 'gff') {
		$self->{gff3}           = 0;
		$self->{gtf}            = 0;
		$self->{typelist}       = '';
		$self->{orphans}        = [];
		$self->{loaded}         = {};
		$self->{duplicate_ids}  = {};
	}
	elsif ($flavor eq 'ucsc') {
		$self->{share}          = 1;
		$self->{source}         = 'UCSC';
		$self->{gene2seqf}      = {};
		$self->{id2count}       = {};
		$self->{counts}         = {};
		$self->{refseqsum}      = {};
		$self->{refseqstat}     = {};
		$self->{kgxref}         = {};
		$self->{ensembldata}    = {};
		$self->{line_count}     = 0;
	}
	
	# process remaining arguments
	if (exists $args{simplify}) {
		$self->simplify( $args{simplify} );
	}
	if (exists $args{do_gene}) {
		$self->do_gene($args{do_gene});
	}
	if (exists $args{do_exon}) {
		$self->do_exon($args{do_exon});
	}
	if (exists $args{do_cds}) {
		$self->do_cds($args{do_cds});
	}
	if (exists $args{do_utr}) {
		$self->do_utr($args{do_utr});
	}
	if (exists $args{do_codon}) {
		$self->do_codon($args{do_codon});
	}
	if (exists $args{version}) {
		$self->{version} = $args{version};
	}
	if (exists $args{source}) {
		$self->source($args{source});
	}
	if (exists $args{class}) {
		my $class = $args{class};
		if (eval "require $class; 1") {
			$self->{sfclass} = $class;
		}
		else {
			croak $@;
		}
	}
	if ($flavor eq 'ucsc') {
		# lots of accessory files for UCSC tables
		if (exists $args{refseqsum}) {
			$self->load_extra_data($args{refseqsum}, 'refseqsum');
		}
		elsif (exists $args{summary}) {
			$self->load_extra_data($args{summary}, 'refseqsum');
		}
		if (exists $args{refseqstat}) {
			$self->load_extra_data($args{refseqstat}, 'refseqstat');
		}
		elsif (exists $args{status}) {
			$self->load_extra_data($args{status}, 'refseqstat');
		}
		if (exists $args{kgxref}) {
			$self->load_extra_data($args{kgxref}, 'kgxref');
		}
		if (exists $args{ensembltogenename}) {
			$self->load_extra_data($args{ensembltogenename}, 'ensembltogene');
		}
		elsif (exists $args{ensname}) {
			$self->load_extra_data($args{ensname}, 'ensembltogene');
		}
		if (exists $args{ensemblsource}) {
			$self->load_extra_data($args{ensemblsource}, 'ensemblsource');
		}
		elsif (exists $args{enssrc}) {
			$self->load_extra_data($args{enssrc}, 'ensemblsource');
		}
	}
	
	# open the file
	if ($args{file}) {
		$self->open_file( $args{file} ) or croak "unable to open file!";
	}
	
	# finished
	return $self;
}

sub do_gene {
	my $self = shift;
	return 0 if ref($self) eq 'Bio::ToolBox::parser::bed';
	if (@_) {
		$self->{'do_gene'} = shift;
	}
	return $self->{'do_gene'};
}	

sub do_exon {
	my $self = shift;
	if (@_) {
		$self->{'do_exon'} = shift;
	}
	return $self->{'do_exon'};
}	

sub do_cds {
	my $self = shift;
	if (@_) {
		$self->{'do_cds'} = shift;
	}
	return $self->{'do_cds'};
}	

sub do_utr {
	my $self = shift;
	if (@_) {
		$self->{'do_utr'} = shift;
	}
	return $self->{'do_utr'};
}	

sub do_codon {
	my $self = shift;
	if (@_) {
		$self->{'do_codon'} = shift;
	}
	return $self->{'do_codon'};
}	

sub do_name {
	my $self = shift;
	return 0 unless ref($self) eq 'Bio::ToolBox::parser::ucsc';
	# does nothing with gff and bed
	if (@_) {
		$self->{'do_name'} = shift;
	}
	return $self->{'do_name'};
}	

sub share {
	my $self = shift;
	return 0 unless ref($self) eq 'Bio::ToolBox::parser::ucsc';
	# does nothing with gff and bed
	if (@_) {
		$self->{'share'} = shift;
	}
	return $self->{'share'};
}	

sub simplify {
	my $self = shift;
	return 0 unless ref($self) eq 'Bio::ToolBox::parser::gff';
	# does nothing with ucsc and bed
	if (defined $_[0]) {
		$self->{simplify} = shift;
	}
	return $self->{simplify};
}

sub source {
	my $self = shift;
	return if ref($self) eq 'Bio::ToolBox::parser::gff';
	if (@_) {
		$self->{'source'} = shift;
	}
	return $self->{'source'};
}

sub version {
	return shift->{version};
}

sub fh {
	my $self = shift;
	if (@_) {
		$self->{'fh'} = shift;
	}
	return $self->{'fh'};
}

sub seq_ids {
	my $self = shift;
	unless (scalar keys %{$self->{seq_ids}}) {
		$self->_get_seq_ids;
	}
	my @s = keys %{$self->{seq_ids}};
	return wantarray ? @s : \@s;
}

sub seq_id_lengths {
	my $self = shift;
	unless (scalar keys %{$self->{seq_ids}}) {
		$self->_get_seq_ids;
	}
	return $self->{seq_ids};
}

sub _get_seq_ids {
	my $self = shift;
	return unless $self->{'eof'};
	foreach (@{ $self->{top_features} }) {
		my $s = $_->seq_id;
		unless (exists $self->{seq_ids}{$s}) {
			$self->{seq_ids}{$s} = 1;
		}
		$self->{seq_ids}{$s} = $_->end if $_->end > $self->{seq_ids}{$s};
	}
}

sub next_top_feature {
	my $self = shift;
	
	# return next item
	if (exists $self->{top_feature_index}) {
		my $i = $self->{top_feature_index};
		if ($i == $self->{last_top_feature_index}) {
			delete $self->{top_feature_index};
			delete $self->{last_top_feature_index};
			return undef;
		}
		$self->{top_feature_index} += 1;
		return $self->{top_features}->[$i];
	}
	# otherwise, it's our first time
	
	# check that we have an open filehandle
	unless ($self->fh) {
		croak("no annotation file loaded to parse!");
	}
	unless ($self->{'eof'}) {
		$self->parse_file or croak "unable to parse file!";
	}
	
	# set up index
	$self->{top_feature_index} = 1;
	$self->{last_top_feature_index} = scalar(@{$self->{top_features}});
	return $self->{top_features}->[0];
}

sub top_features {
	my $self = shift;
	unless ($self->{'eof'}) {
		$self->parse_file;
	}
	my @features = @{ $self->{top_features} };
	return wantarray ? @features : \@features;
}



__END__

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  




