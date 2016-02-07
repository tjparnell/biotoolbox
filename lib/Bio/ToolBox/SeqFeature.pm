package Bio::ToolBox::SeqFeature;
our $VERSION = '1.35';

=head1 NAME

Bio::ToolBox::SeqFeature - Fast, simple SeqFeature implementation

=head1 DESCRIPTION

SeqFeature objects represent functional elements on a genomic or chromosomal 
sequence, such as genes, transcripts, exons, etc. In many cases, especially 
genes, they have a hierarchical structure, typically in this order

    gene
      mRNA or transcript
        exon
        CDS

SeqFeature objects have at a minimum coordinate information, including 
chromosome, start, stop, and strand, and a name or unique identifier. They 
often also have type or source information, which usually follows the 
Sequence Ontology key words.

This is a fast, efficient, simplified SeqFeature implementation that mostly 
implements the L<Bio::SeqFeatureI> API, and could be substituted in for other 
implementations, such L<Bio::SeqFeature::Lite> and L<Bio::SeqFeature::Generic>. 
Unlike the others, however, it inherits no classes or methods and uses an 
unorthodox blessed array to store feature attributes. 

=head1 LIMITATIONS

Because of their underlying array structure, Bio::ToolBox::SeqFeature objects 
should generally not be used as a base class (unless you know the ramifications 
of doing so). The following Bio classes and Interfaces may not work, either because 
they have not been implemented, object structure incompatibility, or simply haven't 
been tested yet for compatibility (it's possible some can be used).

=over 4

=item Bio::AnnotationI

=item Bio::LocationI

=item Bio::RangeI

=item Bio::SeqI

=item Bio::Tools::GFF

=item Bio::DB::SeqFeature::Store

=item Bio::Graphics

=back

=head1 METHODS





=cut

use strict;
use Carp qw(carp cluck croak confess);
use constant {
	SEQID   => 0,
	START   => 1,
	STOP    => 2,
	STRND   => 3,
	NAME    => 4,
	ID      => 5,
	TYPE    => 6,
	SRC     => 7,
	SCORE   => 8,
	PHASE   => 9,	
	SUBF    => 10,
	ATTRB   => 11,
	PARNT   => 12,
};
our $IDCOUNT = 0;

1;

#### Aliases ####
# to maintain compatibility with Bio::SeqFeature::Lite and Bio::SeqFeatureI we 
# put in plenty of aliases to some of the methods
*stop = \&end;
*name = \&display_name;
*id = \&primary_id;
*method = \&primary_tag;
*source = \&source_tag;
*add_segment = \&add_SeqFeature;
*get_all_SeqFeatures = *segments = \&get_SeqFeatures;
*each_tag_value = \&get_tag_values;
*get_all_tags = \&all_tags;
*gff3_string = \&gff_string;


#### METHODS ####
sub new {
	my $class = shift;
	if (ref($class)) {
		$class = ref($class);
	}
	my %args = @_;
	my $self = bless [], $class; 
		# bless ourselves early to take advantage of some methods
		# but otherwise do what we can simply and quickly
	
	# primary options
	$self->[SEQID] = $args{-seq_id} || $args{-seqid} || $args{'-ref'} || 
		$args{chrom} || undef;
	$self->[START] = $args{-start} || undef;
	$self->[STOP] = $args{-end} || $args{-stop} || undef;
	$self->strand($args{-strand}) if exists $args{-strand};
	$self->[NAME] = $args{-display_name} || $args{-name} || undef;
	$self->[ID] = $args{-primary_id} || $args{-id} || undef;
	
	# check orientation
	if (defined $self->[START] and defined $self->[STOP] and 
		$self->[START] > $self->[STOP]
	) {
		# flip the coordinates around
		($self->[START], $self->[STOP]) = ($self->[STOP], $self->[START]);
		$self->[STRND] *= -1;
	}
	
	# additional options
	$args{-type} ||= $args{-primary_tag} || undef;
	if (defined $args{-type}) {
		$self->type($args{-type});
	}
	$args{-source} ||= $args{-source_tag} || undef;
	if (defined $args{-source}) {
		$self->[SRC] = $args{-source};
	}
	if (exists $args{-score}) {
		$self->[SCORE] = $args{-score};
	}
	if (exists $args{-phase}) {
		$self->phase($args{-phase});
	}
	if (exists $args{-attributes} or exists $args{-tags}) {
		$args{-attributes} ||= $args{-tags};
		if (ref($args{-attributes}) eq 'HASH') {
			$self->[ATTRB] = $args{-attributes};
		}
	}
	if (exists $args{-segments}) {
		$self->[SUBF] = [];
		foreach my $s (@{ $args{-segments} }) {
			unless (ref($s) eq $class) {
				croak "segments should be passed as $class objects!";
			}
			push @{$self->[SUBF]}, $s;
		}
	}
	
	return $self; 
}


sub seq_id {
	my $self = shift;
	if (@_) {
		$self->[SEQID] = $_[0];
	}
	my $seqid = defined $self->[SEQID] ? $self->[SEQID] : 
		defined $self->[PARNT] ? $self->[PARNT]->seq_id : undef;
	return $seqid;
}

sub start {
	my $self = shift;
	if (@_) {
		$self->[START] = $_[0];
	}
	return $self->[START];
}

sub end {
	my $self = shift;
	if (@_) {
		$self->[STOP] = $_[0];
	}
	return $self->[STOP];
}

sub strand {
	my $self = shift;
	if (@_) {
		$self->[STRND] = $_[0];
		$self->[STRND] = 1 if $self->[STRND] eq '+';
		$self->[STRND] = -1 if $self->[STRND] eq '-';
		$self->[STRND] = 0 if $self->[STRND] eq '.';
	}
	my $strand = defined $self->[STRND] ? $self->[STRND] : 
		defined $self->[PARNT] ? $self->[PARNT]->strand : 0;
	return $strand;
}

sub display_name {
	my $self = shift;
	if (@_) {
		$self->[NAME] = $_[0];
	}
	return $self->[NAME];
}

sub primary_id {
	my $self = shift;
	if (@_) {
		$self->[ID] = $_[0];
	}
	else {
		# automatically assign a new ID
		$self->[ID] = sprintf("%s%s.%09d", $self->primary_tag, $$, $IDCOUNT++);
	}
	return $self->[ID];
}

sub primary_tag {
	my $self = shift;
	if (@_) {
		$self->[TYPE] = $_[0];
	}
	$self->[TYPE] ||= 'region';
	return $self->[TYPE];
}

sub source_tag {
	my $self = shift;
	if (@_) {
		$self->[SRC] = $_[0];
	}
	my $source = defined $self->[SRC] ? $self->[SRC] : 
		defined $self->[PARNT] ? $self->[PARNT]->source_tag : undef;
	return $source;
}

sub type {
	my $self = shift;
	if (@_) {
		my $type = $_[0];
		if ($type =~ /:/) {
			my ($t, $s) = split /:/, $type;
			$self->[TYPE] = $t;
			$self->[SRC] = $s;
		}
		else {
			$self->[TYPE] = $type;
		}
	}
	$self->[TYPE] ||= undef;
	$self->[SRC] ||= undef;
	if (defined $self->[SRC]) {
		return join(':', $self->[TYPE], $self->[SRC]);
	}
	return $self->[TYPE];
}

sub score {
	my $self = shift;
	if (@_) {
		$self->[SCORE] = $_[0];
	}
	$self->[SCORE] ||= 0;
	return $self->[SCORE];
}

sub phase {
	my $self = shift;
	if (@_) {
		my $p = $_[0];
		unless ($p =~ /^[012\.]$/) {
			cluck "phase must be 0, 1, 2 or .!";
		}
		$self->[PHASE] = $p;
	}
	$self->[PHASE] ||= '.';
	return $self->[PHASE];
}

sub add_SeqFeature {
	my $self = shift;
	$self->[SUBF] ||= [];
	my $count = 0;
	foreach my $s (@_) {
		if (ref($s)) {
			push @{ $self->[SUBF] }, $s;
			$count++;
		}
		else {
			cluck "please use Seqfeature objects when adding sub features!";
		}
	}
	return $count;
}

sub get_SeqFeatures {
	my $self = shift;
	$self->[SUBF] ||= [];
	my @children;
	foreach (@{ $self->[SUBF] }) {
		$_->[PARNT] = $self;
		push @children, $_;
	}
	return wantarray ? @children : \@children;
}

sub add_tag_value {
	my ($self, $key, $value) = @_;
	return unless ($key and $value);
	$self->[ATTRB] ||= {};
	if (exists $self->[ATTRB]->{$key}) {
		if (ref($self->[ATTRB]->{$key}) eq 'ARRAY') {
			push @{ $self->[ATTRB]->{$key} }, $value;
		}
		else {
			my $current = $self->[ATTRB]->{$key};
			$self->[ATTRB]->{$key} = [$current, $value];
		}
	}
	else {
		$self->[ATTRB]->{$key} = $value;
	}
}

sub has_tag {
	my ($self, $key) = @_;
	$self->[ATTRB] ||= {};
	return exists $self->[ATTRB]->{$key};
}

sub get_tag_values {
	my ($self, $key) = @_;
	$self->[ATTRB] ||= {};
	if (exists $self->[ATTRB]->{$key}) {
		if (ref($self->[ATTRB]->{$key}) eq 'ARRAY') {
			return wantarray ? @{ $self->[ATTRB]->{$key} } : $self->[ATTRB]->{$key};
		}
		else {
			return $self->[ATTRB]->{$key};
		}
	}
	else {
		return;
	}
}

sub attributes {
	my $self = shift;
	$self->[ATTRB] ||= {};
	return wantarray ? @{ $self->[ATTRB] } : $self->[ATTRB];
}

sub all_tags {
	my $self = shift;
	$self->[ATTRB] ||= {};
	my @k = keys %{ $self->[ATTRB] };
	return wantarray ? @k : \@k;
}

sub parent {
	my $self = shift;
	return defined $self->[PARNT] ? $self->[PARNT] : undef; 
}




### Range Methods
# Borrowed from Bio::RangeI, but does not do strong/weak strand checks

sub length {
	my $self = shift;
	return $self->end - $self->start + 1;
}

sub overlaps {
	my ($self, $other) = @_;
	return unless ($other and ref($other));
	return unless ($self->seq_id eq $other->seq_id);
	return not (
		$self->start > $other->end or
		$self->end   < $other->start
	);
}

sub contains {
	my ($self, $other) = @_;
	return unless ($other and ref($other));
	return unless ($self->seq_id eq $other->seq_id);
	return (
		$other->start >= $self->start and
		$other->end   <= $self->end
	);
}

sub equals {
	my ($self, $other) = @_;
	return unless ($other and ref($other));
	return unless ($self->seq_id eq $other->seq_id);
	return (
		$other->start == $self->start and
		$other->end   == $self->end
	);
}

sub intersection {
	my ($self, $other) = @_;
	return unless ($other and ref($other));
	return unless ($self->seq_id eq $other->seq_id);
	my ($start, $stop);
	if ($self->start >= $other->start) {
		$start = $self->start;
	}
	else {
		$start = $other->start;
	}
	if ($self->end <= $other->end) {
		$stop = $self->end;
	}
	else {
		$stop = $other->end;
	}
	return if $start > $stop;
	return $self->new(
		-seq_id     => $self->seq_id,
		-start      => $start,
		-end        => $stop,
	);
}

sub union {
	my ($self, $other) = @_;
	return unless ($other and ref($other));
	return unless ($self->seq_id eq $other->seq_id);
	my ($start, $stop);
	if ($self->start <= $other->start) {
		$start = $self->start;
	}
	else {
		$start = $other->start;
	}
	if ($self->end >= $other->end) {
		$stop = $self->end;
	}
	else {
		$stop = $other->end;
	}
	return if $start > $stop;
	return $self->new(
		-seq_id     => $self->seq_id,
		-start      => $start,
		-end        => $stop,
	);
}




### Export methods

sub version {
	my $self = shift;
	if (@_ and $_[0] ne '3') {
		carp "sorry, only GFF version 3 is currently supported!";
	}
	return 3;
}

sub gff_string {
	# exports version 3 GFF, sorry, no mechanism for others....
	my $self = shift;
	my $recurse = shift || 0;
	my $childIDs = shift || undef;
	
	# map all parent-child relationships
	if ($recurse and not defined $childIDs) {
		# we have a top level SeqFeature and we need to go spelunking for IDs
		$childIDs = {};
		foreach my $f ($self->get_SeqFeatures) {
			$f->_spelunk($self->primary_id, $childIDs);
		}
	}
	
	# basics
	my $string = join("\t", (
		$self->seq_id || '.',
		$self->source_tag || '.',
		$self->primary_tag,
		$self->start || 1,
		$self->end || 1,
		defined $self->score ? $self->score : '.',
		$self->strand > 0 ? '+' : $self->strand < 0 ? '-' : '.',
		defined $self->phase ? $self->phase : '.',
	) );
	
	# group attributes
	my $attributes = "Name=" . $self->_encode( $self->display_name ) . 
		'; ID=' . $self->_encode( $self->primary_id );
	if ($childIDs) {
		if (exists $childIDs->{$self->primary_id}) {
			$attributes .= '; Parent=' . join(',', 
				values %{ $childIDs->{$self->primary_id} });
			delete $childIDs->{$self->primary_id};
		}
	}
	foreach my $tag ($self->all_tags) {
		next if $tag eq 'Name';
		next if $tag eq 'ID';
		next if $tag eq 'Parent';
		my $value = $self->get_tag_values($tag);
		if (ref($value) eq 'ARRAY') {
			$value = join(",", map { $self->_encode($_) } @$value);
		}
		else {
			$value = $self->_encode($value);
		}
		$attributes .= "; $tag=$value";
	}
	$string .= "\t$attributes\n";
	
	# recurse
	if ($recurse) {
		foreach my $f ($self->get_SeqFeatures) {
			$string .= $f->gff_string(1, $childIDs);
		}
	}
	return $string;
}

sub _encode {
	my ($self, $value) = @_;
	$value =~ s/([\t\n\r%&\=;, ])/sprintf("%%%X",ord($1))/ge;
	return $value;
}

sub _spelunk {
	my ($self, $parentID, $childIDs);
	$childIDs->{$self->primary_id}{$parentID} += 1;
	foreach my $f ($self->get_SeqFeatures) {
		$f->_spelunk($self->primary_id, $childIDs);
	}
}

sub bed_string {
	my $self = shift;
	my $string = join("\t", (
		$self->seq_id || '.',
		($self->start || 1) - 1,
		$self->end || 1,
		$self->display_name || $self->primary_id,
		defined $self->score ? $self->score : 0,
		$self->strand >= 0 ? '+' : '-',
	) );
	return "$string\n";
}


__END__

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
