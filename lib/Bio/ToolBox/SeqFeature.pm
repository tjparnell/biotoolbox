package Bio::ToolBox::SeqFeature;

use warnings;
use strict;
use Carp qw(carp cluck croak confess);
use constant {
	SEQID => 0,
	START => 1,
	STOP  => 2,
	STRND => 3,
	NAME  => 4,
	ID    => 5,
	TYPE  => 6,
	SRC   => 7,
	SCORE => 8,
	PHASE => 9,
	ATTRB => 10,
	SUBF  => 11,
};

our $VERSION = '1.70';


#### Aliases ####
# to maintain compatibility with Bio::SeqFeature::Lite and Bio::SeqFeatureI we
# put in plenty of aliases to some of the methods
*stop                = \&end;
*name                = \&display_name;
*id                  = \&primary_id;
*method              = \&primary_tag;
*source              = \&source_tag;
*add_segment         = \&add_SeqFeature;
*get_all_SeqFeatures = *segments = \&get_SeqFeatures;
*each_tag_value      = \&get_tag_values;
*get_all_tags        = \&all_tags;
*gff3_string         = \&gff_string;

# avoid once warning errors
*stop                 if 0;
*name                 if 0;
*id                   if 0;
*method               if 0;
*source               if 0;
*add_segment          if 0;
*get_all_SeqFeatures  if 0;
*segments             if 0;
*each_tag_value       if 0;
*get_all_tags         if 0;
*gff3_string          if 0;

#### METHODS ####
sub new {
	my $class = shift;
	if ( ref($class) ) {
		$class = ref($class);
	}
	my %args = @_;
	my $self = bless [], $class;

	# bless ourselves early to take advantage of some methods
	# but otherwise do what we can simply and quickly

	# primary options
	$self->[SEQID] =
		   $args{-seq_id}
		|| $args{-seqid}
		|| $args{'-ref'}
		|| $args{chrom}
		|| undef;
	$self->[START] = int $args{-start} || undef;
	$self->[STOP]  = int $args{-end}   || $args{-stop} || undef;
	$self->strand( $args{-strand} ) if exists $args{-strand};
	$self->[NAME] = $args{-display_name} || $args{-name} || undef;
	$self->[ID]   = $args{-primary_id}   || $args{-id}   || undef;

	# check orientation
	if (    defined $self->[START]
		and defined $self->[STOP]
		and $self->[START] > $self->[STOP] )
	{
		# flip the coordinates around
		( $self->[START], $self->[STOP] ) = ( $self->[STOP], $self->[START] );
		$self->[STRND] *= -1;
	}

	# additional options
	$args{-type} ||= $args{-primary_tag} || undef;
	if ( defined $args{-type} ) {
		$self->type( $args{-type} );
	}
	$args{-source} ||= $args{-source_tag} || undef;
	if ( defined $args{-source} ) {
		$self->[SRC] = $args{-source};
	}
	if ( exists $args{-score} ) {
		$self->[SCORE] = $args{-score};
	}
	if ( exists $args{-phase} ) {
		$self->phase( $args{-phase} );
	}
	if ( exists $args{-attributes} or exists $args{-tags} ) {
		$args{-attributes} ||= $args{-tags};
		if ( ref( $args{-attributes} ) eq 'HASH' ) {
			$self->[ATTRB] = $args{-attributes};
		}
	}
	if ( exists $args{-segments} ) {
		$self->[SUBF] = [];
		foreach my $s ( @{ $args{-segments} } ) {
			unless ( ref($s) eq $class ) {
				croak "segments should be passed as $class objects!";
			}
			push @{ $self->[SUBF] }, $s;
		}
	}

	return $self;
}

sub seq_id {
	my $self = shift;
	if (@_) {
		$self->[SEQID] = $_[0];
	}
	return defined $self->[SEQID] ? $self->[SEQID] : undef;
}

sub start {
	my $self = shift;
	if (@_) {
		$self->[START] = int $_[0];
	}
	return defined $self->[START] ? $self->[START] : undef;
}

sub end {
	my $self = shift;
	if (@_) {
		$self->[STOP] = int $_[0];
	}
	return defined $self->[STOP] ? $self->[STOP] : undef;
}

sub strand {
	my $self = shift;
	if (@_) {
		my $str = $_[0];
		if ( $str eq '+' ) {
			$self->[STRND] = 1;
		}
		elsif ( $str eq '-' ) {
			$self->[STRND] = -1;
		}
		elsif ( $str eq '.' ) {
			$self->[STRND] = 0;
		}
		elsif ( $str =~ /^[\+\-]?1$/ ) {
			$self->[STRND] = $str;
		}
		elsif ( $str eq '0' ) {
			$self->[STRND] = 0;
		}
	}
	return defined $self->[STRND] ? $self->[STRND] : 0;
}

sub display_name {
	my $self = shift;
	if (@_) {
		$self->[NAME] = $_[0];
	}
	return defined $self->[NAME] ? $self->[NAME] : $self->primary_id;
}

sub primary_id {
	my $self = shift;
	if (@_) {
		$self->[ID] = $_[0];
	}
	elsif ( not defined $self->[ID] ) {

		# automatically assign a new ID
		$self->[ID] =
			sprintf( "%s:%d-%d", $self->[SEQID], $self->[START], $self->[STOP] );
	}
	return $self->[ID];
}

sub primary_tag {
	my $self = shift;
	if (@_) {
		$self->[TYPE] = $_[0];
	}
	return defined $self->[TYPE] ? $self->[TYPE] : 'feature';
}

sub source_tag {
	my $self = shift;
	if (@_) {
		$self->[SRC] = $_[0];
	}
	return defined $self->[SRC] ? $self->[SRC] : '';
}

sub type {
	my $self = shift;
	if (@_) {
		my $type = $_[0];
		if ( $type =~ /:/ ) {
			my ( $t, $s ) = split /:/, $type;
			$self->[TYPE] = $t;
			$self->[SRC]  = $s;
		}
		else {
			$self->[TYPE] = $type;
		}
	}
	$self->[TYPE] ||= undef;
	$self->[SRC]  ||= undef;
	if ( defined $self->[SRC] ) {
		return join( ':', $self->[TYPE], $self->[SRC] );
	}
	return $self->[TYPE];
}

sub score {
	my $self = shift;
	if (@_) {
		$self->[SCORE] = $_[0];
	}
	return defined $self->[SCORE] ? $self->[SCORE] : undef;
}

sub phase {
	my $self = shift;
	if (@_) {
		my $p = $_[0];
		unless ( $p =~ /^[012\.]$/ ) {
			cluck "phase must be 0, 1, 2 or .!";
		}
		$self->[PHASE] = $p;
	}
	return defined $self->[PHASE] ? $self->[PHASE] : undef;
}

sub duplicate {
	my $self = shift;
	my $n    = $self->new(
		-seq_id       => $self->[SEQID],
		-start        => $self->[START],
		-end          => $self->[STOP],
		-strand       => $self->strand,
		-type         => $self->primary_tag,
		-source       => $self->source_tag,
		-display_name => $self->name
	);
	$n->score( $self->[SCORE] ) if defined $self->[SCORE];
	$n->phase( $self->[PHASE] ) if defined $self->[PHASE];
	return $n;
}

sub add_SeqFeature {
	my $self = shift;
	$self->[SUBF] ||= [];
	my $count = 0;
	foreach my $s (@_) {
		if ( ref($s) eq ref($self) ) {
			$s->seq_id( $self->seq_id )         unless ( $s->seq_id );
			$s->strand( $self->strand )         unless ( $s->strand );
			$s->source_tag( $self->source_tag ) unless ( $s->source_tag );
			push @{ $self->[SUBF] }, $s;
			$count++;
		}
		else {
			cluck "please use SeqFeature objects when adding sub features!";
		}
	}
	return $count;
}

sub get_SeqFeatures {
	my $self = shift;
	return unless $self->[SUBF];
	my @children;
	foreach ( @{ $self->[SUBF] } ) {
		push @children, $_;
	}
	return wantarray ? @children : \@children;
}

sub delete_SeqFeature {
	my $self = shift;
	return unless $self->[SUBF];
	my $id = shift;
	my $d;
	my $i = 0;
	foreach ( @{ $self->[SUBF] } ) {
		if ( $_->[ID] eq $id ) {
			$d = $i;
			last;
		}
		$i++;
	}
	if ( defined $d ) {
		splice( @{ $self->[SUBF] }, $d, 1 );
		return 1;
	}
	return;
}

sub add_tag_value {
	my ( $self, $key, $value ) = @_;
	return unless ( $key and $value );
	$self->[ATTRB] ||= {};
	if ( exists $self->[ATTRB]->{$key} ) {
		if ( ref( $self->[ATTRB]->{$key} ) eq 'ARRAY' ) {
			push @{ $self->[ATTRB]->{$key} }, $value;
		}
		else {
			my $current = $self->[ATTRB]->{$key};
			$self->[ATTRB]->{$key} = [ $current, $value ];
		}
	}
	else {
		$self->[ATTRB]->{$key} = $value;
	}
}

sub has_tag {
	my ( $self, $key ) = @_;
	return unless $self->[ATTRB];
	return exists $self->[ATTRB]->{$key};
}

sub get_tag_values {
	my ( $self, $key ) = @_;
	return unless $self->[ATTRB];
	if ( exists $self->[ATTRB]->{$key} ) {
		if ( ref( $self->[ATTRB]->{$key} ) eq 'ARRAY' ) {
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
	return unless $self->[ATTRB];
	return wantarray ? %{ $self->[ATTRB] } : $self->[ATTRB];
}

sub all_tags {
	my $self = shift;
	return unless $self->[ATTRB];
	my @k = keys %{ $self->[ATTRB] };
	return wantarray ? @k : \@k;
}

sub remove_tag {
	my ( $self, $key ) = @_;
	return unless $self->[ATTRB];
	if ( exists $self->[ATTRB]->{$key} ) {
		delete $self->[ATTRB]->{$key};
	}
}

### Range Methods
# Borrowed from Bio::RangeI, but does not do strong/weak strand checks

sub length {
	my $self = shift;
	return $self->end - $self->start + 1;
}

sub overlaps {
	my ( $self, $other ) = @_;
	return unless ( $other and ref($other) );
	return unless ( $self->seq_id eq $other->seq_id );
	return not( $self->start > $other->end
		or $self->end < $other->start );
}

sub contains {
	my ( $self, $other ) = @_;
	return unless ( $other and ref($other) );
	return unless ( $self->seq_id eq $other->seq_id );
	return ( $other->start >= $self->start and $other->end <= $self->end );
}

sub equals {
	my ( $self, $other ) = @_;
	return unless ( $other and ref($other) );
	return unless ( $self->seq_id eq $other->seq_id );
	return ( $other->start == $self->start and $other->end == $self->end );
}

sub intersection {
	my ( $self, $other ) = @_;
	return unless ( $other and ref($other) );
	return unless ( $self->seq_id eq $other->seq_id );
	my ( $start, $stop );
	if ( $self->start >= $other->start ) {
		$start = $self->start;
	}
	else {
		$start = $other->start;
	}
	if ( $self->end <= $other->end ) {
		$stop = $self->end;
	}
	else {
		$stop = $other->end;
	}
	return if $start > $stop;
	return $self->new(
		-seq_id => $self->seq_id,
		-start  => $start,
		-end    => $stop,
	);
}

sub union {
	my ( $self, $other ) = @_;
	return unless ( $other and ref($other) );
	return unless ( $self->seq_id eq $other->seq_id );
	my ( $start, $stop );
	if ( $self->start <= $other->start ) {
		$start = $self->start;
	}
	else {
		$start = $other->start;
	}
	if ( $self->end >= $other->end ) {
		$stop = $self->end;
	}
	else {
		$stop = $other->end;
	}
	return if $start > $stop;
	return $self->new(
		-seq_id => $self->seq_id,
		-start  => $start,
		-end    => $stop,
	);
}

sub subtract {
	my ( $self, $other ) = @_;
	return unless ( $other and ref($other) );
	return unless ( $self->seq_id eq $other->seq_id );
	return if $self->equals($other);
	my @pieces;
	if ( $self->overlaps($other) ) {
		my $int = $self->intersection($other);
		if ( $self->start < $int->start ) {
			push @pieces,
				$self->new(
					-seq_id => $self->seq_id,
					-start  => $self->start,
					-end    => $int->start - 1,
				);
		}
		if ( $self->end > $int->end ) {
			push @pieces,
				$self->new(
					-seq_id => $self->seq_id,
					-start  => $int->end + 1,
					-end    => $self->end,
				);
		}
	}
	else {
		push @pieces,
			$self->new(
				-seq_id => $self->seq_id,
				-start  => $self->start,
				-end    => $self->end,
			);
	}
	return wantarray ? @pieces : \@pieces;
}

### Export methods

sub version {
	my $self = shift;
	if ( @_ and $_[0] ne '3' ) {
		carp "sorry, only GFF version 3 is currently supported!";
	}
	return 3;
}

sub gff_string {

	# exports version 3 GFF, sorry, no mechanism for others....
	my $self     = shift;
	my $recurse  = shift || 0;
	my $childIDs = shift || undef;

	# skip myself if we've already printed myself
	return if ( $recurse
		and $childIDs
		and not exists $childIDs->{ $self->primary_id } );

	# map all parent-child relationships
	if ( $recurse and not defined $childIDs ) {

		# we have a top level SeqFeature and we need to go spelunking for IDs
		$childIDs = {};
		foreach my $f ( $self->get_SeqFeatures ) {
			$f->_spelunk( $self->primary_id, $childIDs );
		}
	}

	# basics
	my $string = join(
		"\t",
		(
			$self->seq_id     || '.',
			$self->source_tag || '.',
			$self->primary_tag,
			$self->start || 1,
			$self->end   || 1,
			defined $self->score ? $self->score : '.',
			$self->strand > 0    ? '+' : $self->strand < 0 ? '-' : '.',
			defined $self->phase ? $self->phase : '.',
		)
	);

	# group attributes
	my $attributes;
	my $name = $self->display_name;
	if ($name) {

		# names are optional and may not exist
		$attributes = sprintf "Name=%s;ID=%s", $self->_encode($name),
			$self->_encode( $self->primary_id );
	}
	else {
		# IDs always exist
		$attributes = sprintf "ID=%s", $self->_encode( $self->primary_id );
	}
	if ($childIDs) {
		if ( exists $childIDs->{ $self->primary_id } ) {
			$attributes .= ';Parent='
				. join( ',',
					sort { $a cmp $b }
					map  { $self->_encode($_) }
					keys %{ $childIDs->{ $self->primary_id } } );
			delete $childIDs->{ $self->primary_id };
		}
	}
	foreach my $tag ( $self->all_tags ) {
		next if $tag eq 'Name';
		next if $tag eq 'ID';
		next if $tag eq 'Parent';
		my $value = $self->get_tag_values($tag);
		if ( ref($value) eq 'ARRAY' ) {
			$value = join( ",", map { $self->_encode($_) } @$value );
		}
		else {
			$value = $self->_encode($value);
		}
		$attributes .= ";$tag=$value";
	}
	$string .= "\t$attributes\n";

	# recurse
	if ($recurse) {
		foreach my $f ( $self->get_SeqFeatures ) {
			my $s = $f->gff_string( 1, $childIDs );
			$string .= $s if $s;
		}
	}
	return $string;
}

sub _encode {
	my ( $self, $value ) = @_;
	$value =~ s/([\t\n\r%&\=;, ])/sprintf("%%%X",ord($1))/ge;
	return $value;
}

sub _spelunk {
	my ( $self, $parentID, $childIDs ) = @_;
	$childIDs->{ $self->primary_id }{$parentID} += 1;
	foreach my $f ( $self->get_SeqFeatures ) {
		$f->_spelunk( $self->primary_id, $childIDs );
	}
}

sub bed_string {
	my $self   = shift;
	my $string = join(
		"\t",
		(
			$self->seq_id || '.',
			( $self->start || 1 ) - 1,
			$self->end          || 1,
			$self->display_name || $self->primary_id,
			defined $self->score ? $self->score : 0,
			$self->strand >= 0   ? '+'          : '-',
		)
	);
	return "$string\n";
}

1;

__END__

=head1 NAME

Bio::ToolBox::SeqFeature - Fast, simple SeqFeature implementation

=head1 SYNOPSIS

   # create a transcript
   my $transcript = Bio::ToolBox::SeqFeature->new(
        -seq_id     => chr1,
        -start      => 1001,
        -stop       => 1500,
        -strand     => '+',
   );
   $seqf->primary_tag('mRNA'); # set parameters individually
   
   # create an exon
   my $exon = Bio::ToolBox::SeqFeature->new(
        -start      => 1001,
        -end        => 1200,
        -type       => 'exon',
   );
   
   # associate exon with transcript
   $transcript->add_SeqFeature($exon); 
   my $exon_strand = $exon->strand; # inherits from parent
   
   # walk through subfeatures
   foreach my $f ($transcript->get_all_SeqFeatures) {
   	  printf "%s is a %s\n", $f->display_name, $f->type;
   }
   
   # add attribute
   $transcript->add_tag_value('Status', $status);
   
   # get attribute
   my $value = $transcript->get_tag_values($key);

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
implements the L<Bio::SeqFeatureI> API, and could be substituted for other 
implementations, such L<Bio::SeqFeature::Lite> and L<Bio::SeqFeature::Generic>. 
Unlike the others, however, it inherits no classes or methods and uses an 
unorthodox blessed array to store feature attributes, decreasing memory 
requirements and overhead complexity. 

=head1 METHODS

Refer to the L<Bio::SeqFeatureI> documentation for general implementation and 
ideas, which this module tries to implement.

=head2 Creating new SeqFeature objects

=over 4

=item new

Generate a new SeqFeature object. Pass an array of key =E<gt> value pairs 
with the feature parameters. At a minimum, location information (C<seq_id>, 
C<start>, and C<stop>) should be provided, with name and ID recommended. 
Most of the accession methods may be used as key tags to the C<new> method. 
The following attribute keys are accepted.

=over 4

=item -seq_id

=item -start

=item -end

=item -stop

=item -strand

=item -name

=item -display_name

=item -id

=item -primary_id

=item -type

=item -source

=item -source_tag

=item -primary_tag

=item -score

=item -phase

=item -attributes

Provide an anonymous array of key value pairs representing attribute 
keys and their values.

=item -segments

Provide an anonymous array of SeqFeature objects to add as child 
objects. 

=back

=item duplicate

This method will duplicate a SeqFeature object as a new object, with 
limitations. The C<primary_id> tag, attribute key/values, and subfeatures 
are B<not> duplicated. This basically limits the object to coordinates, 
C<primary_tag>, C<source_tag>, C<display_name>, C<score>, and C<phase>. 
Remaining attributes should be explicitly set by the user.

=back

=head2 Accession methods

These are methods to set and/or retrieve attribute values. Pass a 
single value to set the attribute. The attribute value is always 
returned.

=over 4

=item seq_id

The name of the chromosome or reference sequence. If you are generating 
a new child object, e.g. exon, then seq_id does not need to be provided. 
In these cases, the parent's seq_id is inherited.

=item start

The start coordinate of the feature. SeqFeature objects use the 1-base 
numbering system, following BioPerl convention. Start coordinates are 
always less than the stop coordinate. 

=item end

=item stop

The end coordinate of the feature. Stop coordinates are always greater 
than the start coordinate.

=item strand

The strand that the feature is on. The default value is always unstranded, 
or 0. Any of the following may be supplied: "C<1 0 -1 + . ->". Numeric 
integers 1, 0, and -1 are always returned.

=item source

=item source_tag

A text string representing the source of the feature. This corresponds to 
the second field in a GFF file. The source tag is optional. If not supplied, 
the source tag can be inherited from a parent feature if present. 

=item primary_tag

The type of feature. These usually follow Sequence Ontology terms, but are 
not required. The default value is the generic term "region". Examples 
include gene, mRNA, transcript, exon, CDS, etc. This corresponds to the 
third field in a GFF file.

=item type

A shortcut method which can represent either "primary_tag:source_tag" or, 
if no source_tag is defined, simply "primary_tag".

=item name

=item display_name

A text string representing the name of the feature. The name is not 
required to be a unique value, but generally is. The default name, if 
none is provided, is the C<primary_id>.

=item id

=item primary_id

A text string representing a unique identifier of the feature. If not 
explicitly defined, a unique ID is automatically generated.

=item score

A numeric (integer or floating) value representing the feature. 

=item phase 

An integer (0,1,2) representing the coding frame. Only required for 
CDS features.

=back

=head2 Special Attributes

Special attributes are key value pairs that do not fall under the above 
conventions. It is possible to have more than one value assigned to a 
given key. In a GFF file, this corresponds to the attributes in the 9th 
field, with the exception of special reserved attributes such as C<Name>, 
C<ID>, and C<Parent>.

=over 4

=item add_tag_value

    $seqf->add_tag_value($key, $value);

Sets the special attribute C<$key> to C<$value>. If you have more than one 
value, C<$value> should be an anonymous array of text values. Following 
GFF convention, C<$key> should not comprise of special characters, including 
";,= ".

=item all_tags

=item get_all_tags

Returns an array of all attribute keys.

=item has_tag($key)

Boolean method whether the SeqFeature object contains the attribute.

=item get_tag_values($key)

=item each_tag_value($key)

Returns the value for attribute C<$key>. If multiple values are present, 
it may return an array or array reference.

=item attributes

Returns an array or reference to the key value hash;

=item remove_tag($key)

Deletes the indicated attribute.

=back

=head2 Subfeature Hierarchy

Following Sequence Ontology and GFF conventions, SeqFeatures can have 
subfeatures (children) representing a hierarchical structure, for example 
genes beget transcripts which beget exons. 

Child SeqFeature objects may have more than one parent, for example, shared
exons between alternate transcripts. In which case, only one exon
SeqFeature object is generated, but is added to both transcript objects. 

=over 4

=item add_SeqFeature

=item add_segment

Pass one or more SeqFeature objects to be associated as children.

=item get_SeqFeatures

=item get_all_SeqFeatures

=item segments

Returns an array of all sub SeqFeature objects.

=item delete_SeqFeature($id)

This method will delete a specific child SeqFeature. Pass the method the 
C<primary_id> of the SeqFeature. All SeqFeatures will have a C<primary_id>. 

=back

=head2 Range Methods

These are range methods for comparing one SeqFeature object to another.
They are analogous to L<Bio::RangeI> methods.

They currently do not support strand checks or strand options.

=over 4

=item length

Returns the length of the SeqFeature object.

=item overlaps

	my $overlap_check = $seqf1->overlaps($other_seqf);

Returns a boolean value whether a second SeqFeature object overlaps 
with the self object.

=item contains

	my $check = $seqf1->contains($seqf2);

Returns a boolean value whether the self object completely 
contains a second SeqFeature object.

=item equals

	my $check = $seqf1->equals($seqf2);

Returns a boolean value whether the self object coordinates are 
equivalent to a second SeqFeature object.

=item intersection

	my $check = $seqf1->intersection($seqf2);

Returns a new SeqFeature object representing the intersection or 
overlap area between the self object and a second SeqFeature 
object.

=item union($other)

	my $new_seqf = $seqf1->union($seqf2);

Returns a new SeqFeature object representing the merged interval 
between the self and a second SeqFeature object.

=item subtract

	my $new_seqf = $seqf1->subtract($seqf2);

Returns a new SeqFeature object representing the interval of the 
self object after subtracting a second SeqFeature object. 

=back

=head2 Export Strings

These methods export the SeqFeature object as a text string in the 
specified format. New line characters are included.

=over 4

=item gff_string($recurse)

Exports the SeqFeature object as a GFF3 formatted string. Pass a 
boolean value if you wish to recurse through the hierarchy and 
print subfeatures as a multi-line string. Child-Parent ID attributes 
are smartly handled, including multiple parentage. 

See L<Bio::ToolBox::GeneTools> for methods for printing GTF transcript models.

=item bed_string

Exports the SeqFeature object as a simple BED6 formatted string. 
See L<Bio::ToolBox::GeneTools> for methods for printing BED12 transcript models.

=back

=head1 LIMITATIONS

Because of their underlying array structure, Bio::ToolBox::SeqFeature objects 
should generally not be used as a base class (unless you know the ramifications 
of doing so). The following Bio classes and Interfaces are similar and their API 
was used as a model. However, in most cases they are not likely to work with this 
module because of object structure incompatibility, although this has not been 
explicitly tested. 

=over 4

=item Bio::AnnotationI

=item Bio::LocationI

=item Bio::RangeI

=item Bio::SeqI

=item Bio::Tools::GFF

=item Bio::DB::SeqFeature::Store

=item Bio::Graphics

=back

=head1 SEE ALSO

L<Bio::SeqFeatureI>, L<Bio::SeqFeature::Lite>, L<Bio::Seqfeature::Generic>, 
L<Bio::RangeI>, L<Bio::ToolBox::parser::gff>, L<Bio::ToolBox::parser::ucsc>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
