package Bio::ToolBox::Data::Feature;

use warnings;
use strict;
use English qw(-no_match_vars);
use Carp    qw(carp cluck croak confess);
use Module::Load;
use Bio::ToolBox::db_helper qw(
	get_db_feature
	get_segment_score
	calculate_score
	get_genomic_sequence
	low_level_bam_fetch
);
use Bio::ToolBox::db_helper::constants;

our $VERSION = '1.70';

my $GENETOOL_LOADED = 0;

### Initialization

sub new {

	# this should ONLY be called from Bio::ToolBox::Data* iterators
	my $class = shift;
	my %self  = @_;

	# we trust that new is called properly with data and index values
	return bless \%self, $class;
}

### Set and retrieve values

sub data {
	return shift->{data};
}

sub column_name {
	my ( $self, $column ) = @_;
	return unless defined $column;
	return $self->{data}->name($column);
}

sub feature_type {
	my $self = shift;
	carp 'ERROR: feature_type is a read only method' if @_;
	return $self->{data}->feature_type;
}

sub row_index {
	my $self = shift;
	carp 'ERROR: row_index is a read only method' if @_;
	return $self->{'index'};
}

sub line_number {
	my $self = shift;
	carp 'ERROR: line_number is a read only method' if @_;
	if ( exists $self->{data}->{line_count} ) {
		return $self->{data}->{line_count};
	}
	elsif ( exists $self->{data}->{header_line_count} ) {
		return $self->{data}->{header_line_count} + $self->row_index;
	}
	else {
		return;
	}
}

sub row_values {
	my $self = shift;
	carp 'ERROR: row_values is a read only method' if @_;
	my $row  = $self->{'index'};
	my @data = @{ $self->{data}->{data_table}->[$row] }
		[ 1 .. $#{ $self->{data}->{data_table}->[$row] } ];
	return wantarray ? @data : \@data;
}

sub value {
	my ( $self, $column, $value ) = @_;
	return unless ( defined $column and $column =~ /^\d+/ );
	my $row = $self->{'index'};

	# set a value
	if ( defined $value ) {
		$self->{data}->{data_table}->[$row][$column] = $value;
	}

	my $v = $self->{data}->{data_table}->[$row][$column];
	return length($v) ? $v : '.';    # internal null value, inherited from GFF definition
}

sub seq_id {
	my $self = shift;
	if ( $_[0] ) {
		my $i = $self->{data}->chromo_column;

		# update only if we have an actual column
		if ($i) {
			$self->value( $i, $_[0] );
			$self->{seqid} = $_[0];
		}
		elsif ( exists $self->{feature} ) {
			carp 'ERROR: Unable to update seq_id for parsed SeqFeature objects';
		}
		else {
			carp 'ERROR: No Chromosome column to update!';
		}
	}

	# return cached value
	return $self->{seqid} if exists $self->{seqid};

	# get seq_id value
	my $c;
	my $i = $self->{data}->chromo_column;
	my $j = $self->{data}->coord_column;
	if ($i) {
		$c = $self->value($i);
	}
	elsif ( exists $self->{feature} ) {
		$c = $self->{feature}->seq_id;
	}
	elsif ($j) {
		$self->_extract_coordinate_string($j);
		return $self->{seqid};
	}
	elsif ( $self->feature_type eq 'named' ) {
		my $f = $self->seqfeature;
		if ($f) {
			$c = $f->seq_id;
		}
	}
	$self->{seqid} = $c if defined $c;
	return $c;
}

sub start {
	my $self = shift;
	if ( $_[0] ) {
		my $i = $self->{data}->start_column;

		# update only if we have an actual column
		my $d = $_[0] =~ /^\d+$/ ? 1 : 0;    # looks like an integer
		if ( $i and $d ) {
			if ( $self->{data}->interbase ) {

				# compensate for 0-based, assuming we're always working with 1-based
				$self->value( $i, $_[0] - 1 );
				$self->{start} = $_[0];
			}
			else {
				$self->value( $i, $_[0] );
				$self->{start} = $_[0];
			}
		}
		elsif ( not $d ) {
			carp 'ERROR: Start coordinate value is not an integer';
		}
		elsif ( exists $self->{feature} ) {
			carp 'ERROR: Unable to update Start coordinate for parsed SeqFeature objects';
		}
		else {
			carp 'ERROR: No Start coordinate column to update!';
		}
	}

	# return cached value
	return $self->{start} if exists $self->{start};

	# get start value
	my $s;
	my $i = $self->{data}->start_column;
	my $j = $self->{data}->coord_column;
	if ($i) {

		# collect from table
		$s = int $self->value($i);
		if ( $self->{data}->interbase ) {
			$s += 1;    # compensate for 0-based index
		}
	}
	elsif ( exists $self->{feature} ) {
		$s = $self->{feature}->start;
	}
	elsif ($j) {
		$self->_extract_coordinate_string($j);
		return $self->{start};
	}
	elsif ( $self->feature_type eq 'named' ) {
		my $f = $self->seqfeature;
		if ($f) {
			$s = $f->start;
		}
	}
	$self->{start} = $s if defined $s;
	return $s;
}

*stop = \&end;
*stop if 0;    # avoid once warning

sub end {
	my $self = shift;
	if ( $_[0] ) {
		my $i = $self->{data}->stop_column;

		# update only if we have an actual column
		my $d = $_[0] =~ /^\d+$/ ? 1 : 0;
		if ( $i and $d ) {
			$self->value( $i, $_[0] );
			$self->{end} = $_[0];
		}
		elsif ( not $d ) {
			carp 'ERROR: End coordinate value is not an integer';
		}
		elsif ( exists $self->{feature} ) {
			carp 'ERROR: Unable to update End coordinate for parsed SeqFeature objects';
		}
		else {
			carp 'ERROR: No End coordinate column to update!';
		}
	}

	# return cached value
	return $self->{end} if exists $self->{end};

	# get end value
	my $e;
	my $i = $self->{data}->stop_column;
	my $j = $self->{data}->coord_column;
	if ($i) {
		$e = int $self->value($i);
	}
	elsif ( exists $self->{feature} ) {
		$e = $self->{feature}->end;
	}
	elsif ($j) {
		$self->_extract_coordinate_string($j);
		return $self->{end};
	}
	elsif ( $self->feature_type eq 'named' ) {
		my $f = $self->seqfeature;
		if ($f) {
			$e = $f->end;
		}
	}
	$self->{end} = $e if defined $e;
	return $e;
}

sub strand {
	my $self = shift;
	if ( $_[0] ) {

		# update only if we have an actual column
		my $i = $self->{data}->strand_column;
		if ($i) {
			$self->value( $i, $_[0] );
			$self->{strand} = $self->_strand( $_[0] );
		}
		elsif ( exists $self->{feature} ) {
			carp 'ERROR: Unable to update Strand for parsed SeqFeature objects';
		}
		else {
			carp 'ERROR: No Strand column to update!';
		}
	}

	# return cached value
	return $self->{strand} if exists $self->{strand};

	# get strand value
	my $s = 0;
	if ( my $i = $self->{data}->strand_column ) {
		$s = $self->_strand( $self->value($i) );
	}
	elsif ( exists $self->{feature} ) {
		$s = $self->{feature}->strand;
	}
	elsif ( $self->feature_type eq 'named' ) {
		my $f = $self->seqfeature;
		if ($f) {
			$s = $f->strand;
		}
	}
	$self->{strand} = $s;
	return $s;
}

sub _strand {
	my $self = shift;
	my $str  = shift;
	if ( $str eq '1' or $str eq '-1' or $str eq '0' ) {
		return $str;
	}
	elsif ( $str eq '+' ) {
		return 1;
	}
	elsif ( $str eq '-' ) {
		return -1;
	}
	elsif ( $str eq '.' ) {
		return 0;
	}
	elsif ( $str eq '2' ) {

		# R packages use 1 for + and 2 for -
		return -1;
	}
	elsif ( $str eq '*' ) {

		# R packages use * for .
		return 0;
	}
	else {
		return 0;
	}
}

sub _extract_coordinate_string {
	my ( $self, $i ) = @_;
	if ( $self->value($i) =~ /^ ([\w\.\-]+) : (\d+) (?: \.\. | \-) (\d+) $/x ) {

		# chromosome:start-end or chromosome:start..end
		$self->{seqid} = $1 unless exists $self->{seqid};
		$self->{start} = $2 unless exists $self->{start};
		$self->{end}   = $3 unless exists $self->{end};
	}
	elsif ( $self->value($i) =~ /^([\w\.\-]+) : (\d+) $/x ) {

		# chromosome:start
		$self->{seqid} = $1 unless exists $self->{seqid};
		$self->{start} = $2 unless exists $self->{start};
	}
}

sub peak {
	my $self = shift;
	if ( $self->{data}->format eq 'narrowPeak' ) {
		if ( exists $self->{feature} and $self->{feature}->has_tag('peak') ) {
			return (
				int( $self->{feature}->get_tag_values('peak') ) +
					$self->{feature}->start );
		}
		else {
			# hard coded start and peak columns
			return ( int $self->value(2) + int $self->value(10) + 1 );
		}
	}
	else {
		return $self->midpoint;
	}
}

sub midpoint {
	my $self = shift;
	my $s    = $self->start;
	my $e    = $self->end;
	if ( $s and $e ) {
		return int( ( $s + $e ) / 2 );
	}
	else {
		return undef;
	}
}

*name = \&display_name;
*name if 0;    # avoid once warning

sub display_name {
	my $self = shift;
	if ( $_[0] ) {

		# update only if we have an actual column
		my $i = $self->{data}->name_column;
		if ($i) {
			return $self->value( $i, $_[0] );
		}
		elsif ( exists $self->{feature} ) {
			carp 'ERROR: Unable to update display_name for parsed SeqFeature objects';
		}
		else {
			carp 'ERROR: No Name column to update!';
		}
	}

	# seqfeature
	if ( exists $self->{feature} ) {
		return $self->{feature}->display_name;
	}

	# collect from table
	my $i = $self->{data}->name_column;
	if ($i) {
		return $self->value($i);
	}
	elsif ( my $att = $self->gff_attributes ) {
		return
			   $att->{Name}
			|| $att->{ID}
			|| $att->{transcript_name}
			|| $att->{gene_name}
			|| undef;
	}
}

sub coordinate {
	my $self = shift;
	carp 'ERROR: name is a read only method' if @_;
	my $c = $self->seq_id;
	my $s = $self->start;
	my $e = $self->end;
	if ( $c and $s and $e ) {
		return sprintf "%s:%d-%d", $c, $s, $e;
	}
	elsif ( $c and $s ) {
		return sprintf "%s:%d", $c, $s;
	}
	else {
		return undef;
	}
}

sub type {
	my $self = shift;
	if ( $_[0] ) {

		# update only if we have an actual column
		my $i = $self->{data}->type_column;
		if ($i) {
			return $self->value( $i, $_[0] );
		}
		elsif ( exists $self->{feature} ) {
			carp 'ERROR: Unable to update primary_tag for parsed SeqFeature objects';
		}
		else {
			carp 'ERROR: No Type column to update!';
		}
	}

	# collect from table
	my $i = $self->{data}->type_column;
	if ($i) {
		return $self->value($i);
	}

	# seqfeature
	if ( exists $self->{feature} ) {
		return $self->{feature}->primary_tag;
	}

	# general metadata
	if ( $self->{data}->feature ) {
		return $self->{data}->feature;
	}
	return undef;
}

*id = \&primary_id;
*id if 0;    # avoid once warning

sub primary_id {
	my $self = shift;
	carp 'ERROR: id is a read only method' if @_;
	my $i = $self->{data}->id_column;
	if ($i) {
		my $v = $self->value($i);
		if ( defined $v and $v ne '.' ) {
			return $v;
		}
	}
	return $self->{feature}->primary_id if exists $self->{feature};
	if ( my $att = $self->gff_attributes ) {
		return $att->{ID} || $att->{Name} || $att->{transcript_id};
	}
	return undef;
}

sub length {
	my $self = shift;
	carp 'ERROR: length is a read only method' if @_;
	if ( $self->{data}->vcf ) {

		# special case for vcf files, measure the length of the ALT allele
		return CORE::length( $self->value(5) );
	}
	my $s = $self->start;
	my $e = $self->end;
	if ( defined $s and defined $e ) {
		return $e - $s + 1;
	}
	elsif ( defined $s ) {
		return 1;
	}
	else {
		return undef;
	}
}

sub score {
	my $self = shift;
	my $c    = $self->{data}->score_column;
	return $c ? $self->value($c) : undef;
}

sub attributes {
	my $self = shift;
	return $self->gff_attributes if ( $self->{data}->gff );
	return $self->vcf_attributes if ( $self->{data}->vcf );
	return;
}

sub gff_attributes {
	my $self = shift;
	return unless ( $self->{data}->gff );
	return $self->{attributes} if ( exists $self->{attributes} );
	$self->{attributes} = {};
	foreach my $g ( split( /\s*;\s*/, $self->value(9) ) ) {
		my ( $tag, $value ) = split /\s+|=/, $g;
		next unless ( $tag and $value );

		# unescape URL encoded values, borrowed from Bio::DB::GFF
		$value =~ tr/+/ /;
		$value =~ s/%( [0-9a-fA-F]{2} )/chr hex($1)/xge;
		$self->{attributes}->{$tag} = $value;
	}
	return $self->{attributes};
}

sub vcf_attributes {
	my $self = shift;
	return unless ( $self->{data}->vcf );
	return $self->{attributes} if ( exists $self->{attributes} );
	$self->{attributes} = {};

	# INFO attributes
	my %info;
	if ( $self->{data}->name(8) eq 'INFO' ) {
		%info = map { $_->[0] => defined $_->[1] ? $_->[1] : undef }

			# some tags are simple and have no value, eg SOMATIC
			map { [ split(/=/) ] }
			split( /;/, $self->value(8) );
	}
	$self->{attributes}->{INFO} = \%info;
	$self->{attributes}->{8}    = \%info;

	# Sample attributes
	if ( $self->{data}->number_columns > 9 ) {
		my @formatKeys = split /:/, $self->value(9);
		foreach my $i ( 10 .. $self->{data}->last_column ) {
			my $name       = $self->{data}->name($i);
			my @sampleVals = split /:/, $self->value($i);
			my %sample     = map {
				$formatKeys[$_] => defined $sampleVals[$_] ? $sampleVals[$_] : undef
			} ( 0 .. $#formatKeys );
			$self->{attributes}->{$name} = \%sample;
			$self->{attributes}->{$i}    = \%sample;
		}
	}
	return $self->{attributes};
}

sub rewrite_attributes {
	my $self = shift;
	return $self->rewrite_gff_attributes if ( $self->{data}->gff );
	return $self->rewrite_vcf_attributes if ( $self->{data}->vcf );
	return;
}

sub rewrite_gff_attributes {
	my $self = shift;
	return unless ( $self->{data}->gff );
	return unless exists $self->{attributes};
	my @pairs;    # of key=value items
	if ( exists $self->{attributes}{ID} ) {

		# I assume this does not need to be escaped!
		push @pairs, 'ID=' . $self->{attributes}{ID};
	}
	if ( exists $self->{attributes}{Name} ) {
		my $name = $self->{attributes}{Name};
		$name =~ s/( [\t\n\r%&\=;,\ ] ) /sprintf("%%%X",ord($1))/xge;
		push @pairs, "Name=$name";
	}
	foreach my $key ( sort { $a cmp $b } keys %{ $self->{attributes} } ) {
		next if $key eq 'ID';
		next if $key eq 'Name';
		my $value = $self->{attributes}{$key};
		$key   =~ s/( [\t\n\r%&\=;,\ ] ) /sprintf("%%%X",ord($1))/xge;
		$value =~ s/( [\t\n\r%&\=;,\ ] ) /sprintf("%%%X",ord($1))/xge;
		push @pairs, "$key=$value";
	}
	$self->value( 9, join( '; ', @pairs ) );
	return 1;
}

sub rewrite_vcf_attributes {
	my $self = shift;
	return unless ( $self->{data}->vcf );
	return unless exists $self->{attributes};

	# INFO
	my $info = join(
		';',
		map {
			defined $self->{attributes}->{INFO}{$_}
				? join( '=', $_, $self->{attributes}->{INFO}{$_} )
				: $_
		}
			sort { $a cmp $b }
			keys %{ $self->{attributes}->{INFO} }
	);
	$info ||= '.';    # sometimes we have nothing left
	$self->value( 8, $info );

	# FORMAT
	my @order;
	push @order, 'GT' if exists $self->{attributes}{10}{GT};
	foreach my $key ( sort { $a cmp $b } keys %{ $self->{attributes}{10} } ) {
		next if $key eq 'GT';
		push @order, $key;
	}
	if (@order) {
		$self->value( 9, join( ':', @order ) );
	}
	else {
		$self->value( 9, '.' );
	}

	# SAMPLES
	foreach my $i ( 10 .. $self->{data}->last_column ) {
		if (@order) {
			$self->value( $i, join( ':', map { $self->{attributes}{$i}{$_} } @order ) );
		}
		else {
			$self->value( $i, '.' );
		}
	}
	return 1;
}

### Data collection convenience methods

*feature = \&seqfeature;
*feature if 0;    # avoid once warning

sub seqfeature {
	my $self  = shift;
	my $force = shift || 0;
	carp 'ERROR: feature is a read only method' if @_;
	return $self->{feature}                     if exists $self->{feature};

	# normally this is only for named features in a data table
	# skip this for coordinate features like bed files
	return unless $self->feature_type eq 'named' or $force;

	# retrieve from main Data store
	my $f = $self->{data}->get_seqfeature( $self->{'index'} );
	if ($f) {
		$self->{feature} = $f;
		return $f;
	}

	# retrieve the feature from the database
	return unless $self->{data}->database;
	$f = get_db_feature(
		'db'   => $self->{data}->open_meta_database,
		'id'   => $self->id   || undef,
		'name' => $self->name || undef,
		'type' => $self->type || $self->{data}->feature,
	);
	return unless $f;
	$self->{feature} = $f;
	return $f;
}

sub segment {
	my $self = shift;
	carp 'ERROR: segment is a read only method' if @_;
	return unless $self->{data}->database;
	if ( $self->feature_type eq 'coordinate' ) {
		my $chromo = $self->seq_id;
		my $start  = $self->start;
		my $stop   = $self->end || $start;
		my $db     = $self->{data}->open_meta_database;
		return $db ? $db->segment( $chromo, $start, $stop ) : undef;
	}
	elsif ( $self->feature_type eq 'named' ) {
		my $f = $self->feature;
		return $f ? $f->segment : undef;
	}
	else {
		return undef;
	}
}

sub get_features {
	my $self = shift;
	my %args = @_;
	my $db   = $args{db} || $self->{data}->open_meta_database || undef;
	carp 'ERROR: no database defined to get features!' unless defined $db;
	return                                             unless $db->can('features');

	# convert the argument style for most bioperl db APIs
	my %opts;
	$opts{-seq_id} = $args{chromo} || $self->seq_id;
	$opts{-start}  = $args{start}  || $self->start;
	$opts{-end}    = $args{end}    || $self->end;
	$opts{-type}   = $args{type}   || $self->type;

	return $db->features(%opts);
}

sub get_sequence {
	my $self = shift;
	my %args = @_;
	my $db   = $args{db} || $args{database} || $self->{data}->open_meta_database || undef;

	# this will fail immediately if user doesn't provide valid database

	# get sequence over subfeatures
	$args{subfeature} ||= undef;
	if ( $self->feature_type eq 'named' and $args{subfeature} ) {

		# this is more complicated so we have a dedicated method
		return $self->_get_subfeature_sequence( $db, \%args );
	}

	# get coordinates
	my $seqid  = $args{seq_id} || $args{chromo} || $self->seq_id;
	my $start  = $args{start}  || $self->start;
	my $stop   = $args{stop}   || $args{end} || $self->end;
	my $strand = $self->strand;
	if ( exists $args{strand} ) {

		# user supplied strand, gotta check it
		$strand = $args{strand} =~ /\-|r/i ? -1 : 1;
	}
	if ( exists $args{extend} and $args{extend} ) {
		$start -= $args{extend};
		$start = 1 if $start <= 0;
		$stop += $args{extend};
	}
	return unless ( defined $seqid and defined $start and defined $stop );

	# retrieve and return sequence
	my $seq = get_genomic_sequence( $db, $seqid, $start, $stop );
	if ( $strand == -1 ) {
		$seq =~ tr/gatcGATC/ctagCTAG/;
		$seq = reverse $seq;
	}
	return $seq;
}

sub _get_subfeature_sequence {
	my ( $self, $db, $args ) = @_;

	# get the subfeatures
	my $subfeatures = $self->_get_subfeatures( $args->{subfeature} );
	unless ( @{$subfeatures} ) {
		carp 'ERROR: no subfeatures available! Returning parent sequence!';

		# just return the parent
		undef $args->{subfeature};
		return $self->get_sequence( @{$args} );
	}

	# sort subfeatures
	# this should be done by GeneTools in most cases but just to be sure
	# note that this does NOT merge redundant or overlapping exons!!!!
	my @sorted = map { $_->[0] }
		sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] }
		map { [ $_, $_->start, $_->end ] } @{$subfeatures};

	# collect sequence
	my $sequence;
	foreach my $subf (@sorted) {
		my $seq = get_genomic_sequence( $db, $subf->seq_id, $subf->start, $subf->stop );
		$sequence .= $seq;
	}

	# flip the sequence
	if ( $self->strand == -1 ) {
		$sequence =~ tr/gatcGATC/ctagCTAG/;
		$sequence = reverse $sequence;
	}
	return $sequence;
}

sub _get_subfeatures {
	my $self = shift;
	my $subf = lc shift;

	# load GeneTools
	unless ($GENETOOL_LOADED) {
		load(
			'Bio::ToolBox::GeneTools', qw(get_exons get_cds get_5p_utrs get_3p_utrs
				get_introns)
		);
		if ($EVAL_ERROR) {
			croak "FATAL: missing required modules! $EVAL_ERROR";
		}
		else {
			$GENETOOL_LOADED = 1;
		}
	}

	# feature
	my $feature = $self->seqfeature;
	return unless ($feature);

	# get the subfeatures
	my @subfeatures;
	if ( $subf eq 'exon' ) {
		@subfeatures = get_exons($feature);
	}
	elsif ( $subf eq 'cds' ) {
		@subfeatures = get_cds($feature);
	}
	elsif ( $subf eq '5p_utr' ) {
		@subfeatures = get_5p_utrs($feature);
	}
	elsif ( $subf eq '3p_utr' ) {
		@subfeatures = get_3p_utrs($feature);
	}
	elsif ( $subf eq 'intron' ) {
		@subfeatures = get_introns($feature);
	}
	else {
		croak "ERROR: unrecognized subfeature parameter '$subf'!";
	}

	return \@subfeatures;
}

sub get_score {
	my $self = shift;
	my %args = @_;      # passed arguments to this method

	# verify the dataset for the user, cannot trust whether it has been done or not
	my $db = $args{ddb} || $args{db} || $self->{data}->open_meta_database || undef;
	$args{dataset} = $self->{data}->verify_dataset( $args{dataset}, $db );
	unless ( $args{dataset} ) {
		croak
'FATAL: provided dataset was unrecognized format or otherwise could not be verified!';
	}

	# get positioned scores over subfeatures only
	$args{subfeature} ||= q();
	if ( $self->feature_type eq 'named' and $args{subfeature} ) {

		# this is more complicated so we have a dedicated method
		return $self->_get_subfeature_scores( $db, \%args );
	}

	# build parameter array to pass on to the adapter
	my @params;

	# verify coordinates based on type of feature
	if ( $self->feature_type eq 'coordinate' ) {

		# coordinates are already in the table, use those
		$params[CHR]  = $args{seq_id} || $self->seq_id;
		$params[STRT] = $args{start}  || $self->start;
		$params[STOP] = $args{stop}   || $args{end} || $self->end;
		$params[STR] =
			( exists $args{strand} and defined $args{strand} )
			? $args{strand}
			: $self->strand;
	}
	elsif ( $self->feature_type eq 'named' ) {

		# must retrieve feature from the database first
		my $f = $self->seqfeature;
		return unless $f;
		$params[CHR]  = $args{seq_id} || $f->seq_id;
		$params[STRT] = $args{start}  || $f->start;
		$params[STOP] = $args{stop}   || $args{end} || $f->end;
		$params[STR] =
			( exists $args{strand} and defined $args{strand} )
			? $args{strand}
			: $f->strand;
	}
	else {
		croak
'FATAL: data table does not have identifiable coordinate or feature identification columns for score collection';
	}

	# adjust coordinates as necessary
	if ( exists $args{extend} and $args{extend} ) {
		$params[STRT] -= $args{extend};
		$params[STOP] += $args{extend};
	}

	# check coordinates
	$params[STRT] = 1 if $params[STRT] <= 0;
	if ( $params[STOP] < $params[STRT] ) {

		# coordinates are flipped, reverse strand
		return if ( $params[STOP] <= 0 );
		my $stop = $params[STRT];
		$params[STRT] = $params[STOP];
		$params[STOP] = $stop;
		$params[STR]  = -1;
	}
	return unless ( $params[CHR] and defined $params[STRT] );

	# score attributes
	$params[METH] = $args{'method'}     || 'mean';
	$params[STND] = $args{strandedness} || $args{stranded} || 'all';

	# other parameters
	$params[DB]   = $db;
	$params[RETT] = 0;                # return type should be a calculated value
	$params[DATA] = $args{dataset};

	# get the score
	return get_segment_score(@params);
}

sub _get_subfeature_scores {
	my ( $self, $db, $args ) = @_;

	# get the subfeatures
	my $subfeatures = $self->_get_subfeatures( $args->{subfeature} );
	unless ( @{$subfeatures} ) {

		# no subfeatures, nothing to collect
		return;
	}

	# collect over each subfeature
	my @scores;
	foreach my $exon ( @{$subfeatures} ) {
		my @params;            # parameters to pass on to adapter
		$params[CHR]  = $exon->seq_id;
		$params[STRT] = $exon->start;
		$params[STOP] = $exon->end;
		$params[STR]  = defined $args->{strand} ? $args->{strand} : $exon->strand;
		$params[STND] = $args->{strandedness} || $args->{stranded} || 'all';
		$params[METH] = $args->{method} || 'mean';
		$params[RETT] = 1;     # return type should be an array reference of scores
		$params[DB]   = $db;
		$params[DATA] = $args->{dataset};

		my $exon_scores = get_segment_score(@params);
		push @scores, @{$exon_scores} if defined $exon_scores;
	}

	# combine all the scores based on the requested method
	return calculate_score( $args->{method}, \@scores );
}

sub get_relative_point_position_scores {
	my $self = shift;
	my %args = @_;

	# get the database and verify the dataset
	my $ddb = $args{ddb} || $args{db} || $self->{data}->open_meta_database;
	$args{dataset} = $self->{data}->verify_dataset( $args{dataset}, $ddb );
	unless ( $args{dataset} ) {
		croak
'FATAL: provided dataset was unrecognized format or otherwise could not be verified!';
	}

	# assign some defaults
	$args{strandedness} ||= $args{stranded} || 'all';
	$args{position}     ||= 5;
	$args{coordinate}   ||= undef;
	$args{avoid}        ||= undef;
	$args{'method'}     ||= 'mean';    # in most cases this doesn't do anything
	unless ( $args{extend} ) {
		croak 'FATAL: must provide an extend value!';
	}
	$args{avoid} = undef unless ( $args{db} or $self->{data}->open_meta_database );

	# determine reference coordinate
	unless ( defined $args{coordinate} ) {
		$args{coordinate} = $self->calculate_reference( \%args );
	}

	# build parameter array to pass on to the adapter
	my @params;
	$params[CHR]  = $self->seq_id;
	$params[STRT] = $args{coordinate} - $args{extend};
	$params[STRT] = 1 if $params[STRT] < 1;                                 # sanity check
	$params[STOP] = $args{coordinate} + $args{extend};
	$params[STR]  = defined $args{strand} ? $args{strand} : $self->strand;
	$params[STND] = $args{strandedness};
	$params[METH] = $args{'method'};
	$params[RETT] = 2;      # return type should be a hash reference of positioned scores
	$params[DB]   = $ddb;
	$params[DATA] = $args{dataset};

	# Data collection
	my $pos2data = get_segment_score(@params);

	# Avoid positions
	if ( $args{avoid} ) {
		$self->_avoid_positions( $pos2data, \%args, $params[CHR], $params[STRT],
			$params[STOP] );
	}

	# covert to relative positions
	if ( $args{absolute} ) {

		# do not convert to relative positions
		return wantarray ? %{$pos2data} : $pos2data;
	}
	else {
		# return the collected dataset hash
		return $self->_convert_to_relative_positions( $pos2data,
			$args{coordinate}, $params[STR] );
	}
}

sub get_region_position_scores {
	my $self = shift;
	my %args = @_;

	# get the database and verify the dataset
	my $ddb = $args{ddb} || $args{db} || $self->{data}->open_meta_database;
	$args{dataset} = $self->{data}->verify_dataset( $args{dataset}, $ddb );
	unless ( $args{dataset} ) {
		croak
'FATAL: provided dataset was unrecognized format or otherwise could not be verified!';
	}

	# assign some defaults here, in case we get passed on to subfeature method
	$args{strandedness} ||= $args{stranded} || 'all';
	$args{extend}       ||= 0;
	$args{position}     ||= 5;
	$args{'method'}     ||= 'mean';    # in most cases this doesn't do anything
	$args{avoid} = undef unless ( $args{db} or $self->{data}->open_meta_database );

	# get positioned scores over subfeatures only
	$args{subfeature} ||= q();
	if ( $self->feature_type eq 'named' and $args{subfeature} ) {

		# this is more complicated so we have a dedicated method
		return $self->_get_subfeature_position_scores( \%args, $ddb );
	}

	# Assign coordinates
	# build parameter array to pass on to the adapter
	my @params;
	my $feature = $self->seqfeature || $self;
	$params[CHR]  = $args{chromo} || $args{seq_id} || $feature->seq_id;
	$params[STRT] = $args{start}  || $feature->start;
	$params[STOP] = $args{stop}   || $args{end} || $feature->end;
	$params[STR]  = defined $args{strand} ? $args{strand} : $feature->strand;
	if ( $args{extend} ) {
		$params[STRT] -= $args{extend};
		$params[STOP] += $args{extend};
		$params[STRT] = 1 if $params[STRT] < 1;    # sanity check
	}
	$params[STND] = $args{strandedness};
	$params[METH] = $args{method};
	$params[RETT] = 2;      # return type should be a hash reference of positioned scores
	$params[DB]   = $ddb;
	$params[DATA] = $args{dataset};

	# Data collection
	my $pos2data = get_segment_score(@params);

	# Avoid positions
	if ( $args{avoid} ) {
		$self->_avoid_positions( $pos2data, \%args, $params[CHR], $params[STRT],
			$params[STOP] );
	}

	# covert to relative positions
	if ( $args{absolute} ) {
		return wantarray ? %{$pos2data} : $pos2data;
	}
	else {
		# return data converted to relative positions
		unless ( defined $args{coordinate} ) {
			$args{coordinate} = $self->calculate_reference( \%args );
		}
		return $self->_convert_to_relative_positions( $pos2data,
			$args{coordinate}, $params[STR] );
	}
}

sub _get_subfeature_position_scores {
	my ( $self, $args, $ddb ) = @_;

	# get the subfeatures
	my $subfeatures = $self->_get_subfeatures( $args->{subfeature} );
	unless ( @{$subfeatures} ) {
		carp 'ERROR: no subfeatures available! Returning parent score data!';

		# just return the parent
		undef $args->{subfeature};
		delete $args->{exon} if exists $args->{exon};
		return $self->get_sequence( @{$args} );
	}

	# reset the practical start and stop to the actual subfeatures' final start and stop
	# we can no longer rely on the feature start and stop, consider CDS
	# these subfeatures should already be genomic sorted by GeneTools
	my $practical_start = $subfeatures->[0]->start;
	my $practical_stop  = $subfeatures->[-1]->end;

	# collect over each exon
	# we will adjust the positions of each reported score so that
	# it will appear as if all the exons are adjacent to each other
	# and no introns exist
	my $pos2data    = {};
	my $namecheck   = {};    # to check unique names when using ncount method....
	my $current_end = $practical_start;
	my $adjustment  = 0;
	my $fstrand     = defined $args->{strand} ? $args->{strand} : $self->strand;
	foreach my $exon ( @{$subfeatures} ) {

		my @params;          # parameters to pass on to adapter
		$params[CHR]  = $exon->seq_id;
		$params[STRT] = $exon->start;
		$params[STOP] = $exon->end;
		$params[STR]  = $fstrand;
		$params[STND] = $args->{strandedness};
		$params[METH] = $args->{method};
		$params[RETT] = 2;   # return type should be a hash reference of positioned scores
		$params[DB]   = $ddb;
		$params[DATA] = $args->{dataset};

		# collect scores
		my $exon_scores = get_segment_score(@params);

		# adjust the scores
		$adjustment = $params[STRT] - $current_end;
		$self->_process_exon_scores(
			$exon_scores,  $pos2data,  $adjustment, $params[STRT],
			$params[STOP], $namecheck, $args->{method}
		);

		# reset
		$current_end += $exon->length;
	}

	# collect extensions if requested
	if ( $args->{extend} ) {

		# left side
		my @params;          # parameters to pass on to adapter
		$params[CHR]  = $self->seq_id;
		$params[STRT] = $practical_start - $args->{extend};
		$params[STOP] = $practical_start - 1;
		$params[STR]  = $fstrand;
		$params[STND] = $args->{strandedness};
		$params[METH] = $args->{method};
		$params[RETT] = 2;   # return type should be a hash reference of positioned scores
		$params[DB]   = $ddb;
		$params[DATA] = $args->{dataset};

		my $ext_scores = get_segment_score(@params);

		# no adjustment should be needed
		$self->_process_exon_scores( $ext_scores, $pos2data, 0, $params[STRT],
			$params[STOP], $namecheck, $args->{method} );

		# right side
		# we can reuse our parameter array
		$params[STRT] = $practical_stop + 1;
		$params[STOP] = $practical_stop + $args->{extend};
		$ext_scores   = get_segment_score(@params);

		# the adjustment should be the same as the last exon
		$self->_process_exon_scores(
			$ext_scores,   $pos2data,  $adjustment, $params[STRT],
			$params[STOP], $namecheck, $args->{method}
		);
	}

	# covert to relative positions
	if ( $args->{absolute} ) {
		return wantarray ? %{$pos2data} : $pos2data;
	}
	else {
		# return data converted to relative positions
		# can no longer use original coordinates, but instead the new shifted coordinates
		$args->{practical_start} = $practical_start;
		$args->{practical_stop}  = $current_end;
		my $reference = $self->calculate_reference($args);
		return $self->_convert_to_relative_positions( $pos2data, $reference, $fstrand );
	}
}

sub calculate_reference {

	# my ($self, $args) = @_;
	# I need position
	# optionally alternate start, stop, and/or strand
	my $self = shift;
	my $args;
	if ( scalar @_ == 1 ) {

		# single variable passed
		if ( ref( $_[0] ) eq 'HASH' ) {
			$args = shift @_;
		}
		else {
			# assumed to be the position
			$args = { position => $_[0] };
		}
	}
	else {
		# multiple variables
		$args = {@_};
	}

	# calculate reference coordinate
	my $coordinate;
	my $strand = defined $args->{strand} ? $args->{strand} : $self->strand;
	if ( $args->{position} == 5 ) {
		if ( $strand >= 0 ) {
			$coordinate = $args->{practical_start} || $self->start;
		}
		else {
			$coordinate = $args->{practical_stop} || $self->end;
		}
	}
	elsif ( $args->{position} == 3 ) {
		if ( $strand >= 0 ) {
			$coordinate = $args->{practical_stop} || $self->end;
		}
		else {
			$coordinate = $args->{practical_start} || $self->start;
		}
	}
	elsif ( $args->{position} == 4 ) {

		# strand doesn't matter here
		# but practical coordinates do
		if ( exists $args->{practical_start} ) {
			$coordinate =
				int( ( ( $args->{practical_start} + $args->{practical_stop} ) / 2 ) +
					0.5 );
		}
		else {
			$coordinate = $self->midpoint;
		}
	}
	elsif ( $args->{position} == 10 or $args->{position} == 9 ) {

		# also accept old 0-based index of 9 for narrowPeak
		$coordinate = $self->peak;
	}
	else {
		confess 'FATAL: position must be one of 5, 3, 4, or 10';
	}
	return $coordinate;
}

sub _avoid_positions {
	my ( $self, $pos2data, $args, $seqid, $start, $stop ) = @_;

	# first check the list of avoid types
	if ( ref( $args->{avoid} ) eq 'ARRAY' ) {

		# we have types, presume they're ok
	}
	elsif ( $args->{avoid} eq '1' ) {

		# old style boolean value
		if ( defined $args->{type} ) {
			$args->{avoid} = [ $args->{type} ];
		}
		else {
			# no type provided, we can't avoid that which is not defined!
			# this is an error, but won't complain as we never did before
			$args->{avoid} = $self->type;
		}
	}
	elsif ( $args->{avoid} =~ /w+/i ) {

		# someone passed a string, a feature type perhaps?
		$args->{avoid} = [ $args->{avoid} ];
	}

	### Check for conflicting features
	my $db               = $args->{db} || $self->{data}->open_meta_database;
	my @overlap_features = $self->get_features(
		seq_id => $seqid,
		start  => $start,
		end    => $stop,
		type   => $args->{avoid},
	);

	# get the overlapping features of the same type
	if (@overlap_features) {
		my $primary = $self->primary_id;

		# there are one or more feature of the type in this region
		# one of them is likely the one we're working with
		# but not necessarily - user may be looking outside original feature
		# the others are not what we want and therefore need to be
		# avoided
		foreach my $feat (@overlap_features) {

			# skip the one we want
			next if ( $feat->primary_id eq $primary );

			# now eliminate those scores which overlap this feature
			my $f_start = $feat->start;
			my $f_stop  = $feat->end;
			foreach my $position ( keys %{$pos2data} ) {

				# delete the scored position if it overlaps with
				# the offending feature
				if (    $position >= $f_start
					and $position <= $f_stop )
				{
					delete $pos2data->{$position};
				}
			}
		}
	}
}

sub _convert_to_relative_positions {
	my ( $self, $pos2data, $position, $strand ) = @_;

	my %relative_pos2data;
	if ( $strand >= 0 ) {
		foreach my $p ( keys %{$pos2data} ) {
			$relative_pos2data{ $p - $position } = $pos2data->{$p};
		}
	}
	elsif ( $strand < 0 ) {
		foreach my $p ( keys %{$pos2data} ) {
			$relative_pos2data{ $position - $p } = $pos2data->{$p};
		}
	}
	return wantarray ? %relative_pos2data : \%relative_pos2data;
}

sub _process_exon_scores {
	my ( $self, $exon_scores, $pos2data, $adjustment, $start, $end, $namecheck, $method )
		= @_;

	# ncount method
	if ( $method eq 'ncount' ) {

		# we need to check both names and adjust position
		foreach my $p ( keys %{$exon_scores} ) {
			next unless ( $p >= $start and $p <= $end );
			foreach my $n ( @{ $exon_scores->{$p} } ) {
				if ( exists $namecheck->{$n} ) {
					$namecheck->{$n}++;
					next;
				}
				else {
					$namecheck->{$n} = 1;
					my $a = $p - $adjustment;
					$pos2data->{$a} ||= [];
					push @{ $pos2data->{$a} }, $n;
				}
			}
		}
	}
	else {
		# just adjust scores
		foreach my $p ( keys %{$exon_scores} ) {
			next unless ( $p >= $start and $p <= $end );
			$pos2data->{ $p - $adjustment } = $exon_scores->{$p};
		}
	}
}

sub fetch_alignments {
	my $self = shift;
	my %args = @_;
	$args{db}         ||= $args{dataset} || undef;
	$args{data}       ||= undef;
	$args{callback}   ||= undef;
	$args{subfeature} ||= q();

	# verify - trusting that these are valid, else they will fail lower down in the code
	unless ( $args{db} ) {
		croak 'FATAL: must provide a Bam object database to fetch alignments!';
	}
	unless ( $args{data} and ref( $args{data} ) eq 'HASH' ) {
		croak 'FATAL: must provide a data HASH for the fetch callback!';
	}
	unless ( $args{callback} ) {
		croak 'FATAL: must provide a callback code reference!';
	}

	# array of features to iterate, probably just one or subfeatures
	my @intervals;
	if ( $self->feature_type eq 'named' and $args{subfeature} ) {

		# we have subfeatures to iterate over

		# get the subfeatures
		my $subfeatures = $self->_get_subfeatures( $args{subfeature} );
		if ( @{$subfeatures} ) {
			foreach my $sf ( @{$subfeatures} ) {
				push @intervals, [ $sf->start - 1, $sf->end ];
			}
		}
		else {
			# zere subfeatures? just take the parent then
			push @intervals,
				[
					( $args{start} || $self->start ) - 1,
					$args{stop} || $args{end} || $self->end
				];
		}
	}
	else {
		# take feature as is
		push @intervals,
			[
				( $args{start} || $self->start ) - 1,
				$args{stop} || $args{end} || $self->end
			];
	}

	# get the target id for the chromosome
	# this will fail if the user didn't provide a real bam object!!!
	my ( $tid, undef, undef ) = $args{db}->header->parse_region( $self->seq_id );
	return unless defined $tid;

	# now iterate over the intervals
	foreach my $i (@intervals) {
		$args{data}->{start}  = $i->[0];
		$args{data}->{end}    = $i->[1];
		$args{data}->{strand} = $self->strand;
		low_level_bam_fetch( $args{db}, $tid, $i->[0], $i->[1], $args{callback},
			$args{data} );
	}

	# nothing to return since we're using a data reference
	return 1;
}

### String export

sub bed_string {
	my $self = shift;
	my %args = @_;
	$args{bed} ||= 6;    # number of bed columns
	croak 'FATAL: bed count must be an integer!' unless $args{bed} =~ /^\d+$/;
	croak 'FATAL: bed count must be at least 3!' unless $args{bed} >= 3;

	# coordinate information
	$self->seqfeature;    # retrieve the seqfeature object first
	my $chr   = $args{chromo} || $args{seq_id} || $self->seq_id;
	my $start = $args{start}  || $self->start;
	my $stop =
		   $args{stop}
		|| $args{end}
		|| $self->stop
		|| $start + $self->length - 1
		|| $start;
	if (   $chr eq '.'
		or not CORE::length($chr)
		or $start eq '.'
		or not CORE::length($start) )
	{
		carp sprintf( "ERROR: no valid seq_id or start for data line %d",
			$self->line_number );
		return;
	}
	if ( $start > $stop ) {

		# reversed coordinates? old school way of setting reverse strand
		my $s = $start;
		$start        = $stop;
		$stop         = $s;
		$args{strand} = -1;      # this will override any user provided data?
	}
	$start -= 1;                 # 0-based coordinates
	my $string = "$chr\t$start\t$stop";

	# additional information
	if ( $args{bed} >= 4 ) {
		my $name = $args{name} || $self->name || 'Feature_' . $self->line_number;
		$string .= "\t$name";
	}
	if ( $args{bed} >= 5 ) {
		my $score = exists $args{score} ? $args{score} : $self->score;
		$score = 1 unless defined $score;
		$string .= "\t$score";
	}
	if ( $args{bed} >= 6 ) {
		my $s;
		if ( exists $args{strand} and defined $args{strand} ) {
			$s = $self->_strand( $args{strand} );
		}
		else {
			$s = $self->strand;
		}
		$string .= sprintf( "\t%s", $s == -1 ? '-' : '+' );
	}

	# we could go on with other columns, but there's no guarantee that additional
	# information is available, and we would have to implement user provided data

	# done
	return $string;
}

sub gff_string {
	my $self = shift;
	my %args = @_;

	# coordinate information
	$self->seqfeature;    # retrieve the seqfeature object first
	my $chr   = $args{chromo} || $args{seq_id} || $self->seq_id;
	my $start = $args{start}  || $self->start;
	my $stop =
		   $args{stop}
		|| $args{end}
		|| $self->stop
		|| $start + $self->length - 1
		|| $start;
	if (   $chr eq '.'
		or not CORE::length($chr)
		or $start eq '.'
		or not CORE::length($start) )
	{
		carp sprintf( "ERROR: no valid seq_id or start for data line %d",
			$self->line_number );
		return;
	}
	if ( $start > $stop ) {

		# reversed coordinates? old school way of setting reverse strand
		my $s = $start;
		$start        = $stop;
		$stop         = $s;
		$args{strand} = -1;      # this will override any user provided data?
	}
	my $strand;
	if ( exists $args{strand} and defined $args{strand} ) {
		$strand = $self->_strand( $args{strand} );
	}
	else {
		$strand = $self->strand;
	}
	$strand = $strand == -1 ? '-' : $strand == 1 ? '+' : '.';

	# type information
	my $type = $args{type} || $self->type || undef;
	my ( $source, $primary_tag );
	if ( defined $type and $type =~ /:/ ) {
		( $primary_tag, $source ) = split /:/, $type;
	}
	unless ($source) {
		$source = $args{source} || '.';
	}
	unless ($primary_tag) {
		$primary_tag = $args{primary_tag} || defined $type ? $type : '.';
	}

	# score
	my $score = exists $args{score} ? $args{score} : $self->score;
	$score = '.' unless defined $score;
	my $phase = '.';    # do not even bother!!!!

	# attributes
	my $name       = $args{name} || $self->name || 'Feature_' . $self->line_number;
	my $attributes = "Name=$name";
	my $id         = $args{id} || sprintf( "%08d", $self->line_number );
	$attributes .= ";ID=$id";
	if ( exists $args{attributes} and ref( $args{attributes} ) eq 'ARRAY' ) {
		foreach my $i ( @{ $args{attributes} } ) {
			my $k = $self->{data}->name($i);
			$k =~ s/( [\t\n\r%&\=;,\ ] ) /sprintf("%%%X",ord($1))/xge;
			my $v = $self->value($i);
			$v =~ s/( [\t\n\r%&\=;,\ ] ) /sprintf("%%%X",ord($1))/xge;
			$attributes .= ";$k=$v";
		}
	}

	# done
	my $string = join( "\t",
		$chr,   $source, $primary_tag, $start, $stop,
		$score, $strand, $phase,       $attributes );
	return $string;
}

1;

__END__

=head1 NAME

Bio::ToolBox::Data::Feature - Objects representing rows in a data table

=head1 DESCRIPTION

A Bio::ToolBox::Data::Feature is an object representing a row in the 
data table. Usually, this in turn represents an annotated feature or 
segment in the genome. As such, this object provides convenient 
methods for accessing and manipulating the values in a row, as well as 
methods for working with the represented genomic feature.

This class should NOT be used directly by the user. Rather, Feature 
objects are generated from a Bio::ToolBox::Data::Iterator object 
(generated itself from the L<row_stream|Bio::ToolBox::Data/row_stream> 
function in Bio::ToolBox::Data), or the L<iterate|Bio::ToolBox::Data/iterate> 
function in Bio::ToolBox::Data. Please see the respective documentation 
for more information.

Example of working with a stream object.

	  my $Data = Bio::ToolBox::Data->new(file => $file);
	  
	  # stream method
	  my $stream = $Data->row_stream;
	  while (my $row = $stream->next_row) {
		 # each $row is a Bio::ToolBox::Data::Feature object
		 # representing the row in the data table
		 my $value = $row->value($index);
		 # do something with $value
	  }
	  
	  # iterate method
	  $Data->iterate( sub {
	     my $row = shift;
	     my $number = $row->value($index);
	     my $log_number = log($number);
	     $row->value($index, $log_number);
	  } );


=head1 METHODS

=head2 General information methods

=over 4

=item row_index

Returns the index position of the current data row within the 
data table. Useful for knowing where you are at within the data 
table.

=item feature_type

Returns one of three specific values describing the contents 
of the data table inferred by the presence of specific column 
names. This provides a clue as to whether the table features 
represent genomic regions (defined by coordinate positions) or 
named database features. The return values include:

=over 4

=item * coordinate: Table includes at least chromosome and start

=item * named: Table includes name, type, and/or Primary_ID

=item * unknown: unrecognized

=back

=item column_name

Returns the column name for the given index. 

item data

Returns the parent L<Bio::ToolBox::Data> object, in case you may 
have lost it by going out of scope.

=back

=head2 Methods to access row feature attributes

These methods return the corresponding value, if present in the 
data table, based on the column header name. If the row represents 
a named database object, try calling the L</feature> method first. 
This will retrieve the database SeqFeature object, and the attributes 
can then be retrieved using the methods below or on the actual 
database SeqFeature object.

In cases where there is a table column and a corresponding SeqFeature 
object, for example a start column and a parsed SeqFeature object, the 
table value takes precedence and is returned. You can always obtain the 
SeqFeature's value separately and directly.

These methods do not set attribute values. If you need to change the 
values in a table, use the L</value> method below.

=over 4

=item seq_id

The name of the chromosome the feature is on.

=item start

=item end

=item stop

The coordinates of the feature or segment. Coordinates from known 
0-based file formats, e.g. BED, are returned as 1-based. Coordinates 
must be integers to be returned. Zero or negative start coordinates 
are assumed to be accidents or poor programming and transformed to 1. 
Use the L</value> method if you don't want this to happen.

=item strand

The strand of the feature or segment. Returns -1, 0, or 1. Default is 0, 
or unstranded.

=item midpoint

The calculated midpoint position of the feature.

=item peak

For features in a C<narrowPeak> file, this will report the peak coordinate, 
transformed into a genomic coordinate. 

=item name

=item display_name

The name of the feature.

=item coordinate

Returns a coordinate string formatted as C<seqid:start-stop>. The start
coordinate is converted to 1-based where relevant, in concordance with
HTS tools (I<samtools> and I<tabix>).

=item type

The type of feature. Typically either C<primary_tag> or C<primary_tag:source_tag>. 
In a GFF3 file, this represents columns 3 and 2, respectively. In annotation 
databases such as L<Bio::DB::SeqFeature::Store>, the type is used to restrict 
to one of many different types of features, e.g. gene, mRNA, or exon.

=item id

=item primary_id

Here, this represents the C<primary_ID> in the database. Note that this number 
is generally unique to a specific database, and not portable between databases.

=item length

The length of the feature or segment.

=item score

Returns the value of the Score column, if one is available. Typically 
associated with defined file formats, such as GFF files (6th column), 
BED and related Peak files (5th column), and bedGraph (4th column).

=back

=head2 Accessing and setting values in the row.

=over 4

=item value

  # retrieve a value 
  my $v = $row->value($index);
  # set a value
  $row->value($index, $v + 1);

Returns or sets the value at a specific column index in the 
current data row. Null values return a '.', symbolizing an 
internal null value. 

=item row_values

Returns an array or array reference representing all the values 
in the current data row. 

=back

=head2 Special feature attributes

GFF and VCF files have special attributes in the form of key =E<gt> value pairs. 
These are stored as specially formatted, character-delimited lists in 
certain columns. These methods will parse this information and return as 
a convenient hash reference. The keys and values of this hash may be 
changed, deleted, or added to as desired. To write the changes back to 
the file, use the L</rewrite_attributes> to properly write the attributes 
back to the file with the proper formatting.

=over 4

=item attributes

Generic method that calls either L</gff_attributes> or L</vcf_attributes> 
depending on the data table format. 

=item gff_attributes

Parses the 9th column of GFF files. URL-escaped characters are converted 
back to text. Returns a hash reference of key =E<gt> value pairs. 

=item vcf_attributes

Parses the C<INFO> (8th column) and all sample columns (10th and higher 
columns) in a version 4 VCF file. The Sample columns use the C<FORMAT> 
column (9th column) as keys. The returned hash reference has two levels:
The first level keys are both the column names and index (0-based). The 
second level keys are the individual attribute keys to each value. 
For example:

   my $attr = $row->vcf_attributes;
   
   # access by column name
   my $genotype = $attr->{sample1}{GT};
   my $depth    = $attr->{INFO}{ADP};
   
   # access by 0-based column index 
   my $genotype = $attr->{9}{GT};
   my $depth    = $attr->{7}{ADP}

=item rewrite_attributes

Generic method that either calls L</rewrite_gff_attributes> or 
L</rewrite_vcf_attributes> depending on the data table format.

=item rewrite_gff_attributes

Rewrites the GFF attributes column (the 9th column) based on the 
contents of the attributes hash that was previously generated with 
the L</gff_attributes> method. Useful when you have modified the 
contents of the attributes hash.

=item rewrite_vcf_attributes

Rewrite the VCF attributes for the C<INFO> (8th column), C<FORMAT> (9th 
column), and sample columns (10th and higher columns) based on the 
contents of the attributes hash that was previously generated with 
the L</vcf_attributes> method. Useful when you have modified the 
contents of the attributes hash.

=back

=head2 Convenience Methods to database functions

The next three functions are convenience methods for using the 
attributes in the current data row to interact with databases. 
They are wrappers to methods in the L<Bio::ToolBox::db_helper> 
module.

=over 4

=item seqfeature

=item feature

Returns a SeqFeature object representing the feature or item in 
the current row. If the SeqFeature object is stored in the parent 
C<$Data> object (usually from parsing an annotation file), it is 
immediately returned. Otherwise, the SeqFeature 
object is retrieved from the database using the name and 
type values in the current Data table row. The SeqFeature object 
is requested from the database named in the general metadata. If 
an alternate database is desired, you should change it first using  
the C<$Data>-E<gt>database() method. If the feature name or type is not 
present in the table, then nothing is returned.

See L<Bio::ToolBox::SeqFeature> and L<Bio::SeqFeatureI> for more 
information about working with these objects. See L<Bio::DB::SeqFeature::Store> 
about working with database features.

This method normally only works with "named" feature types in a 
L<Bio::ToolBox::Data> Data table. If your Data table has coordinate 
information, i.e. chromosome, start, and stop columns, then it will 
likely be recognized as a "coordinate" feature_type and not work.

Pass a true value to this method to force the seqfeature lookup. This 
will still require the presence of Name, ID, and/or Type columns to 
perform the database lookup. The L<Bio::ToolBox::Data> method feature() 
is used to determine the type if a Type column is not present.

=item segment

Returns a database Segment object corresponding to the coordinates 
defined in the Data table row. If a named feature and type are 
present instead of coordinates, then the feature is first automatically 
retrieved and a Segment returned based on its coordinates. The 
database named in the general metadata is used to establish the 
Segment object. If a different database is desired, it should be 
changed first using the general L</database> method. 

See L<Bio::DB::SeqFeature::Segment> and L<Bio::RangeI> for more information 
about working with Segment objects.

=item get_features

  my @overlap_features = $row->get_features(type => $type);

Returns seqfeature objects from a database that overlap the Feature 
or interval in the current Data table row. This is essentially a 
convenience wrapper for a Bio::DB style I<features> method using the 
coordinates of the Feature. Optionally pass an array of key value pairs 
to specify alternate coordinates if so desired. Potential keys 
include 

=over 4

=item seq_id

=item start

=item end

=item type

The type of database features to retrieve.

=item db

An alternate database object to collect from.

=back

=item get_sequence

Fetches genomic sequence based on the coordinates of the current seqfeature 
or interval in the current Feature. This requires a database that 
contains the genomic sequence, either the database specified in the 
Data table metadata or an external indexed genomic fasta file. 

If the Feature represents a transcript or gene, then a concatenated 
sequence of the selected subfeatures may be generated and returned. B<Note> 
that redundant or overlapping subfeatures are B<NOT> merged, and 
unexpected results may be obtained.

The sequence is returned as simple string. If the feature is on the reverse 
strand, then the reverse complement sequence is automatically returned. 

Pass an array of key value pairs to specify alternate coordinates if so 
desired. Potential keys include

=over 4

=item subfeature 

Pass a text string representing the type of subfeature from which to collect 
the sequence. Acceptable values include 

=over 4 

=item * exon

=item * cds

=item * 5p_utr

=item * 3p_utr

=item * intron

=back

=item seq_id

=item start

=item end

=item strand

=item extend

Indicate additional basepairs of sequence added to both sides

=item db

The fasta file or database from which to fetch the sequence

=back

=back

=head2 Data collection

The following methods allow for data collection from various 
sources, including bam, bigwig, bigbed, useq, Bio::DB databases, etc. 

=over 4

=item calculate_reference($position)

Calculates and returns the absolute genomic coordinate for a relative 
reference position taking into account feature orientation (strand). 
This is not explicitly data collection, but often used in conjunction
with such. Provide an integer representing the relative position point:

=over 4 

=item * C<3> representing C<3'> end coordinate

=item * C<4> representing mid point coordinate

=item * C<5> representing C<5'> end coordinate

=item * C<9> representing peak summit in F<narrowPeak> formatted files

=back

If necessary, an array or array reference may be provided as an 
alternative parameter with keys including C<position>, C<strand>, 
C<practical_start>, and C<practical_end> if alternate or adjusted 
coordinates should be used instead of the given row feature coordinates.

=item get_score

  my $score = $row->get_score(
       dataset => 'scores.bw',
       method  => 'max',
  );

This method collects a single score over the feature or interval. 
Usually a mathematical or statistical value is employed to derive the 
single score. Pass an array of key value pairs to control data collection.
Keys include the following:

=over 4

=item db

=item ddb

Specify a Bio::DB database from which to collect the data. The default 
value is the database specified in the Data table metadata, if present.
Examples include a L<Bio::DB::SeqFeature::Store> or L<Bio::DB::BigWigSet> 
database.

=item dataset 

Specify the name of the dataset. If a database was specified, then this 
value would be the C<primary_tag> or C<type:source> feature found in the 
database. Otherwise, the name of a data file, such as a bam, bigWig, 
bigBed, or USeq file, would be provided here. This options is required!

=item method

Specify the mathematical or statistical method combining multiple scores 
over the interval into one value. Options include the following:

=over 4

=item * mean

=item * sum

=item * min

=item * max

=item * median

=item * count

Count all overlapping items.

=item * pcount

Precisely count only containing (not overlapping) items.

=item * ncount

Count overlapping unique names only.

=item * range

The difference between minimum and maximum values.

=item * stddev

Standard deviation.

=back

=item strandedness 

Specify what strand from which the data should be taken, with respect 
to the Feature strand. Three options are available. Only really relevant 
for data sources that support strand. 

=over 4

=item * sense

The same strand as the Feature.

=item * antisense

The opposite strand as the Feature.

=item * all

Strand is ignored, all is taken (default).

=back

=item subfeature

Specify the subfeature type from which to collect the scores. Typically 
a SeqFeature object representing a transcript is provided, and the 
indicated subfeatures are collected from object. Pass the name of the 
subfeature to use. Accepted values include the following.

=over 4 

=item * exon

=item * cds

=item * 5p_utr

=item * 3p_utr

=item * intron

=back

=item extend

Specify the number of basepairs that the Data table Feature's 
coordinates should be extended in both directions. Ignored 
when used with the subfeature option.

=item seq_id

=item chromo

=item start

=item end

=item stop

=item strand

Optionally specify zero or more alternate coordinates to use. 
By default, these are obtained from the Data table Feature.

=back
  
=item get_relative_point_position_scores

  while (my $row = $stream->next_row) {
     my $pos2score = $row->get_relative_point_position_scores(
        'ddb'       => '/path/to/BigWigSet/',
        'dataset'   => 'MyData',
        'position'  => 5,
        'extend'    => 1000,
     );
  }

This method collects indexed position scores centered around a 
specific reference point. The returned data is a hash of 
relative positions (example -20, -10, 1, 10, 20) and their score 
values. Pass an array of key value pairs to control data collection.
Keys include the following:

=over 4

=item db

=item ddb

Specify a Bio::DB database from which to collect the data. The default 
value is the database specified in the Data table metadata, if present.
Examples include a L<Bio::DB::SeqFeature::Store> or L<Bio::DB::BigWigSet> 
database.

=item dataset 

Specify the name of the dataset. If a database was specified, then this 
value would be the C<primary_tag> or C<type:source> feature found in the 
database. Otherwise, the name of a data file, such as a bam, bigWig, 
bigBed, or USeq file, would be provided here. This options is required!

=item position

Indicate the position of the reference point relative to the Data table 
Feature. 5 is the 5' coordinate, 3 is the 3' coordinate, and 4 is the 
midpoint (get it? it's between 5 and 3). Default is 5.

=item extend

Indicate the number of base pairs to extend from the reference coordinate. 
This option is required!

=item coordinate

Optionally provide the real chromosomal coordinate as the reference point.

=item absolute 

Boolean option to indicate that the returned hash of positions and scores 
should not be transformed into relative positions but kept as absolute 
chromosomal coordinates.

=item avoid

Provide a C<primary_tag> or C<type:source> database feature type to avoid overlapping 
scores. Each found score is checked for overlapping features and is 
discarded if found to do so. The database should be set to use this.

=item strandedness 

Specify what strand from which the data should be taken, with respect 
to the Feature strand. Three options are available. Only really relevant 
for data sources that support strand. 

=over 4

=item * sense

The same strand as the Feature.

=item * antisense

The opposite strand as the Feature.

=item * all

Strand is ignored, all is taken (default).

=back

=item method

Only required when counting objects.

=over 4

=item * count

Count all overlapping items.

=item * pcount

Precisely count only containing (not overlapping) items.

=item * ncount

Count overlapping unique names only.

=back

=back

=item get_region_position_scores

  while (my $row = $stream->next_row) {
     my $pos2score = $row->get_relative_point_position_scores(
        'ddb'       => '/path/to/BigWigSet/',
        'dataset'   => 'MyData',
        'position'  => 5,
        'extend'    => 1000,
     );
  }

This method collects indexed position scores across a defined 
region or interval. The returned data is a hash of positions and 
their score values. The positions are by default relative to a 
region coordinate, usually to the 5' end. Pass an array of key value 
pairs to control data collection. Keys include the following:

=over 4

=item db

=item ddb

Specify a Bio::DB database from which to collect the data. The default 
value is the database specified in the Data table metadata, if present.
Examples include a L<Bio::DB::SeqFeature::Store> or L<Bio::DB::BigWigSet> 
database.

=item dataset 

Specify the name of the dataset. If a database was specified, then this 
value would be the C<primary_tag> or C<type:source> feature found in the 
database. Otherwise, the name of a data file, such as a bam, bigWig, 
bigBed, or USeq file, would be provided here. This options is required!

=item subfeature

Specify the subfeature type from which to collect the scores. Typically 
a SeqFeature object representing a transcript is provided, and the 
indicated subfeatures are collected from object. When converting to 
relative coordinates, the coordinates will be relative to the length of 
the sum of the subfeatures, i.e. the length of the introns will be ignored.

Pass the name of the subfeature to use. Accepted values include the following.

=over 4 

=item * exon

=item * cds

=item * 5p_utr

=item * 3p_utr

=item * intron

=back

=item extend

Specify the number of basepairs that the Data table Feature's 
coordinates should be extended in both directions. 

=item seq_id

=item chromo

=item start

=item end

=item stop

=item strand

Optionally specify zero or more alternate coordinates to use. 
By default, these are obtained from the Data table Feature.

=item position

Indicate the position of the reference point relative to the Data table 
Feature. 5 is the 5' coordinate, 3 is the 3' coordinate, and 4 is the 
midpoint (get it? it's between 5 and 3). Default is 5.

=item coordinate

Optionally provide the real chromosomal coordinate as the reference point.

=item absolute 

Boolean option to indicate that the returned hash of positions and scores 
should not be transformed into relative positions but kept as absolute 
chromosomal coordinates.

=item avoid

Provide a C<primary_tag> or C<type:source> database feature type to avoid overlapping 
scores. Each found score is checked for overlapping features and is 
discarded if found to do so. The database should be set to use this.

=item strandedness 

Specify what strand from which the data should be taken, with respect 
to the Feature strand. Three options are available. Only really relevant 
for data sources that support strand. 

=over 4

=item * sense

The same strand as the Feature.

=item * antisense

The opposite strand as the Feature.

=item * all

Strand is ignored, all is taken (default).

=back

=item method

Only required when counting objects.

=over 4

=item * count

Count all overlapping items.

=item * pcount

Precisely count only containing (not overlapping) items.

=item * ncount

Count overlapping unique names only.

=back

=back

=item fetch_alignments

  my $sam = $Data->open_database('/path/to/file.bam');
  my $alignment_data = { mapq => [] };
  my $callback = sub {
     my ($a, $data) = @_;
     push @{ $data->{mapq} }, $a->qual;
  };
  while (my $row = $stream->next_row) {
     $row->fetch_alignments(
        'db'        => $sam,
        'data'      => $alignment_data,
        'callback'  => $callback,
     );
  }

This function allows you to iterate over alignments in a Bam file, 
allowing custom information to be collected based on a callback 
code reference that is provided. 

Three parameters are required: C<db>, C<data>, and C<callback>. A 
true value (1) is returned upon success.

=over 4

=item db

Provide an opened, high-level, Bam database object. 

=item callback

Provide a code callback reference to use when iterating over the 
alignments. Two objects are passed to this code function: the 
alignment object and the data structure that is provided. See the 
Bam adapter documentation for details on low-level C<fetch> through 
the Bam index object for details.

=item data

This is a reference to a C<HASH> data object for storing information. 
It is passed to the callback function along with the alignment. Three 
new key =E<gt> value pairs are automatically added: C<start>, C<end>, 
and C<strand>. These correspond to the values for the current queried 
interval. Coordinates are automatically transformed to 0-base coordinate 
system to match low level alignment objects.

=item subfeature

If the feature has subfeatures, such as exons, introns, etc., pass 
the name of the subfeature to restrict iteration only over the 
indicated subfeatures. The C<data> object will inherit the coordinates 
for each subfeatures. Allowed subfeatures include the following:

=over 4 

=item * exon

=item * cds

=item * 5p_utr

=item * 3p_utr

=item * intron

=back

=item start

=item stop

=item end

Provide alternate, custom start and stop coordinates for the row 
feature. Ignored with subfeatures.

=back

=back

=head2 Feature Export

These methods allow the feature to be exported in industry standard 
formats, including the BED format and the GFF format. Both methods 
return a formatted tab-delimited text string suitable for printing to 
file. The string does not include a line ending character.

These methods rely on coordinates being present in the source table. 
If the row feature represents a database item, the L</feature> method 
should be called prior to these methods, allowing the feature to be 
retrieved from the database and coordinates obtained.

=over 4

=item bed_string

Returns a BED formatted string. By default, a 6-element string is 
generated, unless otherwise specified. Pass an array of key values 
to control how the string is generated. The following arguments 
are supported.

=over 4

=item bed

Specify the number of BED elements to include. The number of elements 
correspond to the number of columns in the BED file specification. A 
minimum of 3 (chromosome, start, stop) is required, and maximum of 6 
is allowed (chromosome, start, stop, name, score, strand). 

=item chromo

=item seq_id

Provide a text string of an alternative chromosome or sequence name.

=item start

=item stop

=item end

Provide alternative integers for the start and stop coordinates. 
Note that start values are automatically converted to 0-base 
by subtracting 1.

=item strand

Provide alternate an alternative strand value. 

=item name

Provide an alternate or missing name value to be used as text in the 4th 
column. If no name is provided or available, a default name is generated.

=item score

Provide a numerical value to be included as the score. BED files typically 
use integer values ranging from 1..1000. 

=back

=item gff_string

Returns a GFF3 formatted string. Pass an array of key values 
to control how the string is generated. The following arguments 
are supported.

=over 4

=item chromo

=item seq_id

=item start

=item stop

=item end

=item strand

Provide alternate values from those defined or missing in the current 
row Feature. 

=item source

Provide a text string to be used as the source_tag value in the 2nd 
column. The default value is null ".".

=item primary_tag

Provide a text string to be used as the primary_tag value in the 3rd 
column. The default value is null ".".

=item type

Provide a text string. This can be either a "primary_tag:source_tag" value 
as used by GFF based BioPerl databases, or "primary_tag" alone.

=item score

Provide a numerical value to be included as the score. The default 
value is null ".". 

=item name

Provide alternate or missing name value to be used as the display_name. 
If no name is provided or available, a default name is generated.

=item attributes

Provide an anonymous array reference of one or more row Feature indices 
to be used as GFF attributes. The name of the column is used as the GFF 
attribute key. 

=back

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
