package Bio::ToolBox::Parser::gff;

use warnings;
use strict;
use Carp qw(carp cluck croak confess);
use base 'Bio::ToolBox::Parser';
use Bio::ToolBox::Data;

our $VERSION = '2.02';

my %TYPECOUNTS = ();

sub new {
	my $class = shift;
	return $class->SUPER::new(@_);
}

sub open_file {
	my $self     = shift;
	my $filename = shift || undef;

	# check file
	if ( $filename and $self->file and $filename ne $self->file ) {
		confess 'FATAL: Must open new files with new Parser object!';
	}
	$filename ||= $self->file;
	unless ($filename) {
		carp 'ERROR: No file name passed!';
		return;
	}
	if ( defined $self->{fh} ) {
		return 1;
	}

	# check file format type
	my $filetype = $self->filetype || undef;
	unless ($filetype) {
		( my $flavor, $filetype ) = Bio::ToolBox::Data->taste_file($filename);
		unless ( $flavor eq 'gff' ) {
			confess 'FATAL: File is not a GFF file!!! How did we get here?';
		}
		$self->{filetype} = $filetype;
		if ( $filetype eq 'gtf' ) {
			$self->{gtf} = 1;
		}
		else {
			$self->{gtf} = 0;
		}
	}

	# check type list
	# this must be done prior to opening because it informs parsing subroutines
	# but very limited sampling, not thorough
	my $typelist = Bio::ToolBox::Data->sample_gff_type_list($filename);
	if ( $typelist !~ /\w+/ ) {
		print "GFF file has no evident types!? $filename may not be a valid GFF file\n";
		return;
	}
	$self->{typelist} = $typelist;

	# determine converter subroutine
	if ( $filetype eq 'gtf' ) {
		$self->{convertor_sub} = \&_gtf_to_seqf;

		# double check we have transcript information
		if ( $self->typelist !~ m/transcript | rna/xi ) {

			# we will have to rely on exon and/or cds information to get transcript
			unless ( $self->do_exon or $self->do_cds ) {
				$self->do_exon(1);
			}
		}
		if ( $self->do_gene and $self->typelist !~ /gene/i ) {
			$self->do_exon(1);
		}
	}
	elsif ( $filetype eq 'gff3' ) {
		$self->{convertor_sub} = \&_gff3_to_seqf;
	}
	else {
		$self->{convertor_sub} = \&_gff2_to_seqf;
	}

	# Open filehandle object
	my $fh = Bio::ToolBox::Data->open_to_read_fh($filename)
		or croak "FATAL: cannot open file '$filename'!";
	$self->{fh} = $fh;
	return 1;
}

sub typelist {
	return shift->{typelist} || q();
}

sub next_feature {
	my $self = shift;

	# check that we have an open filehandle
	unless ( $self->fh ) {
		croak 'FATAL: no GFF file loaded to parse!';
	}
	return if $self->{'eof'};

	# look for the next feature line
	while ( my $line = $self->{fh}->getline ) {

		# check first character
		my $firstchar = substr( $line, 0, 1 );

		# skip any comment and pragma lines that we might encounter
		if ( $firstchar eq '#' ) {
			if ( $line =~ /^###$/ ) {

				# a close pragma, all we can do is check for orphans
				# a properly written shouldn't have any orphans, but just in case
				$self->check_orphanage;
				next;
			}
			elsif ( $line =~ m/^\#\#sequence . region/xi ) {

				# sequence region pragma
				my ( $pragma, $seq_id, $start, $stop ) = split /\s+/, $line;
				if ( defined $seq_id and $start =~ /^\d+$/ and $stop =~ /^\d+$/ ) {

					# we're actually only concerned with the stop coordinate
					$self->{seq_ids}{$seq_id} = $stop;
				}
				else {
					print STDERR "malformed sequence-region pragma! $line\n";
				}
				next;
			}
			else {
				# must be some sort of pragma or a comment line, may be useful, keep it
				push @{ $self->{comments} }, $line;
				next;
			}
		}
		elsif ( $firstchar eq '>' ) {

			# fasta header line
			# this is almost always at the end of the file, and rarely is sequence put
			# into GFF files anyway, so let's assume it's the end of the file
			$self->{'eof'} = 1;
			return;
		}

		# line must be a GFF feature
		chomp $line;
		my @fields = split /\t/, $line;
		next unless scalar(@fields) == 9;

		# check the primary_tag and generate the SeqFeature object for known types
		my $type = lc $fields[2];
		if ( $type eq 'cds' ) {
			if ( $self->do_cds ) {
				return &{ $self->{convertor_sub} }( $self, \@fields );
			}
			else {
				next;
			}
		}
		elsif ( $type eq 'exon' ) {
			if ( $self->do_exon ) {
				return &{ $self->{convertor_sub} }( $self, \@fields );
			}
			else {
				next;
			}
		}
		elsif ( $type =~ m/utr | untranslated/x ) {
			if ( $self->do_utr ) {
				return &{ $self->{convertor_sub} }( $self, \@fields );
			}
			else {
				next;
			}
		}
		elsif ( $type =~ /codon/ ) {
			if ( $self->do_codon ) {
				return &{ $self->{convertor_sub} }( $self, \@fields );
			}
			else {
				next;
			}
		}
		elsif ( $type =~ /gene$/ ) {
			if ( $self->do_gene ) {
				return &{ $self->{convertor_sub} }( $self, \@fields );
			}
			else {
				next;
			}
		}
		elsif ( $type =~ m/transcript | rna/x ) {
			return &{ $self->{convertor_sub} }( $self, \@fields );
		}
		elsif ( $type =~ m/(?: chromosome | contig | scaffold )/x ) {

			# gff3 files can record the chromosome as a gff record
			# process this as a region
			$self->{seq_ids}{ $fields[0] } = $fields[4];
			next;
		}
		else {
			# everything else must be some non-standard weird element
			# only process this if the user wants everything
			return &{ $self->{convertor_sub} }( $self, \@fields ) unless $self->simplify;
		}
	}

	# presumably reached the end of the file
	$self->{'eof'} = 1;
	return;
}

*parse_table = \&parse_file;
*parse_table if 0;    # avoid once warning

sub parse_file {
	my $self = shift;

	# check that we have an open filehandle
	unless ( $self->fh ) {
		confess 'FATAL: no file loaded to parse!';
	}
	return 1 if $self->{'eof'};

	# Each line will be processed into a SeqFeature object, and then checked
	# for parentage. Child objects with a Parent tag will be appropriately
	# associated with its parent, or put into an orphanage. Any orphans
	# left will be checked a final time for parents; if a parent can't be
	# found, it will be lost. Features without a parent are assumed to be
	# top-level features.

	printf "  Parsing %s %s format file....\n", $self->simplify ? 'simply' : 'fully',
		$self->filetype;

TOP_FEATURE_LOOP:
	while ( my $feature = $self->next_feature ) {

		### Process the feature
		# add genes, transcripts, and all other parent features
		# to the loaded hash for lookup
		if ( $feature->primary_tag =~ m/(?: gene | rna | transcript )/xi
			or not $feature->has_tag('Parent') )
		{
			# remember this feature as it likely will have children features
			# this ID should be unique in the GFF file
			# otherwise it might be a shared duplicate or a malformed GFF file
			# generally only a concern for top level features
			my $result = $self->_add_loaded_feature($feature);
			if ( $result > 1 ) {
				$self->{duplicate_ids}{ $feature->primary_id } += 1;
			}
		}

		# look for parents and children
		if ( $feature->has_tag('Parent') ) {

			# must be a child
			# there may be more than one parent, per the GFF3 specification
			foreach my $parent_id ( $feature->get_tag_values('Parent') ) {
				my $parent = $self->_get_loaded_feature( $parent_id, qr/\w+/, $feature );
				if ($parent) {

					# associate the child with the parent
					$parent->add_SeqFeature($feature);

					# check boundaries for gtf genes
					# gtf genes may not be explicitly defined so must correct as necessary
					# gff3 files shouldn't have this issue
					if ( $self->{gtf} and $parent->has_tag('autogenerate') ) {
						if ( $feature->start < $parent->start ) {
							$parent->start( $feature->start );
						}
						if ( $feature->end > $parent->end ) {
							$parent->end( $feature->end );
						}

						# check parent's parent too
						if ( $parent->has_tag('Parent') ) {

							# in all likelihood parent is a transcript and there is a
							# gene that probably also needs fixin'
							my ($grandparent_id) = $parent->get_tag_values('Parent');
							my $grandparent = $self->_get_loaded_feature( $grandparent_id,
								qr/gene/, $feature );
							if ( $feature->start < $grandparent->start ) {
								$grandparent->start( $feature->start );
							}
							if ( $feature->end > $grandparent->end ) {
								$grandparent->end( $feature->end );
							}
						}
					}
				}
				else {
					# can't find the parent, maybe not loaded yet?
					# put 'em in the orphanage
					$self->_add_orphan($feature);
				}
			}
		}
		else {
			# must be a parent
			push @{ $self->{top_features} }, $feature;

			# check chromosome length of top features
			my $s = $feature->seq_id;
			unless ( exists $self->{seq_ids}{$s} ) {
				$self->{seq_ids}{$s} = $feature->end;
			}
			$self->{seq_ids}{$s} = $feature->end if $feature->end > $self->{seq_ids}{$s};
		}
	}

	# Finished loading the GFF lines

	# check for orphans
	if ( scalar @{ $self->{orphans} } ) {
		$self->check_orphanage;

		# report
		if ( scalar @{ $self->{orphans} } ) {
			printf
" GFF errors: %d features could not be associated with reported parents!\n",
				scalar @{ $self->{orphans} };
			printf " List: %s\n", join ', ', map { $_->primary_id } @{ $self->{orphans} };
		}
	}

	# report on duplicate IDs
	if ( keys %{ $self->{duplicate_ids} } ) {
		printf " GFF errors: %d IDs were duplicated\n",
			scalar( keys %{ $self->{duplicate_ids} } );
		printf " List: %s\n", join ', ', keys %{ $self->{duplicate_ids} };
	}

	return 1;
}

sub _make_gene_parent {

	# for generating GTF gene parent features
	my ( $self, $fields, $att ) = @_;
	my $gene = $self->{sfclass}->new(
		-seq_id      => $fields->[0],
		-source      => $fields->[1],
		-primary_tag => 'gene',
		-start       => $fields->[3],
		-end         => $fields->[4],
		-strand      => $fields->[6],
		-primary_id  => $att->{'gene_id'},
	);

	# add extra information if possible
	if ( exists $att->{'gene_name'} ) {
		$gene->display_name( $att->{'gene_name'} );
	}
	if ( exists $att->{'gene_biotype'} ) {
		$gene->add_tag_value( 'gene_biotype', $att->{'gene_biotype'} );
	}
	elsif ( exists $att->{'gene_type'} ) {
		$gene->add_tag_value( 'gene_type', $att->{'gene_type'} );
	}
	if ( exists $att->{'gene_source'} ) {
		$gene->add_tag_value( 'gene_source', $att->{'gene_source'} );
	}
	$gene->add_tag_value( 'autogenerate', 1 );    # extra key for auto-generation
	return $gene;
}

sub _make_rna_parent {

	# for generating GTF gene parent features
	my ( $self, $fields, $att ) = @_;
	my $rna = $self->{sfclass}->new(
		-seq_id      => $fields->[0],
		-source      => $fields->[1],
		-primary_tag => 'transcript',
		-start       => $fields->[3],
		-end         => $fields->[4],
		-strand      => $fields->[6],
		-primary_id  => $att->{'transcript_id'},
	);

	# add extra information if possible
	if ( exists $att->{'transcript_name'} ) {
		$rna->display_name( $att->{'transcript_name'} );
	}
	if ( exists $att->{'transcript_biotype'} ) {
		$rna->add_tag_value( 'transcript_biotype', $att->{'transcript_biotype'} );
	}
	elsif ( exists $att->{'transcript_type'} ) {
		$rna->add_tag_value( 'transcript_biotype', $att->{'transcript_type'} );
	}
	if ( exists $att->{'transcript_source'} ) {
		$rna->add_tag_value( 'transcript_source', $att->{'transcript_source'} );
	}
	$rna->add_tag_value( 'autogenerate', 1 );    # extra key for auto-generation
	return $rna;
}

sub _gff3_to_seqf {
	my ( $self, $fields ) = @_;
	my $feature = $self->_gff_to_seqf($fields);

	# process the group tags
	$fields->[8] =~ s/;$//;                      # remove any trailing semi-colon
	my %att = map { split /=/ } split /;\s?/, $fields->[8];

	# add essential attributes
	if ( exists $att{'ID'} ) {
		$feature->primary_id( $att{'ID'} );
	}
	else {
		# generate one
		my $t = $feature->primary_tag;
		$TYPECOUNTS{$t} += 1;
		$feature->primary_id( sprintf( "%s.%d", $t, $TYPECOUNTS{$t} ) );
	}
	if ( exists $att{'Name'} ) {

		# name may be encoded
		my $n = $self->unescape( $att{'Name'} );
		$feature->display_name($n);
	}
	if ( exists $att{'Parent'} ) {

		# always record Parent except for transcripts when genes are not wanted
		unless ( not $self->do_gene and $feature->primary_tag =~ m/rna | transcript/xi ) {
			$feature->add_tag_value( 'Parent', $att{'Parent'} );
		}
	}
	if ( lc $feature->primary_tag eq 'exon' and exists $att{'exon_id'} ) {

		# Ensembl GFF3 stores the exon id but doesn't record it as the ID, why?
		# should not affect parentage as Ensembl doesn't link children
		$feature->primary_id( $att{'exon_id'} );
	}

	# extra attributes as necessary
	unless ( $self->simplify ) {
		foreach my $k (
			grep { !m/(?: Name | ID | Parent | exon_id )/x }
			keys %att
			)
		{
			my $key = $k =~ /%/ ? $self->unescape($k) : $k;
			my $value;
			if ( index( $att{$key}, ',' ) > 0 ) {
				my @a =
					map { index( $_, '%' ) >= 0 ? $self->unescape($_) : $_ }
					split /,/, $att{$key};
				$value = \@a;
			}
			else {
				$value =
					index( $att{$key}, '%' ) >= 0
					? $self->unescape( $att{$key} )
					: $att{$key};
			}
			$feature->add_tag_value( $key, $value );
		}
	}

	return $feature;
}

sub _gtf_to_seqf {
	my ( $self, $fields ) = @_;
	my $feature = $self->_gff_to_seqf($fields);

	# process the group tags
	$fields->[8] =~ s/;$//;    # remove any trailing semi-colon
	my %att;
	foreach ( split /;\s?/, $fields->[8] ) {
		my ( $k, $v ) = split /\s/, $_, 2;
		next unless ( $k and $v );
		$v =~ s/"//g;

		# attribute keys should be unique, but Ensembl re-uses 'tag' key, so tolerate
		if ( exists $att{$k} ) {
			if ( ref( $att{$k} ) eq 'ARRAY' ) {
				push @{ $att{$k} }, $v;
			}
			else {
				my $current = $att{$k};
				$att{$k} = [ $current, $v ];
			}
		}
		else {
			$att{$k} = $v;
		}
	}

	# common attributes
	my $transcript_id = $att{'transcript_id'} || undef;
	my $gene_id       = $att{'gene_id'}       || undef;
	unless ( $gene_id or $transcript_id ) {

		# improperly formatted GTF file without one of these two items, nothing more to do
		return $feature;
	}

	# assign special tags based on the feature type
	my $type = lc $fields->[2];
	if (   $type eq 'exon'
		or $type eq 'cds'
		or $type eq 'utr'
		or $type eq 'five_prime_utr'
		or $type eq 'three_prime_utr'
		or $type eq 'start_codon'
		or $type eq 'stop_codon' )
	{
		$feature->add_tag_value( 'Parent', $transcript_id );

		# exon id if present
		if ( $type eq 'exon' and exists $att{'exon_id'} ) {
			$feature->primary_id( $att{'exon_id'} );
		}
		else {
			# generate an ID
			$TYPECOUNTS{$type} += 1;
			$feature->primary_id( sprintf( "%s.%d", $type, $TYPECOUNTS{$type} ) );
		}

		# check gene parent
		my $gene;
		if ( $self->do_gene and $gene_id ) {
			$gene = $self->_get_loaded_feature( $gene_id, qr/gene/, $feature );
			unless ($gene) {
				$gene = $self->_make_gene_parent( $fields, \%att );
				my $result = $self->_add_loaded_feature($gene);
				if ( $result > 1 ) {
					$self->{duplicate_ids}{$gene_id} += 1;
				}
				push @{ $self->{top_features} }, $gene;
			}
		}

		# check transcript parent
		if ($transcript_id) {
			my $rna = $self->_get_loaded_feature( $transcript_id,
				qr/(?: rna | transcript)/x, $feature );
			unless ($rna) {
				$rna = $self->_make_rna_parent( $fields, \%att );
				my $result = $self->_add_loaded_feature($rna);
				if ( $result > 1 ) {
					$self->{duplicate_ids}{$gene_id} += 1;
				}
				if ($gene) {
					$rna->add_tag_value( 'Parent', $gene_id );
					$gene->add_SeqFeature($rna);
				}
				else {
					# this will become a top feature
					push @{ $self->{top_features} }, $rna;
					if ($gene_id) {
						$rna->add_tag_value( 'gene_id', $gene_id );
					}
				}
			}
		}
	}

	# transcript
	elsif ( $type eq 'transcript' or $type =~ m/rna/ ) {

		# these are sometimes present in GTF files, such as from Ensembl

		# check if this was previously autogenerated
		my $existing = $self->_get_loaded_feature( $transcript_id, $type, $feature );
		if ( $existing and $existing->has_tag('autogenerate') ) {

			# this may happen when lines are out of order
			# rather than replace we just update
			$existing->start( $feature->start );
			$existing->stop( $feature->stop );
			$existing->display_name( $feature->display_name );
			$existing->remove_tag('autogenerate');

			# update additional attributes as necessary
			unless ( $self->simplify ) {
				$self->_add_remaining_gtf_attributes( $existing, \%att );
			}

			# move on to next feature
			return $self->next_feature;
		}

		# add transcript information
		$feature->primary_id($transcript_id);    # this should be present!!!
		if ( exists $att{'transcript_name'} ) {
			$feature->display_name( $att{'transcript_name'} );
		}
		if ( exists $att{'transcript_biotype'} ) {
			$feature->add_tag_value( 'transcript_biotype', $att{'transcript_biotype'} );
		}
		elsif ( exists $att{'transcript_type'} ) {
			$feature->add_tag_value( 'transcript_type', $att{'transcript_type'} );
		}

		# check gene parent
		if ( $self->do_gene and $gene_id ) {
			my $gene = $self->_get_loaded_feature( $gene_id, qr/gene/, $feature );
			unless ($gene) {
				$gene = $self->_make_gene_parent( $fields, \%att );
				my $result = $self->_add_loaded_feature($gene);
				if ( $result > 1 ) {
					$self->{duplicate_ids}{$gene_id} += 1;
				}
				push @{ $self->{top_features} }, $gene;
			}
			$feature->add_tag_value( 'Parent', $gene_id );
		}
		elsif ($gene_id) {

			# not doing gene parents, then at least add the gene_id tag back
			$feature->add_tag_value( 'gene_id', $gene_id );
		}
	}

	# gene
	elsif ( $type eq 'gene' ) {

		# these are sometimes present in GTF files, such as from Ensembl
		# but are not required and often absent

		# check if this was previously autogenerated
		my $existing = $self->_get_loaded_feature( $gene_id, $type, $feature );
		if ( $existing and $existing->has_tag('autogenerate') ) {

			# this may happen when lines are out of order
			# rather than replace we just update
			$existing->start( $feature->start );
			$existing->stop( $feature->stop );
			$existing->display_name( $feature->display_name );
			$existing->remove_tag('autogenerate');

			# update additional attributes as necessary
			unless ( $self->simplify ) {
				$self->_add_remaining_gtf_attributes( $existing, \%att );
			}

			# move on to next feature
			return $self->next_feature;
		}
		else {
			# add gene information
			$feature->primary_id($gene_id);    # this should be present!!!
			if ( exists $att{'gene_name'} ) {
				$feature->display_name( $att{'gene_name'} );
			}
			if ( exists $att{'gene_biotype'} ) {
				$feature->add_tag_value( 'gene_biotype', $att{'gene_biotype'} );
			}
			elsif ( exists $att{'gene_type'} ) {
				$feature->add_tag_value( 'gene_type', $att{'gene_type'} );
			}
		}
	}

	# something else
	else {
		# generate an ID
		$TYPECOUNTS{$type} += 1;
		$feature->primary_id( sprintf( "%s.%d", $type, $TYPECOUNTS{$type} ) );

		# add parent - hope it's already made from CDS or exon features
		if ($transcript_id) {
			$feature->add_tag_value( 'Parent', $transcript_id );
		}
	}

	# store remaining attributes, if any
	unless ( $self->simplify ) {
		$self->_add_remaining_gtf_attributes( $feature, \%att );
	}

	return $feature;
}

sub _add_remaining_gtf_attributes {
	my ( $self, $feature, $att ) = @_;
	foreach my $key ( keys %{$att} ) {
		next if $key eq 'transcript_id';
		next if $key eq 'transcript_name';
		next if $key eq 'gene_id';
		next if $key eq 'gene_name';
		next if $key eq 'exon_id';
		unless ( $feature->has_tag($key) ) {
			$feature->add_tag_value( $key, $att->{$key} );
		}
	}
}

sub _gff2_to_seqf {

	# generic gff1 or gff2 format or poorly defined gff3 or gtf file
	# hope for the best!
	my ( $self, $fields ) = @_;
	my $feature = $self->_gff_to_seqf($fields);

	# process groups
	# we have no uniform method of combining features, so we'll leave the tags
	# as is and hope for the best
	foreach my $g ( split( /\s*;\s*/, $fields->[8] ) ) {
		my ( $tag, $value ) = split /\s+/, $g;
		next unless ( $tag and $value );
		$feature->add_tag_value( $tag, $value );
	}
	return $feature;
}

sub _gff_to_seqf {
	my ( $self, $fields ) = @_;

	# generate the basic SeqFeature
	my $feature = $self->{sfclass}->new(
		-seq_id      => $fields->[0],
		-source      => $fields->[1],
		-primary_tag => $fields->[2],
		-start       => $fields->[3],
		-end         => $fields->[4],
		-strand      => $fields->[6],
	);

	# add more attributes if they're not null
	if ( $fields->[5] ne '.' ) {
		$feature->score( $fields->[5] );
	}
	if ( $fields->[7] ne '.' ) {
		$feature->phase( $fields->[7] );
	}

	# finished
	return $feature;
}

sub unescape {

	# Borrowed unashamedly from bioperl Bio::Tools::GFF
	# which in turn was borrowed from Bio::DB::GFF
	my $self = shift;
	my $v    = shift;
	$v =~ tr/+/ /;
	$v =~ s/%( [0-9a-fA-F]{2} ) /chr hex($1)/xge;
	return $v;
}

sub _add_orphan {
	my ( $self, $feature ) = @_;
	push @{ $self->{orphans} }, $feature;
	return 1;
}

sub check_orphanage {
	my $self = shift;
	return unless scalar @{ $self->{orphans} };

	# go through the list of orphans
	my @reunited;    # list of indices to delete after reuniting orphan with parent
	for ( my $i = 0; $i < scalar @{ $self->{orphans} }; $i++ ) {
		my $orphan  = $self->{orphans}->[$i];
		my $success = 0;

		# find the parent
		foreach my $parent ( $orphan->get_tag_values('Parent') ) {
			my $existing = $self->_get_loaded_feature( $parent, qr/\w+/, $orphan );
			if ($existing) {
				$existing->add_SeqFeature($orphan);
				$success++;
			}
		}

		# delete the orphan from the array if it found it's long lost parent
		push @reunited, $i if $success;
	}

	# clean up the orphanage
	while (@reunited) {
		my $i = pop @reunited;
		splice( @{ $self->{orphans} }, $i, 1 );
	}
}

sub orphans {
	my $self = shift;
	my @orphans;
	foreach ( @{ $self->{orphans} } ) {
		push @orphans, $_;
	}
	return wantarray ? @orphans : \@orphans;
}

sub _get_loaded_feature {
	my ( $self, $id, $type, $feature ) = @_;
	if ( exists $self->{loaded}{$id} ) {
		if ( ref( $self->{loaded}{$id} ) eq 'ARRAY' ) {

			# multiple objects, need to check
			foreach my $f ( @{ $self->{loaded}{$id} } ) {

				# we could do a more stringent intersection test than just overlap
				# but that might miss what we're really after, especially if the
				# parent was autogenerated and needs to be expanded
				# return the first one that matches
				# if there are multiple then there are serious issues
				if ( $f->primary_tag =~ m/$type/ and $feature->overlaps($f) ) {
					return $f;
				}
			}
			return;
		}
		else {
			# assume a single SeqFeature object
			return $self->{loaded}{$id};
		}
	}
	else {
		return;
	}
}

sub _add_loaded_feature {
	my ( $self, $feature ) = @_;
	my $id = $feature->primary_id;
	if ( exists $self->{loaded}{$id} ) {

		# store all of the features as an array
		if ( ref( $self->{loaded}{$id} ) eq 'ARRAY' ) {

			# there's more than two duplicates! pile it on!
			my $check = 1;
			foreach my $f ( @{ $self->{loaded}{$id} } ) {
				$check++ if ( $feature->primary_tag eq $f->primary_tag );
			}
			push @{ $self->{loaded}{$id} }, $feature;
			return $check;
		}
		else {
			my $existing = $self->{loaded}{$id};
			$self->{loaded}{$id} = [ $existing, $feature ];
			if ( $existing->primary_tag eq $feature->primary_tag ) {
				return 2;
			}
			else {
				return 1;
			}
		}
	}
	else {
		# unique ID, so remember it
		$self->{loaded}{$id} = $feature;
		return 1;
	}
	return;
}

1;

__END__

=head1 NAME

Bio::ToolBox::Parser::gff - parse GFF3, GTF, and generic GFF files 

=head1 SYNOPSIS

  use Bio::ToolBox::Parser;
  my $filename = 'file.gff3';
  
  my $Parser = Bio::ToolBox::Parser->new(
  	file    => $filename,
  	do_gene => 1,
  	do_exon => 1,
  ) or die "unable to open gff file!\n";
  # the Parser will taste the file and open the appropriate 
  # subclass parser, gff in this case
  
  while (my $feature = $Parser->next_top_feature() ) {
	# each $feature is a parent SeqFeature object, usually a gene
  	printf "%s:%d-%d\n", $f->seq_id, $f->start, $f->end;
	
	# subfeatures such as transcripts, exons, etc are nested within
	my @children = $feature->get_SeqFeatures();
  }

=head1 DESCRIPTION

This is the GFF specific parser subclass to the L<Bio::ToolBox::Parser>
object, and as such inherits generic methods from the parent. 

This module parses a GFF file into SeqFeature objects. It natively 
handles GFF3, GTF, and general GFF files. 

For both GFF3 and GTF files, fully nested gene models, typically 
gene =E<gt> transcript =E<gt> (exon, CDS, etc), may be built using the appropriate 
attribute tags. For GFF3 files, these include ID and Parent tags; for GTF 
these include C<gene_id> and C<transcript_id> tags. 

For GFF3 files, any feature without a C<Parent> tag is assumed to be a 
parent. Children features referencing a parent feature that has not been 
loaded are considered orphans. Orphans are attempted to be re-associated 
with missing parents after the file is completely parsed. Any orphans left 
may be collected. Files with orphans are considered poorly formatted or 
incomplete and should be fixed. Multiple parentage, for example exons 
shared between different transcripts of the same gene, are fully supported.

Embedded Fasta sequences are ignored.

The SeqFeature objects that are returned are L<Bio::ToolBox::SeqFeature> 
objects (default SeqFeature class). Refer to that documentation for more 
information.

=head1 METHODS

=head2 Initialize and modify the parser. 

In most cases, users should initialize an object using the generic 
L<Bio::ToolBox::Parser> object.

These are class methods to initialize the parser with an annotation file 
and modify the parsing behavior. Most parameters can be set either upon 
initialization or as class methods on the object. Unpredictable behavior 
may occur if you implement these in the midst of parsing a file. 

Do not open subsequent files with the same object. Always create a new 
object to parse a new file

=over 4

=item new

  my $parser = Bio::ToolBox::Parser::gff->new($filename);
  my $parser = Bio::ToolBox::Parser::gff->new(
      file    => 'file.gtf.gz',
      do_gene => 1,
      do_utr  => 1,
  );

Initialize a new gff parser object. Pass a single value (a GFF file name) 
to automatically open (but not parse yet) a file. Alternatively, pass an 
array of key value pairs to control how the file is parsed. Options include 
the following.

=over 4

=item file

Provide a GFF file name to be parsed. It should have a gff, gtf, or gff3 file 
extension. The file may be gzip compressed. 

=item simplify

Pass a boolean value to simplify the SeqFeature objects parsed from the GFF 
file and ignore extraneous attributes.

=item do_gene

Pass a boolean (1 or 0) value to combine multiple transcripts with the same gene 
name under a single gene object. Default is true.

=item do_cds

=item do_exon

=item do_utr

=item do_codon

Pass a boolean (1 or 0) value to parse certain subfeatures. Exon subfeatures 
are always parsed, but C<CDS>, C<five_prime_UTR>, C<three_prime_UTR>, C<stop_codon>, 
and C<start_codon> features may be optionally parsed. Default is false.

=item class

Pass the name of a L<Bio::SeqFeatureI> compliant class that will be used to 
create the SeqFeature objects. The default is to use L<Bio::ToolBox::SeqFeature>, 
which is lighter-weight and consumes less memory. A suitable BioPerl alternative
is L<Bio::SeqFeature::Lite>.

=item simplify

Pass a boolean true value to simplify the attributes of GFF3 and GTF files 
that may have considerable numbers of tags, e.g. Ensembl files. Only 
essential information, including name, ID, and parentage, is retained. 
Useful if you're trying to quickly parse annotation files for basic 
information.

=back

=back

=head2 Access methods

See L<Bio::ToolBox::Parser> for generic methods for accessing the 
features. Below are specific methods to this subclass.

=over 4

=item open_file

Opens a new annotation file in a new object. Normally, this is automatically
done when a Parser object is instantiated with the L<new> method and a file
path was provided. Do not attempt to open a subsequent file with a
pre-existing Parser object; it will fail. 

Pass the path to a GFF annotation file. It may be compressed with gzip.
Success returns 1.

=item next_feature

This will parse the opened file one line at a time, returning the SeqFeature
object for each line of the GFF file. Parent-E<gt>child relationships are
not built. This should only be used for simple files or when
Parent-E<gt>child relationships are not needed.

=item parse_file
=item parse_table

This will parse the opened file entirely into memory, parsing each feature
into SeqFeature objects and assembling into parent-E<gt>child features.
B<NOTE> that for vertebrate genome annotation, this may consume considerable
amount of memory and take a while.

The L<check_orphanage> method is automatically run upon parsing the entire
file.

=item orphans

  my @orphans = $parser->orphans;
  printf "we have %d orphans left over!", scalar @orpans;

This method will return an array of orphan SeqFeature objects that indicated 
they had a parent but said parent could not be found. Typically, this is an 
indication of an incomplete or malformed GFF3 file. Nevertheless, it might 
be a good idea to check this after retrieving all top features.

=item check_orphanage

Method to go through the list of orphan sub features (if any) and attempt to
reunite them with their indicated Parent object if it has been loaded
into memory. This is automatically run when L<parse_file> is called and all
features have been loaded.

=item typelist

Returns a comma-delimited string of the GFF primary tags (column 3) 
observed in the first 1000 lines of an opened file. Useful for 
checking what is in the GFF file. See 
L<Bio::ToolBox::Data::file/sample_gff_type_list>.

=item unescape

This method will unescape special characters in a text string. Certain 
characters, including ";" and "=", are reserved for GFF3 formatting and 
are not allowed, thus requiring them to be escaped.

=back

=head1 SEE ALSO

L<Bio::ToolBox::Parser>, L<Bio::ToolBox::SeqFeature>, L<Bio::Tools::GFF>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

