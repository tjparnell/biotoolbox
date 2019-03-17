package Bio::ToolBox::parser::gff;

our $VERSION = '1.70';

=head1 NAME

Bio::ToolBox::parser::gff - parse GFF3, GTF, and GFF files 

=head1 SYNOPSIS

  use Bio::ToolBox::parser::gff;
  my $filename = 'file.gff3';
  
  my $parser = Bio::ToolBox::parser::gff->new(
  	file    => $filename,
  	do_gene => 1,
  	do_exon => 1,
  ) or die "unable to open gff file!\n";
  
  while (my $feature = $parser->next_top_feature() ) {
	# each $feature is a SeqFeature object
	my @children = $feature->get_SeqFeatures();
  }

=head1 DESCRIPTION

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

Embedded Fasta sequences are ignored, as are most comment and pragma lines.

The SeqFeature objects that are returned are L<Bio::ToolBox::SeqFeature> 
objects. Refer to that documentation for more information.

=head1 METHODS

=head2 Initialize and modify the parser. 

These are class methods to initialize the parser with an annotation file 
and modify the parsing behavior. Most parameters can be set either upon 
initialization or as class methods on the object. Unpredictable behavior 
may occur if you implement these in the midst of parsing a file. 

Do not open subsequent files with the same object. Always create a new 
object to parse a new file

=over 4

=item new

  my $parser = Bio::ToolBox::parser::gff->new($filename);
  my $parser = Bio::ToolBox::parser::gff->new(
      file    => 'file.gtf.gz',
      do_gene => 1,
      do_utr  => 1,
  );

Initialize a new gff parser object. Pass a single value (a GFF file name) 
to open a file. Alternatively, pass an array of key value pairs to control 
how the file is parsed. Options include the following.

=over 4

=item file

Provide a GFF file name to be parsed. It should have a gff, gtf, or gff3 file 
extension. The file may be gzip compressed. 

=item version

Specify the version. Normally this is not needed, as version can be determined 
either from the file extension (in the case of gtf and gff3) or from the 
C<##gff-version> pragma at the top of the file. Acceptable values include 1, 2, 
2.5 (gtf), or 3.

=item class

Pass the name of a L<Bio::SeqFeatureI> compliant class that will be used to 
create the SeqFeature objects. The default is to use L<Bio::ToolBox::SeqFeature>.

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

=back

=item open_file

  $parser->open_file($file) or die "unable to open $file!";

Pass the name of a GFF file to be parsed. The file may optionally be gzipped 
(F<.gz> extension). Do not open a new file when one has already opened a file. 
Create a new object for a new file, or concatenate the GFF files.

=item version

Set or get the GFF version of the current file. Acceptable values include 1, 2, 
2.5 (gtf), or 3. Normally this is determined by file extension or 
C<gff-version> pragma on the first line, and should not need to be set by the 
user in most circumstances.

=item simplify

Pass a boolean true value to simplify the attributes of GFF3 and GTF files 
that may have considerable numbers of tags, e.g. Ensembl files. Only 
essential information, including name, ID, and parentage, is retained. 
Useful if you're trying to quickly parse annotation files for basic 
information.

=back

=head2 Feature retrieval

The following methods parse the GFF file lines into SeqFeature objects. 
It is best if these methods are not mixed; unexpected results may occur. 

=over 4

=item next_top_feature

This method will return a top level parent SeqFeature object 
assembled with child features as sub-features. For example, a gene 
object with mRNA subfeatures, which in turn may have exon and/or CDS 
subfeatures. Child features are assembled based on the existence of 
proper Parent attributes in child features. If no Parent attributes are 
included in the GFF file, then this will behave as L</next_feature>.

Child features (those containing a C<Parent> attribute) 
are associated with the parent feature. A warning will be issued about lost 
children (orphans). Shared subfeatures, for example exons common to 
multiple transcripts, are associated properly with each parent. An opportunity 
to rescue orphans is available using the L</orphans> method.

Note that subfeatures may not necessarily be in ascending genomic order 
when associated with the feature, depending on their order in the GFF3 
file and whether shared subfeatures are present or not. When calling 
subfeatures in your program, you may want to sort the subfeatures. For 
example

  my @subfeatures = map { $_->[0] }
                    sort { $a->[1] <=> $b->[1] }
                    map { [$_, $_->start] }
                    $parent->get_SeqFeatures;

=item top_features

This method will return an array of the top (parent) features defined in 
the GFF file. This is similar to the L</next_top_feature> method except that 
all features are returned at once. 

=item next_feature

This method will return a SeqFeature object representation of 
the next feature (line) in the file. Parent - child relationships are 
NOT assembled; however, undefined parents in a GTF file may still be 
generated, just not returned. 

This method is best used with simple GFF files with no hierarchies 
present. This may be used in a while loop until the end of the file 
is reached. Pragmas are ignored and comment lines and sequence are 
automatically skipped. 

=back

=head2 Other methods

Additional methods for working with the parser object and the parsed 
SeqFeature objects.

=over 4

=item fh

This method returns the L<IO::File> object of the opened GFF file. 

=item parse_file

Parses the entire file into memory. This is automatically called when 
either L</top_features> or L</next_top_feature> is called. 

=item fetch

  my $gene = $parser->fetch($primary_id) or 
     warn "gene $display_name can not be found!";

Fetch a loaded top feature from memory using the C<primary_id> tag, which 
should be unique. Returns the SeqFeature object or C<undef> if not present.
Only useful after </parse_file> is called. 

=item orphans

  my @orphans = $parser->orphans;
  printf "we have %d orphans left over!", scalar @orpans;

This method will return an array of orphan SeqFeature objects that indicated 
they had a parent but said parent could not be found. Typically, this is an 
indication of an incomplete or malformed GFF3 file. Nevertheless, it might 
be a good idea to check this after retrieving all top features.

=item comments

This method will return an array of the comment or pragma lines that may have 
been in the parsed file. These may or may not be useful.

=item typelist

Returns a comma-delimited string of the GFF primary tags (column 3) 
observed in the first 1000 lines of an opened file. Useful for 
checking what is in the GFF file. See 
L<Bio::ToolBox::Data::file/sample_gff_type_list>.

=item from_gff_string

  my $seqfeature = $parser->from_gff_string($string);

This method will parse a single GFF, GTF, or GFF3 formatted string or line 
of text and return a SeqFeature object.

=item unescape

This method will unescape special characters in a text string. Certain 
characters, including ";" and "=", are reserved for GFF3 formatting and 
are not allowed, thus requiring them to be escaped.

=item seq_ids

Returns an array or array reference of the names of the chromosomes or 
reference sequences present in the file. These may be defined by GFF3 
sequence-region pragmas or inferred from the features.

=item seq_id_lengths

  my $seq2len = $parser->seq_id_lengths;
  foreach (keys %$seq2len) {
    printf "chromosome %s is %d bp long\n", $_, $seq2len->{$_};
  }

Returns a hash reference to the chromosomes or reference sequences and 
their corresponding lengths. In this case, the length is either defined 
by the C<sequence-region> pragma or inferred by the greatest end position of 
the top features.

=back

=head1 SEE ALSO

L<Bio::ToolBox::SeqFeature>, L<Bio::ToolBox::parser::ucsc>, L<Bio::Tools::GFF>

=cut

use strict;
use Carp qw(carp cluck croak);
use base 'Bio::ToolBox::Parser'; 
use Bio::ToolBox::Data; 

1;

sub new {
	my $class = shift;
	return $class->SUPER::new(@_);
}


sub open_file {
	my $self = shift;
	
	# check file
	my $filename = shift;
	unless ($filename) {
		cluck("no file name passed!\n");
		return;
	}
	
	# check type list
	my $typelist = Bio::ToolBox::Data->sample_gff_type_list($filename);
	if ($typelist !~ /\w+/) {
		print "GFF file has no evident types!? $filename may not be a valid GFF file\n";
		return;
	}
	$self->{typelist} = $typelist;
	
	# Open filehandle object 
	my $fh = Bio::ToolBox::Data->open_to_read_fh($filename) or
		croak " cannot open file '$filename'!\n";
	
	# check gff version pragma
	my $first = $fh->getline;
	if ($first =~ /^##gff.version\s+([\d\.]+)\s*$/i) {
		# override any version that may have been inferred from the extension
		# based on the assumption that this pragma is correct
		$self->{version} = $1;
	}
	else {
		# no pragma, reopen the file
		$fh->close;
		$fh = Bio::ToolBox::Data->open_to_read_fh($filename);
		# set version based on file type extension????
		if ($filename =~ /\.gtf.*$/i) {
			$self->{version} = '2.5';
			$self->{gtf} = 1;
		}
		elsif ($filename =~ /\.gff3.*$/i) {
			$self->{version} = '3';
			$self->{gff3} = 1;
		}
	}
	$self->{fh} = $fh;
	return 1;
}

sub typelist {
	return shift->{typelist};
}

sub next_feature {
	my $self = shift;
	
	# check that we have an open filehandle
	unless ($self->fh) {
		croak("no GFF file loaded to parse!");
	}
	return if $self->{'eof'};
	
	# check convertor
	unless (defined $self->{convertor_sub}) {
		$self->_assign_gff_converter;
	}
	
	# look for the next feature line
	while (my $line = $self->{fh}->getline) {
		
		# check first character
		my $firstchar = substr($line, 0, 1);
		
		# skip any comment and pragma lines that we might encounter
		if ($firstchar eq '#') {
			if ($line =~ /^###$/) {
				# a close pragma, all we can do is check for orphans
				# a properly written shouldn't have any orphans, but just in case
				$self->check_orphanage;
				next;
			}
			elsif ($line =~ /^##sequence.region/i) {
				# sequence region pragma
				my ($pragma, $seq_id, $start, $stop) = split /\s+/, $line;
				if (defined $seq_id and $start =~ /^\d+$/ and $stop =~ /^\d+$/) {
					# we're actually only concerned with the stop coordinate
					$self->{seq_ids}{$seq_id} = $stop;
				}
				else {
					warn "malformed sequence-region pragma! $line\n";
				}
				next;
			}
			else {
				# must be some sort of pragma or a comment line, may be useful, keep it
				push @{$self->{comments}}, $line;
				next;
			}
		}
		elsif ($firstchar eq '>') {
			# fasta header line
			# this is almost always at the end of the file, and rarely is sequence put 
			# into GFF files anyway, so let's assume it's the end of the file 
			$self->{'eof'} = 1;
			return;
		}
		
		# line must be a GFF feature 
		chomp $line;
		my @fields = split('\t', $line);
		next unless scalar(@fields) == 9;
		
		# check the primary_tag and generate the SeqFeature object for known types
		my $type = lc $fields[2];
		if ($type eq 'cds') {
			if ($self->do_cds) {
				return &{$self->{convertor_sub}}($self, \@fields);
			} else {
				next;
			}
		}
		elsif ($type eq 'exon') {
			if ($self->do_exon) {
				return &{$self->{convertor_sub}}($self, \@fields);
			} else {
				next;
			}
		}
		elsif ($type =~ /utr|untranslated/) {
			if ($self->do_utr) {
				return &{$self->{convertor_sub}}($self, \@fields);
			} else {
				next;
			}
		}
		elsif ($type =~ /codon/) {
			if ($self->do_codon) {
				return &{$self->{convertor_sub}}($self, \@fields);
			} else {
				next;
			}
		}
		elsif ($type =~ /gene$/) {
			if ($self->do_gene) {
				return &{$self->{convertor_sub}}($self, \@fields);
			} else {
				next;
			}
		}
		elsif ($type =~ /transcript|rna/) {
			return &{$self->{convertor_sub}}($self, \@fields);
		}
		elsif ($type =~ /chromosome|contig|scaffold/) {
			# gff3 files can record the chromosome as a gff record
			# process this as a region
			$self->{seq_ids}{$fields[0]} = $fields[4];
			next;
		}
		else {
			# everything else must be some non-standard weird element
			# only process this if the user wants everything
			return &{$self->{convertor_sub}}($self, \@fields) unless $self->simplify;
		}
	}
	
	# presumably reached the end of the file
	$self->{'eof'} = 1;
	return;
}

*parse_table = \&parse_file;

sub parse_file {
	my $self = shift;
	# check that we have an open filehandle
	unless ($self->fh) {
		confess("no file loaded to parse!");
	}
	return 1 if $self->{'eof'};
	
	# Each line will be processed into a SeqFeature object, and then checked 
	# for parentage. Child objects with a Parent tag will be appropriately 
	# associated with its parent, or put into an orphanage. Any orphans 
	# left will be checked a final time for parents; if a parent can't be 
	# found, it will be lost. Features without a parent are assumed to be 
	# top-level features.
	
	printf "  Parsing %s format file....\n", 
		$self->{version} eq '3' ? 'GFF3' : 
		$self->{version} =~ /2\../ ? 'GTF' : 'GFF';
	
	
	TOP_FEATURE_LOOP:
	while (my $feature = $self->next_feature) {
		
		### Process the feature
		if ($feature->primary_tag =~ /gene|rna|transcript/) {
			# remember this feature as it likely will have children features
			
			my $id = $feature->primary_id;
			# remember this feature since we have an ID
			if (exists $self->{loaded}{$id}) {
				# this ID should be unique in the GFF file
				# otherwise it might be a shared duplicate or a malformed GFF file
				# generally only a concern for top level features
				$self->{duplicate_ids}{$id} += 1;
				
				# store all of the features as an array
				if (ref($self->{loaded}{$id}) eq 'ARRAY') {
					# there's more than two duplicates! pile it on!
					push @{ $self->{loaded}{$id} }, $feature; 
				}
				else {
					my $existing = $self->{loaded}{$id};
					$self->{loaded}{$id} = [$existing, $feature];
				}
			} 
			else {
				# unique ID, so remember it
				$self->{loaded}{$id} = $feature;
			}
		}
		# if the feature didn't have an ID, we'll just assume it is
		# a child of another feature, otherwise it may get lost
		
		# look for parents and children
		if ($feature->has_tag('Parent')) {
			# must be a child
			# there may be more than one parent, per the GFF3 specification
			foreach my $parent_id ( $feature->get_tag_values('Parent') ) {
				if (exists $self->{loaded}{$parent_id}) {
					# we've seen this id
					# associate the child with the parent
					my $parent = $self->{loaded}{$parent_id};
					$parent->add_SeqFeature($feature);
					
					# check boundaries for gtf genes
					# gtf genes may not be explicitly defined so must correct as necessary
					# gff3 files won't have this issue
					if ($self->{gtf}) {
						if ($feature->start < $parent->start) {
							$parent->start( $feature->start );
						}
						if ($feature->end > $parent->end) {
							$parent->end( $feature->end );
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
			
			# check chromosome
			my $s = $feature->seq_id;
			unless (exists $self->{seq_ids}{$s}) {
				$self->{seq_ids}{$s} = $feature->end;
			}
			$self->{seq_ids}{$s} = $feature->end if $feature->end > $self->{seq_ids}{$s};
		}
	}
	# Finished loading the GFF lines
	
	# check for orphans
	if (scalar @{ $self->{orphans} }) {
		$self->check_orphanage;
		# report
		if (scalar @{ $self->{orphans} }) {
			printf " GFF errors: %d features could not be associated with reported parents!\n", 
				scalar @{ $self->{orphans} };
		}
	}
	
	# report on duplicate IDs
	if (keys %{ $self->{duplicate_ids} }) {
		printf " GFF errors: %d IDs were duplicated\n", 
			scalar(keys %{ $self->{duplicate_ids} });
	}
	
	return 1;
}


sub _make_gene_parent {
	# for generating GTF gene parent features
	my ($self, $fields, $gene_id) = @_;
	my $gene = $self->{sfclass}->new(
		-seq_id         => $fields->[0],
		-primary_tag    => 'gene',
		-start          => $fields->[3],
		-end            => $fields->[4],
		-strand         => $fields->[6],
		-primary_id     => $gene_id,
	);
	
	if ($fields->[8] =~ /gene_name "([^"]+)";?/) {
		$gene->display_name($1);
	}
	else {
		$gene->display_name($gene_id);
	}
	return $gene;
}


sub _make_rna_parent {
	# for generating GTF gene parent features
	my ($self, $fields, $transcript_id) = @_;
	my $rna = $self->{sfclass}->new(
		-seq_id         => $fields->[0],
		-primary_tag    => 'transcript',
		-start          => $fields->[3],
		-end            => $fields->[4],
		-strand         => $fields->[6],
		-primary_id     => $transcript_id,
	);
	
	if ($fields->[8] =~ /transcript_name "([^"]+)";?/) {
		$rna->display_name($1);
	}
	else {
		$rna->display_name($transcript_id);
	}
	
	# add extra information if possible
	unless ($self->simplify) {
		if ($fields->[8] =~ /transcript_biotype "([^"]+)";?/) {
			$rna->add_tag_value('transcript_biotype', $1);
		}
		elsif ($fields->[8] =~ /transcript_type "([^"]+)";?/) {
			$rna->add_tag_value('transcript_type', $1);
		}
		if ($fields->[8] =~ /transcript_source "([^"]+)";?/) {
			$rna->source($1);
		}
	}
	
	return $rna;
}


sub from_gff_string {
	my ($self, $string) = @_;
	
	# check string
	chomp $string;
	unless ($string) {
		cluck("must pass a string!\n");
		return;
	}
	my @fields = split('\t', $string);
	
	# check convertor
	unless (defined $self->{convertor_sub}) {
		$self->_assign_gff_converter;
	}
	
	# parse appropriately
	return &{$self->{convertor_sub}}($self, \@fields);
}


sub _assign_gff_converter {
	my $self = shift;
	if ($self->{version} eq '3') {
		$self->{gff3} = 1;
		$self->{convertor_sub} = \&_gff3_to_seqf;
	}
	elsif ($self->{version} eq '2.5' or $self->{version} eq '2.2') {
		$self->{gtf} = 1;
		if ($self->simplify) {
			$self->{convertor_sub} = \&_gtf_to_seqf_simple;
		}
		else {
			$self->{convertor_sub} = \&_gtf_to_seqf_full;
		}
		# double check we have transcript information
		if ($self->typelist !~ /transcript|rna/i) {
			# we will have to rely on exon and/or cds information to get transcript 
			unless ($self->do_exon or $self->do_cds) {
				$self->do_exon(1);
			}
		}
		if ($self->do_gene and $self->typelist !~ /gene/i) {
			$self->do_exon(1);
		}
	}
	else {
		$self->{convertor_sub} = \&_gff2_to_seqf;
	}
}


sub _gff3_to_seqf {
	my ($self, $fields) = @_;
	my $group = $fields->[8];
	my $feature = $self->_gff_to_seqf($fields);
	
	# process groups
	foreach my $g (split(/\s*;\s*/, $group)) {
		my ($tag, $value) = split /=/, $g;
		$tag = $self->unescape($tag);
		my @values = map { $self->unescape($_) } split(/,/, $value);
		
		# determine the appropriate attribute based on tag name
		if ($tag eq 'Name') {
			$feature->display_name($values[0]);
		}
		elsif ($tag eq 'ID') {
			$feature->primary_id($values[0]);
		}
		elsif ($tag eq 'Parent') {
			# always record Parent except for transcripts when genes are not wanted
			unless (not $self->do_gene and $feature->primary_tag =~ /rna|transcript/i) {
				foreach (@values) {
					$feature->add_tag_value($tag, $_);
				}
			}
		}
		elsif (lc $tag eq 'exon_id') {
			# ensembl GFF3 store the exon id but doesn't record it as the ID, why?
			$feature->primary_id($values[0]);
		}
		elsif (not $self->{simplify}) {
			foreach (@values) {
				$feature->add_tag_value($tag, $_);
			}
		}
	}
	
	return $feature;
}


sub _gtf_to_seqf_simple {

	my ($self, $fields) = @_;
	my $group = $fields->[8];
	my $feature = $self->_gff_to_seqf($fields);
	
	# extract essential tags
	my ($gene_id, $transcript_id);
	if ($group =~ /gene_id "([^"]+)";?/i) {
		$gene_id = $1;
	}
	if ($group =~ /transcript_id "([^"]+)";?/i) {
		$transcript_id = $1;
	}
	unless ($gene_id or $transcript_id) {
		# improperly formatted GTF file without one of these two items, nothing more to do
		return $feature;
	}
	
	# common subfeatures including exon, CDS, UTR, and codons
	if ($fields->[2] =~ /cds|exon|utr|codon|untranslated/i)  {
		$feature->add_tag_value('Parent', $transcript_id);
		
		# exon id if present
		if ($fields->[2] eq 'exon' and $group =~ /exon_id "([^"]+)";?/i) {
			$feature->primary_id($1);
		}
		
		# check gene parent
		if ($self->do_gene and not exists $self->{loaded}{$gene_id}) {
			my $gene = $self->_make_gene_parent($fields, $gene_id);
			$self->{loaded}{$gene_id} = $gene;
			push @{ $self->{top_features} }, $gene;
		}
		
		# check transcript parent
		unless (exists $self->{loaded}{$transcript_id}) {
			my $rna = $self->_make_rna_parent($fields, $transcript_id);
			$self->{loaded}{$transcript_id} = $rna;
			if ($self->do_gene) {
				$rna->add_tag_value('Parent', $gene_id);
				$self->{loaded}{$gene_id}->add_SeqFeature($rna);
			}
			else {
				push @{ $self->{top_features} }, $rna;
			}
		}
	}
	
	# a transcript feature
	elsif ($fields->[2] =~ /rna|transcript/i) {
		# these are sometimes present in GTF files, such as from Ensembl
		# but are not required and often absent
		$feature->primary_id($transcript_id);
		
		# transcript information
		if ($group =~ /transcript_name "([^"]+)";?/i) {
			$feature->display_name($1);
		}
		
		# check whether parent was loaded and add gene information if not
		if ($self->do_gene) {
			$feature->add_tag_value('Parent', $gene_id);
			if (not exists $self->{loaded}{$gene_id}) {
				my $gene = $self->_make_gene_parent($fields, $gene_id);
				$self->{loaded}{$gene_id} = $gene;
				push @{ $self->{top_features} }, $gene;
			}
		}
	}
	
	# a gene feature
	elsif ($fields->[2] eq 'gene') {
		# these are sometimes present in GTF files, such as from Ensembl
		# but are not required and often absent
		$feature->primary_id($gene_id);
		if ($group =~ /gene_name "([^"]+)";?/i) {
			$feature->display_name($1);
		}
	}
	
	# anything else, like CNS (conserved noncoding sequence) doesn't get any 
	# further special attributes, as far as I can tell
	return $feature;
}


sub _gtf_to_seqf_full {
	
	my ($self, $fields) = @_;
	my $group = $fields->[8];
	my $feature = $self->_gff_to_seqf($fields);
	
	# process the group tags
	my %attributes;
	foreach my $g (split('; ', $group)) { # supposed to be "; " as delimiter
		my ($tag, $value, @bits) = split ' ', $g;
		$value =~ s/[";]//g; # remove the flanking double quotes, assume no internal quotes
		$attributes{$tag} = $value;
	}
	my $transcript_id = $attributes{transcript_id} || undef;
	
	# assign special tags based on the feature type
	if ($fields->[2] =~ /cds|exon|utr|codon|untranslated/i) {
		$feature->add_tag_value('Parent', $attributes{'transcript_id'});
		
		# exon id if present
		if ($fields->[2] eq 'exon' and exists $attributes{'exon_id'}) {
			$feature->primary_id($attributes{'exon_id'});
		}
		
		# check gene parent
		if ($self->do_gene) {
			my $gene_id = $attributes{gene_id} || undef;
			if ($gene_id and not exists $self->{loaded}{$gene_id}) {
				my $gene = $self->_make_gene_parent($fields, $gene_id);
				$self->{loaded}{$gene_id} = $gene;
				push @{ $self->{top_features} }, $gene;
			}
		}
		
		# check transcript parent
		if ($transcript_id and not exists $self->{loaded}{$transcript_id}) {
			my $rna = $self->_make_rna_parent($fields, $transcript_id);
			$self->{loaded}{$transcript_id} = $rna;
			if ($self->do_gene) {
				$rna->add_tag_value('Parent', $attributes{gene_id});
				$self->{loaded}{ $attributes{gene_id} }->add_SeqFeature($rna);
			}
			else {
				push @{ $self->{top_features} }, $rna;
			}
		}
	}
	
	# transcript
	elsif ($fields->[2] =~ /transcript|rna/i) {
		# these are sometimes present in GTF files, such as from Ensembl
		
		# transcript information
		$feature->primary_id($transcript_id);
		delete $attributes{transcript_id};
		if (exists $attributes{transcript_name}) {
			$feature->display_name($attributes{transcript_name});
			delete $attributes{transcript_name};
		}
		
		# check gene parent
		if ($self->do_gene) {
			my $gene_id = $attributes{gene_id} || undef;
			if ($gene_id and not exists $self->{loaded}{$gene_id}) {
				my $gene = $self->_make_gene_parent($fields, $gene_id);
				$self->{loaded}{$gene_id} = $gene;
				push @{ $self->{top_features} }, $gene;
			}
			$feature->add_tag_value('Parent', $gene_id);
		}
	}
	
	# gene
	elsif ($fields->[2] eq 'gene') {
		# these are sometimes present in GTF files, such as from Ensembl
		# but are not required and often absent
		$feature->primary_id($attributes{gene_id});
		delete $attributes{gene_id};
		if (exists $attributes{gene_name}) {
			$feature->display_name($attributes{gene_name});
			delete $attributes{gene_name};
		}
	}
	
	# store remaining attributes, if any
	foreach my $key (keys %attributes) {
		$feature->add_tag_value($key, $attributes{$key});
	}
	return $feature;
}


sub _gff2_to_seqf {
	# generic gff1 or gff2 format or poorly defined gff3 or gtf file
	# hope for the best!
	my ($self, $fields) = @_;
	my $feature = $self->_gff_to_seqf($fields);
	
	# process groups
	# we have no uniform method of combining features, so we'll leave the tags 
	# as is and hope for the best
	foreach my $g (split(/\s*;\s*/, $fields->[8])) {
		my ($tag, $value) = split /\s+/, $g;
		next unless ($tag and $value);
		$feature->add_tag_value($tag, $value);
	}
	return $feature;
}


sub _gff_to_seqf {
	my ($self, $fields) = @_;
	
	# generate the basic SeqFeature
	my $feature = $self->{sfclass}->new(
		-seq_id         => $fields->[0],
		-source         => $fields->[1],
		-primary_tag    => $fields->[2],
		-start          => $fields->[3],
		-end            => $fields->[4],
		-strand         => $fields->[6],
	);
	
	# add more attributes if they're not null
	if ($fields->[5] ne '.') {
		$feature->score($fields->[5]);
	}
	if ($fields->[7] ne '.') {
		$feature->phase($fields->[7]);
	}
	
	# finished
	return $feature;
}


sub unescape {
  # Borrowed unashamedly from bioperl Bio::Tools::GFF
  # which in turn was borrowed from Bio::DB::GFF
  my $self = shift;
  my $v = shift;
  $v =~ tr/+/ /;
  $v =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $v;
}


sub _add_orphan {
	my ($self, $feature) = @_;
	push @{ $self->{orphans} }, $feature;
	return 1;
}


sub check_orphanage {
	my $self = shift;
	return unless scalar @{ $self->{orphans} };
	
	# go through the list of orphans
	my @reunited; # list of indices to delete after reuniting orphan with parent
	for (my $i = 0; $i < scalar @{ $self->{orphans} }; $i++) {
		my $orphan = $self->{orphans}->[$i];
		my $success = 0;
		
		# find the parent
		foreach my $parent ($orphan->get_tag_values('Parent') ) {
			if (exists $self->{loaded}{$parent}) {
				# we have loaded the parent
				# associate each orphan feature with the parent
				$self->{loaded}{$parent}->add_SeqFeature($orphan);
				$success++;
			}
		}
		# delete the orphan from the array if it found it's long lost parent
		push @reunited, $i if $success;
	}
	
	# clean up the orphanage
	while (@reunited) {
		my $i = pop @reunited;
		splice(@{ $self->{orphans} }, $i, 1);
	}
}

sub orphans {
	my $self = shift;
	my @orphans;
	foreach (@{ $self->{orphans} }) {
		push @orphans, $_;
	}
	return wantarray ? @orphans : \@orphans;
}


sub comments {
	my $self = shift;
	my @comments;
	foreach (@{ $self->{comments} }) {
		push @comments, $_;
	}
	return wantarray ? @comments : \@comments;
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

