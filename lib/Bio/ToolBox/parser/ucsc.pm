package Bio::ToolBox::parser::ucsc;
our $VERSION = '1.65';

=head1 NAME

Bio::ToolBox::parser::ucsc - Parser for UCSC genePred, refFlat, etc formats

=head1 SYNOPSIS

  use Bio::ToolBox::parser::ucsc;
  
  ### A simple transcript parser
  my $ucsc = Bio::ToolBox::parser::ucsc->new('file.genePred');
  
  ### A full fledged gene parser
  my $ucsc = Bio::ToolBox::parser::ucsc->new(
        file      => 'ensGene.genePred',
        do_gene   => 1,
        do_cds    => 1,
        do_utr    => 1,
        ensname   => 'ensemblToGene.txt',
        enssrc    => 'ensemblSource.txt',
  );
  
  ### Retrieve one transcript line at a time
  my $transcript = $ucsc->next_feature;
  
  ### Retrieve one assembled gene at a time
  my $gene = $ucsc->next_top_feature;
  
  ### Retrieve array of all assembled genes
  my @genes = $ucsc->top_features;
  
  # Each gene or transcript is a SeqFeatureI compatible object
  printf "gene %s is located at %s:%s-%s\n", 
    $gene->display_name, $gene->seq_id, 
    $gene->start, $gene->end;
  
  # Multiple transcripts can be assembled into a gene
  foreach my $transcript ($gene->get_SeqFeatures) {
    # each transcript has exons
    foreach my $exon ($transcript->get_SeqFeatures) {
      printf "exon is %sbp long\n", $exon->length;
    }
  }
  
  # Features can be printed in GFF3 format
  $gene->version(3);
  print STDOUT $gene->gff_string(1); 
   # the 1 indicates to recurse through all subfeatures
  

=head1 DESCRIPTION

This is a parser for converting UCSC-style gene prediction flat file formats into 
BioPerl-style L<Bio::SeqFeatureI> compliant objects, complete with nested objects 
representing transcripts, exons, CDS, UTRs, start- and stop-codons. Full control 
is available on what to parse, e.g. exons on, CDS and codons off. Additional gene 
information can be added by supplying additional tables of information, such as 
common gene names and descriptions, available from the UCSC repository. 

=head2 Table formats supported

Supported files are tab-delimited text files obtained from UCSC and described 
at L<http://genome.ucsc.edu/FAQ/FAQformat.html#format9>. Formats are identified 
by the number of columns, rather than specific file extensions, column name 
headers, or other metadata. Therefore, unmodified tables should only be used 
for correct parsing. Some errors are reported for incorrect lines. Unadulterated 
files can safely be downloaded from L<http://hgdownload.soe.ucsc.edu/downloads.html>.
Files obtained from the UCSC Table Browser can also be used with caution. Files 
may be gzip compressed.

File formats supported include the following.

=over 4

=item * Gene Prediction (genePred), 10 columns

=item * Gene Prediction with RefSeq gene Name (refFlat), 11 columns

=item * Extended Gene Prediction (genePredExt), 15 columns

=item * Extended Gene Prediction with bin (genePredExt), 16 columns

=item * knownGene table, 12 columns

=back

=head2 Supplemental information

The UCSC gene prediction tables include essential information, but not detailed 
information, such as common gene names, description, protein accession IDs, etc. 
This additional information can be associated with the genes or transcripts during 
parsing if the appropriate tables are supplied. These tables can be obtained from 
the UCSC download site L<http://hgdownload.soe.ucsc.edu/downloads.html>.

Supported tables include the following.

=over 4 

=item * refSeqStatus, for refGene, knownGene, and xenoRefGene tables

=item * refSeqSummary, for refGene, knownGene, and xenoRefGene tables

=item * ensemblToGeneName, for ensGene tables

=item * ensemblSource, for ensGene tables

=item * kgXref, for knownGene tables

=back

=head2 Implementation

For an implementation of this module to generate GFF3 formatted files from UCSC 
data sources, see the L<Bio::ToolBox> script L<ucsc_table2gff3.pl>.

=head1 METHODS

=head2 Initalize the parser object

=over 4

=item new

Initiate a UCSC table parser object. Pass a single value (a table file name) 
to open a table and parse its objects. Alternatively, pass an array of key 
value pairs to control how the table is parsed. Options include the following.

=over 4

=item file

=item table

Provide a file name for a UCSC gene prediction table. The file may be gzip 
compressed. 

=item source

Pass a string to be added as the source tag value of the SeqFeature objects. 
The default value is 'UCSC'. If the file name has a recognizable name, 
such as 'refGene' or 'ensGene', it will be used instead.

=item do_gene

Pass a boolean (1 or 0) value to combine multiple transcripts with the same gene 
name under a single gene object. Default is true.

-item do_exon

=item do_cds

=item do_utr

=item do_codon

Pass a boolean (1 or 0) value to parse certain subfeatures, including exon, 
CDS, five_prime_UTR, three_prime_UTR, stop_codon, and start_codon features. 
Default is false.

=item do_name

Pass a boolean (1 or 0) value to assign names to subfeatures, including exons, 
CDSs, UTRs, and start and stop codons. Default is false.

=item share

Pass a boolean (1 or 0) value to recycle shared subfeatures (exons and UTRs) 
between multiple transcripts of the same gene. This results in reduced 
memory usage, and smaller exported GFF3 files. Default is true. 

=item refseqsum

=item refseqstat

=item kgxref

=item ensembltogene

=item ensemblsource

Pass the appropriate file name for additional information.

=item class

Pass the name of a L<Bio::SeqFeatureI> compliant class that will be used to 
create the SeqFeature objects. The default is to use L<Bio::ToolBox::SeqFeature>.

=back

=back

=head2 Modify the parser object

These methods set or retrieve parameters, and load supplemental files and 
new tables.

=over 4

=item source

=item do_gene

=item do_exon

=item do_cds

=item do_utr

=item do_codon

=item do_name

=item share

These methods retrieve or set parameters to the parsing engine, same as 
the options to the new method.

=item fh

Set or retrieve the file handle of the current table. This module uses 
L<IO::Handle> objects. Be careful manipulating file handles of open tables!

=item open_file

Pass the name of a new table to parse. Existing gene models loaded in 
memory, if any, are discarded. Counts are reset to 0. Supplemental 
tables are not discarded.

=item load_extra_data($file, $type)

	my $file = 'hg19_refSeqSummary.txt.gz';
	my success = $ucsc->load_extra_data($file, 'summary');

Pass two values, the file name of the supplemental file and the type 
of supplemental data. Values can include the following 

=over 4

=item * refseqstatus or status

=item * refseqsummary or summary

=item * kgxref

=item * ensembltogene or ensname

=item * ensemblsource or enssrc

=back

The number of transcripts with information loaded from the supplemental 
data file is returned.

=back

=head2 Feature retrieval

The following methods parse the table lines into SeqFeature objects. 
It is best if methods are not mixed; unexpected results may occur. 

=over 4

=item next_feature

This will read the next line of the table and parse it into a gene or 
transcript object. However, multiple transcripts from the same gene are 
not assembled together under the same gene object. 

=item next_top_feature

This method will return all top features (typically genes), with multiple 
transcripts of the same gene assembled under the same gene object. Transcripts 
are assembled together if they share the same gene name and the transcripts 
overlap. If transcripts share the same gene name but do not overlap, they 
are placed into separate gene objects with the same name but different 
C<primary_id> tags. Calling this method will parse the entire table into 
memory (so that multiple transcripts may be assembled), but only one object 
is returned at a time. Call this method repeatedly using a while loop to 
get all features.

=item top_features

This method is similar to L</next_top_feature>, but instead returns an array 
of all the top features. 

=back

=head2 Other methods

Additional methods for working with the parser object and the parsed 
SeqFeature objects.

=over 4

=item parse_table

Parses the table into memory. If a table wasn't provided using the 
L</new> or L</open_file> methods, then a filename can be passed to this 
method and it will automatically be opened for you. 

=item find_gene

	my $gene = $ucsc->find_gene(
		display_name => 'ABC1',
		primary_id   => 'gene000001',
	);

Pass a gene name, or an array of key =E<gt> values (name, display_name, 
ID, primary_ID, and/or coordinate information), that can be used 
to find a gene already loaded into memory. Only really successful if the 
entire table is loaded into memory. Genes with a matching name are 
confirmed by a matching ID or overlapping coordinates, if available. 
Otherwise the first match is returned.

=item counts

This method will return a hash of the number of genes and RNA types that 
have been parsed.

=item typelist

This method will return a comma-delimited list of the feature types or 
C<primary_tag>s found in the parsed file. Returns a generic list if a 
file has not been parsed.

=item from_ucsc_string

A bare bones method that will convert a tab-delimited text line from a UCSC 
formatted gene table into a SeqFeature object for you. Don't expect alternate 
transcripts to be assembled into genes. 

=item seq_ids

Returns an array or array reference of the names of the chromosomes or 
reference sequences present in the table.

=item seq_id_lengths

Returns a hash reference to the chromosomes or reference sequences and 
their corresponding lengths. In this case, the length is inferred by the 
greatest gene end position.

=back

=head2 Bio::ToolBox::parser::ucsc::builder

This is a private module that is responsible for building SeqFeature 
objects from UCSC table lines. It is not intended for general public use.

=head1 SEE ALSO

L<Bio::ToolBox::SeqFeature>, L<Bio::ToolBox::parser::gff>

=cut

use strict;
use Carp qw(carp cluck croak);
use Bio::ToolBox::Data::file; # only used to get an open filehandle

1;

sub new {
	my $class = shift;
	my $self = {
		'fh'            => undef,
		'version'       => undef,
		'source'        => 'UCSC',
		'top_features'  => [],
		'gene2seqf'     => {},
		'id2count'      => {},
		'counts'        => {},
		'do_gene'       => 1, 
		'do_exon'       => 0,
		'do_cds'        => 0, 
		'do_utr'        => 0, 
		'do_codon'      => 0,
		'do_name'       => 0,
		'share'         => 1, 
		'refseqsum'     => {}, 
		'refseqstat'    => {},
		'kgxref'        => {},
		'ensembldata'   => {},
		'eof'           => 0,
		'line_count'    => 0,
		'sfclass'       => 'Bio::ToolBox::SeqFeature', # default class
		'seq_ids'       => {}, 
	};
	bless $self, $class;
	
	# check for options
	if (@_) {
		if (scalar(@_) == 1) {
			# short and sweet, just a file, we assume
			my $file = shift @_;
			$self->open_file($file);
		}
		else {
			my %options = @_;
			if (exists $options{file} or $options{table}) {
				$options{file} ||= $options{table};
				$self->open_file( $options{file} );
			}
			if (exists $options{do_gene}) {
				$self->do_gene($options{do_gene});
			}
			if (exists $options{do_exon}) {
				$self->do_exon($options{do_exon});
			}
			if (exists $options{do_cds}) {
				$self->do_cds($options{do_cds});
			}
			if (exists $options{do_utr}) {
				$self->do_utr($options{do_utr});
			}
			if (exists $options{do_codon}) {
				$self->do_codon($options{do_codon});
			}
			if (exists $options{do_name}) {
				$self->do_name($options{do_name});
			}
			if (exists $options{share}) {
				$self->share($options{share});
			}
			if (exists $options{source}) {
				$self->source($options{source});
			}
			if (exists $options{refseqsum}) {
				$self->load_extra_data($options{refseqsum}, 'refseqsum');
			}
			elsif (exists $options{summary}) {
				$self->load_extra_data($options{summary}, 'refseqsum');
			}
			if (exists $options{refseqstat}) {
				$self->load_extra_data($options{refseqstat}, 'refseqstat');
			}
			elsif (exists $options{status}) {
				$self->load_extra_data($options{status}, 'refseqstat');
			}
			if (exists $options{kgxref}) {
				$self->load_extra_data($options{kgxref}, 'kgxref');
			}
			if (exists $options{ensembltogenename}) {
				$self->load_extra_data($options{ensembltogenename}, 'ensembltogene');
			}
			elsif (exists $options{ensname}) {
				$self->load_extra_data($options{ensname}, 'ensembltogene');
			}
			if (exists $options{ensemblsource}) {
				$self->load_extra_data($options{ensemblsource}, 'ensemblsource');
			}
			elsif (exists $options{enssrc}) {
				$self->load_extra_data($options{enssrc}, 'ensemblsource');
			}
			if (exists $options{class}) {
				$self->{sfclass} = $options{class};
			}
		}
	}
	
	# done
	return $self;
}

sub version {
	return shift->{version};
}

sub source {
	my $self = shift;
	if (@_) {
		$self->{'source'} = shift;
	}
	return $self->{'source'};
}

sub simplify {
	# this doesn't do anything, for now, but maintain compatibility with gff parser
	return 0;
}

sub do_gene {
	my $self = shift;
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
	if (@_) {
		$self->{'do_name'} = shift;
	}
	return $self->{'do_name'};
}	

sub share {
	my $self = shift;
	if (@_) {
		$self->{'share'} = shift;
	}
	return $self->{'share'};
}	

sub fh {
	my $self = shift;
	if (@_) {
		$self->{'fh'} = shift;
	}
	return $self->{'fh'};
}

sub open_file {
	my $self = shift;
	
	# check file
	my $filename = shift;
	unless ($filename) {
		cluck("no file name passed!");
		return;
	}
	
	# Open filehandle object 
	my $fh = Bio::ToolBox::Data::file->open_to_read_fh($filename) or
		croak " cannot open file '$filename'!\n";
	
	# check number of columns
	my $ncol;
	while (my $line = $fh->getline) {
		next if $line =~ /^#/;
		next unless $line =~ /\w+/;
		my @fields = split "\t", $line;
		$ncol = scalar(@fields);
		last;
	}
	unless ($ncol == 16 or $ncol == 15 or $ncol == 12 or $ncol == 11 or $ncol == 10) {
		carp " File '$filename' doesn't have recognizable number of columns! It has $ncol.";
		return;
	}
	if ($ncol == 10) {
		# turn off gene processing for simple genePred which has no gene names
		$self->do_gene(0);
	}
	
	# reopen file handle
	$fh->close;
	$fh = Bio::ToolBox::Data::file->open_to_read_fh($filename);
	
	# reset source as necessary
	if ($filename =~ /ensgene/i and $self->source eq 'UCSC') {
		$self->source('EnsGene');
	}
	elsif ($filename =~ /xenorefgene/i and $self->source eq 'UCSC') {
		$self->source('xenoRefGene');
	}
	elsif ($filename =~ /refgene/i and $self->source eq 'UCSC') {
		$self->source('refGene');
	}
	elsif ($filename =~ /refseq/i and $self->source eq 'UCSC') {
		$self->source('refSeq');
	}
	elsif ($filename =~ /knowngene/i and $self->source eq 'UCSC') {
		$self->source('knownGene');
	}
	
	# check existing data
	# the reason this parser can handle reading a second file without having to make a 
	# new parser object (like gff and bed) is so that we can potentially recycle the 
	# extra UCSC information 
	if ($self->fh) {
		# close existing
		$self->fh->close;
		# go ahead and clear out existing data
		$self->{'version'}       = undef;
		$self->{'top_features'}  = [];
		$self->{'duplicate_ids'} = {};
		$self->{'gene2seqf'}     = {};
		$self->{'id2count'}      = {};
		$self->{'counts'}        = {};
		$self->{'eof'}           = 0;
		$self->{'line_count'}    = 0;
	}
	
	$self->fh($fh);
	return 1;
}

sub load_extra_data {
	my ($self, $file, $type) = @_;
	unless ($file) {
		cluck "no file name passed!";
		return;
	}
	
	# check the type 
	if ($type =~ /ensembltogene|ensname/i) {
		$type = 'ensembltogene';
	}
	elsif ($type =~ /ensemblsource|enssrc/i) {
		$type = 'ensemblsource';
	}
	elsif ($type =~ /refseqstat|status/i) {
		$type = 'refseqstat';
	}
	elsif ($type =~ /refseqsum|summary/i) {
		$type = 'refseqsum';
	}
	elsif ($type =~ /kgxref/i) {
		$type = 'kgxref';
	}
	else {
		carp "unknown type '$type' to load extra data";
		return;
	}
	
	my $fh = Bio::ToolBox::Data::file->open_to_read_fh($file);
	unless ($fh) {
		carp "unable to open file '$file'! $!";
		return;
	}
	
	# load ensembl data
	my $count = 0;
	if ($type =~ /ensembl/) {
		# we will store gene name in position 0, and source in position 1
		my $index = $type eq 'ensembltogene' ? 0 : 1;
		while (my $line = $fh->getline) {
			
			# process line
			chomp $line;
			next if ($line =~ /^#/);
			my @line_data = split /\t/, $line;
			if (scalar @line_data != 2) {
				carp " file $file doesn't seem right!? Line has " .
					scalar @line_data . " elements!\n";
				return;
			}
			
			# store data into hash
			$self->{'ensembldata'}{ $line_data[0] }->[$index] = $line_data[1];
			$count++;
		}
	}
	
	# load various refSeq data
	else {
		# we just store the line data based on the gene ID
		# each table has different elements
		# here they are for reference
		
		### refSeqStatus table
		# 0	mrnaAcc	RefSeq gene accession name
		# 1	status	Status ('Unknown', 'Reviewed', 'Validated', 'Provisional', 'Predicted', 'Inferred')
		# 2	molecule type ('DNA', 'RNA', 'ds-RNA', 'ds-mRNA', 'ds-rRNA', 'mRNA', 'ms-DNA', 'ms-RNA', 'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'ss-DNA', 'ss-RNA', 'ss-snoRNA', 'tRNA', 'cRNA', 'ss-cRNA', 'ds-cRNA', 'ms-rRNA')	values	molecule type
	
		### refSeqSummary table
		# 0	RefSeq mRNA accession
		# 1	completeness	FullLength ('Unknown', 'Complete5End', 'Complete3End', 'FullLength', 'IncompleteBothEnds', 'Incomplete5End', 'Incomplete3End', 'Partial')	
		# 1	summary	 	text	values	Summary comments
	
		### kgXref table
		# 0	kgID	Known Gene ID
		# 1	mRNA	mRNA ID
		# 2	spID	SWISS-PROT protein Accession number
		# 3	spDisplayID	 SWISS-PROT display ID
		# 4	geneSymbol	Gene Symbol
		# 5	refseq	 RefSeq ID
		# 6	protAcc	 NCBI protein Accession number
		# 7	description	Description
	
		### ensemblToGeneName table
		# 0 Ensembl transcript ID
		# 1 gene name
		
		# load the table
		while (my $line = $fh->getline) {
			chomp $line;
			next if ($line =~ /^#/);
			my @line_data = split /\t/, $line;
	
			# the unique id should be the first element in the array
			# take it off the array, since it doesn't need to be stored there too
			my $id = shift @line_data;
	
			# check for duplicate lines
			if (exists $self->{$type}{$id} ) {
				warn "  $type line for identifier $id exists twice!\n";
				next;
			}
	
			# store data into hash
			$self->{$type}{$id} = [@line_data];
			$count++;
		}
	}
	
	$fh->close;
	return $count;
}

sub typelist {
	my $self = shift;
	my @items;
	foreach my $k (keys %{$self->{counts}}) {
		push @items, $k if $self->{counts}{$k} > 0;
	}
	if (@items) {
		return join(',', @items);
	}
	else {
		# return generic list
		return $self->do_gene ? 'gene,mRNA,ncRNA,exon,CDS' : 'mRNA,ncRNA,exon,CDS';
	}
}

sub next_feature {
	my $self = shift;
	
	# check that we have an open filehandle
	unless ($self->fh) {
		croak("no UCSC file loaded to parse!");
	}
	
	while (my $line = $self->fh->getline) {
		chomp $line;
		if ($line =~ /^#/ or $line !~ /\w+/) {
			$self->{line_count}++;
			next;
		}
		my $builder = Bio::ToolBox::parser::ucsc::builder->new($line, $self);
		$self->{line_count}++;
		unless ($builder) {
			# builder will print its own error message if fails
			warn " unable to parse line number ", $self->{line_count}, "\n";
			next;
		}
		
		# generate the feature from the line
		my $feature;
		if ($self->do_gene) {
			$feature = $builder->build_gene;
		}
		else {
			$feature = $builder->build_transcript;
		}
		
		# return the object, we do this while loop once per valid line
		return $feature;
	}
	
	# presumed end of file
	$self->{'eof'} = 1;
	return;
}

sub next_top_feature {
	my $self = shift;
	unless ($self->{'eof'}) {
		$self->parse_table;
	}
	return shift @{ $self->{top_features} };
}

sub top_features {
	my $self = shift;
	unless ($self->{'eof'}) {
		$self->parse_table;
	}
	my @features = @{ $self->{top_features} };
	return wantarray ? @features : \@features;
}

*parse_file = \&parse_table;

sub parse_table {
	my $self = shift;
	if (@_) {
		$self->open_file(shift) or return;
	}
	unless ($self->fh) {
		carp "must open a file first!";
		return;
	}
	return if ($self->{'eof'});
	
	#### Main Loop
	print "  Parsing UCSC gene table....\n";
	while (my $feature = $self->next_feature) {
		
		# add to gene2seqf hash
		my $gene = $self->find_gene($feature); 
		unless ($gene) {
			# the current feature is not in the hash, so add it
			$self->{gene2seqf}->{ lc $feature->display_name } = [ $feature ];
		}
		
		# check chromosome
		my $s = $feature->seq_id;
		unless (exists $self->{seq_ids}{$s}) {
			$self->{seq_ids}{$s} = $feature->end;
		}
		$self->{seq_ids}{$s} = $feature->end if $feature->end > $self->{seq_ids}{$s};
	}
	
	# add to the top list of features, Schwartzian transform and sort
	# based on the genes found in the gene2seqf hash
	push @{ $self->{top_features} }, 
		map {$_->[2]}
		sort {$a->[0] cmp $b->[0] or $a->[1] <=> $b->[1]}
		map [$_->seq_id, $_->start, $_],
		map @{ $self->{gene2seqf}->{$_} },
		keys %{$self->{gene2seqf}};
	
	return 1;
}

sub find_gene {
	my $self = shift;
	
	# get the name and coordinates from arguments
	my ($name, $id, $chrom, $start, $end, $strand);
	if (scalar @_ == 0) {
		carp "must provide information to find_gene method!";
		return;
	}
	elsif (scalar @_ == 1) {
		$name = $_[0];
	}
	else {
		my %opt = @_;
		$name  = $opt{name} || $opt{display_name} || undef;
		$id    = $opt{id} || $opt{primary_id} || undef;
		$chrom = $opt{chrom} || $opt{seq_id} || undef;
		$start = $opt{start} || undef;
		$end   = $opt{stop} || $opt{end} || undef;
		$strand = $opt{strand} || 0;
	}
	unless ($name) {
		carp "name is required for find_gene!";
		return;
	}
	
	# check if a gene with this name exists
	if (exists $self->{gene2seqf}->{lc $name} ) {
		# we found a matching gene
		
		# pull out the gene seqfeature(s) array reference
		# there may be more than one gene
		my $genes = $self->{gene2seqf}->{ lc $name };
		
		# go through a series of checks to find the appropriate 
		if ($id) {
			foreach my $g (@$genes) {
				if ($g->primary_id eq $id) {
					return $g;
				}
			}
			return; # none of these matched despite having an ID
		}
		if ($chrom and $start and $end) {
			foreach my $g (@$genes) {
				if ( 
					# overlap method borrowed from Bio::RangeI
					($g->strand == $strand) and not (
						$g->start > $end or 
						$g->end < $start
					)
				) {
					# gene and transcript overlap on the same strand
					# we found the intersecting gene
					return $g;
				}
			}
			return; # none of these matched despite having coordinate info
		}
		if (scalar @$genes == 1) {
			# going on trust here that this is the one
			return $genes->[0];
		}
		elsif (scalar @$genes > 1) {
			carp "more than one gene named $name found!";
			return $genes->[0];
		}
		
		# nothing suitable found
		return;
	}
}

sub counts {
	my $self = shift;
	my %counts = %{ $self->{counts} };
	return wantarray ? %counts : \%counts;
}

sub from_ucsc_string {
	my ($self, $string) = @_;
	return unless $string;
	my $builder = Bio::ToolBox::parser::ucsc::builder->new($string, $self);
	return unless $builder;
	if ($self->do_gene) {
		return $builder->build_gene;
	}
	else {
		return $builder->build_transcript;
	}
}

sub seq_ids {
	my $self = shift;
	my @s = keys %{$self->{seq_ids}};
	return wantarray ? @s : \@s;
}

sub seq_id_lengths {
	my $self = shift;
	return $self->{seq_ids};
}



package Bio::ToolBox::parser::ucsc::builder;
use strict;
use Carp qw(carp cluck croak);
our $SFCLASS = ''; # SeqFeature class to use

1;

sub new {
	my ($class, $line, $ucsc) = @_;
	my %self;
	my $format;
	
	# check SeqFeature class
	if ($ucsc->{sfclass} ne $SFCLASS) {
		$SFCLASS = $ucsc->{sfclass};
		eval "require $SFCLASS" or croak $@;
	}
	
	chomp $line;
	my @linedata = split /\t/, $line;
	
	# we're identifying the type of table based on the number of columns
	# may not be the best or accurate method, but it generally works for valid formats
	
	### Extended Gene Prediction Table ###
	if (scalar @linedata == 16) {
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
		
		$format            = 'genePredExt with bin';
		$self{name}        = $linedata[1];
		$self{chrom}       = $linedata[2];
		$self{strand}      = $linedata[3];
		$self{txStart}     = $linedata[4] + 1;
		$self{txEnd}       = $linedata[5];
		$self{cdsStart}    = $linedata[6] + 1;
		$self{cdsEnd}      = $linedata[7];
		$self{exonCount}   = $linedata[8];
		$self{exonStarts}  = $linedata[9];
		$self{exonEnds}    = $linedata[10];
		$self{name2}       = $linedata[12] || undef;
		$self{gene_name}   = $ucsc->{ensembldata}->{ $linedata[1] }->[0] ||
		                     $linedata[12] || undef;
		$self{note}        = $ucsc->{refseqsum}->{ $linedata[1] }->[1] || undef;
		$self{status}      = $ucsc->{refseqstat}->{ $linedata[1] }->[0] || undef;
		$self{completeness} = $ucsc->{refseqsum}->{ $linedata[1] }->[0] || undef;
		if ($linedata[1] =~ /^N[MR]_\d+/) {
			$self{refseq} = $linedata[1];
		}
	}
	### Extended Gene Prediction Table ###
	elsif (scalar @linedata == 15) {
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
		
		$format            = 'genePredExt';
		$self{name}        = $linedata[0];
		$self{chrom}       = $linedata[1];
		$self{strand}      = $linedata[2];
		$self{txStart}     = $linedata[3] + 1;
		$self{txEnd}       = $linedata[4];
		$self{cdsStart}    = $linedata[5] + 1;
		$self{cdsEnd}      = $linedata[6];
		$self{exonCount}   = $linedata[7];
		$self{exonStarts}  = $linedata[8];
		$self{exonEnds}    = $linedata[9];
		$self{name2}       = $linedata[11] || undef;
		$self{gene_name}   = $ucsc->{ensembldata}->{ $linedata[0] }->[0] ||
		                     $linedata[11] || undef;
		$self{note}        = $ucsc->{refseqsum}->{ $linedata[0] }->[1] || undef;
		$self{status}      = $ucsc->{refseqstat}->{ $linedata[0] }->[0] || undef;
		$self{completeness} = $ucsc->{refseqsum}->{ $linedata[0] }->[0] || undef;
		if ($linedata[0] =~ /^N[MR]_\d+/) {
			$self{refseq} = $linedata[0];
		}
	}
	### Known Gene Table ###
	elsif (scalar @linedata == 12) {
		
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
		
		$format            = 'knownGene';
		$self{name}        = $ucsc->{kgxref}->{ $linedata[0] }->[0] || $linedata[0];
		$self{chrom}       = $linedata[1];
		$self{strand}      = $linedata[2];
		$self{txStart}     = $linedata[3] + 1;
		$self{txEnd}       = $linedata[4];
		$self{cdsStart}    = $linedata[5] + 1;
		$self{cdsEnd}      = $linedata[6];
		$self{exonCount}   = $linedata[7];
		$self{exonStarts}  = $linedata[8];
		$self{exonEnds}    = $linedata[9];
		$self{name2}       = $linedata[0];
		$self{gene_name}   = $ucsc->{kgxref}->{ $linedata[0] }->[3] || # geneSymbol
							$ucsc->{kgxref}->{ $linedata[0] }->[0] || # mRNA id
							$ucsc->{kgxref}->{ $linedata[0] }->[4] || # refSeq id
							$linedata[0]; # ugly default
		$self{note}        = $ucsc->{kgxref}->{ $linedata[0] }->[6] || undef;
		$self{refseq}      = $ucsc->{kgxref}->{ $linedata[0] }->[4] || undef;
		$self{status}      = $ucsc->{refseqstat}->{ $self{refseq} }->[0] || undef;
		$self{completeness} = $ucsc->{refseqsum}->{ $self{refseq} }->[0] || undef;
		$self{spid}        = $ucsc->{kgxref}->{ $linedata[0] }->[1] || undef; # SwissProt ID
		$self{spdid}       = $ucsc->{kgxref}->{ $linedata[0] }->[2] || undef; # SwissProt display ID
		$self{protacc}     = $ucsc->{kgxref}->{ $linedata[0] }->[5] || undef; # NCBI protein accession
	}
	### refFlat or Gene Prediction Table ###
	elsif (scalar @linedata == 11) {
		
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
		
		$format            = 'refFlat';
		$self{gene_name}   = $ucsc->{ensembldata}->{ $linedata[1] }->[0] ||
		                     $linedata[0] || undef;
		$self{name2}       = $linedata[0];
		$self{name}        = $linedata[1];
		$self{chrom}       = $linedata[2];
		$self{strand}      = $linedata[3];
		$self{txStart}     = $linedata[4] + 1;
		$self{txEnd}       = $linedata[5];
		$self{cdsStart}    = $linedata[6] + 1;
		$self{cdsEnd}      = $linedata[7];
		$self{exonCount}   = $linedata[8];
		$self{exonStarts}  = $linedata[9];
		$self{exonEnds}    = $linedata[10];
		$self{note}        = $ucsc->{refseqsum}->{ $linedata[1] }->[1] || undef;
		$self{status}      = $ucsc->{refseqstat}->{ $linedata[1] }->[0] || undef;
		$self{completeness} = $ucsc->{refseqsum}->{ $linedata[1] }->[0] || undef;
		if ($linedata[1] =~ /^N[MR]_\d+/) {
			$self{refseq} = $linedata[1];
		}
	}
	### Gene Prediction Table ###
	elsif (scalar @linedata == 10) {
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
		
		$format            = 'genePred';
		$self{name}        = $linedata[0];
		$self{chrom}       = $linedata[1];
		$self{strand}      = $linedata[2];
		$self{txStart}     = $linedata[3] + 1;
		$self{txEnd}       = $linedata[4];
		$self{cdsStart}    = $linedata[5] + 1;
		$self{cdsEnd}      = $linedata[6];
		$self{exonCount}   = $linedata[7];
		$self{exonStarts}  = $linedata[8];
		$self{exonEnds}    = $linedata[9];
		$self{name2}       = $linedata[0]; # re-use transcript name
		$self{gene_name}   = $linedata[0]; # re-use transcript name
		$self{note}        = $ucsc->{refseqsum}->{ $linedata[0] }->[1] || undef;
		$self{status}      = $ucsc->{refseqstat}->{ $linedata[0] }->[0] || undef;
		$self{completeness} = $ucsc->{refseqsum}->{ $linedata[0] }->[0] || undef;
		if ($linedata[0] =~ /^N[MR]_\d+/) {
			$self{refseq} = $linedata[0];
		}
	}
	else {
		# unrecognized line format
		carp "unrecognized format, line has " . scalar @linedata . "elements";
		return;
	}
	
	# verify
	my @errors; 
	push @errors, 'strand'     unless $self{strand}     =~ /^[\+\-]$/;
	push @errors, 'txStart'    unless $self{txStart}    =~ /^\d+$/;
	push @errors, 'txEnd'      unless $self{txEnd}      =~ /^\d+$/;
	push @errors, 'cdsStart'   unless $self{cdsStart}   =~ /^\d+$/;
	push @errors, 'cdsEnd'     unless $self{cdsEnd}     =~ /^\d+$/;
	push @errors, 'exonCount'  unless $self{exonCount}  =~ /^\d+$/;
	push @errors, 'exonStarts' unless $self{exonStarts} =~ /^[\d,]+$/;
	push @errors, 'exonEnds'   unless $self{exonEnds}   =~ /^[\d,]+$/;
	if (@errors) {
		warn "line format for $format has the following errors: @errors";
		return;
	}
	
	# fix values
	$self{strand}     = $self{strand} eq '+' ? 1 : -1;
	$self{exonStarts} = [ map {$_ += 1} ( split ",", $self{exonStarts} ) ];
	$self{exonEnds}   = [ ( split ",", $self{exonEnds} ) ];
	
	# Attempt to identify the transcript type
	my $type = $ucsc->{ensembldata}->{ $self{name} }->[1] || undef;
		# check if we have loaded ensembl source data and use that if available
	if ( $self{cdsStart} - 1 == $self{cdsEnd} ) {
		# there appears to be no coding potential when 
		# txEnd = cdsStart = cdsEnd
		# if you'll look, all of the exon phases should also be -1
		
		if ($type) {
			# we have an ensembl source type, so prefer to use that
			$self{type} = $type;
		}
		
		# otherwise, we may be able to infer some certain 
		# types from the gene name
		elsif ($self{name2} =~ /^mir/i) {
			# a noncoding gene whose name begins with mir is likely a micro RNA
			$self{type} = 'miRNA';
		}
		elsif ($self{name2} =~ /^snr/i) {
			# a noncoding gene whose name begins with snr is likely a snRNA
			$self{type} = 'snRNA';
		}
		elsif ($self{name2} =~ /^sno/i) {
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
		$self{type} = defined $type ? $type : 'mRNA';
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
	return exists $self->{refseq} ? $self->{refseq} : undef;
}

sub note {
	my $self = shift;
	return exists $self->{note} ? $self->{note} : undef;
}

sub status {
	my $self = shift;
	return exists $self->{status} ? $self->{status} : undef;
}

sub completeness {
	my $self = shift;
	return exists $self->{completeness} ? $self->{completeness} : undef;
}

sub ucsc {
	return shift->{ucsc};
}

sub build_gene {
	my $self = shift;
	
	# shortcuts
	my $ucsc = $self->ucsc;
	my $id2count = $ucsc->{id2count};
	my $ensembldata = $ucsc->{ensembldata};
	
	# find a pre-existing gene to update or build a new one
	my $gene = $self->find_gene;
	if ($gene) {
		# update as necessary
		if ( ($self->txStart) < $gene->start) {
			# update the transcription start position
			$gene->start( $self->txStart );
		}
		if ($self->txEnd > $gene->end) {
			# update the transcription stop position
			$gene->end( $self->txEnd );
		}
	}
	else {
		# build a new gene
		$gene = $SFCLASS->new(
			-seq_id        => $self->chrom,
			-source        => $ucsc->source,
			-primary_tag   => 'gene',
			-start         => $self->txStart,
			-end           => $self->txEnd,
			-strand        => $self->strand,
			-phase         => '.',
			-display_name  => $self->gene_name,
		);
		$ucsc->{counts}->{gene} += 1;
	
		# Add a unique primary ID
		my $id = $self->name2;
		if (exists $id2count->{ lc $id }) {
			# we've encountered this gene ID before
		
			# then make name unique by appending the count number
			$id2count->{ lc $id } += 1;
			$id .= '.' . $id2count->{ lc $id };
		}
		else {
			# this is the first transcript with this id
			# set the id counter
			$id2count->{lc $id} = 0;
		}
		$gene->primary_id($id);
		
		# Add an alias
		if ($self->name2 ne $self->gene_name) {
			$gene->add_tag_value('Alias', $self->name2);
		}
	}
	
	# now build the transcript for the gene
	my $transcript = $self->build_transcript($gene);
	$gene->add_SeqFeature($transcript);
	$transcript->add_tag_value('Parent', $gene->primary_id);
	
	# update extra attributes as necessary
	$self->update_attributes($gene);
	
	# finished
	return $gene;
}

sub build_transcript {
	my ($self, $gene) = @_; # gene is not required
	
	# shortcuts
	my $ucsc = $self->ucsc;
	my $id2count = $ucsc->{id2count};
	my $ensembldata = $ucsc->{ensembldata};
	my $counts = $ucsc->{counts};
	
	# Uniqueify the transcript ID and name
	my $id = $self->name;
	if (exists $id2count->{ lc $id } ) {
		# we've encountered this transcript ID before
		
		# now need to make ID unique by appending a number
		$id2count->{ lc $id } += 1;
		$id .= '.' . $id2count->{ lc $id };
	}
	else {
		# this is the first transcript with this id
		$id2count->{lc $id} = 0;
	}
	
	# identify the primary_tag value
	my ($type, $biotype);
	if (exists $ensembldata->{$self->name}) {
		my $t = $ensembldata->{$self->name}->[1] || undef;
		if ($t and $t =~ /protein.coding/i) {
			$type = 'mRNA';
			$biotype = $t;
		}
		elsif ($t and $t =~ /rna|transcript/i) {
			$type = $t;
			$biotype = $t;
		}
		elsif ($t) {
			$type = 'transcript';
			$biotype = $t;
		}
		else {
			$type = $self->type;
		}
	}
	else {
		$type = $self->type;
	}
	
	# Generate the transcript SeqFeature object
	my $transcript = $SFCLASS->new(
		-seq_id        => $self->chrom,
		-source        => $ucsc->source,
		-primary_tag   => $type,
		-start         => $self->txStart,
		-end           => $self->txEnd,
		-strand        => $self->strand,
		-phase         => '.',
		-display_name  => $self->name,
		-primary_id    => $id,
	);
	
	# add gene name as an alias
	if ($self->gene_name ne $self->name2) {
		$transcript->add_tag_value('Alias', $self->gene_name);
	}
	
	# update extra attributes as necessary
	$self->update_attributes($transcript);
	
	# add transcript specific attributes
	if (defined $self->completeness ) {
		$transcript->add_tag_value( 'completeness', $self->completeness );
	}
	if (defined $self->status ) {
		$transcript->add_tag_value( 'status', $self->status );
	}
	if ($biotype) {
		$transcript->add_tag_value( 'biotype', $biotype);
	}
	
	# add the exons
	if ($ucsc->do_exon) {
		$self->add_exons($transcript, $gene);
	}
	
	# add CDS, UTRs, and codons if necessary
	if ( $self->cdsStart - 1 != $self->cdsEnd ) {
		
		if ($ucsc->do_utr) {
			$self->add_utrs($transcript, $gene);
		}
		
		if ($ucsc->do_codon) {
			$self->add_codons($transcript, $gene);
		}
		
		if ($ucsc->do_cds) {
			$self->add_cds($transcript);
		}
	}
	
	# record the type of transcript
	$counts->{$type} += 1;
	
	# transcript is complete
	return $transcript;
}

sub update_attributes {
	my ($self, $seqf) = @_;
	
	# add Note if possible
	if (defined $self->note ) {
		$self->add_unique_attribute($seqf, 'Note', $self->note );
	}
	
	# add refSeq identifier if possible
	if (defined $self->refseq) {
		$self->add_unique_attribute($seqf, 'Dbxref', 'RefSeq:' . $self->refseq);
	}
	
	# add SwissProt identifier if possible
	if (exists $self->{spid} and defined $self->{spid}) {
		$self->add_unique_attribute($seqf, 'Dbxref', 'Swiss-Prot:' . $self->{spid});
	}
	
	# add SwissProt display identifier if possible
	if (exists $self->{spdid} and defined $self->{spdid}) {
		$self->add_unique_attribute($seqf, 'swiss-prot_display_id', $self->{spdid});
	}
	
	# add NCBI protein access identifier if possible
	if (exists $self->{protacc} and defined $self->{protacc}) {
		$self->add_unique_attribute($seqf, 'Dbxref', 'RefSeq:' . $self->{protacc});
	}
}

sub add_unique_attribute {
	my ($self, $seqf, $tag, $value) = @_;
	
	# look for a pre-existing identical tag value
	my $check = 1;
	foreach ($seqf->get_tag_values($tag)) {
		if ($_ eq $value) {
			$check = 0;
			last;
		}
	}
	
	# add it if our value is unique
	$seqf->add_tag_value($tag, $value) if $check;
}

sub add_exons {
	my ($self, $transcript, $gene) = @_;
	my $ucsc = $self->ucsc;
	
	# Add the exons
	EXON_LOOP:
	for (my $i = 0; $i < $self->exonCount; $i++) {
		
		# first look for existing
		if ($ucsc->share and $gene) {
			my $exon = $self->find_existing_subfeature($gene, 'exon', 
				$self->exonStarts->[$i], $self->exonEnds->[$i]);
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
		if ($transcript->strand == 1) {
			# forward strand
			$number = $i;
		}
		else {
			# reverse strand
			$number = abs( $i - $self->exonCount + 1);
		}
		
		# build the exon seqfeature
		my $exon = $SFCLASS->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source,
			-primary_tag   => 'exon',
			-start         => $self->exonStarts->[$i],
			-end           => $self->exonEnds->[$i],
			-strand        => $transcript->strand,
			-primary_id    => $transcript->primary_id . ".exon$number",
		);
		
		# add name if requested
		if ($ucsc->do_name) {
			$exon->display_name( $transcript->display_name . ".exon$number" );
		}
		
		# associate with transcript
		$transcript->add_SeqFeature($exon);
	}
}

sub add_utrs {
	my ($self, $transcript, $gene) = @_;
	my $ucsc = $self->ucsc;
	
	# we will scan each exon and look for a potential utr and build it
	my @utrs;
	UTR_LOOP:
	for (my $i = 0; $i < $self->exonCount; $i++) {
		
		# transform index for reverse strands
		# this will allow numbering from 5'->3'
		my $number; 
		if ($transcript->strand == 1) {
			# forward strand
			$number = $i;
		}
		else {
			# reverse strand
			$number = abs( $i - $self->exonCount + 1);
		}
		
		# identify UTRs
		# we will identify by comparing the cdsStart and cdsStop relative
		# to the exon coordinates
		# the primary tag is determined by the exon strand orientation
		my ($start, $stop, $tag);
		# in case we need to build two UTRs
		my ($start2, $stop2, $tag2);
		
		# Split 5'UTR, CDS, and 3'UTR all on the same exon
		if (
			$self->exonStarts->[$i] < $self->cdsStart
			and
			$self->exonEnds->[$i] > $self->cdsEnd
		) {
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
		elsif (
			$self->exonStarts->[$i] < $self->cdsStart
			and
			$self->exonEnds->[$i] < $self->cdsStart
		) {
			# the exon start/end is entirely before the cdsStart
			$start = $self->exonStarts->[$i];
			$stop  = $self->exonEnds->[$i];
			$tag   = $transcript->strand == 1 ? 'five_prime_UTR' : 'three_prime_UTR';
		}
		
		# Split 5'UTR & CDS on forward, 3'UTR & CDS
		elsif (
			$self->exonStarts->[$i] < $self->cdsStart
			and
			$self->exonEnds->[$i] >= $self->cdsStart
		) {
			# the start/stop codon is in this exon
			# we need to make the UTR out of a portion of this exon 
			$start = $self->exonStarts->[$i];
			$stop  = $self->cdsStart - 1;
			$tag   = $transcript->strand == 1 ? 'five_prime_UTR' : 'three_prime_UTR';
		}
		
		# CDS only
		elsif (
			$self->exonStarts->[$i] >= $self->cdsStart
			and
			$self->exonEnds->[$i] <= $self->cdsEnd
		) {
			# CDS only exon
			next UTR_LOOP;
		}
		
		# Split 3'UTR & CDS on forward, 5'UTR & CDS
		elsif (
			$self->exonStarts->[$i] <= $self->cdsEnd
			and
			$self->exonEnds->[$i] > $self->cdsEnd
		) {
			# the stop/start codon is in this exon
			# we need to make the UTR out of a portion of this exon 
			$start = $self->cdsEnd + 1;
			$stop  = $self->exonEnds->[$i];
			$tag   = $transcript->strand == 1 ? 'three_prime_UTR' : 'five_prime_UTR';
		}
	
		# 3'UTR forward, 5'UTR reverse
		elsif (
			$self->exonStarts->[$i] > $self->cdsEnd
			and
			$self->exonEnds->[$i] > $self->cdsEnd
		) {
			# the exon start/end is entirely after the cdsStop
			# we have a 3'UTR
			$start = $self->exonStarts->[$i];
			$stop  = $self->exonEnds->[$i];
			$tag   = $transcript->strand == 1 ? 'three_prime_UTR' : 'five_prime_UTR';
		}
		
		# Something else?
		else {
			my $warning = "Warning: A malformed UTR that doesn't match known criteria: ";
			$warning .= "cdsStart " . $self->cdsStart;
			$warning .= ", cdsEnd " . $self->cdsEnd;
			$warning .= ", exonStart " . $self->exonStarts->[$i];
			$warning .= ", exonEnd " . $self->exonEnds->[$i];
			warn $warning;
			next UTR_LOOP;
		}
		
		## Generate the UTR objects
		my $utr;
			
		# look for existing utr
		if ($ucsc->share and $gene) {
			$utr = $self->find_existing_subfeature($gene, $tag, $start, $stop); 
		}
			
		# otherwise build the UTR object
		unless ($utr) {
			$utr = $SFCLASS->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-start         => $start,
				-end           => $stop,
				-strand        => $transcript->strand,
				-phase         => '.',
				-primary_tag   => $tag,
				-primary_id    => $transcript->primary_id . ".utr$number",
			);
			$utr->display_name( $transcript->display_name . ".utr$number" ) if 
				$ucsc->do_name;
		}
		
		# store this utr seqfeature in a temporary array
		push @utrs, $utr;
		
		# build a second UTR object as necessary
		if ($start2) {
			my $utr2;
			
			# look for existing utr
			if ($ucsc->share) {
				$utr2 = $self->find_existing_subfeature($gene, $tag2, $start2, $stop2); 
			}
			
			# otherwise build the utr
			unless ($utr2) {
				$utr2 = $SFCLASS->new(
					-seq_id        => $transcript->seq_id,
					-source        => $transcript->source,
					-start         => $start2,
					-end           => $stop2,
					-strand        => $transcript->strand,
					-phase         => '.',
					-primary_tag   => $tag2,
					-primary_id    => $transcript->primary_id . ".utr$number" . "a",
				);
				$utr2->display_name( $transcript->display_name . ".utr$number" . "a" ) 
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
	my ($self, $transcript) = @_;
	
	# we will NOT collapse CDS features since we cannot guarantee that a shared 
	# CDS will have the same phase, since phase is dependent on the translation 
	# start 
	
	# we will scan each exon and look for a potential CDS and build it
	my @cdss;
	my $phase = 0; # initialize CDS phase and keep track as we process CDSs 
	CDS_LOOP:
	for (my $i = 0; $i < $self->exonCount; $i++) {
		
		# transform index for reverse strands
		my $j;
		if ($transcript->strand == 1) {
			# forward strand
			$j = $i;
		}
		else {
			# reverse strand
			# flip the index for exon starts and stops so that we 
			# always progress 5' -> 3' 
			# this ensures the phase is accurate from the start codon
			$j = abs( $i - $self->exonCount + 1);
		}
		
		# identify CDSs
		# we will identify by comparing the cdsStart and cdsStop relative
		# to the exon coordinates
		my ($start, $stop);
		
		# Split 5'UTR, CDS, and 3'UTR all on the same exon
		if (
			$self->exonStarts->[$j] < $self->cdsStart
			and
			$self->exonEnds->[$j] > $self->cdsEnd
		) {
			# exon contains the entire CDS
			$start = $self->cdsStart;
			$stop  = $self->cdsEnd;
		}
		
		# 5'UTR forward, 3'UTR reverse
		elsif (
			$self->exonStarts->[$j] < $self->cdsStart
			and
			$self->exonEnds->[$j] < $self->cdsStart
		) {
			# no CDS in this exon
			next CDS_LOOP;
		}
		
		# Split 5'UTR & CDS on forward, 3'UTR & CDS
		elsif (
			$self->exonStarts->[$j] < $self->cdsStart
			and
			$self->exonEnds->[$j] >= $self->cdsStart
		) {
			# the start/stop codon is in this exon
			# we need to make the CDS out of a portion of this exon 
			$start = $self->cdsStart;
			$stop  = $self->exonEnds->[$j];
		}
		
		# CDS only
		elsif (
			$self->exonStarts->[$j] >= $self->cdsStart
			and
			$self->exonEnds->[$j] <= $self->cdsEnd
		) {
			# entire exon is CDS
			$start = $self->exonStarts->[$j];
			$stop  = $self->exonEnds->[$j];
		}
	
		# Split 3'UTR & CDS on forward, 5'UTR & CDS
		elsif (
			$self->exonStarts->[$j] <= $self->cdsEnd
			and
			$self->exonEnds->[$j] > $self->cdsEnd
		) {
			# the stop/start codon is in this exon
			# we need to make the CDS out of a portion of this exon 
			$start = $self->exonStarts->[$j];
			$stop  = $self->cdsEnd;
		}
	
		# 3'UTR forward, 5'UTR reverse
		elsif (
			$self->exonStarts->[$j] > $self->cdsEnd
			and
			$self->exonEnds->[$j] > $self->cdsEnd
		) {
			# the exon start/end is entirely after the cdsStop
			# we have entirely 5' or 3'UTR, no CDS
			next CDS_LOOP;
		}
		
		# Something else?
		else {
			my $warning = "Warning: A malformed CDS that doesn't match known criteria: ";
			$warning .= "cdsStart " . $self->cdsStart;
			$warning .= ", cdsEnd " . $self->cdsEnd;
			$warning .= ", exonStart " . $self->exonStarts->[$j];
			$warning .= ", exonEnd " . $self->exonEnds->[$j];
			warn $warning;
			next CDS_LOOP;
		}
			
		# build the CDS object
		my $cds = $SFCLASS->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source,
			-start         => $start,
			-end           => $stop,
			-strand        => $transcript->strand,
			-phase         => $phase,
			-primary_tag   => 'CDS',
			-primary_id    => $transcript->primary_id . ".cds$i", 
			-display_name  => $transcript->display_name . ".cds$i",
		);
		# the id and name still use $i for labeling to ensure numbering from 0
		
		# store this utr seqfeature in a temporary array
		push @cdss, $cds;
		
		# reset the phase for the next CDS
			# phase + (3 - (length % 3)), readjust to 0..2 if necessary
			# adapted from Barry Moore's gtf2gff3.pl script
		$phase = $phase + (3 - ( $cds->length % 3) );
		$phase -=3 if $phase > 2;
	}
	
	# associate found UTRs with the transcript
	foreach my $cds (@cdss) {
		$transcript->add_SeqFeature($cds);
	}
}

sub add_codons {
	my ($self, $transcript, $gene) = @_;
	my $ucsc = $self->ucsc;
	
	# generate the start and stop codons
	my ($start_codon, $stop_codon);
	if ($transcript->strand == 1) {
		# forward strand
		
		# share codons if possible
		if ($ucsc->share and $gene) {
			$start_codon = $self->find_existing_subfeature($gene, 'start_codon', 
				$self->cdsStart, $self->cdsStart + 2);
			$stop_codon = $self->find_existing_subfeature($gene, 'stop_codon', 
				$self->cdsEnd - 2, $self->cdsEnd);
		}
		
		# start codon
		unless ($start_codon) {
			$start_codon = $SFCLASS->new(
					-seq_id        => $transcript->seq_id,
					-source        => $transcript->source,
					-primary_tag   => 'start_codon',
					-start         => $self->cdsStart,
					-end           => $self->cdsStart + 2,
					-strand        => 1,
					-phase         => 0,
					-primary_id    => $transcript->primary_id . '.start_codon',
			);
			$start_codon->display_name( $transcript->display_name . '.start_codon' ) if 
				$ucsc->do_name;
		}
		
		# stop codon
		unless ($stop_codon) {
			$stop_codon = $SFCLASS->new(
					-seq_id        => $transcript->seq_id,
					-source        => $transcript->source,
					-primary_tag   => 'stop_codon',
					-start         => $self->cdsEnd - 2,
					-end           => $self->cdsEnd,
					-strand        => 1,
					-phase         => 0,
					-primary_id    => $transcript->primary_id . '.stop_codon',
			);
			$stop_codon->display_name( $transcript->display_name . '.stop_codon' ) if 
				$ucsc->do_name;
		}
	}
	
	else {
		# reverse strand
		
		# share codons if possible
		if ($ucsc->share and $gene) {
			$stop_codon = $self->find_existing_subfeature($gene, 'stop_codon', 
				$self->cdsStart, $self->cdsStart + 2);
			$start_codon = $self->find_existing_subfeature($gene, 'start_codon', 
				$self->cdsEnd - 2, $self->cdsEnd);
		}
		
		# stop codon
		unless ($stop_codon) {
			$stop_codon = $SFCLASS->new(
					-seq_id        => $transcript->seq_id,
					-source        => $transcript->source,
					-primary_tag   => 'stop_codon',
					-start         => $self->cdsStart,
					-end           => $self->cdsStart + 2,
					-strand        => -1,
					-phase         => 0,
					-primary_id    => $transcript->primary_id . '.stop_codon',
			);
			$stop_codon->display_name( $transcript->display_name . '.stop_codon' ) if 
				$ucsc->do_name;
		}
		
		# start codon
		unless ($start_codon) {
			$start_codon = $SFCLASS->new(
					-seq_id        => $transcript->seq_id,
					-source        => $transcript->source,
					-primary_tag   => 'start_codon',
					-start         => $self->cdsEnd - 2,
					-end           => $self->cdsEnd,
					-strand        => -1,
					-phase         => 0,
					-primary_id    => $transcript->primary_id . '.start_codon',
					-display_name  => $transcript->primary_id . '.start_codon',
			);
			$start_codon->display_name( $transcript->display_name . '.start_codon' ) if 
				$ucsc->do_name;
		}
	}
	
	# associate with transcript
	$transcript->add_SeqFeature($start_codon);
	$transcript->add_SeqFeature($stop_codon);
}

sub find_gene {
	my $self = shift;
	
	# check if a gene with this name exists
	if (exists $self->ucsc->{gene2seqf}->{lc $self->gene_name} ) {
		# we found a gene with the same name
		# pull out the gene seqfeature(s) array reference
		# there may be more than one gene
		my $genes = $self->ucsc->{gene2seqf}->{ lc $self->gene_name };
		
		# check that the current transcript intersects with the gene
		# sometimes we can have two separate transcripts with the 
		# same gene name, but located on opposite ends of the chromosome
		# part of a gene family, but unlikely the same gene 200 Mb in 
		# length
		foreach my $g (@$genes) {
			if ( 
				# overlap method borrowed from Bio::RangeI
				($g->strand == $self->strand) and not (
					$g->start > $self->txEnd or 
					$g->end < $self->txStart
				)
			) {
				# gene and transcript overlap on the same strand
				# we found the intersecting gene
				return $g;
			}
		}
	}
	return;
}

sub find_existing_subfeature {
	my ($self, $gene, $type, $start, $stop) = @_;
	
	# we will try to find a pre-existing subfeature at identical coordinates
	foreach my $transcript ($gene->get_SeqFeatures()) {
		# walk through transcripts
		foreach my $subfeature ($transcript->get_SeqFeatures()) {
			# walk through subfeatures of transcripts
			if (
				$subfeature->primary_tag eq $type and
				$subfeature->start == $start and 
				$subfeature->end   == $stop
			) {
				# we found a match
				return $subfeature;
			}
		}
	}
	return;
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

