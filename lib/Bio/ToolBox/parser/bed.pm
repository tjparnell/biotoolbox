package Bio::ToolBox::parser::bed;

use warnings;
use strict;
use Carp qw(carp cluck croak confess);
use base 'Bio::ToolBox::Parser';
use Bio::ToolBox::Data;
use Module::Load;

our $VERSION = '1.70';

sub new {
	my $class = shift;
	return $class->SUPER::new(@_);
}

sub open_file {
	my $self     = shift;
	my $filename = shift || undef;

	# check file
	if ( $filename and $self->file and $filename ne $self->file ) {
		confess 'Must open new files with new Parser object!';
	}
	$filename ||= $self->file;
	unless ($filename) {
		cluck "No file name passed!\n";
		return;
	}
	if ( defined $self->{fh} ) {
		return 1;
	}

	# check file format type
	my $filetype = $self->filetype || undef;
	unless ($filetype) {
		( my $flavor, $filetype ) = Bio::ToolBox::Data->taste_file($filename);
		unless ( $flavor eq 'bed' ) {
			confess 'File is not a BED file!!! How did we get here?';
		}
		$self->{filetype} = $filetype;
	}

	# determine converter subroutine and set column expectations
	if ( $filetype eq 'gappedPeak' ) {
		$self->{convertor_sub} = \&_parse_gappedPeak;
		$self->do_exon(1);    # always parse sub peaks as exons
							  # gappedPeak is essentially bed12 with extra columns
		  # we will use existing code from the ucsc parser to convert bed12 to seqfeatures
		  # we need more object stuff that the ucsc parser expects
		load 'Bio::ToolBox::parser::ucsc::builder';
		$self->{bed}         = 15;
		$self->{id2count}    = {};
		$self->{refseqsum}   = {};
		$self->{refseqstat}  = {};
		$self->{kgxref}      = {};
		$self->{ensembldata} = {};
	}
	elsif ( $filetype eq 'narrowPeak' ) {
		$self->{bed}           = 10;
		$self->{convertor_sub} = \&_parse_narrowPeak;
	}
	elsif ( $filetype eq 'broadPeak' ) {
		$self->{bed}           = 9;
		$self->{convertor_sub} = \&_parse_broadPeak;
	}
	elsif ( $filetype eq 'bedGraph' ) {
		$self->{bed}           = 4;
		$self->{convertor_sub} = \&_parse_bedGraph;
	}
	elsif ( $filetype eq 'bed12' ) {
		$self->{bed}           = 12;
		$self->{convertor_sub} = \&_parse_bed12;

		# we will use existing code from the ucsc parser to convert bed12 to seqfeatures
		# we need more object stuff that the ucsc parser expects
		load 'Bio::ToolBox::parser::ucsc::builder';
		$self->{id2count}    = {};
		$self->{refseqsum}   = {};
		$self->{refseqstat}  = {};
		$self->{kgxref}      = {};
		$self->{ensembldata} = {};
	}
	else {
		# an ordinary bed file
		$self->{bed}           = ( $filetype =~ /bed(\d+)/ )[0];
		$self->{convertor_sub} = \&_parse_bed;
	}

	# Open file
	my $fh = Bio::ToolBox::Data->open_to_read_fh($filename)
		or croak(" Unable to open file '$filename'!");

	# finish
	$self->{fh} = $fh;
	return 1;
}

sub typelist {
	my $self = shift;
	my $ft   = $self->filetype;

	# return generic lists based on what is expected
	if ( $ft eq 'bed12' ) {
		return 'mRNA,ncRNA,exon,CDS';
	}
	elsif ( $ft eq 'narrowPeak' or $ft eq 'broadPeak' ) {
		return 'peak';
	}
	elsif ( $ft eq 'gappedPeak' ) {
		return 'gappedPeak,peak';
	}
	else {
		# return generic
		return 'feature';
	}
}

sub next_feature {
	my $self = shift;

	# check that we have an open filehandle and not finished
	unless ( $self->fh ) {
		croak 'no Bed file loaded to parse!';
	}
	return if $self->{'eof'};

	# loop through the file
	while ( my $line = $self->fh->getline ) {
		if ( $line =~ /^#/ or $line =~ /^(?:track | browser)/x or $line !~ /\w+/ ) {
			push @{ $self->{comments} }, $line;
			$self->{line_count}++;
			next;
		}
		chomp $line;
		my $feature = &{ $self->{convertor_sub} }( $self, $line );
		$self->{line_count}++;
		unless ($feature) {
			printf STDERR "unable to make feature for line %d!\n", $self->{line_count};
			next;
		}

		# return the object, we do this while loop once per valid line
		return $feature;
	}

	# presumed end of file
	$self->{'eof'} = 1;
	$self->fh->close;
	return;
}

*parse_table = \&parse_file;
*parse_table if 0;  # avoid once warning

sub parse_file {
	my $self = shift;

	# check that we have an open filehandle
	unless ( $self->fh ) {
		confess 'no file loaded to parse!';
	}
	return 1 if $self->{'eof'};

	printf "  Parsing %s format file....\n", $self->filetype;

	while ( my $feature = $self->next_feature ) {

		# there are possibly lots and lots of features here
		# we shouldn't have to check IDs for uniqueness, but to maintain compatibility
		# with other parsers and for consistency, we do
		my $id = $feature->primary_id;
		if ( exists $self->{loaded}{$id} ) {
			my $i = 1;
			$id = $feature->primary_id . ".$i";
			while ( exists $self->{loaded}{$id} ) {
				$i++;
				$id = $feature->primary_id . ".$i";
			}

			# we have a unique id, now change it
			$feature->primary_id($id);
		}

		# record the feature
		$self->{loaded}{$id} = $feature;
		push @{ $self->{top_features} }, $feature;

		# check chromosome
		my $s = $feature->seq_id;
		unless ( exists $self->{seq_ids}{$s} ) {
			$self->{seq_ids}{$s} = $feature->end;
		}
		$self->{seq_ids}{$s} = $feature->end if $feature->end > $self->{seq_ids}{$s};
	}
	return 1;
}

sub _parse_narrowPeak {
	my ( $self, $line ) = @_;
	my @data = split /\t/, $line;
	unless ( scalar(@data) == 10 ) {
		croak sprintf( "narrowPeak line %d '%s' doesn't have 10 elements!",
			$self->{line_count}, $line );
	}

	# generate the basic SeqFeature
	my $feature = $self->{sfclass}->new(
		-seq_id      => $data[0],
		-start       => $data[1] + 1,
		-end         => $data[2],
		-name        => $data[3],
		-score       => $data[4],
		-strand      => $data[5],
		-primary_tag => 'peak',
		-primary_id  => sprintf( "%s:%d-%d", $data[0], $data[1], $data[2] ),
	);

	# add extra columns
	$feature->add_tag_value( 'signalValue', $data[6] );
	$feature->add_tag_value( 'pValue',      $data[7] );
	$feature->add_tag_value( 'qValue',      $data[8] );
	$feature->add_tag_value( 'peak',        $data[9] );

	return $feature;
}

sub _parse_broadPeak {
	my ( $self, $line ) = @_;
	my @data = split /\t/, $line;
	unless ( scalar(@data) == 9 ) {
		croak sprintf( "broadPeak line %d '%s' doesn't have 9 elements!",
			$self->{line_count}, $line );
	}

	# generate the basic SeqFeature
	my $feature = $self->{sfclass}->new(
		-seq_id      => $data[0],
		-start       => $data[1] + 1,
		-end         => $data[2],
		-name        => $data[3],
		-score       => $data[4],
		-strand      => $data[5],
		-primary_tag => 'peak',
		-primary_id  => sprintf( "%s:%d-%d", $data[0], $data[1], $data[2] ),
	);

	# add extra columns
	$feature->add_tag_value( 'signalValue', $data[6] );
	$feature->add_tag_value( 'pValue',      $data[7] );
	$feature->add_tag_value( 'qValue',      $data[8] );

	return $feature;
}

sub _parse_bedGraph {
	my ( $self, $line ) = @_;
	my @data = split /\t/, $line;
	unless ( scalar(@data) == 4 ) {
		croak sprintf( "bedGraph line %d '%s' doesn't have 4 elements!",
			$self->{line_count}, $line );
	}

	# generate the basic SeqFeature
	return $self->{sfclass}->new(
		-seq_id     => $data[0],
		-start      => $data[1] + 1,
		-end        => $data[2],
		-score      => $data[3],
		-primary_id => sprintf( "%s:%d-%d", $data[0], $data[1], $data[2] ),
	);
}

sub _parse_bed {
	my ( $self, $line ) = @_;
	my @data = split /\t/, $line;
	unless ( scalar(@data) == $self->{bed} ) {
		croak sprintf( "Bed line %d '%s' doesn't have %d elements!",
			$self->{line_count}, $line, $self->{bed} );
	}

	# generate the basic SeqFeature
	return $self->{sfclass}->new(
		-seq_id     => $data[0],
		-start      => $data[1] + 1,
		-end        => $data[2],
		-name       => $data[3] || undef,
		-score      => $data[4] || undef,
		-strand     => $data[5] || undef,
		-primary_id => sprintf( "%s:%d-%d", $data[0], $data[1], $data[2] ),
	);
}

sub _parse_bed12 {
	my ( $self, $line ) = @_;
	my @data = split /\t/, $line;
	unless ( scalar(@data) == $self->{bed} ) {
		croak sprintf( "Bed line %d '%s' doesn't have %d elements!",
			$self->{line_count}, $line, $self->{bed} );
	}

	# we will take advantage of pre-existing code in the UCSC parser to convert
	# a bed12 line to a fully-fledged processed transcript
	# we just have to go through a genePred format first
	# fortunately, the two are pretty similar in structure

# now convert from bed12
# Chromosome Start End Name Score Strand thickStart thickEnd itemRGB blockCount blockSizes blockStarts
# to genePred
# name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds

	# add missing data
	$data[6]  ||= $data[2];               # cdsStart
	$data[7]  ||= $data[2];               # cdsEnd
	$data[8]  ||= 0;                      # rgb values
	$data[9]  ||= 1;                      # exonCount
	$data[10] ||= $data[7] - $data[6];    # block size
	$data[11] ||= 0;                      # block starts

	# calculate exons
	my @exonSizes  = split /,/, $data[10];
	my @exonStarts = map { $data[1] + $_ } split /,/, $data[11];
	my @exonEnds;
	for ( my $i = 0; $i < $data[9]; $i++ ) {
		push @exonEnds, $exonStarts[$i] + $exonSizes[$i];
	}

	# calculate new genePred elements
	my @new = (
		$data[3],                    # name
		$data[0],                    # chrom
		$data[5],                    # strand
		$data[1],                    # txStart
		$data[2],                    # txStop
		$data[6],                    # cdsStart
		$data[7],                    # cdsEnd
		$data[9],                    # exonCount
		join( ',', @exonStarts ),    # exonStarts
		join( ',', @exonEnds ),      # exonEnds
	);

	# create builder
	my $builder = Bio::ToolBox::parser::ucsc::builder->new( \@new, $self );
	my $feature = $builder->build_transcript;
	$feature->add_tag_value( 'itemRGB', $data[8] );
	$feature->score( $data[4] );
	$feature->primary_id( sprintf( "%s:%d-%d", $data[0], $data[1], $data[2] ) );

	# change the primary ID to match other bed file behavior, not UCSC files'
	return $feature;
}

sub _parse_gappedPeak {
	my ( $self, $line ) = @_;
	my @data = split /\t/, $line;
	unless ( scalar(@data) == 15 ) {
		croak sprintf( "GappedPeak line %d '%s' doesn't have 15 elements!",
			$self->{line_count}, $line );
	}

	# we will take advantage of pre-existing code in the UCSC parser to convert
	# a gappedPeak line into main peak with subpeaks.
	# we just have to go through a genePred format first
	# fortunately, the two are pretty similar in structure

	# calculate exons, er, sub peaks
	my @exonSizes  = split /,/, $data[10];
	my @exonStarts = map { $data[1] + $_ } split /,/, $data[11];
	my @exonEnds;
	for ( my $i = 0; $i < $data[9]; $i++ ) {
		push @exonEnds, $exonStarts[$i] + $exonSizes[$i];
	}

	# calculate new genePred elements
	my @new = (
		$data[3],                    # name
		$data[0],                    # chrom
		$data[5],                    # strand
		$data[1],                    # txStart
		$data[2],                    # txStop
		$data[6],                    # cdsStart
		$data[7],                    # cdsEnd
		$data[9],                    # exonCount
		join( ',', @exonStarts ),    # exonStarts
		join( ',', @exonEnds ),      # exonEnds
	);

	# create builder and process
	my $builder = Bio::ToolBox::parser::ucsc::builder->new( \@new, $self );
	my $feature = $builder->build_transcript;

	# clean up feature and add extra values
	$feature->add_tag_value( 'itemRGB', $data[8] );
	$feature->score( $data[4] );
	$feature->primary_tag('gappedPeak');    # it is not a RNA
	$feature->primary_id( sprintf( "%s:%d-%d", $data[0], $data[1], $data[2] ) );

	# change the primary ID to match other bed file behavior, not UCSC files'
	$feature->add_tag_value( 'signalValue', $data[12] );
	$feature->add_tag_value( 'pValue',      $data[13] );
	$feature->add_tag_value( 'qValue',      $data[14] );
	foreach my $f ( $feature->get_SeqFeatures ) {
		$f->primary_tag('peak');
	}
	return $feature;
}

1;

__END__

=head1 NAME

Bio::ToolBox::parser::bed - Parser for BED-style formats

=head1 SYNOPSIS

  use Bio::ToolBox::Parser;
  my $filename = 'file.bed';
  
  my $Parser = Bio::ToolBox::Parser->new(
  	file    => $filename,
  ) or die "unable to open gff file!\n";
  # the Parser will taste the file and open the appropriate 
  # subclass parser, bed in this case
  
  while (my $feature = $Parser->next_top_feature() ) {
	# each $feature is parent SeqFeature object
  	printf "%s:%d-%d\n", $f->seq_id, $f->start, $f->end;
  }

=head1 DESCRIPTION

This is the BED-style specific parser subclass to the L<Bio::ToolBox::Parser>
object, and as such inherits generic methods from the parent. File formats 
include the following. 

=over 4 

=item Bed

L<Bed|http://genome.ucsc.edu/FAQ/FAQformat.html#format1> files may have 3-12 columns, 
where the first 3-6 columns are basic information about the feature itself, and 
columns 7-12 are usually for defining subfeatures of a transcript model, including 
exons, UTRs (thin portions), and CDS (thick portions) subfeatures. This parser will 
parse these extra fields as appropriate into subfeature SeqFeature objects. Bed files 
are recognized with the file extension F<.bed>.

=item Bedgraph

L<BedGraph|http://genome.ucsc.edu/FAQ/FAQformat.html#format1.8> files are a type of 
wiggle format in Bed format, where the 4th column is a score instead of a name. BedGraph 
files are recognized by the file extension F<.bedgraph> or F<.bdg>.

=item narrowPeak

L<narrowPeak|http://genome.ucsc.edu/FAQ/FAQformat.html#format12> files are a specialized 
Encode variant of bed files with 10 columns (typically denoted as bed6+4), where the 
extra 4 fields represent score attributes to a narrow ChIPSeq peak. These files are 
parsed as a typical bed6 file, and the extra four fields are assigned to SeqFeature 
attribute tags C<signalValue>, C<pValue>, C<qValue>, and C<peak>, respectively. 
NarrowPeak files are recognized by the file extension F<.narrowPeak>. 

=item broadPeak

L<broadPeak|http://genome.ucsc.edu/FAQ/FAQformat.html#format13> files, like narrowPeak, 
are an Encode variant with 9 columns (bed6+3) representing a broad or extended interval 
of ChIP enrichment without a single "peak". The extra three fields are assigned to 
SeqFeature attribute tags C<signalValue>, C<pValue>, and C<qValue>, respectively.
BroadPeak files are recognized by the file extension F<.broadPeak>. 

=back

C<Track> and C<Browser> lines are generally ignored, although a C<track> definition 
line containing a C<type> key will be interpreted if it matches one of the above file 
types. 

=head2 SeqFeature default values

The SeqFeature objects built from the bed file intervals will have some inferred defaults. 

=over 4

=item Coordinate system

SeqFeature objects use the 1-based coordinate system, per the specification of 
L<Bio::SeqFeatureI>, so the 0-based start coordinates of bed files will always be 
parsed into 1-based coordinates.

=item C<display_name>

SeqFeature objects will use the name field (4th column in bed files), if present, as the 
C<display_name>. The SeqFeature object should default to the C<primary_id> if a name was 
not provided.

=item C<primary_id>

It will use a concatenation of the sequence ID, start (original 0-based), and 
stop coordinates as the C<primary_id>, for example 'chr1:0-100'. 

=item C<primary_tag>

Bed files don't have an inherent attribute of feature type (they are all the same 
type), so a default C<primary_tag> is assigned based on the file type. For peak 
files (F<narrowPeak> and F<broadPeak>) this is C<peak>, for F<gappedPeak> this is 
C<gappedPeak> and C<peak> (subfeatures), and for F<bed12> files with transcript models, 
the transcripts will be set to either C<mRNA> or C<ncRNA>, depending on the presence 
of interpreted CDS start and stop (thick coordinates).

=item C<source_tag>

Bed files don't have a concept of a source; default is "".

=item attribute tags

Extra columns in the F<narrowPeak> and F<broadPeak> formats are assigned to attribute tags 
as described above. The C<rgb> values set in bed12 files are also set to an attribute tag.

=back

=head1 METHODS

=head2 Initializing the parser object

In most cases, users should initialize an object using the generic 
L<Bio::ToolBox::Parser> object. 

These are class methods to initialize the parser with an annotation file 
and modify the parsing behavior. Most parameters can be set either upon 
initialization or as class methods on the object. Unpredictable behavior 
may occur if you implement these in the midst of parsing a file. 

Do not open subsequent files with the same object. Always create a new 
object to parse a new file.

=over 4

=item new

  my $parser = Bio::ToolBox::parser::bed->new($filename);
  my $parser = Bio::ToolBox::parser::bed->new(
      file    => 'file.bed',
      do_gene => 1,
      do_cds  => 1,
  );

Initiate a new Bed file parser object. Pass a single value (the bed file name) to 
open the file for parsing. Alternatively, pass an array of key 
value pairs to control how the table is parsed. These options are primarily for 
parsing bed12 files with subfeatures. Options include the following.

=over 4

=item file

Provide the path and file name for a Bed file. The file may be gzip compressed. 

=item source

Pass a string to be added as the source tag value of the SeqFeature objects. 

=item do_exon

=item do_cds

=item do_utr

=item do_codon

For Bed12 formats that represent transcripts, pass a boolean (1 or 0) value to
parse certain subfeatures, including C<exon>, C<CDS>, C<five_prime_UTR>, 
C<three_prime_UTR>, C<stop_codon>, and C<start_codon> features. Default is false.

=item class

Pass the name of a L<Bio::SeqFeatureI> compliant class that will be used to 
create the SeqFeature objects. The default is to use L<Bio::ToolBox::SeqFeature>, 
which is lighter-weight and consumes less memory. A suitable BioPerl alternative
is L<Bio::SeqFeature::Lite>.

=back

=back

=head2 Other methods

Additional methods for working with the parser object and the parsed 
SeqFeature objects.

=over 4

=item typelist

Returns a string representation of the type of SeqFeature types to be encountered in 
the file. Currently this returns generic strings, 'mRNA,ncRNA,exon,CDS' for bed12 
and 'feature' for everything else.

=back

=head1 SEE ALSO

L<Bio::ToolBox::Parser>, L<Bio::ToolBox::SeqFeature>, 
L<Bio::ToolBox::parser::ucsc>, L<Bio::ToolBox::parser::gff>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  



