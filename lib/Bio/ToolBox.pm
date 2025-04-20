package Bio::ToolBox;

use warnings;
use strict;
use Carp qw(cluck);
use Bio::ToolBox::Data;

our $VERSION = '2.02';

sub load_file {
	my $self = shift;
	if ( scalar(@_) == 1 ) {
		return Bio::ToolBox::Data->new( file => $_[0], );
	}
	else {
		return Bio::ToolBox::Data->new(@_);
	}
}

sub parse_file {
	my $self = shift;
	if ( scalar(@_) == 1 ) {
		return Bio::ToolBox::Data->new(
			file       => $_[0],
			parse      => 1,
			simplify   => 1,
			subfeature => 'exon,cds,utr',
		);
	}
	else {
		my %args = @_;
		$args{parse} ||= 1;    # make sure this is present
		return Bio::ToolBox::Data->new(%args);
	}
}

sub new_data {
	my $self = shift;
	if ( scalar(@_) and $_[0] =~ m/^(?: columns | datasets )$/x ) {

		# looks like a correctly formatted list
		return Bio::ToolBox::Data->new(@_);
	}
	else {
		# put provided list into an array
		return Bio::ToolBox::Data->new( columns => [@_] );
	}
}

sub open_file {
	my $self = shift;
	cluck('Bio::ToolBox->open_file has been replaced by read_file');
	return Bio::ToolBox::Data->open_to_read_fh(@_);
}

sub read_file {
	my $self = shift;
	return Bio::ToolBox::Data->open_to_read_fh(@_);
}

sub write_file {
	my $self = shift;
	return Bio::ToolBox::Data->open_to_write_fh(@_);
}

sub open_database {
	my $self = shift;
	return Bio::ToolBox::Data->open_new_database(@_);
}

sub bam_adapter {
	my $self = shift;
	return Bio::ToolBox::Data->bam_adapter(@_);
}

sub big_adapter {
	my $self = shift;
	return Bio::ToolBox::Data->big_adapter(@_);
}

1;

__END__


=head1 NAME

Bio::ToolBox - Tools for querying and analysis of genomic data

=head1 DESCRIPTION

The Bio::ToolBox libraries provide a useful interface for working 
with bioinformatic data. Many bioinformatic data analysis revolves 
around working with tables of information, including lists of 
genomic annotation (genes, promoters, etc.) or defined regions 
of interest (epigenetic enrichment, transcription factor binding 
sites, etc.). This library works with these tables and provides 
a set of common tools for working with them.

=over 4

=item * Opening and saving common tab-delimited text formats

=item * Support for BED, GFF, VCF, narrowPeak files

=item * Scoring intervals and annotation with datasets from microarray or sequencing experiments, including ChIPSeq, RNASeq, and more

=item * ChIPSeq, RNASeq, microarray expression

=item * Support for Bam, BigWig, BigBed, wig, and USeq data formats

=item * Works with any genomic annotation in GTF, GFF3, and UCSC formats

=back

The libraries provide a unified and integrated approach to analyses. 
In many cases, they provide an abstraction layer over a variety of 
different specialized BioPerl and related modules. Instead of 
writing numerous scripts specialized for each data format (wig, 
bigWig, Bam), one script can now work with any data format. 

=head1 LIBRARIES

The libraries and modules are available to extend existing 
scripts or to write your own. 

=over 4

=item L<Bio::ToolBox::Data>

This is the primary library module for working with a table of data, 
either generated as a new list from a database of annotation, or 
opened from a tab-delimited text file, for example a BED file of 
regions. Columns and rows of data may be added, deleted, or manipulated 
with ease. 

Additionally, genomic data may be collected from a wide variety of 
sources using the information in the data table. For example, scoring 
microarray or sequencing data for each interval listed in the data 
table.

This module uses an object-oriented interface. Many of the methods 
and API will be familiar to users of L<Bio::Perl>.

=item L<Bio::ToolBox::Data::Feature>

This is the object class for working with individual rows in a table 
of data. It provides a number of conventions for working with the rows 
in a standard fashion, for example returning the start column value  
regardless of which column it is or whether the table is bed or gff or 
an arbitrary text file. A number of convenience methods are present for 
collecting data from data files. This module is not used directly by the 
user, but its objects are returned when using L<Bio::ToolBox::Data> iterators.

=item L<Bio::ToolBox::Parser>

This is the working base class for parsing annotation files, including
BED and related formats, GFF, GTF, GFF3, and UCSC-derived refFlat, 
genePred, and genePredExt tables. This is designed to slurp an entire 
genome-worth of annotation into memory within a reasonably short amount
of time. Sub-classes include the following.

=over 4

=item L<Bio::ToolBox::Parser::bed>

This parses simple BED formats (3-6 columns), gene-based BED files (12 columns),
ENCODE-style peak formats (narrowPeak, broadPeak, and gappedPeak), and other 
BED-related derivatives. Gene-based BED12 files are parsed into hierarchical
parent and child subfeatures.

=item L<Bio::ToolBox::Parser::gff>

This parses both GTF and GFF3 file formats. Unlike many other GFF parsers 
that work line-by-line only, this maintains parent and child hierarchical 
relationships as parent feature and child subfeatures. To further maintain 
control and reduce unnecessary parsing, unwanted feature types can be 
selectively skipped.

=item L<Bio::ToolBox::Parser::ucsc>

This parses various UCSC file formats, including different refFlat, GenePred, 
and knownGene flavors. Genes, transcripts, and exons are assembled into 
hierarchical child-parent relationships as desired.

=back

=item L<Bio::ToolBox::SeqFeature>

This is a fast, lean, simple object class for representing genomic features. 
It supports, for the most part, the L<Bio::SeqFreatureI> and L<Bio::RangeI> API 
interface without the dependencies. It uses an unorthodox blessed-array object 
structure, which provides measurable improvements in memory consumption and 
speed when loading thousands of annotated SeqFeature objects (think hg19 or hg38 
annotation). 

=item L<Bio::ToolBox::GeneTools>

This is a collection of exportable functions for working with L<Bio::SeqFeatureI> 
compliant objects representing genes and transcripts. It works with objects derived 
from one of the L<"Annotation parsers"> or a L<Bio::DB::SeqFeature::Store> database. 
The functions make hard things easy, such as identifying whether a transcript is 
coding or not (is it encoded in the C<primary_tag> or C<source_tag> or GFF 
attribute or does it have C<CDS> subfeatures?), or identify the alternative exons 
or introns of a multi-transcript gene, or pull out the 5' UTR (which is likely 
not explicitly defined in the table).

=back

=head1 SCRIPTS

The BioToolBox package comes complete with a suite of high-quality production-ready 
scripts ready for a variety of analyses. Look in the scripts folder for details. 
A sampling of what can be done include the following:

=over 4

=item * Annotated feature collection and selection

=item * Data collection and scoring for features

=item * Data file format manipulation and conversion

=item * Low-level processing of sequencing data into customizable wig representation

Scripts have built-in documentation. Execute the script without any options to print 
a synopsis of available options, or add C<--help> to print the full documentation.

=back

=head1 USAGE

This module provides a handful of commonly used convenience methods 
as entry points to working with data files. Most of them use or 
return a L<Bio::ToolBox::Data> object.

=head2 Methods

=over 4

=item load_file

Open a tab-delimited text file as a L<Bio::ToolBox::Data> object. 
Simply pass the file path as a single argument. It assumes the first 
row is the column headers, and comment lines begin with C<#>. 
Compressed files are transparently handled. See the 
L<Bio::ToolBox::Data> C<new> method for more details or options.

  $Data = Bio::ToolBox->load_file('myfile.txt');

=item parse_file

Parse an annotation file, such as BED, GTF, GFF3, UCSC genePred or 
refFlat file, into a L<Bio::ToolBox::Data> table. Each row in the 
resulting table is linked to a parsed SeqFeature gene object. See 
the L<Bio::ToolBox::Data> C<new> method for more details or options.
Default options include parsing subfeatures (exon, cds, and utr) and 
simple GFF attributes.

  $Data = Bio::ToolBox->parse_file('genes.gtf.gz');
    
=item new_data

Generate a new, empty L<Bio::ToolBox::Data> table with the given column 
names. Pass the names of the columns in the new table.

  $Data = Bio::ToolBox->new_data( qw(Name ID Score) );
    
=item read_file

Open a generic file handle for reading. It transparently handles 
compression as necessary. Returns an L<IO::File> object. Pass the 
file path as an argument. 
    
  $fh = Bio::ToolBox->read_file('mydata.txt.gz');
    
=item write_file

Open a generic file handle for writing. It transparently handles 
compression as necessary based on filename extension or passed 
options. It will use the C<pigz> multi-threaded, external, compression
utility if available. See the C<open_to_write_fh> method in 
<Bio::ToolBox::Data::file> for more information.

  $fh = Bio::ToolBox->write_file('mynewdata.txt.gz');

=item open_database

Open a binary database file, including Bam, bigWig, bigBed, Fasta, 
L<Bio::DB::SeqFeature::Store> SQLite file or named MySQL connection, 
USeq file, or any other supported binary or indexed file formats. 
Database type is transparently and automatically checked by looking for 
common file extensions, if present. See the C<open_db_connection> in 
L<Bio::ToolBox::db_helper> for more information.

  $db = Bio::ToolBox->open_database($database);
    
=back

=head1 REPOSITORY

Source code for the Bio::ToolBox package is maintained at 
L<https://github.com/tjparnell/biotoolbox/>. 

Bugs and issues should be submitted at L<https://github.com/tjparnell/biotoolbox/issues>.

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::SeqFeature::Store>, L<Bio::SeqFeatureI>, L<Bio::DB::BigWig>, 
L<Bio::DB::BigBed>, L<Bio::DB::Sam>, L<Bio::DB::HTS>, L<Bio::DB::USeq>, 
L<Bio::ViennaNGS>


=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

=head1 LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. 
