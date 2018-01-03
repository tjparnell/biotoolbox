package Bio::ToolBox;

our $VERSION = '1.53';

1;

=head1 NAME

Bio::ToolBox - Tools for querying and analysis of genomic data

=head1 DESCRIPTION

These libraries provide a useful interface for working with 
bioinformatic data. Many bioinformatic data analysis revolves 
around working with tables of information, including lists of 
genomic annotation (genes, promoters, etc.) or defined regions 
of interest (epigenetic enrichment, transcription factor binding 
sites, etc.). This library works with these tables and provides 
a set of common tools for working with them.

=over 4

=item * Opening and saving common tab-delimited text formats

=item * Support for BED, GFF, VCF, narrowPeak files

=item * Scoring intervals with datasets from microarray and sequencing

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

=item Annotation parsers

Included are two generic parsers for loading an entire genome-worth of 
annotation into memory within a reasonably short amount of time. 

=over 4

=item L<Bio::ToolBox::parser::gff>

This parses both GTF and GFF3 file formats. Unlike many other GFF parsers 
that work line-by-line only, this maintains parent and child hierarchical 
relationships as parent feature and child subfeatures. To further maintain 
control and reduce unnecessary parsing, unwanted feature types can be 
selectively skipped.

=item L<Bio::ToolBox::parser::ucsc>

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
or introns of a multi-transcript gene, or pull out the C<5'> UTR (which is likely 
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
