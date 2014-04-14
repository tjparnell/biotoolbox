package Bio::ToolBox;

our $VERSION = 1.17;

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

=item * Support for BED and GFF files

=item * Scoring intervals with datasets from microarray and sequencing

=item * ChIPSeq, RNASeq, microarray expression, SNP detection

=item * Support for Bam, BigWig, BigBed, wig, and USeq data formats

=item * Intersection with other known annotation

=item * Works with any genomic annotation in GFF3 format

=back

The libraries provide a unified and integrated approach to analyses. 
In many cases, they provide an abstraction layer over a variety of 
different specialized BioPerl and related modules. Instead of 
writing numerous scripts specialized for each data format (wig, 
bigWig, Bam), one script can now work with any data format. 

In many cases, working with genomic annotation in databases assumes 
the use of Bio::DB::SeqFeature::Store formatted databases, available 
with a number of different backend support, including SQLite, MySQL, 
and others. 

=head1 SCRIPTS

The Bio::ToolBox comes complete with an entire suite of high 
quality scripts ready for a wide variety of analyses. 

=over 4

=item * Preparation of databases from public annotation sources

=item * Annotated feature collection and selection

=item * Data collection and scoring

=item * File format manipulation and conversion

=item * Some low-level processing of raw data

=item * Simple analysis and graphing of collected data

=back

=head1 LIBRARIES

The libraries and modules are available to extend existing 
scripts or to write your own. 

There is one primary module, which provides a convenient 
object-oriented interface for working with data files and 
collecting data. It is the one most users will want to work with.

The remaining modules are support and general modules. They 
are collections of exportable subroutines without a convenient 
object-oriented interface. 

=over 4

=item Bio::ToolBox::Data

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
will be familiar to users of BioPerl, from which this module draws 
heavily.

=item Bio::ToolBox::data_helper

This is a helper library for working with the Bio::ToolBox Data 
structure (previously referred colloquially as a tim data structure). 
This is essentially a complex hash of metadata and an array of arrays 
representing the data table.

=item Bio::ToolBox::db_helper

This helper library interacts with databases, including a variety of 
BioPerl-style Bio::DB::* databases such as SeqFeature::Store, Bam, 
BigWig, BigBed, BigWigSet, and USeq. In most cases, unless specifically 
stated otherwise, most database functions assume the use of 
Bio::DB::SeqFeature::Store databases, particularly with regards to 
genomic annotation. The functions are fairly well abstracted, and the 
library will take care of handling the database specifics appropriately.

=item Bio::ToolBox::file_helper

This helper library takes care of file input and output, especially for the 
native Bio::ToolBox Data format, which is just a tab-delimited text table 
with commented metadata lines at the beginning. It transparently handles 
common standard bioinformatic file formats, including BED and GFF, and 
gzip compression.

=item Bio::ToolBox::big_helper

This helper library takes care of converting text versions of wiggle, 
bedGraph, and Bed file formats into UCSC BigWig and BigBed formats.

=back

=head1 EXAMPLE

The following is a simplified example of using the library to create 
a data collection script using an input file of gene identifiers 
(from a database) or a BED file of coordinates. The dataset could be 
a Bam, BigBed, BigWig, or USeq file of genomic data.

  
  use Bio::ToolBox::Data;
  
  my $file = $ARGV[0];
  my $dataset = $ARGV[1];
  unless ($file and $dataset) {
      die "usage: $0 <input_file> <dataset_file>\n";
  }
  
  ### Open a pre-existing file
  my $Data = Bio::ToolBox::Data->new(
        file    => $file,
  );
  
  ### Add new dataset column and metadata
  my $index = $Data->add_column('Data');
  $Data->metadata($index, 'dataset', $dataset);
  $Data->metadata($index, 'method', 'mean');
  $Data->metadata($index, 'stranded', 'sense');
  
  ### Iterate through the Data structure one row at a time
  my $stream = $Data->row_stream;
  while (my $row = $stream->next_row) {
  	  my $value = $row->get_score(
  	      dataset   => $dataset,
  	      method    => 'mean',
  	      stranded  => 'sense',
  	  );
  	  $row->value($index, $value);
  }
  
  ### write the data to file
  my $success = $Data->write_file();
  print "wrote file $success\n"; # file extension will be automatically changed 
  

The example illustrates the simplicity of opening an input file, automatically 
identifying the file format and necessary columns for establishing genomic 
intervals, automatically handling the genomic dataset file from which to 
collect the data, iterating through the table, collecting the score for each 
genomic interval from the dataset, and writing the changed table to a new 
file. Since the input format (BED file) is no longer proper, it will 
automatically change the extension and write a text file complete with metadata.

The script could easily be modified to add or alter functionality. For 
example, use a different method of combining scores for each interval, or 
restricting scoring to a fraction of the defined interval. For a fully-featured 
version of this data collection script, see the BioToolBox script B<get_datasets.pl>.

=head1 REPOSITORY

Source code for the Bio::ToolBox package is maintained at 
L<http://code.google.com/p/biotoolbox/>. Extensive documentation, including How To 
documents, installation guides, and script documentation, can be found there as well.

Bugs and issues should be submitted at L<https://code.google.com/p/biotoolbox/issues/list>.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0. 
