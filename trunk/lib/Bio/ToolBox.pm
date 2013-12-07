package Bio::ToolBox;

our $VERSION = 1.14;

1;

=head1 NAME

Bio::ToolBox - Tools for querying and analysis of genomic data

=head1 DESCRIPTION

These are libraries for the Bio::ToolBox package of bioinformatics 
analysis tools. The original project was focused primarily on the 
generation of functional end-user scripts for bioinformatic data 
conversion, processing, collection, and analysis. The libraries 
were initially simply shared code for these scripts. However, they 
are fully documented, and could be used as the basis for writing 
additional analysis scripts.

In many cases, the Bio::ToolBox scripts provide an abstraction 
layer over a variety of different specialized BioPerl-style 
modules. For example, there is a special emphasis on the collection 
data values for defined genomic coordinate regions, regardless of 
whether the values come from a GFF database, Bam file, BigWig file, 
etc. 

There is also an emphasis on using industry-standard file formats, 
such as BED and GFF3, as well as simple tab-delimited text files 
that may be imported or exported from a variety of other applications. 
Retention of metadata is accomplished by simply prefixing comment 
lines at the beginning of the file.

Currently the libraries are implemented as simple subroutines that 
must be imported into a script's namespace. They are not implemented 
as object-oriented modules, although that may be partially rectified 
in the future.

The libraries are broken into four primary categories.

=over 4

=item Bio::ToolBox::data_helper

This library prepares a Bio::ToolBox data structure (previously referred 
colloquially as a tim data structure) that most of the scripts employ. 
This is essentially a complex hash of metadata and an array of arrays for 
the data table.

This is the library most likely to be turned into friendly object-oriented 
code, as the data structure just needs to be blessed and simple methods 
written.

=item Bio::ToolBox::db_helper

This library works with data libraries, including opening, querying, and 
collecting from a variety of BioPerl-style Bio::DB::* databases, including 
SeqFeature::Store, Bam, BigWig, BigBed, BigWigSet, and USeq. The functions 
are fairly well abstracted, and the library will take care of handling 
the database specifics appropriately.

=item Bio::ToolBox::file_helper

This library takes care of file input and output, especially for the 
Bio::ToolBox data format (previously referred colloquially as the tim data 
format), which is just a tab-delimited text table with commented 
metadata lines at the beginning. It transparently handles common standard 
bioinformatic file formats, including BED and GFF, and gzip compression.

=item Bio::ToolBox::big_helper

This library takes care of converting text versions of wiggle, bedGraph, and 
bed file formats into UCSC BigWig and BigBed formats.

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0. 
