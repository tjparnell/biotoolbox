# NAME

Bio::ToolBox - Tools for querying and analysis of genomic data

# DESCRIPTION

These libraries provide a useful interface for working with 
bioinformatic data. Many bioinformatic data analysis revolves 
around working with tables of information, including lists of 
genomic annotation (genes, promoters, etc.) or defined regions 
of interest (epigenetic enrichment, transcription factor binding 
sites, etc.). This library works with these tables and provides 
a set of common tools for working with them.

- Opening and saving common tab-delimited text formats
- Support for BED, GFF, VCF, narrowPeak files
- Scoring intervals and annotation with datasets from microarray or sequencing
experiments, including ChIPSeq, RNASeq, and more
- Support for Bam, BigWig, BigBed, wig, and USeq data formats
- Works with any genomic annotation in GTF, GFF3, and various UCSC formats

The libraries provide a unified and integrated approach to analyses. 
In many cases, they provide an abstraction layer over a variety of 
different specialized BioPerl and related modules. Instead of 
writing numerous scripts specialized for each data format (wig, 
bigWig, Bam), one script can now work with any data format. 

# INSTALLATION

Installation is simple with the standard Perl incantation.

    perl ./Build.PL
    ./Build
    ./Build test
    ./Build install

Released versions may be obtained though the CPAN repository using 
your favorite package manager. 


## Additional external modules

To make the installation as lean and simple as possible, only the minimal 
additional Perl modules are required, while the remainder are only 
recommended. These can be installed subsequently as necessary as the need 
arises. Most of the database adapters, including those for Bam, BigWig, 
and BigBed, require external C library dependencies that must be compiled 
separately. See the respective modules for installation instructions.

Most scripts should fail gently with warnings about missing modules.

Suggested modules and adapters include the following:

- [Bio::Perl](https://metacpan.org/pod/Bio::Perl)
- [Bio::DB::BigWig](https://metacpan.org/pod/Bio::DB::BigWig)
- [Bio::DB::BigBed](https://metacpan.org/pod/Bio::DB::BigBed)
- [Bio::DB::Sam](https://metacpan.org/pod/Bio::DB::Sam)
- [Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS)
- [Bio::DB::USeq](https://metacpan.org/pod/Bio::DB::USeq)
- [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager)

# LIBRARIES

Several library modules are included in this distribution. These are a highlight 
of the primary user-oriented libraries that are available.

- [Bio::ToolBox::Data](https://metacpan.org/pod/Bio::ToolBox::Data)

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
    and API will be familiar to users of [Bio::Perl](https://metacpan.org/pod/Bio::Perl).

- [Bio::ToolBox::Data::Feature](https://metacpan.org/pod/Bio::ToolBox::Data::Feature)

    This is the object class for working with individual rows in a table 
    of data. It provides a number of conventions for working with the rows 
    in a standard fashion, for example returning the start column value  
    regardless of which column it is or whether the table is bed or gff or 
    an arbitrary text file. A number of convenience methods are present for 
    collecting data from data files. This module is not used directly by the 
    user, but its objects are returned when using [Bio::ToolBox::Data](https://metacpan.org/pod/Bio::ToolBox::Data) iterators.

- Annotation parsers

    Included are two generic parsers for loading an entire genome-worth of 
    annotation into memory within a reasonably short amount of time. 

    - [Bio::ToolBox::parser::gff](https://metacpan.org/pod/Bio::ToolBox::parser::gff)

        This parses both GTF and GFF3 file formats. Unlike many other GFF parsers 
        that work line-by-line only, this maintains parent and child hierarchical 
        relationships as parent feature and child subfeatures. To further maintain 
        control and reduce unnecessary parsing, unwanted feature types can be 
        selectively skipped.

    - [Bio::ToolBox::parser::ucsc](https://metacpan.org/pod/Bio::ToolBox::parser::ucsc)

        This parses various UCSC file formats, including different refFlat, GenePred, 
        and knownGene flavors. Genes, transcripts, and exons are assembled into 
        hierarchical child-parent relationships as desired.

- [Bio::ToolBox::SeqFeature](https://metacpan.org/pod/Bio::ToolBox::SeqFeature)

    This is a fast, lean, simple object class for representing genomic features. 
    It supports, for the most part, the [Bio::SeqFreatureI](https://metacpan.org/pod/Bio::SeqFreatureI) and [Bio::RangeI](https://metacpan.org/pod/Bio::RangeI) API 
    interface without the dependencies. It uses an unorthodox blessed-array object 
    structure, which provides measurable improvements in memory consumption and 
    speed when loading thousands of annotated SeqFeature objects (think hg19 or hg38 
    annotation). 

- [Bio::ToolBox::GeneTools](https://metacpan.org/pod/Bio::ToolBox::GeneTools)

    This is a collection of exportable functions for working with [Bio::SeqFeatureI](https://metacpan.org/pod/Bio::SeqFeatureI) 
    compliant objects representing genes and transcripts. It works with objects derived 
    from one of the ["Annotation parsers"](#annotation-parsers) or a [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store) database. 
    The functions make hard things easy, such as identifying whether a transcript is 
    coding or not (is it encoded in the `primary_tag` or `source_tag` or GFF 
    attribute or does it have `CDS` subfeatures?), or identify the alternative exons 
    or introns of a multi-transcript gene, or pull out the `5'` UTR (which is likely 
    not explicitly defined in the table).

# SCRIPTS

The BioToolBox package comes complete with a suite of high-quality production-ready 
scripts ready for a variety of analyses. Look in the scripts folder for details. 
A sampling of what can be done include the following:

- Annotated feature collection and selection
- Data collection and scoring for features
- Data file format manipulation and conversion
- Low-level processing of sequencing data into customizable wig representation

Scripts have built-in documentation. Execute the script without any options to print 
a synopsis of available options, or add `--help` to print the full documentation.

# CONFIGURATION

There is a small INI-style configuration file, `.biotoolbox.cfg`, written in your 
home directory, which can include paths to helper files and database configurations.

# AUTHOR

Timothy J. Parnell, PhD
Huntsman Cancer Institute
University of Utah
Salt Lake City, UT, 84112

# LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. 
