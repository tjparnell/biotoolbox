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

Basic installation is simple with the standard [Module::Build](https://metacpan.org/pod/Module::Build) 
incantation. This will get you a minimal installation that will work with 
text files (BED, GFF, GTF, etc), but not binary files. 

    perl ./Build.PL
    ./Build
    ./Build test
    ./Build install

To work with binary Bam and BigWig files, see the ["advanced installation"](#advanced-installation) 
below for further guidance. Most scripts should fail gently with warnings about 
if required modules are missing.

Released versions may be obtained though the CPAN repository using your favorite 
package manager. 

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

There is a small INI-style configuration file, `.biotoolbox.cfg`, that may be written 
in your home directory, which can include paths to helper files and database 
configurations.

# ADVANCED INSTALLATION

This is a brief, advanced installation guide for getting a complete installation. I 
recommend using a simple CPAN package manager such as [cpanm](https://metacpan.org/pod/App::cpanminus).

## Locations

For privileged installations (requiring `root` access or `sudo` privilege) you probably 
already know what to do. You can use the `--sudo` or `-S` option to `cpanm`.

For home directory installations using the system perl, you should probably first install 
[local::lib](https://metacpan.org/pod/local::lib), and set the appropriate incantation 
in your `.bash_profile`, `.bashrc`, or other equivalent file as described. For example,

    curl -L https://cpanmin.us | perl - local::lib App::cpanminus \
    && echo 'eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"' >> ~/.bash_profile \
    && . ~/bash_profile

For home directory installations using a newer, modern perl (because many OS Perl 
installations are sadly out of date), please investigate installing your own Perl using 
[PerlBrew](https://perlbrew.pl). The older a Perl installation is, the more stuff has 
to be updated.

## External libraries

There are two external C libraries that are required for reading and writing Bam and 
BigWig files. Note that both [Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS) and 
[Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big) include `INSTALL.pl` scripts within 
their bundles that can compile these external libraries for you in a semi-automated 
control. Proceed here if you wish to have more control on where these are installed.

- [HTSlib](https://github.com/samtools/htslib)

   Follow the directions within for installation. [Version 1.5](https://github.com/samtools/htslib/releases/download/1.5/htslib-1.5.tar.bz2) 
   is known to work well, although newer versions may work as well too. By default, it 
   installs into `/usr/local`, or it may be set to another directory (`$HOME` for example) 
   by adding `--prefix=$HOME` option to the `configure` step. This may also be available 
   via other package managers.

- [libBigWig](https://github.com/dpryan79/libBigWig)

    Follow the directions within for installation. By default, it installs into 
    `/usr/local`. To change to a different location, manually edit the `Makefile`
    to change `prefix` to your desired location, and run `make && make install`.

## Perl modules

The following Perl packages should be explicitly installed. Most of these will 
bring along a number of dependencies (which in turn bring along more dependencies). In 
the end you will install dozens of packages, some of which are needed only for testing. 


- [Bio::Perl](https://metacpan.org/pod/Bio::Perl)

    The Bio::Perl package is a large bundle that brings along a number of extraneous 
    modules and bundled scripts, the vast majority of which is not needed by Bio::ToolBox.
    If you build this manually, or run `cpanm` with the `--interactive` option, you 
    can interactively choose what scripts to include or not include. By default, it 
    installs all additional scripts, potentially cluttering your `bin` folder.

- [Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS)

    It should be able to identify HTSLIB in standard library locations, like `/usr/local` 
    for example, on its own. For non-standard locations, specify the location of the 
    HTSlib path to `Build.PL` using the `--htslib` option. 

- [Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big)

    As with HTSlib, this should be able to identify the libBigWig library in standard 
    locations, but with non-standard locations, you may specify with the `--libbigwig` 
    option to `Build.PL`. 
    
    Note that on MacOS X, the bundle may not be linked properly to the shared library. 
    This will be apparent when the `./Build test` fails dramatically. You will have to 
    manually re-link the bundle to the shared library file with the following command.
    
        install_name_tool -change libBigWig.so /path/to/lib/libBigWig.so blib/arch/auto/Bio/DB/Big/Big.bundle
        
- [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager)

    This is highly recommended to get multi-cpu support. It can get a bit slow otherwise.

- [Set::IntervalTree](https://metacpan.org/pod/Set::IntervalTree)
- [DBD::SQLite](https://metacpan.org/pod/DBD::SQLite)

    If you plan on using BioPerl [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store) 
    databases for annotation then installing SQLite support is suggested. For larger, 
    shared databases, [DBD::mysql](https://metacpan.org/pod/DBD::mysql) is also supported.

An example of installing these with cpanm in your home directory is below.

    cpanm -L $HOME/perl5 Bio::Perl
    cpanm -L $HOME/perl5 --configure-args="--htslib $HOME" Bio::DB::HTS
    cpanm -L $HOME/perl5 --configure-args="--libbigwig $HOME" Bio::DB::Big
    cpanm -L $HOME/perl5 Parallel::ForkManager Set::IntervalTree Bio::ToolBox

## Legacy Perl modules

These are additional legacy Perl modules that are supported, but are not required or 
have been superseded by other modules.

- [Bio::DB::BigWig](https://metacpan.org/pod/Bio::DB::BigWig)
- [Bio::DB::BigBed](https://metacpan.org/pod/Bio::DB::BigBed)
- [Bio::DB::Sam](https://metacpan.org/pod/Bio::DB::Sam)
- [Bio::DB::USeq](https://metacpan.org/pod/Bio::DB::USeq)
- [Bio::Graphics::Wiggle](https://metacpan.org/pod/Bio::Graphics::Wiggle)

# AUTHOR

Timothy J. Parnell, PhD
Huntsman Cancer Institute
University of Utah
Salt Lake City, UT, 84112

# LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. 
