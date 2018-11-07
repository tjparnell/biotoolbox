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
- Works with any genomic annotation in GTF, GFF3, and various UCSC formats, including
refFlat, knownGene, genePred and genePredExt formats

The libraries provide a unified and integrated approach to analyses. 
In many cases, they provide an abstraction layer over a variety of 
different specialized BioPerl and related modules. Instead of 
writing numerous scripts specialized for each data format (wig, 
bigWig, Bam), one script can now work with virtually any data format. 

# INSTALLATION

Basic installation is simple with the standard [Module::Build](https://metacpan.org/pod/Module::Build) 
incantation. This will get you a minimal installation that will work with 
text files (BED, GFF, GTF, etc), but not binary files. 

    perl ./Build.PL
    ./Build
    ./Build test
    ./Build install

To work with binary Bam and BigWig files, see the ["advanced installation"](#advanced-installation) 
below for further guidance. Most scripts should fail gently with warnings 
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

    Included are three generic parsers for loading an entire genome-worth of 
    annotation into memory within a reasonably short amount of time. In each case 
    these load features into a SeqFeature object, 
    [Bio::ToolBox::SeqFeature](https://metacpan.org/pod/Bio::ToolBox::SeqFeature) by 
    default.

    - [Bio::ToolBox::parser::bed](https://metacpan.org/pod/Bio::ToolBox::parser::bed)

        This parses [BED file](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) and 
        related formats, including BED files with 3-12 columns 
        (BED3, BED6, BED12, and in between), [bedGraph](http://genome.ucsc.edu/FAQ/FAQformat.html#format1.8), 
        [narrowPeak](http://genome.ucsc.edu/FAQ/FAQformat.html#format12), and 
        [broadPeak](http://genome.ucsc.edu/FAQ/FAQformat.html#format13). For 
        proper BED12 files, transcripts are parsed with child subfeatures including exon 
        and CDS subfeatures.

    - [Bio::ToolBox::parser::gff](https://metacpan.org/pod/Bio::ToolBox::parser::gff)

        This parses both [GTF](http://mblab.wustl.edu/GTF22.html) and 
        [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) 
        file formats. Unlike most other GFF parsers 
        that work line-by-line only, this maintains parent and child hierarchical 
        relationships as parent feature and child subfeatures. To further maintain 
        control and reduce unnecessary parsing, unwanted feature types can be 
        selectively skipped.

    - [Bio::ToolBox::parser::ucsc](https://metacpan.org/pod/Bio::ToolBox::parser::ucsc)

        This parses various [UCSC file formats](http://genome.ucsc.edu/FAQ/FAQformat.html#format9), 
        including different refFlat, GenePred, and knownGene flavors. 
        Genes, transcripts, and exons are assembled into 
        hierarchical child-parent relationships as desired.

- [Bio::ToolBox::SeqFeature](https://metacpan.org/pod/Bio::ToolBox::SeqFeature)

    This is a fast, lean, simple object class for representing genomic features. 
    It supports, for the most part, the [Bio::SeqFreatureI](https://metacpan.org/pod/Bio::SeqFreatureI) 
    and [Bio::RangeI](https://metacpan.org/pod/Bio::RangeI) API 
    interface without the dependencies. It uses an unorthodox blessed-array object 
    structure, which provides measurable improvements in memory consumption and 
    speed when loading thousands of annotated SeqFeature objects (think hg19 or hg38 
    annotation). 

- [Bio::ToolBox::GeneTools](https://metacpan.org/pod/Bio::ToolBox::GeneTools)

    This is a collection of exportable functions for working with [Bio::SeqFeatureI](https://metacpan.org/pod/Bio::SeqFeatureI) 
    compliant objects representing genes and transcripts. It works with objects derived 
    from one of the ["Annotation parsers"](#annotation-parsers) or a 
    [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store) 
    database. The functions make hard things easy, such as identifying whether a 
    transcript is coding or not (is it encoded in the `primary_tag` or `source_tag` or 
    GFF attribute or does it have `CDS` subfeatures?), or identify the alternative exons 
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

## Examples

The following are just a few examples of highlighted scripts and their usage for 
solutions to common situations. Refer to the scripts' documentation for details on 
options shown and not shown.

- Convert or extract gene annotation

    Gene annotation from [UCSC](http://genome.ucsc.edu) frequently come in UCSC 
    [specific formats](http://genome.ucsc.edu/FAQ/FAQformat.html#format9), including
    refFlat and genePred, whereas gene annotation from [Ensembl](http://ensembl.org) 
    frequently come in GTF or GFF3 formats. At some point, some tool will specifically 
    need annotation in a different format than what you have (since most of 
    bioinformatics seem to be converting from one format to another), 
    [get_features.pl](https://metacpan.org/pod/get_features.pl) can help here.
    
    Simple conversion:
    
         get_features.pl --in file.gff3.gz --refflat
    
    Extract only protein-coding genes, collapsing alternate transcripts into one 
    combined meta-transcript, as a GTF:
    
         get_features.pl --in file.gff3.gz --feature mRNA --collapse --gtf
    
    Pull specific biotype `lincRNA` transcripts from standard chromosomes:
    
         get_features.pl --in file.gff3 --feature transcript --biotype lincRNA \
         --chrskip 'contig|scaffold|unmapped' --gtf
    
    Pull annotation from a [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store)
    SQLite database for use in downstream applications:
    
         get_features.pl --db annotation.sqlite --feature gene --out genes.txt 

- Extract specific gene regions

    Not all parts of genes that you might be interested in are explicitly defined 
    or encoded into a gene annotation table; rather, they are inferred. The 
    [get_gene_regions.pl](https://metacpan.org/pod/get_gene_regions.pl) script can 
    extract these inferred regions.
    
    Extract only the alternate exons from multi-transcript genes:
    
        get_gene_regions.pl --region altExon --in genes.gtf --out altExons
    
    Extract all the transcription start sites of protein coding transcripts, but only 
    report those of alternate transcripts once if they're within 200 bp, and 
    expand the site by 200 bp in both directions, as a bed file:
    
        get_gene_regions.pl --region tss --in genes.gtf --bed --out tss200.bed \
        --feature gene --transcript mRNA --unique --slop 200 --start -200 --end 200 

- Generate wig file representation of a bam file

    Lots of existing programs can generate a coverage file from a bam file, but 
    [bam2wig.pl](https://metacpan.org/pod/bam2wig.pl) will generate wig files in 
    every which way imaginable. Some common scenarios include:
    
    ChIPSeq normalized fragment coverage with empirical fragment size determination (single-end):
    
        bam2wig.pl --extend --shift --model --rpm --bw --in file1.bam --out file1.bw
    
    ChIPSeq normalized fragment coverage with explicit extension, explicit scaling factor, 
    averaged over 3 replicates, excluding repeat elements and certain chromosomes:
    
        bam2wig.pl --extend --extval 200 --bw --out file.bw --rpm --mean \
        --scale 0.0133 --blacklist repeats.bed.gz --chrskip 'chrM|contig|scaffold' \
        --cpu 12 replicate1.bam replicate2.bam replicate3.bam 
    
    ATACSeq cut sites, where cut sites are represented by 50 bp pileup centered over the 
    cut site at the 5' end:
    
        bam2wig.pl --shift --shiftval -25 --extend --extval 50 --in file1.bam
    
    RNASeq stranded coverage, where 2 wig files are generated for each strand, 
    ignoring multi-mappers, with output basename of 'experiment':
    
        bam2wig.pl --span --splice --strand --nosecondary --rpm --bw \
        --in file1.bam --out experiment
    
    Simple point representation of all fragments for counting:
    
        bam2wig.pl --start --bw --in file1.bam --out file1.bw

- Simple data collection

    Sometimes you just need simple, summarized scores over a list of genes or regions. 
    The [get_datasets.pl](https://metacpan.org/pod/get_datasets.pl) script can collect 
    data summarized data from bigWig, bigBed, Bam, USeq, and more.
    
    Collect mean values over Bed intervals:
    
        get_datasets.pl --in file.bed --method mean score1.bw score2.bw score3.bw
    
    In general, [bam2wig.pl](https://metacpan.org/pod/bam2wig.pl) has much better 
    controls for filtering and scoring bam alignments, but in a pinch 
    [get_datasets.pl](https://metacpan.org/pod/get_datasets.pl)
    will also work with bam files. A method of `mean` or `median` will work with 
    raw coverage, while `count`, name count (`ncount`), or precise count (`pcount`) 
    will work with alignment counts.
    
        get_datasets.pl --in file.bed --method ncount file1.bam file2.bam
    
    Collect mean sense RNASeq coverage from a collection of stranded bigWig files in 
    a folder (as shown above) over transcript exons:
    
        get_datasets.pl --in annotation.gtf --feature transcript --subfeature exon \
        --strand sense --ddb /path/to/stranded/bigWigs/ --data experiment --out file.txt
    
- Collect data around a reference point

    Rather than collect a single score for region, you need to collect a range of scores 
    relative to a single reference point, for example a transcription start site.
    The [get_relative_data.pl](https://metacpan.org/pod/get_relative_data.pl) script will 
    collect data in a defined number of windows on both sides of a specified reference 
    point.
    
    Collect mean values in twenty 500 bp windows (+/- 10 kb) from the transcription 
    start site:
    
        get_relative_data.pl --in mRNA.gtf --method mean --win 500 --num 20 \
        --pos 5 --data scores.bw --out mRNA_TSS_10kb.txt

- Collect profile data across an region

    Instead of collecting a single score for a region, you need to generate a profile 
    across a region. The [get_binned_data.pl](https://metacpan.org/pod/get_binned_data.pl)
    script will collect scores in defined fractional windows across a region, setting 
    length to each region an artificial 1000 bp. 
    
    To collect an average profile expression across transcripts, excluding introns, 
    and generate a summary (average profile):
    
        get_binned_data.pl --in mRNA.gtf --method mean --subfeature exon --bins 50 \
        --strand sense --ddb /path/to/stranded/bigWigs/ --data experiment \
        --sum --out mRNA_sense_profile.txt

- Manipulate collected dataset file

    Inevitably, you need to manipulate your data files: add or remove columns, perform 
    mathematical functions on columns, etc. You could open the file in a spreadsheet 
    application, or in an R session, or get creative with `cut`, `sed`, or `awk` commands 
    in the terminal, or just simply run 
    [manipulate_datasets.pl](https://metacpan.org/pod/manipulate_datasets.pl).
    
        manipulate_datasets.pl file.txt
    
    This script is designed to run interactively, allowing you to choose functions via 
    letter options from a menu. Alternatively, functions can be specified on the 
    command line for a programmatic approach in a `bash` script, for example. Note that 
    column indexes are numbered from 0.
    
        manipulate_datasets.pl --in file.txt --func multiply --index 5 --target 0.0133
    

# ADVANCED INSTALLATION

This is a brief, advanced installation guide for getting a complete installation. I 
recommend using a simple CPAN package manager such as [cpanm](https://metacpan.org/pod/App::cpanminus). 
Note that in the following example commands, `cpanm` is given a Perl module name, which 
is used to query [CPAN](https://metacpan.org), but it can also easily take a URL or a 
downloaded archive file. 

Note that `Make` and compilation tools, e.g. `GCC` or `clang`, are required. Most linux 
distributions have these available as an optional install if they're not already 
available. On MacOS, install the Xcode Command Line Tools. 

## Locations

For privileged installations (requiring `root` access or `sudo` privilege) you probably 
already know what to do. You can use the `--sudo` or `-S` option to `cpanm`. Note that 
installing lots of packages in the vendor system perl is generally not recommended, as 
it could interfere with other vital OS functions. It's best to use one of the other two 
methods.

For home directory installations using the system perl, you should first install 
[local::lib](https://metacpan.org/pod/local::lib), and set the appropriate incantation 
in your `.bash_profile`, `.bashrc`, or other equivalent file as described. This can 
also be used for targeted, standalone installations; adjust accordingly. For example,

    curl -L https://cpanmin.us | perl - -l $HOME/perl5 local::lib App::cpanminus \
    && echo 'eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"' >> ~/.bash_profile \
    && . ~/.bash_profile

For home directory installations using a newer, modern perl (because many vendor OS Perl 
installations are sadly out of date), please investigate installing your own Perl using 
[PerlBrew](https://perlbrew.pl), which makes installing newer Perl versions mostly 
painless and hands off; it easily manages multiple side-by-side installations as well as 
multiple `local::lib` installations, in case you want to isolate packages. 

## External libraries

There are two external C libraries that are required for reading Bam and BigWig files. 
Note that both Perl modules [Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS) and 
[Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big) include `INSTALL.pl` scripts within 
their bundles that can compile these external libraries for you in a semi-automated 
control. Proceed here if you wish to have more control over what and where these are 
installed.

- [HTSlib](https://github.com/samtools/htslib)

   Follow the directions within for installation. [Version 1.8](https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2) 
   is known to work well, although newer versions should work too. By default, it 
   installs into `/usr/local`, or it may be set to another directory (`$HOME` for example) 
   by adding `--prefix=$HOME` option to the `configure` step. This may also be available 
   via OS or other package managers.

- [libBigWig](https://github.com/dpryan79/libBigWig)

    Follow the directions within for installation. By default, it installs into 
    `/usr/local`. To change to a different location, manually edit the `Makefile`
    to change `prefix` to your desired location, and run `make && make install`.

## Perl modules

The following Perl packages should be explicitly installed. Most of these will 
bring along a number of dependencies (which in turn bring along more dependencies). In 
the end you will have installed dozens of packages. 


- [Module::Build](https://metacpan.org/pod/Module::Build)

    This may or may not need to be installed, depending on the age of your Perl 
    installation and how much else has been installed. It was part of the standard 
    Perl distribution up until version 5.21. 

- [Bio::Perl](https://metacpan.org/pod/Bio::Perl)

    The Bio::Perl package is a large bundle that brings along a number of extraneous 
    modules and bundled scripts, the vast majority of which is not needed by Bio::ToolBox.
    This is required by Bio::DB::HTS, local annotation SQLite databases, and any of the 
    legacy adapters. 
    
    If you build this manually, or run `cpanm` with the `--interactive` option, you 
    can interactively choose what scripts to include or not include. By default, it 
    installs all additional scripts. None of the included scripts are required by 
    BioToolBox.

- [Bio::DB::HTS](https://metacpan.org/pod/Bio::DB::HTS)

    This provides a perl interface to the HTSlib library for working with Bam files.
    It should be able to identify HTSLIB in standard library locations, such as `/usr/local` 
    for example, on its own. For non-standard locations, specify the location of the 
    HTSlib path to `Build.PL` using the `--htslib` option. 

- [Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big)

    This provides a perl interface to the UCSC-style bigWig and bigBed formats. 
    As with HTSlib, this should be able to identify the libBigWig library in standard 
    locations, but with non-standard locations, you may specify the path with the 
    `--libbigwig` option to `Build.PL`. 
    
    *NOTE*: The distribution from CPAN will install dozens of unnecessary modules for 
    remote URL testing. You may be better off installing from [source](https://github.com/Ensembl/Bio-DB-Big/archive/master.zip).

- [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager)

    This is highly recommended to get multi-cpu support for some of the data collection 
    scripts, which can otherwise get slow with a single thread.

- [Set::IntervalTree](https://metacpan.org/pod/Set::IntervalTree)

    This is necessary for optional functionality for a few scripts.

- [DBD::SQLite](https://metacpan.org/pod/DBD::SQLite)

    If you plan on using BioPerl [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store) 
    databases for annotation, then installing SQLite support is suggested. For larger, 
    shared databases, [DBD::mysql](https://metacpan.org/pod/DBD::mysql) is also supported.

An example of installing these Perl modules with `cpanm` in your home directory is below.
Adjust accordingly.

    cpanm -L $HOME/perl5 Module::Build Bio::Perl
    cpanm -L $HOME/perl5 --configure-args="--htslib $HOME" Bio::DB::HTS
    cpanm -L $HOME/perl5 --configure-args="--libbigwig $HOME" Bio::DB::Big
    cpanm -L $HOME/perl5 Parallel::ForkManager Set::IntervalTree Bio::ToolBox

## External applications

Some programs, notably [bam2wig.pl](https://metacpan.org/pod/distribution/Bio-ToolBox/scripts/bam2wig.pl) 
among others, requires external UCSC utilities for converting wig files to bigWig. You may 
download these from the UCSC Genome Browser utilities section for either 
[linux](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/) or 
[MacOS](http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/). Copy them to your 
`bin` directory in your `PATH`, for example `$HOME/bin`, `$HOME/perl5/bin`, or 
`/usr/local/bin`. Be sure to make them executable by running `chmod +x` on each file.

- wigToBigWig
- bedGraphToBigWig
- bedToBigBed

An example for downloading on linux:

    for name in wigToBigWig bedGraphToBigWig bedToBigBed; \
    do curl http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$name > $HOME/bin/$name \
    && chmod +x $HOME/bin/$name; done;


## Legacy Perl modules

These are additional legacy Perl modules that are supported (for example, if you still 
have a [GBrowse](http://gmod.org/wiki/GBrowse) installation), but are either not required 
or have been superseded by other modules. 

- [Bio::DB::BigWig](https://metacpan.org/pod/Bio::DB::BigWig)
- [Bio::DB::BigBed](https://metacpan.org/pod/Bio::DB::BigBed)
- [Bio::DB::Sam](https://metacpan.org/pod/Bio::DB::Sam)
- [Bio::DB::USeq](https://metacpan.org/pod/Bio::DB::USeq)
- [Bio::Graphics::Wiggle](https://metacpan.org/pod/Bio::Graphics::Wiggle)

## MacOS specific notes

While Macs have a Unix-compatible command-line environment, there are a few issues 
and solutions that I have encountered that may be useful to someone.

- Install XCode command line tools

    You don't need the full blown XCode installed, just the command line tools. 
    Running the following command in terminal will prompt to install them.
    
        xcode-select --install

- Installing your own Perl

    Apple is generally a little slow in updating their Perl compared to the latest 
    available versions, and it is compiled for backwards compatibility with 32-bit 
    `i386` processors (!), at least as of High Sierra (10.13). Some of the errors below will go away if you compile your 
    own Perl, but your success may vary.
    
    While using [PerlBrew](https://perlbrew.pl) generally works well in most cases, I 
    have found the recommendations described [here](https://karl.kornel.us/2015/12/perl-osx-1011-warnings/), 
    with modifications, work well for me. Here is my example for installing version 5.26 
    on a MacOS 10.13 machine.
    
        perlbrew install 5.26.2 --as 5.26 --thread --multi --64all -j 8 \
        -Dccflags="-mmacosx-version=10.13" -Dccdlflags="-mmacosx-version-min=10.13" \
        -Dldflags="-mmacosx-version-min=10.13" -Dlddflags="-mmacosx-version-min=10.13"

    You may adjust options as necessary.

- Linking errors

    When linking Perl modules with XS code (compiled C extensions), especially when using 
    the system Perl, you may see the following errors.
    
        ld: warning: object file was built for newer OSX version (10.13) than being linked (10.4)
    
    This is due to Apple compiling their system Perl with far-reaching backwards 
    compatibility; the Perl binary was compiled for both `i386` and `x86_64`, but in 
    all likelihood your XS was compiled only for `x86_64`. In some cases, this is a 
    harmless error; in other cases, it's a deal breaker. The best solution is to 
    install your own Perl.

- rpath errors

    This is especially notable with the [Bio::DB::Big](https://metacpan.org/pod/Bio::DB::Big) installation as described 
    above, where the `./Build test` fails dramatically because a shared library can not 
    be loaded by the bundle, usually with an error message including this:
    
        Reason: unsafe use of relative rpath libBigWig.so in blib/arch/auto/Bio/DB/Big/Big.bundle with restricted binary
    
    The solution is to manually re-link the bundle to the shared library file with the 
    following command. See this [link](https://stackoverflow.com/questions/33275605/el-capitan-perl-dbd-unsafe-use-of-relative-path) 
    for the source of the  solution.
    
        install_name_tool -change libBigWig.so /path/to/lib/libBigWig.so blib/arch/auto/Bio/DB/Big/Big.bundle

- DB_File errors

    There are [reports of issues](https://github.com/bioperl/bioperl-live/issues/267) 
    regarding certain BioPerl modules that rely on the Berkley database module 
    [DB_File](https://metacpan.org/pod/DB_File). This appears to stem from an issue with 
    the Apple-supplied library in High Sierra as described 
    [here](https://discussions.apple.com/thread/8125401). The best solution is to 
    install your own `berkley-db` library. 
    
    For BioToolBox users, the biggest effect appears to be exceptionally long times 
    during `Build` tests, specifically file `04.DB.t` that uses the in-memory database 
    adapater (maybe 20-30 seconds instead of 1), and excruciatingly long 
    [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store) 
    database builds (possibly days or weeks, I give up). 

# AUTHOR

	Timothy J. Parnell, PhD
	Huntsman Cancer Institute
	University of Utah
	Salt Lake City, UT, 84112

# LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. For details, see the
full text of the license in the file LICENSE.

This package is distributed in the hope that it will be useful, but it
is provided "as is" and without any express or implied warranties. For
details, see the full text of the license in the file LICENSE.
