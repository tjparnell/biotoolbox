# Bio::ToolBox - ucsc_table2gff3

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## ucsc\_table2gff3.pl

A program to convert UCSC gene tables to GFF3 or GTF annotation.

## SYNOPSIS

     ucsc_table2gff3.pl --ftp <text> --db <text>
     
     ucsc_table2gff3.pl [--options] --table <filename>
    
    UCSC database options:
    -f --ftp [refgene|known|all]          specify what tables to retrieve from UCSC
              
    -d --db <text>                        UCSC database name: hg19,hg38,danRer7, etc
    -h --host <text>                      specify UCSC hostname
    
    Input file options:
    -t --table <filename>                 name of table, repeat or comma list
    -k --kgxref <filename>                kgXref file
    -c --chromo <filename>                chromosome file
    
    Conversion options:
    --source <text>                       source text, default UCSC
    --chr   | --nochr         (true)      include chromosomes in output
    --gene  | --nogene        (true)      assemble into genes
    --cds   | --nocds         (true)      include CDS subfeatures
    --utr   | --noutr         (false)     include UTR subfeatures
    --codon | --nocodon       (false)     include start and stop codons
    --share | --noshare       (true)      share subfeatures
    --name  | --noname        (false)     include name
    -g --gtf                              convert to GTF instead of GFF3
    
    General options:
    -z --gz                               compress output
    -v --version                          print version and exit
    -h --help                             show extended documentation

## OPTIONS

The command line flags and descriptions:

### UCSC database options

- --ftp \[refgene|known|all\]

    Request that the current indicated tables and supporting files be 
    downloaded from UCSC via FTP. Four different tables may be downloaded, 
    including _refGene_, and the UCSC _knownGene_ table (if available). 
    Specify all to download all tables. A comma delimited list may also 
    be provided.

- --db &lt;text>

    Specify the genome version database from which to download the requested 
    table files. See [http://genome.ucsc.edu/FAQ/FAQreleases.html](http://genome.ucsc.edu/FAQ/FAQreleases.html) for a 
    current list of available UCSC genomes. Examples included hg19, mm9, and 
    danRer7.

- --host &lt;text>

    Optionally provide the host FTP address for downloading the current 
    gene table files. The default is 'hgdownload.cse.ucsc.edu'.

### Input file options

- --table &lt;filename>

    Provide the name of a UCSC gene or gene prediction table. Tables known 
    to work include the _refGene_, _ensGene_, _xenoRefGene_, and UCSC 
    _knownGene_ tables. Both simple and extended gene prediction tables, as 
    well as refFlat tables are supported. The file may be gzipped. When 
    converting multiple tables, use this option repeatedly for each table. 
    The `--ftp` option is recommended over using this one.

- --kgxref &lt;filename>

    Optionally provide the name of the _kgXref_ file. This file 
    provides additional information for the UCSC _knownGene_ gene table.
    The file may be gzipped.

- --chromo &lt;filename>

    Optionally provide the name of the chromInfo text file. Chromosome 
    and/or scaffold features will then be written at the beginning of the 
    output GFF file (when processing a single table) or written as a 
    separate file (when processing multiple tables). The file may be gzipped.

### Conversion options

- --source &lt;text>

    Optionally provide the text to be used as the GFF source. The default is 
    automatically derived from the source table file name, if recognized, or 
    'UCSC' if not recognized.

- --(no)chr

    When downloading the current gene tables from UCSC using the `--ftp` 
    option, indicate whether (or not) to include the _chromInfo_ table. 
    The default is true. 

- --(no)gene

    Specify whether (or not) to assemble mRNA transcripts into genes. This 
    will create the canonical gene->mRNA->(exon,CDS) heirarchical 
    structure. Otherwise, mRNA transcripts are kept independent. The gene name, 
    when available, are always associated with transcripts through the Alias 
    tag. The default is true.

- --(no)cds

    Specify whether (or not) to include CDS features in the output GFF file. 
    The default is true.

- --(no)utr

    Specify whether (or not) to include three\_prime\_utr and five\_prime\_utr 
    features in the transcript heirarchy. If not defined, the GFF interpreter 
    must infer the UTRs from the CDS and exon features. The default is false.

- --(no)codon

    Specify whether (or not) to include start\_codon and stop\_codon features 
    in the transcript heirarchy. The default is false.

- --(no)share

    Specify whether exons, UTRs, and codons that are common between multiple 
    transcripts of the same gene may be shared in the GFF3. Otherwise, each 
    subfeature will be represented individually. This will reduce the size of 
    the GFF3 file at the expense of increased complexity. If your parser 
    cannot handle multiple parents, set this to --noshare. Due to the 
    possibility of multiple translation start sites, CDS features are never 
    shared. This will have no effect with GTF output. The default is true. 

- --(no)name

    Specify whether you want subfeatures, including exons, CDSs, UTRs, and 
    start and stop codons to have display names. In most cases, this 
    information is not necessary. This will have no effect with GTF output. 
    The default is false.

- --gtf

    Specify that a GTF (version 2.5) format file should be written instead of 
    GFF3. Yes, the name of the program says GFF3, but now we can output GTF 
    too, and changing the name of the program is too late now.

### General options

- --gz

    Specify whether the output file should be compressed with gzip.

- --version

    Print the version number.

- --help

    Display the POD documentation

## DESCRIPTION

This program will convert a UCSC gene or gene prediction table file into a
GFF3 (or optionally GTF) format file. It will build canonical 
gene->transcript->\[exon, CDS, UTR\] heirarchical structures. It will 
attempt to identify non-coding genesas to type using the gene name as inference. 
Various additional informational attributes may also be included with the gene 
and transcriptfeatures, which are derived from supporting table files.

Two table files are currently supported. Gene prediction tables, including 
_refGene_ and UCSC _knownGene_ are supported. Supporting tables include 
_kgXref_. 

Tables obtained from UCSC are typically in the extended GenePrediction 
format, although simple genePrediction and refFlat formats are also 
supported. See [http://genome.ucsc.edu/FAQ/FAQformat.html#format9](http://genome.ucsc.edu/FAQ/FAQformat.html#format9) regarding
UCSC gene prediction table formats. 

The latest table files may be automatically downloaded using FTP from 
UCSC or other host. Since these files are periodically updated, this may 
be the best option. Alternatively, individual files may be specified 
through command line options. Files may be obtained manually through FTP, 
HTTP, or the UCSC Table Browser. 

If provided, chromosome and/or scaffold features will be written as GFF3-style 
sequence-region pragmas (even for GTF files, just in case).

If you need to set up a database using UCSC annotation, you should first 
take a look at the BioToolBox script [db\_setup.pl](https://metacpan.org/pod/db_setup.pl), which provides a 
convenient automated database setup based on UCSC annotation.  

## AUTHOR

    Timothy J. Parnell, PhD
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
