# NAME

get\_gene\_regions.pl

A program to collect specific, often un-annotated, regions from genes.

# SYNOPSIS

get\_gene\_regions.pl \[--options...\] --in &lt;filename> --out &lt;filename>

get\_gene\_regions.pl \[--options...\] --db &lt;text> --out &lt;filename>

    Source data:
    -i --in <filename>            input annotation: GFF3, GTF, genePred, etc
    -d --db <name | filename>     database: name, file.db, or file.sqlite
    
    Feature selection:
    -f --feature <type>           optionally specify gene type or type:source
    -t --transcript               specify the transcript type
         [all|mRNA|ncRNA|snRNA|
         snoRNA|tRNA|rRNA|miRNA|
         lincRNA|misc_RNA]
    -r --region                   specify the gene region to collect
         [tss|tts|cdsStart|cdsStop|
         splice|UTR|exon|
         collapsedExon|altExon|
         uncommonExon|commonExon|
         firstExon|lastExon|intron|
         collapsedIntron|altIntron|
         uncommonIntron|commonIntron|
         firstIntron|lastIntron]
    --gencode                     include only GENCODE tagged genes
    --biotype <regex>             include only specific biotype
    --tsl                         select transcript support level
         [best|best1|best2|best3|
         best4|best5|1|2|3|4|5|NA]
    -u --unique                   select only unique regions
    -l --slop <integer>           duplicate region if within X bp
    -K --chrskip <regex>          skip features from certain chromosomes
    
    Adjustments:
    -b --begin --start integer     specify adjustment to start coordinate
    -e --end --stop integer        specify adjustment to stop coordinate
    
    General options:
    --bed                         output as a bed6 format
    -o --out <filename>              specify output name
    -z --gz                          compress output
    -v --version                     print version and exit
    -h --help

# OPTIONS

The command line flags and descriptions:

## Source data

- --in &lt;filename>

    Provide a gene table or annotation file, including GTF, GFF, GFF3, UCSC 
    refFlat, UCSC genePred or genePredExt, or UCSC knownGene table. Files 
    may be gzipped.

- --db &lt;text>

    Specify the name of a `Bio::DB::SeqFeature::Store` annotation database 
    from which gene or feature annotation may be obtained. Only required if 
    an input gene table is not provided.

## Feature selection

- --feature &lt;type>

    Specify the parental gene feature type (`primary_tag`) or `type:source` when
    using a database. If not specified, a list of available types will be
    presented interactively to the user for selection. This is not relevant for
    GFF3 source files (all gene or transcript features are considered). This is 
    helpful when gene annotation from multiple sources are present in the same 
    database, e.g. refSeq and ensembl sources. More than one feature may be 
    included, either as a comma-delimited list or multiple options.

- --transcript &lt;type>

    Specify the transcript type (usually a gene subfeature) from which to  
    collect the regions. Multiple types may be specified as a comma-delimited 
    list, or 'all' may be specified. If not specified, an interactive list 
    will be presented from which the user may select. Available options include:

         all
         mRNA
         ncRNA
         snRNA
         snoRNA
         tRNA
         rRNA
         miRNA
         lincRNA
         misc_RNA
        

- --region &lt;region>

    Specify the type of region to retrieve. If not specified on the command 
    line, the list is presented interactively to the user for selection. The 
    possibilities are listed below.

        tss           The first base of transcription
        tts           The last base of transcription
        exon          The exons of each transcript
        collapsedExon The exons after collapsing all gene transcripts
        firstExon     The first exon of each transcript
        lastExon      The last exon of each transcript
        altExon       Exons unique to one of several transcripts from a gene
        uncommonExon  Exons shared by 2 or more but not all transcripts
        commonExon    Exons shared by all transcripts from a gene
        intron        Each intron (usually not defined in the GFF3)
        collapsedIntron Introns after collapsing all gene transcripts
        firstIntron   The first intron of each transcript
        lastIntron    The last intron of each transcript
        altIntron     Introns unique to one of several transcripts from a gene
        uncommonIntron Introns shared by 2 or more but not all transcripts
        commonIntron  Introns shared by all transcripts of a gene
        splice        The first and last base of each intron
        UTR           The untranslated regions of each coding transcript
        cdsStart      The first base of the CDS
        cdsStop       The last base of the CDS

- --gencode

    Boolean option to filter transcripts as part of the GENCODE specification. 
    These are marked in Ensembl GTF/GFF3 annotation files as the `tag` attribute 
    with value "basic". Typically, at least one transcript for every gene is 
    marked as part of the GENCODE set. Transcripts not marked as such usually 
    lack sufficient experimental evidence.

- --biotype &lt;regex&lt;gt> 

    Filter transcripts using the `transcript_biotype` or `biotype` 
    GTF/GFF3 attribute, typically found in Ensembl annotation files. Provide 
    a regex compatible string which must match the biotype value to keep the 
    transcripts. For example, to keep specify "miRNA" to keep all micro-RNA 
    transcripts. This works on a subfeature level as well, so that `gene` 
    may be specified as the feature to collect, and only the gene transcripts 
    belonging to the indicating biotype are retained.

- --tsl &lt;level>

    Filter transcripts on the Ensembl GTF/GFF3 attribute 'transcript\_support\_level', 
    which is described at [Ensembl TSL glossary entry](http://uswest.ensembl.org/info/website/glossary.html).
    Provide a level of support to filter. Values include: 

        1       All splice junctions supported by evidence
        2       Transcript flagged as suspect or only support from multiple ESTs
        3       Only support from single EST
        4       Best supporting EST is suspect
        5       No support
        best    Transcripts at the best (lowest) available level are taken
        best1   The word followed by a digit 1-5, indicating any transcript 
                at or better (lower) than the indicated level
        NA      Only transcripts without a level (NA) are retained.

- --unique

    Compare start and stop coordinates of each collected region from 
    each feature and remove duplicate regions. When the --slop option 
    is provided, only the start coordinate plus/minus the slop factor 
    is checked. 

- --slop &lt;integer>

    When identifying unique regions, specify the number of bp to 
    add and subtract to the start position (the slop or fudge factor) 
    of the regions when considering duplicates. Any other region 
    within this window will be considered a duplicate. Useful, for 
    example, when start sites of transcription are not precisely mapped, 
    but not useful with defined introns and exons. This does not take 
    into consideration transcripts from other genes, only the current 
    gene. The default is 0 (no sloppiness).

- --chrskip &lt;regex>

    Exclude features from the output whose sequence ID or chromosome matches 
    the provided regex-compatible string. Expressions should be quoted or 
    properly escaped on the command line. Examples might be 

        'chrM'
        'scaffold.+'
        'chr.+alt|chrUn.+|chr.+_random'

## Adjustments

- --start &lt;integer>
- --begin &lt;integer>
- --stop &lt;integer>
- --end &lt;integer>

    Optionally specify adjustment values to adjust the reported start and 
    end coordinates of the collected regions. A negative value is shifted 
    upstream (5' direction), and a positive value is shifted downstream.
    Adjustments are made relative to the feature's strand, such that 
    a start adjustment will always modify the feature's 5'end, either 
    the feature startpoint or endpoint, depending on its orientation. 

## General options

- --bed

    Automatically convert the output file to a BED file.

- --out &lt;filename>

    Specify the output filename.

- --gz

    Specify whether (or not) the output file should be compressed with gzip.

- --version

    Print the version number.

- --help

    Display this POD documentation.

# DESCRIPTION

This program will collect specific regions from annotated genes and/or 
transcripts. Often these regions are not explicitly defined in the 
source GFF3 annotation, necessitating a script to pull them out. These 
regions include the start and stop sites of transcription, introns, 
the splice sites (both 5' and 3'), exons, the first (5') or last (3') 
exons, or all alternate or common exons of genes with multiple 
transcripts. Importantly, unique regions may only be reported, 
especially important when a single gene may have multiple alternative 
transcripts. A slop factor is included for imprecise annotation.

The program will report the chromosome, start and stop coordinates, 
strand, name, and parent and transcript names for each region 
identified. The reported start and stop sites may be adjusted with 
modifiers. A standard biotoolbox data formatted text file is generated. 
This may be converted into a standard BED or GFF file using the 
appropriate biotoolbox scripts. The file may also be used directly in 
data collection. 

# AUTHOR

    Timothy J. Parnell, PhD
    Howard Hughes Medical Institute
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
