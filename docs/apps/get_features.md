# Bio::ToolBox - get_features

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## get\_features.pl

A program to collect and filter annotated features from source files.

## SYNOPSIS

get\_features.pl --in &lt;filename> --out &lt;filename>

get\_features.pl --db &lt;name> --out &lt;filename>

    Source data:
    -d --db <name | filename>     database: name, file.db, or file.sqlite
    -i --in <filename>            input annotation: GFF3, GTF, genePred, etc
    
    Selection:
    -f --feature <type>           feature: gene, mRNA, transcript, etc
    -u --sub                      include subfeatures (true if gff, gtf, refFlat)
    
    Filter features:
    -l --list <filename>          file of feature IDs to keep
    -K --chrskip <regex>          skip features from certain chromosomes
    -x --exclude <tag=value>      exclude features with specific attribute value
    -n --include <tag=value>      include features with specific attribute value
    --biotype <regex>             include only specific transcript biotype
    --tsl [best|best1|best2|      specify minimum transcript support level 
           best3|best4|best5|       primarily Ensembl annotation 
           1|2|3|4|5|NA]  
    --gencode                     include only GENCODE tagged genes
    
    Adjustments:
    -b --start=<integer>          modify start positions
    -e --stop=<integer>           modify stop positions
    -p --pos [5|m|3|53|p]         relative position from which to modify
    --collapse                    collapse subfeatures from alt transcripts
    
    Report format options:
    -B --bed                      write BED6 (no --sub) or BED12 (--sub) format
    -G --gff                      write GFF3 format
    -g --gtf                      write GTF format
    -r --refflat                  write UCSC refFlat format
    -t --tag <text>               include specific GFF attributes in text output
    --coord                       include coordinates in text output
    
    General options:
    -o --out <filename>           output file name
    --sort                        sort output by genomic coordinates
    -z --gz                       compress output
    -Z --bgz                      bgzip compress output
    -v --version                  print version and exit
    -h --help                     show full documentation

## OPTIONS

The command line flags and descriptions:

### Source data

- --db &lt;text>

    Specify the name of a [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) annotation database 
    from which gene or feature annotation may be derived. A SQLite file 
    or a named relational database may be provided. Used as an alternate 
    to an input file.

- --in &lt;filename>

    Specify the filename of a gene annotation file, including GTF, GFF3, 
    or a UCSC-formatted file including, refFlat, genePred, or knownGene.
    The file may be gzip compressed. Used as an alternate to a database.

### Selection

- --feature &lt;type>

    Provide a feature type to collect. Typically, this corresponds to the 
    GFF type, or 3rd column in a GFF3 or GTF file. Examples include 
    `gene`, `mRNA`, or `transcript`. The default value for input files 
    is '`gene`'. For databases, an interactive list will be presented 
    from which one or more may be chosen.

- --sub

    Optionally include all child subfeatures in the output. For example, 
    transcript, CDS, and/or exon subfeatures of a gene. This option is 
    automatically enabled with GFF, GTF, or refFlat output; it may be 
    turned off with `--nosub`. With BED output, it will force a BED12 
    file to be written. It has no effect with standard text. 

### Filter features

- --list &lt;file>

    Provide a file containing a list of the feature IDs to keep from the 
    input annotation file. The file must have an `ID` or `Primary_ID` 
    column, otherwise an index will be requested interactively from the user. 
    Only those top features whose `ID`, `gene_id`, `transcript_id`, or 
    otherwise SeqFeature `primary_id` exactly matches one from the list 
    will be retained; `names` are not checked. Useful if you need 
    to filter based on external criteria.

- --chrskip &lt;regex>

    Exclude features from the output whose sequence ID or chromosome matches 
    the provided regex-compatible string. Expressions should be quoted or 
    properly escaped on the command line. Examples might be 

        'chrM'
        'scaffold.+'
        'chr.+alt|chrUn.+|chr.+_random'

- --exclude &lt;tag=value>

    Provide a GFF/GTF attribute tag on which to filter out the features 
    matching the indicated value. For example, to filter out protein 
    coding genes using `gene_biotype`, specify "gene\_biotype=protein\_coding". 
    The value is checked by regular expression and can be specified as such.
    This filter does not descend into subfeatures. More than one exclusion 
    tag may be specified with multiple options or a comma-delimited list. 

- --include &lt;tag=value>

    Provide a GFF/GTF attribute tag on which to filter for the features 
    matching the indicated value. For example, to filter for protein 
    coding genes using `gene_biotype`, specify "gene\_biotype=protein\_coding". 
    The value is checked by regular expression and can be specified as such.
    This filter does not descend into subfeatures. More than one inclusion 
    tag may be specified with multiple options or a comma-delimited list. 

- --biotype &lt;regex&lt;gt> 

    Filter transcripts using the `transcript_biotype` or `biotype` 
    GTF/GFF3 attribute, typically found in Ensembl annotation files. Provide 
    a regex compatible string which must match the biotype value to keep the 
    transcripts. For example, to keep specify "miRNA" to keep all micro-RNA 
    transcripts. This works on a subfeature level as well, so that `gene` 
    may be specified as the feature to collect, and only the gene transcripts 
    belonging to the indicating biotype are retained.

- --tsl &lt;level>

    Filter transcripts on the Ensembl GTF/GFF3 attribute `transcript_support_level`, 
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

- --gencode

    Boolean option to filter transcripts as part of the GENCODE specification. 
    These are marked in Ensembl GTF/GFF3 annotation files as the `tag` attribute 
    with value "basic". Typically, at least one transcript for every gene is 
    marked as part of the GENCODE set. Transcripts not marked as such usually 
    lack sufficient experimental evidence.

### Adjustments

- --start=&lt;integer>
- --stop=&lt;integer>

    Optionally specify adjustment values to adjust the reported start and 
    end coordinates of the collected regions. A negative value is shifted 
    upstream (towards the 5 prime end), and a positive value is shifted 
    downstream (towards the 3 prime end). Adjustments are made relative 
    to the indicated position (--pos option, below) based on the feature 
    strand. Adjustments are only allowed when writing standard BED6 or 
    standard text files.

- --pos \[ 5 | m | 3 | 53 | p \]

    Indicate the relative position from which both coordinate adjustments 
    are made. Both start and stop adjustments may be made from the respective 
    5 prime, 3 prime, or middle position as dictated by the feature's strand 
    value. Alternatively, specify '53' to indicate that the start adjustment 
    adjusts the 5 prime end and the stop adjustment adjusts the 3 prime end. 
    If the input file is `narrowPeak` format, adjustments can be made 
    relative to the summit or peak position by specifying 'p' (non-peak files 
    report midpoint). The default is '53'.

- --collapse

    Boolean option to collapse multiple alternate transcripts of a gene into 
    a single artificial transcript, merging overlapping exons and minimizing 
    introns, where appropriate. Genes without alternate transcripts are not 
    collapsed.

### Report format options

- --bed

    With subfeatures enabled, write a BED12 (12-column BED) file. 
    Otherwise, write a standard 6-column BED format file. 

- --gff

    Write a GFF version 3 (GFF3) format output file. Subfeatures are 
    automatically included and coordinate adjustments ignored.

- --gtf

    Write a GTF format (GFF version 2.2 or 2.5) output file. Subfeatures are 
    automatically included and coordinate adjustments ignored.

- --refflat 
- --ucsc

    Write a UCSC-style refFlat format (10 columns) gene annotation table. 
    Subfeatures are automatically included and coordinate adjustments ignored.

- --tag &lt;text>

    When writing a standard text file, optionally include additional 
    GFF/GTF attribute tags. Specify as a comma-delimited list or as 
    separate options.

- --coord

    When writing a standard text file, optionally include the chromosome, 
    start, stop, and strand coordinates. These are automatically included 
    in other formats. This is automatically included when adjusting 
    coordinate positions.

### General options

- --out &lt;filename>

    Specify the output file name. Default is the joined list of features. 

- --sort

    Sort the output file by genomic coordinates. Automatically enabled 
    when compressing with bgzip. This may require more memory.

- --gz

    Specify whether the output file should be compressed with gzip.

- --bgz

    Specify whether the output file should be compressed with block gzip 
    (bgzip) for tabix compatibility.

- --version

    Print the version number.

- --help

    Display this POD documentation.

## DESCRIPTION

This program will extract a list of features from a database or input 
annotation file and write them out to a file. Features may be selected 
using their feature type (the 3rd column in a GFF or GTF file). When 
selecting features from a database, types may be selected interactively 
from a presented list. Features may be filtered based on various 
GFF attributes typically found in Ensembl annotation, including 
`transcript_biotype`, `transcript_support_level`, and GENCODE basic 
tags. Features may also be filtered by chromosome. 

Collected features may be written to a variety of formats, including 
GFF3, GTF, refFlat, simple 6-column BED, or a simple text format. With 
GFF, GTF, and refFlat formats, subfeatures such as exons are automatically 
included (which may also be turned off). With a simple text format, 
the source database or parsed input file are recorded in the header 
metadata for use in subsequent programs. Coordinates may be optionally 
included in the text file, which preempts using parsed features in other 
tools. 

### Coordinate adjustments

Coordinates of the features may be adjusted as desired when writing to 
text or BED file formats. Adjustments may be made relative to either 
the 5 prime, 3 prime, both ends, or the feature midpoint. Positions 
are based on the feature strand. Use the following examples as a guide. 

- upstream 500 bp only

        get_features.pl --start=-500 --stop=-1 --pos 5

- 1 kb total around 5 prime end

        get_features.pl --start=-500 --stop=500 --pos 5

- last 500 bp of feature

        get_features.pl --start=-500 --pos 3

- middle 500 bp of feature

        get_features.pl --start=-250 --stop=250 --pos m

- entire feature plus 1 kb of flanking

        get_features.pl --start=-1000 --stop=1000 --pos 53

Note that positions are always in base coordinates, and the resulting regions 
may be 1 bp longer depending on whether the reference base was included or not.

## AUTHOR

    Timothy J. Parnell, PhD
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
