# NAME

get\_relative\_data.pl

A program to collect data in bins around a relative position.

# SYNOPSIS

get\_relative\_data.pl \[--options\] --in &lt;filename> --out &lt;filename>

get\_relative\_data.pl \[--options\] -i &lt;filename> &lt;data1> &lt;data2...>

    Options for data files:
    -i --in <filename>                  input file: txt bed gff gtf refFlat ucsc
    -o --out <filename>                 optional output file, default overwrite 
    
    Options for new files
    -d --db <name>                      annotation database: mysql sqlite
    -f --feature <type>                 one or more feature types from db or gff
    
    Options for data collection:
    -D --ddb <name|file>                data or BigWigSet database
    -a --data <dataset|filename>        data from which to collect: bw bam etc
    -m --method [mean|median|stddev|    statistical method for collecting data
              min|max|range|sum|count|   default mean
              pcount|ncount]
    -t --strand [all|sense|antisense]   strand of data relative to feature (all)
    --force_strand                      use the specified strand in input file
    --avoid                             avoid neighboring features
    --avtype [type,type,...]            alternative types of feature to avoid
    --long                              collect each window independently
    -r --format <integer>               number of decimal places for numbers
    
    Bin specification:
    -w --win <integer>                  size of windows, default 50 bp
    -n --num <integer>                  number of windows flanking reference, 20
    --up <integer>                        or number of windows upstream
    --down <integer>                      and number of windows downstream
    -p --pos [5|m|3|p]                  reference position, default 5'
    
    Post-processing:
    -U --sum                            generate summary file
    --smooth                            smoothen sparse data
    
    General Options:
    -g --groups                         write columns group index file for plotting
    -z --gz                             compress output file
    -c --cpu <integer>                  number of threads, default 4
    --noparse                           do not parse input file into SeqFeatures
    -v --version                        print version and exit
    -h --help                           show extended documentation

# OPTIONS

The command line flags and descriptions:

## Options for data files

- --in &lt;filename>

    Specify an input file containing either a list of database features or 
    genomic coordinates for which to collect data. Any tab-delimited text 
    file with recognizable headers is supported. Gene annotation file 
    formats are also supported, including bed, gtf, gff3, refFlat, and 
    UCSC native formats such as gene prediction tables are all supported. 
    Gene annotation files will be parsed as sequence features. 
    Files may be gzipped compressed.

- --out &lt;filename>

    Specify the output file name. Required for new files; otherwise, 
    input files will be overwritten unless specified.

## Options for new files

- --db &lt;name | filename>

    Specify the name of a [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) annotation database 
    from which gene or feature annotation may be derived. A database is 
    required for generating new data files with features. This option may 
    skipped when using coordinate information from an input file (e.g. BED 
    file), or when using an existing input file with the database indicated 
    in the metadata.  

- --feature \[type, type:source\]

    Specify the type of feature to map data around. The feature may be 
    listed either as GFF type or GFF type:source. The list 
    of features will be automatically generated from the database. 
    This is only required when an input file is not specified. 

## Options for data collection

- --ddb &lt;name | filename>

    If the data to be collected is from a second database that is separate 
    from the annotation database, provide the name of the data database here. 
    Typically, a second [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) or BigWigSet database 
    is provided here. 

- --data &lt;dataset\_name | filename>

    Provide the name of the dataset to collect the values. If no 
    dataset is specified on the command line, then the program will 
    interactively present a list of datasets from the data database to select. 

    The dataset may be a database file, including bigWig (.bw), 
    bigBed (.bb), or Bam alignment (.bam) files. The files may be local or 
    remote (specified with a http: or ftp: prefix).

    Alternatively, the dataset may be a feature type in a BioPerl 
    [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) or [Bio::DB::BigWigSet](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ABigWigSet) database. Provide 
    either the feature type or `type:source`. 

    More than one datasource may be provided; use multiple data options or list 
    the datasets at the end of the command.

- --method &lt;text>

    Specify the method for combining all of the dataset values within the 
    genomic region of the feature. Accepted values include:

    - mean (default)
    - median
    - sum
    - stddev  Standard deviation of the population (within the region)
    - min
    - max
    - range   Returns difference of max and min
    - count

        Counts the number of overlapping items.

    - pcount (precise count)

        Counts the number of items that precisely fall within the query 
        region. Partially overlapping are not counted.

    - ncount (name count)

        Counts unique names. Useful when spliced alignments overlap more 
        than one exon and you want to avoid double-counting.

- --strand \[sense|antisense|all\]

    Specify whether stranded data should be collected for each of the 
    datasets. Either sense or antisense (relative to the feature) data 
    may be collected. The default value is 'all', indicating all 
    data will be collected.

- --force\_strand

    For features that are not inherently stranded (strand value of 0)
    or that you want to impose a different strand, set this option when
    collecting stranded data. This will reassign the specified strand for
    each feature regardless of its original orientation. This requires the
    presence of a "strand" column in the input data file. This option only
    works with input file lists of database features, not defined genomic
    regions (e.g. BED files). Default is false.

- --avoid

    Indicate whether neighboring features should be avoided when calculating 
    values in a window. After collecting the data, each window is checked for 
    overlapping features; if the window overlaps another feature, no value 
    is reported for that window. This option requires using an annotation 
    database (--db option). This is useful to avoid scoring windows that 
    overlap a neighboring gene, for example. The default is false (return 
    all values regardless of overlap).

- --avtype \[type,type,...\]

    Provide a feature type (primary\_tag or primary\_tag:source) or a 
    comma-delimited list of types to be used when avoiding neighboring 
    features. The default is to avoid features of the same type as that 
    of the query, i.e. collecting data around a feature of type 'gene' 
    will avoid other 'gene' features. This option allows you to avoid 
    other features too, such as 'tRNA' or 'repeat'.

- --long

    Indicate that data should be collected independently for each long 
    window. This may be enabled automatically if the sum of the entire 
    window length passes a predefined threshold. The default for 'short' 
    windows is to collect all of the point data from the dataset first, and 
    then divide the results into the different windows. Datasets consisting 
    of "long" features, for example long alignments, may be counted more 
    than once in long mode when they span multiple windows.

- --format &lt;integer>

    Specify the number of decimal positions to format the collected scores. 
    Default is not to format, often leading to more than the intended 
    significant digits.

## Bin specification

- --win &lt;integer>

    Specify the window size. The default is 50 bp.

- --num &lt;integer>

    Specify the number of windows on either side of the feature position 
    (total number will be 2 x \[num\]). The default is 20, or 1 kb on either 
    side of the reference position if the default window size is used.

- --up &lt;integer>
- --down &lt;integer>

    Alternatively specify the exact number of windows upstream and 
    downstream of the reference position. If only one option is set, 
    then the other option is assumed to be zero. 

- --pos \[5|m|3\]

    Indicate the relative position of the feature around which the 
    data is mapped. Three values are accepted: "5" indicates the 
    5' prime end is used, "3" indicates the 3' end is used, and "m" 
    indicates the middle of the feature is used. The default is to 
    use the 5' end, or the start position of unstranded features. 

## Post-processing

- --(no)sum

    Indicate that the data should be averaged across all features at
    each position, suitable for graphing. A separate text file will 
    be written with the suffix '\_summed' with the averaged data. 
    Default is true (sum).

- --smooth

    Indicate that windows without values should (not) be interpolated
    from neighboring values. The default is false (nosmooth).

## General options

- --groups

    Optionally write a secondary file with the list of column group names and 
    their corresponding dataset group. This can be used to assist in designating 
    the metadata when plotting files, for example in R with pheatmap. The 
    file is named the output basename appended with `.groups.txt`.

- --gz

    Specify whether (or not) the output file should be compressed with gzip.

- --cpu &lt;integer>

    Specify the number of CPU cores to execute in parallel. This requires 
    the installation of Parallel::ForkManager. With support enabled, the 
    default is 2. Disable multi-threaded execution by setting to 1. 

- --noparse

    Prevent input annotation files from being automatically parsed into sequence 
    features. Coordinates will be used as is and new data columns will be appended 
    to the input file. 

- --version

    Print the version number.

- --help

    Display this help

# DESCRIPTION

This program will collect data around a relative coordinate of a genomic 
feature or region. The data is collected in a series of windows flanking the 
feature start (5' position for stranded features), end (3' position), or 
the midpoint position. The number and size of windows are specified via 
command line arguments, or the program will default to 20 windows on both 
sides of the relative position (40 total) of 50 bp size, corresponding to 
2 kb total (+/- 1 kb). Windows without a value may be interpolated 
(smoothed) from neigboring values, if available.

Stranded data may be collected. If the feature does not have an inherent 
strand, one may be specified to enforce stranded collection or a particular 
orientation. 

When features overlap, or the collection windows of one feature overlaps 
with another feature, then data may be ignored and not collected (--avoid).

# EXAMPLES

These are some examples of some common scenarios for collecting data.

- Collect scores in intervals around start

    You want to collect the mean score from a bigWig file in twenty 50 bp 
    intervals flanking the start position of each feature in Bed file.

        get_relative_data.pl --data scores.bw --in input.bed

- Collect scores in intervals around middle

    You want to collect median scores in 20 bp intervals extending 500 bp 
    from the midpoint of each feature.

        get_relative_data.pl --win 20 --num 25 --pos m --data scores.bw --in \
        input.txt

- Collect scores in intervals from annotation database

    You want to collect scores in intervals around the transcription start 
    site of genes in an annotation database, but also avoid intervals that 
    may overlap neighboring genes. You want to collect alignment counts 
    from a Bam file in a stranded fashion. 

        get_relative_data.pl --db annotation --feature gene --avoid --strand \
        sense --method count --data alignments.bam --out gene_tss
          

# AUTHOR

    Timothy J. Parnell, PhD
    Howard Hughes Medical Institute
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
