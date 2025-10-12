# Bio::ToolBox - get\_binned\_data

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## get\_binned\_data.pl

A program to collect data in bins across a list of features.

## SYNOPSIS

    get_binned_data.pl [--options] --in <filename> --out <filename>
     
    get_binned_data.pl [--options] -i <filename> <data1> <data2...>
     
     Options for data files:
     -i --in <filename>                  input file: txt bed gff gtf refFlat ucsc
     -o --out <filename>                 optional output file, default overwrite 
     
     Options for new files:
     -d --db <name>                      annotation database: mysql sqlite
     -f --feature <type>                 one or more feature types from db or gff
     
     Options for data collection:
     -D --ddb <name|file>                data or BigWigSet database
     -a --data <dataset|filename>        data from which to collect: bw bam etc
     -m --method [mean|median|stddev|    statistical method for collecting data
           min|max|range|sum|count|      default mean
           pcount|ncount]
     -t --strand [all|sense|antisense]   strand of data relative to feature (all)
     -u --subfeature [exon|cds|          collect over gene subfeatures 
           5p_utr|3p_utr] 
     --long                              collect each window independently
     -r --format <integer>               number of decimal places for numbers
     --mapq <integer>                    minimum map quality of counted alignments
     
     Bin specification:
     -b --bins <integer>                 number of bins feature is divided (10)
     -x --ext <integer>                  number of extended bind outside feature
     -X --extsize <integer>              size of extended bins
     --min <integer>                     minimum size of feature to divide
     
     Post-processing:
     -U --sum                            generate summary file
     --smooth                            smoothen sparse data
     
     General options:
     -g --groups                         write columns group index file for plotting
     -z --gz                             compress output file
     -c --cpu <integer>                  number of threads, default 4
     --noparse                           do not parse input file into SeqFeatures
     -v --version                        print version and exit
     -h --help                           show extended documentation

## OPTIONS

The command line flags and descriptions:

### Options for data files

- --in &lt;filename>

    Specify an input file containing either a list of database features or 
    genomic coordinates for which to collect data. Any tab-delimited text 
    file with recognizable headers is supported. Gene annotation file 
    formats are also supported, including bed, gtf, gff3, refFlat, and 
    UCSC native formats such as gene prediction tables are all supported. 
    Gene annotation files will be parsed as sequence features. 
    Files may be gzipped compressed.

- --out &lt;filename>

    Specify the output file name. Default is to overwrite the input text 
    file. Required if generating a new file from a database.

### Options for new files

- --db &lt;name>

    Specify the name of a [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) annotation database 
    from which gene or feature annotation may be collected rather than providing 
    an input file. The name may be that of a MySQL database or a SQLite file. 

- --feature &lt;type | type:source | alias>,...

    Specify the type of feature from which to collect values. This is required 
    only for new feature tables. Three types of values may be passed: the 
    feature type, feature type and source expressed as 'type:source'. 
    More than one feature may be included as a comma-delimited list (no spaces). 

### Options for data collection

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

- --strand \[all|sense|antisense\]

    Specify whether stranded data should be collected. Three values are 
    allowed: all datasets should be collected (default), only sense 
    datasets, or only antisense datasets should be collected.

- --subfeature \[ exon | cds | 5p\_utr | 3p\_utr \]

    Optionally specify the type of subfeature to collect from, rather than 
    the entire gene. If the parent feature is gene and the subfeature is exon, 
    then all transcripts of the gene will be collapsed. The other subfeatures 
    (cds, 5p\_utr, and 3p\_utr) will not work with gene features but only with 
    coding mRNA transcripts. Note that the long option is incompatible. 
    Default is null. 

- --exons

    Legacy option for specifying --subfeature exon.

- --long

    Indicate that data should be collected independently for each long 
    window. This may be enabled automatically if the sum of the entire 
    window length passes a predefined threshold. The default for 'short' 
    windows is to collect all of the point data from the dataset first, and 
    then divide the results into the different windows. Datasets consisting 
    of "long" features, for example long alignments, may be counted more 
    than once in long mode when they span multiple windows. Not compatible 
    when subfeatures are enabled.

- --format &lt;integer>

    Specify the number of decimal positions to format the collected scores. 
    Default is not to format, often leading to more than the intended 
    significant digits.

- --mapq &lt;integer>

	Specify the minimum mapping quality of alignments to be considered when
	counting from a Bam file. Default is 0, which will include all alignments,
	including multi-mapping (typically MAPQ of 0). Set to an integer in range
	of 0..255. Only affects count methods, including `count`, `ncount`, and
	`pcount`. Other methods involving coverage, e.g. `mean`, do not filter
	alignments.

### Bin specification

- --bins &lt;integer>

    Specify the number of bins that will be generated over the length 
    of the feature. The size of the feature is a percentage of the 
    feature length. The default number is 10, which results in bins of 
    size equal to 10% of the feature length. 

- --ext &lt;integer>

    Specify the number of extended bins on either side of the feature. 
    The bins are of the same size as determined by the feature 
    length and the --bins value. The default is 0. 

- --extsize &lt;integer>

    Specify the exact bin size in bp of the extended bins rather than
    using a percentage of feature of length.

- --min &lt;integer>

    Specify the minimum feature size to be averaged. Features with a
    length below this value will not be skipped (all bins will have
    null values). This is to avoid having bin sizes below the average 
    microarray tiling distance. The default is undefined (no limit).

### Post-processing

- --sum

    Indicate that the data should be averaged across all features at
    each position, suitable for graphing. A separate text file will be
    written with the suffix '\_summed' with the averaged data. The default 
    is false.

- --smooth

    Indicate that windows without values should (not) be interpolated
    from neighboring values. The default is false.

### General options

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
    default is 4. Disable multi-threaded execution by setting to 1. 

- --noparse

    Prevent input annotation files from being automatically parsed into sequence 
    features. Coordinates will be used as is and new data columns will be appended 
    to the input file. 

- --version

    Print the version number.

- --help

    This help text.

## DESCRIPTION

This program will collect data across a gene or feature body into numerous 
percentile bins. It is used to determine if there is a spatial 
distribution preference for the dataset over gene bodies. The number 
of bins may be specified as a command argument (default 10). Additionally, 
extra bins may be extended on either side of the gene (default 0 on either 
side). The bin size is determined as a percentage of gene length.

## EXAMPLES

These are some examples of some common scenarios for collecting data.

- Collect scores in intervals

    You want to collect the mean score from a bigWig file in 10% intervals 
    across each feature in a Bed file.

        get_binned_data.pl --data scores.bw --in input.bed

- Collect scores in intervals plus extended regions

    You want to collect the maximum score in 5% intervals across each each 
    feature as well as five 100 bp intervals outside of each interval.

        get_binned_data.pl --bins 20 --method max --ext 5 --extsize 100 --data \
        scores.bw --in input.txt

- Collect scores in intervals for genes

    You want to collect stranded alignment counts from a Bam file for genes 
    in an annotation database.

        get_binned_data.pl --db annotation --feature gene --strand sense \
        --method count --data alignments.bam --out gene_profile --sum
        

## AUTHOR

    Timothy J. Parnell, PhD
    Howard Hughes Medical Institute
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
