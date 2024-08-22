# Bio::ToolBox

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## correlate\_position\_data.pl

A script to calculate correlations between two datasets along the length of a feature.

## SYNOPSIS

correlate\_position\_data.pl \[--options\] &lt;filename>

    Options for data files:
    -i --in <filename>               input file: txt bed etc
    -o --out <filename>              optional output file, default overwrite 
    -d --db <name>                   alternate annotation database
    
    Options for data sources
    -D --ddb <name|file>             data or BigWigSet database
    -r --ref <dataset|filename>      reference data: bw, name, etc
    -t --test <dataset|filename>     test data: bw, name, etc
    
    Options for correlating data
    --pval                           calculate P-value by ANOVA
    --shift                          determine optimal shift to match datasets
    --radius <integer>               radius in bp around reference point to calculate
    -p --pos [5|m|3]                 reference point to measure correlation (m)
    --norm [rank|sum]                normalization method between datasets
    --force_strand                   force an alternate strand
    
    General options:
    -c --cpu <interger>              number of threads (4)
    -z --gz                          compress output with gz
    -v --version                     print version and exit
    -h --help                        show extended documentation

## OPTIONS

The command line flags and descriptions:

### Options for data files

- --in &lt;filename>

    Specify the input file of features. The features may be comprised of 
    name and type, or chromosome, start, and stop. Strand information is 
    optional, but is respected if included. A feature list may be 
    generated with the program [get\_features.pl](https://metacpan.org/pod/get_features.pl). The file may be 
    compressed with gzip.

- --out &lt;filename>

    Specify the output filename. By default it rewrites the input file.

- --db &lt;name | filename>

    Specify the name of a [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) annotation database 
    from which gene or feature annotation may be derived. A database is 
    required for generating new data files with features. This option may 
    skipped when using coordinate information from an input file (e.g. BED 
    file), or when using an existing input file with the database indicated 
    in the metadata. 

### Options for data sources

- --ddb &lt;name | filename>

    If the data to be collected is from a second database that is separate 
    from the annotation database, provide the name of the data database here. 
    Typically, a second [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) or BigWigSet database 
    is provided here. 

- --ref &lt;type | filename>
- --test &lt;type | filename>

    Define both the reference and test datasets with which to compare and 
    correlate. These may be GFF type or name in a database or BigWigSet, or 
    they may be a BigWig or even Bam file. Both options are required. If 
    not provided, they may be interactively chosen from the database.

### Options for correlating data

- --pval

    Perform an ANOVA analysis between the test and reference data sets and 
    report a P-value. By default, this performs a dependent, parametric 
    ANOVA. This requires the [Statistic::ANOVA](https://metacpan.org/pod/Statistic%3A%3AANOVA) module to be installed. 
    Please refer to the module documentation for details. If your needs 
    require a change to the test, you may edit the parameters at the top 
    of this script. For convenience, the P-values are reported as -Log10(P) 
    transformed values. The default is false.

- --shift

    Optionally specify whether an optimal shift should be calculated that 
    would result in a better Pearson correlation value. The default is 
    false.

- --radius &lt;integer>

    Define the radius in basepairs around a reference point to determine 
    the window size for the correlation analysis. This value is required 
    when calculating an optimal shift (--shift option). The default is to 
    take the length of the feature as the window for calculating the 
    correlation.  

- --pos \[5|m|3\]

    Indicate the relative position of the feature to be used as the 
    reference point around which the window (determined by the radius 
    value) for collecting data will be centered. Three values are 
    accepted: "5" indicates the 5' prime end is used, "3" indicates the 
    3' end is used, and "m" indicates the middle of the feature is used. 
    The default is to use the midpoint. 

- --norm \[rank|sum\]

    Optionally define a method of normalizing the scores between the 
    reference and test data sets prior to calculating the correlation. 
    Two methods are currently supported: "rank" converts all values 
    to rank values (the mean rank is reported for identical values) 
    and essentially calculating a Spearman's rank correlation, while 
    "sum" scales all values so that the absolute sums are identical. 
    Normalization occurs after missing or zero values are interpolated. 
    The default is no normalization.

- --force\_strand

    If enabled, a strand orientation will be enforced when determining the 
    optimal shift. This does not affect the correlation calculation, only 
    the direction of the reported shift. This requires the presence of a 
    data column in the input file with strand information. The default is 
    no enforcement of strand.

### General options

- --gz

    Specify whether (or not) the output file should be compressed with gzip.

- --version

    Print the version number.

- --help

    Display this POD documentation.

## DESCRIPTION

This program will calculate statistics between the positioned scores of
two different datasets over a window from an annotated feature or
chromosomal segment. These statistics will help determine whether the
positions or distribution of scores across the window vary or underwent
a positional shift between a test and a reference dataset. For example,
if the enrichment of nucleosome signal from a ChIP experiment shifts in
genomic position, indicating a change in nucleosome position. 

Two statistics may be calculated. First, it will calculate a a Pearson
linear correlation coefficient (r value) between the datasets (default). 
Additionally, an ANOVA analysis may be performed between the datasets and 
generate a P-value. 

By default, the correlation is determined between the data points 
collected over the entire length of the feature. Alternatively, a 
radius and reference point (default is midpoint) may be provided 
that sets the window for collecting scores and calculating a correlation.

In general, to ensure a more reliable Pearson value, fragment ChIP or 
nucleosome coverage should be used rather than point (start or midpoint) 
data, as it will give more reliable results. Fragment coverage is more 
akin to smoothened data and gives better results than interpolated point 
data. 

Normalized read-depth data should be used when possible. If necessary, 
Values can be normalized using one of two methods. The values may be 
converted to rank positions (compare to Kendall's tau), or scaled such 
that the absolute sum values are equal (for example, when working with 
sequence tag read counts).

In addition to calculating a correlation coefficient, an optimal shift 
may also be calculated. This essentially shifts the data, 1 bp at a time, 
in order to identify a shift that would produce a higher correlation. In 
other words, what amount of movement to the left or right would make the 
test data look like the reference data? The window is shifted from -2 
radius to +2 radius relative to the reference point, and the highest 
correlation is reported along with the shift value that generated it. 

## AUTHOR

    Timothy J. Parnell, PhD
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
