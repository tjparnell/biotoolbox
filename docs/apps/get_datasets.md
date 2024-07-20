# NAME

get\_datasets.pl

A program to collect data for a list of features

# SYNOPSIS

get\_datasets.pl \[--options...\] &lt;filename>

get\_datasets.pl \[--options...\] --in &lt;filename> &lt;data1> &lt;data2...>

    Options for data files:
    -i --in <filename>                  input file: txt bed gff gtf refFlat ucsc
    -o --out <filename>                 optional output file, default overwrite 
    
    Options for new files:
    -d --db <name>                      annotation database: mysql sqlite
    -f --feature <type>                 one or more feature types from db or gff
    
    Options for feature "genome":
    --win <integer>                     size of windows across genome (500 bp)
    --step <integer>                    step size of windows across genome
    --chrskip <regex>                   regular expression to skip chromosomes
    --blacklist <filename>              file of intervals to skip (bed, gff, txt)
    --prefix <text>                     prefix text for naming windows
    
    Options for data collection:
    -D --ddb <name|file>                data or BigWigSet database
    -a --data <dataset|filename>        data from which to collect: bw bam etc
    -m --method [mean|median|stddev|    statistical method for collecting data
              min|max|range|sum|count|   default mean
              pcount|ncount]
    -t --strand [all|sense|antisense]   strand of data relative to feature (all)
    -u --subfeature [exon|cds|          collect over gene subfeatures 
          5p_utr|3p_utr|intron] 
    --force_strand                      use the specified strand in input file
    --fpkm [region|genome]              calculate FPKM using which total count
    --tpm                               calculate TPM values
    -r --format <integer>               number of decimal places for numbers
    --discard <number>                  discard features whose sum below threshold
    
    Adjustments to features:
    -x --extend <integer>               extend the feature in both directions
    -b --begin --start <integer>        adjust relative start coordinate
    -e --end --stop <integer>           adjust relative stop coordinate
    -p --pos [5|m|3|53|p]               relative position to adjust (default 5')
    --fstart=<decimal>                  adjust fractional start
    --fstop=<decimal>                   adjust fractional stop
    --limit <integer>                   minimum size to take fractional window
    
    General options:
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

    Specify the output file name. Required for new feature tables; optional for 
    current files. If this is argument is not specified then the input file is 
    overwritten.

## Options for new files

- --db &lt;name | filename>

    Specify the name of a `Bio::DB::SeqFeature::Store` annotation database 
    from which gene or feature annotation may be derived. A database is 
    required for generating new data files with features. This option may 
    skipped when using coordinate information from an input file (e.g. BED 
    file), or when using an existing input file with the database indicated 
    in the metadata.  

- --feature &lt;type | type:source | alias>,...

    Specify the type of feature from which to collect values. This is required 
    only for new feature tables. Three types of values may be passed: the 
    feature type, feature type and source expressed as 'type:source', or an 
    alias to one or more feature types. More than one feature may be included 
    as a comma-delimited list (no spaces). 

## Options for feature "genome"

- --feature genome

    To collect genomic intervals (or regions) simply specify 'genome' as 
    the feature type.

- --win &lt;integer>

    When generating a new genome interval list (feature type 'genome'), 
    optionally specify the window size.  

- --step &lt;integer>

    Optionally indicate the step size when generating a new list of intervals 
    across the genome. The default is equal to the window size.

- --chrskip &lt;regex>

    Provide a regular expression to skip certain chromosomes. Perl-based 
    regular expressions are employed. Expressions should be quoted or 
    properly escaped on the command line. Examples might be 

        'chrM'
        'scaffold.+'
        'chr.+alt|chrUn.+|chr.+_random'

- --blacklist &lt;file>

    Provide a file of genomic intervals to avoid. Examples might include 
    multi-copy repetitive elements, ribosomal RNA, or heterochromatic regions.
    The file should be any text file interpretable by [Bio::ToolBox::Data](https://metacpan.org/pod/Bio%3A%3AToolBox%3A%3AData) 
    with chromosome, start, and stop coordinates, including BED and GFF formats.

- --prefix &lt;text>

    Provide a text string to prefix the name of generated genomic windows. 
    Names will be appended with an incrementing, unformatted digit. 

## Options for data collection

- --ddb &lt;name>

    If the data to be collected is from a second database that is separate 
    from the annotation database, provide the name of the data database here. 
    Typically, a second [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) or BigWigSet database 
    is provided here. 

- --data &lt;type1,type2,type3&type4,...>
- --data &lt;file1,...>
- --data none

    Provide the name of the dataset to collect the values. Use this argument 
    repeatedly for each dataset to be collected. Two or more datasets may be
    merged into one by delimiting with an ampersand "&" (no spaces!). If no 
    dataset is specified on the command line, then the program will 
    interactively present a list of datasets from the database to select. 

    The dataset may be a feature type in a BioPerl [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) 
    or [Bio::DB::BigWigSet](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ABigWigSet) database. Provide either the feature type or 
    type:source. The feature may point to another data file whose path is 
    stored in the feature's attribute tag (for example a binary 
    Bio::Graphics::Wiggle `.wib` file, a bigWig file, or Bam file), or the 
    features' scores may be used in data collection.

    Alternatively, the dataset may be a database file, including bigWig (.bw), 
    bigBed (.bb), useq (.useq), or Bam alignment (.bam) files. The files may 
    be local or remote (specified with a http: or ftp: prefix).

    Note that counting Bam alignments is _very_ limited. All supplementary, 
    secondary, and marked duplicates are ignored. One end of the alignment 
    must be within the feature (or both ends with method `pcount`). Splices 
    and indels are ignored. Paired-end alignment strand is inferred from the
    first read. In some cases, pre-filtering the alignments or converting 
    to another format (bigWig or bigBed) may be preferable.

    To force the program to simply write out the list of collected features 
    without collecting data, provide the dataset name of "none".

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

- --strand \[all | sense | antisense\]

    Specify whether stranded data should be collected for each of the 
    datasets. Either sense or antisense (relative to the feature) data 
    may be collected. Note that strand is not supported with some 
    data files, including bigWig files (unless specified through a GFF3 feature 
    attribute or BigWigSet database) and Bam files (score coverage
    is not but count is). The default value is 'all', indicating all data 
    will be collected.  

- --force\_strand

    For features that are not inherently stranded (strand value of 0)
    or that you want to impose a different strand, set this option when
    collecting stranded data. This will reassign the specified strand for
    each feature regardless of its original orientation. This requires the
    presence of a "strand" column in the input data file. This option only
    works with input file lists of database features, not defined genomic
    regions (e.g. BED files). Default is false.

- --subfeature \[ exon | cds | 5p\_utr | 3p\_utr | intron \]

    Optionally specify the type of subfeature to collect from, rather than 
    the entire gene. If the parent feature is gene and the subfeature is exon
    or intron, then all transcripts of the gene will be collapsed. The other 
    subfeatures (cds, 5p\_utr, and 3p\_utr) will not work with gene features but 
    only with coding mRNA transcripts. Note that the options extend, start, stop, 
    fstart, and fstop are ignored. Default is null. 

- --exons

    Legacy option for specifying --subfeature exon.

- --fpkm \[region|genome\]

    Optionally indicate that counts should be converted to Fragments Per 
    Kilobase per Million mapped (FPKM). This is a method for normalizing 
    sequence read depth and is used with Bam (or optionally bigBed) files. 
    Two methods exist for normalizing: 

    - region 

        Uses the sum of counts over all input regions examined and ignores 
        non-counted reads. This is the traditional method of calculating FPKM 
        and should be used preferentially with genes.

    - genome

        Uses the sum of all reads across the genome, regardless of whether 
        it was counted in an input region or not. This might be used when a 
        more global normalization is needed.

    The region method is best used with RNASeq data and a complete gene 
    annotation table. The genome method is best used with partial annotation 
    tables or other Seq types, such as ChIPSeq. This option can only be used 
    with one of the count methods (count, ncount, pcount). A sum method may be 
    cautiously allowed if, for example, using bigWig point data. The FPKM values 
    are appended as additional columns in the output table.

- --tpm

    Calculate Transcripts Per Million, a normalization method analogous to FPKM 
    but less biased to sequencing depth. This uses explicitly the counts collected 
    over the input regions, and not the entire genome.

- --format &lt;integer>

    Specify the number of decimal positions to format the collected scores. 
    Default is not to format, often leading to more than the intended 
    significant digits.

- --discard &lt;number>

    Discard features whose sum of newly collected datasets are less than the 
    indicated value. This is intended as a time-saving feature when collecting 
    alignment counts over features or genomic windows, where some features are 
    expected to return a zero count. Note that this only checks the datasets 
    that were newly collected. For more advanced filtering, see 
    [manipulate\_datasets.pl](https://metacpan.org/pod/manipulate_datasets.pl).

## Adjustments to features

- --extend &lt;integer>

    Optionally specify the bp extension that will be added to both sides of the 
    feature's region.

- --start &lt;integer>
- --stop &lt;integer>
- --begin &lt;integer>
- --end &lt;integer>

    Optionally specify adjustment values to adjust the region to collect values 
    relative to the feature position defined by the `--pos` option (default is 
    the 5' position). A negative value is shifted upstream (5' direction), 
    and a positive value is shifted downstream. Adjustments are always made 
    relative to the feature's strand. Default value is 0 (no change).

- --pos \[5|m|3|53|p\]

    Indicate the relative position of the feature with which the 
    data is collected when combined with the "start" and "stop" or "fstart" 
    and "fstop" options. Five values are accepted: \`5\` indicates the 
    5' prime end is used, \`3\` indicates the 3' end is used, \`m\` 
    indicates the middle of the feature is used, \`p\` indicates a peak 
    summit position is used (narrowPeak input only), and \`53\` indicates that 
    both ends are modified, i.e. start modifies start and end modifies end 
    (strand relative). The default is to use the 5' end, or the start 
    position of unstranded features. 

- --fstart=&lt;number>
- --fstop=&lt;number>

    Optionally specify the fractional start and stop position of the region to 
    collect values as a function of the feature's length and relative to the 
    specified feature position defined by the `--pos` option (default is 5'). The 
    fraction should be presented as a decimal number, e.g. 0.25. Prefix a 
    negative sign to specify an upstream position. Default values are 0 (fstart)
    and 1 (fstop), or no change. 

- --limit &lt;integer>

    Optionally specify the minimum size limit for subfractionating a feature's 
    region. Used in combination with fstart and fstop to prevent taking a 
    subregion from a region too small to support it. The default is 10 bp.

## General options

- --gz

    Indicate whether the output file should (not) be compressed by gzip. 
    If compressed, the extension `.gz` is appended to the filename. If a compressed 
    file is opened, the compression status is preserved unless specified otherwise.

- --cpu &lt;integer>

    Specify the number of CPU cores to execute in parallel. This requires 
    the installation of [Parallel::ForkManager](https://metacpan.org/pod/Parallel%3A%3AForkManager). With support enabled, the 
    default is 4. Disable multi-threaded execution by setting to 1. 

- --noparse

    Prevent input annotation files from being automatically parsed into sequence 
    features. Coordinates will be used as is and new data columns will be appended 
    to the input file. 

- --version

    Print the version number.

- --help

    Display the POD documentation for this program.

# DESCRIPTION

This program will collect dataset values from a variety of sources, including 
features in a BioPerl Bio::DB::SeqFeature::Store database, binary wig files 
`.wib` loaded in a database using Bio::Graphics::Wiggle, bigWig files, 
bigBed files, Bam alignment files, or a Bio::DB::BigWigSet database. 

The values are collected for a list of known database features (genes, 
transcripts, etc.) or genomic regions (defined by chromosome, start, and 
stop). The list may be provided as an input file or generated as a new 
list from a database. Output data files may be reloaded for additional 
data collection.

At each feature or interval, multiple data points within the genomic segment 
are combined statistically and reported as a single value for the feature. 
The method for combining datapoints may be specified; the default method is 
the mean of all datapoints.

The coordinates of the features may be adjusted in numerous ways, including 
specifying a specific relative start and stop, a fractional start and stop, 
an extension to both start and stop, and specifying the relative position 
(5' or 3' or midpoint).

Stranded data may be collected, if the dataset supports stranded information. 
Also, two or more datasets may be combined and treated as one. Note that 
collecting stranded data may significantly slow down data collection.

# EXAMPLES

These are some examples of some common scenarios for collecting data.

- Simple mean scores

    You want to collect the mean score from a bigWig file for each feature 
    in a BED file of intervals.

        get_datasets.pl --in input.bed --data scores.bw

- Collect normalized counts

    You want to collect normalized read counts from multiple Bam files 
    for each feature in a BED file. This will count alignment names (safe for 
    paired-end alignments) over the intervals, and transform to Fragments (Reads) 
    Per Million, a depth-normalizing function based on the total number of 
    fragments counted in each dataset. 

        get_datasets.pl --in input.bed --method ncount --fpkm region *.bam 

- Collect stranded RNASeq data

    You have stranded RNASeq data, and you would like to determine the 
    expression level for all genes from an annotation file. Use the 
    `ncount` method to count alignment names to avoid double counting 
    alignments split over multiple exons.

        get_datasets.pl --in annotation.gtf --feature transcript --subfeature exon \
        --strand sense --method ncount --out expression.txt *.bam

- Restrict to specific region

    You have ChIPSeq enrichment scores in a bigWig file and you now want 
    to score just the transcription start site of known transcripts in a 
    [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) annotation database. Here you will 
    restrict to 500 bp flanking the TSS.

        get_datasets.pl --db annotation.sqlite --feature mRNA --start=-500 \
        --stop=500 --pos 5 --data scores.bw --out tss_scores.txt

- Avoid first and last 1 Kb of each interval

        get_datasets.pl --in file.bed --start=1000 \
        --stop=-1000 --pos 53 --data scores.bw --out file_scores.txt

- Count intervals

    You have identified all possible transcription factor binding sites in 
    the genome and put them in a bigBed file. Now you want to count how 
    many exist in each upstream region of each gene.

        get_datasets.pl --db annotation.gtf --feature gene --start=-5000 \
        --stop=0 --data tfbs.bb --method count --out tfbs_sums.txt

- Many datasets at once

    While you can provide multiple data source files as a space-delimited 
    list at the end of the command, you can also treat a folder or directory 
    of bigWig files as a special database, known as a BigWigSet. Each file 
    becomes a database feature, and you can interactively choose one or more 
    from which to collect. Each dataset is appended to the input file as a new 
    column. Provide the folder as a data database `--ddb` option.

        get_datasets.pl --in input.txt --ddb /path/to/bigwigset

- Stranded BigWig data

    You can generate stranded RNASeq coverage from a Bam file using the 
    BioToolBox script bam2wig.pl, which yields rnaseq\_f.bw and rnaseq\_r.bw 
    files. These are automatically interpreted as stranded datasets in a 
    BigWigSet folder context; see above.

        get_datasets.pl --in input.txt --strand sense \
        --ddb /path/to/rnaseq/bigwigset 

- Binned coverage across the genome

    You are interested in sequencing depth across the genome to look for 
    depleted regions. You count reads in 1 kb intervals across the genome.

        get_datasets.pl --db genome.fasta --feature genome --win 1000 \
        --data alignments.bam --method count --out coverage.txt

- Middle of feature

    You are interested in the maximum score in the central 50% of each 
    feature.

        get_datasets.pl --in input.txt --fstart=0.25 --fstop=0.75 \
        --data scores.bw 
        

# AUTHOR

    Timothy J. Parnell, PhD
    Howard Hughes Medical Institute
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
