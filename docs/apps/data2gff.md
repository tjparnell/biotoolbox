# Bio::ToolBox - data2gff

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## data2gff.pl

A program to convert a generic data file to GFF format.

## SYNOPSIS

data2gff.pl \[--options...\] &lt;filename>

    File options:
    -i --in <filename>                    input file: txt
    -o --out <filename>                   output file name
    -H --noheader                         input file has no header row
    -0 --zero                             file is in 0-based coordinate system
    
    Column indices:
    -a --ask                              interactive selection of columns
    -c --chr <index>                      chromosome column
    -b --begin --start <index>            start coordinate column
    -e --end --stop <index>               stop coordinate column
    -s --score <index>                    score column
    -t --strand <index>                   strand column
    -n --name <text | index>              name column or base name text
    -d --id <index>                       primary ID column
    -g --tags <index,index,...>           zero or more columns for tag attributes
    -r --source <text | index>            source column or text
    -y --type <text | index>              type column or text
    
    General options:
    --unique                              make IDs unique
    --sort                                sort output by genomic coordinates
    -z --gz                               compress output file
    -Z --bgz                              bgzip compress output file
    -v --version                          print version and exit
    -h --help                             show extended documentation

## OPTIONS

The command line flags and descriptions:

### File options

- --in &lt;filename>

    Specify an input file containing either a list of database features or 
    genomic coordinates for which to convert to GFF format. The file should be a 
    tab-delimited text file, one row per feature, with columns representing 
    feature identifiers, attributes, coordinates, and/or data values. Files may 
    be gzipped compressed.

- --out &lt;filename>

    Optionally specify the name of of the output file. The default is to use 
    the input file base name. The '.gff' extension is automatically
    added if required.

- --noheader

    The input file does not have column headers, often found with UCSC 
    derived annotation data tables. 

- --zero

    Input file is in interbase or 0-based coordinates. This should be 
    automatically detected for most known file formats, e.g. BED.

### Column indices

- --ask

    Indicate that the program should interactively ask for column indices or
    text strings for the GFF attributes, including coordinates, source, type, 
    etc. It will present a list of the column names to choose from. Enter 
    nothing for non-relevant columns or to accept default values.

- --chr &lt;column\_index>

    The index of the dataset in the data table to be used 
    as the chromosome or sequence ID column in the gff data.

- --start &lt;column\_index>
- --begin &lt;column\_index>

    The index of the dataset in the data table to be used 
    as the start position column in the gff data.

- --stop &lt;column\_index>
- --end &lt;column\_index>

    The index of the dataset in the data table to be used 
    as the stop or end position column in the gff data.

- --score &lt;column\_index>

    The index of the dataset in the data table to be used 
    as the score column in the gff data.

- --strand &lt;column\_index>

    The index of the dataset in the data table to be used
    for strand information. Accepted values might include
    any of the following "+, -, 1, -1, 0, .".

- --name &lt;text | column\_index>

    Enter either the text that will be shared name among 
    all the features, or the index of the dataset in the data 
    table to be used as the name of each gff feature. This 
    information will be used in the 'group' column.

- --id &lt;column\_index>

    The index of the dataset in the data table to be used
    as the unique ID of each gff feature. This information
    will be used in the 'group' column of GFF v3 files 
    only. The default is to automatically generate a 
    unique identifier.

- --tags &lt;column\_indices>

    Provide a comma delimited list of column indices that contain 
    values to be included as group tags in the GFF features. The 
    key will be the column name.

- --source &lt;text | column\_index>

    Enter either a text string or a column index representing the 
    GFF source that should be used for the features. The default is 
    'data'.

- --type &lt;text | column\_index>

    Enter either a text string or a column index representing the 
    GFF 'type' or 'method' that should be used for the features. If 
    not defined, it will use the column name for either 
    the 'score' or 'name' column, if defined. As a last resort, it 
    will use the most creative method of 'Experiment'.

### General options

- --unique

    Indicate whether the feature names should be made unique. A count 
    number is appended to the name of subsequent features to make them 
    unique. Using a base text string for the name will automatically 
    generate unique names.

- --sort

    Sort the output file by genomic coordinates. Automatically enabled 
    when compressing with bgzip. 

- --gz

    Indicate whether the output file should be compressed with gzip.

- --bgz

    Specify whether the output file should be compressed with block gzip 
    (bgzip) for tabix compatibility.

- --version

    Print the version number.

- --help

    Display the POD documentation

## DESCRIPTION

This program will convert a data file into a GFF version 3 formatted text file. 
Only simple conversions are performed, where each data line is converted 
to a single feature. Complex features with parent-child relationships (such 
as genes) should be converted with something more advanced.

## AUTHOR

    Timothy J. Parnell, PhD
    Howard Hughes Medical Institute
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
