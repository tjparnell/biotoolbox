# NAME

pull\_features.pl

A program to pull out a specific list of data rows from a data file.

# SYNOPSIS

pull\_features.pl --data &lt;filename> --list &lt;filename> --out &lt;filename>

    File options:
    -d --data <filename>          Source of all data rows or features
    -l --list <filename>          List of specific row or feature names
    -o --out <filename>           Output file, or basename for group files
    
    Column index options:
    -x --dindex <index>           Data column index of row name to lookup
    -X --lindex <index>           List column index of row name to lookup
    -g --gindex <index>           Group column index for lookup
    
    Output options:
    -r --order [list | data]      Order of items in output based on
    -U --sum                      Generate a summary file
    --sumonly                     Skip output, just make a summary file
    --start <integer>             First data column to make a summary file
    --stopi <integer>             Last data column to make a summary file
    --log                         Summarized data is in log2 space
    
    General options:
    -v --version                  print version and exit
    -h --help                     show full documentation

# OPTIONS

The command line flags and descriptions:

## File options

- --data

    Specify a tab-delimited text file as the data source file. One of 
    the columns in the input file should contain the identifiers to be 
    used in the lookup. The file may be gzipped.

- --list

    Specify the name of a text file containing the list of feature 
    names or identifiers to pull. The file may be a single column or 
    tab-delimited multi-column file with column headers. A .kgg file 
    from a Cluster k-means analysis may be used.

- --out

    Specify the output file name. 

## Column index options

- --dindex &lt;integer>
- --lindex &lt;integer>

    Specify the index numbers of the columns in the data and list 
    files, respectively, containing the identifiers to match features. 
    If not specified, then the program will attempt to identify  
    appropriate matching columns with the same header name. If none 
    are specified, the user must select interactively from a list of 
    available column names. 

- --gindex &lt;integer>

    Specify the group column from the list file. This allows the data 
    file to be split into multiple output group files. A column named 
    'group' will automatically be identified. A .kgg file will 
    automatically use the Cluster column as the group index.

## Output options

- --order \[list | data\]

    Optionally specify the order of features in the output file. Two 
    options are available. Specify 'list' to match the order of features 
    in the list file. Or specify 'data' to match the order of features 
    in the data file. The default is list.

- --sum

    Indicate that the pulled data should be averaged across all 
    features at each position, suitable for graphing. A separate text 
    file with '\_summed' appended to the filename will be written.

- --sumonly

    Indicate that only a summary file should be written, and that the 
    pulled data file should be skipped. Useful if you're just after 
    the summary for graphing purposes.

- --starti &lt;integer>

    When re-summarizing the pulled data, indicate the start column 
    index that begins the range of datasets to summarize. Defaults 
    to the leftmost column without a standard feature description
    name.

- --stopi &lt;integer>

    When re-summarizing the pulled data, indicate the stop column
    index the ends the range of datasets to summarize. Defaults
    to the last or rightmost column.

- --log

    The data is in log2 space. Only necessary when re-summarizing the
    pulled data.

## General options

- --version

    Print the version number.

- --help

    Display this POD documentation.

# DESCRIPTION

Given a list of requested unique feature identifiers, this program will 
pull out those features (rows) from a datafile and write a new file. This 
program compares in function to a popular spreadsheet VLOOKUP command. 
The list is provided as a separate text file, either as a single column 
file or a multi-column tab-delimited from which one column is selected. 
All rows from the source data file that match an identifier in the list 
will be written to the new file. The order of the features in the output 
file may match either the list file or the data file. 

If the list file has a second group column, then the rows for each group 
will be written to separate files, with the output file name appended with 
the group identifier. Use the gindex option to specify the group column.

The program will also accept a Cluster gene file (with .kgg extension) 
as a list file with group information, where the clusters are the groups. 

The program will optionally regenerate a summed data file, in which values 
in the specified data columns are averaged and written out as rows in a 
separate data file. Compare this function to the summary option in the 
biotoolbox scripts [get\_relative\_data.pl](https://metacpan.org/pod/get_relative_data.pl) or [average\_gene.pl](https://metacpan.org/pod/average_gene.pl).

# AUTHOR

    Timothy J. Parnell, PhD
    Howard Hughes Medical Institute
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
