# Bio::ToolBox - split\_data\_file

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## split\_data\_file.pl

A program to split a data file by rows based on common data values.

## SYNOPSIS

split\_data\_file.pl \[--options\] &lt;filename>

    File options:
    -i --in <filename>                (txt bed gff gtf vcf refFlat ucsc etc)
    -p --prefix <text>                output file prefix (input basename)
    -H --noheader                     input file has no headers
    
    Splitting options:
    -x --index <column_index>         column with values to split upon
    -t --tag <text>                   use VCF/GFF attribute
    -m --max <integer>                maximum number of items per output file
    
    General options:
    -z --gz                           compress output file
    -v --version                      print version and exit
    -h --help                         show extended documentation

## OPTIONS

The command line flags and descriptions:

### File options

- --in &lt;filename>

    Specify the file name of a data file. It must be a tab-delimited text file. 
    The file may be compressed with gzip.

- --prefix &lt;text>

    Optionally provide a filename prefix for the output files. The default 
    prefix is the input filename base name. If no prefix is desired, using 
    just the values as filenames, then set the prefix to 'none'.

- --noheader

    Indicate that the input file has no column header line, and that dummy 
    headers will be provided. Not necessary for BED, GFF, or recognized UCSC 
    file formats.

### Splitting options

- --index &lt;column\_index>

    Provide the index number of the column or dataset containing the values 
    used to split the file. If not specified, then the index is requested 
    from the user in an interactive mode.

- --tag &lt;text>

    Provide the attribute tag name that contains the values to split the 
    file. Attributes are supported by GFF and VCF files. If splitting a 
    VCF file, please also provide the column index. The INFO column is 
    index 8, and sample columns begin at index 10.

- --max &lt;integer>

    Optionally specify the maximum number of data lines to write to each 
    file. Each group of specific value data is written to one or more files. 
    Enter as an integer; underscores may be used as thousands separator, e.g. 
    100\_000. 

### General options

- --gz

    Indicate whether the output files should be compressed 
    with gzip. Default behavior is to preserve the compression 
    status of the input file.

- --version

    Print the version number.

- --help

    Display the POD documentation

## DESCRIPTION

This program will split a data file into multiple files based on common 
values in the data table. All rows with the same value will be 
written into the same file. A good example is chromosome, where all 
data points for a given chromosome will be written to a separate file, 
resulting in multiple files representing each chromosome found in the 
original file. The column containing the values to split and group 
should be indicated; if the column is not sepcified, it may be 
selected interactively from a list of column headers. 

This program can also split files based on an attribute tag in GFF or 
VCF files. Attributes are often specially formatted delimited key value 
pairs associated with each feature in the file. Provide the name of the 
attribute tag to split the file. Since attributes may vary based on 
the feature type, an interactive list is not supplied from which to 
choose the attribute.

If the max argument is set, then each group will be written to one or 
more files, with each file having no more than the indicated maximum 
number of data lines. This is useful to keep the file size reasonable, 
especially when processing the files further and free memory is 
constrained. A reasonable limit may be 100K or 1M lines.

The resulting files will be named using the basename of the input file, 
appended with the unique group value (for example, the chromosome name)
demarcated with a #. If a maximum line limit is set, then the file part 
number is appended to the basename, padded with zeros to three digits 
(to assist in sorting). Each file will have duplicated and preserved 
metadata. The original file is preserved.

This program is intended as the complement to 'join\_data\_files.pl'.

## AUTHOR

    Timothy J. Parnell, PhD
    Howard Hughes Medical Institute
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
