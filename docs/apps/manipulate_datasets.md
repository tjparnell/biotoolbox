# Bio::ToolBox

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## manipulate\_datasets.pl

A progam to manipulate tab-delimited data files.

## SYNOPSIS

manipulate\_datasets.pl \[--options ...\] &lt;filename> 

    File options:
    -i --in <filename>                input data file
    -o --out <filename>               output file, default overwrite
    -H --noheader                     input file has no header row
    
    Non-interactive functions:
    -f --func [ reorder | delete | rename | new | number | concatenate | 
                split | coordinate | sort | gsort | null | duplicate | 
                above | below | specific | keep | addname | cnull | 
                absolute | minimum | maximum | log | delog | format | pr | 
                add | subtract | multiply | divide | combine | scale | 
                zscore | ratio | diff | normdiff | center | rewrite | 
                export | treeview | summary | stat ]
    -x --index <integers>             column index to work on
    
    Operation options:
    -n --exp --num <integer>          numerator column index for ratio
    -d --con --den <integer>          denominator column index for ratio
    -t --target <text> or <number>    target value for certain functions
    --place [r | n]                   replace column contents or new column
    --(no)zero                        include zero in certain functions
    --dir [i | d]                     sort order: increase or decrease
    --name <text>                     name of new column
    --log                             values are in log scale
    
    General Options:
    -z --gz                           compress output file
    -Z --bgz                          bgzip compress output file
    -v --version                      print version and exit
    -h --help                         show extended documentation

## OPTIONS

The command line flags and descriptions:

### File options

- --in &lt;filename>

    Provide the name of an input file. The file must be a tab-delimited text file,
    with each row specifiying a genomic feature (gene) or region, and each column
    representing identifying information or a dataset value. The first row should
    include the name for each column, e.g. the name of the database dataset from
    which the values were collected.

    If the file was generated using a BioToolBox script, then each column may have
    metadata stored in a header comment line at the beginning of the file.
    Manipulations on the data in a column dataset will be added to this metadata and
    recorded upon file write.

    The input file may be compressed using the gzip program and recognized with 
    the extension '.gz'.  

- --out &lt;filename>

    Optionally provide an alternative output file name. If no name is provided, 
    then the input file will be overwritten with a new file. Appropriate 
    extensions will be appended as necessary.

- --noheader

    Indicate that the input file has no column header line, and that dummy 
    headers will be provided. Not necessary for BED, GFF, or recognized UCSC 
    file formats.

### Non-interactive functions

- --func &lt;function>

    The program is designed to be run interactively. However, single manipulations 
    may be performed on single datasets by specifying a function name and any 
    other required options. These functions include the following.

    **reorder** **delete** **rename** **new** **number** **concatenate**
    **split** **coordinate** **sort** **gsort** **null** **duplicate** **above**
    **below** **specific** **keep** **cnull** **absolute** **minimum**
    **maximum** **log** **delog** **format** **pr** **add** **subtract**
    **multiply** **divide** **combine** **scale** **zscore** **ratio** **diff**
    **normdiff** **center** **rewrite** **export** **treeview** **summary** **stat**

    Refer to the FUNCTIONS section for details.

- --index &lt;integers>

    Provide the index number of the column(s) on which to perform the 
    function(s). Multiple indices may also be specified using a comma delimited 
    list without spaces. Ranges of contiguous indices may be specified using a 
    dash between the start and stop. For example, '1,2,5-7,9' would indicate 
    datasets '1, 2, 5, 6, 7, and 9'. 

### Operation options

- --exp &lt;integer>
- --num &lt;integer>

    Specify the index number to be used for the experiment or numerator 
    column with the 'ratio' or 'difference' functions. Both flags are aliases
    for the same thing.

- --con &lt;integer>
- --den &lt;integer>

    Specify the index number to be used for the control or denominator 
    column with the 'ratio' or 'difference' functions. Both flags are aliases
    for the same thing.

- --target &lt;string> or &lt;number>

    Specify the target value when using various functions. This is a catch-all 
    option for a number of functions. Please refer to the function description 
    for more information.

- --place \[r | n\]

    Specify the placement of a transformed column. Two values are accepted ('r' 
    or 'n'):

        - (r)eplace the original column with new values
        - add as a (n)ew column

    Defaults to new placement when executed automatically using the --func 
    option, or prompts the user when executed interactively.

- --(no)zero

    Specify that zero values should or should not be included when 
    calculating statistics or converting null values on a column. The default is 
    undefined, meaning the program may ask the user interactively whether to 
    include zero values or not.

- --dir \[i | d\]

    Specify the direction of a sort: 

        - (i)ncreasing
        - (d)ecreasing
        

- --name &lt;string>

    Specify a new column name when re-naming a column using the rename function 
    from the command line. Also, when generating a new column using a defined 
    function (--func &lt;function>) from the command line, the new column will use 
    this name.

- --log 

    Indicate whether the data is (not) in log2 space. This is required to ensure 
    accurate mathematical calculations in some manipulations. This is not necessary 
    when the log status is appropriately recorded in the dataset metadata.

### General options

- --gz 

    Indicate whether the output file should be gzip compressed. The compression 
    status of the input file will be preserved if overwriting.

- --bgz

    Specify whether the output file should be compressed with block gzip 
    (bgzip) for tabix compatibility.

- --version

    Print the version number.

- --help

    Display the POD documentation using perldoc. 

## DESCRIPTION

This program allows some common mathematical and other manipulations on one
or more columns in a datafile. The program is designed as a simple
replacement for common manipulations performed in a full featured
spreadsheet program, e.g. Excel, particularly with datasets too large to be
loaded, all in a conveniant command line program. The program is designed
to be operated primarily as an interactive program, allowing for multiple
manipulations to be performed. Alternatively, single manipulations may be
performed as specified using command line options. As such, the program can
be called in shell scripts.

Note that the datafile is loaded entirely in memory. For extremely large 
datafiles, e.g. binned genomic data, it may be best to first split the 
file into chunks (use `split_data_file.pl`), perform the manipulations, 
and recombine the file (use `join_data_file.pl`). This could be done 
through a simple shell script.

The program keeps track of the number of manipulations performed, and if 
any are performed, will write out to file the changed data. Unless an 
output file name is provided, it will overwrite the input file (NO backup is
made!).

## FUNCTIONS

This is a list of the functions available for manipulating columns. These may 
be selected interactively from the main menu (note the case sensitivity!), 
or specified on the command line using the --func option.

- **stat** (menu option **t**)

    Print some basic statistics for a column, including mean, 
    median, standard deviation, min, and max. If 0 values are present,
    indicate whether to include them (y or n)

- **reorder** (menu option **R**)

    The column may be rewritten in a different order. The new order 
    is requested as a string of index numbers in the desired order. 
    Also, a column may be deleted by skipping its number or duplicated
    by including it twice.

- **delete** (menu option **D**)

    One or more column may be selected for deletion.

- **rename** (menu option **n**)

    Assign a new name to a column. For automatic execution, use the --name 
    option to specify the new name. Also, for any automatically executed 
    function (using the --func option) that generates a new column, the 
    column's new name may be explicitly defined with this option.

- **number** (menu option **b**)

    Assign a number to each datapoint (or line), incrementing from 1 
    to the end. The numbered column will be inserted after the requested 
    column index.

- **concatenate** (menu option **C**)

    Concatenate the values from two or more columns into a single new 
    column. The character used to join the values may be specified 
    interactively or by the command line option --target (default is '\_' 
    in automatic execution mode). The new column is appended at the end.

- **split** (menu option **T**)

    Split a column into two or more new columns using a specified character 
    as the delimiter. The character may be specified interactively or 
    with the --target command line option (default is '\_' in automatic 
    execution mode). The new columns are appended at the end. If the 
    number of split items are not equal amongst the rows, absent values 
    are appended with null values.

- **coordinate** (menu option **O**)

    Generate a coordinate string from the chromosome, start, and, if 
    present, stop coordinate values as a new column. The string will 
    have the format "chr:start-stop" or "chr:start". This is useful 
    in making unique identifiers or working with genome browsers.

- **sort** (menu option **o**)

    The entire data table is sorted by a specific column, or by the 
    mean of a list of columns if more than one is provided. The first
    datapoint is checked for the presence of letters, and the data 
    table is then sorted either asciibetically or numerically. The 
    direction of sort, (i)ncreasing or (d)ecreasing, is requested. 

- **gsort** (menu option **g**)

    The entire data table is sorted by increasing genomic position, 
    first by chromosome then by start position. These columns must exist 
    and have recognizable names (e.g. 'chromo', 'chromosome', 'start').

- **null** (menu option **N**)

    Delete rows that contain a null value in one or more 
    columns. Some of the other functions may not work properly if
    a non-value is present. If 0 values are present, indicate whether
    to toss them (y or n). This may also be specified as a command line 
    option using the --except flag.

- **duplicate** (menu option **P**)

    Delete rows with duplicate values. One or more columns may be 
    selected to search for duplicate values. Values are simply concatenated 
    when multiple columns are selected. Rows with duplicated values are 
    deleted, always leaving the first row.

- **above** (menu option **A**)

    Delete rows with values that are above a certain threshold value. 
    One or more columns may be selected to test values for the 
    threshold. The threshold value may be requested interactively or 
    specified with the --target option.

- **below** (menu option **B**)

    Delete rows with values that are below a certain threshold value. 
    One or more columns may be selected to test values for the 
    threshold. The threshold value may be requested interactively or 
    specified with the --target option.

- **specific** (menu option **S**)

    Delete rows with values that contain a specific value, either text 
    or number. One or more columns may be selected to check for values. 
    The specific values may be selected interactively from a list or 
    specified with the --target option.

- **keep** (menu option **K**)

    Keep only those rows with values that contain a specific value, 
    either text or number. One or more columns may be selected to check 
    for values. The specific values may be selected interactively from a 
    list or specified with the --target option.

- **addname** (menu option **M**)

    Add or update the name of each feature or row. If the data table 
    already has a Name column, the value will be updated. Otherwise a 
    new column will be added. The name will be a text prefix followed 
    by an integer (row index). The prefix may be defined by setting the 
    \--target option, interactively provided by the user, or taken from 
    the general table feature metadata.

- **cnull** (menu option **U**)

    Convert null values to a specific value. One or more columns may 
    be selected to convert null values. The new value may be requested 
    interactively or defined with the --target option.  

- **absolute** (menu option **G**)

    Convert signed values to their absolute value equivalents. One or 
    more columns may be selected to convert.

- **minimum** (menu option **I**)

    Reset datapoints whose values are less than a specified minimum 
    value to the minimum value. One or more columns may be selected 
    to reset values to the minimum. The minimum value may be requested 
    interactively or specified with the --target option. 

- **maximum** (menu option **X**)

    Reset datapoints whose values are greater than a specified maximum 
    value to the maximum value. One or more columns may be selected 
    to reset values to the maximum. The maximum value may be requested 
    interactively or specified with the --target option. 

- **add** (menu option **a**)

    Add a value to a column. A real number may be supplied, or the words
    'mean', 'median', or 'sum' may be entered as a proxy for those statistical
    values of the column. The column may either be replaced or added
    as a new one. For automatic execution, specify the number using the
    \--target option.

- **subtract** (menu option **u**)

    Subtract a value from a column. A real number may be supplied, or the words
    'mean', 'median', or 'sum' may be entered as a proxy for those statistical
    values of the column. The column may either be replaced or added
    as a new one. For automatic execution, specify the number using the
    \--target option.

- **multiply** (menu option **y**)

    Multiply a column by a value. A real number may be supplied, or the words
    'mean', 'median', or 'sum' may be entered as a proxy for those statistical
    values of the column. The column may either be replaced or added
    as a new one. For automatic execution, specify the number using the
    \--target option.

- **divide** (menu option **v**)

    Divide a column by a value. A real number may be supplied, or the words
    'mean', 'median', or 'sum' may be entered as a proxy for those statistical
    values of the column. The column may either be replaced or added
    as a new one. For automatic execution, specify the number using the
    \--target option.

- **scale** (menu option **s**)

    A column may be a median scaled as a means of normalization 
    with other columns. The current median of the column requested is
    presented, and a new median target is requested. The column may 
    either be replaced with the median scaled values or added as a new 
    column. For automatic execution, specify the new median target 
    with the --target option.

- **pr** (menu option **p**)

    A column may be converted to a percentile rank, whereby the
    data values are sorted in ascending order and assigned a new value 
    from 0..1 based on its rank in the sorted order. The column may 
    either be replaced with the percentile rank or added as a new
    column. The original order of the column is maintained.

- **zscore** (menu option **Z**)

    Generate a Z-score or standard score for each value in a column. The
    Z-score is the number of standard deviations the value is away from
    the column's mean, such that the new mean is 0 and the standard 
    deviation is 1. Provides a simple method of normalizing columns
    with disparate dynamic ranges.

- **log** (menu option **l**)

    A column may be converted to log values. The column may either 
    be replaced with the log values or added as a new column. Use 
    the --target option to specify the base (usually 2 or 10).

- **delog** (menu option **L**)

    A column that is currently in log space may be converted back to
    normal numbers. The column may either be replaced with the 
    new values or added as a new column. Use the --target option to 
    specify the base (usually 2 or 10). The base may be obtained from the 
    metadata.

- **format** (menu option **f**)

    Format the numbers of a column to a given number of decimal places. 
    An integer must be provided. The column may either be replaced or 
    added as a new column. For automatic execution, use the --target 
    option to specify the number decimal places.

- **combine** (menu option **c**)

    Mathematically combine the data values in two or more columns. The 
    methods for combining the values include mean, median, min, max, 
    stdev, or sum. The method may be specified on the command line 
    using the --target option. The combined data values are added as a 
    new column.

- **ratio** (menu option **r**)

    A ratio may be generated between two columns. The experiment and 
    control columns are requested and the ratio is added as a new
    column. For log2 numbers, the control is subtracted from the
    experiment. The log2 status is checked in the metadata for the 
    specified columns, or may be specified as a command line option, or
    asked of the user.

- **diff** (menu option **d**)

    A simple difference is generated between two existing columns. The 
    values in the 'control' column are simply subtracted from the 
    values in the 'experimental' column and recorded as a new column.
    For enumerated columns (e.g. tag counts from Next Generation 
    Sequencing), the columns should be subsampled to equalize the sums 
    of the two columns. The indices for the experimental and control columns 
    may either requested from the user or supplied by the --exp and 
    \--con command line options. 

- **normdiff** (menu option **z**)

    A normalized difference is generated between two existing columns. 
    The difference between 'control' and 'experimental' column values 
    is divided by the square root of the sum (an approximation of the 
    standard deviation). This is supposed to yield fewer false positives
    than a simple difference (see Nix et al, BMC Bioinformatics, 2008).
    For enumerated datasets (e.g. tag counts from Next Generation 
    Sequencing), the datasets should be subsampled to equalize the sums 
    of the two datasets. The indices for the experimental and control columns 
    may either requested from the user or supplied by the --exp and 
    \--con command line options. 

- **center** (menu option **e**)

    Center normalize the datapoints in a row by subtracting the mean or
    median of the datapoints. The range of columns is requested or 
    provided by the --index option. Old values are replaced by new 
    values. This is useful for visualizing data as a heat map, for example.

- **new** (menu option **w**)

    Generate a new column which contains an identical value for 
    each datapoint (row). The value may be either requested interactively or 
    supplied using the --target option. This function may be useful for 
    assigning a common value to all of the data points before joining the 
    data file with another.

- **summary** (menu option **y**)

    Write out a summary of collected windowed data file, in which the mean 
    for each of the data columns is calculated, transposed (columns become 
    rows), and written to a new data file. This is essentially identical to 
    the summary function from the biotoolbox analysis scripts 
    [map\_relative\_data.pl](https://metacpan.org/pod/map_relative_data.pl) and [pull\_features.pl](https://metacpan.org/pod/pull_features.pl). It assumes that each 
    column has start and stop metadata. The program will automatically 
    identify available columns to summarize based on their name. In 
    interactive mode, it will request the contiguous range of start and 
    ending columns to summarize. The contiguous columns may also be 
    indicated using the --index option. The method of summarizing the 
    data can be specified interactively or with the --target option. 
    Methods include 'mean' (default), 'median', or 'trimmean', where
    the top and bottom 1% of values are discarded and a mean determined
    from the remaining 98% of values. By default, a new file using the 
    input file base name appended with '\_&lt;method>\_summary' is written, or
    a filename may be specified using the --out option.

- **export** (menu option **x**)

    Export the data into a simple tab-delimited text file that contains no 
    metadata or header information. Non-values '.' are converted to  
    true nulls. If an output file name is specified using the --outfile 
    option, it will be used. Otherwise, a possible filename will be 
    suggested based on the input file name. If any modifications are 
    made to the data structure, a normal data file will still be written. 
    Note that this could overwrite the exported file if the output file name
    was specified on the command line, as both file write subroutines will 
    use the same name!

- **treeview** (menu option **i**)

    Export the data to the CDT format compatible with both Treeview and 
    Cluster programs for visualizing and/or generating clusters. Specify the 
    columns containing a unique name and the columns to be analyzed (e.g. 
    \--index &lt;name>,&lt;start-stop>). Extraneous columns are removed. 
    Additional manipulations on the columns may be performed prior to 
    exporting. These may be chosen interactively or using the codes 
    listed below and specified using the --target option.

        su - decreasing sort by sum of row values
        sm - decreasing sort by mean of row values
        cg - median center features (rows)
        cd - median center datasets (columns)
        zd - convert columns to Z-scores
        pd - convert columns to percentile ranks
        L2 - convert values to log2
        L10 - convert values to log10
        n0 - convert nulls to 0.0

    A simple Cluster data text file is written (default file name 
    "&lt;basename>.cdt"), but without the GWEIGHT column or EWEIGHT row. The 
    original file will not be rewritten.

- **rewrite** (menu option **W**)

    Force the data file contents to be re-written. Useful if you want to 
    write an intermediate file during numerous interactive manipulations. 
    Consider this as a 'Save as...'.

## AUTHOR

    Timothy J. Parnell, PhD
    Howard Hughes Medical Institute
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
