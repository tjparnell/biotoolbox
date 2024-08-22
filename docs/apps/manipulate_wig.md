# Bio::ToolBox

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## manipulate\_wig.pl

A progam to manipulate wiggle files.

## SYNOPSIS

manipulate\_wig.pl \[options\] -i &lt;file1.wig> -o &lt;file1.out.wig>

    File Options: 
    -i --in <file>            Input file. Accepts 'stdin'.
    -o --out <file>           Output file. Accepts 'stdout'.
    
    Selection functions:
    -k --skip <regex>         Skip lines where chromosomes match regex 
    -y --apply <regex>        Only apply manipulations to matching chromosomes 
    
    Manipulation functions (in order of execution):
    -u --null                 Convert null, NA, N/A, NaN, inf values to 0
    -d --delog [2|10]         Delog values of given base 
    -b --abs                  Convert to the absolute value 
    -m --mult <float>         Multiply score by the given value
    -a --add <float>          Add the given value to the score
    -l --log [2|10]           Convert to log2 or log10. 
    -n --min <float>          Set the minimum score
    -x --max <float>          Set the maximum score
    -p --place <int>          Format score to decimal positions
    -z --zero                 Discard lines with zero values

    BigWig support:
    --chromo <file>           Chromosome sizes file for writing bigWig
    --db <file>               Indexed file to obtain chromosome info
    --bw2w <path>             Path to UCSC bigWigToWig utility
    --w2bw <path>             Path to UCSC wigToBigWig utility
    
    General functions:
    -t --stats                Calculate statistics 
    -v --version              print version and exit
    -h --help                 show extended documentation

## OPTIONS

The command line flags and descriptions:

### File options

- --in &lt;file>

    Specify the input wig file. All three formats, variableStep, fixedStep, and 
    bedGraph, are supported. Files may be gzipped. BigWig files are supported, 
    so long as the UCSC bigWigToWig utility is available. Alternatively, the input 
    may be read from standard input by specifying 'stdin' as the file name. 

- --out &lt;file>

    Specify the output wig file. The output format will be the same format as the
    input. The file may be gzipped by appending `.gz` to the name. BigWig files are
    supported, so long as the UCSC wigToBigWig utility is available and a chromosome
    file is provided. Alternatively, the output may be sent to standard output by
    specifying 'stdout' as the file name. 

### Selection functions

- --skip &lt;regex>

    Selectively skip (discard) lines corresponding to certain chromosomes that 
    match the provided regular expression. For example, skip the 
    mitochondrial and random contigs, use "chrM|chrUn|random".

- --apply &lt;regex>

    Selectively apply manipulation functions to certain chromosomes that match 
    provided regular expression, leaving remaining lines untouched. For example, 
    to apply a normalization to the X chromosome, use 'chrX'.

### Manipulation functions

- --null

    Convert lines with a score of `null`, `NA`, `N/A`, `NaN`, or `inf` to 
    a value of 0. 

- --delog \[2|10\]

    Convert lines from log space in the indicated base.

- --abs

    Convert line scores to absolute values.

- --mult &lt;float>

    Multiply line scores by the indicated value.

- --add &lt;float>

    Add the indicated value to each line score.

- --log \[2|10\]

    Convert the line score to a log equivalent in the indicated base space.

- --min &lt;float>

    Set the minimum floor score. Any score below the indicated value 
    will be set to the indicated value.

- --max &lt;float>

    Set the maximum ceiling score. Any score above the indicated value 
    will be set to the indicated value.

- --place &lt;integer>

    Format the score value to the indicated number of decimal positions.

- --zero 

    Discard lines with a score value of zero.

### BigWig support

- chromo &lt;file>

    When writing to a bigWig output file, provide a chromosome sizes text 
    file for use with the `wigToBigWig` utility. Alternatively, use a 
    database file, below.

- db &lt;file>

    When writing to a bigWig output file, provide an indexed database file, 
    such as another bigWig file, Bam, indexed Fasta, etc, for automatically 
    generating a chromosome sizes text file to use with the `wigToBigWig` 
    utility. If a bigWig input file was specified, it will be conveniently 
    substituted as a database. **Note** that the `--skip` option will be 
    applied to the generated chromosome file.

- bw2w &lt;path>

    If the UCSC `bigWigToWig` utility is not in your environment `PATH`, 
    provide the path with this option.

- w2bw &lt;path>

    If the UCSC `wigToBigWig` utility is not in your environment `PATH`, 
    provide the path with this option.

### General functions

- --stats

    Calculate the statistics across the genome at base pair resolution. 
    Statistics are calculated after any processing indicated above are 
    performed. Only applied chromosomes are calculated. Results are 
    printed to STDOUT, or STDERR if standard out is used as an output.

- --version

    Print the version number of the program and exit.

- --help

    Display the POD documentation using perldoc. 

## DESCRIPTION

A program to manipulate the score value of wig files. This will process all 
forms of text based wig files, including fixedStep, variableStep, and bedGraph. 
Files may be gzip compressed. BigWig files are also transparently supported as 
both input and output, provided that the appropriate UCSC utility files are 
available.

**NOTE:** More than one option may be specified! The options above are the order 
in which the score is manipulated. If they are not in the order you want, you 
may have to pipe to sequential instances. Use 'stdin' and 'stdout' for filenames.

## AUTHOR

    Timothy J. Parnell, PhD
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
