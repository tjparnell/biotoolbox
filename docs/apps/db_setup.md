# Bio::ToolBox - db_setup

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## db\_setup.pl

A program to setup a SeqFeature::Store SQLite database from UCSC data

## SYNOPSIS

db\_setup.pl \[--options...\] <UCSC database>

    Options:
    -d --db <UCSC database>
    -p --path </path/to/store/database> 
    -t --table [refGene|knownGene|all]
    -k --keep
    -V --verbose
    -v --version
    -h --help

## OPTIONS

The command line flags and descriptions:

- --db <UCSC database>

    Provide the short UCSC database name for the species and version you want 
    to use. See [http://genome.ucsc.edu/FAQ/FAQreleases.html](http://genome.ucsc.edu/FAQ/FAQreleases.html) for a current 
    list of available UCSC genomes. Examples include hg19, mm10, danRer7, 
    sacCer3, etc.

- --path &lt;/path/to/store/database>

    Specify the optional path to store the SQLite database file. The default 
    path is the `~/lib`.

- --table \[refGene|knownGene|all\]

    Provide one or more UCSC tables to load into the database. They may be 
    specified as comma-delimited list (no spaces) or as separate, repeated 
    arguments. The default is refGene and knownGene (if available).

- --keep

    Keep the downloaded UCSC tables and converted GFF3 files. Default is to 
    delete them.

- --verbose

    Show realtime database loading progress.

- --version

    Print the version number.

- --help

    Display this POD documentation.

## DESCRIPTION

This program will simplify the task of generating an annotation database. You 
provide the short name of the UCSC database for the species and genome version 
you are interested in, and the script will automatically download gene annotation 
and build a _Bio::DB::SeqFeature::Store_ database for use with BioToolBox 
scripts. 

## AUTHOR

    Timothy J. Parnell, PhD
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
