# Bio::ToolBox - db_types

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## db\_types.pl

A program to print out the available feature types in a database.

## SYNOPSIS

db\_types.pl &lt;database>

    Options:
    --db <database>
    --version
    --help
    

## OPTIONS

The command line flags and descriptions:

- --db &lt;database>

    Specify the name of a [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) annotation database 
    from which gene or feature annotation may be derived. If not specified, 
    then a list of known databases in the BioToolBox configuration file 
    `.biotoolbox.cfg` will be presented as a list to the user.

- --version

    Print the version number.

- --help

    Display this POD documentation.

## DESCRIPTION

This program will print a list of all of the known feature types present 
in a [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) database. The types are organized into 
groups by their source tag.

BigWigSet databases, comprised of a directory of BigWig files and a 
metadata file, are also supported.

## AUTHOR

    Timothy J. Parnell, PhD
    Howard Hughes Medical Institute
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
