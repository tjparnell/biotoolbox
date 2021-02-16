# NAME

Bio::ToolBox - Tools for querying and analysis of genomic data

# DESCRIPTION

This package provides a number of Perl modules and scripts for working 
with common bioinformatic data. Many bioinformatic data analysis revolves 
around working with tables of information, including lists of 
genomic annotation (genes, promoters, etc.) or defined regions 
of interest (epigenetic enrichment, transcription factor binding 
sites, etc.). This library works with these tables and provides 
a set of common tools for working with them.

- Opening and saving common tab-delimited text formats
- Support for BED, GFF, VCF, narrowPeak files
- Scoring intervals and annotation with datasets from microarray or sequencing
experiments, including ChIPSeq, RNASeq, and more
- Support for Bam, BigWig, BigBed, wig, and USeq data formats
- Works with any genomic annotation in GTF, GFF3, and various UCSC formats, including
refFlat, knownGene, genePred and genePredExt formats

The libraries provide a unified and integrated approach to analyses. 
In many cases, they provide an abstraction layer over a variety of 
different specialized [Bio::Perl](https://metacpan.org/pod/Bio::Perl) 
and related modules. Instead of writing numerous scripts specialized for 
each data format (wig, bigWig, Bam), one script can now work with virtually 
any data format. 

# INSTALLATION

[Released versions](https://metacpan.org/pod/Bio::ToolBox) can be installed 
from [CPAN](https://metacpan.org) using your favorite installer. For example,
using [CPAN Minus](https://metacpan.org/pod/App::cpanminus) 

    cpanm Bio::ToolBox

Manual installation is simple with the standard [Module::Build](https://metacpan.org/pod/Module::Build) 
incantation. 

    perl ./Build.PL
    ./Build
    ./Build test
    ./Build install

In either case, this will get you a minimal installation that will work with 
text files (BED, GFF, GTF, etc), but not binary files. To work with binary Bam and 
BigWig files, two additional [external libraries](docs/AdvancedInstallation.md#external-libraries)
must also be compiled and installed; This is not hard, and you likely already have 
one (maybe both) installed on your system. Most scripts should fail gently with 
warnings if required modules are missing.

For step-by-step instructions to get a complete installation, see the 
[Advanced Installation guide](docs/AdvancedInstallation.md).

## Docker

For those so inclined, a [dockerfile](Docker/ReadMe.md) is provided for your convenience 
in building a Docker image following the advanced installation guide.

# LIBRARIES

Several user-oriented library modules are included in this distribution for 
working with bioinformatic data. They provide a foundation for the included 
analysis scripts, and can be used for custom coding projects. See 
[Libraries](docs/Libraries.md) for more information.

# SCRIPTS

The BioToolBox package comes complete with a suite of high-quality production-ready 
scripts ready for a variety of analyses. A sampling of what can be done include 
the following:

- Annotated feature collection and selection
- Data collection and scoring for features
- Data file format manipulation and conversion
- Low-level processing of sequencing data into customizable wig representation

Scripts have built-in documentation. Execute the script without any options to print 
a synopsis of available options, or add `--help` to print the full documentation.

See [Script examples](docs/Scripts.md) for more information.

# AUTHOR

	Timothy J. Parnell, PhD
	Huntsman Cancer Institute
	University of Utah
	Salt Lake City, UT, 84112

# LICENSE

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. For details, see the
full text of the license in the file LICENSE.

This package is distributed in the hope that it will be useful, but it
is provided "as is" and without any express or implied warranties. For
details, see the full text of the license in the file LICENSE.
