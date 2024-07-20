# Bio::ToolBox

## Tools for querying and analysis of genomic data

This package provides a number of Perl modules and applications for working 
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

## Applications

This package includes a number of production-quality, well-documented, command line 
applications for the advanced and novice bioinformatic analyst. These are 
described on the [Applications](Applications.md) page, with links therein to 
help documentation. 

For some common scenario and real-world examples, see the [Examples](Examples.md)
page. Most of these applications are used by the author on a weekly basis in production.

For an example of using this package as a basis for an epigenetic pipeline, see the
[Multi-Replica Macs ChIPSeq](https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq)
package.


## Code Libraries

The Bio::ToolBox package includes a number of Perl modules with extensive 
documentation that can be used in custom scripts. See the
[Included Libraries](Libraries.md) page for a brief explanation of what is available.


## Installation

[Released versions](https://metacpan.org/pod/Bio::ToolBox) can be installed 
from [CPAN](https://metacpan.org) using your favorite Perl package manger. 
For example, using [CPAN Minus](https://metacpan.org/pod/App::cpanminus) 

    cpanm Bio::ToolBox

Manual installation is simple with the standard 
[Module::Build](https://metacpan.org/pod/Module::Build) incantation. 

    perl ./Build.PL
    ./Build
    ./Build test
    ./Build install

In either case, this will get you a minimal installation that will work with
text files (BED, GFF, GTF, etc), but not binary files. To work with binary Bam
and BigWig files, two additional
[external libraries](AdvancedInstallation.md#external-libraries)
must also be compiled and installed; This is not hard, and you likely already
have one (possibly both) installed on your system. Most scripts should fail
gently with warnings if required modules are missing. These modules are
indicated as `recommended` rather than `required` in `Build.PL` script.

For a step-by-step instructions to get a complete installation, see the 
[Advanced Installation guide](AdvancedInstallation.md).

For MacOS-specific installation issues, see the [MacOS notes](MacOSNotes.md) page.
Even though MacOS is unix, it has its own idiosyncrasies that I've discovered 
myself, and are written here in hope that they may be helpful to someone somewhere.


## Frequently Asked Questions

OK, maybe not frequently, but the [FAQ](FAQ.md) page is basically a list of 
questions, tips, and hints that only the author knows but that may be 
useful for others to know.


## License

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0. For details, see the
full text of the license in the
[LICENSE](https://github.com/tjparnell/biotoolbox/blob/master/LICENSE).

This package is distributed in the hope that it will be useful, but it
is provided “as is” and without any express or implied warranties. For
details, see the full text of the license in the file LICENSE.



	
