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

# DOCUMENTATION

See the [documentation page](https://tjparnell.github.io/biotoolbox) for full details.

# INSTALLATION

[Released versions](https://metacpan.org/pod/Bio::ToolBox) can be installed 
from CPAN using your favorite Perl installer.

Manual installation is simple with the standard `Module::Build` incantation. 

    perl ./Build.PL
    ./Build
    ./Build test
    ./Build install

In either case, this will get you a minimal installation that will work with 
text files (BED, GFF, GTF, etc), but not binary files. To work with binary files,
including Bam and BigWig files, see the details in the
[Installation Guide](https://tjparnell.github.io/biotoolbox/AdvancedInstallation.html).

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
