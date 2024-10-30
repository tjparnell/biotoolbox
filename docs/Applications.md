# Bio::ToolBox - Applications

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## Included Applications

There are 22 production-quality script applications included in the package. These
are listed below, grouped by general function. Also see the 
[application examples](Examples.md) page.

Applications have built-in documentation. Execute the script without any options
to print a synopsis of available options, or add `-h` or `--help` to print the full
documentation.


## Data conversion

These are applications for converting from one bioinformatic file type to another.

### bam2wig.pl

[bam2wig](apps/bam2wig.md) will generate coverage or point data representations
of alignments from a bam file in virtually every single possible way imaginable. 
It can generate fixedStep, variableStep, and bedGraph formats in either text or
bigWig file formats. Cross-strand correlation shift models may also be determined
for single-end alignments. Common scenarios and settings are described in 
[recommended settings](apps/bam2wig.md#RECOMMENDED_SETTINGS).

### data2bed.pl

[data2bed](apps/data2bed.md) is a general purpose application for converting any 
text file into a standard bed format, assuming that chromosome coordinate
information is available. 

### data2wig.pl

[data2wig](apps/data2wig.md) is a general purpose application for converting any 
text file into a wiggle file format. One score column may be chosen, or
multiple score columns may be provided and combined using a simple mathematical
method. Support for writing directly to bigWig files is included.

### data2fasta.pl

[data2fasta](apps/data2fasta.md) is a general purpose application for converting any 
text file into a fasta file. Sequence may either be in the source file, or extracted
from a provided genomic fasta file.

### data2gff.pl

[data2gff](apps/data2gff.md) is a general purpose application for converting any 
text file into a standard GFF3 format. Columns may be specified for attribute
tags. 

### ucsc_table2gff3.pl

[ucsc_table2gff3](apps/ucsc_table2gff3.md) will convert a limited set of file
table formats for genes from UCSC into more common GFF3 format.



## Feature annotation

### get_features.pl

[get_features](apps/get_features.md) will collect features from an annotation
file or local SeqFeature database. Usually a subset of features based on 
filtering criteria or provided list will be retrieved. For example, collecting 
all gene features matching a specific biotype or attribute quality from a genome
annotation GFF3 or GTF file. Features may be written out in the original format 
or converted to other formats. Subfeatures, such as CDS or exon features, are
included as appropriate.

### get_gene_regions.pl

[get_gene_regions](apps/get_gene_regions.md) will collect specific gene regions
from an annotation file that are not typically annotated but rather inferred.
These include, for example, introns, alternate or common exons or introns, first
or last exons or introns, untranslated regions (UTRs), etc. 

### get_feature_info.pl

[get_feature_info](apps/get_feature_info) will collect additional information 
from an annotation file for list of features. Often these may be key=value 
attributes embedded in a GTF or GFF3 annotation file and difficult to extract 
otherwise.


## Data collection

These are applications for collecting genomic data from datasets, typically 
scoring annotation features with genomic data provided in bigWig or bam file
formats. Due to technical and performance reasons, bigWig files are preferred
over bam files, although the latter have limited supported; see the
[FAQ](FAQ.md) for details.

### get_datasets.pl

[get_datasets](apps/get_datasets.md) is a general purpose data collection
application for collecting a single numerical value for a list of genomic
features. For example, collecting the mean alignment coverage over genes. A
variety of input annotation and data formats are supported, as well as methods
of collection. A variety of examples and
[usage settings](apps/get_datasets.md#EXAMPLES) are provided. 

### get_binned_data.pl

[get_binned_data](apps/get_binned_data.md) will divide the provided genomic
features into bins and collect a single numerical value for each bin. In this
manner, a profile of scores across the length of the genomic feature may be
generated. An option is available to generate a summary file of the column (bin)
means.

### get_relative_data.pl

[get_relative_data](apps/get_relative_data.md) will collect genomic data in 
bins around a fixed reference point of the provided genomic features. For 
example, coverage around the Transcription Start Site (the 5' end of a gene)
may be collected. An option is available to generate a summary file of the column
(bin) means.

### correlate_position_data.pl

[correlate_position_data](apps/correlate_position_data.md) will calculate a 
correlation between two datasets, e.g. bigWig files of scores, along the 
length of a genomic feature. This can help to determine whether a distribution 
of scores along the genomic feature has shifted between the scores. For example,
in chromatin biology, whether nucleosomal occupancy has shifted across mapped
nucleosomal positions upon transcriptional activation.


## Data manipulation

These applications allow one to manipulate collected datasets in tab-delimited
text files on the command line, without having to invoke complex `awk` or `sed`
commands or round-tripping through a GUI spreadsheet application. Keep in mind
that files generated by this package often have metadata or comment lines
prefixed by a `#` character, which can stymy or otherwise become lost using
other tools.

### manipulate_datasets.pl

[manipulate_datasets](apps/manipulate_datasets.md) is an interactive, menu-driven
application for performing all sorts of column, row, and value manipulations. 
Single functions may be specified on the command line.

### manipulate_wig.pl

[manipulate_wig](apps/manipulate_wig.md) will perform various numeric
transformations on the scores in a text wiggle or bigWig file. 


## File manipulation

These applications work on columns or rows of one or more tab-delimited text
files. Metadata or comment lines (prefixed by a `#` character) are correctly
handled.

### merge_datasets.pl

[merge_datasets](apps/merge_datasets.md) will copy the columns from two or more
tab-delimited text files into a single file. If row numbers are unequal or are
otherwise in a different order, a specified column may be used as a lookup
identifier to match rows. 

### split_data_file.pl

[split_data_file](apps/split_data_file.md) will split a tab-delimited text file
by rows, either by a specific value in a specified column or by maximum allowable
rows per file. 

### join_data_file.pl

[join_data_file](apps/join_data_file.md) will concatenate the rows of multiple
tab-delimited text files into a single file.

### pull_features.pl

[pull_features](apps/pull_features.md) will take a list of identifiers in a list
file and pull the corresponding rows from a data file and write a new file. 
This is useful to subset a dataset. An option is available to generate a new
summary file (mean of all columns). 


## Database 

These applications are used to work with a local 
[Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store)
annotation database. See the [FAQ](FAQ.md) for details.

### db_setup.pl

[db_setup](apps/db_setup.pl) will assist in setting up a local SQLite database,
particularly with UCSC-derived annotation.

### db_types.pl

[db_types](apps/db_types.md) will list the available feature types in a local
database.

### get_intersecting_features.pl

[get_intersecting_features](apps/get_intersecting_features.md) will search a local
database and identify the intersecting features that overlap the query features or
coordinates. Useful for annotating novel genomic intervals of interest.

