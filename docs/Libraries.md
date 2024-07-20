# LIBRARIES

Several library modules are included in this distribution. The following is a brief
description of the primary user-oriented libraries that are available. Links will 
take you to the full documentation on MetaCPAN.

- [Bio::ToolBox::Data](https://metacpan.org/pod/Bio::ToolBox::Data)

    This is the primary library module for working with a table of data, 
    either generated as a new list from a database of annotation, or 
    opened from a tab-delimited text file, for example a BED file of 
    regions. Columns and rows of data may be added, deleted, or manipulated 
    with ease. 

    Additionally, genomic data may be collected from a wide variety of 
    sources using the information in the data table. For example, scoring 
    microarray or sequencing data for each interval listed in the data 
    table.

    This module uses an object-oriented interface. Many of the methods 
    and API will be familiar to users of [Bio::Perl](https://metacpan.org/pod/Bio::Perl).

- [Bio::ToolBox::Data::Feature](https://metacpan.org/pod/Bio::ToolBox::Data::Feature)

    This is the object class for working with individual rows in a table 
    of data. It provides a number of conventions for working with the rows 
    in a standard fashion, for example returning the start column value  
    regardless of which column it is or whether the table is bed or gff or 
    an arbitrary text file. A number of convenience methods are present for 
    collecting data from data files. This module is not used directly by the 
    user, but its objects are returned when using
    [Bio::ToolBox::Data](https://metacpan.org/pod/Bio::ToolBox::Data) iterators.

- [Bio::ToolBox::Parser](https://metacpan.org/pod/Bio::ToolBox::Parser)

    This is the base class for working with common annotation formats. By
    default, these will slurp an entire genomic annotation file into memory
    in a reasonably short amount of time, including parsing gene annotation
    into a nested, hierarchical structure. To minimize memory usage, the 
    [Bio::ToolBox::SeqFeature](https://metacpan.org/pod/Bio::ToolBox::SeqFeature)
    is used by default. There are three format-specific parsers.

    - [Bio::ToolBox::Parser::bed](https://metacpan.org/pod/Bio::ToolBox::Parser::bed)

        This parses [BED file](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) and 
        related formats, including BED files with 3-12 columns 
        (BED3, BED6, BED12, and in between), [bedGraph](http://genome.ucsc.edu/FAQ/FAQformat.html#format1.8), 
        [narrowPeak](http://genome.ucsc.edu/FAQ/FAQformat.html#format12), and 
        [broadPeak](http://genome.ucsc.edu/FAQ/FAQformat.html#format13). For 
        proper BED12 files, transcripts are parsed with child subfeatures including exon 
        and CDS subfeatures.

    - [Bio::ToolBox::Parser::gff](https://metacpan.org/pod/Bio::ToolBox::Parser::gff)

        This parses both [GTF](http://mblab.wustl.edu/GTF22.html) and 
        [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) 
        file formats. Unlike most other GFF parsers 
        that work line-by-line only, this maintains parent and child hierarchical 
        relationships as parent feature and child subfeatures. To further maintain 
        control and reduce unnecessary parsing, unwanted feature types can be 
        selectively skipped.

    - [Bio::ToolBox::Parser::ucsc](https://metacpan.org/pod/Bio::ToolBox::Parser::ucsc)

        This parses various [UCSC file formats](http://genome.ucsc.edu/FAQ/FAQformat.html#format9), 
        including different refFlat, GenePred, and knownGene flavors. 
        Genes, transcripts, and exons are assembled into 
        hierarchical child-parent relationships as desired.

- [Bio::ToolBox::SeqFeature](https://metacpan.org/pod/Bio::ToolBox::SeqFeature)

    This is a fast, lean, simple object class for representing genomic features. 
    It supports, for the most part, the [Bio::SeqFreatureI](https://metacpan.org/pod/Bio::SeqFreatureI) 
    and [Bio::RangeI](https://metacpan.org/pod/Bio::RangeI) API 
    interface without the dependencies. It uses an unorthodox blessed-array object 
    structure, which provides measurable improvements in memory consumption and 
    speed when loading thousands of annotated SeqFeature objects (think hg19 or hg38 
    annotation). 

- [Bio::ToolBox::GeneTools](https://metacpan.org/pod/Bio::ToolBox::GeneTools)

    This is a collection of exportable functions for working with
    [Bio::SeqFeatureI](https://metacpan.org/pod/Bio::SeqFeatureI) 
    compliant objects representing genes and transcripts. It works with objects derived 
    from one of the ["Annotation parsers"](#annotation-parsers) or a 
    [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store) 
    database. The functions make hard things easy, such as identifying whether a 
    transcript is coding or not (is it encoded in the `primary_tag` or `source_tag` or 
    GFF attribute or does it have `CDS` subfeatures?), or identify the alternative exons 
    or introns of a multi-transcript gene, or pull out the `5'` UTR (which is likely 
    not explicitly defined in the table).

