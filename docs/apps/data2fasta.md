# NAME

data2fasta.pl

A program to retrieve sequences from a list of features

# SYNOPSIS

data2fasta.pl \[--options...\] &lt;filename>

    File Options:
    -i --in <filename>                input file: txt, gff, bed, ucsc, vcf, etc
    -o --out <filename>               output file name
    
    Database:
    -d --db <name|fasta>              annotation database with sequence or fasta
    
    Feature selection:
    -f --feature <text>               feature when parsing gff3, gtf, or ucsc input
    -u --subfeature [exon|cds|        collect over subfeatures 
          5p_utr|3p_utr] 
    
    Column indices:
    -n --name --id <index>            name or ID column
    -s --seq <index>                  column with sequence
    -c --chr <index>                  chromosome column
    -b --begin --start <index>        start coordinate column
    -e --end --stop <index>           stop coordinate column
    -t --strand <index>               strand column
    -x --extend <integer>             extend coordinates in both directions
    --desc <index>                    description column
    
    Fasta output options:
    --cat                             concatenate all sequences into one
    --pad <integer>                   pad concatenated sequences with Ns
    
    General options:
    -z --gz                           compress output fasta file
    -v --version                      print version and exit
    -h --help                         show extended documentation

# OPTIONS

The command line flags and descriptions:

## File options

- --in &lt;filename>

    Specify the input data file. The file may be a tab-delimited text file 
    with coordinate columns for fetching genomic sequence. Alternatively it 
    may be an annotation file such as GFF3, GTF, refFlat, genePred, etc, 
    in which case sequence may be selected from certain features and 
    subfeatures, for example mRNA and CDS. **Note** that no collapsing of 
    redundant or overlapping subfeatures is performed; see [get\_features.pl](https://metacpan.org/pod/get_features.pl). 
    Finally, text files with sequence in a column, for example oligo sequences, 
    may be used, skipping the need for database sequence retrieval.
    The file may be compressed with gzip.

- --out &lt;filename>

    Specify the output filename. By default it uses the input file basename.

## Database

- --db &lt;name|fasta>

    Provide the name of an uncompressed Fasta file (multi-fasta is ok) or 
    directory containing multiple fasta files representing the genomic 
    sequence. If Bam file support is available, then the fasta will be 
    indexed and searched with a fasta index .fai file. If not, then the 
    fasta can by indexed by the older [Bio::DB::Fasta](https://metacpan.org/pod/Bio%3A%3ADB%3A%3AFasta) adapter, which 
    also supports a directory of multiple fasta files. If the index is 
    not present, then the parent directory must be writeable.
    Alternatively, the name of a [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) 
    annotation database that contains genomic sequence may also be provided. 
    The database name may be obtained from the input file metadata. 
    Required only if collecting sequence from genomic coordinates.

## Feature selection

- --feature &lt;text>

    When parsing a gene annotation file such as a GFF3, GTF, or UCSC format 
    file, provide a feature type to select features if desired.

- --subfeature \[exon|cds|5p\_utr|3p\_utr\]

    When collecting from subfeatures, indicate the subfeature type from 
    list available. No merging of overlapping or redundant subfeatures 
    is performed here. See [get\_features.pl](https://metacpan.org/pod/get_features.pl).

## Column indices

- --name --id &lt;column\_index>

    Optionally specify the index for the name or ID column. It may be 
    automatically determined from the column header.

- --seq &lt;column\_index>

    Optionally specify the index for the sequence column. It may be 
    automatically determined from the column header.

- --chr &lt;column\_index>

    Optionally specify the index for the chromosome column. It may be 
    automatically determined from the column header.

- --start &lt;column\_index>
- --begin &lt;column\_index>

    Optionally specify the index for the start position column. It may be 
    automatically determined from the column header.

- --stop &lt;column\_index>
- --end &lt;column\_index>

    Optionally specify the index for the stop position column. It may be 
    automatically determined from the column header.

- --strand &lt;column\_index>

    Optionally specify the index for the strand column. It may be 
    automatically determined from the column header.

- --extend &lt;integer>

    Optionally provide the number of extra base pairs to extend the start 
    and stop positions. This will then include the given number of base 
    pairs of flanking sequence from the database. This only applies when 
    sequence is obtained from the database.

- --desc &lt;column\_index>

    Optionally specify the index of the description column. It may be 
    automatically determined from the column header.

## Fasta output options

- --cat

    Optionally indicate that all of the sequences should be concatenated 
    into a single Fasta sequence. The default is to write a multi-fasta 
    file with separate sequences.

- --pad &lt;integer>

    When concatenating sequences into a single Fasta sequence, optionally 
    indicate the number of 'N' bases to insert between the individual 
    sequences. The default is zero.

## General options

- --gz

    Specify whether (or not) the output file should be compressed with gzip.

- --version

    Print the version number.

- --help

    Display this POD documentation.

# DESCRIPTION

This program will take a tab-delimited text file (BED file, 
for example) and generate either a multi-sequence fasta file containing the 
sequences of each feature defined in the input file, or optionally a single 
concatenated fasta file. If concatenating, the individual sequences may be 
padded with the given number of 'N' bases. 

This program has two modes. If the name and sequence is already present in 
the file, it will generate the fasta file directly from the file content.

Alternatively, if only genomic position information (chromosome, start, 
stop, and optionally strand) is present in the file, then the sequence will 
be retrieved from a database. Multiple database adapters are supported for 
indexing genomic fastas, including the [Bio::DB::HTS](https://metacpan.org/pod/Bio%3A%3ADB%3A%3AHTS) package, the 
[Bio::DB::Sam](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASam) package, or the BioPerl [Bio::DB::Fasta](https://metacpan.org/pod/Bio%3A%3ADB%3A%3AFasta) adapter. Annotation 
databases such as [Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio%3A%3ADB%3A%3ASeqFeature%3A%3AStore) are also supported.  
If strand information is provided, then the sequence reverse complement 
is returned for reverse strand coordinates.

# AUTHOR

    Timothy J. Parnell, PhD
    Dept of Oncological Sciences
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
