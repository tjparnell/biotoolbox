# NAME

bam2wig.pl

A program to convert Bam alignments into a wig representation file.

# SYNOPSIS

bam2wig.pl \[--options...\] &lt;file.bam>

bam2wig.pl --extend --rpm --mean --out file --bw file1.bam file2.bam

    Required options:
     -i --in <filename.bam>        repeat if multiple bams, or comma-delimited list
    
    Reporting options (pick one):
     -s --start                    record at 5' position
     -d --mid                      record at midpoint of alignment or pair
     -a --span                     record across entire alignment or pair
     -e --extend                   extend alignment (record predicted fragment)
     --cspan                       record a span centered on midpoint
     --smartcov                    record paired coverage without overlaps, splices
     --ends                        record paired endpoints
     --coverage                    raw alignment coverage
    
    Alignment reporting options:
     -l --splice                   split alignment at N splices
     -t --strand                   record separate strands as two wig files
     --flip                        flip the strands for convenience
     
    Paired-end alignments:
     -p --pe                       process paired-end alignments, both are checked
     -P --fastpe                   process paired-end alignments, only F are checked
     --minsize <integer>           minimum allowed insertion size (30)
     --maxsize <integer>           maximum allowed insertion size (600)
     --first                       only process paired first read (0x40) as single-end
     --second                      only process paired second read (0x80) as single-end
     
    Alignment filtering options:
     -K --chrskip <regex>          regular expression to skip chromosomes
     -B --blacklist <file>         interval file of regions to skip (bed, gff, txt)
     -q --qual <integer>           minimum mapping quality (0)          
     -S --nosecondary              ignore secondary (0x100) alignments (false)
     -D --noduplicate              ignore duplicate (0x400) alignments (false)
     -U --nosupplementary          ignore supplementary (0x800) alignments (false)
     --intron <integer>            maximum allowed gap (intron) size in bp (none)
     
     Shift options:
     -I --shift                    shift reads in the 3' direction
     -x --extval <integer>         explicit extension size in bp (default is to calculate)
     -H --shiftval <integer>       explicit shift value in bp (default is to calculate) 
     --chrom <integer>             number of chromosomes to sample (4)
     --minr <float>                minimum pearson correlation to calculate shift (0.5)
     --zmin <float>                minimum z-score from average to test peak for shift (3)
     --zmax <float>                maximum z-score from average to test peak for shift (10)
     -M --model                    write peak shift model file for graphing
     
    Score options:
     -r --rpm                      scale depth to Reads Per Million mapped
     -m --mean                     average multiple bams (default is addition)
     --scale <float>               explicit scaling factor, repeat for each bam file
     --fraction                    assign fractional counts to multi-mapped alignments
     --splfrac                     assign fractional count to each spliced segment
     --format <integer>            number of decimal positions (4)
     --chrnorm <float>             use chromosome-specific normalization factor
     --chrapply <regex>            regular expression to apply chromosome-specific factor
    
    Wig format:
     --bin <integer>               bin size for span or extend mode (10)
     --bdg                         bedGraph, default for span and extend at bin 1
     --fix                         fixedStep, default for bin > 1
     --var                         varStep, default for start, mid
     --nozero                      do not write zero score intervals in bedGraph
     
    Output options:
     -o --out <filename>           output file name, default is bam file basename
     -b --bw                       convert to bigWig format (supports bdg, fix, var)
     --bwapp /path/to/wigToBigWig  path to external converter (default searches \$PATH)
     -z --gz                       gzip compress text output 
     
    General options:
     -c --cpu <integer>            number of parallel processes (4)
     --temp <directory>            directory to write temporary files (output path)
     -V --verbose                  report additional information
     -v --version                  print version information
     -h --help                     show full documentation

# OPTIONS

The command line flags and descriptions:

## Input

- --in &lt;filename>

    Specify the input Bam alignment file. More than one file may be 
    specified, either with repeated options, a comma-delimited list, 
    or simply appended to the command. Bam files will be automatically 
    indexed if necessary.

## Reporting Options

- --start

    Specify that the 5' position should be recorded in the wig file.

- --mid

    Specify that the midpoint of the alignment (single-end) or fragment 
    (paired-end) will be recorded in the wig file.

- --span

    Specify that the entire span of the alignment (single-end) or 
    fragment (paired-end) will be recorded in the wig file. 

- --extend

    Specify that the alignment should be extended in the 3' direction 
    and that the entire length of the extension be recorded in the wig 
    file. The extension may be defined by the user or empirically 
    determined.

- --cspan

    Specify that a defined span centered at the alignment (single-end) 
    or fragment (paired-end) midpoint will be recorded in the wig file.
    The span is defined by the extension value.

- --smartcov

    Smart alignment coverage of paired-end alignments without 
    double-counting overlaps or recording gaps (intron splices). 

- --ends

    Record both endpoints of paired-end fragments, i.e. the outermost 
    or 5' ends of properly paired fragments. This may be useful with 
    ATAC-Seq, Cut&Run-Seq, or other cleavage experiments where you want 
    to record the locations of cutting yet retain the ability to filter 
    paired-end fragment sizes.

- --coverage

    Specify that the raw alignment coverage be calculated and reported 
    in the wig file. This utilizes a special low-level operation and 
    precludes any alignment filtering or post-normalization methods. 
    Counting overlapping bases in paired-end alignments are dependent on 
    the bam adapter (older versions would double-count).

- --position \[start|mid|span|extend|cspan|coverage\]

    Legacy option for supporting previous versions of bam2wig. 

## Alignment reporting options

- --splice

    Indicate that the bam file contains alignments with splices, such as 
    from RNASeq experiments. Alignments will be split on cigar N operations 
    and each sub fragment will be recorded. This only works with single-end 
    alignments, and is disabled for paired-end reads (just treat as single-end). 
    Only start and span recording options are supported.

- --strand

    Indicate that separate wig files should be written for each strand. 
    The output file basename is appended with either '\_f' or '\_r' for 
    both files. Strand for paired-end alignments are determined by the 
    strand of the first read.

- --flip

    Flip the strand of the output files when generating stranded wig files. 
    Do this when RNA-Seq alignments map to the opposite strand of the 
    coding sequence, depending on the library preparation method. 

## Paired-end alignments

- --pe

    The Bam file consists of paired-end alignments, and only properly 
    mapped pairs of alignments will be counted. Properly mapped pairs 
    include FR reads on the same chromosome, and not FF, RR, RF, or 
    pairs aligning to separate chromosomes. Both alignments are required 
    to be present before the pair is counted. The default is to treat 
    all alignments as single-end.

- --fastpe

    The Bam file consists of paired-end alignments, but to increase processing 
    time and be more tolerant of weird pairings, only the forward alignment is 
    required and considered; all reverse alignments are ignored. The default is 
    to treat all alignments as single-end.

- --minsize &lt;integer>

    Specify the minimum paired-end fragment size in bp to accept for recording. 
    Default is 30 bp.

- --maxsize &lt;integer>

    Specify the maximum paired-end fragment size in bp to accept for recording. 
    Default is 600 bp.

- --first

    Take only the first read of a pair, indicated by flag 0x40, and record as 
    a single-end alignment. No test of insert size or proper pair status is 
    made.

- --second

    Take only the second read of a pair, indicated by flag 0x80, and record as 
    a single-end alignment. No test of insert size or proper pair status is 
    made.

## Alignment filtering options:

- --qual &lt;integer>

    Set a minimum mapping quality score of alignments to count. The mapping 
    quality is a range from 0-255, with higher numbers indicating lower 
    probability of a mapping error. Multi-mapping alignments often have a 
    map quality of 0. The default is 0 (accept everything).

- --nosecondary

    Boolean flag to skip secondary alignments, indicated by the 
    alignment bit flag 0x100. Secondary alignments typically represent 
    alternative mapping locations, or multi-mapping events. By default,  
    secondary alignments are included. 

- ---noduplicate

    Boolean flag to skip duplicate alignments, indicated by the 
    alignment bit flag 0x400. Duplicates alignments may represent a PCR or 
    optical duplication. By default, duplicate alignments are included. 

- --nosupplementary

    Boolean flag to skip supplementary alignments, indicated by 
    the alignment bit flag 0x800. Supplementary alignments are typically 
    associated with chimeric fragments. By default, supplementary alignments 
    are included.

- --chrskip &lt;regex>

    Provide a regular expression to skip certain chromosomes. Perl-based 
    regular expressions are employed. Expressions should be quoted or 
    properly escaped on the command line. Examples might be 

        'chrM'
        'scaffold.+'
        'chr.+alt|chrUn.+|chr.+_random'

- --blacklist &lt;file>

    Provide a file of genomic intervals from which to exclude alignments. 
    Examples might include repeats, ribosomal RNA, or heterochromatic regions.
    The file should be any text file interpretable by [Bio::ToolBox::Data](https://metacpan.org/pod/Bio%3A%3AToolBox%3A%3AData) 
    with chromosome, start, and stop coordinates, including BED and GFF formats.
    Note that this only excludes overlapping alignments, and does not include 
    extended alignments.

- --intron &lt;integer>

    Provide a positive integer as the maximum intron size allowed in an alignment 
    when splitting on splices. If an N operation in the CIGAR string exceeds this 
    limit, the alignment is skipped. Default is 0 (no filtering).

## Shift options

- --shift

    Specify that the positions of the alignment should be shifted towards 
    the 3' end. Useful for ChIP-Seq applications, where only the ends of 
    the fragments are counted and often seen as separated discrete peaks 
    on opposite strands flanking the true target site. This option is 
    disabled with paired-end and spliced reads (where it is not needed). 

- --shiftval &lt;integer>

    Provide the value in bp that the recorded position should be shifted. 
    The value should be 1/2 the average length of the library insert size.
    The default is to automatically and empirically determine the 
    appropriate shift value using cross-strand correlation (recommended). 

- --extval &lt;integer>

    Manually set the length for reads to be extended. By default, the shift 
    value is determined empirically and extension is set to 2X the shift 
    value. This is also used for the cspan mode.

- --chrom &lt;integer>

    Indicate the number of sequences or chromosomes to sample when 
    empirically determining the shift value. The reference sequences 
    listed in the Bam file header are taken in order of decreasing 
    length, and one or more are taken as a representative sample of 
    the genome. The default value is 4. 

- --minr &lt;float>

    Provide the minimum Pearson correlation value to accept a shift 
    value when empirically determining the shift value. Enter a decimal value 
    between 0 and 1. Higher values are more stringent. The default 
    is 0.5.

- --zmin &lt;float>

    Specify the minimum z-score (or number of standard deviations) from 
    the chromosomal mean depth to test for a peak shift. Increase this 
    number to test for strong robust peaks, which give a better estimations 
    of the shift value. Default is 3.

- --zmax &lt;float> 

    Specify the maximum z-score (or number of standard deviations) from 
    the chromosomal mean depth to test for a peak shift. This excludes 
    erroneous peaks due to repetitive sequence alignments with high coverage. 
    Increase this number to include more robust peaks that can give a 
    better estimation of the shift value. Default is 10.

- --model

    Indicate that the shift model profile data should be written to 
    file for examination. The average profile, including for each 
    sampled chromosome, are reported for the forward and reverse strands, 
    as  well as the shifted profile. A standard text file is generated 
    using the output base name. The default is to not write the model 
    shift data.

## Score Options

- --rpm

    Convert the data to Reads (or Fragments) Per Million mapped. This is useful 
    for comparing read coverage between different datasets. The default is 
    no RPM conversion. 

- --scale &lt;float>

    Optionally provide your own scaling factor. This will be multiplied with 
    every position when generating the wig file. This may be combined with the 
    rpm factor as well. When combining multiple bam files, either a single scale 
    factor may be supplied for all files, or individual scale factors may be 
    supplied for each bam file. If supplying multiple, use the option multiple 
    times or give a comma-delimited list. The values should be in the same order 
    as the bam files. 

- --mean
- --separate

    When processing multiple bam files, this option will take the mean or average 
    across all bam files. Without this option, the bam files are simply added. 
    When combined with the rpm option, each bam file will be scaled separately 
    before taking the average.  

- --fraction

    Indicate that multi-mapping alignments should be given fractional counts 
    instead of full counts. The number of alignments is determined using the 
    NH alignment tag. If a read has 10 alignments, then each alignment is 
    given a count of 0.1. 

- --splfrac

    Indicate that spliced segments should be given a fractional count. This 
    allows a count to be assigned to each spliced segment while avoiding 
    double-counting. Best used with RNASeq spliced point data (--start or 
    \--mid); not recommended for --span.

- --format &lt;integer>

    Indicate the number of decimal postions reported in the wig file. This 
    is only applicable when rpm, scale, or fraction options are provided. 
    The default value is 4 decimal positions.

- --chrnorm &lt;float>

    Apply a normalization factor to the counts on specific chromosomes only. 
    Usually this is to normalize, for example, variable copy-number chromosomes, 
    such as transfected vector sequences, or haploid sex chromosomes. 

- --chrapply &lt;regex>

    Specify the Perl-based regular expression to match the chromosome(s) to 
    apply the specific normalization factor. For example, 'chrX$' to specify 
    the X chromosome only.

## Wig format

- --bin &lt;integer>

    Specify the bin size in bp for the output wig file. In general, specifying 
    a larger bin size will decrease the run time and memory requirements in 
    exchange for loss of resolution. The default for span, center span, or 
    extend mode is 10 bp; all other modes is 1 bp. 

- --bdg

    Specify that the output wig format is a bedGraph-style wig format. This is 
    the default format for extend, span, and cspan modes of operation.

- --fix

    Specify that the output wig format is in fixedStep wig format. This is the 
    default format for coverage mode of operation.

- --var

    Specify that the output wig format is in variableStep wig format. This is 
    the default format for start and midpoint modes of operation.

- --nozero

    When writing bedGraph format, skip (do not write) intervals with a value of 
    zero. Does not apply to fixedStep or variableStep formats.

## Output Options

- --out &lt;filename>

    Specify the output base filename. An appropriate extension will be 
    added automatically. By default it uses the base name of the 
    input file.

- --bw

    Specify whether or not the wig file should be further converted into 
    an indexed, compressed, binary BigWig file. The default is false.

- --bwapp /path/to/wigToBigWig

    Optionally specify the full path to the UCSC _wigToBigWig_ conversion 
    utility. The application path may be set in the `.biotoolbox.cfg` file 
    or found in the default environment `$PATH`, which makes this option 
    mostly unnecessary. 

- --gz

    Specify whether (or not) the output text file should be compressed with 
    gzip. Disable with `--nogz`. Does not apply to bigWig format.

## General options

- --cpu &lt;integer>

    Specify the number of parallel instances to run simultaneously. This requires 
    the installation of [Parallel::ForkManager](https://metacpan.org/pod/Parallel%3A%3AForkManager). With support enabled, the 
    default is 4. Disable multi-threaded execution by setting to 1. 

- --temp &lt;directory>

    Optionally specify an alternate temporary directory path where the temporary 
    files will be written. The default is the specified output file path, or the 
    current directory. Temporary files will always be written in a subdirectory of 
    the path specified with the template "bam2wigTEMP\_XXXX".

- --verbose

    Print extra informational statements during processing. The default is false.

- --version

    Print the version number.

- --help

    Display this POD documentation.

# DESCRIPTION

This program will enumerate aligned sequence tags and generate a wig, 
or optionally BigWig, file. Alignments may be counted and recorded 
in several different ways. Strict enumeration may be performed and 
recorded at either the alignment's start or midpoint position. 
Alternatively, either the alignment or fragment may be recorded 
across its span. Finally, a basic unstranded, unshifted, and 
non-transformed alignment coverage may be generated. 

Both paired-end and single-end alignments may be counted. Alignments 
with splices (e.g. RNA-Seq) may be counted singly or separately. 
Alignment counts may be separated by strand, facilitating analysis of 
RNA-Seq experiments. 

For ChIP-Seq experiments, the alignment position may be shifted 
in the 3 prime direction. This effectively merges the separate peaks 
(representing the ends of the enriched fragments) on each strand 
into a single peak centered over the target locus. Alternatively, 
the entire predicted fragment may be recorded across its span. 
This extended method of recording infers the mean size of the 
library fragments, thereby emulating the coverage of paired-end 
sequencing using single-end sequence data. The shift value is 
empirically determined from the sequencing data or 
provided by the user. If requested, the shift model profile may be 
written to file. 

The output wig file may be either a variableStep, fixedStep, or 
bedGraph format. The wig file may be further converted into a 
compressed, indexed, binary bigWig format, dependent on the 
availability of the appropriate conversion utilities. 

# RECOMMENDED SETTINGS

The type of wig file to generate for your Bam sequencing file can vary 
depending on your particular experimental application. Here are a few 
common sequencing applications and my recommended settings for generating 
the wig or bigWig file.

- Straight coverage

    To generate a straight-forward coverage map, similar to what most genome 
    browsers display when using a Bam file as source. **NOTE** that this mode 
    is pure raw coverage, and does not include any filtering methods. The other 
    modes allow alignment filtering.

        bam2wig.pl --coverage --in <bamfile>

- Smart paired-end coverage

    When you have paired-end alignments and need explicit alignment coverage
    without double-counting overlaps (as would occur if you counted as
    single-end span) or uncovered insertion (as would occur if you counted as 
    paired-end span) and not counting gaps (e.g. intron splices, as would occur
    with span mode), use the smart paired-end coverage mode. This properly
    assembles coverage from paired-end alignments taking into account overlaps
    and gaps.

        bam2wig --smartcov --in <bamfile>

- Single-end ChIP-Seq

    When sequencing Chromatin Immuno-Precipitation products, one generally 
    performs a 3 prime shift adjustment to center the fragment's end reads 
    over the predicted center and putative target. To adjust the positions 
    of tag count peaks, let the program empirically determine the shift 
    value from the sequence data (recommended). Otherwise, if you know 
    the mean size of your ChIP eluate fragments, you can use the --shiftval 
    option. 

    To evaluate the empirically determined shift value, be sure to include 
    the --model option to examine the profiles of stranded and shifted read 
    counts and the distribution of cross-strand correlations.

    Depending on your downstream applications and/or preferences, you 
    can record strict enumeration (start positions) or coverage (extend 
    position).

    Finally, to compare ChIP-Seq alignments from multiple experiments, 
    convert your reads to Reads Per Million Mapped, which will help to 
    normalize read counts.

        bam2wig.pl --start --shift --model --rpm --in <bamfile>
        
        bam2wig.pl --extend --model --rpm --in <bamfile>

- Paired-end ChIP-Seq

    If both ends of the ChIP eluate fragments are sequenced, then we do not 
    need to calculate a shift value. Instead, we will simply count at the 
    midpoint of each properly-mapped sequence pair, or record the defined 
    fragment span. 

        bam2wig.pl --mid --pe --rpm --in <bamfile>
        
        bam2wig.pl --span --pe --rpm --in <bamfile>

- Unstranded RNA-Seq

    With RNA-Sequencing, we may be interested in either coverage (generating 
    a transcriptome map) or simple tag counts (differential gene expression), 
    so we can count in one of two ways. 

    To compare RNA-Seq data from different experiments, convert the read 
    counts to Reads Per Million Mapped, which will help to normalize read 
    counts.

        bam2wig --span --splice --rpm --in <bamfile>
        
        bam2wig --mid --rpm --in <bamfile>

- Stranded, single-end RNA-Seq

    If the library was generated in such a way as to preserve strand, then 
    we can separate the counts based on the strand of the alignment. Note 
    that the reported strand may be accurate or flipped, depending upon 
    whether first-strand or second-strand synthesized cDNA was sequenced, 
    and whether your aligner took this into account. Check the Bam 
    alignments in a genome browser to confirm the orientation relative to 
    coding sequences. If alignments are opposite to the direction of 
    transcription, you can include the --flip option to switch the output.

        bam2wig --span ---splice --strand --rpm --in <bamfile>

        bam2wig --pos mid --strand --rpm --in <bamfile>
        

- Paired-end RNA-Seq

    Use the smart paired-end coverage mode to properly record paired-end 
    alignments with splice junctions. 

        bam2wig --smartcov --strand --rpm --in <bamfile>
        

# TEXT REPRESENTATION OF RECORDING ALIGNMENTS

To help users visualize how this program records alignments in a wig 
file, drawn below are 10 alignments, five forward and five reverse. 
They may be interpreted as either single-end or paired-end. Drawn 
below are the numbers that would be recorded in a wig file for various 
parameter settings. Note that alignments are not drawn to scale and 
are drawn for visualization purposes only. Values of X represent 10.

- Alignments

        ....>>>>>>.....................................<<<<<<.............
        .....>>>>>>..................................<<<<<<...............
        ........>>>>>>.......................................<<<<<<.......
        ........>>>>>>.........................................<<<<<<.....
        ..........>>>>>>............................................<<<<<<

- Starts

        ....11..2.1.......................................1.1.....1.1....1

- Midpoints

        ......11..2.1..................................1.1.....1.1....1...

- Stranded Starts

        F...11..2.1.......................................................
        R.................................................1.1.....1.1....1

- Span (Coverage)

        ....122244433311.............................112222111122221211111

- Mid Span (extend value 2)

        ......121.2211.................................1111....1111...11..

- Stranded Span

        F...122244433311..................................................
        R............................................112222111122221211111

- Shifted Starts (shift value 26)

        ........................1.1...11121.1..1..........................

- Shifted Span (shift value 26)

        ...................11222211112344365544411........................

- Extend (extend value 52)

        12223445789999XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX999887665321111

- Paired-End Midpoints

        ............................1...111..1............................

- Paired-End Mid span (extend value 6)

        ..........................111123333432111.........................

- Paired-End Span

        ....12224455555555555555555555555555555555555555555443333332211111

# AUTHOR

    Timothy J. Parnell, PhD
    Huntsman Cancer Institute
    University of Utah
    Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
