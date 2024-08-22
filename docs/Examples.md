# Bio::ToolBox

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## Application Examples

The BioToolBox package comes complete with a suite of high-quality, production-ready 
applications (Perl scripts) ready for a variety of analyses. A sampling of what 
can be done include the following:

- Annotated feature collection and selection
- Data collection and scoring for features
- Data file format manipulation and conversion
- Low-level processing of sequencing data into customizable wig representation

The [Applications](Applications.md) page describes each one.

The following are just a few examples of highlighted scripts and their usage for 
solutions to common situations. Refer to the applications' documentation for details
on options shown and not shown.

### Convert or extract gene annotation

Gene annotation from [UCSC](http://genome.ucsc.edu) frequently come in UCSC 
[specific formats](http://genome.ucsc.edu/FAQ/FAQformat.html#format9), including
refFlat and genePred, whereas gene annotation from [Ensembl](http://ensembl.org) 
frequently come in GTF or GFF3 formats. At some point, some tool will specifically 
need annotation in a different format than what you have (since most of 
bioinformatics seem to be converting from one format to another), 
[get_features.pl](apps/get_features.md) can help here.

Simple conversion:

	 get_features.pl --in file.gff3.gz --refflat

Extract only protein-coding genes, collapsing alternate transcripts into one 
combined meta-transcript, as a GTF:

	 get_features.pl --in file.gff3.gz --feature mRNA --collapse --gtf

Pull specific biotype `lincRNA` transcripts from standard chromosomes:

	 get_features.pl --in file.gff3 --feature transcript --biotype lincRNA \
	 --chrskip 'contig|scaffold|unmapped' --gtf

Storing annotation in a SQL database isn't used frequently any more (at least by
the author), but support is still present using the Bio::Perl adapter
[Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store).
To pull annotation from a local SQLite database file for use in downstream 
applications:

	 get_features.pl --db annotation.sqlite --feature gene --out genes.txt 

### Extract specific gene regions

Not all parts of genes that you might be interested in are explicitly defined 
or encoded into a gene annotation table; rather, they are inferred. The 
[get_gene_regions.pl](apps/get_gene_regions.md) script can 
extract these inferred regions.

Extract only the alternate exons from multi-transcript genes:

	get_gene_regions.pl --region altExon --in genes.gtf --out altExons

Extract all the transcription start sites of protein coding transcripts, but only 
report those of alternate transcripts once if they're within 200 bp, and 
expand the site by 200 bp in both directions, as a bed file:

	get_gene_regions.pl --region tss --in genes.gtf --bed --out tss200.bed \
	--feature gene --transcript mRNA --unique --slop 200 --start -200 --end 200 

### Generate wig file representation of a bam file

Lots of existing programs can generate a coverage file from a bam file, but 
[bam2wig.pl](apps/bam2wig.md) will generate wig files in 
every which way imaginable. Some common scenarios include:

ChIPSeq normalized fragment coverage with empirical fragment size determination
(single-end):

	bam2wig.pl --extend --shift --model --rpm --bw --in file1.bam --out file1.bw

ChIPSeq normalized fragment coverage with explicit extension, explicit scaling factor, 
averaged over 3 replicates, excluding repeat elements and certain chromosomes:

	bam2wig.pl --extend --extval 200 --bw --out file.bw --rpm --mean \
	--scale 0.0133 --blacklist repeats.bed.gz --chrskip 'chrM|contig|scaffold' \
	--cpu 12 replicate1.bam replicate2.bam replicate3.bam 

ATACSeq cut sites, where cut sites are represented by 50 bp pileup centered over the 
cut site at the 5' end:

	bam2wig.pl --shift --shiftval -25 --extend --extval 50 --in file1.bam

RNASeq stranded coverage, where 2 wig files are generated for each strand, 
ignoring multi-mappers, with output basename of 'experiment':

	bam2wig.pl --span --splice --strand --nosecondary --rpm --bw \
	--in file1.bam --out experiment

Simple point representation of all fragments for counting:

	bam2wig.pl --start --bw --in file1.bam --out file1.bw

### Simple data collection

Sometimes you just need simple, summarized scores over a list of genes or regions. 
The [get_datasets.pl](apps/get_datasets.md) script can collect 
data summarized data from bigWig, bigBed, Bam, USeq, and more.

Collect mean values over Bed intervals:

	get_datasets.pl --in file.bed --method mean score1.bw score2.bw score3.bw

In general, [bam2wig.pl](apps/bam2wig.md) has much better 
controls for filtering and scoring bam alignments, but in a pinch 
[get_datasets.pl](apps/get_datasets.md)
will also work with bam files. A method of `mean` or `median` will work with 
raw coverage, while `count`, name count (`ncount`), or precise count (`pcount`) 
will work with alignment counts.

	get_datasets.pl --in file.bed --method ncount file1.bam file2.bam

Collect mean sense RNASeq coverage from a collection of stranded bigWig files in 
a folder (as shown above) over transcript exons:

	get_datasets.pl --in annotation.gtf --feature transcript --subfeature exon \
	--strand sense --ddb /path/to/stranded/bigWigs/ --data experiment --out file.txt
    
### Collect data around a reference point

Rather than collect a single score for region, you need to collect a range of scores 
relative to a single reference point, for example a transcription start site.
The [get_relative_data.pl](apps/get_relative_data.md) script will 
collect data in a defined number of windows on both sides of a specified reference 
point.

Collect mean values in twenty 500 bp windows (+/- 10 kb) from the transcription 
start site:

	get_relative_data.pl --in mRNA.gtf --method mean --win 500 --num 20 \
	--pos 5 --data scores.bw --out mRNA_TSS_10kb.txt

### Collect profile data across an region

Instead of collecting a single score for a region, you need to generate a profile 
across a region. The [get_binned_data.pl](apps/get_binned_data.md)
script will collect scores in defined fractional windows across a region, setting 
length to each region an artificial 1000 bp. 

To collect an average profile expression across transcripts, excluding introns, 
and generate a summary (average profile):

	get_binned_data.pl --in mRNA.gtf --method mean --subfeature exon --bins 50 \
	--strand sense --ddb /path/to/stranded/bigWigs/ --data experiment \
	--sum --out mRNA_sense_profile.txt

### Manipulate collected dataset file

Inevitably, you need to manipulate your data files: add or remove columns, perform 
mathematical functions on columns, etc. You could open the file in a spreadsheet 
application, or in an R session, or get creative with `cut`, `sed`, or `awk` commands 
in the terminal, or just simply run 
[manipulate_datasets.pl](apps/manipulate_datasets.md).

	manipulate_datasets.pl file.txt

This script is designed to run interactively, allowing you to choose functions via 
letter options from a menu. Press `m` for the menu, `q` to save the file and quit, 
and `Q` to just quit without saving. 

Alternatively, single functions can be specified on the command line for a
programmatic approach in a `bash` script, for example. 

	manipulate_datasets.pl --in file.txt --func multiply --index 5 --target 0.0133

Be careful when specifying functions on the command line: the program will overwrite
the input file unless an output file is specified with `--out`. Also, manipulations
to columns of data are automatically written to new columns unless specified otherwise
with `--place n`.


