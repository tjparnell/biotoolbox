# Bio::ToolBox - FAQ

|[Home](ReadMe.md)|[Install](AdvancedInstallation.md)|[Libraries](Libraries.md)|[Applications](Applications.md)|[Examples](Examples.md)|[FAQ](FAQ.md)|

## Frequently Asked Questions

These may or may not have actually been asked, but it's a collection of hints that the 
programmer understands but a casual user might not, as well as rationale.

- Programs don't recognize a UCSC gene table (refFlat, knownGene, genePred, etc)

	UCSC doesn't have official file extensions, and their downloads page just 
	have `.txt.gz` extensions. Furthermore, they don't have proper column headers. Downloads 
	from the table browser will stick a header line, prefixed with a `#` but no space
	between it and the first word. I have work arounds to detect those headers, but 
	what about the files from the download page?
	
	Programs that are designed to potentially interpret a gene table, such as
	[get_datasets](apps/get_datasets.md), will "taste" a file for potential UCSC
	formats based on the number of columns present. Other programs, like
	[manipulate_datasets](apps/manipulate_datasets.md), don't bother and will take
	the first line as the header regardless. 
	
	Some programs accept a `--noheader` flag, and it will insert dummy column headers.
	
	Otherwise, you can help yourself by changing the extension from `.txt` to something 
	more descriptive, like `.refflat`, `.genepred`, `.knowngene`, or even the most 
	generic `.ucsc`. Don't forget the `.gz` if it's compressed.

- What is the difference between Start and Start0?

	`Start` represents the starting coordinate in 1-base coordinate. `Start0` represents 
	the starting coordinate in 0-base coordinate. It's just my way of differentiating 
	between the two.
	
	Many annotation formats come in two flavors of coordinate system: 1-base system
	(counting each nucleotide in a sequence starting at 1) or 0-base (or interbase) system
	(counting between bases, hence starting at 0). The GFF family of annotation file
	formats (including GTF and GFF3) use 1-base. The UCSC family of annotation formats
	(BED, refFlat, genePred, etc) use 0-base. SAM files are 1-based, but binary BAM files
	are internally 0-based, while VCF files are 1-based. In other words, every format is
	different. The [BioPerl](https://bioperl.org) libraries, of which much of BioToolBox
	was initially based on, uses 1-base for everything. BioToolBox inherently transforms
	0-based coordinates to 1-base formats internally, at least when it is aware of what
	the file is using, hence the purpose of naming columns differently.

- Why do so many programs reference a database and how do I use one?

	In the early days of BioToolBox, much of the analysis was based on
	[BioPerl](https://bioperl.org) databases, notably
	[Bio::DB::SeqFeature::Store](https://metacpan.org/pod/Bio::DB::SeqFeature::Store),
	where annotation as well as datasets (microarray values) were stored. These were SQL
	databases, backed by either MySQL or SQLite. These are still supported, although less
	so as annotation files can now be parsed on the fly or datasets stored in bigWig or
	Bam databases. 
	
	For annotation, working with a database can be arguably faster, especially when 
	working with an annotation set over and over again. Use the BioPerl script, 
	`bp_seqfeature_load.pl`, to load GFF3 annotation into a Bio::DB::SeqFeature::Store 
	database. Unless you're working in a multi-user environment, I recommend using a 
	SQLite file backend as it is much simpler. For MySQL environments, create a simple, 
	read-only account with a generic password and put it in the `.biotoolbox.cfg` file 
	in your home directory. This is a plain-text file, which is why the account should 
	be read-only. Numerous databases and their paths, connections, and/or credentials 
	can be stored here.

- What is the `.biotoolbox.cfg` file?

	A small INI-style text file used to store various general configurations. It was
	initially a place to store database credentials. It can also be used to store the
	paths of various external utility applications, such as `wigToBigWig`, when they are
	not in your environment `PATH`. It is mostly legacy, and no longer needed in most
	cases. See the
	[Bio::ToolBox::db_helper::config](https://metacpan.org/pod/Bio::ToolBox::db_helper::config) 
	documentation.

- What is a BigWigSet, and how do I use one?

	A BigWigSet is a collection of bigWig files in a folder, with which you can associate 
	metadata by using a `metadata.txt` file, including but not limited to sample name, 
	source, strand, etc. It was a concept developed by Lincoln Stein, implemented in his 
	[Bio::DB::BigWigSet](https://metacpan.org/pod/Bio::DB::BigWigSet) adapter, for use with 
	the [GBrowse](http://gmod.org/wiki/GBrowse) genome viewer. It has been 
	re-implemented in 
	[Bio::ToolBox::db_helper::big](https://metacpan.org/pod/Bio::ToolBox::db_helper::big)
	because of its usefulness. 
	
	A simple BigWigSet, consisting of one or more bigWig files in a directory, can be
	specified by providing the path to the directory as an argument to the `--ddb` option
	in data collection scripts. The `metadata.txt` file is optional. If you prepare one,
	it is a simple INI-style text file, with block headers consisting of the file names.
	For example, 
	    
	    [file1.bw]
		name = mydata
		type = ChIPSeq
	
		[file2.bw]
		name = mydata2
		type = ChIPSeq

- Why do you prefer bigWig over bam files?

	As currently implemented, the bam adapter for use in data collection scripts is
	pretty limited in terms of alignment filtering and manipulating and recording
	alignment positions. The [bam2wig](apps/bam2wig.md) script, on the other hand,
	has extensive filtering and reporting capabilities. The bigWig adapters are
	generally much faster than bam adapters in collecting data, and you can always
	visualize exactly what you are collecting in a genome browser. If your workflow
	requires complicated alignment filtering and counting without a bigWig conversion
	step, there are other tools out there....

- How do you make the programs run faster?

	Perl is relatively easy to code, but it's by no means a speed demon. When possible, 
	I fork the main process into child processes to split the work into pieces for 
	parallel execution, which helps considerably. Some things, however, just can't be 
	split, e.g. parsing a GFF file. For iterating over bam files, each chromosome 
	reference is split into a separate fork (no love for uni-chromosomal bacteria). It's 
	not whole-integer multiples faster, but it's better than single-thread. 
	
	For data collection, give more cpu cores when possible. The data collection
	scripts self-limits the number of forks based on how big the input data file is.
	The OS will generally take care of disk I/O and concurrency. If you're constantly
	loading large annotation files, consider making them smaller (take just what you
	need with [get_features](apps/get_features.md) or
	[get_gene_regions](apps/get_gene_regions.md)), use simpler annotation formats
	(BED and refFlat are fastest), or import them into an annotation SQLite database
	file (see above).
	
	There is an (unofficial) [effort](http://perl11.org) to make Perl faster by tweaking 
	the internals. BioToolBox will install under [cperl](http://perl11.org/cperl/), 
	although getting the prerequisites installed is not a trivial task. The speed gain 
	is modest.

- Why do you fork instead of using threads?

	It's easier? Forking a child process is less complicated, as memory is automatically 
	set up as read-only, copy-on-write, without the headaches of initializing variables as 
	shareable. More importantly, however, is that this package relies heavily on 
	external libraries and XS compiled-C code, which are not natively thread-safe. 
	Hence the reason for always re-opening database handles (SQL, Bam, bigWig, etc) upon 
	forking.

- Why do your column indexes begin at 0 instead of 1?

	**UPDATE**: With version 2 release of Bio::ToolBox, all column indexes are now
	indexed at 1 instead of 0. 
	
	Inexperience, laziness, and legacy code. In Perl (and many programming languages), 
	an array of variables are indexed beginning at 0. Since data tables are loaded as 
	an array (rows) of arrays (columns), the first column is 0. When I first started as 
	an (inexperienced) programmer, it was too hard (or I was too lazy?) to always shift 
	an index by one to make it appear as if the first column began at position 1. I 
	could do so now, but it might break (some?) legacy code. So far, no one has 
	complained.

- Why Perl and not Python?

	This project traces back all the way to 2005. At that time, the
	[BioPerl](https://bioperl.org) project was already in full swing, with ready-to-use
	modules for working with annotation, sequences, and storing data (annotation,
	microarray scores, etc) in databases with fully-featured APIs. The BioPython and
	BioJava equivalents paled in comparison at the time, so it was a "no brainer" decision
	to go with Perl. I'm too invested to start over. 


