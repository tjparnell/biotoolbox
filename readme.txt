README for Tim's BioToolBox

DESCRIPTION

These are a collection of various perl script programs for processing, 
converting, analyzing, and manipulating genomic data and/or features. 
These programs include the following capabilities
	
	* Tools to process raw Agilent microarray files
	* Tools to process next generation sequencing BAM files
	* Tools to prepare bioinformatic data (microarray & sequencing) into 
	  data formats (GFF v.3, wig, bigWig, bigBed, bam) for loading into a 
	  Bio::DB::SeqFeature::Store database
	* Tools for collecting bioinformatic data from a 
	  Bio::DB::SeqFeature::Store database relative to any genomic feature 
	  within the database
	* Tools for manipulating analyzing collected data, including 
	  generating graphs
	* Tools for converting data formats between popular formats


INSTALLATION

The programs have been developed and used on Mac OS X versions 10.5 and 10.6. The programs should work on any sane unix-like environment; Microsoft Windows compatability is not tested nor guaranteed. A working Perl environment is, of course, also required, and usually present in virtually every unix-like operating system distribution. Perl version 5.10.0 have been tested. 64-bit executable is desireable but not required, if only because some of the data files can get particularly large and I occasionally have had out-of-memory errors with 32-bit Perl.

Several Perl modules are required for these programs to work. Ideally I would have Build Module to check and help install these, but until then, here is a simple list of module distributions that I use. Note that dependencies of these, of course, should also be installed. All can be obtained through CPAN, which I heartily recommend using.
	
	* Bioperl
	* Bio::Graphics
	* Statistics::Lite
	* Statisitics::Descriptive
	* GD::Graph
	* GD::Graph::smoothlines
Optional, but highly recommended
	* Bio-SamTools
	* Bio-BigFile
	
In addition, some other libraries need to be installed which are required for the above to work.
	
	* samtools (http://samtools.sourceforge.net), required for Bio-SamTools, but also the executables for generally working with bam files
	* gd2, (http://www.boutell.com/gd/) required for the GD modules. On a Mac OS X system, this can be installed using Fink. Be sure that Fink is configured to use the same architecture as default Perl is using (64 bit on 10.6).
	* mysql (http://www.mysql.com) The Bio::DB::SeqFeature::Store database is configured to use a mysql database backend (although others could be used). On Mac OS X (the client version, not Server), this can also be installed through Fink.
	
In addition to the required modules, other programs that I've found to be quite useful and/or indespensible to my work are also highly suggested.

	* GBrowse (http://gmod.org) for visualizing the genomic data, particularly from Bio::DB::GFF or Bio::DB::SeqFeature::Store databases.
	* Cluster (http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/cluster/) for generating heirarchical and k-means clusters
	* Treeview (http://jtreeview.sourceforge.net) for visualizing the clusters

Since I work primarily on Mac OS X, I have compiled a a HowTo for setting up the system to work with these and other bioinformatics tools. It can be found at http://code.google.com/p/biotoolbox/wiki/SetupForMacOSX. Linux users may be able to glean useful information from that document as well.

USAGE

The programs are distributed as executable perl scripts. The programs may 
be called directly from the bin directory, or it may be added to your 
environment $PATH for easy execution. The programs will likely die prematurely 
if a required module is not found. Note that the tim_*_helper modules are found in the lib directory, and the perl scripts should be able to find them automatically if you leave the directory structure as is.

There is configuration file that may be customized for your particular installation. The default file is lib/tim_db_helper.cfg. It is a simple INI-style file that is used to set up database connection profiles, feature aliases, helper application locations, etc. The file is intended to be edited by users. There are three locations where the file may be stored: 1) in the lib directory along with the library modules, 2) in the user's home root directory, or 3) a custom location defined by the environment variable 'TIM_DB_HELPER'.

All programs are designed to be run from the command line. All require command line options. Executing the program without any options will present a synopsis of the options that are available. Most programs also have a --help option, which will display detailed information about the program and execution (usually by displaying the internal POD). The options are given in the long format (--help, for example), but may be shortened to single letters if the first letter is unique (-h, for example).

All of the programs can be called from a bash or another perl script. Some programs, manipulate_datasets.pl for example, also provide an interactive interface to allow for spontaneous work or when the exact index number or name of the dataset in the file or database is not immediately known.

Many of the programs are designed to input and output a file format referred to in colloquial tems as a 'tim data format'. This is merely a tabbed-delimited text format (unix line endings), where the rows represent genomic features, bins, etc. and the columns represent descriptive information and data. The first line in the table are the column headings. At the beginning of the file are metadata lines, prefixed by a #, which describe the content in each column using a simple key=value system. This metadata is useful to describe how and where the data was obtained, and what format it is in, for example log2. This frees up having to cram all the metadata into the filename. The metadata lines are disposable in most cases, and can be safely deleted before importing the file into another program, such as Excel. More detailed information can be found in the POD documentation of tim_file_helper.pm or online.

ONLINE DOCUMENTATION

Setting up Mac OS X for bioinformatics
http://code.google.com/p/biotoolbox/wiki/SetupForMacOSX

Setting up a Bio::DB::SeqFeature::Store database
http://code.google.com/p/biotoolbox/wiki/LoadingDatabase

Example Data collection
http://code.google.com/p/biotoolbox/wiki/ExampleFeatureDataCollection

Explanation of the text data report format
http://code.google.com/p/biotoolbox/wiki/TimDataFormat

Mapping SNPs
http://code.google.com/p/biotoolbox/wiki/MappingSNPs

Up to data list of programs
http://code.google.com/p/biotoolbox/wiki/ProgramList


