README for Tim's BioToolBox

DESCRIPTION

These are a collection of various perl script programs for processing, 
converting, analyzing, and manipulating genomic data and/or features. 
These programs include the following capabilities
	
	* Tools for converting data formats between popular formats, including
	  GFF3, wig, bed, BAM, bigWig, bigBed, Bar, and USeq
	* Collecting data relative to any annotated genomic feature
	* Manipulation and analysis of collected data
	* Generating simple graphs of collected data
	* Processing raw Agilent microarray files
	* Simple processing of Next Generation Sequencing alignment files
	* Precise mapping of nucleosomes from paired-end sequencing



INSTALLATION

Platform 
The programs have been developed and used on Mac OS X versions 10.5 and
10.6. The programs should work on any sane unix-like environment, including
Linux; Microsoft Windows compatability is not tested nor guaranteed.

Perl 
A working Perl environment is, of course, also required, and usually
present in virtually every unix-like operating system distribution. Perl
version 5.10.0 has been tested. 64-bit executable is desireable but not
required, if only because some of the data files can get particularly
large, leading to out-of-memory errors with 32-bit Perl.

No install 
There is no installation script or compilation required. The scripts are
designed to be run as is from the scripts directory. Biotoolbox-specific
libraries are included in the lib directory and should be automatically
found if the directory structure is maintained.

Dependencies 
Several different Perl modules are required for specific programs to work.
To check for these dependencies, run the script 'check_dependencies.pl'
located in the biotoolbox root directory. The script will assist in
installing or upgrading missing or out-of-date modules through CPAN, with
only one exception. Depending upon your installation, you may need to
execute with administrative privelages: run 'sudo check_dependencies.pl'
and enter the password.

Note that not all of modules are required. The best way to find out if a
dependency is absolutely required is to try running the desired script.
Perl will appropriately complain if it can't find the module.

Mac OS X 
Since I work primarily on Mac OS X, I have compiled a a HowTo for
setting up the system to work with these and other bioinformatics tools. It
can be found at http://code.google.com/p/biotoolbox/wiki/SetupForMacOSX.
Linux users may be able to glean useful information from that document as
well.



USAGE

Configuration 
There is a configuration file that may be customized for your particular
installation. The default file is lib/biotoolbox.cfg. It is a simple
INI-style file that is used to set up database connection profiles, feature
aliases, helper application locations, etc. The file is intended to be
edited by users. More documentation can be found within the file itself.

There are three locations where the file may be stored: 1) in the lib
directory along with the library modules, 2) in the user's home root
directory, or 3) a custom location defined by the environment variable
'BIOTOOLBOX'.

Execution 
All biotoolbox scripts are designed to be run from the command line or
executed from another script. Some programs, for example
manipulate_datasets.pl, also provide an interactive interface to allow for
spontaneous work or when the exact index number or name of the dataset in
the file or database is not immediately known.

Help 
All scripts require command line options for execution. Executing the
program without any options will present a synopsis of the options that are
available. Most programs also have a --help option, which will display
detailed information about the program and execution (usually by displaying
the internal POD). The options are given in the long format (--help, for
example), but may be shortened to single letters if the first letter is
unique (-h, for example).

File Formats 
Many of the programs are designed to input and output a file format
referred to in colloquial tems as a 'tim data format'. This is merely a
tabbed-delimited text format (unix line endings), where the rows represent
genomic features, bins, etc. and the columns represent descriptive
information and data. The first line in the table are the column headings.

At the beginning of the file are metadata lines, prefixed by a #, which
describe the content in each column using a simple key=value system. This
metadata is useful to describe how and where the data was obtained, and
what format it is in, for example log2. This frees up having to cram all
the metadata into the filename. The metadata lines are dispensible in most
cases, and can be safely deleted before importing the file into another
program, such as Excel. More detailed information can be found in the POD
documentation of tim_file_helper.pm or online.



ONLINE DOCUMENTATION

Setting up Mac OS X for bioinformatics
http://code.google.com/p/biotoolbox/wiki/SetupForMacOSX

Setting up a Bio::DB::SeqFeature::Store database
http://code.google.com/p/biotoolbox/wiki/LoadingDatabase

Description of supported data file formats
http://code.google.com/p/biotoolbox/wiki/DataFileFormats

Example Data collection
http://code.google.com/p/biotoolbox/wiki/ExampleFeatureDataCollection

Explanation of the text data report format
http://code.google.com/p/biotoolbox/wiki/TimDataFormat

Mapping SNPs 
http://code.google.com/p/biotoolbox/wiki/MappingSNPs

Up to data list of programs
http://code.google.com/p/biotoolbox/wiki/ProgramList


