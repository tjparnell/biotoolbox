README for Tim's BioToolBox

DESCRIPTION

These are a collection of various perl script programs for processing, 
converting, analyzing, and manipulating genomic data and/or features. 
These programs include the following capabilities
	
	* Tools for converting data formats between popular formats, including
	  GFF3, wig, bed, BAM, bigWig, bigBed, and Bar
	* Collecting data relative to any annotated genomic feature
	* Manipulation and analysis of collected data
	* Generating simple graphs of collected data
	* Processing raw microarray files
	* Simple processing of Next Generation Sequencing alignment files
	* Precise mapping of nucleosomes from paired-end sequencing



INSTALLATION

*Platform 
The programs have been developed and used on Mac OS X versions 10.5 through
10.7. The programs should work on any sane unix-like environment, including
Linux; Microsoft Windows compatability is not tested nor guaranteed.

*Perl 
A working Perl environment is, of course, also required, and usually
present in virtually every unix-like operating system distribution. Perl
versions 5.10 and 5.12 have been tested. 64-bit executable is desireable, but 
not required, if only because some of the data files can get particularly
large, leading to out-of-memory errors with 32-bit Perl.

*Perl Module Dependencies 
Several different Perl modules are required for specific programs to work, 
most notably BioPerl, among a few others.

To check for these dependencies, run the script 'check_dependencies.pl'
located in the biotoolbox root directory. The script will helpfully assist 
in installing or upgrading missing or out-of-date modules through CPAN, with
only one exception. Depending upon your installation, you may need to
execute with administrative privelages: run 'sudo check_dependencies.pl'
and enter the password.

Note that not all of the modules are required. The best way to find out if 
a dependency is required is to try running the desired script. Perl will 
appropriately complain if it can't find the module.

*No install 
There is no installation script or compilation required. The scripts are
designed to be run as is from the scripts directory. Biotoolbox-specific
libraries are included in the lib directory and should be automatically
found if the directory structure is maintained. The biotoolbox libraries 
were not intended for general consumption (they do not have a friendly 
object-oriented API), but they are fully documented in case someone dares 
to use them.

For those who need additional guidance in setting up their system, there 
is an online HowTo manual located here: 
http://code.google.com/p/biotoolbox/wiki/BioToolBoxSetUp



USAGE

*Configuration 
There is a configuration file that may be customized for your particular
installation. The default file is biotoolbox.cfg. It is a simple
INI-style file that is used to set up database connection profiles, feature
aliases, helper application locations, etc. The file is intended to be
edited by users. More documentation can be found within the file itself.

There are three locations where the file may be stored: 1) in the biotoolbox 
root directory, 2) in the user's home root directory, or 3) a custom 
location defined by the environment variable 'BIOTOOLBOX'.

*Execution 
All biotoolbox scripts are designed to be run from the command line or
executed from another script. Some programs, for example
manipulate_datasets.pl, also provide an interactive interface to allow for
spontaneous work or when the exact index number or name of the dataset in
the file or database is not immediately known.

*Help 
All scripts require command line options for execution. Executing the
program without any options will present a synopsis of the options that are
available. Most programs also have a --help option, which will display
detailed information about the program and execution (usually by displaying
the internal POD). The options are given in the long format (--help, for
example), but may be shortened to single letters if the first letter is
unique (-h, for example).

*File Formats 
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




UPDATES AND FIXES

The biotoolbox scripts are under continual development by the author. The 
latest distribution may be downloaded from the project website at 
http://code.google.com/p/biotoolbox/. 

The latest bugfixes and updates may also be obtained through SVN.

Please contact the author for bugs. Feature requests are also accepted, 
within time constraints. Contact information is at the project website.




ONLINE DOCUMENTATION

Setting up a computer for BioToolBox
http://code.google.com/p/biotoolbox/wiki/BioToolBoxSetUp

Setting up Mac OS X for bioinformatics
http://code.google.com/p/biotoolbox/wiki/SetupForMacOSX

Up to data list of BioToolBox programs
http://code.google.com/p/biotoolbox/wiki/ProgramList

Setting up a Bio::DB::SeqFeature::Store database
http://code.google.com/p/biotoolbox/wiki/LoadingDatabase

Description of supported data file formats
http://code.google.com/p/biotoolbox/wiki/DataFileFormats

Example Data collection
http://code.google.com/p/biotoolbox/wiki/ExampleFeatureDataCollection

Example working with Next Generation Sequencing Data
http://code.google.com/p/biotoolbox/wiki/ExampleNextGenerationSequencing

Explanation of the text data report format
http://code.google.com/p/biotoolbox/wiki/TimDataFormat

Mapping SNPs 
http://code.google.com/p/biotoolbox/wiki/MappingSNPs


