package tim_file_helper;

### modules
require Exporter;
use strict;
use Carp;
use File::Basename qw(fileparse);
use IO::File;
use IO::Zlib;
use Storable qw(store_fd fd_retrieve store retrieve);
use Statistics::Lite qw(mean min);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	generate_tim_data_structure
	verify_data_structure
	find_column_index
);


### Variables
# Export
our @ISA = qw(Exporter);
our @EXPORT = qw(
);
our @EXPORT_OK = qw(
	load_tim_data_file 
	open_tim_data_file 
	write_tim_data_file 
	open_to_read_fh
	open_to_write_fh
	convert_genome_data_2_gff_data 
	convert_and_write_to_gff_file
	write_summary_data
);

# Count for unique GFF3 ID creation
our $GFF3_ID_COUNT = 0; 

# List of acceptable filename extensions
	# include gzipped versions, but list uncompressed versions first
	# need to escape the periods so that they match periods and not 
	# any character - apparently fileparse() uses a regex
our @SUFFIX_LIST = qw(
	\.store
	\.store\.gz
	\.txt
	\.txt\.gz
	\.gff
	\.gff\.gz
	\.gtf
	\.gtf\.gz
	\.gff3
	\.gff3\.gz
	\.bed
	\.bed\.gz
	\.sgr
	\.sgr\.gz
	\.kgg
); 


### The True Statement
1; 



### Load new version data table from file

sub load_tim_data_file {
	# this subroutine will load a tim data file entirely into memory
	
	# retrieve file name as passed argument and check it
	my $filename = shift; 
	$filename = _check_file($filename);
	unless ($filename) {
		carp " file '$filename' does not exist!\n";
		return;
	}
	
	# split filename into its base components
	my ($basename, $path, $extension) = fileparse($filename, @SUFFIX_LIST);
	
	# The resulting data structure
	my $inputdata_ref;
	
	# Open the file based on the file type and load datastructure
	if ($extension =~ /store/i) {
		# the file is a binary Storable.pm file
		
		if ($extension eq '.store.gz') {
			# the stored file is compressed
			# need to filter through gunzip first
			# we're opening a standard globtype filehandle for compatibility
			# with Storable.pm
			
			open *FILE, "gunzip -c $filename |" or croak "unable to open gunzip -c $filename";
			$inputdata_ref = fd_retrieve(\*FILE);
			close FILE;
		}
		else {
			# an uncompressed stored file
			
			$inputdata_ref = retrieve($filename);
		}
		
		# check file
		unless (verify_data_structure($inputdata_ref) ) {
			carp "badly formatted data file!";
			return;
		}
		
		# record current filename data
		$inputdata_ref->{'filename'} = $filename; # the original filename
		$inputdata_ref->{'basename'} = $basename; # filename basename
		$inputdata_ref->{'extension'} = $extension; # the filename extension
	}
	
	#
	else {
		# a text file that must be parsed into a data structure
		
		# open the file and parse the metadata
		# this will return the metadata hash and an open filehandle
		my $fh; # open filehandle for processing the data table
		($fh, $inputdata_ref) = open_tim_data_file($filename);
		unless ($fh) {
			return;
		}
		
		# prepare the data table
		my @datatable;
		
		# put the column names into the data table
		push @datatable, $inputdata_ref->{'column_names'};
		delete $inputdata_ref->{'column_names'}; # we no longer need this
		
		# load the data table
		while (my $line = $fh->getline) {		
			chomp $line;
			
			# the current file position should be at the beginning of the
			# data table information
			
			# simply read each line in the file, explode the line into an 
			# anonymous array and push it to the data_table array
			my @linedata = split /\t/, $line;
			
			# convert null values to internal '.'
			for (my $i = 0; $i < $inputdata_ref->{'number_columns'}; $i++ ) {
				if (!defined $linedata[$i]) {
					# a null value to convert
					$linedata[$i] = '.';
				}
				if ($linedata[$i] =~ /^n\/?a$/i) {
					# value matches na or n/a, a null value
					$linedata[$i] = '.';
				}
			}
			
			# check the number of elements
			if (scalar @linedata != $inputdata_ref->{'number_columns'} ) {
				carp "line $. has wrong number of columns! " .
					scalar(@linedata) . " instead of " . 
					$inputdata_ref->{'number_columns'};
			}
			
			# store the array
			push @datatable, [ @linedata ];
		}
		
		# associate the data table with the data hash
		$inputdata_ref->{'data_table'} = \@datatable;
		
		# record the index number of the last data row
		$inputdata_ref->{'last_row'} = scalar @datatable - 1;
		
		# completed loading the file
		$fh->close;
		
	}
	
	# finished
	return $inputdata_ref;
}







#### Open a tim data file, process the metadata, and return the open filehandle

sub open_tim_data_file {
	# This subroutine will open a tim data text file, process the headers and
	# extract the metadata, load the metadata into a data hash, then return
	# the metadata and the open filehandle for further processing.
	
	# This may be useful in avoiding perl malloc errors when trying to load 
	# a humungous data file into memory. In that case, the file may processed 
	# line by line.
	
	# get filename
	my $file = shift;
	my $filename = _check_file($file);
	unless ($filename) {
		carp " file '$file' does not exist!\n";
		return;
	}
	
	# split filename into its base components
	my ($basename, $path, $extension) = fileparse($filename, @SUFFIX_LIST);
	
	# check for Storable files
	if ($extension =~ /store/i) {
		carp " stored files should be opened with tim_file_helper::load_tim_data_file";
		return;
	}
	
	# open the file
	my $fh = open_to_read_fh($filename);
	unless (defined $fh) {
		return;
	}
	
	
	# generate the data hash to store the file metadata into
	my $inputdata = generate_tim_data_structure(undef);
	unless ($inputdata) {
		carp " cannot generate tim data structure!\n";
		return;
	}
	$inputdata->{'filename'} = $filename; # the original filename
	$inputdata->{'basename'} = $basename; # filename basename
	$inputdata->{'extension'} = $extension; # the filename extension
	$inputdata->{'path'} = $path; # the original path
	$inputdata->{'column_names'} = []; # array for the column names
	$inputdata->{'headers'} = undef; # boolean to indicate column headers present
	
	# read and parse the file
	# we will ONLY parse the header lines prefixed with a #, as well as the 
	# first data row which contains the column names, except for gff files
	# which don't have column names
	my $header_line_count = 0;
	PARSE_HEADER_LOOP:
	while (my $line = $fh->getline) {		
		# we are not chomping the line here because of possible side effects
		# with UCSC tables where we have to count elements in the first line
		# and potential header lines, and the first line has a null value at 
		# the end
		
		# Parse the datafile metadata headers
		
		# the generating program
		if ($line =~ m/^# Program (.+)$/) {
			$inputdata->{'program'} = $1;
			chomp $inputdata->{'program'}; # in case it comes through the grep
			$header_line_count++;
		}
		
		# the source database
		elsif ($line =~ m/^# Database (.+)$/) {
			$inputdata->{'db'} = $1;
			chomp $inputdata->{'db'};
			$header_line_count++;
		}
		
		# the type of feature in this datafile
		elsif ($line =~ m/^# Feature (.+)$/) {
			$inputdata->{'feature'} = $1;
			chomp $inputdata->{'feature'};
			$header_line_count++;
		}
		
		# column or dataset specific information
		elsif ($line =~ m/^# Column/) {
			my %temphash; # a temporary hash to put the column metadata into
			my $index; # to record the column index
			
			# strip the Column metadata identifier
			my $metadataline = $line; # to avoid manipulating $line
			chomp $metadataline;
			$metadataline =~ s/^# Column_\d+ //; 
			# this follows the syntax from write_tim_data_file()
			# the digits represent the column index + 1
			
			# break up the column metadata
			foreach (split /;/, $metadataline) {
				my ($key, $value) = split /=/;
				if ($key eq 'index') {
					# record the index value
					$index = $value;
				}
				# store the key & value
				$temphash{$key} = $value;
			}
			
			# store the column metadata hash into the main data hash
			# use the index as the key
			if (defined $index) {
				# index was defined in the file's metadata
				# check for pre-existing main metadata hash for this column
				# for example, gff files will have standard hashes created above
				# or possibly a badly formatted previous metadata line
				if (exists $inputdata->{$index}) {
					# we will simply overwrite the previous metadata hash
					# harsh, I know, but what to do?
					# if it was canned metadata for a gff file, that's ok
					$inputdata->{$index} = \%temphash;
				}
				else {
					# metadata hash doesn't exist, so we will add it
					$inputdata->{$index} = \%temphash;
					# also update the number of columns
					$inputdata->{'number_columns'} += 1;
				}
			} 
			else {
				# index was not defined in the file's metadata
				# a poorly formed metadata line, complain to user
				carp "badly formatted column metadata line at line $.!";
				
				# attempt a recovery anyway
				if (%temphash) {
					# we have loaded something into the temporary hash from
					# the file's metadata line
					
					# assign a new index number
					my $new_index = $inputdata->{'number_columns'};
					$temphash{'index'} = $new_index;
					
					# check for pre-existing main metadata hash
					if (exists $inputdata->{$new_index}) {
						# we will simply overwrite the previous metadata hash
						$inputdata->{$new_index} = \%temphash;
					}
					else {
						# metadata hash doesn't exist, so add it
						$inputdata->{$new_index} = \%temphash;
						$inputdata->{'number_columns'} += 1;
					}
				}
			}
			$header_line_count++;
		}
		
		# gff version header
		elsif ($line =~ /^##gff-version\s+([\d\.]+)$/) {
			# store the gff version in the hash
			# this may or may not be present in the gff file, but want to keep
			# it if it is
			$inputdata->{'gff'} = $1;
			chomp $inputdata->{'gff'};
			$header_line_count++;
		}
		
		# any other nonstandard header
		elsif ($line =~ /^#/) {
			# store in an anonymous array in the inputdata hash
			push @{ $inputdata->{'other'} }, $line;
			$header_line_count++;
		}
		
		# the remainder is the data table itself
		else {
			# the first row in the data table are (usually) the column names 
			# we only want the names, not the rest of the table
			
			# specific file formats have implicit predefined column formats
			# these file formats do NOT have column headers
			# we will first check for those file formats and process accordingly
			
			# working with a gff file
			if ($extension =~ /g[tf]f/i) {
				# gff files have nine defined columns
				# there are different specifications and variants:
				# gff (v.1), gff v.2, gff v.2.5 (aka gtf), gff v.3 (gff3)
				# however, the columns are the same in all versions
				# more info on gff can be found here
				# http://gmod.org/wiki/GFF3
				
				# set values
				$inputdata->{'number_columns'} = 9; # supposed to be 9 columns
				if ($inputdata->{'gff'} == 0) { 
					# try and determine type based on extension (yeah, right)
					if ($extension =~ /gtf/) {
						$inputdata->{'gff'} = 2.5;
					}
					elsif ($extension =~ /gff3/) {
						$inputdata->{'gff'} = 3;
					}
					else {
						# likely version 2, never really encounter v. 1
						$inputdata->{'gff'} = 2;
					}
				}
				
				# a list of gff column names
				my @gff_names = qw(
					Chromosome
					Source
					Type
					Start
					Stop
					Score
					Strand
					Phase
					Group
				);
				
				# set the metadata for the each column
					# some of these may already be defined if there was a 
					# column metadata specific column in the file
				for (my $i = 0; $i < 9; $i++) {
					# loop for each column
					# set name unless it already has one from metadata
					$inputdata->{$i}{'name'} = $gff_names[$i] unless 
						exists $inputdata->{$i}{'name'};
					$inputdata->{$i}{'index'} = $i;
					# assign the name to the column header
					$inputdata->{'column_names'}->[$i] = 
						$inputdata->{$i}{'name'};
				}
				
				# set headers flag to false
				$inputdata->{'headers'} = 0;
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
			
			# working with a bed file
			elsif ($extension =~ /bed/i) {
				# bed files have a loose format
				# they require a minimum of 3 columns, and have a max of 12
				# 3, 6, and 12 column files are the most common
				# there are also something called paired bed files floating around
				# The official details and specifications may be found at 
				# http://genome.ucsc.edu/FAQ/FAQformat#format1
				
				# first determine the number of columns we're working
				my @elements = split /\s/, $line;
					# normally tab-delimited, but the specs are not explicit
				my $column_count = scalar @elements;
				$inputdata->{'number_columns'} = $column_count; 
				
				# Define the columns and metadata
				if ($column_count >= 3) {
					# the first three columns are required: chrom start end
					
					# set column names
					my @bed_names = qw(
						Chromosome
						Start
						End
					);
					
					# set the metadata for each column
						# some of these may already be defined if there was a 
						# column metadata specific column in the file
					for (my $i = 0; $i < 3; $i++) {
						# loop for each column
						# set name unless it already has one from metadata
						$inputdata->{$i}{'name'} = $bed_names[$i] unless 
							exists $inputdata->{$i}{'name'};
						$inputdata->{$i}{'index'} = $i;
						# assign the name to the column header
						$inputdata->{'column_names'}->[$i] = 
							$inputdata->{$i}{'name'};
					}
				}
				
				if ($column_count >= 4) {
					# name of the bed line feature
					
					# column metadata
					$inputdata->{3}{'name'} = 'Name' unless 
						exists $inputdata->{3}{'name'};
					$inputdata->{3}{'index'} = 3;
					
					# column header name
					$inputdata->{'column_names'}->[3] = 
						$inputdata->{3}{'name'};
				}	
				
				if ($column_count >= 5) {
					# score of the bed line feature
					
					# column metadata
					$inputdata->{4}{'name'} = 'Score' unless 
						exists $inputdata->{4}{'name'};
					$inputdata->{4}{'index'} = 4;
					
					# column header name
					$inputdata->{'column_names'}->[4] = 
						$inputdata->{4}{'name'};
				}	
				
				if ($column_count >= 6) {
					# strand of the bed line feature
					
					# column metadata
					$inputdata->{5}{'name'} = 'Strand' unless 
						exists $inputdata->{5}{'name'};
					$inputdata->{5}{'index'} = 5;
					
					# column header name
					$inputdata->{'column_names'}->[5] = 
						$inputdata->{5}{'name'};
				}	
				
				if ($column_count >= 7) {
					# start position for block (exon)
					
					# column metadata
					$inputdata->{6}{'name'} = 'thickStart' unless 
						exists $inputdata->{6}{'name'};
					$inputdata->{6}{'index'} = 6;
					
					# column header name
					$inputdata->{'column_names'}->[6] = 
						$inputdata->{6}{'name'};
				}	
				
				if ($column_count >= 8) {
					# end position for block (exon)
					
					# column metadata
					$inputdata->{7}{'name'} = 'thickEnd' unless 
						exists $inputdata->{7}{'name'};
					$inputdata->{7}{'index'} = 7;
					
					# column header name
					$inputdata->{'column_names'}->[7] = 
						$inputdata->{7}{'name'};
				}	
				
				if ($column_count >= 9) {
					# RGB value of bed feature
					
					# column metadata
					$inputdata->{8}{'name'} = 'itemRGB' unless 
						exists $inputdata->{8}{'name'};
					$inputdata->{8}{'index'} = 8;
					
					# column header name
					$inputdata->{'column_names'}->[8] = 
						$inputdata->{8}{'name'};
				}	
				
				if ($column_count >= 10) {
					# The number of blocks (exons)
					
					# column metadata
					$inputdata->{9}{'name'} = 'blockCount' unless 
						exists $inputdata->{9}{'name'};
					$inputdata->{9}{'index'} = 9;
					
					# column header name
					$inputdata->{'column_names'}->[9] = 
						$inputdata->{9}{'name'};
				}	
				
				if ($column_count >= 11) {
					# The size of the blocks (exons)
					
					# column metadata
					$inputdata->{10}{'name'} = 'blockSizes' unless 
						exists $inputdata->{10}{'name'};
					$inputdata->{10}{'index'} = 10;
					
					# column header name
					$inputdata->{'column_names'}->[10] = 
						$inputdata->{10}{'name'};
				}	
				
				if ($column_count >= 12) {
					# The start positions of the blocks (exons)
					
					# column metadata
					$inputdata->{11}{'name'} = 'blockStarts' unless 
						exists $inputdata->{11}{'name'};
					$inputdata->{11}{'index'} = 11;
					
					# column header name
					$inputdata->{'column_names'}->[11] = 
						$inputdata->{11}{'name'};
				}	
				
				if ($column_count > 12) {
					# why would there be extra columns in here!!??
					
					warn " BED file '$filename' has too many columns! Bad formatting?\n";
					
					# process anyway
					for (my $i = 11; $i < $column_count; $i++) {
						# column metadata
						$inputdata->{$i}{'name'} = "Column_$i" unless 
							exists $inputdata->{$i}{'name'};
						$inputdata->{$i}{'index'} = $i;
						
						# column header name
						$inputdata->{'column_names'}->[$i] = 
							$inputdata->{$i}{'name'};
					}
				}
				
				else {
					die " BED file '$filename' doesn't have at least 3 columns!\n";
				}
				
				# set headers flag to false
				$inputdata->{'headers'} = 0;
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
			
			# a SGR file
			elsif ($extension =~ /sgr/i) {
				# a sgr file contains three columns: chromo, position, score
				# this is a very simple file format, useful in exporting and
				# importing to binary BAR files used in T2, USeq, and IGB
				
				# set values
				$inputdata->{'number_columns'} = 3; # supposed to be 3 columns
				
				# set column names
				push @{ $inputdata->{'column_names'} }, qw(
					Chromo
					Start
					Score
				);
				
				# set the metadata for the each column
					# some of these may already be defined if there was a 
					# column metadata specific column in the file
				$inputdata->{0} = {
					'name'  => 'Chromo',
					'index' => 0,
				} unless exists $inputdata->{0};
				$inputdata->{1} = {
					'name'  => 'Start',
					'index' => 1,
				} unless exists $inputdata->{1};
				$inputdata->{2} = {
					'name'  => 'Score',
					'index' => 2,
				} unless exists $inputdata->{2};
				
				# set headers flag to false
				$inputdata->{'headers'} = 0;
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
			
			# Data table files from UCSC Table Browser
				# these will have one comment line marked with #
				# that really contains the column headers
			elsif (
				# !exists $inputdata->{0} and
				scalar @{ $inputdata->{'other'} } >= 1 and
				scalar(split /\t/, $inputdata->{'other'}->[-1]) == 
					scalar(split /\t/, $line)
			) {
				# lots of requirements, but this checks that (column 
				# metadata has not already been loaded), there is at least one 
				# unknown comment line, and that the last comment line is 
				# splitable and has the same number of elements as the first 
				# data line
				
				# process the real header line
				my $header_line = pop @{ $inputdata->{'other'} };
				chomp $header_line;
				
				# generate the metadata
				my $i = 0;
				foreach (split /\t/, $header_line) {
					$inputdata->{$i} = { 
						'name'   => $_,
						'index'  => $i,
					};
					push @{ $inputdata->{'column_names'} }, $_;
					$i++;
				}
				$inputdata->{'number_columns'} = $i;
				
				# we do not count the current line as a header
				
				# set headers flag to true
				$inputdata->{'headers'} = 1;
				
				last PARSE_HEADER_LOOP;
			}
			
			# all other file formats, including tim data files
			else {
				# we have not yet parsed the row of data column names
				# we will do so now
				
				my @namelist = split /\t/, $line;
				chomp $namelist[-1];
				
				# we will define the columns based on
				for my $i (0..$#namelist) {
					# confirm that a file metadata exists for this column
					if (exists $inputdata->{$i}) {
						unless ($namelist[$i] eq $inputdata->{$i}->{'name'}) {
							carp "metadata and header names for column $i do not match!";
							# set the name to match the actual column name
							$inputdata->{$i}->{'name'} = $namelist[$i];
						}
					} 
					
					# otherwise be nice and generate it here
					else {
						$inputdata->{$i} = {
							'name'  => $namelist[$i],
							'index' => $i,
						};
					}
				}
				
				# check the number of columns
				if (scalar @namelist != $inputdata->{'number_columns'} ) {
					# adjust to match actual content
					$inputdata->{'number_columns'} = scalar @namelist;
				}
				
				# put the column names in the metadata
				push @{ $inputdata->{'column_names'} }, @namelist;
				
				# count as a header line
				$header_line_count++;
				
				# set headers flag to true
				$inputdata->{'headers'} = 1;
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
		}
		
	}
	
	
	# close and re-open the file
		# I tried using seek functions - but they don't work with binary gzip 
		# files, and I can't get the seek function to return the same position
		# as simply advancing through the file like below
		# so I'll just do it the old way and close/open and advance
	$fh->close;
	$fh = open_to_read_fh($filename);

	# advance to the first data line position
	for (1 .. $header_line_count) {
		my $line = $fh->getline;
	}
	
	# now return the advanced filehandle and the parsed metadata
	return ($fh, $inputdata);
	
}






### Write out a data file from the data hash


sub write_tim_data_file {
	
	# collect passed arguments
	my $argument_ref = shift;
	unless ($argument_ref) {
		carp "no arguments passed!";
		return;
	}
	my $datahash_ref = $argument_ref->{'data'} || undef;
	my $filename = $argument_ref->{'filename'} || undef;
	my $gz = $argument_ref->{'gz'} || undef;
	my $format = $argument_ref->{'format'} || undef;
	
	unless (defined $datahash_ref) {
		# we need data to write
		carp "no data to write!\n";
		return;
	}
	unless (verify_data_structure($datahash_ref) ) {
		carp "bad data structure!";
		return;
	}
	
	# determine filename
	unless ($filename) {
		# we need a filename to write
		
		# check for a filename in the metadata
		if (exists $datahash_ref->{'filename'}) {
			# re-use the original file name 
			$filename = $datahash_ref->{'filename'};
		}
		else {
			# complain about no file name
			carp "no filename given!\n";
			return;
		}
	}
	
	# split filename into its base components
	my ($name, $path, $extension) = fileparse($filename, @SUFFIX_LIST);
	
	# Adjust filename extension if necessary
	unless ($extension) {
			
		# a binary store file
		if (defined $format and $format eq 'store') {
			$extension = '.store';
		}
		
		# a gff file
		elsif (
			$datahash_ref->{'gff'} and 
			$datahash_ref->{'number_columns'} == 9
		) {
			$extension = '.gff';
		} 
		
		# original file had an extension, re-use it
		elsif (exists $datahash_ref->{'extension'}) {
			$extension = $datahash_ref->{'extension'};
		}
		
		# normal data text file
		else {
			$extension = '.txt';
		}
	}
	
	# determine format 
	unless ($format) {
		if (defined $argument_ref->{'simple'}) {
			# an old method of specifying simple
			$format = 'simple';
		}
		elsif ($extension) {
			# check extension from the parsed filename, if present
			if ($extension =~ /store/i) {
				# filename suggests the Storable format
				$format = 'store';
			}
			elsif ($extension =~ /sgr/) {
				# sgr is simple format, no headers 
				$format = 'simple';
			}
			elsif ($extension =~ /bed/) {
				# bed is simple format, no headers 
				$format = 'simple';
			}
			elsif ($extension =~ /g[tf]f/) {
				# gff has no headers, but allows comment lines
				# this has a special exception later when writing text
				unless (
					exists $datahash_ref->{'gff'} and 
					$datahash_ref->{'gff'} > 0
				) {
					# set this to true
					$datahash_ref->{'gff'} = 2; # default version?
				}
				$format = 'text'; #
			}
			else {
				# everything else is text
				$format = 'text';
			}
		}
		else {
			# somehow we got this far without defining? use default text
			$format = 'text';
		}
	}
	
	# check zip status if necessary
	unless (defined $gz) {
		# look at filename extension as a clue
		# in case we're overwriting the input file, keep the zip status
		if ($extension =~ m/\.gz$/i) {
			$gz = 1;
		}
		else {
			$gz = 0; # default
		}
	}
	
	# adjust gzip extension as necessary
	if ($gz) {
		# add gz extension if necessary
		unless ($extension =~ m/\.gz$/i) {
			$extension .= '.gz';
		}
	}
	else {
		$extension =~ s/\.gz$//i; # strip gz extension if present
	}
	
	# generate the new filename
	my $newname = $path . $name . $extension;
	
	
	# Write file based on the format
	if ($format eq 'store') {
		# Storable binary data format using architecture independent format
		
		if ($gz) {
			# compress the stored file
			# we're opening the file using a standard globtype file handle
			# to make it compatible with Storable.pm
			
			open *FILE, "| gzip >$newname " or croak "unable to open gzip $newname";
			store_fd($datahash_ref, \*FILE);
			close FILE;
		}
		else {
			# write an uncompressed store file
			
			store($datahash_ref, $newname);
		}
		
		
	}
	elsif ($format eq 'text' or $format eq 'simple') {
		# A text format
	
		# Open file for writing
		my $fh = open_to_write_fh($newname, $gz);
		unless (defined $fh) { 
			return;
		}
		
		
		# Write the headers
		if ($format eq 'text') {
			# default text format has metadata headers
			# 'simple format
			
			# write gff statement if gff format
			if ($datahash_ref->{'gff'}) {
				# write gff statement
				print {$fh} "##gff-version $datahash_ref->{gff}\n";
			}
			
			# Write the primary headers
			if ($datahash_ref->{'program'}) {
				# write program header if present
				print {$fh} '# Program ' . $datahash_ref->{'program'} . "\n";
			}
			if ($datahash_ref->{'db'}) {
				# write database header if present
				print {$fh} '# Database ' . $datahash_ref->{'db'} . "\n";
			}
			if ($datahash_ref->{'feature'}) {
				# write feature header if present
				print {$fh} '# Feature ' . $datahash_ref->{'feature'} . "\n";
			}
			
			# Write the miscellaneous headers
			foreach ( @{ $datahash_ref->{'other'} } ) {
				# write remaining miscellaneous header lines if present
				unless (/\n$/s) {
					# append newline if not present
					$_ .= "\n";
				}
				# check for comment character at beginning
				if (/^# /) {
					print {$fh} $_;
				}
				else {
					print {$fh} "# " . $_;
				}
			}
		
			# Write the column metadata headers
			for (my $i = 0; $i < $datahash_ref->{'number_columns'}; $i++) {
				# each column metadata in the hash is referenced by the column's
				# index number as the key
				# we will take each index one at a time in increasing order
				
				# some files do not need or tolerate metadata lines, for those 
				# known files the metadata lines will be skipped
				
				# these column metadata lines do not need to be written if they
				# only have two values, presumably name and index, for files 
				# that don't normally have column headers, e.g. gff
				if ($datahash_ref->{'extension'} =~ /sgr|kgg/i) {
					# these do not need metadata
					next;
				}
				elsif (
					scalar( keys %{ $datahash_ref->{$i} } ) == 2 and
					$datahash_ref->{'extension'} =~ /gtf|gff|bed/i
				) {
					# only two keys are present and it's a gff or bed file
					# these are standard metadata keys (name, index) and do 
					# not need to be written
					# more than two keys suggests an interesting metadata key
					# so it will be written
					next;
				}
				
				# we will put each key=value pair into @pairs, listed asciibetically
				my @pairs; # an array of the key value pairs from the metadata hash
				# put name and index first
				push @pairs, 'name=' . $datahash_ref->{$i}{'name'};
				push @pairs, 'index=' . $datahash_ref->{$i}{'index'};
				# put remainder in alphabetical order
				foreach (sort {$a cmp $b} keys %{ $datahash_ref->{$i} } ) {
					next if $_ eq 'name';
					next if $_ eq 'index';
					push @pairs,  $_ . '=' . $datahash_ref->{$i}{$_};
				}
				
				# Finally write the header line, joining the pairs with a semi-colon
				# into a single string.
				# The column identifier is comprised of the word 'Column' and the
				# the index number + 1 (shift to 1-based counting) joined by '_'.
				print {$fh} "# Column_", $i+1, " ", join(";", @pairs), "\n";
			}
		}
		
		
		# Write the table column headers
			# we'll check the boolean value of whether headers existed
		if (exists $datahash_ref->{'headers'}) {
			if ($datahash_ref->{'headers'}) {
				# table headers existed in original source file, 
				# and they should be written again
				print {$fh} 
					join("\t", @{ $datahash_ref->{'data_table'}[0] }), "\n";
			}
		}
		else {
			# doesn't look like the datahash came from a loaded file
			# print headers dependent on file type
			unless ($extension =~ /g[tf]f|sgr|bed/) {
				# these special files do not use headers
				# everything else writes the header
				print {$fh} 
					join("\t", @{ $datahash_ref->{'data_table'}[0] }), "\n";
			}
		}
			
		
		# Write the data table
		if ($format eq 'simple') {
			
			# the simple format will strip the non-value '.' from the table
			for (my $i = 1; $i <= $datahash_ref->{'last_row'}; $i++) {
				# we will step though the data_table array one row at a time
				# convert the non-value '.' to undefined
				# and print using a tab-delimited format
				my @linedata;
				foreach ( @{ $datahash_ref->{'data_table'}[$i] }) {
					if ($_ eq '.') {
						push @linedata, q{}; # an undefined value
					} else {
						push @linedata, $_;
					}
				}
				print {$fh} join("\t", @linedata) . "\n";
			}
		}
		
		else {
			# normal data files
			for (my $i = 1; $i <= $datahash_ref->{'last_row'}; $i++) {
				# we will step though the data_table array one row at a time
				# we will join each row's array of elements into a string to print
				# using a tab-delimited format
				print {$fh} 
					join("\t", @{ $datahash_ref->{'data_table'}[$i] }), "\n";
			}
		}
		
		# done writing
		$fh->close;
	}
	
	else {
		# some unrecognized file format was passed on
		carp "unrecognized file format!\n";
		return;
	}
	
	# if we made it this far, it should've been a success!
	# return the new file name as indication of success
	return $newname;
}






#### Open a file for reading

sub open_to_read_fh {
	# a simple subroutine to open a filehandle for reading a file
	
	# check filename
	my $file = shift;
	my $filename = _check_file($file);
	unless ($filename) {
		carp " file '$file' does not exist!\n";
		return;
	}
	
	
	# check for binary store file
	if ($filename =~ /\.store(?:\.gz)?$/) {
		carp " binary store files should not be opened with this subroutine!\n";
		return;
	}
	

	# Open filehandle object as appropriate
	my $fh; # filehandle
	if ($filename =~ /\.gz$/i) {
		# the file is compressed with gzip
		$fh = IO::Zlib->new;
	} 
	else {
		# the file is uncompressed and space hogging
		$fh = IO::File->new;
	}
	
	# Open file and return
	if ($fh->open($filename, "r") ) {
		return $fh;
	}
	else {
		carp "unable to open file '$filename': " . $fh->error . "\n";
		return;
	}
}






#### Open a file for writing

sub open_to_write_fh {
	# A simple subroutine to open a filehandle for writing a file
	
	my ($filename, $gz, $append) = @_;
	
	# check filename
	unless ($filename) {
		carp " no filename to write!";
		return;
	}
	
	# check zip status if necessary
	unless (defined $gz) {
		# look at filename extension as a clue
		# in case we're overwriting the input file, keep the zip status
		if ($filename =~ m/\.gz$/i) {
			$gz = 1;
		}
		else {
			$gz = 0; # default
		}
	}
	
	# check file append mode
	unless (defined $append) {
		# default is not to append
		$append = 0;
	}
	
	# determine write mode
	my $mode;
	if ($gz and $append) {
		# append a gzip file
		$mode = 'ab';
	}
	elsif (!$gz and $append) {
		# append a normal file
		$mode = 'a';
	}
	elsif ($gz and !$append) {
		# write a new gzip file
		$mode = 'wb';
	}
	else {
		# write a new normal file
		$mode = 'w';
	}
	
	
	# Generate appropriate filehandle object and name
	my $fh;
	if ($gz) {
		# write a space-saving compressed file
		
		# add gz extension if necessary
		unless ($filename =~ m/\.gz$/i) {
			$filename .= '.gz';
		}
		
		$fh = new IO::Zlib;
	}
	else {
		# write a normal space-hogging file
		
		# strip gz extension if present
		$filename =~ s/\.gz$//i; 
		
		$fh = new IO::File;
	}
	
	# Open file for writing and return
	if ($fh->open($filename, $mode) ) {
		return $fh;
	}
	else {
		carp " unable to open file '$filename': " . $fh->error . "\n";
		return;
	}
}





#### Convert a tim data table into GFF format

sub convert_genome_data_2_gff_data {
	# a subroutine to convert the data table format from genomic bins or
	# windows to a gff format for writing as a gff file
	
	### Establish general gff variables
	
	# get passed arguments
	my $arg_ref = shift;
	unless ($arg_ref) {
		carp "no arguments passed!";
		return;
	}
	my $input_data_ref = $arg_ref->{'data'};
	unless (verify_data_structure($input_data_ref) ) {
		carp "bad data structure!";
		return;
	}
	
	# chromosome
	my $chr_index;
	if (exists $arg_ref->{'chromo'} and $arg_ref->{'chromo'} =~ /^\d+$/) {
		$chr_index = $arg_ref->{'chromo'};
	}
	else {
		$chr_index = find_column_index($input_data_ref, '^chr|seq|refseq');
	}
		
	# start position
	my $start_index;
	if (exists $arg_ref->{'start'} and $arg_ref->{'start'} =~ /^\d+$/) {
		$start_index = $arg_ref->{'start'};
	}
	else {
		$start_index = find_column_index($input_data_ref, 'start');
	}
		
	# stop position
	my $stop_index;
	if (exists $arg_ref->{'stop'} and $arg_ref->{'stop'} =~ /^\d+$/) {
		$stop_index = $arg_ref->{'stop'};
	}
	else {
		$stop_index = find_column_index($input_data_ref, 'stop|end');
	}
	
	# check that we have coordinates
	unless ( defined $chr_index ) {
		carp " unable to identify chromosome index!";
		return;
	}
	unless ( defined $start_index ) {
		carp " unable to identify start index!";
		return;
	}
	
	# score
	my $score_index;
	if (exists $arg_ref->{'score'} and $arg_ref->{'score'} =~ /^\d+$/) {
		$score_index = $arg_ref->{'score'};
	}
	
	# name
	my $name;
	my $name_index;
	if (exists $arg_ref->{'name'} and $arg_ref->{'name'} ne q()) {
		if ($arg_ref->{'name'} =~ /^\d+$/) {
			# name is a single digit, most likely a index
			$name_index = $arg_ref->{'name'};
		}
		else {
			# name is likely a string
			$name = $arg_ref->{'name'};
		}
	}
	
	# strand
	my $strand_index;
	if ( exists $arg_ref->{'strand'} and $arg_ref->{'strand'} =~ /^\d+$/ ) {
		$strand_index = $arg_ref->{'strand'};
	}
	
	# set gff version, default is 3
	my $gff_version = $arg_ref->{'version'} || 3;
	
	
	# get array of tag indices
	my @tag_indices;
	if (exists $arg_ref->{'tags'}) {
		@tag_indices = @{ $arg_ref->{'tags'} };
	}
	
	# identify the unique ID index
	my $id_index;
	if (exists $arg_ref->{'id'} and $arg_ref->{'id'} =~ /^\d+$/ ) {
		$id_index = $arg_ref->{'id'} ;
	}
	
	# reference to the data table
	my $data_table_ref = $input_data_ref->{'data_table'};
	
	
	### Identify default values
	
	# set default gff data variables
	my $source;
	if (exists $arg_ref->{'source'} and $arg_ref->{'source'} ne q() ) {
		# defined in passed arguments
		$source = $arg_ref->{'source'};
	}
	else {
		# the default is data
		$source = 'data';
	}
	my ($method, $method_index);
	if (exists $arg_ref->{'method'} and $arg_ref->{'method'} ne q() ) {
		# defined in passed arguments
		if ($arg_ref->{'method'} =~ /^\d+$/) {
			# the method looks like a single digit, most likely an index value
			$method_index = $arg_ref->{'method'};
		}
		else {
			# explicit method string
			$method = $arg_ref->{'method'};
		}
	}
	elsif (exists $arg_ref->{'type'} and $arg_ref->{'type'} ne q() ) {
		# defined in passed arguments, alternate name
		if ($arg_ref->{'type'} =~ /^\d+$/) {
			# the method looks like a single digit, most likely an index value
			$method_index = $arg_ref->{'type'};
		}
		else {
			# explicit method string
			$method = $arg_ref->{'type'};
		}
	}
	elsif (defined $name_index) {
		# the name of the dataset for the features' name
		$method = $input_data_ref->{$name_index}{'name'};
	}
	elsif (defined $score_index) {
		# the name of the dataset for the score
		$method = $input_data_ref->{$score_index}{'name'};
	}
	else {
		$method = 'Experiment';
	}
	
	# fix method and source if necessary
	# replace any whitespace or dashes with underscores
	$source =~ s/[\s\-]/_/g;
	$method =~ s/[\s\-]/_/g if defined $method;
	
	
	### Other processing
	# convert the start postion to 1-based from 0-based
	my $convert_zero_base = $arg_ref->{'zero'} || 0;
	
	
	### Reorganize the data table
		# NOTE: this will destroy any information that may be here and not
		# included in the gff data
		# since we're working with referenced data, you better hope that you
		# don't want this data later....
	
	# relabel the data table headers
	$data_table_ref->[0] = [ 
		qw( RefSeq Source Method Start Stop Score Strand Phase Group) 
	];
	
	# re-write the data table
	for my $row (1..$input_data_ref->{'last_row'}) {
		
		# collect coordinate information
		my $refseq = $data_table_ref->[$row][$chr_index];
		my $start = $data_table_ref->[$row][$start_index];
		my $stop;
		if (defined $stop_index) {
			$stop = $data_table_ref->[$row][$stop_index];
		}
		else {
			$stop = $start;
		}
		if ($convert_zero_base) {
			# coordinates are 0-based, shift the start postion
			$start += 1;
		}
		if ($arg_ref->{'midpoint'} and $start != $stop) {
			# if the midpoint is requested, then assign the midpoint to both
			# start and stop
			my $position = sprintf "%.0f", ( ($start + $stop) / 2 );
			$start = $position;
			$stop = $position;
		}
		
		# collect strand information
		my $strand;
		if (defined $strand_index) {
			my $value = $data_table_ref->[$row][$strand_index];
			if ($value =~ m/\A [f \+ 1 w]/xi) {
				# forward, plus, one, watson
				$strand = '+';
			}
			elsif ($value =~ m/\A [r \- c]/xi) {
				# reverse, minus, crick
				$strand = '-';
			}
			elsif ($value =~ m/\A [0 \.]/xi) {
				# zero, period
				$strand = '.';
			}
			else {
				# unidentified, assume it's non-stranded
				$strand = '.';
			}
		}
		else {
			# no strand information
			$strand = '.';
		}
		
		# collect gff method/type name
		my $gff_method;
		if (defined $method_index) {
			# variable method, index was defined
			$gff_method = $data_table_ref->[$row][$method_index];
		}
		else {
			# explicit method
			$gff_method = $method;
		}
		
		# collect score information
		my $score;
		if (defined $score_index) {
			$score = $data_table_ref->[$row][$score_index];
		}
		else {
			$score = '.';
		}
		
		# collect group information
		my $group;
		if ($gff_version == 3) {
			
			# define and record GFF ID
			if (defined $id_index) {
				# this assumes that the $id_index values are all unique
				# user's responsibility to fix it otherwise
				$group = 'ID=' . $data_table_ref->[$row][$id_index] . ';';
			}
			else {
				# ID is not really needed unless we have subfeatures
				# therefore we no longer automatically include ID
# 				if (defined $name) {
# 					# a name string was explicitly defined
# 					# increment the GFF ID count to make unique
# 					$GFF3_ID_COUNT++;
# 					# Record ID
# 					$group = 'ID=' . $name . '.' . $GFF3_ID_COUNT . ';';
# 				}
# 				else {
# 					# use the method as the name
# 					# increment the GFF ID count to make unique
# 					$GFF3_ID_COUNT++;
# 					# Record ID 
# 					$group = 'ID=' . $gff_method . '.' . $GFF3_ID_COUNT . ';';
# 				}
			}
			
			# define and record the GFF Name
			if (defined $name_index) {
				# a name is provided for each feature
				$group .= 'Name=' . $data_table_ref->[$row][$name_index];
			}
			elsif (defined $name) {
				# a name string was explicitly defined
				$group .= "Name=$name";
			}
			else {
				# use the method as the name
				$group .= "Name=$gff_method";
			}
		}
		else { 
			# gff_version 2
			if (defined $name_index) {
				$group = "$gff_method \"" . $data_table_ref->[$row][$name_index] . "\"";
			}
			else {
				$group = "Experiment $gff_method";
			}
		}
		
		# add group tag information if present
		foreach (@tag_indices) {
			unless ($data_table_ref->[$row][$_] eq '.') {
				# no tag if null value
				$group .= ';' . lc($input_data_ref->{$_}{name}) . '=' . 
					$data_table_ref->[$row][$_];
			}
		}
		
		# rewrite in gff format
		$data_table_ref->[$row] = [ (
			$refseq,
			$source, 
			$gff_method,
			$start,
			$stop,
			$score,
			$strand, 
			'.', # phase, this generally isn't used
			$group
		) ];
	}
	
	
	
	### Reorganize metadata
	# there may be some useful metadata in the current data hash that 
	# will be pertinant to the re-generated gff data
	# we need to keep this metadata, toss the rest, and re-write new
	
	# from the tim_db_helper, get_new_genome_list() only really has useful
	# metadata from the start column, index 1
	# also keep any metadata from the score and name columns, if defined
	
	# keep some current metadata
	my $start_metadata_ref = $input_data_ref->{$start_index};
	my $score_metadata_ref; # new empty hashes
	my $group_metadata_ref;
	if (defined $score_index) {
		$score_metadata_ref = $input_data_ref->{$score_index};
	}
	if (defined $name_index) {
		$group_metadata_ref = $input_data_ref->{$name_index};
	}
	
	# delete old metadata
	for (my $i = 0; $i < $input_data_ref->{'number_columns'}; $i++) {
		# delete the existing metadata hashes
		# they will be replaced with new ones
		delete $input_data_ref->{$i};
	}
	
	# define new metadata
	$input_data_ref->{0} = {
		'name'  => 'RefSeq',
		'index' => 0,
	};
	$input_data_ref->{1} = {
		'name'  => 'Source',
		'index' => 1,
	};
	$input_data_ref->{2} = {
		'name'  => 'Method',
		'index' => 2,
	};
	$input_data_ref->{3} = $start_metadata_ref;
	$input_data_ref->{3}{'name'} = 'Start';
	$input_data_ref->{3}{'index'} = 3;
	$input_data_ref->{4} = {
		'name'  => 'Stop',
		'index' => 4,
	};
	$input_data_ref->{5} = $score_metadata_ref;
	$input_data_ref->{5}{'name'} = 'Score';
	$input_data_ref->{5}{'index'} = 5;
	$input_data_ref->{6} = {
		'name'  => 'Strand',
		'index' => 6,
	};
	$input_data_ref->{7} = {
		'name'  => 'Phase',
		'index' => 7,
	};
	$input_data_ref->{8} = $group_metadata_ref;
	$input_data_ref->{8}{'name'} = 'Group';
	$input_data_ref->{8}{'index'} = 8;
	
	# reset the number of columns
	$input_data_ref->{'number_columns'} = 9;
	
	# set the gff metadata to write a gff file
	$input_data_ref->{'gff'} = $gff_version;
	
	# set headers to false
	$input_data_ref->{'headers'} = 0;
	
	# success
	return 1;
}


#### Export a tim data table to GFF file

sub convert_and_write_to_gff_file {
	# a subroutine to export the data table format from genomic bins or
	# windows to a gff file
	
	## Establish general variables
	
	# reset GFF3 counter
	$GFF3_ID_COUNT = 0;
	
	# get passed arguments
	my $arg_ref = shift;
	unless ($arg_ref) {
		carp "no arguments passed!";
		return;
	}
	
	# basics
	my $input_data_ref = $arg_ref->{'data'};
	unless ($input_data_ref) {
		carp "no data structure passed!";
		return;
	}
	unless (verify_data_structure($input_data_ref) ) {
		carp "bad data structure!";
		return;
	}
	# reference to the data table
	my $data_table_ref = $input_data_ref->{'data_table'};
	
	# chromosome
	my $chr_index;
	if (exists $arg_ref->{'chromo'} and $arg_ref->{'chromo'} =~ /^\d+$/) {
		$chr_index = $arg_ref->{'chromo'};
	}
	else {
		$chr_index = find_column_index($input_data_ref, '^chr|seq|refseq');
	}
		
	# start position
	my $start_index;
	if (exists $arg_ref->{'start'} and $arg_ref->{'start'} =~ /^\d+$/) {
		$start_index = $arg_ref->{'start'};
	}
	else {
		$start_index = find_column_index($input_data_ref, 'start');
	}
		
	# stop position
	my $stop_index;
	if (exists $arg_ref->{'stop'} and $arg_ref->{'stop'} =~ /^\d+$/) {
		$stop_index = $arg_ref->{'stop'};
	}
	else {
		$stop_index = find_column_index($input_data_ref, 'stop|end');
	}
	
	# check that we have coordinates
	unless ( defined $chr_index ) {
		carp " unable to identify chromosome index!";
		return;
	}
	unless ( defined $start_index ) {
		carp " unable to identify start index!";
		return;
	}
	
	# score
	my $score_index;
	if (exists $arg_ref->{'score'} and $arg_ref->{'score'} =~ /^\d+$/) {
		$score_index = $arg_ref->{'score'};
	}
	
	# name
	my $name;
	my $name_index;
	if (exists $arg_ref->{'name'} and $arg_ref->{'name'} ne q()) {
		if ($arg_ref->{'name'} =~ /^\d+$/) {
			# name is a single digit, most likely a index
			$name_index = $arg_ref->{'name'};
		}
		else {
			# name is likely a string
			$name = $arg_ref->{'name'};
		}
	}
	
	# strand
	my $strand_index;
	if ( exists $arg_ref->{'strand'} and $arg_ref->{'strand'} =~ /^\d+$/ ) {
		$strand_index = $arg_ref->{'strand'};
	}
	
	# GFF file version, default is 3
	my $gff_version = $arg_ref->{'version'} || 3;
	
	# get array of tag indices
	my @tag_indices;
	if (exists $arg_ref->{'tags'}) {
		@tag_indices = @{ $arg_ref->{'tags'} };
	}
	
	# identify the unique ID index
	my $id_index;
	if (exists $arg_ref->{'id'} and $arg_ref->{'id'} =~ /^\d+$/ ) {
		$id_index = $arg_ref->{'id'} ;
	}
	
	
	
	## Set default gff data variables
	my $source;
	if (exists $arg_ref->{'source'} and $arg_ref->{'source'} ne q() ) {
		# defined in passed arguments
		$source = $arg_ref->{'source'};
	}
	else {
		# the default is data
		$source = 'data';
	}
	
	my ($method, $method_index);
	if (exists $arg_ref->{'method'} and $arg_ref->{'method'} ne q() ) {
		# defined in passed arguments
		if ($arg_ref->{'method'} =~ /^\d+$/) {
			# the method looks like a single digit, most likely an index value
			$method_index = $arg_ref->{'method'};
		}
		else {
			# explicit method string
			$method = $arg_ref->{'method'};
		}
	}
	elsif (exists $arg_ref->{'type'} and $arg_ref->{'type'} ne q() ) {
		# defined in passed arguments, alternate name
		if ($arg_ref->{'type'} =~ /^\d+$/) {
			# the method looks like a single digit, most likely an index value
			$method_index = $arg_ref->{'type'};
		}
		else {
			# explicit method string
			$method = $arg_ref->{'type'};
		}
	}
	elsif (defined $name) {
		# the name provided
		$method = $name;
	}
	elsif (defined $name_index) {
		# the name of the dataset for the features' name
		$method = $input_data_ref->{$name_index}{'name'};
	}
	elsif (defined $score_index) {
		# the name of the dataset for the score
		$method = $input_data_ref->{$score_index}{'name'};
	}
	else {
		$method = 'Experiment';
	}
	# fix method and source if necessary
	# replace any whitespace or dashes with underscores
	$source =~ s/[\s\-]/_/g;
	$method =~ s/[\s\-]/_/g if defined $method;
	
	
	## Open output file
	# get the filename
	my $filename;
	if ( $arg_ref->{'filename'} ne q() ) {
		$filename = $arg_ref->{'filename'};
		# remove unnecessary extensions
		$filename =~ s/\.gz$//;
		$filename =~ s/\.(txt|store)$//;
		unless ($filename =~ /\.gff$/) {
			# add extension if necessary
			$filename .= '.gff';
		}
	}
	elsif (defined $name) {
		# specific name provided
		$filename = $name . '.gff';
	}
	elsif (defined $method and $method ne 'Experiment') {
		# use the method name, so long as it is not the default Experiment
		$filename = $method . '.gff';
	}
	elsif (defined $input_data_ref->{'basename'}) {
		# use the base file name for lack of a better name
		$filename = $input_data_ref->{'basename'} . '.gff';
	}
	else {
		# what, still no name!!!????
		$filename = 'your_stupid_output_gff_file_with_no_name.gff';
	}
	if ($gff_version == 3) {
		$filename .= '3'; # make extension gff3
	}
	
	# open the file for writing 
	my $gz = $arg_ref->{'gz'} || 0;
	my $output_gff = open_to_write_fh($filename, $gz);
	
	# write basic headers
	print {$output_gff} "##gff_version $gff_version\n";
	if (exists $input_data_ref->{'filename'}) {
		# record the original file name for reference
		print {$output_gff} "# Exported from file '", 
			$input_data_ref->{'filename'}, "'\n";
	}
	
	# Write the column metadata headers
		# write the metadata lines only if there is useful information
		# and only for the pertinent columns (chr, start, score)
		# we will check the relavent columns for extra information beyond
		# that of name and index
		# we will write the metadata then for that dataset that is being used
		# substituting the column name and index appropriate for a gff file
		
		# check the chromosome metadata
		if (scalar( keys %{ $input_data_ref->{$chr_index} } ) > 2) {
			# chromosome has extra keys of info
			print {$output_gff} "# Column_1 ";
			my @pairs;
			foreach (sort {$a cmp $b} keys %{ $input_data_ref->{$chr_index} } ) {
				if ($_ eq 'index') {
					push @pairs, "index=0";
				}
				elsif ($_ eq 'name') {
					push @pairs, "name=RefSeq";
				}
				else {
					push @pairs,  $_ . '=' . $input_data_ref->{$chr_index}{$_};
				}
			}
			print {$output_gff} join(";", @pairs), "\n";
		}
		
		# check the start metadata
		if (scalar( keys %{ $input_data_ref->{$start_index} } ) > 2) {
			# start has extra keys of info
			print {$output_gff} "# Column_4 ";
			my @pairs;
			foreach (sort {$a cmp $b} keys %{ $input_data_ref->{$start_index} } ) {
				if ($_ eq 'index') {
					push @pairs, "index=3";
				}
				elsif ($_ eq 'name') {
					push @pairs, "name=Start";
				}
				else {
					push @pairs,  $_ . '=' . $input_data_ref->{$start_index}{$_};
				}
			}
			print {$output_gff} join(";", @pairs), "\n";
		}
		
		# check the score metadata
		if (
			defined $score_index and
			scalar( keys %{ $input_data_ref->{$score_index} } ) > 2
		) {
			# score has extra keys of info
			print {$output_gff} "# Column_6 ";
			my @pairs;
			foreach (sort {$a cmp $b} keys %{ $input_data_ref->{$score_index} } ) {
				if ($_ eq 'index') {
					push @pairs, "index=5";
				}
				elsif ($_ eq 'name') {
					push @pairs, "name=Score";
				}
				else {
					push @pairs,  $_ . '=' . $input_data_ref->{$score_index}{$_};
				}
			}
			print {$output_gff} join(";", @pairs), "\n";
		}
		
		# check the name metadata
		if (
			defined $name_index and
			scalar( keys %{ $input_data_ref->{$name_index} } ) > 2
		) {
			# score has extra keys of info
			print {$output_gff} "# Column_9 ";
			my @pairs;
			foreach (sort {$a cmp $b} keys %{ $input_data_ref->{$name_index} } ) {
				if ($_ eq 'index') {
					push @pairs, "index=8";
				}
				elsif ($_ eq 'name') {
					push @pairs, "name=Group";
				}
				else {
					push @pairs,  $_ . '=' . $input_data_ref->{$name_index}{$_};
				}
			}
			print {$output_gff} join(";", @pairs), "\n";
		}
		
	
			
	# Write the gff features
	for my $row (1..$input_data_ref->{'last_row'}) {
		
		# collect coordinate information
		my $refseq = $data_table_ref->[$row][$chr_index];
		my $start = $data_table_ref->[$row][$start_index];
		my $stop;
		if (defined $stop_index) {
			$stop = $data_table_ref->[$row][$stop_index];
		}
		else {
			$stop = $start;
		}
		if ($arg_ref->{'midpoint'} and $stop != $stop) {
			# if the midpoint is requested, then assign the midpoint to both
			# start and stop
			my $position = sprintf "%.0f", ($start + $stop)/2;
			$start = $position;
			$stop = $position;
		}
		
		# collect score information
		my $score;
		if (defined $score_index) {
			$score = $data_table_ref->[$row][$score_index];
		}
		else {
			$score = '.';
		}
		
		# collect strand information
		my $strand;
		if (defined $strand_index) {
			my $value = $data_table_ref->[$row][$strand_index];
			if ($value =~ m/\A [f \+ 1 w]/xi) {
				# forward, plus, one, watson
				$strand = '+';
			}
			elsif ($value =~ m/\A [r \- c]/xi) {
				# reverse, minus, crick
				$strand = '-';
			}
			elsif ($value =~ m/\A [0 \.]/xi) {
				# zero, period
				$strand = '.';
			}
			else {
				# unidentified, assume it's non-stranded
				$strand = '.';
			}
		}
		else {
			# no strand information
			$strand = '.';
		}
		
		# collect gff method/type name
		my $gff_method;
		if (defined $method_index) {
			# variable method, index was defined
			$gff_method = $data_table_ref->[$row][$method_index];
		}
		else {
			# explicit method
			$gff_method = $method;
		}
		
		
		# collect group information based on version
		my $group;
		if ($gff_version == 3) {
			
			# define and record GFF ID
			if (defined $id_index) {
				# this assumes that the $id_index values are all unique
				# user's responsibility to fix it otherwise
				$group = 'ID=' . $data_table_ref->[$row][$id_index] . ';';
			}
			else {
				# ID is not really needed unless we have subfeatures
				# therefore we no longer automatically include ID
# 				if (defined $name) {
# 					# a name string was explicitly defined
# 					# increment the GFF ID count to make unique
# 					$GFF3_ID_COUNT++;
# 					# Record ID
# 					$group = 'ID=' . $name . '.' . $GFF3_ID_COUNT . ';';
# 				}
# 				else {
# 					# use the method as the name
# 					# increment the GFF ID count to make unique
# 					$GFF3_ID_COUNT++;
# 					# Record ID 
# 					$group = 'ID=' . $gff_method . '.' . $GFF3_ID_COUNT . ';';
# 				}
			}
			
			# define and record the GFF Name
			if (defined $name_index) {
				# a name is provided for each feature
				$group .= 'Name=' . $data_table_ref->[$row][$name_index];
			}
			elsif (defined $name) {
				# a name string was explicitly defined
				$group .= "Name=$name";
			}
			else {
				# use the method as the name
				$group .= "Name=$gff_method";
			}
		}
		else {
			# gff version 2
			if (defined $name) {
				$group = "$gff_method \"$name\"";
			}
			elsif (defined $name_index) {
				$group = "$gff_method \"" . $data_table_ref->[$row][$name_index] . "\"";
			}
			else {
				# really generic
				$group = "Experiment \"$gff_method\"";
			}
		}
		
		# add group tag information if present
		foreach (@tag_indices) {
			unless ($data_table_ref->[$row][$_] eq '.') {
				# no tag if null value
				$group .= ';' . lc($input_data_ref->{$_}{name}) . '=' . 
					$data_table_ref->[$row][$_];
			}
		}
		
		# Write gff feature
		print {$output_gff} join("\t", (
			$refseq,
			$source, 
			$gff_method,
			$start,
			$stop,
			$score,
			$strand, 
			'.', # phase, this generally isn't used
			$group
		) ), "\n";
		
	}
	
	# success
	$output_gff->close;
	return $filename;
}



### Generate a summary data file

sub write_summary_data {
	
	# Collect passed arguments
	my $argument_ref = shift;
	unless ($argument_ref) {
		carp "no arguments passed!";
		return;
	}
	my $datahash_ref =   $argument_ref->{'data'} || undef;
	my $outfile =        $argument_ref->{'filename'} || undef;
	my $dataset =        $argument_ref->{'dataset'} || undef;
	my $startcolumn =    $argument_ref->{'startcolumn'} || undef;
	my $endcolumn =      $argument_ref->{'endcolumn'} ||
	                     $argument_ref->{'stopcolumn'} || undef;
	my $log =            $argument_ref->{'log'} || undef;
	
	
	
	# Check required values
	unless (defined $datahash_ref) {
		carp "no data structure passed!\n";
	}
	unless (verify_data_structure($datahash_ref) ) {
		carp "bad data structure!";
		return;
	}
	unless (defined $outfile) {
		if (defined $datahash_ref->{'basename'}) {
			# use the opened file's basename if it exists
			$outfile = $datahash_ref->{'basename'};
		}
		else {
			carp "no filename passed to write_summary_data!\n";
			return;
		}
	}
	
	
	
	# Auto calculate missing arguments
	unless (defined $startcolumn) {
		# we will attempt to determine where to start the summing process
		# we will skip those columns with common feature-description names
		
		my @acceptable_indices; # array of acceptable indices
		my %skip = (
			# these are example names generated from subs in tim_db_helper.pm
			'SystematicName'  => 1,
			'Name'            => 1,
			'Alias'           => 1,
			'Aliases'         => 1,
			'Type'            => 1,
			'Class'           => 1,
			'GeneClass'       => 1,
			'Chromosome'      => 1,
			'Start'           => 1,
			'Stop'            => 1,
			'Gene'            => 1,
		);
		
		# walk through the dataset names
		for (my $i = 0; $i < $datahash_ref->{'number_columns'}; $i++) {
			unless (exists $skip{ $datahash_ref->{$i}{'name'} } ) {
				# if the dataset name is not in the hash of skip names
				push @acceptable_indices, $i;
			}
		}
		
		# The start column should be the lowest index number in the list of
		# acceptable_indices array.
		# Assuming, of course, that feature descriptor datasets (columns) are
		# leftmost only.
		$startcolumn = min(@acceptable_indices);
	}
	unless (defined $endcolumn) {
		# take the last or rightmost column
		$endcolumn = $datahash_ref->{'number_columns'} - 1;
	}
	unless ($dataset) {
		# the original dataset name (i.e. the name of the dataset in the 
		# database from which the column's data was derived) should be the same 
		# in all the columns
		$dataset = $datahash_ref->{$startcolumn}{'dataset'};
	}
	unless (defined $log) {
		# the log flag should be set in the column metadata and should be the
		# same in all
		$log = $datahash_ref->{$startcolumn}{'log2'};
	}
		
	
	
	# Prepare array to store the summed data
	my $summed_feature = $datahash_ref->{'feature'} . '_averaged_windows';
	my $summed_data = generate_tim_data_structure(
		$summed_feature, 
		qw(
			Window
			Midpoint
			Score
		)
	);
	$summed_data->{'db'} = $datahash_ref->{'database'};
	$summed_data->{0}{'number_features'} = $datahash_ref->{'last_row'};
	$summed_data->{2}{'log2'} = $log;
	$summed_data->{2}{'dataset'} = $dataset;
	
	
	# Collect summarized data
	# we will walk through the columns and take an average for each 
	# column or window
	for (
		my $column = $startcolumn;
		$column <= $endcolumn;
		$column++
	) { 
		
		# determine the midpoint position of the window
		my $midpoint = mean(
			# this assumes the column metadata has start and stop
			$datahash_ref->{$column}{'start'},	
			$datahash_ref->{$column}{'stop'},	
		) or undef; 
		
		
		# collect the values in the column
		my @values;
		for my $row (1..$datahash_ref->{'last_row'}) {
			my $value = $datahash_ref->{'data_table'}->[$row][$column];
			unless ($value eq '.') { 
				push @values, $value;
			}
		}
		
		# adjust if log value
		if ($log) {
			@values = map { 2 ** $_ } @values;
		}
		
		# determine mean value
		my $window_mean = mean(@values);
		if ($log) { 
			$window_mean = log($window_mean) / log(2);
		}
		
		# push to summed output
		push @{ $summed_data->{'data_table'} }, [ (
			$datahash_ref->{$column}{'name'}, 
			$midpoint, 
			$window_mean
		) ];
		$summed_data->{'last_row'} += 1;
	}
	
	# Write summed data
	$outfile =~ s/\.txt(\.gz)?$//; # strip any .txt or .gz extensions if present
	my $written_file = write_tim_data_file( {
		'data'      => $summed_data,
		'filename'  => $outfile . '_summed',
	} );
	
	# Return
	if ($written_file) {
		return $written_file;
	}
	else {
		return;
	}
}





### Internal subroutine to check for file existance
sub _check_file {
	my $filename = shift;
	
	# check for file existance
	if (-e $filename) {
		# confirmed full filename and path
		return $filename;
	}
	else {
		# file name is either incomplete or non-existent
		
		# first try adding some common file extensions in case those are missing
		my $new_filename;
		foreach my $ext (@SUFFIX_LIST) {
			$ext =~ s/\\//g; # the list is designed as a regular expression list
							# with the periods escaped for File::Basename
							# but we need to un-escape them for this purpose
			if (-e $filename . $ext) {
				$new_filename = $filename . $ext;
				last;
			}
		}
		
		# return
		if (defined $new_filename) {
			# found it
			return $new_filename;
		}
		else {
			# we've now tried all the extension combinations
			# and failed to find it
			# return nothing
			return;
		}
	}

}





__END__

=head1 NAME

tim_file_helper

=head1 DESCRIPTION

These are general file helper subroutines to work with data text files, 
primarily opening, loading, and writing. Specifically, it is designed 
to work with I<tim data text files>, which is a generic tab delimited
text file of rows and columns along with special metadata column headers.
While it best uses this I<tim data format>, it will really read any tab 
delimited text file. Special file formats used in bioinformatics, including 
GFF and BED files, are automatically recognized by their file extension and 
appropriate metadata added. 

Files opened using these subroutines are stored in a specific complex data 
structure described below. This format allows for data access as well as 
records metadata about each column (dataset) and the file in general. This
metadata helps preserve a "history" of the dataset: where it came from, how
it was collected, and how it was processed.

Additional subroutines are also present for general processing and output of
this data structure.

The I<tim data file format> is described below, and following that a 
description of the data structure.

=head1 FORMAT OF TIM DATA TEXT FILE

The tim data file format is not indicated by a special file extension. 
Rather, a generic '.txt' extension is used to preserve functionality with
other text processing programs. The file is essentially a simple tab 
delimited text file representing rows (lines) and columns (demarcated by the
tabs). 

What makes it unique are the metadata header lines, each prefixed by a '# '.
These metadata lines describe the data within the table with regards to its
type, source, methodology, history, and processing. The metadata is designed
to be read by both human and computer. Opening files without this metadata 
will result in basic default metadata assigned to each column. Special files
recognized by their extension (e.g. GFF or BED) will have appropriate 
metadata assigned.

The specific metadata lines that are specifically recognized are listed below.

=over 4

=item Feature

The Feature describes the types of features represented on each row in the 
data table. These can include gene, transcript, genome, etc.

=item Database

The name of the database used in generation of the feature table. This 
is often also the database used in collecting the data, unless the dataset
metadata specifies otherwise.

=item Program

The name of the program generating the data table and file. It usually 
includes the whole path of the executable.

=item Column

The next header lines include column specific metadata. Each column 
will have a separate header line, specified initially by the word 
'Column', followed by an underscore and the column number (1-based). 
Following this is a series of 'key=value' pairs separated by ';'. 
Spaces are generally not allowed. Obviously '=' or ';' are not 
allowed or they will interfere with the parsing. The metadata 
describes how and where the data was collected. Additionally, any 
modifications performed on the data are also recorded here. The only 
two keys that are required are the first two, 'name' and 'index'. If
the file being read does not contain metadata, then it will be generated
with these two basic keys.

=back

A list of standard column header keys is below, but is not exhaustive. 

=over

=item name

The name of the column. This should be identical to the table header.

=item index

The 0-based index number of the column.

=item database

Included if different from the main database indicated above.

=item window

The size of the window for genome datasets

=item step

The step size of the window for genome datasets

=item dataset

The name of the dataset(s) from which data is collected. Comma delimited.

=item start

The starting point for the feature in collecting values

=item stop

The stopping point of the feature in collecting values

=item extend

The extension of the region in collecting values

=item strand

The strandedness of the data collecte. Values include 'sense',
'antisense', or 'none'

=item method

The method of collecting values

=item log2

boolean indicating the values are in log2 space or not

=back


Finally, the data table follows the metadata. The table consists of 
tab-delimited data. The same number of fields should be present in each 
row. Each row represents a genomic feature or landmark, and each column 
contains either identifying information or a collected dataset. 
The first row will always contain the column names, except in special
cases such as the GFF format where the columns are strictly defined.
The column name should be the same as defined in the column's metadata.
When loading GFF files, the header names and metadata are automatically
generated for conveniance. 


=head1 BINARY REPRESENTATION OF TIM DATA FILE

Alternatively, the internal memory data structure may be written and 
read directly as a binary file using the Storable.pm module, negating 
the need for both metadata and data parsing. This may significantly 
decrease both load and write times owing to the faster compiled C code. 
This is evident particularly with very large datasets, at the cost of 
slightly larger files. 

Binary stored files are recognized by the file extension '.store'. They 
are stored in an architecture dependent manner. Be careful when migrating 
between different machine architectures, and when in doubt, use the text 
format. As with text files, the files may be compressed using gzip, 
and recognized as such with the '.gz' extension. 

Many of the programs which utilize this module may not directly offer the 
choice of writing a binary file. By explicitly including the '.store' 
extension in the output filename, the write modules within should 
then write a binary file.


=head1 USAGE

Call the module at the beginning of your perl script and it will import 
both read and write modules.
  
  use tim_db_helper;
  

The specific usage for each subroutine is detailed below.


=over


=item load_tim_data_file()

This is a newer, updated file loader and parser for tim's data files. It will
completely parse and load the file contents into the described data structure 
in memory. Files with metadata lines (described in tim data format) will 
have the metadata lines loaded. Files without metadata lines will have basic 
metadata (column name and index) automatically generated. The first 
non-header line should contain the column (dataset) name. Recognized file 
formats without headers, including GFF, BED, and SGR, will have the columns 
automatically named.

This subroutine uses the open_tim_data_file() subroutine and completes the 
loading of the file into memory.

Pass the module the filename. The file may be compressed with gzip, recognized
by the .gz extension.

The subroutine will return a scalar reference to the hash, described above. 
Failure to read or parse the file will return an empty value.


Example:
	
	my $filename = 'my_data.txt.gz';
	my $data_ref = load_tim_data_file($filename);
	
=item open_tim_data_file()

This is a file opener and metadata parser for data files, including tim's 
data formatted files and other recognized data formats (gff, bed, sgr). It 
will open the file, parse the metadata, and return an open file handle 
ready for reading. It will NOT load the entire file contents into memory. 
This is to allow for processing those gigantic data files that will break 
Perl with malloc errors. 

The subroutine will open the file, parse the header lines (marked with
a # prefix) into a metadata hash as described above, parse the data column 
names (the first row in the table), set the file pointer to the first row of
data in the table, and return the open file handle along with a scalar 
reference to the metadata hash. The calling program may then process the file 
through the filehandle line by line as appropriate.

The data column names may be found in an array in the data hash under the 
key 'column_names';

Pass the module the filename. The file may be compressed with gzip, recognized
by the .gz extension.

The subroutine will return two items: a scalar reference to the file handle,
and a scalar reference to the data hash, described as above. The file handle
is an IO::Handle object and may be manipulated as such.
Failure to read or parse the file will return an empty value.

Example:
	
	my $filename = 'my_data.txt.gz';
	my ($fh, $metadata_ref) = open_tim_data_file($filename);
	while (my $line = $fh->getline) {
		...
	}
	$fh->close;



=item write_tim_data_file()

This subroutine will write out a data file formatted for tim's data files. 
Either text or binary files may be written. Please refer to L<FORMAT OF 
TIM DATA TEXT FILE> and L<BINARY REPRESENTATION OF TIM DATA FILE> for more 
information regarding the file format. If the 'gff' key is true in the data 
hash, then a gff file will be written.

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  data     => A scalar reference to the tim data structure ad described
              in B<tim_data_helper.pm>. 
  Optional: 
  filename => A scalar value containing the name of the file to 
              write. This value is required for new data files and 
              optional for overwriting existing files (the filename 
              stored in the metadata is used). Appropriate extensions 
              are added (e.g, .store, .txt, .gz, etc) as neccessary. 
  format   => A string to indicate the file format to be written.
              Acceptable values include 'text', 'store', and 'simple'.
              Text files are text in nature, include all metadata, and
              usually have '.txt' extensions. Store files are binary
              Storable.pm files representing all data and metadata,
              and have '.store' extensions. Simple files are
              tab-delimited text files without metadata, useful for
              exporting data. If the format is not specified, the
              extension of the passed filename will be used as a
              guide. The default behavior is to write standard text
              files.
  gz       => A boolean value (1 or 0) indicating whether the file 
              should be written through a gzip filter to compress. If 
              this value is undefined, then the file name is checked 
              for the presence of the '.gz' extension and the value 
              set appropriately. Default is false.
  simple   => A boolean value (1 or 0) indicating whether a simple 
              tab-delimited text data file should be written. This is 
              an old alias for setting 'format' to 'simple'.

The subroutine will return true if the write was successful, otherwise it will
return undef. The true value is the name of the file written, including any 
changes to the extension if necessary. 

Note that by explicitly providing the filename extension, some of these 
options may be set without providing the arguments to the subroutine. 
Specifically, using '.store' vs '.txt' will set the file format, and 
adding '.gz' will compress the file. The arguments always take 
precendence over the filename extensions, however.

Example

	my $filename = 'my_data.txt.gz';
	my $data_ref = load_tim_data_file($filename);
	...
	my $success_write = write_tim_data_file( {
		'data'     => $data_ref,
		'filename' => $filename,
		'format'   => 'simple',
	} );
	if ($success_write) {
		print "wrote $success_write!";
	}


=item open_to_read_fh()

This subroutine will open a file for reading. If the passed filename has
a '.gz' extension, it will appropriately open the file through a gunzip 
filter.

Pass the subroutine the filename. It will return a scalar reference to the
open filehandle. The filehandle is an IO::Handle object and may be manipulated
as such.

Example
	
	my $filename = 'my_data.txt.gz';
	my $fh = open_to_read_fh($filename);
	while (my $line = $fh->getline) {
		# do something
	}
	$fh->close;
	

=item open_to_write_fh()

This subroutine will open a file for writing. If the passed filename has
a '.gz' extension, it will appropriately open the file through a gzip 
filter.

Pass the subroutine three values: the filename, a boolean value indicating
whether the file should be compressed with gzip, and a boolean value 
indicating that the file should be appended. The gzip and append values are
optional. The compression status may be determined automatically by the 
presence or absence of the passed filename extension; the default is no 
compression. The default is also to write a new file and not to append.

If gzip compression is requested, but the filename does not have a '.gz' 
extension, it will be automatically added. However, the change in file name 
is not passed back to the originating program; beware!

The subroutine will return a scalar reference to the open filehandle. The 
filehandle is an IO::Handle object and may be manipulated as such.

Example
	
	my $filename = 'my_data.txt.gz';
	my $gz = 1; # compress output file with gzip
	my $fh = open_to_write_fh($filename, $gz);
	# write to new compressed file
	print {$fh} "something interesting\n";
	$fh->close;
	

=item convert_genome_data_2_gff_data()

This subroutine will convert an existing data hash structure as described above
and convert it to a defined gff data structure, i.e. one that has the nine 
defined columns. Once converted, a gff data file may then be written using the
write_tim_data_file() subroutine. To convert and write the gff file in one 
step, see the following subroutine, convert_and_write_gff_file();

NOTE: This method is DESTRUCTIVE!!!!
Since the data table will be completely reorganized, any extraneous data 
in the data table will be discarded. Since referenced data is being 
used, any data loss may be significant and unexpected. A normal data file
should be written first to preserve extraneous data, and the conversion to
gff data be the last operation done.

Since the gff data structure requires genomic coordinates, this data must be 
present as identifiable datasets in the data table and metadata. It looks 
specifically for datasets labeled 'Chromosome', 'Start', and 'Stop' or 'End'. 
Failure to identify these datasets will simply return nothing. A dataset 
generated with get_new_genome_list() in 'tim_db_helper.pm' will generate these
datasets. 

The subroutine must be passed a reference to an anonymous hash with the 
arguments. The keys include

  Required:
  data     => A scalar reference to the data hash. The data hash 
              should be as described in this module.
  Optional: 
  chromo   => The index of the column in the data table that contains
              the chromosome or reference name. By default it 
              searches for the first column with a name that begins 
              with 'chr' or 'refseq' or 'seq'.
  start    => The index of the column with the start position. By 
              default it searches for the first column with a name 
              that contains 'start'.
  stop     => The index of the column with the stop position. By 
              default it searches for the first column with a name 
              that contains 'stop' or 'end'.
  score    => The index of the dataset in the data table to be used 
              as the score column in the gff data.
  name     => The name to be used for the GFF features. Pass either 
              the index of the dataset in the data table that 
              contains the unique name for each gff feature, or a 
              text string to be used as the name for all of the 
              features. This information will be used in the 
              'group' column.
  strand   => The index of the dataset in the data table to be used
              for strand information. Accepted values might include
              any of the following 'f(orward), r(everse), w(atson),
              c(rick), +, -, 1, -1, 0, .).
  source   => A scalar value to be used as the text in the 'source' 
              column. Default is 'data'.
  type     => A scalar value to be used as the text in the 'method'
              or 'type' column. If not defined, it will use the 
              name of the dataset used for either the 'score' or 
              'name' column, if defined. As a last resort, it 
              will use the most creative method of 'Experiment'.
              Alternatively, if a single integer is passed, then 
              it is assumed to be the column index for a unique 
              method value for each feature.
  method   => Alias for "type".
  midpoint => A boolean (1 or 0) value to indicate whether the 
              midpoint between the actual 'start' and 'stop' values
              should be used instead of the actual values. Default 
              is false.
  zero     => The coordinates are 0-based (interbase). Convert to 
              1-based format (bioperl conventions).
  tags     => Provide an anonymous array of indices to be added as 
              tags in the Group field of the GFF feature. The tag's 
              key will be the column's name. As many tags may be 
              added as desired.
  id       => Provide the index of the column containing unique 
              values which will be used in generating the GFF ID 
              in v.3 GFF files. If not provided, the ID is 
              automatically generated from the name.
  version  => The GFF version (2 or 3) to be written. The default is 
              version 3.

The subroutine will return true if the conversion was successful, otherwise it
will return nothing.

Example

	my $data_ref = load_tim_data_file($filename);
	...
	my $success = convert_genome_data_2_gff_data( {
		'data'     => $data_ref,
		'score'    => 3,
		'midpoint' => 1,
	} );
	if ($success) {
		# write a gff file
		my $success_write = write_tim_data_file( {
			'data'     => $data_ref,
			'filename' => $filename,
		} );
		if ($success_write) {
			print "wrote $success_write!";
		}
	}


=item convert_and_write_to_gff_file()

This subroutine will convert a tim data structure as described above into 
GFF format and write the file. It will preserve the current data structure 
and convert the data on the fly as the file is written, unlike the 
destructive subroutine convert_genome_data_2_gff_data(). 

Either a v.2 or v.3 GFF file may be written. The only metadata written 
is the original data's filename (if present) and any dataset (column) 
metadata that contains more than the basics (name and index).

Since the gff data structure requires genomic coordinates, this data must be 
present as identifiable datasets in the data table and metadata. It looks 
specifically for datasets labeled 'Chromosome', 'Start', and optionally 
'Stop' or 'End'. Failure to identify these datasets will simply return 
nothing. A dataset generated with get_new_genome_list() in 
'tim_db_helper.pm' will generate these coordinate datasets. 

If successful, the subroutine will return the name of the output gff file
written.

The subroutine must be passed a reference to an anonymous hash with the 
arguments. The keys include

  Required:
  data     => A scalar reference to the data hash. The data hash 
              should be as described in this module.
  Optional: 
  filename => The name of the output GFF file. If not specified, 
              the default value is, in order, the method, name of 
              the indicated 'name' dataset, name of the indicated 
              'score' dataset, or the originating file basename.
  version  => The version of GFF file to write. Acceptable values 
              include '2' or '3'. For v.3 GFF files, unique ID 
              values will be auto generated, unless provided with a 
              'name' dataset index. Default is to write v.3 files.
  chromo   => The index of the column in the data table that contains
              the chromosome or reference name. By default it 
              searches for the first column with a name that begins 
              with 'chr' or 'refseq' or 'seq'.
  start    => The index of the column with the start position. By 
              default it searches for the first column with a name 
              that contains 'start'.
  stop     => The index of the column with the stop position. By 
              default it searches for the first column with a name 
              that contains 'stop' or 'end'.
  score    => The index of the dataset in the data table to be used 
              as the score column in the gff data.
  name     => The name to be used for the GFF features. Pass either 
              the index of the dataset in the data table that 
              contains the unique name for each gff feature, or a 
              text string to be used as the name for all of the 
              features. This information will be used in the 
              'group' column.
  strand   => The index of the dataset in the data table to be used
              for strand information. Accepted values might include
              any of the following 'f(orward), r(everse), w(atson),
              c(rick), +, -, 1, -1, 0, .).
  source   => A scalar value to be used as the text in the 'source' 
              column. Default is 'data'.
  method   => A scalar value to be used as the text in the "method"
              or "type" column. If not defined, it will use the 
              name of the dataset used for either the 'score' or 
              'name' column, if defined. As a last resort, it 
              will use the most creative method of 'Experiment'.
              Alternatively, if a single integer is passed, then 
              it is assumed to be the column index for a unique 
              method value for each feature.
  type     => Alias for "method".
  midpoint => A boolean (1 or 0) value to indicate whether the 
              midpoint between the actual 'start' and 'stop' values
              should be used instead of the actual values. Default 
              is false.
  tags     => Provide an anonymous array of indices to be added as 
              tags in the Group field of the GFF feature. The tag's 
              key will be the column's name. As many tags may be 
              added as desired.
  id       => Provide the index of the column containing unique 
              values which will be used in generating the GFF ID 
              in v.3 GFF files. If not provided, the ID is 
              automatically generated from the name.

Example

	my $data_ref = load_tim_data_file($filename);
	...
	my $success = convert_and_write_to_gff_file( {
		'data'     => $data_ref,
		'score'    => 3,
		'midpoint' => 1,
		'filename' => "$filename.gff",
		'version'  => 2,
	} );
	if ($success) {
		print "wrote file '$success'!";
	}
	

=item write_summary_data()

This subroutine will summarize the data in a data file, generating mean values
for all the values in each dataset (column), and writing an output file with
the summarized data. This is useful for data collected in windows across a 
feature, for example, microarray data values across the body of genes, and 
then generating a composite or average gene occupancy.

The output file is a tim data tab-delimited file as described above with three
columns: The Name of the window, the Midpoint of the window (calculated as the
mean of the start and stop points for the window), and the mean value. The 
table is essentially rotated 90º from the original table; the averages of each
column dataset becomes rows of data.

Pass the subroutine an anonymous hash of arguments. These include:

  Required:
  data        => A scalar reference to the data hash. The data hash 
                 should be as described in this module.
  filename    => The base filename for the file. This will be 
                 appended with '_summed' to differentiate from the 
                 original data file. This may be automatically  
                 obtained from the metadata of an opened file if 
                 not specified, otherwise it will not work.
  Optional: 
  startcolumn => The index of the beginning dataset containing the 
                 data to summarized. This may be automatically 
                 calculated by taking the leftmost column without
                 a known feature-description name (using examples 
                 from tim_db_helper.pm).
  stopcolumn  => The index of the last dataset containing the 
                 data to summarized. This may be automatically 
                 calculated by taking the rightmost column. 
  dataset     => The name of the original dataset used in 
                 collecting the data. It may be obtained from the 
                 metadata for the startcolumn.
  log         => The data is in log2 space. It may be obtained 
                 from the metadata for the startcolumn.

Example

	my $main_data_ref = load_tim_data_file($filename);
	...
	my $summary_success = write_summary_data( {
		'data'         => $main_data_ref,
		'filename'     => $outfile,
		'startcolumn'  => 4,
	} );


=back

=head1 INTERNAL SUBROUTINES

These are internally used subroutines and are not exported for general usage.

=over

=item _check_file

This subroutine confirms the existance of a passed filename. If not 
immediately found, it will attempt to append common file extensions 
and verifiy its existence. This allows the user to pass only the base 
file name and not worry about missing the extension. This may be useful 
in shell scripts.

=back



=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  

