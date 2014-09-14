package Bio::ToolBox::file_helper;

### modules
require Exporter;
use strict;
use Carp qw(carp cluck croak confess);
use File::Basename qw(fileparse);
use IO::File;
use Statistics::Lite qw(mean min);
use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure
	verify_data_structure
	find_column_index
);
our $VERSION = '1.20';


### Variables
# Export
our @ISA = qw(Exporter);
our @EXPORT = qw(
);
our @EXPORT_OK = qw(
	load_tim_data_file 
	parse_filename
	open_tim_data_file 
	write_tim_data_file 
	open_to_read_fh
	open_to_write_fh
	convert_genome_data_2_gff_data 
	convert_and_write_to_gff_file
	write_summary_data
);

# List of acceptable filename extensions
	# include gzipped versions, but list uncompressed versions first
	# need to escape the periods so that they match periods and not 
	# any character - apparently fileparse() uses a regex
our @SUFFIX_LIST = qw(
	\.txt
	\.txt\.gz
	\.txt\.bz2
	\.gff
	\.gff\.gz
	\.gtf
	\.gtf\.gz
	\.gff3
	\.gff3\.gz
	\.bed
	\.bed\.gz
	\.bdg
	\.bdg\.gz
	\.bedgraph
	\.bedgraph.gz
	\.sgr
	\.sgr\.gz
	\.kgg
	\.cdt
	\.vcf
	\.vcf\.gz
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
	
	# open the file and parse the metadata
	# this will return the metadata hash and an open filehandle
	my ($fh, $inputdata) = open_tim_data_file($filename);
	unless ($fh) {
		return;
	}
	
	# prepare the data table
	my @datatable;
	
	# put the column names into the data table
	push @datatable, $inputdata->{'column_names'};
	delete $inputdata->{'column_names'}; # we no longer need this
	
	# determine strand column 
		# should be true for BED and GFF files or any other file with strand
	my $strand_i = find_column_index($inputdata, '^strand$');
	my $plusminus_count = 0;
	
	# load the data table
	while (my $line = $fh->getline) {		
		
		# the current file position should be at the beginning of the
		# data table information
		
		# skip comment lines
		if ($line =~ /^#/) {
			push @{ $inputdata->{'other'} }, $line;
			next;
		}
		
		# no real line, just empty space
		if ($line !~ m/\w+/) {
			next;
		}
		
		# simply read each line in the file, explode the line into an 
		# anonymous array and push it to the data_table array
		my @linedata = split /\t/, $line;
		
		# check the number of elements
		if (scalar @linedata != $inputdata->{'number_columns'} ) {
			if ($line =~ /\r/ and $line !~ /\n/) {
				croak "File '$filename' does not appear to have unix line endings!\n" . 
					" Please convert to unix-style line endings and try again\n";
			}
			else {
				carp "File '$filename' is inconsistent! line $. has ", scalar(@linedata),
					" columns instead of expected ", $inputdata->{'number_columns'}, "\n";
				return;
			}
		}
		
		# chomp the last element
		# we do this here to ensure the tab split above gets all of the values
		# otherwise trailing null values aren't included in @linedata
		# be sure to handle both newlines and carriage returns
		$linedata[-1] =~ s/[\r\n]+$//;
		
		# convert null values to internal '.'
		for (my $i = 0; $i < $inputdata->{'number_columns'}; $i++ ) {
			if (!defined $linedata[$i]) {
				# not defined position in the array?
				$linedata[$i] = '.';
			}
			elsif ($linedata[$i] eq '') {
				# a null value
				$linedata[$i] = '.';
			}
			elsif ($linedata[$i] =~ /^n\/?a$/i) {
				# value matches na or n/a, a null value
				$linedata[$i] = '.';
			}
		}
		
		# convert 0-based starts to 1-base for BED source files
		if ($inputdata->{'bed'}) {
			# add 1 to each start position
			$linedata[1] += 1;
		}
		
		# convert strand information to integers
		if (defined $strand_i) {
			# convert any interpretable value to signed value
			if ($linedata[$strand_i] eq '+') {
				# just a simple plus
				$linedata[$strand_i] = 1;
				$plusminus_count++;
			}
			elsif ($linedata[$strand_i] eq '-') {
				# simple minus
				$linedata[$strand_i] = -1;
				$plusminus_count++;
			}
			elsif ($linedata[$strand_i] eq '.') {
				# unstranded GFF format, not BED
				$linedata[$strand_i] = 0;
				$plusminus_count++;
			}
			# otherwise assume bioperl convention -1, 0, 1
			# if it is not, then hope for the best
			# I am dropping support for the ancient forward, watson, reverse, crick
			# who uses those anyway?????
			# some old bioperl scripts may still support it
		}
		
		# store the array
		push @datatable, [ @linedata ];
	}
	
	# record metadata for base number and strand style
	if ($inputdata->{'bed'}) {
		$inputdata->{1}{'base'} = 1;
	}
	if (defined $strand_i) {
		# add an internal use only metadata value to indicate we've changed strand values
		if ($plusminus_count == (scalar @datatable - 1)) {
			# all the loaded lines are plusminus style
			$inputdata->{$strand_i}{'strand_style'} = 'plusminus';
			if (exists $inputdata->{$strand_i}{'AUTO'}) {
				# update automatically generated metadata
				$inputdata->{$strand_i}{'AUTO'}++;
			}
		}
	}
	
		
	# associate the data table with the data hash
	$inputdata->{'data_table'} = \@datatable;
	
	# record the index number of the last data row
	$inputdata->{'last_row'} = scalar @datatable - 1;
	
	# completed loading the file
	$fh->close;
	
	# verify the structure
	verify_data_structure($inputdata);
	
	# finished
	return $inputdata;
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
	my ($basename, $path, $extension) = parse_filename($filename);
	
	# open the file
	my $fh = open_to_read_fh($filename);
	unless (defined $fh) {
		return;
	}
	
	
	# generate the data hash to store the file metadata into
	my $inputdata = generate_tim_data_structure(undef);
	unless ($inputdata) {
		cluck " cannot generate tim data structure!\n";
		return;
	}
	$inputdata->{'filename'} = $filename; # the original filename
	$inputdata->{'basename'} = $basename; # filename basename
	$inputdata->{'extension'} = $extension; # the filename extension
	$inputdata->{'path'} = $path; # the original path
	$inputdata->{'column_names'} = []; # array for the column names
	$inputdata->{'headers'} = q(); # boolean to indicate column headers present
	$inputdata->{'program'} = q(); # clear program, this will be read from file
	
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
		
		# check for Mac-style return line endings
		if ($line =~ /\r/ and $line !~ /\n/) {
			croak "File '$filename' does not appear to have unix line endings!\n" . 
				" Please convert to unix-style line endings and try again\n";
		}
		
		# Parse the datafile metadata headers
		
		# no real line, just empty space
		if ($line !~ m/\w+/) {
			$header_line_count++;
			next;
		}
		
		# the generating program
		elsif ($line =~ m/^# Program (.+)$/) {
			$inputdata->{'program'} = $1;
			$inputdata->{'program'} =~ s/[\r\n]+$//;
			$header_line_count++;
		}
		
		# the source database
		elsif ($line =~ m/^# Database (.+)$/) {
			$inputdata->{'db'} = $1;
			$inputdata->{'db'} =~ s/[\r\n]+$//;
			$header_line_count++;
		}
		
		# the type of feature in this datafile
		elsif ($line =~ m/^# Feature (.+)$/) {
			$inputdata->{'feature'} = $1;
			$inputdata->{'feature'} =~ s/[\r\n]+$//;
			$header_line_count++;
		}
		
		# column or dataset specific information
		elsif ($line =~ m/^# Column_(\d+)/) {
			# the column number will become the index number
			# it should be 0-based
				# but it may be 1-based
				# if the metadata includes an index key, then this number is 
				# 1-based and the index key value is the real 0-based number
			my $index = -1; # not a real index number
			$index    = $1; # capture this from grep above
			
			# strip the Column metadata identifier
			my $metadataline = $line; # to avoid manipulating $line
			$metadataline =~ s/[\r\n]+$//;
			$metadataline =~ s/^# Column_\d+ //; 
			
			# break up the column metadata
			my %temphash; # a temporary hash to put the column metadata into
			foreach (split /;/, $metadataline) {
				my ($key, $value) = split /=/;
				if ($key eq 'index') {
					if ($index != $value) {
						# the index value obtained from the Column number
						# above is not correct - it was 1-based and we 
						# need 0-based
						# the value from the metadata index key will be 
						# correct, so we will use that
						$index = $value;
					}
				}
				# store the key & value
				$temphash{$key} = $value;
			}
			
			# create a index metadata key if not already present
			# the rest of biotoolbox may expect this to be present
			unless (exists $temphash{'index'}) {
				$temphash{'index'} = $index;
			}
			
			# store the column metadata hash into the main data hash
			# use the index as the key
			if ($index >= 0) {
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
			$inputdata->{'gff'} =~ s/[\r\n]+$//;
			$header_line_count++;
		}
		
		# any other nonstandard header
		elsif ($line =~ /^#/) {
			# store in an anonymous array in the inputdata hash
			push @{ $inputdata->{'other'} }, $line;
			$header_line_count++;
		}
		
		# a track line 
		elsif ($line =~ /^track\s+/i) {
			# common with wig, bed, or bedgraph files for use with
			# the UCSC genome browser
			# treat as a comment line, there's not that much useful info
			# in here, possibly the name, that's it
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
			
			# Data tables with a commented header line
				# including data table files from UCSC Table Browser
				# also VCF files as well
				# these will have one comment line marked with #
				# that really contains the column headers
			if ( _commented_header_line($inputdata, $line) ) {
				# lots of requirements, but this checks that (column 
				# metadata has not already been loaded), there is at least one 
				# unknown comment line, and that the last comment line is 
				# splitable and has the same number of elements as the first 
				# data line
				
				# process the real header line
				my $header_line = pop @{ $inputdata->{'other'} };
				$header_line =~ s/[\r\n]+$//;
				
				# generate the metadata
				my $i = 0;
				foreach (split /\t/, $header_line) {
					$inputdata->{$i} = { 
						'name'   => $_,
						'index'  => $i,
						'AUTO'   => 3,
					} unless exists $inputdata->{$i};
					
					push @{ $inputdata->{'column_names'} }, $_;
					$i++;
				}
				$inputdata->{'number_columns'} = $i;
				
				# we do not count the current line as a header
				
				# set headers flag to true
				$inputdata->{'headers'} = 1;
				
				# check for special formatted files
					# these include GFF and BED, which don't normally have 
					# header lines, even commented ones, but you never know 
					# what crazy users like Sue will feed these programs
					# we set the gff or bed flag
				if ($extension =~ /gtf|gff3?/) {
					# set the gff version based on the extension
					unless ($inputdata->{'gff'}) {
						$inputdata->{'gff'} = 
							$extension =~ /gtf/  ? 2.5 :
							$extension =~ /gff3/ ? 3   :
							2;
					}
					
					# set the feature type
					unless (defined $inputdata->{'feature'}) {
						$inputdata->{'feature'} = 'region';
					}
				}
				elsif ($extension =~ /bed|bdg|bedgraph/) {
					$inputdata->{'bed'} = $inputdata->{'number_columns'};
					
					# set the feature type
					unless (defined $inputdata->{'feature'}) {
						$inputdata->{'feature'} = 'region';
					}
				}
				
				last PARSE_HEADER_LOOP;
			}
			
			# a GFF file
			elsif ($extension =~ /g[tf]f/i) {
				# gff files have nine defined columns
				# there are different specifications and variants:
				# gff (v.1), gff v.2, gff v.2.5 (aka gtf), gff v.3 (gff3)
				# however, the columns are the same in all versions
				# more info on gff can be found here
				# http://gmod.org/wiki/GFF3
				
				# set column number
				$inputdata->{'number_columns'} = 9; # supposed to be 9 columns
				
				# set the gff version based on the extension
				unless ($inputdata->{'gff'}) {
					$inputdata->{'gff'} = 
						$extension =~ /gtf/  ? 2.5 :
						$extension =~ /gff3/ ? 3   :
						2;
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
					# set metadata unless it's already loaded
					unless (exists $inputdata->{$i}) {
						$inputdata->{$i}{'name'}  = $gff_names[$i];
						$inputdata->{$i}{'index'} = $i;
						$inputdata->{$i}{'AUTO'}  = 3;
					}
					# assign the name to the column header
					$inputdata->{'column_names'}->[$i] = 
						$inputdata->{$i}{'name'};
				}
				
				# set headers flag to false
				$inputdata->{'headers'} = 0;
				
				# set the feature type
				unless (defined $inputdata->{'feature'}) {
					$inputdata->{'feature'} = 'region';
				}
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
			
			# a Bed or BedGraph file
			elsif ($extension =~ /bdg|bed/i) {
				# bed files have a loose format
				# they require a minimum of 3 columns, and have a max of 12
				# 3, 6, and 12 column files are the most common
				# there are also something called paired bed files floating around
				# The official details and specifications may be found at 
				# http://genome.ucsc.edu/FAQ/FAQformat#format1
				
				# a special type of bed file is the bedgraph, using 
				# either a bdg, bedgraph, or simply bed extension
				# these only have four columns, no more, no less
				# the fourth column is score, not name
				
				# first check for a commented header line, in case someone like 
				# Sue puts one in (geez)
				
				# first determine the number of columns we're working
				my @elements = split /\s+/, $line;
					# normally tab-delimited, but the specs are not explicit
				my $column_count = scalar @elements;
				$inputdata->{'number_columns'} = $column_count; 
				$inputdata->{'bed'} = $column_count;
				
				
				
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
						unless (exists $inputdata->{$i}) {
							$inputdata->{$i}{'name'}  = $bed_names[$i];
							$inputdata->{$i}{'index'} = $i;
							$inputdata->{$i}{'AUTO'}  = 3;
						}
						# assign the name to the column header
						$inputdata->{'column_names'}->[$i] = 
							$inputdata->{$i}{'name'};
					}
					
					# add additional columns if necessary
					if ($column_count >= 4) {
						# this is a special column
						# it may be either a name or a score
						# determine this by the file extension, not an 
						# exact science
						
						if ($extension =~ /bdg|graph/i) {
							# a bedgraph file
							# this is the score column
							
							# column metadata
							unless (exists $inputdata->{3}) {
								$inputdata->{3}{'name'}  = 'Score';
								$inputdata->{3}{'index'} = 3;
								$inputdata->{3}{'AUTO'}  = 3;
							}
							
							# column header name
							$inputdata->{'column_names'}->[3] = 
								$inputdata->{3}{'name'};
						}
						
						else {
							# a plain old bed file
							# this is the name column
							
							# column metadata
							unless (exists $inputdata->{3}) {
								$inputdata->{3}{'name'}  = 'Name';
								$inputdata->{3}{'index'} = 3;
								$inputdata->{3}{'AUTO'}  = 3;
							}
							
							# column header name
							$inputdata->{'column_names'}->[3] = 
								$inputdata->{3}{'name'};
						}
					}	
					
					if ($column_count >= 5) {
						# score of the bed line feature
						
						# column metadata
						unless (exists $inputdata->{4}) {
							$inputdata->{4}{'name'}  = 'Score';
							$inputdata->{4}{'index'} = 4;
							$inputdata->{4}{'AUTO'}  = 3;
						}
						
						# column header name
						$inputdata->{'column_names'}->[4] = 
							$inputdata->{4}{'name'};
					}	
					
					if ($column_count >= 6) {
						# strand of the bed line feature
						
						# column metadata
						unless (exists $inputdata->{5}) {
							$inputdata->{5}{'name'}  = 'Strand';
							$inputdata->{5}{'index'} = 5;
							$inputdata->{5}{'AUTO'}  = 3;
						}
						
						# column header name
						$inputdata->{'column_names'}->[5] = 
							$inputdata->{5}{'name'};
					}	
					
					if ($column_count >= 7) {
						# start position for block (exon)
						
						# column metadata
						unless (exists $inputdata->{6}) {
							$inputdata->{6}{'name'}  = 'thickStart';
							$inputdata->{6}{'index'} = 6;
							$inputdata->{6}{'AUTO'}  = 3;
						}
						
						# column header name
						$inputdata->{'column_names'}->[6] = 
							$inputdata->{6}{'name'};
					}	
					
					if ($column_count >= 8) {
						# end position for block (exon)
						
						# column metadata
						unless (exists $inputdata->{7}) {
							$inputdata->{7}{'name'}  = 'thickEnd';
							$inputdata->{7}{'index'} = 7;
							$inputdata->{7}{'AUTO'}  = 3;
						}
						
						# column header name
						$inputdata->{'column_names'}->[7] = 
							$inputdata->{7}{'name'};
					}	
					
					if ($column_count >= 9) {
						# RGB value of bed feature
						
						# column metadata
						unless (exists $inputdata->{8}) {
							$inputdata->{8}{'name'}  = 'itemRGB';
							$inputdata->{8}{'index'} = 8;
							$inputdata->{8}{'AUTO'}  = 3;
						}
						
						# column header name
						$inputdata->{'column_names'}->[8] = 
							$inputdata->{8}{'name'};
					}	
					
					if ($column_count >= 10) {
						# The number of blocks (exons)
						
						# column metadata
						unless (exists $inputdata->{9}) {
							$inputdata->{9}{'name'}  = 'blockCount';
							$inputdata->{9}{'index'} = 9;
							$inputdata->{9}{'AUTO'}  = 3;
						}
						
						# column header name
						$inputdata->{'column_names'}->[9] = 
							$inputdata->{9}{'name'};
					}	
					
					if ($column_count >= 11) {
						# The size of the blocks (exons)
						
						# column metadata
						unless (exists $inputdata->{10}) {
							$inputdata->{10}{'name'}  = 'blockSizes';
							$inputdata->{10}{'index'} = 10;
							$inputdata->{10}{'AUTO'}  = 3;
						}
						
						# column header name
						$inputdata->{'column_names'}->[10] = 
							$inputdata->{10}{'name'};
					}	
					
					if ($column_count >= 12) {
						# The start positions of the blocks (exons)
						
						# column metadata
						unless (exists $inputdata->{11}) {
							$inputdata->{11}{'name'}  = 'blockStarts';
							$inputdata->{11}{'index'} = 11;
							$inputdata->{11}{'AUTO'}  = 3;
						}
						
						# column header name
						$inputdata->{'column_names'}->[11] = 
							$inputdata->{11}{'name'};
					}	
					
					if ($column_count > 12) {
						# why would there be extra columns in here!!??
						
						carp " BED file '$filename' has too many columns! Bad formatting?\n";
						
						# process anyway
						for (my $i = 11; $i < $column_count; $i++) {
							# column metadata
							unless (exists $inputdata->{$i}) {
								$inputdata->{$i}{'name'}  = "Column_$i";
								$inputdata->{$i}{'index'} = $i;
								$inputdata->{$i}{'AUTO'}  = 3;
							}
							
							# column header name
							$inputdata->{'column_names'}->[$i] = 
								$inputdata->{$i}{'name'};
						}
					}
				}
				
				# less than 3 columns!???
				else {
					carp " BED file '$filename' doesn't have at least 3 columns!\n";
					return;
				}
				
				# set the feature type
				unless (defined $inputdata->{'feature'}) {
					$inputdata->{'feature'} = 'region';
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
					'AUTO'  => 3,
				} unless exists $inputdata->{0};
				$inputdata->{1} = {
					'name'  => 'Start',
					'index' => 1,
					'AUTO'  => 3,
				} unless exists $inputdata->{1};
				$inputdata->{2} = {
					'name'  => 'Score',
					'index' => 2,
					'AUTO'  => 3,
				} unless exists $inputdata->{2};
				
				# set headers flag to false
				$inputdata->{'headers'} = 0;
				
				# set the feature type
				unless (defined $inputdata->{'feature'}) {
					$inputdata->{'feature'} = 'region';
				}
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
			
			# all other file formats, including tim data files
			else {
				# we have not yet parsed the row of data column names
				# we will do so now
				
				my @namelist = split /\t/, $line;
				$namelist[-1] =~ s/[\r\n]+$//;
				
				# we will define the columns based on
				for my $i (0..$#namelist) {
					# confirm that a file metadata exists for this column
					if (exists $inputdata->{$i}) {
						unless ($namelist[$i] eq $inputdata->{$i}->{'name'}) {
							cluck "metadata and header names for column $i do not match!";
							# set the name to match the actual column name
							$inputdata->{$i}->{'name'} = $namelist[$i];
						}
					} 
					
					# otherwise be nice and generate it here
					else {
						$inputdata->{$i} = {
							'name'  => $namelist[$i],
							'index' => $i,
							'AUTO'  => 3,
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



### Parse the filename using the list suffix list

sub parse_filename {
	my $filename = shift;
	my ($basename, $path, $extension) = fileparse($filename, @SUFFIX_LIST);
	return ($basename, $path, $extension);
}






### Write out a data file from the data hash


sub write_tim_data_file {
	
	# collect passed arguments
	my %args = @_; 
	$args{'data'}     ||= undef;
	$args{'filename'} ||= undef;
	$args{'format'}   ||= undef;
	unless (exists $args{'gz'}) {$args{'gz'} = undef} 
		# this is a boolean value, need to be cognizant of 0
		# this will be checked below
	
	# check the data
	my $data = $args{'data'};
	unless ($data) {
		# we need data to write
		cluck "no data to write!\n";
		return;
	}
	unless (verify_data_structure($data) ) {
		cluck "bad data structure!";
		return;
	}
	
	# determine filename
	unless ($args{'filename'}) {
		# we need a filename to write
		
		# check for a filename in the metadata
		if (exists $data->{'filename'}) {
			# re-use the original file name 
			$args{'filename'} = $data->{'filename'};
		}
		else {
			# complain about no file name
			cluck "no filename given!\n";
			return;
		}
	}
	
	# split filename into its base components
	my ($name, $path, $extension) = fileparse($args{'filename'}, @SUFFIX_LIST);
	
	# Adjust filename extension if necessary
	if ($extension) {
		# we have an extension
		# make sure it makes sense
		
		# GFF file
		if ($extension =~ /gff/i) {
			unless ($data->{'gff'}) {
				# it's not set as a gff data
				# let's set it to true and see if it passes verification
				$data->{'gff'} = 3; # default
				if ( 
					verify_data_structure($data) and
					$data->{'gff'}
				) {
					# keep the gff extension, it seems  good
				}
				else {
					# it's not good, set it back to text
					warn " re-setting extension from $extension to .txt\n";
					$extension =~ s/gff3?/txt/i;
				}
			}
		}
		
		# BED file
		elsif ($extension =~ /bed|bdg/i) {
			unless ($data->{'bed'}) {
				# it's not set as a bed data
				# let's set it to true and see if it passes verification
				$data->{'bed'} = 1; # a fake true
				if ( 
					verify_data_structure($data) and
					$data->{'bed'}
				) {
					# keep the bed extension, it seems  good
				}
				else {
					# if it's not BED data, we don't use the extension
					# change it to text
					warn " re-setting extension from $extension to .txt\n";
					$extension =~ s/bed|bdg/txt/i;
				}
			}
		}
		
		# SGR file
		elsif ($extension =~ /sgr/i) {
			if (
				exists $data->{'extension'} and
				$data->{'extension'} =~ /sgr/i
			) {
				# if the original file extension was sgr
				# then it likely passed verification above
				# so we will keep it
			}
			else {
				# original file was not SGR
				# let's pretend it was and see if still passes 
				# verification
				# the sgr verification relies on the recorded extension
				$data->{'extension'} = '.sgr';
				
				# re-verify the structure
				verify_data_structure($data);
				
				# we'll take the extension the verification sets
				if ($data->{'extension'} =~ /txt/) {
					warn " re-setting extension from $extension to .txt\n";
				}
				$extension = $data->{'extension'};
			}
		}
		
	}
	
	else {
		# no extension was available
		# try and determine one
				
		# a gff file
		if ($data->{'gff'}) {
			$extension = '.gff';
			
			# for GFF3 files only
			if ($data->{'gff'} == 3) {
				$extension .= '3';
			}
		} 
		
		# a bed file
		elsif ($data->{'bed'}) {
			# check whether it's a real bed or bedgraph file
			if (
				$data->{'number_columns'} == 4 and 
				$data->{3}{'name'} =~ /score/i
			) {
				# a bedgraph file
				$extension = '.bdg';
			}
			else {
				# a regular bed file
				$extension = '.bed';
			}
		}
		
		# original file had an extension, re-use it
		elsif (exists $data->{'extension'}) {
			
			# structure is gff
			if ($data->{'extension'} =~ /gff/i) {
				
				# check to see that we still have a valid gff structure
				if ($data->{'gff'}) {
					$extension = $data->{'extension'};
				}
				else {
					# not valid gff, write a txt file
					$extension = '.txt';
				}
			}
			
			# structure is bed
			elsif ($data->{'extension'} =~ /bed|bdg/i) {
				# check to see that we still have a valid bed structure
				if ($data->{'bed'}) {
					$extension = $data->{'extension'};
				}
				else {
					# not valid bed, write a txt file
					$extension = '.txt';
				}
			}
			
			# unstructured format
			else {
				$extension = $data->{'extension'};
			}
		}
		
		# normal data text file
		else {
			$extension = '.txt';
		}
	}
	
	# determine format 
	unless ($args{'format'}) {
		if (defined $args{'simple'}) {
			# an old method of specifying simple
			$args{'format'} = 'simple';
		}
		elsif ($extension) {
			# check extension from the parsed filename, if present
			if ($extension =~ /sgr|cdt/i) {
				# sgr is simple format, no headers 
				$args{'format'} = 'simple';
			}
			else {
				# everything else is text
				$args{'format'} = 'text';
			}
		}
		else {
			# somehow we got this far without defining? use default text
			$args{'format'} = 'text';
		}
	}
	
	# check zip status if necessary
	unless (defined $args{'gz'}) {
		# look at filename extension as a clue
		# in case we're overwriting the input file, keep the zip status
		if ($extension =~ m/\.gz$/i) {
			$args{'gz'} = 1;
		}
		else {
			$args{'gz'} = 0; # default
		}
	}
	
	# adjust gzip extension as necessary
	if ($args{'gz'} and $extension !~ m/\.gz$/i) {
		$extension .= '.gz';
	}
	elsif (not $args{'gz'} and $extension =~ /\.gz$/i) {
		$extension =~ s/\.gz$//i;
	}
	
	# check filename length
	# assuming a maximum of 256, at least on Mac with HFS+, don't know about Linux
	if (length($name . $extension) > 255) {
		my $limit = 253 - length($extension);
		$name = substr($name, 0, $limit) . '..';
		warn " filename too long! Truncating to $limit characters\n";
	}
	
	# generate the new filename
	my $newname = $path . $name . $extension;
	
	
	# Convert base to interbase coordinates if necessary
	if ($extension =~ /\.bed|bdg/i and $data->{'bed'} > 0) {
		# we are writing a confirmed bed file 
		if (
			exists $data->{1}{'base'} and 
			$data->{1}{'base'} == 1
		) {
			# the start coordinates are in base format
			# need to convert back to interbase
			for (my $row = 1; $row <= $data->{'last_row'}; $row++) {
				# subtract 1 to each start position
				$data->{'data_table'}->[$row][1] -= 1;
			}
			delete $data->{1}{'base'};
		}
	}
	
	
	# Convert strand information
	my $strand_i = find_column_index($data, '^strand$');
	if (
		defined $strand_i and
		exists $data->{$strand_i}{'strand_style'} and 
		$data->{$strand_i}{'strand_style'} eq 'plusminus'
	) {
		# strand information was originally BED and GFF style +,.,-
		# then convert back to that format before writing
		for (my $row = 1; $row <= $data->{'last_row'}; $row++) {
			if ($data->{'data_table'}->[$row][$strand_i] == 1) {
				$data->{'data_table'}->[$row][$strand_i] = '+';
			}
			elsif ($data->{'data_table'}->[$row][$strand_i] == -1) {
				$data->{'data_table'}->[$row][$strand_i] = '-';
			}
			elsif ($data->{'data_table'}->[$row][$strand_i] == 0) {
				$data->{'data_table'}->[$row][$strand_i] = '.';
			}
		}
	}
	
	
	# Open file for writing
	my $fh = open_to_write_fh($newname, $args{'gz'});
	unless (defined $fh) { 
		return;
	}
	
	
	# Write the headers
	if ($args{'format'} eq 'text') {
		# default text format has metadata headers
		# 'simple format
		
		# write gff statement if gff format
		if ($data->{'gff'}) {
			# write gff statement
			$fh->print("##gff-version $data->{gff}\n");
		}
		
		# Write the primary headers
		unless ($extension =~ m/gff|bed|bdg|sgr|kgg|cdt/i) {
			# we only write these for text files, not gff or bed files
			
			if ($data->{'program'}) {
				# write program header if present
				$fh->print('# Program ' . $data->{'program'} . "\n");
			}
			if ($data->{'db'}) {
				# write database header if present
				$fh->print('# Database ' . $data->{'db'} . "\n");
			}
			if ($data->{'feature'}) {
				# write feature header if present
				$fh->print('# Feature ' . $data->{'feature'} . "\n");
			}
		}
		
		# Write the miscellaneous headers
		foreach ( @{ $data->{'other'} } ) {
			# write remaining miscellaneous header lines if present
			# we do this for all files
			
			unless (/\n$/s) {
				# append newline if not present
				$_ .= "\n";
			}
			# check for comment character at beginning
			if (/^#/) {
				$fh->print($_);
			}
			else {
				$fh->print("# " . $_);
			}
		}
	
		# Write the column metadata headers
		for (my $i = 0; $i < $data->{'number_columns'}; $i++) {
			# each column metadata in the hash is referenced by the column's
			# index number as the key
			# we will take each index one at a time in increasing order
			
			# some files do not need or tolerate metadata lines, for those 
			# known files the metadata lines will be skipped
			
			# these column metadata lines do not need to be written if they
			# only have two values, presumably name and index, for files 
			# that don't normally have column headers, e.g. gff
			if ($extension =~ /sgr|kgg|cdt/i) {
				# these do not need metadata
				next;
			}
			elsif (
				exists $data->{$i}{'AUTO'} and
				scalar( keys %{ $data->{$i} } ) == 
					$data->{$i}{'AUTO'}
			) {
				# some of the metadata values were autogenerated and 
				# we have the same number of keys as were autogenerated
				# no need to write these
				next;
			}
			elsif (
				$extension =~ m/gff|bed|bdg/i and
				scalar( keys %{ $data->{$i} } ) == 2
			) {
				# only two metadata keys exist, name and index
				# GFF and BED files do not these to be written
				# so skip
				next;
			}
			
			# we will put each key=value pair into @pairs, listed asciibetically
			my @pairs; # an array of the key value pairs from the metadata hash
			# put name first
			# we are no longer writing the index number
			push @pairs, 'name=' . $data->{$i}{'name'};
			# put remainder in alphabetical order
			foreach (sort {$a cmp $b} keys %{ $data->{$i} } ) {
				next if $_ eq 'name'; # already written
				next if $_ eq 'index'; # internal use only
				next if $_ eq 'AUTO'; # internal use only
				next if $_ eq 'strand_style'; # internal use only
				push @pairs,  $_ . '=' . $data->{$i}{$_};
			}
			
			# Finally write the header line, joining the pairs with a 
			# semi-colon into a single string.
			# The column identifier is comprised of the word 'Column' 
			# and the index number joined by '_'.
			$fh->print("# Column_$i ", join(";", @pairs), "\n");
		}
	}
	
	
	# Write the table column headers
	if (
		$data->{'headers'} and
		$data->{'gff'} == 0 and
		$data->{'bed'} == 0 and
		$extension !~ /sgr/i
	) {
		# table headers existed in original source file, 
		# and this is not a GFF, BED, or SGR file,
		# therefore headers should be written
		$fh->print( 
			join("\t", @{ $data->{'data_table'}[0] }), "\n");
	}
		
	
	# Write the data table
	if ($args{'format'} eq 'simple') {
		
		# the simple format will strip the non-value '.' from the table
		for (my $i = 1; $i <= $data->{'last_row'}; $i++) {
			# we will step though the data_table array one row at a time
			# convert the non-value '.' to undefined
			# and print using a tab-delimited format
			my @linedata;
			foreach ( @{ $data->{'data_table'}[$i] }) {
				if ($_ eq '.') {
					push @linedata, q{}; # an undefined value
				} else {
					push @linedata, $_;
				}
			}
			$fh->print(join("\t", @linedata) . "\n");
		}
	}
	
	else {
		# normal data files
		for (my $i = 1; $i <= $data->{'last_row'}; $i++) {
			# we will step though the data_table array one row at a time
			# we will join each row's array of elements into a string to print
			# using a tab-delimited format
			$fh->print( 
				join("\t", @{ $data->{'data_table'}[$i] }), "\n");
		}
	}
	
	# done writing
	$fh->close;
	
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
	
	
	# Open filehandle object as appropriate
	my $fh; # filehandle
	if ($filename =~ /\.gz$/i) {
		# the file is compressed with gzip
		$fh = IO::File->new("gzip -dc $filename |") or 
			carp "unable to read '$filename' $!\n";
	} 
	elsif ($filename =~ /\.bz2$/i) {
		# the file is compressed with bzip2
		$fh = IO::File->new("bzip2 -dc $filename |") or 
			carp "unable to read '$filename' $!\n";
	} 
	else {
		# the file is uncompressed and space hogging
		$fh = IO::File->new($filename, 'r') or 
			carp "unable to read '$filename' $!\n";
	}
	return $fh if defined $fh;	
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
	
	# check filename length
	# assuming a maximum of 256, at least on Mac with HFS+, don't know about Linux
	my $name = fileparse($filename);
	if (length $name > 255) {
		carp " filename is too long! please shorten\n";
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
	
	# add gz extension if necessary
	if ($gz and $filename !~ m/\.gz$/i) {
		$filename .= '.gz';
	}
	
	
	# Generate appropriate filehandle object
	my $fh;
	if (not $gz and not $append) {
		$fh = IO::File->new($filename, 'w') or 
			carp "cannot write to file '$filename' $!\n";
	}
	elsif ($gz and !$append) {
		$fh = IO::File->new("| gzip >$filename") or 
			carp "cannot write to compressed file '$filename' $!\n";
	}
	elsif (not $gz and $append) {
		$fh = IO::File->new(">> $filename") or 
			carp "cannot append to file '$filename' $!\n";
	}
	elsif ($gz and $append) {
		$fh = IO::File->new("| gzip >>$filename") or 
			carp "cannot append to compressed file '$filename' $!\n";
	}
	return $fh if defined $fh;
}





#### Convert a tim data table into GFF format

sub convert_genome_data_2_gff_data {
	# a subroutine to convert the data table format from genomic bins or
	# windows to a gff format for writing as a gff file
	
	# get passed arguments
	my %args = @_; 
	unless (%args) {
		cluck "no arguments passed!";
		return;
	}
	
	# check data structure
	$args{'data'} ||= undef;
	my $data = $args{'data'};
	unless (verify_data_structure($data) ) {
		cluck "bad data structure!";
		return;
	}
	
	
	### Establish general gff variables
	
	# chromosome
	my $chr_index;
	if (
		exists $args{'chromo'} and 
		$args{'chromo'} =~ /^\d+$/ and
		exists $data->{ $args{'chromo'} }
	) {
		$chr_index = $args{'chromo'};
	}
	else {
		$chr_index = find_column_index($data, '^chr|seq|refseq');
	}
		
	# start position
	my $start_index;
	if (
		exists $args{'start'} and 
		$args{'start'} =~ /^\d+$/ and
		exists $data->{ $args{'start'} }
	) {
		$start_index = $args{'start'};
	}
	else {
		$start_index = find_column_index($data, 'start');
	}
		
	# stop position
	my $stop_index;
	if (
		exists $args{'stop'} and 
		$args{'stop'} =~ /^\d+$/ and
		exists $data->{ $args{'stop'} }
	) {
		$stop_index = $args{'stop'};
	}
	else {
		$stop_index = find_column_index($data, 'stop|end');
	}
	
	
	# check that we have required coordinates
	unless ( defined $chr_index ) {
		cluck " unable to identify chromosome index!";
		return;
	}
	unless ( defined $start_index ) {
		cluck " unable to identify start index!";
		return;
	}
	
	# score
	my $score_index;
	if (
		exists $args{'score'} and 
		$args{'score'} =~ /^\d+$/ and
		exists $data->{ $args{'score'} }
	) {
		$score_index = $args{'score'};
	}
	
	# name
	my $name;
	my $name_index;
	if (exists $args{'name'} and $args{'name'} ne q()) {
		if (
			$args{'name'} =~ /^\d+$/ and
			exists $data->{ $args{'name'} }
		) {
			# name is a single digit, most likely a index
			$name_index = $args{'name'};
		}
		else {
			# name is likely a string
			$name = _escape( $args{'name'} );
		}
	}
	
	# strand
	my $strand_index;
	if (
		exists $args{'strand'} and 
		$args{'strand'} =~ /^\d+$/ and
		exists $data->{ $args{'strand'} }
	) {
		$strand_index = $args{'strand'};
	}
	
	# set gff version, default is 3
	my $gff_version = $args{'version'} || 3;
	
	
	# get array of tag indices
	my @tag_indices;
	if (exists $args{'tags'}) {
		@tag_indices = @{ $args{'tags'} };
	}
	
	# identify the unique ID index
	my $id_index;
	if (
		exists $args{'id'} and 
		$args{'id'} =~ /^\d+$/ and
		exists $data->{ $args{'id'} }
	) {
		$id_index = $args{'id'} ;
	}
	
	# reference to the data table
	my $data_table = $data->{'data_table'};
	
	
	### Identify default values
	
	# gff source tag
	my ($source, $source_index);
	if (exists $args{'source'} and $args{'source'} ne q() ) {
		# defined in passed arguments
		if (
			$args{'source'} =~ /^\d+$/ and 
			exists $data->{ $args{'source'} }
		) {
			# looks like an index
			$source_index = $args{'source'};
		}
		else {
			# a text string
			$source = $args{'source'};
		}
	}
	else {
		# the default is data
		$source = 'data';
	}
	
	# gff method or type column
	my ($method, $method_index);
	if (exists $args{'method'} and $args{'method'} ne q() ) {
		# defined in passed arguments
		if (
			$args{'method'} =~ /^\d+$/ and
			exists $data->{ $args{'method'} }
		) {
			# the method looks like a single digit, most likely an index value
			$method_index = $args{'method'};
		}
		else {
			# explicit method string
			$method = $args{'method'};
		}
	}
	elsif (exists $args{'type'} and $args{'type'} ne q() ) {
		# defined in passed arguments, alternate name
		if (
			$args{'type'} =~ /^\d+$/ and
			exists $data->{ $args{'type'} }
		) {
			# the method looks like a single digit, most likely an index value
			$method_index = $args{'type'};
		}
		else {
			# explicit method string
			$method = $args{'type'};
		}
	}
	elsif (defined $name) {
		# the name provided
		$method = $name;
	}
	elsif (defined $name_index) {
		# the name of the dataset for the features' name
		$method = $data->{$name_index}{'name'};
	}
	elsif (defined $score_index) {
		# the name of the dataset for the score
		$method = $data->{$score_index}{'name'};
	}
	else {
		$method = 'Experiment';
	}
	
	# fix method and source if necessary
	# replace any whitespace or dashes with underscores
	$source =~ s/[\s\-]/_/g if defined $source;
	$method =~ s/[\s\-]/_/g if defined $method;
	
	
	### Other processing
	# convert the start postion to 1-based from 0-based
	my $convert_zero_base = $args{'zero'} || 0;
	
	
	### Reorganize the data table
		# NOTE: this will destroy any information that may be here and not
		# included in the gff data
		# since we're working with referenced data, you better hope that you
		# don't want this data later....
	
	# relabel the data table headers
	$data_table->[0] = [ 
		qw( Chromosome Source Type Start Stop Score Strand Phase Group) 
	];
	
	# re-write the data table
	for my $row (1..$data->{'last_row'}) {
		
		# collect coordinate information
		my $refseq = $data_table->[$row][$chr_index];
		my $start = $data_table->[$row][$start_index];
		my $stop;
		if (defined $stop_index) {
			$stop = $data_table->[$row][$stop_index];
		}
		else {
			$stop = $start;
		}
		if ($convert_zero_base) {
			# coordinates are 0-based, shift the start postion
			$start += 1;
		}
		if ($args{'midpoint'} and $start != $stop) {
			# if the midpoint is requested, then assign the midpoint to both
			# start and stop
			my $position = sprintf "%.0f", ( ($start + $stop) / 2 );
			$start = $position;
			$stop = $position;
		}
		
		# collect strand information
		my $strand;
		if (defined $strand_index) {
			my $value = $data_table->[$row][$strand_index];
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
			$gff_method = $data_table->[$row][$method_index];
		}
		else {
			# explicit method
			$gff_method = $method;
		}
		
		# collect source tag
		my $gff_source;
		if (defined $source_index) {
			# variable source tag
			$gff_source = $data_table->[$row][$source_index];
		}
		else {
			# static source tag
			$gff_source = $source;
		}
		
		# collect score information
		my $score;
		if (defined $score_index) {
			$score = $data_table->[$row][$score_index];
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
				$group = 'ID=' . $data_table->[$row][$id_index] . ';';
			}
			
			# define and record the GFF Name
			if (defined $name_index) {
				# a name is provided for each feature
				$group .= 'Name=' . _escape( 
					$data_table->[$row][$name_index] );
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
				$group = "$gff_method \"" . $data_table->[$row][$name_index] . "\"";
			}
			else {
				$group = "Experiment $gff_method";
			}
		}
		
		# add group tag information if present
		foreach (@tag_indices) {
			unless ($data_table->[$row][$_] eq '.') {
				# no tag if null value
				$group .= ';' . lc($data->{$_}{name}) . '=' . 
					_escape( $data_table->[$row][$_] );
			}
		}
		
		# rewrite in gff format
		$data_table->[$row] = [ (
			$refseq,
			$gff_source, 
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
	
	# from the Bio::ToolBox::db_helper, get_new_genome_list() only really has useful
	# metadata from the start column, index 1
	# also keep any metadata from the score and name columns, if defined
	
	# keep some current metadata
	my $start_metadata_ref = $data->{$start_index};
	my $score_metadata_ref; # new empty hashes
	my $group_metadata_ref;
	if (defined $score_index) {
		$score_metadata_ref = $data->{$score_index};
	}
	if (defined $name_index) {
		$group_metadata_ref = $data->{$name_index};
	}
	
	# delete old metadata
	for (my $i = 0; $i < $data->{'number_columns'}; $i++) {
		# delete the existing metadata hashes
		# they will be replaced with new ones
		delete $data->{$i};
	}
	
	# define new metadata
	$data->{0} = {
		'name'  => 'Chromosome',
		'index' => 0,
		'AUTO'  => 3,
	};
	$data->{1} = {
		'name'  => 'Source',
		'index' => 1,
		'AUTO'  => 3,
	};
	$data->{2} = {
		'name'  => 'Type',
		'index' => 2,
		'AUTO'  => 3,
	};
	$data->{3} = $start_metadata_ref;
	$data->{3}{'name'} = 'Start';
	$data->{3}{'index'} = 3;
	if (keys %{ $data->{3} } == 2) {
		$data->{3}{'AUTO'} = 3;
	}
	$data->{4} = {
		'name'  => 'Stop',
		'index' => 4,
		'AUTO'  => 3,
	};
	$data->{5} = $score_metadata_ref;
	$data->{5}{'name'} = 'Score';
	$data->{5}{'index'} = 5;
	if (keys %{ $data->{5} } == 2) {
		$data->{5}{'AUTO'} = 3;
	}
	$data->{6} = {
		'name'  => 'Strand',
		'index' => 6,
		'AUTO'  => 3,
	};
	$data->{7} = {
		'name'  => 'Phase',
		'index' => 7,
		'AUTO'  => 3,
	};
	$data->{8} = $group_metadata_ref;
	$data->{8}{'name'} = 'Group';
	$data->{8}{'index'} = 8;
	if (keys %{ $data->{8} } == 2) {
		$data->{8}{'AUTO'} = 3;
	}
	
	# reset the number of columns
	$data->{'number_columns'} = 9;
	
	# set the gff metadata to write a gff file
	$data->{'gff'} = $gff_version;
	
	# reset feature
	$data->{'feature'} = 'region';
	
	# set headers to false
	$data->{'headers'} = 0;
	
	# success
	return 1;
}


#### Export a tim data table to GFF file

sub convert_and_write_to_gff_file {
	# a subroutine to export the data table format from genomic bins or
	# windows to a gff file
	
	# get passed arguments
	my %args = @_; 
	unless (%args) {
		cluck "no arguments passed!";
		return;
	}
	
	# check data structure
	$args{'data'} ||= undef;
	my $data = $args{'data'};
	unless (verify_data_structure($data) ) {
		cluck "bad data structure!";
		return;
	}
	my $data_table = $data->{'data_table'};
	
	
	## Establish general variables
	
	# chromosome
	my $chr_index;
	if (
		exists $args{'chromo'} and 
		$args{'chromo'} =~ /^\d+$/ and
		exists $data->{ $args{'chromo'} }
	) {
		$chr_index = $args{'chromo'};
	}
	else {
		$chr_index = find_column_index($data, '^chr|seq|refseq');
	}
		
	# start position
	my $start_index;
	if (
		exists $args{'start'} and 
		$args{'start'} =~ /^\d+$/ and
		exists $data->{ $args{'start'} }
	) {
		$start_index = $args{'start'};
	}
	else {
		$start_index = find_column_index($data, 'start');
	}
		
	# stop position
	my $stop_index;
	if (
		exists $args{'stop'} and 
		$args{'stop'} =~ /^\d+$/ and
		exists $data->{ $args{'stop'} }
	) {
		$stop_index = $args{'stop'};
	}
	else {
		$stop_index = find_column_index($data, 'stop|end');
	}
	
	# check that we have coordinates
	unless ( defined $chr_index ) {
		cluck " unable to identify chromosome index!";
		return;
	}
	unless ( defined $start_index ) {
		cluck " unable to identify start index!";
		return;
	}
	
	# score
	my $score_index;
	if (
		exists $args{'score'} and 
		$args{'score'} =~ /^\d+$/ and
		exists $data->{ $args{'score'} }
	) {
		$score_index = $args{'score'};
	}
	
	# name
	my $name;
	my $name_index;
	if (exists $args{'name'} and $args{'name'} ne q()) {
		if (
			$args{'name'} =~ /^\d+$/ and
			exists $data->{ $args{'name'} }
		) {
			# name is a single digit, most likely a index
			$name_index = $args{'name'};
		}
		else {
			# name is likely a string
			$name = _escape( $args{'name'} );
		}
	}
	
	# strand
	my $strand_index;
	if (
		exists $args{'strand'} and 
		$args{'strand'} =~ /^\d+$/ and
		exists $data->{ $args{'strand'} }
	) {
		$strand_index = $args{'strand'};
	}
	
	# GFF file version, default is 3
	my $gff_version = $args{'version'} || 3;
	
	# get array of tag indices
	my @tag_indices;
	if (exists $args{'tags'}) {
		@tag_indices = @{ $args{'tags'} };
	}
	
	# identify the unique ID index
	my $id_index;
	if (
		exists $args{'id'} and 
		$args{'id'} =~ /^\d+$/ and
		exists $data->{ $args{'id'} }
	) {
		$id_index = $args{'id'} ;
	}
	
	
	
	## Set default gff data variables
	
	# gff source tag
	my ($source, $source_index);
	if (exists $args{'source'} and $args{'source'} ne q() ) {
		# defined in passed arguments
		if (
			$args{'source'} =~ /^\d+$/ and 
			exists $data->{ $args{'source'} }
		) {
			# looks like an index
			$source_index = $args{'source'};
		}
		else {
			# a text string
			$source = $args{'source'};
		}
	}
	else {
		# the default is data
		$source = 'data';
	}
	
	# gff method or type column
	my ($method, $method_index);
	if (exists $args{'method'} and $args{'method'} ne q() ) {
		# defined in passed arguments
		if (
			$args{'method'} =~ /^\d+$/ and
			exists $data->{ $args{'method'} }
		) {
			# the method looks like a single digit, most likely an index value
			$method_index = $args{'method'};
		}
		else {
			# explicit method string
			$method = $args{'method'};
		}
	}
	elsif (exists $args{'type'} and $args{'type'} ne q() ) {
		# defined in passed arguments, alternate name
		if (
			$args{'type'} =~ /^\d+$/ and
			exists $data->{ $args{'type'} }
		) {
			# the method looks like a single digit, most likely an index value
			$method_index = $args{'type'};
		}
		else {
			# explicit method string
			$method = $args{'type'};
		}
	}
	elsif (defined $name) {
		# the name provided
		$method = $name;
	}
	elsif (defined $name_index) {
		# the name of the dataset for the features' name
		$method = $data->{$name_index}{'name'};
	}
	elsif (defined $score_index) {
		# the name of the dataset for the score
		$method = $data->{$score_index}{'name'};
	}
	else {
		$method = 'Experiment';
	}
	# fix method and source if necessary
	# replace any whitespace or dashes with underscores
	$source =~ s/[\s\-]/_/g if defined $source;
	$method =~ s/[\s\-]/_/g if defined $method;
	
	
	## Open output file
	# get the filename
	my $filename;
	if ( $args{'filename'} ne q() ) {
		$filename = $args{'filename'};
		# remove unnecessary extensions
		$filename =~ s/\.gz$//;
		$filename =~ s/\.txt$//;
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
	elsif (defined $data->{'basename'}) {
		# use the base file name for lack of a better name
		$filename = $data->{'basename'} . '.gff';
	}
	else {
		# what, still no name!!!????
		$filename = 'your_stupid_output_gff_file_with_no_name.gff';
	}
	if ($gff_version == 3) {
		$filename .= '3'; # make extension gff3
	}
	
	# open the file for writing 
	my $gz = $args{'gz'};
	my $output_gff = open_to_write_fh($filename, $gz);
	
	# write basic headers
	print {$output_gff} "##gff-version $gff_version\n";
	if (exists $data->{'filename'}) {
		# record the original file name for reference
		print {$output_gff} "# Exported from file '", 
			$data->{'filename'}, "'\n";
	}
	
	
	### Write the column metadata headers
	# write the metadata lines only if there is useful information
	# and only for the pertinent columns (chr, start, score)
	# we will check the relavent columns for extra information beyond
	# that of name and index
	# we will write the metadata then for that dataset that is being used
	# substituting the column name and index appropriate for a gff file
	
	# check the chromosome metadata
	if (scalar( keys %{ $data->{$chr_index} } ) > 2) {
		# chromosome has extra keys of info
		print {$output_gff} "# Column_0 ";
		my @pairs;
		foreach (sort {$a cmp $b} keys %{ $data->{$chr_index} } ) {
			if ($_ eq 'index') {
				next;
			}
			elsif ($_ eq 'name') {
				push @pairs, "name=Chromosome";
			}
			else {
				push @pairs,  $_ . '=' . $data->{$chr_index}{$_};
			}
		}
		print {$output_gff} join(";", @pairs), "\n";
	}
	
	# check the start metadata
	if (scalar( keys %{ $data->{$start_index} } ) > 2) {
		# start has extra keys of info
		print {$output_gff} "# Column_3 ";
		my @pairs;
		foreach (sort {$a cmp $b} keys %{ $data->{$start_index} } ) {
			if ($_ eq 'index') {
				next;
			}
			elsif ($_ eq 'name') {
				push @pairs, "name=Start";
			}
			else {
				push @pairs,  $_ . '=' . $data->{$start_index}{$_};
			}
		}
		print {$output_gff} join(";", @pairs), "\n";
	}
	
	# check the score metadata
	if (
		defined $score_index and
		scalar( keys %{ $data->{$score_index} } ) > 2
	) {
		# score has extra keys of info
		print {$output_gff} "# Column_5 ";
		my @pairs;
		foreach (sort {$a cmp $b} keys %{ $data->{$score_index} } ) {
			if ($_ eq 'index') {
				next;
			}
			elsif ($_ eq 'name') {
				push @pairs, "name=Score";
			}
			else {
				push @pairs,  $_ . '=' . $data->{$score_index}{$_};
			}
		}
		print {$output_gff} join(";", @pairs), "\n";
	}
	
	# check the name metadata
	if (
		defined $name_index and
		scalar( keys %{ $data->{$name_index} } ) > 2
	) {
		# score has extra keys of info
		print {$output_gff} "# Column_8 ";
		my @pairs;
		foreach (sort {$a cmp $b} keys %{ $data->{$name_index} } ) {
			if ($_ eq 'index') {
				next;
			}
			elsif ($_ eq 'name') {
				push @pairs, "name=Group";
			}
			else {
				push @pairs,  $_ . '=' . $data->{$name_index}{$_};
			}
		}
		print {$output_gff} join(";", @pairs), "\n";
	}
	
			
	
	### Write the gff features
	for my $row (1..$data->{'last_row'}) {
		
		# collect coordinate information
		my $refseq = $data_table->[$row][$chr_index];
		my $start = $data_table->[$row][$start_index];
		my $stop;
		if (defined $stop_index) {
			$stop = $data_table->[$row][$stop_index];
		}
		else {
			$stop = $start;
		}
		if ($args{'midpoint'} and $stop != $stop) {
			# if the midpoint is requested, then assign the midpoint to both
			# start and stop
			my $position = sprintf "%.0f", ($start + $stop)/2;
			$start = $position;
			$stop = $position;
		}
		
		# collect score information
		my $score;
		if (defined $score_index) {
			$score = $data_table->[$row][$score_index];
		}
		else {
			$score = '.';
		}
		
		# collect strand information
		my $strand;
		if (defined $strand_index) {
			my $value = $data_table->[$row][$strand_index];
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
			$gff_method = $data_table->[$row][$method_index];
		}
		else {
			# explicit method
			$gff_method = $method;
		}
		
		# collect source tag
		my $gff_source;
		if (defined $source_index) {
			# variable source tag
			$gff_source = $data_table->[$row][$source_index];
		}
		else {
			# static source tag
			$gff_source = $source;
		}
		
		# collect group information based on version
		my $group;
		if ($gff_version == 3) {
			
			# define and record GFF ID
			if (defined $id_index) {
				# this assumes that the $id_index values are all unique
				# user's responsibility to fix it otherwise
				$group = 'ID=' . $data_table->[$row][$id_index] . ';';
			}
			
			# define and record the GFF Name
			if (defined $name_index) {
				# a name is provided for each feature
				$group .= 'Name=' . _escape( 
					$data_table->[$row][$name_index] );
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
				$group = "$gff_method \"" . $data_table->[$row][$name_index] . "\"";
			}
			else {
				# really generic
				$group = "Experiment \"$gff_method\"";
			}
		}
		
		# add group tag information if present
		foreach (@tag_indices) {
			unless ($data_table->[$row][$_] eq '.') {
				# no tag if null value
				$group .= ';' . lc($data->{$_}{name}) . '=' . 
					_escape( $data_table->[$row][$_] );
			}
		}
		
		# Write gff feature
		print {$output_gff} join("\t", (
			$refseq,
			$gff_source, 
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
	my %args = @_; 
	unless (%args) {
		cluck "no arguments passed!";
		return;
	}
	
	# parameters
	my $outfile =        $args{'filename'}    || undef;
	my $dataset =        $args{'dataset'}     || undef;
	my $startcolumn =    $args{'startcolumn'} || undef;
	my $endcolumn =      $args{'endcolumn'}   || $args{'stopcolumn'}  || undef;
	my $log =            $args{'log'}         || undef;
	
	
	
	# Check required values
	my $data = $args{'data'} ||= undef;
	unless (defined $data) {
		cluck "no data structure passed!\n";
		return;
	}
	unless (verify_data_structure($data) ) {
		cluck "bad data structure!";
		return;
	}
	unless (defined $outfile) {
		if (defined $data->{'basename'}) {
			# use the opened file's filename if it exists
			# prepend the path if it exists
			# the extension will be added later
			$outfile = $data->{'path'} . $data->{'basename'};
		}
		else {
			cluck "no filename passed to write_summary_data!\n";
			return;
		}
	}
	
	
	
	# Auto calculate missing arguments
	unless (defined $startcolumn) {
		# we will attempt to determine where to start the summing process
		# we will skip those columns with common feature-description names
		
		my @acceptable_indices; # array of acceptable indices
		my %skip = (
			# these are example names generated from subs in Bio::ToolBox::db_helper
			'systematicname'  => 1,
			'name'            => 1,
			'id'              => 1,
			'alias'           => 1,
			'aliases'         => 1,
			'type'            => 1,
			'class'           => 1,
			'geneclass'       => 1,
			'chromosome'      => 1,
			'chromo'          => 1,
			'seq_id'          => 1,
			'seqid'           => 1,
			'start'           => 1,
			'stop'            => 1,
			'end'             => 1,
			'gene'            => 1,
			'strand'          => 1,
			'length'          => 1,
			'primary_id'      => 1,
		);
		
		# walk through the dataset names
		for (my $i = 0; $i < $data->{'number_columns'}; $i++) {
			unless (exists $skip{ lc $data->{$i}{'name'} } ) {
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
		$endcolumn = $data->{'number_columns'} - 1;
	}
	unless ($dataset) {
		# the original dataset name (i.e. the name of the dataset in the 
		# database from which the column's data was derived) should be the same 
		# in all the columns
		$dataset = $data->{$startcolumn}{'dataset'} || 'data_scores';
	}
	unless (defined $log) {
		# the log flag should be set in the column metadata and should be the
		# same in all
		$log = $data->{$startcolumn}{'log2'} || 0;
	}
	
	# Prepare score column name
		# we will use the basename of the output file name to make it 
		# easier in downstream applications
	my ($data_name, undef, undef) = fileparse($outfile, @SUFFIX_LIST);
	
	# Prepare array to store the summed data
	my $summed_data = generate_tim_data_structure(
		'averaged_windows', 
		'Window',
		'Midpoint',
		$data_name
	);
	$summed_data->{'db'} = $data->{'database'};
	$summed_data->{0}{'number_features'} = $data->{'last_row'};
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
		my $midpoint = int mean(
			# this assumes the column metadata has start and stop
			$data->{$column}{'start'},	
			$data->{$column}{'stop'},	
		) or undef; 
		
		
		# collect the values in the column
		my @values;
		for my $row (1..$data->{'last_row'}) {
			my $value = $data->{'data_table'}->[$row][$column];
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
			$data->{$column}{'name'}, 
			$midpoint, 
			$window_mean
		) ];
		$summed_data->{'last_row'} += 1;
	}
	
	# Write summed data
	$outfile =~ s/\.txt(\.gz)?$//; # strip any .txt or .gz extensions if present
	my $written_file = write_tim_data_file(
		'data'      => $summed_data,
		'filename'  => $outfile . '_summed',
		'gz'        => 0,
	);
	
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


### Internal subroutine to escape special characters for GFF3 files
sub _escape {
	my $string = shift;
	# this magic incantation was borrowed from Bio::Tools::GFF
	$string =~ s/([\t\n\r%&\=;, ])/sprintf("%%%X",ord($1))/ge;
	return $string;
}


### Internal subroutine to check if a comment line contains headers
sub _commented_header_line {
	my ($inputdata, $line) = @_;
	
	# prepare arrays from the other lines and current line
	my @commentdata;
	if ( scalar @{ $inputdata->{'other'} } >= 1 ) {
		# take the last line in the other array
		@commentdata = split /\t/, $inputdata->{'other'}->[-1];
	}
	my @linedata = split /\t/, $line;
	
	# check if the counts are equal
	if (scalar @commentdata == scalar @linedata) {
		return 1;
	}
	else {
		return 0;
	}
}



__END__

=head1 NAME

Bio::ToolBox::file_helper

=head1 SYNOPSIS

  use Bio::ToolBox::file_helper qw(
    load_tim_data_file
    write_tim_data_file
    open_to_read_fh
    open_to_write_fh
  );
  
  my $input_data = load_tim_data_file($file) or die "can't open file!";
  
  my ($fh, $metadata) = open_tim_data_file($file) or die;
  
  while (my $line = $fh->getline) {
    ...
  }
  
  my $input_fh = open_to_read_fh($file);
  
  my $output_fh = open_to_write_fh($file);
  
  my $output_fh = open_to_write_fh($file, $gz, $append);
  
  my $success = write_tim_data_file(
    'data'       => $data,
    'filename'   => $file,
    'gz'         => $gz,
  );

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
'Column', followed by an underscore and the column number (0-based). 
Following this is a series of 'key=value' pairs separated by ';'. 
Spaces are generally not allowed. Obviously '=' or ';' are not 
allowed or they will interfere with the parsing. The metadata 
describes how and where the data was collected. Additionally, any 
modifications performed on the data are also recorded here. The only 
key that is required is 'name'. 
If the file being read does not contain metadata, then it will be auto 
generated with basic metadata.

=back

A list of standard column header keys is below, but is not exhaustive. 

=over

=item name

The name of the column. This should be identical to the table header.

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

=head1 USAGE

Call the module at the beginning of your perl script and pass a list of the 
desired modules to import. None are imported by default.
  
  use Bio::ToolBox::db_helper qw(load_tim_data_file write_tim_data_file);
  
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

BED and BedGraph style files, recognized by .bed or .bdg file extensions, 
have their start coordinate adjusted by +1 to convert from 0-based interbase 
numbering system to 1-based numbering format, the convention used by BioPerl. 
A metadata attribute is applied informing the user of the change. When writing 
a valid Bed or BedGraph file, converted start positions are changed back to 
interbase format.

Strand information is parsed from recognizable symbols, including "+, -, 1, 
-1, f, r, w, c, 0, .",  to the BioPerl convention of 1, 0, and -1. Valid 
BED and GFF files are changed back when writing these files. 

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
Please refer to L<FORMAT OF TIM DATA TEXT FILE> for more 
information regarding the file format. If the 'gff' key is true in the data 
hash, then a gff file will be written.

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  data     => A scalar reference to the tim data structure ad described
              in C<Bio::ToolBox::data_helper>. 
  Optional: 
  filename => A scalar value containing the name of the file to 
              write. This value is required for new data files and 
              optional for overwriting existing files (the filename 
              stored in the metadata is used). Appropriate extensions 
              are added (e.g, .txt, .gz, etc) as neccessary. 
  format   => A string to indicate the file format to be written.
              Acceptable values include 'text', and 'simple'.
              Text files are text in nature, include all metadata, and
              usually have '.txt' extensions. Simple files are
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
The arguments always take precendence over the filename extensions, however.

Example

	my $filename = 'my_data.txt.gz';
	my $data_ref = load_tim_data_file($filename);
	...
	my $success_write = write_tim_data_file(
		'data'     => $data_ref,
		'filename' => $filename,
		'format'   => 'simple',
	);
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
	$fh->print("something interesting\n");
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
generated with get_new_genome_list() in Bio::ToolBox::db_helper will generate these
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
  source   => A scalar value representing either the index of the 
              column containing values, or a text string to 
              be used as the GFF source value. Default is 'data'.
  type     => A scalar value representing either the index of the 
              column containing values, or a text string to 
              be used as the GFF type or method value. If not 
              defined, it will use the column name of the dataset 
              used for either the 'score' or 'name' column, if 
              defined. As a last resort, it will use the most 
              creative method of 'Experiment'.
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
	my $success = convert_genome_data_2_gff_data(
		'data'     => $data_ref,
		'score'    => 3,
		'midpoint' => 1,
	);
	if ($success) {
		# write a gff file
		my $success_write = write_tim_data_file(
			'data'     => $data_ref,
			'filename' => $filename,
		);
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
Bio::ToolBox::db_helper will generate these coordinate datasets. 

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
  source   => A scalar value representing either the index of the 
              column containing values, or a text string to 
              be used as the GFF source value. Default is 'data'.
  type     => A scalar value representing either the index of the 
              column containing values, or a text string to 
              be used as the GFF type or method value. If not 
              defined, it will use the column name of the dataset 
              used for either the 'score' or 'name' column, if 
              defined. As a last resort, it will use the most 
              creative method of 'Experiment'.
  method   => Alias for "type".
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
	my $success = convert_and_write_to_gff_file(
		'data'     => $data_ref,
		'score'    => 3,
		'midpoint' => 1,
		'filename' => "$filename.gff",
		'version'  => 2,
	);
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
                 from Bio::ToolBox::db_helper).
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
	my $summary_success = write_summary_data(
		'data'         => $main_data_ref,
		'filename'     => $outfile,
		'startcolumn'  => 4,
	);


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
