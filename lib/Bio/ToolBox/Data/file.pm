package Bio::ToolBox::Data::file;
our $VERSION = 1.27;

=head1 NAME

Bio::ToolBox::Data::file - File functions to Bio:ToolBox::Data family

=head1 DESCRIPTION

File methods for reading and writing data files for both L<Bio::ToolBox::Data> 
and L<Bio::ToolBox::Data::Stream> objects. This module should not be used 
directly. See the respective modules for more information.

=cut

use strict;
use Carp qw(carp cluck croak confess);
require Exporter;
use File::Basename qw(fileparse);
use IO::File;
use Statistics::Lite qw(mean min);
our $VERSION = 1.27;


### Variables
# column names
our $std_col_names = _set_standard_column_names();

# Export
our @ISA = qw(Exporter);
our @EXPORT = qw();
our @EXPORT_OK = qw(
	open_to_read_fh
	open_to_write_fh
	check_file
);

# List of acceptable filename extensions
our $suffix = qr/\.(?:txt|gff3?|gtf|bed|bdg|bedgraph|sgr|kgg|cdt|vcf|narrowpeak|broadpeak|reff?lat|genepred|ucsc)(?:\.gz)?/i;


### The True Statement
1; 



### Load new version data table from file

sub load_file {
	my ($self, $filename) = @_;
	
	# open the file and load metadata
	$self->parse_filename($filename);
	my $fh = $self->open_to_read_fh($self->filename);
	return unless $fh;
	$self->parse_headers($fh);
	$self->{data_table}->[0] eq $data->{'column_names'}; 
	
	# set metadata for converting 0-based starts to 1-based
	$self->{'0based_starts'} = [];
	if ($data->{'ucsc'} or $data->{'bed'}) {
		# same thing for ucsc files
		# adjust both transcription and coding start
		#### We should be doing thickStart and blockStarts too ####
		# but I'm not sure how useful this really would be, plus adds complexity 
		# that will slow file loading down - it's already pretty complicated
		foreach my $name (qw(start txStart cdsStart peak)) {
			my $c = $self->find_column($name);
			next unless defined $c;
			next if $self->metadata($c, 'base');
			push @{ $self->{'0based_starts'} }, $c;
		}
	}
	
	# internal strand plus/minus count for determining conversion
	$self->{plusminus_count} = 0;
	
	# Load the data table
	while (my $line = $fh->getline) {		
		# the current file position should be at the beginning of the
		# data table information
		
		# skip comment and empty lines
		if ($line =~ /^#/) {
			$self->add_comment($line);
			next;
		}
		next if $line !~ m/\w+/;
		
		# process the line
		$self->add_data_line($line);
	}
	
	# update metadata as necessary
	if ($self->{plusminus_count}) {
		# we have converted strand information
		if (exists $data->{$strand_i}{'strand_style'}) {
			# just in case, reset it to plusminus
			$data->{$strand_i}{'strand_style'} = 'plusminus';
		}
		else {
			$data->{$strand_i}{'strand_style'} = 'plusminus';
			if (exists $data->{$strand_i}{'AUTO'}) {
				# update automatically generated metadata
				$data->{$strand_i}{'AUTO'}++;
			}
		}
	}
	else {
		# no plusminus count, make sure metadata was not set automatically
		# but this is technically a problem
		if (exists $data->{$strand_i}{'strand_style'}) {
			carp "File format is suspicious! format suggests plus/minus strand format" . 
				" but none was found!?\n";
			delete $data->{$strand_i}{'strand_style'};
			$data->{$strand_i}{'AUTO'}--;
		}
	}
	foreach my $s (@starts) {
		# each column of 0-based start that has been updated
		$data->{$s}{'base'} = 1;
		if (exists $data->{$s}{'AUTO'}) {
			# update automatically generated metadata
			$data->{$s}{'AUTO'}++;
		}
	}
	delete $self->{'0based_starts'};
	delete $self->{plusminus_count};
	
	# record the index number of the last data row
	$self->{'last_row'} = scalar @datatable - 1;
	
	# completed loading the file
	$fh->close;
	
	# verify the structure
	return unless $self->verify;
	
	# finished
	return 1;
}







sub parse_headers {
	my $self = shift;
	my $fh = shift;
	my $noheader = shift || 0; # boolean to indicate no headers are present
	unless (ref $fh =~ /^IO/) {
		confess " must pass an open IO::Handle compatible object!\n";
	}
	
	$data->{'filename'} = $filename; # the original filename
	$data->{'basename'} = $basename; # filename basename
	$data->{'extension'} = $extension; # the filename extension
	$data->{'path'} = $path; # the original path
	$data->{'column_names'} = []; # array for the column names
	$data->{'headers'} = q(); # boolean to indicate column headers present
	$data->{'program'} = q(); # clear program, this will be read from file
	
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
			die "File '$filename' does not appear to have unix line endings!\n" . 
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
			$data->{'program'} = $1;
			$data->{'program'} =~ s/[\r\n]+$//;
			$header_line_count++;
		}
		
		# the source database
		elsif ($line =~ m/^# Database (.+)$/) {
			$data->{'db'} = $1;
			$data->{'db'} =~ s/[\r\n]+$//;
			$header_line_count++;
		}
		
		# the type of feature in this datafile
		elsif ($line =~ m/^# Feature (.+)$/) {
			$data->{'feature'} = $1;
			$data->{'feature'} =~ s/[\r\n]+$//;
			$header_line_count++;
		}
		
		# column or dataset specific information
		elsif ($line =~ m/^# Column_(\d+)/) {
			# the column number will become the index number
			my $index = $1; # capture this from grep above
			
			# load metadata
			_column_metadata($data, $line, $index);
			
			$header_line_count++;
		}
		
		# gff version header
		elsif ($line =~ /^##gff-version\s+([\d\.]+)$/) {
			# store the gff version in the hash
			# this may or may not be present in the gff file, but want to keep
			# it if it is
			$data->{'gff'} = $1;
			$data->{'gff'} =~ s/[\r\n]+$//;
			$header_line_count++;
		}
		
		# any other nonstandard header
		elsif ($line =~ /^#/) {
			# store in an anonymous array in the inputdata hash
			push @{ $data->{'other'} }, $line;
			$header_line_count++;
		}
		
		# a track or browser line 
		elsif ($line =~ /^(?:track|browser)\s+/i) {
			# common with wig, bed, or bedgraph files for use with
			# the UCSC genome browser
			# treat as a comment line, there's not that much useful info
			# in here, possibly the name, that's it
			push @{ $data->{'other'} }, $line;
			$header_line_count++;
		}
		
		# the remainder is the data table itself
		else {
			# the first row in the data table are (usually) the column names 
			# we only want the names, not the rest of the table
			
			# specific file formats have implicit predefined column formats
			# these file formats do NOT have column headers
			# we will first check for those file formats and process accordingly
			
			### Data tables with a commented header line 
			if ( _commented_header_line($data, $line) ) {
				# including data table files from UCSC Table Browser
				# also VCF files as well
				# these will have one comment line marked with #
				# that really contains the column headers
				
				# lots of requirements, but this checks that (column 
				# metadata has not already been loaded), there is at least one 
				# unknown comment line, and that the last comment line is 
				# splitable and has the same number of elements as the first 
				# data line
				
				# process the real header line
				my $header_line = pop @{ $data->{'other'} };
				$header_line =~ s/[\r\n]+$//;
				_commented_header_line_metadata($data, $header_line);
				
				last PARSE_HEADER_LOOP;
			}
			
			
			### a GFF file
			elsif ($extension =~ /g[tf]f/i) {
				# gff files have nine defined columns
				# there are different specifications and variants:
				# gff (v.1), gff v.2, gff v.2.5 (aka gtf), gff v.3 (gff3)
				# however, the columns are the same in all versions
				# more info on gff can be found here
				# http://gmod.org/wiki/GFF3
				
				# generate metadata for gff files
				add_gff_metadata($data);
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
			
			
			### a Bed or BedGraph file
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
				
				# generate metadata
				add_bed_metadata($data, $line);
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
			
			
			### a peak file
			elsif ($extension =~ /peak/i) {
				# three different types of peak files are available
				# see http://genome.ucsc.edu/FAQ/FAQformat.html
				
				# generate metadata
				add_peak_metadata($data, $line);
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
			
			
			### a UCSC gene table
			elsif ($extension =~ /ref+lat|genepred|ucsc/i) {
				# these are tricky, as we will try to determine contents by counting
				# the columns, which may not be accurate
				# not only that, but this presumes the extension even makes it recognizable
				# see http://genome.ucsc.edu/FAQ/FAQformat.html#format9 for details
				# also biotoolbox script ucsc_table2gff3.pl
				
				# generate metadata
				add_ucsc_metadata($data, $line);
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
			
			
			### a SGR file
			elsif ($extension =~ /sgr/i) {
				# a sgr file contains three columns: chromo, position, score
				# this is a very simple file format, useful in exporting and
				# importing to binary BAR files used in T2, USeq, and IGB
				
				# generate metadata
				add_sgr_metadata($data, $line);
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
			
			
			### standard text file with headers, i.e. everything else
			else {
				# we have not yet parsed the row of data column names
				# we will do so now
				
				# generate metadata
				add_standard_metadata($data, $line);
				
				# count as a header line
				$header_line_count++;
				
				# end this loop
				last PARSE_HEADER_LOOP;
			}
		}
		
	}
	
	# No header was requested
	if ($noheader) {
		# user indicated that no header was present
		# which means that what we used as a header is actually the first data row
		
		# fix the column names
		for (my $i = 0; $i < $data->{'number_columns'}; $i++) {
			my $name = $data->{$i}{'name'};
			$data->{$i}{'name'} = "Column$i ($name)";
			$data->{'column_names'}->[$i] = $data->{$i}{'name'};
			$data->{$i}{'AUTO'} = 3;
		}
		
		# move the counter back one line
		$header_line_count -= 1;
		
		# set headers flag to false
		$data->{'headers'} = 0; 
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
	return ($fh, $data);
	
}


sub add_data_line {
	my ($self, $line) = @_;
	
	# do not chomp the line yet, just split into an array
	my @linedata = split /\t/, $line;
	
	# check the number of elements
	if (scalar @linedata != $self->number_columns ) {
		if ($line =~ /\r/ and $line !~ /\n/) {
			die "File does not appear to have unix line endings!\n" . 
				" Please convert to unix-style line endings and try again\n";
		}
		warn "Number of columns in line is inconsistent\n";
		# we will verify the table after loading all the lines to verify 
	}
	
	# chomp the last element
	# we do this here to ensure the tab split above gets all of the values
	# otherwise trailing null values aren't included in @linedata
	# be sure to handle both newlines and carriage returns
	$linedata[-1] =~ s/[\r\n]+$//;
	
	# convert null values to internal '.'
	for (my $i = 0; $i < $self->number_columns; $i++) {
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
	
	# convert strand as necessary
	my $strand_i = $self->strand_column;
	if (defined $strand_i) {
		# convert any interpretable value to signed value
		if ($linedata[$strand_i] eq '+') {
			# just a simple plus
			$linedata[$strand_i] = 1;
			$self->{plusminus_count}++;
		}
		elsif ($linedata[$strand_i] eq '-') {
			# simple minus
			$linedata[$strand_i] = -1;
			$self->{plusminus_count}++;
		}
		elsif ($linedata[$strand_i] eq '.') {
			# unstranded GFF format, not BED
			$linedata[$strand_i] = 0;
			$self->{plusminus_count}++;
		}
		# otherwise assume bioperl convention -1, 0, 1
		# if it is not, then hope for the best
		# I am dropping support for the ancient forward, watson, reverse, crick
		# who uses those anyway?????
		# some old bioperl scripts may still support it
	}
	
	# adjust start positions
	foreach my $s (@{ $self->{'0based_starts'} }) {
		$linedata[$s] += 1;
	}
	
	# add the line of data
	push @{ $self->{data_table} }, \@linedata;
	return 1;
}




### Parse the filename using the list suffix list

sub parse_filename {
	my $self = shift;
	my $filename = $self->check_file(shift @_);
	my ($basename, $path, $extension) = fileparse($filename, $suffix);
	$self->{filename}  = $filename;
	$self->{basename}  = $basename;
	$self->{path}      = $path;
	$self->{extension} = $extension;
}






### Write out a data file from the data hash


sub write_data_file {
	
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
	my ($name, $path, $extension) = fileparse($args{'filename'}, $suffix);
	
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
		elsif ($extension =~ /bed|bdg|peak/i) {
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
			elsif ($data->{'extension'} =~ /bed|bdg|peak/i) {
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
	if ($extension =~ /\.bed|bdg|peak/i and $data->{'bed'} > 0) {
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
			if (exists $data->{1}{'AUTO'}) {
				$data->{1}{'AUTO'} -= 1;
			}
		}
	}
	if ($data->{'ucsc'} > 0) {
		my $t = find_column_index($data, 'txStart');
		my $c = find_column_index($data, 'cdsStart');
		
		# check if we need to convert these, probably do, but just in case
		my ($do_t, $do_c);
		$do_t = 1 if (exists $data->{$t}{'base'} and $data->{$t}{'base'} == 1);
		$do_t = 1 if (exists $data->{$c}{'base'} and $data->{$c}{'base'} == 1);
		
		# convert back to interbase
		for (my $row = 1; $row <= $data->{'last_row'}; $row++) {
			# subtract 1 to each start position
			$data->{'data_table'}->[$row][$t] -= 1 if $do_t;
			$data->{'data_table'}->[$row][$c] -= 1 if $do_c;
		}
		if ($do_t) {
			delete $data->{$t}{'base'};
			$data->{$t}{'AUTO'} -= 1 if exists $data->{$t}{'AUTO'};
		}
		if ($do_c) {
			delete $data->{$c}{'base'};
			$data->{$c}{'AUTO'} -= 1 if exists $data->{$c}{'AUTO'};
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
		unless ($extension =~ m/gff|bed|bdg|sgr|kgg|cdt|peak/i) {
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
				$extension =~ m/gff|bed|bdg|peak/i and
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
	my ($self, $file) = @_;
	my $filename = $self->check_file($file);
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
	my ($self, $filename, $gz, $append) = @_;
	
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
	my ($data_name, undef, undef) = fileparse($outfile, $suffix);
	
	# Prepare array to store the summed data
	my $summed_data = generate_data_structure(
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
	my $written_file = write_data_file(
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





### Subroutine to check for file existance
sub check_file {
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
		foreach my $ext (qw(.gz .txt .txt.gz .bed .bed.gz)) {
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
	my ($data, $line) = @_;
	
	# prepare arrays from the other lines and current line
	my @commentdata;
	if ( scalar @{ $data->{'other'} } >= 1 ) {
		# take the last line in the other array
		@commentdata = split /\t/, $data->{'other'}->[-1];
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


### Internal subroutine to process metadata for standard columns
sub _column_metadata {
	my ($data, $line, $index) = @_;
	
	# strip the Column metadata identifier
	$line =~ s/[\r\n]+$//;
	$line =~ s/^# Column_\d+ //; 
	
	# break up the column metadata
	my %temphash; # a temporary hash to put the column metadata into
	foreach (split /;/, $line) {
		my ($key, $value) = split /=/;
		if ($key eq 'index') {
			if ($index != $value) {
				# the value from the metadata index key should be 
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
	if (exists $data->{$index}) {
		# we will simply overwrite the previous metadata hash
		# harsh, I know, but what to do?
		# if it was canned metadata for a gff file, that's ok
		warn "Warning: more than one metadata line exists for index $index!\n";
		$data->{$index} = \%temphash;
	}
	else {
		# metadata hash doesn't exist, so we will add it
		$data->{$index} = \%temphash;
		# also update the number of columns
		$data->{'number_columns'} += 1;
	}
}

### Internal subroutine to generate metadata for commented header line
sub _commented_header_line_metadata {
	my ($data, $header_line) = @_;
	
	# generate the metadata
	my $i = 0;
	foreach (split /\t/, $header_line) {
		$data->{$i} = { 
			'name'   => $_,
			'index'  => $i,
			'AUTO'   => 3,
		} unless exists $data->{$i};
		
		push @{ $data->{'column_names'} }, $_;
		$i++;
	}
	$data->{'number_columns'} = $i;
	
	# we do not count the current line as a header
	
	# set headers flag to true
	$data->{'headers'} = 1;
}


### Subroutine to generate metadata for gff files
sub add_gff_metadata {
	my $data = shift;
	
	# set column number
	$data->{'number_columns'} = 9; # supposed to be 9 columns
	
	# set the gff version based on the extension
	unless ($data->{'gff'}) {
		$data->{'gff'} = 
			$data->{extension} =~ /gtf/  ? 2.5 :
			$data->{extension} =~ /gff3/ ? 3   :
			2;
	}
	
	# set the metadata for the each column
		# some of these may already be defined if there was a 
		# column metadata specific column in the file
	for (my $i = 0; $i < 9; $i++) {
		# loop for each column
		# set metadata unless it's already loaded
		unless (exists $data->{$i}) {
			$data->{$i}{'name'}  = $std_col_names->{gff}->[$i];
			$data->{$i}{'index'} = $i;
			$data->{$i}{'AUTO'}  = 3;
		}
		# assign the name to the column header
		$data->{'column_names'}->[$i] = 
			$data->{$i}{'name'};
	}
	
	# set strand style
	$data->{6}{'strand_style'} = 'plusminus';
	$data->{6}{'AUTO'}++;
	
	# set headers flag to false
	$data->{'headers'} = 0;
	
	# set the feature type
	unless (defined $data->{'feature'}) {
		$data->{'feature'} = 'region';
	}
}


### Subroutine to generate metadata for BED files
sub add_bed_metadata {
	my ($data, $line) = @_;

	# first determine the number of columns we're working
	my @elements = split /\s+/, $line;
		# normally tab-delimited, but the specs are not explicit
	my $column_count = scalar @elements;
	$data->{'number_columns'} = $column_count; 
	$data->{'bed'} = $column_count;
	
	# malformed bed file
	if ($column_count < 3) {
		warn " file $data->{filename} appears to be malformed! Extension suggests " . 
			" BED file but only has $column_count columns!\n";
	}
	
	# set column names
	my $bed_names = $std_col_names->{bed12};
	
	# bedGraph files
	if ($data->{extension} =~ /bdg|graph/i) {
		# fourth column is Score
		$bed_names = $std_col_names->{bdg};
	}
	
	# set the metadata for each column
		# some of these may already be defined if there was a 
		# column metadata specific column in the file
	for (my $i = 0; $i < $column_count; $i++) {
		# loop for each column
		# set name unless it already has one from metadata
		unless (exists $data->{$i}) {
			$data->{$i}{'name'}  = $bed_names->[$i] || 'extraColumn';
			$data->{$i}{'index'} = $i;
			$data->{$i}{'AUTO'}  = 3;
		}
		# assign the name to the column header
		$data->{'column_names'}->[$i] = 
			$data->{$i}{'name'};
	}
	
	# set strand style
	if ($column_count >= 6) {
		$data->{5}{'strand_style'} = 'plusminus';
		$data->{5}{'AUTO'}++;
	}
	
	# set the feature type
	unless (defined $data->{'feature'}) {
		$data->{'feature'} = 'region';
	}
	
	# set headers flag to false
	$data->{'headers'} = 0;
}


### Subroutine to generate metadata for broadpeak and narrowpeak files
sub add_peak_metadata {
	my ($data, $line) = @_;

	my @elements = split /\s+/, $line;
		# normally tab-delimited, but the specs are not explicit
	my $column_count = scalar @elements;
	$data->{'number_columns'} = $column_count; 
	$data->{'bed'} = $column_count;
		# technically this is called bed6+4 or bed6+3, but for our 
		# purposes here, we will stick to column count to avoid breaking
	
	# determine type of peak file based on extension
	# narrowPeak
	if ($data->{extension} =~ /narrowpeak/i) {
		# set the metadata for each column
			# some of these may already be defined if there was a 
			# column metadata specific column in the file
		for (my $i = 0; $i < $column_count; $i++) {
			# loop for each column
			# set name unless it already has one from metadata
			unless (exists $data->{$i}) {
				$data->{$i}{'name'} = $std_col_names->{narrowpeak}->[$i] || 'extraColumn';
				$data->{$i}{'index'} = $i;
				$data->{$i}{'AUTO'}  = 3;
			}
			# assign the name to the column header
			$data->{'column_names'}->[$i] = 
				$data->{$i}{'name'};
		}
	}
	
	# broadPeak
	elsif ($data->{extension} =~ /broadpeak/i) {
		# set the metadata for each column
			# some of these may already be defined if there was a 
			# column metadata specific column in the file
		for (my $i = 0; $i < $column_count; $i++) {
			# loop for each column
			# set name unless it already has one from metadata
			unless (exists $data->{$i}) {
				$data->{$i}{'name'} = $std_col_names->{broadpeak}->[$i] || 'extraColumn';
				$data->{$i}{'index'} = $i;
				$data->{$i}{'AUTO'}  = 3;
			}
			# assign the name to the column header
			$data->{'column_names'}->[$i] = 
				$data->{$i}{'name'};
		}
	}
	
	# some other type of peak
	else {
		# set the metadata for each column
			# some of these may already be defined if there was a 
			# column metadata specific column in the file
		for (my $i = 0; $i < $column_count; $i++) {
			# loop for each column
			# set name unless it already has one from metadata
			unless (exists $data->{$i}) {
				$data->{$i}{'name'}  = $std_col_names->{bed6}->[$i] || 'extraColumn';
				$data->{$i}{'index'} = $i;
				$data->{$i}{'AUTO'}  = 3;
			}
			# assign the name to the column header
			$data->{'column_names'}->[$i] = 
				$data->{$i}{'name'};
		}
	}
	
	# set strand style
	if ($column_count >= 6) {
		$data->{5}{'strand_style'} = 'plusminus';
		$data->{5}{'AUTO'}++;
	}
	
	# set the feature type
	unless (defined $data->{'feature'}) {
		$data->{'feature'} = 'region';
	}
	
	# set headers flag to false
	$data->{'headers'} = 0;
}


### Subroutine to generate metadata for various UCSC gene files
sub add_ucsc_metadata {
	my ($data, $line) = @_;
	
	my @elements = split /\s+/, $line;
		# normally tab-delimited, but the specs are not explicit
	my $column_count = scalar @elements;
	$data->{'number_columns'} = $column_count; 
	$data->{'ucsc'} = $column_count;
	
	
	my $names;
	if ($column_count == 16) {
		# an extended gene prediction table, e.g. refGene, ensGene
		# as downloaded from the UCSC Table Browser or FTP site
		$names = $std_col_names->{ucsc16};
	}
	elsif ($column_count == 15) {
		# an extended gene prediction table, e.g. refGene, ensGene
		# without the bin value
		$names = $std_col_names->{ucsc15};
	}
	elsif ($column_count == 12) {
		# a known gene prediction table
		$names = $std_col_names->{ucsc12};
	}
	elsif ($column_count == 11) {
		# a refFlat gene prediction table
		$names = $std_col_names->{ucsc11};
	}
	elsif ($column_count == 10) {
		# a refFlat gene prediction table
		$names = $std_col_names->{ucsc10};
	}
	
	# assign the column names and metadata
	for (my $i = 0; $i < $column_count; $i++) {
		# loop for each column
		# set name unless it already has one from metadata
		unless (exists $data->{$i}) {
			$data->{$i}{'name'}  = $names->[$i] || 'extraColumn';
			$data->{$i}{'index'} = $i;
			$data->{$i}{'AUTO'}  = 3;
		}
		# assign the name to the column header
		$data->{'column_names'}->[$i] = 
			$data->{$i}{'name'};
	}
	
	# set strand style
	my $strand_i = find_column_index($data, 'strand');
	$data->{$strand_i}{'strand_style'} = 'plusminus';
	$data->{$strand_i}{'AUTO'}++;
	
	# set the feature type
	unless (defined $data->{'feature'}) {
		$data->{'feature'} = 'gene';
	}
	
	# set headers flag to false
	$data->{'headers'} = 0;
}


### Subroutine to generate metadata for SGR files
sub add_sgr_metadata {
	my ($data, $line) = @_;
	
	# set column metadata
	for (my $i = 0; $i < 3; $i++) {
		# loop for each column
		# set name unless it already has one from metadata
		unless (exists $data->{$i}) {
			$data->{$i}{'name'}  = $std_col_names->{sgr}->[$i] || 'extraColumn';
			$data->{$i}{'index'} = $i;
			$data->{$i}{'AUTO'}  = 3;
		}
		# assign the name to the column header
		$data->{'column_names'}->[$i] = 
			$data->{$i}{'name'};
	}
	$data->{'number_columns'} = 3; 


	# set headers flag to false
	$data->{'headers'} = 0;
	
	# set the feature type
	unless (defined $data->{'feature'}) {
		$data->{'feature'} = 'region';
	}
}


### Internal subroutine to generate metadata for standard files
sub add_standard_metadata {
	my ($data, $line) = @_;
	
	my @namelist = split /\t/, $line;
	$namelist[-1] =~ s/[\r\n]+$//;
	
	# we will define the columns based on
	for my $i (0..$#namelist) {
		# confirm that a file metadata exists for this column
		if (exists $data->{$i}) {
			unless ($namelist[$i] eq $data->{$i}->{'name'}) {
				cluck "metadata and header names for column $i do not match!";
				# set the name to match the actual column name
				$data->{$i}->{'name'} = $namelist[$i];
			}
		} 
		
		# otherwise be nice and generate it here
		else {
			$data->{$i} = {
				'name'  => $namelist[$i],
				'index' => $i,
				'AUTO'  => 3,
			};
		}
	}
	
	# check the number of columns
	if (scalar @namelist != $data->{'number_columns'} ) {
		# adjust to match actual content
		$data->{'number_columns'} = scalar @namelist;
	}
	
	# put the column names in the metadata
	push @{ $data->{'column_names'} }, @namelist;
	
	# set headers flag to true
	$data->{'headers'} = 1;
}



### Internal subroutine to generate hash of standard file format column names
sub _set_standard_column_names {
	my %names;
	
	$names{gff} = [qw(Chromosome Source Type Start Stop Score Strand Phase Group)];
	$names{bed12} = [qw(Chromosome Start End Name Score Strand  
		thickStart thickEnd itemRGB blockCount blockSizes blockStarts)];
	$names{bed6} = [qw(Chromosome Start End Name Score Strand)];
	$names{bdg} = [qw(Chromosome Start End Score)];
	$names{narrowpeak} = [qw(Chromosome Start End Name Score Strand signalValue 
		pValue qValue peak)];
	$names{broadpeak} = [qw(Chromosome Start End Name Score Strand signalValue 
		pValue qValue)];
	$names{sgr} = [qw(Chromo Start Score)];
	$names{ucsc16} = [qw(bin name chrom strand txStart txEnd cdsStart cdsEnd 
		exonCount exonStarts exonEnds score name2 cdsStartSt cdsEndStat exonFrames)];
	$names{ucsc15} = [qw(name chrom strand txStart txEnd cdsStart cdsEnd 
		exonCount exonStarts exonEnds score name2 cdsStartSt cdsEndStat exonFrames)];
	$names{ucsc12} = [qw(name chrom strand txStart txEnd cdsStart cdsEnd 
		exonCount exonStarts exonEnds proteinID alignID)];
	$names{ucsc11} = [qw(geneName transcriptName chrom strand txStart txEnd cdsStart 
		cdsEnd exonCount exonStarts exonEnds)];
	$names{ucsc10} = [qw(name chrom strand txStart txEnd cdsStart cdsEnd exonCount 
		exonStarts exonEnds)];
	$names{refflat} = $names{ucsc10};
	$names{genepred} = $names{ucsc11};
	$names{genepredext} = $names{ucsc15};
	$names{knowngene} = $names{ucsc12};
	
	return \%names;
}


__END__

=head1 NAME

Bio::ToolBox::file_helper

=head1 DESCRIPTION

These are subroutines for providing file IO for the L<Bio::ToolBox::Data> 
data structure. In other words, they are not object methods, but rather 
exportable subroutines for reading and writing data into the complex data 
structure underlying the L<Bio::ToolBox::Data> object. These subroutines 
and data structures predate the L<Bio::ToolBox::Data> object model and 
exist for backwards compatibility. End-users are strongly encouraged to 
use the L<Bio::ToolBox::Data> API. 

These file IO methods work with any generic tab-delimited text file 
of rows and columns. It also properly handles comment, metadata, and 
column-specific metadata custom to Bio::ToolBox programs.
Special file formats used in bioinformatics, including for example
GFF and BED files, are automatically recognized by their file extension and 
appropriate metadata added. 

Files opened using these subroutines are stored in a specific complex data 
structure described below. This format allows for data access as well as 
records metadata about each column (dataset) and the file in general. This
metadata helps preserve a "history" of the dataset: where it came from, how
it was collected, and how it was processed.

Additional subroutines are also present for general processing and output of
this data structure.

The data file format is described below, and following that a 
description of the data structure.

=head1 FORMAT OF BIOTOOLBOX DATA TEXT FILE

The BioToolBox data file format is not indicated by a special file extension. 
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
  
  use Bio::ToolBox::db_helper qw(load_data_file write_data_file);
  
The specific usage for each subroutine is detailed below.

=over

=item load_data_file()

This is a newer, updated file loader and parser for BioToolBox data files. It will
completely parse and load the file contents into the described data structure 
in memory. Files with metadata lines (described in BioToolBox data format) will 
have the metadata lines loaded. Files without metadata lines will have basic 
metadata (column name and index) automatically generated. The first 
non-header line should contain the column (dataset) name. Recognized file 
formats without headers, including GFF, BED, and SGR, will have the columns 
automatically named.

This subroutine uses the open_data_file() subroutine and completes the 
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
	my $data_ref = load_data_file($filename);
	
=item open_data_file()

This is a file opener and metadata parser for data files, including BioToolBox 
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
is an L<IO::Handle> object and may be manipulated as such.
Failure to read or parse the file will return an empty value.

Example:
	
	my $filename = 'my_data.txt.gz';
	my ($fh, $metadata_ref) = open_data_file($filename);
	while (my $line = $fh->getline) {
		...
	}
	$fh->close;



=item write_data_file()

This subroutine will write out a data file formatted for BioToolBox data files. 
Please refer to L<FORMAT OF BIOTOOLBOX DATA TEXT FILE> for more 
information regarding the file format. If the 'gff' key is true in the data 
hash, then a gff file will be written.

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  data     => A scalar reference to the data structure ad described
              in L<Bio::ToolBox::data_helper>. 
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
	my $data_ref = load_data_file($filename);
	...
	my $success_write = write_data_file(
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
	

=item write_summary_data()

This subroutine will summarize the data in a data file, generating mean values
for all the values in each dataset (column), and writing an output file with
the summarized data. This is useful for data collected in windows across a 
feature, for example, microarray data values across the body of genes, and 
then generating a composite or average gene occupancy.

The output file is a BioToolBox data tab-delimited file as described above with three
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

	my $main_data_ref = load_data_file($filename);
	...
	my $summary_success = write_summary_data(
		'data'         => $main_data_ref,
		'filename'     => $outfile,
		'startcolumn'  => 4,
	);


=item check_file

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
it under the terms of the Artistic License 2.0.  
