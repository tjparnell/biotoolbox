package Bio::ToolBox::Data::file;
our $VERSION = '1.62';

=head1 NAME

Bio::ToolBox::Data::file - File functions to Bio:ToolBox::Data family

=head1 DESCRIPTION

File methods for reading and writing data files for both L<Bio::ToolBox::Data> 
and L<Bio::ToolBox::Data::Stream> objects. This module should not be used 
directly. See the respective modules for more information.

=cut

use strict;
use Carp qw(carp cluck croak confess);
use File::Basename qw(fileparse);
use File::Which;
use IO::File;

# List of acceptable filename extensions
our $SUFFIX = qr/\.(?:txt|gff3?|gtf|bed|bdg|bedgraph|sgr|kgg|cdt|vcf|narrowpeak|broadpeak|reff?lat|genepred|ucsc|maf)(?:\.gz)?/i;

# gzip application
our $gzip_app;
our $bgzip_app;

### The True Statement
1; 



### Load new version data table from file

sub load_file {
	my $self = shift;
	my $file = shift;
	my $noheader = shift || 0; # may not be present
	
	# check that we have an empty table
	if ($self->last_row != 0 or $self->number_columns != 0 or $self->filename) {
		carp "Cannot load file onto an existing data table!";
		return;
	}
	
	# open the file and load metadata
	my $filename = $self->check_file($file);
	unless ($filename) {
		print " file '$file' does not exist!\n";
		return;
	}
	$self->add_file_metadata($filename);
	$self->open_to_read_fh or return;
	$self->parse_headers($noheader);
	
	# Load the data table
	while (my $line = $self->{fh}->getline) {		
		# the current file position should be at the beginning of the
		# data table information
		next if $line !~ m/\w+/;
		
		# skip comment and empty lines
		if (substr($line,0,1) eq '#') {
			$self->add_comment($line);
			next;
		}
		
		# process the line
		$self->add_data_line($line);
	}
	
	# completed loading the file
	$self->close_fh;
	delete $self->{fh};
	
	# verify the structure
	return unless $self->verify;
	
	# finished
	return 1;
}

sub taste_file {
	my $self = shift;
	my $file = shift;
	my $filename = $self->check_file($file) or return;
	my $Taste = $self->new;
	$Taste->add_file_metadata($filename);
	$Taste->open_to_read_fh or return;
	$Taste->parse_headers;
	
	# check existing metadata
	if ($Taste->gff) {
		return 'gff';
	}
	elsif ($Taste->bed == 12) {
		return 'bed';
	}
	elsif ($Taste->ucsc) {
		return 'ucsc';
	}
	
	# load first 10 data lines
	$Taste->{data_table}->[0] = $Taste->{'column_names'}; # set table header names
	for (my $i = 1; $i <= 10; $i++) {
		my $line = $Taste->fh->getline;
		next if $line !~ m/\w+/;
		next if $line =~ /^#/;
		$Taste->add_data_line($line);
	}
	$Taste->close_fh;
	
	# check if the number of columns match a known format, then verify
	my $number = $Taste->number_columns;
	if ($number == 9) {
		# possibly a GFF file
		$Taste->gff(2);
		$Taste->verify(1);
		return 'gff' if $Taste->gff == 2;
	}
	elsif ($number == 10) {
		# possibly a genePred file
		$Taste->ucsc(10);
		$Taste->verify(1);
		return 'ucsc' if $Taste->ucsc == 10;
		$Taste->add_ucsc_metadata(10,1); # force metadata
		$Taste->verify(1);
		return 'ucsc' if $Taste->ucsc == 10;
	}
	elsif ($number == 11) {
		# possibly a refFlat file
		$Taste->ucsc(11);
		$Taste->verify(1);
		return 'ucsc' if $Taste->ucsc == 11;
		$Taste->add_ucsc_metadata(11,1); # force metadata
		$Taste->verify(1);
		return 'ucsc' if $Taste->ucsc == 11;
	}
	elsif ($number == 12) {
		# possibly a knownGene or BED12 file
		$Taste->ucsc(12);
		$Taste->verify(1);
		return 'ucsc' if $Taste->ucsc == 12;
		$Taste->bed(12);
		$Taste->verify(1);
		return 'bed' if $Taste->bed == 12;
		$Taste->add_ucsc_metadata(12,1); # force metadata
		$Taste->verify(1);
		return 'ucsc' if $Taste->ucsc == 12;
		$Taste->add_bed_metadata(12,1); # force metadata
		$Taste->verify(1);
		return 'bed' if $Taste->bed == 12;
	}
	elsif ($number == 15) {
		# possibly a genePredExt file
		$Taste->ucsc(15);
		$Taste->verify(1);
		return 'ucsc' if $Taste->ucsc == 15;
		$Taste->add_ucsc_metadata(15,1); # force metadata
		$Taste->verify(1);
		return 'ucsc' if $Taste->ucsc == 15;
	}
	elsif ($number == 16) {
		# possibly a genePredExt file
		$Taste->ucsc(16);
		$Taste->verify(1);
		return 'ucsc' if $Taste->ucsc == 16; 
		$Taste->add_ucsc_metadata(16,1); # force metadata
		$Taste->verify(1);
		return 'ucsc' if $Taste->ucsc == 16;
	}
	return;
}

sub sample_gff_type_list {
	my ($self, $file) = @_;
	return unless $file =~ /\.g[tf]f\d?(?:\.gz)?$/i; # assume extension is accurate
	my $fh = $self->open_to_read_fh($file) or die "can't open $file!\n";
	my %types;
	my $count = 0; 
	while ($count < 1000) {
		my $line = $fh->getline or last;
		next if $line !~ m/\w+/;
		next if $line =~ /^#/;
		my @fields = split('\t', $line);
		$types{ $fields[2] } += 1;
		$count++;
	}
	$fh->close;
	return join(',', keys %types);
}


sub parse_headers {
	my $self = shift;
	my $noheader = shift || 0; # boolean to indicate no headers are present
	
	# filehandle
	my $fh;
	if (exists $self->{fh} and $self->{fh}) {
		$fh = $self->{fh};
	}
	elsif ($self->filename) {
		$fh = $self->open_to_read_fh($self->filename);
	}
	else {
		confess " file metadata and/or open filehandle must be set before parsing headers!";
	}
	
	# check that we have an empty table
	if ($self->last_row != 0 or $self->number_columns != 0) {
		cluck "Cannot parse file headers onto an existing data table!";
		return;
	}	
	
	# read and parse the file
	# we will ONLY parse the header lines prefixed with a #, as well as the 
	# first data row which contains the column names
	$self->program(undef); # reset this to blank, it will be filled by file metadata
	my $header_line_count = 0;
	while (my $line = $fh->getline) {		
		# we are not chomping the line here because of possible side effects
		# with UCSC tables where we have to count elements in the first line
		# and potential header lines, and the first line has a null value at 
		# the end
		
		# check for Mac-style return line endings
		if ($line =~ /\r/ and $line !~ /\n/) {
			my $filename = $self->filename;
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
			my $p = $1;
			$p =~ s/[\r\n]+$//;
			$self->program($p);
			$header_line_count++;
		}
		
		# the source database
		elsif ($line =~ m/^# Database (.+)$/) {
			my $d = $1;
			$d =~ s/[\r\n]+$//;
			$self->database($d);
			$header_line_count++;
		}
		
		# the type of feature in this datafile
		elsif ($line =~ m/^# Feature (.+)$/) {
			my $f = $1;
			$f =~ s/[\r\n]+$//;
			$self->feature($f);
			$header_line_count++;
		}
		
		# column or dataset specific information
		elsif ($line =~ m/^# Column_(\d+)/) {
			# the column number will become the index number
			my $index = $1; 
			$self->add_column_metadata($line, $index);
			$header_line_count++;
		}
		
		# gff version header
		elsif ($line =~ /^##gff-version\s+([\d\.]+)$/) {
			# store the gff version in the hash
			# this may or may not be present in the gff file, but want to keep
			# it if it is
			my $g = $1;
			$g =~ s/[\r\n]+$//;
			$self->gff($g);
			$header_line_count++;
		}
		
		# VCF version header
		elsif ($line =~ /^##fileformat=VCFv([\d\.]+)$/) {
			# store the VCF version in the hash
			# this may or may not be present in the gff file, but want to keep
			# it if it is
			my $v = $1;
			$v =~ s/[\r\n]+$//;
			$self->vcf($v);
			$self->add_comment($line); # so that it is written properly
			$header_line_count++;
		}
		
		# any other nonstandard header
		elsif ($line =~ /^#/) {
			$self->add_comment($line);
			$header_line_count++;
		}
		
		# a track or browser line 
		elsif ($line =~ /^(?:track|browser)\s+/i) {
			# common with wig, bed, or bedgraph files for use with UCSC genome browser
			# treat as a comment line, there's not that much useful info here
			$self->add_comment($line);
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
			if ( $self->_commented_header_line($line) ) {
				# these will have one comment line marked with #
				# that really contains the column headers
				
				# process the real header line
				$self->add_standard_metadata( pop @{ $self->{'comments'} } );
			}
			# we will continue here in case the commented header line was part of a 
			# formatted file type, which will be checked by extension below
			
			### a GFF file
			if ($self->extension =~ /g[tf]f/i) {
				$self->add_gff_metadata;
			}
			
			### a Bed or BedGraph file
			elsif ($self->extension =~ /bdg|bed/i) {
				my $count = scalar(split '\t', $line);
				$self->add_bed_metadata($count);
			}
			
			### a peak file
			elsif ($self->extension =~ /peak/i) {
				my $count = scalar(split '\t', $line);
				$self->add_peak_metadata($count);
			}
			
			### a UCSC gene table
			elsif ($self->extension =~ /ref+lat|genepred|ucsc/i) {
				my $count = scalar(split '\t', $line);
				$self->add_ucsc_metadata($count);
			}
			
			### a SGR file
			elsif ($self->extension =~ /sgr/i) {
				$self->add_sgr_metadata;
			}
			
			### standard text file with headers, i.e. everything else
			unless ($self->number_columns) {
				# we have not yet parsed the row of data column names
				# we will do so now
				$self->add_standard_metadata($line);
				
				# count as a header line
				$header_line_count++;
			}
		}
		last if $self->number_columns != 0;
	}
	
	# if we didn't find columns, check that it wasn't actually commented
	# for example, an empty vcf file
	if ($self->number_columns == 0 and 
		scalar(@{ $self->{comments} }) and 
		$self->{comments}->[-1] =~ /\t/
	) {
		# process the real header line
		$self->add_standard_metadata( pop @{ $self->{'comments'} } );
	}
		
	
	# No header was requested
	if ($noheader and not $self->bed and not $self->gff and 
		not $self->vcf and not $self->ucsc
	) {
		# which means that what we used as a header is actually the first data row
		
		# fix the column names
		for (my $i = 0; $i < $self->number_columns; $i++) {
			my $name = $self->name($i);
			$self->name($i, "Column$i ($name)");
			$self->{$i}{'AUTO'} = 3;
		}
		# adjust metadata
		$header_line_count -= 1;
		$self->{'headers'} = -1; # special case, we never write headers here 
	}
	
	# Header sanity check
		# some exported file formats, such as from R, do not include a proper 
		# header for the first column, as these are assumed to be row names
		# this will result in incorrectly parsed files where the last columns 
		# will be merged into one column with an internal tab - not good
		# need to handle these
	my $nextline = $fh->getline;
	if ($nextline) {
		my @nextdata = split '\t', $nextline;
		if (scalar(@nextdata) - 1 == $self->number_columns) {
			# whoops! we caught a off-by-one discrepancy between header and data row
			my $old_last = $self->last_column;
			$fh->close; # having this open complicates changing columns....
			
			# add a new "column" (just metadata for now) and move it to the beginning
			my $i = $self->add_column('Column1');
			$self->reorder_column($i, 0 .. $old_last);
		}
	}
	
	# re-open the file
		# I tried using seek functions - but they don't work with binary gzip 
		# files, and I can't get the seek function to return the same position
		# as simply advancing through the file like below
		# so I'll just do it the old way and close/open and advance
	$fh->close if $fh; # may have been closed from the sanity check above
	$fh = $self->open_to_read_fh;
	for (1 .. $header_line_count) {
		my $line = $fh->getline;
	}
	$self->{header_line_count} = $header_line_count;
	$self->{fh} = $fh;
	return 1;
}


sub add_data_line {
	my ($self, $line) = @_;
	
	# do not chomp the line yet, just split into an array
	my @linedata = split '\t', $line, $self->{number_columns};
	
	# chomp the last element
	# we do this here to ensure the tab split above gets all of the values
	# otherwise trailing null values aren't included in @linedata
	# be sure to handle both newlines and carriage returns
	$linedata[-1] =~ s/[\r\n]+$//;
	
	# check for extra remaining tabs
	if (index($linedata[-1], "\t") != -1) {
		die "FILE INCONSISTENCY ERRORS! line has additional tabs (columns) than headers!\n $line";
	}
		
	# add the line of data
	push @{ $self->{data_table} }, \@linedata;
	$self->{last_row} += 1;
	return 1;
}


### Parse the filename using the list suffix list
sub add_file_metadata {
	my ($self, $filename) = @_;
	confess "no valid filename!" unless defined $filename;
	my ($basename, $path, $extension) = fileparse($filename, $SUFFIX);
	unless ($extension) {
		# look for a nonstandard extension, allowing for .gz extension
		if ($filename =~ /(\.\w+(?:\.gz)?)$/i) {
			$extension = $1;
			$basename =~ s/$extension\Z//;
		}
	}
	$self->{filename}  = $filename;
	$self->{basename}  = $basename;
	$self->{path}      = $path;
	$self->{extension} = $extension;
}


### Write out a data file
sub write_file {
	my $self = shift;
	
	# collect passed arguments
	my %args;
	if (scalar(@_) == 1) {
		$args{'filename'} = $_[0];
	}
	else {
		%args = @_;
	}
	$args{'filename'} ||= $args{'file'} || undef;
	$args{'format'}   ||= undef;
	unless (exists $args{'gz'}) {$args{'gz'} = undef} 
	
	# check the data
	unless ($self->verify) {
		cluck "bad data structure!";
		return;
	}
	
	# check filename
	unless ($args{'filename'} or $self->filename) {
		cluck "no filename given!\n";
		return;
	}
	
	# split filename into its base components
	my ($name, $path, $extension) = 
		fileparse($args{'filename'} || $self->filename, $SUFFIX);
	
	# Adjust filename extension if necessary
	if ($extension =~ /(g[tf]f)/i) {
		if (not $self->gff) {
			# let's set it to true and see if it passes verification
			$self->{'gff'} = $extension =~ /gtf/i ? 2.5 : 3; # default
			unless ($self->verify and $self->gff) {
				warn " GFF structure changed, re-setting extension from $extension to .txt\n";
				$extension =~ s/g[tf]f3?/txt/i;
			}
		}
	}
	elsif ($extension =~ /bedgraph|bed|bdg|narrowpeak|broadpeak/i) {
		if (not $self->bed) {
			# let's set it to true and see if it passes verification
			$self->{'bed'} = 1; # a fake true
			unless ($self->verify and $self->bed) {
				warn " BED structure changed, re-setting extension from $extension to .txt\n";
				$extension = $extension =~ /gz$/i ? '.txt.gz' : '.txt';
			}
		}
	}
	elsif ($extension =~ /vcf/i) {
		if (not $self->vcf) {
			# let's set it to true and see if it passes verification
			$self->{'vcf'} = 1; # a fake true
			unless ($self->verify and $self->vcf) {
				warn " VCF structure changed, re-setting extension from $extension to .txt\n";
				$extension = $extension =~ /gz$/i ? '.txt.gz' : '.txt';
			}
		}
	}
	elsif ($extension =~ /sgr/i) {
		unless ($self->{'extension'} =~ /sgr/i) {
			# original file was not SGR
			# let's pretend it was and see if still passes 
			# the sgr verification relies on the recorded extension
			$self->{'extension'} = '.sgr';
			$self->verify;
			if ($self->extension =~ /txt/) {
				warn " SGR structure changed, re-setting extension from $extension to .txt\n";
			}
			$extension = $self->{'extension'};
		}
	}
	elsif ($extension =~ /reff?lat|genepred|ucsc/i) {
		if ($self->ucsc != $self->number_columns) {
			# it's not set as a ucsc data
			# let's set it to true and see if it passes verification
			$self->ucsc($self->number_columns);
			unless ($self->verify and $self->ucsc) {
				warn " UCSC structure changed, re-setting extension from $extension to .txt\n";
				$extension = $extension =~ /gz$/i ? '.txt.gz' : '.txt';
			}
		}
	}
	elsif ($extension =~ /txt/i) {
		# plain old text file, sounds good to me
		# make sure headers are enabled
		$self->{'headers'} = 1;
	}
	elsif (not $extension) {
		# no extension was available
		# try and determine one from metadata
			
		if ($self->gff) {
			$extension = $self->gff == 3 ? '.gff3' : $self->gff == 2.5 ? '.gtf' : '.gff';
		} 
		elsif ($self->bed) {
			if (
				$self->number_columns == 4 and 
				$self->name(3) =~ /score/i
			) {
				$extension = '.bdg'; # a bedGraph file
			}
			else {
				$extension = '.bed'; # a regular bed file
			}
		}
		elsif ($self->ucsc) {
			# use a generic ucsc format, don't bother to customize it
			$extension = '.ucsc';
		}
		elsif ($self->vcf) {
			$extension = '.vcf';
		}
		elsif ($name =~ /(\.\w{3}(?:\.gz)?)$/i) {
			# a non-standard 3 letter file extension
			# anything else might be construed as part of the filename, so run the 
			# risk of adding a default extension below
			$extension = $1;
			$name =~ s/$extension\Z//;
		}
		elsif ($self->extension) {
			# original file had an extension, re-use it if appropriate
			# why wouldn't this get picked up above???? probably old cruft, 
			# or a non-standard or unknown file extension
			# leave it in for the time being, shouldn't hurt anything
			if ($self->extension =~ /g[tf]f/i) {
				$extension = $self->gff ? $self->extension : '.txt';
			}
			elsif ($self->extension =~ /bed|bdg|peak/i) {
				$extension = $self->bed ? $self->extension : '.txt';
			}
			else {
				# an unstructured format
				$extension = $self->extension;
			}
		}
		else {
			# normal data text file
			$extension = '.txt';
		}
	}
	# otherwise the extension must be good, hope for the best
	
	# determine format 
	# this is an arcane specification of whether we want a "simple" no metadata 
	# format, or an ordinary text format that may or may not have metadata
	# it's currently not hurting much, so leave it in for now?
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
		if ($extension =~ m/\.vcf\.gz/i) {
			# vcf requires bgzip
			$args{'gz'} = 2;
		}
		elsif ($extension =~ m/\.gz$/i) {
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
	# don't even get me started on Windows NTFS path length limitation
		if (length($name . $extension) > 255) {
		my $limit = 253 - length($extension);
		$name = substr($name, 0, $limit) . '..';
		warn " filename too long! Truncating to $limit characters\n";
	}
	
	# generate the new filename
	my $newname = $path . $name . $extension;
	
	
	# Convert strand information
	my $strand_i = $self->strand_column;
	if (defined $strand_i and ($self->gff or $self->bed or $self->ucsc) ) {
		# convert to +/-/. nomenclature as necessary
		if ($self->gff) {
			for my $row (1 .. $self->last_row) {
				my $s = $self->{'data_table'}->[$row][$strand_i];
				if ($s =~ /\d/) {
					$s = $s == 1 ? '+' : $s == -1 ? '-' : '.';
				}
				$self->{'data_table'}->[$row][$strand_i] = $s;
			}
		}
		elsif ($self->bed or $self->ucsc) {
			for my $row (1 .. $self->last_row) {
				my $s = $self->{'data_table'}->[$row][$strand_i];
				if ($s =~ /\d/) {
					$s = $s >= 0 ? '+' : '-';
				}
				$self->{'data_table'}->[$row][$strand_i] = $s;
			}
		}
	}
	
	
	# Open file for writing
	my $fh = $self->open_to_write_fh($newname, $args{'gz'});
	return unless defined $fh;
	
	
	# Write the headers
	if ($args{'format'} eq 'text') {
		# default text format has metadata headers
		
		# write gff statement if gff format
		if ($self->gff) {
			$fh->print('##gff-version ' . $self->gff . "\n");
		}
		
		# Write the primary headers
		unless (
			$self->gff or $self->bed or $self->ucsc or $self->vcf or
			$extension =~ m/sgr|kgg|cdt|peak/i
		) {
			# we only write these for normal text files, not defined format files
			
			if ($self->program) {
				$fh->print('# Program ' . $self->program . "\n");
			}
			if ($self->database) {
				$fh->print('# Database ' . $self->database . "\n");
			}
			if ($self->feature) {
				$fh->print('# Feature ' . $self->feature . "\n");
			}
		}
		
		# Write the miscellaneous headers
		foreach ( @{ $self->{'comments'} } ) {
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
		for (my $i = 0; $i < $self->number_columns; $i++) {
			# each column metadata in the hash is referenced by the column's
			# index number as the key
			# we will take each index one at a time in increasing order
			
			# some files do not need or tolerate metadata lines, for those 
			# known files the metadata lines will be skipped
			
			# these column metadata lines do not need to be written if they
			# only have two values, presumably name and index, for files 
			# that don't normally have column headers, e.g. gff
			if (
				exists $self->{$i}{'AUTO'} and
				scalar( keys %{ $self->{$i} } ) == $self->{$i}{'AUTO'}
			) {
				# some of the metadata values were autogenerated and 
				# we have the same number of keys as were autogenerated
				# no need to write these
				next;
			}
			elsif (scalar( keys %{ $self->{$i} } ) == 2) {
				# only two metadata keys exist, name and index
				# these are so simple it's not worth writing them
				next;
			}
			elsif ($extension =~ /sgr|kgg|cdt/i or $self->ucsc or $self->vcf) {
				# these do not support metadata lines
				next;
			}
			
			# we will put each key=value pair into @pairs, listed asciibetically
			my @pairs; # an array of the key value pairs from the metadata hash
			# put name first
			# we are no longer writing the index number
			push @pairs, 'name=' . $self->{$i}{'name'};
			# put remainder in alphabetical order
			foreach (sort {$a cmp $b} keys %{ $self->{$i} } ) {
				next if $_ eq 'name'; # already written
				next if $_ eq 'index'; # internal use only
				next if $_ eq 'AUTO'; # internal use only
				push @pairs,  $_ . '=' . $self->{$i}{$_};
			}
			
			# Finally write the header line, joining the pairs with a 
			# semi-colon into a single string.
			# The column identifier is comprised of the word 'Column' 
			# and the index number joined by '_'.
			$fh->print("# Column_$i ", join(";", @pairs), "\n");
		}
	}
	
	
	# Write the table column headers
	if ($self->{'headers'} == 1) {
		$fh->printf("%s\n", join("\t", @{ $self->{'data_table'}[0] }));
	}
		
	
	# Write the data table
	if ($args{'format'} eq 'simple') {
		
		# the simple format will strip the non-value '.' from the table
		for (my $i = 1; $i <= $self->last_row; $i++) {
			# we will step though the data_table array one row at a time
			# convert the non-value '.' to undefined
			# and print using a tab-delimited format
			my @linedata;
			foreach ( @{ $self->{'data_table'}[$i] }) {
				if ($_ eq '.') {
					push @linedata, undef;
				} else {
					push @linedata, $_;
				}
			}
			$fh->printf("%s\n", join("\t", @linedata));
		}
	}
	
	else {
		# normal data files
		for (my $i = 1; $i <= $self->last_row; $i++) {
			# we will step though the data_table array one row at a time
			# we will join each row's array of elements into a string to print
			# using a tab-delimited format
			$fh->printf("%s\n", join("\t", @{ $self->{'data_table'}[$i] }));
		}
	}
	
	# done writing
	$fh->close;
	
	# if we made it this far, it should've been a success!
	# return the new file name as indication of success
	return $newname;
}

sub save {
	return shift->write_file(@_);
}


#### Open a file for reading
sub open_to_read_fh {
	my $self = shift;
	my $file = shift || undef;
	my $obj = ref($self) =~ /^Bio::ToolBox/ ? 1 : 0;
	
	# check file
	if (not $file and $obj) {
		$file = $self->filename || undef;
	}
	unless ($file) {
		carp " no filename provided or associated with object!";
		return;
	}
	
	# Open filehandle object as appropriate
	my $fh; 
	if ($file =~ /\.gz$/i) {
		# the file is compressed with gzip
		$fh = IO::File->new("gzip -dc $file |") or 
			carp "unable to read '$file' $!\n";
	} 
	elsif ($file =~ /\.bz2$/i) {
		# the file is compressed with bzip2
		$fh = IO::File->new("bzip2 -dc $file |") or 
			carp "unable to read '$file' $!\n";
	} 
	else {
		# the file is uncompressed and space hogging
		$fh = IO::File->new($file, 'r') or 
			carp "unable to read '$file' $!\n";
	}
	
	if ($obj) {
		$self->{fh} = $fh;
	}
	return $fh;	
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
		if ($filename =~ m/\.vcf(\.gz)?$/i) {
			$gz = 2; # bgzip
		}
		if ($filename =~ m/\.gz$/i) {
			$gz = 1; # regular gzip
		}
		else {
			$gz = 0; # default
		}
	}
	
	# gzip compression application
	if ($gz == 1 and not $gzip_app) {
		# use parallel gzip if possible
		# this is stored in a global variable so we only have to look once
		$gzip_app = which('pigz');
		if ($gzip_app) {
			$gzip_app .= ' -p 3'; # use a conservative 3 processes, plus perl, so 4 total
		}
		else {
			# default is the standard gzip application
			# should be available in any application
			$gzip_app = which('gzip');
		}
		unless ($gzip_app) {
			croak "no gzip application in PATH to open compressed file handle output!\n";
		}
	}
	elsif ($gz == 2 and not $bgzip_app) {
		# use parallel bgzip if possible
		# this is stored in a global variable so we only have to look once
		$bgzip_app = which('bgzip');
		if ($bgzip_app) {
			# I'm going to assume this is a recent bgzip with multi-threading
			# use 3 threads, same as with pigz
			$bgzip_app .= ' -@ 3 -c';
		}
		unless ($bgzip_app) {
			croak "no bgzip application in PATH to open compressed file handle output!\n";
		}
	}
	my $gzipper = $gz == 1 ? $gzip_app : $gz == 2 ? $bgzip_app : undef;
	
	# check file append mode
	unless (defined $append) {
		# default is not to append
		$append = 0;
	}
	
	
	# Generate appropriate filehandle object
	my $fh;
	if (not $gzipper and not $append) {
		$fh = IO::File->new($filename, 'w') or 
			carp "cannot write to file '$filename' $!\n";
	}
	elsif ($gzipper and !$append) {
		$fh = IO::File->new("| $gzipper >$filename") or 
			carp "cannot write to compressed file '$filename' $!\n";
	}
	elsif (not $gzipper and $append) {
		$fh = IO::File->new(">> $filename") or 
			carp "cannot append to file '$filename' $!\n";
	}
	elsif ($gzipper and $append) {
		$fh = IO::File->new("| $gzipper >>$filename") or 
			carp "cannot append to compressed file '$filename' $!\n";
	}
	return $fh if defined $fh;
}


### Subroutine to check for file existance
sub check_file {
	my ($self, $filename) = @_;
	
	# check for file existance
	if (-e $filename) {
		# confirmed full filename and path
		return $filename;
	}
	else {
		# file name is either incomplete or non-existent
		# try adding some common file extensions in case those are missing
		my $new_filename;
		foreach my $ext (qw(.gz .txt .txt.gz .bed .bed.gz)) {
			if (-e $filename . $ext) {
				$new_filename = $filename . $ext;
				last;
			}
		}
		return $new_filename;
	}
}


sub fh {
	my $self = shift;
	return $self->{fh} if exists $self->{fh};
	return;
}

sub close_fh {
	my $self = shift;
	$self->{fh}->close if (exists $self->{fh} and $self->{fh});
}


### Internal subroutine to check if a comment line contains headers
sub _commented_header_line {
	my ($data, $line) = @_;
	
	# prepare arrays from the other lines and current line
	my @commentdata;
	if ( scalar @{ $data->{'comments'} } >= 1 ) {
		# take the last line in the other array
		@commentdata = split '\t', $data->{'comments'}->[-1];
	}
	my @linedata = split '\t', $line;
	
	# check if the counts are equal
	if (scalar @commentdata == scalar @linedata) {
		return 1;
	}
	else {
		return 0;
	}
}


### Internal subroutine to process metadata for standard columns
sub add_column_metadata {
	my ($data, $line, $index) = @_;
	
	# strip the Column metadata identifier
	$line =~ s/[\r\n]+$//;
	$line =~ s/^# Column_\d+ //; 
	
	# break up the column metadata
	my %temphash; # a temporary hash to put the column metadata into
	foreach (split ';', $line) {
		my ($key, $value) = split '=';
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
	}
	return 1;
}

### Subroutine to generate metadata for gff files
	# gff files have nine defined columns
	# there are different specifications and variants:
	# gff (v.1), gff v.2, gff v.2.5 (aka gtf), gff v.3 (gff3)
	# however, the columns are the same in all versions
	# more info on gff can be found http://gmod.org/wiki/GFF3
sub add_gff_metadata {
	my $self = shift;
	my $version = shift || undef;
	my $force = shift || 0;
	
	# set the gff version based on the extension if it isn't already
	unless ($self->gff) {
		$self->{gff} = defined $version ? $version :
			$self->extension =~ /gtf/  ? 2.5 :
			$self->extension =~ /gff3/ ? 3   :
			2;
	}
	
	# set the metadata for the each column
		# some of these may already be defined if there was a 
		# column metadata specific column in the file
	my $column_names = $self->standard_column_names('gff');
	for (my $i = 0; $i < 9; $i++) {
		# loop for each column
		# set metadata unless it's already loaded
		if ($force or not exists $self->{$i}) {
			$self->{$i}{'name'}  = $column_names->[$i];
			$self->{$i}{'index'} = $i;
			$self->{$i}{'AUTO'}  = 3;
		}
		# assign the name to the column header
		if ($force or not defined $self->{'column_names'}->[$i]) {
			$self->{'column_names'}->[$i] = $self->{$i}{'name'};
		}
	}
	$self->{data_table}->[0] = $self->{'column_names'};
	
	# set column number always to 9
	if ($force or $self->{'number_columns'} == 0) { 
		$self->{'number_columns'} = 9;
	}
	
	# set headers flag to false
	$self->{'headers'} = 0 unless $self->{0}{'name'} =~ /^#/;
	
	# set the feature type
	unless (defined $self->{'feature'}) {
		$self->{'feature'} = 'region';
	}
	return 1;
}


### Subroutine to generate metadata for BED files
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
sub add_bed_metadata {
	my ($self, $column_count) = @_;
	my $force = shift || 0;

	$self->{'number_columns'} = $column_count; 
	$self->{'bed'} = $column_count;
	
	# set column names
	my $bed_names = $self->standard_column_names(
		# try to get this from the file extension if available
		$self->extension =~ /bdg|graph/i ? 'bdg' : 'bed12'
	);
	
	# set the metadata for each column
		# some of these may already be defined if there was a 
		# column metadata specific column in the file
	for (my $i = 0; $i < $column_count; $i++) {
		# loop for each column
		# set name unless it already has one from metadata
		if ($force or not exists $self->{$i}) {
			$self->{$i}{'name'}  = $bed_names->[$i] || 'extraColumn';
			$self->{$i}{'index'} = $i;
			$self->{$i}{'AUTO'}  = 3;
		}
		# assign the name to the column header
		if ($force or not defined $self->{'column_names'}->[$i]) {
			$self->{'column_names'}->[$i] = $self->{$i}{'name'};
		}
	}
	$self->{data_table}->[0] = $self->{'column_names'};

	# set the feature type
	unless (defined $self->{'feature'}) {
		$self->{'feature'} = 'region';
	}
	
	# set headers flag to false
	$self->{'headers'} = 0 unless $self->{0}{'name'} =~ /^#/;
	return 1;
}


### Subroutine to generate metadata for broadpeak and narrowpeak files
	# three different types of peak files are available
	# see http://genome.ucsc.edu/FAQ/FAQformat.html
sub add_peak_metadata {
	my ($self, $column_count) = @_;
	my $force = shift || 0;

	$self->{'number_columns'} = $column_count; 
	$self->{'bed'} = $column_count;
		# technically this is called bed6+4 or bed6+3, but for our 
		# purposes here, we will stick to column count to avoid breaking stuff
	
	# column names determined by extension
	my $column_names = $self->standard_column_names(
		$self->extension =~ /narrowpeak/i ? 'narrowpeak' :
		$self->{extension} =~ /broadpeak/i ? 'broadpeak' : 'bed6'
	);
	
	# add metadata
	for (my $i = 0; $i < $column_count; $i++) {
		if ($force or not exists $self->{$i}) {
			$self->{$i}{'name'} = $column_names->[$i] || 'extraColumn';
			$self->{$i}{'index'} = $i;
			$self->{$i}{'AUTO'}  = 3;
		}
		# assign the name to the column header
		if ($force or not defined $self->{'column_names'}->[$i]) {
			$self->{'column_names'}->[$i] = $self->{$i}{'name'};
		}
	}
	$self->{data_table}->[0] = $self->{'column_names'};
	
	# set the feature type
	unless (defined $self->{'feature'}) {
		$self->{'feature'} = 'region';
	}
	
	# set headers flag to false
	$self->{'headers'} = 0 unless $self->{0}{'name'} =~ /^#/;
	return 1;
}


### Subroutine to generate metadata for various UCSC gene files
	# these are tricky, as we will try to determine contents by counting
	# the columns, which may not be accurate
	# not only that, but this presumes the extension even makes it recognizable
	# see http://genome.ucsc.edu/FAQ/FAQformat.html#format9 for details
	# also biotoolbox script ucsc_table2gff3.pl
sub add_ucsc_metadata {
	my ($self, $column_count) = @_;
	my $force = shift || 0;
	
	$self->{'number_columns'} = $column_count; 
	$self->{'ucsc'} = $column_count;
	
	# set names based on column count
	my $column_names = $self->standard_column_names('ucsc' . $column_count);
	
	# assign the column names and metadata
	for (my $i = 0; $i < $column_count; $i++) {
		# loop for each column
		# set name unless it already has one from metadata
		if ($force or not exists $self->{$i}) {
			$self->{$i}{'name'}  = $column_names->[$i] || 'extraColumn';
			$self->{$i}{'index'} = $i;
			$self->{$i}{'AUTO'}  = 3;
		}
		# assign the name to the column header
		if ($force or not defined $self->{'column_names'}->[$i]) {
			$self->{'column_names'}->[$i] = $self->{$i}{'name'}; 
		}
	}
	$self->{data_table}->[0] = $self->{'column_names'};
	
	# set the feature type
	unless (defined $self->{'feature'}) {
		$self->{'feature'} = 'gene';
	}
	
	# set headers flag to false
	$self->{'headers'} = 0 unless $self->{0}{'name'} =~ /^#/;
	return 1;
}


### Subroutine to generate metadata for SGR files
	# a sgr file contains three columns: chromo, position, score
	# this is a very simple file format, useful in exporting and
	# importing to binary BAR files used in T2, USeq, and IGB
sub add_sgr_metadata {
	my $self = shift;
	my $force = shift || 0;
	
	# set column metadata
	my $column_names = $self->standard_column_names('sgr');
	for (my $i = 0; $i < 3; $i++) {
		# loop for each column
		# set name unless it already has one from metadata
		if ($force or not exists $self->{$i}) {
			$self->{$i}{'name'}  = $column_names->[$i] || 'extraColumn';
			$self->{$i}{'index'} = $i;
			$self->{$i}{'AUTO'}  = 3;
		}
		# assign the name to the column header
		if ($force or not defined $self->{'column_names'}->[$i]) {
			$self->{'column_names'}->[$i] = $self->{$i}{'name'}; 
		}
	}
	$self->{data_table}->[0] = $self->{'column_names'};
	$self->{'number_columns'} = 3; 


	# set headers flag to false
	$self->{'headers'} = 0 unless $self->{0}{'name'} =~ /^#/;
	
	# set the feature type
	unless (defined $self->{'feature'}) {
		$self->{'feature'} = 'region';
	}
	return 1;
}


### Internal subroutine to generate metadata for standard files
sub add_standard_metadata {
	my ($self, $line) = @_;
	
	my @namelist = split '\t', $line;
	$namelist[-1] =~ s/[\r\n]+$//;
	
	# we will define the columns based on
	for my $i (0..$#namelist) {
		# make up a name if one doesn't exist
		$namelist[$i] ||= "Column_$i";
		
		# confirm that a file metadata exists for this column
		if (exists $self->{$i}) {
			unless ($namelist[$i] eq $self->{$i}->{'name'}) {
				warn "metadata and header names for column $i do not match!";
				# set the name to match the actual column name
				$self->{$i}->{'name'} = $namelist[$i];
			}
		} 
		
		# otherwise be nice and generate it here
		else {
			$self->{$i} = {
				'name'  => $namelist[$i],
				'index' => $i,
				'AUTO'  => 3,
			};
		}
	}
	
	# check the number of columns
	if (scalar @namelist != $self->{'number_columns'} ) {
		# adjust to match actual content
		$self->{'number_columns'} = scalar @namelist;
	}
	
	# put the column names in the metadata
	$self->{'column_names'} = \@namelist;
	$self->{data_table}->[0] = $self->{'column_names'};
	
	# set headers flag to true
	$self->{'headers'} = 1;
	return 1;
}


### Internal subroutine to generate hash of standard file format column names
sub standard_column_names {
	my ($self, $type) = @_;
	
	if ($type eq 'gff') {
		return [qw(Chromosome Source Type Start Stop Score Strand Phase Group)];
	}
	elsif ($type eq 'bed12') {
		return [qw(Chromosome Start0 End Name Score Strand  
			thickStart0 thickEnd itemRGB blockCount blockSizes blockStarts0)];
	}
	elsif ($type eq 'bed6') {
		return [qw(Chromosome Start0 End Name Score Strand)];
	}
	elsif ($type eq 'bdg') {
		return [qw(Chromosome Start0 End Score)];
	}
	elsif ($type eq 'narrowpeak') {
		return [qw(Chromosome Start0 End Name Score Strand signalValue 
			pValue qValue peak)];
	}
	elsif ($type eq 'broadpeak') {
		return [qw(Chromosome Start0 End Name Score Strand signalValue 
			pValue qValue)];
	}
	elsif ($type eq 'sgr') {
		return [qw(Chromo Start Score)];
	}
	elsif ($type eq 'ucsc16') {
		return [qw(bin name chrom strand txStart0 txEnd cdsStart0 cdsEnd exonCount 
			exonStarts0 exonEnds score name2 cdsStartSt cdsEndStat exonFrames)];
	}
	elsif ($type eq 'ucsc15' or $type eq 'genepredext') {
		return [qw(name chrom strand txStart0 txEnd cdsStart0 cdsEnd exonCount 
			exonStarts0 exonEnds score name2 cdsStartSt cdsEndStat exonFrames)];
	}
	elsif ($type eq 'ucsc12' or $type eq 'knowngene') {
		return [qw(name chrom strand txStart0 txEnd cdsStart0 cdsEnd exonCount 
			exonStarts0 exonEnds proteinID alignID)];
	}
	elsif ($type eq 'ucsc11' or $type eq 'refflat') {
		return [qw(geneName transcriptName chrom strand txStart0 txEnd cdsStart0 
			cdsEnd exonCount exonStarts0 exonEnds)];
	}
	elsif ($type eq 'ucsc10' or $type eq 'genepred') {
		return [qw(name chrom strand txStart0 txEnd cdsStart0 cdsEnd exonCount 
			exonStarts exonEnds)];
	}
	else {
		confess "unrecognized standard column name format '$type'!";
	}
}


__END__

=head1 DESCRIPTION

These are methods for providing file IO for the L<Bio::ToolBox::Data> 
data structure. These file IO methods work with any generic tab-delimited 
text file of rows and columns. It also properly handles comment, metadata, 
and column-specific metadata custom to L<Bio::ToolBox> programs.
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

=head1 RECOGNIZED FILE FORMATS

L<Bio::ToolBox> will recognize a number of standard bioinformatic file 
formats, almost all of which are recognized by their extension. Recognition 
is NOT guaranteed if an alternate file extension is used!!!!

These formats include

=over 4

=item BED 

These include file extensions F<.bed>, F<.bedgraph>, and F<.bdg>.
Bed files must have 3-12 columns. BedGraph files must have 4 columns.

=item GFF

These include file extensions F<.gff>, F<.gff3>, and F<.gtf>. 
The specific format may also be recognized by the C<gff-version> pragma. 
These files must have 9 columns.

=item UCSC tables

These include file extensions F<.refFlat>, F<.genePred>, and F<.ucsc>. In 
some cases, a simple F<.txt> can also be recognized if the file matches the 
expected file structure. Different formats are typically recognized by the 
number of columns, and can include simple refFlat, gene prediction, extended 
gene prediction, and known Gene tables. The Bin column may or may not be present.

=item Peak files

These include file extensions F<.narrowPeak> and F<.broadPeak>. 
These are special "BED6+4" file formats. 

=item CDT

These include file extension F<.cdt>. 
Cluster data files used with Cluster 3.0 and Treeview.

=item SGR

Rare file format of chromosome, position, score. File extension F<.sgr>.

=item TEXT

Almost any tab-delimited text file with a F<.txt> extension can be loaded.

=item Compressed files

File extension F<.gz> and F<.bz2> are recognized as compressed files.
Compressed files are usually read through an external decompression 
program. All of the above formats can be loaded as compressed files.


=back

=head1 DEFAULT BIO::TOOLBOX DATA TEXT FILE FORMAT

When not writing to a defined format, e.g. BED or GFF, a L<Bio::ToolBox> 
Data structure is written as a simple tab-delimited text file, with the 
first line being the column header names. Such files are easily parsed 
by other programs. 

If additional metadata is included in the Data object, then these are 
written as comment lines, prefixed by a "# ", before the table. Metadata 
can describe the data within the table with regards to its type, source, 
methodology, history, and processing. The metadata is designed to be read 
by both human and computer. Opening files without this metadata 
will result in basic default metadata assigned to each column. 

Some common metadata lines that are specifically recognized are listed below.

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
modifications performed on the data are also recorded here. 

A list of common column metadata keys is shown. 

=over 4

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

The strandedness of the data collected. Values include 'sense',
'antisense', or 'none'

=item method

The method of collecting values

=item log2

boolean indicating the values are in log2 space or not

=back

=back

=head1 USER METHODS REFERENCE

These methods are generally available to L<Bio::ToolBox::Data> objects 
and can be used by the user.

=over

=item load_file

This will load a file into a new, empty Data table. This function is 
called automatically when a filename is provided to the new() function. 
The existence of the file is first checked (appending common missing 
extensions as necessary), metadata and column headers processed and/or 
generated from default settings, the content loaded into the table, and 
the structure verified. Error messages may be printed if the structure or 
format is inconsistent or doesn't match the expected format, e.g a file 
with a F<.bed> extension doesn't match the UCSC specification.
Pass the name of the filename.

=item taste_file

Tastes, or checks, a file for a certain flavor, or known gene file formats. 
This is based on file extension, metadata headers, and/or file content 
in the first 10 lines or so. Returns a string based on the file format.
Values include gff, bed, ucsc, or undefined. Useful for determining if 
the file represents a known gene table format that lacks a defined file 
extension, e.g. UCSC formats. Pass the path to the file to check.

=item sample_gff_type_list

Checks the different types of features available in a GFF formatted file. 
It will temporarily open the file, read the first 1000 lines or so, and 
compile a list of the values in the 3rd column of the GFF file. It will 
return a comma-delimited string of these values upon success, suitable 
for regular expression checking. Pass the name of the GFF file to check.

=item add_file_metadata

Add or update the file metadata to a L<Bio::ToolBox::Data> object. 
This will automatically parse the path, basename, and recognized file extension.
Pass the file name.

=item write_file

=item save

This method will write out a L<Bio::ToolBox::Data> structure to file. 
Zero or more values may be passed to the method.

Pass no values, and the filename stored in the metadata will be used in 
writing the file, effectively overwriting itself. No filename will generate 
an error. 

Pass a single value representing the filename to write. The current 
working directory is assumed if no path is provided in the filename.

Pass an array of key =E<gt> values for fine control of the write process. 
Keys include the following: 

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
  gz       => A value (2, 1, or 0) indicating whether the file 
              should be written through a gzip filter to compress. If 
              this value is undefined, then the file name is checked 
              for the presence of the '.gz' extension and the value 
              set appropriately. Default is false (no compression).
              Set to 1 to use ordinary gzip, or set to 2 to use block 
              gzip (bgzip) compression for tabix compatibility.
  simple   => A boolean value (1 or 0) indicating whether a simple 
              tab-delimited text data file should be written. This is 
              an old alias for setting 'format' to 'simple'.

The method will return the real name of the file written if the write was 
successful. The filename may be modified slightly as necessary, for example 
append or change the file extension to match the specified file format.

=item open_to_read_fh

This subroutine will open a file for reading. If the passed filename has
a F<.gz> extension, it will appropriately open the file through a gunzip 
filter.

Pass the subroutine the filename. It will return a scalar reference to the
open filehandle. The filehandle is an L<IO::Handle> object and may be manipulated
as such.

Example

	my $filename = 'my_data.txt.gz';
	my $fh = Bio::ToolBox::Data::file->open_to_read_fh($filename);
	while (my $line = $fh->getline) {
		# do something
	}
	$fh->close;
	

=item open_to_write_fh

This subroutine will open a file for writing. If the passed filename has
a F<.gz> extension, it will appropriately open the file through a gzip 
filter.

Pass the subroutine three values: the filename, a boolean value indicating
whether the file should be compressed with gzip, and a boolean value 
indicating that the file should be appended. The gzip and append values are
optional. The compression status may be determined automatically by the 
presence or absence of the passed filename extension; the default is no 
compression. The default is also to write a new file and not to append.

If gzip compression is requested, but the filename does not have a F<.gz> 
extension, it will be automatically added. However, the change in file name 
is not passed back to the originating program; beware!

The subroutine will return a scalar reference to the open filehandle. The 
filehandle is an L<IO::Handle> object and may be manipulated as such.

Example

	my $filename = 'my_data.txt.gz';
	my $gz = 1; # compress output file with gzip
	my $fh = Bio::ToolBox::Data::file->open_to_write_fh($filename, $gz);
	# write to new compressed file
	$fh->print("something interesting\n");
	$fh->close;

=back

=head1 OTHER METHODS

These methods are used internally by L<Bio::ToolBox::Core> and other objects 
are not recommended for use by general users. 

=over 4

=item parse_headers

This will determine the file format, parse any metadata lines that may 
be present, add metadata and inferred column names for known file formats, 
and determine the table column header names. This is automatically called 
by L</load_file>, and generally need not be called.

Pass a true boolean option if there were no headers in the file.

=item add_data_line

Parses a text line from the file into a Data table row. Pass the text line.

=item check_file

This subroutine confirms the existance of a passed filename. If not 
immediately found, it will attempt to append common file extensions 
and verify its existence. This allows the user to pass only the base 
file name and not worry about missing the extension. This may be useful 
in shell scripts. Pass the file name.

=item add_column_metadata

Parse a column metadata line from a file into a Data structure.

=item add_gff_metadata

Add default column metadata for a GFF file. 
Specify which GFF version. A second boolean value can be passed to 
force the method.

=item add_bed_metadata

Add default column metadata for a BED file. Specify the number of BED 
columns. Pass a second boolean to force the method.

=item add_peak_metadata

Add default column metadata for a narrowPeak or broadPeak file. 
Specify the number of columns. Pass a second boolean to force the method.

=item add_ucsc_metadata

Add default column metadata for a UCSC refFlat or genePred file. 
Specify the number of columns to define the format. 
Pass a second boolean to force the method.

=item add_sgr_metadata

Add default column metadata for a SGR file. Pass a boolean to force the method.

=item add_standard_metadata

Add default column metadata for a generic file. Pass the text line 
containing the tab-delimited column headers.

=item standard_column_names

Returns an anonymous array of standard file format column header names. 
Pass a value representing the file format. Values include gff, bed12, 
bed6, bdg, narrowpeak, broadpeak, sgr, ucsc16, ucsc15, genepredext, 
ucsc12, knowngene, ucsc11, genepred, ucsc10, refflat.

=back

=head1 SEE ALSO

L<Bio::ToolBox::Data>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
