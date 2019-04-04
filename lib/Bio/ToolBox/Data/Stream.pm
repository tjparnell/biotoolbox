package Bio::ToolBox::Data::Stream;
our $VERSION = '1.66';

=head1 NAME

Bio::ToolBox::Data::Stream - Read, Write, and Manipulate Data File Line by Line

=head1 SYNOPSIS

  use Bio::ToolBox::Data;
  
  ### Open a pre-existing file
  my $Stream = Bio::ToolBox::Data->new(
        in      => 'regions.bed',
        stream  => 1,
  );
  
  # or directly
  my $Stream = Bio::ToolBox::Data::Stream->new(
        in      => 'regions.bed',
  );
  
  ### Open a new file for writing
  my $Stream = Bio::ToolBox::Data::Stream->new(
        out     => 'output.txt',
        columns => [qw(chromosome start stop name)],
  );
  
  
  ### Working line by line
  while (my $line = $Stream->next_line) {
  	  # get the positional information from the file data
  	  # assuming that the input file had these identifiable columns
  	  # each line is Bio::ToolBox::Data::Feature item
  	  my $seq_id = $line->seq_id;
  	  my $start  = $line->start;
  	  my $stop   = $line->end;
  	  
  	  # change values
  	  $line->value(1, 100); # index, new value
  }
  
  
  ### Working with two file streams
  my $inStream = Bio::ToolBox::Data::Stream->new(
        file    => 'regions.bed',
  );
  my $outStream = $inStream->duplicate('regions_ext100.bed');
  my $sc = $inStream->start_column;
  my $ec = $inStream->end_column;
  while (my $line = $inStream->next_line) {
      # adjust positions by 100 bp
      my $s = $line->start;
      my $e = $line->end;
      $line->value($sc, $s - 100);
      $line->value($ec, $e + 100);
      $outStream->write_row($line);
  }
  
  
  ### Finishing
  # close your file handles when you are done
  $Stream->close_fh;

=head1 DESCRIPTION

This module works similarly to the L<Bio::ToolBox::Data> object, except that 
rows are read from a file handle rather than a memory structure. This 
allows very large files to be read, manipulated, and even written without 
slurping the entire contents into a memory.

For an introduction to the L<Bio::ToolBox::Data> object and methods, refer to 
its documentation and the L<Bio::ToolBox::Data::Feature> documentation. 

Typically, manipulations are only performed on one row at a time, not on an 
entire table. Therefore, large scale table manipulations, such as sorting, is 
not possible. 

A typical workflow consists of opening two Stream objects, one for reading and 
one for writing. Rows are read, one at a time, from the read Stream, manipulated 
as necessary, and then written to the write Stream. Each row is passed as a 
L<Bio::ToolBox::Data::Feature> object. It can be manipulated as such, or the 
corresponding values may be dumped as an array. Working with the row data 
as an array is required when adding or deleting columns, since these manipulations 
are not allowed with a Feature object. The write Stream can then be passed 
either the Feature object or the array of values to be written.


=head1 METHODS

=head2 Initializing the structure

A new Bio::ToolBox::Data::Stream object may be generated directly, or indirectly 
through the L<Bio::ToolBox::Data> module.

=over 4

=item new

	my $Stream = Bio::ToolBox::Data::Stream->new(
	   in           => $filename,
	);
	my $Stream = Bio::ToolBox::Data->new(
	   stream       => 1,
	   in           => $filename,
	);

Options to the new function are listed below. Streams are inherently either 
read or write mode, determined by the mode given through the options.

=over 4

=item in

Provide the path of the file to open for reading. File types are 
recognized by the extension, and compressed files (.gz) are supported. File 
types supported include all those listed in L<Bio::ToolBox::file_helper>. 

=item out

Provide the path of the file to open for writing. No check is made 
for pre-existing files; if it exists it will be overwritten! A new data 
object is prepared, therefore column names must be provided. 

=item noheader

Boolean option indicating that the input file does not have file headers, 
in which case dummy headers are provided. This is not necessary for 
defined file types that don't normally have file headers, such as 
BED, GFF, or UCSC files. Ignored for output files.

=item columns

	my $Stream = Bio::ToolBox::Data::Stream->new(
	   out      => $filename,
	   columns  => [qw(Column1 Column2 ...)],
	);

When a new file is written, provide the names of the columns as an 
anonymous array. If no columns are provided, then a completely empty 
data structure is made. Columns must be added with the add_column() 
method below.

=item gff

When writing a GFF file, provide a GFF version. When this is given, the 
nine standard column names and metadata are automatically provided based 
on the file format specification. Note that the column names are not 
actually written in the file, but are maintained for internal use. 
Acceptable versions include 1, 2, 2.5 (GTF), and 3 (GFF3).

=item bed

When writing a BED file, provide the number of bed columns that the file 
will have. When this is given, the standard column names and metadata 
will be automatically provided based on the standard file format 
specification. Note that column names are not actually written to the file, 
but are maintained for internal use. Acceptable values are integers from 
3 to 12. 

=item ucsc

When writing a UCSC-style file format, provide the number of bed columns 
that the file will have. When this is given, the standard column names and 
metadata will be automatically provided based on the file format specification. 
Note that column names are not actually written to the file, but are maintained 
for internal use. Acceptable values include 10 (refFlat without gene names), 
11 (refFlat with gene names), 12 (knownGene gene prediction table), and 15 
(an extended gene prediction or genePredExt table).

=item gz

Boolean value to change the compression status of the output file. If 
overwriting an input file, the default is maintain the compression status, 
otherwise no compression. Pass a 0 for no compression, 1 for standard 
gzip compression, or 2 for block gzip (bgzip) compression for tabix 
compatibility.

=back

=item duplicate

   my $Out_Stream = $Stream->duplicate($new_filename);

For an opened-to-read Stream object, you may duplicate the object as a new 
opened-to_write Stream object that maintains the same columns and metadata. 
A new different filename must be provided. 

=back

=head2 General Metadata

There is a variety of general metadata regarding the Data structure that 
is available. 

The following methods may be used to access or set these 
metadata properties. Note that metadata is only written at the beginning 
of the file, and so must be set prior to iterating through the file.

=over 4

=item feature

Returns or sets the name of the features used to collect 
the list of features. The actual feature types are listed 
in the table, so this metadata is merely descriptive.

=item feature_type

Returns one of three specific values describing the contents 
of the data table inferred by the presence of specific column 
names. This provides a clue as to whether the table features 
represent genomic regions (defined by coordinate positions) or 
named database features. The return values include:

=over 4

=item coordinate: Table includes at least chromosome and start

=item named: Table includes name, type, and/or Primary_ID

=item unknown: unrecognized

=back

=item program

Returns or sets the name of the program generating the list.

=item database

Returns or sets the name or path of the database from which the 
features were derived.

=item gff

Returns or sets the version of loaded GFF files. Supported versions 
included 1, 2, 2.5 (GTF), and 3.

=item bed

Returns or sets the BED file version. Here, the BED version is simply 
the number of columns.

=item ucsc

Returns or sets the UCSC file format version. Here, the version is 
simply the number of columns. Supported versions include 10 (gene 
prediction), 11 (refFlat, or gene prediction with gene name), 12 
(knownGene table), 15 (extended gene prediction), or 16 (extended 
gene prediction with bin).

=item vcf

Returns or sets the VCF file version number. VCF support is limited.

=back

=head2 File information

These methods provide information about the file from which the 
data table was loaded. This does not include parsed annotation tables.

=over 4

=item filename

=item path

=item basename

=item extension

Returns the filename, full path, basename, and extension of 
the filename. Concatenating the last three values will reconstitute 
the first original filename.

=item add_file_metadata

  $Data->add_file_metadata('/path/to/file.txt');

Add filename metadata. This will automatically parse the path, 
basename, and recognized extension from the passed filename and 
set the appropriate metadata attributes.

=back

=head2 Comments

Comments are the other commented lines from a text file (lines 
beginning with a #) that were not parsed as metadata.

=over 4

=item comments

Returns a copy of the array containing commented lines.

=item add_comment

Appends the text string to the comment array.

=item delete_comment

Deletes a comment. Provide the array index of the comment to 
delete. If an index is not provided, ALL comments will be deleted!

=item vcf_headers

For VCF files, this will partially parse the VCF headers into a 
hash structure that can be queried or manipulated. Each header 
line is parsed for the primary key, being the first word after the 
## prefix, e.g. INFO, FORMAT, FILTER, contig, etc. For the simple 
values, they are stored as the value. For complex entries, such as 
with INFO and FORMAT, a second level hash is created with the ID 
extracted and used as the second level key. The value is always the 
always the remainder of the string.

For example, the following would be a simple parsed vcf header in 
code representation.

  $vcf_header = {
     FORMAT => {
        GT = q(ID=GT,Number=1,Type=String,Description="Genotype"),
        AD = q(ID=AD,Number=.,Type=Integer,Description="ref,alt Allelic depths"),
     },
     fileDate => 20150715,
  }

=item rewrite_vcf_headers

If you have altered the vcf headers exported by the vcf_headers() 
method, then this method will rewrite the hash structure as new 
comment lines. Do this prior to writing the new file stream
or else you will lose your changed VCF header metadata.

=back

=head2 Column Metadata

Information about the columns may be accessed. This includes the 
names of the column and shortcuts to specific identifiable columns, 
such as name and coordinates. In addition, each column may have 
additional metadata. Each metadata is a series of key =E<gt> 
value pairs. The minimum keys are 'index' (the 0-based index 
of the column) and 'name' (the column header name). Additional 
keys and values may be queried or set as appropriate. When the 
file is written, these are stored as commented metadata lines at 
the beginning of the file. Setting metadata is futile after 
reading or writing has begun.

=over 4

=item list_columns

Returns an array or array reference of the column names 
in ascending (left to right) order.

=item number_columns

Returns the number of columns in the Data table. 

=item last_column

Returns the array index of the last (rightmost) column in the 
Data table.

=item name

  $Stream->name($index, $new_name);
  my $name = $Stream->name($i);

Convenient method to return the name of the column given the 
index number. A column may also be renamed by passing a new name.

=item metadata

  $Stream->metadata($index, $key, $new_value);
  my $value = $Stream->metadata($index, $key)

Returns or sets the metadata value for a specific $key for a 
specific column $index.

This may also be used to add a new metadata key. Simply provide 
the name of a new $key that is not present

If no key is provided, then a hash or hash reference is returned 
representing the entire metadata for that column.

=item copy_metadata

  $Stream->copy_metadata($source, $target);

This method will copy the metadata (everything except name and 
index) between the source column and target column. Returns 1 if 
successful.  

=item delete_metadata

  $Stream->delete_metadata($index, $key);

Deletes a column-specific metadata $key and value for a specific 
column $index. If a $key is not provided, then all metadata keys 
for that index will be deleted.

=item find_column

  my $i = $Stream->find_column('Gene');
  my $i = $Stream->find_column('^Gene$')

Searches the column names for the specified column name. This 
employs a case-insensitive grep search, so simple substitutions 
may be made.

=item chromo_column

=item start_column

=item stop_column

=item strand_column

=item name_column

=item type_column

=item id_column

These methods will return the identified column best matching 
the description. Returns C<undef> if that column is not present. 
These use the L</find_column> method with a predefined list of 
aliases.

=back

=head2 Modifying Columns

These methods allow modification to the number and order of the 
columns in a Stream object. These methods can only be employed 
prior to opening a file handle for writing, i.e. before the first 
L</write_row> method is called. This enables one, for example, to 
duplicate a read-only Stream object to create a write-only Stream, 
add or delete columns, and then begin the row iteration.

=over 4

=item add_column

  my $i = $Stream->add_column($name);

Appends a new column at the rightmost position (highest 
index). It adds the column header name and creates a 
new column metadata hash. Pass a text string representing 
the new column name. It returns the new column index if 
successful.

=item copy_column

  my $j = $Stream->copy_column($i);

This will copy a column, appending the duplicate column at 
the rightmost position (highest index). It will duplicate 
column metadata as well. It will return the new index 
position.

=item delete_column

Deletes one or more specified columns. Any remaining 
columns rightwards will have their indices shifted 
down appropriately. If you had identified one of the 
shifted columns, you may need to re-find or calculate 
its new index.

=item reorder_column

  $Data->reorder_column($c,$b,$a,$a);

Reorders columns into the specified order. Provide the 
new desired order of indices. Columns could be duplicated 
or deleted using this method. The columns will adopt their 
new index numbers.

=back

=head2 Row Data Access

Once a file Stream object has been opened, and metadata and/or 
columns adjusted as necessary, then the file contents can be 
iterated through, one row at a time. This is typically a one-way 
direction. If you need to go back or start over, the easiest thing 
to do is re-open the file as a new Stream object. 

There are two main methods, L</next_row> for reading and L</write_row> 
for writing. They cannot and should not be used on the same Stream 
object.

=over 4

=item next_row

=item next_line

=item read_line

This method reads the next line in the file handle and returns a 
L<Bio::ToolBox::Data::Feature> object. This object represents the 
values in the current file row. 

Note that strand values and 0-based start coordinates are automatically 
converted to BioPerl conventions if required by the file type.

=item add_row

=item add_line

=item write_row

=item write_line

  $Data->add_row(\@values);
  $Data->add_row($Row); # Bio::ToolBox::Data::Feature object

This method writes a new row or line to a file handle. The first 
time this method is called the file handle is automatically opened for 
writing. Up to this point, columns may be manipulated. After this point, 
columns cannot be adjusted (otherwise the file structure becomes 
inconsistent).

This method may be implemented in one of three ways, based on the type 
data that is passed. 

=over 4

=item * A Feature object

A Feature object representing a row from another L<Bio::ToolBox::Data> 
data table or Stream. The values from this object will be automatically 
obtained. Modified strand and 0-based coordinates may be adjusted back 
as necessary.

=item * An array reference of values

Pass an array reference of values. The number of elements should match the 
number of expected columns. The values will be automatically joined using tabs. 
This implementation should be used if you using values from another Stream 
and the number of columns have been modified.

Manipulation of strand and 0-based starts may be performed if the 
metadata indicates this should be done.

=item * A string

Pass a text string. This assumes the column values are already tab 
concatenated. A new line character is appended if one is not included. 
No data manipulation (strand or 0-based starts) or sanity checking of the 
required number of columns is performed. Use with caution!

=back

=item iterate

    $Stream->iterate( sub {
       my $row = shift;
       my $number = $row->value($index);
       my $log_number = log($number);
       $row->value($index, $log_number);
    } );

A convenience method that will process a code reference for every line 
in the file. Pass a subroutine or code reference. The subroutine will 
receive the line as a L<Bio::ToolBox::Data::Feature> object, just as with 
the L</read_line> method. 

=back

=head2 File Handle methods 

The below methods work with the file handle. When you are finished with 
a Stream, you should be kind and close the file handle properly.

=over 4

=item mode

Returns the write mode of the Stream object. Read-only objects 
return false (0) and write-only Stream objects return true (1).

=item close_fh

Closes the file handle.

=item fh

Returns the L<IO::File> compatible file handle object representing 
the file handle. Use with caution.

=back

=head1 SEE ALSO

L<Bio::ToolBox::Data>, L<Bio::ToolBox::Data::Feature>

=cut

use strict;
use Carp qw(carp cluck croak confess);
use base 'Bio::ToolBox::Data::core';
use Bio::ToolBox::Data::Feature;

1;


#### Initialize ####

sub new {
	my $class = shift;
	my %args  = @_;
	
	# file arguments
	$args{in}  ||= $args{file} || undef;
	$args{out} ||= undef;
	unless ($args{in} or $args{out}) {
		cluck "a filename must be specified with 'in' or 'out' argument keys!\n";
		return;
	}
	if (defined $args{in} and defined $args{out}) {
		cluck "cannot define both 'in' and 'out' arguments!\n";
		return;
	} 
	$args{noheader} ||= 0;
	
	# prepare object
	my $self = $class->SUPER::new();
	
	# open an existing file for reading
	if ($args{in}) {
		
		# check and open file
		my $filename = $self->check_file($args{in});
		unless ($filename) {
			carp sprintf "file '%s' does not exist!", $args{in};
			return;
		}
		$self->add_file_metadata($filename);
		$self->open_to_read_fh or return;
		$self->{mode} = 0; # read mode
		
		# parse column headers
		$self->parse_headers($args{noheader});
		$self->{line_count} = $self->{header_line_count};
		
		# push a dummy row, this will get tossed when the first next_row() is called
		$self->{data_table}->[1] = $self->{'column_names'}; 
	}
	
	# prepare to write to a new stream
	elsif ($args{out}) {
		
		# add file name information
		$self->add_file_metadata($args{out});
		
		# we will not open the file handle quite yet in case the user 
		# wants to modify metadata
		$self->{mode} = 1; # set to write mode
		$self->{fh} = undef;
		
		# get names of columns user may have passed
		my @columns;
		if (exists $args{columns}) {
			@columns = @{ $args{columns} };
		}
		elsif (exists $args{datasets}) {
			@columns = @{ $args{datasets} };
		}
		
		# add the column names
		if (@columns) {
			foreach my $c (@columns) {
				$self->add_column($c);
			}
		}
		elsif (exists $args{gff} and $args{gff}) {
			# use standard names for the number of columns indicated
			# we trust that the user knows the subtle difference between gff versions
			$self->add_gff_metadata($args{gff});
			unless ($self->extension =~ /g[tf]f/) {
				$self->{extension} = $args{gff} == 2.5 ? '.gtf' : 
					$args{gff} == 3 ? '.gff3' : '.gff';
			}
		}
		elsif (exists $args{bed} and $args{bed}) {
			# use standard names for the number of columns indicated
			unless ($args{bed} =~ /^\d{1,2}$/ and $args{bed} >= 3) {
				carp "bed parameter must be an integer 3-12!";
				return;
			}	
			$self->add_bed_metadata($args{bed});
			unless ($self->extension =~ /bed|peak/) {
				$self->{extension} = '.bed';
			}
		}
		elsif (exists $args{ucsc} and $args{ucsc}) {
			# a ucsc format such as refFlat, genePred, or genePredExt
			my $u = $self->add_ucsc_metadata($args{ucsc});
			unless ($u) {
				carp "unrecognized number of columns for ucsc format!";
				return;
			};
			unless ($self->extension =~ /ucsc|ref+lat|genepred/) {
				$self->{extension} = '.ucsc';
			}
		}
		# else it will be an empty object with no columns
		
		# append gz if necessary
		if (exists $args{gz} and $args{gz} and $self->extension !~ /gz$/) {
			$self->{extension} .= '.gz';
		}
		
		# rebuild the filename after modifying the extension
		$self->{filename} = $self->{path} . $self->{basename} . $self->{extension};
		
		# add feature
		$args{feature} ||= $args{features} || undef;
		$self->feature($args{feature}) unless $self->feature;
	}
	
	return $self;
}


sub duplicate {
	my ($self, $filename) = @_;
	unless ($filename) {
		carp "a new filename must be provided!";
		return;
	}
	if ($filename eq $self->filename) {
		carp "provided filename is not unique from that in metadata!";
		return;
	}
	
	# duplicate the data structure
	my $columns = $self->list_columns;
	my $Dup = $self->new(
		'out' => $filename, 
		'columns' => $columns,
	) or return;
	
	# copy the metadata
	for (my $i = 0; $i < $self->number_columns; $i++) {
		# column metadata
		my %md = $self->metadata($i);
		$Dup->{$i} = \%md;
	}
	foreach (qw(feature program db bed gff vcf ucsc headers)) {
		# various keys
		$Dup->{$_} = $self->{$_};
	}
	my @comments = $self->comments;
	push @{$Dup->{comments}}, @comments;
	
	return $Dup;
}



### Column manipulation

sub add_column {
	my ($self, $name) = @_;
	return unless $name;
	unless ($self->mode) {
		cluck "We have a read-only Stream object, cannot add columns";
		return;
	}
	if (defined $self->{fh}) {
		# Stream file handle is opened
		cluck "Cannot modify columns when a Stream file handle is opened!";
		return;
	}
	
	my $column = $self->number_columns;
	$self->{$column} = {
		'name'      => $name,
		'index'     => $column,
	};
	$self->{data_table}->[0][$column] = $name;
	$self->{number_columns}++;
	delete $self->{column_indices} if exists $self->{column_indices};
	if ($self->gff or $self->bed or $self->ucsc or $self->vcf) {
		# check if we maintain integrity, at least insofar what we test
		$self->verify(1); # silence so user doesn't get these messages
	}
	return $column;
}

sub copy_column {
	my $self = shift;
	unless ($self->mode) {
		confess "We have a read-only Stream object, cannot add columns";
	}
	if (defined $self->{fh}) {
		# Stream file handle is opened
		confess "Cannot modify columns when a Stream file handle is opened!";
	}
	my $index = shift;
	return unless defined $index;
	
	my $new_index = $self->add_column( $self->name($index) );
	$self->copy_metadata($index, $new_index);
	return $new_index;
}



#### Row Access ####

*next_line = *read_line = \&next_row;

sub next_row {
	my $self = shift;
	if ($self->{mode}) {
		confess "Stream object is write-only! cannot read";
	}
	
	# read and add the next line in the file
	my $line = $self->{fh}->getline or return;
	$self->{line_count}++;
	if (substr($line,0,1) eq '#') {
		# we shouldn't have internal comment lines, but just in case....
		# could be a gff3 pragma
		$self->add_comment($line);
		return $self->next_row;
	}
	
	# add the current line to the data table as row 1
	pop @{ $self->{data_table} }; # remove the old line
	$self->add_data_line($line);
	
	# return the feature
	return Bio::ToolBox::Data::Feature->new(
		'data'      => $self,
		'index'     => 1, 
	);	
}


*add_row = *add_line = *write_line = \&write_row;

sub write_row {
	my $self = shift;
	my $data = shift;
	unless ($self->{mode}) {
		confess "Stream object is read-only! cannot write";
	}
	
	# open the file handle if it hasn't been opened yet
	unless (defined $self->{fh}) {
		# we first write a standard empty data file with metadata and headers
		my $newfile = $self->write_file($self->filename);
		unless ($newfile) {
			die "unable to write file!";
		}
		
		# just in case the filename is changed when writing the file
		if ($newfile ne $self->filename) {
			$self->add_file_metadata($newfile);
		}
		
		# then we re-open the file for appending
		my $fh = $self->open_to_write_fh($newfile, undef, 1) or 
			die "unable to append to file $newfile!";
		$self->{fh} = $fh;
	}
	
	# identify what kind of data we are dealing with
	my $data_ref = ref $data;
	if ($data_ref eq 'Bio::ToolBox::Data::Feature') {
		# user passed a Feature object
		$self->{fh}->print( join("\t", ($data->row_values)), "\n" );
	}
	elsif ($data_ref eq 'ARRAY') {
		# user passed an array of values
		$self->{fh}->print( join("\t", @$data), "\n");
	}
	else {
		# assume the passed data is a string
		# make sure it has a newline
		unless ($data =~ /\n$/) {
			$data .= "\n";
		}
		$self->{fh}->print($data);
	}
	return 1;
}

sub iterate {
	my $self = shift;
	my $code = shift;
	unless (ref $code eq 'CODE') {
		cluck "iterate_function() method requires a code reference!";
		return;
	}
	while (my $row = $self->next_row) {
		&$code($row);
	}
	return 1;
}




#### File handle ####

sub mode {
	my $self = shift;
	return $self->{mode};
}

sub DESTROY {
	my $self = shift;
	$self->close_fh;
}


####################################################

__END__

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
