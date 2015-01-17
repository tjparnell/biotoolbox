package Bio::ToolBox::Data::Stream;
our $VERSION = 1.24;

=head1 NAME

Bio::ToolBox::Data::Stream - Read, Write, and Manipulate Data File Line by Line

=head1 SYNOPSIS
  
  use Bio::ToolBox::Data;
  
  ### Open a pre-existing file
  my $Stream = Bio::ToolBox::Data->new(
        file    => 'regions.bed',
        stream  => 1,
  );
  
  # or directly
  my $Stream = Bio::ToolBox::Data::Stream->new(
        file    => 'regions.bed',
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

=over 4

=item new()

Create a new Bio::ToolBox::Data::Stream object. For simplicity, a new file 
may also be opened using the L<Bio::ToolBox::Data> new function.
	
	my $Stream = Bio::ToolBox::Data->new(
	   stream       => 1,
	   file         => $filename,
	);

Options to the new function are listed below. Streams are inherently either 
read or write mode, determined by the mode given through the options.

=over 4

=item file =E<gt> $filename

Provide the path and name of the file to open. File types are recognized by 
the extension, and compressed files (.gz) are supported. File types supported 
include all those listed in L<Bio::ToolBox::file_helper>. Files are 
checked for existence. Existing files are assumed to be read, and non-existent 
files are assumed to be written, unless otherwise specified by the mode 
option. This option is required.

=item overwrite =E<GT> boolean

If a file exists and you wish to overwrite, pass this option with a true value.

=item columns =E<gt> [qw(Column1 Column2 ...)]

When a new file is written, provide the names of the columns as an 
anonymous array. 

=back

=item duplicate($filename)

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

=over

=item feature()

=item feature($text)

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

=item program($name)

Returns or sets the name of the program generating the list.

=item database($name)

Returns or sets the name or path of the database from which the 
features were derived.

=back

The following methods may be used to access metadata only.

=over

=item gff

=item bed

Returns the GFF version number or the number of BED columns 
indicating that the Data structure is properly formatted as 
such. A value of 0 means they are not formatted as such.

=back

=head2 File information

=over 4

=item filename($text)

Returns or sets the filename for the Data structure. If you set 
a new filename, the path, basename, and extension are 
automatically derived for you. If a path was not provided, 
the current working directory is assumed. 

=item path

=item basename

=item extension

Returns the full path, basename, and extension of the filename. 
Concatenating these three values will reconstitute the 
original filename.

=back

=head2 Comments

Comments are the other commented lines from a text file (lines 
beginning with a #) that were not parsed as metadata.

=over 4

=item comments

Returns a copy of the array containing commented lines.

=item add_comment($text)

Appends the text string to the comment array.

=item delete_comment

=item delete_comment($index)

Deletes a comment. Provide the array index of the comment to 
delete. If an index is not provided, ALL comments will be deleted!

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

=item name($index)

Convenient method to return the name of the column given the 
index number.

=item metadata($index, $key)

=item metadata($index, $key, $new_value)

Returns or sets the metadata value for a specific $key for a 
specific column $index.

This may also be used to add a new metadata key. Simply provide 
the name of a new $key that is not present

If no key is provided, then a hash or hash reference is returned 
representing the entire metadata for that column.

=item find_column($name)

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
These use the find_column() method with a predefined list of 
aliases.

=back

=head2 Modifying Columns

These methods allow modification to the number and order of the 
columns in a Stream object. These methods can only be employed 
prior to opening a file handle for writing, i.e. before the first 
write_row() method is called. This enables one, for example, to 
duplicate a read-only Stream object to create a write-only Stream, 
add or delete columns, and then begin the row iteration.

=over

=item add_column($name)

Appends a new column at the rightmost position (highest 
index). It adds the column header name and creates a 
new column metadata hash. Pass a text string representing 
the new column name. It returns the new column index if 
successful.

=item copy_column($index)

This will copy a column, appending the duplicate column at 
the rightmost position (highest index). It will duplicate 
column metadata as well. It will return the new index 
position.

=item delete_column($index1, $index2, ...)

Deletes one or more specified columns. Any remaining 
columns rightwards will have their indices shifted 
down appropriately. If you had identified one of the 
shifted columns, you may need to re-find or calculate 
its new index.

=item reorder_column($index1,  $index, ...)

Reorders columns into the specified order. Provide the 
new desired order of indices. Columns could be duplicated 
or deleted using this method. The columns will adopt their 
new index numbers.

=head2 Row Data Access

Once a file Stream object has been opened, and metadata and/or 
columns adjusted as necessary, then the file contents can be 
iterated through, one row at a time. This is typically a one-way 
direction. If you need to go back or start over, the easiest thing 
to do is re-open the file as a new Stream object. 

There are two main methods, next_row() for reading and write_row() 
for writing. They cannot and should not be used on the same Stream 
object.

=over 4

=item next_row()

This method reads the next line in the file handle and returns a 
L<Bio::ToolBox::Data::Feature> object. This object represents the 
values in the current file row. 

Note that strand values and 0-based start coordinates are automatically 
converted to BioPerl conventions if required by the file type.

=item add_row()

=item write_row()

This method writes a new row or line to a file handle. The first 
time this method is called the file handle is automatically opened for 
writing. Up to this point, columns may be manipulated. After this point, 
columns cannot be adjusted (otherwise the file structure becomes 
inconsistent).

This method may be implemented in one of three ways, based on the type 
data that is passed. 

=over 4

=item A <Bio::ToolBox::Data::Feature> object

A Feature object representing a row from another <Bio::ToolBox::Data> 
data table or Stream. The values from this object will be automatically 
obtained. B<Note:> Only pass this object if the number and names of the columns 
are identical between read and write Streams, otherwise very strange 
things may happen! If you modify the number of columns, then use the second 
approach below. Modified strand and 0-based coordinates may be adjusted back 
as necessary.

=item An array of values

Pass an array of values. The number of elements should match the number 
of expected columns. The values will be automatically joined using tabs. 
This implementation should be used if you using values from another Stream 
and the number of columns have been modified.

Manipulation of strand and 0-based starts may be performed if the 
metadata indicates this should be done.

=item A string

Pass a text string. This assumes the values are already concatenated. 
A new line character is appended if one is not included. No data 
manipulation (strand or 0-based starts) or sanity checking of the 
required number of columns is performed. Use with caution!

=back

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

=cut


use strict;
use Carp qw(carp cluck croak confess);

use Bio::ToolBox::data_helper qw(
	generate_data_structure
	find_column_index
);
use Bio::ToolBox::file_helper qw(
	open_data_file
	write_data_file
	open_to_write_fh
	parse_filename
	process_data_line
	check_file
);
use Bio::ToolBox::Data::common;
use Bio::ToolBox::Data::Feature;

1;



#### Initialize ####

sub new {
	my $class = shift;
	my %args  = @_;
	
	$args{features} ||= $args{feature} || 'feature';
	$args{filename} ||= $args{file} || undef;
	$args{overwrite} ||= 0;
	
	unless (defined $args{filename}) {
		carp "a filename must be provided!";
		return;
	}
	
	# prepare
	my $data;
	my $fh;
	my $mode;
	
	# check if the file exists
	my $filename = check_file($args{filename});
	
	# default behavior is to open an existing file for reading
	if ($filename and not $args{overwrite}) {
		($fh, $data) = open_data_file($filename);
		unless ($data) {
			croak "unable to read file $filename!";
		}
		$mode = 0; # read mode
		
		# add the column headers
		push @{ $data->{data_table} }, $data->{'column_names'};
		delete $data->{'column_names'}; # we no longer need this
		
		# look for potential start columns to convert from 0-based
		my @starts;
		if ($data->{'ucsc'} or $data->{'bed'}) {
			foreach my $name (qw(start txStart cdsStart peak)) {
				my $c = find_column_index($data, $name);
				next unless defined $c;
				next if (exists $data->{$c}{'base'} and $data->{$c}{'base'} == 1);
				push @starts, $c;
			}
		}
		$data->{column_starts} = \@starts;
		
		# potential strand column that may need to be converted
		my $strand_i = find_column_index($data, '^strand$');
		if (defined $strand_i) {
			$data->{strand_check} = $strand_i;
		}
	}
	
	# alternate behavior if file does not exist is to create an empty stream
	# and prepare it for writing
	else {
		my @datasets;
		if (exists $args{datasets}) {
			@datasets = @{ $args{datasets} };
		}
		elsif (exists $args{columns}) {
			@datasets = @{ $args{columns} };
		}
		unless (@datasets) {
			carp "no column names provided!";
			return;
		}
		my $feature = $args{features} || 'feature';
		$data = generate_data_structure($feature, @datasets);
		
		# add file name information
		my ($basename, $path, $extension) = parse_filename($filename);
		$data->{filename}  = $filename;
		$data->{basename}  = $basename;
		$data->{path}      = $path;
		$data->{extension} = $extension;
		
		# add extra information for checking strand and start information
		$data->{strand_check} = undef;
		$data->{column_starts} = [];
		
		# we will not open the file handle quite yet in case the user 
		# wants to modify metadata
		$mode = 1; # write mode
	}
	
	$data->{fh} = $fh;
	$data->{mode} = $mode;
	return bless $data, $class;
}


sub duplicate {
	my $self = shift;
	my $filename = shift;
	unless ($filename) {
		carp "a new filename must be provided!";
		return;
	}
	if ($filename eq $self->filename) {
		carp "provided filename is not unique from that in metadata!";
		return;
	}
	
	# duplicate the data structure
	my $data = generate_data_structure( $self->feature, $self->list_columns );
	
	# copy the metadata
	for (my $i = 0; $i < $self->number_columns; $i++) {
		# column metadata
		my %md = $self->metadata($i);
		$data->{$i} = \%md;
	}
	foreach (qw(program db bed gff ucsc headers)) {
		# various keys
		$data->{$_} = $self->{$_};
	}
	$data->{column_starts} = [ @{ $self->{column_starts} } ];
	my @comments = $self->comments;
	push @{$data->{other}}, @comments;
	
	# add file name information
	undef $data->{fh};
	$data->{mode} = 1;
	my ($basename, $path, $extension) = parse_filename($filename);
	$data->{filename}  = $filename;
	$data->{basename}  = $basename;
	$data->{path}      = $path;
	$data->{extension} = $extension;
	
	return bless $data, 'Bio::ToolBox::Data::Stream';
}



### Column manipulation

# see also Bio::ToolBox::Data::common for imported methods

sub add_column {
	my ($self, $name) = @_;
	return unless $name;
	if (defined $self->{fh}) {
		# Stream file handle is opened
		cluck "Cannot modify columns when a Stream file handle is opened!";
		return;
	}
	
	my $column = $self->{number_columns};
	$self->{$column} = {
		'name'      => $name,
		'index'     => $column,
	};
	$self->{data_table}->[0][$column] = $name;
	$self->{number_columns}++;
	delete $self->{column_indices} if exists $self->{column_indices};
	return $column;
}

sub copy_column {
	my $self = shift;
	if (defined $self->{fh}) {
		# Stream file handle is opened
		cluck "Cannot modify columns when a Stream file handle is opened!";
		return;
	}
	my $index = shift;
	return unless defined $index;
	
	my $new_index = $self->add_column( $self->name($index) );
	$self->copy_metadata($index, $new_index);
	return $new_index;
}



#### Row Access ####

sub next_row {
	my $self = shift;
	if ($self->mode) {
		cluck "Stream object is write-only! cannot read";
		return;
	}
	
	# process the data line, converting strand and 0-based coordinates as necessary
	my $line = $self->{fh}->getline;
	return unless $line;
	my ($linedata, $plusminus) = process_data_line($line, $self->number_columns, 
		$self->strand_column, @{ $self->{column_starts} });
	
	# update the data table
	# it will always be row 1
	$self->{data_table}->[1] = $linedata;
	
	# update metadata as necessary
	if (defined $self->{strand_check} and $plusminus) {
		$self->metadata($self->strand_column, 'strand_style', 'plusminus');
		my $auto = $self->metadata($self->strand_column, 'AUTO');
		if ($auto) {
			# update automatically generated metadata
			$self->metadata($self->strand_column, 'AUTO', $auto++);
		}
		undef $self->{strand_check}; # we no longer need to check this
	}
	
	# return the feature
	return Bio::ToolBox::Data::Feature->new(
		'data'      => $self,
		'index'     => 1, 
	);	
}


sub add_row {
	shift->write_row(@_);
}


sub write_row {
	my $self = shift;
	my $data = shift;
	unless ($self->mode) {
		cluck "Stream object is read-only! cannot write";
		return;
	}
	
	# open the file handle if it hasn't been opened yet
	unless (defined $self->{fh}) {
		# we first write a standard empty data file with metadata and headers
		my $filename = write_data_file(
			data     => $self,
			filename => $self->filename,
		);
		unless ($filename) {
			die "unable to write file!";
		}
		
		# just in case the filename is changed when writing the file
		$self->filename($filename);
		
		# then we re-open the file for appending
		my $fh = open_to_write_fh($filename, undef, 1) or 
			die "unable to append to file $filename!";
		$self->{fh} = $fh;
	}
	
	# identify what kind of data we are dealing with
	my $data_ref = ref $data;
	my @values;
	if ($data_ref eq 'Bio::ToolBox::Data::Feature') {
		# user passed a Feature object
		
		# get the values
		@values = $data->row_values;
		
		# check strand information
		unless (defined $self->{strand_check}) {
			# we do this check only once, presuming that if it's true for the first 
			# one, it's probably true for all
			# too expensive to do every time we perform a write, right?
			# must dig back into the original Stream object, which may or may not be self
			if ($data->{data}->strand_column) {
				if ($data->{data}->metadata( $data->{data}->strand_column, 'strand_style') 
					eq 'plusminus'
				) {
					# we do need to convert strand back
					$self->{strand_check} = 1;
				}
				else {
					# we should not convert strand back
					$self->{strand_check} = 0;
				}
			}
			
			# also check start coordinates while we're at it
			if (scalar @{ $data->{data}->{column_starts} } > 0) {
				$self->{column_starts} = [ @{ $data->{data}->{column_starts} } ];
			}
			else {
				$self->{column_starts} = [];
			}
		}
	}
	elsif ($data_ref eq 'ARRAY') {
		# user passed an array of values
		@values = @$data;
		
		# check strand information
		unless (defined $self->{strand_check}) {
			# we do this check only once, presuming that if it's true for the first 
			# one, it's probably true for all
			# too expensive to do every time we perform a write, right?
			# use the metadata from the current Stream object which may not be accurate
			if ($self->strand_column) {
				if (
					$self->metadata($self->strand_column, 'strand_style') eq 'plusminus'
				) {
					# we need to convert strand back
					$self->{strand_check} = 1;
				}
				else {
					# we do not need to convert strand back
					$self->{strand_check} = 0;
				}
			}
		}
	}
	
	# write the values as necessary
	if (@values) {
		# we have an array of values to write
		# make adjustments as necessary
		if ($self->{strand_check}) {
			# strand adjustment
			my $i = $data->{data}->strand_column;
			if ($values[$i] >= 0) {
				$values[$i] = '+';
			}
			else {
				$values[$i] = '-';
			}
		}
		if ( @{$self->{column_starts}} ) {
			# 0-based start coordinate adjustment
			foreach my $c ( @{$self->{column_starts}} ) {
				$values[$c] -= 1;
			}
		}
		
		# write the values
		$self->{fh}->print( join("\t", @values), "\n" );
	}
	
	else {
		# no values given
		# assume the passed data is a string
		
		# make sure it has a newline
		unless ($data =~ /\n$/) {
			$data .= "\n";
		}
		$self->{fh}->print($data);
	}
}



#### File handle ####

sub mode {
	my $self = shift;
	return $self->{mode};
}

sub fh {
	my $self = shift;
	return $self->{fh};
}

sub close_fh {
	my $self = shift;
	$self->{fh}->close;
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
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
