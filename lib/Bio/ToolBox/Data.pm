package Bio::ToolBox::Data;
our $VERSION = 1.24;

=head1 NAME

Bio::ToolBox::Data - Reading, writing, and manipulating data structure

=head1 SYNOPSIS
  
  use Bio::ToolBox::Data;
  
  ### Create new gene list from database
  my $Data = Bio::ToolBox::Data->new(
        db      => 'hg19',
        feature => 'gene:ensGene',
  );
  
  my $Data = Bio::ToolBox::Data->new(
        db      => 'hg19',
        feature => 'genome',
        win     => 1000,
        step    => 1000,
  );
  
  
  ### Open a pre-existing file
  my $Data = Bio::ToolBox::Data->new(
        file    => 'coordinates.bed',
  );
  
  
  ### Get a specific value
  my $value = $Data->value($row, $column);
  
  
  ### Replace or add a value
  $Data->value($row, $column, $new_value);
  
  
  ### Iterate through a Data structure one row at a time
  my $stream = $Data->row_stream;
  while (my $row = $stream->next_row) {
  	  # get the positional information from the file data
  	  # assuming that the input file had these identifiable columns
  	  my $seq_id = $row->seq_id;
  	  my $start  = $row->start;
  	  my $stop   = $row->end;
  	  
  	  # generate a Bio::Seq object from the database using 
  	  # these coordinates 
  	  my $region = $db->segment($seq_id, $start, $stop);
  	  
  	  # modify a row value
  	  my $value = $row->value($column);
  	  my $new_value = $value + 1;
  	  $row->value($column, $new_value);
  }
  
  
  ### write the data to file
  my $success = $Data->write_file(
       filename     => 'new_data.txt',
       gz           => 1,
  );
  print "wrote new file $success\n"; # file is new_data.txt.gz
  
=head1 DESCRIPTION

This module works with the primary Bio::ToolBox Data structure. Simply, it 
is a complex data structure representing a tabbed-delimited table (array 
of arrays), with plenty of options for metadata. Many common bioinformatic 
file formats are simply tabbed-delimited text files (think BED and GFF). 
Each row is a feature or genomic interval, and each column is a piece of 
information about that feature, such as name, type, and/or coordinates. 
We can append to that file additional columns of information, perhaps 
scores from genomic data sets. We can record metadata regarding how 
and where we obtained that data. Finally, we can write the updated 
table to a new file.

=head1 METHODS

=head2 Initializing the structure

=over 4

=item new()

Initialize a new Data structure. This generally requires options, 
provided as an array of key =E<gt> values. A new list of features 
may be obtained from an annotation database or an existing file 
may be loaded. If you do not pass any options, a new empty 
structure will be generated for you to populate. 

These are the options available.

=over 4

=item file =E<gt> $filename

Provide the path and name to an existing tabbed-delimited text 
file. BED and GFF files and their variants are accepted. Except 
for structured files, e.g. BED and GFF, the first line is 
assumed to be column header names. Commented lines (beginning 
with #) are parsed as metadata. The files may be compressed 
(gzip or bzip2).

=item feature =E<gt> $type

=item feature =E<gt> "$type:$source"

=item feature =E<gt> 'genome'

For de novo lists from an annotation database, provide the GFF 
type or type:source (columns 3 and 2) for collection. A comma 
delimited string may be accepted (not an array). 

For a list of genomic intervals across the genome, specify a 
feature of 'genome'.

=item db =E<gt> $name

=item db =E<gt> $path

=item db =E<gt> $database_object

Provide the name of the database from which to collect the 
features. It may be a short name, whereupon it is checked in 
the Bio::ToolBox configuration file C<.biotoolbox.cfg> for 
connection information. Alternatively, a path to a database 
file or directory may be given. 

If you already have an opened Bio::DB::SeqFeature::Store database 
object, you can simply pass that. See Bio::ToolBox::db_helper for 
more information. However, this in general should be discouraged, 
since the name of the database will not be properly recorded when 
saving to file. It is better to simply pass the name of database 
again; multiple connections to the same database are smartly handled 
in the background.

=item win =E<gt> $integer

=item step =E<gt> $integer

If generating a list of genomic intervals, optionally provide 
the window and step values. Default values are defined in 
the Bio::ToolBox configuration file C<.biotoolbox.cfg>.

=item columns =E<gt> [qw(Column1 Column2 ...)]

=item datasets =E<gt> [qw(Column1 Column2 ...)]

When no file is given or database given to search, then a new, 
empty Data object is returned. In this case, you may optionally 
provide the column names in advance as an anonymous array. You 
may also optionally provide a general feature name, if desired.

=back

If successful, the method will return a Bio::ToolBox::Data object.

=back

=head2 General Metadata

There is a variety of general metadata regarding the Data 
structure. 

The following methods may be used to access or set these 
metadata properties.  

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

=head2 The Data table

The Data table is the array of arrays containing all of the 
actual information. 

=over 4

=item list_columns

Returns an array or array reference of the column names 
in ascending (left to right) order.

=item number_columns

Returns the number of columns in the Data table. 

=item last_row

Returns the array index number of the last row. 
Since the header row is index 0, this is also the 
number of actual content rows.

=item column_values($index)

Returns an array or array reference representing the values 
in the specified column. This includes the column header as 
the first element. Pass the method the column index.

=item add_column($name)

=item add_column(\@column_data)

Appends a new column to the Data table at the 
rightmost position (highest index). It adds the column 
header name and creates a new column metadata hash. 
Pass the method one of two possibilities. Pass a text 
string representing the new column name, in which case 
no data will be added to the column. Alternatively, pass 
an array reference, and the contents of the array will 
become the column data. If the Data table already has 
rows, then the passed array reference must have the same 
number of elements.

It returns the new column index if successful.

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

=item add_row

=item add_row(\@values)

Add a new row of data values to the end of the Data table. 
Optionally provide a reference to an array of values to 
put in the row. The array is filled up with C<undef> for 
missing values, and excess values are dropped.

=item delete_row($row1, $row2, ...)

Deletes one or more specified rows. Rows are spliced out 
highest to lowest index to avoid issues. Be very careful 
deleting rows while simultaneously iterating through the 
table!

=item row_values($row)

Returns a copy of an array for the specified row index. 
Modifying this returned array does not migrate back to the 
Data table; Use the value method below instead.

=item value($row, $column)

=item value($row, $column, $new_value)

Returns or sets the value at a specific row or column index.
Index positions are 0-based (header row is index 0). 

=back

=head2 Column Metadata

Each column has metadata. Each metadata is a series of key =E<gt> 
value pairs. The minimum keys are 'index' (the 0-based index 
of the column) and 'name' (the column header name). Additional 
keys and values may be queried or set as appropriate. When the 
file is written, these are stored as commented metadata lines at 
the beginning of the file.

=over 4

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

=item copy_metadata($source, $target)

This method will copy the metadata (everything except name and 
index) between the source column and target column. Returns 1 if 
successful.  

=item delete_metadata($index, $key);

Deletes a column-specific metadata $key and value for a specific 
column $index. If a $key is not provided, then all metadata keys 
for that index will be deleted.

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

=head2 Efficient Data Access

Most of the time we need to iterate over the Data table, one row 
at a time, collecting data or processing information. These methods 
simplify the process.

=over 4

=item iterate(sub {})

This method will process a code reference on every row in the data 
table. Pass a subroutine or code reference. The subroutine will 
receive the row as a Bio::ToolBox::Data::Feature object. With this 
object, you can retrieve values, set values, and add new values. 
For example
    
    $Data->iterate( sub {
       my $row = shift;
       my $number = $row->value($index);
       my $log_number = log($number);
       $row->value($index, $log_number);
    } );

=item row_stream()

This returns an Bio::ToolBox::Data::Iterator object, which has one 
method, next_row(). Call this method repeatedly until it returns 
C<undef> to work through each row of data.

Users of the Bio::DB family of database adaptors may recognize the 
analogy to the seq_stream() method.

=item next_row()

Called from a Bio::ToolBox::Data::Iterator object, it returns a 
Bio::ToolBox::Data::Feature object. This object represents the 
values in the current Data table row.

An example using the iterator is shown below.
  
  my $stream = $Data->row_stream;
  while (my $row = $stream->next_row) {
     # each $row is a Bio::ToolBox::Data::Feature object
     # representing the row in the data table
     my $value = $row->value($index);
     # do something with $value
  }

=back

=head2 Data Table Functions

These methods alter the Data table en masse. 

=over 4

=item verify()

This method will verify the Data structure, including the metadata and the 
Data table. It ensures that the table has the correct number of rows and 
columns as described in the metadata, and that each column has the basic 
metadata. 

If the Data structure is marked as a GFF or BED structure, then the table 
is checked that the structure matches the proper format. If not, for 
example when additional columns have been added, then the GFF or BED value 
is set to null. 

This method is automatically called prior to writing the Data table to file.

=item sort_data($index, $direction);

This method will sort the Data table by the values in the indicated column. 
It will automatically determine whether the contents of the column are 
numbers or alphanumeric, and will sort accordingly, either numerically or 
asciibetically. The first non-null value in the column is used to determine. 
The sort may fail if the values are not consistent. Pass a second optional 
value to indicate the direction of the sort. The value should be either 
'i' for 'increasing' or 'd' for 'decreasing'. The default order is 
increasing. 

=item gsort_data

This method will sort the Data table by increasing genomic coordinates. It 
requires the presence of chromosome and start (or position) columns, 
identified by their column names. These are automatically identified. 
Failure to find these columns mean a failure to sort the table. Chromosome 
names are sorted first by their digits (e.g. chr2 before chr10), and then 
alphanumerically. Base coordinates are sorted by increasing value. 
Identical positions are kept in their original order.

=item splice_data($current_part, $total_parts)

This method will splice the Data table into $total_parts number of pieces, 
retaining the $current_part piece. The other parts are discarded. This 
method is intended to be used when a program is forked into separate 
processes, allowing each child process to work on a subset of the original 
Data table. 

Two values are passed to the method. The first is the current part number, 
1-based. The second value is the total number of parts that the table 
should be divided, corresponding to the number of concurrent processes. 
One easy approach to forking is to use L<Parallel::ForkManager>. The 
example below shows how to fork into four concurrent processes.
	
	my $Data = Bio::ToolBox::Data->new(file => $file);
	my $pm = Parallel::ForkManager->new(4);
	for my $i (1..4) {
		$pm->start and next;
		### in child ###
		$Data->splice_data($i, 4);
		# do something with this portion
		# then save to a temporary unique file
		$Data->save("$file_$i");
		$pm->finish;
	}
	$pm->wait_all_children;
	# reload children files
	$Data->reload_children(glob "$file_*");

Since each forked child process is separate from their parent process, 
their contents must be reloaded into the current Data object. The 
L<Parallel::ForkManager> documentation recommends going through a disk 
file intermediate. Therefore, write each child Data object to file using 
a unique name. Once all children have been reaped, they can be reloaded 
into the current Data object using the reload_children() method.

Remember that if you fork your script into child processes, any database 
connections must be re-opened; they are typically not clone safe. If you 
have an existing database connection by using the open_database() method, 
it should be automatically re-opened for you when you use the splice_data() 
method, but you will need to call open_database() again in the child 
process to obtain the new database object.

=item reload_children(@children_files)

Discards the current data table in memory and reloads two or more files 
written from forked children processes. Provide the name of the child 
files in the order you want them loaded. The files will be automatically 
deleted if successfully loaded. Returns the number of lines reloaded on 
success.

=back

=head2 File Functions

When you are finished modifying the Data table, it may then be written out 
as a tabbed-delimited text file. If the format corresponds to a valide BED or 
GFF file, then it may be written in that format. 

Several functions are available for writing the Data table, exporting to a 
compatible GFF file format, or writing a summary of the Data table.

=over 4

=item save()

=item write_file()

=item write_file($filename)

=item write_file(%options)

These methods will write the Data structure out to file. It will 
be first verified as to proper structure. Opened BED and GFF files 
are checked to see if their structure is maintained. If so, they 
are written in the same format; if not, they are written as regular 
tab-delimited text files. You may pass additional options.

=over 4

=item filename =E<gt> $filename

Optionally pass a new filename. Required for new objects; previous 
opened files may be overwritten if a new name is not provided. If 
necessary, the file extension may be changed; for example, BED files 
that no longer match the defined format lose the .bed and gain a .txt 
extension. Compression may or add or strip .gz as appropriate. If 
a path is not provided, the current working directory is used.

=item gz =E<gt> boolean

Change the compression status of the output file. The default is to 
maintain the status of the original opened file.

=back

If the file save is successful, it will return the full path and 
name of the saved file, complete with any changes to the file extension.

=item summary_file(%options)

Write a separate file summarizing columns of data (mean values). 
The mean value of each column becomes a row value, and each column 
header becomes a row identifier (i.e. the table is transposed). The 
best use of this is to summarize the mean profile of windowed data 
collected across a feature. See the Bio::ToolBox scripts 
L<get_relative_data.pl> and L<average_gene.pl> as examples. 
You may pass these options. 

=over 4

=item filename =E<gt> $filename

Pass an optional new filename. The default is to take the basename 
and append "_summed" to it.

=item startcolumn =E<gt> $index

=item stopcolumn =E<gt> $index

Provide the starting and ending columns to summarize. The default 
start is the leftmost column without a recognized standard name. 
The default ending column is the last rightmost column. Indexes are 
0-based.

=back

If successful, it will return the name of the file saved.

=back

=head2 Verifying Datasets

When working with row Features and collecting scores, the dataset 
from which you are collecting must be verified prior to collection. 
This ensures that the proper database adaptor is available and loaded, 
and that the dataset is correctly specified (otherwise nothing would be 
collected). This verification is normally performed transparently when 
you call L<get_score|Bio::ToolBox::Data::Feature/get_score>() or 
L<get_position_scores|Bio::ToolBox::Data::Feature/get_position_scores>().
However, datasets may be explicitly verified prior to calling the score 
methods. 

=over 4

=item verify_dataset($dataset)

=item verify_dataset($dataset, $database)

Pass the name of the dataset (GFF type or type:source) for a GFF3-based 
database, e.g. <Bio::DB::SeqFeature::Store>, or path and file name for a 
data file, e.g. Bam, BigWig, BigBed, or USeq file. If a separate database 
is being used, pass the name or opened database object as a second 
parameter. For more advance options, see 
L<Bio::ToolBox::db_helper/verify_or_request_feature_types>. 

The name of the verified dataset, with a prefix if necessary, is returned.

=back

=cut

use strict;
use Carp qw(carp cluck croak confess);

use Bio::ToolBox::data_helper qw(
	generate_data_structure
	verify_data_structure
	sort_data_structure
	gsort_data_structure
	splice_data_structure
);
use Bio::ToolBox::db_helper qw(
	get_new_feature_list 
	get_new_genome_list 
);
use Bio::ToolBox::file_helper qw(
	load_data_file
	write_data_file
	write_summary_data
);
use Bio::ToolBox::Data::common;
use Bio::ToolBox::utility;


1;



### Initialize

sub new {
	my $class = shift;
	my %args  = @_;
	$args{features} ||= $args{feature} || undef;
	
	# check for stream
	if (exists $args{stream} and $args{stream}) {
		my $ok;
		eval {require Bio::ToolBox::Data::Stream; $ok = 1;};
		if ($ok) {
			return Bio::ToolBox::Data::Stream->new(@_);
		}
		else {
			croak "Cannot load Bio::ToolBox::Data::Stream!";
		}
	}
	
	my $data;
	if (exists $args{file}) {
		# load from file
		# extensions and file existence will be taken care of
		$data = load_data_file($args{file}) or 
			croak "Cannot load file $args{file}!\n";
	}
	elsif (exists $args{db} and exists $args{features} ) {
		# generate new list
		if ($args{features} eq 'genome') {
			$data = get_new_genome_list(%args) or 
				croak "Cannot generate new genome list!\n";
		}
		else {
			$data = get_new_feature_list(%args) or 
				croak "Cannot generate new feature list!\n";
		}
	}
	else {
		# a new empty structure
		my @datasets;
		if (exists $args{datasets}) {
			@datasets = @{ $args{datasets} };
		}
		elsif (exists $args{columns}) {
			@datasets = @{ $args{columns} };
		}
		my $feature = $args{features} || 'feature';
		$data = generate_data_structure($feature, @datasets);
	}
	
	return bless $data, $class;
}



### File functions

# we don't have file read or open, see the new method instead

sub write_file {
	my $self = shift;
	my %args;
	if (scalar @_ == 1) {
		$args{filename} = $_[0];
	}
	else {
		%args = @_;
	}
	$args{data} = $self;
	return write_data_file(%args);
}

sub save {
	return shift->write_file(@_);
}

sub summary_file {
	my $self = shift;
	my %args = @_;
	$args{data} = $self;
	return write_summary_data(%args);
}

sub write_gff {
	my $self = shift;
	my %args = @_;
	$args{data} = $self;
	return convert_and_write_to_gff_file(%args);
}

sub reload_children {
	my $self = shift;
	my @files = @_;
	return unless @files;
	
	# prepare Stream
	my $ok;
	eval {require Bio::ToolBox::Data::Stream; $ok = 1;};
	unless ($ok) {
		carp "unable to load Bio::ToolBox::Data::Stream! $@";
		return;
	}
	
	# open first stream
	my $first = shift @files;
	my $Stream = Bio::ToolBox::Data::Stream->new(in => $first);
	unless ($Stream) {
		carp "unable to load first child file $first";
		return;
	}
	
	# destroy current data table and metadata
	$self->{data_table} = [];
	for (my $i = 0; $i < $self->number_columns; $i++) {
		delete $self->{$i};
	}
	
	# copy the metadata
	foreach (qw(program feature db bed gff ucsc headers number_columns last_row)) {
		# various keys
		$self->{$_} = $Stream->{$_};
	}
	for (my $i = 0; $i < $Stream->number_columns; $i++) {
		# column metadata
		my %md = $Stream->metadata($i);
		$self->{$i} = \%md;
	}
	my @comments = $Stream->comments;
	push @{$self->{other}}, @comments;
	
	# first row column headers
	$self->{data_table}->[0] = [ @{ $Stream->{data_table}->[0] } ];
	
	# load the data
	while (my $row = $Stream->next_row) {
		my $a = $row->row_values;
		$self->add_row($a);
	}
	$Stream->close_fh;
	
	# load remaining files
	foreach my $file (@files) {
		my $Stream = Bio::ToolBox::Data::Stream->new(in => $file);
		if ($Stream->number_columns != $self->number_columns) {
			confess "fatal error: child file $file has a different number of columns!";
		}
		while (my $row = $Stream->next_row) {
			my $a = $row->row_values;
			$self->add_row($a);
		}
		$Stream->close_fh;
	}
	$self->verify;
	
	# clean up
	unlink($first, @files);
	return $self->last_row;
}




### Data structure manipulation

sub verify {
	my $self = shift;
	return verify_data_structure($self);
}

sub sort_data {
	my $self = shift;
	my $index = shift;
	my $direction = shift || 'i';
	return unless exists $self->{$index}{name};
	return sort_data_structure($self, $index, $direction);
}

sub gsort_data {
	my $self = shift;
	return gsort_data_structure($self, $self->chromo_column, $self->start_column);
}

sub splice_data {
	my $self = shift;
	splice_data_structure($self, @_);
	if (exists $self->{db_connection}) {
		# re-open a new un-cached database connection
		$self->open_database(1);
	}
	return 1;
}

sub convert_gff {
	my $self = shift;
	my %args = @_;
	$args{data} = $self;
	return convert_genome_data_2_gff_data(%args);
}



### Column manipulation

sub column_values {
	my ($self, $column) = @_;
	return unless defined $column;
	return unless exists $self->{$column}{name};
	my @values = map {$self->value($_, $column)} (0 .. $self->last_row);
	return wantarray ? @values : \@values;
}

sub add_column {
	my ($self, $name) = @_;
	return unless $name;
	my $column = $self->{number_columns};
	
	# check for array of column data
	if (ref $name eq 'ARRAY') {
		if ($self->last_row > 1) {
			# table has existing data beyond column headers
			if ($self->last_row == (scalar @$name - 1)) {
				# same number of elements, add it the table
				$self->{$column} = {
					'name'  => $name->[0],
					'index' => $column,    
				};
				for (my $r = 0; $r <= $self->last_row; $r++) {
					$self->{data_table}->[$r][$column] = $name->[$r];
				}
			}
			else {
				# different number of elements
				cluck "array has different number of elements than Data table!\n"; 
				return;
			}
		}
		else {
			# table has no existing data
			$self->{$column} = {
				'name'  => $name->[0],
				'index' => $column,    
			};
			for (my $i = 0; $i < scalar @$name; $i++) {
				$self->value($i, $column, $name->[$i]);
			}
			$self->{last_row} = scalar @$name - 1;
			$self->{headers} = 1; # boolean to indicate the table now has headers
		}
	}
	elsif (ref $name eq '') {
		# just a name
		$self->{$column} = {
			'name'      => $name,
			'index'     => $column,
		};
		$self->{data_table}->[0][$column] = $name;
	}
	else {
		cluck "must pass a scalar value or array reference";
		return;
	}
	
	$self->{number_columns}++;
	delete $self->{column_indices} if exists $self->{column_indices};
	return $column;
}

sub copy_column {
	my $self = shift;
	my $index = shift;
	return unless defined $index;
	my $data = $self->column_values($index);
	my $new_index = $self->add_column($data);
	$self->copy_metadata($index, $new_index);
	return $new_index;
}



### Rows and Data access

sub add_row {
	my $self = shift;
	my @row_data;
	if ($_[0] and ref $_[0] eq 'ARRAY') {
		@row_data = @{ $_[0] };
	}
	else {
		@row_data = map {'.'} (1 .. $self->{number_columns});
	}
	if (scalar @row_data > $self->{number_columns}) {
		cluck "row added has more elements than columns in data structure!\n"; 
		splice @row_data, 0, $self->{number_columns};
	}
	until (scalar @row_data == $self->{number_columns}) {
		push @row_data, '.';
	}
	my $row_index = $self->{last_row} + 1;
	$self->{data_table}->[$row_index] = \@row_data;
	$self->{last_row}++;
	return $row_index;
}

sub delete_row {
	my $self = shift;
	my @deleted = sort {$b <=> $a} @_;
	while (@deleted) {
		my $d = shift @deleted;
		splice( @{ $self->{data_table} }, $d, 1);
		$self->{last_row}--;
	}
	return 1;
}

sub row_values {
	my ($self, $row)  = @_;
	my @data = @{ $self->{data_table}->[$row] };
	return wantarray ? @data : \@data;
}

sub value {
	my ($self, $row, $column, $value) = @_;
	return unless (defined $row and defined $column);
	if (defined $value) {
		$self->{data_table}->[$row][$column] = $value;
	}
	return $self->{data_table}->[$row][$column];
}

sub row_stream {
	my $self = shift;
	return Bio::ToolBox::Data::Iterator->new($self);
}

sub iterate {
	my $self = shift;
	my $code = shift;
	unless (ref $code eq 'CODE') {
		cluck "iterate_function() method requires a code reference!";
		return;
	}
	my $stream = $self->row_stream;
	while (my $row = $stream->next_row) {
		&$code($row);
	}
	return 1;
}


####################################################

package Bio::ToolBox::Data::Iterator;
use Bio::ToolBox::Data::Feature;

sub new {
	my ($class, $data) = @_;
	my %iterator = (
		'index'     => 1,
		'data'      => $data,
	);
	return bless \%iterator, $class;
}

sub next_row {
	my $self = shift;
	return if $self->{'index'} > $self->{data}->{last_row}; # no more
	my $i = $self->{'index'};
	$self->{'index'}++;
	return Bio::ToolBox::Data::Feature->new(
		'data'      => $self->{data},
		'index'     => $i,
	);	
}

sub row_index {
	my $self = shift;
	return $self->{'index'};
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
