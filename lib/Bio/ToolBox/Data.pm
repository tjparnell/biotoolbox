package Bio::ToolBox::Data;
our $VERSION = 1.14;

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
  
=head1 PREFACE

This is an object-oriented interface to the Bio::ToolBox Data structure. 
Most of the remaining Bio::ToolBox libraries are collections of 
exported subroutines, and the data structures required a lot of 
manual manipulation (and redundant code - sigh). They were written 
before I learned to fully appreciate the benefits of OO-code. Hence, 
this module is an attempt to right the wrongs of my early practices.

Many of the provided scripts that accompany the Bio::ToolBox distribution 
do not use the this OO interface. They can be cryptic, obtuse, and hard 
to follow. New scripts should follow this interface instead.

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
may be obtained from an annotation database, an existing file 
may be loaded, or a new empty structure may be generated. 

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
file or directory may be given. If you already have an opened 
Bio::DB::SeqFeature::Store database object, you can simply 
pass that. See Bio::ToolBox::db_helper for more information.

=item win =E<gt> $integer

=item step =E<gt> $integer

If generating a list of genomic intervals, optionally provide 
the window and step values. Default values are defined in 
the Bio::ToolBox configuration file C<.biotoolbox.cfg>.

=back

If successful, the method will return a Bio::ToolBox::Data object.

=back

=head2 General Metadata

There is a variety of general metadata regarding the Data 
structure. 

The following methods may be used to access or set these 
metadata properties.  

=over

=item feature($text)

Returns or sets the name of the features used to collect 
the list of features. The actual feature types are listed 
in the table, so this is metadata is merely descriptive.

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
actual information. It has some metadata as well.

=over 4

=item number_columns

Returns the number of columns in the Data table. 

=item last_row

Returns the array index number of the last row. 
Since the header row is index 0, this is also the 
number of actual content rows.

=item add_column($name)

Appends a new empty column to the Data table at the 
rightmost position (highest index). It adds the column 
header name and creates a new column metadata hash. 
But it doesn't actually fill values in every row. 
It returns the new column index.

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

=item metadata($index, $key)

=item metadata($index, $key, $new_value)

Returns or sets the metadata value for a specific key for a 
specific column index.

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

=back

=head2 Bio::ToolBox::Data::Feature Methods

These are methods for working with the current data row generated 
using the next_row() method from a row_stream iterator.

=over 4

=item seq_id

=item start

=item end

=item strand

=item name

=item type

=item id

These methods return the corresponding appropriate value, if 
present. These rely on the corresponding find_column methods.

=item value($index)

=item value($index, $new_value)

Returns or sets the value at a specific column index in the 
current data row.

=back

The next three functions are convenience methods for using the 
attributes in the current data row to interact with databases. 
They are wrappers to methods in the Bio::ToolBox::db_helper 
module.

=over 4

=item feature

Returns a SeqFeature object from the database using the name and 
type values in the current Data table row. The SeqFeature object 
is requested from the database named in the general metadata. If 
an alternate database is desired, you should change it first using  
the general database() method. If the feature name or type is not 
present in the table, then nothing is returned.

=item segment

Returns a database Segment object corresponding to the coordinates 
defined in the Data table row. The database named in the general 
metadata is used to establish the Segment object. If a different 
database is desired, it should be changed first using the 
general database() method. 

=item get_score(%args)

This is a convenience method for the 
Bio::ToolBox::db_helper::get_chromo_region_score() method. It 
will return a single score value for the region defined by the 
coordinates or typed named feature in the current data row. If 
the Data table has coordinates, then those will be automatically 
used. If the Data table has typed named features, then the 
coordinates will automatically be looked up for you by requesting 
a SeqFeature object from the database.

The name of the dataset from which to collect the data must be 
provided. This may be a GFF type in a SeqFeature database, a 
BigWig member in a BigWigSet database, or a path to a BigWig, 
BigBed, Bam, or USeq file. Additional parameters may also be 
specified; please see the Bio::ToolBox::db_helper::
get_chromo_region_score() method for full details.

Here is an example of collecting mean values from a BigWig 
and adding the scores to the Data table.
  
  my $index = $Data->add_column('MyData');
  my $stream = $Data->row_stream;
  while (my $row = $stream->next_row) {
     my $score = $row->get_score(
        'method'    => 'mean',
        'dataset'   => '/path/to/MyData.bw',
     );
     $row->value($index, $score);
  }

=item get_position_scores(%args)

This is a convenience method for the Bio::ToolBox::db_helper::
get_region_dataset_hash() method. It will return a hash of 
positions =E<gt> scores over the region defined by the 
coordinates or typed named feature in the current data row. 
The coordinates for the interrogated region will be 
automatically provided.

Just like the get_score() method, the dataset from which to 
collect the scores must be provided, along with any other 
optional arguments. See the documentation for the 
Bio::ToolBox::db_helper::get_region_dataset_hash() method 
for more details.

Here is an example for collecting positioned scores around 
the 5 prime end of a feature from a BigWigSet directory.
  
  my $stream = $Data->row_stream;
  while (my $row = $stream->next_row) {
     my %position2score = $row->get_position_scores(
        'ddb'       => '/path/to/BigWigSet/',
        'dataset'   => 'MyData',
        'position'  => 5,
        'start'     => -500,
        'stop'      => 500,
     )
     # do something with %position2score
  }

=back

=head2 Data Table Functions

These methods alter the Data table en masse. 

=over 4

=item verify

This method will verify the Data structure, including the metadata and the 
Data table. It ensures that the table has the correct number of rows and 
columns as described in the metadata, and that each column has the basic 
metadata. 

If the Data structure is marked as a GFF or BED structure, then the table 
is checked that the structure matches the proper format. If not, for 
example when additional columns have been added, then the GFF or BED value 
is set to null. 

This method is automatically called prior to writing the Data table to file.

=item splice_data($current_part, $total_parts)

This method will splice the Data table into $total_parts number of pieces, 
retaining the $current_part piece. The other parts are discarded. This 
method is intended to be used when a program is forked into separate 
processes, allowing each child process to work on a subset of the original 
Data table. 

Two values are passed to the method. The first is the current part number, 
1-based. The second value is the total number of parts that the table 
should be divided, corresponding to the number of concurrent processes. 
For example, to fork the program into four concurrent processes.
	
	my $Data = Bio::ToolBox::Data->new(file => $file);
	my $pm = Parallel::ForkManager->new(4);
	for my $i (1..4) {
		$pm->start and next;
		### in child
		$Data->splice_data($i, 4);
		$db = $Data->open_database; # a clone-safe new db object
		# do something with this portion
		$Data->save('filename' => "file#$i");
		$pm->finish;
	}
	$pm->wait_all_children;

There is no convenient method for merging the modified contents of the 
table from each child process back into the original Data table, as 
each child is essentially isolated from the parent. The Parallel::ForkManager 
documentation recommends going through a disk file intermediate. See the 
accompanying BioToolBox script F<join_data_file.pl> for concatenating Data 
table files together.

Remember that if you fork your script into child processes, any database 
connections must be re-opened; they are typically not clone safe. If you 
have an existing database connection by using the open_database() method, 
it should be automatically re-opened for you when you use the splice_data() 
method, but you will need to call open_database() again in the child 
process to obtain the new database object.

=item convert_gff(%options)

This method will irreversibly convert the Data table into a GFF format.
Table columns will be added, deleted, reordered, and renamed as necessary 
to generate the GFF structure. An array of options should be passed to 
control the conversion step. 

=over 4

=item version =E<gt> <2|3>

Provide the GFF version. The default is version 3.

=item chromo =E<gt> $index

=item start =E<gt> $index

=item stop =E<gt> $index

=item strand =E<gt> $index

Provide the column indices for the appropriate columns. These should  
be automatically identified from the column header names. Indices 
are 0-based.

=item score =E<gt> $index

Provide the index column name for whatever score column. This is 
not automatically determined. 

=item source =E<gt> $index|$text

=item type =E<gt> $index|$text

=item name =E<gt> $index|$text

Provide either a column index (0-based) or a text name to be used 
for all the features. Integers between 0 and the rightmost column 
index are presumed to be an index; everything else is taken as text.

=item tag =E<gt> \@indices

Provide an array reference of column indices to be used for GFF tags.

=item id =E<gt> $index

Provide a column index of unique values to be used for GFF3 ID tag.

=item midpoint =E<gt> <boolean>

Flag to use the midpoint instead of actual start and stop coordinates.

=back

=back

=head2 Data Table File Functions

When you are finished modifying the Data table, it may then be written out 
as a tabbed-delimited text file. If the format corresponds to a valide BED or 
GFF file, then it may be written in that format. 

Several functions are available for writing the Data table, exporting to a 
compatible GFF file format, or writing a summary of the Data table.

=over 4

=item write_file()

=item save()

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

=item summary_file()

Write a separate file summarizing columns of data (mean values). 
The mean value of each column becomes a row value, and each column 
header becomes a row identifier (i.e. the table is transposed). The 
best use of this is to summarize the mean profile of windowed data 
collected across a feature. See the Bio::ToolBox scripts 
C<get_relative_data.pl> and C<average_gene.pl> as an example. 
You may pass options. 

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

=item write_gff()

This will write out the existing data in GFF format. A number of 
options may be passed to control the conversion.

=over 4

=item filename =E<gt> $filename

Optionally pass the filename to save. A suitable default will be 
generated if not provided.

=item version =E<gt> <2|3>

Provide the GFF version. The default is version 3.

=item chromo =E<gt> $index

=item start =E<gt> $index

=item stop =E<gt> $index

=item strand =E<gt> $index

Provide the column indices for the appropriate columns. These should  
be automatically identified from the column header names. Indices 
are 0-based.

=item score =E<gt> $index

Provide the index column name for whatever score column. This is 
not automatically determined. 

=item source =E<gt> $index|$text

=item type =E<gt> $index|$text

=item name =E<gt> $index|$text

Provide either a column index (0-based) or a text name to be used 
for all the features. Integers between 0 and the rightmost column 
index are presumed to be an index; everything else is taken as text.

=item tag =E<gt> \@indices

Provide an array reference of column indices to be used for GFF tags.

=item id =E<gt> $index

Provide a column index of unique values to be used for GFF3 ID tag.

=item midpoint =E<gt> <boolean>

Flag to use the midpoint instead of actual start and stop coordinates.

=back

=back

=cut

use strict;
use Carp qw(carp cluck croak confess);

use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure
	verify_data_structure
	splice_data_structure
	index_data_table
	find_column_index
);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	write_tim_data_file
	convert_genome_data_2_gff_data 
	convert_and_write_to_gff_file
	write_summary_data
	parse_filename
);
use Bio::ToolBox::db_helper qw(
	get_new_feature_list 
	get_new_genome_list 
	open_db_connection
	verify_or_request_feature_types
);

1;



### Initialize

sub new {
	my $class = shift;
	my %args  = @_;
	$args{features} ||= $args{feature} || 'feature';
	
	my $data;
	if (exists $args{file} and -e $args{file}) {
		# load from file
		$data = load_tim_data_file($args{file}) or 
			croak "Cannot load file $args{file}!\n";
	}
	elsif (exists $args{db} and exists $args{features}) {
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
		$data = generate_tim_data_structure($args{features}, @datasets);
	}
	
	return bless $data, $class;
}



### File functions

# we don't have file read or open, see the new method instead

sub write_file {
	my $self = shift;
	my %args = @_;
	$args{data} = $self;
	return write_tim_data_file(%args);
}

sub save {
	my $self = shift;
	return $self->write_file(@_);
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





### Data structure manipulation

sub verify {
	my $self = shift;
	return verify_data_structure($self);
}

sub splice_data {
	my $self = shift;
	splice_data_structure($self, @_);
	if ($self->{db} and exists $self->{db_connection}) {
		# re-open a new un-cached database connection
		my $db = open_db_connection($self->{db}, 1);
		$self->{db_connection} = $db if $db;
	}
	return 1;
}

sub convert_gff {
	my $self = shift;
	my %args = @_;
	$args{data} = $self;
	return convert_genome_data_2_gff_data(%args);
}




### Metadata

sub feature {
	my $self = shift;
	if (@_) {
		$self->{feature} = shift;
	}
	return $self->{feature};
}

sub program {
	my $self = shift;
	if (@_) {
		$self->{program} = shift;
	}
	return $self->{program};
}

sub database {
	my $self = shift;
	if (@_) {
		$self->{db} = shift;
		if (exists $self->{db_connection}) {
			my $db = open_db_connection($self->{db});
			$self->{db_connection} = $db if $db;
		}
	}
	return $self->{db};
}

sub open_database {
	my $self = shift;
	return unless $self->{db};
	if (exists $self->{db_connection}) {
		return $self->{db_connection};
	}
	else {
		my $db = open_db_connection($self->{db});
		return unless $db;
		$self->{db_connection} = $db;
		return $db;
	}
}

sub gff {
	my $self = shift;
	return $self->{gff};
}

sub bed {
	my $self = shift;
	return $self->{bed};
}

sub number_columns {
	my $self = shift;
	return $self->{number_columns};
}

sub last_row {
	my $self = shift;
	return $self->{last_row};
}

sub filename {
	my $self = shift;
	if (@_) {
		my $filename = shift;
		my ($basename, $path, $extension) = parse_filename($filename);
		$self->{filename}  = $filename;
		$self->{basename}  = $basename;
		$self->{path}      = $path;
		$self->{extension} = $extension;
	}
	return $self->{filename};
}

sub basename {
	my $self = shift;
	return $self->{basename};
}

sub path {
	my $self = shift;
	return $self->{path};
}

sub extension {
	my $self = shift;
	return $self->{extension};
}

sub comments {
	my $self = shift;
	my @comments = @{ $self->{other} };
	return @comments;
}

sub add_comment {
	my $self = shift;
	my $comment = shift or return;
	push @{ $self->{other} }, $comment;
	return 1;
}

sub delete_comment {
	my $self = shift;
	my $index = shift;
	if (defined $index) {
		eval {splice @{$self->{other}}, $index, 1};
	}
	else {
		$self->{other} = [];
	}
}

sub metadata {
	my $self = shift;
	my ($index, $key, $value) = @_;
	return unless defined $index;
	return unless exists $self->{$index};
	if ($key and $value) { 
		# we are setting a new value
		$self->{$index}{$key} = $value;
		return 1;
	}
	elsif ($key and not defined $value) {
		# retrieve a value
		return unless exists $self->{$index}{$key};
		return $self->{$index}{$key};
	}
	else {
		my %hash = %{ $self->{$index} };
		return wantarray ? %hash : \%hash;
	}
}




### Columns or Datasets

sub find_column {
	my ($self, $name) = @_;
	return unless $name;
	return find_column_index($self, $name);
}

sub add_column {
	my ($self, $name) = @_;
	return unless $name;
	my $column = $self->{number_columns};
	$self->{$column} = {
		'name'      => $name,
		'index'     => $column,
	};
	$self->{data_table}->[0][$column] = $name;
	$self->{number_columns}++;
	return $column;
}

sub delete_column {
	my $self = shift;
	
	my @deletion_list = sort {$a <=> $b} @_;
	my @retain_list; 
	for (my $i = 0; $i < $self->number_columns; $i++) {
		# compare each current index with the first one in the list of 
		# deleted indices. if it matches, delete. if not, keep
		if ( $i == $deletion_list[0] ) {
			# this particular index should be deleted
			shift @deletion_list;
		}
		else {
			# this particular index should be kept
			push @retain_list, $i;
		}
	}
	return $self->reorder_column(@retain_list);
}

sub reorder_column {
	my $self = shift;
	my @order = @_;
	
	# reorder data table
	for (my $row = 0; $row <= $self->last_row; $row++) {
		my @old = $self->row_values($row);
		my @new = map { $old[$_] } @order;
		splice( @{ $self->{data_table} }, $row, 1, \@new);
	}
	
	# reorder metadata
	my %old_metadata;
	for (my $i = 0; $i < $self->number_columns; $i++) {
		# copy the metadata info hash into a temporary hash
		$old_metadata{$i} = $self->{$i};
		delete $self->{$i}; # delete original
	}
	for (my $i = 0; $i < scalar(@order); $i++) {
		# now copy back from the old_metadata into the main data hash
		# using the new index number in the @order array
		$self->{$i} = { %{ $old_metadata{ $order[$i] } } };
		# assign new index number
		$self->{$i}{'index'} = $i;
	}
	$self->{'number_columns'} = scalar @order;
	return 1;
}

sub _find_column_indices {
	my $self = shift;
	my $name   = find_column_index($self, '^name');
	my $type   = find_column_index($self, '^type|class');
	my $id     = find_column_index($self, '^primary_id');
	my $chromo = find_column_index($self, '^chr|seq|ref|ref.?seq');
	my $start  = find_column_index($self, '^start|position');
	my $stop   = find_column_index($self, '^stop|end');
	my $strand = find_column_index($self, '^strand');
	$self->{column_indices} = {
		'name'      => $name,
		'type'      => $type,
		'id'        => $id,
		'seq_id'    => $chromo,
		'chromo'    => $chromo,
		'start'     => $start,
		'stop'      => $stop,
		'end'       => $stop,
		'strand'    => $strand,
	};
	return 1;
}

sub chromo_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{chromo};
}

sub start_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{start};
}

sub stop_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{stop};
}

sub strand_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{strand};
}

sub name_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{name};
}

sub type_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{type};
}

sub id_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{id};
}

sub verify_dataset {
	my $self = shift;
	my $dataset = shift;
	my $database = shift; # name or object?
	return unless $dataset;
	$database ||= $self->open_database;
	if (exists $self->{verfied_dataset}{$dataset}) {
		return $self->{verfied_dataset}{$dataset};
	}
	else {
		if ($dataset =~ /^(?:file|http|ftp)/) {
			# local or remote file already verified?
			$self->{verfied_dataset}{$dataset} = $dataset;
			return $dataset;
		}
		my ($verified) = verify_or_request_feature_types(
			# normally returns an array of verified features, we're only checking one
			db      => $database,
			feature => $dataset,
		);
		if ($verified) {
			$self->{verfied_dataset}{$dataset} = $verified;
			return $verified;
		}
	}
	return;
}



### Rows and Data access

sub add_row {
	my $self = shift;
	my $row_data = shift;
	if (scalar @$row_data > $self->{number_columns}) {
		cluck "row added has more elements than columns in data strucure!\n"; 
		splice @$row_data, 0, $self->{number_columns};
	}
	until (scalar @$row_data == $self->{number_columns}) {
		push @$row_data, undef;
	}
	my $row = $self->{last_row} + 1;
	$self->{data_table}->[$row] = $row_data;
	$self->{last_row}++;
	return 1;
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





####################################################

package Bio::ToolBox::Data::Iterator;

sub new {
	my ($class, $data) = @_;
	my $chromo = $data->chromo_column;
	my $start  = $data->start_column;
	my $type   = $data->type_column;
	my $name   = $data->name_column;
	my $feature_type;
	if (defined $chromo and defined $start) {
		$feature_type = 'coordinate';
	}
	elsif (defined $type and defined $name) {
		$feature_type = 'named';
	}
	my %iterator = (
		'index'     => 1,
		'data'      => $data,
		'last'      => $data->{last_row},
		'feature'   => $feature_type,
	);
	return bless \%iterator, $class;
}

sub next_row {
	my $self = shift;
	return if $self->{'index'} > $self->{'last'}; # no more
	my $i = $self->{'index'};
	$self->{'index'}++;
	return Bio::ToolBox::Data::Feature->new(
		'data'      => $self->{data},
		'index'     => $i,
		'feature'   => $self->{feature},
	);	
}

sub row_index {
	my $self = shift;
	return $self->{'index'};
}



####################################################

package Bio::ToolBox::Data::Feature;

use Carp;
use Bio::ToolBox::db_helper qw(
	get_feature
	get_chromo_region_score
	get_region_dataset_hash
);

sub new {
	my $class = shift;
	my %self = @_;
	return bless \%self, $class;
}


### Set and retrieve values

sub row_index {
	my $self = shift;
	return $self->{'index'};
}

sub row_values {
	my $self  = shift;
	my $row = $self->{'index'};
	my @data = @{ $self->{data}->{data_table}->[$row] };
	return wantarray ? @data : \@data;
}

sub value {
	my ($self, $column, $value) = @_;
	return unless defined $column;
	my $row = $self->{'index'};
	
	if (defined $value) {
		# set a value
		$self->{data}->{data_table}->[$row][$column] = $value;
	}
	return $self->{data}->{data_table}->[$row][$column];
}

sub seq_id {
	my $self = shift;
	my $i = $self->{data}->chromo_column;
	return defined $i ? $self->value($i) : undef;
}

sub start {
	my $self = shift;
	my $i = $self->{data}->start_column;
	return defined $i ? $self->value($i) : undef;
}

sub end {
	my $self = shift;
	my $i = $self->{data}->stop_column;
	return defined $i ? $self->value($i) : undef;
}

sub stop {
	return shift->end;
}

sub strand {
	my $self = shift;
	my $i = $self->{data}->strand_column;
	return defined $i ? $self->value($i) : 0; # default is no strand
}

sub name {
	my $self = shift;
	my $i = $self->{data}->name_column;
	return defined $i ? $self->value($i) : undef;
}

sub type {
	my $self = shift;
	my $i = $self->{data}->type_column;
	return defined $i ? $self->value($i) : undef;
}

sub id {
	my $self = shift;
	my $i = $self->{data}->id_column;
	return defined $i ? $self->value($i) : undef;
}


### Data collection convenience methods

sub feature {
	my $self = shift;
	my $id   = $self->id;
	my $name = $self->name;
	my $type = $self->type;
	return unless ($name and $type);
	return unless $self->{data}->database;
	return get_feature(
		'db'    => $self->{data}->open_database,
		'id'    => $self->id,
		'name'  => $self->name,
		'type'  => $self->type,
	);
}

sub segment {
	my $self   = shift;
	my $chromo = $self->seq_id;
	my $start  = $self->start;
	my $stop   = $self->end;
	return unless ($chromo and defined $start and defined $stop);
	return unless $self->{data}->database;
	my $db = $self->{data}->open_database;
	return $db ? $db->segment($chromo, $start, $stop) : undef;
}

sub get_score {
	my $self = shift;
	my %args = @_;
	
	if ($self->{feature} eq 'named') {
		# we must get the coordinates ourselves via lookup
		my $f = $self->feature;
		return unless $f;
		$args{chromo} ||= $f->seq_id;
		$args{start}  ||= $f->start;
		$args{stop}   ||= $f->end;
		$args{strand} ||= $f->strand;
	}
	elsif ($self->{feature} eq 'coordinate') {
		$args{chromo} ||= $self->seq_id;
		$args{start}  ||= $self->start;
		$args{stop}   ||= $self->end;
		$args{strand} ||= $self->strand;
	}
	else {
		# die otherwise we will have this error every time
		croak "data table does not have coordinates or feature attributes for score collection\n";
	}
	
	# verify the dataset for the user, cannot trust whether it has been done or not
	my $db = $args{ddb} || $args{db} || $self->{data}->open_database || undef;
	$args{dataset} = $self->{data}->verify_dataset($args{dataset}, $db);
	unless ($args{dataset}) {
		croak "provided dataset was unrecognized format or otherwise could not be verified!\n";
	}
	
	$args{db} ||= $self->{data}->open_database;
	
	return get_chromo_region_score(%args);
}

sub get_position_scores {
	my $self = shift;
	my %args = @_;
	
	if ($self->{feature} eq 'named') {
		# we must get the coordinates ourselves via lookup
		$args{id}     ||= $self->id;
		$args{name}   ||= $self->name;
		$args{type}   ||= $self->type;
	}
	elsif ($self->{feature} eq 'coordinate') {
		$args{chromo} ||= $self->seq_id;
		$args{start}  ||= $self->start;
		$args{stop}   ||= $self->end;
		$args{strand} ||= $self->strand;
	}
	else {
		# die otherwise we will have this error every time
		croak "data table does not have coordinates or feature attributes for score collection\n";
	}
	
	# verify the dataset for the user, cannot trust whether it has been done or not
	my $db = $args{ddb} || $args{db} || $self->{data}->open_database;
	$args{dataset} = $self->{data}->verify_dataset($args{dataset}, $db);
	unless ($args{dataset}) {
		croak "provided dataset was unrecognized format or otherwise could not be verified!\n";
	}
	
	$args{db} ||= $self->{data}->open_database;
	unless ($args{db} or $args{ddb}) {
		# make sure we provide a database
		if ($self->{feature} eq 'coordinate' and $args{dataset} =~ /^(?:file|http|ftp)/) {
			$args{ddb} = $args{dataset};
		}
	}
	
	return get_region_dataset_hash(%args);
}


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
