package Bio::ToolBox::Data;
our $VERSION = '1.69';

=head1 NAME

Bio::ToolBox::Data - Reading, writing, and manipulating data structure

=head1 SYNOPSIS

  use Bio::ToolBox::Data;
  
  ### Open a pre-existing file
  # a data table with same columns as input file
  my $Data = Bio::ToolBox::Data->new(
        file    => 'coordinates.bed',
  );
  
  ### Parse a GTF, GFF3, refFlat, or genePred gene table
  # data table with names and references to fully parsed 
  # SeqFeature transcript objects with exon SeqFeature subfeatures
  my $Data = Bio::ToolBox::Data->new(
  	    file    => 'annotation.gtf.gz',
  	    parse   => 1,
  	    feature => 'transcript',
  	    subfeature => 'exon'
  );
  
  ### New gene list from a Bio::DB::SeqFeature::Store database
  # data table with name and reference ID to database SeqFeature objects
  my $Data = Bio::ToolBox::Data->new(
        db      => 'hg19.sqlite',
        feature => 'gene:ensGene',
  );
    
  
  ### Get a specific value
  my $value = $Data->value($row, $column);
  
  ### Replace or add a value
  $Data->value($row, $column, $new_value);
  
  ### Add a column
  my $new_index = $Data->add_column('Data2');
  
  ### Find a column
  my $old_index = $Data->find_column('Data1');
  
  ### Return a specific row as an object
  my $row = $Data->get_row($i); # Bio::ToolBox::Data::Feature
  
  ### Iterate through a Data structure one row at a time
  my $stream = $Data->row_stream;
  while (my $row = $stream->next_row) {
  	  # each row is returned as Bio::ToolBox::Data::Feature object
  	  # get the positional information from the file data
  	  # assuming that the input file had these identifiable columns
  	  my $seq_id = $row->seq_id;
  	  my $start  = $row->start;
  	  my $stop   = $row->end;
  	  
  	  # work with the referenced SeqFeature object
  	  my $seqf = $row->seqfeature;
  	  
  	  # generate a new Bio::Seq object from the database using 
  	  # these coordinates 
  	  my $region = $db->segment($seq_id, $start, $stop);
  	  
  	  # modify a row value
  	  my $value = $row->value($old_index);
  	  my $new_value = $value + 1;
  	  $row->value($new_index, $new_value);
  }
  
  ### Iterate through a Data table with a code reference
  $Data->iterate(\&my_code);
  
  
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
file formats are simply tabbed-delimited text files (think BED, GFF, VCF). 
Each row is a feature or genomic interval, and each column is a piece of 
information about that feature, such as name, type, and/or coordinates. 
We can append to that file additional columns of information, perhaps 
scores from genomic data sets. We can record metadata regarding how 
and where we obtained that data. Finally, we can write the updated 
table to a new file.

=head1 METHODS

=head2 Initializing the structure

=over 4

=item new

Initialize a new Data structure. This generally requires options, 
provided as an array of key =E<gt> values. A new list of features 
may be obtained from an annotation database or an existing file 
may be loaded. If you do not pass any options, a new empty 
structure will be generated for you to populate. 

These are the options available.

=over 4

=item file

=item in

  my $Data = Bio::ToolBox::Data->new(file => $file);

Provide the path and name to an existing tabbed-delimited text 
file from which to load the contents. This is a shortcut to the 
load_file() method. See that method for more details.

=item stream

  my $Data = Bio::ToolBox::Data->new(file => $file, stream => 1);

Boolean option indicating that the file should be opened as a file  
stream. A L<Bio::ToolBox::Data::Stream> object will be returned. This 
is a convenience method. 

=item noheader

  my $Data = Bio::ToolBox::Data->new(file => $file, noheader => 1);

Boolean option indicating that the file does not have file headers, 
in which case dummy headers are provided. This is not necessary for 
defined file types that don't normally have file headers, such as 
BED, GFF, or UCSC files.

=item parse

  my $Data = Bio::ToolBox::Data->new(file => $file, parse => 1);

Boolean option indicating that a gene annotation table or file should 
be parsed into SeqFeature objects and a general table of names and IDs 
representing those objects be generated. The annotation file may 
be specified in one of two ways: Through the file option above, 
or in the database metadata of an existing table file representing 
previously parsed objects.

=item db

  my $Data = Bio::ToolBox::Data->new(db => 'hg19', feature => 'gene');

Provide the name of the database from which to collect the 
features. It may be a short name, whereupon it is checked in 
the L<Bio::ToolBox> configuration file F<.biotoolbox.cfg> for 
connection information. Alternatively, a path to a database 
file or directory may be given. 

If you already have an opened L<Bio::DB::SeqFeature::Store> database 
object, you can simply pass that. See L<Bio::ToolBox::db_helper> for 
more information. However, this in general should be discouraged, 
since the name of the database will not be properly recorded when 
saving to file. It is better to simply pass the name of database 
again; multiple connections to the same database are smartly handled 
in the background.

=item feature

  my $Data = Bio::ToolBox::Data->new(file => $filename, feature => 'gene');

For de novo lists from an annotation database, provide the GFF 
type or type:source (columns 3 and 2) for collection. A comma 
delimited string may be accepted for databases. For parsed files, 
only a simple string is accepted.

For a list of genomic intervals across the genome, specify a 
feature of 'genome' with a database object.

=item subfeature

When parsing annotation files such as GTF, one or more subfeature 
types, e.g. C<exon> or C<utr>, may be specified as a comma-delimited 
string. This ensures that the subfeatures will be parsed into 
SeqFeature objects. Otherwise, only the top level features will be 
parsed. This expedites parsing by skipping unwanted features.

=item win

=item step

  my $Data = Bio::ToolBox::Data->new(db => $dbname, win => $integer);

If generating a list of genomic intervals, optionally provide 
the window and step values. The default is 500 bp.

=item chrskip

Provide a regular expression compatible or C<qr> string for skipping or 
excluding specific or classes of chromosomes, for example the mitochondrial 
chromosome or unmapped contigs. This works with both feature collection 
and genomic intervals. The default is to take all chromosomes.

=back

When no file is given or database given to search, then a new, 
empty Data object is returned. In this case, you may optionally 
specify the names of the columns or indicate a specific file 
format structure to assign column names. The following options can 
then be provided.

=over 4

=item columns

=item datasets

  my $Data = Bio::ToolBox::Data->new(columns => [qw(Column1 Column2 ...)] );

Provide the column names in advance as an anonymous array. 

=item gff

Pass the GFF version of the file: 1, 2, 2.5 (GTF), or 3.

=item bed

Pass the number of BED columns (3-12).

=item ucsc 

Pass the number of columns to indicate the type of UCSC format. These 
include 10 (refFlat without gene names), 11 (refFlat with gene names), 
12 (knownGene gene prediction table), and 15 
(an extended gene prediction or genePredExt table).

=back

If successful, the method will return a Bio::ToolBox::Data object.

=item duplicate

This will create a new Data object containing the same column 
headers and metadata, but lacking the table content, i.e. no 
rows of data. File name metadata, if present in the original, is 
not preserved. The purpose here, for example, is to allow one 
to selectively copy rows from one Data object to another.

=item parse_table

  $Data->parse_table($file)
  $Data->parse_table( {
         file => $file, 
         feature => 'gene',
         subfeature => 'exon',
         chrskip => 'chrM|contig',
  } );

This will parse a gene annotation table into SeqFeature objects. 
If this is called from an empty Data object, then the table will 
be filled with the SeqFeature object names and IDs. If this is 
called from a non-empty Data object, then the table's contents 
will be associated with the SeqFeature objects using their name and 
ID. The stored SeqFeature objects can be retrieved using the 
L</get_seqfeature> method.

Pass the method a single argument. This may be either a simple 
scalar to a filename to be parsed, or a reference to hash of 
one or more argument options. Possible options include:

=over 4

=item file

The file name of the gene table to be parsed. This may 
be a GFF, GFF3, GTF, or any of the UCSC gene table formats. 
These will be parsed using Bio::ToolBox::parser::* adapters.

=item feature

A regular expression compatible string or C<qr> string to match 
the top features C<primary_tag> to keep. The C<source_tag> is not 
checked. The default is 'gene', although a transcript type could 
alternatively be specified (in which case transcripts are not 
parsed in gene features).

=item subfeature

A regular expression compatible string or C<qr> string to match 
any sub features C<primary_tag> to parse. The C<source_tag> is not checked.
Typically these include exon, CDS, or UTR. The default is nothing.

=item chrskip

A regular expression compatible string or C<qr> string to match 
chromosomes to be skipped or excluded. Any feature with a matching 
C<seq_id> chromosome name will be skipped.

=back

=back

=head2 General Metadata

There is a variety of general metadata regarding the Data 
structure. The following methods may be used to access or set these 
metadata properties. Some of these are stored as comment lines at 
the beginning of the file, and will be read or set when the file is 
loaded.

=over

=item feature

  $Data->feature('something');
  my $feature = $Data->feature;

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

=item * coordinate

Table includes at least chromosome and start columns.

=item * named

Table includes name, type, and/or Primary_ID, possibly 
referring to named database features.

=item * unknown

Table is unrecognized format.

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

=head2 Metadata comments

Comments are any other commented lines from a text file (lines 
beginning with a #) that were not parsed as general metadata.

=over 4

=item comments

Returns a copy of the array containing commented lines. Each 
comment line becomes an element in the array.

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
comment lines. Do this prior to writing or saving the Data sturcture
or else you will lose your changed VCF header metadata.

=back

=head2 The Data table

The Data table is the array of arrays containing all of the 
actual information. Rows and columns are indexed using 0-based 
indexing as with all Perl arrays. Row 0 is always the column 
header row containing the column names, regardless whether an 
actual header name existed in the original file format (e.g. 
BED or GFF formats). Any individual table "cell" can be 
specified as C<[$row][$column]>. 

=over 4

=item list_columns

Returns an array or array reference of the column names 
in ascending (left to right) order.

=item number_columns

Returns the number of columns in the Data table. 

=item last_column

Returns the array index number for the last (right most) 
column. This number is always 1 less than the value 
returned by number_columns() due to 0-based indexing.

=item last_row

Returns the array index number of the last row. 
Since the header row is index 0, this is also the 
number of actual content rows.

=item column_values

  my $values = $Data->column_values($i);

Returns an array or array reference representing the values 
in the specified column. This includes the column header as 
the first element. Pass the method the column index.

=item add_column

  # add empty column with name
  my $i = $Data->add_column($name);
  
  # add entire column
  my $i = $Data->add_column(\@array);

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

=item copy_column

  my $j = $Data->copy_column($i);

This will copy a column, appending the duplicate column at 
the rightmost position (highest index). It will duplicate 
column metadata as well. It will return the new index 
position.

=item delete_column

  $Data->delete_column($i, $j);

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

=item add_row

  $Data->add_row(\@values);
  $Data->add_row($Row); # Bio::ToolBox::Data::Feature object

Add a new row of data values to the end of the Data table. 
Optionally provide either a reference to an array of values 
to put in the row, or pass a L<Bio::ToolBox::Data::Feature> 
Row object, such as one obtained from another Data object. 
If the number of columns do not match, the array is filled 
up with null values for missing columns, or excess values 
are dropped.

=item get_row

  $Data->get_row($i); # Bio::ToolBox::Data::Feature object

Returns the specified row as a L<Bio::ToolBox::Data::Feature> 
object.

=item delete_row

  $Data->delete_row($i, $j);

Deletes one or more specified rows. Rows are spliced out 
highest to lowest index to avoid issues. Be very careful 
deleting rows while simultaneously iterating through the 
table!

=item row_values

  my $row_values = $Data->row_values($i);

Returns a copy of an array for the specified row index. 
Modifying this returned array does not migrate back to the 
Data table; Use the L</value> method below instead.

=item value

  my $value = $Data->value($row, $column);
  $Data->value($row, $column, $new_value);

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

=item name

  $Data->name($index, $new_name);
  my $name = $Data->name($i);

Convenient method to return the name of the column given the 
index number. A column may also be renamed by passing a new name.

=item metadata

  $Data->metadata($index, $key, $new_value);
  my $value = $Data->metadata($index, $key)

Returns or sets the metadata value for a specific $key for a 
specific column $index.

This may also be used to add a new metadata key. Simply provide 
the name of a new $key that is not present

If no key is provided, then a hash or hash reference is returned 
representing the entire metadata for that column.

=item copy_metadata

  $Data->copy_metadata($source, $target);

This method will copy the metadata (everything except name and 
index) between the source column and target column. Returns 1 if 
successful.  

=item delete_metadata

  $Data->delete_metadata($index, $key);

Deletes a column-specific metadata $key and value for a specific 
column $index. If a $key is not provided, then all metadata keys 
for that index will be deleted.

=item find_column

  my $i = $Data->find_column('Gene');
  my $i = $Data->find_column('^Gene$')

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

=head2 Working with databases

These are methods for working primarily with annotation databases.

=over 4

=item open_database

This is wrapper method that tries to do the right thing and passes 
on to either L</open_meta_database> or L</open_new_database> methods. 
Basically a legacy method for L</open_meta_database>.

=item open_meta_database

Open the database that is listed in the metadata. Returns the 
database connection. Pass a true value to force a new database 
connection to be opened, rather than returning a cached connection 
object (useful when forking).

=item open_new_database

Convenience method for opening a second or new database that is 
not specified in the metadata, useful for data collection. This 
is a shortcut to L<Bio::ToolBox::db_helper/open_db_connection>.
Pass the database name.

=back

=head2 SeqFeature Objects

SeqFeature objects corresponding to data rows can be stored in the Data 
object. This can be useful if the SeqFeature object is not readily 
available from a database or is processor intensive in generating or 
parsing. Note that storing large numbers of objects will increase memory 
usage. 

SeqFeature objects are usually either L<Bio::DB::SeqFeature>, 
L<Bio::SeqFeature::Lite>, or L<Bio::DB::SeqFeature> objects, depending 
upon their source. More information can obtained from the L<Bio::SeqFeatureI> 
abstract API.

=over 4

=item store_seqfeature

  $Data->store_seqfeature($row_index, $seqfeature);

Stores the SeqFeature object for the given row index. Only one SeqFeature 
object can be stored per row.

=item get_seqfeature

  my $feature = $Data->get_seqfeature($row_index);

Retrieves the SeqFeature object for the given row index.

=item delete_seqfeature

  $Data->store_seqfeature($row_index);

Removes the SeqFeature object for the given row index. 

=item collapse_gene_transcripts

  my $success = $Data->collapse_gene_transcripts;

This method will iterate through a Data table and collapse multiple alternative 
transcript subfeatures of stored gene SeqFeature objects in the table. Exons 
of multiple transcripts will be merged, maximizing exon size and minimizing 
intron size. Genes with only one transcript will not be affected. Stored 
SeqFeature objects that are do not have a C<primary_tag> of "gene" are silently 
skipped. Refer to the L<Bio::ToolBox::GeneTools/collapse_transcripts> method 
for more details. The number of rows successfully collapsed is returned. 

=item add_transcript_length

  my $length_index = $Data->add_transcript_length;
  my $length_index = $Data->add_transcript_length('cds');

This method will generate a new column in the Data table representing the 
length of a transcript, i.e. the sum of the length of subfeatures for 
the stored SeqFeature object in the Data table. The default subfeature is 
C<exon>; however, alternative subfeature types may be passed to the method 
and used. These include C<cds>, C<5p_utr>, and C<3p_utr> for CDS, the 5C<'> UTR, 
and the 3C<'> UTR, respectively. See the corresponding transcript length 
methods in L<Bio::ToolBox::GeneTools> for more information. If a length 
is not calculated, for example the feature C<primary_tag> is a "gene", 
then the simple length of the feature is recorded. 

The name of the new column is one of "Merged_Transcript_Length" for exon, 
"Transcript_CDS_Length", "Transcript_5p_UTR_Length", or "Transcript_3p_UTR_Length". 
The index of the new length column is returned. 

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

=item sort_data

  $Data->sort_data($index, 'i'); # increasing sort

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

=item splice_data

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

This method will splice the Data table into C<$total_parts> number of pieces, 
retaining the C<$current_part> piece. The other parts are discarded. This 
method is intended to be used when a program is forked into separate 
processes, allowing each child process to work on a subset of the original 
Data table. 

Two values are passed to the method. The first is the current part number, 
1-based. The second value is the total number of parts that the table 
should be divided, corresponding to the number of concurrent processes. 
One easy approach to forking is to use L<Parallel::ForkManager>. The 
above example shows how to fork into four concurrent processes.

Since each forked child process is separate from their parent process, 
their contents must be reloaded into the current Data object. The 
L<Parallel::ForkManager> documentation recommends going through a disk 
file intermediate. Therefore, write each child Data object to file using 
a unique name. Once all children have been reaped, they can be reloaded 
into the current Data object using the L</reload_children> method.

Remember that if you fork your script into child processes, any database 
connections must be re-opened; they are typically not clone safe. If you 
have an existing database connection by using the L</open_database> method, 
it should be automatically re-opened for you when you use the L</splice_data> 
method, but you will need to call L</open_database> again in the child 
process to obtain the new database object.

=item reload_children

Discards the current data table in memory and reloads two or more files 
written from forked children processes. Provide the name of the child 
files in the order you want them loaded. The files will be automatically 
deleted if successfully loaded. Returns the number of lines reloaded on 
success.

=back

=head2 File Functions

The Data table may be read in from a file or written out as a file. In 
all cases, it is a tab-delimited text file, whether as an ad hoc table 
or a specific bioinformatic format, e.g. BED, GFF, etc. Multiple common 
file formats are supported. Column headers are assumed, except in those 
cases where it is not, e.g. BED, GFF, etc. Metadata may be included as 
commented lines at the beginning of the file, prefixed with a C<#> symbol.
Reading and writing gzip compressed files is fully supported.

=over 4

=item load_file

  $Data->load_file($filename);

This will load a file into a new, empty Data table. This function is 
called automatically when a filename is provided to the L</new> function. 
The existence of the file is first checked (appending common missing 
extensions as necessary), metadata and column headers processed and/or 
generated from default settings, the content loaded into the table, and 
the structure verified. Error messages may be printed if the structure or 
format is inconsistent or doesn't match the expected format, e.g a file 
with a F<.bed> extension doesn't match the UCSC specification.
Pass the name of the filename.

=item taste_file

  my $flavor = $Data->taste_file($filename);
  # returns gff, bed, ucsc, or undef

Tastes, or checks, a file for a certain flavor, or known gene file formats. 
This is based on file extension, metadata headers, and/or file content 
in the first 10 lines or so. Returns a string based on the file format.
Values include gff, bed, ucsc, or undefined. Useful for determining if 
the file represents a known gene table format that lacks a defined file 
extension, e.g. UCSC formats.

=item save

=item write_file

  my $success = $Data->save;
  my $success = $Data->save('my_file.txt');
  my $success = $Data->save(filename => $file, gz => 1);
  print "file $success was saved!\n";

Pass the file name to be written. If no file name is passed, then 
the filename and path stored in the metadata are used, if present.

These methods will write the Data structure out to file. It will 
be first verified as to proper structure. Opened BED and GFF files 
are checked to see if their structure is maintained. If so, they 
are written in the same format; if not, they are written as regular 
tab-delimited text files. 

You may pass additional options.

=over 4

=item filename

Optionally pass a new filename. Required for new objects; previous 
opened files may be overwritten if a new name is not provided. If 
necessary, the file extension may be changed; for example, BED files 
that no longer match the defined format lose the .bed and gain a .txt 
extension. Compression may or add or strip .gz as appropriate. If 
a path is not provided, the current working directory is used.

=item gz

Boolean value to change the compression status of the output file. If 
overwriting an input file, the default is maintain the compression status, 
otherwise no compression. Pass a 0 for no compression, 1 for standard 
gzip compression, or 2 for block gzip (bgzip) compression for tabix 
compatibility.

=back

If the file save is successful, it will return the full path and 
name of the saved file, complete with any changes to the file extension.

=item summary_file

Write a separate file summarizing columns of data (mean values). 
The mean value of each column becomes a row value, and each column 
header becomes a row identifier (i.e. the table is transposed). The 
best use of this is to summarize the mean profile of windowed data 
collected across a feature. See the Bio::ToolBox scripts 
L<get_relative_data.pl> and L<get_binned_data.pl> as examples. 
For data from L<get_binned_data.pl> where the columns are expressed 
as percentile bins, the reported midpoint column is automatically 
converted based on a length of 1000 bp.

You may pass these options. They are optional.

=over 4

=item filename

Pass an optional new filename. The default is to take the basename 
and append "_summed" to it.

=item startcolumn

=item stopcolumn

Provide the starting and ending columns to summarize. The default 
start is the leftmost column without a recognized standard name. 
The default ending column is the last rightmost column. Indexes are 
0-based.

=item dataset

Pass a string that is the name of the dataset. This could be collected 
from the metadata, if present. This will become the name of the score 
column if defined.

=back

The name of the summarized column is either the provided dataset name, 
the defined basename in the metadata of the Data structure, or a generic 
name. If successful, it will return the name of the file saved.

=back

=head2 Verifying Datasets

When working with row Features and collecting scores, the dataset 
from which you are collecting must be verified prior to collection. 
This ensures that the proper database adaptor is available and loaded, 
and that the dataset is correctly specified (otherwise nothing would be 
collected). This verification is normally performed transparently when 
you call L<get_score|Bio::ToolBox::Data::Feature/get_score> or 
L<get_position_scores|Bio::ToolBox::Data::Feature/get_position_scores>.
However, datasets may be explicitly verified prior to calling the score 
methods. 

=over 4

=item verify_dataset

 my $dataset = $Data->verify_dataset($dataset, $database);

Pass the name of the dataset (GFF type or type:source) for a GFF3-based 
database, e.g. <Bio::DB::SeqFeature::Store>, or path and file name for a 
data file, e.g. Bam, BigWig, BigBed, or USeq file. If a separate database 
is being used, pass the name or opened database object as a second 
parameter. For more advance options, see 
L<Bio::ToolBox::db_helper/verify_or_request_feature_types>. 

The name of the verified dataset, with a prefix if necessary, is returned.

=back

=head2 Efficient Data Access

Most of the time we need to iterate over the Data table, one row 
at a time, collecting data or processing information. These methods 
simplify the process.

=over 4

=item iterate

    $Data->iterate( sub {
       my $row = shift;
       my $number = $row->value($index);
       my $log_number = log($number);
       $row->value($index, $log_number);
    } );

This method will process a code reference on every row in the data 
table. Pass a subroutine or code reference. The subroutine will 
receive the row as a L<Bio::ToolBox::Data::Feature> object. With this 
object, you can retrieve values, set values, and add new values. 

=item row_stream

This returns an C<Bio::ToolBox::Data::Iterator> object, which has one 
method, C<next_row()>. Call this method repeatedly until it returns 
C<undef> to work through each row of data.

Users of the C<Bio::DB> family of database adaptors may recognize the 
analogy to the C<seq_stream()> method.

=item next_row

  my $stream = $Data->row_stream;
  while (my $row = $stream->next_row) {
     # each $row is a Bio::ToolBox::Data::Feature object
     # representing the row in the data table
     my $value = $row->value($index);
     # do something with $value
  }

Called from a C<Bio::ToolBox::Data::Iterator> object, it returns a 
L<Bio::ToolBox::Data::Feature> row object. If SeqFeature objects are 
associated with the row, perhaps from a parsed input annotation file, 
then they are automatically associated with the row object. (Previous 
versions required separately calling the seqfeature() row method to 
perform this.)

=back

=head1 SEE ALSO

L<Bio::ToolBox::Data::Feature>, L<Bio::ToolBox::SeqFeature>, L<Bio::DB::Sam>,
L<Bio::DB::HTS>, L<Bio::DB::BigWig>, L<Bio::DB::BigBed>, L<Bio::DB::USeq>, 
L<Bio::DB::SeqFeature::Store>, L<Bio::Perl>

=cut

use strict;
use Carp qw(carp cluck croak confess);
use List::Util qw(sum0);
use base 'Bio::ToolBox::Data::core';
use Bio::ToolBox::db_helper qw(
	get_new_feature_list  
	get_new_genome_list
	get_db_feature
);
use Bio::ToolBox::utility qw(simplify_dataset_name sane_chromo_sort);
use Module::Load;

1;


### Initialize

sub new {
	my $class = shift;
	my %args  = @_;
	if (ref($class)) {
		$class = ref($class);
	}
	
	# check for important arguments
	$args{features} ||= $args{feature} || 'gene';
	$args{stream} ||= $args{Stream} || 0;
	$args{file} ||= $args{in} || undef;
	$args{parse} ||= 0;
	$args{noheader} ||= 0;
	
	# check for stream
	if ($args{stream}) {
		$class = "Bio::ToolBox::Data::Stream";
		load($class);
		return $class->new(@_);
	}
	
	# initialize
	my $self = $class->SUPER::new();
	
	# prepare a new table based on the provided arguments
	if ($args{file} and $args{parse}) {
		# parse from file
		$args{subfeature} ||= '';
		unless ( $self->parse_table(\%args) ) {
			my $l = $self->load_file($args{file});
			return unless $l;
			if ($self->database =~ /^Parsed:(.+)$/ and $self->feature_type eq 'named') {
				# looks like the loaded file was from a previously parsed table
				# let's try this again
				# this may die if it doesn't work
				$args{file} = $1;
				$args{feature} = $self->feature;
				$self->parse_table(\%args); 
			}
		}
	}
	elsif ($args{file}) {
		# load from file
		unless ( $self->load_file($args{file}, $args{noheader}) ) {
			return;
		}
	}
	elsif (exists $args{db} and $args{features}) {
		# generate new list
		$self->feature($args{features});
		$self->database($args{db});
		$args{data} = $self;
		my $result;
		if ($args{features} eq 'genome') {
			# create a genomic interval list
			# first must parse any exclusion list provided as a new Data object
			$args{blacklist} ||= $args{exclude} || undef;
			if (defined $args{blacklist}) {
				my $exclusion_Data = $self->new(
					file => $args{blacklist},
					parse => 0 # we don't need to parse
				);
				if ($exclusion_Data and $exclusion_Data->feature_type eq 'coordinate') {
					printf "   Loaded %d exclusion list items\n", 
						$exclusion_Data->number_rows;
					delete $args{blacklist};
					delete $args{exclude};
					$args{exclude} = $exclusion_Data;
				}
				else {
					carp " Cannot not load exclusion coordinate list file $args{blacklist}!";
					return;
				}
			}
			$result = get_new_genome_list(%args);
		}
		else {
			$result = get_new_feature_list(%args);
		}
		unless ($result) {
			carp " Cannot generate new $args{features} list!";
			return;
		}
	}
	else {
		# a new empty structure
		
		# check to see if user provided column names
		$args{columns} ||= $args{datasets} || undef;
		if (defined $args{columns}) {
			 foreach my $c ( @{ $args{columns} } ) {
			 	$self->add_column($c);
			 }
			 $self->{feature} = $args{feature} if exists $args{feature};
		}
		
		# or possibly a specified format structure 
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
	}
	
	return $self;
}

sub duplicate {
	my $self = shift;
	
	# duplicate the data structure
	my $columns = $self->list_columns;
	my $Dupe = $self->new(
		'columns' => $columns,
	) or return;
	
	# copy the metadata
	for (my $i = 0; $i < $self->number_columns; $i++) {
		# column metadata
		my %md = $self->metadata($i);
		$Dupe->{$i} = \%md;
	}
	foreach (qw(feature program db bed gff ucsc vcf headers)) {
		# various keys
		$Dupe->{$_} = $self->{$_};
	}
	my @comments = $self->comments;
	push @{$Dupe->{comments}}, @comments;
	
	return $Dupe;
}

sub parse_table {
	my $self = shift;
	my $args = shift;
	my ($file, $feature, $subfeature, $simplify);
	if (ref $args) {
		$file = $args->{file} || '';
		$feature = $args->{feature} || '';
		$subfeature = $args->{subfeature} || '';
		$simplify = (exists $args->{simplify} and defined $args->{simplify}) ? 
			$args->{simplify} : 1; # default is to simplify
	}
	else {
		# no hash reference, assume just a file name
		$file = $args;
		undef $args;
		$feature = undef;
		$subfeature = '';
		$simplify = 1;
	}
	unless ($file) {
		carp "no annotation file provided to parse!";
		return;
	}
	
	# the file format determines the parser class
	my $flavor = $self->taste_file($file) or return;
	my $class = 'Bio::ToolBox::parser::' . $flavor;
	eval {load $class};
	if ($@) {
		carp "unable to load $class! cannot parse $flavor!";
		return;
	}
	
	# open parser
	my $parser = $class->new() or return;
	$parser->open_file($file) or return;
	my $typelist = $parser->typelist;
	
	# set feature based on the type list from the parser
	unless ($feature) {
		if ($typelist =~ /gene/i) {
			$feature = 'gene';
		}
		elsif ($typelist eq 'region') {
			$feature = 'region';
		}
		else {
			$feature = 'rna'; # generic RNA
		}
	}
	
	# set parser parameters
	$parser->simplify($simplify);
	if ($subfeature) {
		$parser->do_exon(1) if $subfeature =~ /exon/i;
		$parser->do_cds(1) if $subfeature =~ /cds/i;
		$parser->do_utr(1) if $subfeature =~ /utr|untranslated/i;
		$parser->do_codon(1) if $subfeature =~/codon/i;
	}
	if ($feature =~ /gene$/i) {
		$parser->do_gene(1);
	}
	else {
		$parser->do_gene(0);
	}
	my $mrna_check = 0;
	if (lc($feature) eq 'mrna' and $parser->typelist !~ /mrna/i and not $self->last_row) {
		# user requested mRNA for a new data file but it's not present in the type list
		# look for it the hard way by parsing CDS too - sigh
		load('Bio::ToolBox::GeneTools', 'is_coding');
		$parser->do_cds(1);
		$mrna_check = 1;
	}
	
	# parse the table
	$parser->parse_file or return;
	
	# store the SeqFeature objects
	if ($self->last_row > 0) {
		# we already have a table, presumably representing the features
		my $count = 0;
		$self->iterate( sub {
			my $row = shift;
			my $f = $parser->find_gene(
				name  => $row->name,
				id    => $row->id,
			);
			if ($f) {
				$self->store_seqfeature($row->row_index, $f);
				$count++;
			}
		} );
		unless ($count == $self->last_row) {
			die <<PARSEFAIL;
Not all features in the input file could be matched to a corresponding SeqFeature 
object in the annotation file $file.
Double check your input and annotation files. You can create a new table by just 
providing your annotation file.
PARSEFAIL
		}
	}
	else {
		# create a new table
		$self->add_column('Primary_ID');
		$self->add_column('Name');
		
		# check for chromosome exclude
		my $chr_exclude;
		if ($args) {
			$chr_exclude = $args->{chrskip} || undef;
		}
		
		# fill table with features
		while (my $f = $parser->next_top_feature) {
			if ($chr_exclude) {
				next if $f->seq_id =~ /$chr_exclude/i;
			}
			my $type = $f->type;
			if ($f->type =~ /$feature/i or ($mrna_check and is_coding($f)) ) {
				my $index = $self->add_row([ $f->primary_id, $f->display_name ]);
				$self->store_seqfeature($index, $f);
			}
		}
		unless ($self->last_row) {
			printf " Zero '%s' features found!\n Check your feature or try generic features like gene, mRNA, or transcript\n",
				$feature;
		}
		$self->database("Parsed:$file");
		$self->feature($feature);
		$self->add_comment("Chromosomes excluded: $chr_exclude") if $chr_exclude;
		
		# add input parsed file metadata
		$self->add_file_metadata($file);
		# but delete some stuff, just want basename
		undef $self->{extension};
		undef $self->{filename};
		undef $self->{path};
	}
	return 1;
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
	my $name_ref = ref $name;
	if ($name_ref eq 'ARRAY') {
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
	elsif ($name_ref eq 'Bio::DB::GFF::Typename') {
		# a Typename object that was selected from a SeqFeature::Store database
		$self->{$column} = {
			'name'      => $name->asString,
			'index'     => $column,
		};
		$self->{data_table}->[0][$column] = $name->asString;
	}
	elsif ($name_ref eq '') {
		# just a name
		$self->{$column} = {
			'name'      => $name,
			'index'     => $column,
		};
		$self->{data_table}->[0][$column] = $name;
	}
	else {
		cluck "unrecognized reference '$name_ref'! pass a scalar value or array reference";
		return;
	}
	
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
	elsif ($_[0] and ref $_[0] eq 'Bio::ToolBox::Data::Feature') {
		@row_data = $_[0]->row_values;
	}
	elsif ($_[0] and $_[0] =~ /\t/) {
		@row_data = split /\t/, $_[0];
	}
	else {
		@row_data = map {'.'} (1 .. $self->{number_columns});
	}
	if (scalar @row_data > $self->{number_columns}) {
		cluck "row added has more elements than table columns! truncating row elements\n"; 
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

sub get_row {
	my ($self, $row_number) = @_;
	return unless defined $self->{data_table}->[$row_number];
	my @options = (
		'data'      => $self,
		'index'     => $row_number,
	);
	if (
		exists $self->{SeqFeatureObjects} and
		defined $self->{SeqFeatureObjects}->[$row_number] 
	) {
		push @options, 'feature', $self->{SeqFeatureObjects}->[$row_number];
	}
	return Bio::ToolBox::Data::Feature->new(@options);	
}

sub delete_row {
	my $self = shift;
	my @deleted = sort {$b <=> $a} @_;
	while (@deleted) {
		my $d = shift @deleted;
		splice( @{ $self->{data_table} }, $d, 1);
		$self->{last_row}--;
		if (exists $self->{SeqFeatureObjects}) {
			splice( @{ $self->{SeqFeatureObjects} }, $d, 1);
		}
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
		confess "iterate_function() method requires a code reference!";
	}
	my $stream = $self->row_stream;
	while (my $row = $stream->next_row) {
		&$code($row);
	}
	return 1;
}


#### Stored SeqFeature manipulation

sub store_seqfeature {
	my ($self, $row_i, $seqfeature) = @_;
	unless (defined $row_i and ref($seqfeature)) {
		confess "must provide a row index and SeqFeature object!";
	}
	confess "invalid row index" unless ($row_i <= $self->last_row);
	$self->{SeqFeatureObjects} ||= [];
	$self->{SeqFeatureObjects}->[$row_i] = $seqfeature;
	return 1;
}

sub delete_seqfeature {
	my ($self, $row_i) = @_;
	confess "invalid row index" unless ($row_i <= $self->last_row);
	return unless $self->{SeqFeatureObjects};
	undef $self->{SeqFeatureObjects}->[$row_i];
}

sub collapse_gene_transcripts {
	my $self = shift;
	unless ($self->feature_type eq 'named') {
		carp "Table does not contain named features!";
		return;
	}
	
	# load module
	my $class = "Bio::ToolBox::GeneTools";
	eval {load $class, qw(collapse_transcripts)};
	if ($@) {
		carp "unable to load $class! cannot collapse transcripts!";
		return;
	}
	
	# collapse the transcripts
	my $success = 0;
	if (exists $self->{SeqFeatureObjects}) {
		# we should have stored SeqFeature objects, probably from parsed table
		for (my $i = 1; $i <= $self->last_row; $i++) {
			my $feature = $self->{SeqFeatureObjects}->[$i] or next;
			my $collSeqFeat = collapse_transcripts($feature) or next;
			$success += $self->store_seqfeature($i, $collSeqFeat);
		}
	}
	else {
		# no stored SeqFeature objects, probably names pointing to a database
		# we will have to fetch the feature from a database
		my $db = $self->open_meta_database(1) or  # force open a new db connection
			confess "No SeqFeature objects stored and no database connection!";
		my $name_i = $self->name_column;
		my $id_i = $self->id_column;
		my $type_i = $self->type_column;
		for (my $i = 1; $i <= $self->last_row; $i++) {
			my $feature = get_db_feature(
				db    => $db,
				name  => $self->value($i, $name_i) || undef,
				type  => $self->value($i, $type_i) || undef,
				id    => $self->value($i, $id_i) || undef,
			) or next;
			my $collSeqFeat = collapse_transcripts($feature) or next;
			$success += $self->store_seqfeature($i, $collSeqFeat);
		}
	}
	
	return $success;
}

sub add_transcript_length {
	my $self = shift;
	my $subfeature = shift || 'exon';
	unless (exists $self->{SeqFeatureObjects}) {
		carp "no SeqFeature objects stored for collapsing!";
		return;
	}
	
	# load module
	eval {load "Bio::ToolBox::GeneTools", qw(get_transcript_length get_transcript_cds_length
		get_transcript_5p_utr_length get_transcript_3p_utr_length)};
	if ($@) {
		carp "unable to load Bio::ToolBox::GeneTools! cannot collapse transcripts!";
		return;
	}
	
	# determine name and calculation method
	my $length_calculator;
	my $length_name;
	if ($subfeature eq 'exon') {
		$length_calculator = \&get_transcript_length;
		$length_name= 'Merged_Transcript_Length';
	}
	elsif ($subfeature eq 'cds') {
		$length_calculator = \&get_transcript_cds_length;
		$length_name= 'Transcript_CDS_Length';
	}
	elsif ($subfeature eq '5p_utr') {
		$length_calculator = \&get_transcript_5p_utr_length;
		$length_name= 'Transcript_5p_UTR_Length';
	}
	elsif ($subfeature eq '3p_utr') {
		$length_calculator = \&get_transcript_3p_utr_length;
		$length_name= 'Transcript_3p_UTR_Length';
	}
	else {
		carp "unrecognized subfeature type '$subfeature'! No length calculated";
		return;
	}
	
	# add new column and calculate
	my $length_i = $self->add_column($length_name);
	if (exists $self->{SeqFeatureObjects}) {
		# we should have stored SeqFeature objects, probably from parsed table
		for (my $i = 1; $i <= $self->last_row; $i++) {
			my $feature = $self->{SeqFeatureObjects}->[$i] or next;
			my $length = &$length_calculator($feature);
			$self->{data_table}->[$i][$length_i] = $length || $feature->length;
		}
	}
	else {
		# no stored SeqFeature objects, probably names pointing to a database
		# we will have to fetch the feature from a database
		my $db = $self->open_meta_database(1) or  # force open a new db connection
			confess "No SeqFeature objects stored and no database connection!";
		my $name_i = $self->name_column;
		my $id_i = $self->id_column;
		my $type_i = $self->type_column;
		for (my $i = 1; $i <= $self->last_row; $i++) {
			my $feature = get_db_feature(
				db    => $db,
				name  => $self->value($i, $name_i) || undef,
				type  => $self->value($i, $type_i) || undef,
				id    => $self->value($i, $id_i) || undef,
			) or next;
			my $length = &$length_calculator($feature);
			$self->{data_table}->[$i][$length_i] = $length || $feature->length;
			$self->store_seqfeature($i, $feature);
				# to store or not to store? May explode memory, but could save from 
				# expensive db calls later
		}
	}
	
	return $length_i;
}



### Data structure manipulation

sub sort_data {
	my $self = shift;
	my $index = shift;
	my $direction = shift || 'i';
	
	# confirm options
	return unless exists $self->{$index}{name};
	unless ($direction =~ /^[id]/i) {
		carp "unrecognized sort order '$direction'! Must be i or d";
		return;
	}
	
	# Sample the dataset values
	# this will be used to guess the sort method, below
	my $example; # an example of the dataset
	foreach (my $i = 1; $i <= $self->last_row; $i++) {
		# we want to avoid null values, so keep trying
		# null being . or any variation of N/A, NaN, inf
		my $v = $self->value($i, $index);
		if (defined $v and $v !~ /^(?:\.|n\/?a|nan|\-?inf)$/i) {
			# a non-null value, take it
			$example = $v;
			last;
		} 
	}
	
	# Determine sort method, either numerical or alphabetical
	my $sortmethod; 
	if ($example =~ /[a-z]/i) { 
		# there are detectable letters
		$sortmethod = 'ascii';
	} 
	elsif ($example =~ /^\-?\d+\.?\d*$/) {
		# there are only digits, allowing for minus sign and a decimal point
		# I don't think this allows for exponents, though
		$sortmethod = 'numeric';
	} 
	else { 
		# unable to determine (probably alphanumeric), sort asciibetical
		$sortmethod = 'ascii';
	}
	
	# Re-order the datasets
	# Directly sorting the @data array is proving difficult. It keeps giving me
	# a segmentation fault. So I'm using a different approach by copying the 
	# @data_table into a temporary hash.
		# put data_table array into a temporary hash
		# the hash key will the be dataset value, 
		# the hash value will be the reference the row data
	my %datahash;
	
	# reorder numerically
	if ($sortmethod eq 'numeric') {
		for my $row (1..$self->last_row) {
			my $value = $self->value($row, $index); 
			
			# check to see whether this value exists or not
			while (exists $datahash{$value}) {
				# add a really small number to bump it up and make it unique
				# this, of course, presumes that none of the dataset values
				# are really this small - this may be an entirely bad 
				# assumption!!!!! I suppose we could somehow calculate an 
				# appropriate value.... nah.
				# don't worry, we're only modifying the value used for sorting,
				# not the actual value
				$value += 0.00000001; 
			}
			# store the row data reference
			$datahash{$value} = $self->{data_table}->[$row]; 
		}
		
		# re-fill the array based on the sort direction
		if ($direction =~ /^i/i) { 
			# increasing sort
			my $i = 1; # keep track of the row, skip the header
			foreach (sort {$a <=> $b} keys %datahash) {
				# put back the reference to the anonymous array of row data
				$self->{data_table}->[$i] = $datahash{$_};
				$i++; # increment for next row
			}
		} 
		else { 
			# decreasing sort
			my $i = 1; # keep track of the row, skip the header
			foreach (sort {$b <=> $a} keys %datahash) {
				# put back the reference to the anonymous array of row data
				$self->{data_table}->[$i] = $datahash{$_};
				$i++; # increment for next row
			}
		}
		
		# summary prompt
		printf " Data table sorted numerically by the contents of %s\n",
			$self->name($index);
	} 
	
	# reorder asciibetically
	elsif ($sortmethod eq 'ascii') {
		for my $row (1..$self->last_row) {
			# get the value to sort by
			my $value = $self->value($row, $index); 
			if (exists $datahash{$value}) { 
				# not unique
				my $n = 1;
				my $lookup = $value . sprintf("03%d", $n);
				# we'll try to make a unique value by appending 
				# a number to the original value
				while (exists $datahash{$lookup}) {
					# keep bumping up the number till it's unique
					$n++;
					$lookup = $value . sprintf("03%d", $n);
				}
				$datahash{$lookup} = $self->{data_table}->[$row];
			} 
			else {
				# unique
				$datahash{$value} = $self->{data_table}->[$row];
			}
		}
		
		# re-fill the array based on the sort direction
		if ($direction eq 'i' or $direction eq 'I') { 
			# increasing
			my $i = 1; # keep track of the row
			foreach (sort {$a cmp $b} keys %datahash) {
				# put back the reference to the anonymous array of row data
				$self->{data_table}->[$i] = $datahash{$_};
				$i++; # increment for next row
			}
		} 
		elsif ($direction eq 'd' or $direction eq 'D') { 
			# decreasing
			my $i = 1; # keep track of the row
			foreach (sort {$b cmp $a} keys %datahash) {
				# put back the reference to the anonymous array of row data
				$self->{data_table}->[$i] = $datahash{$_};
				$i++; # increment for next row
			}
		}
		
		# summary prompt
		printf " Data table sorted asciibetically by the contents of '%s'\n",
			$self->name($index);
	}
	
	return 1;
}

sub gsort_data {
	my $self = shift;
	
	# identify indices
	unless ($self->feature_type eq 'coordinate') {
		carp "no chromosome and start/position columns to sort!\n";
		return;
	}
	my $chromo_i = $self->chromo_column;
	my $start_i  = $self->start_column;
	my $stop_i   = $self->stop_column || $start_i;
	
	# collect the data by chromosome: start and end
	my %chrom2row;
	for my $row (1 .. $self->last_row) {
		my $c = $self->{data_table}->[$row]->[$chromo_i];
		$chrom2row{$c} ||= [];
		push @{ $chrom2row{$c} }, [
			$self->{data_table}->[$row]->[$start_i], 
			$self->{data_table}->[$row]->[$stop_i], 
			$self->{data_table}->[$row]
		];
	}
	
	# get sane chromosome order
	my @chroms = sane_chromo_sort(keys %chrom2row); 
	
	# sort the table
	my @sorted_data;
	push @sorted_data, $self->{data_table}->[0];
	if ($self->gff) {
		# sort by increasing start and decreasing end, i.e. decreasing length
		# this should mostly handle putting genes/transcripts before exons
		foreach my $chr (@chroms) {
			push @sorted_data, 
				map { $_->[2] }
				sort { $a->[0] <=> $b->[0] or $b->[1] <=> $a->[1] }
				@{ $chrom2row{$chr} };
		}
	}
	else {
		# sort by increasing start followed by increasing end, i.e. increasing length
		foreach my $chr (@chroms) {
			push @sorted_data, 
				map { $_->[2] }
				sort { $a->[0] <=> $b->[0] or $a->[1] <=> $b->[1] }
				@{ $chrom2row{$chr} };
		}
	}
	$self->{data_table} = \@sorted_data;
	return 1;
}

sub splice_data {
	my ($self, $part, $total_parts) = @_;
	
	unless ($part and $total_parts) {
		confess "ordinal part and total number of parts not passed\n";
	}
	my $part_length = int($self->last_row / $total_parts);
	
	# check for SeqFeatureObjects array
	if (exists $self->{SeqFeatureObjects}) {
		# it needs to be the same length as the data table, it should be
		while (scalar @{$self->{SeqFeatureObjects}} < scalar @{$self->{data_table}}) {
			push @{$self->{SeqFeatureObjects}}, undef;
		}
	}
	
	# splicing based on which part we do 
	if ($part == 1) {
		# remove all but the first part
		splice( 
			@{$self->{'data_table'}}, 
			$part_length + 1 
		);
		if (exists $self->{SeqFeatureObjects}) {
			splice( 
				@{$self->{SeqFeatureObjects}}, 
				$part_length + 1 
			);
		}
	}
	elsif ($part == $total_parts) {
		# remove all but the last part
		splice( 
			@{$self->{'data_table'}}, 
			1,
			$part_length * ($total_parts - 1) 
		);
		if (exists $self->{SeqFeatureObjects}) {
			splice( 
				@{$self->{SeqFeatureObjects}}, 
				1,
				$part_length * ($total_parts - 1) 
			);
		}
	}
	else {
		# splicing in the middle requires two rounds
		
		# remove the last parts
		splice( 
			@{$self->{'data_table'}}, 
			($part * $part_length) + 1
		);
		
		# remove the first parts
		splice( 
			@{$self->{'data_table'}}, 
			1,
			$part_length * ($part - 1) 
		);
		
		if (exists $self->{SeqFeatureObjects}) {
			splice( 
				@{$self->{SeqFeatureObjects}}, 
				($part * $part_length) + 1
			);
			splice( 
				@{$self->{SeqFeatureObjects}}, 
				1,
				$part_length * ($part - 1) 
			);
		}
	}
	
	# update last row metadata
	$self->{'last_row'} = scalar(@{$self->{'data_table'}}) - 1;
	
	# re-open a new un-cached database connection
	if (exists $self->{db_connection}) {
		delete $self->{db_connection};
	}
	return 1;
}



### File functions

sub reload_children {
	my $self = shift;
	my @files = @_;
	return unless @files;
	
	# prepare Stream
	my $class = "Bio::ToolBox::Data::Stream";
	eval {load $class};
	if ($@) {
		carp "unable to load $class! can't reload children!";
		return;
	}
	
	# open first stream
	my $first = shift @files;
	my $Stream = $class->new(in => $first);
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
		$self->add_row($row);
	}
	$Stream->close_fh;
	
	# load remaining files
	foreach my $file (@files) {
		my $Stream = $class->new(in => $file);
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




### Export summary data file
sub summary_file {
	my $self = shift;
	
	# Collect passed arguments
	my %args = @_; 
	
	# parameters
	# either one or more datasets can be summarized
	my $outfile = $args{'filename'} || undef;
	my @datasets;
	if ($args{dataset} and ref $args{dataset} eq 'ARRAY') {
		@datasets = @{ $args{dataset} };
	}
	elsif ($args{dataset} and ref $args{dataset} eq 'SCALAR') {
		push @datasets, $args{dataset};
	}
	my @startcolumns;
	if ($args{startcolumn} and ref $args{startcolumn} eq 'ARRAY') {
		@startcolumns = @{ $args{startcolumn} };
	}
	elsif ($args{startcolumn} and ref $args{startcolumn} eq 'SCALAR') {
		push @startcolumns, $args{startcolumn};
	}
	my @endcolumns;
	$args{endcolumn} ||= $args{stopcolumn};
	if ($args{endcolumn} and ref $args{endcolumn} eq 'ARRAY') {
		@endcolumns = @{ $args{endcolumn} };
	}
	elsif ($args{endcolumn} and ref $args{endcolumn} eq 'SCALAR') {
		push @endcolumns, $args{endcolumn};
	}
	
	
	# Check required values
	unless ($self->verify) {
		cluck "bad data structure!";
		return;
	}
	unless (defined $outfile) {
		if ($self->basename) {
			# use the opened file's filename if it exists
			# prepend the path if it exists
			# the extension will be added later
			$outfile = $self->path . $self->basename;
		}
		else {
			cluck "no filename passed to write_summary_data!\n";
			return;
		}
	}
	
	# Identify possible dataset columns
	my %possibles;
	my %skip = map {$_ => 1} qw (systematicname name id alias aliases type class 
			geneclass chromosome chromo seq_id seqid start stop end gene strand 
			length primary_id merged_transcript_length transcript_cds_length 
			transcript_5p_utr_length transcript_3p_utr_length);
	# walk through the dataset names
	$possibles{unknown} = [];
	for (my $i = 0; $i < $self->number_columns; $i++) {
		if (not exists $skip{ lc $self->{$i}{'name'} }) {
			my $d = $self->metadata($i, 'dataset') || undef;
			if (defined $d) {
				# we have what appears to be a dataset column
				$possibles{$d} ||= [];
				push @{$possibles{$d}}, $i;
			}
			else {
				# still an unknown possibility
				push @{$possibles{unknown}}, $i;
			}
		}
	}
	
	# check datasets
	unless (@datasets) {
		if (scalar keys %possibles > 1) {
			# we will always have the unknown category, so anything more than one 
			# means we found legitimate dataset columns
			delete $possibles{unknown};
		}
		@datasets = sort {$a cmp $b} keys %possibles;
	}
	
	# check starts
	if (scalar @startcolumns != scalar @datasets) {
		@startcolumns = (); # ignore what we were given?
		foreach my $d (@datasets) {
			# take the first column with this dataset
			push @startcolumns, $possibles{$d}->[0];
		}
	}
	
	# check stops
	if (scalar @endcolumns != scalar @datasets) {
		@endcolumns = (); # ignore what we were given?
		foreach my $d (@datasets) {
			# take the last column with this dataset
			push @endcolumns, $possibles{$d}->[-1];
		}
	}
	
	
	# Prepare Data object to store the summed data
	my $summed_data = $self->new(
		feature => 'averaged_windows', 
		columns => ['Window','Midpoint'],
	);

	# Go through each dataset
	foreach my $d (0 .. $#datasets) {
		
		# Prepare score column name
		my $data_name = simplify_dataset_name($datasets[$d]);
		
		# add column
		my $i = $summed_data->add_column($data_name);
		$summed_data->metadata($i, 'dataset', $datasets[$d]);
		
		# tag for remembering we're working with percentile bins
		my $do_percentile = 0;
		
		# remember the row
		my $row = 1;
		
		# Collect summarized data
		for (
			my $column = $startcolumns[$d];
			$column <= $endcolumns[$d];
			$column++
		) { 
		
			# determine the midpoint position of the window
			# this assumes the column metadata has start and stop
			my $midpoint = int(sum0($self->metadata($column, 'start'), 
				$self->metadata($column, 'stop')) / 2); 
		
			# convert midpoint to fraction of 1000 for plotting if necessary
			if (substr($self->name($column), -1) eq '%') {
				$midpoint *= 10; # midpoint * 0.01 * 1000 bp
				$do_percentile++;
			}
			if ($do_percentile and substr($self->name($column), -2) eq 'bp') {
				# working on the extension after the percentile bins
				$midpoint += 1000;
			}
		
			# collect the values in the column
			my @values;
			for my $row (1..$self->last_row) {
				my $v = $self->value($row, $column);
				push @values, $v eq '.' ? 0 : $v;  # treat nulls as zero
			}
		
			# adjust if log value
			my $log = $self->metadata($column, 'log2') || 0;
			if ($log) {
				@values = map { 2 ** $_ } @values;
			}
		
			# determine mean value
			my $window_mean = sum0(@values) / scalar(@values);
			if ($log) { 
				$window_mean = log($window_mean) / log(2);
			}
		
			# push to summed output
			if ($d == 0) {
				# this is the first dataset, so we need to add a row
				$summed_data->add_row( [ $self->{$column}{'name'}, $midpoint, $window_mean ] );
			}
			else {
				# we're summarizing multiple datasets, we already have name midpoint
				# first do sanity check
				if ($summed_data->value($row, 1) != $midpoint) {
					carp("unable to summarize multiple datasets with nonequal columns of data!");
					return;
				}
				$summed_data->value($row, $i, $window_mean);
			}
			$row++;
		}
	}
	
	
	
	# Write summed data
	$outfile =~ s/\.txt(\.gz)?$//i; # strip any .txt or .gz extensions if present
	my $written_file = $summed_data->write_file(
		'filename'  => $outfile . '_summary.txt',
		'gz'        => 0,
	);
	return $written_file;
}





####################################################

package Bio::ToolBox::Data::Iterator;
use Bio::ToolBox::Data::Feature;
use Carp;

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
	my @options = (
		'data'      => $self->{data},
		'index'     => $i,
	);
	if (
		exists $self->{data}->{SeqFeatureObjects} and
		defined $self->{data}->{SeqFeatureObjects}->[$i] 
	) {
		push @options, 'feature', $self->{data}->{SeqFeatureObjects}->[$i];
	}
	return Bio::ToolBox::Data::Feature->new(@options);	
}

sub row_index {
	my $self = shift;
	carp "row_index is a read only method" if @_;
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
