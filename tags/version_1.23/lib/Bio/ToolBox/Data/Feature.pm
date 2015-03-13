package Bio::ToolBox::Data::Feature;
our $VERSION = 1.23;

=head1 NAME

Bio::ToolBox::Data::Feature - Objects representing rows in a data table

=head1 DESCRIPTION

A Bio::ToolBox::Data::Feature is an object representing a row in the 
data table. Usually, this in turn represents an annotated feature or 
segment in the genome. As such, this object provides convenient 
methods for accessing and manipulating the values in a row, as well as 
methods for working with the represented genomic feature.

This class should not be used directly by the user. Rather, Feature 
objects are generated from a Bio::ToolBox::Data::Iterator object 
(generated itself from the row_stream() function in Bio::ToolBox::Data), 
or the iterate() function in Bio::ToolBox::Data. Please see the 
documentation for Bio::ToolBox::Data for more information.

Example of working with a stream object.
	
	  my $Data = Bio::ToolBox::Data->new(file => $file);
	  
	  # stream method
	  my $stream = $Data->row_stream;
	  while (my $row = $stream->next_row) {
		 # each $row is a Bio::ToolBox::Data::Feature object
		 # representing the row in the data table
		 my $value = $row->value($index);
		 # do something with $value
	  }
	  
	  # iterate method
	  $Data->iterate( sub {
	     my $row = shift;
	     my $number = $row->value($index);
	     my $log_number = log($number);
	     $row->value($index, $log_number);
	  } );


=head1 METHODS

=head2 General information methods

=over 4

=item row_index

Returns the index position of the current data row within the 
data table. Useful for knowing where you are at within the data 
table.

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

=back

=head2 Methods to access row feature attributes

These methods return the corresponding value, if present in the 
data table, based on the column header name. If the row represents 
a named database object, try calling the feature() method first. 
This will retrieve the database SeqFeature object, and the attributes 
can then be retrieved using the methods below or on the actual 
database SeqFeature object.

These methods do not set attribute values. If you need to change the 
values in a table, use the value() method below.

=over 4

=item seq_id

The name of the chromosome the feature is on.

=item start

=item end

=item stop

The coordinates of the feature or segment. Coordinates from known 
0-based file formats, e.g. BED, are returned as 1-based. Coordinates 
must be integers to be returned. Zero or negative start coordinates 
are assumed to be accidents or poor programming and transformed to 1. 
Zero or negative stop coordinates are just ignored. Use the value() 
method if you don't want this to happen.

=item strand

The strand of the feature or segment. Returns -1, 0, or 1. Default is 0, 
or unstranded.

=item name

=item display_name

The name of the feature.

=item type

The type of feature. Typically either primary_tag or primary_tag:source_tag. 
In a GFF3 file, this represents columns 3 and 2, respectively. In annotation 
databases such as Bio::DB::SeqFeature::Store, the type is used to restrict 
to one of many different types of features, e.g. gene, mRNA, or exon.

=item id

Here, this represents the primary_ID in the database. Note that this 
number is unique to a specific database, and not portable between databases.

=item length

The length of the feature or segment.

=back

=head2 Accessing and setting values in the row.

=over 4

=item value($index)

=item value($index, $new_value)

Returns or sets the value at a specific column index in the 
current data row.

=item row_values

Returns an array or array reference representing all the values 
in the current data row. 

=back

=head2 Convenience Methods to database functions

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
the $Data-E<gt>database() method. If the feature name or type is not 
present in the table, then nothing is returned.

See Bio::DB::SeqFeature and Bio::SeqFeatureI for more information 
about working with these objects.

=item segment

Returns a database Segment object corresponding to the coordinates 
defined in the Data table row. If a named feature and type are 
present instead of coordinates, then the feature is first automatically 
retrieved and a Segment returned based on its coordinates. The 
database named in the general metadata is used to establish the 
Segment object. If a different database is desired, it should be 
changed first using the general database() method. 

See Bio::DB::SeqFeature::Segment and Bio::RangeI for more information 
about working with Segment objects.

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

If you wish to override coordinates that are present in the 
Data table, for example to extend or shift the given coordinates 
by some amount, then simply pass the new start and end 
coordinates as options to this method.

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

If you wish to override coordinates that are present in the 
Data table, for example to extend or shift the given coordinates 
by some amount, then simply pass the new start and end 
coordinates as options to this method.

Here is an example for collecting positioned scores around 
the 5 prime end of a feature from a BigWigSet directory.
  
  my $stream = $Data->row_stream;
  while (my $row = $stream->next_row) {
     my %position2score = $row->get_position_scores(
        'ddb'       => '/path/to/BigWigSet/',
        'dataset'   => 'MyData',
        'position'  => 5,
     )
     # do something with %position2score
  }

=back

=cut

use strict;
use Carp;
use Bio::ToolBox::db_helper qw(
	get_feature
	get_chromo_region_score
	get_region_dataset_hash
);




### Initialization

# this should only be called from Bio::ToolBox::Data::Iterator
sub new {
	my $class = shift;
	my %self = @_;
	return unless exists $self{data}; # this must be a Bio::ToolBox::Data object
	return bless \%self, $class;
}


### Set and retrieve values

sub feature_type {
	my $self = shift;
	return $self->{data}->feature_type;
}

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
	if (defined $i) {
		my $v = $self->value($i);
		return (defined $v and $v ne '.') ? $v : undef;
	}
	return $self->{feature}->seq_id if exists $self->{feature};
	return undef;
}

sub start {
	my $self = shift;
	my $i = $self->{data}->start_column;
	my $v = $self->value($i) if defined $i;
	$v = $self->{feature}->start if (not defined $v and exists $self->{feature});
	return undef unless ($v =~ /^\-?\d+$/);
	return $v < 1 ? 1 : $v;
}

sub end {
	my $self = shift;
	my $i = $self->{data}->stop_column;
	my $v = $self->value($i) if defined $i;
	$v = $self->{feature}->end if (not defined $v and exists $self->{feature});
	return undef unless ($v =~ /^\-?\d+$/);
	return $v < 1 ? undef : $v;
}

sub stop {
	return shift->end;
}

sub strand {
	my $self = shift;
	my $i = $self->{data}->strand_column;
	return $self->value($i) if defined $i;
	return $self->{feature}->strand if exists $self->{feature};
	return 0;
}

sub name {
	my $self = shift;
	my $i = $self->{data}->name_column;
	return $self->value($i) if defined $i;
	return $self->{feature}->display_name if exists $self->{feature};
	return undef;
}

sub display_name {
	return shift->name;
}

sub type {
	my $self = shift;
	my $i = $self->{data}->type_column;
	return $self->value($i) if defined $i;
	return $self->{feature}->primary_tag if exists $self->{feature};
	return undef;
}

sub id {
	my $self = shift;
	my $i = $self->{data}->id_column;
	return $self->value($i) if defined $i;
	return $self->{feature}->primary_id if exists $self->{feature};
	return undef;
}

sub length {
	my $self = shift;
	my $s = $self->start;
	my $e = $self->end;
	if (defined $s and defined $e) {
		return $e - $s + 1;
	}
	else {
		return undef;
	}
}

### Data collection convenience methods

sub feature {
	my $self = shift;
	return $self->{feature} if exists $self->{feature};
	return undef unless $self->{data}->database;
	
	# retrieve the feature from the database
	my $id   = $self->id;
	my $name = $self->name;
	my $type = $self->type || $self->{data}->feature;
	return undef unless ($id or ($name and $type));
	my $f = get_feature(
		'db'    => $self->{data}->open_database,
		'id'    => $id,
		'name'  => $name, # we can handle "name; alias" lists later
		'type'  => $type,
	);
	$self->{feature} = $f if $f;
	return $f;
}

sub segment {
	my $self   = shift;
	return unless $self->{data}->database;
	if ($self->feature_type eq 'coordinate') {
		my $chromo = $self->seq_id;
		my $start  = $self->start;
		my $stop   = $self->end || $start;
		my $db = $self->{data}->open_database;
		return $db ? $db->segment($chromo, $start, $stop) : undef;
	}
	elsif ($self->feature_type eq 'named') {
		my $f = $self->feature;
		return $f ? $f->segment : undef;
	}
	else {
		return undef;
	}
}

sub get_score {
	my $self = shift;
	my %args = @_;
	
	# verify coordinates based on type of feature
	if ($self->feature_type eq 'coordinate') {
		# coordinates are already in the table, use those
		$args{chromo} ||= $self->seq_id;
		$args{start}  ||= $self->start;
		$args{stop}   ||= $self->end;
	}
	elsif ($self->feature_type eq 'named') {
		# must retrieve feature from the database first
		my $f = $self->feature;
		return unless $f;
		$args{chromo} ||= $f->seq_id;
		$args{start}  ||= $f->start;
		$args{stop}   ||= $f->end;
	}
	else {
		croak "data table does not have identifiable coordinate or feature identification columns for score collection";
	}
	unless (exists $args{strand} and defined $args{strand}) {
		$args{strand} = $self->strand; 
	}
# 	unless ($args{chromo} and defined $args{start}) {
# 		return;
# 	}
	
	# verify the dataset for the user, cannot trust whether it has been done or not
	my $db = $args{ddb} || $args{db} || $self->{data}->open_database || undef;
	$args{dataset} = $self->{data}->verify_dataset($args{dataset}, $db);
	unless ($args{dataset}) {
		croak "provided dataset was unrecognized format or otherwise could not be verified!";
	}
	
	$args{db} ||= $self->{data}->open_database;
	
	return get_chromo_region_score(%args);
}

sub get_position_scores {
	my $self = shift;
	my %args = @_;
	
	# the get_region_dataset_hash() method can handle both db features and coordinates
	# therefore, this method will also handle both
	# set arguments based on the feature type.
	if ($self->feature_type eq 'named') {
		# confirm that we have appropriate 
		$args{id}     ||= $self->id;
		$args{name}   ||= $self->name;
		$args{type}   ||= $self->type;
		# do NOT assign strand here, it will be assigned in db_helper
	}
	elsif ($self->feature_type eq 'coordinate') {
		$args{chromo} ||= $self->seq_id;
		$args{start}  ||= $self->start;
		$args{stop}   ||= $self->end;
		unless (exists $args{strand} and defined $args{strand}) {
			$args{strand} = $self->strand;
		}
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
