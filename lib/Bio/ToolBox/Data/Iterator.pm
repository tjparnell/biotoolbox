package Bio::ToolBox::Data::Iterator;

use warnings;
use strict;
use Carp;
use Bio::ToolBox::Data::Feature;

our $VERSION = 1.70;

sub new {
	my ( $class, $data ) = @_;
	my %iterator = (
		'index' => 1,
		'data'  => $data,
	);
	return bless \%iterator, $class;
}

sub next_row {
	my $self = shift;
	return if $self->{'index'} > $self->{data}->{last_row};    # no more
	my $i = $self->{'index'};
	$self->{'index'}++;
	my @options = (
		'data'  => $self->{data},
		'index' => $i,
	);
	if ( exists $self->{data}->{SeqFeatureObjects}
		and defined $self->{data}->{SeqFeatureObjects}->[$i] )
	{
		push @options, 'feature', $self->{data}->{SeqFeatureObjects}->[$i];
	}
	return Bio::ToolBox::Data::Feature->new(@options);
}

sub row_index {
	my $self = shift;
	carp 'row_index is a read only method' if @_;
	return $self->{'index'};
}

1;

=head1 NAME

Bio::ToolBox::Data::Iterator - Class for iterating through Data tables

=head1 SYNOPSIS

  # open a Bio::ToolBox::Data object
  # initiate an iterator stream
  my $stream = $Data->row_stream;
  while (my $row = $stream->next_row) {
     # each $row is a Bio::ToolBox::Data::Feature object
     # representing the row in the data table
     my $value = $row->value($index);
     # do something with $value
  }

=head1 DESCRIPTION

This is an iteration object for iterating through the rows of a
L<Bio::ToolBox::Data> object table.

This should not be created directly by end-users. Rather, see
L<Bio::ToolBox::Data\row_stream> for details.

=head1 METHODS

There is essentially only one method for end-users, C<next_row>.

=over 4

=item new

Pass a valid  L<Bio::ToolBox::Data> object.

=item next_row

This returns the next row in a Data table as a
L<Bio::ToolBox::Data::Feature> row object. If SeqFeature objects are
associated with the row, perhaps from a parsed input annotation file,
then they are automatically associated with the row object. (Previous
versions required separately calling the seqfeature() row method to
perform this.)

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.


