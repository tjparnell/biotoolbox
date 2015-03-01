package Bio::ToolBox::Data::common;
our $VERSION = 1.25;

=head1 NAME

Bio::ToolBox::Data::common - common methods for Data objects

=head1 DESCRIPTION

Common methods for metadata and manipulation in a L<Bio::ToolBox::Data> 
data table and L<Bio::ToolBox::Data::Stream> file stream. This module 
should not be used directly. See the respective modules for more information.

=cut

use strict;
require Exporter;
use Carp qw(carp cluck croak confess);
use Bio::ToolBox::data_helper qw(find_column_index);
use Bio::ToolBox::file_helper qw(parse_filename);
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
);
our @ISA = qw(Exporter);
our @EXPORT = qw(open_database verify_dataset delete_column reorder_column
	feature feature_type program database gff bed number_columns 
	last_row filename basename path extension comments add_comment delete_comment
	list_columns name metadata copy_metadata delete_metadata 
	find_column _find_column_indices chromo_column start_column stop_column 
	strand_column name_column type_column id_column);


#### Database methods #####

sub open_database {
	my $self = shift;
	my $force = shift || 0;
	return unless $self->{db};
	if (exists $self->{db_connection}) {
		return $self->{db_connection} unless $force;
	}
	my $db = open_db_connection($self->{db}, $force);
	return unless $db;
	$self->{db_connection} = $db;
	return $db;
}

sub verify_dataset {
	my $self = shift;
	my $dataset = shift;
	my $database = shift; # name or object?
	return unless $dataset;
	if (exists $self->{verfied_dataset}{$dataset}) {
		return $self->{verfied_dataset}{$dataset};
	}
	else {
		if ($dataset =~ /^(?:file|http|ftp)/) {
			# local or remote file already verified?
			$self->{verfied_dataset}{$dataset} = $dataset;
			return $dataset;
		}
		$database ||= $self->open_database;
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



#### Column Manipulation ####

sub delete_column {
	my $self = shift;
	if (defined $self->{fh}) {
		# Stream file handle is opened
		cluck "Cannot modify columns when a Stream file handle is opened!";
		return;
	}
	
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
	if (defined $self->{fh}) {
		# Stream file handle is opened
		cluck "Cannot modify columns when a Stream file handle is opened!";
		return;
	}
	
	# reorder data table
	my @order = @_;
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
	delete $self->{column_indices} if exists $self->{column_indices};
	return 1;
}



#### General Metadata ####

sub feature {
	my $self = shift;
	if (@_) {
		$self->{feature} = shift;
	}
	return $self->{feature};
}

sub feature_type {
	my $self = shift;
	if (exists $self->{feature_type}) {
		return $self->{feature_type};
	}
	my $feature_type;
	if (defined $self->chromo_column and defined $self->start_column) {
		$feature_type = 'coordinate';
	}
	elsif (defined $self->id_column or 
		( defined $self->type_column and defined $self->name_column )
	) {
		$feature_type = 'named';
	}
	else {
		$feature_type = 'unknown';
	}
	$self->{feature_type} = $feature_type;
	return $feature_type;
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

sub gff {
	my $self = shift;
	if ($_[0] and $_[0] =~ /^[123]/) {
		$self->{gff} = $_[0];
	}
	return $self->{gff};
}

sub bed {
	my $self = shift;
	if ($_[0] and $_[0] =~ /^\d+/) {
		$self->{bed} = $_[0];
	}
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
			# parse_filename is exported from Bio::ToolBox::file_helper
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



#### General Comments ####

sub comments {
	my $self = shift;
	my @comments = @{ $self->{other} };
	foreach (@comments) {s/[\r\n]+//g}
	# comments are not chomped when loading
	# side effect of dealing with rare commented header lines with null values at end
	return @comments;
}

sub add_comment {
	my $self = shift;
	my $comment = shift or return;
	# comment is not required to be prefixed with "# ", it will be added when saving
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



#### Column Metadata ####

sub list_columns {
	my $self = shift;
	my @list;
	for (my $i = 0; $i < $self->number_columns; $i++) {
		push @list, $self->{$i}{'name'};
	}
	return wantarray ? @list : \@list;
}

sub name {
	my $self = shift;
	my ($index, $new_name) = @_;
	return unless defined $index;
	return unless exists $self->{$index}{name};
	if (defined $new_name) {
		$self->{$index}{name} = $new_name;
		if (exists $self->{data_table}) {
			$self->{data_table}->[0][$index] = $new_name;
		}
		elsif (exists $self->{column_names}) {
			$self->{column_names}->[$index] = $new_name;
		}
	}
	return $self->{$index}{name};
}

sub metadata {
	my $self = shift;
	my ($index, $key, $value) = @_;
	return unless defined $index;
	return unless exists $self->{$index};
	if ($key and $key eq 'name') {
		return $self->name($index, $value);
	}
	if ($key and defined $value) { 
		# we are setting a new value
		$self->{$index}{$key} = $value;
		return $value;
	}
	elsif ($key and not defined $value) {
		if (exists $self->{$index}{$key}) {
			# retrieve a value
			return $self->{$index}{$key};
		}
		else {
			# set a new empty key
			$self->{$index}{$key} = q();
			return 1;
		}
	}
	else {
		my %hash = %{ $self->{$index} };
		return wantarray ? %hash : \%hash;
	}
}

sub delete_metadata {
	my $self = shift;
	my ($index, $key) = @_;
	return unless defined $index;
	if (defined $key) {
		if (exists $self->{$index}{$key}) {
			return delete $self->{$index}{$key};
		}
	}
	else {
		# user wants to delete the metadata
		# but we need to keep the basics name and index
		foreach my $key (keys %{ $self->{$index} }) {
			next if $key eq 'name';
			next if $key eq 'index';
			delete $self->{$index}{$key};
		}
	}
}

sub copy_metadata {
	my ($self, $source, $target) = @_;
	return unless (exists $self->{$source}{name} and exists $self->{$target}{name});
	my $md = $self->metadata($source);
	delete $md->{name};
	delete $md->{'index'};
	delete $md->{'AUTO'} if exists $md->{'AUTO'}; # presume this is no longer auto index
	foreach (keys %$md) {
		$self->{$target}{$_} = $md->{$_};
	}
	return 1;
}



#### Column Indices ####

sub find_column {
	my ($self, $name) = @_;
	return unless $name;
	return find_column_index($self, $name);
}

sub _find_column_indices {
	my $self = shift;
	# these are hard coded index name regex to accomodate different possibilities
	# these do not include parentheses for grouping
	# non-capturing parentheses will be added later in the sub for proper 
	# anchoring and grouping - long story why, don't ask
	my $name   = find_column_index($self, '^name|geneName|transcriptName|geneid|id|alias');
	my $type   = find_column_index($self, '^type|class|primary_tag');
	my $id     = find_column_index($self, '^primary_id');
	my $chromo = find_column_index($self, '^chr|seq|ref|ref.?seq');
	my $start  = find_column_index($self, '^start|position|pos|txStart$');
	my $stop   = find_column_index($self, '^stop|end|txEnd');
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

sub end_column {
	return shift->stop_column;
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

