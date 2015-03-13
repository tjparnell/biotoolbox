package Bio::ToolBox::Data::meta;
our $VERSION = 1.20;

=head1 NAME

Bio::ToolBox::Data::meta - metadata for Data tables

=head1 DESCRIPTION

General methods for describing the metadata in a Bio::ToolBox::Data 
data table. This module should not be used directly. See 
Bio::ToolBox::Data for more information.

=cut

use Bio::ToolBox::data_helper qw(find_column_index);
use Bio::ToolBox::file_helper qw(parse_filename);
use strict;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(feature feature_type program database gff bed number_columns 
	last_row filename basename path extension comments 
	list_columns name metadata find_column _find_column_indices chromo_column  
	start_column stop_column strand_column name_column type_column id_column);




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

sub comments {
	my $self = shift;
	my @comments = map {s/[\r\n]+//g} @{ $self->{other} };
	# comments are not chomped when loading
	# side effect of dealing with rare commented header lines with null values at end
	return @comments;
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

