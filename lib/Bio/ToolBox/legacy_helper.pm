package Bio::ToolBox::legacy_helper;
our $VERSION = '1.30';

=head1 NAME

Bio::ToolBox::legacy_helper - exported methods to support legacy API

=head1 DESCRIPTION

These are legacy methods that used to be provided by Bio::ToolBox::data_helper 
and Bio::ToolBox::file_helper, but have now been superseded by the object 
oriented API of L<Bio::ToolBox::Data>. All new scripts should use the 
L<Bio::ToolBox::Data> API and NOT these methods. This module will go away 
in the near future. (Actually, it will just be moved to the <Bio::ToolBox::Extra> 
distribution, where most legacy programs live out any remaining usefulness.)

=cut

use strict;
require Exporter;
use Carp qw(carp cluck croak confess);
use Bio::ToolBox::Data;
our @ISA = qw(Exporter);
our @EXPORT = qw(
	generate_data_structure
	verify_data_structure
	find_column_index
	open_data_file 
	load_data_file
	write_data_file 
	open_to_read_fh
	open_to_write_fh
	write_summary_data
	check_file
);
our $CLASS = 'Bio::ToolBox::Data';

1;

sub generate_data_structure {
	# Collect the feature
	my $feature = shift;
	
	# Collect the array of dataset headers
	my @datasets = @_;
	
	# Initialize the hash structure
	my $Data = $CLASS->new(
		feature     => $feature,
		datasets    => \@datasets,
	);
	
	# the old libraries used a data structure remarkably similar
	# to the Data object, so even though we are returning a blessed object, 
	# the old programs will access the underlying data structure just as before and 
	# should work in the same manner, more or less
	return $Data;
}

sub verify_data_structure {
	my $Data = shift;
	
	unless (ref($Data) eq $CLASS) {
		# try an impromptu blessing and hope this works!
		bless($Data, $CLASS);
	}
	return $Data->verify;
}

sub find_column_index {
	my ($Data, $name) = @_;
	
	unless (ref($Data) eq $CLASS) {
		# try an impromptu blessing and hope this works!
		bless($Data, $CLASS);
	}
	return $Data->find_column($name);
}

sub open_data_file {
	my $file = shift;
	my $Data = $CLASS->new();
	my $filename = $Data->check_file($file);
	$Data->add_file_metadata($filename);
	my $fh = $Data->open_to_read_fh or return;
	$Data->parse_headers($fh);
	return ($Data->{fh}, $Data);
}

sub load_data_file {
	my $file = shift;
	my $Data = $CLASS->new(file => $file);
	return $Data;
}

sub write_data_file {
	my %args = @_; 
	$args{'data'}     ||= undef;
	$args{'filename'} ||= undef;
	$args{'format'}   ||= undef;
	unless (exists $args{'gz'}) {$args{'gz'} = undef} 
	my $Data = $args{data};
	return unless $Data;
	unless (ref($Data) eq $CLASS) {
		# try an impromptu blessing and hope this works!
		bless($Data, $CLASS);
	}
	return $Data->write_file(
		filename    => $args{'filename'},
		'format'    => $args{'format'},
		'gz'        => $args{gz},
	);
}

sub open_to_read_fh {
	my $filename = shift;
	return unless defined $filename;
	return $CLASS->open_to_read_fh($filename);
}

sub open_to_write_fh {
	my ($filename, $gz, $append) = @_;
	return unless defined $filename;
	return $CLASS->open_to_read_fh($filename, $gz, $append);
}

sub write_summary_data {
	my %args = @_; 
	my $Data = $args{data} || undef;
	return unless defined $Data;
	unless (ref($Data) eq $CLASS) {
		# try an impromptu blessing and hope this works!
		bless($Data, $CLASS);
	}
	delete $args{data};
	return $Data->summary_file(%args);
}

sub check_file {
	my $filename = shift;
	return $CLASS->check_file($filename);
}

__END__

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

