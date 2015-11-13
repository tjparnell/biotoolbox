package Bio::ToolBox::Gene::DB;
our $VERSION = '1.35';

=head1 NAME

Bio::ToolBox::Gene::DB - Database for gene SeqFeatures

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=over 4

=back

=cut

use strict;
use Carp qw(carp cluck croak confess);
use Archive::Zip qw( :ERROR_CODES );
use Storable qw(nfreeze thaw);

### The True Statement
1; 

sub new {
	my $class = shift;
	my %args = @_;
	my %self = {};
	
	
	
	return bless \%self, $class;
}


__END__

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
