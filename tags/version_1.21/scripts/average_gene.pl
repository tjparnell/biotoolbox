#!/usr/bin/perl

# documentation at end of file

use strict;
use FindBin qw($Bin);

print <<MESSAGE;
####################        Warning         ####################
The original 'average_gene.pl' script has been renamed 
'get_binned_data.pl' to reflect a more accurate script function.
Please use 'get_binned_data.pl' in the future. No functionality 
has changed.

Executing get_binned_data.pl for you in 5 seconds.
################################################################
MESSAGE

sleep 5;

exec ("$Bin/get_binned_data.pl", @ARGV) or 
	die " Could not execute $Bin/get_binned_data.pl!";

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
