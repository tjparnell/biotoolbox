package Bio::ToolBox::utility;
our $VERSION = '1.69';

=head1 NAME

Bio::ToolBox::utility - common utility functions for Bio::ToolBox

=head1 DESCRIPTION

These are general subroutines that don't fit in with the other modules.

=head1 REGULAR SUBROUTINES

The following subroutines are automatically exported when you use this module.

=over 4

=item parse_list

	my $index_request = '1,2,5-7';
	my @indices = parse_list($index_request); # returns [1,2,5,6,7]

This subroutine parses a scalar value into a list of values. The scalar is 
a text string of numbers (usually column or dataset indices) delimited by 
commas and/or including a range. For example, a string "1,2,5-7" would become 
an array of [1,2,5,6,7].

Pass the module the scalar string.

It will return the array of numbers.

=item format_with_commas

	my $count = '4327908475';
	printf " The final count was %s\n", format_with_commas($count);

This subroutine process a large number (e.g. 4327908475) into a human-friendly 
version with commas delimiting the thousands (4,327,908,475).

Pass the module a scalar string with a number value.

It will return a scalar value containing the formatted number.

=item ask_user_for_index

	my @answers = ask_user_for_index($Data, 'Please enter 2 or more columns   ');

This subroutine will present the list of column names from a L<Bio::ToolBox::Data> 
structure along with their numeric indexes to the user and prompt for one 
or more to be selected and entered. The function is smart enough to only print 
the list once (if it hasn't changed) so as not to annoy the user with repeated 
lists of header names when used more than once. A text prompt should be provided, 
or a generic one is used. The list of indices are validated, and a warning printed 
for invalid responses. The responses are then returned as a single value or array, 
depending on context.

=item simplify_dataset_name

	my $simple_name = simplify_dataset_name($dataset);

This subroutine will take a dataset name and simplify it. Dataset names may 
often be file names of data files, such as Bam and bigWig files. These may 
include a C<file:>, C<http:>, or C<ftp:> prefix, one or more directory paths, 
and one or more file name extensions. Additionally, more than one dataset 
may be combined, for example two stranded bigWig files, with an ampersand. 
This function will safely remove the prefix, directories, and everything after 
the first period. 





=back

=cut

use strict;
use Carp;
use File::Spec;
require Exporter;
use Bio::ToolBox::Data::file;


### Variables
# Export
our @ISA = qw(Exporter);
our @EXPORT = qw(
	parse_list
	format_with_commas
	ask_user_for_index
	simplify_dataset_name
);
our $DATA_COLNAMES  = undef;
our $DATA_FILENAME  = undef;

### The True Statement
1; 



#################   The Subroutines   ###################

### Parse string into list
sub parse_list {
	# this subroutine will parse a string into an array
	# it is designed for a string of numbers delimited by commas
	# a range of numbers may be specified using a dash
	# hence 1,2,5-7 would become an array of 1,2,5,6,7
	
	my $string = shift;
	return unless defined $string;
	if ($string =~ /[^\d,\-\s\&]/) {
		carp " the string contains characters that can't be parsed\n";
		return;
	}
	my @list;
	foreach (split /[,\s+]/, $string) {
		# check for a range
		if (/\-/) { 
			my ($start, $stop) = split /\-/;
			# add each item in the range to the list
			for (my $i = $start; $i <= $stop; $i++) {
				push @list, $i;
			}
			next;
		} 
		else {
			# either an ordinary number or an "&"ed list of numbers 
			push @list, $_;
		}
	}
	return @list;
}



### Format a number into readable comma-delimited by thousands number
sub format_with_commas {
	# for formatting a number with commas
	my $number = shift;
	if ($number =~ /[^\d,\-\.]/) {
		carp " the string contains characters that can't be parsed\n";
		return $number;
	}
	
	# check for decimals
	my ($integers, $decimals);
	if ($number =~ /^\-?(\d+)\.(\d+)$/) {
		$integers = $1;
		$decimals = $2;
	}
	else {
		$integers = $number;
	}
	
	# format
	my @digits = split //, $integers;
	my @formatted;
	while (@digits) {
		if (@digits > 3) {
			unshift @formatted, pop @digits;
			unshift @formatted, pop @digits;
			unshift @formatted, pop @digits;
			unshift @formatted, ',';
		}
		else {
			while (@digits) {
				unshift @formatted, pop @digits;
			}
		}
	}
	
	# finished
	my $final = join("", @formatted);
	$final .= $decimals if defined $decimals;
	return $final;
}


sub ask_user_for_index {
	my $Data = shift;
	my $line = shift || ' Enter the desired column index   ';
	unless (ref($Data) =~ /Bio::ToolBox::Data/) {
		carp "Must pass a Bio::ToolBox::Data object!\n";
		return;
	}
	
	# print column header names only if we have not done so before
	unless (
		# we use filename and column number as indicators 
		$Data->filename eq $DATA_FILENAME and 
		join(";", $Data->list_columns) eq $DATA_COLNAMES
	) {
		print " These are the columns in the file\n";
		my $i = 0;
		foreach ($Data->list_columns) {
			print "  $i\t$_\n";
			$i++;
		}
		# remember for next time
		$DATA_FILENAME = $Data->filename;
		$DATA_COLNAMES = join(";", $Data->list_columns);
	}
	print $line;
	
	# get response
	my $response = <STDIN>;
	chomp $response;
	my @indices = parse_list($response); 
	
	# verify
	my @good;
	foreach (@indices) {
		if ($Data->name($_)) {
			push @good, $_;
		}
		else {
			print "  $_ is not a valid index!\n";
		}
	}
	return wantarray ? @good : $good[0];
}


sub simplify_dataset_name {
	my $dataset = shift;
	my $new_name;
	
	# strip any file prefix
	$dataset =~ s/^(?:file|http|ftp):\/*//;
	
	if ($dataset =~ /&/) {
		# a combination dataset
		foreach (split /&/, $dataset) {
			my $n = simplify_dataset_name($_);
			if ($new_name) {
				$new_name .= '&' . $n;
			}
			else {
				$new_name = $n;
			}
		}
	}
	elsif ($dataset =~ m|/|) {
		# appears to have paths
		my (undef, undef, $file_name) = File::Spec->splitpath($dataset);
		$file_name =~ s/^([\w\d\-\_]+)\..+$/$1/i; # take everything up to first .
		$new_name = $file_name;
	}
	else {
		# strip everything after first period, like above
		$dataset =~ s/^([\w\d\-\_]+)\..+$/$1/i; # take everything up to first .
		$new_name = $dataset;
	}
	return $new_name;
}


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

