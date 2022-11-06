package Bio::ToolBox::utility;

use warnings;
use strict;
use Carp;
use File::Spec;
use IO::Prompt::Tiny qw(prompt);
require Exporter;

our $VERSION = '1.70';

### Variables
# Export
our @ISA    = qw(Exporter);

our @EXPORT_OK = qw(
	parse_list
	format_with_commas
	ask_user_for_index
	simplify_dataset_name
	sane_chromo_sort
);

my $DATA_COLNAMES = undef;
my $DATA_FILENAME = undef;

### Parse string into list
sub parse_list {

	# this subroutine will parse a string into an array
	# it is designed for a string of numbers delimited by commas
	# a range of numbers may be specified using a dash
	# hence 1,2,5-7 would become an array of 1,2,5,6,7

	my $string = shift;
	return unless defined $string;
	if ( $string =~ /[^\d,\-\s\&]/x ) {
		carp ' the string contains characters that cannot be parsed';
		return;
	}
	my @list;
	foreach ( split /[,\s+]/, $string ) {

		# check for a range
		if (/\-/) {
			my ( $start, $stop ) = split /\-/;

			# add each item in the range to the list
			for ( my $i = $start; $i <= $stop; $i++ ) {
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

	# check number
	my ( $integers, $decimals, $sign );
	if ( $number =~ m/^ (\-)? (\d+) \. (\d+) $/x ) {
		$sign     = $1;
		$integers = $2;
		$decimals = $3;
	}
	elsif ( $number =~ m/^ (\-)? (\d+) $/x ) {
		$sign     = $1;
		$integers = $2;
	}
	else {
		carp ' the string contains characters that cannot be parsed';
		return $number;
	}

	# format
	my @digits = split //, $integers;
	my @formatted;
	while (@digits) {
		if ( @digits > 3 ) {
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
	my $final = $sign ? $sign : q();
	$final .= join( q(), @formatted );
	$final .= '.' . $decimals if defined $decimals;
	return $final;
}

sub ask_user_for_index {
	my $Data = shift;
	my $line = shift || ' Enter the desired column index   ';
	unless ( ref($Data) =~ /Bio::ToolBox::Data/x ) {
		carp 'Must pass a Bio::ToolBox::Data object!';
		return;
	}

	# print column header names only if we have not done so before
	unless (
		# we use filename and column number as indicators
		$Data->filename eq $DATA_FILENAME
		and join( ';', $Data->list_columns ) eq $DATA_COLNAMES
		)
	{
		print " These are the columns in the file\n";
		my $i = 0;
		foreach ( $Data->list_columns ) {
			print "  $i\t$_\n";
			$i++;
		}

		# remember for next time
		$DATA_FILENAME = $Data->filename;
		$DATA_COLNAMES = join( ';', $Data->list_columns );
	}
	print $line;

	# get response
	my $response = prompt($line);
	my @indices = parse_list($response);

	# verify
	my @good;
	foreach (@indices) {
		if ( $Data->name($_) ) {
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
	$dataset =~ s/^(?: file | http | ftp ):\/*//x;

	if ( $dataset =~ /&/ ) {

		# a combination dataset
		foreach ( split /&/, $dataset ) {
			my $n = simplify_dataset_name($_);
			if ($new_name) {
				$new_name .= '&' . $n;
			}
			else {
				$new_name = $n;
			}
		}
	}
	else {
		# a single dataset
		# this could be either a file name or an entry in a BioPerl or BigWigSet database
		# remove any possible paths
		( undef, undef, $new_name ) = File::Spec->splitpath($dataset);

		# remove any known file extensions
		$new_name =~
s/\. (?: bw | bam | bb | useq | bigwig | bigbed | g[tf]f3? | cram | wig | bdg | bedgraph ) (?:\.gz)? $//xi;

		# remove common non-useful stuff
		# trying to imagine all sorts of possible things
		$new_name =~
s/[_\.\-] (?: sort | sorted | dedup | dedupe | deduplicated | rmdup | mkdup | markdup | dup | unique | filt | filtered ) \b //xgi;
		$new_name =~
s/[_\.\-] (?: coverage | rpm | ext\d* | extend\d* | log2fe | log\d+ | qvalue | fragment | count | lambda_control | fe | fold.?enrichment | ratio | log\d*ratio ) \b //xgi;
	}
	return $new_name;
}

sub sane_chromo_sort {
	my @chroms = @_;
	return unless scalar @chroms;

	# let's try and sort in some kind of rational order
	my @numeric;
	my @romanic;
	my @mixed;
	my @alphic;
	my @sex;
	my @mito;
	foreach my $c (@chroms) {

		my $name;
		if ( ref($c) eq 'ARRAY' ) {
			$name = $c->[0];
		}
		else {
			$name = $c;
		}

		# identify the type of chromosome name to sort
		if ( $name =~ m/^ (?:chr)? ( [wxyz] ) $/xi ) {

			# sex chromosomes
			push @sex, [ $1, $c ];
		}
		elsif ( $name =~ m/^ (?:chr)? (?: m | mt | mito ) (?:dna)? $/xi ) {

			# mitochondrial
			push @mito, [ $name, $c ];
		}
		elsif ( $name =~ m/^ (?:chr)? (\d+) $/xi ) {

			# standard numeric chromosome
			push @numeric, [ $1, $c ];
		}
		elsif ( $name =~ m/^ (?:chr)? ( [IVX]+ ) $/x ) {

			# Roman numerals - silly Saccharomyces cerevisiae
			push @romanic, [ $1, $c ];
		}
		elsif ( $name =~ m/^ ( [a-zA-Z_\-\.]+ ) (\d+)/x ) {

			# presumed contigs and such?
			push @mixed, [ $1, $2, $name, $c ];
		}
		else {
			# everything else
			push @alphic, [ $name, $c ];
		}
	}

	# check romanic
	if ( scalar @romanic ) {

		# looks like we have romanic chromosomes
		if ( scalar @sex ) {

			# probably caught up chrX, unlikely WYZ
			my @x = grep { $sex[$_]->[0] =~ m/^X$/ } ( 0 .. $#sex );
			foreach ( reverse @x ) {

				# I'm assuming and hoping there's only one chrX found
				# but reverse the list, just in case - assuming grep returns in order
				push @romanic, ( splice( @sex, $_, 1 ) );
			}
		}
		if ( scalar @numeric ) {

			# well, shoot, this is weird, mix of both numeric and romanic chromosomes?
			# just merge romanic with alphic and hope for the best
			push @alphic, @romanic;
		}
		else {
			# convert the romanic to numeric
			while (@romanic) {
				my $r = shift @romanic;
				my $c = $r->[0];
				$c =~ s/IV/4/;
				$c =~ s/IX/9/;
				$c =~ s/V/5/;
				$c =~ s/I/1/g;
				my $n = 0;
				foreach ( split m//, $c ) {
					if ( $_ eq 'X' ) {
						$n += 10;
					}
					else {
						$n += $_;
					}
				}
				push @numeric, [ $n, $r->[1] ];
			}
		}
	}

	# sort
	my @sorted;
	push @sorted, map { $_->[1] } sort { $a->[0] <=> $b->[0] } @numeric;
	push @sorted, map { $_->[1] } sort { $a->[0] cmp $b->[0] } @sex;
	push @sorted, map { $_->[1] } sort { $a->[0] cmp $b->[0] } @mito;
	push @sorted,
		map { $_->[3] }
		sort { $a->[0] cmp $b->[0] or $a->[1] <=> $b->[1] or $a->[2] cmp $b->[2] } @mixed;
	push @sorted, map { $_->[1] } sort { $a->[0] cmp $b->[0] } @alphic;

	return @sorted;
}

1;

__END__

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

=item sane_chromo_sort

    my @chromo = $db->seq_ids;
    my @sorted = sane_chromo_sort(@chromo);

This subroutine will take a list of chromosome or sequence identifiers and sort
them into a reasonably sane order: standard numeric identifiers first (numeric
order), sex chromosomes (alphabetical), mitochondrial, names with text and
numbers (text first alphabetically, then numbers numerically) for contigs and
such, and finally anything else (aciibetically). Any 'chr' prefix is ignored.
Roman numerals are properly handled numerically. 

The provided list may be a list of SCALAR values (chromosome names) or ARRAY 
references, with the first element assumed to be the name, e.g. 
C<[$name, $length]>. 

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

