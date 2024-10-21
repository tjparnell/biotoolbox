#!/usr/bin/perl

# documentation at end of file

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case bundling);
use Bio::ToolBox::Data::Stream;
use Bio::ToolBox::utility qw(ask_user_for_index format_with_commas);

our $VERSION = '2.01';

print "\n This script will split a data file by features\n\n";

### Quick help
unless (@ARGV) {

	# when no command line options are present
	# print SYNOPSIS
	pod2usage(
		{
			'-verbose' => 0,
			'-exitval' => 1,
		}
	);
}

### Get command line options and initialize values
my ( $infile, $index, $tag, $max, $gz, $prefix, $noheader, $help, $print_version, );

# Command line options
GetOptions(
	'i|in=s'        => \$infile,           # specify the input data file
	'x|index|col=i' => \$index,            # index for the column to use for splitting
	't|tag=s'       => \$tag,              # attribute tag name
	'm|max=i'       => \$max,              # maximum number of lines per file
	'p|prefix=s'    => \$prefix,           # output file prefix
	'H|noheader'    => \$noheader,         # file has no headers
	'z|gz!'         => \$gz,               # compress output files
	'h|help'        => \$help,             # request help
	'v|version'     => \$print_version,    # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {

	# print entire POD
	pod2usage(
		{
			'-verbose' => 2,
			'-exitval' => 1,
		}
	);
}

# Print version
if ($print_version) {
	print " Biotoolbox script split_data_file.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

### Check for required values
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		print STDERR " FATAL: no input file! use --help for more information\n";
		exit 1;
	}
}
unless ( defined $gz ) {
	if ( $infile =~ /\.gz$/ ) {

		# input file is compressed, so keep it that way
		$gz = 1;
	}
	else {
		$gz = 0;
	}
}
$noheader = defined $noheader ? $noheader : 0;

### Load Input file
my $Input = Bio::ToolBox::Data::Stream->new( in => $infile, noheader => $noheader )
	or die "Unable to open input file!\n";

# Identify the column
unless ( defined $index or defined $tag ) {
	$index = ask_user_for_index( $Input,
		"  Enter the column index number containing the values to split by   " );
	unless ( defined $index ) {
		print STDERR " FATAL: Must provide a valid index!\n";
		exit 1;
	}
}
if ($tag) {
	unless ( $Input->gff or $Input->vcf ) {
		print STDERR
			" FATAL: Input file must be in GFF or VCF format to use attribute tags!";
		exit 1;
	}
	if ( $Input->vcf and not defined $index ) {
		print STDERR
			" FATAL: Please provide a column index for accessing VCF attributes.\n"
			. " The INFO column is index 8, and sample columns begin\n"
			. " at index 10.\n";
		exit 1;
	}
	elsif ( $Input->gff ) {
		$index = 8;
	}
}

### Split the file
printf " Splitting file by elements in column %s%s...\n",
	$Input->name($index),
	$tag ? ', attribute tag $tag' : q();
my %out_files;    # a hash of the file names written
                  # we can't assume that all the data elements we're splitting on are
                  # contiguous in the file
                  # if they're not, then we would be simply re-writing over the
                  # previous block
                  # also, we're enforcing a maximum number of lines per file
                  # so we'll remember the files we've written, and re-open that file
                  # to write the next block of data
my $split_count = 0;

while ( my $row = $Input->next_row ) {

	# Get the check value
	my $check;
	if ($tag) {
		my $attrib = $row->attributes;
		$check = $attrib->{$tag} || $attrib->{$index}{$tag} || undef;
	}
	else {
		$check = $row->value($index);
	}
	unless ( exists $out_files{$check}{'stream'} ) {
		request_new_file_name($check);
	}

	# write the row
	$out_files{$check}{'stream'}->add_row($row);
	$out_files{$check}{'number'} += 1;
	$out_files{$check}{'total'}  += 1;

	# Check the number of lines collected, close if necessary
	if ( defined $max and $out_files{$check}{'number'} == $max ) {

		# we've reached the maximum number of data lines for this current data
		$out_files{$check}{'stream'}->close_fh;
		delete $out_files{$check}{'stream'};
	}
}

### Finish
# Properly close out all file handles
$Input->close_fh;
foreach my $value ( keys %out_files ) {
	$out_files{$value}{'stream'}->close_fh if exists $out_files{$value}{'stream'};
}

# report
print " Split '$infile' into $split_count files\n";
foreach my $value ( sort { $a cmp $b } keys %out_files ) {
	printf "  wrote %s lines in %d file%s for '$value'\n",
		format_with_commas( $out_files{$value}{total} ), $out_files{$value}{parts},
		$out_files{$value}{parts} > 1 ? 's' : q();
}

sub request_new_file_name {

	# calculate a new file name based on the current check value and part number
	my $value          = shift;
	my $filename_value = $value;
	$filename_value =~ s/[\: \| \\ \/ \+ \* \? \# \( \) \[ \] \{ \} \  ]+/_/xg;

	# replace unsafe characters

	my $file;
	if ( $prefix and $prefix eq 'none' ) {
		$file = $Input->path . $filename_value;
	}
	elsif ($prefix) {
		$file = $prefix . '#' . $filename_value;
	}
	else {
		$file = $Input->path . $Input->basename . '#' . $filename_value;
	}

	# add the file part number, if we're working with maximum line files
	# padded for proper sorting
	if ( defined $max ) {
		if ( defined $out_files{$value}{'parts'} ) {
			$out_files{$value}{'parts'} += 1;    # increment
			$file .= '_' . sprintf( "%03d", $out_files{$value}{'parts'} );
		}
		else {
			$out_files{$value}{'parts'} = 1;     # initial
			$file .= '_' . sprintf( "%03d", $out_files{$value}{'parts'} );
		}
	}
	else {
		# only 1 part is necessary
		$out_files{$value}{'parts'} = 1;
	}

	# finish the file name
	$file .= $Input->extension;
	$out_files{$value}{'number'} = 0;

	# open an output Stream
	if ( exists $out_files{$value}{'stream'} ) {

		# an open stream, close it
		$out_files{$value}{'stream'}->close_fh;
	}
	my $Stream = $Input->duplicate($file);
	$out_files{$value}{'stream'} = $Stream;

	# check the total
	unless ( exists $out_files{$value}{'total'} ) {
		$out_files{$value}{'total'} = 0;
	}

	# keept track of the number of files opened
	$split_count++;
}

__END__

=head1 NAME

split_data_file.pl

A program to split a data file by rows based on common data values.

=head1 SYNOPSIS

split_data_file.pl [--options] <filename>
  
  File options:
  -i --in <filename>                (txt bed gff gtf vcf refFlat ucsc etc)
  -p --prefix <text>                output file prefix (input basename)
  -H --noheader                     input file has no headers
  
  Splitting options:
  -x --index <column_index>         column with values to split upon
  -t --tag <text>                   use VCF/GFF attribute
  -m --max <integer>                maximum number of items per output file
  
  General options:
  -z --gz                           compress output file
  -v --version                      print version and exit
  -h --help                         show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 File options

=over 4

=item --in E<lt>filenameE<gt>

Specify the file name of a data file. It must be a tab-delimited text file. 
The file may be compressed with gzip.

=item --prefix E<lt>textE<gt>

Optionally provide a filename prefix for the output files. The default 
prefix is the input filename base name. If no prefix is desired, using 
just the values as filenames, then set the prefix to 'none'.

=item --noheader

Indicate that the input file has no column header line, and that dummy 
headers will be provided. Not necessary for BED, GFF, or recognized UCSC 
file formats.

=back

=head2 Splitting options

=over 4

=item --index E<lt>column_indexE<gt>

Provide the index number of the column or dataset containing the values 
used to split the file. If not specified, then the index is requested 
from the user in an interactive mode.

=item --tag E<lt>textE<gt>

Provide the attribute tag name that contains the values to split the 
file. Attributes are supported by GFF and VCF files. If splitting a 
VCF file, please also provide the column index. The INFO column is 
index 8, and sample columns begin at index 10.

=item --max E<lt>integerE<gt>

Optionally specify the maximum number of data lines to write to each 
file. Each group of specific value data is written to one or more files. 
Enter as an integer; underscores may be used as thousands separator, e.g. 
100_000. 

=back

=head2 General options

=over 4

=item --gz

Indicate whether the output files should be compressed 
with gzip. Default behavior is to preserve the compression 
status of the input file.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will split a data file into multiple files based on common 
values in the data table. All rows with the same value will be 
written into the same file. A good example is chromosome, where all 
data points for a given chromosome will be written to a separate file, 
resulting in multiple files representing each chromosome found in the 
original file. The column containing the values to split and group 
should be indicated; if the column is not sepcified, it may be 
selected interactively from a list of column headers. 

This program can also split files based on an attribute tag in GFF or 
VCF files. Attributes are often specially formatted delimited key value 
pairs associated with each feature in the file. Provide the name of the 
attribute tag to split the file. Since attributes may vary based on 
the feature type, an interactive list is not supplied from which to 
choose the attribute.

If the max argument is set, then each group will be written to one or 
more files, with each file having no more than the indicated maximum 
number of data lines. This is useful to keep the file size reasonable, 
especially when processing the files further and free memory is 
constrained. A reasonable limit may be 100K or 1M lines.

The resulting files will be named using the basename of the input file, 
appended with the unique group value (for example, the chromosome name)
demarcated with a #. If a maximum line limit is set, then the file part 
number is appended to the basename, padded with zeros to three digits 
(to assist in sorting). Each file will have duplicated and preserved 
metadata. The original file is preserved.

This program is intended as the complement to 'join_data_files.pl'.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
