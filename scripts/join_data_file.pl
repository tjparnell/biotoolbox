#!/usr/bin/perl

# documentation at end of file

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use IO::Prompt::Tiny qw(prompt);
use Bio::ToolBox::Data::Stream;
use Bio::ToolBox::utility qw(format_with_commas);

our $VERSION = '1.70';

print "\n This script will concatenate two or more data files\n\n";

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
my ( $outfile, $gz, $help, $print_version, );

# Command line options
GetOptions(
	'o|out=s'   => \$outfile,          # specify the input data file
	'z|gz!'     => \$gz,               # compress output files
	'h|help'    => \$help,             # request help
	'v|version' => \$print_version,    # print the version
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
	print " Biotoolbox script join_data_file.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

### Check for required values
unless ( scalar @ARGV > 1 ) {
	print " WARNING: Only one file was provided to join!\n";
	exit;
}

### Load first file
my $first_file = shift @ARGV;
my $first_data = Bio::ToolBox::Data::Stream->new( in => $first_file )
	or die "Unable to open first file '$first_file'!\n";

### Prepare output file name
# get outfile name
unless ($outfile) {

	# if the file was a set generated by split_data_file.pl
	# then it may use the # as a demarcation symbol for the basename and split value
	# look for it and regenerate the original basename
	if ( $first_data->basename =~ /^(.+) \# \w+ $/x ) {
		$outfile = $first_data->path . $1 . $first_data->extension;
	}
	else {
		# ask the user for input, what else to do?
		my $p       = ' please enter the output file name: ';
		my $default = File::Spec->catfile( $first_data->path, 
			sprintf( "%s.joined%s", $first_data->basename, $first_data->extension ) );
		$outfile = prompt( $p, $default );
	}
}

# check extension
unless ( defined $gz ) {
	if ( $outfile =~ /\.gz$/ ) {
		$gz = 1;
	}
	elsif ( $first_file =~ /\.gz$/ ) {

		# first input file is compressed, so keep it that way
		$gz = 1;
	}

	# otherwise, keep it undefined
}

### Begin writing output joined file
my $Output = $first_data->duplicate($outfile)
	or die " unable to write output file!\n";
my $line_count = 0;
while ( my $row = $first_data->next_row ) {
	$Output->add_row($row);
	$line_count++;
}
$first_data->close_fh;
printf " merged %s lines from %s\n", format_with_commas($line_count),
	$first_data->filename;

### Now write the remaining files
foreach my $file (@ARGV) {

	# open the file
	my $Data = Bio::ToolBox::Data::Stream->new( in => $file )
		or die "\n Unable to open file '$file'! Unable to proceed!\n";

	# check that file extension matches
	if ( $Output->extension ne $Data->extension ) {

		# double-check we don't have gz extensions confounding our comparison
		my $outext = $Output->extension;
		$outext =~ s/\.gz$//i;
		my $inext = $Data->extension;
		$inext =~ s/\.gz$//i;
		if ( $outext ne $inext ) {
			printf( "\n WARNING: File extensions do not match! Compare %s with %s!\n",
				$Output->extension, $Data->extension );
		}
	}

	# check for equal number of columns
	unless ( $Output->number_columns == $Data->number_columns ) {
		printf STDERR 
"\n FATAL: Column number mismatch!\n %s has %d columns instead of %d!\n Unable to proceed!\n",
			$file, $Data->number_columns, $Output->number_columns;
		exit 1;
	}

	# check first and last column names
	my @output_names = $Output->list_columns;
	my @data_names   = $Data->list_columns;
	if ( join(q( ), @output_names) ne join(q( ), @data_names) ) {
		print "\n  WARNING! Column header names don't match!!\n";
		for my $i ( 0 .. $#output_names ) {
			if ( $output_names[$i] ne $data_names[$i] ) {
				printf "  compare index $i, '%s' with '%s'\n      ", $output_names[$i],
					$data_names[$i];
			}
		}
	}

	# continue writing the file
	# for the sake of speed over data checking, we're taking raw input lines
	# and writing directly to the output
	# write to the raw output filehandle if we can get it
	# going through Data::Feature objects is considerably slower and not necessary
	my $in_fh  = $Data->fh;
	my $out_fh = $Output->fh;    # only get a filehandle if we've started writing the file
	if ($out_fh) {
		while ( my $line = $in_fh->getline ) {
			$out_fh->print($line);
			$line_count++;
		}
	}
	else {
		# not as fast but fast enough
		while ( my $line = $in_fh->getline ) {
			$Output->write_row($line);
			$line_count++;
		}
	}
	printf " merged %s lines from %s\n", format_with_commas($line_count), $Data->filename;
}

### Finish
$Output->close_fh;
print " Wrote combined file '$outfile'\n";

__END__

=head1 NAME

join_data_file.pl

A program to join two or more data files and concatenate rows.

=head1 SYNOPSIS

join_data_file.pl [--options] <file1> <file2> ...
  
  Options:
  -o --out <filename>       provide output file name, default file1
  -z --gz                   compress output
  -v --version              print version and exit
  -h --help                 show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --out E<lt>filenameE<gt>

Provide the name of the output file. If the input files were 
split using 'split_data_file.pl', then the original base name 
may be reconstituted. Otherwise, the user will be asked for 
an output file name.

=item --gz

Indicate whether the output files should be compressed 
with gzip. Default behavior is to preserve the compression 
status of the first input file.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will join two or or more data files, essentially concatanating
the files but intelligently dealing with the metadata and column headers. 
Checks are made to ensure that the number of columns in the subsequent files 
match the first file.

The program will not merge datasets from multiple files; see 
the program 'merge_datasets.pl' for that.

This program is intended as the complement to 'split_data_files.pl'.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
