#!/usr/bin/perl

# A script to join two or more tim data files

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper qw(
	open_tim_data_file
	write_tim_data_file
	open_to_write_fh
);

print "\n This script will join two or more data files\n\n";


### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Get command line options and initialize values
my (
	$outfile,
	$gz,
	$help
);

# Command line options
GetOptions( 
	'out=s'     => \$outfile, # specify the input data file
	'gz!'       => \$gz, # compress output files
	'help'      => \$help # request help
);

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}




### Check for required values
unless (scalar @ARGV > 1) {
	die "  OOPS! Two or more data files must be given!\n use $0 --help\n";
}
unless (defined $gz) {
	if ($ARGV[0] =~ /\.gz$/) {
		# first input file is compressed, so keep it that way
		$gz = 1;
	}
	else {
		$gz = 0;
	}
}




### Load first file
my $first_file = shift @ARGV;
print " Joining file '$first_file'...  ";
my ($in_fh, $metadata_ref) = open_tim_data_file($first_file);
unless ($in_fh) {
	die "Unable to open first file '$first_file'!\n";
}

# generate data table
$metadata_ref->{'data_table'} = [];

# add column headers
push @{ $metadata_ref->{'data_table'} }, $metadata_ref->{'column_names'};





### Begin writing file
unless ($outfile) {
	# if the file was a set generated by split_data_file.pl
	# then it may use the # as a demarcation symbol for the basename and split
	# value
	# look for it and regenerate the original basename
	if ($metadata_ref->{'basename'} =~ /^(.+)\#\w+$/) {
		$outfile = $1;
	}
	else {
		# ask the user for input, what else to do?
		print " please enter the output file name   ";
		$outfile = <STDIN>;
		chomp $outfile;
	}
}

# First write the metadata file
my $new_outfile = write_tim_data_file( {
	'data'     => $metadata_ref,
	'filename' => $outfile,
	'gz'       => $gz,
} );
unless ($new_outfile) {
	die " unable to write output file!\n";
}

# Now reopen for appended writing
my $out_fh = open_to_write_fh($new_outfile, $gz, 1);

# Continue writing the first file
my $line_count = 0;
while (my $line = $in_fh->getline) {
	print {$out_fh} $line;
	$line_count++;
}
$in_fh->close;
print "$line_count data lines merged\n";




### Now write the remaining files
foreach my $file (@ARGV) {
	
	print " Joining file '$file'...  ";
	
	# open the file
	my ($file_fh, $file_data_ref) = open_tim_data_file($file);
	unless ($file_fh) {
		die " Unable to open file '$file'! Unable to proceed!\n";
	}
	
	# check that metadata matches
	unless (
		$file_data_ref->{number_columns} == 
			$metadata_ref->{number_columns} 
	) {
		die " Number of file columns don't match! Unable to proceed!\n";
	}
	unless (
		$file_data_ref->{column_names}->[0] eq 
			$metadata_ref->{column_names}->[0]
		and
		$file_data_ref->{column_names}->[-1] eq 
			$metadata_ref->{column_names}->[-1]
	) {
		die " Column headers don't match! Unable to proceed!\n";
	}
	
	# continue writing the file
	while (my $line = $file_fh->getline) {
		print {$out_fh} $line;
		$line_count++;
	}
	$file_fh->close;
	$file_fh = undef;
	print "$line_count data lines merged\n";
}




### Finish
$out_fh->close;
print " Wrote combined file '$new_outfile'\n";




__END__

=head1 NAME

join_data_file.pl

=head1 SYNOPSIS

join_data_file.pl [--options] <file1> <file2> ...
  
  --out <filename>
  --(no)gz
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --out <filename>

Provide the name of the output file. If the input files were 
split using 'split_data_file.pl', then the original base name 
may be reconstituted. Otherwise, the user will be asked for 
an output file name.

=item --(no)gz

Indicate whether the output files should be compressed 
with gzip. Default behavior is to preserve the compression 
status of the first input file.

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












