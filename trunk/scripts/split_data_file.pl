#!/usr/bin/env perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::ToolBox::file_helper qw(
	open_tim_data_file
	write_tim_data_file
	open_to_write_fh
);
my $VERSION = '1.14';

print "\n This script will split a data file by features\n\n";


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
	$infile, 
	$index,
	$max,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'        => \$infile, # specify the input data file
	'index|col=i' => \$index, # index for the column to use for splitting
	'max=i'       => \$max, # maximum number of lines per file
	'gz!'         => \$gz, # compress output files
	'help'        => \$help, # request help
	'version'     => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script split_data_file.pl, version $VERSION\n\n";
	exit;
}





### Check for required values
unless ($infile) {
	$infile = shift @ARGV or
		die "  OOPS! No source data file specified! \n use $0 --help\n";
}
unless (defined $gz) {
	if ($infile =~ /\.gz$/) {
		# input file is compressed, so keep it that way
		$gz = 1;
	}
	else {
		$gz = 0;
	}
}



### Load file
my ($in_fh, $metadata_ref) = open_tim_data_file($infile);
unless ($in_fh) {
	die "Unable to open data table!\n";
}

# add column headers
push @{ $metadata_ref->{'data_table'} }, $metadata_ref->{'column_names'};
$metadata_ref->{'last_row'} = 0;



### Identify the column
# Request the index from the user if necessary
unless (defined $index) {
	
	# present a list to the user
	print "\n  This is the list of columns in the data file\n";
	for (my $i = 0; $i < $metadata_ref->{'number_columns'}; $i++) {
		print "   $i\t", $metadata_ref->{$i}{'name'}, "\n";
	}
	
	# get user response
	print "  Enter the column index number containing the values to split by   ";
	$index = <STDIN>;
	chomp $index;
}

# check index
unless ($index =~ /^\d+$/ and exists $metadata_ref->{$index}) {
	die " unknown column index!\n";
}



### Split the file
print " Splitting file by elements in column '$metadata_ref->{$index}{name}'...\n";
my %written_files; # a hash of the file names written
	# we can't assume that all the data elements we're splitting on are 
	# contiguous in the file
	# if they're not, then we would be simply re-writing over the 
	# previous block
	# also, we're enforcing a maximum number of lines per file
	# so we'll remember the files we've written, and re-open that file 
	# to write the next block of data
my $previous_value;
my $split_count = 0;
while (my $line = $in_fh->getline) {
	
	# Collect line data and the check value
	chomp $line;
	my @data = split /\t/, $line;
	my $check_value = $data[$index];
	
	# For the first line only
	unless (defined $previous_value) {
		$previous_value = $check_value;
	}
	
	# Determine whether to write or proceed to next line
	if ($check_value eq $previous_value ) {
		# the same value, so keep in same array
		
		push @{ $metadata_ref->{'data_table'} }, [ @data ];
		$metadata_ref->{'last_row'} += 1;
	}
	else {
		# different value, new data section
		
		# write the current data out to file(s)
		write_current_data_to_file_part($previous_value);
			# this should automatically clear the previous table data
		
		# now add the current row of data
		push @{ $metadata_ref->{'data_table'} }, [ @data ];
		$metadata_ref->{'last_row'} += 1;
		
		# reset
		$previous_value = $check_value;
	}
	
	# Check the number of lines collected, write if necessary
	if (defined $max and $metadata_ref->{'last_row'} == $max) {
		# we've reached the maximum number of data lines for this current 
		# data block
		
		# need to force a data write
		write_current_data_to_file_part($previous_value);
	}
}



### Finish

# Final write 
write_current_data_to_file_part($previous_value);

# Properly close out all file handles
$in_fh->close;
foreach my $value (keys %written_files) {
	$written_files{$value}{'fh'}->close;
}


# report
print " Split '$infile' into $split_count files\n";
foreach my $value (sort {$a cmp $b} keys %written_files) {
	print "  wrote $written_files{$value}{total} lines in " . 
		"$written_files{$value}{parts} files for '$value'\n";
}




sub write_current_data_to_file_part {

	# get the current value we're working with
	my $value = shift;
	my $last_row = $metadata_ref->{'last_row'};
	
	# open the file and write
		# check for a pre-existing file to be added, or start a new one
	if (defined $written_files{$value}{'file'}) {
		# we have a current file that is partially written
		my $file = $written_files{$value}{'file'};
		my $out_fh = $written_files{$value}{'fh'};
		
		# begin writing the data
		# determine how many lines we can still write to this file
		
		# maximum is defined but we'll have to do two writes
		if (
			defined $max and
			($max - $written_files{$value}{'number'}) < $last_row
		) {
			# we'll have to do two writes
			# finish up current, then write the remainder
			my $limit = $max - $written_files{$value}{'number'};
			
			# write the lines up to the current limit
			for (my $row = 1; $row <= $limit; $row++) {
				$out_fh->print(
					join("\t", @{ $metadata_ref->{'data_table'}->[$row] } ) .
				"\n");
			}
			$written_files{$value}{'total'} += $limit;
			
			# clear the table contents of the written lines
			splice( @{ $metadata_ref->{'data_table'} }, 1, $limit );
			$metadata_ref->{'last_row'} = 
				scalar @{ $metadata_ref->{'data_table'} } - 1;
			
			# we're finished with this file
			$written_files{$value}{'file'}   = undef;
			$written_files{$value}{'fh'}     = undef;
			$written_files{$value}{'number'} = 0;
			
			# now write the remainder
			write_current_data_to_file_part($value);
		}
		
		# maximum is either not defined, or there is enough room left to write
		else {
			# we can write with impunity
			
			# write the lines
			for (my $row = 1; $row <= $last_row; $row++) {
				$out_fh->print(
					join("\t", @{ $metadata_ref->{'data_table'}->[$row] } ) .
				"\n");
				
				# keep track of the number of lines
				$written_files{$value}{'number'} += 1;
				$written_files{$value}{'total'} += 1;
			}
			
			# clear the table contents
			splice( @{ $metadata_ref->{'data_table'} }, 1, $last_row );
			$metadata_ref->{'last_row'} = 0;
		}
	}
	
	else {
		# put in a request for a file name
		request_new_file_name($value);
		
		# write the file
		# this should be within the maximum line limit, so we should be safe
		my $success = write_tim_data_file(
			'data'     => $metadata_ref,
			'filename' => $written_files{$value}{'file'},
			'gz'       => $gz,
		);
		if ($success) {
			# record the number of lines written
			$written_files{$value}{'number'} += $last_row;
			$written_files{$value}{'total'} += $last_row;
			
			# update the file name in case it was changed by the write method
			$written_files{$value}{'file'} = $success;
		}
		else {
			warn "   unable to write $last_row lines! data lost!\n";
		}
		
		# clear the table contents
		splice( @{ $metadata_ref->{'data_table'} }, 1, $last_row );
		$metadata_ref->{'last_row'} = 0;
		
		# check whether we've filled up the file
		if (defined $max and $last_row == $max) {
			$written_files{$value}{'file'} = undef;
			$written_files{$value}{'number'} = 0;
		}
		else {
			# reopen the file for future writing
			$written_files{$value}{'fh'} = open_to_write_fh(
				$written_files{$value}{'file'},
				$gz,
				1
			);
		}
	}
	
}

sub request_new_file_name {
	# calculate a new file name based on the current check value and part number
	my $value = shift;
	my $filename_value = $value;
	$filename_value =~ s/[\:\|\\\/\+\*\?\# ]+/_/g; # replace unsafe characters
	my $file = $metadata_ref->{'path'} . $metadata_ref->{'basename'} . 
		'#' . $filename_value;
	
	# add the file part number, if we're working with maximum line files
	# padded for proper sorting
	if (defined $max) {
		if (defined $written_files{$value}{'parts'}) {
			$written_files{$value}{'parts'} += 1; # increment
			$file .= '_' . sprintf("%03d", $written_files{$value}{'parts'});
		}
		else {
			$written_files{$value}{'parts'} = 1; # initial
			$file .= '_' . sprintf("%03d", $written_files{$value}{'parts'});
		}
	}
	else {
		# only 1 part is necessary
		$written_files{$value}{'parts'} = 1;
	}
	
	# finish the file name
	$file .= $metadata_ref->{'extension'};
	$written_files{$value}{'file'}   = $file;
	$written_files{$value}{'number'} = 0;
	
	# check the total
	unless (exists $written_files{$value}{'total'}) {
		$written_files{$value}{'total'}  = 0;
	}
	
	# keept track of the number of files opened
	$split_count++;
}



__END__

=head1 NAME

split_data_file.pl

A script to split a data file by rows based on common data values.

=head1 SYNOPSIS

split_data_file.pl [--options] <filename>
  
  Options:
  --in <filename>
  --index <column_index>
  --max <integer>
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a data file. It must be a tab-delimited text file,
preferably in the tim data format as described in Bio::ToolBox::file_helper, 
although any format should work. The file may be compressed with gzip.

=item --index <column_index>

=item --col <column_index>

Provide the index number of the column or dataset containing the values 
used to split the file. If not specified, then the index is requested 
from the user in an interactive mode.

=item --max <integer>

Optionally specify the maximum number of data lines to write to each 
file. Each group of specific value data is written to one or more files. 
Enter as an integer; underscores may be used as thousands separator, e.g. 
100_000. 

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
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
