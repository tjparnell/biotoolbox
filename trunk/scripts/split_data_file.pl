#!/usr/bin/perl

# A script to split a tim data file based on common unique data values

use strict;
use Pod::Usage;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper;

print "\n This script will split a data file\n\n";


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
	$gz,
	$help
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # specify the input data file
	'col=i'     => \$index, # index for the column to use for splitting
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
unless ($infile) {
	$infile = shift @ARGV or
		die "  OOPS! No source data file specified! \n use $0 --help\n";
}
unless (defined $index) {
	# default is to use the first column
	$index = 0;
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
if ($infile =~ /\.store(?:\.gz)$/) {
	die "Unable to split a binary store file! convert to text first\n";
}
my ($in_fh, $metadata_ref) = open_tim_data_file($infile);
unless ($in_fh) {
	die "Unable to open data table!\n";
}

# generate data table
$metadata_ref->{'data_table'} = [];

# add column headers
push @{ $metadata_ref->{'data_table'} }, $metadata_ref->{'column_names'};





### Split the file
print " Splitting file by elements in column '$metadata_ref->{$index}{name}'...\n";
my %written_files; # a hash of the file names written
	# we can't assume that all the data elements we're splitting on are 
	# contiguous in the file
	# if they're not, then we would be simply re-writing over the 
	# previous block
	# instead, we'll remember the files we've written, and re-open that file 
	# to write the next block of data
my $previous_value;
my $split_count = 0;
my $line_count = 0;
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
		$line_count++;
	}
	else {
		# different value, new data section
		
		# update last_row index
		$metadata_ref->{'last_row'} = 
			scalar( @{ $metadata_ref->{'data_table'} } ) - 1;
		
		# write the current data to file
		my $outfile = $metadata_ref->{'basename'} . '#' . $previous_value;
		if (exists $written_files{$outfile}) {
			# this set of data is part of another block of identical data
			my $out_fh = open_to_write_fh(
				$written_files{$outfile},
				$gz,
				1
			) or warn "   unable to re-open file! data lost!\n";
			
			# add to the data table
			for (my $row = 1; $row <= $metadata_ref->{'last_row'}; $row++ ) {
				# print each data row, skipping the header
				print {$out_fh} join("\t", 
					@{ $metadata_ref->{'data_table'}->[$row] } ) . "\n";
			}
			
			# finish
			$out_fh->close;
			print "   wrote $line_count additional lines to file '" . 
				$written_files{$outfile} . "'\n";
		}
		
		else {
			# write a new file for this block of data
			my $success_write = write_tim_data_file( {
				'data'     => $metadata_ref,
				'filename' => $outfile,
				'gz'       => $gz,
			} );
			if ($success_write) {
				print "   wrote $line_count lines to file '$success_write'\n";
				# remember the file
				$written_files{$outfile} = $success_write;
					# these file names may not be identical due to extensions
			}
			else {
				warn "   unable to write $line_count lines!\n";
			}
			$split_count++;
		}
		
		# regenerate the data table
		$metadata_ref->{'data_table'} = [];
		push @{ $metadata_ref->{'data_table'} }, $metadata_ref->{'column_names'};
		$metadata_ref->{'last_row'} = undef;
		
		# now add the current row of data
		push @{ $metadata_ref->{'data_table'} }, [ @data ];
		
		# reset
		$previous_value = $check_value;
		$line_count = 1;
	}
}

### Finish
$in_fh->close;

# Final write 
{
	# update last_row index
	$metadata_ref->{'last_row'} = 
		scalar( @{ $metadata_ref->{'data_table'} } ) - 1;
	
	my $outfile = $metadata_ref->{'basename'} . '#' . $previous_value;
	if (exists $written_files{$outfile}) {
		# this set of data is part of another block of identical data
		my $out_fh = open_to_write_fh(
			$written_files{$outfile},
			$gz,
			1
		) or warn "   unable to re-open file! data lost!\n";
		
		# add to the data table
		for (my $row = 1; $row <= $metadata_ref->{'last_row'}; $row++ ) {
			# print each data row, skipping the header
			print {$out_fh} join("\t", 
				@{ $metadata_ref->{'data_table'}->[$row] } ) . "\n";
		}
		
		# finish
		$out_fh->close;
		print "   wrote $line_count additional lines to file '" . 
			$written_files{$outfile} . "'\n";
	}
	
	else {
		# write a new file for this block of data
		my $success_write = write_tim_data_file( {
			'data'     => $metadata_ref,
			'filename' => $outfile,
			'gz'       => $gz,
		} );
		if ($success_write) {
			print "   wrote $line_count lines to file '$success_write'\n";
			# remember the file
			$written_files{$outfile} = $success_write;
				# these file names may not be identical due to extensions
		}
		else {
			warn "   unable to write $line_count lines!\n";
		}
		$split_count++;
	}
	
}

print " Split '$infile' into $split_count files\n";




__END__

=head1 NAME

split_data_file.pl

=head1 SYNOPSIS

split_data_file.pl --col <column_index> <filename>
  
  --in <filename>
  --col <column_index>
  --(no)gz
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a data file. It must be a tab-delimited text file,
preferably in the tim data format as described in 'tim_file_helper.pm', 
although any format should work. The file may be compressed with gzip.

=item --col <column_index>

Provide the index number of the column or dataset 
containing the values used to split the file. Default 
index is 0 (first column).

=item --(no)gz

Indicate whether the output files should be compressed 
with gzip. Default behavior is to preserve the compression 
status of the input file.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will split a data file into multiple files based on common 
unique values in the data table. All rows with the same value will be 
written into the same file. A good example is chromosome, where all 
data points for a given chromosome will be written to a separate file, 
resulting in multiple files representing each chromosome found in the 
original file. The column containing the values to split and group 
should be indicated; the default is to take the first column. 

NOTE: The program assumes that the file is sorted based on this column; 
splitting an unsorted file will likely yield numerous, small, and 
overwritten files with likely data loss.

The resulting files will be named using the basename of the input file, 
appended with the unique group value (for example, the chromosome name)
demarcated with a #. Each file will have duplicated and preserved 
metadata. The original file is preserved.

This program is intended as the complement to 'join_data_files.pl'.



=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112












