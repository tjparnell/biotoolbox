#!/usr/bin/perl

# A program to merge two or more dataset files


use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	generate_tim_data_structure
	find_column_index
);
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
);
my $VERSION = '1.8.5';

print "\n A progam to merge datasets from two files\n";


# Print help
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
	$lookup,
	$automatic,
	$outfile,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'lookup!'   => \$lookup, # force merging by value lookup
	'auto!'     => \$automatic, # select columns automatically
	'out=s'     => \$outfile, # name of output file 
	'gz!'       => \$gz, # compress output
	'help'      => \$help, # request help
	'version'   => \$print_version, # print the version
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
	print " Biotoolbox script merge_datasets.pl, version $VERSION\n\n";
	exit;
}




### Check for requirements and set defaults

# input files
unless (scalar @ARGV >= 2) {
	die " At least two input files are required!\n";
}

# automatic mode
if ($automatic) {
	unless ($outfile) {
		die " Must provide an output file name in automatic mode!\n";
	}
}

# Get the alphabet-number conversion hashes
my %letter_of = get_letters(); # convert numbers into letters
my %number_of = get_numbers(); # to convert letters into numbers

# Set up output
my $output_data_ref; # the reference scalar for the output data structure

# name of lookup column to be used for all files
my $lookup_name;


### Process and merge the files

# exactly two filenames are provided
if (scalar @ARGV == 2) {
	my $input_data1_ref = read_file(shift @ARGV);
	my $input_data2_ref = read_file(shift @ARGV);
	merge_two_datasets($input_data1_ref, $input_data2_ref);
}

# more than two filenames
else {
	# merge the first two
	my $input_data1_ref = read_file(shift @ARGV);
	my $input_data2_ref = read_file(shift @ARGV);
	merge_two_datasets($input_data1_ref, $input_data2_ref);
	$input_data1_ref = undef;
	$input_data2_ref = undef;
	
	# merge the subsequent files
	foreach (@ARGV) {
		my $input_data_ref = read_file($_);
		add_datasets($input_data_ref);
	}
}




### Offer to re-name the datasets
unless ($automatic) {
	print " Do you wish to rename the headers? y or n   ";
	my $response = <STDIN>;
	chomp $response;
	
	if ($response eq 'y' or $response eq 'Y') { 
		# we will be re-naming headers
		rename_dataset_names();
	}
}



### Print the output

# Request the output file name
# default is to simply overwrite file1
unless ($outfile) {
	print " Enter the output file name [", $output_data_ref->{'filename'}, "] ";
	$outfile = <STDIN>;
	chomp $outfile;
	if ($outfile eq '') {
		# use file 1 as the default file name
		$outfile = $output_data_ref->{'filename'};
	} 
}

# write the file
my $file_written = write_tim_data_file( {
	'data'      => $output_data_ref,
	'filename'  => $outfile,
	'gz'        => $gz,
} );
if ($file_written) {
	print " Wrote file '$file_written'\n";
}
else {
	print " No file written!\n";
}





########################   Subroutines   ###################################

### Read the input data file
sub read_file {
	# subroutine to load the time data file and print basic stats about it
	my $filename = shift;
	
	# load the file data
	my $file_data_ref = load_tim_data_file($filename);
	unless ($file_data_ref) {
		die " file '$filename' not loaded!";
	}
	
	# print the results
	print "\n Loaded '$filename' with ", $file_data_ref->{'last_row'}, 
		" data rows and ", $file_data_ref->{'number_columns'}, " columns\n";
	
	# return
	return $file_data_ref;
}


### Merge two datasets together
sub merge_two_datasets {
	my ($input_data1_ref, $input_data2_ref) = @_;
	
	# Check the input data
	my $check = check_data_tables($input_data1_ref, $input_data2_ref);
	if ($check and !$lookup) {
		print " Files have non-equal numbers of data rows! Enabling lookup\n";
	}
	
	
	# Merge by lookup values
	if ($lookup or $check) {
		# we need to merge by lookup values
		merge_two_datasets_by_lookup($input_data1_ref, $input_data2_ref);
		return;
	}
	
	
	# Otherwise Merge the two datasets blindly
	# determine the new order
	my @order;
	if ($automatic) {
		# automatic selection
		@order = automatically_determine_order(
			$input_data1_ref, $input_data2_ref);
		
	}
	else {
		# manual selection from user
		@order = request_new_order($input_data1_ref, $input_data2_ref);
	}
	
	# Initialize the output data structure if necessary
	$output_data_ref = initialize_output_data_structure($input_data1_ref); 

	# assign the datasets to the output data in requested order
	foreach my $request (@order) {
		
		# determine the current column index we're working with
		my $column = $output_data_ref->{'number_columns'};
		
		# add the dataset from the appropriate input file
		if ($request =~ /^\d+$/) {
			# a digit indicates a dataset from file1
			
			# print dataset name in automatic mode
			if ($automatic) {
				print "  Merging column " . 
					$input_data1_ref->{$request}{'name'} . "\n";
			}
			
			# copy the dataset
			for my $row (0 .. $input_data1_ref->{'last_row'}) {
				$output_data_ref->{'data_table'}->[$row][$column] = 
					$input_data1_ref->{'data_table'}->[$row][$request];
			}
			
			# copy the metadata
			copy_metadata($input_data1_ref, $request);
		} 
		elsif ($request =~ /^[a-z]+$/i) {
			# a letter indicates a dataset from file2
			
			# first convert back to number
			my $number = $number_of{$request};
			
			# print dataset name in automatic mode
			if ($automatic) {
				print "  Merging column " . 
					$input_data2_ref->{$number}{'name'} . "\n";
			}
			
			# copy the dataset
			for my $row (0 .. $input_data2_ref->{'last_row'}) {
				$output_data_ref->{'data_table'}->[$row][$column] = 
					$input_data2_ref->{'data_table'}->[$row][$number];
			}
			
			# copy the metadata
			copy_metadata($input_data2_ref, $number);
		} 
		else {
			die " unrecognized symbol '$request' in request! nothing done!\n";
		}
	}
}



### merge two datasets by lookup value
sub merge_two_datasets_by_lookup {
	my ($input_data1_ref, $input_data2_ref) = @_;
	
	# determine lookup indices
	my ($lookup_i1, $lookup_i2) = 
		request_lookup_indices($input_data1_ref, $input_data2_ref);
		
	# determine order
	my @order;
	if ($automatic) {
		# automatic selection
		@order = automatically_determine_order(
			$input_data1_ref, $input_data2_ref);
		
	}
	else {
		# manual selection from user
		@order = request_new_order($input_data1_ref, $input_data2_ref);
	}
	
	# rearrange as necessary to make the first data structure dominant
	if ($order[0] =~ /[a-z]+/i) {
		# it looks like the user wants the second file to be dominant
		# the dominant file is where we'll be taking all of the values
		
		# switch the references
		($input_data1_ref, $input_data2_ref) = 
			($input_data2_ref, $input_data1_ref);
		
		# switch the lookup indices
		($lookup_i1, $lookup_i2) = ($lookup_i2, $lookup_i1);
		
		# switch the numbers and letters
		map {
			if (/\d+/) {
				$_ = $letter_of{$_};
			}
			else {
				$_ = $number_of{$_};
			}
		} @order;
	}
	
	# Index the second dataset
		# we'll be putting the lookup values into an index hash
		# where the lookup value is the key and the row number is the value
	$input_data2_ref->{'index'} = {};
	for (my $row = 1; $row <= $input_data2_ref->{'last_row'}; $row++) {
		my $key = $input_data2_ref->{'data_table'}->[$row][$lookup_i2];
		if (exists $input_data2_ref->{'index'}{$key}) {
			# value is not unique
			warn " lookup value '$key' in file " . 
				$input_data2_ref->{'filename'} . 
				", row $row is a duplicate!\n" . 
				" Using the first occurence value\n";
		}
		else {
			# value is ok
			$input_data2_ref->{'index'}{$key} = $row;
		}
	}
	
	# Initialize the output data structure if necessary
	$output_data_ref = initialize_output_data_structure($input_data1_ref); 
	
	# assign the datasets to the output data in requested order
	foreach my $request (@order) {
		
		# determine the current column index we're working with
		my $column = $output_data_ref->{'number_columns'};
		
		# add the dataset from the appropriate input file
		if ($request =~ /\d+/) {
			# a digit indicates a dataset from file1
			# we're assuming that file1 is dominant, we copy all the 
			# rows from file1 into output, no lookup required
			
			# print dataset name in automatic mode
			if ($automatic) {
				print "  Merging column " . 
					$input_data1_ref->{$request}{'name'} . "\n";
			}
			
			# copy the dataset including header
			for my $row (0 .. $input_data1_ref->{'last_row'}) {
				$output_data_ref->{'data_table'}->[$row][$column] = 
					$input_data1_ref->{'data_table'}->[$row][$request];
				
			}
			
			# copy the metadata
			copy_metadata($input_data1_ref, $request);
		} 
		
		elsif ($request =~ /[a-z]+/i) {
			# a letter indicates a dataset from file2
			# we will have to perform the lookup here
			
			# first convert back to number
			my $request_number = $number_of{$request};
			
			# print dataset name in automatic mode
			if ($automatic) {
				print "  Merging column " . 
					$input_data2_ref->{$request_number}{'name'} . "\n";
			}
			
			# copy the header
			$output_data_ref->{'data_table'}->[0][$column] = 
				$input_data2_ref->{'data_table'}->[0][$request_number];
			
			# copy the dataset
			for my $row (1 .. $input_data1_ref->{'last_row'}) {
				# identify the appropriate row in file2 by lookup value
				my $lookup = $input_data1_ref->{'data_table'}->[$row][$lookup_i1];
				my $row2 = $input_data2_ref->{'index'}{$lookup} || undef;
				
				# copy the appropriate value
				if (defined $row2) {
					$output_data_ref->{'data_table'}->[$row][$column] = 
						$input_data2_ref->{'data_table'}->[$row2][$request_number];
				}
				else {
					$output_data_ref->{'data_table'}->[$row][$column] = '.';
				}
			}
			
			# copy the metadata
			copy_metadata($input_data2_ref, $request_number);
		} 
		
		else {
			die " unrecognized  symbol '$request' in request! nothing done!\n";
		}
	}
	
}



### Add datasets from one data file to the output data structure
sub add_datasets {
	my $data_ref = shift;
	
	# Check the input data
	my $check = check_data_tables($output_data_ref, $data_ref);
	
	
	# Determine if we need to do lookup
	my $lookup_i;
	if ($check or $lookup) {
		# we need to merge by lookup values
		
		# check whether we have a lookup column defined
		unless ($lookup_name) {
			# we have to do this by treating as two new datasets
			merge_two_datasets_by_lookup($output_data_ref, $data_ref);
			return;
		}
		
		# proceed with single lookup
		$lookup_i = request_lookup_indices($data_ref);
		
		# index the new data reference
			# we'll be putting the lookup values into an index hash
			# where the lookup value is the key and the row number is the value
		$data_ref->{'index'} = {};
		for (my $row = 1; $row <= $data_ref->{'last_row'}; $row++) {
			my $key = $data_ref->{'data_table'}->[$row][$lookup_i];
			if (exists $data_ref->{'index'}{$key}) {
				# value is not unique
				warn " lookup value '$key' in file " . 
					$data_ref->{'filename'} . 
					", row $row is a duplicate!\n" . 
					" Using the first occurence value\n";
			}
			else {
				# value is ok
				$data_ref->{'index'}{$key} = $row;
			}
		}
		
	}
	
	
	# Merge the new datasets with current output
	# determine the new order
	my @order;
	if ($automatic) {
		# automatic selection
		@order = automatically_determine_order($data_ref);
		
	}
	else {
		# manual selection from user
		@order = request_new_order($data_ref);
	}
	
	# assign the datasets to the output data in requested order
	foreach my $request (@order) {
		
		# determine the current column index we're working with
		my $column = $output_data_ref->{'number_columns'};
		
		# add the dataset from the appropriate input file
		if ($request =~ /^\d+$/) {
			
			# print dataset name in automatic mode
			if ($automatic) {
				print "  Merging column " . 
					$data_ref->{$request}{'name'} . "\n";
			}
			
			# copy the dataset by lookup or blindly
			if (defined $lookup_i) {
				# merging by lookup
				
				# copy the header
				$output_data_ref->{'data_table'}->[0][$column] = 
					$data_ref->{'data_table'}->[0][$request];
				
				# copy the dataset
				for my $row (1 .. $output_data_ref->{'last_row'}) {
					# identify the appropriate row in file2 by lookup value
					my $lookup = 
						$output_data_ref->{'data_table'}->[$row][$lookup_i];
					my $row2 = $data_ref->{'index'}{$lookup} || undef;
					
					# copy the appropriate value
					if (defined $row2) {
						$output_data_ref->{'data_table'}->[$row][$column] = 
							$data_ref->{'data_table'}->[$row2][$request];
					}
					else {
						$output_data_ref->{'data_table'}->[$row][$column] = '.';
					}
				}
			}
			else {
				# merging blindly
				
				for my $row (0 .. $data_ref->{'last_row'}) {
					$output_data_ref->{'data_table'}->[$row][$column] = 
						$data_ref->{'data_table'}->[$row][$request];
				}
			}
			
			# copy the metadata
			copy_metadata($data_ref, $request);
		} 
		else {
			die " unrecognized  symbol '$request' in request! nothing done!\n";
		}
	}
}


### Request the new order
sub request_new_order {
	my @data_refs = @_;
	# Either one or two data file references may be passed
	# The dataset names will be printed listed by their index. Numbers will be
	# used for the first data file, letters for the second, if present.
	
	# Print the column headers
	if (scalar @data_refs == 1) {
		# Only one data file
		# we'll be using only numbers to list the datasets
		print_datasets($data_refs[0], 'number');
	}
	elsif (scalar @data_refs == 2) {
		# Two data files
		# use numbers for the first one
		print_datasets($data_refs[0], 'number');
		# use letters for the second one
		print_datasets($data_refs[1], 'letter');
		
	}
		
	# Request response from user
	print " Enter the columns' indices in the desired final order.\n" .
		" Separate indices with commas or specify a range (start - stop).\n   ";
	my $response = <STDIN>;
	chomp $response;
	$response =~ s/\s+//g;
	
	# Parse response
	my @order = parse_list($response); 
	
	# done
	print " using order: ", join(", ", @order), "\n";
	return @order;
}


### Request lookup indices for two datasets
sub request_lookup_indices {
	my ($data1, $data2) = @_;
	# One or two data references may be passed
	# We'll ask for the lookup value indices for each
	
	
	# Determine the lookup index values
	my $index1;
	my $index2;
	
	# Determine if we need one or two
	if (!defined $data2 and $lookup_name) {
		# we already have a lookup name
		# just need to find same column for one dataset
		
		my $index1 = find_column_index($data1, "^$lookup_name\$");
		unless (defined $index1) {
			die " Cannot find lookup column with name '$lookup_name'" . 
				" in file " . $data1->{'filename'} . "\n";
		}
		
		# print the found column name
		print "  using column $index1 (", 
			$data1->{$index1}{'name'}, 
			") as lookup index for file ", 
			$data1->{'filename'}, "\n";
		return $index1;
	}
	
	# First try some known column identifiers we could use automatically
	foreach my $name (qw(name id transcript gene)) {
		
		# identify possibilities
		my $possible_index1 = find_column_index($data1, "^$name\$");
		my $possible_index2 = find_column_index($data2, "^$name\$");
		
		# check if something was found
		if (defined $possible_index1 and defined $possible_index2) {
			
			# assign
			$index1 = $possible_index1;
			$index2 = $possible_index2;
			$lookup_name = $name; # for future lookups
			
			# report
			print "  using column $index1 (", 
				$data1->{$index1}{'name'}, 
				") as lookup index for file ", 
				$data1->{'filename'}, "\n";
			print "  using column $index2 (", 
				$data2->{$index2}{'name'}, 
				") as lookup index for file ", 
				$data2->{'filename'}, "\n";
		}
	}
	
	unless (defined $index1 and defined $index2) {
		# Automatic identification didn't work
		# must bother the user for help
		if ($automatic) {
			die " Unable to identify appropriate lookup columns automatically!\n" .
				" Please execute interactively to identify lookup columns\n";
		}
		
		print " Unable to identify appropriate lookup columns automatically!\n";
		
		# Print the index headers
		# use numbers for the first one
		print_datasets($data1, 'number');
		# use letters for the second one
		print_datasets($data2, 'letter');
	
		# Request first index responses from user
		print " Enter the unique identifier index for lookup in the first file   ";
		$index1 = <STDIN>;
		chomp $index1;
		unless ($index1 =~ /^\d+$/ and exists $data1->{$index1}) {
			# check that it's valid
			die " unknown index value!\n";
		}
		
		# Request second index responses from user
		print " Enter the unique identifier index for lookup in the second file   ";
		$index2 = <STDIN>;
		chomp $index2;
		unless ($index2 =~ /^[a-z]+$/i) { 
			# check that it's valid
			die " unknown index value!\n";
		}
		$index2 = $number_of{$index2}; # convert to a number
		unless (exists $data1->{$index2}) {
			# check that it's valid
			die " unknown index value!\n";
		}
		
	}
	
	# done
	return $index1, $index2;
}




### Automatically determine the order
sub automatically_determine_order {
	
	# get the two data structures
	# we could have two, or only one
	my ($data1, $data2, $number);
	if (scalar @_ == 1) {
		# comparing new data with the output data we are building
		$data1 = $output_data_ref;
		$data2 = shift @_;
		$number = 1;
	}
	else {
		# two new datasets 
		$data1 = shift @_;
		$data2 = shift @_;
		$number = 2;
	}
	
	# generate quick hash of names to exclude from data1
	my %exclude;
	for (my $i = 0; $i < $data1->{'number_columns'}; $i++) {
		$exclude{ $data1->{$i}{'name'} } = 1;
	}
	
	# generate the order
	my @order;
	
	# add the first dataset if provided
	if ($number == 2) {
		# we will automatically take all of the columns from the first data
		for (my $i = 0; $i < $data1->{'number_columns'}; $i++) {
			push @order, $i;
		}
	}
	
	# automatically identify those columns in second data not present in
	# the first one to take
	for (my $i = 0; $i < $data2->{'number_columns'}; $i++) {
		if ($data2->{$i}{'name'} eq 'Score') {
			# name is score, take it
			if ($number == 1) {
				# one data file only, push number
				push @order, $i;
			}
			else {
				# two data files, push letter
				push @order, $letter_of{$i};
			}
		}
		
		elsif (exists $exclude{ $data2->{$i}{'name'} } ) {
			# non-unique column, skip
			next;
		}
		
		else {
			# must be a unique column, take it
			if ($number == 1) {
				# one dataset only, push number
				push @order, $i;
			}
			else {
				# two datasets, push letter
				push @order, $letter_of{$i};
			}
		}
	}
	
	# done
	return @order;
}


### Parse a string into a list
sub parse_list {
	my $string = shift;
	
	# Parse the string
	my @list;
	foreach my $item (split /,/, $string) {
		
		# comma-delimited list of items
		# may contain ranges of incremental items, check for those
		
		if ($item =~ /\-/) {
			# item is a range, specified as 'start - stop'
			
			if ($item =~ /^(\d+)\-(\d+)$/) {
				# range contains numbers, so must be from file1
				for (my $i = $1; $i <= $2; $i++) {
					# we will loop through from specified start to stop
					# add each number to the order
					push @list, $i;
				}
			} 
			
			elsif ($item =~ /^([a-z]+)\-([a-z]+)$/) {
				# range does not contain numbers, so must be from file2
				for (my $i = $number_of{$1}; $i <= $number_of{$2}; $i++) {
					# we will loop through from specified start to stop
					# converting each letter to the corresponding number for 
					# use in the for loop
					# then convert the number in $i back to a letter to add
					# to the order 
					push @list, $letter_of{$i};
				}
			} 
			
			else {
				# unrecognizable
				die " unrecognizable range order!\n";
			}
		} 
		
		else {
			# item is not a range but a single number/letter
			# add it to the order
			push @list, $item;
		}
	
	}
	
	return @list;
}



### Print the dataset indices and names for the passed data file 
sub print_datasets {
	
	# pass the data structure reference and the index type ('number' or 'letter')
	my ($data_ref, $index_type) = @_;
	
	# array to be used in the default order
	my @order;
	
	# print the dataset names for this datafile
	print " These are the headers in file '" . $data_ref->{'filename'} . "'\n";
	foreach (my $i = 0; $i < $data_ref->{'number_columns'}; $i++) {
		# walk through each dataset (column) in file
		
		# print the dataset name and it's index
		# and record the index in the order array for the default order
		if ($index_type eq 'number') {
			print '  ', $i, "\t", $data_ref->{$i}{'name'}, "\n";
			push @order, $i;
		}
		else {
			# use letters instead of numbers as the index key
			my $letter = $letter_of{$i};
			print '  ', $letter, "\t", $data_ref->{$i}{'name'}, "\n";
			push @order, $letter;
		}
	}
	
	return @order;
}




### Check the data tables for similarity
sub check_data_tables {
	my ($input_data1_ref, $input_data2_ref) = @_;
	
	# Check the feature types
	if ( 
		$input_data1_ref->{'feature'} and 
		$input_data2_ref->{'feature'} and
		$input_data1_ref->{'feature'} ne $input_data2_ref->{'feature'} 
	) {
		# Each file has feature type defined, but they're not the same
		# probably really don't want to combine these
		print " WARNING! The metadata feature types for both files don't match!!!\n";
		print "   Continuing anyway...\n";
	}
	
	# check line numbers
	my $status = 0;
	if ( $input_data1_ref->{'last_row'} != $input_data2_ref->{'last_row'} ) {
		# the number of rows in each data table don't equal
		# we will need to do this by lookup
		$status = 1;
	}
	return $status;
}



### Initialize the output data structure
sub initialize_output_data_structure {
	my $data_ref = shift;
	
	# generate brand new output data structure
	my $output_data = generate_tim_data_structure(
		$data_ref->{'feature'}, # re-use the same feature as file1
	);
	
	# we'll re-use the values from file1 
	$output_data->{'program'}   = $data_ref->{'program'};
	$output_data->{'db'}        = $data_ref->{'db'};
	$output_data->{'filename'}  = $data_ref->{'filename'}; # borrow file name
	$output_data->{'last_row'}  = $data_ref->{'last_row'}; # should be identical

	return $output_data;
}



### Copy the dataset metadata
sub copy_metadata {
	
	# collect arguments
	my ($data_ref, $dataset_index) = @_;
	
	# determine current dataset index to copy into
	my $current_index = $output_data_ref->{'number_columns'};
	
	# copy the metadata
	$output_data_ref->{$current_index} = { %{ $data_ref->{$dataset_index} } };
	
	# check if should rename the dataset
	if ($automatic and $data_ref->{$dataset_index}{'name'} eq 'Score') {
		# only in automatic mode and the dataset is a generic Score
		# no opportunity to rename interactively
		# use file basename appended with Score
		$output_data_ref->{$current_index}{'name'} = 
			$data_ref->{'basename'} . '_Score';
		$output_data_ref->{'data_table'}->[0][$current_index] = 
			$output_data_ref->{$current_index}{'name'};
	}
	
	# set the original file name
	$output_data_ref->{$current_index}{'original_file'} = 
		$data_ref->{'filename'};
	
	# reset the index
	$output_data_ref->{$current_index}{'index'} = $current_index;
	
	# increment the number of columns
	$output_data_ref->{'number_columns'} += 1;
}
	


### Re-name the dataset names
sub rename_dataset_names {
	print " For each header, type a new name or push enter to accept the current name\n";
	
	for (my $i = 0; $i < $output_data_ref->{'number_columns'}; $i++) {
		# walk through the list of columns
		
		# print the current name and it's originating file name
		print '  ', $output_data_ref->{$i}{'original_file'}, ': ', 
			$output_data_ref->{$i}{'name'}, ' ->  ';
		
		# request user input
		my $new_name = <STDIN>;
		chomp $new_name;
		
		# rename as requested
		if ($new_name) {
			# user entered a new name
			# assign new names in both metadata and column header
			$output_data_ref->{$i}{'name'} = $new_name;
			$output_data_ref->{'data_table'}->[0][$i] = $new_name;
		}
	}
}





### Generate the lookup hash for numbers to letters
sub get_letters {
	my %hash = (
		0 => 'a',
		1 => 'b',
		2 => 'c', 
		3 => 'd',
		4 => 'e',
		5 => 'f',
		6 => 'g',
		7 => 'h',
		8 => 'i',
		9 => 'j',
		10 => 'k',
		11 => 'l',
		12 => 'm',
		13 => 'n',
		14 => 'o',
		15 => 'p',
		16 => 'q',
		17 => 'r',
		18 => 's',
		19 => 't',
		20 => 'u',
		21 => 'v',
		22 => 'w',
		23 => 'x',
		24 => 'y',
		25 => 'z',
		26 => 'aa',
		27 => 'bb',
		28 => 'cc',
		29 => 'dd',
		30 => 'ee',
		31 => 'ff',
		32 => 'gg',
		33 => 'hh',
		34 => 'ii',
		35 => 'jj',
		36 => 'kk',
		37 => 'll',
		38 => 'mm',
		39 => 'nn',
		40 => 'oo',
		41 => 'pp',
		42 => 'qq',
		43 => 'rr',
		44 => 'ss',
		45 => 'tt',
		46 => 'uu',
		47 => 'vv',
		48 => 'ww',
		49 => 'xx',
		50 => 'yy',
		51 => 'zz',
		52 => 'aaa',
		53 => 'bbb',
		54 => 'ccc',
		55 => 'ddd',
		56 => 'eee',
		57 => 'fff',
		58 => 'ggg',
		59 => 'hhh',
		60 => 'iii',
		61 => 'jjj',
		62 => 'kkk',
		63 => 'lll',
		64 => 'mmm',
		65 => 'nnn',
		66 => 'ooo',
		67 => 'ppp',
		68 => 'qqq',
		69 => 'rrr',
		70 => 'sss',
		71 => 'ttt',
		72 => 'uuu',
		73 => 'vvv',
		74 => 'www',
		75 => 'xxx',
		76 => 'yyy',
		77 => 'zzz',
	);
	return %hash;
}


### Generate the lookup hash for letters to numbers
sub get_numbers {
	my %hash = (
		'a' => 0,
		'b' => 1,
		'c' => 2, 
		'd' => 3,
		'e' => 4,
		'f' => 5,
		'g' => 6,
		'h' => 7,
		'i' => 8,
		'j' => 9,
		'k' => 10,
		'l' => 11,
		'm' => 12,
		'n' => 13,
		'o' => 14,
		'p' => 15,
		'q' => 16,
		'r' => 17,
		's' => 18,
		't' => 19,
		'u' => 20,
		'v' => 21,
		'w' => 22,
		'x' => 23,
		'y' => 24,
		'z' => 25,
		'aa' => 26,
		'bb' => 27,
		'cc' => 28,
		'dd' => 29,
		'ee' => 30,
		'ff' => 31,
		'gg' => 32,
		'hh' => 33,
		'ii' => 34,
		'jj' => 35,
		'kk' => 36,
		'll' => 37,
		'mm' => 38,
		'nn' => 39,
		'oo' => 40,
		'pp' => 41,
		'qq' => 42,
		'rr' => 43,
		'ss' => 44,
		'tt' => 45,
		'uu' => 46,
		'vv' => 47,
		'ww' => 48,
		'xx' => 49,
		'yy' => 50,
		'zz' => 51,
		'aaa' => 52,
		'bbb' => 53,
		'ccc' => 54,
		'ddd' => 55,
		'eee' => 56,
		'fff' => 57,
		'ggg' => 58,
		'hhh' => 59,
		'iii' => 60,
		'jjj' => 61,
		'kkk' => 62,
		'lll' => 63,
		'mmm' => 64,
		'nnn' => 65,
		'ooo' => 66,
		'ppp' => 67,
		'qqq' => 68,
		'rrr' => 69,
		'sss' => 70,
		'ttt' => 71,
		'uuu' => 72,
		'vvv' => 73,
		'www' => 74,
		'xxx' => 75,
		'yyy' => 76,
		'zzz' => 77,
	);
	return %hash;
}

__END__

=head1 NAME

merge_datasets.pl

=head1 SYNOPSIS

merge_datasets.pl [--options...] <file1> <file2> ...
  
  Options:
  --lookup
  --auto
  --out <filename> 
  --(no)gz
  --version
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --lookup

Force the program to merge data by using lookup values in each file. 
Enable this option if the data rows are not sorted identically.
This will be done automatically if the number of data rows are not 
equal between the files.

=item --auto

Execute in automatic mode, where all columns from the first data file 
are retained, and only uniquely named or Score columns from subsequent 
data files are merged. If lookup is enabled, then an appropriate 
lookup column will be chosen.

=item --out <filename>

Specify the output filename. By default it uses the first file name.
Required in automatic mode.

=item --(no)gz

Specify whether (or not) the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will merge two or more tab-delimited data files into one file. 
Datasets or columns from each file are merged together into an output file. 
Columns are appended to the end (after the rightmost column). 

By default, the program is run in an interactive mode allowing the columns 
to be chosen from a presented list. Alternatively, it may be run in 
automatic mode, where uniquely named datasets from subsequent files are 
appended to the first file. Score columns from specific formatted files 
(BED, BedGraph, GFF) are also automatically taken. 

The program blindly assumes that rows (features) are equivalent in all of the 
datasets if there are an equal number of data rows. However, if there are an 
unequal number of data rows, or the user forces by using the --lookup option, 
then dataset values are looked up first using specified lookup values before 
merging (compare with Excel VLOOKUP function). In this case, the dataset 
lookup indices from each file are identified automatically. Potential 
lookup columns include 'Name', 'ID', 'Transcript', or 'Gene'. Failing that, 
they are chosen interactively.

When a lookup is performed, the first index in the order determines which 
file is dominant, meaning that all rows from that file are included, and only 
the rows that match by the lookup value are included from the second file. Null 
values are recorded when no match is found in the second file.

After merging in interactive mode, an opportunity for interactively 
renaming the dataset names is presented.

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
