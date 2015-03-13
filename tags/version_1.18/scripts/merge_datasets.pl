#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure
	find_column_index
);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	write_tim_data_file
);
my $VERSION = '1.17';

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
	$manual,
	$user_lookup_name,
	$use_coordinate,
	$outfile,
	$gz,
	$help,
	$print_version,
);
my @order_requests;

# Command line options
GetOptions( 
	'lookup!'   => \$lookup, # force merging by value lookup
	'auto!'     => \$automatic, # select columns automatically
	'manual!'   => \$manual, # always run interactively
	'lookupname|lun=s' => \$user_lookup_name, # alternate lookup column name
	'coordinate!' => \$use_coordinate, # use coordinates for lookup
	'index=s'   => \@order_requests, # determine order in advance
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

# user provided index list
my $user_list = @order_requests ? 1 : 0;

# automatic mode
if ($automatic and $manual) {
	die " AAGH! I can't do both manual and automatic at the same time!\n";
}
if ($automatic or $user_list) {
	unless ($outfile) {
		die " Must provide an output file name in automatic mode!\n";
	}
}

# Get the alphabet-number conversion hashes
my ($letter_of, $number_of) = get_number_letters(); # convert numbers into letters

# Set up output
my $output_data_ref; # the reference scalar for the output data structure

# name of lookup column to be used for all files
my $lookup_name;
my $output_lookup_i;
if ($use_coordinate) {
	$lookup_name = 'Coordinate';
}


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
	undef $input_data1_ref;
	undef $input_data2_ref;
	
	# merge the subsequent files
	foreach (@ARGV) {
		my $input_data_ref = read_file($_);
		add_datasets($input_data_ref);
	}
}

# clean up coordinate column
if (
	$lookup_name =~ /^coordinate$/i and 
	$output_data_ref->{0}{'name'} =~ /^coordinate$/i
) {
	# delete the coordinate metadata
	for my $i (1 .. $output_data_ref->{'number_columns'}-1) {
		$output_data_ref->{$i-1} = $output_data_ref->{$i};
	}
	$output_data_ref->{'number_columns'}--;
	
	# delete the table column
	for my $row (0 .. $output_data_ref->{'last_row'}) {
		shift @{ $output_data_ref->{'data_table'}->[$row] };
	}
} 



### Offer to re-name the datasets
unless ($automatic or $user_list) {
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
my $file_written = write_tim_data_file(
	'data'      => $output_data_ref,
	'filename'  => $outfile,
	'gz'        => $gz,
);
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
	
	# delete any pre-existing original_file metadata
	# we'll be writing new
	for (my $i = 0; $i < $file_data_ref->{'number_columns'}; $i++) {
		if (exists $file_data_ref->{$i}{'original_file'}) {
			delete $file_data_ref->{$i}{'original_file'};
		}
	}
	
	# add coordinates if necessary
	if ($file_data_ref->{'bed'} or $file_data_ref->{'gff'} or $use_coordinate) {
		
		# refuse if user does not want to use coordinates
		if (defined $use_coordinate and $use_coordinate == 0) {
			return $file_data_ref;
		}
		
		# identify coordinate columns
		my $coord_i = $file_data_ref->{'number_columns'};
		my $chr_i   = find_column_index($file_data_ref, '^chr|seq|ref|ref.?seq');
		my $start_i = find_column_index($file_data_ref, '^start|position');
		my $stop_i  = find_column_index($file_data_ref, '^stop|end');
		unless (defined $chr_i and defined $start_i) {
			# cannot add coordinate column, do without ?
			warn " cannot generate coordinates for file\n";
			return $file_data_ref;
		}
		
		# add new metadata
		$file_data_ref->{$coord_i} = {
			'name'      => 'Coordinate',
			'index'     => $coord_i,
		};
		$file_data_ref->{'data_table'}->[0][$coord_i] = 'Coordinate';
		$file_data_ref->{'number_columns'}++;
		
		# generate coordinates
		if (defined $stop_i) {
			# merge chromosome:start-stop
			for my $row (1 .. $file_data_ref->{'last_row'}) {
				$file_data_ref->{'data_table'}->[$row][$coord_i] = join("", 
					$file_data_ref->{'data_table'}->[$row][$chr_i], ':', 
					$file_data_ref->{'data_table'}->[$row][$start_i], '-',
					$file_data_ref->{'data_table'}->[$row][$stop_i]
				);
			}
		}
		else {
			# merge chromosome:start
			for my $row (1 .. $file_data_ref->{'last_row'}) {
				$file_data_ref->{'data_table'}->[$row][$coord_i] = join("", 
					$file_data_ref->{'data_table'}->[$row][$chr_i], ':', 
					$file_data_ref->{'data_table'}->[$row][$start_i]
				);
			}
		}
	}
	
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
			my $number = $number_of->{$request};
			
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
	
	# add coordinate column to output if necessary
	if ($lookup_name =~ /^coordinate$/i) {
		
		# check if we taking from file 1 or 2
		if ($order[0] =~ /[a-z]+/i) {
			unshift @order, $lookup_i2;
		}
		else {
			unshift @order, $lookup_i1;
		} 
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
				$_ = $letter_of->{$_};
			}
			else {
				$_ = $number_of->{$_};
			}
		} @order;
	}
	
	# Index the second dataset
		# we'll be putting the lookup values into an index hash
		# where the lookup value is the key and the row number is the value
	$input_data2_ref->{'index'} = {};
	my $index_warning = 0;
	for (my $row = 1; $row <= $input_data2_ref->{'last_row'}; $row++) {
		my $key = $input_data2_ref->{'data_table'}->[$row][$lookup_i2];
		if (exists $input_data2_ref->{'index'}{$key}) {
			# value is not unique
			$index_warning++;
		}
		else {
			# value is ok
			$input_data2_ref->{'index'}{$key} = $row;
		}
	}
	if ($index_warning) {
		warn " Warning: $index_warning rows had two or more duplicate lookup values\n" . 
			"  for column $lookup_i2 in file " . $input_data2_ref->{'filename'} . 
			"\n  Only the first occurence was used\n";
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
			my $request_number = $number_of->{$request};
			
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
		
		# make sure we remember the lookup_column_index in the new output
		# if it is necessary
		if ($lookup_name and not defined $output_lookup_i) {
			# we have a lookup name, but the index for this lookup in
			# the output data structure is not yet defined
			
			# check if the current output index is it
			if ($output_data_ref->{$column}{'name'} =~ /\A $lookup_name \Z/xi) {
				# this current column matches the lookup name
				# so we will use it
				$output_lookup_i = $column;
			}
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
		if ($manual or not $lookup_name) {
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
				
				# check that we have the output lookup index
				unless (defined $output_lookup_i) {
					die " The lookup index for column '$lookup_name' is " .
						"not defined in the output!\n" . 
						" Please ensure you include column '$lookup_name' in" .
						" the output file.\n";
				}
				
				# copy the header
				$output_data_ref->{'data_table'}->[0][$column] = 
					$data_ref->{'data_table'}->[0][$request];
				
				# copy the dataset
				for my $row (1 .. $output_data_ref->{'last_row'}) {
					# identify the appropriate row in file2 by lookup value
					my $lookup = 
						$output_data_ref->{'data_table'}->[$row][$output_lookup_i];
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
	my $list;
	
	# Check if user supplied the order by command line
	if (@order_requests) {
		$list = shift @order_requests;
	}
	
	# else request interactively
	else {
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
		$list = <STDIN>;
		chomp $list;
		$list =~ s/\s+//g;
	}
	
	# Parse response
	my @order = parse_list($list); 
	
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
	
	
	# Automatic methods of determining the lookup index
	unless ($manual) {
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
		# don't forget to add the user-requested lookup name
		my @name_list = qw(coordinate name id transcript gene);
		if ($user_lookup_name) {
			unshift @name_list, $user_lookup_name;
		}
		foreach my $name (@name_list) {
			
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
				
				# don't go through remaining list
				last;
			}
		}
	}
	
	unless (defined $index1 and defined $index2) {
		# Automatic identification didn't work
		# must bother the user for help
		if ($automatic) {
			die " Unable to identify appropriate lookup columns automatically!\n" .
				" Please execute interactively to identify lookup columns\n";
		}
		
		print " Unable to identify appropriate lookup columns automatically!\n"
			unless $manual;
		
		# Print the index headers
		# use numbers for the first one
		print_datasets($data1, 'number', 1);
		# use letters for the second one
		print_datasets($data2, 'letter', 1);
	
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
		$index2 = $number_of->{$index2}; # convert to a number
		unless (exists $data2->{$index2}) {
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
				push @order, $letter_of->{$i};
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
				push @order, $letter_of->{$i};
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
				for (my $i = $number_of->{$1}; $i <= $number_of->{$2}; $i++) {
					# we will loop through from specified start to stop
					# converting each letter to the corresponding number for 
					# use in the for loop
					# then convert the number in $i back to a letter to add
					# to the order 
					push @list, $letter_of->{$i};
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
	my ($data_ref, $index_type, $include_coordinate) = @_;
	
	# array to be used in the default order
	my @order;
	
	# print the dataset names for this datafile
	print " These are the headers in file '" . $data_ref->{'filename'} . "'\n";
	
	# print Numbers
	if ($index_type eq 'number') {
		foreach (my $i = 0; $i < $data_ref->{'number_columns'}; $i++) {
			# skip the coordinate
			next if ($data_ref->{$i}{'name'} eq 'Coordinate' and not $include_coordinate);
			
			# print the dataset name and it's index
			# and record the index in the order array for the default order
			print '  ', $i, "\t", $data_ref->{$i}{'name'}, "\n";
			push @order, $i;
		}
	}
	
	# print letters
	else {
		foreach (my $i = 0; $i < $data_ref->{'number_columns'}; $i++) {
			# skip the coordinate
			next if ($data_ref->{$i}{'name'} eq 'Coordinate' and not $include_coordinate);
			
			# use letters instead of numbers as the index key
			my $letter = $letter_of->{$i};
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
	$output_data->{'last_row'}  = $data_ref->{'last_row'}; # should be identical
	$output_data->{'headers'}   = $data_ref->{'headers'};
	
	# assign the filename
	if ($outfile) {
		# use the new output file name
		$output_data->{'filename'}  = $outfile;
	}
	else {
		# borrow the first file name
		$output_data->{'filename'}  = $data_ref->{'filename'}; # borrow file name
	}
	
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
	unless (exists $output_data_ref->{$current_index}{'original_file'}) {
		$output_data_ref->{$current_index}{'original_file'} = 
			$data_ref->{'filename'};
	}
	
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



### Generate lookup hashes for converting numbers to letters and back
sub get_number_letters {
	
	# hashes
	my %n2l; # numbers to letters
	my %l2n; # letters to numbers
	my %lookup = (
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
	);
	
	# fill up the lookup hashes
	my $first;
	for (my $i = 0; $i < 702; $i++) {
		# this gives range from a to zz
		
		# check the first letter
		if ( $i % 26 == 0) {
			$first = $lookup{ ($i / 26) - 1 };
		}
		
		# generate the letter
		# two letters [null..z][a..z]
		my $letter = $first . $lookup{ $i % 26 };
		
		# store in hashes
		$n2l{ $i } = $letter;
		$l2n{ $letter } = $i;
	}
	
	return (\%n2l, \%l2n);
}


__END__

=head1 NAME

merge_datasets.pl

A program to merge two or more data files by appending columns.

=head1 SYNOPSIS

merge_datasets.pl [--options...] <file1> <file2> ...
  
  Options:
  --lookup
  --auto
  --manual
  --index <number,letter,range>
  --lookupname | --lun <text>
  --coordinate
  --out <filename> 
  --gz
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

=item --manual

Execute in full manual mode, where columns for lookup and merge are 
chosen interactively from lists.

=item --index <number,letter,range>

For advanced users, provide the list of indices to merge from both 
datasets. The format is the same as in the interactive mode. Column 
indices from the first file are represented by numbers, the second 
file by letters. Values are comma-delimited, and ranges may be 
provided as "start-stop". To indicate indices from subsequent files, 
provide separate --index options for each subsequent file. The default 
is to run the program interactively.

=item --lookupname <text>

=item --lun <text>

Provide an alternate column name to identify the columns automatically 
in the input files containing the lookup values when performing the 
lookup. Each file should have the same lookup column name. Default 
values include 'Name', 'ID', 'Transcript', or 'Gene'.

=item --coordinate

Force using genomic coordinates when performing a lookup. This effectively 
merges the chromosome:start-stop coordinates into a single string for 
lookup matching. The Coordinates column is temporary. This is automatically 
enabled when working with BED or GFF files that may or may not have unique 
names but will have unique coordinates. This may be disabled, for example 
to force using a BED feature name as the lookup value, by specifying 
--nocoordinate.

=item --out <filename>

Specify the output filename. By default it uses the first file name.
Required in automatic mode.

=item --gz

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

The columns may be specified using one of three methods. By default, 
the program is run in an interactive mode allowing the columns to be 
chosen from a presented list for each input file provided. Second, the 
columns may be specified manually on the command line using one or 
more --index options. Third, the program may be executed in full 
automatic mode, where uniquely named datasets from subsequent files 
are automatically appended to the first file. Score columns from 
specific formatted files (BED, BedGraph, GFF) are also automatically 
taken. 

The program blindly assumes that rows (features) are equivalent and in the 
same order for all of the datasets. However, if there are an 
unequal number of data rows, or the user forces by using the --lookup option, 
then dataset values are looked up first using specified lookup values before 
merging (compare with Excel VLOOKUP function). In this case, the dataset 
lookup indices from each file are identified automatically. Potential 
lookup columns include 'Name', 'ID', 'Transcript', 'Gene', or any user 
provided name (using the --lookupname option). Files with genomic coordinates, 
including GFF and BED, may use a temporary coordinate string derived from the 
features coordinates. If a lookup column can not be identified automatically, 
then they are chosen interactively.

When a lookup is performed, the first index in the order determines which 
file is dominant, meaning that all rows from that file are included, and only 
the rows that match by the lookup value are included from the second file. Null 
values are recorded when no match is found in the second file.

After merging in interactive mode only, an opportunity for interactively 
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
