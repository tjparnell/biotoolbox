#!/usr/bin/perl

# documentation at end of file

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use IO::Prompt::Tiny qw(prompt);
use Bio::ToolBox::Data;

our $VERSION = '2.02';

print "\n A progam to merge datasets from two files\n";

# Print help
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
my (
	$lookup,           $automatic,      $manual,
	$user_lookup_name, $use_coordinate, $outfile,
	$gz,               $help,           $print_version,
);
my @order_requests;

# Command line options
GetOptions(
	'l|lookup!'          => \$lookup,              # force merging by value lookup
	'a|auto!'            => \$automatic,           # select columns automatically
	'm|manual!'          => \$manual,              # always run interactively
	'n|lookupname|lun=s' => \$user_lookup_name,    # alternate lookup column name
	'c|coordinate!'      => \$use_coordinate,      # use coordinates for lookup
	'x|index=s'          => \@order_requests,      # determine order in advance
	'o|out=s'            => \$outfile,             # name of output file
	'z|gz!'              => \$gz,                  # compress output
	'h|help'             => \$help,                # request help
	'v|version'          => \$print_version,       # print the version
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
	print " Biotoolbox script merge_datasets.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

### Check for requirements and set defaults

# input files
unless ( scalar @ARGV >= 2 ) {
	print STDERR " FATAL: At least two input files are required!\n";
	exit 1;
}

# user provided index list
my $user_list = @order_requests ? 1 : 0;

# automatic mode
if ( $automatic and $manual ) {
	print STDERR " FATAL: Cannot specify both manual and automatic at the same time!\n";
	exit 1;
}
if ( $automatic or $user_list ) {
	unless ($outfile) {
		print STDERR " FATAL: Must provide an output file name in automatic mode!\n";
		exit 1;
	}
}

# Get the alphabet-number conversion hashes
my ( $letter_of, $number_of ) = get_number_letters();    # convert numbers into letters

# Set up output
my $output_data;    # the reference scalar for the output data structure

# name of lookup column to be used for all files
my $lookup_name;
my $output_lookup_i;
if ($use_coordinate) {
	$lookup_name = 'MergeDatasetCoordinate';
}

### Process and merge the files

# exactly two filenames are provided
if ( scalar @ARGV == 2 ) {
	my $input_data1 = read_file( shift @ARGV );
	my $input_data2 = read_file( shift @ARGV );
	merge_two_datasets( $input_data1, $input_data2 );
}

# more than two filenames
else {
	# merge the first two
	my $input_data1 = read_file( shift @ARGV );
	my $input_data2 = read_file( shift @ARGV );
	merge_two_datasets( $input_data1, $input_data2 );

	# undefine these data sources, hopefully to release memory
	undef $input_data1;
	undef $input_data2;

	# merge the subsequent files
	while (@ARGV) {
		my $input_data = read_file( shift @ARGV );
		add_datasets($input_data);
	}
}

# clean up coordinate column if it was generated
if ( my $c = $output_data->find_column('MergeDatasetCoordinate') ) {
	$output_data->delete_column($c);
}

### Offer to re-name the datasets
unless ( $automatic or $user_list ) {
	my $p        = ' Do you wish to rename the headers? y/n (default n):  ';
	my $response = prompt( $p, 'n' );

	if ( lc $response eq 'y' ) {
		rename_dataset_names();
	}
}

# clean up metadata we added but shouldn't keep
# this was put in here for remembering original file names for renaming purposes
for my $i ( 1 .. $output_data->number_columns ) {
	$output_data->delete_metadata( $i, 'original_file' );
}

### Print the output
# Request the output file name
# default is to simply overwrite file1
unless ($outfile) {
	my $p = sprintf " Enter the output file name [%s]: ", $output_data->{'filename'};
	$outfile = prompt( $p, $output_data->{'filename'} );
}

# write the file
my $file_written = $output_data->write_file(
	'filename' => $outfile,
	'gz'       => $gz,
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
	my $Data = Bio::ToolBox::Data->new( file => $filename );
	unless ($Data) {
		die " file '$filename' not loaded!";
	}

	# print the results
	printf "\n Loaded '$filename' with %s rows and %s columns\n",
		$Data->last_row, $Data->number_columns;

	# delete any pre-existing original_file metadata
	# we'll be writing new
	for my $i ( 1 .. $Data->number_columns ) {
		$Data->delete_metadata( $i, 'original_file' );
	}

	# add coordinates if necessary
	if ( $Data->bed or $Data->gff or $use_coordinate ) {

		# refuse if user does not want to use coordinates
		if ( defined $use_coordinate and $use_coordinate == 0 ) {
			return $Data;
		}

		# identify coordinate columns
		if ( ( defined $Data->chromo_column and defined $Data->start_column )
			or defined $Data->coord_column )
		{
			my @columns;
			push @columns, $Data->chromo_column;
			push @columns, $Data->start_column;
			push @columns, $Data->stop_column;
			push @columns, $Data->coord_column;
			printf "  generating coordinate lookup string from columns %s\n",
				join( ', ', grep {defined} @columns );
		}
		else {
			print
"  cannot generate coordinates for file!! No chromosome, start, or coordinate columns\n";
			return $Data;
		}

		# add new column
		my $coord_i = $Data->add_column('MergeDatasetCoordinate');

		# generate two coordinate strings for lookup
		# one should be 1-based, the other 0-based, hopefully we will match on one
		$Data->iterate(
			sub {
				my $row = shift;
				my $c   = $row->seq_id;
				my $e   = $row->end;
				my $v;
				if ( defined $e ) {
					$v = sprintf "%s:%d-%d,%s:%d-%d", $c, $row->start, $e,
						$c, $row->start - 1, $e;
				}
				else {
					$v = sprintf "%s:%d,%s:%d", $c, $row->start, $c, $row->start - 1;
				}
				$row->value( $coord_i, $v );
			}
		);
	}

	# return
	return $Data;
}

### Merge two datasets together
sub merge_two_datasets {
	my ( $input_data1, $input_data2 ) = @_;

	# Check the input data
	my $check = check_data_tables( $input_data1, $input_data2 );
	if ( $check and not $lookup ) {
		print "\n Files have non-equal numbers of data rows! Enabling lookup\n";
	}

	# Merge by lookup values
	if ( $lookup or $check ) {

		# we need to merge by lookup values
		return merge_two_datasets_by_lookup( $input_data1, $input_data2 );
	}

	# Otherwise Merge the two datasets blindly
	# determine the new order
	my @order;
	if ($automatic) {

		# automatic selection
		@order = automatically_determine_order( $input_data1, $input_data2 );
	}
	else {
		# manual selection from user
		@order = request_new_order( $input_data1, $input_data2 );
	}

	# Initialize the output data structure if necessary
	$output_data = initialize_output_data_structure($input_data1);

	# add comment lines from second files, first was take care of during initialization
	foreach ( $input_data2->comments ) {
		$output_data->add_comment($_);
	}

	# assign the datasets to the output data in requested order
	foreach my $request (@order) {

		# add the dataset from the appropriate input file
		if ( $request =~ /^\d+$/ ) {

			# a digit indicates a dataset from file1

			# print dataset name in automatic mode
			if ($automatic) {
				printf "  Merging column %s\n", $input_data1->name($request);
			}

			# copy the dataset
			my $values = $input_data1->column_values($request);
			my $index  = $output_data->add_column($values);

			# copy the metadata
			copy_metadata( $input_data1, $request, $index );
		}
		elsif ( $request =~ /^[a-z]+$/i ) {

			# a letter indicates a dataset from file2

			# first convert back to number
			my $number = $number_of->{$request};

			# print dataset name in automatic mode
			if ($automatic) {
				printf "  Merging column %s\n", $input_data2->name($number);
			}

			# copy the dataset
			my $values = $input_data2->column_values($number);
			my $index  = $output_data->add_column($values);

			# copy the metadata
			copy_metadata( $input_data2, $number, $index );
		}
		else {
			# we should never get here
			die " unrecognized symbol '$request' in request! nothing done!\n";
		}
	}
}

### merge two datasets by lookup value
sub merge_two_datasets_by_lookup {
	my ( $input_data1, $input_data2 ) = @_;

	# determine lookup indices
	my ( $lookup_i1, $lookup_i2 ) = request_lookup_indices( $input_data1, $input_data2 );

	# determine order
	my @order;
	if ($automatic) {

		# automatic selection
		@order = automatically_determine_order( $input_data1, $input_data2 );
	}
	else {
		# manual selection from user
		@order = request_new_order( $input_data1, $input_data2 );
	}

	# add coordinate column to output if necessary
	if ( $lookup_name
		and ( $lookup_name eq 'MergeDatasetCoordinate' and not $automatic ) )
	{

		# check if we taking from file 1 or 2
		if ( $order[0] =~ /[a-z]+/i ) {
			if ( $input_data2->name($lookup_i2) eq 'MergeDatasetCoordinate' ) {
				unshift @order, $letter_of->{$lookup_i2};
			}
		}
		else {
			if ( $input_data1->name($lookup_i1) eq 'MergeDatasetCoordinate' ) {
				unshift @order, $lookup_i1;
			}
		}
	}

	# rearrange as necessary to make the first data structure dominant
	if ( $order[0] =~ /[a-z]+/ ) {

		# it looks like the user wants the second file to be dominant
		# the dominant file is where we'll be taking all of the values

		# switch the references
		( $input_data1, $input_data2 ) = ( $input_data2, $input_data1 );

		# switch the lookup indices
		( $lookup_i1, $lookup_i2 ) = ( $lookup_i2, $lookup_i1 );

		# switch the numbers and letters
		foreach (@order) {
			if (/\d+/) {
				$_ = $letter_of->{$_};
			}
			else {
				$_ = $number_of->{$_};
			}
		}
	}

	# Index the second dataset
	# we'll be putting the lookup values into an index hash
	# where the lookup value is the key and the row number is the value
	my $index = index_dataset( $input_data2, $lookup_i2 );

	# Initialize the output data structure if necessary
	$output_data = initialize_output_data_structure($input_data1);

	# add comment lines from second files, first was take care of during initialization
	foreach ( $input_data2->comments ) {
		$output_data->add_comment($_);
	}

	# assign the datasets to the output data in requested order
	foreach my $request (@order) {

		# add the dataset from the appropriate input file based on the type of request
		my $column;

		# a digit indicates a dataset from file1
		if ( $request =~ /\d+/ ) {

			# we're assuming that file1 is dominant, we copy all the
			# rows from file1 into output, no lookup required

			# print dataset name in automatic mode
			if ($automatic) {
				printf "  Merging column %s\n", $input_data1->name($request);
			}

			# copy the dataset including header
			my $values = $input_data1->column_values($request);
			$column = $output_data->add_column($values);

			# copy the metadata
			copy_metadata( $input_data1, $request, $column );
		}

		# a letter indicates a dataset from file2
		elsif ( $request =~ /[a-z]+/i ) {

			# we will have to perform the lookup here

			# first convert back to number
			my $request_number = $number_of->{$request};

			# print dataset name in automatic mode
			if ($automatic) {
				printf "  Merging column %s\n", $input_data2->name($request_number);
			}

			# add new empty column
			$column = $output_data->add_column( $input_data2->name($request_number) );

			# copy the dataset via lookup process
			if ( $input_data1->name($lookup_i1) eq 'MergeDatasetCoordinate' ) {
				foreach my $r1 ( 1 .. $input_data1->last_row ) {

					# identify the appropriate row in file2 by lookup value
					# using both primary and alternate coordinate string
					my ( $k1, $k2 ) = split /,/, $input_data1->value( $r1, $lookup_i1 );
					my $r2 = $index->{$k1} || $index->{$k2} || undef;

					# copy the appropriate value
					if ( defined $r2 ) {
						$output_data->value( $r1, $column,
							$input_data2->value( $r2, $request_number ) );
					}
					else {
						$output_data->value( $r1, $column, '.' );    # null value
					}
				}
			}
			else {
				foreach my $r1 ( 1 .. $input_data1->last_row ) {

					# identify the appropriate row in file2 by lookup value
					my $lookup_val = $input_data1->value( $r1, $lookup_i1 );
					my $r2         = $index->{$lookup_val} || undef;

					# copy the appropriate value
					if ( defined $r2 ) {
						$output_data->value( $r1, $column,
							$input_data2->value( $r2, $request_number ) );
					}
					else {
						$output_data->value( $r1, $column, '.' );    # null value
					}
				}
			}

			# copy the metadata
			copy_metadata( $input_data2, $request_number, $column );
		}

		# something else
		else {
			print " WARNING: unknown request '$request'. skipping this column\n";
			next;
		}

		# make sure we remember the lookup_column_index in the new output
		# if it is necessary
		if ( $lookup_name and not defined $output_lookup_i ) {

			# we have a lookup name, but the index for this lookup in
			# the output data structure is not yet defined

			# check if the current output index is it
			if ( $output_data->name($column) =~ /\A $lookup_name \Z/xi ) {

				# this current column matches the lookup name
				# so we will use it
				$output_lookup_i = $column;
			}
		}
	}
}

### Add datasets from one data file to the output data structure
sub add_datasets {
	my $data = shift;

	# Check the input data
	my $check = check_data_tables( $output_data, $data );

	# Determine if we need to do lookup
	my $lookup_i;
	my $index;    # lookup hash for the current data table
	if ( $check or $lookup ) {

		# we need to merge by lookup values

		# check whether we need to do the lookup manually or automatically
		if ( $manual or not $lookup_name ) {

		  # we don't have an automatic lookup name defined or user wants to do it manually
		  # so we do this by treating as two new datasets
			merge_two_datasets_by_lookup( $output_data, $data );
			return;
		}

		# index the new data reference
		# we'll be putting the lookup values into an index hash
		# where the lookup value is the key and the row number is the value
		$lookup_i = request_lookup_indices($data);
		$index    = index_dataset( $data, $lookup_i );
	}

	# add comment lines from second files, first was take care of during initialization
	foreach ( $data->comments ) {
		$output_data->add_comment($_);
	}

	# Merge the new datasets with current output
	# determine the new order
	my @order;
	if ($automatic) {

		# automatic selection
		@order = automatically_determine_order($data);

	}
	else {
		# manual selection from user
		@order = request_new_order($data);
	}

	# assign the datasets to the output data in requested order
	foreach my $request (@order) {

		# check the request
		if ( $request !~ /^\d+$/ ) {
			print " WARNING: unknown request '$request'. skipping this column\n";
			next;
		}

		# print dataset name in automatic mode
		if ($automatic) {
			printf "  Merging column %s\n", $data->name($request);
		}

		# copy the dataset by lookup or blindly
		my $column;
		if ( defined $lookup_i ) {

			# merging by lookup

			# check that we have the output lookup index
			unless ( defined $output_lookup_i ) {
				print <<MESSAGE;
 FATAL: The lookup index for column '$lookup_name' is not defined in the output!
 Please ensure you include column '$lookup_name' in the output file.
MESSAGE
				exit 1;
			}

			# add new empty column
			$column = $output_data->add_column( $data->name($request) );

			# copy the dataset via lookup process
			if ( $output_data->name($output_lookup_i) eq 'MergeDatasetCoordinate' ) {
				$output_data->iterate(
					sub {
						my $row = shift;

						# lookup the row using either primary or alternate coordinates
						# then replace if a value is found
						my ( $key1, $key2 ) = split /,/,
							$row->value($output_lookup_i);
						my $r = $index->{$key1} || $index->{$key2} || undef;
						if ($r) {
							$row->value( $column, $data->value( $r, $request ) );
						}
						else {
							$row->value( $column, '.' );    # null value
						}
					}
				);
			}
			else {
				$output_data->iterate(
					sub {
						my $row = shift;

					   # identify the appropriate row in source data table by lookup value
						my $r = $index->{ $row->value($output_lookup_i) } || undef;
						if ($r) {
							$row->value( $column, $data->value( $r, $request ) );
						}
						else {
							$row->value( $column, '.' );    # null value
						}
					}
				);
			}
		}

		# merging blindly
		else {
			# copy the dataset including header
			my $values = $data->column_values($request);
			$column = $output_data->add_column($values);
		}

		# copy the metadata
		copy_metadata( $data, $request, $column );
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
		if ( scalar @data_refs == 1 ) {

			# Only one data file
			# we'll be using only numbers to list the datasets
			print_datasets( $data_refs[0], 'number' );
		}
		elsif ( scalar @data_refs == 2 ) {

			# Two data files
			# use numbers for the first one
			print_datasets( $data_refs[0], 'number' );

			# use letters for the second one
			print_datasets( $data_refs[1], 'letter' );

		}

		# Request response from user
		my $p = "\n Enter the columns' indices in the desired final order.\n"
			. " Separate indices with commas or specify a range (start - stop):   ";
		$list = prompt($p);
		$list =~ s/\s+//g;
	}

	# Parse response
	my @order = parse_list($list);
	while ( scalar @order == 0 ) {

		# an empty order is returned if the list was unparseable
		# allow user to try again
		my $p = " Unable to parse your list into valid indices. Please try again\n   ";
		$list = prompt($p);
		$list =~ s/\s+//g;
		@order = parse_list($list);
	}

	# done
	printf " using order: %s\n", join( ', ', @order );
	return @order;
}

### Request lookup indices for two datasets
sub request_lookup_indices {
	my ( $data1, $data2 ) = @_;

	# One or two data references may be passed
	# We'll ask for the lookup value indices for each

	# Determine the lookup index values
	my $index1;
	my $index2;

	# Automatic methods of determining the lookup index
	if ( not $manual ) {

		# Only one data table that we're adding, find the lookup column
		if ( not defined $data2 and $lookup_name ) {

			# we already have a lookup name
			# just need to find same column for one dataset
			$index1 = $data1->find_column("^$lookup_name\$");
			if ( defined $index1 ) {
				printf "\n  using column $index1 (%s) as lookup index for file '%s'\n",
					$data1->name($index1), $data1->filename;
				return $index1;
			}
			else {
				printf "\n Cannot find lookup column '%s' in file '%s'!!!\n",
					$lookup_name, $data1->filename;

				# try some other possibilities?
				if ( $lookup_name eq 'MergeDatasetCoordinate' ) {
					foreach my $name (qw(coordinate primary_id)) {
						$index1 = $data1->find_column("^$name\$");
						if ( defined $index1 ) {
							printf " Trying column $index1 (%s) as lookup index\n",
								$data1->name($index1);
							return $index1;
						}
					}
				}
				exit 1;
			}
		}

		# First try some known column identifiers we could use automatically
		# don't forget to add the user-requested lookup name
		my @name_list;
		if ($use_coordinate) {
			@name_list = qw(MergeDatasetCoordinate coordinate id);
		}
		else {
			@name_list = qw(name primary_id gene.?id transcript.?id gene transcript);
		}
		if ($user_lookup_name) {
			unshift @name_list, $user_lookup_name;
		}
		foreach my $name (@name_list) {

			# identify possibilities
			my $possible_index1 = $data1->find_column("^$name\$");
			my $possible_index2 = $data2->find_column("^$name\$");

			# check if something was found
			if ( defined $possible_index1 and defined $possible_index2 ) {

				# assign
				$index1      = $possible_index1;
				$index2      = $possible_index2;
				$lookup_name = $name;              # for future lookups

				# report
				printf "\n  using column $index1 (%s) as lookup index for file '%s'\n",
					$data1->name($index1), $data1->filename;
				printf "  using column $index2 (%s) as lookup index for file '%s'\n",
					$data2->name($index2), $data2->filename;

				# don't go through remaining list
				last;
			}
		}
	}

	if ( not defined $index1 and not defined $index2 ) {

		# Automatic identification didn't work
		# must bother the user for help
		if ($automatic) {
			print <<MESSAGE;

 Unable to identify appropriate lookup columns automatically!
 Use --coordinate option if coordinates should be used for lookup.
 Otherwise, please execute interactively to identify lookup columns.
MESSAGE
			exit 1;
		}

		unless ($manual) {
			print <<MESSAGE;

 Unable to identify appropriate lookup columns automatically!
 If using coordinates, re-execute with the --coordinate option.
MESSAGE
		}

		# Print the index headers
		# use numbers for the first one
		print_datasets( $data1, 'number', 1 );

		# use letters for the second one
		print_datasets( $data2, 'letter', 1 );

		# Request first index responses from user
		my $p = ' Enter the unique identifier index for lookup in the first file:  ';
		while ( not defined $index1 ) {
			$index1 = prompt($p);

			# check that it's valid
			unless ( $index1 =~ /^\d+$/ and exists $data1->{$index1} ) {
				print "  unknown index value! Try again\n";
				undef $index1;
			}
		}

		# Request second index responses from user
		$p = ' Enter the unique identifier index for lookup in the second file:  ';
		while ( not defined $index2 ) {
			$index2 = prompt($p);

			# check that it's valid
			unless ( $index2 =~ /^[a-z]+$/i ) {
				print "  index value must be a letter! Try again\n";
				undef $index2;
				next;
			}
			$index2 = $number_of->{$index2};    # convert to a number
												# check that it's valid
			unless ( exists $data2->{$index2} ) {
				print "  unknown index value! Try again\n";
				undef $index2;
			}
		}
	}

	# done
	return $index1, $index2;
}

### Automatically determine the order
sub automatically_determine_order {

	# get the two data structures
	# we could have two, or only one
	my ( $data1, $data2, $number );
	if ( scalar @_ == 1 ) {

		# comparing new data with the output data we are building
		$data1  = $output_data;
		$data2  = shift @_;
		$number = 1;
	}
	else {
		# two new datasets
		$data1  = shift @_;
		$data2  = shift @_;
		$number = 2;
	}

	# generate quick hash of names to exclude from data1
	my %exclude;
	for my $i ( 1 .. $data1->number_columns ) {
		$exclude{ $data1->name($i) } = 1;
	}

	# generate the order
	my @order;

	# add the first dataset if provided
	if ( $number == 2 ) {

		# we will automatically take all of the columns from the first data
		for my $i ( 1 .. $data1->number_columns ) {
			push @order, $i;
		}
	}

	# automatically identify those columns in second data not present in
	# the first one to take
	for my $i ( 1 .. $data2->number_columns ) {

		# name is score, take it
		if ( $data2->name($i) eq 'Score' ) {

			# take number or letter depending on number of data files
			if ( $number == 1 ) {
				push @order, $i;
			}
			else {
				# two data files, push letter
				push @order, $letter_of->{$i};
			}
		}

		elsif ( exists $exclude{ $data2->name($i) } ) {

			# non-unique column, skip
			next;
		}

		else {
			# must be a unique column, take it
			# take number or letter depending on number of data files
			if ( $number == 1 ) {
				push @order, $i;
			}
			else {
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
	foreach my $item ( split /,/, $string ) {

		# comma-delimited list of items
		# may contain ranges of incremental items, check for those

		if ( $item =~ /\-/ ) {

			# item is a range, specified as 'start - stop'

			if ( $item =~ /^(\d+) \- (\d+)$/x ) {

				# range contains numbers, so must be from file1
				for ( my $i = $1; $i <= $2; $i++ ) {

					# we will loop through from specified start to stop
					# add each number to the order
					push @list, $i;
				}
			}

			elsif ( $item =~ /^ ([a-z]+) \- ([a-z]+) $/x ) {

				# range does not contain numbers, so must be from file2
				for ( my $i = $number_of->{$1}; $i <= $number_of->{$2}; $i++ ) {

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
				return;
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
	my ( $data, $index_type, $include_coordinate ) = @_;

	# array to be used in the default order
	my @order;

	# print the dataset names for this datafile
	printf "\n These are the headers in file '%s'\n", $data->filename;

	# print Numbers
	if ( $index_type eq 'number' ) {
		for my $i ( 1 .. $data->number_columns ) {

			# skip the coordinate
			if ( $data->name($i) eq 'MergeDatasetCoordinate'
				and not $include_coordinate )
			{
				next;
			}

			# print the dataset name and it's index
			# and record the index in the order array for the default order
			printf " %4s  %s\n", $i, $data->name($i);
			push @order, $i;
		}
	}

	# print letters
	else {
		for my $i ( 1 .. $data->number_columns ) {

			# skip the coordinate
			if ( $data->name($i) eq 'MergeDatasetCoordinate'
				and not $include_coordinate )
			{
				next;
			}

			# use letters instead of numbers as the index key
			my $letter = $letter_of->{$i};
			printf " %4s  %s\n", $letter, $data->name($i);
			push @order, $letter;
		}
	}
	return @order;
}

### Check the data tables for similarity
sub check_data_tables {
	my ( $input_data1, $input_data2 ) = @_;

	# Check the feature types
	if (    $input_data1->feature
		and $input_data2->feature
		and $input_data1->feature ne $input_data2->feature )
	{
		# Each file has feature type defined, but they're not the same
		# probably really don't want to combine these
		print " WARNING! The metadata feature types for both files don't match!!!\n";
		print "   Continuing anyway...\n";
	}

	# check line numbers
	my $status = 0;
	if ( $input_data1->last_row != $input_data2->last_row ) {

		# the number of rows in each data table don't equal
		# we will need to do this by lookup
		$status = 1;
	}
	return $status;
}

### Index the values in a data table
sub index_dataset {
	my ( $data, $lookup_i ) = @_;

	# load up the lookup index hash for the current data table
	my %index;
	my $index_warning = 0;
	if ( $data->name($lookup_i) eq 'MergeDatasetCoordinate' ) {
		$data->iterate(
			sub {
				my $row = shift;
				my ( $key1, $key2 ) = split /,/, $row->value($lookup_i);
				die "no primary,alternate coordinates!"
					unless ( defined $key1 and defined $key2 );
				if ( exists $index{$key1} ) {
					$index_warning++;
				}
				else {
					$index{$key1} = $row->row_index;
				}
				if ( exists $index{$key2} ) {
					$index_warning++;
				}
				else {
					$index{$key2} = $row->row_index;
				}
			}
		);
	}
	else {
		$data->iterate(
			sub {
				my $row = shift;
				my $key = $row->value($lookup_i);
				if ( exists $index{$key} ) {
					$index_warning++;
				}
				else {
					$index{$key} = $row->row_index;
				}
			}
		);
	}
	if ($index_warning) {
		my $f = $data->filename;
		my $n = $data->name($lookup_i);
		print <<MESSAGE;
  WARNING: $index_warning rows had two or more duplicate lookup values
  for column '$n' (index $lookup_i) in file '$f'.
  Only the first occurence was used.
MESSAGE

	}
	return \%index;
}

### Initialize the output data structure
sub initialize_output_data_structure {
	my $Data1 = shift;

	# generate brand new output data structure
	my $output_data1 = Bio::ToolBox::Data->new(
		feature => $Data1->feature,    # re-use the same feature as file1
	);

	# we'll re-use the values from file1
	$output_data1->program( $Data1->program );    # force overwrite value
	$output_data1->database( $Data1->database );
	foreach ( $Data1->comments ) {
		$output_data1->add_comment($_);
	}

	# assign the filename
	if ($outfile) {

		# use the new output file name
		$output_data1->add_file_metadata($outfile);
	}
	else {
		# borrow the first file name
		$output_data1->add_file_metadata( $Data1->filename );    # borrow file name
	}

	return $output_data1;
}

### Copy the dataset metadata
sub copy_metadata {

	# collect arguments
	my ( $data, $request, $index ) = @_;

	# copy the metadata
	my %md = $data->metadata($request);
	foreach my $k ( keys %md ) {
		next if $k eq 'index';
		next if $k eq 'name';
		$output_data->metadata( $index, $k, $md{$k} );
	}

	# check if should rename the dataset
	if ( $automatic and $data->name($request) eq 'Score' ) {

		# only in automatic mode and the dataset is a generic Score
		# no opportunity to rename interactively
		# use file basename appended with Score
		$output_data->name( $index, $data->basename . '_Score' );
	}

	# add original filename
	$output_data->metadata( $index, 'original_file', $data->filename );
}

### Re-name the dataset names
sub rename_dataset_names {
	print " For each header, type a new name or push enter to accept the current name\n";

	for my $i ( 1 .. $output_data->number_columns ) {

		# walk through the list of columns

		# request user input
		# print the current name and it's originating file name
		my $p = sprintf "  %s: %s ->  ", $output_data->metadata( $i, 'original_file' ),
			$output_data->name($i);
		my $new_name = prompt( $p, q() );

		# rename as requested
		if ($new_name) {
			$output_data->name( $i, $new_name );
		}
	}
}

### Generate lookup hashes for converting numbers to letters and back
sub get_number_letters {

	# hashes
	my %n2l;    # numbers to letters
	my %l2n;    # letters to numbers
	my %lookup = (
		0  => 'a',
		1  => 'b',
		2  => 'c',
		3  => 'd',
		4  => 'e',
		5  => 'f',
		6  => 'g',
		7  => 'h',
		8  => 'i',
		9  => 'j',
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
	for my $i ( 0 .. 701 ) {

		# this gives range from a to zz

		# check the first letter
		if ( $i % 26 == 0 ) {
			$first = $lookup{ ( $i / 26 ) - 1 };
		}

		# generate the letter
		# two letters [null..z][a..z]
		my $letter;
		if ( defined $first ) {
			$letter = sprintf "%s%s", $first, $lookup{ $i % 26 };
		}
		else {
			$letter = $lookup{ $i % 26 };
		}

		# store in hashes, converting 0-base numbers to 1-base numbers
		$n2l{ $i + 1 } = $letter;
		$l2n{$letter} = $i + 1;
	}

	return ( \%n2l, \%l2n );
}

__END__

=head1 NAME

merge_datasets.pl

A program to merge columns of two or more data files.

=head1 SYNOPSIS

merge_datasets.pl [--options...] <file1> <file2> ...
  
  Merge options:
  -l --lookup                      lookup rows by specific value
  -a --auto                        automatically select unique columns
  -m --manual                      manually select lookup and each data column
  -x --index <number,letter,range> specify columns: 1st file numbers, 2nd letters 
  -n --lookupname <text>           name of lookup column
  -c --coordinate                  use coordinate string as lookup
  
  General options:
  -o --out <filename>              default first file
  -z --gz                          compress output
  -v --version                     print version and exit
  -h --help                        show full documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 Merge options

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

=item --lookupname E<lt>textE<gt>

=item --lun E<lt>textE<gt>

Provide an alternate column name to identify the columns automatically 
in the input files containing the lookup values when performing the 
lookup. Each file should have the same lookup column name. Default 
values include 'Name', 'ID', 'Transcript', or 'Gene'.

=item --coordinate

Use genomic coordinates when performing a lookup. Both 1-base and 0-base
coordinate files can be safely and effectively matched. The generated
coordinates column is temporary. This is automatically enabled when working
with coordinate files, such as BED and GFF. This may be disabled by
specifying --nocoordinate, or simply specifying a column lookup name.

=back

=head2 General options

=over 4

=item --out E<lt>filenameE<gt>

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
it under the terms of the Artistic License 2.0.  
