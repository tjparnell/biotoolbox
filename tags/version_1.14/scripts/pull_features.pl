#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::data_helper qw(find_column_index);
use Bio::ToolBox::file_helper qw(
	open_tim_data_file
	load_tim_data_file
	write_tim_data_file
	write_summary_data
);
my $VERSION = '1.14';

print "\n A script to pull out specific features from a data file\n";

### Quick help
unless (@ARGV) { # when no command line options are present
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Get command line options
my (
	$datafile, 
	$listfile,
	$outfile,
	$data_index,
	$list_index,
	$order,
	$sum,
	$startcolumn,
	$stopcolumn,
	$log,
	$help,
	$print_version,
);
GetOptions( 
	'data=s'     => \$datafile, # the input data file
	'list=s'     => \$listfile, # the list file
	'out=s'      => \$outfile, # the new output file name
	'dindex=i'   => \$data_index, # index to look up in the data file
	'lindex=i'   => \$list_index, # index of look up values in list file
	'order=s'    => \$order, # the order to keep the values
	'sum'        => \$sum, # flag to re-sum the pulled values
	'starti=i'   => \$startcolumn, # index of column to start summarizing
	'stopi=i'    => \$stopcolumn, # index of column to stop summarizing
	'log!'       => \$log, # values are in log, respect log status
	'help'       => \$help, # flag to print help
	'version'    => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script pull_features.pl, version $VERSION\n\n";
	exit;
}



### Check for required values

unless (defined $datafile) {
	die " no input data file specified!\n";
}

unless (defined $listfile) {
	die " no list data specified!\n";
}

unless (defined $outfile) {
	die " no output data file name given!\n";
}

my $list_order;
if ($order) {
	if ($order eq 'list') {
		$list_order = 1;
	}
	elsif ($order eq 'data') {
		$list_order = 0;
	}
	else {
		die " unrecognized order request '$order'! Enter list or data\n";
	}
}
else {
	# default is the list order
	$list_order = 1;
}



### Open files
# file handles and metadata
my $list_data = load_tim_data_file($listfile) or
	die " unable to open list file!\n";

my ($data_fh, $data_md) = open_tim_data_file($datafile) or
	die " unable to open data file!\n";



### Determine indices
identify_indices();



### Load the list of specified values
print " Collecting lookup values from file '$listfile'...\n";
my ($requests, $pulled) = collect_request_list();
	# these are global references to two data hashes
	# the first is for a lookup for the feature requests and identify the group number
	# the second is for storing all the found pulled data


### Pull out the desired features
print " Pulling features...\n";
my ($found_count, $notfound_count) = pull_requested_features();
print "  $found_count features were found and pulled\n";
if ($notfound_count > 0) {
	print "  $notfound_count features were not found\n";
}



### Write the output files
if ($found_count) {
	# first generate output structures
	generate_output_data_structures();
	# then
	write_files();
}
else {
	print "  Nothing found! Nothing to write!\n";
}





########################   Subroutines   ###################################

### Identify the indices
sub identify_indices {
	
	# first check whether we have a simple list file
	if ($list_data->{'number_columns'} == 1) {
		# simple one-column file
		# the answer is obvious
		print "  list file only has 1 column, using it\n";
		$list_index = 0;
		# be VERY careful, though, it may not have a column header name
	}
	elsif ($listfile =~ /\.kgg$/i) {
		# a cluster gene file
		# the answer is obvious
		print "  using .kgg file as list\n";
		$list_index = 0;
		# KGG files have column headers
	}
		
	
	# look for the corresponding list index if data was specified
	if (defined $data_index and !defined $list_index) {
		# we have the data index but need the list index
		
		# get the column header name
		my $lookup = $data_md->{$data_index}{'name'};
		
		# find it in the list
		my $possible = find_column_index($list_data, "^$lookup\$");
		
		# check
		if (defined $possible) {
			# found something
			$list_index = $possible;
			print "  found column '", $list_data->{$list_index}{'name'}, 
				"', using list index $list_index\n";
			return;
		}
	}
	
	# look for the corresponding data index if list was specified
	elsif (!defined $data_index and defined $list_index) {
		# we have the list index but need the data index
		
		# get the column header name
		my $lookup = $list_data->{$list_index}{'name'};
		
		# find it in the list
		my $possible = find_column_index($data_md, "^$lookup\$");
		
		# check
		if (defined $possible) {
			# found something
			$data_index = $possible;
			print "  found column '", $data_md->{$data_index}{'name'}, 
				"', using data index $data_index\n";
			return;
		}
		else {
			# did not find something
			# is it possible it is a simple list of names?
			if ($list_data->{'number_columns'} == 1) {
				# we likely only have a simple list of features
				# adjust the data table accordingly by adding a name
				# this will avoid not finding the first element in the output file
				# if this is not true, then the user will just get an error that 
				# one feature cannot be found.....
				unshift @{ $list_data->{'data_table'} }, 'Name';
				$list_data->{0}{'name'} = 'Name';
				$list_data->{'last_row'}++;
			}
		}
	}
	
	# neither was specified
	elsif (!defined $data_index and !defined $list_index) {
		
		# try some known column identifiers
		foreach my $name (qw(name id transcript gene)) {
			
			# identify possibilities
			my $possible_data_i = find_column_index($data_md, "^$name\$");
			my $possible_list_i = find_column_index($list_data, "^$name\$");
			
			# check if something was found
			if (defined $possible_data_i and defined $possible_list_i) {
				
				# assign
				$data_index = $possible_data_i;
				$list_index = $possible_list_i;
				
				# report
				print " found identical columns\n";
				print "  using list column '", $list_data->{$list_index}{'name'}, 
					"', index $list_index\n";
				print "  using data column '", $data_md->{$data_index}{'name'}, 
					"', index $data_index\n";
				return;
			}
		}
	}
	
	# End automatic guessing of index numbers, ask the user
	foreach ( ['list', $list_data], ['data', $data_md] ) {
		# a complicated foreach loop to look for either one
		
		# the type and metadata structure
		my ($type, $metadata) = @$_;
		
		# skip the ones we already know
		next if ($type eq 'list' and defined $list_index);
		next if ($type eq 'data' and defined $data_index);
		
		# print the headers
		print "\n These are the columns in the $type file.\n";
		for (my $i = 0; $i < $metadata->{'number_columns'}; $i++) {
			print "   $i\t$metadata->{$i}{name}\n";
		}
	
		# process the answer
		print " Enter the unique identifier lookup column index from the $type file    ";
		my $answer = <STDIN>;
		chomp $answer;
		
		# check answer
		if (exists $metadata->{$answer}) {
			# answer appears to be a column index
			if ($type eq 'list') {
				$list_index = $answer;
			}
			else {
				$data_index = $answer;
			}
		}
		else {
			die " Invalid response!\n";
		}
	}
}



### Subroutine to collect list values from a file
sub collect_request_list {
	
	my %requests; # the identifier values to look up
	# request{ unique_id } = group#
		# for KGG lists where we are splitting each of the groups into 
		# separate files, we need to know how many and which ones
		# can't trust whether all groups are in the KGG file or just a few
	my %pulled;
	# pulled_data{ group# } -> { unique_id } = [ line_data ]
	# still need a list of the order:
	# pulled_data{ group# } -> { 'feature_order' } = [unique_id,...]
	
	
	# Check whether we're working with a Cluster .kgg file
	if ($listfile =~ /\.kgg$/i) {
		for my $row (1 .. $list_data->{'last_row'}) {
			
			my $id    = $list_data->{'data_table'}->[$row][0];
			my $group = $list_data->{'data_table'}->[$row][1];
			
			# store the identifier in the requests hash
			# the gene identifier is the key, the cluster group number is the value
			$requests{$id} = $group;
			
			# prepare the pulled data hash
			$pulled{$group}->{$id} = [];
			
			# record the ID for use in writing in the file list order
			if (exists $pulled{$group}{'list_order'} ) {
				push @{ $pulled{$group}{'list_order'} }, $id;
			}
			else {
				$pulled{$group}{'list_order'} = [ ($id) ];
			}
		}
	}
	
	# otherwise a simple text file
	else {
		# prepare the order array
		# using group ID of 0
		$pulled{0}{'feature_order'} = [];
		
		# collect the identifiers
		for my $row (1 .. $list_data->{'last_row'}) {
			my $id = $list_data->{'data_table'}->[$row][$list_index];
			
			$requests{$id} = 0;
			$pulled{0}{$id} = [];
			push @{ $pulled{0}{'list_order'} }, $id;
		}
	}
	
	# prepare the dump array for keeping the features in data file order 
	# instead of list file order
	foreach my $group (keys %pulled) {
		$pulled{$group}{'data_dump'} = [];
	}
	
	return (\%requests, \%pulled);
}



### Subroutine to generate the output data structures
sub generate_output_data_structures {
	
	# we will simply duplicate the current data metadata structure 
	# for each output data structure
	# there may be more than one output data structure, primarily with 
	# kgg source files where each group will be put into a new file
	
	# we will store the output data structure in the %pulled hash structure
	# under the key 'output_data' under the appropriate group id
	
	# only one output data structure
	if (scalar keys %{$pulled} == 1) {
		
		# copy
		my %new_md = %{ $data_md };
		
		# update the file name data
		my $newfile = $outfile;
		my $ext = $data_md->{'extension'};
		unless ($newfile =~ m/$ext\Z/) {
			# add the current extension if necessary
			$newfile .= $ext;
		}
		$new_md{'filename'} = $newfile;
		$new_md{'basename'} = undef;
		$new_md{'path'} = undef;
		
		# add the data table
		$new_md{'data_table'} = [];
		# add the column headers
		push @{ $new_md{'data_table'} }, [ @{ $data_md->{'column_names'} } ];
		
		# store the structure
		$pulled->{0}{'output_data'} = \%new_md;
	}
	
	# multiple data structures
	else {
		foreach my $group (keys %$pulled) {
			# copy
			my %new_md = %{ $data_md };
			
			# update the file name data
			my $newfile = $outfile;
			my $ext = $data_md->{'extension'};
			$newfile =~ s/$ext\Z//; # remove the extension if present
			$new_md{'filename'} = $newfile . '_g' . $group . $ext;
			$new_md{'basename'} = undef;
			$new_md{'path'} = undef;
			
			# add the data table
			$new_md{'data_table'} = [];
			# add the column headers
			push @{ $new_md{'data_table'} }, [ @{ $data_md->{'column_names'} } ];
			
			# store the structure
		$pulled->{$group}{'output_data'} = \%new_md;
		}
	}
}



### Subroutine to pull the requested features
sub pull_requested_features {
	
	# we will now walk through the data file and pull out those lines 
	# which match the feature identifier
	my $found = 0;
	my $notfound = 0;
	while (my $line = $data_fh->getline) {
		$line =~ s/[\r\n]+$//; # chomp
		my @d = split /\t/, $line;
		
		# check if the identifier exists in our requests hash
		if (exists $requests->{ $d[$data_index] } ) {
			# found
			my $id = $d[$data_index];
			my $group = $requests->{ $d[$data_index] };
						
			# record the line data
			if ($list_order) {
				# store by ID so that it can be sorted later
				push @{ $pulled->{$group}{$id} }, \@d;
			}
			else {
				# just dump it in the same order as the data file
				push @{ $pulled->{$group}{'data_dump'} }, \@d;
			}
			$found++;
		}
		else {
			# not found
			$notfound++;
		}
	}
	
	return ($found, $notfound);
}



### Subroutine to write the output files
sub write_files {
	# Write the files
	
	foreach my $group (keys %$pulled) {
		
		# de-reference the output data for ease
		my $data = $pulled->{$group}{'output_data'};
		
		# copy the pulled data lines into the output data structure
		if ($list_order) {
			# we will use the list order in the output file
			foreach my $id (@{ $pulled->{$group}{'list_order'} }) {
				# this is a list of the feature ids in the same order as the input 
				# list file 
				push @{ $data->{'data_table'} }, @{ $pulled->{$group}{$id} };
			}
		}
		else {
			# we will use the data file order in the output file
			# just copy the contents over
			push @{ $data->{'data_table'} }, @{ $pulled->{$group}{'data_dump'} };
		}
		$data->{'last_row'} = scalar @{ $data->{'data_table'} } - 1;
		
		# write the file
		my $write_results = write_tim_data_file(
			'data'      => $data,
			# no file name, using the filename recorded in the out data
		);
		if ($write_results) {
			print "  Wrote new datafile '$write_results'\n";
		}
		else {
			print "  Unable to write datafile '$write_results'!!!\n";
		}
		
		# Summarize the pulled data
		if ($sum) {
			print " Generating final summed data...\n";
			my $sumfile = write_summary_data(
				'data'         => $data,
				'filename'     => $data->{'filename'}, # must provide
				'startcolumn'  => $startcolumn,
				'endcolumn'    => $stopcolumn,
				'log'          => $log,
			);
			if ($sumfile) {
				print "  Wrote summary file '$sumfile'\n";
			}
			else {
				print "  Unable to write summary file!\n";
			}
		}
	}
}



__END__

=head1 NAME

pull_features.pl

A script to pull out a specific list of data rows from a data file.

=head1 SYNOPSIS

pull_features.pl --data <filename> --list <filename> --out <filename>
  
  Options:
  --data <filename>
  --list <filename>
  --out <filename>
  --dindex <integer>
  --lindex <integer>
  --order [list | data]
  --sum
  --starti <integer>
  --stopi <integer>
  --log
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --data

Specify a tab-delimited text file as the data source file. One of 
the columns in the input file should contain the identifiers to be 
used in the lookup. The file may be gzipped.

=item --list

Specify the name of a text file containing the list of feature 
names or identifiers to pull. The file may be a single column or 
tab-delimited multi-column file with column headers. A .kgg file 
from a Cluster k-means analysis may be used.

=item --out

Specify the output file name. 

=item --dindex <integer>

=item --lindex <integer>

Specify the index numbers of the columns in the data and list 
files, respectively, containing the identifiers to match features. 
If not specified, then the program will attempt to identify  
appropriate matching columns with the same header name. If none 
are specified, the user must select interactively from a list of 
available column names. 

=item --order [list | data]

Optionally specify the order of features in the output file. Two 
options are available. Specify 'list' to match the order of features 
in the list file. Or specify 'data' to match the order of features 
in the data file. The default is list.

=item --sum

Indicate that the pulled data should be averaged across all 
features at each position, suitable for graphing. A separate text 
file with '_summed' appended to the filename will be written.

=item --starti <integer>

When re-summarizing the pulled data, indicate the start column 
index that begins the range of datasets to summarize. Defaults 
to the leftmost column without a standard feature description
name.

=item --stopi <integer>

When re-summarizing the pulled data, indicate the stop column
index the ends the range of datasets to summarize. Defaults
to the last or rightmost column.

=item --log

The data is in log2 space. Only necessary when re-summarizing the
pulled data.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

Given a list of requested unique feature identifiers, this program will 
pull out those features (rows) from a datafile and write a new file. This 
program compares in function to a popular spreadsheet VLOOKUP command. 
The list is provided as a separate text file, either as a single column 
file or a multi-column tab-delimited from which one column is selected. 
All rows from the source data file that match an identifier in the list 
will be written to the new file. The order of the features in the output 
file may match either the list file or the data file. 

The program will also accept a Cluster gene file (with .kgg extension) 
as a list file. In this case, all of the genes for each cluster are 
written into separate files, with the output file name appended with the 
cluster number. 

The program will optionally regenerate a summed data file, in which values 
in the specified data columns are averaged and written out as rows in a 
separate data file. Compare this function to the summary option in the 
biotoolbox scripts L<get_relative_data.pl> or L<average_gene.pl>.
 
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
