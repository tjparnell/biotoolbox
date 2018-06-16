#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use Bio::ToolBox::Data;
use Bio::ToolBox::utility;
my $VERSION =  '1.60';

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
	$group_index,
	$order,
	$sum,
	$sum_only,
	$startcolumn,
	$stopcolumn,
	$log,
	$help,
	$print_version,
);
GetOptions( 
	'd|data=s'     => \$datafile, # the input data file
	'l|list=s'     => \$listfile, # the list file
	'o|out=s'      => \$outfile, # the new output file name
	'x|dindex=i'   => \$data_index, # index to look up in the data file
	'X|lindex=i'   => \$list_index, # index of look up values in list file
	'g|gindex=i'   => \$group_index, # index of group in list file
	'r|order=s'    => \$order, # the order to keep the values
	'U|sum!'       => \$sum, # flag to re-sum the pulled values
	'sumonly!'     => \$sum_only, # only save the summary file
	'starti=i'     => \$startcolumn, # index of column to start summarizing
	'stopi=i'      => \$stopcolumn, # index of column to stop summarizing
	'log!'         => \$log, # values are in log, respect log status
	'h|help'       => \$help, # flag to print help
	'v|version'    => \$print_version, # print the version
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
	print " Biotoolbox script pull_features.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
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

$sum = 1 if $sum_only;

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
my $List = Bio::ToolBox::Data->new(file => $listfile) or
	die " unable to open list file!\n";

my $Data = Bio::ToolBox::Data->new(file => $datafile) or
	die " unable to open data file!\n";



### Determine indices



### Load the list of specified values
print " Collecting lookup values from file '$listfile'...\n";
identify_indices();
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
	if ($List->number_columns == 1) {
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
		$group_index = 1;
		# KGG files have column headers
	}
		
	
	# look for the corresponding list index if data was specified
	if (defined $data_index and !defined $list_index) {
		# we have the data index but need the list index
		
		# get the column header name
		my $lookup = $Data->name($data_index);
		
		# find it in the list
		my $possible = $List->find_column("^$lookup\$");
		
		# check
		if (defined $possible) {
			# found something
			$list_index = $possible;
			printf "  found column '%s', using list index $list_index\n", 
				$List->name($list_index);
			return;
		}
	}
	
	# look for the corresponding data index if list was specified
	elsif (!defined $data_index and defined $list_index) {
		# we have the list index but need the data index
		
		# get the column header name
		my $lookup = $List->name($list_index);
		
		# find it in the list
		my $possible = $Data->find_column("^$lookup\$");
		
		# check
		if (defined $possible) {
			# found something
			$data_index = $possible;
			printf "  found column '%s', using data index $data_index\n", 
				$Data->name($data_index);
			return;
		}
		else {
			# did not find something
			# is it possible it is a simple list of names?
			if ($List->number_columns == 1) {
				# we likely only have a simple list of features
				# adjust the data table accordingly by adding a name
				# this will avoid not finding the first element in the output file
				# if this is not true, then the user will just get an error that 
				# one feature cannot be found.....
				#### WE ARE MESSING WITH THE INTERNALS OF THE OBJECT HERE
				#### DONT DO THIS!!!!!!!
				unshift @{ $List->{'data_table'} }, 'Name';
				$List->{0}{'name'} = 'Name';
				$List->{'last_row'}++;
			}
		}
	}
	
	# neither was specified
	elsif (!defined $data_index and !defined $list_index) {
		
		# just try the built-in name column index
		# this uses a few common names
		$data_index = $Data->name_column;
		$list_index = $List->name_column;
		
		if (defined $data_index and defined $list_index) {
			# report
			printf "  using list column '%s', index $list_index\n", 
				$List->name($list_index);
			printf "  using data column '%s', index $data_index\n", 
				$Data->name($data_index);
			return;
		}
		else {
			# nope, do not have a match, forget our guesses
			undef $data_index;
			undef $list_index;
		}
	}
	
	# check for group number
	if ($List->number_columns > 1 and not defined $group_index) {
		my $i = $List->find_column('group');
		if (defined $i) {
			$group_index = $i;
		}
		# do we ask for a group or not????? probably not.... keep original functionality
	}
	
	# End automatic guessing of index numbers, ask the user
	unless (defined $list_index) {
		$list_index = ask_user_for_index($List, 
			" Enter the unique identifier lookup column index from the list file    ");
	}
	unless (defined $data_index) {
		$data_index = ask_user_for_index($Data, 
			" Enter the unique identifier lookup column index from the data file    ");
	}
	printf " We are using list lookup index $list_index, %s, and data lookup index $data_index, %s\n",
		$List->name($list_index), $Data->name($data_index);
}



### Subroutine to collect list values from a file
sub collect_request_list {
	
	my %requests; # the identifier values to look up
	# request{ unique_id } = group#
		# for KGG lists where we are splitting each of the groups into 
		# separate files, we need to know how many and which ones
		# can't trust whether all groups are in the KGG file or just a few
	my %pulled;
		# hash pulled{ group# } -> { unique_id } = [ line_data ]
		# still need a list of the order:
		# hash pulled{ group# } -> { 'feature_order' } = [unique_id,...]
	
	
	# check if we have multiple groups to work with
	if (defined $group_index) {
		$List->iterate( sub {
			my $row = shift;
			my $id    = $row->value($list_index);
			my $group = $row->value($group_index);
			
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
		} );
	}
	
	# otherwise a simple one group list file
	else {
		# prepare the order array
		# using group ID of 0
		$pulled{0}{'feature_order'} = [];
		
		# collect the identifiers
		$List->iterate( sub {
			my $row = shift;
			my $id = $row->value($list_index);
			
			$requests{$id} = 0;
			$pulled{0}{$id} = [];
			push @{ $pulled{0}{'list_order'} }, $id;
		} );
	}
	
	# prepare the dump array for keeping the features in data file order 
	# instead of list file order
	foreach my $group (keys %pulled) {
		$pulled{$group}{'data_dump'} = [];
	}
	
	return (\%requests, \%pulled);
}



### Subroutine to pull the requested features
sub pull_requested_features {
	
	# we will now walk through the data file and pull out those lines 
	# which match the feature identifier
	my $found = 0;
	my $notfound = 0;
	$Data->iterate( sub {
		my $row = shift;
		
		# check if the identifier exists in our requests hash
		if (exists $requests->{ $row->value($data_index) } ) {
			# found
			my $id = $row->value($data_index);
			my $group = $requests->{$id};
						
			# record the line data
			if ($list_order) {
				# store by ID so that it can be sorted later
				push @{ $pulled->{$group}{$id} }, $row;
			}
			else {
				# just dump it in the same order as the data file
				push @{ $pulled->{$group}{'data_dump'} }, $row;
			}
			$found++;
		}
		else {
			# not found
			$notfound++;
		}
	} );
	
	return ($found, $notfound);
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
		
		# generate new file name data
		my $newfile = $outfile;
		my $ext = $Data->extension;
		unless ($newfile =~ m/$ext\Z/) {
			# add the current extension if necessary
			$newfile .= $ext;
		}
		
		# duplicate and store away
		my $New = $Data->duplicate;
		$New->add_file_metadata($newfile);
		$pulled->{0}{'output_data'} = $New;
	}
	
	# multiple data structures
	else {
		foreach my $group (keys %$pulled) {
			
			# generate new file name
			my $newfile = $outfile;
			my $ext = $Data->extension;
			$newfile =~ s/$ext\Z//; # remove the extension if present
			$newfile .= '_g' . $group . $ext;
			
			# duplicate and store away
			my $New = $Data->duplicate;
			$New->add_file_metadata($newfile);
			$pulled->{$group}{'output_data'} = $New;
		}
	}
}



### Subroutine to write the output files
sub write_files {
	# Write the files
	
	foreach my $group (keys %$pulled) {
		
		# de-reference the output data for ease
		my $group_Data = $pulled->{$group}{'output_data'};
		
		# copy the pulled data lines into the output data structure
		if ($list_order) {
			# we will use the list order in the output file
			foreach my $id (@{ $pulled->{$group}{'list_order'} }) {
				# this is a list of the feature ids in the same order as the input 
				# list file 
				foreach ( @{ $pulled->{$group}{$id} } ) {
					$group_Data->add_row($_);
				}
			}
		}
		else {
			# we will use the data file order in the output file
			# just copy the contents over
			foreach ( @{ $pulled->{$group}{'data_dump'} } ) {
				$group_Data->add_row($_);
			}
		}
		
		# write the file
		if (not $sum_only) {
			my $write_results = $group_Data->save;
				# no file name, using the filename recorded in the out data
			if ($write_results) {
				print " Wrote new datafile '$write_results'\n";
			}
			else {
				print " Unable to write datafile '$write_results'!!!\n";
			}
		}
		
		# Summarize the pulled data
		if ($sum) {
			my $sumfile = $group_Data->summary_file(
				'startcolumn'  => $startcolumn,
				'endcolumn'    => $stopcolumn,
				'log'          => $log,
			);
			if ($sumfile) {
				print " Wrote summary file '$sumfile'\n";
			}
			else {
				print " Unable to write summary file!\n";
			}
		}
	}
}



__END__

=head1 NAME

pull_features.pl

A program to pull out a specific list of data rows from a data file.

=head1 SYNOPSIS

pull_features.pl --data <filename> --list <filename> --out <filename>
  
  File options:
  -d --data <filename>          Source of all data rows or features
  -l --list <filename>          List of specific row or feature names
  -o --out <filename>           Output file, or basename for group files
  
  Column index options:
  -x --dindex <index>           Data column index of row name to lookup
  -X --lindex <index>           List column index of row name to lookup
  -g --gindex <index>           Group column index for lookup
  
  Output options:
  -r --order [list | data]      Order of items in output based on
  -U --sum                      Generate a summary file
  --sumonly                     Skip output, just make a summary file
  --start <integer>             First data column to make a summary file
  --stopi <integer>             Last data column to make a summary file
  --log                         Summarized data is in log2 space
  
  General options:
  -v --version                  print version and exit
  -h --help                     show full documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 File options

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

=back

=head2 Column index options

=over 4

=item --dindex E<lt>integerE<gt>

=item --lindex E<lt>integerE<gt>

Specify the index numbers of the columns in the data and list 
files, respectively, containing the identifiers to match features. 
If not specified, then the program will attempt to identify  
appropriate matching columns with the same header name. If none 
are specified, the user must select interactively from a list of 
available column names. 

=item --gindex E<lt>integerE<gt>

Specify the group column from the list file. This allows the data 
file to be split into multiple output group files. A column named 
'group' will automatically be identified. A .kgg file will 
automatically use the Cluster column as the group index.

=back

=head2 Output options

=over 4

=item --order [list | data]

Optionally specify the order of features in the output file. Two 
options are available. Specify 'list' to match the order of features 
in the list file. Or specify 'data' to match the order of features 
in the data file. The default is list.

=item --sum

Indicate that the pulled data should be averaged across all 
features at each position, suitable for graphing. A separate text 
file with '_summed' appended to the filename will be written.

=item --sumonly

Indicate that only a summary file should be written, and that the 
pulled data file should be skipped. Useful if you're just after 
the summary for graphing purposes.

=item --starti E<lt>integerE<gt>

When re-summarizing the pulled data, indicate the start column 
index that begins the range of datasets to summarize. Defaults 
to the leftmost column without a standard feature description
name.

=item --stopi E<lt>integerE<gt>

When re-summarizing the pulled data, indicate the stop column
index the ends the range of datasets to summarize. Defaults
to the last or rightmost column.

=item --log

The data is in log2 space. Only necessary when re-summarizing the
pulled data.

=back

=head2 General options

=over 4

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

If the list file has a second group column, then the rows for each group 
will be written to separate files, with the output file name appended with 
the group identifier. Use the gindex option to specify the group column.

The program will also accept a Cluster gene file (with .kgg extension) 
as a list file with group information, where the clusters are the groups. 

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
it under the terms of the Artistic License 2.0.  
