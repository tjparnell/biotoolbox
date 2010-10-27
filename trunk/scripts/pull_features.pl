#!/usr/bin/perl

# A script to pull out a specific subset or list of features or lines of data 
# from a data file. Compare to Excel's VLOOKUP command, only faster.

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
	write_summary_data
);

print "\n This script will pull out specific features from a data file\n";

### Quick help
unless (@ARGV) { # when no command line options are present
	print " Usage for $0
  Required:
  --data <filename>
  --list <filename>
  --out <filename>
  Optional:
  --dindex <integer>
  --lindex <integer>
  --sum
  --starti <integer>
  --stopi <integer>
  --log
  --help\n";
	exit;
}



### Get command line options
my (
	$datafile, 
	$listfile,
	$outfile,
	$sum,
	$startcolumn,
	$stopcolumn,
	$log,
	$help, 
);
my ($data_index, $list_index) = (-0.5, -0.5); # default 'null' values
	# I can't use 0 or a real null that could be interpreted as 0
	# because 0 may be a valid index!
GetOptions( 
	'data=s'     => \$datafile, # the input data file
	'list=s'     => \$listfile, # the list file
	'out=s'      => \$outfile, # the new output file name
	'dindex=i'   => \$data_index, # index to look up in the data file
	'lindex=i'   => \$list_index, # index of look up values in list file
	'sum'        => \$sum, # flag to re-sum the pulled values
	'starti=i'   => \$startcolumn, # index of column to start summarizing
	'stopi=i'    => \$stopcolumn, # index of column to stop summarizing
	'log!'       => \$log, # values are in log, respect log status
	'help'       => \$help, # flag to print help
) or die " a command is wrong! please use --help for more info\n";

if ($help) {
	print_help();
	exit;
}


### Check for required values

unless (defined $datafile) {
	die " no input data file specified!\n";
}

unless (defined $outfile) {
	die " no output data file name given!\n";
}

unless (defined $listfile) {
	die " no list data specified!\n";
}



### Load datafile
print " Loading data file '$datafile'...\n";
my $main_data_ref = load_tim_data_file($datafile);
unless ($main_data_ref) {
	die " No data loaded!\n";
}



### Load the list of specified values

# we will load the list of gene names into a hash
my %request;
print " Collecting lookup values from file '$listfile'...\n";
%request = collect_list_from_file();

unless (%request) {
	die " list of lookup values is empty!\n";
}
#print " found ". keys(%request) . " lookup values\n";




### Determine the dataset index
if ($data_index == -0.5) { # the default 'null' value
	# we must ask the user
	
	# print the headers
	print "\n There are multiple columns in the data file.\n";
	for (my $i = 0; $i < $main_data_ref->{'number_columns'}; $i++) {
		print "   $i\t$main_data_ref->{$i}{name}\n";
	}
	
	# process the answer
	print " Enter the number of the column with the gene names to match   ";
	my $answer = <STDIN>;
	chomp $answer;
	
	# check answer and return
	if (exists $main_data_ref->{$answer}) {
		# answer appears to be a column index
		$data_index = $answer;
	}
	else {
		die " Invalid response!\n";
	}
}




### Pull out the desired features
print " Pulling features...\n";
my @new_data_table; # the new data table with the pulled features

# We will walk through the original data table and lookup the value in the
# request hash. If it exists, copy the feature to the new data table.

# first copy the header row
push @new_data_table, [ @{ $main_data_ref->{'data_table'}->[0] } ];

# then the rest of the rows
foreach my $row ( @{ $main_data_ref->{'data_table'} } ) {
	my $lookup_value = $row->[$data_index];
	if (exists $request{ $lookup_value }) {
		# this row is in the requested list
		# copy to new table
		push @new_data_table, [ @{$row} ];
		# record the find
		$request{ $lookup_value } += 1;
	}
}

# Assign the new data table to the main data structure, replacing the old
$main_data_ref->{'data_table'} = \@new_data_table;
# re-calculate the last row index
$main_data_ref->{'last_row'} = scalar @new_data_table - 1;

# calculate the number of finds
my $found = 0;
my $multiples = 0;
my $notfound = 0;
foreach (keys %request) {
	if ($request{$_} == 0) {
		# not found
		$notfound++;
	}
	elsif ($request{$_} == 1) {
		# found once
		$found++;
	}
	elsif ($request{$_} > 1) {
		# found more than once
		$multiples++;
	}
}



### Write out the new data file
if ($found > 0) {
	# features were found, so let's write the new file
	
	# Print summary
	print "  $found features were found and pulled\n";
	if ($notfound > 0) {
		print "  $notfound features were not found\n";
	}
	if ($multiples > 0) {
		print "  $multiples features were found more than once!!!!\n";
	}
	
	# Write the file
	my $write_results = write_tim_data_file( {
		'data'      => $main_data_ref,
		'filename'  => $outfile,
	} );
	# report write results
	if ($write_results) {
		print "  Wrote new datafile '$outfile'\n";
	}
	else {
		print "  Unable to write datafile '$outfile'!!!\n";
	}
	
	# Summarize the pulled data
	if ($sum) {
		print " Generating final summed data...\n";
		my $sumfile = write_summary_data( {
			'data'         => $main_data_ref,
			'filename'     => $outfile,
			'startcolumn'  => $startcolumn,
			'endcolumn'    => $stopcolumn,
			'log'          => $log,
		} );
		if ($sumfile) {
			print "  Wrote summary file '$sumfile'\n";
		}
		else {
			print "  Unable to write summary file!\n";
		}
	}
	

}
else {
	print "  Nothing found! Nothing to write!\n";
}



sub print_help { # subroutine to print the online help documentation
	print "
 Command line options for $0:

 Given a list of features, this program will pull out those features (rows) 
 from a datafile (compare to Microsoft Excel's VLOOKUP command). The list 
 may be provided as either a separate text file or may be grabbed from the 
 clipbard (note this doesn't work through remote connections!). It will write
 a new data file containing only those features it found. 
 
 The list file may be a simple text file containing the feature names of the 
 features to pull. The first line should be the column header name and should
 be identical to the column header name in the data file used in the lookup.
 If more than one column is present in the list of lookup names, the index 
 will be requested interactively from the user, or be specified as a command
 line argument. The same format applies to data in the clipboard; it simply
 skips the step of writing a separate list file.
 
 The command line arguments include:
 
  Required:
  --data <filename>
  --list <filename> OR --paste
  --out <filename>
  Optional:
  --dindex <integer>
  --lindex <integer>
  --sum
  --starti <integer>
  --stopi <integer>
  --log
  --help
 
 The program writes out a standard tim data format file (see tim_file_helper).

  The command line flags and descriptions:
  Required:
  --data     Specify a tim data formatted input file of genes.
  --out      Specify the output file name. 
  --list     Specify the name of a text file containing the feature names
             or values to look up. The file must contain a column header 
             that matches a column header name in the data file. Multiple 
             columns may be present.
  Optional:
  --dindex   Specify the index number of the column in the data file 
             containing the data to look up and match. Automatically
             determined by matching the header names.
  --lindex   Specify the index number of the column in the list file (or 
             clipboard contents) containing the values to look up and match 
             if more than one column is present. Defaults to interactively 
             asking the user.
  --sum      Indicate that the pulled data should be averaged across all 
             features at each position, suitable for graphing. A separate text 
             file with '_summed' appended to the filename will be written.
  --starti   When re-summarizing the pulled data, indicate the start column 
             index that begins the range of datasets to summarize. Defaults 
             to the leftmost column without a standard feature description
             name.
  --stopi    When re-summarizing the pulled data, indicate the stop column
             index the ends the range of datasets to summarize. Defaults
             to the last or rightmost column.
  --log      The data is in log2 space. Only necessary when re-summarizing the
             pulled data.
  --help     This help text

";
}









### Subroutine to collect list values from a file
sub collect_list_from_file {
	
	# load the file
	my $list_data_ref = load_tim_data_file($listfile);
	unless ($list_data_ref) {
		die " No list file loaded!\n";
	}
	#print " there are $list_data_ref->{last_row} values in list '$listfile'\n";
	
	my @list; # the array to put the final list of features into	
	
	# Check for whether the list file is a Cluster .kgg file
	if ($listfile =~ /\.kgg$/i) {
		# It is a .kgg file
		# This file is a simple two column tab-delimited file generated by 
		# Cluster. The first column contains the gene names. The second 
		# column is the group number.
		
		# we'll use the list_index as the the group number to pull out 
		if ($list_index == -0.5) {
			# we'll need to ask first
			
			# collect the groups and the numbers
			my %groups;
			for (my $row = 1; $row <= $list_data_ref->{'last_row'}; $row++) {
				$groups{ $list_data_ref->{'data_table'}->[$row][1] } += 1;
			}
			
			# ask user
			print "\n These are the group numbers in '$listfile':\n";
			foreach (sort {$a <=> $b} keys %groups) {
				print "   Group $_\thas $groups{$_} genes\n";
			}
			print " Enter the group number to use   ";
			my $answer = <STDIN>;
			chomp $answer;
			if (exists $groups{$answer}) {
				$list_index = $answer;
			}
			else {
				die " unkown response!\n";
			}
		}
		
		# now pull out the list of group genes
		for (my $row = 1; $row <= $list_data_ref->{'last_row'}; $row++) {
			if ($list_data_ref->{'data_table'}->[$row][1] == $list_index) {
				push @list, $list_data_ref->{'data_table'}->[$row][0];
				push @list, 0; # this will become the value in the request hash
			}
		}
	}
	
	else {
		# an ordinary text file
		
		# take the list from the appropriate column
		# first determine which column
		if ($list_data_ref->{'number_columns'} == 1) {
			# there is only one column in the file
			# that must be it!
			$list_index == 0;
		}
		
		elsif ($list_index == -0.5) { 
			# default 'null' value
			# since there are multiple columns, and the index wasn't defined,
			# we'll ask the user which column
			# this assumes that the first line is the column headers
			
			# print the headers
			print "\n There are multiple columns in the list.\n";
			for (my $i = 0; $i < $list_data_ref->{'number_columns'}; $i++) {
				print "   $i\t$list_data_ref->{$i}{name}\n";
			}
			
			# process the answer
			print " Enter the number of the column with the gene names to use   ";
			my $answer = <STDIN>;
			chomp $answer;
			
			# check answer and return
			if (exists $list_data_ref->{$answer}) {
				# answer appears to be a column index
				$list_index = $answer;
			}
			else {
				die " Invalid response!\n";
			}
			
		}
				
		# then walk through and collect the values from the appropriate column
		for (my $row = 1; $row < $list_data_ref->{'last_row'}; $row++) {
			# add the value in the $list_index column to the list array
			push @list, $list_data_ref->{'data_table'}->[$row][$list_index];
			push @list, 0;
		}
		
	}
	
	# return the list
	return @list;
}









