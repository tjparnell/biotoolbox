#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	find_column_index
);
use tim_file_helper qw(
	open_tim_data_file
	write_tim_data_file
	write_summary_data
);
my $VERSION = '1.10';

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




### Open files
# file handles and metadata
my ($list_fh, $list_md) = open_tim_data_file($listfile) or
	die " unable to open list file!\n";

my ($data_fh, $data_md) = open_tim_data_file($datafile) or
	die " unable to open data file!\n";



### Determine indices
identify_indices();



### Load the list of specified values
print " Collecting lookup values from file '$listfile'...\n";
my ($requests, $number_output_files) = collect_request_list();



### Open output files
# this is an array of one or more output data structures
my @out_datas = generate_output_data_structures($number_output_files);



### Pull out the desired features
print " Pulling features...\n";
my ($found_count, $notfound_count) = pull_requested_features();
print "  $found_count features were found and pulled\n";
if ($notfound_count > 0) {
	print "  $notfound_count features were not found\n";
}



### Write the output files
if ($found_count) {
	write_files();
}
else {
	print "  Nothing found! Nothing to write!\n";
}





########################   Subroutines   ###################################

### Identify the indices
sub identify_indices {
	
	# first check whether we have a simple list file
	if ($list_md->{'number_columns'} == 1) {
		# simple one-column file
		# the answer is obvious
		print "  list file only has 1 column, using it\n";
		$list_index = 0;
		return;
	}
	elsif ($listfile =~ /\.kgg$/i) {
		# a cluster gene file
		# the answer is obvious
		print "  using .kgg file as list\n";
		$list_index = 0;
		return;
	}
		
	
	# look for the corresponding list index if data was specified
	if (defined $data_index and !defined $list_index) {
		# we have the data index but need the list index
		
		# get the column header name
		my $lookup = $data_md->{$data_index}{'name'};
		
		# find it in the list
		my $possible = find_column_index($list_md, "^$lookup\$");
		
		# check
		if (defined $possible) {
			# found something
			$list_index = $possible;
			print "  found list column '", $list_md->{$list_index}{'name'}, 
				"', using index $list_index\n";
			return;
		}
	}
	
	# look for the corresponding data index if list was specified
	elsif (!defined $data_index and defined $list_index) {
		# we have the list index but need the data index
		
		# get the column header name
		my $lookup = $list_md->{$list_index}{'name'};
		
		# find it in the list
		my $possible = find_column_index($data_md, "^$lookup\$");
		
		# check
		if (defined $possible) {
			# found something
			$data_index = $possible;
			print "  found data column '", $data_md->{$data_index}{'name'}, 
				"', using index $data_index\n";
			return;
		}
	}
	
	# neither was specified
	elsif (!defined $data_index and !defined $list_index) {
		
		# try some known column identifiers
		foreach my $name (qw(name id transcript gene)) {
			
			# identify possibilities
			my $possible_data_i = find_column_index($data_md, "^$name\$");
			my $possible_list_i = find_column_index($list_md, "^$name\$");
			
			# check if something was found
			if (defined $possible_data_i and defined $possible_list_i) {
				
				# assign
				$data_index = $possible_data_i;
				$list_index = $possible_list_i;
				
				# report
				print " found identical columns\n";
				print "  using list column '", $list_md->{$list_index}{'name'}, 
					"', index $list_index\n";
				print "  using data column '", $data_md->{$data_index}{'name'}, 
					"', index $data_index\n";
				return;
			}
		}
	}
	
	# Give up and ask the user
	foreach ( ['list', $list_md], ['data', $data_md] ) {
		
		# the type and metadata structure
		my $type = $_->[0];
		my $metadata = $_->[1];
		
		# print the headers
		print "\n These are the columns in the $type file.\n";
		for (my $i = 0; $i < $metadata->{'number_columns'}; $i++) {
			print "   $i\t$metadata->{$i}{name}\n";
		}
	
		# process the answer
		print " Enter the $type column index to use    ";
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
	my %number_of_clusters; # the number of clusters (output files)
	
	# Check whether we're working with a Cluster .kgg file
	if ($listfile =~ /\.kgg$/i) {
		while (my $line = $list_fh->getline) {
			chomp $line;
			my @d = split /\t/, $line;
			
			# the gene identifier is the key, the cluster number is the value
			$requests{ $d[0] } = $d[1];
			
			# record the number clusters
			$number_of_clusters{ $d[1] } += 1;
		}
	}
	
	# otherwise a simple text file
	else {
		# collect the identifiers
		while (my $line = $list_fh->getline) {
			chomp $line;
			my @d = split /\t/, $line;
			$requests{ $d[ $list_index ] } = 0;
		}
		
		# we have just one cluster
		$number_of_clusters{0} = 1; 
	}
	
	return (\%requests, scalar(keys %number_of_clusters) );
}



### Subroutine to generate the output data structures
sub generate_output_data_structures {
	my $number = shift;
	
	my @outputs; # array of output structures
		# we will simply duplicate the current data metadata structure 
		# for each output data structure
		# there may be more than one output data structure, primarily with 
		# kgg source files where each group will be put into a new file
	
	# only one output data structure
	if ($number == 1) {
		
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
		push @outputs, \%new_md;
	}
	
	# multiple data structures
	else {
		for my $i (0 .. $number-1) {
			# copy
			my %new_md = %{ $data_md };
			
			# update the file name data
			my $newfile = $outfile;
			my $ext = $data_md->{'extension'};
			$newfile =~ s/$ext\Z//; # remove the extension if present
			$new_md{'filename'} = $newfile . '_g' . $i . $ext;
			$new_md{'basename'} = undef;
			$new_md{'path'} = undef;
			
			# add the data table
			$new_md{'data_table'} = [];
			# add the column headers
			push @{ $new_md{'data_table'} }, [ @{ $data_md->{'column_names'} } ];
			
			# store the structure
			push @outputs, \%new_md;
		}
	}
	
	return @outputs;
}



### Subroutine to pull the requested features
sub pull_requested_features {
	
	# we will now walk through the data file and pull out those lines 
	# which match the feature identifier
	my $found = 0;
	my $notfound = 0;
	while (my $line = $data_fh->getline) {
		chomp $line;
		my @d = split /\t/, $line;
		
		# check if the identifier exists in our requests hash
		if (exists $requests->{ $d[$data_index] } ) {
			# found
			
			# determine which output data table to store it
			my $out = $out_datas[ $requests->{ $d[$data_index] } ];
			
			# record the line data
			push @{ $out->{'data_table'} }, \@d;
			$found++;
		}
		else {
			# not found
			$notfound++;
		}
	}
	
	# re-calculate the last row index
	foreach my $out (@out_datas) {
		$out->{'last_row'} = scalar(@{ $out->{'data_table'} }) - 1;
	}
	
	return ($found, $notfound);
}



### Subroutine to write the output files
sub write_files {
	# Write the files
	
	foreach my $out (@out_datas) {
		
		# check if there's any data in here
		next unless $out->{'last_row'} > 0; # nothing but the header
		
		# write the file
		my $write_results = write_tim_data_file(
			'data'      => $out,
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
				'data'         => $out,
				'filename'     => $out->{'filename'}, # must provide
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

Given a list of requested feature identifiers, this program will pull 
out those features (rows) from a datafile and write a new file. This 
program compares in function to a popular spreadsheet VLOOKUP command. 
The list is provided as a separate text file, either as a single column 
file or a multi-column tab-delimited from which one column is selected. 
All rows from the source data file that match an identifier in the list 
will be written to the new file. 

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
