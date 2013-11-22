#!/usr/bin/perl

# A script to graph histogram plots for one or two microarray data sets

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use GD::Graph::lines;
use GD::Graph::bars;
use Statistics::Lite qw(mean max);
use Statistics::Descriptive;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	parse_list
);
use tim_file_helper qw(
	load_tim_data_file
);

print "\n This script will plot histograms of value frequencies\n\n";

### Quick help
unless (@ARGV) { # when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Get command line options
my (
	$infile, 
	$dataset,
	$out,
	$binnumber,
	$binsize,
	$lines,
	$start,
	$max,
	$directory,
	$help
);
GetOptions( 
	'in=s'        => \$infile, # the input file
	'index=s'     => \$dataset, # columns (datasets) to plot
	'out=s'       => \$out, # output file name
	'bins=i'      => \$binnumber, # the number of bins to put the data into
	'size=f'      => \$binsize, # the size of each bin
	'min=f'       => \$start, # the starting value to calculate the bins
	'max=f'       => \$max, # maximum value for x-axis
	'lines!'      => \$lines, # indicate whether graph should be a linegraph
	'dir=s'       => \$directory, # optional name of the graph directory
	'help'        => \$help, # flag to print help
);

if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}



### Check required and default values
unless ($infile) {
	$infile = shift @ARGV or
		die " no input file! use --help for more information\n";
}

unless (defined $start) {
	$start = 0; # default start of 0
}
unless (defined $binsize) {
	if (defined $binnumber and defined $max) {
		$binsize = ($max - $start) / $binnumber;
	}
	else {
		die " need to specify at least two of bins, binsize, or max! see help\n";
	}
}
unless (defined $binnumber) {
	if (defined $binsize and defined $max) {
		$binnumber = int( ($max - $start) / $binsize);
	}
	else {
		die " need to specify at least two of bins, binsize, or max! see help\n";
	}
}
unless (defined $max) {
	if (defined $binsize and defined $binnumber) {
		$max = $start + ($binnumber * $binsize);
	}
	else {
		die " need to specify at least two of bins, binsize, or max! see help\n";
	}
}
		








####### Main ###########

### Load the file
# load the file using subroutine from tim_db_helper.pm
print " Loading data from file $infile....\n";
my $main_data_ref = load_tim_data_file($infile);
unless ($main_data_ref) {
	die " No data loaded!\n";
}
my $data_table_ref = $main_data_ref->{'data_table'};

# load the dataset names into hashes
my (%dataset_by_name, %dataset_by_id); # hashes for name and id
for (my $i = 0; $i < $main_data_ref->{'number_columns'}; $i++) {
	my $name = $main_data_ref->{$i}{'name'};
	if (
		# check column header names for gene or window attribute information
		$name =~ /name/i or 
		$name =~ /class/i or
		$name =~ /alias/i or
		$name =~ /^chr/i or
		$name =~ /start/i or
		$name =~ /stop/i or 
		$name =~ /type/i
	) { 
		# skip on to the next header
		next; 
	} 
	else { 
		# record the data set name
		$dataset_by_id{$i} = $name;
		$dataset_by_name{$name} = $i;
	}
}	

# determine the bins for the frequency distribution
my @bins; # an array for the bins
for (my $i = ($binsize + $start); $i < $max; $i += $binsize) {
	push @bins, $i;
	#print "$i ";
}

# Prepare output directory
unless ($directory) {
	$directory = $main_data_ref->{'basename'} . '_graphs';
}
unless (-e "$directory") {
	mkdir $directory or die "Can't create directory $directory\n";
}




### Get the list of datasets to pairwise compare

# A list of datasets was provided as a file
if ($dataset) {
	my @list = parse_list($dataset);
	graph_designated_datasets(@list);
} 

# Interactive session
else {
	graph_datasets_interactively();
}

print " Finished!\n";

	





########################### Subroutines #######################################




## subroutine to graph all the data sets in the loaded data file
sub graph_all_datasets {
	
	# put the dataset IDs into sorted array
	my @list = sort {$a <=> $b} keys %dataset_by_id; 
	
	foreach (@list) {
		graph_one($_);
	}
}


## subroutine to graph all the data sets in the loaded data file
sub graph_designated_datasets {
	my @columns = @_;
	foreach (@columns) {
		# we're plotting two datasets
		if (/&/) { 
			my ($one, $two) = split /&/;
			if ( 
				(exists $dataset_by_id{$one}) and 
				(exists $dataset_by_id{$two}) 
			) {
				graph_two("$dataset_by_id{$one}", "$dataset_by_id{$two}");
			} 
			else {
				print " Column numbers '$_' is not valid\n";
			}
		} 
		
		# we're plotting only one dataset
		else { 
			if (exists $dataset_by_id{$_} ) {
				graph_one($dataset_by_id{$_});
			} 
			else {
				print " Column number '$_' is not valid\n";
			}
		}
	}
}


## Subroutine to interact with user and ask for data set pairs to graph sequentially
sub graph_datasets_interactively {
	
	# present the dataset name list
	print " These are the data sets in the file $infile\n";
	foreach (sort {$a <=> $b} keys %dataset_by_id) {
		print "   $_\t$dataset_by_id{$_}\n";
	}
	
	# collect user answer
	print " Enter the numbers dataset to graph  ";
	my $answer = <STDIN>; # get user answer
	chomp $answer;
	
	# this loop will keep going until no dataset (undefined) is returned
	while ($answer) {
		if ($answer eq 'q') {
			last;
		}
		if ($answer =~ /,|&/) { # two datasets are requested
			$answer =~ s/\s*//g; 
			my ($one, $two) = split /[&,]/, $answer;
			if ( 
				(exists $dataset_by_id{$one}) and 
				(exists $dataset_by_id{$two}) 
			) {
				graph_two("$dataset_by_id{$one}", "$dataset_by_id{$two}");
			} 
			else {
				print " One of these numbers is not valid\n";
			}
			
		} else { # only one dataset is requested
			if (exists $dataset_by_id{$answer} ) {
				graph_one($dataset_by_id{$answer});
			} 
			else {
				print " One of these numbers is not valid\n";
			}
		}
		print "\n Enter the next set of datasets or enter 'q' or nothing to exit  ";
		$answer = <STDIN>;
		chomp $answer;
	}
}


## Graph one dataset
sub graph_one {
	
	# Retrieve the values for plotting
	my $name = shift;
	print " Preparing graph for $name....\n";
	my $index = $dataset_by_name{$name};
	my @values;
	for my $i (1..$main_data_ref->{'last_row'}) { # walk through the data file
		my $y = $data_table_ref->[$i][$index];
		# only take numerical data
		# must have a numeric value from both datasets, otherwise skip
		if ($y ne '.') {
			push @values, $y; # put into the values array
		}
	}
	#print "  found " . scalar @values . " useable values\n";
	
# 	# Add the maximum value to @bins to ensure we count everyone
# 	if (max(@values) > max(@bins)) {
# 		push @bins, max(@values);
# 	}
	
	# Determine the data frequency
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@values);
	my %frequency = $stat->frequency_distribution(\@bins);
	my @yvalue; # an array of arrays for the graph data
	foreach (sort {$a <=> $b} keys %frequency) {
		push @yvalue, $frequency{$_};
		#print "   x $_ , y $frequency{$_}\n"; # print values to check
	}
	my @data = ( [@bins], [@yvalue] );
	
	# Now (finally) prepare the graph
	pop @bins; # remove the maximum that we added
	my $title = "Distribution of $out $name";
	if ($lines) {
		graph_this_as_lines($name, undef, $title, \@data);
	} 
	else {
		graph_this_as_bars($name, undef, $title, \@data);
	}
}	




## Graph two datasets
sub graph_two {
	
	# Retrieve the values for plotting
	my ($name1, $name2) = @_;
	print " Preparing graph for $name1 and $name2....\n";
	my $index1 = $dataset_by_name{$name1};
	my $index2 = $dataset_by_name{$name2};
	my (@values1, @values2);
	for my $i (1..$main_data_ref->{'last_row'}) { 
		# walk through the data file
		my $value1 = $data_table_ref->[$i][$index1];
		my $value2 = $data_table_ref->[$i][$index2];
		# only take numerical data
		# must have a numeric value from both datasets, otherwise skip
		if ($value1 ne '.') {
			push @values1, $value1; # put into the values array
		}
		if ($value2 ne '.') {
			push @values2, $value2; # put into the values array
		}
	}
	
# 	# Add the maximum value to @bins to count everyone
# 	my $max1 = max(@values1);
# 	my $max2 = max(@values2);
# 	if ($max1 >= $max2) {
# 		push @bins, $max1;
# 	} 
# 	else {
# 		push @bins, $max2;
# 	}
	
	# Determine the data frequency
	# we have to first determine which dataset has the biggest max value
	# this will determine which dataset we will use first
	# we will then use the bins specified from the this first dataset in 
	# processing the second dataset
	my (@yvalue1, @yvalue2); # an array of arrays for the graph data
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@values1);
	my %frequency = $stat->frequency_distribution(\@bins);
	foreach (sort {$a <=> $b} keys %frequency) {
		push @yvalue1, $frequency{$_};
	}
	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@values2);
	%frequency = $stat->frequency_distribution(\@bins);
	foreach (sort {$a <=> $b} keys %frequency) {
		push @yvalue2, $frequency{$_};
	}
	
	my @data = ( \@bins, \@yvalue1, \@yvalue2 );
	
	# Now (finally) prepare the graph
	pop @bins; # remove the maximum that we added
	my $title = "Distributions of $out $name1 & $name2";
	if ($lines) {
		graph_this_as_lines($name1, $name2, $title, \@data);
	} 
	else {
		graph_this_as_bars($name1, $name2, $title, \@data);
	}
}



# Write a line graph
sub graph_this_as_lines { 	
	my ($name1, $name2, $title, $dataref) = @_;
	my $name; # generic name
	
	my $graph = GD::Graph::lines->new(600,400);
	$graph->set(
		title			=> $title,
		x_label         => $name2 ? 'values' : "$name1 value", 
		y_label         => 'number',
		x_max_value     => $max,
		x_min_value     => $start,
		x_number_format	=> "%.1f",
		x_tick_number	=> scalar @bins,
		x_label_skip    => 4,
		y_long_ticks	=> 1,
		transparent		=> 0,
		
	) or warn $graph->error;
	
	# options for two datasets
	if ($name2) {
		$graph->set_legend($name1, $name2);
		$graph->set(legend_placement => 'BC'); # bottom center
		$name = "$name1\_&_$name2";
	} 
	else {
		$name = $name1;
	}
	
	# Write the graph file
	my $gd = $graph->plot($dataref) or warn $graph->error;
	my $filename;
	if ($out) {
		$filename = $out;
		$filename =~ s/\.\w*$//;
		$filename = $filename . "_$name";
	} 
	else {
		$filename = "distribution_$name";
	}
	$filename = File::Spec->catfile($directory, $filename);
	$filename = check_file_uniqueness($filename, 'png');
	
	open IMAGE, ">$filename" or die " Can't open output '$filename'!\n";
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	print "wrote line graph '$filename'\n";
}





# Write a bar graph
sub graph_this_as_bars { 	
	my ($name1, $name2, $title, $dataref) = @_;
	my $name; # generic name
	
	my $graph = GD::Graph::bars->new(600,400);
	$graph->set(
		title			=> $title,
		x_label         => $name2 ? 'values' : "$name1 value",
		y_label         => 'number',
		x_number_format	=> "%.1f",
		#x_tick_number	=> scalar @bins,
		bar_spacing     => 2,
		x_label_skip    => 4,
		y_long_ticks	=> 1,
		transparent		=> 0,
		
	) or warn $graph->error;
	
	# options for two datasets
	if ($name2) {
		$graph->set_legend($name1, $name2);
		$graph->set(legend_placement => 'BC'); # bottom center
		$name = "$name1\_&_$name2";
	} 
	else {
		$name = $name1;
	}
	
	# Write the graph file
	my $gd = $graph->plot($dataref) or warn $graph->error;
	my $filename;
	if ($out) {
		$filename = $out;
		$filename =~ s/\.\w*$//;
		$filename = $filename . "_$name";
	} 
	else {
		$filename = "distribution_$name";
	}
	$filename = File::Spec->catfile($directory, $filename);
	$filename = check_file_uniqueness($filename, 'png');
	
	open IMAGE, ">$filename" or die " Can't open output file '$filename'!\n";
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	print "wrote bar graph '$filename'\n";
}



## Make a unique filename
sub check_file_uniqueness {
	my ($filename, $extension) = @_;
	my $number = 1;
	
	# check whether the file is unique
	if (-e "$filename\.$extension") {
		# file already exists, need to make it unique
		my $test = $filename . '_' . $number . '.' . $extension;
		while (-e $test) {
			$number++;
			$test = $filename . '_' . $number . '.' . $extension;
		}
		# found a unique file name
		return $test;
	}
	else {
		# filename is good
		return "$filename\.$extension";
	}
}




__END__

=head1 NAME

graph_histogram.pl

A script to graph a histogram of a dataset of values

=head1 SYNOPSIS

graph_histogram.pl --bins <integer> --size <number> <filename> 
   
   --in <filename>
   --index <column_index>
   --bins <integer>
   --size <number>
   --min <number>
   --max <number>
   --lines
   --out <base_filename>
   --dir <output_directory>
   --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a previously generated feature dataset.
The tim data format is preferable, although any other tab-delimited text 
data formats may be usable. See the file description in C<tim_db_helper.pm>.

=item --index <column_index>

Specify the column number(s) corresponding to the dataset(s) in
the file to graph. Number is 0-based index. Each dataset should be 
demarcated by a comma. A range of indices may also be specified using 
a dash to demarcate the beginning and end of the inclusive range. 
Two datasets may also be graphed together; these indices should be joined 
with an ampersand. For example, "2,4-6,5&6" will individually graph 
datasets 2, 4, 5, 6, and a combination 5 and 6 graph.

If no dataset indices are specified, then they may be chosen 
interactively from a list.

=item --bins <integer>

Specify the number of bins or partitions into which the data will be 
grouped. This argument is optional if --max and --size are provided.

=item --size <number>

Specify the size of each bin or partition. A decimal number may be 
provided. This argument is optional if --bins and --max are provided.

=item --min <number>

Optionally indicate the minimum value of the bins. When generating 
the list of bins, this is used as the starting value. Default is 0. 
A negative number may be provided using the format --min=-1. 

=item --max <number>

Specify the maximum bin value. This argument is optional if --bins 
and --size are provided.

=item --lines

Optionally specify a line graph to be generated instead of the 
default vertical bar graph.   

=item --out

Optionally specify the output filename prefix. The default value is 
"distribution_".

=item --dir

Optionally specify the name of the target directory to place the 
graphs. The default value is the basename of the input file 
appended with "_graphs".

=item --help

Print this help documenation

=back

=head1 DESCRIPTION

This program will generate PNG graphic files representing the histogram 
of the values in one or two datasets. The size of each bin or partition 
must be provided, as well as either the number of bins or the maximum 
bin value. The resulting files are written to a subdirectory named after 
the input file. The files are named after the dataset name (column 
header) with a prefix. 

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















