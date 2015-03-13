#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use Statistics::Descriptive;
use Bio::ToolBox::data_helper qw(parse_list);
use Bio::ToolBox::file_helper qw(load_tim_data_file);
my $gd_ok;
eval {
	require GD::Graph::lines; 
	require GD::Graph::bars; 
	$gd_ok = 1;
};
my $VERSION = '1.15';

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
	$y_max,
	$y_ticks,
	$x_skip,
	$x_offset,
	$x_format,
	$directory,
	$help,
	$print_version,
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
	'ymax=i'      => \$y_max, # maximum value on y axis
	'yticks=i'    => \$y_ticks, # number of ticks on y axis
	'skip=i'      => \$x_skip, # number of ticks to skip on x axis
	'offset=i'    => \$x_offset, # skip number of x axis ticks before labeling
	'format=i'    => \$x_format, # format decimal numbers of x axis
	'dir=s'       => \$directory, # optional name of the graph directory
	'help'        => \$help, # flag to print help
	'version'     => \$print_version, # print the version
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
	print " Biotoolbox script graph_histogram.pl, version $VERSION\n\n";
	exit;
}




### Check required and default values
unless ($gd_ok) {
	die "Module GD::Graph must be installed to run this script.\n";
}

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
unless (defined $x_skip) {
	$x_skip = 4;
}
unless (defined $x_format) {
	$x_format = 0;
	# set it to equal the number of decimals in binsize
	if ($binsize =~ /\.(\d+)$/) {
		$x_format = length $1;
	}
}
unless (defined $x_offset) {
	$x_offset = 0;
}
unless (defined $y_ticks) {
	$y_ticks = 4;
}







####### Main ###########

### Load the file
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
		$name =~ /^(?:name|id|alias|chromosome|start|stop|end|strand|type|class)$/i
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
my $format = $binsize =~ /\.(\d+)$/ ? length $1 : $x_format;
for my $i (1 .. $binnumber) {
	push @bins, sprintf "%.$format" . "f", $start + ($i * $binsize);
}

# Prepare output directory
unless ($directory) {
	$directory = $main_data_ref->{'path'} . $main_data_ref->{'basename'} . '_graphs';
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
		unless (
			$y eq '.' or 
			$y < $start or
			$y > $max
		) {
			push @values, $y; # put into the values array
		}
	}
	#print "  found " . scalar @values . " useable values\n";
	
	# Determine the data frequency
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@values);
	my %frequency = $stat->frequency_distribution(\@bins);
	my @yvalue; # an array of arrays for the graph data
	foreach (sort {$a <=> $b} keys %frequency) {
		push @yvalue, $frequency{$_};
		# print "   x $_ , y $frequency{$_}\n"; # print values to check
	}
	my $string = '%.' . $x_format . 'f';
	my @xlabels = map { sprintf $string, $_ } @bins;
	my @data = ( [@xlabels], [@yvalue] );
	
	# Now (finally) prepare the graph
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
		unless (
			$value1 eq '.' or
			$value1 < $start or
			$value1 > $max
		) {
			push @values1, $value1; # put into the values array
		}
		unless (
			$value2 eq '.' or
			$value2 < $start or
			$value2 > $max
		) {
			push @values2, $value2; # put into the values array
		}
	}
	
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
	
	my $string = '%.' . $x_format . 'f';
	my @xlabels = map { sprintf $string, $_ } @bins;
	my @data = ( \@xlabels, \@yvalue1, \@yvalue2 );
	
	# Now (finally) prepare the graph
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
		x_tick_number	=> scalar @bins,
	) or warn $graph->error;
	
	# set axes
	set_graph_axes($graph);
	
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
		bar_spacing     => 2,
	) or warn $graph->error;
	
	# set axes
	set_graph_axes($graph);
	
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


## Set some generic axis options regardless of graph type
sub set_graph_axes {
	my $graph = shift;
	
	# set ticks and label number format
	$graph->set(
		x_label_skip    => $x_skip,
		x_tick_offset   => $x_offset,
		x_number_format	=> '%.' . $x_format . 'f', # "%.2f"
		y_tick_number	=> $y_ticks,
		y_long_ticks	=> 1,
		transparent		=> 0,
	) or warn $graph->error;
	
	# set y axis maximum
	if ($y_max) {
		$graph->set(y_max_value => $y_max) or warn $graph->error;
	}
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

A script to graph a histogram (bar or line) of one or more datasets.

=head1 SYNOPSIS

graph_histogram.pl --bins <integer> --size <number> <filename> 

graph_histogram.pl --bins <integer> --max <number> <filename> 

graph_histogram.pl --size <number> --max <number> <filename> 
   
  Options:
  --in <filename>
  --index <column_index>
  --bins <integer>
  --size <number>
  --min <number>
  --max <number>
  --ymax <integer>
  --yticks <integer>
  --skip <integer>
  --offset <integer>
  --format <integer>
  --lines
  --out <base_filename>
  --dir <output_directory>
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. The 
first row should be column headers. Text files generated by other 
B<BioToolBox> scripts are acceptable. Files may be gzipped compressed.

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

=item --ymax <integer>

Specify the maximum Y axis value. The default is automatically determined.

=item --yticks <integer>

Specify explicitly the number of major ticks for the Y axes. 
The default is 4.

=item --skip <integer>

Specify the ordinal number of X axis major ticks to label. This 
avoids overlapping labels. The default is 4 (every 4th tick is labeled).

=item --offset <integer>

Specify the number of X axis ticks to skip at the beginning before starting 
to label them. This may help in adjusting the look of the graph. The 
default is 0.

=item --format <integer>

Specify the number of decimal places the X axis labels should be formatted. 
The default is the number of decimal places in the bin size parameter.

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

=item --version

Print the version number.

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
