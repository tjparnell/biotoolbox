#!/usr/bin/perl

# A script to graph histogram plots for one or two microarray data sets

use strict;
use Getopt::Long;
use GD::Graph::lines;
use GD::Graph::bars;
use Statistics::Lite qw(mean max);
use Statistics::Descriptive;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper;

print "\n This script will plot histograms of value frequencies\n\n";

### Get command line options
my (
	$infile, 
	$out,
	$binnumber,
	$binsize,
	$lines,
	$all,
	$max,
	$help
);
my @columns; # an array of the columns (datasets) within in the file to plot
			 # should be 0-based, may have two datasets comma delimited
my $directory = 'graphs'; # default value unless otherwise specified
my $start = 0; # the starting value to calculate the bins
GetOptions( 
	'in=s'       => \$infile, # the input file
	'col=s'      => \@columns, # columns (datasets) to plot
	'out=s'      => \$out, # output file name
	'bins=i'     => \$binnumber, # the number of bins to put the data into
	'binsize=f'  => \$binsize, # the size of each bin
	'start=i'    => \$start, # the starting value to calculate the bins
	'lines'      => \$lines, # indicate whether graph should be a linegraph
	'all'        => \$all, # flag to plot all data sets pairwise
	'max=i'      => \$max, # maximum value for x-axis
	'dir=s'      => \$directory, # optional name of the graph directory
	'help'       => \$help, # flag to print help
);

if ($help) {
	print_help();
	exit;
}

unless ($infile) {die " Missing input file!\n Please use --help for more information\n"}
unless ($binnumber and $binsize) {die " Missing bin number and bin size!\n Please use --help for more information\n"}

unless ($max) {
	# default value is calculated
	$max = $start + ($binnumber * $binsize);
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
		$name =~ /stop/i
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
unless (-e "$directory") {
	mkdir $directory or die "Can't create directory $directory\n";
}




### Get the list of datasets to pairwise compare

# A list of datasets was provided as a file
if (@columns) {
	graph_designated_datasets();
} 

# Interactive session
else {
	graph_datasets_interactively();
}

	
##### Subroutines ##########
## subroutine to print the online help documentation
sub print_help {
	print "
 Command line options for $0

 Put more detailed info in here. Sometime soon. You lazy pickle.
  
  The command line flags and descriptions:
  --in       Specify the file name for a gene data list containing the
             microarray data. Generate using the script 'get_datasets.pl'
  --col      Specify the column number(s) corresponding to the dataset(s) in
             the file to graph. Number is 0-based index. Use repeatedly for 
             each dataset to graph. To graph two datasets together on one plot, 
             use a comma between the numbers. If the column is not specified, 
             then the program willask interactively which datasets to plot.
  --bins     Specify the number of bins the data will be grouped into
  --binsize  Specify the size of each bin. Note that any value above 
             (bins * binsize) will not be counted.
  --start    Optionally indicate the starting value that the bins should be
             generated. Default is 0.
  --lines    Optionally specify a line graph instead of default bar graph    
  --max      Optionally specify the maximum X-axis value. Default is calculated
             from (bins * binsize).
  --out      Optionally specify the output base filename
  --dir      Specify an optional name for the output subdirectory name 
  --help     This help text

";
}



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
	foreach (@columns) {
		# we're plotting two datasets
		if (/,/) { 
			s/\s*//g; 
			my ($one, $two) = split /,/;
			if ( 
				(exists $dataset_by_id{$one}) and 
				(exists $dataset_by_id{$two}) 
			) {
				graph_two("$dataset_by_id{$one}", "$dataset_by_id{$two}");
			} 
			else {
				print " Column numbers $_ is not valid\n";
			}
		} 
		
		# we're plotting only one dataset
		else { 
			if (exists $dataset_by_id{$_} ) {
				graph_one($dataset_by_id{$_});
			} 
			else {
				print " Column number $_ is not valid\n";
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
		if ($answer =~ /,/) { # two datasets are requested
			$answer =~ s/\s*//g; 
			my ($one, $two) = split /,/, $answer;
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
		print "\n Enter the next set of datasets or push return to exit  ";
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
	
	# Add the maximum value to @bins to count everyone
	my $max_value = max(@values);
	push @bins, $max_value;
	
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
		graph_this_as_lines($name, "", $title, \@data);
	} 
	else {
		graph_this_as_bars($name, "", $title, \@data);
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
	
	# Add the maximum value to @bins to count everyone
	my $max1 = max(@values1);
	my $max2 = max(@values2);
	if ($max1 >= $max2) {
		push @bins, $max1;
	} 
	else {
		push @bins, $max2;
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
		x_label         => 'values',
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
		$filename = $filename . "_$name" . ".png";
	} 
	else {
		$filename = "distribution_$name" . '.png';
	}
	open IMAGE, ">$directory/$filename" or die " Can't open output file!\n";
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	print "wrote line graph $filename in directory $directory\n";
}

# Write a bar graph
sub graph_this_as_bars { 	
	my ($name1, $name2, $title, $dataref) = @_;
	my $name; # generic name
	
	my $graph = GD::Graph::bars->new(600,400);
	$graph->set(
		title			=> $title,
		x_label         => 'values',
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
		$filename = $filename . "_$name" . ".png";
	} 
	else {
		$filename = "distribution_$name" . '.png';
	}
	open IMAGE, ">$directory/$filename" or die " Can't open output file!\n";
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	print "wrote bar graph $filename in directory $directory\n";
}




