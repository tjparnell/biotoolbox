#!/usr/bin/perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
use Statistics::Descriptive;
use Bio::ToolBox::Data;
use Bio::ToolBox::utility;
my $gd_ok;
eval {
	require GD;
	require GD::Graph::lines; # for the line type graph
	require GD::Graph::mixed; # for the scatter plot
	$gd_ok = 1;
};
my $gd_smooth;
eval {
	require GD::Graph::smoothlines; 
	$gd_smooth = 1;
};
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};
my $VERSION = '1.20';

print "\n This script will graph correlation plots for two data sets\n\n";

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
	$infile, 
	$type, 
	$index,
	$all, 
	$norm, 
	$moving_average,
	$min, 
	$max, 
	$x_min,
	$x_max,
	$y_min,
	$y_max,
	$x_ticks,
	$y_ticks,
	$ticks,
	$x_format,
	$y_format,
	$format,
	$places,
	$directory,
	$out,
	$numbers,
	$cpu,
	$help,
	$print_version,
);
my @pairs; # an array of pairs
GetOptions( 
	'in=s'      => \$infile, # the input file
	'type=s'    => \$type, # the type of graph to generate
	'index=s'   => \$index, # a list of x,y pairs to plot
	'pair=s'    => \@pairs, # the x,y pairs to plot
	'all'       => \$all, # flag to plot all data sets pairwise
	'norm!'     => \$norm, # flag to skip percentile normalizing
	'ma=s'      => \$moving_average, # window, setp values for moving average
	'min=f'     => \$min, # mininum axis coordinate
	'max=f'     => \$max, # maximum axis coordinate
	'xmin=f'    => \$x_min, # minimum x axis coordinate
	'ymin=f'    => \$y_min, # minimum y axis coordinate
	'xmax=f'    => \$x_max, # maximum x axis coordinate
	'ymax=f'    => \$y_max, # maximum y axis coordinate
	'ticks=i'   => \$ticks, # the number of major axes ticks
	'xticks=i'  => \$x_ticks, # the number o major x axis ticks
	'yticks=i'  => \$y_ticks, # the number o major y axis ticks
	'format=i'  => \$format, # number of places to format tick labels
	'xformat=i' => \$x_format, # number of places to format x axis tick labels
	'yformat=i' => \$y_format, # number of places to format y axis tick labels
	'dir=s'     => \$directory, # optional name of the graph directory
	'out=s'     => \$out, # output file name
	'numbers'   => \$numbers, # print the graph numbers in addition to the graph
	'cpu=i'     => \$cpu, # number of parallel threads to execute
	'help'      => \$help, # flag to print help
	'version'   => \$print_version, # print the version
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
	print " Biotoolbox script graph_data.pl, version $VERSION\n\n";
	exit;
}



### check requirements
unless ($gd_ok) {
	die <<GD_WARNING
Modules GD and GD::Graph failed to load, either because they are not installed or they 
are missing an external dependency (libgd). Please install these to run this script.
GD_WARNING
}

unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		die " Missing input file!\n";
	}
}


### set defaults
if ($type) {
	unless ($type =~ /scatter|line|smooth/i ) {
		die " Unknown graph type '$type'!\n Please use --help for more information\n";
	}
}
else {
	$type = 'scatter';
	print " Using default graph type of 'scatter'\n";
}
if ($type =~ /smooth/i) {
	die "Perl module GD::Graph::smoothlines must be installed to graph smooth plots\n"
		unless $gd_smooth;
}
unless (defined $norm) {
	# default is no normalization (percent rank)
	$norm = 0;
}
if (defined $min) {
	# assign general minimum value to specific axes
	unless (defined $x_min) {
		$x_min = $min;
	}
	unless (defined $y_min) {
		$y_min = $min;
	}
}
if (defined $max) {
	# assign general maximum value to specific axes
	unless (defined $x_max) {
		$x_max = $max;
	}
	unless (defined $y_max) {
		$y_max = $max;
	}
}
unless (defined $x_ticks) {
	$x_ticks = defined $ticks ? $ticks : 4;
}
unless (defined $y_ticks) {
	$y_ticks = defined $ticks ? $ticks : 4;
}
unless (defined $x_format) {
	$x_format = defined $format ? $format : 0;
}
unless (defined $y_format) {
	$y_format = defined $format ? $format : 0;
}
if ($moving_average) {
	unless ( $moving_average =~ /^(\d+),(\d+)$/ and $1 >= $2 ) {
		die " moving average values is not (well) defined!\n Use --help for more information\n";
	}
}
if ($parallel) {
	# conservatively enable 2 cores
	$cpu ||= 2;
}
else {
	# disable cores
	print " disabling parallel CPU execution, no support present\n" if $cpu;
	$cpu = 0;
}



####### Main ###########

### Prepare global variables and set up for execution

print " Loading data from file $infile....\n";
my $Data = Bio::ToolBox::Data->new(file => $infile) or
	die " No data loaded!\n";

# load the dataset names into hashes
my %dataset_by_id; # hashes for name and id
my $i = 0;
foreach my $name ($Data->list_columns) {
	
	# check column header names for gene or window attribute information
	# these won't be used for graph generation, so we'll skip them
	next if $name =~ /^(?:name|id|class|type|alias|probe|chr|
		chromo|chromosome|seq|sequence|refseq|contig|scaffold|start|stop|end|mid|
		midpoint|strand|primary_id)$/xi;
	
	# record the data set name
	$dataset_by_id{$i} = $name;
	$i++;
}	


# Prepare output directory
unless ($directory) {
	# directory may be specified from a command line argument, otherwise
	# generate default directory from input file name
	$directory = $Data->path . $Data->basename . '_graphs';
}
unless (-e "$directory") {
	mkdir $directory or die "Can't create directory '$directory'\n";
}


# Set up correlation values Data object
my $statfile = File::Spec->catfile($directory, $Data->basename . "_stats.txt");
my $stat_Data;
if (-e $statfile) {
	# if the file exists, append to it
	$stat_Data = Bio::ToolBox::Data->new(file => $statfile) or 
		die " unable to write '$statfile'!\n";
} 
else {
	# open a new file
	$stat_Data->Bio::ToolBox::Data->new(
		feature => 'correlations',
		columns => [qw(
			SourceFile
			GraphFile
			GraphType
			Normalized
			X_Data
			Y_Data
			Intercept
			Slope
			Pearson_coefficient
			R^2
		)],
	); 
}






### Get the list of datasets to pairwise compare

if ($all) {
	# All datasets in the input file will be compared to each other
	print " Comparing all data set pairwise combinations....\n";
	graph_all_datasets();
}
elsif ($index or @pairs) {
	# A list of dataset pairs was provided upon execution
	
	# provided as a string in the index option
	if ($index) {
		push @pairs, split /,/, $index;
	}
	
	# walk through each 
	my @to_do;
	foreach my $pair (@pairs) {
		my ($x, $y) = split /,|&/, $pair;
		if ( 
			(exists $dataset_by_id{$x}) and 
			(exists $dataset_by_id{$y}) 
		) {
			push @to_do, [$x, $y];
		} 
		else {
			print " One of these numbers ($x, $y) is not valid\n";
		}
	}
	
	# graph the pairs
	graph_provided_datasets(@to_do);
}
else {
	# Interactive session
	graph_datasets_interactively();
}



# Write out the linear regression statistics
$stat_Data->write_file(
	filename  => $statfile,
	gz        => 0,
	simple    => 1, # no metadata
);


### The end of the main program	
exit;



##### Subroutines ##########



## subroutine to graph all the data sets in the loaded data file
sub graph_all_datasets {
	
	# a list of all the datasets to work with
	my @list;
	if ($index) {
		# user provided a list to work with
		@list = parse_list($index);
	}
	else {
		# put the dataset IDs into sorted array
		@list = sort {$a <=> $b} keys %dataset_by_id; 
	}
	
	# Now process through the list 
	# the first dataset will be x, and the subsequent will by y
	# all of the ys will be done for a given x, then x is incremented and the ys repeated
	# any y less than the current x will be skipped to avoid repeats
	my @to_do;
	for (my $x = 0; $x < (scalar @list - 1); $x++) {
		for (my $y = 1; $y < (scalar @list); $y++) {
			if ($y <= $x) {next}; # skip pairs we've already done
			push @to_do, [$list[$x], $list[$y]];
		}
	}
	
	# graph the datasets
	graph_provided_datasets(@to_do);
}


## Graph the dataset pairs provided by the user
sub graph_provided_datasets {
	my @to_do = @_;
	# the to_do list is an array of [x,y] arrays
	
	
	# We can graph either in parallel or serially
	
	# Parallel Execution
	if ($cpu > 1 and scalar(@to_do) >= $cpu) {
		
		# prepare ForkManager
		print " Forking into $cpu children for parallel graph generation\n";
		my $pm = Parallel::ForkManager->new($cpu);
		$pm->run_on_finish( sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result) = @_;
			$stat_Data->add_row($result) if $exit_code == 1;
		});
		
		# run in parallel
		foreach (@to_do) {
			$pm->start and next;
			
			# in child
			my $result = graph_this($_->[0], $_->[1]);
			if ($result) {
				# return exit code of 1 means success
				$pm->finish(1, \$result); 
			}
			else {
				# exit code of 0 means failure, no correlation to report
				print " Failed to generate graph for ", $dataset_by_id{$_->[0]}, 
					" and ", $dataset_by_id{$_->[1]}, "\n";
				$pm->finish(0);
			}
		}
		$pm->wait_all_children;
	}
	
	# Serial Execution
	else {
		foreach (@to_do) {
			my $result = graph_this($_->[0], $_->[1]); 
			if ($result) {
				$stat_Data->add_row($result);
			}
			else {
				print " Failed to generate graph for ", $dataset_by_id{$_->[0]}, 
					" and ", $dataset_by_id{$_->[1]}, "\n";
			}
		}
	}
}


## Subroutine to interact with user and ask for data set pairs to graph sequentially
sub graph_datasets_interactively {
	
	# print messages 
	
	# present the dataset name list
	print " These are the data sets in the file $infile\n";
	foreach (sort {$a <=> $b} keys %dataset_by_id) {
		print "   $_\t$dataset_by_id{$_}\n";
	}
	
	# collect user input
	print " Enter the numbers for the X and Y datasets separated by a comma  ";
	my $answer = <STDIN>; # get user answer
	chomp $answer;
	
	# this loop will keep going until no dataset (undefined) is returned
	while ($answer) {
		$answer =~ s/\s*//g;
		my ($x, $y) = split /,|&/, $answer;
		if ( 
			(exists $dataset_by_id{$x}) and 
			(exists $dataset_by_id{$y}) 
		) {
			# generate the graph, and return the correlation stats as a string
			my $result = graph_this($x, $y);
			if ($result) {
				$stat_Data->add_row($result);
			}
			else {
				print " Failed to generate graph for ", $dataset_by_id{$x}, 
					" and ", $dataset_by_id{$y}, "\n";
			}
		} 
		else {
			print " One of these numbers is not valid\n";
		}
		
		# get ready for next
		print "\n Enter the next set of datasets or push return to exit  ";
		$answer = <STDIN>;
		chomp $answer;
	}
}


## Subroutine to graph two given data sets
sub graph_this {
	
	# Retrieve the x and y values for plotting
	my ($xid, $yid) = @_;
	
	# get the name of the datasets
	my $xname = $dataset_by_id{$xid}; 
	my $yname = $dataset_by_id{$yid};
	
	# collect the values into separate x and y values
	my (@xvalues, @yvalues);
	$Data->iterate( sub { 
		my $row = shift;
		my $x = $row->value($xid);
		my $y = $row->value($yid);
		# only take numerical data
		# must have a numeric value from both datasets, otherwise skip
		if ( ($x ne '.') and ($y ne '.') ) {
			push @xvalues, $x; 
			push @yvalues, $y;
		}
	});
	
	# Sort the data points
	sort_data(\@xvalues, \@yvalues);
	
	# Smooth by moving average
	if ($moving_average) {
		smooth_data(\@xvalues, \@yvalues);
	}
	
	# Normalize values to a percent ranke (0..1) if requested
	if ($norm) {
		normalize_this(\@xvalues);
		normalize_this(\@yvalues);
	}
	
	# Determine graph type and plot accordingly
	# the correlation statistics will be returned as a string
	if ($type eq 'scatter') {
		return graph_scatterplot($xname, $yname, \@xvalues, \@yvalues);
	} 
	elsif ($type eq 'line') {
		return graph_line_plot($xname, $yname, \@xvalues, \@yvalues);
	}
	elsif ($type eq 'smooth') {
		return graph_smoothed_line_plot($xname, $yname, \@xvalues, \@yvalues);
	}
}



## Subroutine to sort data from lowest to highest
sub sort_data {
	
	# references to the x and y datasets
	my ($xref, $yref) = @_;
	
	# Store the X and Y values temporarily into a sorting hash
	my %sorthash;
	for (my $i = 0; $i < scalar @{$xref}; $i++) {
		# put both the X and Y values into a hash that is keyed by a unique 
		# number - the value's original position in the array
		# the values are in a simple anonymous array
		$sorthash{$i} = [ $xref->[$i], $yref->[$i] ];
	}
	
	# replace the original values with the sorted values
	my $i = 0;
	foreach (sort { 
		$sorthash{$a}->[0] <=> $sorthash{$b}->[0] or
		$a <=> $b
	} keys %sorthash) {
		# sort first by increasing X values, second by the original array 
		# postion 
		# put the sorted X values and corresponding Y value back into the 
		# original array
		# since the number of values aren't changing, this shouldn't be 
		# a problem
		$xref->[$i] = $sorthash{$_}->[0];
		$yref->[$i] = $sorthash{$_}->[1];
		$i++;
	}
}



## Subroutine to smooth data by moving average
sub smooth_data {
	
	# references to the x and y datasets
	my ($xref, $yref) = @_;
	
	# copy the data to new temporary arrays
	my @xdata = @{ $xref };
	my @ydata = @{ $yref };
	
	# empty the original x and y dataset arrays
	# they will be refilled with the smoothened data
	@{ $xref } = ();
	@{ $yref } = ();
	
	# the data should already be sorted
	
	# Determine the moving average values from command line argument
	my ($window, $step) = split /,/, $moving_average;
	my $halfstep = sprintf "%.0f", ($window / 2);
	
	# Smooth the Y data by taking a moving average
	my $stat = Statistics::Descriptive::Sparse->new();
	for (
		my $position = $halfstep; 
		$position < ( scalar(@xdata) - $halfstep) ; 
		$position += $step
	) {
		# The x value will be the midpoint of the window, and the y value will
		# be the average of all values within the window.
		# The x value and corresponding window will be incremented by the step
		# value until the end of the array is reached.
		
		# collect the window values from the sorted y dataset
		my @window; 
		for (my $n = $position - $halfstep; $n <= $position + $halfstep; $n++) {
			push @window, $ydata[$n];
		}
		
		# calculate mean
		$stat->add_data(@window);
		my $window_mean = $stat->mean();
		$stat->clear();
		
		# record the new x and y values
		push @{ $yref }, $window_mean;
		push @{ $xref }, $xdata[$position];
	}
}



## Subroutine to normalize a simple array of data and scale it to 0..1
sub normalize_this {
	my $aref = shift; # the reference to the array to be normalized
	my %qvalues; # a hash for the quantile lookup
	
	# put values into quantile lookup hash
	my $i = 1;
	foreach (sort {$a <=> $b} @{ $aref } ) {
		# for the quantile hash, the key will be the original data value
		# and the hash value will be the rank in the ordered data array
		# the rank is stored in anonymous array since there may be more 
		# than one identical data value and we want the ranks for all of 
		# them
		push @{ $qvalues{$_} }, $i;
		$i++;
	}
	
	# scale to 0..1
	my $total = scalar @{ $aref };
	my $scale_factor = 1 / $total;
	
	# replace the current value with the scaled quantile value
	my $stat = Statistics::Descriptive::Sparse->new();
	for (my $i = 0; $i < $total; $i++) {
		# usually there is only one rank value for each data value
		# sometimes, there are two or more identical data values within the 
		# the dataset, in which case we will simply take the mean rank
		
		my $final_rank;
		if (scalar @{ $qvalues{$aref->[$i]} } == 1) {
			# there's only value
			$final_rank = @{ $qvalues{$aref->[$i]} }[0];
		}
		else {
			# there's more than one value, must take mean
			$stat->add_data( @{ $qvalues{$aref->[$i]} } );
			$final_rank = $stat->mean();
			$stat->clear();
		}
		
		# the input value will be replaced with the adjusted rank value
		$aref->[$i] = $final_rank * $scale_factor; # the new rank on a scale of 0..1
	}
}



## Subroutine to graph a scatterplot and regression line of the data
sub graph_scatterplot {
	# the passed values
	my ($xname, $yname, $xref, $yref) = @_;
	
	# calculate statistics
	my ($q, $m, $r, $rsquared) = get_stats($xref, $yref);
	my $r_formatted = sprintf "%.2f", $r;
	my $rsquared_formatted = sprintf "%.2f", $rsquared;

	
	# Fill the data array for the graphing engine
	# data array will be comprised of three arrays, first is x, second is y (for
	# the scatter plot of points), third is calculated y (for the linear 
	# regression line)
	my @data = ( [ ], [ ], [ ] ); # initialize
	for (my $i = 0; $i < scalar @{ $xref }; $i++) {
		my $x = $xref->[$i]; # current x value
		push @{ $data[0] }, $x; # x values
		push @{ $data[1] }, $yref->[$i]; # real y values
		push @{ $data[2] }, ($x * $m) + $q; # calculated y value
	}
	
	# Now prepare the graph
	my $title = "$xname vs. $yname (r $r_formatted)";
	my $graph = GD::Graph::mixed->new(600,600);
	$graph->set(
		types			=> [qw(points lines)],
		x_label			=> $xname,
		y_label			=> $yname,
		title			=> $title,
		transparent		=> 0,
		markers			=> [1],
		marker_size		=> 1,
		line_types		=> [1],
		line_width		=> 1,
		dclrs			=> [qw(lblue red)],
	) or warn $graph->error;
	
	# set axes
	set_graph_axes($graph);
	
	# Generate graph file name
	my $filename = $xname . '_and_' . $yname;
	if ($out) {
		# add output prefix if requested
		$filename = $out . '_' . $filename;
	}
	$filename = File::Spec->catfile($directory, $filename);
	$filename = check_file_uniqueness($filename, 'png');
	
	# Write the graph file
	my $gd = $graph->plot(\@data) or warn $graph->error;
	open IMAGE, ">$filename" or do {
		warn " Can't open output file '$filename'!\n";
		return;
	};
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	
	# write the graph numbers if requested
	if ($numbers) {
		$filename =~ s/png$/txt/; # change the extension
		my $fh = open_to_write_fh($filename) or 
			warn " Can't write numbers output file '$filename'!\n";
		if ($fh) {
			$fh->print("X_Values\tY_Values\tLinear_Regression_Y_Values\n");
			for (my $i = 0; $i < scalar @{ $data[0] }; $i++) {
				$fh->print("$data[0][$i]\t$data[1][$i]\t$data[2][$i]\n");
			}
			$fh->close;
		}
	}
	
	# return stats
	print " Generated graph for $xname vs $yname\n" . 
		"  Pearson correlation is $r_formatted, R^2 is $rsquared_formatted\n";
	return [ (
		$infile,
		$filename,
		'scatter',
		$norm,
		$xname,
		$yname,
		$q,
		$m,
		$r,
		$rsquared
	) ];
}




## Subroutine to graph a smoothed line of the data
sub graph_line_plot {
	# the passed values
	my ($xname, $yname, $xref, $yref) = @_;
	
	# calculate statistics
	my ($q, $m, $r, $rsquared) = get_stats($xref, $yref);
	my $r_formatted = sprintf "%.2f", $r;
	my $rsquared_formatted = sprintf "%.2f", $rsquared;

	
	# Put the smoothed values into the data array for the graphing engine
	my @data = ($xref, $yref);
	
	# Now (finally) prepare the graph
	my $xtitle = $xname;
	my $ytitle = $yname;
	if ($moving_average) {
		$ytitle .= ' (smoothed)'; # if smoothed by moving average
	}
	my $title = "$xtitle vs. $ytitle (r $r_formatted)";
	my $graph = GD::Graph::lines->new(600,600);
	$graph->set(
		x_label			=> $xtitle,
		y_label			=> $ytitle,
		title			=> $title,
		y_long_ticks	=> 1,
		transparent		=> 0,
	) or warn $graph->error;
	
	# set axes
	set_graph_axes($graph);
	
	# Generate graph file name
	my $filename = $xname . '_and_' . $yname;
	if ($out) {
		# add output prefix if requested
		$filename = $out . '_' . $filename;
	}
	$filename = File::Spec->catfile($directory, $filename);
	$filename = check_file_uniqueness($filename, 'png');
	
	# Write the graph file
	my $gd = $graph->plot(\@data) or warn $graph->error;
	open IMAGE, ">$filename" or do {
		warn " Can't open output file '$filename'!\n";
		return;
	};
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	
	# write the graph numbers if requested
	if ($numbers) {
		$filename =~ s/png$/txt/; # change the extension
		my $fh = open_to_write_fh($filename) or 
			warn " Can't write numbers output file '$filename'!\n";
		if ($fh) {
			$fh->print("X_Values\tY_Values\n");
			for (my $i = 0; $i < scalar @{ $data[0] }; $i++) {
				$fh->print("$data[0][$i]\t$data[1][$i]\n");
			}
			$fh->close;
		}
	}
	
	# return stats
	print " Generated graph for $xname vs $yname\n" . 
		"  Pearson correlation is $r_formatted, R^2 is $rsquared_formatted\n";
	return [ (
		$infile,
		$filename,
		'scatter',
		$norm,
		$xname,
		$yname,
		$q,
		$m,
		$r,
		$rsquared
	) ];
}


## Subroutine to graph a smoothed line of the data
sub graph_smoothed_line_plot {
	# the passed values
	my ($xname, $yname, $xref, $yref) = @_;
	
	
	# calculate statistics
	my ($q, $m, $r, $rsquared) = get_stats($xref, $yref);
	my $r_formatted = sprintf "%.2f", $r;
	my $rsquared_formatted = sprintf "%.2f", $rsquared;

	
	# Put the smoothed values into the data array for the graphing engine
	my @data = ($xref, $yref);
	
	# Now (finally) prepare the graph
	my $xtitle = $xname;
	my $ytitle = $yname;
	if ($moving_average) {
		$ytitle .= ' (smoothed)'; # if smoothed by moving average
	}
	my $title = "$xtitle vs. $ytitle (r $r_formatted)";
	my $graph = GD::Graph::smoothlines->new(600,600);
	$graph->set(
		x_label			=> $xtitle,
		y_label			=> $ytitle,
		title			=> $title,
		y_long_ticks	=> 1,
		transparent		=> 0,
	) or warn $graph->error;
	
	# set axes
	set_graph_axes($graph);
	
	# Generate graph file name
	my $filename = $xname . '_and_' . $yname;
	if ($out) {
		# add output prefix if requested
		$filename = $out . '_' . $filename;
	}
	$filename = File::Spec->catfile($directory, $filename);
	$filename = check_file_uniqueness($filename, 'png');
	
	# Write the graph file
	my $gd = $graph->plot(\@data) or warn $graph->error;
	open IMAGE, ">$filename" or do {
		warn " Can't open output file '$filename'!\n";
		return;
	};
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	
	# write the graph numbers if requested
	if ($numbers) {
		$filename =~ s/png$/txt/; # change the extension
		my $fh = open_to_write_fh($filename) or 
			warn " Can't write numbers output file '$filename'!\n";
		if ($fh) {
			$fh->print("X_Values\tY_Values\n");
			for (my $i = 0; $i < scalar @{ $data[0] }; $i++) {
				$fh->print("$data[0][$i]\t$data[1][$i]\n");
			}
			$fh->close;
		}
	}
	
	# return stats
	print " Generated graph for $xname vs $yname\n" . 
		"  Pearson correlation is $r_formatted, R^2 is $rsquared_formatted\n";
	return [ (
		$infile,
		$filename,
		'scatter',
		$norm,
		$xname,
		$yname,
		$q,
		$m,
		$r,
		$rsquared
	) ];
}


## Set some generic axis options regardless of graph type
sub set_graph_axes {
	my $graph = shift;
	
	# set ticks and label number format
	$graph->set(
		x_tick_number	=> $x_ticks,
		y_tick_number	=> $y_ticks,
		x_number_format	=> '%.' . $x_format . 'f', # "%.2f",
		y_number_format => '%.' . $y_format . 'f',
		x_label_position => 0.5,
	) or warn $graph->error;
	
	
	# set min max values on the graph
	if ($norm) { 
		# explicitly set the axis values if the data was normalized
		$graph->set(
			y_min_value		=> 0,
			y_max_value		=> 1,
			x_min_value		=> 0,
			x_max_value		=> 1
		) or warn $graph->error;
	} 
	else {
		# set each x,y minimum and maximum value if user defined
		if (defined $x_min) {
			$graph->set(x_min_value => $x_min) or warn $graph->error;
		}
		if (defined $y_min) {
			$graph->set(y_min_value => $y_min) or warn $graph->error;
		}
		if (defined $x_max) {
			$graph->set(x_max_value => $x_max) or warn $graph->error;
		}
		if (defined $y_max) {
			$graph->set(y_max_value => $y_max) or warn $graph->error;
		}
	}
	# otherwise we let it calculate automatic values
	
	# the default tiny font is too small for 800x600 graphic
	$graph->set_x_axis_font(GD::Font->Small) or warn $graph->error;
	$graph->set_y_axis_font(GD::Font->Small) or warn $graph->error;
}



## Determine statistics and perform linear regression on the data
sub get_stats {
	
	my ($xref, $yref) = @_;
	
	# initialize variables
	my $q; # intercept
	my $m; # slope, satisfies equation $y = $m*$x + $q
	my $r; # the pearson correlation coefficient between x and y
	my $rms; # root mean square, I won't be using this
	my $rsquared; # the coefficient of determination
	
	# calculate
	my $stat = Statistics::Descriptive::Full->new(); # initialize
	$stat->add_data($yref);
	($q, $m, $r, $rms) = $stat->least_squares_fit( @{ $xref } );
	$rsquared = $r ** 2;
	
	# return
	return ($q, $m, $r, $rsquared);
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

graph_data.pl

A script to graph XY line or dot plots between data sets.

=head1 SYNOPSIS

graph_data.pl [--options] <filename>
  
  Options:
  --in <filename>
  --type [scatter | line | smooth]
  --pair <X_index>,<Y_index>
  --index <X_index&Y_index,...>
  --all
  --ma <window>,<step>
  --norm
  --min=<value>
  --xmin=<value>
  --ymin=<value>
  --max=<value>
  --xmax=<value>
  --ymax=<value>
  --ticks <integer>
  --xticks <integer>
  --yticks <integer>
  --format <integer>
  --xformat <integer>
  --yformat <integer>
  --out <base_filename>
  --dir <foldername>
  --cpu <integer>
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

=item --type [scatter | line | smooth]

Indicate the type of graph to plot. Three different values are 
accepted. B<Scatter> graphs will plot all pairwise X,Y values and 
draw a linear regression line through them. B<Line> graphs will 
plot a continuous line connecting all pairwise X,Y data values. 
With noisy data, the plot will look the best if the values are 
first smoothed using a moving average (--ma option). Finally, a 
B<smooth> graph is the same as a line graph, except the data 
values are smoothed using a bezier curve function. Note that the 
bezier smoothing function is not equivalent or as effective 
as a moving average. The default value is a scatter plot.

=item --pair <X_index>,<Y_index>

Specify the two datasets to plot together. Use the datasets'
index (0-based) expressed as 'X,Y' with no spaces. Use the 
option repeatedly to plot multiple graphs. If no datasets are 
set, then the lists may be selected interactively from a list.

=item --index <X_index&Y_index,...>

An alternative method of specifying the datasets as a 
comma-delimited list of X&Y datasets, where the X and Y 
indices are demarcated by an ampersand. If no datasets are 
set, then the lists may be selected interactively from a list.

=item --all

Indicate that all available datasets in the input file should be
plotted together. Redundant graphs are skipped, e.g. Y,X versus X,Y.
If you wish to graph only a subset of datasets, provide a list 
and/or range using the --index option.

=item --norm

Datasets should (not) be normalized by converting to percentile 
rank values (0..1). This is helpful when the two datasets are 
not in similar scales. Default is false.

=item --ma <window>,<step>

Specify the values to smooth the data by moving average. 
Express the values as 'window,step', no spaces. All Y values 
within the window are averaged together and plotted 
against the window midpoint X value. The window position is 
then incremented by the step size. The step size should be 
equal or less than the window size. Both values must be 
real integers. This data manipulation is extremely useful and 
recommended for noisy datasets.

=item --min=<value>

=item --xmin=<value>

=item --ymin=<value>

Specify explicitly the minimum values for either the X or Y axes. 
Both may be set independently or to the same value with the --min 
option. The default is automatically calculated.

=item --max=<value>

=item --xmax=<value>

=item --ymax=<value>

Specify explicitly the maximum values for either the X or Y axes. 
Both may be set independently or to the same value with the --max 
option. The default is automatically calculated.

=item --ticks <integer>

=item --xticks <integer>

=item --yticks <integer>

Specify explicitly the number of major ticks for either the X or Y axes. 
Both may be set independently or to the same value with the --ticks 
option. The default is 4.

=item --format <integer>

=item --xformat <integer>

=item --yformat <integer>

Specify explicitly the number of decimal places to format the labels 
for the major axes' ticks. Both may be set independently or to the same 
value with the --format option. The default is 0.

=item --out <base_filename>

Optionally specify the output filename prefix.

=item --dir <foldername>

Specify an optional name for the output subdirectory name. Default
is the input filename base with '_graphs' appended.

=item --cpu <integer>

Specify the number of CPU cores to execute in parallel. This requires 
the installation of Parallel::ForkManager. With support enabled, the 
default is 2. Disable multi-threaded execution by setting to 1. 
Parallel execution is only applicable when a list of datasets are 
provided or the --all option is enabled; interactive execution is 
performed serially.

=item --version

Print the version number.

=item --help

Display this help as a POD.

=back

=head1 DESCRIPTION

This program will graph pairwise data sets against each other. This is useful 
for determining correlations between data sets. The graphs are generated as 
PNG images and written to a subdirectory. 

Three types of graphs may be generated, specified by the --type argument.

A scatter plot will plot all pairwise values between two datasets as a 
point on an X,Y graph. A line representing the linear regression of the 
two datasets is plotted above the points. 

A line plot will plot the pairwise values between two datasets as a 
continuous line. 

A smooth plot is a line plot of pairwise values between two datasets and 
smoothed by a bezier curve.

The datasets to be plotted may be specified either on the command line 
using the --pair or --index arguments. If not specified, the program 
defaults to an interactive mode and the user may repeatedly select the 
datasets from a list of available datasets in the data file. 

The data in the datasets may be manipulated in several ways prior to plotting.
The data may be converted to a percentile rank, smoothed by a moving average, 
constrained to minimum and maximum values, etc. 

If the graph doesn't look like you expect, and you are not normalizing by 
converting to a percent rank (--norm), try explicitly setting the --min and 
--max values. The GD::Graph module tries its best at setting these 
automatically, but sometimes does funny things.

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
