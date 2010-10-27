#!/usr/bin/perl

# A script to graph correlation plots for two microarray data sets

use strict;
use Pod::Usage;
use Getopt::Long;
eval {
	use GD::Graph::lines; # for the line type graph
	use GD::Graph::mixed; # for the scatter plot
};
eval {
	use GD::Graph::smoothlines; # for bezier smoothed line graph
};
use Statistics::Descriptive;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper qw(
	load_tim_data_file
);

print "\n This script will graph correlation plots for two microarry data sets\n\n";

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
	$log,
	$directory,
	$numbers,
	$help, 
);
my @pairs; # an array of pairs
GetOptions( 
	'in=s'    => \$infile, # the input file
	'type=s'  => \$type, # the type of graph to generate
	'index=s' => \$index, # a list of x,y pairs to plot
	'pair=s'  => \@pairs, # the x,y pairs to plot
	'all'     => \$all, # flag to plot all data sets pairwise
	'norm!'   => \$norm, # flag to skip percentile normalizing
	'ma=s'    => \$moving_average, # window, setp values for moving average
	'min=i'   => \$min, # mininum axis coordinate
	'max=i'   => \$max, # maximum axis coordinate
	'log!'    => \$log, # values are in log, respect log status
	'dir=s'   => \$directory, # optional name of the graph directory
	'numbers' => \$numbers, # print the graph numbers in addition to the graph
	'help'    => \$help, # flag to print help
);

if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		die " Missing input file!\n";
	}
}
if ($type) {
	# check the requested type
	my %accepted_types = (
		# a quickie hash to look up valid graph types, add more here
		'scatter'  => 1,
		'line'     => 1,
		'smooth'   => 1,
	);
	unless (exists $accepted_types{$type} ) {
		die " Unknown graph type '$type'!\n Please use --help for more information\n";
	}
}
else {
	$type = 'scatter';
	print " Using default graph type of 'scatter'\n";
}
unless (defined $norm) {
	# default is no normalization (percent rank)
	$norm = 0;
}
unless (defined $log) {
	# default is no log
	$log = 0;
}





####### Main ###########

### Prepare global variables and set up for execution

# load the file using subroutine from tim_file_helper.pm
print " Loading data from file $infile....\n";
my $main_data_ref = load_tim_data_file($infile);
unless ($main_data_ref) {
	die " No data loaded!\n";
}
my $data_table_ref = $main_data_ref->{'data_table'};

# load the dataset names into hashes
my %dataset_by_id; # hashes for name and id
for (my $i = 0; $i < $main_data_ref->{'number_columns'}; $i++) {
	my $name = $main_data_ref->{$i}{'name'};
	if (
		# check column header names for gene or window attribute information
		# these won't be used for graph generation, so we'll skip them
		$name =~ /name/i or 
		$name =~ /class/i or
		$name =~ /alias/i or
		$name =~ /^chr/i or
		$name =~ /start/i or
		$name =~ /stop/i or
		$name =~ /end/i or
		$name =~ /^probe/i
	) { 
		# skip on to the next header
		next; 
	} 
	else { 
		# record the data set name
		$dataset_by_id{$i} = $name;
	}
}	


# Prepare output directory
unless ($directory) {
	# directory may be specified from a command line argument, otherwise
	# generate default directory from input file name
	$directory = $main_data_ref->{'basename'} . '_graphs';
}
unless (-e "$directory") {
	mkdir $directory or die "Can't create directory '$directory'\n";
}


# Set up an output array of the correlation values
my @correlation_output; 






### Get the list of datasets to pairwise compare

# A list of dataset pairs was provided upon execution
if ($index) {
	my @list = split /,/, $index;
	push @pairs, @list;
}
if (@pairs) {
	foreach my $pair (@pairs) {
		my ($x, $y) = split /,|&/, $pair;
		if ( 
			(exists $dataset_by_id{$x}) and 
			(exists $dataset_by_id{$y}) 
		) {
			graph_this($x, $y);
		} 
		else {
			print " One of these numbers ($x, $y) is not valid\n";
		}
	}
}

# All datasets in the input file will be compared to each other
elsif ($all) {
	print " Comparing all data set pairwise combinations....\n";
	graph_all_datasets();
}

# Interactive session
else {
	graph_datasets_interactively();
}



### Write out the linear regression statistics

# open file
my $statfile = $main_data_ref->{'basename'} . "_stats.txt";
if (-e "$directory/$statfile") {
	# if the file exists, append to it
	open STATFILE, ">>$directory/$statfile" or warn " unable to write $directory/$statfile\n";
} 
else {
	# open a new file
	open STATFILE, ">$directory/$statfile" or warn " unable to write $directory/$statfile\n";
	print STATFILE join("\t", qw(
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
	) ) . "\n"; 
}

# print to file
foreach (@correlation_output) {
	print STATFILE "$_\n";
}
close STATFILE;


### The end of the main program	




##### Subroutines ##########



## subroutine to graph all the data sets in the loaded data file
sub graph_all_datasets {
	
	# put the dataset IDs into sorted array
	my @list = sort {$a <=> $b} keys %dataset_by_id; 
	
	# Now process through the list 
	# the first dataset will be x, and the subsequent will by y
	# all of the ys will be done for a given x, then x is incremented and the ys repeated
	# any y less than the current x will be skipped to avoid repeats
	for (my $x = 0; $x < (scalar @list - 1); $x++) {
		for (my $y = 1; $y < (scalar @list); $y++) {
			if ($y <= $x) {next}; # skip pairs we've already done
			graph_this($list[$x], $list[$y]);
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
			graph_this($x, $y);
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
	print " Preparing graph for $xname vs. $yname....\n";
	
	# collect the values
	my (@xvalues, @yvalues);
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) { 
		# walk through the data file
		my $x = $data_table_ref->[$row][$xid];
		my $y = $data_table_ref->[$row][$yid];
		# only take numerical data
		# must have a numeric value from both datasets, otherwise skip
		if ( ($x ne '.') and ($y ne '.') ) {
			push @xvalues, $x; # put into the values array
			push @yvalues, $y;
		}
	}
	
	# check log status
	my $changed_x_log_status = 0;
	my $changed_y_log_status = 0;
	if ($log) {
		# we need to work with log values
		# we will obey any log2 flags present in the metadata first
		# failing that, we will blindly assume the data needs to be transformed
		
		# first check the metadata for the x dataset
		if (exists $main_data_ref->{$xid}{'log2'}) {
			# x dataset has log2 flag
			if ($main_data_ref->{$xid}{'log2'} == 1) {
				# it is in log2 space, so de-log
				@xvalues = map { 2 ** $_ } @xvalues;
				$changed_x_log_status = 1;
			}
		} else {
			# x dataset has no log2 flag
			# then we will blindly assume it needs to be transformed
			@xvalues = map { 2 ** $_ } @xvalues;
			$changed_x_log_status = 1;
		}
		# then check the metadata for the y dataset
		if (exists $main_data_ref->{$yid}{'log2'}) {
			# y dataset has log2 flag
			if ($main_data_ref->{$yid}{'log2'} == 1) {
				# it is in log2 space, so de-log
				@yvalues = map { 2 ** $_ } @yvalues;
				$changed_y_log_status = 1;
			}
		} else {
			# y dataset has no log2 flag
			# then we will blindly assume it needs to be transformed
			@yvalues = map { 2 ** $_ } @yvalues;
			$changed_y_log_status = 1;
		}
	}
	
	
	# Sort the data points
	sort_data(\@xvalues, \@yvalues);
	
	
	# Smooth by moving average
	if ($moving_average) {
		# use subroutine
		smooth_data(\@xvalues, \@yvalues);
	}
	
	
	# Either normalize or convert back to log
	if ($norm) {
		# Normalize values to a percent rank (0..1)
		# they do not need to be re-logged
		normalize_this(\@xvalues);
		normalize_this(\@yvalues);
	}
	elsif ($changed_x_log_status or $changed_y_log_status) {
		# Values were originally log2, change them back before plotting
		if ($changed_x_log_status) {
			@xvalues = map { log($_) / log(2) } @xvalues;
		}
		if ($changed_y_log_status) {
			@yvalues = map { log($_) / log(2) } @yvalues;
		}
	}
	
	
	# Determine graph type and plot accordingly
	if ($type eq 'scatter') {
		graph_scatterplot($xname, $yname, \@xvalues, \@yvalues);
	} 
	elsif ($type eq 'line') {
		graph_line_plot($xname, $yname, \@xvalues, \@yvalues);
	}
	elsif ($type eq 'smooth') {
		graph_smoothed_line_plot($xname, $yname, \@xvalues, \@yvalues);
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
	unless (
		defined $window and      # both values are defined
		defined $step and 
		$window =~ /^\d+$/ and   # both values are numbers
		$step =~ /^\d+$/ and 
		$step <= $window         # step is less or equal to window
	){
		die " moving average values is not (well) defined!\n Use --help for more information\n";
	}
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
		
		markers			=> [1],
		marker_size		=> 1,
		line_types		=> [1],
		line_width		=> 1,
		
		x_tick_number	=> 4,
		y_tick_number	=> 4,
		x_number_format	=> "%.2f",
		y_number_format => "%.2f",
		x_label_position => 0.5,
		
		transparent		=> 0,
		dclrs			=> [qw(lblue red)],
	
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
	elsif (defined $min or defined $max) {
		if (defined $min) {
			$graph->set(
				y_min_value => $min,
				x_min_value => $min,
			) or warn $graph->error;
		}
		if (defined $max) {
			$graph->set(
				y_max_value => $max,
				x_max_value => $max,
			) or warn $graph->error;
		}
	}
	# otherwise we let it calculate automatic values
		
	# Write the graph file
	my $gd = $graph->plot(\@data) or warn $graph->error;
	my $filename = $xname . '_and_' . $yname . '.png';
	my $filenumber = 1; # an incremental number to make a unique file name
	while (-e "$directory/$filename") {
		$filename = $xname . '_and_' . $yname . '_' . $filenumber . '.png';
		$filenumber++;
	}
	open IMAGE, ">$directory/$filename" or die " Can't open output image file!\n";
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	
	# record stats
	print "  Pearson correlation is $r_formatted, R^2 is $rsquared_formatted\n";
	push @correlation_output, join("\t", (
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
	) );
	
	# write the graph numbers if requested
	if ($numbers) {
		$filename =~ s/png$/txt/; # change the extension
		open NUMBERS, ">$directory/$filename" or die " Can't open output text file!\n";
		print NUMBERS "X_Values\tY_Values\tLinear_Regression_Y_Values\n";
		for (my $i = 0; $i < scalar @{ $data[0] }; $i++) {
			print NUMBERS "$data[0][$i]\t$data[1][$i]\t$data[2][$i]\n";
		}
		close NUMBERS;
	}
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
		
		x_tick_number	=> 4,
		y_tick_number	=> 4,
		y_long_ticks	=> 1,
		x_number_format	=> "%.2f",
		y_number_format => "%.2f",
		x_label_position => 0.5,
		
		transparent		=> 0,
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
	elsif (defined $min or defined $max) {
		if (defined $min) {
			$graph->set(
				y_min_value => $min,
				x_min_value => $min,
			) or warn $graph->error;
		}
		if (defined $max) {
			$graph->set(
				y_max_value => $max,
				x_max_value => $max,
			) or warn $graph->error;
		}
	}
	# otherwise we let it calculate automatic values
	
	# Write the graph file
	my $gd = $graph->plot(\@data) or warn $graph->error;
	my $filename = "$xname" . '_and_' . "$yname" . '.png';
	my $filenumber = 1; # an incremental number to make a unique file name
	while (-e "$directory/$filename") {
		$filename = "$xname" . '_and_' . "$yname" . '_' . $filenumber . '.png';
		$filenumber++;
	}
	open IMAGE, ">$directory/$filename" or die " Can't open output file!\n";
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	
	# record stats
	print "  Pearson correlation is $r_formatted, R^2 is $rsquared_formatted\n";
	push @correlation_output, join("\t", (
		$infile,
		$filename,
		'line',
		$norm,
		$xname,
		$yname . ' (smoothed)',
		$q,
		$m,
		$r,
		$rsquared
	) );
	
	# write the graph numbers if requested
	if ($numbers) {
		$filename =~ s/png$/txt/; # change the extension
		open NUMBERS, ">$directory/$filename" or die " Can't open output text file!\n";
		print NUMBERS "X_Values\tY_Values\n";
		for (my $i = 0; $i < scalar @{ $data[0] }; $i++) {
			print NUMBERS "$data[0][$i]\t$data[1][$i]\n";
		}
		close NUMBERS;
	}
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
		
		x_tick_number	=> 4,
		y_tick_number	=> 4,
		y_long_ticks	=> 1,
		x_number_format	=> "%.2f",
		y_number_format => "%.2f",
		x_label_position => 0.5,
		
		transparent		=> 0,
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
	elsif (defined $min or defined $max) {
		if (defined $min) {
			$graph->set(
				y_min_value => $min,
				x_min_value => $min,
			) or warn $graph->error;
		}
		if (defined $max) {
			$graph->set(
				y_max_value => $max,
				x_max_value => $max,
			) or warn $graph->error;
		}
	}
	# otherwise we let it calculate automatic values
	
	# Write the graph file
	my $gd = $graph->plot(\@data) or warn $graph->error;
	my $filename = "$xname" . '_and_' . "$yname" . '.png';
	my $filenumber = 1; # an incremental number to make a unique file name
	while (-e "$directory/$filename") {
		$filename = "$xname" . '_and_' . "$yname" . '_' . $filenumber . '.png';
		$filenumber++;
	}
	open IMAGE, ">$directory/$filename" or die " Can't open output file!\n";
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	
	# record stats
	print "  Pearson correlation is $r_formatted, R^2 is $rsquared_formatted\n";
	push @correlation_output, join("\t", (
		$infile,
		$filename,
		'line',
		$norm,
		$xname,
		$yname . ' (smoothed)',
		$q,
		$m,
		$r,
		$rsquared
	) );
	
	# write the graph numbers if requested
	if ($numbers) {
		$filename =~ s/png$/txt/; # change the extension
		open NUMBERS, ">$directory/$filename" or die " Can't open output text file!\n";
		print NUMBERS "X_Values\tY_Values\n";
		for (my $i = 0; $i < scalar @{ $data[0] }; $i++) {
			print NUMBERS "$data[0][$i]\t$data[1][$i]\n";
		}
		close NUMBERS;
	}
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







__END__

=head1 NAME

graph_data.pl

=head1 SYNOPSIS
 
  graph_data.pl [--options] <filename>
  
  Options:
  --in <filename>
  --type [scatter | line | smooth]
  --pair <X_index>,<Y_index>
  --index <X_index&Y_index,...>
  --all
  --ma <window>,<step>
  --(no)norm
  --(no)log
  --min=<value>
  --max=<value>
  --dir <foldername>
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4


=item --in <filename>

Specify the input data file. A tim data formatted file is 
preferred, although any tabbed delimited text file where 
datasets are represented as columns will suffice. The file 
may be gzip'ed.

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

=item --(no)log

Indicate whether dataset values are in log2 space or not. If set 
to true and the log2 status is indicated in the metadata, then 
the metadata status is preserved. Default is false (and metadata 
ignored).

=item --(no)norm

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

=item --min=<value>, --max=<value>

Specify explicitly the minimum and/or maximum axis value(s). 
Default is calculated. Sometimes setting this improves 
the appearance of the graph. They may be set independently.

=item --dir <foldername>

Specify an optional name for the output subdirectory name. Default
is the input filename base with '_graphs' appended.

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
converted from log2 values, constrained to minimum and maximum values, etc. 
See options for details.

If the graph doesn't look like you expect, and you are not normalizing by 
converting to a percent rank (--norm), try explicitly setting the --min and 
--max values. The GD::Graph module tries its best at setting these 
automatically, but sometimes does funny things, particularly with log2 data 
that spans 0.


=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112






