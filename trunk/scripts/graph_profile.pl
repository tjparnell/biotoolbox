#!/usr/bin/perl

# A script to graph bezier-curve smoothed profile plots for summed genomic 
# data, e.g. average gene, TSS data, etc.

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(mean median);
use GD::Graph::smoothlines; # for bezier smoothed line graph
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper;

print "\n This script will graph profile plots of genomic data\n\n";

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
	$all, 
	$data, 
	$center, 
	$min, 
	$max, 
	$x_index,
	$log,
	$directory,
	$help, 
);
GetOptions( 
	'in=s'    => \$infile, # the input file
	'index=s' => \$data, # a list of datasets to graph
	'all'     => \$all, # flag to plot all data sets individually
	'cen!'    => \$center, # flag to center normalize the datasets
	'min=f'   => \$min, # mininum y axis coordinate
	'max=f'   => \$max, # maximum y axis coordinate
	'x=i'     => \$x_index, # index of the X-axis values
	'log!'    => \$log, # values are in log, respect log status
	'dir=s'   => \$directory, # optional name of the graph directory
	'help'    => \$help, # flag to print help
);

if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}


##### Check required and default variables

unless ($infile) {
	if (@ARGV) {
		# something left over from the commandline, must be infile
		$infile = shift @ARGV;
	}
	else {
		die " Missing input file!\n Please use --help for more information\n";
	}
}
unless (defined $center) {
	# default is no median center normalization
	$center = 0;
}
unless (defined $log) {
	# default is no log
	$log = 0;
}





####### Main ###########

### Load the file
# load the file using subroutine from tim_file_helper.pm
print " Loading data from file $infile....\n";
my $main_data_ref = load_tim_data_file($infile);
unless ($main_data_ref) {
	die " No data loaded!\n";
}
my $data_table_ref = $main_data_ref->{'data_table'};

# load the dataset names into hashes
my %dataset_by_id; # hashes for name and id
my ($name_i, $midpoint_i); # index numbers for the window name and midpoint
for (my $i = 0; $i < $main_data_ref->{'number_columns'}; $i++) {
	# check column header names for gene or window attribute information
	my $name = $main_data_ref->{$i}{'name'}; # name of the dataset
	# identify the indices for the window or midpoint
	if ($name =~ /window/i) {
		$name_i = $i;
	}
	elsif ($name =~ /midpoint/i) {
		$midpoint_i = $i;
	}
	else { 
		# record the data set name
		$dataset_by_id{$i} = $name;
	}
}	


# Prepare output directory
unless ($directory) {
	$directory = $main_data_ref->{'basename'} . '_graphs';
}
unless (-e "$directory") {
	mkdir $directory or die "Can't create directory $directory\n";
}


### Get the list of datasets to pairwise compare

# A list of dataset pairs was provided upon execution
if ($data) {
	my @datasets = _parse_list($data);
	foreach my $dataset (@datasets) {
		
		# check for multiple datasets to be graphed together
		if ($dataset =~ /&/) {
			my @list = split /&/, $dataset;
			my $check = -1; # assume all are correct initially
			foreach (@list) {
				unless (exists $dataset_by_id{$_} ) {
					$check = $_; # remember the bad one
					last; # one at a time, please
				}
			}
			
			# graph the datasets
			if ($check == -1) {
				# all dataset indexes are valid
				graph_this(@list);
			}
			else {
				print " The index '$check' is not valid! Nothing graphed\n";
			}
		}
		
		# single dataset
		else {
			if (exists $dataset_by_id{$dataset}) {
				graph_this($dataset);
			}
			else {
				print " The index '$dataset' is not valid! Nothing graphed\n";
			}
		}
		
	}
}

# All datasets in the input file will be compared to each other
elsif ($all) {
	print " Graphing all datasets individually....\n";
	foreach (sort { $a <=> $b } keys %dataset_by_id) {
		# graph each dataset index alone in incremental order
		graph_this($_);
	}
}

# Interactive session
else {
	graph_datasets_interactively();
}


### The end



##### Subroutines ##########


## Subroutine to interact with user and ask for data set pairs to graph sequentially
sub graph_datasets_interactively {
	
	# print messages 
	
	# present the dataset name list
	print " These are the data sets in the file $infile\n";
	foreach (sort {$a <=> $b} keys %dataset_by_id) {
		print "   $_\t$dataset_by_id{$_}\n";
	}
	
	# collect user input
	print " Enter the dataset indices separated by a comma  ";
	my $answer = <STDIN>; # get user answer
	chomp $answer;
	
	# this loop will keep going until no dataset (undefined) is returned
	while ($answer) {
		$answer =~ s/\s*//g;
		my @datasets = _parse_list($answer);
		
		# validate the indices
		my $check = -1; # assume all are correct initially
		foreach (@datasets) {
			unless (exists $dataset_by_id{$_} ) {
				$check = $_; # remember the bad one
				last; # one at a time, please
			}
		}
		
		# graph the datasets
		if ($check == -1) {
			# all dataset indexes are valid
			graph_this(@datasets);
		}
		else {
			print " The index '$check' is not valid! Nothing graphed\n";
		}
		
		# get ready for next
		print "\n Enter the next set of datasets or push return to exit  ";
		$answer = <STDIN>;
		chomp $answer;
	}
}


## Subroutine to graph two given data sets
sub graph_this {
	
	# Retrieve the dataset indices for plotting
	my @datasets = @_;
	
	# Get the name of the datasets
	# join the dataset names together into a string
	my $graph_name = join " & ", map { $dataset_by_id{$_} } @datasets;
	print " Preparing graph for $graph_name....\n";
	
	# Collect the values
	my @graph_data; # a complex array of arrays
	
	# first collect the x values, which may be either user defined, the 
	# window midpoint, or the window name
	unless (defined $x_index) {
		if (defined $midpoint_i) {
			# we'll use the midpoint preferentially, for now
			$x_index = $midpoint_i;
		} 
		elsif (defined $name_i) {
			# if we didn't find that, use the window name
			$x_index = $name_i;
		}
		else {
			die "WTH!? no x index found!\n";
		}
	}
	push @graph_data, []; # first empty sub-array for the x values
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		# we'll walk through the data table
		# and push the appropriate x values into the graph data first sub-array
		push @{ $graph_data[0] }, $data_table_ref->[$row][$x_index];
	}
	
	# next collect the dataset values
	foreach my $dataset (@datasets) {
		push @graph_data, []; # new empty sub-array for each y data
		for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
			# we'll walk through the data table
			# and push the requested y values into the graph data last sub-array
			my $value = $data_table_ref->[$row][$dataset];
			if ($value ne '.') {
				# a real value
				push @{ $graph_data[-1] }, $value;
			}
			else {
				# a null value
				push @{ $graph_data[-1] }, undef;
			}
		}
	}
	
	
	# Check log status
	if ($log) {
		# we need to work with log values
		# we will obey any log2 flags present in the metadata first
		# failing that, we will blindly assume the data needs to be transformed
		
		foreach my $ydata (1..$#datasets) {
			# we'll walk through datasets array of requested indices
			# using the index number in the datasets array rather than the 
			# the index number of the main data table
			# this way we can convert to the index number of the corresponding 
			# dataset in the graph_data array
			# skipping the x dataset, of course
			
			# first check the metadata for the x dataset
			if (exists $main_data_ref->{ $datasets[$ydata] }{'log2'}) {
				# x dataset has log2 flag
				if ($main_data_ref->{ $datasets[$ydata] }{'log2'} == 1) {
					# it is in log2 space, so de-log
					
					@{ $graph_data[ $ydata + 1 ] } = 
						map { 2 ** $_ } @{ $graph_data[ $ydata + 1 ] };
						# the corresponding @graph_data index is $ydata + 1, 
						# since we're skipping the x dataset
						# convert all values in that dataset from log2 space
				}
			} 
			else {
				# dataset has no log2 flag
				# then we will blindly assume it needs to be transformed
				
				@{ $graph_data[ $ydata + 1 ] } = 
					map { 2 ** $_ } @{ $graph_data[ $ydata + 1 ] };
					# the corresponding @graph_data index is $ydata + 1, 
					# since we're skipping the x dataset
					# convert all values in that dataset from log2 space
			}
			
		}
	}	
	
	
	# Median center the y dataset values if requested
	if ($center) {
		# We will center each Y dataset
		
		foreach my $ydata (1..$#graph_data) {
			# using the index of the graph data sub_arrays
			# skipping the first sub-array, the x values
			
			# determine the median value of this dataset
			my $median_value = median( @{ $graph_data[$ydata] } );
			
			# subtract the median value from each value
			@{ $graph_data[$ydata] } = 
				map { $_ - $median_value } @{ $graph_data[$ydata] };
		}
	}
	
	
	
	
	####  Prepare the graph  ####
	
	# Initialize the graph
	my $graph = GD::Graph::smoothlines->new(800,600);
	$graph->set(
		'x_label'           => $main_data_ref->{$x_index}{'name'},
		'title'             => 'Profile' . $main_data_ref->{'feature'},
		'transparent'       => 0, # no transparency
		'y_tick_number'     => 8,
		'y_long_ticks'      => 0,
		'zero_axis'         => 1, 
		'y_number_format'   => "%.2f",
		'x_label_position'  => 0.5,
		'x_label_skip'      => 2,
		'line_width'        => 2,
		
	) or warn $graph->error;
	
	# Set the legend
	$graph->set_legend( map { $dataset_by_id{$_} } @datasets ) 
		or warn $graph->error;
	$graph->set('legend_placement' => 'BR') or warn $graph->error; # BottomRight
	
	# Set fonts
	# the default tiny font is too small for 800x600 graphic
	# possibilities: gdTinyFont gdSmallFont gdMediumBoldFont gdLargeFont gdGiantFont
	$graph->set_legend_font(GD::gdMediumBoldFont) or warn $graph->error;
	$graph->set_x_label_font(GD::gdLargeFont) or warn $graph->error;
	$graph->set_title_font(GD::gdGiantFont) or warn $graph->error;
	$graph->set_x_axis_font(GD::gdSmallFont) or warn $graph->error;
	$graph->set_y_axis_font(GD::gdSmallFont) or warn $graph->error;
	
	
	# Set min max values on the graph if explicitly defined
	if (defined $min) {
		$graph->set( 'y_min_value' => $min ) or warn $graph->error;
	}
	if (defined $max) {
		$graph->set( 'y_max_value' => $max ) or warn $graph->error;
	}
	# otherwise we let it calculate automatic values
	
	
	# Write the graph file
	my $gd = $graph->plot(\@graph_data) or warn $graph->error;
	my $filename =  $graph_name . '_profile.png';
	my $filenumber = 1; # an incremental number to make a unique file name
	while (-e "$directory/$filename") {
		# it's possible that re-running the program on the same datasets
		# will result in the file name and overwriting pre-existing files.
		# to avoid this, we will check for the existence of a file with the 
		# proposed filename and append a unique number if necessary
		$filename = $graph_name . '_profile' . $filenumber . '.png';
		$filenumber++;
	}
	open IMAGE, ">$directory/$filename" or die " Can't open output file!\n";
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
	
	print " wrote file '$filename'\n";
}



sub _parse_list {
	# this subroutine will parse a string into an array
	# it is designed for a string of numbers delimited by commas
	# a range of numbers may be specified using a dash
	# hence 1,2,5-7 would become an array of 1,2,5,6,7
	
	my $string = shift;
	$string =~ s/\s+//g; 
	my @list;
	foreach (split /,/, $string) {
		# check for a range
		if (/\-/) { 
			my ($start, $stop) = split /\-/;
			# add each item in the range to the list
			for (my $i = $start; $i <= $stop; $i++) {
				push @list, $i;
			}
			next;
		} 
		else {
			# ordinary number
			push @list, $_;
		}
	}
	return @list;
}





__END__

=head1 NAME

graph_profile.pl

A script to graph values along a specific X-axis

=head1 SYNOPSIS

graph_profile.pl <filename> 
   
   --in <filename>
   --index <index1,index2,...>
   --all
   --(no)cen
   --(no)log
   --min <integer>
   --max <integer>
   --x <integer>
   --dir <foldername>
   --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a previously generated feature dataset.
The tim data format is preferable, although any other tab-delimited text 
data formats may be usable. See the file description in L<tim_db_helper.pm>.

=item --index <index>

Specify the column number(s) corresponding to the dataset(s) in
the file to graph. Number is 0-based index. Each dataset should be 
demarcated by a comma. A range of indices may also be specified using 
a dash to demarcate the beginning and end of the inclusive range. 
Multiple datasets to be graphed together should be joined with an ampersand. 
For example, "2,4-6,5&6" will individually graph datasets 2, 4, 5, 6, 
and a combination 5 and 6 graph.

If no dataset indices are specified, then they may be chosen 
interactively from a list.

=item --(no)log

Dataset values are (not) in log2 space, or status should be respected 
if indicated in the file metadata.

=item --cen

Datasets should be median centered prior to graphing. Useful when 
graphing multiple datasets together when they have different 
medians.

=item --min <integer>, --max <integer>

Specify the minimum and/or maximum values for the Y-axis. The default 
values are automatically determined from the dataset.

=item --x <index>

Specify the index of the X-axis dataset. Unless specified, the program 
automatically uses the columns labeled 'Midpoint' or 'Window', if 
present. 

=item --dir

Optionally specify the name of the target directory to place the 
graphs. The default value is the basename of the input file 
appended with "_graphs".

=item --help

Print this help documenation

=back

=head1 DESCRIPTION

This program will generate PNG graphic files representing a profile of 
values plotted along a specific X-axis. For example, plotting values 
along genomic coordinates or relative positions along a feature. The 
X-axis values are static and each dataset is plotted against it. One or 
more datasets may be plotted on a single graph, each in a different 
color, with a legend printed beneath the graph. The graph is a simple 
Bezier-smoothed line graph.

The resulting files are written to a subdirectory named after 
the input file. The files are named after the dataset name (column 
header) with a prefix. 

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112





