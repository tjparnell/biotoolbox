#!/usr/bin/perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long;
use Statistics::Lite qw(:all);
use Bio::ToolBox::Data;
use Bio::ToolBox::utility;
my $VERSION = 1.23;

print "\n A tool for manipulating datasets in data files\n";



### Quick help
unless (@ARGV) { # when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Get command line options and initialize values
my ( # command line option variables
	$infile, 
	$outfile, 
	$function, 
	$opt_index,
	$opt_numerator, 
	$opt_denominator, 
	$opt_target, 
	$opt_placement,
	$opt_exception, 
	$opt_zero, 
	$opt_direction, 
	$opt_name, 
	$opt_log, 
	$gz, 
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # name of input file
	'out=s'     => \$outfile, # name of new output file 
	'func=s'    => \$function, # name of the function to  perform
	'index=s'   => \$opt_index, # index number(s) of the dataset to work on
	'exp|num=i' => \$opt_numerator, # index number of numerator dataset
	'con|den=i' => \$opt_denominator, # index number of denominator dataset
	'target=s'  => \$opt_target, # target
	'place=s'   => \$opt_placement, # placement of transformed dataset
	'except=s'  => \$opt_exception, # old argument exception to deal with 0 values
	'zero!'     => \$opt_zero, # include 0 values
	'dir=s'     => \$opt_direction, # sort order
	'name=s'    => \$opt_name, # new dataset name
	'log!'      => \$opt_log, # data values are in log2 space
	'gz!'       => \$gz, # write gzipped data file
	'help'      => \$help, # request help
	'version'   => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Get file name
unless ($infile) {
	$infile = shift @ARGV;
}


### Print help if requested
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script manipulate_datasets.pl, version $VERSION\n\n";
	exit;
}


### Check for required values
unless ($infile) {
	die "No file name supplied!";
}



### Load file
my $Data = Bio::ToolBox::Data->new(file => $infile) 
	or die "no data loaded from file '$infile'!";
printf "    Loaded '$infile' with %s data rows and %s columns\n\n", 
	format_with_commas($Data->last_row), $Data->number_columns;



### Initialize more variables
# parse the requested index string in to an array of indices
my @opt_indices;
if (defined $opt_index) {
	@opt_indices = parse_list($opt_index) or 
		die " requested index '$opt_index' cannot be parsed!\n";
}

# set zero usage
if ($opt_exception =~ /^y/i) {
	$opt_zero = 1;
}

# initialize modification count
my $modification = 0; 
	# the function subroutines return a value indicating success (the data
	# table has been modified) or failure (or at least no changes has been
	# written to the data array)
	# the $modification variable stores the sum of these return values and 
	# hence indicates whether a new file needs to be written or not

# set menu hashes
my %letter_to_function = _get_letter_to_function_hash();
	# these are hashes that convert a single letter to the function name

# set the function to subroutine hash
my %function_to_subroutine = _get_function_to_subroutine_hash();
	# this hash will convert the function name to the actual subroutine that 
	# does the work


### Determine execution mode
if ($function) {
	# Function was requested upon execution
	automatic_execution();
	write_and_quit_function();
} 
else {
	# otherwise interactive execution
	interactive_execution();
	write_and_quit_function();
}




################################################################################
##############             Main Program Subroutines               ##############
################################################################################






sub print_menu {
	print <<MENU;
These are the functions available:
 
Column manipulation
  R  (R)eorder columns in a different order
  D  (D)elete a column
  n  Re(n)ame a column
  w  Generate a ne(w) column with a single identical value
  b  Num(b)er the rows
  C  (C)oncatenate columns
  T  spli(T) column into new columns
  O  Generate c(O)ordinate string
Row manipulation
  o  S(o)rt all rows by a specific column
  g  (g)enomic sort by chromosome, start
  N  Delete rows with (N)ull values
  P  Delete rows with du(P)licate values, keeping first
  A  Delete rows with values (A)bove threshold
  B  Delete rows with values (B)elow threshold
  S  Delete rows with (S)pecific values
  K  (K)eep only rows with specific values
Conversions
  U  Convert n(U)ll values to a specific value
  G  Convert si(G)ned values to an absolute value
  I  Set a m(I)nimum value
  X  Set a ma(X)imum value
  l  (l)og convert the column
  L  De-(L)og the column
  f  (f)ormat decimal numbers in a column
  p  (p)ercentile rank convert a column
Mathematical manipulation
  a  (a)dd a specific value to a column
  u  S(u)btract a specific value from a column
  y  Multipl(y) a specific value with a column
  v  Di(v)ide a column by a specific value 
  c  (c)ombine columns with math operation
  s  Median (s)cale a column
  Z  Generate (Z)-score values of a column
  r  Generate a (r)atio between two columns
  d  Generate a (d)ifference between two columns
  z  Generate a normali(z)ed difference between two columns
  e  Median c(e)nter row values
File
  W  Re(W)rite the file
  x  E(x)port into a simple tab-delimited text file
  i  Export a .cdt file for Treev(i)ew or Cluster analysis
  Y  Write out a mean summar(Y) profile of the data
Other
  t  Print S(t)atistics on a column
  V  (V)iew table contents
  h  (h)elp
  q  (q)uit, saving changes if necessary
  Q  (Q)uit without saving changes
MENU
	# unlisted option: print this (m)enu
	# unused letters: E F H jJ k M 
	return; # return 0, nothing done
}


sub print_online_help {
	# display FUNCTIONS from POD for online help
	# upon quitting perldoc it should return to the program
	pod2usage( {
		'-verbose'  => 99,
		'-sections' => 'FUNCTIONS',
		'-exitval'  => 'NOEXIT',
	} );
	return;
}



sub automatic_execution {
	# use global variables defined by command-line arguments to automatically
	# execute the manipulations
	
	# print options as a record for when I keep the output
	print " Executing automatic manipulation with function '$function'\n";
	my %option_values = (
		'index'      => $opt_index, 
		'experiment' => $opt_numerator, 
		'control'    => $opt_denominator, 
		'target'     => $opt_target, 
		'place'      => $opt_placement, 
		'include zero' => $opt_zero, 
		'direction'  => $opt_direction, 
		'name'       => $opt_name, 
		'log'        => $opt_log, 
	);
	foreach my $value (sort {$a cmp $b} keys %option_values) {
		if (defined $option_values{$value}) {
			print "   using '$option_values{$value}' for option $value\n";
		}
	}
	
	# set default placement to new
		# I hate having a shell script stall asking the user for placement
	unless (defined $opt_placement) {
		$opt_placement = 'n';
	}
	
	# execute
	if (exists $function_to_subroutine{$function} ) {
		# the function is recognized as legitimate
		$modification += &{ $function_to_subroutine{$function} };
	}
	else {
		die "unknown function '$function'!";
	}
}


sub interactive_execution {
	# interact with the user to perform an unlimited number of executions
	
	# Ask for the function
	print " Functions are chosen by symbol. Use [m] for menu or [h] for help\n";
	print " Enter the symbol for the function you would like to perform   ";
	my $request = <STDIN>;
	chomp $request;
	while (1) {
		if (exists $letter_to_function{$request} ) {
			# first check that the letter corresponds to a function
			
			# perform the function
			$modification += &{ 
				$function_to_subroutine{
					$letter_to_function{$request}
				}
			};
			
			# prepare for the next function
			print " Which function is next? [m]enu, [h]elp, just [Q]uit, save & [q]uit  ";
			$request = <STDIN>;
			chomp $request;
		}
		
		else {
			# unrecognized command
			print " unrecognized command. [m]enu, [h]elp, just [Q]uit, save & [q]uit  ";
			$request = <STDIN>;
			chomp $request;
		}
	}
	return;
}



################################################################################
##############               Function Subroutines                 ##############
################################################################################


sub write_and_quit_function {
	# write out file as necessary and quit
	
	if ($modification > 0) { 
		# a value greater than 0 indicates that changes to the data array have been
		# made and that we need to write an output file
	
		# write the file
		my $write_results = $Data->write_file(
			'filename'  => $outfile,
			'gz'        => $gz,
		);
	
		# report write results
		print " $modification manipulations performed\n";
		if ($write_results) {
			print " Wrote datafile $write_results\n";
		}
		else {
			print " Failed to write datafile\n";
		}
	} 
	else {
		# no need to write output
		print " No changes written\n";
	}
	exit;
}


sub quit_function {
	# do not write file
	print " No changes written\n";
	exit;
}


sub print_statistics_function {
	# print simple statistics for the dataset
	
	# request dataset(s)
	my @indices = _request_indices(
		" Enter one or more column index numbers to calculate statistics  "
	);
	unless (@indices) {
		warn " unknown index number(s).\n";
		return;
	}
	
	# get statistics and print
	foreach my $index (@indices) {
		my %statdata = _get_statistics_hash($index);
		unless (%statdata) { 
			warn " unable to get statistics for column index $index!\n"; 
			return;
		}

		# print the metadata and the calculated statistics
		printf "  Statistics for column $index, '%s'\n", $Data->name($index);
		my %metadata = $Data->metadata($index);
		foreach my $key (keys %metadata ) {
			# print the metadata
			if ($key eq 'name') {
				next; # skip the name, it's been done
			}
			elsif ($key eq 'index') {
				next; # skip the index, it's been done
			}
			elsif ($key eq 'AUTO') {
				next; # skip the name, it's been done
			}
			else {
				printf "   $key => '%s'\n", $metadata{$key};
			}
		}
		print "   count    = $statdata{count}\n" ;
		print "   mean     = $statdata{mean}\n", ;
		print "   median   = $statdata{median}\n", ;
		print "   std dev  = $statdata{stddevp}\n", ;
		print "   min      = $statdata{min}\n", ;
		print "   max      = $statdata{max}\n";
		print "   mode     = $statdata{mode}\n";
		print "   sum      = $statdata{sum}\n";
	}
	
	# we're returning 0 because the data file has not been modified
	return; 
}


sub reorder_function {
	# this subroutine will re-order the datasets (columns) in the file
	
	# determine the new order for the columns
	my @order; # array for the new order
	my $sub_request = 0; # a flag to indicate a request from another subroutine
	if (@_) {
		# the new order may be passed from another subroutine
		@order = @_;
		$sub_request = 1; # set to true
	}
	else {
		# otherwise request from user
		my $line = 
			" Enter the indices in the desired order." . 
			" Indices may skipped or duplicated.\n" . 
			" Enter as comma delimited list, and/or range (start - stop)\n   ";
		@order = _request_indices($line);
	}
	unless (@order) {
		warn " No order! Nothing done!\n";
		return;
	}
	
	# re-order the data columns
	$Data->reorder_column(@order);
	
	# completion
	if ($sub_request) {
		# the re-ordering was requested by another subroutine
		# suppress the completion statement
		return 1;
	}
	else {
		# explicit user request, print completion statement
		print " re-ordered data as '" . join(", ", @order) . "\n";
		return 1;
	}
}


sub delete_function {
	# this will delete datasets (columns)
	
	# request the datasets to be deleted
	my @deletion_list;
	if (@_) {
		# the deletion may be passed on from another subroutine
		@deletion_list = @_;
	}
	else {
		# otherwise request list from user
		@deletion_list = _request_indices(
		" Enter one or more column index numbers to be deleted.\n  "
		);
	}
	unless (@deletion_list) {
		warn " No list for deletion! Nothing done!\n";
		return;
	}
	
	$Data->delete_column(@deletion_list);
	print " datasets '" . join(", ", @deletion_list) . "' deleted\n";
	return 1;
}


sub rename_function {
	# this subroutine will re-name a dataset name
	
	# determine or request dataset index and newname
	my ($index, $newname);
	if (scalar @_ == 2) {
		# passed from internal subroutine
		$index   = $_[0];
		$newname = $_[1];
	}
	else {
		# request from user
		
		# index
		$index = _request_index(
			" Enter the index number of the column to rename  "
		);
		if ($index == -1) {
			warn " unknown index number. nothing done\n";
			return;
		}
		
		# name
		if ($function and $opt_name) {
			# new name is specified from the command line during automatic execution
			# use this global value
			$newname = $opt_name;
		}
		else {
			# request a new name from the user
			print " Enter a new name...  ";
			$newname = <STDIN>;
			chomp $newname;
		}
	}
	
	# assign new name
	my $oldname = $Data->name($index);
	$Data->name($index, $newname);
	print " $oldname re-named to $newname\n";
	return 1;
}


sub concatenate_function {
	# this subroutine will concatenate two or more column values into a new column
	
	# request datasets
	my @indices = _request_indices(
			" Enter two or more column index numbers to concatenate  "
	);
	unless (@indices) {
		warn " unknown index number(s). nothing done\n";
		return;
	}
	unless (scalar @indices >= 2) {
		warn " at least two columns are required to concatenate! nothing done\n";
		return;
	}
	
	# identify joining character
	my $glue;
	if ($opt_target) {
		$glue = $opt_target;
	}
	elsif ($function) {
		# automatic execution, don't bother user, default is _
		print " Joining values with character '_'\n";
		$glue = '_';
	}
	else {
		print " Enter a character to join values, return for nothing    ";
		$glue = <STDIN>;
		chomp $glue;
	}
	
	# generate new column
	my $new_name;
	if ($function and $opt_name) {
		# automatic execution and new name was specifically given 
		$new_name = $opt_name;
	}
	else {
		$new_name = join('_', map {$Data->name($_)} @indices);
		if (length($new_name) > 30) {
			# I don't like long names!
			$new_name = 'cat_' . join(',', @indices);
		}
	}
	my $new_position = $Data->add_column($new_name);
	
	# concatenate values
	$Data->iterate( sub {
		my $row = shift;
		my $cat = join($glue, map {$row->value($_)} @indices);
		$row->value($new_position, $cat);
	} );
	
	printf " Values from %s concatenated into new column '$new_name'\n", 
		join(", ", map {$Data->name($_)} @indices);
	return 1;
}


sub split_function {
	# split a column into multiple new columns
	# this is tricky as we must know how many columns to create
	
	# Request dataset
	my $index = _request_index(
		" Enter the index number of a column to split  ");
	if ($index == -1) {
		warn " unknown index number. nothing done\n";
		return;
	}
	
	# identify joining character
	my $glue;
	if ($opt_target) {
		$glue = $opt_target;
	}
	elsif ($function) {
		# automatic execution, don't bother user, default is _
		print " Splitting on default character '_'\n";
		$glue = '_';
	}
	else {
		print " Enter the character to split the values, return for nothing    ";
		$glue = <STDIN>;
		chomp $glue;
	}
	
	# add new column
	my $new_position = $Data->add_column( $Data->name($index) . '_array' );
	
	# split values
	my $number = 1; # maximum number of split values found
	$Data->iterate( sub {
		my $row = shift;
		my @values = split($glue, $row->value($index));
		if (scalar @values > $number) {
			$number = scalar @values;
		}
		$row->value($new_position, \@values); # store array reference for now
	} );
	
	# add more columns as necessary
	my @positions;
	for my $i (1 .. $number) {
		my $p = $Data->add_column( $Data->name($index) . '_' . $i );
		push @positions, $p;
	}
	
	# now store the split values from the array ref into the new columns
	$Data->iterate( sub {
		my $row = shift;
		my $values = $row->value($new_position); # an array ref
		my $i = 0;
		foreach my $p (@positions) {
			$row->value($p, $values->[$i] || '.');
			$i++;
		}
	} );
	
	# delete the array reference Column
	$Data->delete_column($new_position);
	
	printf " Split column %s into $number new columns\n", $Data->name($index);
	return 1;
}


sub coordinate_function {
	# this subroutine will generate a coordinate string from coordinate values
	
	# identify the coordinates
	my $chr_i   = $Data->chromo_column;
	my $start_i = $Data->start_column;
	my $stop_i  = $Data->stop_column;
	unless (defined $chr_i and defined $start_i) {
		# cannot add coordinate column, do without ?
		warn " cannot generate coordinates, no chromosome or start column found\n";
		return;
	}
	
	# the new index position 
	my $new_position = $Data->add_column('Coordinate');
	
	# generate coordinates
	if (defined $stop_i) {
		# we have a stop coordinate to use
		$Data->iterate( sub {
			my $row = shift;
			my $coord = join("", $row->seq_id, ':', $row->start, '-', $row->end);
			$row->value($new_position, $coord);
		} );
	}
	else {
		# we don't have a stop coordinate to use
		$Data->iterate( sub {
			my $row = shift;
			my $coord = join(':', $row->seq_id, $row->start);
			$row->value($new_position, $coord);
		} );
	}
	
	print " Coordinate string generated as new column $new_position\n";
	return 1;
}


sub median_scale_function {
	# this subroutine will median scale a dataset
	
	# request datasets
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter one or more column index numbers to median scale  "
		);
	}
	unless (@indices) {
		warn " unknown index number(s). nothing done\n";
		return;
	}
	
	# Where to put new values?
	my $placement = _request_placement();
	
	# Obtain the target median value
	my $target;
	if (defined $opt_target) {
		# use the command line specified target
		$target = $opt_target;
	}
	else {
		# request target from user
		print " Enter the new median target  ";
		$target = <STDIN>;
		chomp $target;
	}
	
	# Work through the requested datasets
	my @datasets_modified; # a list of which datasets were modified
	INDEX_LOOP: foreach my $index (@indices) {
		
		# Retrieve values and calculate median
		my %statdata = _get_statistics_hash($index, 'n');
		unless (%statdata) { 
			warn " unable to get statistics for dataset " . 
				$Data->name($index) . ", index $index!\n"; 
			next INDEX_LOOP;
		}
		print " The median value for dataset " . 
			$Data->name($index) . " is $statdata{'median'}\n";
	
		# Calculate correction value
		my $correction_value = $target / $statdata{median};
	
		# Replace values
		$index = _prepare_new_destination($index, '_scaled') if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			next if $row->value($index) eq '.'; # null value, nothing to do
			my $v = $correction_value * $row->value($index);
			$row->value($index, $v);
		} );
			
		# annotate metadata
		$Data->metadata($index, 'median_scaled', $target);
		
		# results
		push @datasets_modified, $Data->name($index);
	}
	
	# report results
	if (@datasets_modified) {
		printf " %s were median scaled to $target\n", 
			join(", ", @datasets_modified);
	}
	return scalar(@datasets_modified);
}


sub percentile_rank_function {
	# this subroutine will convert a dataset into a percentile rank

	# request datasets
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter one or more column index numbers to convert to percentile rank  "
		);
	}
	unless (@indices) {
		warn " unknown index number(s). nothing done\n";
		return;
	}
	
	# Where to put new values?
	my $placement = _request_placement();	
	
	# Process each index request
	my @datasets_modified; # a list of which datasets were modified
	INDEX_LOOP: foreach my $index (@indices) {
		
		# Calculate percent rank of values
			# remove null values
		my @values = grep {!/^\.$/} $Data->column_values($index);
		my $total = scalar @values;
		my %percentrank;
		my $n = 1;
		foreach (sort { $a <=> $b } @values) {
			# sort by increasing hash values, not hash keys
			# percentrank is key value (index) divided by total
			$percentrank{$_} = $n/$total;
			$n++;
		}
		
		# Replace the contents with the calculated percent rank
		$index = _prepare_new_destination($index, '_pr') if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			next if $row->value($index eq '.');
			$row->value($index, $percentrank{ $row->value($index) });
		} );
		
		# update metadata
		$Data->metadata($index, 'converted', 'percent_rank');
		
		# done
		push @datasets_modified, $Data->name($index);
	}	
	
	# report results
	if (@datasets_modified) {
		printf " %s were converted to percent rank\n", 
			join(", ", @datasets_modified);
	}
	return scalar(@datasets_modified);
}


sub zscore_function {
	# this subroutine will generate a z-score for each value in a dataset

	# identify the datasets to convert
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter one or more column index numbers to convert to z-scores  "
		);
	}
	unless (@indices) {
		warn " Unknown columns. Nothing done.\n";
		return;
	}
	
	# Where to put new values?
	my $placement = _request_placement();	
	
	# Process each index request
	my @datasets_modified; # a list of which datasets were modified
	foreach my $index (@indices) {
		
		# generate statistics on the dataset
		my %statdata = _get_statistics_hash($index, 'y');
		unless (%statdata) {
			warn " unable to generate statistics for index $index! skipping\n";
			next;
		}
		
		# Replace the current values
		$index = _prepare_new_destination($index, '_Zscore') if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			next if $row->value($index) eq '.';
			my $v = ($row->value($index) - $statdata{'mean'}) / $statdata{'stddevp'};
			$row->value($index, $v);
		} );
		
		# update metadata
		$Data->metadata($index, 'converted', 'Z-score');
		
		# done
		push @datasets_modified, $Data->name($index);
	}
	
	# report results
	if (@datasets_modified) {
		printf " %s were converted to Z-scores\n", 
			join(", ", @datasets_modified);
	}
	return scalar(@datasets_modified);
}




sub sort_function {
	# This will sort the entire data table by the values in one dataset
	
	# Request dataset
	my $index;
	if (@_) {
		# from another subroutine
		$index = shift @_;
	}
	else {
		$index = _request_index(
		" Enter the index number of a column to sort by  ");
	}
	if ($index == -1) {
		warn " unknown index number. nothing done\n";
		return;
	}
	
	# Ask the sort direction
	my $direction;
	if ($opt_direction) {
		# direction was specified on the command line
		$direction = $opt_direction;
	}
	else {
		# otherwise ask the user for the direction
		print " Sort by (i)ncreasing or (d)ecreasing order?  ";
		$direction = <STDIN>;
		chomp $direction;
		unless ($direction =~ /^i|d$/i) {
			warn " unknown order; nothing done\n";
			return;
		}
	}
	
	# sort
	$Data->sort_data($index, $direction);
	
	# remove any pre-existing sorted metadata since no longer valid
	for (my $i = 0; $i < $Data->number_columns; $i++) {
		$Data->delete_metadata($i, 'sorted');
	}
	
	# annotate metadata
	if ($direction =~ /i/i) {
		$Data->metadata($index, 'sorted', "increasing")
			unless $Data->metadata($index, 'AUTO'); # internal flag to not accept metadata
	}
	else {
		$Data->metadata($index, 'sorted', "decreasing")
			unless $Data->metadata($index, 'AUTO'); # internal flag to not accept metadata
	}
	
	return 1;
}



sub genomic_sort_function {
	# This will sort the entire data table by chromosome and start position
	
	my $result = $Data->gsort_data;
	unless ($result) {
		print " Data table not sorted\n";
		return;
	}
	
	# remove any pre-existing sorted metadata since no longer valid
	for (my $i = 0; $i < $Data->number_columns; $i++) {
		$Data->delete_metadata($i, 'sorted');
	}
	
	# annotate metadata
	my $chr_i = $Data->chromo_column;
	my $start_i = $Data->start_column;
	$Data->metadata($chr_i, 'sorted', 'genomic') unless 
		$Data->metadata($chr_i, 'AUTO');
	$Data->metadata($start_i, 'sorted', 'genomic') unless 
		$Data->metadata($start_i, 'AUTO');
	
	print " Data table is sorted by genomic order\n";
	return 1;
}



sub toss_nulls_function {
	# Toss out datapoints (lines) that have a non-value in the specified dataset
	
	# generate the list of datasets to check
	my @order = _request_indices(
		" Enter one or more column index numbers to check for non-values\n   "); 
	unless (@order) {
		warn " No valid columns! Nothing done!\n";
		return;
	}
	
	# Collection exception rule from commandline
	my $zero = $opt_zero;
	
	# Identify those rows that need to be deleted
	my @todelete;
	$Data->iterate( sub {
		my $row = shift;
		my $check = 0; 
		foreach my $i (@order) {
			my $v = $row->value($i);
			if ($v eq '.') {
				$check++;
			} 
			elsif (not defined $v) {
				$check++;
			} 
			elsif ($v == 0) {
				# we have a 0 value, what to do?
				if (not defined $zero and $function) {
					# running automatically, do not both user
					$zero = 0;
				}
				elsif (not defined $zero and not defined $function) {
					# ask the user if we haven't already
					print " Also toss values of 0? y or n  ";
					my $answer = <STDIN>;
					if ($answer =~ /^y/i) {
						$zero = 1;
					}
					else {
						$zero = 0;
					}
				}
				if ($zero) { 
					# Yes, toss the 0 values
					$check++;
				}
			}
		}
		# mark for deletion if the row fails the check
		push @todelete, $row->row_index if $check;
	} );
	
	# Delete
	$Data->delete_row(@todelete);
	
	# update metadata
	foreach my $index (@order) {
		$Data->metadata($index, 'deleted_non_value_features',  scalar(@todelete))
			unless $Data->metadata($index, 'AUTO');
	}
	
	# report
	printf " %s rows with null values in %s were deleted.\n", scalar(@todelete), 
		join(', ', map {$Data->name($_)} @order);
	printf " %s rows are remaining\n", $Data->last_row;
	return 1;
}






sub toss_duplicates_function {
	# Toss out datapoints (lines) that have a duplicate values
	
	# generate the list of datasets to check
	my @order = _request_indices(
		" Enter one or more column index numbers to check for duplicates\n   "); 
	unless (@order) {
		warn " No valid columns! Nothing done!\n";
		return;
	}
	
	# initialize variables
	my %values2check;
	my @todelete;
	
	# check values
	$Data->iterate( sub {
		my $row = shift;
		# we will simply concatenate all values to check for duplicity
		my $value = join('_', map {$row->value($_)} @order);
		if (exists $values2check{$value}) {
			# yes, it exists, mark for destruction
			$values2check{$value}++;
			push @todelete, $row->row_index;
		}
		else {
			# nope, it's good
			$values2check{$value} = 1;
		}
	} );
	
	# delete
	$Data->delete_row(@todelete);
	
	# update metadata
	foreach my $index (@order) {
		$Data->metadata($index, 'deleted_duplicate_features',  scalar(@todelete))
			unless $Data->metadata($index, 'AUTO');
	}
	
	# print result
	printf " %s rows with duplicate values in %s were deleted.\n", 
		scalar(@todelete), join(', ', map {$Data->name($_)} @order);
	printf " %s rows are remaining\n", $Data->last_row;
	return 1;
}

sub toss_above_threshold_function {
	# Toss out datapoints (lines) whose value is above a certain threshold
	return toss_threshold_function('above');
}


sub toss_below_threshold_function {
	# Toss out datapoints (lines) whose value is below a certain threshold
	return toss_threshold_function('below');
}


sub toss_threshold_function {
	# Toss out lines whose specified datapoint value is above or below a 
	# specified threshold
	
	# get the direction for excluding datapoints, above or below
	my $direction = shift;
	
	# generate the list of datasets to check
	my @order = _request_indices(
		" Enter one or more column index numbers to toss values $direction a value\n   "); 
	unless (@order) {
		warn " No valid datasets! Nothing done!\n";
		return;
	}
	
	# identify the threshold
	my $threshold;
	if (defined $opt_target) {
		# specified on the command line
		$threshold = $opt_target;
	}
	else {
		# interactively ask the user
		if ($direction eq 'above') {
			print " Toss values that exceed this value:  ";
		}
		else {
			print " Toss values that are below this value:  ";
		}
		$threshold = <STDIN>;
		chomp $threshold;
	}
	
	# Check values
	my @todelete;
	if ($direction eq 'above') {
		$Data->iterate( sub {
			my $row = shift;
			my $check = 0; 
			foreach my $i (@order) {
				my $v = $row->value($i);
				next unless defined $v;
				next if $v eq '.';
				$check++ if $v > $threshold;
			}
			# mark for deletion if the row fails the check
			push @todelete, $row->row_index if $check;
		} );
	}
	elsif ($direction eq 'below') {
		$Data->iterate( sub {
			my $row = shift;
			my $check = 0; 
			foreach my $i (@order) {
				my $v = $row->value($i);
				next unless defined $v;
				next if $v eq '.';
				$check++ if $v < $threshold;
			}
			# mark for deletion if the row fails the check
			push @todelete, $row->row_index if $check;
		} );
	}
	
	# Delete
	$Data->delete_row(@todelete);
	
	# update metadata
	foreach my $index (@order) {
		$Data->metadata($index, "deleted_$direction\_$threshold", scalar(@todelete) ) 
			unless $Data->metadata($index, 'AUTO');
	}
	
	# report
	printf " %s rows with values $direction $threshold in %s were deleted.\n", 
		scalar(@todelete), join(', ', map {$Data->name($_)} @order);
	printf " %s rows are remaining\n", $Data->last_row;
	return 1;
}


sub toss_specific_values_function {
	return do_specific_values_function(1);
}


sub keep_specific_values_function {
	return do_specific_values_function(0);
}


sub do_specific_values_function {
	# toss specific values
	my $toss = shift; # boolean to toss (1) or keep (0) 
	
	# generate the list of datasets to check
	my @list = _request_indices(
		" Enter one or more column index numbers to check for specific values\n   "); 
	unless (@list) {
		warn " No valid columns! Nothing done!\n";
		return;
	}
	
	# determine values
	my %wanted;
	if ($opt_target) {
		# list provided
		%wanted = map {$_ => 0} split(/,/, $opt_target);
	}
	else {
		# generate potential values, hope there aren't too many
		my %possibilities;
		$Data->iterate( sub {
			my $row = shift;
			foreach my $l (@list) {
				$possibilities{ $row->value($l) } += 1;
			}
		} );
		
		# present list to user
		print " These are the values (occurrences) in the indicated columns\n";
		my %lookup;
		my $i = 1;
		foreach (sort { $a cmp $b } keys %possibilities) {
			# sort asciibetically by name
			printf "   $i\t$_ (%s)\n", $possibilities{$_};
			$lookup{$i} = $_;
			$i++;
		}
		printf " Enter the number(s) corresponding to the values to %s. Enter as a \n" . 
			" comma-delimited list or range.    ", $toss ? 'toss' : 'keep';
		my $response = <STDIN>;
		chomp $response;
		foreach (parse_list($response)) {
			if (exists $lookup{$_}) {
				$wanted{ $lookup{$_} } = 0;
			}
		}
	}
	unless (%wanted) {
		warn " No specific values provided! Nothing done!\n";
		return;
	}
	
	# Identify rows to delete
	my @todelete;
	if ($toss) {
		# we are tossing lines that contain the specific value
		$Data->iterate( sub {
			my $row = shift;
			foreach my $l (@list) {
				if (exists $wanted{ $row->value($l) }) {
					push @todelete, $row->row_index;
					last;
				}
			}
		} );
	}
	else {
		# we are keeping only those lines that contain the specific value
		$Data->iterate( sub {
			my $row = shift;
			my $check = 1; # default is to delete this row, so check starts true
			foreach my $l (@list) {
				if (exists $wanted{ $row->value($l) }) {
					$check = 0;
					last;
				}
			}
			push @todelete, $row->row_index if $check;
		} );
	}
	
	# Delete
	$Data->delete_row(@todelete);
	
	# update metadata
	foreach my $index (@list) {
		if ($toss) {
			$Data->metadata($index, "deleted_specific_values", join(',', keys %wanted) ) 
				unless $Data->metadata($index, 'AUTO');
		}
		else {
			$Data->metadata($index, "kept_specific_values", join(',', keys %wanted) ) 
				unless $Data->metadata($index, 'AUTO');
		}
	}
	
	# report
	printf " %s rows with specific values in %s were deleted.\n", scalar(@todelete), 
		join(', ', map {$Data->name($_)} @list);
	printf " %s data lines are remaining\n", $Data->last_row;
	return 1;
}



sub convert_nulls_function {
	# Convert null values to something else
	
	# identify the datasets to check
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter one or more dataset index numbers to convert null values  "
		);
	}
	unless (@indices) {
		warn " no valid indices. Nothing done.\n";
		return;
	}
	
	# request replacement value
	my $new_value;
	if (defined $opt_target) {
		# command line option
		$new_value = $opt_target;
	}
	else {
		# interactively ask the user
		print " Enter the new value to convert nulls to  ";
		$new_value = <STDIN>;
		chomp $new_value;
	}
	
	# check zero status
	my $zero;
	if (defined $opt_zero) {
		# command line option
		$zero = $opt_zero;
	}
	
	# request placement
	my $placement = _request_placement();
	
	## Process the datasets and subtract their values
	my @datasets_modified; # a list of which datasets were modified
	my $total_count = 0; # total number of resets done
	foreach my $index (@indices) {
		
		# number of resets we do for this index
		my $count  = 0; 
		
		# reset values
		$index = _prepare_new_destination($index, '_convert_nulls') if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			my $v = $row->value($index);
			if (not defined $v) {
				$row->value($index, $new_value);
				$count++;
			}
			elsif ($v eq '.') {
				$row->value($index, $new_value);
				$count++;
			}
			elsif ($v == 0) {
				# zero value, what to do?
				if (not defined $zero and $function) {
					# running automatically, do not both user
					$zero = 0;
				}
				elsif (not defined $zero and not defined $function) {
					# wasn't defined on the command line, running interactively, 
					# so stop the program and ask the user
					print " Include 0 values to convert? y or n  ";
					my $answer = <STDIN>;
					$zero = ($answer =~ /^y/i) ? 1 : 0;
					# remember for next time, chances are user may still want this
					# value again in the future
					$opt_zero = $zero;
				}
				
				if ($zero) {
					$row->value($index, $new_value);
					$count++;
				}
			}
		} );
		
		# update metadata
		if ($count) {
			$Data->metadata($index, 'null_value', $new_value) unless 
				$Data->metadata($index, 'AUTO');
			$total_count += $count;
			push @datasets_modified, $Data->name($index);
		}
	}
		
	
	# report results
	if (@datasets_modified) {
		printf " $total_count null values were converted for %s\n", 
			join(", ", @datasets_modified);
	}
	return scalar(@datasets_modified);
}



sub convert_absolute_function {
	# Convert signed values to their absolute value
	
	# identify the datasets to check
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter one or more dataset index numbers to make absolute  "
		);
	}
	unless (@indices) {
		warn " no valid indices. Nothing done.\n";
		return;
	}
	
	# request placement
	my $placement = _request_placement();
	
	## Process the datasets and subtract their values
	my @datasets_modified; # a list of which datasets were modified
	my $total_count = 0; # total number of conversions done
	my $total_failed = 0;
	foreach my $index (@indices) {
		
		# number of resets we do
		my $count  = 0; 
		my $failed = 0;
		
		# reset minimum values
		$index = _prepare_new_destination($index, '_absolute') if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			my $v = $row->value($index);
			next if $v eq '.';
			my $new_value;
			eval { $new_value = abs($v) };
			if (defined $new_value) {
				$row->value($index, $new_value);
				$count++;
			}
			else {
				$failed++;
			}
		} );
		
		# update metadata
			
		# results
		if ($count) {
			$Data->metadata($index, 'convert', 'absolute') unless 
				$Data->metadata($index, 'AUTO');
			$total_count += $count;
			push @datasets_modified, $Data->name($index);
		}
		else {
			$total_failed += $failed;
		} 
	}
	
	# report results
	if (@datasets_modified) {
		printf " $total_count values were converted to absolute values for %s\n", 
			join(", ", @datasets_modified);
	}
	if ($total_failed) {
		print " $total_failed values could not be converted\n";
	}
	return scalar(@datasets_modified);
}




sub minimum_function {
	# Set a minimum value 
	
	# request datasets
	my @indices = _request_indices(
		" Enter one or more dataset index numbers to reset minimum values  ");
	unless (@indices) {
		warn " no valid indices. nothing done\n";
		return;
	}
	
	# request value
	my $value;
	if (defined $opt_target) {
		# command line option
		$value = $opt_target;
	}
	else {
		# interactively ask the user
		print " Enter the minimum value to accept  ";
		$value = <STDIN>;
		chomp $value;
	}
	
	# request placement
	my $placement = _request_placement();
	
	## Process the datasets and subtract their values
	my @datasets_modified; # a list of which datasets were modified
	my $total_count = 0; # total number of conversions done
	foreach my $index (@indices) {
		
		# number of resets we do
		my $count  = 0; 
		
		# reset minimum values
		$index = _prepare_new_destination($index, '_minimum_reset') if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			my $v = $row->value($index);
			next if $v eq '.';
			if ($v < $value) {
				$row->value($index, $value);
				$count++;
			}
		} );
		
		# results
		if ($count) {
			$Data->metadata($index, 'minimum_value', $value) unless 
				$Data->metadata($index, 'AUTO');
			$total_count += $count;
			push @datasets_modified, $Data->name($index);
		}
	}
	
	# report results
	if (@datasets_modified) {
		printf " $total_count values were reset to a minimum value for %s\n", 
			join(", ", @datasets_modified);
	}
	return scalar(@datasets_modified);
}




sub maximum_function {
	# Set a maximum value
	
	# request datasets
	my @indices = _request_indices(
		" Enter one or more dataset index numbers to reset maximum values  ");
	unless (@indices) {
		warn " no valid indices. nothing done\n";
		return;
	}
	
	# request value
	my $value;
	if (defined $opt_target) {
		# command line option
		$value = $opt_target;
	}
	else {
		# interactively ask the user
		print " Enter the maximum value to accept  ";
		$value = <STDIN>;
		chomp $value;
	}
	
	# request placement
	my $placement = _request_placement();
	
	
	## Process the datasets and subtract their values
	my @datasets_modified; # a list of which datasets were modified
	my $total_count = 0; # total number of conversions done
	foreach my $index (@indices) {
		
		# number of resets we do
		my $count  = 0; 
		
		# reset minimum values
		$index = _prepare_new_destination($index, '_maximum_reset') if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			my $v = $row->value($index);
			next if $v eq '.';
			if ($v < $value) {
				$row->value($index, $value);
				$count++;
			}
		} );
		
		# results
		if ($count) {
			$Data->metadata($index, 'maximum_value', $value) unless 
				$Data->metadata($index, 'AUTO');
			$total_count += $count;
			push @datasets_modified, $Data->name($index);
		}
	}
	
	# report results
	if (@datasets_modified) {
		printf " $total_count values were reset to a maximum value for %s\n", 
			join(", ", @datasets_modified);
	}
	return scalar(@datasets_modified);
}




sub log_function {
	# this subroutine will convert dataset values to log2 space
	
	# request datasets
	my @indices;
	my $base;
	if (@_) {
		# provided from an internal subroutine
		($base, @indices) = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter one or more dataset index numbers to convert to log  "
		);
		if (defined $opt_target) {
			# specified on the command line
			$base = $opt_target;
		}
		else {
			# interactively ask the user
			print " What log base to use? [2 10]:  ";
			$base = <STDIN>;
			chomp $base;
		}
		
	}
	unless (@indices) {
		warn " unknown index number(s). nothing done\n";
		return;
	}
	unless ($base =~ /^\d+$/) {
		warn " unrecognized base number. nothing done\n";
		return;
	}
	my $factor = log($base);
	
	# request placement
	my $placement = _request_placement();
	
	# process each index request
	my @datasets_modified; # a list of which datasets were modified
	my $total_count = 0; # total number of conversions done
	my $total_failed = 0;
	foreach my $index (@indices) {
		
		# check the current metadata status
		my $check = $Data->metadata($index, 'log') || $Data->metadata($index, 'log2') || 0;
		if ($check != 0) {
			warn " dataset $Data->name($index) metadata " .
				"reports it is currently in log scale. Continue? y/n\n";
			my $response = <STDIN>;
			next if $response =~ /n/i;
		}
		
		# Placement dictates method
		my $count = 0; # conversion count
		my $failed = 0;
		
		# perform log conversions
		$index = _prepare_new_destination($index, "_log$base") if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			my $v = $row->value($index);
			# check the value contents and process appropriately
			if ($v == 0) { 
				# cannot take log of 0, change to null
				$row->value($index, '.'); 
				$failed++;
			} 
			elsif ($v eq '.') {
				# a null value, do nothing
				$failed++;
			} 
			else {
				my $new_value = log($v) / $factor;
				$row->value($index, $new_value);
				$count++;
			}
		} );
			
		# update metadata
		$Data->metadata($index, 'log', $base) unless 
			$Data->metadata($index, 'AUTO');
		
		# results
		if ($count) {
			$total_count += $count;
			push @datasets_modified, $Data->name($index);
		}
		$total_failed += $failed;
	}
	
	# report results
	if (@datasets_modified) {
		printf " $total_count values were converted to log$base for %s\n", 
			join(", ", @datasets_modified);
	}
	if ($total_failed) {
		print " $total_failed values could not be converted\n";
	}
	return scalar(@datasets_modified);
}



sub delog_function {
	# this subroutine will convert a dataset from log to normal numbers
	
	# request datasets
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter one or more dataset index numbers to convert from log2  "
		);
	}
	unless (@indices) {
		warn " unknown index number(s). nothing done\n";
		return;
	}
	
	# request placement
	my $placement = _request_placement();
	
	# process each index request
	my $base;
	my @datasets_modified; # a list of which datasets were modified
	my $total_count = 0; # total number of conversions done
	my $total_failed = 0;
	foreach my $index (@indices) {
		
		# check the log metadata status
		$base ||= $Data->metadata($index, 'log') || 0;
		unless ($base) {
			$base = 2 if ($Data->metadata($index, 'log2'));
			$base = 10 if ($Data->metadata($index, 'log10'));
		}
		if ($base == 0) {
			if (defined $opt_target) {
				$base = $opt_target;
			}
			else {
				print " What log base is the data in? [2 10]:  ";
				$base = <STDIN>;
				chomp $base;
			}
		}
		unless ($base =~ /^\d+$/) {
			warn " Unrecognized base integer '$base'. Nothing done.\n";
			return scalar(@datasets_modified);
		}
		
		# Placement dictates method
		my $count = 0; # conversion count
		my $failed = 0;
		$index = _prepare_new_destination($index, "_delog$base") if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			my $v = $row->value($index);
			# check the value contents and process appropriately
			if ($v eq '.') {
				# a null value, do nothing
				$failed++;
			} 
			else {
				my $new_value = $base ** $v;
				$row->value($index, $new_value);
				$count++;
			}
		} );
		
		# update metadata
		$Data->metadata($index, 'log', 0) unless 
			$Data->metadata($index, 'AUTO');
		$Data->delete_metadata($index, 'log2');
		$Data->delete_metadata($index, 'log10');
			
		# results
		if ($count) {
			$total_count += $count;
			push @datasets_modified, $Data->name($index);
		}
		$total_failed += $failed;
	}
	
	# report results
	if (@datasets_modified) {
		printf " $total_count values were converted from log$base for %s\n", 
			join(", ", @datasets_modified);
	}
	if ($total_failed) {
		print " $total_failed values could not be converted\n";
	}
	return scalar(@datasets_modified);
}



sub format_function {
	# this subroutine will format the numbers in a dataset using sprintf
	
	# request datasets
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter one or more dataset index numbers to format  "
		);
	}
	unless (@indices) {
		warn " unknown index number(s). nothing done\n";
		return;
	}
	
	# ask for the number of decimal positions to format to
	my $positions;
	if (defined $opt_target) {
		# specified on the command line
		$positions = $opt_target;
	}
	else {
		# interactively ask the user
		print " Format the numbers to how many decimal positions?  ";
		$positions = <STDIN>;
		chomp $positions;
	}
	unless ($positions =~ /^\d+$/) {
		warn " Unknown number of positions; formatting NOT done\n";
		return;
	}
	my $format_string = '%.' . $positions . 'f';
	
	# ask for placement
	my $placement = _request_placement();
	
	# format each index request
	my @datasets_modified; # a list of which datasets were modified
	foreach my $index (@indices) {
		$index = _prepare_new_destination($index, '_formatted') if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			my $v = $row->value($index);
			if ($v ne '.') {
				$row->value($index, sprintf($format_string, $v));
			}
		});
	
		$Data->metadata($index, 'formatted', $positions);
		push @datasets_modified, $Data->name($index);
	}
	
	# report results
	if (@datasets_modified) {
		printf " formatted values to $positions decimal positions for %s\n", 
			join(", ", @datasets_modified);
	}
	return scalar(@datasets_modified);
}


sub combine_function {
	# mathematically combine two or more datasets
	
	# request datasets
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter two or more dataset index numbers to mathematically combine  "
		);
	}
	unless (@indices) {
		warn " unknown index number(s). nothing done\n";
		return;
	}
	unless (scalar @indices >= 2) {
		warn " two or more indices are required to combine. nothing done.\n";
		return;
	}
	
	# request value
	my $method;
	if (defined $opt_target) {
		# command line option
		$method = $opt_target;
	}
	else {
		# interactively ask the user
		print " Enter the method: [mean median min max stdev sum]  ";
		$method = <STDIN>;
		chomp $method;
	}
	my $combine_method = 
		$method eq 'mean'   ? \&mean : 
		$method eq 'median' ? \&median :
		$method eq 'min'    ? \&min :
		$method eq 'max'    ? \&max :
		$method eq 'stdev'  ? \&stddevp : 
		$method eq 'sum'    ? \&sum :
		undef;
	unless ($combine_method) {
		warn " unknown method. nothing done\n";
		return;
	}
	
	# generate new column
	my $new_name;
	if ($function and $opt_name) {
		# automatic execution and new name was specifically given 
		$new_name = $opt_name;
	}
	else {
		$new_name = $method;
	}
	my $new_position = $Data->add_column($new_name);
		
	# combine datasets
	my $failure_count = 0; # number of failures
	$Data->iterate( sub {
		my $row = shift;
		
		# collect data, throwing out null values
		my @data = grep {!/^\.$/} map {$row->value($_)} @indices;
		
		# combine the data
		my $v;
		if (@data) {
			$v = &{$combine_method}(@data);
		}
		else {
			$v = '.';
			$failure_count++;
		}
		$row->value($new_position, $v);
	} );
	
	# generate new metadata
	$Data->metadata($new_position, 'datasets', join(',', map {$Data->name($_)} @indices));
	$Data->metadata($new_position, 'combine', $method);
	
	# finish
	printf " combined datasets %s by $method\n", 
		join(', ', map { $Data->name($_) } @indices);
	print " $failure_count data rows could not be combined\n" if $failure_count;
	return 1;
}




sub ratio_function {
	# Generate a new ratio between two datasets
	
	# request dataset
	my ($numerator, $denominator);
	if (defined $opt_numerator and defined $opt_denominator) {
		if ( _validate_index_list($opt_numerator, $opt_denominator) ) {
			$numerator = $opt_numerator;
			$denominator = $opt_denominator;
		}
		else {
			warn " unknown index number(s); nothing done\n";
			return;
		}	
	}
	else {
		($numerator, $denominator) = _request_indices(
			" Enter two dataset index numbers for the ratio as " 
			. "'numerator, denominator'\n  ");
		unless (defined $numerator and defined $denominator) {
			warn " unknown index number(s); nothing done\n";
			return;
		}
	}
	
	# check if log2 numbers
	# log2 ratios performed by subtraction, regular number ratios by division
	my $log;
	if (defined $opt_log) {
		# if unable to identify the log2 status in metadata
		# use the command line option variable if set
		$log = $opt_log;
	}
	else {
		# set log status based on mutual metadata flags
		my $numerator_log = $Data->metadata($numerator, 'log2') || 
			$Data->metadata($numerator, 'log');
		my $denominator_log = $Data->metadata($denominator, 'log2') || 
			$Data->metadata($denominator, 'log');
		
		if ($numerator_log and $denominator_log) { 
			if ($numerator_log == $denominator_log) {
				# both set, both equal, great!
				$log = $numerator_log;
			}
			else {
				warn " Columns appear to have different log status!\n" . 
					"  numerator is $numerator_log and denominator is " . 
					" $denominator_log\n  Nothing done.";
				return;
			}
		}
		elsif ($numerator_log or $denominator_log) { 
			warn " Columns appear to have different log status!\n" . 
				"  numerator is $numerator_log and denominator is " . 
				" $denominator_log\n  Nothing done.";
			return;
		}
		else {
			# looks like neither is log state
			$log = 0;
		}
	}
	$log = 2 if $log == 1; # assume log2 status if simply true
		
	# generate new column
	my $new_name;
	if ($function and $opt_name) {
		# automatic execution and new name was specifically given 
		$new_name = $opt_name;
	}
	else {
		$new_name = $Data->name($numerator) . '_' .
			$Data->name($denominator) . '_ratio';
	}
	my $new_position = $Data->add_column($new_name);
		
	# generate ratio
	my $failure_count = 0; # number of failures
	$Data->iterate( sub {
		my $row = shift;
		my $n = $row->value($numerator);
		my $d = $row->value($denominator);
		
		# add new value
		if ($n eq '.' or $d eq '.') {
			# either value is null
			$row->value($new_position, '.');
			$failure_count++;
		}
		elsif ($d == 0 and !$log) {
			# denominator is 0, avoid div by 0 errors, assign null
			$row->value($new_position, '.');
			$failure_count++;
		}
		else {
			# both values good, calculate ratio
			if ($log) {
				# perform a subtraction with log values
				$row->value($new_position, ($n - $d) );
			} 
			else {
				# perform a division with non-log values
				$row->value($new_position, ($n / $d) );
			}
		}
	} );
	
	# annotate the new medadata hash
	$Data->metadata($new_position, 'method', 'ratio');
	$Data->metadata($new_position, 'log', $log);
	$Data->metadata($new_position, 'denominator', $Data->name($denominator));
	$Data->metadata($new_position, 'numerator', $Data->name($numerator));
	
	# print conclusion
	printf " ratio between %s and %s\n generated as new dataset\n", 
		$Data->name($numerator), $Data->name($denominator);
	if ($failure_count) {
		print " $failure_count datapoints could not generate ratios\n";
	}
	return 1;
}

sub difference_function {
	# Generate a new difference between two datasets
	
	# Determine whether the difference should be normalized
	# normalization divides the difference by the square root of the sum, an
	# estimation of standard deviation
	my $normalization = shift; # pass a true value to normalize
	
	# Request datasets
	my ($experiment_index, $control_index); 
	if (defined $opt_numerator and defined $opt_denominator) {
		if ( _validate_index_list($opt_numerator, $opt_denominator) ) {
			# defined as command line options
			$experiment_index = $opt_numerator;
			$control_index = $opt_denominator;
		}
		else {
			warn " unknown index number(s); nothing done\n";
			return;
		}	
	}
	else {
		# request from user
		my $line;
		if ($normalization) {
			$line = " Enter two dataset index numbers to generate the normalized difference\n" 
				. " as 'experiment, control'\n  ";
		}
		else {
			$line = " Enter two dataset index numbers to generate the difference\n" 
				. " as 'experiment, control'\n  ";
		}
		($experiment_index, $control_index) = _request_indices($line);
		unless (defined $experiment_index and defined $control_index) {
			warn " unknown index number(s); nothing done\n";
			return;
		}
	}
	
	# Generate new column
	my $new_name;
	if ($function and $opt_name) {
		# this was an automatically executed function
		# and a new name was specified on the command line
		$new_name = $opt_name;
	} 
	elsif ($normalization) {
		$new_name = $Data->name($experiment_index) . '_' .
			$Data->name($control_index) . '_normdiff';
	}
	else {
		$new_name = $Data->name($experiment_index) . '_' .
			$Data->name($control_index) . '_diff';
	}
	my $new_position = $Data->add_column($new_name);
		
	# Generate difference
	my $failure_count = 0; # number of failures
	$Data->iterate( sub {
		my $row = shift;
		
		# calculate difference
		if ($row->value($experiment_index) eq '.' or 
			$row->value($control_index) eq '.'
		) {
			# can't do anything with nulls, new value will be null
			$row->value($new_position, '.');
			$failure_count++;
		}
		else {
			my $diff = $row->value($experiment_index) - $row->value($control_index);
			if ($normalization) {
				# determine a normalized difference value
				my $sum = $row->value($experiment_index) + $row->value($control_index);
				if ($sum == 0) {
					# avoid pesky divide-by-0 errors, the difference is also 0
					$row->value($new_position, 0);
				}
				else {
					$row->value($new_position, ($diff / sqrt($sum)) );
				}
			}
			else {
				# determine a straight difference value
				$row->value($new_position, $diff);
			}
		}
	} );
	
	
	# Annotate the new medadata hash
	$Data->metadata($new_position, 'method', 
		$normalization ? 'normalized_difference' : 'difference');
	$Data->metadata($new_position, 'experiment', $Data->name($experiment_index) ); 
	$Data->metadata($new_position, 'control', $Data->name($control_index) ); 
	
	# Print conclusion
	printf "%s difference between %s and %s\n generated as new column\n", 
		$normalization ? ' normalized' : '', $Data->name($experiment_index), 
		$Data->name($control_index);
	if ($failure_count) {
		print " $failure_count rows could not generate ratios\n";
	}
	
	return 1;
}


sub normalized_difference_function {
	# Generate a normalized difference between two datasets
	
	# this actually calls the difference_function subroutine
		# pass a true value to force normalization
	if ( difference_function(1) ) {
		return 1;
	}
	else {
		warn " unable to generate a normalized difference\n";
		return;
	}
}



sub add_function {
	# add a specific value to dataset values
	return math_function('add', @_);
}


sub subtract_function {
	# subtract a specific value from dataset values
	return math_function('subtract', @_);
}


sub multiply_function {
	# multiply dataset values by a specific value
	return math_function('multiply', @_);
}


sub divide_function {
	# divide dataset values by a specific value
	return math_function('divide', @_);
}


sub math_function {
	# General function to perform simple math on each datapoint value 
	# in a dataset by a specific value
	
	# determine the mathematical function
		# function is one of add, subtract, multiply, divide
	my $math = shift;
	my $calculate; # a reference to the appropriate function
	if ($math eq 'add') {
		$calculate = sub {
			return $_[0] + $_[1];
		};
	}
	elsif ($math eq 'subtract') {
		$calculate = sub {
			return $_[0] - $_[1];
		};
	}
	elsif ($math eq 'multiply') {
		$calculate = sub {
			return $_[0] * $_[1];
		};
	}
	elsif ($math eq 'divide') {
		$calculate = sub {
			return $_[0] / $_[1];
		};
	}

	
	# request dataset indices
	my $line = " Enter one or more dataset index numbers to $math  ";
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices($line);
	}
	unless (@indices) {
		warn " no valid indices. nothing done\n";
		return;
	}
	
	# request the value to perform the mathematical function 
	my $request_value;
	if (defined $opt_target) {
		# command line option
		$request_value = $opt_target;
	}
	else {
		# interactively ask the user
		print " Enter the numeric value, 'mean', 'median', or 'sum' to $math  ";
		$request_value = <STDIN>;
		chomp $request_value;
	}
	
	
	# request placement
	my $placement = _request_placement();
	
	
	# generate past tense verb
	my $mathed = $math;
	$mathed =~ s/^(multipl)y$/$1ied/;
	$mathed =~ s/^(add)$/$1ed/;
	$mathed =~ s/^(divide)$/$1d/;
	$mathed =~ s/^(subtract)$/$1ed/;
	
	## Process the datasets and subtract their values
	my $dataset_modification_count = 0; # a count of how many processed
	foreach my $index (@indices) {
		
		# determine the actual value to divide by
		my $value; # the actual number to divide by
		if ($request_value =~ /median|mean|sum/i) { 
			
			# collect dataset statistics
			my %stathash = _get_statistics_hash($index);
			unless (%stathash) { 
				warn " unable to get statistics! nothing done\n";
				return;
			}
			
			# assign the appropriate value
			if ($request_value eq 'median') {
				$value = $stathash{'median'};
				printf "  median value for %s is $value\n", $Data->name($index);
			}
			elsif ($request_value eq 'mean') {
				$value = $stathash{'mean'};
				printf "  mean value for %s is $value\n", $Data->name($index);
			}
			elsif ($request_value eq 'sum') {
				$value = $stathash{'sum'};
				printf "  sum value for %s is $value\n", $Data->name($index);
			}
		} 
		elsif ($request_value =~ /^\-?\d+(?:\.\d+)?(?:[eE]\-?\d+)?$/) { 
			# a numeric value, allowing for negatives, decimals, exponents
			$value = $request_value;
		} 
		else {
			warn " unrecognized value '$request_value'; nothing done\n";
			return;
		}
		
		
		# generate subtraction product
		my $failed_count = 0; # failed count  
		$index = _prepare_new_destination($index, "_$mathed\_$value") if $placement =~ /^n/i;
		$Data->iterate( sub {
			my $row = shift;
			if ($row->value($index) eq '.') {
				$failed_count++;
			}
			else {
				$row->value($index, &$calculate( $row->value($index), $value) );
			}
		} );	
		# print conclusion
		printf " dataset %s  was $mathed by $value\n", $Data->name($index);
		if ($failed_count) {
			print " $failed_count datapoints could not be $mathed\n";
		}
		$dataset_modification_count++;
	}
	
	# done
	return $dataset_modification_count;
}



sub number_function {
	# This subroutine will number the datapoints or lines
	
	# request dataset
	my $line = " The numbers will be entered as a dataset in front of " . 
		"which current dataset?\n Enter nothing to place the numbers " .
		"at the end.    ";
	my $index = _request_index($line);
	# a return of -1 indicates the numbers will be put at the end
	
	# generate new column
	# we'll put the numbers at the end for now
	my $new_name;
	if ($function and $opt_name) {
		$new_name = $opt_name;
	}
	else {
		$new_name = 'Numbers';
	}
	my $new_position = $Data->add_column($new_name);
	
	# put in numbers for each line
	$Data->iterate( sub {
		my $row = shift;
		$row->value($new_position, $row->row_index);
	} );
	
	
	# Reorder the columns
	if ($index >= 0) {
		my @new_order;
		for (my $i = 0; $i < $new_position; $i++) {
			# we will walk through the list of index numbers sequentially
			# put each index into the new_order array
			# when we encounter the requested index to put the new dataset, 
			# we'll put the new_position, and then the requested index
			if ($i == $index) {
				# we want to put the new Numbers dataset before this index
				push @new_order, $new_position;
				push @new_order, $i;
			}
			else {
				# otherwise, put the current index in the order
				push @new_order, $i;
			}
		}
		$Data->reorder_column(@new_order);
	} 
	
	print " Lines numbered at column $index\n";
	return 1;
}



sub write_summary_function {
	# this will write out a summary file of the data
	
	# determine indices to summarize
	my ($startcolumn, $stopcolumn);
	if ($function) {
		# running under automatic mode
		# check if user supplied indices
		if (scalar @opt_indices >= 2) {
			# assume contiguous indices, use the first and last one
			# it may also be simply the start and stop indices
			$startcolumn = $opt_indices[0];
			$stopcolumn  = $opt_indices[-1];
		}
		# otherwise the summary module will automatically deduce the columns
	}
	else {
		# request indices only when running interactively and not automatically
		($startcolumn, $stopcolumn) = _request_indices(
			" Enter the starting and/or ending indices of the datasets to summarize\n". 
			" Or nothing for automatic detection     "
		);
	}
	
	# write the summary
	my $sumfile = $Data->summary_file(
		'filename'     => $outfile,
		'startcolumn'  => $startcolumn,
		'endcolumn'    => $stopcolumn,
	);
	
	# report outcome
	if ($sumfile) {
		print " wrote summary file '$sumfile'\n";
	}
	else {
		print " unable to write summary file!\n";
	}
	
	# since no changes have been made to the data structure, return 0
	return;
}


sub export_function {
	# this will export the file into an even simpler text format
	
	# generate a possible new name based on the input name
	my $possible_name = $Data->path . $Data->basename . '_out.txt';
	
	# determine the export file name
	my $outfilename;
	if (defined $outfile) {
		# use the outfile name specified on the command line
		# note that this could be overwritten later if $modification > 0
		# but to allow for automated execution, we can't ask the user
		# for verification
		$outfilename = $outfile;
	}
	elsif ($function) {
		# automatic execution, don't ask the user
		$outfilename = $possible_name;
	}
	else {
		# ask for new filename
		print " Provide an exported file name [$possible_name]  ";
		my $answer = <STDIN>;
		chomp $answer;
		if ($answer eq '') {
			# user wants to use the suggested name
			$outfilename = $possible_name;
		} 
		else {
			# user requested own name
			$outfilename = $answer;
		}
	}
		
	
	# write the file
	my $write_results = $Data->write_file(
		'filename'  => $outfilename,
		'gz'        => $gz,
		'simple'    => 1,
	);
	
	# report write results
	if ($write_results) {
		print " Exported data to simple text file '$write_results'\n";
	}
	else {
		print " Unable to export data to file!\n";
	}
	
	# since no changes have been made to the data structure, return
	return;
}



sub export_treeview_function {
	# this is a specialized function to export a datafile into a format 
	# compatible with the Treeview program
	
	print " Exporting CDT file for Treeview and Cluster analysis\n";
	
	# First check for previous modifications
	if ($modification and not $function) {
		print " There are existing unsaved changes to the data. Do you want to\n";
		print " save these first before making required, irreversible changes? ".
			"y/n  ";
		my $answer = <STDIN>;
		chomp $answer;
		if ($answer =~ /^y/i) {
			rewrite_function();
		}
	}
	
	# Set options for placement for subsequent manipulations
	$opt_placement = 'r';
	
	
	### Get user information for processing
	# Identify the dataset columns
	# we need one unique name column, and a range of data columns
	# the rest will be deleted
	my @datasets = _request_indices(
		" Enter one unique name column, followed by a range of data columns  "
	);
	unless (@datasets) {
		warn " Unknown datasets. Nothing done.\n";
		return;
	}
	
	# Identify the manipulations requested
	my @manipulations;
	if ($function) {
		# automatic function, use the command line target option
		@manipulations = split /,/, $opt_target;
	}
	else {
		# ask the user
		print <<LIST;
Available dataset manipulations\n";
  su - decreasing sort by sum of row values
  sm - decreasing sort by mean of row values
  cg - median center features (genes)\n";
  cd - median center datasets\n";
  zd - convert dataset to Z-scores\n";
  pd - convert dataset to percentile rank\n";
  L2 - convert dataset to log2\n";
  n0 - convert null values to 0\n";
Enter the manipulation(s) in order of desired execution  
LIST
		my $answer = <STDIN>;
		chomp $answer;
		@manipulations = split /[,\s]+/, $answer;
	}
	
	
	### First, delete extraneous datasets or columns
	# the CDT format for Treeview expects a unique ID and NAME column
	# we will duplicate the first column
	unshift @datasets, $datasets[0];
	
	# perform a reordering of the columns
	$Data->reorder_column(@datasets);
	
	# rename the first two columns
	$Data->name(0, 'ID');
	$Data->name(1, 'NAME');
	
	# we now have just the columns we want
	# reset the dataset indices to what we currently have
	# name and ID should be index 0 and 1
	@datasets = (2 .. $Data->number_columns - 1);
	
	
	### Second, perform dataset manipulations
	foreach (@manipulations) {
		if (/^su$/i) {
			# decreasing sort by sum of row values
			$opt_target = 'sum';
			combine_function(@datasets);
			my $i = $Data->number_columns - 1;
			$opt_direction = 'd';
			sort_function($i);
			$Data->delete_column($i);
		}
		if (/^sm$/i) {
			# decreasing sort by sum of row values
			$opt_target = 'mean';
			combine_function(@datasets);
			my $i = $Data->number_columns - 1;
			$opt_direction = 'd';
			sort_function($i);
			$Data->delete_column($i);
		}
		elsif (/^cg$/i) {
			# Median center features
			print " median centering features....\n";
			center_function(@datasets);
		}
		elsif (/^cd$/i) {
			# Median center datasets
			print " median centering datasets....\n";
			$opt_target = 'median';
			$opt_zero = 1;
			subtract_function(@datasets);
		}
		elsif (/^zd$/i) {
			# Z-score convert dataset
			print " converting datasets to Z-scores....\n";
			zscore_function(@datasets);
		}
		elsif (/^pd$/i) {
			# convert dataset to percentile rank
			print " converting datasets to percentile ranks....\n";
			percentile_rank_function(@datasets);
		}
		elsif (/^l2$/i) {
			# convert dataset to log2 values
			print " converting datasets to log2 values....\n";
			log2_function(@datasets);
		}
		elsif (/^n0$/i) {
			# convert nulls to 0
			print " converting null values to 0.0....\n";
			$opt_target = '0.0';
			$opt_zero   = 1 unless defined $opt_zero; # convert 0s too
			convert_nulls_function(@datasets);
		}
		else {
			warn " unkown manipulation '$_'!\n";
		}
	}
	
	
	### Third, export a simple file
	if ($outfile) {
		unless ($outfile =~ /\.cdt$/i) {
			# make sure it has .cdt extension
			$outfile .= '.cdt';
		}
	}
	else {
		# generate file name
		$outfile = $Data->path . $Data->basename . '.cdt';
	}
	$gz = 0;
	export_function();
	
	
	### Finally, reset the modification status 
	# We don't to record these changes upon exit as they are specific to 
	# exporting for treeview, any pre-existing changes should have been 
	# saved earlier
	$modification = 0;
	return;
}



sub center_function {
	# this will center the datapoints in each row
	
	# identify the datasets to be centered
	my @datasets;
	if (@_) {
		# provided from an internal subroutine
		@datasets = @_;
	}
	else {
		# otherwise request from user
		@datasets = _request_indices(
			" Enter the range of dataset indices to be centered  "
		);
	}
	unless (@datasets) {
		warn " Unknown datasets. Nothing done.\n";
		return;
	}
	
	# Center the data
	# We will determine the median value for all the datapoints in the 
	# indicated datasets for each feature (row). Then the median value will
	# be subtracted from each value. The median value of the adjusted 
	# datapoints should then be effectively 0.
	$Data->iterate( sub {
		my $row = shift;
		
		# collect the datapoint values in each dataset for the feature  
		my @values;
		foreach my $d (@datasets) {
			push @values, $row->value($d) if $row->value($d) ne '.';
		}
		next unless (@values);
		
		# determine median value
		my $median_value = median(@values);
		
		# adjust the datapoints
		foreach my $d (@datasets) {
			next if $row->value($d) eq '.';
			$row->value($d, ($row->value($d) - $median_value) );
		}
	} );
	
	# annotate metadata
	foreach my $d (@datasets) {
		$Data->metadata($d, 'centered', 'median');
	}
	
	# print conclusion
	printf " Datasets %s  median centered\n", join(", ", @datasets);
	return 1;
}


sub new_column_function {
	# this will generate a new dataset
	
	# request column name
	my $name;
	if (defined $opt_name) {
		# command line option
		$name = $opt_name;
	}
	elsif ($function) {
		$name = 'New_Column';
	}
	else {
		# interactively ask the user
		print " Enter the name for the new column   ";
		$name = <STDIN>;
		chomp $name;
	}
	
	# request value
	my $value;
	if (defined $opt_target) {
		# command line option
		$value = $opt_target;
	}
	elsif ($function) {
		$value = '.'; # null
	}
	else {
		# interactively ask the user
		print " Enter the common value to be assigned in the new column   ";
		$value = <STDIN>;
		chomp $value;
	}
	
	# generate the new dataset
	my $new_position = $Data->add_column($name);
	$Data->iterate( sub {
		shift->value($new_position, $value);
	} );
	
	# done
	print " Added new dataset '$name' at index $new_position with value '$value'\n";
	return 1;
}



sub rewrite_function {
	
	# check output file name
	my $rewrite_filename;
	if (defined $outfile) {
		# use the defined outfile name given
		$rewrite_filename = $outfile;
	} 
	else {
		# ask the user for a new name
		print " Enter a new file name [$infile]  ";
		my $answer = <STDIN>;
		chomp $answer;
		if ($answer) {
			$rewrite_filename = $answer;
		}
		else {
			# accept default of using the input file name
			$rewrite_filename = $infile;
		}
	}
	
	# write the file
	my $write_results = $Data->write_file(
		'filename'  => $rewrite_filename,
		'gz'        => $gz,
	);
	
	# report write results
	if ($write_results) {
		print " $modification manipulations performed\n Wrote datafile $write_results\n";
	}
	else {
		print " unable to re-write to file '$rewrite_filename'!\n";
	}
	
	return;
}


sub view_function {
	
	my $start_row = 0;
	my $stop_row = 10;
	
	# print the contents
	while ($stop_row) {
		
		# print the table contents
		$stop_row = $Data->last_row if $stop_row > $Data->last_row;
		for (my $i = $start_row; $i <= $stop_row; $i++) {
			my $text = join("  ", $Data->row_values($i));
			print "$text\n";
		}
		
		# continue or return the next menu response
		print "\nPress Return for next 10 lines, or enter a (m)enu command   ";
		my $answer = <STDIN>;
		chomp $answer;
		if ($answer) {
			if (exists $letter_to_function{$answer}) {
				$modification += &{ 
					$function_to_subroutine{
						$letter_to_function{$answer}
					}
				};
			}
			else {
				print "unkown response\n";
			}
			return;
		}
		else {
			# print next 10 lines
			$start_row += 10;
			$stop_row += 10;
		}
	}
	return;
}





################################################################################
##############               Internal Subroutines                 ##############
################################################################################


sub _get_letter_to_function_hash {
	# this hash converts the one-letter menu key to the function name
	# the key is the menu letter
	# the value is the function name
	return (
		't' => "stat",
		'R' => "reorder",
		'D' => "delete",
		'n' => "rename",
		'b' => "number",
		'C' => "concatenate",
		'T' => "split",
		'O' => "coordinate",
		'o' => "sort",
		'g' => "gsort",
		'N' => "null",
		'P' => "duplicate",
		'A' => "above",
		'B' => "below",
		'S' => "specific",
		'K' => "keep",
		'U' => "cnull",
		'G' => "absolute",
		'I' => "minimum",
		'X' => "maximum",
		'a' => "add",
		'u' => "subtract",
		'y' => "multiply",
		'v' => "divide",
		's' => "scale",
		'p' => "pr",
		'Z' => "zscore",
		'l' => "log",
		'L' => "delog",
		'f' => "format",
		'c' => "combine",
		'r' => "ratio",
		'd' => "diff",
		'z' => "normdiff",
		'e' => "center",
		'w' => "new",
		'Y' => "summary",
		'x' => "export",
		'W' => "rewrite",
		'i' => "treeview",
		'h' => "help",
		'V' => "view",
		'q' => "write_quit",
		'Q' => "quit",
		'm' => "menu",
	);
}

sub _get_function_to_subroutine_hash {
	# this hash converts the function name to the actual subroutine for the function
	# the key is the function name
	# the value is a scalar reference to the subroutine
	return (
		'stat'        => \&print_statistics_function,
		'reorder'     => \&reorder_function,
		'delete'      => \&delete_function,
		'rename'      => \&rename_function,
		'number'      => \&number_function,
		'concatenate' => \&concatenate_function,
		'split'       => \&split_function,
		'coordinate'  => \&coordinate_function,
		'sort'        => \&sort_function,
		'gsort'       => \&genomic_sort_function,
		'null'        => \&toss_nulls_function,
		'duplicate'   => \&toss_duplicates_function,
		'above'       => \&toss_above_threshold_function,
		'below'       => \&toss_below_threshold_function,
		'specific'    => \&toss_specific_values_function,
		'keep'        => \&keep_specific_values_function,
		'cnull'       => \&convert_nulls_function,
		'absolute'    => \&convert_absolute_function,
		'minimum'     => \&minimum_function,
		'maximum'     => \&maximum_function,
		'add'         => \&add_function,
		'subtract'    => \&subtract_function,
		'multiply'    => \&multiply_function,
		'divide'      => \&divide_function,
		'scale'       => \&median_scale_function,
		'pr'          => \&percentile_rank_function,
		'zscore'      => \&zscore_function,
		'log'         => \&log_function,
		'log2'        => \&log_function, # holdover from previous
		'delog'       => \&delog_function,
		'delog2'      => \&delog_function,
		'format'      => \&format_function,
		'combine'     => \&combine_function,
		'ratio'       => \&ratio_function,
		'diff'        => \&difference_function,
		'normdiff'    => \&normalized_difference_function,
		'center'      => \&center_function,
		'new'         => \&new_column_function,
		'summary'     => \&write_summary_function,
		'export'      => \&export_function,
		'rewrite'     => \&rewrite_function,
		'treeview'    => \&export_treeview_function,
		'view'        => \&view_function,
		'help'        => \&print_online_help,
		'menu'        => \&print_menu,
		'write_quit'  => \&write_and_quit_function,
		'quit'        => \&quit_function,
	);
}


sub _print_headers {
	# this subroutine will print the names of the datasets from the column 
	# headers. It will also generate a hash of the dataset names with their 
	# index numbers as keys
	print " These are the datasets in the file\n";
	my $i = 0;
	foreach ($Data->list_columns) {
		print "\t$i\t$_\n";
		$i++;
	}
}



sub _request_index {
	# this subroutine will determine which dataset index to use
	
	# if index is specified on the command line, that will be used
	# alternatively, it will ask the user which dataset to process.	
	# it will return the index number
	
	my $line = shift; # the custom request line to give the user
	my $index;
	
	if (@opt_indices) {
		# index array is specified on the command line
		# use the first element in the global index array
		$index = $opt_indices[0];
		unless ( _validate_index_list($index) )  {
			# return an error value
			$index = -1;
		}
	}
	
	else {
		# request interactively from the user
		
		_print_headers();
		print $line; # print the user prompt
		$index = <STDIN>;
		chomp $index;
		unless ( _validate_index_list($index) )  {
			# return an error value
			$index = -1;
		}
	}
	
	return $index;
}


sub _request_indices {
	# this subroutine will determine which datasets are to be used
	# if the indices are specified on the command line, those will be used
	# alternatively, it will ask the user for the indices interactively
	
	my $line = shift; # the custom request line to give the user
	my @indices;
	
	# get list of indices
	if (@opt_indices) {
		# index array is specified on the command line
		@indices = @opt_indices;
	}
	else {
		# request interactively from the user
		_print_headers();
		print $line; # print the user prompt
		my $response = <STDIN>;
		chomp $response;
		@indices = parse_list($response); 
	}
	
	# check the list of indices
	unless ( _validate_index_list(@indices) ) {
		return;
	}
	
	return @indices;
}


sub _request_placement {
	# this subroutine will determine where to put a new dataset
	# replace current or generate new
	
	my $placement;
	if (defined $opt_placement) {
		# use the command line specified placement
		$placement = $opt_placement;
	}
	else {
		# request placement from user
		print " (r)eplace dataset or generate (n)ew dataset?  ";
		$placement = <STDIN>;
		chomp $placement;
	}
	
	return $placement;
}



sub _prepare_new_destination {
	# User requested new column placement
	my ($index, $suffix) = @_;
	
	my $new_index = $Data->copy_column($index);
	my $new_name;
	if ($function and $opt_name) {
		# automatic execution and new name was specifically given 
		$new_name = $opt_name;
	}
	else {
		$new_name = $Data->name($index) . $suffix;
	}
	$Data->name($new_index, $new_name);
		
	return $new_index;
}



sub _validate_index_list {
	# this subroutine will ensure that the list of index numbers represent 
	# actual columns, i.e., is each number < than the number of datasets in 
	# the data table
	
	foreach my $num (@_) {
		# check that each number represents a metadata hash 
		# and presumably dataset
		unless ($num =~ /^\d+$/ and $num <= $Data->number_columns) {
			warn " requested index number '$num' is not valid!\n";
			return;
		}
	}
	return 1;
}



sub _get_statistics_hash {
	# internal subroutine to get statistics
	
	my ($index, $zero) = @_;
	# the index number (column) of the dataset
	# the exception rule to work with 0 numbers
	unless (defined $zero) {
		# use global command-line specified value if present
		$zero = $opt_zero;
	}
	
	# collect the values in the dataset
	my @invalues = $Data->column_values($index);
	shift @invalues; # remove header
	my @goodvalues;
	foreach my $v (@invalues) {
		next if $v eq '.';
		if ($v == 0) {
			# we need to determine whether we can accept 0 values
			unless (defined $zero) {
				# wasn't defined on the command line, so stop the program and ask the user
				print " Include 0 values in the statistics? y or n  ";
				my $answer = <STDIN>;
				if ($answer =~ /^y/i) {
					$zero = 1;
				}
				else {
					$zero = 0;
				}
				
				# remember for next time, chances are user may still want this
				# value again in the future
				$opt_zero = $zero;
			}
			next unless $zero;
		} 
		push @goodvalues, $v;
	}
	
	# check that we have values
	unless (@goodvalues) {
		warn " no valid values collected from dataset index '$index'!";
		return;
	}
	
	return statshash(@goodvalues);
}




__END__

=head1 NAME

manipulate_datasets.pl

A progam to manipulate tab-delimited data files.

=head1 SYNOPSIS

manipulate_datasets.pl [--options ...] <input_filename> 

  Options:
  --in <input_filename>
  --func [ reorder | delete | rename | new | number | concatenate | split | 
           sort | gsort | null | duplicate | above | below | specific | keep
           coordinate | cnull | absolute | minimum | maximum | log | delog | 
           format | pr | add | subtract | multiply | divide | combine | scale | 
           zscore | ratio | diff | normdiff | center | rewrite | export | 
           treeview | summary | stat ]
  --index <integers>
  --exp | --num <integer>
  --con | --den <integer>
  --target <text> or <number>
  --place [r | n]
  --(no)zero
  --dir [i | d]
  --name <text>
  --out <filename>
  --log
  --gz
  --version
  --help
  --doc

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <input_filename>

Provide the name of an input file. The file must be a tab-delimited text file,
with each row specifiying a genomic feature (gene) or region, and each column
representing identifying information or a dataset value. The first row should
include the name for each column, e.g. the name of the database dataset from
which the values were collected.

If the file was generated using a BioToolBox script, then each column may have
metadata stored in a header comment line at the beginning of the file.
Manipulations on the data in a column dataset will be added to this metadata and
recorded upon file write.

The input file may be compressed using the gzip program and recognized with 
the extension '.gz'.  

=item --func <function>

The program is designed to be run interactively. However, single manipulations 
may be performed on single datasets by specifying a function name and any 
other required options. These functions include the following.
  
B<reorder> B<delete> B<rename> B<new> B<number> B<concatenate> B<split> B<coordinate>
B<sort> B<gsort> B<null> B<duplicate> B<above> B<below> B<specific> B<keep>
B<cnull> B<absolute> B<minimum> B<maximum> B<log> B<delog> B<format> B<pr>
B<add> B<subtract> B<multiply> B<divide> B<combine> B<scale> B<zscore> B<ratio> B<diff> B<normdiff> B<center>
B<rewrite> B<export> B<treeview> B<summary> B<stat>
  
Refer to the FUNCTIONS section for details.

=item --index <integers>

Provide the 0-based index number of the column(s) on which to perform the 
function(s). Multiple indices may also be specified using a comma delimited 
list without spaces. Ranges of contiguous indices may be specified using a 
dash between the start and stop. For example, '1,2,5-7,9' would indicate 
datasets '1, 2, 5, 6, 7, and 9'. 

=item --exp <integer>

=item --num <integer>

Specify the 0-based index number to be used for the experiment or numerator 
column with the 'ratio' or 'difference' functions. Both flags are aliases
for the same thing.

=item --con <integer>

=item --den <integer>

Specify the 0-based index number to be used for the control or denominator 
column with the 'ratio' or 'difference' functions. Both flags are aliases
for the same thing.

=item --target <string> or <number>

Specify the target value when using various functions. This is a catch-all 
option for a number of functions. Please refer to the function description 
for more information.

=item --place [r | n]

Specify the placement of a transformed column. Two values are accepted ('r' 
or 'n'):
  
  - (r)eplace the original column with new values
  - add as a (n)ew column

Defaults to new placement when executed automatically using the --func 
option, or prompts the user when executed interactively.

=item --(no)zero

Specify that zero values should or should not be included when 
calculating statistics or converting null values on a column. The default is 
undefined, meaning the program may ask the user interactively whether to 
include zero values or not.

=item --dir [i | d]

Specify the direction of a sort: 
  
  - (i)ncreasing
  - (d)ecreasing
  
=item --name <string>

Specify a new column name when re-naming a column using the rename function 
from the command line. Also, when generating a new column using a defined 
function (--func <function>) from the command line, the new column will use 
this name.

=item --out <filename>

Optionally provide an alternative output file name. If no name is provided, 
then the input file will be overwritten with a new file. Appropriate 
extensions will be appended as necessary.

=item --log 

Indicate whether the data is (not) in log2 space. This is required to ensure 
accurate mathematical calculations in some manipulations. This is not necessary 
when the log status is appropriately recorded in the dataset metadata.

=item --gz 

Indicate whether the output file should (not) be compressed. The appropriate extension will be 
added. If this option is not specified, then the compression status of the input file will be 
preserved.

=item --version

Print the version number.

=item --help

Display the POD documentation using perldoc. 

=back

=head1 DESCRIPTION

This program allows some common mathematical and other manipulations on one
or more columns in a datafile. The program is designed as a simple
replacement for common manipulations performed in a full featured
spreadsheet program, e.g. Excel, particularly with datasets too large to be
loaded, all in a conveniant command line program. The program is designed
to be operated primarily as an interactive program, allowing for multiple
manipulations to be performed. Alternatively, single manipulations may be
performed as specified using command line options. As such, the program can
be called in shell scripts.

Note that the datafile is loaded entirely in memory. For extremely large 
datafiles, e.g. binned genomic data, it may be best to first split the 
file into chunks (use C<split_data_file.pl>), perform the manipulations, 
and recombine the file (use C<join_data_file.pl>). This could be done 
through a simple shell script.

The program keeps track of the number of manipulations performed, and if 
any are performed, will write out to file the changed data. Unless an 
output file name is provided, it will overwrite the input file (NO backup is
made!).

=head1 FUNCTIONS

This is a list of the functions available for manipulating columns. These may 
be selected interactively from the main menu (note the case sensitivity!), 
or specified on the command line using the --func option.

=over 4

=item B<stat> (menu option B<t>)

Print some basic statistics for a column, including mean, 
median, standard deviation, min, and max. If 0 values are present,
indicate whether to include them (y or n)

=item B<reorder> (menu option B<R>)

The column may be rewritten in a different order. The new order 
is requested as a string of index numbers in the desired order. 
Also, a column may be deleted by skipping its number or duplicated
by including it twice.

=item B<delete> (menu option B<D>)

One or more column may be selected for deletion.

=item B<rename> (menu option B<n>)

Assign a new name to a column. For automatic execution, use the --name 
option to specify the new name. Also, for any automatically executed 
function (using the --func option) that generates a new column, the 
column's new name may be explicitly defined with this option.

=item B<number> (menu option B<b>)

Assign a number to each datapoint (or line), incrementing from 1 
to the end. The numbered column will be inserted after the requested 
column index.

=item B<concatenate> (menu option B<C>)

Concatenate the values from two or more columns into a single new 
column. The character used to join the values may be specified 
interactively or by the command line option --target (default is '_' 
in automatic execution mode). The new column is appended at the end.

=item B<split> (menu option B<T>)

Split a column into two or more new columns using a specified character 
as the delimiter. The character may be specified interactively or 
with the --target command line option (default is '_' in automatic 
execution mode). The new columns are appended at the end. If the 
number of split items are not equal amongst the rows, absent values 
are appended with null values.

=item B<coordinate> (menu option B<O>)

Generate a coordinate string from the chromosome, start, and, if 
present, stop coordinate values as a new column. The string will 
have the format "chr:start-stop" or "chr:start". This is useful 
in making unique identifiers or working with genome browsers.

=item B<sort> (menu option B<o>)

The entire data table is sorted by a specific column. The first
datapoint is checked for the presence of letters, and the data 
table is then sorted either asciibetically or numerically. If the 
sort method cannot be automatically determined, it will ask. The 
direction of sort, (i)ncreasing or (d)ecreasing, is requested. 

=item B<gsort> (menu option B<g>)

The entire data table is sorted by increasing genomic position, 
first by chromosome then by start position. These columns must exist 
and have recognizable names (e.g. 'chromo', 'chromosome', 'start').

=item B<null> (menu option B<N>)

Delete rows that contain a null value in one or more 
columns. Some of the other functions may not work properly if
a non-value is present. If 0 values are present, indicate whether
to toss them (y or n). This may also be specified as a command line 
option using the --except flag.

=item B<duplicate> (menu option B<P>)

Delete rows with duplicate values. One or more columns may be 
selected to search for duplicate values. Values are simply concatenated 
when multiple columns are selected. Rows with duplicated values are 
deleted, always leaving the first row.

=item B<above> (menu option B<A>)

Delete rows with values that are above a certain threshold value. 
One or more columns may be selected to test values for the 
threshold. The threshold value may be requested interactively or 
specified with the --target option.

=item B<below> (menu option B<B>)

Delete rows with values that are below a certain threshold value. 
One or more columns may be selected to test values for the 
threshold. The threshold value may be requested interactively or 
specified with the --target option.

=item B<specific> (menu option B<S>)

Delete rows with values that contain a specific value, either text 
or number. One or more columns may be selected to check for values. 
The specific values may be selected interactively from a list or 
specified with the --target option.

=item B<keep> (menu option B<K>)

Keep only those rows with values that contain a specific value, 
either text or number. One or more columns may be selected to check 
for values. The specific values may be selected interactively from a 
list or specified with the --target option.

=item B<cnull> (menu option B<U>)

Convert null values to a specific value. One or more columns may 
be selected to convert null values. The new value may be requested 
interactively or defined with the --target option.  

=item B<absolute> (menu option B<G>)

Convert signed values to their absolute value equivalents. One or 
more columns may be selected to convert.

=item B<minimum> (menu option B<I>)

Reset datapoints whose values are less than a specified minimum 
value to the minimum value. One or more columns may be selected 
to reset values to the minimum. The minimum value may be requested 
interactively or specified with the --target option. 

=item B<maximum> (menu option B<X>)

Reset datapoints whose values are greater than a specified maximum 
value to the maximum value. One or more columns may be selected 
to reset values to the maximum. The maximum value may be requested 
interactively or specified with the --target option. 

=item B<add> (menu option B<a>)

Add a value to a column. A real number may be supplied, or the words
'mean', 'median', or 'sum' may be entered as a proxy for those statistical
values of the column. The column may either be replaced or added
as a new one. For automatic execution, specify the number using the
--target option.

=item B<subtract> (menu option B<u>)

Subtract a value from a column. A real number may be supplied, or the words
'mean', 'median', or 'sum' may be entered as a proxy for those statistical
values of the column. The column may either be replaced or added
as a new one. For automatic execution, specify the number using the
--target option.

=item B<multiply> (menu option B<y>)

Multiply a column by a value. A real number may be supplied, or the words
'mean', 'median', or 'sum' may be entered as a proxy for those statistical
values of the column. The column may either be replaced or added
as a new one. For automatic execution, specify the number using the
--target option.

=item B<divide> (menu option B<v>)

Divide a column by a value. A real number may be supplied, or the words
'mean', 'median', or 'sum' may be entered as a proxy for those statistical
values of the column. The column may either be replaced or added
as a new one. For automatic execution, specify the number using the
--target option.

=item B<scale> (menu option B<s>)

A column may be a median scaled as a means of normalization 
with other columns. The current median of the column requested is
presented, and a new median target is requested. The column may 
either be replaced with the median scaled values or added as a new 
column. For automatic execution, specify the new median target 
with the --target option.

=item B<pr> (menu option B<p>)

A column may be converted to a percentile rank, whereby the
data values are sorted in ascending order and assigned a new value 
from 0..1 based on its rank in the sorted order. The column may 
either be replaced with the percentile rank or added as a new
column. The original order of the column is maintained.

=item B<zscore> (menu option B<Z>)

Generate a Z-score or standard score for each value in a column. The
Z-score is the number of standard deviations the value is away from
the column's mean, such that the new mean is 0 and the standard 
deviation is 1. Provides a simple method of normalizing columns
with disparate dynamic ranges.

=item B<log> (menu option B<l>)

A column may be converted to log values. The column may either 
be replaced with the log values or added as a new column. Use 
the --target option to specify the base (usually 2 or 10).

=item B<delog> (menu option B<L>)

A column that is currently in log space may be converted back to
normal numbers. The column may either be replaced with the 
new values or added as a new column. Use the --target option to 
specify the base (usually 2 or 10). The base may be obtained from the 
metadata.

=item B<format> (menu option B<f>)

Format the numbers of a column to a given number of decimal places. 
An integer must be provided. The column may either be replaced or 
added as a new column. For automatic execution, use the --target 
option to specify the number decimal places.

=item B<combine> (menu option B<c>)

Mathematically combine the data values in two or more columns. The 
methods for combining the values include mean, median, min, max, 
stdev, or sum. The method may be specified on the command line 
using the --target option. The combined data values are added as a 
new column.

=item B<ratio> (menu option B<r>)

A ratio may be generated between two columns. The experiment and 
control columns are requested and the ratio is added as a new
column. For log2 numbers, the control is subtracted from the
experiment. The log2 status is checked in the metadata for the 
specified columns, or may be specified as a command line option, or
asked of the user.

=item B<diff> (menu option B<d>)

A simple difference is generated between two existing columns. The 
values in the 'control' column are simply subtracted from the 
values in the 'experimental' column and recorded as a new column.
For enumerated columns (e.g. tag counts from Next Generation 
Sequencing), the columns should be subsampled to equalize the sums 
of the two columns. The indices for the experimental and control columns 
may either requested from the user or supplied by the --exp and 
--con command line options. 

=item B<normdiff> (menu option B<z>)

A normalized difference is generated between two existing columns. 
The difference between 'control' and 'experimental' column values 
is divided by the square root of the sum (an approximation of the 
standard deviation). This is supposed to yield fewer false positives
than a simple difference (see Nix et al, BMC Bioinformatics, 2008).
For enumerated datasets (e.g. tag counts from Next Generation 
Sequencing), the datasets should be subsampled to equalize the sums 
of the two datasets. The indices for the experimental and control columns 
may either requested from the user or supplied by the --exp and 
--con command line options. 

=item B<center> (menu option B<e>)

Center normalize the datapoints in a row by subtracting the mean or
median of the datapoints. The range of columns is requested or 
provided by the --index option. Old values are replaced by new 
values. This is useful for visualizing data as a heat map, for example.

=item B<new> (menu option B<w>)

Generate a new column which contains an identical value for 
each datapoint (row). The value may be either requested interactively or 
supplied using the --target option. This function may be useful for 
assigning a common value to all of the data points before joining the 
data file with another.

=item B<summary> (menu option B<y>)

Write out a summary of collected windowed data file, in which the mean 
for each of the data columns is calculated, transposed (columns become 
rows), and written to a new data file. This is essentially identical to 
the summary function from the biotoolbox analysis scripts 
L<map_relative_data.pl> and L<pull_features.pl>. It assumes that each 
column has start and stop metadata. The program will automatically 
identify available columns to summarize based on their name. In 
interactive mode, it will request the contiguous range of start and 
ending columns to summarize. The contiguous columns may also be 
indicated using the --index option. By default, a new file using the 
input file base name appended with '_summary' is written, or a 
filename may be specified using the --out option.

=item B<export> (menu option B<x>)

Export the data into a simple tab-delimited text file that contains no 
metadata or header information. Non-values '.' are converted to  
true nulls. If an output file name is specified using the --outfile 
option, it will be used. Otherwise, a possible filename will be 
suggested based on the input file name. If any modifications are 
made to the data structure, a normal data file will still be written. 
Note that this could overwrite the exported file if the output file name
was specified on the command line, as both file write subroutines will 
use the same name!

=item B<treeview> (menu option B<i>)

Export the data to the CDT format compatible with both Treeview and 
Cluster programs for visualizing and/or generating clusters. Specify the 
columns containing a unique name and the columns to be analyzed (e.g. 
--index <name>,<start-stop>). Extraneous columns are removed. 
Additional manipulations on the columns may be performed prior to 
exporting. These may be chosen interactively or using the codes 
listed below and specified using the --target option.
  
  su - decreasing sort by sum of row values
  sm - decreasing sort by mean of row values
  cg - median center features (rows)
  cd - median center datasets (columns)
  zd - convert columns to Z-scores
  pd - convert columns to percentile ranks
  L2 - convert values to log2
  n0 - convert nulls to 0.0

A simple Cluster data text file is written (default file name 
"<basename>.cdt"), but without the GWEIGHT column or EWEIGHT row. The 
original file will not be rewritten.

=item B<rewrite> (menu option B<W>)

Force the data file contents to be re-written. Useful if you want to 
write an intermediate file during numerous interactive manipulations. 
Consider this as a 'Save as...'.

=back

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
