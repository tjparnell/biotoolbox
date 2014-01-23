#!/usr/bin/env perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long;
use Statistics::Lite qw(:all);
use Bio::ToolBox::data_helper qw(
	find_column_index
	parse_list
);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	write_tim_data_file
	write_summary_data
);
my $VERSION = '1.14';

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
my $main_data_ref = load_tim_data_file($infile) 
	or die "no data loaded from file '$infile'!";
my $data_table_ref = $main_data_ref->{'data_table'}; 
	# an easier reference to the data table
$infile = $main_data_ref->{'filename'};
	# update input file name in case a missing extension was added
print "    Loaded '$infile' with ", $main_data_ref->{'last_row'}, 
	" data rows and ", $main_data_ref->{'number_columns'}, " columns\n\n";

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
} 
else {
	# otherwise interactive execution
	interactive_execution();
}



### Write output
if ($modification > 0) { 
	# a value greater than 0 indicates that changes to the data array have been
	# made and that we need to write an output file
	
	# check output file name
	unless (defined $outfile) {
		# we will overwrite the input file
		$outfile = $infile;
	}
	
	# write the file
	my $write_results = write_tim_data_file(
		'data'      => $main_data_ref,
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




################################################################################
##############             Main Program Subroutines               ##############
################################################################################



sub print_menu {
	print " These are the functions available:\n" .
		"  t  Print S(t)atistics on a dataset\n" .
		"  R  (R)eorder columns in a different order\n" .
		"  D  (D)elete a dataset\n" .
		"  n  Re(n)ame a dataset\n" .
		"  b  Num(b)er the datapoints\n" .
		"  C  Generate (C)oordinate string\n" .
		"  o  S(o)rt all datasets by a specific dataset\n" .
		"  g  (g)enomic sort by chromosome, start\n" .
		"  N  Toss data lines with (N)ull values\n" .
		"  P  Toss data lines with du(P)licate values\n" .
		"  A  Toss data lines with values (A)bove threshold\n" .
		"  B  Toss data lines with values (B)elow threshold\n" .
		"  U  Convert n(U)ll values to a specific value\n" .
		"  L  Convert signed values to an abso(L)ute value\n" .
		"  I  Set a m(I)nimum value\n" .
		"  X  Set a ma(X)imum value\n" .
		"  a  (a)dd a specific value to a dataset\n" .
		"  u  S(u)btract a specific value from a dataset\n" .
		"  y  Multipl(y) a specific value with a dataset\n" .
		"  v  Di(v)ide a dataset by a specific value\n" . 
		"  s  Median (s)cale a dataset\n" .
		"  p  (p)ercentile rank convert a dataset\n" .
		"  Z  Generate (Z)-score values of a dataset\n" .
		"  l  (l)og2 convert the dataset\n" .
		"  2  De-log(2) convert the dataset\n" .
		"  f  (f)ormat decimal numbers in a dataset\n" .
		"  c  (c)ombine datasets together\n" .
		"  r  Generate a (r)atio between two datasets\n" .
		"  d  Generate a (d)ifference between two datasets\n" .
		"  z  Generate a normali(z)ed difference between two datasets\n" .
		"  su (su)bsample a dataset\n" .
		"  si Convert data to (si)gned data according to strand\n" .
		"  st Merge (st)randed datasets into one\n" .
		"  e  C(e)nter normalize feature datapoints \n" .
		"  w  Generate a ne(w) column with a single identical value\n" .
		"  Y  Write out a summar(Y) file of the data\n" .
		"  x  E(x)port into a simple tab-delimited text file\n" .
		"  W  Re(W)rite the file\n" .
		"  T  Export for (T)reeview or Cluster analysis\n" .
		"  h  (h)elp\n" .
		"  q  (q)uit, saving changes if necessary\n" .
		"  Q  (Q)uit without saving changes\n"
		#  m  print this (m)enu
	;
	# unused letters: E F G H jJ kK L O S V Y 
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
	my $response = <STDIN>;
	chomp $response;
	while (1) {
		if (exists $letter_to_function{$response} ) {
			# first check that the letter corresponds to a function
			
			# if the response is quit
			if ($response eq 'Q') {
				# quit without saving changes
				$modification = 0; # pretend we never made changes
				return;
			}
			elsif ($response eq 'q') {
				# quit, saving changes if necessary
				return;
			}
			
			# perform the function
			$modification += &{ 
				$function_to_subroutine{
					$letter_to_function{$response}
				}
			};
			
			# prepare for the next function
			print " Which function is next? [m]enu, [h]elp, just [Q]uit, save & [q]uit  ";
			$response = <STDIN>;
			chomp $response;
		}
		
		else {
			# unrecognized command
			print " unrecognized command. [m]enu, [h]elp, just [Q]uit, save & [q]uit  ";
			$response = <STDIN>;
			chomp $response;
		}
	}
}



################################################################################
##############               Function Subroutines                 ##############
################################################################################



sub print_statistics_function {
	# print simple statistics for the dataset
	
	# request dataset(s)
	my @indices = _request_indices(
		" Enter one or more dataset index numbers to calculate statistics  "
	);
	unless (@indices) {
		warn " unknown index number(s).\n";
		return;
	}
	
	# determine what to do with zero values
	my $zero = $opt_zero; 
	
	# get statistics and print
	foreach my $index (@indices) {
		my %statdata = _get_statistics_hash($index, $zero);
		unless (%statdata) { 
			warn " unable to get statistics for dataset index $index!\n"; 
			return;
		}

		# print the metadata and the calculated statistics
		printf "  Statistics for dataset '%s'\n", $main_data_ref->{$index}{'name'};
		print "   index => $index\n";
		foreach my $key (keys %{ $main_data_ref->{$index} } ) {
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
				printf "   $key => '%s'\n", $main_data_ref->{$index}{$key};
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
	for my $i (0..$main_data_ref->{'last_row'}) {
		# we will work one row at a time
		# including the header row (index 0)
		# assign the old row of data values to a temporary array @old
		my @old = @{ $data_table_ref->[$i] };
		
		my @new;
		foreach (@order) {
			# walk through the new order, assign the old value to the new array
			push @new, $old[$_];
		}
		
		# splice the new row of data values back into the data table, replacing the old
		splice(@{ $data_table_ref }, $i, 1, \@new);
	}
	
	# re-order the metadata hashes
	my %old_metadata;
	for (my $i = 0; $i < $main_data_ref->{'number_columns'}; $i++) {
		# copy the metadata info hash into a temporary hash
		$old_metadata{$i} = $main_data_ref->{$i};
		delete $main_data_ref->{$i}; # delete original
	}
	for (my $i = 0; $i < scalar(@order); $i++) {
		# now copy back from the old_metadata into the main data hash
		# using the new index number in the @order array
		$main_data_ref->{$i} = { %{ $old_metadata{ $order[$i] } } };
			# the above line is ugly syntax, but basically we're creating 
			# a new anonymous hash with the same data contents as the
			# old one. We're referencing the old one by the actual number in
			# the @order array (referenced itself by position)
			# this is to accomodate dataset duplication
		# assign new index number
		$main_data_ref->{$i}{'index'} = $i;
	}
	
	# assign a new number_columns value
	$main_data_ref->{'number_columns'} = scalar @order;
	
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
		" Enter one or more dataset index numbers to be deleted.\n  "
		);
	}
	unless (@deletion_list) {
		warn " No list for deletion! Nothing done!\n";
		return;
	}
	
	# rather than writing a whole new deletion subroutine, we will simply turn the 
	# request for deletion(s) into a re-ordering request, and pass it on to the 
	# re-ordering subroutine. The result will be the same.
	
	# Determine the new order
	@deletion_list = sort {$a <=> $b} @deletion_list;
	my @deleted_names; # an array of the deleted dataset names
	my @order; # the final order of the retained datasets
	for (my $i = 0; $i < $main_data_ref->{'number_columns'}; $i++) {
		# compare each current index with the first one in the list of 
		# deleted indices. if it matches, delete. if not, keep
		
		if ( $i == $deletion_list[0] ) {
			# this particular index should be deleted
			# record the name
			push @deleted_names, $main_data_ref->{$i}{'name'};
			# remove from the deletion_list
			shift @deletion_list;
		}
		
		else {
			# this particular index should be kept
			push @order, $i;
		}
	}
	
	# Perform reordering (deletion)
	my $deletion_success = reorder_function(@order);
	
	# completion
	if ($deletion_success) {
		print " datasets '" . join(", ", @deleted_names) . "' deleted\n";
		return 1;
	}
	else {
		print " nothing deleted!\n";
		return;
	}
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
			" Enter the index number of the dataset to rename  "
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
	my $oldname = $main_data_ref->{$index}{'name'};
	$main_data_ref->{$index}{'name'} = $newname; # assign metadata
	$data_table_ref->[0][$index] = $newname; # assign the column header
	print " $oldname re-named to $newname\n";
	return 1;
}


sub coordinate_function {
	# this subroutine will generate a coordinate string from coordinate values
	
	# identify the coordinates
	my $chr_i   = find_column_index($main_data_ref, '^chr|seq|ref|ref.?seq');
	my $start_i = find_column_index($main_data_ref, '^start|position');
	my $stop_i  = find_column_index($main_data_ref, '^stop|end');
	unless (defined $chr_i and defined $start_i) {
		# cannot add coordinate column, do without ?
		warn " cannot generate coordinates, no chromosome or start column found\n";
		return;
	}
	
	# the new index position is equivalent to the number of columns
	my $new_position = $main_data_ref->{'number_columns'};
	
	# generate coordinates
	if (defined $stop_i) {
		# merge chromosome:start-stop
		for my $row (1 .. $main_data_ref->{'last_row'}) {
			$data_table_ref->[$row][$new_position] = join("", 
				$data_table_ref->[$row][$chr_i], ':', 
				$data_table_ref->[$row][$start_i], '-',
				$data_table_ref->[$row][$stop_i]
			);
		}
	}
	else {
		# merge chromosome:start
		for my $row (1 .. $main_data_ref->{'last_row'}) {
			$data_table_ref->[$row][$new_position] = join("", 
				$data_table_ref->[$row][$chr_i], ':', 
				$data_table_ref->[$row][$start_i]
			);
		}
	}
	
	# generate new metadata
	$main_data_ref->{$new_position} = {
		'name'   => 'Coordinate',
		'index'  => $new_position,
	};
	$main_data_ref->{'number_columns'}++;
	$data_table_ref->[0][$new_position] = 'Coordinate';
		
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
			" Enter one or more dataset index numbers to median scale  "
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
		
		# Check the log2 metadata status
		if (exists $main_data_ref->{$index}{'log2'} ) {
			if ( $main_data_ref->{$index}{'log2'} == 1 ) {
				warn " dataset $main_data_ref->{$index}{'name'} metadata " .
					"reports it is currently in log2 scale.\n" . 
					" You should de-log prior median scaling. Continue with " . 
					"median scaling? y/n\n";
				my $response = <STDIN>;
				next INDEX_LOOP if $response =~ /n/i;
			}
		}
	
		# Retrieve values and calculate median
		my %statdata = _get_statistics_hash($index, 'n');
		unless (%statdata) { 
			warn " unable to get statistics for dataset " . 
				$main_data_ref->{$index}{'name'} . ", index $index!\n"; 
			next INDEX_LOOP;
		}
		print " The median value for dataset " . 
			$main_data_ref->{$index}{'name'} . " is $statdata{'median'}\n";
	
		# Calculate correction value
		my $correction_value = $target / $statdata{median};
	
		# Calculate all of the new values and place accordingly
		if ($placement eq 'r' or $placement eq 'R') {
			# Replace the current value
		
			for my $i (1..$main_data_ref->{'last_row'}) {
				if ( $data_table_ref->[$i][$index] eq '.') {
					# null value, nothing to do
					next;
				}
				else {
					$data_table_ref->[$i][$index] = 
						$correction_value * $data_table_ref->[$i][$index];
				}
			}
		
			# annotate metadata
			$main_data_ref->{$index}{'median_scaled'} = $target;
		
			# results
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		} 
	
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Place as a new dataset
		
			# the new index position is equivalent to the number of columns
			my $new_position = $main_data_ref->{'number_columns'};
		
			# calculate new values
			for my $i (1..$main_data_ref->{'last_row'}) {
				if ( $data_table_ref->[$i][$index] eq '.') {
					# null value, nothing to do
					$data_table_ref->[$i][$new_position] = '.';
				}
				else {
					$data_table_ref->[$i][$new_position] = 
						$correction_value * $data_table_ref->[$i][$index];
				}
			}
		
			# copy the medadata hash and annotate
			my $new_name;
			if ($function and $opt_name) {
				# automatic execution and new name was specifically given 
				$new_name = $opt_name;
			}
			else {
				$new_name = $main_data_ref->{$index}{'name'} . '_scaled';
			}
			_generate_new_metadata(
				$index,
				$new_position,
				'median_scaled',
				$target,
				$new_name,
			);
		
			# results
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		} 
	
		else {
			warn " median scaling NOT done; unknown placement request\n";
			return;
		}
	}
	# report results
	if (@datasets_modified) {
		my $string = $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= " median scaled to $target";
		$string .= $placement =~ /n/i ? " and written as new datasets\n" : "\n";
		print $string;
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
			" Enter one or more dataset index numbers to convert to percentile rank  "
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
		
		# Retrieve values from the specified dataset
		my %values; # these will be put in a hash
		for my $i (1..$main_data_ref->{'last_row'}) {
			my $value = $data_table_ref->[$i][$index];
			if ($value eq '.' or $value eq '') {
				warn 
					" There are null values in dataset index $index! " .
					"Please toss them before proceeding.\n skipping\n";
				next INDEX_LOOP;
			} else {
				$values{$i} = $value;
			}
		}
		
		# Calculate percent rank
		my %percentrank;
		my $n = 1;
		my $total = scalar keys %values;
		foreach (sort { $values{$a} <=> $values{$b} } keys %values) {
			# sort by increasing hash values, not hash keys
			# percentrank is key value (index) divided by total
			$percentrank{$_} = $n/$total;
			$n++;
		}
		
		# Put new values back in
		if ($placement eq 'r' or $placement eq 'R') {
			
			# Replace the contents of the original dataset
			for my $i (1..$main_data_ref->{'last_row'}) {
				$data_table_ref->[$i][$index] = $percentrank{$i};
			}
			
			# update metadata
			$main_data_ref->{$index}{'converted'} = 'percent_rank';
			
			# done
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		} 
		
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Generate a new dataset
			
			# the new index position is equivalent to the number of columns
			my $new_position = $main_data_ref->{'number_columns'};
			
			# copy values
			for my $i (1..$main_data_ref->{'last_row'}) {
				$data_table_ref->[$i][$new_position] = $percentrank{$i};
			}
			
			# copy the medadata hash and annotate
			my $new_name;
			if ($function and $opt_name) {
				# automatic execution and new name was specifically given 
				$new_name = $opt_name;
			}
			else {
				$new_name = $main_data_ref->{$index}{'name'} . '_pr';
			}
			_generate_new_metadata(
				$index,
				$new_position,
				'converted',
				'percent_rank',
				$new_name,
			);
			
			# done
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		} 
		
		else {
			warn " percent rank conversion NOT done; unknown placement request\n";
			return;
		}
	}	
	
	# report results
	if (@datasets_modified) {
		my $string = $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= " converted to percent rank";
		$string .= $placement =~ /n/i ? " and written as new datasets\n" : "\n";
		print $string;
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
			" Enter one or more dataset index numbers to convert to z-scores  "
		);
	}
	unless (@indices) {
		warn " Unknown datasets. Nothing done.\n";
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
		
		# Put new values back in
		if ($placement eq 'r' or $placement eq 'R') {
			
			# Replace the current values
			for my $row (1 .. $main_data_ref->{'last_row'}) {
				next if $data_table_ref->[$row][$index] eq '.';
				$data_table_ref->[$row][$index] = 
					( $data_table_ref->[$row][$index] - $statdata{'mean'} ) / 
					$statdata{'stddevp'};
			}
			
			# update metadata
			$main_data_ref->{$index}{'converted'} = 'Z-score';
			
			# done
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		}
		
		# Generate a new dataset
		elsif ($placement eq 'n' or $placement eq 'N') {
			
			# generate the z scores into a new dataset
			my $new_position = $main_data_ref->{'number_columns'};
			for my $row (1 .. $main_data_ref->{'last_row'}) {
				if ($data_table_ref->[$row][$index] eq '.') {
					$data_table_ref->[$row][$new_position] = '.';
					next;
				}
				$data_table_ref->[$row][$new_position] = 
					( $data_table_ref->[$row][$index] - $statdata{'mean'} ) / 
					$statdata{'stddevp'};
			}
		
			# copy the medadata hash and annotate
			my $new_name = $main_data_ref->{$index}{'name'} . '_Zscore';
			_generate_new_metadata(
				$index,
				$new_position,
				'converted',
				'Z-score',
				$new_name,
			);
			
			# done
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		}
		
		else {
			warn " Z-score conversion NOT done; unknown placement request\n";
			return;
		}
	}
	
	# report results
	if (@datasets_modified) {
		my $string = $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= " converted to Z-scores";
		$string .= $placement =~ /n/i ? " and written as new datasets\n" : "\n";
		print $string;
	}
	return scalar(@datasets_modified);
}




sub sort_function {
	# This will sort the entire data table by the values in one dataset
	
	# Request dataset
	my $index = _request_index(
		" Enter the index number of a numeric dataset to sort by  ");
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
	
	# Sample the dataset values
	# this will be used to guess the sort method, below
	my $example; # an example of the dataset
	{ 
		# get an example value of the dataset
		my $i = 1;
		while ($example eq undef) {
			# we want to avoid a non-value '.', so keep trying
			if ($data_table_ref->[$i][$index] ne '.') {
				# a non-null value, take it
				$example = $data_table_ref->[$i][$index];
			} 
			else {
				# a null value, proceed to next one
				$i++;
			}
		}
	}
	
	# Determine sort method, either numerical or alphabetical
	my $sortmethod; 
	if ($example =~ /[a-z]/i) { 
		# there are detectable letters
		$sortmethod = 'ascii';
	} 
	elsif ($example =~ /^\-?\d+\.?\d*$/) {
		# there are only digits, allowing for minus sign and a decimal point
		$sortmethod = 'numeric';
	} 
	else { 
		# unable to determine (probably alphanumeric), ask the user
		print " What is the sorting method? (a)scii or (n)umeric?   ";
		my $answer = <STDIN>;
		chomp $answer;
		if ($answer eq 'a') {
			$sortmethod = 'ascii';
		} 
		elsif ($answer eq 'n') {
			$sortmethod = 'numeric';
		} 
		else {
			warn " Unknown sort method. Nothing done.\n";
			return;
		}
	}
	
	
	# Remove the table header
	# this keeps the header out of the sorting process
	my $header = shift @{ $data_table_ref }; 
	# calculate our own temporary last_row index, since the main data value
	# is not valid because we moved the header out
	my $last_row = scalar @{ $data_table_ref } - 1;
	
	# Re-order the datasets
	# Directly sorting the @data array is proving difficult. It keeps giving me
	# a segmentation fault. So I'm using a different approach by copying the 
	# @data_table into a temporary hash.
		# put data_table array into a temporary hash
		# the hash key will the be dataset value, 
		# the hash value will be the reference the row data
	my %datahash;
	
	# reorder numerically
	if ($sortmethod eq 'numeric') {
		print " Sorting '", $main_data_ref->{$index}{'name'}, "' numerically...\n";
		for my $row (0..$last_row) {
			
			# get the value to sort by
			my $value = $data_table_ref->[$row][$index]; 
			
			# check for alphabet characters
			if ($value =~ /[abcdf-z]+/i) { 
				# check for any letter except for e, which may represent
				# exponent values in scientific notation
				warn "  Unable to numeric sort with alphabet characters in dataset!\n";
				return;
			}
			# check to see whether this value exists or not
			while (exists $datahash{$value}) {
				# add a really small number to bump it up and make it unique
				# this, of course, presumes that none of the dataset values
				# are really this small - this may be an entirely bad 
				# assumption!!!!! I suppose we could somehow calculate an 
				# appropriate value.... nah.
				# don't worry, we're only modifying the value used for sorting,
				# not the actual value
				$value += 0.00000001; 
			}
			
			# store the row data reference
			$datahash{$value} = $data_table_ref->[$row]; 
		}
		
		# re-fill the array based on the sort direction
		if ($direction eq 'i' or $direction eq 'I') { 
			# increasing sort
			my $i = 0; # keep track of the row
			foreach (sort {$a <=> $b} keys %datahash) {
				# put back the reference to the anonymous array of row data
				$data_table_ref->[$i] = $datahash{$_};
				$i++; # increment for next row
			}
		} 
		
		elsif ($direction eq 'd' or $direction eq 'D') { 
			# decreasing sort
			my $i = 0; # keep track of the row
			foreach (sort {$b <=> $a} keys %datahash) {
				# put back the reference to the anonymous array of row data
				$data_table_ref->[$i] = $datahash{$_};
				$i++; # increment for next row
			}
		}
		
		# restore the table header
		unshift @{ $data_table_ref }, $header;
		
		# summary prompt
		print " Data table sorted numerically by the contents of " .
			$main_data_ref->{$index}{'name'} . "\n";
		
	} 
	
	# reorder asciibetically
	elsif ($sortmethod eq 'ascii') {
		print " Sorting '", $main_data_ref->{$index}{'name'}, "' asciibetically...\n";
		for my $row (0..$last_row) {
			
			# get the value to sort by
			my $value = $data_table_ref->[$row][$index]; 
			
			# check to see if this is a unique value
			if (exists $datahash{$value}) { 
				# not unique
				my $n = 1;
				my $lookup = $value . sprintf("03%d", $n);
				# we'll try to make a unique value by appending 
				# a number to the original value
				while (exists $datahash{$lookup}) {
					# keep bumping up the number till it's unique
					$n++;
					$lookup = $value . sprintf("03%d", $n);
				}
				$datahash{$lookup} = $data_table_ref->[$row];
			} 
			else {
				# unique
				$datahash{$value} = $data_table_ref->[$row];
			}
		}
		
		# re-fill the array based on the sort direction
		if ($direction eq 'i' or $direction eq 'I') { 
			# increasing
			my $i = 0; # keep track of the row
			foreach (sort {$a cmp $b} keys %datahash) {
				# put back the reference to the anonymous array of row data
				$data_table_ref->[$i] = $datahash{$_};
				$i++; # increment for next row
			}
		} 
		
		elsif ($direction eq 'd' or $direction eq 'D') { 
			# decreasing
			my $i = 0; # keep track of the row
			foreach (sort {$b cmp $a} keys %datahash) {
				# put back the reference to the anonymous array of row data
				$data_table_ref->[$i] = $datahash{$_};
				$i++; # increment for next row
			}
		}
		
		# restore the table header
		unshift @{ $data_table_ref }, $header;
		
		# summary prompt
		print " Data table sorted asciibetically by the contents of '" .
			$main_data_ref->{$index}{'name'} . "'\n";
	}
	
	# remove any pre-existing sorted metadata since no longer valid
	for (my $i = 0; $i < $main_data_ref->{'number_columns'}; $i++) {
		if (exists $main_data_ref->{$i}{'sorted'}) {
			delete $main_data_ref->{$i}{'sorted'};
		}
	}
	
	# annotate metadata
	if ($direction =~ /i/i) {
		$main_data_ref->{$index}{'sorted'} = $sortmethod . "_increasing";
	}
	else {
		$main_data_ref->{$index}{'sorted'} = $sortmethod . "_decreasing";
	}
	
	return 1;
}



sub genomic_sort_function {
	# This will sort the entire data table by chromosome and start position
	
	# attempt to automatically identify the chromo and start indices
	my $chromo_i = find_column_index($main_data_ref, '^chr|seq|refseq');
	my $start_i = find_column_index($main_data_ref, '^start|position');
	
	# if unable to auto-identify columns, request from user
	unless (defined $chromo_i and defined $start_i) {
		my $line = " Please enter the index numbers for chromosome and start," .
			"\n separated by a comma   ";
		($chromo_i, $start_i) = _request_indices($line);
		unless (defined $chromo_i and defined $start_i) {
			warn " unknown index! nothing done\n";
			return;
		}
	}
	
	# load the data into a temporary hash
	# Directly sorting the @data array is proving difficult. It keeps giving me
	# a segmentation fault. So I'm using a different approach by copying the 
	# data_table into temporary hashes.
	# 
	# The datalines will be put into a hash of hashes: The first key will be 
	# the chromosome name, the second hash will be the start value.
	# 
	# To deal with some chromosomes that don't have numbers (e.g. chrM), we'll
	# use two separate hashes: one is for numbers, the other for strings
	# when it comes time to sort, we'll put the numbers first, then strings
	
	my %num_datahash;
	my %str_datahash;
	for my $row (1 .. $main_data_ref->{'last_row'}) { 
		
		my $startvalue = $data_table_ref->[$row][$start_i];
		
		# check for alphabet characters
		if ($startvalue =~ /[a-z]+/i) { 
			warn "  Unable to numeric sort with alphabet characters in start data!\n";
			return;
		}
		
		# put the dataline into the appropriate temporary hash
		if ($data_table_ref->[$row][$chromo_i] =~ /^(?:chr)?(\d+)$/) {
			# dealing with a numeric chromosome name
			# restricting to either chr2 or just 2 but not 2-micron
			my $chromovalue = $1;
			while (exists $num_datahash{$chromovalue}{$startvalue}) { 
				# if another item already exists at this location
				# add a really small number to bump it up and make it unique
				$startvalue += 0.001; 
			}
			$num_datahash{$chromovalue}{$startvalue} = $data_table_ref->[$row];
		} 
		else {
			# dealing with a non-numeric chromosome name
			my $chromovalue = $data_table_ref->[$row][$chromo_i];
			# use the entire chromosome name as key
			while (exists $str_datahash{$chromovalue}{$startvalue}) { 
				# if another item already exists at this location
				# add a really small number to bump it up and make it unique
				$startvalue += 0.001; 
			}
			$str_datahash{$chromovalue}{$startvalue} = $data_table_ref->[$row];
		}
	}
	
	
	# Now re-load the data array with sorted data
	# put the numeric chromosome data back first
	my $i = 1; # keep track of the row
	foreach my $chromovalue (sort {$a <=> $b} keys %num_datahash) {
		# first, numeric sort on increasing chromosome number
		foreach my $startvalue (
			sort {$a <=> $b} keys %{ $num_datahash{$chromovalue} } 
		) {
			# second, numeric sort on increasing position value
			$data_table_ref->[$i] = $num_datahash{$chromovalue}{$startvalue};
			$i++; # increment for next row
		}
	}
	# next put the string chromosome data back
	foreach my $chromovalue (sort {$a cmp $b} keys %str_datahash) {
		# first, ascii sort on increasing chromosome name
		foreach my $startvalue (
			sort {$a <=> $b} keys %{ $str_datahash{$chromovalue} } 
		) {
			# second, numeric sort on increasing position value
			$data_table_ref->[$i] = $str_datahash{$chromovalue}{$startvalue};
			$i++; # increment for next row
		}
	}
	
	# remove any pre-existing sorted metadata since no longer valid
	for (my $i = 0; $i < $main_data_ref->{'number_columns'}; $i++) {
		if (exists $main_data_ref->{$i}{'sorted'}) {
			delete $main_data_ref->{$i}{'sorted'};
		}
	}
	
	# annotate metadata
	$main_data_ref->{$chromo_i}{'sorted'} = 'genomic';
	$main_data_ref->{$start_i}{'sorted'} = 'genomic';
	
	print " Data table is sorted by genomic order\n";
	return 1;
}



sub toss_nulls_function {
	# Toss out datapoints (lines) that have a non-value in the specified dataset
	
	# generate the list of datasets to check
	my @order = _request_indices(
		" Enter one or more dataset index numbers to check for non-values\n   "); 
	unless (@order) {
		warn " No valid datasets! Nothing done!\n";
		return;
	}
	
	# Collection exception rule from commandline
	my $zero = $opt_zero;
	
	# Begin checking and tossing
	# We will walk through each row in the data table checking for non-values.
	# Any row that has a non-value will be skipped and ignored. Any row that 
	# does NOT have a non-value, i.e. good, will be copied into a NEW 
	# data table. At the end, this new data table will replace the old one.
	my @new_data_table;
	
	# Copy the header over
	push @new_data_table, $data_table_ref->[0];
	
	# Walk through the temp data array
	my $tosscount = 0; # count how many we tossed
	for my $row (1..$main_data_ref->{'last_row'}) {
		my $check = 0; # we will use a check variable
		
		# check for non-values in each of the requested datasets
		foreach my $index (@order) {
			if ($data_table_ref->[$row][$index] eq '.') {
				# my standard null value
				$check++;
			} 
			elsif ($data_table_ref->[$row][$index] eq undef) {
				# a true null value, these should've been converted to '.'
				$check++;
			} 
			elsif ($data_table_ref->[$row][$index] == 0) {
				# we have a 0 value, what to do?
				unless (defined $zero) {
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
		
		# decide whether to keep or toss
		if ($check > 0) { 
			# this data line fails the check
			# there is at least 1 non-value present in the requested columns
			# we will skip it
			$tosscount++;
		} 
		else { 
			# the check is good
			# copy the row data reference to the new array
			push @new_data_table, $data_table_ref->[$row];
		}
	}
	
	# Re-assign the new data table to the main data array
	$main_data_ref->{'data_table'} = \@new_data_table;
	$data_table_ref = $main_data_ref->{'data_table'};
	
	# re-calculate the last row index
	$main_data_ref->{'last_row'} = scalar @{ $data_table_ref } - 1;
	
	# update metadata
	foreach my $index (@order) {
		$main_data_ref->{$index}{'tossed'} = $tosscount . '_non_value_features';
	}
	
	# report
	if (scalar @order == 1) {
		# single dataset checked
		print " $tosscount data lines were removed that had null values in dataset "
			. $main_data_ref->{ $order[0] }{'name'} . "\n";
	}
	else {
		# generate a string of the searched data set names
		my $names = join ", ", map { $main_data_ref->{$_}{'name'} } @order;
		print " $tosscount data lines were removed that had null values in " .
			"datasets '$names'\n";
	}
	print " ", $main_data_ref->{'last_row'}, " data lines are remaining\n";
	return 1;
}






sub toss_duplicates_function {
	# Toss out datapoints (lines) that have a duplicate values
	
	# generate the list of datasets to check
	my @order = _request_indices(
		" Enter one or more dataset index numbers to check for duplicates\n   "); 
	unless (@order) {
		warn " No valid datasets! Nothing done!\n";
		return;
	}
	
	# initialize variables
	my %values2check;
	my @rows2toss;
	
	# check values
	foreach my $row (1 .. $main_data_ref->{'last_row'}) {
		
		# generate a value by concatanating all requested values
		my $value;
		foreach my $i (@order) {
			$value .= $data_table_ref->[$row][$i];
		}
		
		# check the value
		if (exists $values2check{$value}) {
			# yes, it exists, mark for destruction
			$values2check{$value} += 1;
			push @rows2toss, $row;
		}
		else {
			# nope, it's good
			$values2check{$value} = 1;
		}
	}
	
	# delete the duplicate rows
	foreach my $row (sort {$b <=> $a} @rows2toss) {
		# we sort from the bottom up so that we can find accurate row indexes
		unless (splice @{$data_table_ref}, $row, 1) {
			die " unrecoverable error! splicing data table failed!\n";
		}
	}
	
	# re-calculate the last row index
	$main_data_ref->{'last_row'} = scalar @{ $data_table_ref } - 1;
	
	# update metadata
	foreach my $index (@order) {
		$main_data_ref->{$index}{'tossed'} = scalar(@rows2toss) . 
			'_duplicate_features';
	}
	
	# print result
	print " tossed " . scalar(@rows2toss) . " duplicate features found in " . 
		"datasets " . join(', ', map { $main_data_ref->{$_}{'name'} } @order) . 
		"\n";
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
		" Enter one or more dataset index numbers to toss values $direction a value\n   "); 
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
	
	# Begin checking and tossing
	# We will walk through each row in the data table checking the values.
	# Any row whose datapoints do satisfy the test, i.e. bad, will be skipped. 
	# Any row whose datapoints do NOT satisfy the test, i.e. good, will be 
	# copied into a NEW data table. 
	# At the end, this new data table will replace the old one.
	my @new_data_table;
	
	# Copy the header over
	push @new_data_table, $data_table_ref->[0];
	
	# Walk through the temp data array
	my $tosscount = 0; # count how many we tossed
	for my $row (1..$main_data_ref->{'last_row'}) {
		my $check = 0; # we will use a check variable
		
		# check for non-values in each of the requested datasets
		if ($direction eq 'above') {
			# looking for values that exceed threshold
			
			foreach my $index (@order) {
				if ($data_table_ref->[$row][$index] eq '.') {
					# my standard null value, can't verify
					next;
				} 
				elsif ($data_table_ref->[$row][$index] > $threshold) {
					# this value exceeds the threshold
					$check++;
				} 
			}
		}
		
		else {
			# looking for values that are below the threshold
			
			foreach my $index (@order) {
				if ($data_table_ref->[$row][$index] eq '.') {
					# my standard null value, can't verify
					next;
				} 
				elsif ($data_table_ref->[$row][$index] < $threshold) {
					# this value exceeds the threshold
					$check++;
				} 
			}
		}
		
		# decide whether to keep or toss
		if ($check > 0) { 
			# this data line fails the check
			# there is at least 1 value present in the requested columns 
			# that exceeds the threshold value
			# we will skip it
			$tosscount++;
		} 
		else { 
			# the check is good
			# copy the row data reference to the new array
			push @new_data_table, $data_table_ref->[$row];
		}
	}
	
	# Re-assign the new data table to the main data array
	$main_data_ref->{'data_table'} = \@new_data_table;
	$data_table_ref = $main_data_ref->{'data_table'};
	
	# re-calculate the last row index
	$main_data_ref->{'last_row'} = scalar @{ $data_table_ref } - 1;
	
	# update metadata
	foreach my $index (@order) {
		$main_data_ref->{$index}{'tossed'} = $tosscount . '_lines_' . $direction 
			. '_threshold_' . $threshold;
	}
	
	# report
	if (scalar @order == 1) {
		# single dataset checked
		print " $tosscount data lines were removed that had values $direction threshold in "
			. $main_data_ref->{ $order[0] }{'name'} . "\n";
	}
	else {
		# generate a string of the searched data set names
		my $names = join ", ", map { $main_data_ref->{$_}{'name'} } @order;
		print " $tosscount data lines were removed that had values $direction " .
			"threshold in datasets '$names'\n";
	}
	print " ", $main_data_ref->{'last_row'}, " data lines are remaining\n";
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
	
	# request value
	my $value;
	if (defined $opt_target) {
		# command line option
		$value = $opt_target;
	}
	else {
		# interactively ask the user
		print " Enter the new value to convert nulls to  ";
		$value = <STDIN>;
		chomp $value;
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
		
		# reset minimum values
		if ($placement eq 'r' or $placement eq 'R') {
			# Replace the contents of the original dataset
			
			for my $i (1..$main_data_ref->{'last_row'}) {
				# check for valid numbers
				if ($data_table_ref->[$i][$index] eq '.') {
					# null value, need to change
					# change it in situ
					$data_table_ref->[$i][$index] = $value;
					$count++;
				}
				elsif ($data_table_ref->[$i][$index] == 0) {
					# zero value, what to do?
					unless (defined $zero) {
						# wasn't defined on the command line, so stop the program and ask the user
						print " Include 0 values to convert? y or n  ";
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
					
					if ($zero) {
						# change it in situ
						$data_table_ref->[$i][$index] = $value;
						$count++;
					}
				}
			}
			
			# update metadata
			$main_data_ref->{$index}{'null_value'} = $value;
			
			# results
			if ($count) {
				$total_count += $count;
				push @datasets_modified, $main_data_ref->{$index}{'name'};
			}
		} 
		
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Generate a new dataset
			
			# the new index position is equivalent to the number of columns
			my $new_position = $main_data_ref->{'number_columns'};
			
			# calculate new values
			for my $i (1..$main_data_ref->{'last_row'}) {
				# check for null values
				if ($data_table_ref->[$i][$index] eq '.') {
					# null value, need to change
					$data_table_ref->[$i][$new_position] = $value;
					$count++;
				} 
				elsif ($data_table_ref->[$i][$index] == 0) {
					# zero value, what to do?
					unless (defined $zero) {
						# wasn't defined on the command line, so stop the program and ask the user
						print " Include 0 values to convert? y or n  ";
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
					
					if ($zero) {
						# change it in situ
						$data_table_ref->[$i][$new_position] = $value;
						$count++;
					}
				}
				else {
					# acceptable
					$data_table_ref->[$i][$new_position] = 
						$data_table_ref->[$i][$index];
				}
			}
			
			# copy the medadata hash and annotate
			my $new_name;
			if ($function and $opt_name) {
				# automatic execution and new name was specifically given 
				$new_name = $opt_name;
			}
			else {
				$new_name = $main_data_ref->{$index}{'name'} . "_convert_nulls";
			}
			_generate_new_metadata(
				$index,
				$new_position,
				'null_value',
				$value,
				$new_name,
			);
			
			# results
			$total_count += $count;
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		} 
		
		else {
			warn " null values NOT changed; unknown placement request\n";
			return;
		}
	}
	
	# report results
	if (@datasets_modified) {
		my $string = " $total_count null values were converted for";
		$string .= $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= $placement =~ /n/i ? " and generated as new datasets\n" : "\n";
		print $string;
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
		if ($placement eq 'r' or $placement eq 'R') {
			# Replace the contents of the original dataset
			
			for my $i (1..$main_data_ref->{'last_row'}) {
				# check for valid numbers
				if ($data_table_ref->[$i][$index] eq '.') {
					# null value, cannot change
					next;
				}
				else {
					# change it in situ
					my $new_value;
					eval { $new_value = abs( $data_table_ref->[$i][$index] ) };
					if (defined $new_value) {
						$data_table_ref->[$i][$index] = $new_value;
						$count++;
					}
					else {
						$failed++;
					}
				} 
			}
			
			# update metadata
			$main_data_ref->{$index}{'convert'} = 'absolute';
			
			# results
			if ($count) {
				$total_count += $count;
				push @datasets_modified, $main_data_ref->{$index}{'name'};
			}
			$total_failed += $failed;
		} 
		
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Generate a new dataset
			
			# the new index position is equivalent to the number of columns
			my $new_position = $main_data_ref->{'number_columns'};
			
			# calculate new values
			for my $i (1..$main_data_ref->{'last_row'}) {
				# check for null values
				if ($data_table_ref->[$i][$index] eq '.') {
					# null value, cannot change
					$data_table_ref->[$i][$new_position] = '.';
				} 
				else {
					# acceptable
					my $new_value;
					eval { $new_value = abs( $data_table_ref->[$i][$index] ) };
					if (defined $new_value) {
						$data_table_ref->[$i][$new_position] = $new_value;
						$count++;
					}
					else {
						$data_table_ref->[$i][$new_position] = 
							$data_table_ref->[$i][$index];
						$failed++;
					}
				}
			}
			
			# copy the medadata hash and annotate
			my $new_name;
			if ($function and $opt_name) {
				# automatic execution and new name was specifically given 
				$new_name = $opt_name;
			}
			else {
				$new_name = $main_data_ref->{$index}{'name'} . "_absolute";
			}
			_generate_new_metadata(
				$index,
				$new_position,
				'convert',
				'absolute',
				$new_name,
			);
			
			# results
			$total_count += $count;
			$total_failed += $failed;
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		} 
		
		else {
			warn " values NOT changed; unknown placement request\n";
			return;
		}
	}
	
	# report results
	if (@datasets_modified) {
		my $string = " $total_count values were converted to absolute values for";
		$string .= $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= $placement =~ /n/i ? " and generated as new datasets\n" : "\n";
		print $string;
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
		if ($placement eq 'r' or $placement eq 'R') {
			# Replace the contents of the original dataset
			
			for my $i (1..$main_data_ref->{'last_row'}) {
				# check for valid numbers
				if ($data_table_ref->[$i][$index] eq '.') {
					# null value, nothing to do
					next;
				} 
				elsif ($data_table_ref->[$i][$index] < $value) {
					# current value below minimum value
					# change it in situ
					$data_table_ref->[$i][$index] = $value;
					$count++;
				}
			}
			
			# update metadata
			$main_data_ref->{$index}{'minimum_value'} = $value;
			
			# results
			if ($count) {
				$total_count += $count;
				push @datasets_modified, $main_data_ref->{$index}{'name'};
			}
		} 
		
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Generate a new dataset
			
			# the new index position is equivalent to the number of columns
			my $new_position = $main_data_ref->{'number_columns'};
			
			# calculate new values
			for my $i (1..$main_data_ref->{'last_row'}) {
				# check for valid numbers
				if ($data_table_ref->[$i][$index] eq '.') {
					# null value, nothing to do
					$data_table_ref->[$i][$new_position] = '.';
				} 
				elsif ($data_table_ref->[$i][$index] < $value) {
					# current value below minimum value
					$data_table_ref->[$i][$new_position] = $value;
					$count++;
				} 
				else {
					# acceptable
					$data_table_ref->[$i][$new_position] = 
						$data_table_ref->[$i][$index];
				}
			}
			
			# copy the medadata hash and annotate
			my $new_name;
			if ($function and $opt_name) {
				# automatic execution and new name was specifically given 
				$new_name = $opt_name;
			}
			else {
				$new_name = $main_data_ref->{$index}{'name'} . "_minimum_reset";
			}
			_generate_new_metadata(
				$index,
				$new_position,
				'minimum_value',
				$value,
				$new_name,
			);
			
			# results
			$total_count += $count;
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		} 
		
		else {
			warn " minimum value NOT reset; unknown placement request\n";
			return;
		}
	}
	
	# report results
	if (@datasets_modified) {
		my $string = " $total_count values were reset to a minimum value for";
		$string .= $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= $placement =~ /n/i ? " and generated as new datasets\n" : "\n";
		print $string;
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
		
		# reset maximum values
		if ($placement eq 'r' or $placement eq 'R') {
			# Replace the contents of the original dataset
			
			for my $i (1..$main_data_ref->{'last_row'}) {
				# check for valid numbers
				if ($data_table_ref->[$i][$index] eq '.') {
					# null value, nothing to do
					next;
				} 
				elsif ($data_table_ref->[$i][$index] > $value) {
					# current value below maximum value
					# change it in situ
					$data_table_ref->[$i][$index] = $value;
					$count++;
				}
			}
			
			# update metadata
			$main_data_ref->{$index}{'maximum_value'} = $value;
			
			# results
			if ($count) {
				$total_count += $count;
				push @datasets_modified, $main_data_ref->{$index}{'name'};
			}
		} 
		
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Generate a new dataset
			
			# the new index position is equivalent to the number of columns
			my $new_position = $main_data_ref->{'number_columns'};
			
			# calculate new values
			for my $i (1..$main_data_ref->{'last_row'}) {
				# check for valid numbers
				if ($data_table_ref->[$i][$index] eq '.') {
					# null value, nothing to do
					$data_table_ref->[$i][$new_position] = '.';
				} 
				elsif ($data_table_ref->[$i][$index] > $value) {
					# current value below maximum value
					$data_table_ref->[$i][$new_position] = $value;
					$count++;
				} 
				else {
					# acceptable
					$data_table_ref->[$i][$new_position] = 
						$data_table_ref->[$i][$index];
				}
			}
			
			# copy the medadata hash and annotate
			my $new_name;
			if ($function and $opt_name) {
				# automatic execution and new name was specifically given 
				$new_name = $opt_name;
			}
			else {
				$new_name = $main_data_ref->{$index}{'name'} . "_maximum_reset";
			}
			_generate_new_metadata(
				$index,
				$new_position,
				'maximum_value',
				$value,
				$new_name,
			);
			
			# results
			$total_count += $count;
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		} 
		
		else {
			warn " maximum value NOT reset; unknown placement request\n";
			return;
		}
	}
	
	# report results
	if (@datasets_modified) {
		my $string = " $total_count values were reset to a maximum value for";
		$string .= $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= $placement =~ /n/i ? " and generated as new datasets\n" : "\n";
		print $string;
	}
	return scalar(@datasets_modified);
}




sub log2_function {
	# this subroutine will convert dataset values to log2 space
	
	# request datasets
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter one or more dataset index numbers to convert to log2  "
		);
	}
	unless (@indices) {
		warn " unknown index number(s). nothing done\n";
		return;
	}
	
	# request placement
	my $placement = _request_placement();
	
	# process each index request
	my @datasets_modified; # a list of which datasets were modified
	my $total_count = 0; # total number of conversions done
	my $total_failed = 0;
	LOG2_LOOP: foreach my $index (@indices) {
		
		# check the log2 metadata status
		if (exists $main_data_ref->{$index}{'log2'} ) {
			if ( $main_data_ref->{$index}{'log2'} == 1 ) {
				warn " dataset $main_data_ref->{$index}{'name'} metadata " .
					"reports it is currently in log2 scale. Continue? y/n\n";
				my $response = <STDIN>;
				next LOG2_LOOP if $response =~ /n/i;
			}
		}
		
		# Placement dictates method
		my $count = 0; # conversion count
		my $failed = 0;
		if ($placement eq 'r' or $placement eq 'R') { 
			# Replace the current dataset
			
			# perform log2 conversion
			for my $i (1..$main_data_ref->{'last_row'}) {
				# walk through each value in the table
				
				# check the value contents and process appropriately
				if ($data_table_ref->[$i][$index] == 0) { 
					# cannot take log of 0
					$data_table_ref->[$i][$index] = '.'; 
					# change to null value
					$failed++;
				} 
				elsif ($data_table_ref->[$i][$index] eq '.') {
					# a null value, do nothing
					$failed++;
				} 
				else {
					# a numeric value, calculate the log2 value
					$data_table_ref->[$i][$index] = 
						log($data_table_ref->[$i][$index]) / log(2);
					$count++;
				}
			}
			
			# update metadata
			$main_data_ref->{$index}{'log2'} = 1;
			
			# results
			if ($count) {
				$total_count += $count;
				push @datasets_modified, $main_data_ref->{$index}{'name'};
			}
			$total_failed += $failed;
		}
		
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Generate a new dataset
			
			# the new index position is equivalent to the number of columns
			my $new_index = $main_data_ref->{'number_columns'};
			
			# perform log2 conversion
			for my $i (1..$main_data_ref->{'last_row'}) {
				# walk through each value in the table
				
				# check the value contents and process appropriately
				if ($data_table_ref->[$i][$index] == 0) { 
					# cannot take log of 0
					$data_table_ref->[$i][$new_index] = '.'; 
					# change to null value
					$failed++;
				} 
				elsif ($data_table_ref->[$i][$index] eq '.') {
					# a null value
					$data_table_ref->[$i][$new_index] = '.';
					$failed++;
				} 
				else {
					# a numeric value, calculate the log2 value
					$data_table_ref->[$i][$new_index] = 
						log($data_table_ref->[$i][$index]) / log(2);
					$count++;
				}
			}
			
			# annotate new metadata
			my $new_name = $main_data_ref->{$index}{'name'} . '_log2';
			_generate_new_metadata(
				$index,
				$new_index,
				'log2',
				1,
				$new_name,
			);
			
			# results
			$total_count += $count;
			$total_failed += $failed;
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		}
		
		else {
			# Unknown placement
			warn " log2 conversion NOT done; unknown placement request\n";
			return; # can't proceed with any index
		}
	}
	
	# report results
	if (@datasets_modified) {
		my $string = " $total_count values were converted to log2 scale for";
		$string .= $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= $placement =~ /n/i ? " and generated as new datasets\n" : "\n";
		print $string;
	}
	if ($total_failed) {
		print " $total_failed values could not be converted\n";
	}
	return scalar(@datasets_modified);
}



sub delog2_function {
	# this subroutine will convert a dataset from log2 to normal base10 numbers
	
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
	my @datasets_modified; # a list of which datasets were modified
	my $total_count = 0; # total number of conversions done
	my $total_failed = 0;
	LOG2_LOOP: foreach my $index (@indices) {
		
		# check the log2 metadata status
		if (exists $main_data_ref->{$index}{'log2'} ) {
			unless ( $main_data_ref->{$index}{'log2'} == 1 ) {
				warn " dataset $main_data_ref->{$index}{'name'} metadata " .
					"reports it is not in log2 scale. Continue? y/n\n";
				my $response = <STDIN>;
				next LOG2_LOOP if $response =~ /n/i;
			}
		}
		
		# Placement dictates method
		my $count = 0; # conversion count
		my $failed = 0;
		if ($placement eq 'r' or $placement eq 'R') { 
			# Replace the current dataset
			
			# perform log2 conversion
			for my $i (1..$main_data_ref->{'last_row'}) {
				# walk through each value in the table
				
				# check the value contents and process appropriately
				if ($data_table_ref->[$i][$index] eq '.') {
					# a null value, do nothing
					$failed++;
				} 
				else {
					# a numeric value, de-log2 the value
					$data_table_ref->[$i][$index] = 
						2 ** $data_table_ref->[$i][$index];
					$count++;
				}
			}
			
			# update metadata
			$main_data_ref->{$index}{'log2'} = 0;
			
			# results
			if ($count) {
				$total_count += $count;
				push @datasets_modified, $main_data_ref->{$index}{'name'};
			}
			$total_failed += $failed;
		}
		
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Generate a new dataset
			
			# the new index position is equivalent to the number of columns
			my $new_index = $main_data_ref->{'number_columns'};
			
			# perform log2 conversion
			for my $i (1..$main_data_ref->{'last_row'}) {
				# walk through each value in the table
				
				# check the value contents and process appropriately
				if ($data_table_ref->[$i][$index] eq '.') {
					# a null value
					$data_table_ref->[$i][$new_index] = '.';
					$failed++;
				} 
				else {
					# a numeric value, calculate the log2 value
					$data_table_ref->[$i][$new_index] = 
						2 ** $data_table_ref->[$i][$index];
					$count++;
				}
			}
			
			# annotate new metadata
			my $new_name = $main_data_ref->{$index}{'name'} . '_delog2';
			_generate_new_metadata(
				$index,
				$new_index,
				'log2',
				0,
				$new_name,
			);
			
			# results
			$total_count += $count;
			$total_failed += $failed;
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		}
		
		else {
			# Unknown placement
			warn " log2 de-conversion NOT done; unknown placement request\n";
			return; # can't proceed with any index
		}
	}
	
	# report results
	if (@datasets_modified) {
		my $string = " $total_count values were converted from log2 scale for";
		$string .= $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= $placement =~ /n/i ? " and generated as new datasets\n" : "\n";
		print $string;
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
		if ($placement eq 'r' or $placement eq 'R') {
			# Replace the contents of the original dataset
		
			for my $i (1..$main_data_ref->{'last_row'}) {
				if ($data_table_ref->[$i][$index] ne '.') {
					$data_table_ref->[$i][$index] = 
						sprintf $format_string, $data_table_ref->[$i][$index];
				}
			}
		
			# update metadata
			$main_data_ref->{$index}{'formatted'} = $positions;
		
			# results
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		} 
	
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Generate a new dataset
		
			# the new index position is equivalent to the number of columns
			my $new_position = $main_data_ref->{'number_columns'};
		
			# calculate new values
			for my $i (1..$main_data_ref->{'last_row'}) {
				if ($data_table_ref->[$i][$index] eq '.') {
					# working with a null value
					$data_table_ref->[$i][$new_position] = '.';
				} 
			
				else {
					# working with a real value
					# let's hope it is a number that may formatted
					$data_table_ref->[$i][$new_position] = 
						sprintf $format_string, $data_table_ref->[$i][$index];
				}
			}
		
			# copy the medadata hash and annotate
			my $new_name = $main_data_ref->{$index}{'name'} . '_formatted';
			_generate_new_metadata(
				$index,
				$new_position,
				'formatted',
				$positions,
				$new_name,
			);
		
			# results
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		} 
	
		else {
			warn " formatting not done; unknown placement request\n";
			return;
		}
	}
	
	# report results
	if (@datasets_modified) {
		my $string = $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= " formatted to $positions decimal positions";
		$string .= $placement =~ /n/i ? " and written as new datasets\n" : "\n";
		print $string;
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
	my %method2sub = (
		'mean'    => \&mean,
		'median'  => \&median,
		'min'     => \&min,
		'max'     => \&max,
		'stdev'   => \&stddevp,
		'sum'     => \&sum,
	);
	unless (exists $method2sub{$method}) {
		warn " unknown method. nothing done\n";
		return;
	}
	
	# identify log status
		# although we won't use this in the calculation....
	my $log;
	if (exists $main_data_ref->{ $indices[0] }{'log2'}) {
		$log = $main_data_ref->{ $indices[0] }{'log2'};
		foreach my $index (@indices) {
			if ($main_data_ref->{$index}{'log2'} != $log) {
				warn " unmatched log status between datasets! nothing done\n";
				return;
			}
		}	
	}
	
	# the new index position is equivalent to the number of columns
	my $new_position = $main_data_ref->{'number_columns'};
		
	# combine datasets
	my $failure_count = 0; # number of failures
	for my $row (1..$main_data_ref->{'last_row'}) {
		
		# collect row data
		my @data;
		foreach my $index (@indices) {
			if ($data_table_ref->[$row][$index] ne '.') {
				# only if not a null value
				push @data, $data_table_ref->[$row][$index];
			}
		}
		
		# combine the data
		if (@data) {
			# we have at least one datapoint, combine
			$data_table_ref->[$row][$new_position] = 
				&{ $method2sub{$method} }(@data);
		}
		else {
			# no datapoints, record null
			$data_table_ref->[$row][$new_position] = '.';
			$failure_count++;
		}
	}
	
	# generate new metadata
	my $new_name;
	if ($function and $opt_name) {
		# automatic execution and new name was specifically given 
		$new_name = $opt_name;
	}
	else {
		$new_name = "$method";
		foreach my $index (@indices) {
			$new_name .= '_' . $main_data_ref->{$index}{name};
		}
	}
	$main_data_ref->{$new_position}{'name'} = $new_name;
	$data_table_ref->[0][$new_position] = $new_name;
	$main_data_ref->{$new_position}{'index'} = $new_position;
	$main_data_ref->{$new_position}{'combine_method'} = $method;
	$main_data_ref->{$new_position}{'log2'} = $log if defined $log;
	$main_data_ref->{$new_position}{'datasets'} = join(',', 
		map { $main_data_ref->{$_}{'name'} } @indices);
	$main_data_ref->{'number_columns'} += 1;
	
	# finish
	print " combined datasets " . join(', ', 
		map { $main_data_ref->{$_}{'name'} } @indices) . " by $method" .
		" and generated a new dataset\n";
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
	
	elsif (
		exists $main_data_ref->{$numerator}{'log2'} and
		exists $main_data_ref->{$denominator}{'log2'}
	) {
		# great! both datasets have metadata log2 status flags
		# set log status based on mutual metadata flags
		if (
			$main_data_ref->{$numerator}{'log2'} == 1 and
			$main_data_ref->{$denominator}{'log2'} == 1
		) {
			# both datasets have the log2 flag set to true
			$log = 1;
		}
		elsif (
			$main_data_ref->{$numerator}{'log2'} == 0 and
			$main_data_ref->{$denominator}{'log2'} == 0
		) {
			# both datasets have the log2 flag set to false
			$log = 0;
		}
		elsif (
			$main_data_ref->{$numerator}{'log2'} != 
			$main_data_ref->{$denominator}{'log2'}
		) {
			# both datasets' log2 flags are not equal
			warn " datasets have different log status! nothing done";
			return;
		}
	}
	
	else {
		# when all else fails, ask the user
		print " Are these dataset values in log2 space? y/n   ";
		my $response = <STDIN>;
		chomp $response;
		if ($response =~ /^y/i) {
			$log = 1;
		} elsif ($response =~ /^n/i) {
			$log = 0;
		} else {
			warn " unknown response; nothing done\n";
			return;
		}
	}
	
	# the new index position is equivalent to the number of columns
	my $new_position = $main_data_ref->{'number_columns'};
		
	# generate ratio
	my $failure_count = 0; # number of failures
	for my $i (1..$main_data_ref->{'last_row'}) {
		
		# either value is null
		if (
			$data_table_ref->[$i][$numerator] eq '.' or 
			$data_table_ref->[$i][$denominator] eq '.'
		) {
			$data_table_ref->[$i][$new_position] = '.';
			$failure_count++;
		}
		
		# denominator is 0
		elsif ($data_table_ref->[$i][$denominator] == 0 and $log == 0) {
			# check if denominator is 0
			# we'll keep it simple, assign it a null value
			$data_table_ref->[$i][$new_position] = '.';
			$failure_count++;
		}
		
		# both values are good, determine ratio
		else {
			if ($log) {
				# perform a subtraction with log values
				$data_table_ref->[$i][$new_position] = 
					( $data_table_ref->[$i][$numerator] - 
					$data_table_ref->[$i][$denominator] );
			} 
			else {
				# perform a division with non-log values
				$data_table_ref->[$i][$new_position] =  
					( $data_table_ref->[$i][$numerator] / 
					$data_table_ref->[$i][$denominator] );
			}
		}
	}
	
	# annotate the new medadata hash
	my $new_name;
	if ($function and $opt_name) {
		# automatic execution and new name was specifically given 
		$new_name = $opt_name;
	}
	else {
		$new_name = $main_data_ref->{$numerator}{'name'} . '_' .
			$main_data_ref->{$denominator}{'name'} . '_ratio';
	}
	$main_data_ref->{$new_position}{'name'} = $new_name;
	$data_table_ref->[0][$new_position] = $new_name;
	$main_data_ref->{$new_position}{'index'} = $new_position;
	$main_data_ref->{$new_position}{'method'} = 'ratio';
	$main_data_ref->{$new_position}{'log2'} = $log;
	$main_data_ref->{$new_position}{'experiment'} = 
		$main_data_ref->{$numerator}{'name'};
	$main_data_ref->{$new_position}{'control'} = 
		$main_data_ref->{$denominator}{'name'};
	$main_data_ref->{'number_columns'} += 1;
	
	# print conclusion
	print " ratio between $main_data_ref->{$numerator}{'name'} and " . 
		"$main_data_ref->{$denominator}{'name'}\n generated as new dataset\n";
	if ($failure_count > 0) {
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
	
	# Check for log2 status
# 	if (
# 		(
# 			exists $main_data_ref->{$experiment_index}{'log2'} and
# 			$main_data_ref->{$experiment_index}{'log2'} == 1
# 		)
# 		or
# 		(
# 			exists $main_data_ref->{$control_index}{'log2'} and
# 			$main_data_ref->{$control_index}{'log2'} == 1
# 		)
# 	) {
# 		# one or both datasets are in the log2 space
# 		# we cannot proceed with log2 datasets
# 		warn " one or both datasets are in log2 space; cannot proceed\n";
# 		return;
# 	}

	
	# the new index position is equivalent to the number of columns
	my $new_position = $main_data_ref->{'number_columns'};
		
	# Generate difference
	my $failure_count = 0; # number of failures
	for my $i (1..$main_data_ref->{'last_row'}) {
		
		# either value is null
		if (
			$data_table_ref->[$i][$experiment_index] eq '.' or 
			$data_table_ref->[$i][$control_index] eq '.'
		) {
			$data_table_ref->[$i][$new_position] = '.';
			$failure_count++;
		}
		
		# both values are non-null, determine difference
		else {
			my $diff = $data_table_ref->[$i][$experiment_index] - 
				$data_table_ref->[$i][$control_index];
			
			if ($normalization) {
				# determine a normalized difference value
				my $sum = $data_table_ref->[$i][$experiment_index] +
					$data_table_ref->[$i][$control_index];
				if ($sum == 0) {
					# avoid pesky divide-by-0 errors, the difference is also 0
					$data_table_ref->[$i][$new_position] = 0;
				}
				else {
					$data_table_ref->[$i][$new_position] = $diff / sqrt($sum);
				}
			}
			else {
				# determine a straight difference value
				$data_table_ref->[$i][$new_position] = $diff;
			}
		}
	}
	
	# Generate the dataset's new name
	my $new_name;
	if ($function and $opt_name) {
		# this was an automatically executed function
		# and a new name was specified on the command line
		$new_name = $opt_name;
	} 
	elsif ($normalization) {
		$new_name = $main_data_ref->{$experiment_index}{'name'} . '_' .
			$main_data_ref->{$control_index}{'name'} . '_normdiff';
	}
	else {
		$new_name = $main_data_ref->{$experiment_index}{'name'} . '_' .
			$main_data_ref->{$control_index}{'name'} . '_diff';
	}
	
	# Determine the new method for the dataset's metadata
	my $new_method;
	if ($normalization) {
		$new_method = 'normalized_difference';
	}
	else {
		$new_method = 'difference';
	}
	
	# Annotate the new medadata hash
	$main_data_ref->{$new_position}{'name'} = $new_name;
	$data_table_ref->[0][$new_position] = $new_name;
	$main_data_ref->{$new_position}{'index'} = $new_position;
	$main_data_ref->{$new_position}{'method'} = $new_method;
	$main_data_ref->{$new_position}{'log2'} = 0;
	$main_data_ref->{$new_position}{'experiment'} = 
		$main_data_ref->{$experiment_index}{'name'};
	$main_data_ref->{$new_position}{'control'} = 
		$main_data_ref->{$control_index}{'name'};
	$main_data_ref->{'number_columns'} += 1;
	
	# Print conclusion
	print " normalized" if $normalization;
	print " difference between '$main_data_ref->{$experiment_index}{name}' and " . 
		"'$main_data_ref->{$control_index}{name}'\n generated as new dataset\n";
	if ($failure_count > 0) {
		print " $failure_count datapoints could not generate ratios\n";
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
				print "  median value for " . 
					$main_data_ref->{$index}{'name'} . " is $value\n";
			}
			elsif ($request_value eq 'mean') {
				$value = $stathash{'mean'};
				print "  mean value for " . 
					$main_data_ref->{$index}{'name'} . " is $value\n";
			}
			elsif ($request_value eq 'sum') {
				$value = $stathash{'sum'};
				print "  sum value for " . 
					$main_data_ref->{$index}{'name'} . " is $value\n";
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
		my $count = 0; # failed count  
		if ($placement eq 'r' or $placement eq 'R') {
			# Replace the contents of the original dataset
			
			for my $i (1..$main_data_ref->{'last_row'}) {
				# check for valid numbers
				if ($data_table_ref->[$i][$index] eq '.') {
					$count++;
					next;
				} else {
					$data_table_ref->[$i][$index] = 
						&$calculate($data_table_ref->[$i][$index], $value);
				}
			}
			
			# update metadata
			$main_data_ref->{$index}{$math} = $value;
			
			# print conclusion
			print " dataset $main_data_ref->{$index}{'name'} was $mathed "
				. "by value '$value'\n";
			$dataset_modification_count++;
		} 
		
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Generate a new dataset
			
			# the new index position is equivalent to the number of columns
			my $new_position = $main_data_ref->{'number_columns'};
			
			# calculate new values
			for my $i (1..$main_data_ref->{'last_row'}) {
				# check for valid numbers
				if ($data_table_ref->[$i][$index] eq '.') {
					$data_table_ref->[$i][$new_position] = '.';
					$count++;
				} else {
					$data_table_ref->[$i][$new_position] = 
						&$calculate($data_table_ref->[$i][$index], $value);
				}
			}
			
			# copy the medadata hash and annotate
			my $new_name;
			if ($function and $opt_name) {
				# automatic execution and new name was specifically given 
				$new_name = $opt_name;
			}
			else {
				$new_name = $main_data_ref->{$index}{'name'} . "_$mathed\_$value";
			}
			_generate_new_metadata(
				$index,
				$new_position,
				$math,
				$value,
				$new_name,
			);
			
			# print conclusion
			print " dataset '$main_data_ref->{$index}{'name'}' was $mathed by value " 
				. "'$value' and generated as a new dataset\n";
			$dataset_modification_count++;
		} 
		
		else {
			warn " $math NOT done; unknown placement request\n";
			return;
		}
		
		if ($count > 0) {
			print " $count datapoints could not be $mathed\n";
		}
	}
	
	# done
	return $dataset_modification_count;
}




sub convert_strand_to_sign {
	# this subroutine will convert one dataset to a signed dataset according to 
	# strand information
	
	# request datasets
	my @indices;
	if (@_) {
		# provided from an internal subroutine
		@indices = @_;
	}
	else {
		# otherwise request from user
		@indices = _request_indices(
			" Enter one or more dataset index numbers to convert to signed data  "
		);
	}
	unless (@indices) {
		warn " unknown index number(s). nothing done\n";
		return;
	}
	
	# identify the strand index
	my $strand_i = find_column_index($main_data_ref, 'strand|direction');
	unless (defined $strand_i) {
		# can't find it? ask the user
		$strand_i = _request_index(
			" Enter the index for the strand information column   ");
		if ($strand_i == -1) {
			warn " unknown strand index number. nothing done\n";
			return;
		}
	}
	
	# request placement
	my $placement = _request_placement();
	
	# proceed with the conversions
	my @datasets_modified; # a list of which datasets were modified
	my $total_count = 0; # total number of conversions done
	foreach my $index (@indices) {
		my $count = 0;
		
		if ($placement eq 'r' or $placement eq 'R') { 
			# Replace the current dataset
		
			# perform sign conversion
			for my $row (1..$main_data_ref->{'last_row'}) {
				# walk through each value in the table
			
				# check the value contents and process appropriately
				if ($data_table_ref->[$row][$index] == 0) { 
					# zero is inherently unsigned, skip
					next;
				} 
				elsif ($data_table_ref->[$row][$index] eq '.') {
					# a null value, do nothing
					next;
				} 
				else {
					# a presumed numeric value, change the sign
					if ($data_table_ref->[$row][$strand_i] =~ /^\-|r|c/i) {
						# looks like it is reverse: (minus), (r)everse, (c)rick
						# then prepend data value with a (minus)
						$data_table_ref->[$row][$index] = 
							-($data_table_ref->[$row][$index]);
						$count++;
					}
					elsif ($data_table_ref->[$row][$strand_i] !~ /^\+|1|f|w|0|\./i) {
						warn " unrecognized strand symbol for data row $row!\n";
						# do nothing for these
					}
					else {
						# do nothing for forward or no strand
					}
				}
			}
		
			# update metadata
			$main_data_ref->{$index}{'strand'} = 'signed';
		
			# results
			if ($count) {
				$total_count += $count;
				push @datasets_modified, $main_data_ref->{$index}{'name'};
			}
		}
	
		elsif ($placement eq 'n' or $placement eq 'N') {
			# Generate a new dataset

			# the new index position is equivalent to the number of columns
			my $new_index = $main_data_ref->{'number_columns'};
		
			# perform sign conversion
			for my $row (1..$main_data_ref->{'last_row'}) {
				# walk through each value in the table
			
				# check the value contents and process appropriately
				if ($data_table_ref->[$row][$index] == 0) { 
					# zero is inherently unsigned, skip
					$data_table_ref->[$row][$new_index] = 0;
				} 
				elsif ($data_table_ref->[$row][$index] eq '.') {
					# a null value, do nothing
					$data_table_ref->[$row][$new_index] = '.';
				} 
				else {
					# a presumed numeric value, change the sign
					if ($data_table_ref->[$row][$strand_i] =~ /^\-|r|c/i) {
						# looks like it is reverse: (minus), (r)everse, (c)rick
						# then prepend data value with a (minus)
						$data_table_ref->[$row][$new_index] = 
							-($data_table_ref->[$row][$index]);
						$count++;
					}
					elsif ($data_table_ref->[$row][$strand_i] =~ /^\+|1|f|w|0|\./i) {
						# forward or no strand, simply copy as is
						$data_table_ref->[$row][$new_index] = 
							$data_table_ref->[$row][$index];
					}
					else {
						warn " unrecognized strand symbol for data row $row!\n";
						# go ahead and copy the data
						$data_table_ref->[$row][$new_index] = 
							$data_table_ref->[$row][$index];
					}
				}
			}
		
			# annotate new metadata
			my $new_name = $main_data_ref->{$index}{'name'} . '_signed';
			_generate_new_metadata(
				$index,
				$new_index,
				'strand',
				'signed',
				$new_name,
			);
		
			# results
			$total_count += $count;
			push @datasets_modified, $main_data_ref->{$index}{'name'};
		}
		else {
			# Unknown placement
			warn " sign conversion NOT done; unknown placement request\n";
			return; # can't proceed with any index
		}
	}
	
	# report results
	if (@datasets_modified) {
		my $string = " $total_count values reverse strand values' signs were changed for";
		$string .= $#datasets_modified ? " datasets " : " dataset ";
		$string .= join(", ", @datasets_modified);
		$string .= $placement =~ /n/i ? " and generated as new datasets\n" : "\n";
		print $string;
	}
	return scalar(@datasets_modified);
}



sub merge_stranded_data {
	# this subroutine will merge two datasets representing forward and reverse
	# into one signed dataset

	# request datasets
	my @indices = _request_indices(
		" Enter the two indices for the forward and reverse strand datasets to merge   "
	);
	unless (@indices) {
		warn " unknown index number(s). nothing done\n";
		return;
	}
	unless (scalar @indices == 2) {
		warn " two indices must be provided! nothing done\n";
		return;
	}
	my ($f_index, $r_index) = @indices;
	
	# the new index position is equivalent to the number of columns
	my $new = $main_data_ref->{'number_columns'};
	
	# calculate the numbers for each line
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		# check all the possibilities
		
		if (
			$data_table_ref->[$row][$f_index] =~ /^[0|\.]$/ and
			$data_table_ref->[$row][$r_index] =~ /^[0|\.]$/
		) {
			# both datasets either have zero or no data
			# use forward dataset as the default
			$data_table_ref->[$row][$new] = $data_table_ref->[$row][$f_index];
		}
		
		elsif (
			$data_table_ref->[$row][$f_index] =~ /^[0|\.]$/ and
			$data_table_ref->[$row][$r_index] !~ /^[0|\.]$/
		) {
			# forward dataset only has zero or no data
			# use reverse dataset, reversing sign
			$data_table_ref->[$row][$new] = 
				-($data_table_ref->[$row][$r_index]);
		}
		
		elsif (
			$data_table_ref->[$row][$f_index] !~ /^[0|\.]$/ and
			$data_table_ref->[$row][$r_index] =~ /^[0|\.]$/
		) {
			# reverse dataset only has zero or no data
			# use forward dataset 
			$data_table_ref->[$row][$new] = $data_table_ref->[$row][$f_index];
		}
		
		else {
			# both datasets have data
			# perform simple subtraction
			$data_table_ref->[$row][$new] = $data_table_ref->[$row][$f_index] - 
				$data_table_ref->[$row][$r_index];
		}
	}
	
	# annotate new metadata
	my $new_name = $main_data_ref->{$f_index}{'name'} . '_&_' .
		$main_data_ref->{$r_index}{'name'} . '_merged';
	_generate_new_metadata(
		$f_index,
		$new,
		'strand',
		'signed',
		$new_name,
	);
		
	# finished
	print " $main_data_ref->{$f_index}{name} and $main_data_ref->{$r_index}{name}" .
		" stranded datasets were merged as signed data and" .
		" recorded as a new dataset\n";
	return 1;
}



sub number_function {
	# This subroutine will number the datapoints or lines
	
	# request dataset
	my $line = " The numbers will be entered as a dataset in front of " . 
		"which current dataset?\n Enter nothing to place the numbers " .
		"at the end.    ";
	my $index = _request_index($line);
	# a return of -1 indicates the numbers will be put at the end
	
	# number the datasets
	# we'll put the numbers at the end for now
	
	# the new index position is equivalent to the number of columns
	my $new_position = $main_data_ref->{'number_columns'};
	
	# calculate the numbers for each line
	for (my $i = 1; $i <= $main_data_ref->{'last_row'}; $i++) {
		# the line number is put at the end of the row
		$data_table_ref->[$i][$new_position] = $i;
	}
	
	# generate metadata
	my $new_name;
	if ($function and $opt_name) {
		$new_name = $opt_name;
	}
	else {
		$new_name = 'Numbers';
	}
	$main_data_ref->{$new_position}{'name'} = $new_name;
	$main_data_ref->{$new_position}{'index'} = $new_position;
	$data_table_ref->[0][$new_position] = 'Numbers';
	$main_data_ref->{'number_columns'} += 1;
	
	
	# move the dataset if necessary
	if ($index >= 0) {
		
		# first determine the new order for the datasets
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
		
		# then request to reorder the datasets
		my $success = reorder_function(@new_order);
		
		# print report
		if ($success) {
			print " Data points numbered and inserted before $index\n";
		}
		else {
			print " Data points were numbered but could not be moved from the end position\n";
		}
	} 
	
	else {
		# No re-ordering is requested
		print " Data points numbered and inserted at the end\n";
	}
		
	
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
	my $sumfile = write_summary_data(
		'data'         => $main_data_ref,
		'filename'     => $outfile,
		'startcolumn'  => $startcolumn,
		'endcolumn'    => $stopcolumn,
		'log'          => $opt_log,
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
	my $possible_name = $main_data_ref->{'basename'} . '_out.txt';
	
	# determine the export file name
	my $outfilename;
	if (defined $outfile) {
		# use the outfile name specified on the command line
		# note that this could be overwritten later if $modification > 0
		# but to allow for automated execution, we can't ask the user
		# for verification
		print " using export file name '$outfile'\n";
		$outfilename = $outfile;
	}
	elsif ($function) {
		# automatic execution, don't ask the user
		$outfilename = $possible_name;
	}
	else {
		# ask for new filename
		print " Exported file name? [$possible_name]  ";
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
	my $write_results = write_tim_data_file(
		'data'      => $main_data_ref,
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
	
	print " Preparing export for Treeview and Cluster analysis\n";
	
	# First check for previous modifications
	if ($modification) {
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
		print " Available dataset manipulations\n";
		print "   cg - median center features (genes)\n";
		print "   cd - median center datasets\n";
		print "   zd - convert dataset to Z-scores\n";
		print "   pd - convert dataset to percentile rank\n";
		print "   L2 - convert dataset to log2\n";
		print "   n0 - convert null values to 0\n";
		print " Enter the manipulation(s) in order of desired execution   ";
		my $answer = <STDIN>;
		chomp $answer;
		@manipulations = split /[,\s]+/, $answer;
	}
	
	
	### First, delete extraneous datasets or columns
	# the CDT format for Treeview expects a unique ID and NAME column
	# we will duplicate the first column
	unshift @datasets, $datasets[0];
	
	# perform a reordering of the columns
	reorder_function(@datasets);
	
	# rename the first two columns
	rename_function(0, 'ID');
	rename_function(1, 'NAME');
	
	# we now have just the columns we want
	# reset the dataset indices to what we currently have
	# name and ID should be index 0 and 1
	@datasets = (2 .. $main_data_ref->{'number_columns'} -1);
	
	
	### Second, perform dataset manipulations
	foreach (@manipulations) {
		if (/^cg$/i) {
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
		$outfile = $main_data_ref->{'path'} . 
			$main_data_ref->{'basename'} . '.cdt';
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
	for my $i (1..$main_data_ref->{'last_row'}) {
		
		# collect the datapoint values in each dataset for the feature  
		my @values;
		foreach my $dataset (@datasets) {
			if ($data_table_ref->[$i][$dataset] ne '.') {
				push @values, $data_table_ref->[$i][$dataset];
			}
		}
		unless (@values) {
			# no datapoints to work with
			next;
		}
		
		# determine median value
		my $median_value = median(@values);
		
		# adjust the datapoints
		foreach my $dataset (@datasets) {
			if ($data_table_ref->[$i][$dataset] ne '.') {
				$data_table_ref->[$i][$dataset] = 
					$data_table_ref->[$i][$dataset] - $median_value;
			}
		}
	}
	
	# annotate metadata
	foreach my $dataset (@datasets) {
		$main_data_ref->{$dataset}{'centered'} = 'median';
	}
	
	# print conclusion
	print " Datasets " . join(", ", @datasets) . " median centered\n";
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
	else {
		# interactively ask the user
		print " Enter the common value to be assigned in the new column   ";
		$value = <STDIN>;
		chomp $value;
	}
	
	# generate the new dataset metadata
	# the new index position is equivalent to the number of columns
	my $new_position = $main_data_ref->{'number_columns'};
	$main_data_ref->{$new_position} = {
		'name'    => $name,
		'index'   => $new_position,
	};
	$data_table_ref->[0][$new_position] = $name;
	$main_data_ref->{'number_columns'} += 1;
	
	# fill in the dataset
	for my $row (1..$main_data_ref->{'last_row'}) {
		$data_table_ref->[$row][$new_position] = $value;
	}
	
	# done
	print " Added new dataset '$name' at index $new_position with value '$value'\n";
	return 1;
}



sub subsample_function {
	# this will subsample the datasets to a target sum
	
	# identify the datasets to be subsampled
	my @datasets = _request_indices(
		" Enter the range of dataset indices to be subsampled  "
	);
	unless (@datasets) {
		warn " Unknown datasets. Nothing done.\n";
		return;
	}
	
	# identify the target sum
	my $target;
	if (defined $opt_target) {
		# use the command line specified target
		$target = $opt_target;
	}
	else {
		# request target from user
		print " Enter the specific target sum, or the index of a dataset to target  ";
		$target = <STDIN>;
		chomp $target;
		
	}
	
	# check whether the target is specific or a dataset index
	if ($target < $main_data_ref->{'number_columns'}) {
		# the number looks like it is a dataset index
		
		# validate just in case, then calculate the 
		if ( _validate_index_list($target) ) {
			
			# first check for log2
			if (
				exists $main_data_ref->{$target}{'log2'} and 
				$main_data_ref->{$target}{'log2'} == 1
			) {
				warn " target index is in log2 space! Nothing done.\n";
				return;
			}
			
			# determine the statistics for this dataset
			my %target_stats = _get_statistics_hash($target, 'y');
			unless (%target_stats) {
				warn " unable to get statistics on target dataset index; nothing done\n";
				return;
			}
			
			# set the actual target to the sum of this dataset
			$target = $target_stats{'sum'};
		}
		else {
			warn " nothing done\n";
			return;
		}
	}
	
	# Begin subsampling datasets
	my $subsample_count = 0; # count of successful subsamplings
	foreach my $dataset (@datasets) {
		
		# check for log2 space
		if (
			exists $main_data_ref->{$target}{'log2'} and 
			$main_data_ref->{$target}{'log2'} == 1
		) {
			warn " dataset is in log2 space! Nothing done.\n";
			next;
		}
		
		# get current stats
		my %dataset_stats = _get_statistics_hash($dataset, 'y');
		unless (%dataset_stats) {
			warn " unable to get statistics on dataset index $dataset!\n";
			next;
		}
		my $starting_sum = $dataset_stats{'sum'};
		
		# subsample
		my $new_index = _subsample_dataset($dataset, $starting_sum, $target);
		if ($new_index) {
			print " generated new subsampled dataset at index $new_index\n";
		}
		else {
			print " subsampling failed\n";
		}
		
		$subsample_count++;
	}
	
	# done
	return $subsample_count;
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
	my $write_results = write_tim_data_file(
		'data'      => $main_data_ref,
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
		'C' => "coordinate",
		'o' => "sort",
		'g' => "gsort",
		'N' => "null",
		'P' => "duplicate",
		'A' => "above",
		'B' => "below",
		'U' => "cnull",
		'L' => "absolute",
		'I' => "minimum",
		'X' => "maximum",
		'a' => "add",
		'u' => "subtract",
		'y' => "multiply",
		'v' => "divide",
		's' => "scale",
		'p' => "pr",
		'Z' => "zscore",
		'l' => "log2",
		'2' => "delog2",
		'f' => "format",
		'c' => "combine",
		'su' => "subsample",
		'r' => "ratio",
		'd' => "diff",
		'z' => "normdiff",
		'si' => "strandsign",
		'st' => "mergestrand",
		'e' => "center",
		'w' => "new",
		'Y' => "summary",
		'x' => "export",
		'W' => "rewrite",
		'T' => "treeview",
		'h' => "help",
		'q' => "quit",
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
		'coordinate'  => \&coordinate_function,
		'sort'        => \&sort_function,
		'gsort'       => \&genomic_sort_function,
		'null'        => \&toss_nulls_function,
		'duplicate'   => \&toss_duplicates_function,
		'above'       => \&toss_above_threshold_function,
		'below'       => \&toss_below_threshold_function,
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
		'log2'        => \&log2_function,
		'delog2'      => \&delog2_function,
		'format'      => \&format_function,
		'combine'     => \&combine_function,
		'subsample'   => \&subsample_function,
		'ratio'       => \&ratio_function,
		'diff'        => \&difference_function,
		'normdiff'    => \&normalized_difference_function,
		'strandsign'  => \&convert_strand_to_sign,
		'mergestrand' => \&merge_stranded_data,
		'center'      => \&center_function,
		'new'         => \&new_column_function,
		'summary'     => \&write_summary_function,
		'export'      => \&export_function,
		'rewrite'     => \&rewrite_function,
		'treeview'    => \&export_treeview_function,
		'help'        => \&print_online_help,
		'menu'        => \&print_menu,
	);
}


sub _print_headers {
	# this subroutine will print the names of the datasets from the column 
	# headers. It will also generate a hash of the dataset names with their 
	# index numbers as keys
	print " These are the datasets in the file\n";
	for (my $i = 0; $i < $main_data_ref->{'number_columns'}; $i++) {
		my $name = $main_data_ref->{$i}{'name'};
		print "\t$i\t$name\n";
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




sub _validate_index_list {
	# this subroutine will ensure that the list of index numbers represent 
	# actual columns, i.e., is each number < than the number of datasets in 
	# the data table
	
	foreach my $number (@_) {
		# check that each number represents a metadata hash 
		# and presumably dataset
		if (
			$number =~ /^\d+$/ and 
			exists $main_data_ref->{$number}
		) {
			# some integer number represented by a metadata hash is good
			# nothing to do, it is ok
		}
		else {
			# any other number is not good
			warn " requested index number '$number' doesn't represent a dataset!\n";
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
	my @values;
	for my $i ( 1..$main_data_ref->{'last_row'} ) {
		my $value = $data_table_ref->[$i][$index];
		if ($value eq '.') {
			# an internally represented null value
			next;
		} 
		elsif ($value == 0) {
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
			if ($zero) {
				push @values, $value;
			}
			# otherwise, we skip the 0 value and move on
		} 
		else {
			push @values, $value;
		}
	}
	
	# check that we have values
	unless (@values) {
		warn " no valid values collected from dataset index '$index'!";
		return;
	}
	
	my %statdata = statshash(@values);
	return %statdata;
}

sub _generate_new_metadata {
	# subroutine to copy a pre-existing metadata hash into a new metadata
	# hash and update it
	
	# collect arguments
	my ($old_index, $new_index, $key, $value, $name) = @_;
	
	# determine if new name was explicitly given during automatic execution
	if ($function and $opt_name) {
		$name = $opt_name;
	}
	
	# generate a new anonymous metadata array based on the old one	
	$main_data_ref->{$new_index} = { %{ $main_data_ref->{$old_index} } };
	
	# update the new key = value pair for the function
	$main_data_ref->{$new_index}{$key} = $value;
	
	# set and update the new name
	$main_data_ref->{$new_index}{'name'} = $name;
	$data_table_ref->[0][$new_index] = $name;
	
	# update index
	$main_data_ref->{$new_index}{'index'} = $new_index;
	
	# reset the number of columns
	$main_data_ref->{'number_columns'} += 1;
	
	return 1;
}

sub _subsample_dataset { 
	# a subroutine to sub sample a dataset
	# pass the index of the dataset, the current sum, and the target sum
	my ($dataset, $current_sum, $target_sum) = @_;
	
	# Copy the dataset first
	# we'll be using the reorder_function to duplicate
	my @new_order; 
	for (my $i = 0; $i < $main_data_ref->{'number_columns'}; $i++) {
		# keep all the current datasets
		push @new_order, $i;
	}
	push @new_order, $dataset; # duplicate the dataset at the end
	unless (reorder_function(@new_order) ) {
		warn " unable to duplicate dataset during subsampling!\n";
		return;
	}
	my $new_index = $main_data_ref->{'number_columns'} - 1; 
		# index of the duplicated dataset, which is at the end
	
	# Perform subsampling 
	my $difference = $current_sum - $target_sum;
	print " sub sampling dataset '" . $main_data_ref->{$dataset}{'name'} . 
		"' by $difference....\n";
	my $number = $main_data_ref->{'last_row'} - 1; # number to limit rand
	my $count = 0;
	my $limit = 50 * $number; # if it's truly random we should hit it 
	while ($difference) {
		# we need to match the sums, so we will subtract the sum difference
		# from the duplicate dataset one count at a time
		# this will be done from random datapoints across the dataset
		
		# pick a datapoint at random
		my $row = int rand($number);
		$row += 1; # to avoid the title line (row 0) and include the last row
		
		# check for valid value
		if ($data_table_ref->[$row][$new_index] eq '.') {
			# a null value, keep going
		}
		elsif ($data_table_ref->[$row][$new_index] > 0) {
			# a valid non-zero count
			# subtract one from the value
			$data_table_ref->[$row][$new_index] -= 1;
			
			# reduce difference to proceed to next round
			$difference -= 1;
		}
		else {
			# what else is in here? nothing we can do
		}
		
		# limit to avoid infinite loops
		if ($count >= $limit) {
			warn " reached limit of $limit iterations! Unable to sub sample!\n";
			return;
		}
		$count++;
	}
	
	# update metadata
	$main_data_ref->{$new_index}{'name'} .= '_subsampled';
	$data_table_ref->[0][$new_index] .= '_subsampled';
	$main_data_ref->{$new_index}{'subsample_target'} = $target_sum;
	
	# finished, return the index of the subsampled dataset
	return $new_index;
}



__END__

=head1 NAME

manipulate_datasets.pl

A progam to manipulate tab-delimited data files.

=head1 SYNOPSIS

manipulate_datasets.pl [--options ...] <input_filename> 

  Options:
  --in <input_filename>
  --func [ stat | reorder | delete | rename | number | coordinate | sort |
          gsort | null | duplicate | above | below | cnull | absolute | 
          minimum | maximum | add | subtract | multiply | divide | scale | 
          pr | zscore | log2 | delog2 | format | combine | subsample | 
          ratio | diff | normdiff | strandsign | mergestrand | center | 
          new | summary | export | rewrite | treeview ]
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
  
  stat
  reorder
  delete
  rename
  number
  coordinate
  sort
  gsort
  null
  duplicate
  above
  below
  cnull
  absolute
  minimum
  maximum
  add
  subtract
  multiply
  divide
  scale
  pr
  zscore
  log2
  delog2
  format
  combine
  subsample
  ratio
  diff
  normdiff
  strandsign
  mergestrand
  center
  new
  summary
  export
  rewrite
  treeview

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

=item B<stat> (menu option 't')

Print some basic statistics for a column, including mean, 
median, standard deviation, min, and max. If 0 values are present,
indicate whether to include them (y or n)

=item B<reorder> (menu option 'R')

The column may be rewritten in a different order. The new order 
is requested as a string of index numbers in the desired order. 
Also, a column may be deleted by skipping its number or duplicated
by including it twice.

=item B<delete> (menu option 'D')

One or more column may be selected for deletion.

=item B<rename> (menu option 'n')

Assign a new name to a column. For automatic execution, use the --name 
option to specify the new name. Also, for any automatically executed 
function (using the --func option) that generates a new column, the 
column's new name may be explicitly defined with this option.

=item B<number> (menu option 'b')

Assign a number to each datapoint (or line), incrementing from 1 
to the end. The numbered column will be inserted after the requested 
column index.

=item B<coordinate> (menu option 'C')

Generate a coordinate string from the chromosome, start, and, if 
present, stop coordinate values as a new column. The string will 
have the format "chr:start-stop" or "chr:start". This is useful 
in making unique identifiers or working with genome browsers.

=item B<sort> (menu option 'o')

The entire data table is sorted by a specific column. The first
datapoint is checked for the presence of letters, and the data 
table is then sorted either asciibetically or numerically. If the 
sort method cannot be automatically determined, it will ask. The 
direction of sort, (i)ncreasing or (d)ecreasing, is requested. 

=item B<gsort> (menu option 'g')

The entire data table is sorted by increasing genomic position, 
first by chromosome then by start position. These columns must exist 
and have recognizable names (e.g. 'chromo', 'chromosome', 'start').

=item B<null> (menu option 'N')

Toss datapoints (rows) that contain a null value in one or more 
columns. Some of the other functions may not work properly if
a non-value is present. If 0 values are present, indicate whether
to toss them (y or n). This may also be specified as a command line 
option using the --except flag.

=item B<duplicate> (menu option 'P')

Toss datapoints with duplicate values. One or more columns may be 
selected to search for duplicate values. Values are simply concatenated 
when multiple columns are selected. Rows with duplicated values are 
deleted, always leaving the first row.

=item B<above> (menu option 'A')

Toss datapoints with values that are above a certain threshold value. 
One or more columns may be selected to test values for the 
threshold. The threshold value may be requested interactively or 
specified with the --target option.

=item B<below> (menu option 'B')

Toss datapoints with values that are below a certain threshold value. 
One or more columns may be selected to test values for the 
threshold. The threshold value may be requested interactively or 
specified with the --target option.

=item B<cnull> (menu option 'U')

Convert null values to a specific value. One or more columns may 
be selected to convert null values. The new value may be requested 
interactively or defined with the --target option.  

=item B<absolute> (menu option 'L')

Convert signed values to their absolute value equivalents. One or 
more columns may be selected to convert.

=item B<minimum> (menu option 'I')

Reset datapoints whose values are less than a specified minimum 
value to the minimum value. One or more columns may be selected 
to reset values to the minimum. The minimum value may be requested 
interactively or specified with the --target option. 

=item B<maximum> (menu option 'X')

Reset datapoints whose values are greater than a specified maximum 
value to the maximum value. One or more columns may be selected 
to reset values to the maximum. The maximum value may be requested 
interactively or specified with the --target option. 

=item B<add> (menu option 'a')

Add a value to a column. A real number may be supplied, or the words
'mean', 'median', or 'sum' may be entered as a proxy for those statistical
values of the column. The column may either be replaced or added
as a new one. For automatic execution, specify the number using the
--target option.

=item B<subtract> (menu option 'u')

Subtract a value from a column. A real number may be supplied, or the words
'mean', 'median', or 'sum' may be entered as a proxy for those statistical
values of the column. The column may either be replaced or added
as a new one. For automatic execution, specify the number using the
--target option.

=item B<multiply> (menu option 'y')

Multiply a column by a value. A real number may be supplied, or the words
'mean', 'median', or 'sum' may be entered as a proxy for those statistical
values of the column. The column may either be replaced or added
as a new one. For automatic execution, specify the number using the
--target option.

=item B<divide> (menu option 'v')

Divide a column by a value. A real number may be supplied, or the words
'mean', 'median', or 'sum' may be entered as a proxy for those statistical
values of the column. The column may either be replaced or added
as a new one. For automatic execution, specify the number using the
--target option.

=item B<scale> (menu option 's')

A column may be a median scaled as a means of normalization 
with other columns. The current median of the column requested is
presented, and a new median target is requested. The column may 
either be replaced with the median scaled values or added as a new 
column. For automatic execution, specify the new median target 
with the --target option.

=item B<pr> (menu option 'p')

A column may be converted to a percentile rank, whereby the
data values are sorted in ascending order and assigned a new value 
from 0..1 based on its rank in the sorted order. The column may 
either be replaced with the percentile rank or added as a new
column. The original order of the column is maintained.

=item B<zscore> (menu option 'Z')

Generate a Z-score or standard score for each value in a column. The
Z-score is the number of standard deviations the value is away from
the column's mean, such that the new mean is 0 and the standard 
deviation is 1. Provides a simple method of normalizing columns
with disparate dynamic ranges.

=item B<log2> (menu option 'l')

A column may be converted to log base 2 values. The column
may either be replaced with the log2 values ar added as a new 
column.

=item B<delog2> (menu option '2')

A column that is currently in log2 space may be converted back to
normal base10 numbers. The column may either be replaced with the 
new values or added as a new column.

=item B<format> (menu option 'f')

Format the numbers of a column to a given number of decimal places. 
An integer must be provided. The column may either be replaced or 
added as a new column. For automatic execution, use the --target 
option to specify the number decimal places.

=item B<combine> (menu option 'c')

Mathematically combine the data values in two or more columns. The 
methods for combining the values include mean, median, min, max, 
stdev, or sum. The method may be specified on the command line 
using the --target option. The combined data values are added as a 
new column.

=item B<subsample> (menu option 'su')

Subsample an enumerated dataset. Datapoints within a column are chosen 
randomly and the values reduced by one until the target sum is reached. 
This assumes an enumerated dataset, e.g. tag counts from
Next-gen sequencing, where the sum of tags from two or more columns
must be normalized to each other. This function assumes that the data
are positive integers; log2 columns will not be used. Provide a
target sum to reach using the --target option. If this number is less
than or equal to the number of columns in the file, then it is
assumed to be an index number and the sum of that column will be used
as the target. The subsampled dataset is always added as new column.

=item B<ratio> (menu option 'r')

A ratio may be generated between two columns. The experiment and 
control columns are requested and the ratio is added as a new
column. For log2 numbers, the control is subtracted from the
experiment. The log2 status is checked in the metadata for the 
specified columns, or may be specified as a command line option, or
asked of the user.

=item B<diff> (menu option 'd')

A simple difference is generated between two existing columns. The 
values in the 'control' column are simply subtracted from the 
values in the 'experimental' column and recorded as a new column.
For enumerated columns (e.g. tag counts from Next Generation 
Sequencing), the columns should be subsampled to equalize the sums 
of the two columns. The indices for the experimental and control columns 
may either requested from the user or supplied by the --exp and 
--con command line options. 

=item B<normdiff> (menu option 'z')

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

=item B<strandsign> (menu option 'si')

Convert a column's values to signed data according to strand. Forward 
strand data is positive, and reverse strand is negative. This function 
may not be appropriate with logarithmic or other datasets that include 
negative values. Provide the index of a single column to convert. A 
second column should provide the strand information and be labeled 
with a name that includes either 'strand' or 'direction'. Strand 
information may include 1, -1, +, -, f, r, w, c, ., or 0. The stranded 
data may overwrite the data or written to a new column.

=item B<mergestrand> (menu option 'st')

Merge two stranded columns into a single new column, with the forward 
dataset represented as a positive value, and the reverse dataset as a 
negative value. Datapoints which contain values on both strands are 
recorded as a simple difference (forward - reverse). This function may 
not be appropriate with logarithmic or other datasets that include 
negative values. Provide the indices of the two columns as forward, 
reverse. Metadata is not checked for validity. 

=item B<center> (menu option 'e')

Center normalize the datapoints in a row by subtracting the mean or
median of the datapoints. The range of columns is requested or 
provided by the --index option. Old values are replaced by new 
values. This is useful for visualizing data as a heat map, for example.

=item B<new> (menu option 'w')

Generate a new column which contains an identical value for 
each datapoint (row). The value may be either requested interactively or 
supplied using the --target option. This function may be useful for 
assigning a common value to all of the data points before joining the 
data file with another.

=item B<summary> (menu option 'y')

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

=item B<export> (menu option 'x')

Export the data into a simple tab-delimited text file that contains no 
metadata or header information. Non-values '.' are converted to  
true nulls. If an output file name is specified using the --outfile 
option, it will be used. Otherwise, a possible filename will be 
suggested based on the input file name. If any modifications are 
made to the data structure, a normal data file will still be written. 
Note that this could overwrite the exported file if the output file name
was specified on the command line, as both file write subroutines will 
use the same name!

=item B<treeview> (menu option 'T')

Export the data to the CDT format compatible with both Treeview and 
Cluster programs for visualizing and/or generating clusters. Specify the 
columns containing a unique name and the columns to be analyzed (e.g. 
--index <name>,<start-stop>). Extraneous columns are removed. 
Additional manipulations on the columns may be performed prior to 
exporting. These may be chosen interactively or using the codes 
listed below and specified using the --target option.
  
  cg - median center features (rows)
  cd - median center datasets (columns)
  zd - convert columns to Z-scores
  pd - convert columns to percentile ranks
  L2 - convert values to log2
  n0 - convert nulls to 0.0

A simple Cluster data text file is written (default file name 
"<basename>.cdt"), but without the GWEIGHT column or EWEIGHT row. The 
original file will not be rewritten.

=item B<rewrite> (menu option 'W')

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
