#!/usr/bin/perl
 
# A script to process Agilent microarrays.
# This will take as many Agilent feature extraction files as you give it
# It will optionally quantile normalize the data, average multiple probe features,
# optionally combine the red and green values of the replicate files, and perform
# linear regression to determine Pearson correlation coefficients between the replicates.
# It will output the data in a tab delimited text file listed by probe name.

use strict;
use Getopt::Long;
use Statistics::Lite qw(mean median);
use Statistics::Descriptive;
# use GD::Graph::points;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper qw(
	write_tim_data_file
);
my $VERSION = '1.0.0';

print "\n A script to process Agilent microarray files\n\n";

### Quick help
unless (@ARGV) { # when no command line options are present
	print " Usage for $0
  Required:
  --in filename or --inlist filename
  --out filename
  --ratio [rg, gr, g, r, mix] or --inlist filename
  Optional:
  --median integer
  --nlist filename
  --separate
  --(no)norm
  --(no)log
  --processed
  --help\n";
	exit;
}



### Get command line options
my (
	$in_file_list,
	$outfile,
	$ch_ratio,
	$median_target,
	$norm_list,
	$separate,
	$quant,
	$log,
	$processed,
	$help,
	$print_version,
);
my @in_files;
GetOptions( 
	'in=s'       => \@in_files, # array of input files
	'inlist=s'   => \$in_file_list, # a recipe list of input files
	'out=s'      => \$outfile, # output file name
	'ratio=s'    => \$ch_ratio, # channel ratio
	'median=i'   => \$median_target, # median scale target value
	'nlist=s'    => \$norm_list, # filename for list of control probes to normalize
	'separate'   => \$separate, # quantile normalize experiment and control separately
	'norm!'      => \$quant, # do not quantile normalize
	'log!'       => \$log, # convert to log2
	'processed'  => \$processed, # use processed values rather than raw
	'help'       => \$help, # print help
	'version'    => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

if ($help) {
	help(); 
	exit;
}
# Print version
if ($print_version) {
	print " Biotoolbox script process_agilent.pl, version $VERSION\n\n";
	exit;
}

unless (@in_files or $in_file_list) {die " No input file name(s) specified! use --help\n"};
unless ($outfile) {die " No output file name specified! use --help\n"};
unless ($ch_ratio or $in_file_list) {die " No channel ratio specified! use --help\n"};
unless (defined $quant) {$quant = 1}
unless (defined $log) {$log = 0};


### Define global data containing variables
my %channellist; # a hash to record the channel ratios: first key=file name, second key = red or green, value=code for experiment (e) or control (c) or none (n)
my %evalues; # a hash of arrays to keep the experimental values: key=feature_number, value=array of red values
my %cvalues; # a hash of arrays to keep the control values: key=feature_number, value=array of green values
			  # the value arrays will represent the replicate values from each input file
my %probes; # a hash of arrays to keep the probe names and corresponding feature numbers: key=probe_id, value=feature_number(s)
my @meanslist; # keep an array of arrays of the probe mean experimental & control values for linear regression analysis
my %output_data; # a complex hash for the output data
my @report; # output array for reporting statistical analysis
my $e_norm_target; # the value for normalizing experimental values
my $c_norm_target; # the value for normalizing control values


### Prepare output report
my $dt = localtime(time); # get the date and time
push @report, "Processing Agilent files on $dt\n\n";



### Determine list of input files and experimental and control channels
# the list of input files may be provided directly by command line arguments
# or provided in a recipe list in a file specified by a command line argument

if ($in_file_list) {
	# load the recipe list file
	# this will fill up the in_file array
	# the recipe also includes the channel ratio directions
	load_infile_list($in_file_list);
}
else {
	# input file names were specified on the command line
	# determine the channel ratio direction for the files
	determine_channel_ratio(); 
}

@in_files = sort {$a cmp $b} (@in_files); # sort the filenames, just so they're in neat little order
# then collect the names of the replicates and their file origins
my @e_name_list; # an array for the names of the experiment data set files 
my @c_name_list;
foreach my $file (@in_files) {
	# we want to record the file name and channel for each experiement and control set
	# using the information in the channellist hash
	
	if ($channellist{$file}{'red'} eq 'e') {
		push @e_name_list, "$file:red";
	} elsif ($channellist{$file}{'red'} eq 'c') {
		push @c_name_list, "$file:red";
	}
	if ($channellist{$file}{'green'} eq 'e') {
		push @e_name_list, "$file:green";
	} elsif ($channellist{$file}{'green'} eq 'c') {
		push @c_name_list, "$file:green";
	}
}
unless (@c_name_list) {
	# if we only have experimental data and no control data
	# then set the separate flag to process the experimental only
	$separate = 1;
}



### Load the data from the Agilent files
load_my_agilent_files();



### Perform quantile normalization on the data
if ($quant) {
	print " Performing quantile normalization....\n";
	perform_normalization();
}



### Perform median scaling
if ($norm_list) {
	print " Normalizing values to the probe list in $norm_list....\n";
	get_normalization_probes();
}
if ($median_target or $norm_list) {
	print " Median scaling and/or normalizing values....\n";
	perform_scaling_and_normalization();
}

## testing output
#  {
#  	open OUTFILE, ">check_output2.txt";
#  	foreach (sort {$a <=> $b} keys %evalues) {
#  		print OUTFILE "$_ exp: " . join(", ", @{ $evalues{$_} }) . " con: " . join(", ", @{ $cvalues{$_} }) . "\n";
#  	}
#  	close OUTFILE;
#  }


### Preparing the data sets for output
if (@c_name_list) { # there is both experimental and control data values
	print " Combining replicate data sets for both experimental and control data....\n";
	push @report, "\nCombining the replicate data sets for both experimental and control data\n";
	output_combined_records();
} else { # there is only experimental data values (single channel)
	print " Combining replicate data sets for experimental data....\n";
	push @report, "\nCombining the replicate data sets for experimental data\n";
	output_single_records();
}
if ($log) {
	push @report, "Converting combined data sets to log base 2\n";
	print " Converting combined data sets to log base 2\n";
	# this is actually done when preparing the data sets for output above
}


### Linear regression analysis
print " Performing linear regression analysis....\n";
# calculate statistics on the combined experimental & control values
if (%cvalues) {
	perform_regression_on_means();
}
# calculate statistics pairwise among replicates
perform_regression_on_replicates();

# 
# ### Generating MA plot
# if (%cvalues) {
# 	# generate a MA plot between experimental and control values
# 	# this really isn't all that useful....
# 	print " Generating MA plot....\n";
# 	generate_ma_plot();
# }


### Write output
$outfile =~ s/\.txt$//; # strip extension if present
# write main output
my $written_file = write_tim_data_file( {
	# we will write a tim data file
	# appropriate extensions and compression should be taken care of
	'data'     => \%output_data,
	'filename' => $outfile,
} );
if ($written_file) {
	print " Wrote data file '$written_file'\n";
}
else {
	print " unable to write data file!\n";
}
push @report, "\nWrote output file to $written_file\n\n\n";

# write the report file
open FILE, ">$outfile\_report.txt";
print FILE @report;
close FILE;


print "I'm finished!!!!\n";



##########  Subroutines  ############
################################################################################

#### Help 
sub help { # subroutine to print the online help documentation
	print "

This program will process raw Agilent feature extraction files. 

It will take as many Agilent files as you give it. It will quantile normalize
the data, calculate the experimental/control ratio, and report the normalized 
data and ratio for each experimental probe on the array in a text file. The 
experimental and control channels may be quantile normalized either together or
separately. The data may be median scaled to a desired value and converted to 
log2 values. A subset of the probes on the array may be used as normalization 
controls. 

The program will determine correlation and r-squared values between datasets 
and report them. It will write out a tab-delimited data file consisting of 
normalized data values for experiment and control, and the experiment/control
ratio for each probeID. It will also write out a separate report text file.

Command line options for process_agilent.pl
  --in       Specify an input file name. Use repeatedly for each replicate 
             file to be specified
  --out      Specify the output file name
  --ratio    Specify the experiment/control ratio. Acceptable values include
             - rg     use red/green ( Cy5/Cy3 )
             - gr     use green/red ( Cy3/Cy5 )
             - g      green ( Cy3 ) only, e.g. HybMap arrays
             - r      red ( Cy5 ) only
             - mix    for dye swapped or complex experiments; interactively 
                      choose which channels are which for each input file
  --nlist    Optionally specify the name of a text file containing probe names
             to normalize the array against. One probe ID per line in the file.
  --median   Optionally specify a number to median scale values 
  --separate Quantile normalize experiment and control channels separately
  --nonorm   Do NOT quantile normalize
  --log      Convert values to log base 2
  --help     This help text
  ";
}



#### Loading an infile recipe list
sub load_infile_list {
	# this will load a processing recipe file
	# a tab delimited text file giving the names of the input files and the 
	# the usage for the red and green channels
	# one row per file
	# line consists of file_name, red_channel, green_channel
	# each channel must be one of three values: (e)xperiment, (c)ontrol, or (n)one
	my $listfile = shift;
	open LISTFILE, $listfile or die("unable to read recipe file '$listfile'!\n");
	while (my $line = <LISTFILE>) {
		if ($line =~ /^#/) {next} # skip any comment lines
		chomp $line;
		my ($file, $red, $green) = split /\t/, $line;
		if (-e $file) {
			# the input file exists
			
			# process the red channel request
			if ($red =~ /^([ecn])/) {
				# red channel matches acceptable usage
				# remember the first letter
				my $channel = $1;
				if (exists $channellist{$file}{'red'}) {
					# this file has been listed before
					if ($channellist{$file}{'red'} eq 'n' and 
						$channel ne 'n'
					) {
						# the channel was previously listed as not being used
						# and now is going to be used
						$channellist{$file}{'red'} = $channel;
					}
					elsif ($channellist{$file}{'red'} ne 'n' and 
						$channel eq 'n'
					) {
						# the channel was previously listed for use
						# and current use is none
						# this is ok, nothing to do
					}
					else {
						die "recipe file '$listfile' has errors! Input file " .
							"'$file' listed twice!\n";
					}
				}
				else {
					# this file has not been listed before
					$channellist{$file}{'red'} = $channel;
				}
				
			}
			else {
				die "recipe file '$listfile' has unrecognizable red channel value for file '$file'!\n";
			}
			
			# process the green channel request
			if ($green =~ /^([ecn])/) {
				# green channel matches acceptable usage
				# remember the first letter
				my $channel = $1;
				if (exists $channellist{$file}{'green'}) {
					# this file has been listed before
					if ($channellist{$file}{'green'} eq 'n' and 
						$channel ne 'n'
					) {
						# the channel was previously listed as not being used
						# and now is going to be used
						$channellist{$file}{'green'} = $channel;
					}
					elsif ($channellist{$file}{'green'} ne 'n' and 
						$channel eq 'n'
					) {
						# the channel was previously listed for use
						# and current use is none
						# this is ok, nothing to do
					}
					else {
						die "recipe file '$listfile' has errors! Input file " .
							"'$file' listed twice!\n";
					}
				}
				else {
					# this file has not been listed before
					$channellist{$file}{'green'} = $channel;
				}
			}
			else {
				die "recipe file '$listfile' has unrecognizable green channel value for file '$file'!\n";
			}
		}
		
		else {
			# the file does not exist!
			die "input file '$file' cannot be found!\n";
		}
	}
	close LISTFILE;
	
	# add the list of input files to the input file array
	push @in_files, keys %channellist;
}





#### Determining the ratio direction (red/green or green/red or mixed)
sub determine_channel_ratio {
	# we're recording the channel list code in a hash of hash
	# the first key is the file name
	# the second key is the channel, red or green
	# the values are letters representing (e)xperiment, (c)ontrol, or (n)one
	
	if ($ch_ratio eq 'rg') {
		foreach (@in_files) {
			$channellist{$_}{'red'} = 'e';
			$channellist{$_}{'green'} = 'c';
		}
		push @report, "Using red channel for experiment and green channel for control\n";
	} 
	elsif ($ch_ratio eq 'gr') {
		foreach (@in_files) {
			$channellist{$_}{'red'} = 'c';
			$channellist{$_}{'green'} = 'e';
		}
		push @report, "Using green channel for experiment and red channel for control\n";
	} 
	elsif ($ch_ratio eq 'g') {
		foreach (@in_files) {
			$channellist{$_}{'red'} = 'n';
			$channellist{$_}{'green'} = 'e';
		}
		push @report, "Using green channel only for experiment\n";
	} 
	elsif ($ch_ratio eq 'r') {
		foreach (@in_files) {
			$channellist{$_}{'red'} = 'e';
			$channellist{$_}{'green'} = 'n';
		}
		push @report, "Using red channel only for experiment\n";
	} 
	elsif ($ch_ratio eq 'mix') {
		# an interactive session to request the ratio for each input file
		push @report, "Using different experiment/control channel ratios for each input file:\n";
		print "\n  You requested a mix of channel ratios. These are your input files.\n";
		print "  For each file and channel, type 'e' for experiment, 'c' for control,\n  or 'n' for none\n\n";
		my %human_readable = ( # a quickie hash for letter abbreviations and human readable names
				'e' => 'experimental',
				'c' => 'control',
				'n' => 'none'
		);
		foreach (@in_files) {
			# print the name of each file, then print channel name and request a 
			# response for each channel
			print "\t$_\n\t  red  ";
			my $r_answer = <STDIN>;
			chomp $r_answer;
			until ( exists $human_readable{$r_answer} ) {
				print "\t  unknown response, try again   ";
				$r_answer = <STDIN>;
				chomp $r_answer;
			}
			print "\t  green  ";
			my $g_answer = <STDIN>;
			chomp $g_answer;
			until ( exists $human_readable{$g_answer} ) {
				print "\t  unknown response, try again   ";
				$g_answer = <STDIN>;
				chomp $g_answer;
			}
			
			$channellist{$_}{'red'} = $r_answer;
			$channellist{$_}{'green'} = $g_answer;
			push @report, "  File $_ red channel is $human_readable{$r_answer} and green channel is $human_readable{$g_answer}\n";
			
		}
	} else {
		die "unknown channel ratio value $ch_ratio!\n";
	}
}




#### Loading Files
sub load_my_agilent_files {
	foreach my $myfile (@in_files) {
		print " Reading input file $myfile....\n";
		open INFILE, "$myfile" or die("unable to open file $myfile!!!\n");
		
		# column indices values
		my ($featurenum_i, $controltype_i, $probename_i, $gsig_i, $rsig_i); 
		my $feature_count = 0;
		
		# walk through the file
		FILELOOP: while (my $line = <INFILE>) {
			
			# Determine what line we're working with based on whether column
			# indices have been identified yet or not
			
			if (defined $featurenum_i) {
				# column indices have been defined
				# now working with data lines
				my @data = split /\t/, $line;
				if ($data[$controltype_i] == 0) { # only take the experimental probe values
					
					# Need to record the featurenumber for each probe name
					push @{ $probes{$data[$probename_i]} }, $data[$featurenum_i];
					
					# collect red channel values
					if ($channellist{$myfile}{'red'} eq 'e') { 
						push @{ $evalues{$data[$featurenum_i]} },  $data[$rsig_i];
					} 
					elsif ($channellist{$myfile}{'red'} eq 'c') {
						push @{ $cvalues{$data[$featurenum_i]} },  $data[$rsig_i];
					}
					
					# collect green channel values
					if ($channellist{$myfile}{'green'} eq 'e') { 
						push @{ $evalues{$data[$featurenum_i]} },  $data[$gsig_i];
					} 
					elsif ($channellist{$myfile}{'green'} eq 'c') {
						push @{ $cvalues{$data[$featurenum_i]} },  $data[$gsig_i];
					}
					
					# count features
					$feature_count++;
				}
			}
			
			else {
				# feature column indices not defined yet
				
				# we are looking for the FEATURES line
				if ( $line =~ /^FEATURES/ ) {
					# we have the column header line
					# now need to identify the indices for specific columns
					
					# load the column headers into a quick temp hash
					my @data = split /\t/, $line;
					my %name_lookup;
					for (my $i = 0; $i < scalar @data; $i++) {
						$name_lookup{ $data[$i] } = $i;
					}
					
					# assign the index numbers
					$featurenum_i = $name_lookup{'FeatureNum'} or die 
						" unable to find 'FeautureNum' column!";
					$controltype_i = $name_lookup{'ControlType'} or die 
						" unable to find 'ControlType' column!";
					$probename_i = $name_lookup{'ProbeName'} or die 
						" unable to find 'ProbeName' column!";
					if ($processed) {
						# use processed signal
						$gsig_i = $name_lookup{'gProcessedSignal'} or die 
							" unable to find 'gProcessedSignal' column!";
						$rsig_i = $name_lookup{'rProcessedSignal'} or die 
							" unable to find 'rProcessedSignal' column!";
					}
					else{
						# use the raw signal, default
						$gsig_i = $name_lookup{'gMeanSignal'} or die 
							" unable to find 'gMeanSignal' column!";
						$rsig_i = $name_lookup{'rMeanSignal'} or die 
							" unable to find 'rMeanSignal' column!";
					}
					
					# move on to next line
					next FILELOOP;
				}
				
				else {
					# not the FEATURES line we're looking for
					next FILELOOP;
				}
			}
		}
		close INFILE;
		
		push @report, "Read file $myfile and collected $feature_count features\n";
	}
}



#### Perform normalization
sub perform_normalization {
	
	# Normalize experiment and control values separately
	if ($separate) {
		push @report, "\nQuantile normalizing experimental and control data sets separately\n";
		
		# prepare the data for normalization
		my @e_inputvalues; # an array of arrays for experiment
		my @c_inputvalues; # an array of arrays for control
		if (%cvalues) { # we have both experimental and control data
			my $i = 0;
			foreach (sort {$a <=> $b} keys %evalues) { # sort to keep the spots in order in the array
				push @e_inputvalues, [ @{ $evalues{$_} } ]; # push the experimental values first
				# the %cvalues hash should be sorted in the same manner
				push @c_inputvalues, [ @{ $cvalues{$_} } ]; # then push control values to the control array
				$i++;
			}
		} else { # we only have experimental data
			my $i = 0;
			foreach (sort {$a <=> $b} keys %evalues) { # sort to keep the spots in order in the array
				push @e_inputvalues, [ @{ $evalues{$_} } ]; # push the experimental values first
				$i++;
			}
		}
		
		# perform the normalization
		my @e_outputvalues = normalize(\@e_inputvalues);
		my @c_outputvalues;
		if (@c_inputvalues) { # only process the control array values if present
			@c_outputvalues = normalize(\@c_inputvalues);
		}
		
		# put the normalized values back
		if (%cvalues) { # we have both experimental and control data
			my $i = 0;
			foreach (sort {$a <=> $b} keys %evalues) {
				$evalues{$_} = [ @{ $e_outputvalues[$i] } ]; # 
				# the %cvalues hash should be sorted in the same manner
				$cvalues{$_} = [ @{ $c_outputvalues[$i] } ]; # 
				$i++;
			}
		} else { # we only have experimental data
			my $i = 0;
			foreach (sort {$a <=> $b} keys %evalues) {
				$evalues{$_} = [ @{ $e_outputvalues[$i] } ]; # 
				$i++;
			}
		}
		
		
		
	# Normalize experiment and control values together
	} else {
		push @report, "\nQuantile normalizing experimental and control data sets together\n";
		
		# prepare the data for normalization
		my @inputvalues; # an array of arrays
		my $i = 0;
		foreach (sort {$a <=> $b} keys %evalues) { # sort to keep the spots in order in the array
			my @data;
			push @data, @{ $evalues{$_} }; # push the experimental values first
			# the %cvalues hash should be sorted in the same manner
			push @data, @{ $cvalues{$_} }; # then push control
			$inputvalues[$i] = \@data; 
			$i++;
		}
		
		# perform the normalization
		my @outputvalues = normalize(\@inputvalues);
		
		# put the normalized values back
		$i = 0;
		my $e_end = scalar( @e_name_list ) - 1; # determine the array indices endpoints to correspond to experiment and control values
		my $c_start = $e_end + 1;
		my $c_end = scalar( @{ $outputvalues[0] } ) - 1;
		foreach (sort {$a <=> $b} keys %evalues) {
			$evalues{$_} = [ @{ $outputvalues[$i] }[0..$e_end] ]; # take the experiment values
			# the %cvalues hash should be sorted in the same manner
			$cvalues{$_} = [ @{ $outputvalues[$i] }[$c_start..$c_end] ]; # then take the control values 
			$i++;
		}
	}
}



### Quantile normalization subroutine
sub normalize {
	my $input = shift; # reference to input array
	# determine the number of replicate lists (columns) of data
	my $column_number = scalar @{ $input->[0] };
	# determine the number of data elements (rows) of data; 
	# this assumes all rows have the same number of columns - in theory it should
	my $row_number = scalar @$input;
	my %quantile_list; # the quantile list as a hash
	
	# Collect the data for each list
	for (my $i = 0; $i < $column_number; $i++) { # step through each replicate list one at time
		my @current_list; # put the values in each list into an array (don't need to remember order)
		for (0..($row_number-1)) {
			push @current_list, $input->[$_][$i];
		}
		
		# now need to sort the values in increasing order and put into quantile hash
		my $rank = 0;
		foreach (sort {$a <=> $b} @current_list) {
			# the quantile hash is keyed by the index number of the data elements (rows)
			# for each column (data set) we will be sorting the values in increasing order
			# and adding to a cumuluative sum for that rank index
			$quantile_list{$rank} += $_;
			$rank++;
		}
	}
	
	# Once data is collected, take the mean value of all the quantile values
	foreach (keys %quantile_list) {
		# we're dividing the cumulative sum for each rank index by the number of datasets
		# to generate the mean quantile value
		$quantile_list{$_} = $quantile_list{$_} / $column_number;
	}
	
	# Now replace initial values with the mean quantile value
	my @output;
	for (my $col = 0; $col < $column_number; $col++) { # for each data set
		my %startvalues; # temporary hash for the starting input value
		my %endvalues; # temporary hash of the ending value
		
		# Retrieve the data for the current list and put into starting values hash
		for my $row (0..($row_number-1)) {
			# recording the position index as key and original value as the hash value
			$startvalues{$row} = $input->[$row][$col];
		}
		
		# Determine the quantile value equivalent for the starting value
		my $index = 0;
		foreach (sort { $startvalues{$a} <=> $startvalues{$b} } keys %startvalues) {
			# sort by the hash value but return the hash key in $_
			# we are doing an increasing sort of the original starting value and 
			# getting its rank position, which may not necessarily be in order
			# then lookup the quantile value for this rank position
			$endvalues{$_} = $quantile_list{$index};
			# we are putting it in an output hash, the key is the original rank
			# position and the value is the quantile value
			# the quantile value is looked up by its own increasing rank position, $index
			$index++;
		}
		
		# Replace the starting values with the found quantile value
		for (0..($row_number-1)) {
			$output[$_][$col] = $endvalues{$_};
		}
				
	}
	
	return @output;
}



#### Determine the median scale from external list of normalization probes
sub get_normalization_probes {	
	# load the file
	open INFILE, $norm_list or die("unable to open median list file $norm_list!!!\n");
	my @mprobelist = <INFILE>; # median probe list array
	close INFILE;
	
	
	# collect the normalization probe values
	my @e_norm_values; # an array of the experimental values for each normalization probe
	my @c_norm_values; # an array of the control values for each normalization probe
	my $mprobenumber = 0;
	foreach my $myprobe (@mprobelist) {
		if ($myprobe =~ /^#/) {next} # skip comment lines;
		chomp $myprobe;
		foreach ( @{ $probes{$myprobe} } ) {
			# remember that each probe may be represented more than once on the array
			# we will be collecting feature values
			# calculate the mean value across all replicates for this probe feature
			my $value = &mean( @{ $evalues{$_} } );
			# then add the mean value to the list of normalization probes
			push @e_norm_values, $value;
			# repeat for the control values
			if (%cvalues) { # process control array only if it exists
				$value = &mean( @{ $cvalues{$_} } );
				push @c_norm_values, $value;
			}
		}
		$mprobenumber++; # count how many probes there are
	}
	
	# calculate the normalization value, which is the median of the normalization probe values
	$e_norm_target = &median(@e_norm_values); # global normalization value
	if (@c_norm_values) { # process control array only if it exists
		$c_norm_target = &median(@c_norm_values);
	}
	
	# report
	push @report, " Loaded $mprobenumber probes from file $norm_list to use in normalization\n";
	push @report, "   The experimental normalization value is $e_norm_target\n";
	if ($c_norm_target) {
		push @report, "   The control normalization value is $c_norm_target\n";
	}
}




#### Perform median scaling
sub perform_scaling_and_normalization {	
	# if there is an external normalization list, we will scale the median of
	# those probes to the desired median target specified by the user, and then
	# apply it to the entire data list
	
	# if there is no external normalization list, we will scale the median of 
	# the entire data list to the desired median target specified by the user, 
	# and then apply it to the entire data list
	
	if ($median_target) {
		push @report, "Median scaling data values to target of $median_target\n";
	}
	
	# determine the experimental scaling factor
	for (my $i = 0; $i < scalar(@e_name_list); $i++) {
		my $e_factor; # the scaling factor
		
		# both median scaling target and normalization
		if ($norm_list and $median_target) {
			$e_factor = $median_target / $e_norm_target;
			push @report, "  Median scaling and normalization factor for experiment $e_name_list[$i] is $e_factor\n";
		
		# normalization only
		} elsif ($norm_list and !$median_target) {
			$e_factor = 1 / $e_norm_target;
			push @report, "  Normalization factor for experiment $e_name_list[$i] is $e_factor\n";
			
		# median scaling only
		} elsif (!$norm_list and $median_target) {
			my @e_list;
			foreach my $n (keys %evalues) { # grab all values for a specific dataset
				push @e_list, $evalues{$n}[$i]; 
			}
			my $e_median = median(@e_list);
			$e_factor = $median_target / $e_median; # the scaling factor
			push @report, "  Median scaling factor for experiment $e_name_list[$i] is $e_factor\n";
		}
		
		# scale the experimental values
		foreach my $n (keys %evalues) {
			$evalues{$n}[$i] = $evalues{$n}[$i] * $e_factor;
		}
	}
	
	
	# determine the control scaling factor
	if (@c_name_list) {
		for (my $i = 0; $i < scalar(@c_name_list); $i++) {
			my $c_factor; # the scaling factor
			
			# both median scaling target and normalization
			if ($norm_list and $median_target) {
				$c_factor = $median_target / $c_norm_target;
				push @report, "  Median scaling and normalization factor for control $c_name_list[$i] is $c_factor\n";
			
			# normalization only
			} elsif ($norm_list and !$median_target) {
				$c_factor = 1 / $c_norm_target;
				push @report, "  Normalization factor for control $c_name_list[$i] is $c_factor\n";
				
			# median scaling only
			} elsif (!$norm_list and $median_target) {
				my @c_list;
				foreach my $n (keys %cvalues) { # grab all values for a specific dataset
					push @c_list, $cvalues{$n}[$i]; 
				}
				my $c_median = median(@c_list);
				$c_factor = $median_target / $c_median; # the scaling factor
				push @report, "  Median scaling factor for control $c_name_list[$i] is $c_factor\n";
			}
			
			# scale the control values
			foreach my $n (keys %cvalues) {
				$cvalues{$n}[$i] = $cvalues{$n}[$i] * $c_factor;
			}
		}
	}
}


#### Output combined two channel records
sub output_combined_records {
	# we will be generating an output data hash in the format for 
	# tim data text files, see tim_file_helper.pm
	# the main data hash is %output_data
	my @data_table; # an array of the processed microarray values
	
	# prepare the general metadata
	$output_data{'program'} = $0;
	$output_data{'feature'} = 'microarray_probes';
	$output_data{'gff'} = 0;
	$output_data{'number_columns'} = 4;
	
	# prepare the column metadata
	$output_data{0} = {
		'name'   => 'ProbeID',
		'index'  => 0,
	};
	$output_data{1} = {
		'name'   => 'Experiment',
		'index'  => 1,
		'source' => join(",", @e_name_list),
		'log2'   => $log,
	};
	$output_data{2} = {
		'name'   => 'Control',
		'index'  => 2,
		'source' => join(",", @c_name_list),
		'log2'   => $log,
	};
	$output_data{3} = {
		'name'   => 'Ratio',
		'index'  => 3,
		'log2'   => $log,
	};
	if ($norm_list) {
		# normalization list
		$output_data{1}{'normalization_list'} = $norm_list;
		$output_data{2}{'normalization_list'} = $norm_list;
	}
	if ($median_target) {
		# target for median scaling
		$output_data{1}{'median_target'} = $median_target;
		$output_data{2}{'median_target'} = $median_target;
	}
	if ($quant) {
		# quantile normalization performed
		if ($separate) {
			# each channel separately
			$output_data{1}{'quantile_normalize'} = 'separate';
			$output_data{2}{'quantile_normalize'} = 'separate';
		}
		else {
			# both channels together
			$output_data{1}{'quantile_normalize'} = 'together';
			$output_data{2}{'quantile_normalize'} = 'together';
		}
	}
	else {
		# no quantile normalization done
		$output_data{1}{'quantile_normalize'} = 'none';
		$output_data{2}{'quantile_normalize'} = 'none';
	}
	
	# prepare the column headers
	push @data_table, [ qw(ProbeID Experiment Control Ratio) ];
	
	# collect and combine the feature data
	foreach my $probe (sort {$a cmp $b} keys %probes) {
		my @e_list; # an array for all the feature values for a given probe
		my @c_list; # 
		foreach ( @{ $probes{$probe} } ) { # each probe may have more than one feature (spot)
			# collect the experiment & control values for all the features
			push @c_list, @{ $cvalues{$_} }; 
			push @e_list, @{ $evalues{$_} };
		}
		unless (@e_list) {die "no experimental feature values for probe $probe!\n"}
		unless (@c_list) {die "no control feature values for probe $probe!\n"}
		my $e_mean = &mean(@e_list); # calculate mean from all the features
		my $c_mean = &mean(@c_list);
		
		# copy the mean values into a global array for linear regression later
		push @{ $meanslist[0] }, $e_mean; # experiment is first sub-array
		push @{ $meanslist[1] }, $c_mean; # control is second sub-array
		
		# calculate ratio
		my $ratio;
		if ($c_mean == 0 or $c_mean eq '') {
			#print " probe $probe control is $c_mean with values @c_list\n";
			$ratio = '.';
		} 
		else {
			$ratio = $e_mean / $c_mean;
		}
		if ($log) { 
			# convert to log base 2 numbers
			$e_mean = ( log($e_mean) / log(2) );
			$c_mean = ( log($c_mean) / log(2) );
			unless ($ratio eq '.') {
				$ratio = ( log($ratio) / log(2) );
			}
		}
		push @data_table, [ ($probe, $e_mean, $c_mean, $ratio) ];
	}
	
	# set the data_table into the data structure
	$output_data{'data_table'} = \@data_table;
}



#### Output single channel records
sub output_single_records {
	# when there is no control data (%cvalues is empty)
	# we will be generating an output data hash in the format for 
	# tim data text files, see tim_file_helper.pm
	# the main data hash is %output_data
	my @data_table; # an array of the processed microarray values
	
	# prepare the general metadata
	$output_data{'program'} = $0;
	$output_data{'feature'} = 'microarray_probes';
	$output_data{'gff'} = 0;
	$output_data{'number_columns'} = 2;
	
	# prepare the column metadata
	$output_data{0} = {
		'name'   => 'ProbeID',
		'index'  => 0,
	};
	$output_data{1} = {
		'name'   => 'Experiment',
		'index'  => 1,
		'source' => join(",", @e_name_list),
		'log2'   => $log,
	};
	if ($norm_list) {
		# normalization list
		$output_data{1}{'normalization_list'} = $norm_list;
	}
	if ($median_target) {
		# target for median scaling
		$output_data{1}{'median_target'} = $median_target;
	}
	if ($quant) {
		# quantile normalization done
		$output_data{1}{'quantile_normalize'} = 'yes';
	}
	else {
		# no quantile normalization done
		$output_data{1}{'quantile_normalize'} = 'none';
	}
	
	
	
	# prepare the header file
	push @data_table, [ qw(ProbeID Experiment) ];
	
	# collect and combine the feature data
	foreach my $probe (sort {$a cmp $b} keys %probes) {
		my @e_list; # an array for all the feature values for a given probe
		foreach ( @{ $probes{$probe} } ) { # each probe may have more than one feature (spot)
			push @e_list, @{ $evalues{$_} }; # collect the experiment & control values for all the features
		}
		my $e_mean = &mean(@e_list); # calculate mean from all the features
		if ($log) { # convert to log base 2 numbers
			$e_mean = ( log($e_mean) / log(2) );
		}
		push @data_table, [ ($probe, $e_mean) ];
	}
	
	# set the data_table into the data structure
	$output_data{'data_table'} = \@data_table;
}



#### Perform regression on combined data means
sub perform_regression_on_means {
	my $stat = Statistics::Descriptive::Full->new(); # initialize
	# pass the reference to the experiment means list, this will by the y data
	$stat->add_data(\@{ $meanslist[0] });
	# the control means list will be the x data
	my ($q, $m, $r, $rms) = $stat->least_squares_fit(@{ $meanslist[1] });
	# intercept, slope, Pearson r, and root-mean-square are returned
	push @report, "\nLinear regression analysis:\n";
	push @report, " Experimental_Mean vs Control_mean\n";
	push @report, "   Intercept is " . sprintf("%.3f", $q) . ", slope is ". sprintf("%.3f", $m) . "\n";
	push @report, "   Pearson Correlation coefficient (r) is " . sprintf("%.3f", $r) . "\n";
	push @report, "   Coefficient of Determination (R^2) is ". sprintf("%.3f", $r ** 2) . "\n\n";
}



#### Perform regression on many replicates
sub perform_regression_on_replicates {
	
	# we will be walking systematically through the filelist and replicates
	# the first data set will be x, and will be compared to the subsequent data set, y
	# the y data set will then be incremented until all data sets have been compared 
	# the x data set will then be incremented by 1, and the process repeated
	
	# experimental data
	for (my $x = 0; $x < (scalar @e_name_list - 1); $x++) { # starting the x data set at 0
		for (my $y = 1; $y < scalar @e_name_list; $y++) { # starting the y data set at 1
			if ($y <= $x) {next}; # skip pairs we've already done
			
			# collect the x and y values
			my @e_xvalues; # temporay array to put the replicate values into for regression analysis 
			my @e_yvalues;
			foreach my $feature (keys %evalues) {
				push @e_xvalues, $evalues{$feature}[$x];
				push @e_yvalues, $evalues{$feature}[$y];
			}
			
			# calculate statistics for experimental values
			my $stat = Statistics::Descriptive::Full->new();
			$stat->add_data(@e_yvalues);
			my ($q, $m, $r, $rms) = $stat->least_squares_fit(@e_xvalues);
			# intercept, slope, Pearson r, and root-mean-square are returned
			
			# write out to report
			push @report, " Experimental $e_name_list[$x] vs $e_name_list[$y] values\n";
			push @report, "   Intercept is " . sprintf("%.3f", $q) . ", slope is ". sprintf("%.3f", $m) . "\n";
			push @report, "   Pearson Correlation coefficient (r) is " . sprintf("%.3f", $r) . "\n";
			push @report, "   Coefficient of Determination (R^2) is ". sprintf("%.3f", $r ** 2) . "\n\n";
		}
	}
	
	# control data
	if (@c_name_list) {
		for (my $x = 0; $x < (scalar @c_name_list - 1); $x++) { # starting the x data set at 0
			for (my $y = 1; $y < scalar @c_name_list; $y++) { # starting the y data set at 1
				if ($y <= $x) {next}; # skip pairs we've already done

				# collect the x and y values
				my @c_xvalues; # temporay array to put the replicate values into for regression analysis 
				my @c_yvalues;
				foreach my $feature (keys %cvalues) {
					push @c_xvalues, $cvalues{$feature}[$x];
					push @c_yvalues, $cvalues{$feature}[$y];
				}
				
				# calculate statistics for experimental values
				my $stat = Statistics::Descriptive::Full->new();
				$stat->add_data(@c_yvalues);
				my ($q, $m, $r, $rms) = $stat->least_squares_fit(@c_xvalues);
				# intercept, slope, Pearson r, and root-mean-square are returned
				
				# write out to report
				push @report, " Control $c_name_list[$x] vs $c_name_list[$y] values\n";
				push @report, "   Intercept is " . sprintf("%.3f", $q) . ", slope is ". sprintf("%.3f", $m) . "\n";
				push @report, "   Pearson Correlation coefficient (r) is " . sprintf("%.3f", $r) . "\n";
				push @report, "   Coefficient of Determination (R^2) is ". sprintf("%.3f", $r ** 2) . "\n\n";
			}
		}
	}
}



#### Generate a MA plot
sub generate_ma_plot {
	# M is the intensity ratio, defined log2R - log2G, plotted on y axis
	# A is the average intensity, defined as (log2R + log2G)/2, plotted on x axis
	
	# generate the data
	my $total = scalar @{ $meanslist[1] };
	my @data; # an array for the graphing data
	for my $i (0..($total - 1)) {
		# we will use the final experiment and control data to generate M and A
		# use @meanslist, experiment in subarray 0, control in subarray 1
		my $a = ( ( log($meanslist[0][$i]) / log(2) ) + ( log($meanslist[1][$i]) / log(2) ) ) / 2;
		my $m = ( log($meanslist[0][$i]) / log(2) ) - ( log($meanslist[1][$i]) / log(2) );
		push @{ $data[0] }, $a;
		push @{ $data[1] }, $m;
	}
	
	# prepare the graph
		my $graph = GD::Graph::points->new(800,600);
		$graph->set(
			x_label			=> 'average intensity',
			y_label			=> 'intensity ratio (exp - con)',
			title			=> "MA Plot for $outfile",
			markers			=> [1],
			marker_size		=> 1,
			x_tick_number	=> 4,
			y_tick_number	=> 4,
			x_number_format	=> "%.2f",
			y_number_format => "%.2f",
			x_label_position => 0.5,
			transparent		=> 0,
		) or warn $graph->error;
	
	# write the image file
	my $gd = $graph->plot(\@data) or warn $graph->error;
	my $filename = $outfile . '.png';
	open IMAGE, ">$filename" or warn " Can't open image output file!\n";
	binmode IMAGE;
	print IMAGE $gd->png;
	close IMAGE;
}

