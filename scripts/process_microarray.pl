#!/usr/bin/perl

# documentation at end of file
 
use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(mean median);
use Statistics::Descriptive;
use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure
	find_column_index
);
use Bio::ToolBox::file_helper qw(
	open_to_read_fh
	open_to_write_fh
	open_tim_data_file
	write_tim_data_file
);
use constant LOG2 => log(2);
my $VERSION = '1.18';

print "\n A script to process microarray files\n\n";

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
	$recipe,
	$outfile,
	$ch_ratio,
	$median_target,
	$norm_list,
	$separate,
	$quant,
	$log,
	$processed,
	$gz,
	$help,
	$print_version,
);
my @in_files;
my @channels;
GetOptions( 
	'in=s'       => \@in_files, # array of input files
	'recipe=s'   => \$recipe, # a recipe list of input files
	'out=s'      => \$outfile, # output file name
	'channel=s'  => \@channels, # channel designations
	'median=i'   => \$median_target, # median scale target value
	'nlist=s'    => \$norm_list, # filename for list of control probes to normalize
	'separate'   => \$separate, # quantile normalize experiment and control separately
	'norm!'      => \$quant, # do not quantile normalize
	'log!'       => \$log, # convert to log2
	'processed'  => \$processed, # use processed values rather than raw
	'gz!'        => \$gz, # compress output
	'help'       => \$help, # print help
	'version'    => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script process_microarray.pl, version $VERSION\n\n";
	exit;
}




### Check for requirements and assign default values
unless (@in_files or $recipe) {
	die " No input file name(s) or recipe file specified! use --help\n";
}
if (scalar @in_files == 1 and $in_files[0] =~ /,/) {
	# comma delimited list of infiles
	my $list = shift @in_files;
	@in_files = split /,/, $list;
}
if (scalar @in_files > 0 and scalar @channels == 0) {
	die " No channel designations specified! use --help\n";
}
if (scalar @channels == 1 and $channels[0] =~ /,/) {
	# comma delimited list of infiles
	my $list = shift @channels;
	@channels = split /,/, $list;
}
if (scalar @in_files != scalar @channels) {
	die " Not an equal number of input files and channel designations!\n";
}

unless ($outfile) {
	if ($recipe) {
		$outfile = $recipe;
		$outfile =~ s/\.txt$//;
		$outfile .= '_output';
	}
	else {
		die " No output file name specified! use --help\n";
	}
}
unless (defined $quant) {
	$quant = 1;
}
unless (defined $log) {
	$log = 1;
}




### Define global data containing variables
my %evalues; 
	# a hash of arrays to keep the experimental values: 
	# key        = unique_id, 
	# value      = array of red values
my %cvalues; 
	# a hash of arrays to keep the control values: 
	# key        = unique_id, 
	# value      = array of green values
my %probes; 
	# a hash of arrays to keep the probe names and corresponding 
	# feature numbers: 
	# key        = probe_id, 
	# second key = feature number
	# value      = unique_id
my $unique_id = 1; 
	# unique_id for the probe IDs
my %data_sets;
	# a hash to remember the data sources and
	# key       = experiment or control
	# value     = array of data source filenames and channels
my @report; 
	# output array for reporting results and analyses






### Load the experimental and control 
	# the list of input files may be provided directly by command line arguments
	# or provided in a recipe list in a file specified by a command line argument
load_microarray_data();
unless (%evalues) {
	die " Nothing loaded for experiment data!\n";
}

# check that all probes loaded have equal counts of values
check_microarray_value_counts();




### Perform quantile normalization on the data
if ($quant) {
	print " Quantile normalizing values....\n";
	perform_quantile_normalization();
}




### Perform median scaling and probe normalization
# if a list of control probes to be normalized against was provided, 
# then it will be normalized here
if ($norm_list) {
	print " Normalizing to list of control probes....\n";
}
if ($median_target) {
	print " Median scaling values to $median_target....\n";
}
perform_probe_normalization_median_scaling();




### Preparing the data sets for output
my $output;
if (%cvalues) { 
	# there is both experimental and control data values
	print " Combining replicate data sets for both experimental and control data....\n";
	push @report, "\nCombining the replicate data sets for both experimental and control data\n";
	$output = output_combined_records();
} 
else { 
	# there is only experimental data values (single channel)
	print " Combining replicate data sets for experimental data....\n";
	push @report, "\nCombining the replicate data sets for experimental data\n";
	$output = output_single_records();
}

if ($log) {
	push @report, "Output values converted to log base 2\n";
	print " Output values converted to log base 2\n";
	# this is actually done when preparing the data sets for output above
}




### Linear regression analysis
print " Performing linear regression analysis....\n";
push @report, "\nLinear regression analysis:\n";
push @report, join("\t", qw(
	Group
	DataSetX
	DataSetY
	Intercept
	Slope
	PearsonCorrelation(r)
	CoefficientDetermination(r^2)
) ) . "\n";

# calculate statistics on the combined experimental & control values
if (%cvalues) {
	perform_regression_on_experiment_control($output);
}

# calculate statistics pairwise among replicates
perform_regression_on_replicates('experiment', \%evalues);
if (%cvalues) {
	perform_regression_on_replicates('control', \%cvalues);
}




### Write output
$outfile =~ s/\.txt$//; # strip extension if present
# write main output
my $written_file = write_tim_data_file(
	# we will write a tim data file
	'data'     => $output,
	'filename' => $outfile,
	'gz'       => $gz,
);
if ($written_file) {
	print " Wrote data file '$written_file'\n";
	push @report, "\nWrote output file to $written_file\n\n\n";
}
else {
	print " unable to write data file!\n";
}

# write the report file
my $report_fh = open_to_write_fh("$outfile\_report.txt");
$report_fh->print(@report);
$report_fh->close;
print "Wrote report file to $outfile\_report.txt\n";

print "Finished!\n";



############################   Subroutines   ###################################

#### Load the microarray data files
sub load_microarray_data {
	
	# load the contents of the recipe file into the in_files and channels 
	# desigination arrays
	if ($recipe) {
		# this will load a processing recipe file
		# a space delimited text file giving the names of the input files 
		# and the the designation usage for the red and green channels
		# one row per file
		# each line consists of file_name, red_channel designation, 
		# green_channel designation
		# each designation must be one of three values: (e)xperiment, (c)ontrol, 
		# or (n)one
		
		my $fh = open_to_read_fh($recipe) or 
			die " unable to open recipe file '$recipe'!\n";
		
		# process the file
		while (my $line = $fh->getline) {
			next if ($line =~ /^#/); # skip any comment lines
			next unless ($line =~ /\w/i); # skip blank lines
			chomp $line;
			my ($file, $channel) = split /\s+/, $line, 2;
			if ($file and $channel) {
				# defined file and channel designation
				# push into appropriate arrays
				push @in_files, $file;
				push @channels, $channel;
			}
			else {
				# bad structure
				warn "bad recipe file line format!\n";
			}
		}
		$fh->close;
	}
	
	# load the arrays by guessing the appropriate type
	while (@in_files) {
		my $file    = shift @in_files;
		my $channel = shift @channels;
		
		if ($file =~ /\.txt(?:\.gz)?$/i and $channel =~ /^[ecn]{2}$/) {
			# a .txt file and two channel designation
			# we presume an Agilent file
			load_agilent_file($file, $channel);
		}
		elsif ($file =~ /\.gpr(?:\.gz)?$/i and $channel =~ /^[ecn]{2}$/) {
			# a .gpx file and two channel designation
			# we presume a GenePix file
			load_genepix_file($file, $channel);
		}
		elsif ($file =~ /\.pair(?:\.gz)?$/i and $channel =~ /^[ecn]{1}$/) {
			# a .pair file and one channel designation
			# we presume a NimbleGen file
			load_nimblegen_file($file, $channel);
		}
		else {
			warn "unrecognized file extension or channel designation!\n" . 
				"  file '$file' and channel '$channel'\n";
		}
		
	}
	
	# summary
	push @report, "Loaded " . scalar(keys %probes) . " oligo probes" . 
		" represented by " . scalar(keys %evalues) . " microarray features\n";
}




#### Load Agilent Files
sub load_agilent_file {
	my ($file, $channel) = @_;
	
	# open the file
	print " Reading file $file....\n";
	my $fh = open_to_read_fh($file) or 
		warn " unable to open file!\n";
	return unless $fh;
	
	# column indices values
	my ($featurenum_i, $controltype_i, $probename_i, $gsig_i, $rsig_i); 
	
	# assign the separate channel designations
	my ($red, $green) = split //, $channel, 2;
	
	# walk through the file
	my $feature_count = 0;
	while (my $line = $fh->getline) {
		
		# Agilent files have multiple header lines and a relatively lax format
		# Plus the columns have changed over different software versions
		# Because of the complex header format, we can't use my standard 
		# find_column_index() method, so we'll do things manually <sigh!>
		
		# process header lines
		if (!defined $probename_i) {
			# feature column indices not defined yet
			
			# we are looking for the FEATURES line
			if ( $line =~ /^FEATURES/ ) {
				# we finally have the column header line
				# now need to identify the indices for specific columns
				
				# load the column headers into a quick temp hash
				my @data = split /\t/, $line;
				my %name_lookup;
				for (my $i = 0; $i < scalar @data; $i++) {
					$name_lookup{ $data[$i] } = $i;
				}
				
				# assign the index numbers
				$controltype_i = $name_lookup{'ControlType'} or die 
					" unable to find 'ControlType' column!";
				$featurenum_i = $name_lookup{'FeatureNum'} or die 
					" unable to find 'FeautureNum' column!";
				$probename_i = $name_lookup{'ProbeName'} or die 
					" unable to find 'ProbeName' column!";
				if ($processed) {
					# use processed signal
					$gsig_i = $name_lookup{'gProcessedSignal'} or die 
						" unable to find 'gProcessedSignal' column!";
					$rsig_i = $name_lookup{'rProcessedSignal'} or die 
						" unable to find 'rProcessedSignal' column!";
				}
				else {
					# use the raw signal, default
					$gsig_i = $name_lookup{'gMeanSignal'} or die 
						" unable to find 'gMeanSignal' column!";
					$rsig_i = $name_lookup{'rMeanSignal'} or die 
						" unable to find 'rMeanSignal' column!";
				}
			}
			
			else {
				# not the FEATURES line we're looking for
				next;
			}
		}
		
		# processing data lines
		else {
			# column indices have been defined
			# now just working with data lines
			my @data   = split /\t/, $line;
			my $probe  = $data[$probename_i];
			my $number = $data[$featurenum_i];
			
			# only take the experimental probe values
			next if ($data[$controltype_i] != 0);
				
			# Need to identify the lookup unique_id for this probe
			my $lookup;
			if (exists $probes{$probe} ) {
				# we've seen this probe before
				
				# check if it is the same microarray feature 
				if (exists $probes{$probe}{$number}) {
					# same feature, this value likely coming from a new dataset
					$lookup = $probes{$probe}{$number};
				}
				else {
					# new microarray feature for the same probe
					# add new id
					$probes{$probe}{$number} = $unique_id;
					$lookup = $unique_id;
					$unique_id++; # prepare for the next one
				}
			}
			else {
				# new microarray feature for this probe
				# add new id
				$probes{$probe}{$number} = $unique_id;
				$lookup = $unique_id;
				$unique_id++; # prepare for the next one
			}
			
			# collect red channel values
			if ($red =~ /^e/i) { 
				push @{ $evalues{$lookup} },  $data[$rsig_i];
			} 
			elsif ($red =~ /^c/i) {
				push @{ $cvalues{$lookup} },  $data[$rsig_i];
			}
			
			# collect green channel values
			if ($green =~ /^e/i) { 
				push @{ $evalues{$lookup} },  $data[$gsig_i];
			} 
			elsif ($green =~ /^c/i) {
				push @{ $cvalues{$lookup} },  $data[$gsig_i];
			}
			
			# increment counter
			$feature_count++;
		}
		
	}
	
	# done
	$fh->close;
	
	# print reports and remember the dataset names for later
	if ($feature_count) {
		push @report, "Read file $file\n Collected $feature_count features\n";
		
		# red channel reports
		if ($red eq 'e') {
			push @report, "  Cy5 red channel collected into experiment data array\n";
			push @{ $data_sets{'experiment'} }, "$file:red";
		}
		elsif ($red eq 'c') {
			push @report, "  Cy5 red channel collected into control data array\n";
			push @{ $data_sets{'control'} }, "$file:red";
		}
		elsif ($red eq 'n') {
			push @report, "  Cy5 red channel not collected\n";
		}
		
		# green channel reports
		if ($green eq 'e') {
			push @report, "  Cy3 green channel collected into experiment data array\n";
			push @{ $data_sets{'experiment'} }, "$file:green";
		}
		elsif ($green eq 'c') {
			push @report, "  Cy3 green channel collected into control data array\n";
			push @{ $data_sets{'control'} }, "$file:green";
		}
		elsif ($green eq 'n') {
			push @report, "  Cy3 green channel not collected\n";
		}
	}
	else {
		warn " Failed to import features from file!\n";
	}
}



#### Load GenePix files
sub load_genepix_file {
	my ($file, $channel) = @_;
	
	# open the file
	print " Reading file $file....\n";
	my $fh = open_to_read_fh($file) or 
		warn " unable to open file!\n";
	return unless $fh;
	
	# column indices values
	my ($probename_i, $gsig_i, $rsig_i); 
		# GenePix files don't have a column of unique numbers
		# like Agilent and Nimblegen, just the probe_id
	
	# assign the separate channel designations
	my ($red, $green) = split //, $channel, 2;
	
	# walk through the file
	my $feature_count = 0;
	my $row = 1;
	while (my $line = $fh->getline) {
		
		# GenePix files have multiple header lines and a relatively lax format
		# Plus the columns have changed over different software versions
		# Because of the complex header format, we can't use my standard 
		# find_column_index() method, so we'll do things manually <sigh!>
		
		# process header lines
		if (!defined $probename_i) {
			# feature column indices not defined yet
			
			# we are looking for the FEATURES line
			if ( $line =~ /column.+row.+name.+id/i ) {
				# we finally have the column header line
				# assuming this hasn't changed
				# now need to identify the indices for specific columns
				
				# load the column headers into a quick temp hash
				my @data = split /\t/, $line;
				my %name_lookup;
				for (my $i = 0; $i < scalar @data; $i++) {
					my $value = $data[$i];
					$value =~ s/"//g; # strip quotations that may be present
					$name_lookup{$value} = $i;
				}
				
				# assign the index numbers
				$probename_i = $name_lookup{'ID'} or 
					die " unable to find 'ID' column!";
				if ($processed) {
					# use processed signal
					# this feature mean - background mean
					$gsig_i = $name_lookup{'F532 Mean - B532'} or die 
						" unable to find 'F532 Mean - B532' column!";
					$rsig_i = $name_lookup{'F635 Mean - B635'} or die 
						" unable to find 'F635 Mean - B635' column!";
				}
				else {
					# use the raw signal, default
					$gsig_i = $name_lookup{'F532 Mean'} or die 
						" unable to find 'F532 Mean' column!";
					$rsig_i = $name_lookup{'F635 Mean'} or die 
						" unable to find 'F635 Mean' column!";
				}
				
				$row++;
			}
			
			else {
				# not the FEATURES line we're looking for
				$row++;
				next;
			}
		}
		
		# processing data lines
		else {
			# column indices have been defined
			# now just working with data lines
			my @data   = split /\t/, $line;
			my $probe  = $data[$probename_i];
			
			# generate a feature number
			# this is simply the row number, since the gpx file format doesn't 
			# have unique feature numbers, ID is good enough for them but 
			# not for my program.... :-P
			my $number = $row;
			
			# Need to identify the lookup unique_id for this probe
			my $lookup;
			if (exists $probes{$probe} ) {
				# we've seen this probe before
				
				# check if it is the same microarray feature 
				if (exists $probes{$probe}{$number}) {
					# same feature, this value likely coming from a new dataset
					$lookup = $probes{$probe}{$number};
				}
				else {
					# new microarray feature for the same probe
					# add new id
					$probes{$probe}{$number} = $unique_id;
					$lookup = $unique_id;
					$unique_id++; # prepare for the next one
				}
			}
			else {
				# new microarray feature for this probe
				# add new id
				$probes{$probe}{$number} = $unique_id;
				$lookup = $unique_id;
				$unique_id++; # prepare for the next one
			}
			
			# collect red channel values
			if ($red =~ /^e/i) { 
				push @{ $evalues{$lookup} },  $data[$rsig_i];
			} 
			elsif ($red =~ /^c/i) {
				push @{ $cvalues{$lookup} },  $data[$rsig_i];
			}
			
			# collect green channel values
			if ($green =~ /^e/i) { 
				push @{ $evalues{$lookup} },  $data[$gsig_i];
			} 
			elsif ($green =~ /^c/i) {
				push @{ $cvalues{$lookup} },  $data[$gsig_i];
			}
			
			# increment counter
			$row++;
			$feature_count++;
		}
		
	}
	
	# done
	$fh->close;
	
	# print reports and remember the dataset names for later
	if ($feature_count) {
		push @report, "Read file $file\n Collected $feature_count features\n";
		
		# red channel reports
		if ($red eq 'e') {
			push @report, "  Cy5 red channel collected into experiment data array\n";
			push @{ $data_sets{'experiment'} }, "$file:red";
		}
		elsif ($red eq 'c') {
			push @report, "  Cy5 red channel collected into control data array\n";
			push @{ $data_sets{'control'} }, "$file:red";
		}
		elsif ($red eq 'n') {
			push @report, "  Cy5 red channel not collected\n";
		}
		
		# green channel reports
		if ($green eq 'e') {
			push @report, "  Cy3 green channel collected into experiment data array\n";
			push @{ $data_sets{'experiment'} }, "$file:green";
		}
		elsif ($green eq 'c') {
			push @report, "  Cy3 green channel collected into control data array\n";
			push @{ $data_sets{'control'} }, "$file:green";
		}
		elsif ($green eq 'n') {
			push @report, "  Cy3 green channel not collected\n";
		}
	}
	else {
		warn " Failed to import features from file!\n";
	}
}



#### Load NimbleGen files
sub load_nimblegen_file {
	my ($file, $channel) = @_;
	
	# open the file
	print " Reading file $file....\n";
	my ($fh, $metadata) = open_tim_data_file($file) or 
		warn " unable to open file!\n";
	return unless $fh;
	
	# determine column indices
	my $featurenum_i = find_column_index($metadata, '^MATCH_INDEX') or 
		die " unable to find MATCH_INDEX column index! Is this a Nimblegen .pair file?\n";
	my $probename_i = find_column_index($metadata, '^PROBE_ID') or 
		die " unable to find PROBE_ID column index! Is this a Nimblegen .pair file?\n";
	my $signal_i = find_column_index($metadata, '^PM') or 
		die " unable to find PM column index! Is this a Nimblegen .pair file?\n";
	
	# walk through the file
	my $feature_count = 0;
	while (my $line = $fh->getline) {
		
		# get data
		chomp $line;
		my @data   = split /\t/, $line;
		my $probe  = $data[$probename_i];
		my $number = $data[$featurenum_i];
		
		# Need to identify the lookup unique_id for this probe
		# I don't think we have duplicate probes here on Nimblegen slides
		# but, just in case, we'll use the same strategy as with Agilent files
		# and use both PROBE_ID and MATCH_INDEX (a unique number) to uniquely 
		# identify probes, especially across multi-slide sets
		my $lookup;
		if (exists $probes{$probe} ) {
			# we've seen this probe before
			
			# check if it is the same microarray feature 
			if (exists $probes{$probe}{$number}) {
				# same feature, this value likely coming from a new dataset
				$lookup = $probes{$probe}{$number};
			}
			else {
				# new microarray feature for the same probe
				# add new id
				$probes{$probe}{$number} = $unique_id;
				$lookup = $unique_id;
				$unique_id++; # prepare for the next one
			}
		}
		else {
			# new microarray feature for this probe
			# add new id
			$probes{$probe}{$number} = $unique_id;
			$lookup = $unique_id;
			$unique_id++; # prepare for the next one
		}
		
		# collect the signal value
		# it looks like Nimblegen slide allows for both Perfect Match (PM) and 
		# MisMatch (MM) oligos. There is a column for both. We're only 
		# interested in the PM values
		if ($channel =~ /^e/i) { 
			push @{ $evalues{$lookup} },  $data[$signal_i];
		} 
		elsif ($channel =~ /^c/i) {
			push @{ $cvalues{$lookup} },  $data[$signal_i];
		}
		elsif ($channel =~ /^n/i) {
			warn " why bother loading a file that won't be used!!??";
			return;
		}
		# increment counter
		$feature_count++;
	}	
	
	# done
	$fh->close;
	
	# print reports and remember the dataset names for later
	if ($feature_count) {
		push @report, "Read file $file\n Collected $feature_count features\n";
		if ($channel eq 'e') {
			push @report, "  PM values collected into experiment data array\n";
			push @{ $data_sets{'experiment'} }, "$file";
		}
		elsif ($channel eq 'c') {
			push @report, "  PM values collected into control data array\n";
			push @{ $data_sets{'control'} }, "$file";
		}
	}
	else {
		warn " Failed to import features from file!\n";
	}
}



### check the number of values for each probe
sub check_microarray_value_counts {
	
	# need to verify the number of values for each microarray feature loaded
	# they should all the be same, or the quantile normalization won't work
	
	# check the experiment values
	my %evalue_counts;
	foreach my $id (keys %evalues) {
		my $count = scalar @{ $evalues{$id} };
		$evalue_counts{$count} += 1;
	}
	
	# evaluate counts
	if (scalar( keys %evalue_counts ) > 1) {
		warn " There is a problem with the number of values loaded for each\n" .
			" microarray feature in the experiment group! There should be\n" . 
			" an equal number for each! These are the counts I found:\n";
		foreach my $count (sort {$a <=> $b} keys %evalue_counts) {
			warn "   $count\t$evalue_counts{$count}\n";
		}
		warn " I am unable to continue with processing!\n";
		
		# dump the values
		warn " writing out the experiment values for your evaluation...\n";
		my $dump = open_to_write_fh("$outfile\_experiment_values.txt");
		$dump->print(join("\t", qw(Probe FeatureNumber Values...)) . "\n");
		foreach my $probe (sort {$a cmp $b} keys %probes) {
			foreach my $number (sort {$a <=> $b} keys %{ $probes{$probe} } ) {
				$dump->print( 
					join("\t", 
						$probe,
						$number,
						@{ $evalues{ $probes{$probe}{$number} } },
					) 
				. "\n");
			}
		}
		warn " wrote dump file $outfile\_experiment_values.txt\n";
		exit;
	}	
	
	# check the control values if necessary
	if (%cvalues) {
		my %cvalue_counts;
		
		# count them
		foreach my $id (keys %cvalues) {
			my $count = scalar @{ $cvalues{$id} };
			$cvalue_counts{$count} += 1;
		}
		
		# evaluate
		if (scalar( keys %cvalue_counts ) > 1) {
			warn " There is a problem with the number of values loaded for each\n" .
				" microarray feature in the control group! There should be\n" . 
				" an equal number for each! These are the counts I found:\n";
			foreach my $count (sort {$a <=> $b} keys %cvalue_counts) {
				warn "   $count\t$cvalue_counts{$count}\n";
			}
			warn " I am unable to continue with processing!\n";
			
			# dump the values
			warn " writing out the control values for your evaluation...\n";
			my $dump = open_to_write_fh("$outfile\_control_values.txt");
			$dump->print(join("\t", qw(Probe FeatureNumber Values...)) . "\n");
			foreach my $probe (sort {$a cmp $b} keys %probes) {
				foreach my $number (sort {$a <=> $b} keys %{ $probes{$probe} } ) {
					$dump->print( 
						join("\t", 
							$probe,
							$number,
							@{ $cvalues{ $probes{$probe}{$number} } },
						) 
					. "\n");
				}
			}
			warn " wrote dump file $outfile\_control_values.txt\n";
			exit;
		}	
	}
	
}



#### Perform normalization
sub perform_quantile_normalization {
	
	# Normalize experiment and control values separately
	if ($separate) {
		push @report, "\nQuantile normalizing experimental and control data sets separately\n";
		
		# prepare the data for normalization
		
		# an array of arrays for experiment and control values
		my @e_inputvalues; 
		my @c_inputvalues; 
		
		# we have both experimental and control data
		if (%cvalues) { 
			
			foreach my $id (sort {$a <=> $b} keys %evalues) { 
				# sort to keep the spots in order in the array
				
				# push the experimental values first
				push @e_inputvalues, [ @{ $evalues{$id} } ]; 
				
				# then push control values to the control array
				# the %cvalues hash would be sorted in the same manner
				if (exists $cvalues{$id}) {
					push @c_inputvalues, [ @{ $cvalues{$id} } ]; 
				}
			}
		} 
		
		# we only have experimental data
		else { 
			foreach my $id (sort {$a <=> $b} keys %evalues) { 
				# sort to keep the spots in order in the array
				
				# push the experimental values first
				push @e_inputvalues, [ @{ $evalues{$id} } ]; 
			}
		}
		
		# perform the experiment normalization
		my (@e_outputvalues, @c_outputvalues);
		if (scalar @{ $e_inputvalues[0] } > 1) {
			# only normalize if there is more than one dataset
			@e_outputvalues = normalize(\@e_inputvalues);
		}
		else {
			push @report, " quantile normalization for experiment not performed, only one dataset\n";
		}
		
		# perform the control normalization
		if (scalar @{ $c_inputvalues[0] } > 1) {
			# only normalize if there is more than one dataset
			@c_outputvalues = normalize(\@c_inputvalues);
		}
		else {
			push @report, " quantile normalization for control not performed, only one dataset\n";
		}
		
		# put the normalized values back
		if (@e_outputvalues) { 
			# we have normalized experimental data
			my $i = 0;
			foreach my $id (sort {$a <=> $b} keys %evalues) {
				# we're sorting in the same manner, so we can take the 
				# elements of the output array sequentially since it will 
				# be the same
				$evalues{$id} = [ @{ $e_outputvalues[$i] } ]; # 
				$i++;
			}
		}
		if (@c_outputvalues) { 
			# we have normalized control data
			my $i = 0;
			foreach my $id (sort {$a <=> $b} keys %cvalues) {
				# we're sorting in the same manner, so we can take the 
				# elements of the output array sequentially since it will 
				# be the same
				$cvalues{$id} = [ @{ $c_outputvalues[$i] } ]; # 
				$i++;
			}
		} 
	} 
	
	# Normalize experiment and control values together
	else {
		push @report, "\nQuantile normalizing experimental and control data sets together\n";
		
		# prepare the data for normalization
		my @inputvalues; # an array of arrays
		my $experiment_count = 0; # number of experiment data sets
		my $control_count    = 0; # number of control data sets
		foreach my $id (sort {$a <=> $b} keys %evalues) { 
			# sort to keep the spots in order in the array
			my @data; # temporary array to put both experimental and 
					  # control values, to become the secondary data array in 
					  # @inputvalues
			
			# push the experimental values first
			push @data, (@{ $evalues{$id} }); 
			
			# then push control
			# the %cvalues hash keys would be sorted in the same manner
			if (exists $cvalues{$id}) {
				push @data, (@{ $cvalues{$id} });
			}
			
			# push to the input set
			push @inputvalues, \@data; 
			
			# record the number of data sets in experiment and control
			unless ($experiment_count) {
				# this assumes that we have the same data set count for all
				# probes - why wouldn't we?
				$experiment_count = scalar @{ $evalues{$id} };
				$control_count    = scalar @{ $cvalues{$id} }; 
			}
		}
		
		# perform the normalization, if we have more than one dataset
		unless (scalar @{ $inputvalues[0] } > 1) {
			push @report, " Quantile normalization not necessary, only one dataset!\n";
			return;
		}
		my @outputvalues = normalize(\@inputvalues);
		
		# determine the array indices endpoints to correspond to 
		# experiment and control values
		my $e_end   = $experiment_count - 1; 
		my $c_start = $experiment_count;
		my $c_end   = $experiment_count + $control_count - 1;
		
		# put the normalized values back
		my $i = 0;
		foreach (sort {$a <=> $b} keys %evalues) {
			# we're sorting in the same manner, so we can take the 
			# elements of the output array sequentially since it will 
			# be the same
			
			# take the experiment values
			$evalues{$_} = [ @{ $outputvalues[$i] }[0..$e_end] ]; 
			
			# then take the control values 
			# the %cvalues hash should be sorted in the same manner
			$cvalues{$_} = [ @{ $outputvalues[$i] }[$c_start..$c_end] ]; 
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
	for (my $i = 0; $i < $column_number; $i++) { 
		# step through each replicate list one at time
		my @current_list; 
		
		# put the values in each list into an array (don't need to remember order)
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
	for (my $col = 0; $col < $column_number; $col++) { 
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
	
	# open the file
	my $fh = open_to_read_fh($norm_list) or 
		die("unable to open normalization probe list file '$norm_list'!!!\n");
	
	# initialize variables
	my @e_norm_values; 
	my @c_norm_values; 
	my $norm_probe_count = 0;
	my ($e_norm_target, $c_norm_target);
	
	# collect the normalization probe values
	while (my $line = $fh->getline) {
		next if $line =~ /^#/; # skip comment lines;
		next unless $line =~ /\w+/i; # skip empty lines
		chomp $line;
		# each line should consist of a probe name
		# if not, we simply won't find it in the probe hash
		
		if (exists $probes{$line}) {
			# remember that each probe may be represented by more than one  
			# oligo feature on the the array
			# and then we may have multiple data sets (microarray 
			# hybridizations) for each slide
			# we will calculate the mean value 
			# first across all data sets for a given microarray feature,
			# then across all the microarray features for a given probe
			
			# take the mean of each data set value for each microarray 
			# probe feature
			my @feature_values = map { 
				mean( @{ $evalues{$_} } )  
			} values %{ $probes{$line} };
			
			# then take the mean of all of the duplicate microarray probe 
			# feature positions
			push @e_norm_values, mean(@feature_values);
			
			# repeat for the control values as necessary
			if (%cvalues) { 
				# take the mean of each data set value for each microarray 
				# probe feature
				my @feature_values = map { 
						mean( @{ $cvalues{$_} } )  
					} values %{ $probes{$line} }; 
				
				# then take the mean of all of the replicate microarray probe 
				# positions
				push @c_norm_values, mean(@feature_values);
			}
			
			# count how many probes there are
			$norm_probe_count++; 
		}
	}
	$fh->close;
	
	# calculate the normalization value, 
	# which is the median of the all the normalization probe values
	$e_norm_target = median(@e_norm_values); 
	if (@c_norm_values) { 
		# process control array only if it exists
		$c_norm_target = median(@c_norm_values);
	}
	
	# report
	push @report, " Loaded $norm_probe_count probes from file '$norm_list' to use in normalization\n";
	push @report, "   The experimental normalization value is $e_norm_target\n";
	if ($c_norm_target) {
		push @report, "   The control normalization value is $c_norm_target\n";
	}
	
	# return
	return ($e_norm_target, $c_norm_target);
}




#### Perform median scaling
sub perform_probe_normalization_median_scaling {	
	# if there is an external normalization list, we will scale the median of
	# those probes to the desired median target specified by the user, and then
	# apply it to the entire data list
	
	# if there is no external normalization list, we will scale the median of 
	# the entire data list to the desired median target specified by the user, 
	# and then apply it to the entire data list
	
	# Get list of probes to which the array will be normalized
	# and calculate the normalization factor to apply
	my ($e_norm_target, $c_norm_target);
	if ($norm_list) {
		($e_norm_target, $c_norm_target) = get_normalization_probes();
	}
	
	if ($median_target) {
		push @report, " Median scaling data values to target of $median_target\n";
	}
	
	
	# determine the experimental scaling factor
	for my $data_index (0..(scalar @{ $evalues{1} } - 1) ) {
		# work through each experimental data set
		# calculate a different factor for each data set, based on whether
		# we are median scaling, probe normalizing, or both
		# and apply that 
		
		# the scaling factor
		my $e_factor; 
		
		# both median scaling target and normalization
		if ($norm_list and $median_target) {
			$e_factor = $median_target / $e_norm_target;
			push @report, "  Median scaling and normalization factor for" . 
				" experiment data index $data_index is $e_factor\n";
		}
		
		# normalization only
		elsif ($norm_list and !$median_target) {
			$e_factor = 1 / $e_norm_target;
			push @report, "  Normalization factor for experiment data index" . 
				" $data_index is $e_factor\n";
		}
		
		# median scaling only
		elsif (!$norm_list and $median_target) {
			
			# grab all values for a specific dataset
			my @list;
			foreach my $n (keys %evalues) { 
				push @list, $evalues{$n}->[$data_index]; 
			}
			
			# calculate that the median for that dataset
			my $e_median = median(@list);
			
			# calculate the scaling factor
			$e_factor = $median_target / $e_median; 
			push @report, "  Median scaling factor for experiment data index" . 
				" $data_index is $e_factor\n";
		}
		
		# scale the experimental values if necessary
		if ($e_factor) {
			foreach my $probe (keys %evalues) {
				$evalues{$probe}->[$data_index] = 
					$evalues{$probe}->[$data_index] * $e_factor;
			}
		}
	}
	
	
	# determine the control scaling factor
	if (%cvalues) {
		for my $data_index (0..(scalar @{ $cvalues{1} } - 1) ) {
			# work through each controlal data set
			# calculate a different factor for each data set, based on whether
			# we are median scaling, probe normalizing, or both
			# and apply that 
			
			# the scaling factor
			my $c_factor; 
			
			# both median scaling target and normalization
			if ($norm_list and $median_target) {
				$c_factor = $median_target / $c_norm_target;
				push @report, "  Median scaling and normalization factor for" .
					" control data set index $data_index is $c_factor\n";
			}
			
			# normalization only
			elsif ($norm_list and !$median_target) {
				$c_factor = 1 / $c_norm_target;
				push @report, "  Normalization factor for control data index" .
					" $data_index is $c_factor\n";
			}
			
			# median scaling only
			elsif (!$norm_list and $median_target) {
				
				# grab all values for a specific dataset
				my @list;
				foreach my $n (keys %cvalues) { 
					push @list, $cvalues{$n}->[$data_index]; 
				}
				
				# calculate that the median for that dataset
				my $e_median = median(@list);
				
				# calculate the scaling factor
				$c_factor = $median_target / $e_median; 
				push @report, "  Median scaling factor for control data index" .
					" $data_index is $c_factor\n";
			}
			
			# scale the control values if necessary
			if ($c_factor) {
				foreach my $probe (keys %cvalues) {
					$cvalues{$probe}->[$data_index] = 
						$cvalues{$probe}->[$data_index] * $c_factor;
				}
			}
		}
	}
	# done with normalizing
}


#### Output combined two channel records
sub output_combined_records {
	
	# prepare output data structure
	my $data = generate_tim_data_structure( qw(
		microarray_probes
		ProbeID
		Experiment
		Control
		Ratio
	) );
	
	# Add metadata
	$data->{'program'}      = $0;
	
	# Add sources
	$data->{1}{'source'}    = join(",", @{ $data_sets{'experiment'} } );
	$data->{2}{'source'}    = join(",", @{ $data_sets{'control'} } );
	
	# Add log2 status
	$data->{1}{'log2'}      = $log;
	$data->{2}{'log2'}      = $log;
	$data->{3}{'log2'}      = $log;
	
	# Add normalization list
	if ($norm_list) {
		# normalization list
		$data->{1}{'normalization_list'} = $norm_list;
		$data->{2}{'normalization_list'} = $norm_list;
	}
	
	# Add median target
	if ($median_target) {
		# target for median scaling
		$data->{1}{'median_target'} = $median_target;
		$data->{2}{'median_target'} = $median_target;
	}
	
	# Add quantile normalization
	if ($quant) {
		# quantile normalization performed
		if ($separate) {
			# each channel separately
			$data->{1}{'quantile_normalize'} = 'separate';
			$data->{2}{'quantile_normalize'} = 'separate';
		}
		else {
			# both channels together
			$data->{1}{'quantile_normalize'} = 'together';
			$data->{2}{'quantile_normalize'} = 'together';
		}
	}
	else {
		# no quantile normalization done
		$data->{1}{'quantile_normalize'} = 'none';
		$data->{2}{'quantile_normalize'} = 'none';
	}
	
	# collect and combine the feature data
	foreach my $probe (sort {$a cmp $b} keys %probes) {
		
		# Collect Experimental and Control values
			# first take the mean of the microarray data set values for each 
			# probe feature on the microarray, then take the mean of all 
			# the probe features (there may only be one)
		my $experiment = mean(
			map { 
				mean( @{ $evalues{$_} } )  
			} values %{ $probes{$probe} }
		);
		my $control = mean(
			map { 
				mean( @{ $cvalues{$_} } )  
			} values %{ $probes{$probe} }
		);
		
		# calculate ratio
		my $ratio;
		if ($control == 0 or $control eq '') {
			$ratio = '.';
		} 
		else {
			$ratio = $experiment / $control;
		}
		
		# calculate log2 ratios if requested
		if ($log) { 
			$experiment = ( log($experiment) / LOG2 );
			$control    = ( log($control) / LOG2 );
			unless ($ratio eq '.') {
				$ratio  = ( log($ratio) / LOG2 );
			}
		}
		
		# store the data
		push @{ $data->{'data_table'} }, 
			[ ($probe, $experiment, $control, $ratio) ];
	}
	
	# Update the row count
	$data->{'last_row'} = scalar @{ $data->{'data_table'} } - 1;
	
	# done
	return $data;
}



#### Output single channel records
sub output_single_records {
	
	# prepare output data structure
	my $data = generate_tim_data_structure( qw(
		microarray_probes
		ProbeID
		Experiment
	) );
	
	# Add metadata
	$data->{'program'}      = $0;
	
	# Add sources
	$data->{1}{'source'}    = join(",", @{ $data_sets{'experiment'} } );
	
	# Add log2 status
	$data->{1}{'log2'}      = $log;
	
	# Add normalization list
	if ($norm_list) {
		# normalization list
		$data->{1}{'normalization_list'} = $norm_list;
	}
	
	# Add median target
	if ($median_target) {
		# target for median scaling
		$data->{1}{'median_target'} = $median_target;
	}
	
	# Add quantile normalization
	if ($quant) {
		# quantile normalization performed
		$data->{1}{'quantile_normalize'} = 'yes';
	}
	else {
		# no quantile normalization done
		$data->{1}{'quantile_normalize'} = 'none';
	}
	
	# collect and combine the feature data
	foreach my $probe (sort {$a cmp $b} keys %probes) {
		
		# Collect Experimental
			# first take the mean of the microarray data set values for each 
			# probe feature on the microarray, then take the mean of all 
			# the probe features (there may only be one)
		my $experiment = mean(
			map { 
				mean( @{ $evalues{$_} } )  
			} values %{ $probes{$probe} }
		);
		
		# calculate log2 ratios if requested
		if ($log) { 
			$experiment = ( log($experiment) / LOG2 );
		}
		
		# store the data
		push @{ $data->{'data_table'} }, 
			[ ($probe, $experiment) ];
	}
	
	# Update the row count
	$data->{'last_row'} = scalar @{ $data->{'data_table'} } - 1;
	
	# done
	return $data;
}



#### Perform regression on combined data means
sub perform_regression_on_experiment_control {
	
	my $data = shift;
	
	# initialize
	my $stat = Statistics::Descriptive::Full->new(); 
	my @experiment;
	my @control;
	
	# collect the data
	for my $i (1 .. $data->{'last_row'}) {
		# hard coded indices because we just made the table
		push @experiment, $data->{'data_table'}->[$i][1];
		push @control,    $data->{'data_table'}->[$i][2];
	}
	
	# calculate statistics
	$stat->add_data(@experiment);
	my ($q, $m, $r, $rms) = $stat->least_squares_fit(@control);
		# intercept, slope, Pearson r, and root-mean-square are returned
	
	# report
	push @report, join("\t", 
		(
			'Final',
			'Control',
			'Experiment',
			sprintf("%.3f", $q),
			sprintf("%.3f", $m),
			sprintf("%.3f", $r),
			sprintf("%.3f", $r ** 2),
		)
	) . "\n";
}



#### Perform regression on many replicates
sub perform_regression_on_replicates {
	
	# get the passed data set names and values hash references
	my ($group, $data_values) = @_;
	
	# we will be walking systematically through the name list
	# the first data set will be x, and will be compared to the subsequent data set, y
	# the y data set will then be incremented until all data sets have been compared 
	# the x data set will then be incremented by 1, and the process repeated
	my $number = scalar @{ $data_values->{1} };
	return unless ($number > 1);
	
	# compare the replicates
	for (my $x = 0; $x < $number - 1; $x++) { 
		# starting the x data set at 0
		for (my $y = 1; $y < $number; $y++) { 
			# starting the y data set at 1 to avoid 0
			
			# skip pairs we've already done
			next if ($y <= $x); 
			
			# collect the x and y values into arrays
			my @xvalues; 
			my @yvalues;
			foreach my $id (keys %{$data_values}) {
				push @xvalues, 
					$data_values->{$id}->[$x];
				push @yvalues, 
					$data_values->{$id}->[$y];
			}
			
			# calculate statistics for experimental values
			my $stat = Statistics::Descriptive::Full->new();
			$stat->add_data(@yvalues);
			my ($q, $m, $r, $rms) = $stat->least_squares_fit(@xvalues);
			# intercept, slope, Pearson r, and root-mean-square are returned
			
			# write out to report
			push @report, join("\t", 
				(
					$group,
					$x,
					$y,
					sprintf("%.3f", $q),
					sprintf("%.3f", $m),
					sprintf("%.3f", $r),
					sprintf("%.3f", $r ** 2),
				) 
			) . "\n";
		}
	}
}


__END__

=head1 NAME

process_microarray.pl

A script to quantile normalize two or more raw microarray data sets.

=head1 SYNOPSIS

process_microarray.pl --in <filename1>,<filename2> --channel <redgreen> 
--out <filename>

process_microarray.pl --recipe <file>
  
  Options:
  --in <filename>
  --recipe <filename>
  --out <filename>
  --channel <redgreen> [e | c | n]
  --median <integer>
  --nlist <filename>
  --(no)norm
  --separate
  --(no)log
  --processed
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input microarray file(s). Either provide a single comma-delimited 
list, or use this argument repeatedly for each file. Agilent text files, 
Nimblegen .pair files, or GenePix .gpr files are accepted. The file(s) may be 
compressed with gzip.

=item --recipe <filename>

Specify a text file containing the list of input files to process, along 
with the channel assignment for experiment and control microarray values. 
See below for the format.

=item --out <filename>

Specify the output filename. By default it uses the basename of the recipe 
filename, if available.

=item --channel [e|c|n]

Specify (e)xperiment, (c)ontrol or (n)one designation for each channel 
of each input file. For two-channel microarrays, provide two designations 
for red (Cy5) and green (Cy3) channels, respectively. For one-channel 
microarrays, provide only one designation per file. There must be one 
designation for each specified file. Use this argument repeatedly for 
each file, or provide a single comma-delimited list. 

=item --median <integer>

Optionally specify the target value to which the microarray probes will 
be median scaled. The default value is no scaling. 

=item --nlist <filename>

Optionally specify the filename of a text file containing a list of control 
probe ID names to which all the rest of the probes will be normalized 
against. The control probe median value will be used as the scaling target. 

=item --(no)norm

Optionally indicate whether to quantile normalize (or not) the values  
between microarray sets. The default is true (quantile normalize values).

=item --separate

Optionally indicate that experiment and control microarray data sets 
should be quantile normalized separately and not together. The default 
is false (quantile everything together).

=item --(no)log

Optionally indicate that the experiment, control, and ratio values 
should (not) be converted to a log2 ratio. The default is true.

=item --processed

Optionally indicate that the Processed values in an Agilent text file 
should be used. The default is to use the Raw values.

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will process raw microarray data files and return normalized 
values for each named probe. It will quantile normalize replicate datasets, 
median scale the values, and calculate an experiment/control log2 ratio. It 
supports Agilent Feature Extraction files, GenePix files, and NimbleGen 
pair files.

It will accept as many input files as you provide, although quantile 
normalization requires at least two replicates. Experiment and control data 
sets may be quantile normalized together, or separately. Experiment and 
control data sets are not limited to the same data file (for two-channel 
microarray hybridizations), giving creative flexibility in assigning 
data files and microarray hybridizations. Control data sets may also be 
excluded, in which case all provided data sources are treated as experiment.

Multi-slide hybridizations are easily taken into account if there are not 
duplicate probes across the slides with identical feature numbers. Each 
microarray probe is assigned a unique number based on probe name and 
feature number, and each one must have the same number of data set or 
replicate values. Otherwise, normalization cannot proceed. If you get this 
error message, then you may have to filter the raw files first, or process 
separately.

If a reference control set of oligo probes are included in the microarray 
design, their probe names may be supplied as an external normalization file, 
in which case the median value of those probes are used as the target for 
median scaling of the entire array. 

The program will determine correlation and r-squared values between datasets 
and report them. For multi-slide sets, it is best if all slides in a set are 
listed sequentially for this to work best. 

The program will write out a tab-delimited data file consisting of 
normalized data values for experiment and control, and the experiment/control
ratio for each probeID. It will also write out a separate report text file.

=head1 RECIPE FILE

To facilitate complex arrangements of data files, hybridizations, channel 
swapping, and differing numbers of experiment and control data sets, a 
recipe file may be generated. This is a simple, two-column (whitespace 
delimited) text file. 

The first column is the data source path and filename.

The second column is one or two letter abbreviation specifying the 
channel and designation, following the nomenclature under the channel 
option above. Use one letter for one-channel (red or green) data files, 
or two letters for two-channel (red and green) data files. Use e for 
experiment, c for control, or n for none. For two-channel files, specify 
red first, then green.

Blank lines and comment lines (prefix #) are ignored.

An example is below
 
 # two-channel file, red is experiment, green is control
 path/to/file.txt	ec
 # two-channel file, red is experiment, green is not used
 path/to/file.txt	en
 
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
