#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename qw(fileparse);
use Statistics::Lite qw(mean median sum);
use Bio::ToolBox::data_helper qw(index_data_table);
use Bio::ToolBox::db_helper qw(
	get_new_genome_list
);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	write_tim_data_file
	open_to_read_fh
	convert_and_write_to_gff_file
);
eval {
	# check for bam support
	require Bio::ToolBox::db_helper::bam;
	Bio::ToolBox::db_helper::bam->import;
};
my $VERSION = '1.15';


print "\n This script will generate genomic binned data\n\n";


### Quick help
unless (@ARGV) {
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Get command line options and initialize values
my (
	$datafile, 
	$infile, 
	$new, 
	$dbname,
	$outfile,
	$win,
	$step,
	$method, 
	$paired,
	$span,
	$midpoint,
	$shift,
	$interpolate,
	$log,
	$gffout,
	$gz,
	$help,
	$print_version,
);


# Command line options
GetOptions( 
	'data=s'      => \$datafile, # the source data files
	'in=s'        => \$infile, # specify a pre-existing genomic bin data file
	'new'         => \$new, # generate a new genomic bin data file
	'db=s'        => \$dbname, # specify the database name
	'out=s'       => \$outfile, # name of output file 
	'win=i'       => \$win, # the window size of the genomic bins
	'step=i'      => \$step, # the step size for moving the window
	'method=s'    => \$method, # the method of binning the data
	'paired'      => \$paired, # bam data is paired
	'span'        => \$span, # assign value across entire feature span
	'midpoint'    => \$midpoint, # assign value at the feature midpoint
	'shift=i'     => \$shift, # 3' shift value
	'interpolate' => \$interpolate, # interpolate missing data
	'log!'        => \$log, # values in log2 space
	'gff'         => \$gffout, # output a gff file
	'gz!'         => \$gz, # compress output
	'help'        => \$help, # request help
	'version'     => \$print_version, # print the version
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
	print " Biotoolbox script bin_genomic_data.pl, version $VERSION\n\n";
	exit;
}



### Check for required values
unless ($datafile) {
	die "  OOPS! No source data file specified! \n use --help\n";
}
unless ($infile or $new) {
	die "  OOPS! Must specify bin data file or use new! \n use --help\n";
}
if ($new) {
	unless ($dbname and $outfile and $win) {
		die "  OOPS! missing parameters for new bin data! \n use --help\n";
	}
}
unless ($method =~ /enumerate|count|mean|median|sum/) {
	die "  OOPS! unknown method '$method' specified! \n use --help\n";
}
if ($method eq 'enumerate') {
	# count is just an alias, make internally consistent
	$method = 'count';
}
if (defined $log) {
	if ($log and $method eq 'count') {
		die "  OOPS! method 'count' must be used with non-log data!\n";
	}
}
else {
	$log = 0;
}
unless (defined $span) {$span = 0}
unless (defined $midpoint) {$midpoint = 0}
	# without these two options set we'll just use the read's start point





### Prepare genomic bins
my ($main_data, $column, $data_type) = prepare_data_structure(); 
	# return the main data structure, column index, and data type

# simple reference to the data table
my $table_ref = $main_data->{'data_table'};




### Collecting data into the genomic bins

# Determine the feature processing method
my $process_feature;
if ($span) {
	# collect data across the entire span of the feature
	$process_feature = \&process_feature_span;
}
elsif ($midpoint) {
	# process the feature midpoint
	$process_feature = \&process_feature_midpoint;
}
else {
	# default to process the feature startpoint
	$process_feature = \&process_feature_startpoint;
}

# Walk through the source data file
	# how we collect is based on the file type
	# determine file type based on extension
if ($data_type =~ /gff/) {
	get_gff_data();
} 
elsif ($data_type =~ /bed/) {
	get_bed_data();
} 
elsif ($data_type =~ /sgr/) {
	get_sgr_data();
} 
elsif ($data_type =~ /bam/) {
	get_bam_data();
} 
elsif ($data_type =~ /wig|bdg/) {
	get_wig_data();
} 
else {
	die " unknown source file type!\n";
}


# Secondary processing
secondary_processing();


# Interpolate missing binned data from neighbors
if ($interpolate) {
	interpolate_data();
}




### Write output
# set the output file name
unless ($outfile) {
	$outfile = $infile;
}

# write file according to file type
my $write_success;
if ($gffout) {
	# the file should be gff format
	# we will need to convert the data table into gff format
	print " Writing gff format...\n";
	$write_success = convert_and_write_to_gff_file(
			'data'     => $main_data,
			'score'    => $column,
			'source'   => $datafile,
			'method'   => $main_data->{$column}{'name'},
			'filename' => $outfile,
			'version'  => 2,
			'gz'       => $gz,
	);
}
else {
	$write_success = write_tim_data_file(
		'data'      => $main_data,
		'filename'  => $outfile,
		'gz'        => $gz,
	);
}

# conclusion
if ($write_success) {
	print " Wrote file $write_success\n";
}
else {
	warn " Unable to write output file!\n";
}






#### Subroutines ###############################################################


sub prepare_data_structure {
	my $main_data;
	
	# prepare a new data file
	if ($new) {
		print " Generating a new set of genomic bins....\n";
		# generate new data table
		$main_data = get_new_genome_list(
				'db'        => $dbname, 
				'win'       => $win, 
				'step'      => $step,
		);
		unless ($main_data) {
			die "Unable to generate a new data table!\n";
		}
		# set program metadata
		$main_data->{'program'} = $0;
	} 
	
	# otherwise load a pre-existing file
	else {
		print " Loading data file $infile....\n";
		$main_data = load_tim_data_file($infile);
		unless ($main_data) {
			die "Unable to load data file!\n";
		}
		unless ($main_data->{'feature'} eq 'genome') {
			die " The input file was not generated with a genome data set\n";
		}
	}
	
	# index the data table
	print " Indexing data table....\n";
	index_data_table($main_data) or 
		die " unable to index the data table!\n";
	
	# Determine name for the dataset
		# we will simply use the basename of the data file
		# also get the extension for recognizing the file type to use below
	my ($name, $path, $ext) = fileparse($datafile, 
		qw(
			\.gff
			\.gff\.gz
			\.gff3
			\.gff3\.gz
			\.bed
			\.bed\.gz
			\.sgr
			\.sgr\.gz
			\.wig
			\.wig\.gz
			\.bdg
			\.bdg\.gz
			\.bam
		) 
	);
	
	# Set metadata
	my $column = $main_data->{'number_columns'};
	$main_data->{$column} = {
			'datafile' => $datafile,
			'log2'     => $log,
			'method'   => $method,
			'name'     => $name,
			'index'    => $column,
	};
	if ($span) {
		$main_data->{$column}{'span'} = 1;
	}
	elsif ($midpoint) {
		$main_data->{$column}{'reference_point'} = 'midpoint';
	}
	else {
		$main_data->{$column}{'reference_point'} = 'start';
	}
	if ($shift) {
		$main_data->{$column}{'shift'} = $shift;
	}
	
	# insert column header name
	$main_data->{'data_table'}->[0][$column] = $name; 
	
	# update number of columns
	$main_data->{'number_columns'} += 1;
		
	# finished
	# return data structure, the new column index, and the data file extension
	return ($main_data, $column, $ext);
}


sub get_gff_data {
	print " Collecting features from GFF file '$datafile'....\n";
	
	# open file io object
	my $fh = open_to_read_fh($datafile);
	unless ($fh) {
		die "unable to open '$datafile': $!";
	} 
	
	# collect the features
	my $bin_count = 0;
	my $feature_count = 0;
	while (my $line = $fh->getline) {
		if ($line =~ /^#/) {next} # skip comment lines
		chomp $line;
		my ($chromo, $start, $stop, $score) = (split /\t/, $line)[0,3,4,5];
		
		# process the feature
		$bin_count += &{$process_feature}(
			$chromo,
			$start,
			$stop,
			$score
		);
		$feature_count++;
	}
	
	$fh->close;
	print "  $feature_count features loaded, a total of $bin_count bins were modified\n";
}



sub get_bed_data {
	print " Collecting features from BED file '$datafile'....\n";
	
	# open file io object
	my $fh = open_to_read_fh($datafile);
	unless ($fh) {
		die "unable to open '$datafile': $!";
	} 
	
	# collect the features
	my $bin_count = 0;
	my $feature_count = 0;
	while (my $line = $fh->getline) {
		if ($line =~ /^#/) {next} # skip comment lines
		chomp $line;
		my @data = split /\t/, $line;
		
		# interbase conversion of start
		$data[1] += 1;
		
		# the score column is optional
		my $score;
		if (scalar @data >= 5) {
			$score = $data[4];
		}
		
		# process the feature
		$bin_count += &{$process_feature}(
			$data[0],
			$data[1],
			$data[2],
			$score
		);
		$feature_count++;
	}
	
	$fh->close;
	print "  $feature_count features loaded, a total of $bin_count bins were modified\n";
}



sub get_sgr_data {
	print " Collecting features from sgr file '$datafile'....\n";
	
	# open file io object
	my $fh = open_to_read_fh($datafile);
	unless ($fh) {
		die "unable to open '$datafile': $!";
	} 
	
	# collect the features
	my $bin_count = 0;
	my $feature_count = 0;
	while (my $line = $fh->getline) {
		if ($line =~ /^#/) {next} # skip comment lines
		chomp $line;
		my ($chromo, $start, $score) = split /\t/, $line;
		
		# process the feature
		$bin_count += &{$process_feature}(
			$chromo,
			$start,
			$start,
			$score
		);
		$feature_count++;
	}
	
	$fh->close;
	print "  $feature_count features loaded, a total of $bin_count bins were modified\n";
}



sub get_bam_data {
	print " Collecting features from BAM file '$datafile'....\n";
	
	# open bam io object
	unless (exists &open_bam_db) {
		die " unable to load Bam file support! Is Bio::DB::Sam installed?\n"; 
	}
	my $sam = open_bam_db($datafile) or 
		die " unable to open BAM file '$datafile'!\n";
	
	# initialize counts
	my $bin_count = 0;
	my $feature_count = 0;
		
	# collect bam features in segments across the genome
	for my $tid (0 .. $sam->n_targets - 1) {
		# Because I am using the high level interface and collecting the 
		# read data as Bio::DB::Alignment objects, doing this across the 
		# entire genome takes a staggering amount of memory, even using 
		# the iterator (it has to keep it somewhere!)
		# To get around this, we'll be collecting features in small 
		# segments.
		# Furthermore, we'll be skipping any segments not represented 
		# in the index
		
		# check chromosome
		my $refseq = $sam->target_name($tid);
		print "    for chromosome $refseq....\n";
		unless (exists $main_data->{index}{$refseq} ) {
			# if it's not in the index, there's no need to check it
			next;
		}
		
		# walk through chromosome by segment
		my $seq_length = $sam->target_len($tid);
		for (my $pos = 1; $pos <= $seq_length; $pos += 20000) {
			# a window size of 20 kb is completely arbitrary
			
			# determine the segment end
			my $end_pos = $pos + 20500 - 1;
				# we make the segment end a little bigger than the normal 
				# end of the segment (20000) to allow for all of the right 
				# ends to be recovered.
			if ($end_pos > $seq_length) {
				$end_pos = $seq_length;
			}
			
			# check that this window is in the index
			my $pos_check = int( $pos / $main_data->{index_increment} );
			my $end_pos_check = int( $end_pos / 
				$main_data->{index_increment} );
			if (
				$main_data->{index}{$refseq}{$pos_check} or
				$main_data->{index}{$refseq}{$end_pos_check}
			) {
				# at least a portion of the segment is represented in the 
				# index, so go ahead and proceed
				
				# collect the features, either matched pairs or single matches
				if ($paired) {
					# features are read pairs
					my $iterator = $sam->features(
							-type     => 'read_pair',
							-iterator => 1,
							-seq_id   => $refseq,
							-start    => $pos,
							-end      => $end_pos,
					#) or next; # move on if no features found
					);
					unless ($iterator) {
						next;
					}
					
					
					# collect the read pairs
					while (my $pair = $iterator->next_seq) {
						
						# get individual reads
						my ($left, $right) = $pair->get_SeqFeatures;
						unless (defined $right) {
							# the right read may be undefined because it is 
							# outside of the defined segment, i.e. its 
							# start is gt $end_pos
							# therefore we'll skip these and hope we catch 
							# them later
							next;
						}
						
						# skip unmapped
						if ($left->unmapped or $left->munmapped) {next} 
						
						# get coordinates
						my $start = $left->start;
						my $end = $right->end;
						
						# skip pairs that begin outside of our current 
						# coordinates, or the right end begins beyone the 
						# normal window end, 
						# this should ensure we don't count some pairs twice
						if ($start < $pos) {next}
						if ($start > ($pos + 20000 - 1) ) {next}
						
						# process the feature
						$bin_count += &{$process_feature}(
							$refseq,
							$start,
							$end,
							'.'
						);
						$feature_count++;
					}
				}
				else {
					# features are single reads
					my $iterator = $sam->features(
							-type     => 'match',
							-iterator => 1,
							-seq_id   => $refseq,
							-start    => $pos,
							-end      => $end_pos,
					) or next;
					
					# collect the alignments
					while (my $alignment = $iterator->next_seq) {
						# skip unaligned reads
						if ($alignment->unmapped) {next} 
						
						# get coordinates
						my $start = $alignment->start;
						my $end = $alignment->end;
						
						# skip pairs that begin outside of our current 
						# normal coordinates, this should ensure we don't count
						# some pairs twice
						if ($start < $pos) {next}
						if ($start > ($pos + 20000 - 1) ) {next}
						
						# shift if requested
						if ($shift) {
							# shift the read's position by $shift bp towards 3' end
							my $position;
							if ($alignment->strand > 0) {
								# forward strand
								$position = $start + $shift;
							}
							else {
								# reverse strand
								$position = $end - $shift;
							}
							
							$bin_count += &{$process_feature}(
								$refseq,
								$position,
								$position,
								'.'
							);
						}	
						
						else {
							$bin_count += &{$process_feature}(
								$refseq,
								$start,
								$end,
								'.'
							);
						}
						$feature_count++;
					}
				}
			
			}
			
		}
		
	}
	
	print "  $feature_count features loaded, a total of $bin_count bins were modified\n";
}



sub get_wig_data {
	print " Collecting features from wig file '$datafile'....\n";
	
	# open file io object
	my $fh = open_to_read_fh($datafile);
	unless ($fh) {
		die "unable to open '$datafile': $!";
	} 
	
	# initialize variables
	my $bin_count = 0;
	my $feature_count = 0;
	# reusuable variables
	my $chromo;
	my $fixstart;
	my $step;
	my $span = 1; # default of 1 bp
	
	# collect the data
	while (my $line = $fh->getline) {
		my $start; # specific line variables
		my $stop;
		my $score;
		
		# The wiggle file can have 3 different formats: BED format, variable step, 
		# and fixed step. We need to determine whether each line is a definition
		# line or a data line, based on the line's contents and/or number of 
		# elements. The definition lines will fill the reusable variables above
		# and help in filling out the specific variables.
		
		## check the line's contents
		$line =~ s/[\r\n]+$//; # strip all line endings
		my @data = split /\s+/, $line;
		
		# a track line
		if ($data[0] =~ /track/i) {
			# not much useable information in here for us
		}
		
		# a variable step definition line
		elsif ($data[0] =~ /^variablestep$/i) { 
			foreach (@data) {
				if (/chrom=(\w+)/) {$chromo = $1}
				if (/span=(\w+)/)  {$span = $1}
			}
			next;
		} 
		
		# a fixed step definition line
		elsif ($data[0] =~ /^fixedstep$/i) { 
			foreach (@data) {
				if (/chrom=(\w+)/) {$chromo = $1}
				if (/span=(\w+)/)  {$span = $1}
				if (/start=(\w+)/) {$fixstart = $1}
				if (/step=(\w+)/)  {$step = $1}
			}
			next;
		} 
		
		# a BED data line
		elsif (scalar @data == 4) {
			$chromo = $data[0];
			$start = $data[1] + 1; # the BED line alone uses 0-based indexing
			$stop = $data[2];
			$score = $data[3];
		} 
		
		# a variable step data line
		elsif (scalar @data == 2) { 
			unless ($chromo) { 
				die "Bad formatting! variable step data but chromosome not defined!\n";
			}
			$start = $data[0];
			$stop = $start + $span - 1;
			$score = $data[1];
		} 
		
		# a fixed step data line
		elsif (scalar @data == 1) { 
			unless ($chromo) { 
				die "Bad formatting! fixed step data but chromosome not defined!\n";
			}
			unless ($fixstart) { 
				die "Bad formatting! fixed step data but start not defined!\n";
			}
			unless ($step) { 
				die "Bad formatting! fixed step data but step not defined!\n";
			}
			$start = $fixstart;
			$fixstart += $step; # prepare for next round
			$stop = $start + $span - 1;
			$score = $data[0];
		}
		
		# process the feature
		$bin_count += &{$process_feature}(
			$chromo,
			$start,
			$stop,
			$score
		);
		$feature_count++;
	}
	
	$fh->close;
	print "  $feature_count features loaded, a total of $bin_count bins were modified\n";
	

}



sub process_feature_startpoint {
	
	# collect the feature elements
	my ($chromo, $start, $stop, $score) = @_;
	
	return process_feature($chromo, $start, $score);
}



sub process_feature_midpoint {
	
	# collect the feature elements
	my ($chromo, $start, $stop, $score) = @_;
	
	# determine the feature's midpoint as the reference position
	my $position;
	if ($start == $stop) { 
		# already have a midpoint
		$position = $start;
	} else { 
		# calculate the midpoint
		$position = int( ($start + $stop)/2 );
	}
	
	return process_feature($chromo, $position, $score);
}



sub process_feature_span {
	
	# collect the feature elements
	my ($chromo, $start, $stop, $score) = @_;
	
	# Process at each position along the feature
	my $found = 0;
	for (
		my $position = $start; 
		$position <= $stop; 
		$position += $main_data->{1}{'step'}
	) {
		# process 
		$found += process_feature($chromo, $position, $score);
	}
	
	return $found;
}


sub process_feature {
	
	# collect the feature elements
	my ($chromo, $position, $score) = @_;
	
	# calculate the index value
		# this is used to quickly find the point in the data table 
		# where the feature's coodinates may be found
	my $index_value = int( $position / $main_data->{'index_increment'} );
	my $starting_row = $main_data->{'index'}->{$chromo}{$index_value} 
		|| undef;
	unless (defined $starting_row) {
		# no index, return 0
		return 0;
	} 
	
	# calculate the last reasonable position
	my $last = $position + $main_data->{1}{'win'} + 1;
	
	# Process according to the method
	my $found = 0;
	if ($method eq 'count') {
		# enumeration is easy, we simply add the count
		for (
			my $row = $starting_row; 
			$row < $main_data->{'last_row'}; 
			$row++
		) {
			# check each row to see whether our position fits here
			if (
				# midpoint falls within this window
				$table_ref->[$row][0] eq $chromo and
				$table_ref->[$row][1] <= $position and
				$table_ref->[$row][2] >= $position
			) {
				$table_ref->[$row][$column] += 1;
				$found++;
			}
			elsif ($table_ref->[$row][0] ne $chromo) {
				# looks like we've moved off this chromosome
				last;
			}
			elsif ($table_ref->[$row][1] > $last) {
				# we've moved past the region of interest
				last;
			}
		}
	}
	else {
		# other methods require that we collect all the values first
		# before actually calculating
		# therefore we'll be adding the scores to each bin
		# scores will be kept as a string to keep the data_table simpler
		# we'll split them later
		for (
			my $row = $starting_row; 
			$row < $main_data->{'last_row'}; 
			$row++
		) {
			# check each row to see whether our position fits here
			if (
				$table_ref->[$row][0] eq $chromo and
				$table_ref->[$row][1] <= $position and
				$table_ref->[$row][2] >= $position
			) {
				if (defined $table_ref->[$row][$column]) {
					$table_ref->[$row][$column] .= ";$score";
				}
				else {
					$table_ref->[$row][$column] = $score;
				}
				$found++;
			}
			elsif ($table_ref->[$row][0] ne $chromo) {
				# looks like we've moved off this chromosome
				last;
			}
			elsif ($table_ref->[$row][1] > $last) {
				# we've moved past the region of interest
				last;
			}
		}
	}
	
	return $found;
}



sub secondary_processing {
	
	# mean
	if ($method eq 'mean') {
		print " Calculating genomic bin mean values....\n";
		# we have only collected the raw data
		# now need to calculate the mean values
		for (my $row = 1; $row <= $main_data->{'last_row'}; $row++) {
			if (defined $table_ref->[$row][$column]) {
				# we have raw data, split and determine mean
				$table_ref->[$row][$column] = _calculate_mean(
					split(/;/, $table_ref->[$row][$column]) );
			}
			else {
				# no data
				$table_ref->[$row][$column] = '.';
			}
		}
	}
	
	# median
	elsif ($method eq 'median') {
		print " Calculating genomic bin median values....\n";
		# we have only collected the raw data
		# now need to calculate the median values
		for (my $row = 1; $row <= $main_data->{'last_row'}; $row++) {
			if (defined $table_ref->[$row][$column]) {
				# we have raw data, split and determine median
				$table_ref->[$row][$column] = _calculate_median(
					split(/;/, $table_ref->[$row][$column]) );
			}
			else {
				# no data
				$table_ref->[$row][$column] = '.';
			}
		}
	}
	
	# sum
	elsif ($method eq 'sum') {
		print " Calculating genomic bin sums....\n";
		# we have only collected the raw data
		# now need to sum the values
		for (my $row = 1; $row <= $main_data->{'last_row'}; $row++) {
			if (defined $table_ref->[$row][$column]) {
				# we have raw data, split and determine sum
				$table_ref->[$row][$column] = _calculate_sum(
					split(/;/, $table_ref->[$row][$column]) );
			}
			else {
				# no data
				$table_ref->[$row][$column] = 0;
			}
		}
	}

	# count
	elsif ($method eq 'count') {
		# make sure we have a value in each bin
		for (my $row = 1; $row <= $main_data->{'last_row'}; $row++) {
			unless (defined $table_ref->[$row][$column]) {
				# no count was here
				# then place a zero in here
				$table_ref->[$row][$column] = 0;
			}
		}
	}
}



sub interpolate_data {
	print " Interpolating missing data....\n";
	
	# count how many bins are interpolated
	my $interpcount = 0; 
	
	# identify the first data row
	my $row = 1;
	while ($row < $main_data->{'last_row'}) {
		if ($table_ref->[$row][$column] eq '.') {
			$row++;
		}
		else {
			last;
		}
	}
	
	# now start interpolating missing data
	INTERPOLATE_LOOP:
	while ($row < $main_data->{'last_row'}) {
		
		# check for null values
		if ($table_ref->[$row][$column] eq '.') {
			# we have a null value
			
			# now we need to look for the next non-value
			my $next_data = $row + 1;
			while ($next_data < $main_data->{'last_row'}) {
				
				# check for same chromosome
				if ($table_ref->[$next_data][0] ne $table_ref->[$row][0]) {
					# we've come to the next chromosome
					# that means we won't find another valid data point to 
					# interpolate from, so need to move on
					$row = $next_data;
					undef $next_data;
					last;
				}
				
				if ($table_ref->[$next_data][$column] eq '.') {
					# still a null value
					$next_data++;
				}
				else {
					# we have a value
					last;
				}
			}
			if (!$next_data) {
				# this means we moved on to the next chromosome
				# need to find the first data point for this chromosome
				while ($row < $main_data->{'last_row'}) {
					if ($table_ref->[$row][$column] eq '.') {
						$row++;
					}
					else {
						last;
					}
				}
				next INTERPOLATE_LOOP;
			}
			if ($next_data == $main_data->{'last_row'}) {
				# that's it, we're done
				last INTERPOLATE_LOOP;
			}
			
			
			# calculate the fractional value
			my $number = $next_data - ($row - 1); 
				# the number of rows between between good datapoints that need
				# to be interpolated
			if ($number == 0) {
				warn " encounted 0 again at row $row and next row $next_data!\n";
				$row++;
				next INTERPOLATE_LOOP;
			}
			my $fractional_value = 
				( $table_ref->[$next_data][$column] - 
					$table_ref->[$row - 1][$column] ) / $number;
			
			# replace all the null values
			for (my $i = $row; $i < $next_data; $i++) {
				
				# calculate the new fractional value for this row
				my $new_value = $table_ref->[$row - 1][$column] +
					( $fractional_value * ($i - ($row - 1)) );

				# replace
				$table_ref->[$i][$column] = $new_value;
				$interpcount++;
			}
			
			# reset the row value
			$row = $next_data + 1;
		}
		
		else {
			# we have a value, move on
			$row++;
		}
		
	}
	
	print "  Interpolated $interpcount bins\n";
}



sub _calculate_mean {
	# a subroutine to calculate the mean
	# takes into account the log status
	if ($log) { 
		# working with log2 numbers
		my @data = map { 2**$_ } @_; # de-log
		return log( mean(@data) ) / log(2);
	}
	else {
		return mean(@_);
	}
}


sub _calculate_median {
	# a subroutine to calculate the median
	# takes into account the log status
	if ($log) { 
		# working with log2 numbers
		my @data = map { 2**$_ } @_; # de-log
		return log( median(@data) ) / log(2);
	}
	else {
		return median(@_);
	}
}


sub _calculate_sum {
	# a subroutine to calculate the sum
	# takes into account the log status
	if ($log) { 
		# working with log2 numbers
		my @data = map { 2**$_ } @_; # de-log
		return log( sum(@data) ) / log(2);
	}
	else {
		return sum(@_);
	}
}



__END__

=head1 NAME

bin_genomic_data.pl

A script to bin genomic data into windows.

=head1 SYNOPSIS
 
 bin_genomic_data.pl --data <filename> --in <filename> --method [method]
 
 bin_genomic_data.pl --data <filename> --new --win <integer> --db <db_name>
     --out <filename> --method [method]
 
  Options:
  --data <filename>
  --in <filename>
  --new
  --method [count | mean | median |sum]
  --out <filename>
  --db <database_name>
  --win <integer>
  --step <integer>
  --paired
  --span
  --midpoint
  --shift <integer>
  --interpolate
  --log
  --gff
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --data <filename>

Specify the source genomic data file containing the datapoints to be 
binned. Supported files include .gff, .bed, , .wig, .sgr, and .bam files. 
Text files may be gzipped.

=item --in <filename>

Specify the filename of an existing data table of genomic bins or intervals. 
It may be generated using this program or B<get_datasets.pl>. A Bed file 
may be acceptable. Required unless C<--new> is provided. The file may be 
gzip compressed. 

=item --new

Indicate that a new data file should be generated. Required unless C<--in> 
is used.

=item --method [count | mean | median |sum]

Indicate which method should be used to collect the binned data. 
Acceptable values include:
 
 - count         Count the number of unique features in the bin.
 - mean          Take the mean score of features in the bin.
 - median        Take the median score of features in the bin.
 - sum           Sum the scores of features in the bin

=item --out <filename>

Specify the output file name. Required if generating a new binned 
data file, otherwise if an input file was specified it will be 
overwritten.

=item --db <database_name>

Specify a database. This is used only to retrieve chromosome 
information to generate genomic bins. Required if --new is used.

=item --win <integer>

Specify the genomic bin window size in bp. Required if --new is used.

=item --step <integer>

Specify the step size for moving the window when generating bins. 
The default is to equal the window size.

=item --paired

Indicate that the data bam file consists of paired-end alignments 
rather than single-end alignments.

=item --span

Assign feature value (or count) at each genomic bin across the 
feature. This is dependent on the metadata step value when the 
genomic bins were generated.

=item --midpoint

Assign the the feature value (or count) at the feature's midpoint 
rather than the default feature's start point.

=item --shift <integer>

For stranded single-end alignments from a bam data file, the start 
position may be shifted in the 3' direction by the indicated 
number of bp. This is to account for ChIP-Seq data where the peak 
of tag counts is offset from the actual center of the sequenced 
fragments. Use a shift value of 1/2 the mean fragment length of the 
sequencing library.

=item --interpolate

Optionally interpolate missing bin values from flanking bins. 

=item --log

Flag to indicate that source data is log2, and to calculate 
accordingly and report as log2 data.

=item --gff

Write a new gff output file instead of a normal tim data file.

=item --gz

Compress the output file through gzip.

=item --version

Print the version number.

=item --help

Display the POD documentation.

=back

=head1 DESCRIPTION

This program will convert genomic data into genomic binned data. The genomic 
bins are generated from the indicated database as segments of the genome of 
the specified window size. The step size may be optionally specified, or 
defaults to the window size. The source data may be collected from multiple 
sources, including GFF, BED, BAM, WIG, or SGR files. The data may be 
collected in one of two ways.

First, the program can enumerate (or count) features which overlap with a 
genomic bin. Sequence tag (BAM file) alignments (both single- and 
paired-end), BED, and GFF features may be counted. The feature start point, 
or optionally the midpoint, is used to record the feature. Alternatively, the 
counts may be enumerated across the entire span of the feature when using 
small (single bp) bins.

Second, the program can statistically combine the scores (values) of features
that overlap each genomic bin. Feature scores may be collected from GFF, BED, 
WIG, or SGR files. Either the mean, median, or sum value can be recorded 
for the genomic bin. 

Genomic bins lacking a score may be interpolated from neighboring bins. Up to
four consecutive non-value bins may interpolated from neighboring bins that
contain values.

The program reads/writes a data file compatable with the 'get_datasets.pl' 
program. It can optionally write a new GFF file based on the genomic binned
data.

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
