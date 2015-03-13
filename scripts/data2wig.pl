#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(mean median sum max);
use Bio::ToolBox::data_helper qw(
	find_column_index
);
use Bio::ToolBox::big_helper qw(wig_to_bigwig_conversion);
use Bio::ToolBox::file_helper qw(
	open_tim_data_file
	open_to_write_fh
);
my $VERSION = '1.15';

print "\n This script will export a data file to a wig file\n\n";


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
	$infile, 
	$outfile,
	$step,
	$bedgraph,
	$step_size,
	$span,
	$chr_index,
	$start_index,
	$stop_index,
	$score_index,
	$track_name,
	$use_track,
	$midpoint,
	$interbase,
	$format,
	$method,
	$log2,
	$bigwig,
	$bw_app_path,
	$database,
	$chromo_file,
	$keep,
	$gz,
	$help,
	$print_version,
);


# Command line options
GetOptions( 
	'in=s'      => \$infile, # name of input file
	'out=s'     => \$outfile, # name of output gff file 
	'step=s'    => \$step, # wig step method
	'bed|bdg!'  => \$bedgraph, # write a bedgraph file
	'size=i'    => \$step_size, # wig step size
	'span=i'    => \$span, # the wig span size
	'chr=i'     => \$chr_index, # index for the chromosome column
	'start|pos=i' => \$start_index, # index for the start column
	'stop|end=i'=> \$stop_index, # index for the stop column
	'index|score=i' => \$score_index, # index for the score column
	'name=s'    => \$track_name, # name string for the track
	'track!'    => \$use_track, # boolean to include a track line
	'mid!'      => \$midpoint, # boolean to use the midpoint
	'zero|inter!' => \$interbase, # shift from interbase
	'format=i'  => \$format, # format output to indicated number of places
	'method=s'  => \$method, # method for combining duplicate values
	'log!'      => \$log2, # data is in log2 format
	'bigwig|bw' => \$bigwig, # generate a binary bigwig file
	'db=s'      => \$database, # database for bigwig file generation
	'chromof=s' => \$chromo_file, # name of a chromosome file
	'bwapp=s'   => \$bw_app_path, # path to wigToBigWig utility
	'keep!'     => \$keep, # keep the wig file after converting to bigWig
	'gz!'       => \$gz, # boolean to compress output file
	'help'      => \$help, # request help
	'version'   => \$print_version, # print the version
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
	print " Biotoolbox script data2wig.pl, version $VERSION\n\n";
	exit;
}



### Check for required or default values
unless ($infile) {
	$infile = shift @ARGV or
		die "  OOPS! No source data file specified! \n use $0 --help\n";
}
unless (defined $use_track) {
	# default is to write a track
	# this will be changed later if we're writing a bigwig
	$use_track = 1;
}



### Load input file
my ($in_fh, $metadata_ref) = open_tim_data_file($infile);
unless ($in_fh) {
	die "Unable to open data table!\n";
}



### Check and/or ask for specific options

identify_indices();

check_track_name();

check_step();

check_log2();

set_bigwig_options();

my $method_sub = set_method_sub();



### Open output file
unless ($outfile) {
	# automatically generate output file name based on track name
	$outfile = $track_name;
}
unless ($outfile =~ /\.(?:wig|bdg|bedgraph)(?:\.gz)?$/i) {
	# add extension
	$outfile .= $bedgraph ? '.bdg' : '.wig';
}
my $out_fh = open_to_write_fh($outfile, $gz) or 
	die " unable to open output file '$outfile' for writing!\n";

# write track line
if ($use_track) {
	print {$out_fh} "track type=wiggle_0 name=$track_name\n";
}



### Start the conversion 
print " converting '" . $metadata_ref->{$score_index}{'name'} . "'....\n";
if ($bedgraph) {
	convert_to_bedgraph();
}
elsif ($step eq 'fixed') {
	convert_to_fixedStep();
}
elsif ($step eq 'variable') {
	convert_to_variableStep();
}



# close files
$in_fh->close;
$out_fh->close;



### Finish Up
if ($bigwig) {
	# requested to continue and generate a binary bigwig file
	print " temporary wig file '$outfile' generated\n";
	print " converting to bigwig file....\n";
	convert_to_bigwig();
}
else {
	# no big wig file needed, we're finished
	print " finished! wrote file '$outfile'\n";
}




############ Subroutines ###############

sub identify_indices {
	
	# automatically identify the indices based on file type if possible
	
	# gff
	if ( $metadata_ref->{'gff'} ) {
		# these indices are assumptions based on proper GFF format
		$chr_index   = 0 unless defined $chr_index;
		$start_index = 3 unless defined $start_index;
		$stop_index  = 4 unless defined $stop_index;
		$score_index = 5 unless defined $score_index;
	}
	
	# bedgraph
	elsif ( 
		$metadata_ref->{'bed'} == 4 and
		$metadata_ref->{'extension'} =~ /graph|bdg/
	) {
		# a bedgraph format
		$chr_index   = 0 unless defined $chr_index;
		$start_index = 1 unless defined $start_index;
		$stop_index  = 2 unless defined $stop_index;
		$score_index = 3 unless defined $score_index;
		
		# automatically set interbase
		$interbase = 1;
	}
	
	# traditional bed
	elsif ( 
		$metadata_ref->{'bed'} >=5
	) {
		# a bed format, using the score column
		$chr_index   = 0 unless defined $chr_index;
		$start_index = 1 unless defined $start_index;
		$stop_index  = 2 unless defined $stop_index;
		$score_index = 4 unless defined $score_index;
		
		# automatically set interbase
		$interbase = 1;
	}
	
	# sgr
	elsif ( $metadata_ref->{'extension'} =~ /sgr/ ) {
		# SGR format only has start and score
		$chr_index   = 0 unless defined $chr_index;
		$start_index = 1 unless defined $start_index;
		$score_index = 2 unless defined $score_index;
	}
	
	# non-standard text file or tim data format text file
	else {
		# we will automatically look for the coordinate columns
		
		# chromosome
		unless (defined $chr_index) {
			$chr_index = find_column_index($metadata_ref, '^chr|seq|ref');
			
			# this is a required index
			unless (defined $chr_index) {
				die " No chromosome or sequence ID column" . 
						" found in the file!!!\n  Please specify with" . 
						" the --chr option\n";
			}
		}
		
		# start
		unless (defined $start_index) {
			$start_index = find_column_index($metadata_ref, '^start|pos');
			
			# this is a required index
			unless (defined $start_index) {
				# still no start? how about a midpoint? position?
				# keep trying
				$start_index = find_column_index($metadata_ref, 
					'midpoint|mid|position|index');
				
				# fail
				unless (defined $start_index) {
					# still nothing found, fail
					die " No start, midpoint, or position coordinate column" . 
						" found in the file!!!\n  Please specify with" . 
						" the --start option\n";
				}
			}
		}
		
		# stop
		unless (defined $stop_index) {
			# this is not absolutely required
			$stop_index = find_column_index($metadata_ref, '^stop|end');
			
			if ($midpoint and not defined $stop_index) {
				die " Stop index is required to calculate midpoint!\n";
			}
		}
		
		# score
		unless (defined $score_index) {
			# first look for a generic score index
			$score_index = find_column_index($metadata_ref, '^score$');
			
			unless (defined $score_index) {
				# ask the user for help if we can't find it
				# print the column names
				print " These are the column names in the datafile\n";
				for (my $i = 0; $i < $metadata_ref->{'number_columns'}; $i++) {
					print "   $i\t", $metadata_ref->{$i}{'name'}, "\n";
				}
			
				# ask for the score index
				print " Enter the index for the feature score column  ";
				$score_index = <STDIN>;
				chomp $score_index;
			}
		}
	}
	
	# finished with indices
}


sub check_track_name {
	# determine what the track name will be
	
	unless (defined $track_name) {
		if ( 
			$metadata_ref->{'gff'} or
			$metadata_ref->{'bed'} or
			$metadata_ref->{'extension'} =~ /sgr/
		) {
			# if it is a gff/sgr/bed file, use the base file name
			# I have a practice of using the GFF type as the file name
			# this is easier than attempting to read the type column from the gff file
			$track_name = $metadata_ref->{'basename'};
		}
		else {
			# use the name of the score column
			$track_name = $metadata_ref->{$score_index}{'name'};
		}
	}
}


sub check_step {
	# In my biotoolbox scripts that generate genomic bins or windows
	# the win(dow) and step values are recorded as metadata under 
	# the start column
	# Therefore, we will check for this metadata to confirm and/or 
	# determine the step and span parameters
	
	
	# check step 
	if ($bedgraph) {
		# write a bedgraph file
		$step = 'bed';
		$use_track = 0; # no track data
		return;
	}
	elsif ($step eq 'bed') {
		# write a bedgraph file
		$bedgraph = 1;
		$use_track = 0; # no track data
		return;
	}
	elsif ($step eq 'variable') {
		# this is ok, we can work with it
	}
	elsif ($step eq 'fixed') {
		# double check that the data file supports this
		# assign the step size as necessary
		if (exists $metadata_ref->{$start_index}{'step'} ) {
			
			if (defined $step_size) {
				if ($step_size != $metadata_ref->{$start_index}{'step'}) {
					warn " Requested step size $step_size does not match" .
						" metadata step size " . 
						$metadata_ref->{$start_index}{'step'} . "!\n" .
						"  Accurate wig file is not guaranteed!!!\n";
				}
			}
			
			else {
				# define it from the metadata
				$step_size = $metadata_ref->{$start_index}{'step'};
			}
		}
		
		elsif (!defined $step_size) {
			warn " Fixed step size not defined by user or metadata! Using 'variableStep'\n";
			$step = 'variable';
		}
	}
	else {
		# attempt to determine automatically
		if ( exists $metadata_ref->{$start_index}{'step'} ) {
			
			# set step size
			$step = 'fixed';
			$step_size = $metadata_ref->{$start_index}{'step'};
			print " Automatically generating 'fixedStep' wig with " . 
				"step of $step_size bp\n";
		}
		else {
			print " Automatically generating 'variableStep' wig\n";
			$step = 'variable';
		}
	}
	
	# check span
	if ($span) {
		# user set it, confirm with metadata if possible
		if (
			exists $metadata_ref->{$start_index}{'win'} and
			$metadata_ref->{$start_index}{'win'} != $span
		) {
			# the requested span and metadata window size do not match
			print " Requested span size $span does not match" .
						" metadata window size " . 
						$metadata_ref->{$start_index}{'win'} . "!\n" .
						"  Accurate wig file is not guaranteed!!!\n";
		}
	}
	else {
		# attempt to determine automatically 
		if (
			exists $metadata_ref->{$start_index}{'win'} and
			not $midpoint
		) {
			# set the span equal to the window size
			$span = $metadata_ref->{$start_index}{'win'};
			print " Automatically setting span to $span bp\n";
		}
		else {
			# default span is 1 bp
			$span = 1;
		}
	}
	
	# confirm span
	if ($midpoint and $span > 1) {
		# cannot use span parameter if the midpoint is requested
		$span = 1;
		print " Reverting span to 1 bp when midpoint is requested\n";
	}
	if ($step eq 'fixed' and $span and $span > $step_size) {
		# span greater than the step will confuse wig parsers, 
		# including conversion to bigWig files
		# the wig format definition cannot have > 1 value per position
		print " Span of $span bp is greater than step of $step_size bp!\n" . 
			"   Wig cannot have overlapping positions. Switching to midpoint\n";
		$span = 1;
		$midpoint = 1;
	}
}


sub check_log2 {
	unless (defined $log2) {
		# check the metadata for the score dataset
		if (exists $metadata_ref->{$score_index}{'log2'}) {
			$log2 = $metadata_ref->{$score_index}{'log2'};
		}
		else {
			# otherwise default is false
			$log2 = 0;
		}
	}
}


sub set_bigwig_options {
	if ($bigwig) {
		# if we're generating bigwig file, no track is needed
		$use_track = 0;
		
		# force no compression
		$gz = 0;
		
		# check that we have a source for chromosome info
		unless ($database or $chromo_file) {
			if (exists $metadata_ref->{db}) {
				$database = $metadata_ref->{db};
			}
			else {
				die " No database name or chromosome file provided for generating bigwig file!\n";
			}
		}
	}
}


sub set_method_sub {
	
	# for combining values from duplicate positions we need a method
	if ($method eq 'mean') {
		return \&mean;
	}
	elsif ($method eq 'median') {
		return \&median;
	}
	elsif ($method eq 'sum') {
		return \&sum;
	}
	elsif ($method eq 'max') {
		return \&max;
	}
}


sub convert_to_fixedStep {
	
	# keep track of current chromosome name and length
	my $current_chr; # current chromosome
	
	# walk through the data file
	while (my $line = $in_fh->getline) {
		chomp $line;
		my @data = split /\t/, $line;
		
		# adjust for interbase 0-base coordinates
		if ($interbase) {
			$data[$start_index] += 1;
		}
		
		# skip negative or zero coordinates
		next if $data[$start_index] <= 0;
		
		# write definition line if necessary
		if ($data[$chr_index] ne $current_chr) {
			# new chromosome, new definition line
			
			# need to determine start position first
			my $start = calculate_position(@data);
			
			# print definition line
			my $definition = 'fixedStep chrom=' . $data[$chr_index] .
				 " start=$start step=$step_size span=$span";
			$out_fh->print("$definition\n");
				
			
			# reset the current chromosome
			$current_chr = $data[$chr_index];
		}
		
		
		# adjust score formatting as requested
		my $score;
		if ($data[$score_index] eq '.') {
			# an internal null value
			# treat as 0 if fixedStep
			$score = 0;
		}
		else {
			$score = $data[$score_index];
		}
		if (defined $format) {
			# format if requested
			$score = format_score($score);
		}
		
		# check we haven't gone over the chromosome end
		# we assume the last interval would accurately record the chromosome end
		if ( defined $stop_index ) {
			if ( ($data[$start_index] + $step_size - 1) > $data[$stop_index]) {
				# size is too big
				# we will not write this last data point
				# we may have some data loss at the end of the chromosome
				#warn " $current_chr clipped at $data[$stop_index]\n";
				next;
			}
		}
		
		# write fixed data line
		$out_fh->print("$score\n");
	}
}


sub convert_to_variableStep {
	my $current_chr; # current chromosome
	my $previous_pos; # previous position to avoid duplicates in wig file
	my @scores; # reusable array for putting multiple data points in
	while (my $line = $in_fh->getline) {
		chomp $line;
		my @data = split /\t/, $line;
		
		# write definition line if necessary
		if ($data[$chr_index] ne $current_chr) {
			
			# first check and write orphan scores
			# this might happen if there was only one score on the entire chr
			if (@scores) {
				if (scalar @scores == 1) {
					# print the one score
					$out_fh->print("$previous_pos $scores[0]\n");
				}
				else {
					# more than one score
					
					# combine the scores if possible
					if ($method_sub) {
						my $new_score;
						if ($log2) {
							# convert to log2 first
							@scores = map {2 ** $_} @scores;
							# combine and convert back to log2
							$new_score = log( &{$method_sub}(@scores) ) / log(2);
						}
						else {
							$new_score = &{$method_sub}(@scores);
						}
						
						# print the combined score
						$out_fh->print("$previous_pos $new_score\n");
					}
					else {
						die " there are " . scalar(@scores) . " scores for " . 
							"position $current_chr $previous_pos!!!!\n" .
							" Please define a combination method!!! see help\n";
					}
				}
			}
			
			# new chromosome, new definition line
			my $definition = 'variableStep chrom=' . $data[$chr_index] . 
				" span=$span";
			$out_fh->print("$definition\n");
			
			# reset the current chromosome
			$current_chr = $data[$chr_index];
			$previous_pos = undef;
			@scores = ();
		}
		
		
		# adjust for interbase 0-base coordinates
		if ($interbase) {
			$data[$start_index] += 1;
		}
		
		# collect the score
		my $score;
		if ($data[$score_index] eq '.') {
			# internal null value, skip these
			next;
		}
		if (defined $format) {
			# format if requested
			$score = format_score( $data[$score_index] );
		}
		else {
			# no formatting, take as is
			$score = $data[$score_index];
		}
		
		
		# calculate the position that we will use
		my $position = calculate_position(@data);
			
		# skip negative or zero coordinates
		next if $position <= 0;
		
		# check for duplicate positions and write appropriately
		if (!defined $previous_pos) {
			# new chromosome
			$previous_pos = $position;
			push @scores, $score;
		}
		
		elsif ($position == $previous_pos) {
			# same position, add to the score list
			push @scores, $score;
		}
		
		elsif ($position < $previous_pos) {
			# error!!! unsorted file!!!!
			die " file is not sorted by increasing position!\n   " . 
				"chromosome $current_chr, compare current position " . 
				"$position\n    versus previous position $previous_pos\n";
		}
		
		else {
			# we have moved on to the next position
			# now print the previous scores
			if (scalar @scores == 1) {
				# print the one score
				$out_fh->print("$previous_pos $scores[0]\n");
			}
			else {
				# more than one score
				
				# combine the scores if possible
				if ($method_sub) {
					my $new_score;
					if ($log2) {
						# convert to log2 first
						@scores = map {2 ** $_} @scores;
						# combine and convert back to log2
						$new_score = log( &{$method_sub}(@scores) ) / log(2);
					}
					else {
						$new_score = &{$method_sub}(@scores);
					}
					
					# print the combined score
					$out_fh->print("$previous_pos $new_score\n");
				}
				else {
					die " there are " . scalar(@scores) . " scores for " . 
						"position $current_chr $previous_pos!!!!\n" .
						" Please define a combination method!!! see help\n";
				}
			}
			
			# reset for next
			$previous_pos = $position;
			@scores = ($score);
		}
	}

}


sub convert_to_bedgraph {
	
	# check for indices
	unless (defined $chr_index and defined $start_index and defined $stop_index) {
		die " One or more indices for chromosome, start, or stop is not defined" . 
			" or found!\n Unable to write a bedgraph file!\n";
	}
	
	# variables to check for overlap
	my $current_chr; # current chromosome
	my $previous_pos; # previous position to avoid overlap
	while (my $line = $in_fh->getline) {
		chomp $line;
		my @data = split /\t/, $line;
		
		# adjust start position
		unless ($interbase) {
			$data[$start_index]--;
		}
		
		# check coordinates
		if (defined $previous_pos and defined $current_chr) {
			
			# check if on the same chromosome
			if ($current_chr eq $data[$chr_index]) {
				# check for overlap
				if ($data[$start_index] < $previous_pos) {
					die " There are overlapping intervals or the file is not sorted by" .
						" coordinates!\n Compare $data[$chr_index]:$data[$start_index]" . 
						" with previous stop position $previous_pos\n";
				}
				# otherwise it is ok
				$previous_pos = $data[$stop_index];
			}
			else {
				# new chromosome
				$current_chr = $data[$chr_index];
				$previous_pos = $data[$stop_index];
			}
		}
		else {
			# define the current 
			$current_chr = $data[$chr_index];
			$previous_pos = $data[$stop_index];
		}
		
		# collect the score
		my $score;
		if ($data[$score_index] eq '.') {
			# internal null value, skip these
			next;
		}
		if (defined $format) {
			# format if requested
			$score = format_score( $data[$score_index] );
		}
		else {
			# no formatting, take as is
			$score = $data[$score_index];
		}
		
		# write the feature line
		$out_fh->print(join("\t", $data[$chr_index], $data[$start_index], 
			$data[$stop_index], $score), "\n");
	}
}


sub calculate_position {
	my @data = @_;
	my $position;
	
	if ($midpoint) {
		# user requested to use the midpoint
		if ( 
			defined $stop_index and 
			$data[$start_index] != $data[$stop_index] 
		) {
			# not same position, so calculate midpoint
			$position = int(  
				( ( $data[$start_index] + $data[$stop_index] ) / 2 ) + 0.5 );
		}
		else {
			# same position
			$position = $data[$start_index];
		}
	}
	else {
		# otherwise use the start position
		$position = $data[$start_index];
	}
	
	return $position;
}


sub format_score {
	my $score = shift;
	
	# format the score value to the indicated number of spaces
	if ($format == 0) {
		# no decimal places
		return sprintf( "%.0f", $score);
	}
	elsif ($format == 1) {
		# 1 decimal place
		return sprintf( "%.1f", $score);
	}
	elsif ($format == 2) {
		# 2 decimal places
		return sprintf( "%.2f", $score);
	}
	elsif ($format == 3) {
		# 3 decimal places
		return sprintf( "%.3f", $score);
	}
}


sub convert_to_bigwig {
	
	# perform the conversion
	my $bw_file = wig_to_bigwig_conversion(
			'wig'       => $outfile,
			'db'        => $database,
			'chromo'    => $chromo_file,
			'bwapppath' => $bw_app_path,
	);

	
	# confirm
	if ($bw_file) {
		print " bigwig file '$bw_file' generated\n";
		unlink $outfile unless $keep; # remove the wig file
	}
	else {
		die " bigwig file not generated! see standard error\n";
	}
}





__END__

=head1 NAME

data2wig.pl

A script to convert a generic data file into a wig file.

=head1 SYNOPSIS

data2wig.pl [--options...] <filename> 
  
  Options:
  --in <filename>
  --out <filename> 
  --step [fixed | variable | bed]
  --bed | --bdg
  --size <integer>
  --span <integer>
  --index | --score <column_index>
  --chr <column_index>
  --start | --pos <column_index>
  --stop | --end <column_index>
  --name <text>
  --(no)track
  --mid
  --inter | --zero
  --format [0 | 1 | 2 | 3]
  --method [mean | median | sum | max]
  --log
  --bigwig | --bw
  --chromof <filename>
  --db <database>
  --bwapp </path/to/wigToBigWig>
  --keep
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. Genome 
coordinates are required. The first row should be column headers. Text 
files generated by other B<BioToolBox> scripts are acceptable. Files may 
be gzipped compressed.

=item --out <filename>

Optionally specify the name of of the output file. The track name is 
used as default. The '.wig' extension is automatically added if required.

=item --step [fixed | variable | bed]

The type of step progression for the wig file. Three wig formats are 
available: 
  - fixedStep: where data points are positioned at equal distances 
        along the chromosome
  - variableStep: where data points are variably positioned along 
        the chromosome. 
  - bed (bedGraph): where scores are associated with intervals 
        defined by start and stop coordinates.
The fixedStep wig file has one column of data (score), the variableStep 
wig file has two columns (position and score), and the bedGraph has 
four columns of data (chromosome, start, stop, score). If the option 
is not defined, then the format is automatically determined from the 
metadata of the file.

=item --bed

=item --bdg

Convenience option to specify a bedGraph file should be written. Same as 
specifying --step=bed.

=item --size <integer>

Optionally define the step size in bp for 'fixedStep' wig file. This 
value is automatically determined from the table's metadata, if available. 
If the C<--step> option is explicitly defined as 'fixed', then the step 
size may also be explicitly defined. If this value is not explicitly
defined or automatically determined, the variableStep format is used by
default.

=item --span <integer>

Optionally indicate the size of the region in bp to which the data value 
should be assigned. The same size is assigned to all data values in the 
wig file. This is useful, for example, with microarray data where all of 
the oligo probes are the same length and you wish to assign the value 
across the oligo rather than the midpoint. The default is inherently 1 bp. 

=item --index <column_index>

=item --score <column_index>

Indicate the column index (0-based) of the dataset in the data table 
to be used for the score. If a GFF file is used as input, the score column is 
automatically selected. If not defined as an option, then the program will
interactively ask the user for the column index from a list of available
columns.

=item --chr <column_index>

Optionally specify the column index (0-based) of the chromosome or 
sequence identifier. This is required to generate the wig file. It may be 
identified automatically from the column header names.

=item --start <column_index>

=item --pos <column_index>

Optionally specify the column index (0-based) of the start or chromosome 
position. This is required to generate the wig file. It may be 
identified automatically from the column header names.

=item --start <column_index>

=item --end <column_index>

Optionally specify the column index (0-based) of the stop or end 
position. It may be identified automatically from the column header names.

=item --name <text>

The name of the track defined in the wig file. The default is to use 
the name of the chosen score column, or, if the input file is a GFF file, 
the base name of the input file. 

=item --(no)track

Do (not) include the track line at the beginning of the wig file. Wig 
files normally require a track line, but if you will be converting to 
the binary bigwig format, the converter requires no track line. Why it 
can't simply ignore the line is beyond me. This option is automatically 
set to false when the C<--bigwig> option is enabled.

=item --mid

A boolean value to indicate whether the 
midpoint between the actual 'start' and 'stop' values
should be used. The default is to use only the 'start' position. 

=item --zero

=item --inter

Source data is in interbase coordinate (0-base) system. Shift the 
start position to base coordinate (1-base) system. Wig files are by 
definition 1-based. This is automatically enabled when converting 
from Bed or BedGraph files. Default is false.

=item --format [0 | 1 | 2 | 3]

Indicate the number of decimal places the score value should
be formatted. Acceptable values include 0, 1, 2, or 3 places.
The default is to not format the score value.

=item --method [mean | median | sum | max]

Define the method used to combine multiple data values at a single 
position. Wig files do not tolerate multiple identical positions.

=item --log

If multiple data values need to be combined at a single identical 
position, indicate whether the data is in log2 space or not. This 
affects the mathematics behind the combination method.

=item --bigwig

=item --bw

Indicate that a binary BigWig file should be generated instead of 
a text wiggle file. A .wig file is first generated, then converted to 
a .bw file, and then the .wig file is removed.

=item --chromof <filename>

When converting to a BigWig file, provide a two-column tab-delimited 
text file containing the chromosome names and their lengths in bp. 
Alternatively, provide a name of a database, below.

=item --db <database>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
or other indexed data file, e.g. Bam or bigWig file, from which chromosome 
length information may be obtained. For more information about using databases, 
see L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. It 
may be supplied from the input file metadata.

=item --bwapp </path/to/wigToBigWig>

Specify the path to the UCSC wigToBigWig or bedGraphToBigWig conversion 
utility. The default is to first check the BioToolBox configuration 
file C<biotoolbox.cfg> for the application path. Failing that, it will 
search the default environment path for the utility. If found, it will 
automatically execute the utility to convert the wig file.

=item --keep

Keep the wig or bedGraph file after converting to a bigWig file. The 
default is to delete the file if the bigWig conversion is successful.

=item --gz

A boolean value to indicate whether the output wiggle 
file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will convert any tab-delimited data text file into a wiggle 
formatted text file. This requires that the file contains not only the 
scores bu also chromosomal coordinates, i.e. chromosome, start, and 
(optionally) stop. The program should automatically detect these 
columns (if appropriately labeled) or they can be specified. An option 
exists to use the midpoint of a region, e.g. microarray probe.

The wig file format is specified by documentation supporting the UCSC 
Genome Browser and detailed here: http://genome.ucsc.edu/goldenPath/help/wiggle.html.
Three formats are supported, 'fixedStep', 'variableStep', and 'bedGraph'. 
The format may be requested or determined empirically from the input file 
metadata. Genomic bin files generated with C<BioToolBox> scripts record 
the window and step values in the metadata, which are used to determine 
the span and step wig values, respectively. The variableStep format 
is otherwise generated by default. The span is, by default, 1 bp.

Wiggle files cannot tolerate multiple datapoints at the same identical 
position, e.g. multiple microarray probes matching a repetitive sequence. 
An option exists to mathematically combine these positions into one value.

Strand is not inherently supported in wig files. If you have stranded data, 
they should be split into separate files. The C<BioToolBox> script 
C<split_data_file.pl> can be used for this purpose.

A binary BigWig file may also be further generated from the  
text wiggle file. The binary format is preferential to the text version 
for a variety of reasons, including fast, random access and no loss in 
data value precision. More information can be found at this location:
http://genome.ucsc.edu/goldenPath/help/bigWig.html. Conversion requires 
BigWig file support, supplied by the external C<wigToBigWig> or 
C<bedGraphToBigWig> utility available from UCSC.

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
