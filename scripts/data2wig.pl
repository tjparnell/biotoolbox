#!/usr/bin/perl

# A script to convert a generic data file into a wig file
# this presumes it has chromosomal coordinates to convert

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(mean median sum max);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	find_column_index
);
# use tim_db_helper has moved down below and is loaded on demand
use tim_file_helper qw(
	open_tim_data_file
	open_to_write_fh
);
use tim_db_helper::config;
eval {
	# check for bigwig file conversion support
	require tim_db_helper::bigwig;
	tim_db_helper::bigwig->import;
};

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
	$step_size,
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
	$gz,
	$help
);


# Command line options
GetOptions( 
	'in=s'      => \$infile, # name of input file
	'out=s'     => \$outfile, # name of output gff file 
	'step=s'    => \$step, # wig step method
	'size=i'    => \$step_size, # wig step size
	'chr=i'     => \$chr_index, # index for the chromosome column
	'start=i'   => \$start_index, # index for the start column
	'stop|end=i'=> \$stop_index, # index for the stop column
	'score=i'   => \$score_index, # index for the score column
	'name=s'    => \$track_name, # name string for the track
	'track!'    => \$use_track, # boolean to include a track line
	'mid!'      => \$midpoint, # boolean to use the midpoint
	'inter!'    => \$interbase, # shift from interbase
	'format=i'  => \$format, # format output to indicated number of places
	'method=s'  => \$method, # method for combining duplicate values
	'log!'      => \$log2, # data is in log2 format
	'bigwig|bw' => \$bigwig, # generate a binary bigwig file
	'db=s'      => \$database, # database for bigwig file generation
	'chromof=s' => \$chromo_file, # name of a chromosome file
	'bwapp=s'   => \$bw_app_path, # path to wigToBigWig utility
	'gz!'       => \$gz, # boolean to compress output file
	'help'      => \$help # request help
);

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
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
unless ($outfile =~ /\.wig$/i) {
	# add extension
	$outfile .= '.wig';
}
my $out_fh = open_to_write_fh($outfile, $gz) or 
	die " unable to open output file '$outfile' for writing!\n";

# write track line
if ($use_track) {
	print {$out_fh} "track type=wiggle_0 name=$track_name\n";
}



### Start the conversion 
print " converting '" . $metadata_ref->{$score_index}{'name'} . "'....\n";
if ($step eq 'fixed') {
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
		# if it is a gff file, it will be the score index
		$chr_index   = 0;
		$start_index = 3;
		$stop_index  = 4;
		$score_index = 5;
	}
	
	# bed
	elsif ( 
		$metadata_ref->{'extension'} =~ /bed/ and
		$metadata_ref->{'number_columns'} >=5
	) {
		# this is the default bed format, not a bedgraph
		$chr_index   = 0;
		$start_index = 1;
		$stop_index  = 2;
		$score_index = 4;
	}
	
	# sgr
	elsif ( $metadata_ref->{'extension'} =~ /sgr/ ) {
		# if it is a sgr file, it will be the score index
		$chr_index   = 0;
		$start_index = 1;
		$score_index = 2;
	}
	
	# non-standard text file or tim data format text file
	else {
		# we will automatically look for the coordinate columns
		$chr_index    = find_column_index($metadata_ref, '^chr|seq|refseq');
		$start_index  = find_column_index($metadata_ref, 'start');
		$stop_index   = find_column_index($metadata_ref, 'stop|end');
		
		# check that we have the required coordinates
		unless (defined $start_index) {
			# no start? how about a midpoint? position?
			$start_index = find_column_index($metadata_ref, 'midpoint|mid|position');
			unless (defined $start_index) {
				# still nothing found, fail
				die " No start, midpoint, or position coordinate column found in the file!\n";
			}
		}
		unless (defined $chr_index) {
			die " No chromosome or sequence ID column found in the file!\n";
		}
		
		# identify the score index
		unless (defined $score_index) {
			# ask the user for help if it wasn't provided
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
	
	# finished with indices
}


sub check_track_name {
	# determine what the track name will be
	
	unless (defined $track_name) {
		if ( 
			$metadata_ref->{'gff'} or
			$metadata_ref->{'extension'} =~ /sgr/ or
			$metadata_ref->{'extension'} =~ /bed/
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
	# check step 
	if ($step eq 'variable') {
		# this is ok, we can work with it
	}
	elsif ($step eq 'fixed') {
		# double check that the data file supports this
		# assign the step size as necessary
		if ( exists $metadata_ref->{$start_index}{'step'} ) {
			$step_size = $metadata_ref->{$start_index}{'step'};
		}
		else {
			warn " no step size indicated for 'fixedStep' wig file! Using 'variableStep'\n";
			$step = 'variable';
		}
	}
	else {
		# attempt to determine automatically
		if ( exists $metadata_ref->{$start_index}{'step'} ) {
			print " Automatically generating 'fixedStep' wig....\n";
			$step = 'fixed';
			$step_size = $metadata_ref->{$start_index}{'step'};
		}
		else {
			print " Automatically generating 'variableStep' wig....\n";
			$step = 'variable';
		}
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
	my $current_chr; # current chromosome
	while (my $line = $in_fh->getline) {
		chomp $line;
		my @data = split /\t/, $line;
		
		# adjust for interbase 0-base coordinates
		if ($interbase) {
			$data[$start_index] += 1;
		}
		
		# write definition line if necessary
		if ($data[$chr_index] ne $current_chr) {
			# new chromosome, new definition line
			
			# need to determine start position first
			my $start = calculate_position(@data);
			
			# print definition line
			print {$out_fh} 'fixedStep chrom=' . $data[$chr_index] . ' start=' .
				$start . ' step=' . $step_size . "\n";
				
			
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
		
		# write fixed data line
		print {$out_fh} "$score\n";
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
			# new chromosome, new definition line
			print {$out_fh} 'variableStep chrom=' . $data[$chr_index] . "\n";
			
			# reset the current chromosome
			$current_chr = $data[$chr_index];
			$previous_pos = undef;
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
				"chromosome $current_chr, compare $position versus $previous_pos\n";
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
	
	# check that bigwig conversion is supported
	unless (exists &wig_to_bigwig_conversion) {
		warn "\n  Support for converting to bigwig format is not available\n" . 
			"  Please convert manually. See documentation for more info\n";
		print " finished\n";
		exit;
	}
	
	# open database connection if necessary
	my $db;
	if ($database) {
		eval {
			use tim_db_helper qw(open_db_connection);
		};
		if ($@) {
			warn " unable to load tim_db_helper! Is BioPerl installed?\n";
		}
		else {
			$db = open_db_connection($database);
		}
	}
	
	# find wigToBigWig utility
	unless ($bw_app_path) {
		# check for an entry in the configuration file
		$bw_app_path = $TIM_CONFIG->param('applications.wigToBigWig') || 
			undef;
	}
	unless ($bw_app_path) {
		# next check the system path
		$bw_app_path = `which wigToBigWig` || undef;
	}
			
	# perform the conversion
	my $bw_file = wig_to_bigwig_conversion( {
			'wig'       => $outfile,
			'db'        => $db,
			'chromo'    => $chromo_file,
			'bwapppath' => $bw_app_path,
	} );

	
	# confirm
	if ($bw_file) {
		print " bigwig file '$bw_file' generated\n";
		unlink $outfile; # remove the wig file
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
  --step [fixed | variable]
  --size <integer>
  --score <column_index>
  --name <text>
  --(no)track
  --(no)mid
  --inter
  --format [0 | 1 | 2 | 3]
  --method [mean | median | sum | max]
  --(no)log
  --bigwig|--bw
  --chromof <filename>
  --db <database>
  --(no)gz
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a data file, with or without the --in flag. 
The file may be any tab-delimited text file (preferably in the tim data 
format as described in C<tim_file_helper.pm>), GFF, SGR, or BED file. 
Recognizeable genome coordinate columns should be present, including 
chromosome, start, and stop. Data files collected using the 'genome' 
windows feature are ideal. The file may be compressed with gzip.

=item --out <filename>

Optionally specify the name of of the output file. The track name is 
used as default. The '.wig' extension is automatically added if required.

=item --step [fixed | variable]

The type of step progression for the wig file. Two wig formats are available:
'fixedStep' where data points are positioned at equal distances along the 
chromosome, and 'variableStep' where data points are not equally spaced 
along the chromosome. The 'fixedStep' wig file has one column of data 
values (score), while the 'variableStep' wig file has two columns
(position and score). If the option is not defined, then the format is 
automatically determined from the metadata of the file.

=item --size <integer>

Optionally define the step size in bp for 'fixedStep' wig file. This 
value is automatically determined from the table's metadata, if available. 
If the --step option is explicitly defined as 'fixed', then the step size 
may also be explicitly defined. If this value is not explicitly
defined or automatically determined, the variableStep format is used by
default.

=item --score <column_index>

Indicate the column index (0-based) of the dataset in the data table 
to be used for the score. If a GFF file is used as input, the score column is 
automatically selected. If not defined as an option, then the program will
interactively ask the user for the column index from a list of available
columns.

=item --name <text>

The name of the track defined in the wig file. The default is to use 
the name of the chosen score column, or, if the input file is a GFF file, 
the base name of the input file. 

=item --(no)track

Do (not) include the track line at the beginning of the wig file. Wig 
files normally require a track line, but if you will be converting to 
the binary bigwig format, the converter requires no track line. Why it 
can't simply ignore the line is beyond me. This option is automatically 
set to false when the --bigwig option is enabled.

=item --(no)mid

A boolean value to indicate whether the 
midpoint between the actual 'start' and 'stop' values
should be used. The default is to use only the 'start' position. 

=item --inter

Source data is in interbase coordinate (0-base) system. Shift the 
start position to base coordinate (1-base) system. Wig files are by 
definition 1-based. Default is false.

=item --format [0 | 1 | 2 | 3]

Indicate the number of decimal places the score value should
be formatted. Acceptable values include 0, 1, 2, or 3 places.
The default is to not format the score value.

=item --method [mean | median | sum | max]

Define the method used to combine multiple data values at a single 
position. Wig files do not tolerate multiple identical positions.

=item --(no)log

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

Specify the database from which chromosome lengths can be derived when 
generating a bigwig file. This option is only required when generating 
bigwig files. It may also be supplied from the metadata in the source 
data file.

=item --(no)gz

A boolean value to indicate whether the output wiggle 
file should be compressed with gzip.

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
Genome Browser and detailed at this location:
http://genome.ucsc.edu/goldenPath/help/wiggle.html
Two formats are supported, 'fixedStep' and 'variableStep'. 

Wiggle files cannot tolerate multiple datapoints at the same identical 
position, e.g. multiple microarray probes matching a repetitive sequence. 
An option exists to mathematically combine these positions into one value.

A binary BigWig file may also be further generated from the  
text wiggle file. The binary format is preferential to the text version 
for a variety of reasons, including fast, random access and no loss in 
data value precision. More information can be found at this location:
http://genome.ucsc.edu/goldenPath/help/bigWig.html. Conversion requires 
BigWig file support, supplied by the biotoolbox module 
C<tim_db_helper::bigwig>. 


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











