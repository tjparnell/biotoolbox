#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(mean median sum max);
use Bio::ToolBox::Data;
use Bio::ToolBox::utility qw(ask_user_for_index);
use Bio::ToolBox::big_helper qw(wig_to_bigwig_conversion);
my $VERSION =  '1.33';

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
	$use_track = $bigwig ? 0 : 1;
}



### Load input file
my $Input = Bio::ToolBox::Data->new(file => $infile) or
	die "Unable to open file '$infile'!\n";



### Check and/or ask for specific options

check_indices();

check_track_name();

check_step();

set_bigwig_options() if $bigwig;

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
my $out_fh = Bio::ToolBox::Data->open_to_write_fh($outfile, $gz) or 
	die " unable to open output file '$outfile' for writing!\n";

# write track line
if ($use_track) {
	print {$out_fh} "track type=wiggle_0 name=$track_name\n";
}



### Start the conversion 
printf " converting '%s'....\n", $Input->name($score_index);
# sort by genomic coordinates just in case
$Input->gsort_data if $Input->feature_type eq 'coordinate';
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

sub check_indices {
	
	# check coordinates
	unless (defined $chr_index or defined $Input->chromo_column) {
		$chr_index = ask_user_for_index($Input, 
			" Enter the index for the chromosome column  ");
		unless (defined $chr_index) {
			die " No identifiable chromosome column index!\n";
		}
	}
	unless (defined $start_index or defined $Input->start_column) {
		$start_index = ask_user_for_index($Input, 
			" Enter the index for the start or position column  ");
		unless (defined $start_index) {
			die " No identifiable start column index!\n";
		}
	}
	# stop column is optional	
	
	# score
	unless (defined $score_index) {
		# first look for a generic score index
		$score_index = ask_user_for_index($Input, 
			" Enter the index for the score column  ");
		unless (defined $score_index) {
			die " No identifiable score column index!\n";
		}
	}
}


sub check_track_name {
	# determine what the track name will be
	$track_name = $Input->name($score_index);
	if ($track_name =~/^score$/i) {
		# some sort of bed or gff file standard score column
		$track_name = $Input->basename;
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
		if (defined $Input->metadata($start_index, 'step') ) {
			if (defined $step_size) {
				if ($step_size != $Input->metadata($start_index, 'step')) {
					die " Requested step size $step_size does not match" .
						" metadata step size!!!\n";
				}
			}
			else {
				# define it from the metadata
				$step_size = $Input->metadata($start_index, 'step');
			}
		}
		elsif (!defined $step_size) {
			warn " Fixed step size not defined by user or metadata! Using 'variableStep'\n";
			$step = 'variable';
		}
	}
	else {
		# attempt to determine automatically
		if (defined $Input->metadata($start_index, 'step') ) {
			# set step size
			$step = 'fixed';
			$step_size = $Input->metadata($start_index, 'step');
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
		if (defined $Input->metadata($start_index, 'win') and
			$Input->metadata($start_index, 'win') != $span
		) {
			# the requested span and metadata window size do not match
			die " Requested span size $span does not match metadata window size!!!\n";
		}
	}
	else {
		# attempt to determine automatically 
		if (defined $Input->metadata($start_index, 'win') and not $midpoint) {
			# set the span equal to the window size
			$span = $Input->metadata($start_index, 'win');
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


sub set_bigwig_options {
	# if we're generating bigwig file, no track is needed
	$use_track = 0;
	
	# force no compression
	$gz = 0;
	
	# check that we have a source for chromosome info
	unless ($database or $chromo_file) {
		$database = $Input->database or
		die " No database name or chromosome file provided for generating bigwig file!\n";
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
	else {
		# default is mean
		return \&mean;
	}
}


sub convert_to_fixedStep {
	
	# keep track of current chromosome name and length
	my $current_chr; # current chromosome
	
	# walk through the data file
	my $iterator = $Input->row_stream;
	while (my $row = $iterator->next_row) {
		# coordinates
		my $chromosome = defined $chr_index ? $row->value($chr_index) : $row->seq_id;
		my $start = calculate_position($row);
		next if $start <= 0; # skip negative or zero coordinates
		
		# write definition line if necessary
		if ($chromosome ne $current_chr) {
			# new chromosome, new definition line
			my $def = sprintf "fixedStep chrom=%s start=%d step=%d span=%d\n",
				 $chromosome, $start, $step_size, $span;
			$out_fh->print($def);
			
			# reset the current chromosome
			$current_chr = $chromosome;
		}
		
		# adjust score formatting as requested
		my $score = $row->value($score_index);
		$score = 0 if $score eq '.';
		if (defined $format) {
			# format if requested
			$score = format_score($score);
		}
		
		# write fixed data line
		$out_fh->print("$score\n");
	}
}


sub convert_to_variableStep {
	my $current_chr; # current chromosome
	my $previous_pos = 0; # previous position to avoid duplicates in wig file
	my @scores; # reusable array for putting multiple data points in
	my $iterator = $Input->row_stream;
	while (my $row = $iterator->next_row) {
		# coordinates
		my $chromosome = defined $chr_index ? $row->value($chr_index) : $row->seq_id;
		my $start = calculate_position($row);
		next if $start <= 0; # skip negative or zero coordinates
		
		# write definition line if necessary
		if ($chromosome ne $current_chr) {
			# first check and write orphan scores
			# this might happen if there was only one score on the entire chr
			if (@scores) {
				if (scalar @scores == 1) {
					# print the one score
					$out_fh->print("$previous_pos $scores[0]\n");
				}
				else {
					# more than one score, combine them
					my $score = &{$method_sub}(@scores);
					$out_fh->print("$previous_pos $score\n");
				}
				@scores = ();
			}
			
			# print new definition line and reset for next
			$out_fh->print("variableStep chrom=$chromosome span=$span\n");
			$current_chr = $chromosome;
			$previous_pos = $start;
		}
		
		# collect the score
		my $score = $row->value($score_index);
		next if ($score eq '' or $score eq '.'); # skip null values
		if (defined $format) {
			# format if requested
			$score = format_score($score);
		}
		
		# check for duplicate positions and write appropriately
		if ($start == $previous_pos) {
			# same position, add to the score list
			push @scores, $score;
		}
		else {
			# we have moved on to the next position
			# now print the previous scores
			if (scalar @scores == 1) {
				# print the one score
				$out_fh->print("$previous_pos $scores[0]\n");
			}
			elsif (scalar @scores > 1) {
				# more than one score
				my $score = &{$method_sub}(@scores);
				$out_fh->print("$previous_pos $score\n");
			}
			
			if ($start < $previous_pos) {
				# print warning that output wig will not be sorted
				warn " warning! output wig will not be sorted correctly! chromosome $chromosome: $start > $previous_pos\n"; 
			}
			
			# reset for next
			$previous_pos = $start;
			@scores = ($score);
		}
	}
}


sub convert_to_bedgraph {
	
	# variables to check for overlap
	my $current_chr; # current chromosome
	my $previous_pos; # previous position to avoid overlap
	my $iterator = $Input->row_stream;
	while (my $row = $iterator->next_row) {
		# coordinates
		my $chromosome = defined $chr_index ? $row->value($chr_index) : $row->seq_id;
		my $start = defined $start_index ? $row->value($start_index) : $row->start;
		my $stop  = defined $stop_index ? $row->value($stop_index) : $row->stop;
		
		# adjust start position
		unless ($interbase) {
			$start--;
		}
		
		# check coordinates
		if (defined $previous_pos and defined $current_chr) {
			# check if on the same chromosome
			if ($current_chr eq $chromosome) {
				# check for overlap
				if ($start < $previous_pos) {
					die " There are overlapping intervals or the file is not sorted by" .
						" coordinates!\n Compare $chromosome:$start" . 
						" with previous stop position $previous_pos\n";
				}
				# otherwise it is ok
				$previous_pos = $stop;
			}
			else {
				# new chromosome
				$current_chr = $chromosome;
				$previous_pos = $stop;
			}
		}
		else {
			$current_chr = $chromosome;
			$previous_pos = $stop;
		}
		
		# collect the score
		my $score = $row->value($score_index);
		next if $score eq '.'; # skip null values
		if (defined $format) {
			# format if requested
			$score = format_score($score);
		}
		
		# write the feature line
		$out_fh->print(join("\t", $chromosome, $start, $stop, $score), "\n");
	}
}


sub calculate_position {
	my $row = shift;
	my $start;
	if (defined $start_index) {
		$start = $row->value($start_index);
		$start += 1 if $interbase;
	}
	else {
		$start = $row->start; # this should handle interbase automatically
	}
	if ($midpoint) {
		# calculate midpoint and return
		my $end = defined $stop_index ? $row->value($stop_index) : $row->end;
		$end ||= $start; # in case no end was defined
		return $start == $end ? $start : int( ( ($start + $end) / 2) + 0.5) ;
	}
	else {
		return $start;
	}
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
definition 1-based. This is automatically handled for most input  
files. Default is false.

=item --format [0 | 1 | 2 | 3]

Indicate the number of decimal places the score value should
be formatted. Acceptable values include 0, 1, 2, or 3 places.
The default is to not format the score value.

=item --method [mean | median | sum | max]

Define the method used to combine multiple data values at a single 
position. Wig files do not tolerate multiple identical positions. 
Default is mean.

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
it under the terms of the Artistic License 2.0.  
