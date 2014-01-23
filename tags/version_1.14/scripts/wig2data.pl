#!/usr/bin/env perl

# documentation at end of file 

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::data_helper qw(generate_tim_data_structure);
use Bio::ToolBox::file_helper qw(
	write_tim_data_file
	open_to_read_fh
	open_to_write_fh
	convert_genome_data_2_gff_data
);
my $VERSION = '1.14';

print "\n This program will convert wiggle files to a tabbed text file\n\n";

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
	$gff,
	$type,
	$source,
	$places,
	$midpoint,
	$version,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the solexa data file
	'out=s'     => \$outfile, # name of output file 
	'gff!'      => \$gff, # write a gff file
	'type=s'    => \$type, # the name of the data, goes into the type field of GFF
	'source=s'  => \$source, # the source of the data, goes into the source field of GFF
	'format=i'  => \$places, # indicate number of decimal places 
	'midpoint!' => \$midpoint, # use midpoint instead of start and stop
	'version=i' => \$version, # the gff version
	'gz!'       => \$gz, # compress output
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
	print " Biotoolbox script wig2data.pl, version $VERSION\n\n";
	exit;
}




### Check for requirements
unless ($infile) {
	$infile = shift @ARGV or
		die "  OOPS! No source data file specified! \n use $0 --help\n";
}
unless ($source) {
	# default is none
	$source = '.';
}
unless (defined $gz) {
	$gz = 0;
}
unless ($version) {
	$version = 3;
}



### Open the file 
print " Converting file '$infile'...\n";
my $in_fh = open_to_read_fh($infile) or 
	die " unable to open file!\n";




### Do the conversion

# prepare output data structure

# initialize variables
my $out_data_ref = initialize_data_structure();
my (
	# reusuable variables
	$refseq,
	$fixstart,
	$step,
	$span,
	$out_fh,
); 
my $count = 0;

# main loop
while (my $line = $in_fh->getline) {
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
		# but we can use the track name as the type
		foreach (@data) {
			if (/name=(.+)/) {
				$type = $1;
				$type =~ s/ /_/g; # convert spaces to underscores
				$type =~ s/"//g; # remove quotation marks
				last;
			}
		}
		unless (defined $type) {
			# default is to use the base filename
			$type = $infile;
			$type =~ s/\.wig(?:\.gz)$//i;
		}
	}
	
	# a variable step definition line
	elsif ($data[0] =~ /^variablestep$/i) { 
		foreach (@data) {
			if (/chrom=(\w+)/) {$refseq = $1}
			if (/span=(\w+)/) {$span = $1}
		}
		next;
		
	} 
	
	# a fixed step definition line
	elsif ($data[0] =~ /^fixedstep$/i) { 
		foreach (@data) {
			if (/chrom=(\w+)/) {$refseq = $1}
			if (/span=(\w+)/) {$span = $1}
			if (/start=(\w+)/) {$fixstart = $1}
			if (/step=(\w+)/) {$step = $1}
		}
		next;
	} 
	
	# a BED data line
	elsif (scalar @data == 4) {
		$refseq = $data[0];
		$start = $data[1] + 1; # the BED line alone uses 0-based indexing
		$stop = $data[2];
		$score = $data[3];
	} 
	
	# a variable step data line
	elsif (scalar @data == 2) { 
		unless ($refseq) { 
			die "Bad formatting! variable step data but chromosome not defined!\n";
		}
		$start = $data[0];
		if ($span) {
			$stop = $start + $span;
		} 
		else {
			$stop = $start;
		}
		$score = $data[1];
	} 
	
	# a fixed step data line
	elsif (scalar @data == 1) { 
		unless ($refseq) { 
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
		if ($span) {
			$stop = $start + $span;
		} 
		else {
			$stop = $start;
		}
		$score = $data[0];
	}
	
	# format the score value
	my $fscore;
	if (defined $places) {
		if ($places == 0) {
			$fscore = sprintf "%.0f", $score;
		} 
		elsif ($places == 1) {
			$fscore = sprintf "%.1f", $score;
		} 
		elsif ($places == 2) {
			$fscore = sprintf "%.2f", $score;
		} 
		elsif ($places == 3) {
			$fscore = sprintf "%.3f", $score;
		}
	} 
	else {
		$fscore = $score;
	}
	
	# add the gff data to the data structure
	push @{ $out_data_ref->{'data_table'} }, [
		$refseq,
		$start,
		$stop,
		$fscore
	];
	$count++;
	
	# temporarily write output
	if ($count == 50000) {
		
		# progressively write out the converted data
		write_progressive_data();
		
		# regenerate the output data table
		$out_data_ref = initialize_data_structure();
		
		# reset count
		$count = 0;
	}
	
}


### Write final output
write_progressive_data();




### Finish
$in_fh->close;
$out_fh->close;
print " wrote file '$outfile'\n";






########################   Subroutines   ###################################


sub initialize_data_structure {
	my $data = generate_tim_data_structure(
		'wig_data',
		qw(
			Chromo
			Start
			Stop
			Score
		)
	) or die " unable to generate tim data structure!\n";
	
	# add metadata
	$data->{3}{'original_file'} = $infile;
	if (defined $places) {
		$out_data_ref->{3}{'formatted'} = $places;
	}
	
	# finished
	return $data;
}




sub write_progressive_data {
	# a subroutine to progressively write out the converted data
	
	# update last line
	$out_data_ref->{'last_row'} = scalar @{ $out_data_ref->{'data_table'} } -1;
	
	# convert to gff if requested
	if ($gff) {
		convert_genome_data_2_gff_data(
			'data'     => $out_data_ref,
			'score'    => 3,
			'source'   => $source,
			'type'     => $type,
			'midpoint' => $midpoint,
			'version'  => $version,
		) or die " Unable to convert to GFF format!\n";
	}
	
	# check for filename
	unless ($outfile) {
		# default is to use the type
		$outfile = $type;
	}
	
	# check for file handle
	if (defined $out_fh) {
		# the output file has been opened and partially written
		# we now only need to write the data table portion and not the 
		# metadata
		
		for my $row (1 .. $out_data_ref->{'last_row'}) {
			print {$out_fh} join(
					"\t", @{ $out_data_ref->{'data_table'}->[$row] }
				), "\n";
				 
		}
	}
	else {
		# we will need to open the output file to write if it's not 
		# opened yet
		
		# rather than generating new code for writing the gff file,
		# we will simply use the write_tim_data_file sub
		my $new_outfile = write_tim_data_file(
			'data'      => $out_data_ref,
			'filename'  => $outfile,
			'gz'        => $gz,
		);
		if ($new_outfile) {
			# success
			# reassign the name
			$outfile = $new_outfile;
		}
		else {
			die " unable to write output file!\n";
		}
		
		# but now we will have to reopen the file for appended writing
		$out_fh = open_to_write_fh($outfile, $gz, 1);
	}

}




__END__

=head1 NAME

wig2data.pl

A script to convert a text wiggle file to a tab-delimited text file. 

=head1 SYNOPSIS

wig2data.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --out <filename> 
  --gff
  --type <text>
  --source <text>
  --format [0,1,2,3]
  --midpoint
  --version [2,3]
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a wig file. The file may be compressed with gzip.

=item --out <filename>

Specify the output filename. By default it uses the GFF type as the 
basename.

=item --gff

Indicate whether the output file should be in GFF format. If false, a 
standard tim data tab delimited text file will be written with four 
columns: chromosome, start, stop, and score. The default value is false.

=item --type <text>

Specify the text string to be used as the GFF feature 'type' or 
'method' (the 3rd column). By default it uses the name specified in the 
track line; otherwise, it uses the basename of the input wig file.

=item --source <text>

Specify the text string to be used as the GFF feature 'source' 
(the 2nd column). The default value is none.

=item --format [0,1,2,3]

Specify the number of decimal places to which the wig file will be 
formatted. Default is no formatting.

=item --midpoint

Specify whether (or not) a midpoint position should be calculated 
between the start and stop positions and be used in the output GFF 
file. This only pertains to BED style wig files where both the 
start and stop positions are reported and stepped wig files where 
a span value is specified. The default value is false.

=item --version [2,3]

Specify the GFF version. The default is version 3.

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will convert a wiggle data file into a tabbed delimited text 
data file. The data file will have four columns: chromosome, start, stop, 
and score. Alternatively, a GFF file may be written, in which case the 
GFF source and type values should be specified.

Wiggle files are used with the UCSC Genome Browser and can have multiple
formats, including BED (also referred to as bedgraph, variable step, and
fixed step. This program can convert all three formats. Improperly
formatted wig files may cause the program to die. More information about
wig files may be obtained from here:
http://genome.ucsc.edu/goldenPath/help/wiggle.html

The start position will be the coordinate listed in the wig file. If a span
value is indicated in the wiggle file, the GFF stop will equal start plus
span; otherwise stop will equal the start value. The score value may be
formatted to the indicated number of decimal places.
 
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
