#!/usr/bin/perl

# documentation at end of file

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(mean stddevp median min max);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	write_tim_data_file
	convert_and_write_to_gff_file
);
eval {
	# check for bam support
	require Bio::ToolBox::db_helper::bam;
	Bio::ToolBox::db_helper::bam->import;
};
my $VERSION = '1.15';


print "\n A script to get exact nucleosome fragment sizes from a Bam file\n\n";

### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Command line options
my (
	$infile, 
	$bamfile, 
	$minsize, 
	$maxsize, 
	$winsize,
	$AT_ends,
	$gff,
	$type,
	$source,
	$outfile,
	$help,
	$print_version,
);
GetOptions( 
	'in=s'       => \$infile, # the input nucleosome data file
	'bam=s'      => \$bamfile, # name of bam file 
	'min=i'      => \$minsize, # the minimum cutoff size for paired-read segments
	'max=i'      => \$maxsize, # the maximum cutoff size for paired-read segments
	'win=i'      => \$winsize, # the window size for looking for nucleosome mids
	'at!'        => \$AT_ends, # discard non-AT ends
	'gff!'       => \$gff, # boolean to write a GFF
	'type=s'     => \$type, # the GFF type
	'source=s'   => \$source, # the GFF source
	'out=s'      => \$outfile, # output file name
	'help'       => \$help, # request help
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
	print " Biotoolbox script get_actual_nuc_sizes.pl, version $VERSION\n\n";
	exit;
}


### Check for required values
unless ($infile and $bamfile) {
	die "  Both input nucleosome data and BAM files must be specified!\n";
}




### Open files

# Input data file
print " opening files....\n";
my $main_data_ref = load_tim_data_file($infile) or 
	die "  Unable to load data file!\n";
unless (
	$main_data_ref->{'feature'} eq 'nucleosome' and
	$main_data_ref->{'program'} =~ /map_nucleosomes/
) {
	# we require output from the map_nucleosomes.pl program
	die "  The data file '$infile' doesn't have nucleosome features!";
}

# BAM file
unless (exists &open_bam_db) {
	die " unable to load Bam file support! Is Bio::DB::Sam installed?\n"; 
}
my $sam = open_bam_db($infile) or die " unable to open bam file '$infile'!\n";




### Collecting nucleosome sizes

# set up new metadata and get indices
my ($count_i, $mean_i, $std_i) = generate_new_metadata($main_data_ref);

# collect nucleosome data
print " collecting nucleosome sizes....\n";
collect_nuc_sizes(
	$main_data_ref, $count_i, $mean_i, $std_i);

# generate final statistics
print_final_stats($main_data_ref, $mean_i);




### Write output

# write data file
unless ($outfile) {
	# overwrite the input file
	$outfile = $infile;
}
my $success = write_tim_data_file(
	'data'      => $main_data_ref,
	'filename'  => $outfile,
);
if ($success) {
	print " Wrote data file '$success'\n";
}
else {
	print " Unable to write data file!\n";
}

# write GFF file
if ($gff) {
	# write if requested
	
	# filename
	$outfile =~ s/\.txt/.gff/;
	
	# check required values
	unless ($type) {
		# default GFF3 compatible SO term
		$type = 'histone_binding_site';
	}
	unless ($source) {
		# set default source, the name of this program
		$source = 'map_nucleosomes.pl';
	}
	
	# write
	$success = convert_and_write_to_gff_file(
		'data'     => $main_data_ref,
		'filename' => $outfile,
		'version'  => 3,
		'score'    => $count_i,
		'name'     => 4,
		'type'     => $type,
		'source'   => $source,
		'tag'      => [6, $std_i], # fuzziness, stdev of length
	);
	if ($success) {
		print " Wrote GFF file '$success'\n";
	}
	else {
		print " Unable to write GFF file!\n";
	}
	
}

### The End





############################# Subroutines #####################################

sub generate_new_metadata {
	my $data_ref = shift;
	
	# count column
	my $count_i = $data_ref->{'number_columns'};
	$data_ref->{$count_i} = {
		'index'     => $count_i,
		'name'      => 'nuc_count',
		'bamfile'   => $bamfile,
	};
	$data_ref->{'data_table'}->[0][$count_i] = 'nuc_count';
	$data_ref->{'number_columns'} += 1;
	
	# mean column
	my $mean_i = $data_ref->{'number_columns'};
	$data_ref->{$mean_i} = {
		'index'     => $mean_i,
		'name'      => 'nuc_size_mean',
		'bamfile'   => $bamfile,
	};
	$data_ref->{'data_table'}->[0][$mean_i] = 'nuc_size_mean';
	$data_ref->{'number_columns'} += 1;
	
	# std column
	my $std_i = $data_ref->{'number_columns'};
	$data_ref->{$std_i} = {
		'index'     => $std_i,
		'name'      => 'nuc_size_std',
		'bamfile'   => $bamfile,
	};
	$data_ref->{'data_table'}->[0][$std_i] = 'nuc_size_std';
	$data_ref->{'number_columns'} += 1;
		
	
	# add extra metadata
	if ($minsize) {
		$data_ref->{$mean_i}{'min_size'} = $minsize;
	}
	if ($maxsize) {
		$data_ref->{$mean_i}{'max_size'} = $maxsize;
	}
	if ($AT_ends) {
		$data_ref->{$count_i}{'AT_ends'} = 'skipped';
	}
	if ($winsize) {
		$data_ref->{$count_i}{'window_size'} = $winsize;
	}
	else {
		$data_ref->{$count_i}{'window_size'} = 'fuzziness_score';
	}
	
	return ($count_i, $mean_i, $std_i);
}



sub collect_nuc_sizes {
	my ($data_ref, $count_i, $mean_i, $std_i) = @_;
	my $table = $data_ref->{'data_table'};
	
	# walk through the nucleosomes
	for my $row (1 .. $data_ref->{'last_row'} ) {
		
		# we are assuming we are working with the output from map_nucleosomes.pl
		# therefore the indices are hard coded
		# this will have to change to make it general purpose
		
		# determine the search window endpoints
		# this is centered around the determined midpoint
		# this may be expressly defined or we'll use the fuzziness value
		my $midpoint = $table->[$row][3];
		my ($start, $end);
		if ($winsize) {
			# window size explicitly defined
			$start = $midpoint - $winsize;
			$end = $midpoint + $winsize;
		}
		else {
			# we'll use the fuzziness value
			$start = $midpoint - $table->[$row][6];
			$end = $midpoint + $table->[$row][6];
		}
		
		# determine the search parameters
			# we can't actually use the start and end points we just defined
			# otherwise we won't find anything
			# limitation of bam module in that both paired reads must be within  
			# the window in order to be found, hence we specify a nice big fat 
			# window to search
		my $search_start = $start > 150 ? $start - 150 : 1, # no negative starts
		my $chr_length = $sam->length( $table->[$row][0] );
		my $search_end = ($end + 150) < $chr_length ? ($end + 150) : $chr_length;
		
		# collect the features
		my @pairs = $sam->features(
			-seq_id     => $table->[$row][0],
			-start      => $search_start,
			-end        => $search_end,
			-type       => 'read_pair',
		);
		
		# walk through the features
		my @sizes;
		foreach my $pair (@pairs) {
			
			# get individual reads
			my ($left, $right) = $pair->get_SeqFeatures;
			unless (defined $left and defined $right) {
				# one or both alignments are not defined
				# does this mean one of the alignments is unaligned?
				next;
			}
			unless ($left->proper_pair) {
				# the pair of reads are not properly mapped
				# is this redundant?
				next;
			}
			
			# get and check size
			my $size = $pair->length;
			if (defined $minsize and $size < $minsize) {
				next;
			}
			if (defined $maxsize and $size > $maxsize) {
				next;
			}
			
			# check AT ends if requested
			if ($AT_ends) {
				# both ends must be either A or T to continue, otherwise report
				my $leftseq = $left->query->dna; 
				my $rightseq = $right->query->dna; 
					# this should be reverse-complemented
				unless ($leftseq =~ /^[aAtT]/ and $rightseq =~ /[aAtT]$/) {
					# one of the ends is not correct
					next;
				}
			}
			
			# check whether inside the window
			my $pair_midpoint = int( ($pair->start + $pair->end) / 2);
			if ($pair_midpoint >= $start and $pair_midpoint <= $end) {
				# the midpoint is within our window
				# include this size
				push @sizes, $size;
			}
		}
		
		if (@sizes) {
			# calculate the statistics
			$table->[$row][$count_i] = scalar @sizes;
			my $mean = sprintf "%.0f", mean(@sizes);
			$table->[$row][$mean_i] = $mean;
			$table->[$row][$std_i] = sprintf "%.0f", stddevp(@sizes);
			
			# determine new coordinates
			# we'll use the mean size of the fragment, centered around midpoint
			# we're updating the original start, stop coordinates
			# these are hard coded since we're assuming output from map_nucleosomes
			$table->[$row][1] = $midpoint - int($mean/2);
			$table->[$row][2] = $midpoint + int($mean/2);
		}
		else {
			warn " no pairs found for nucleosome $table->[$row][4]!?\n";
			$table->[$row][$count_i] = 0;
			$table->[$row][$mean_i] = '.';
			$table->[$row][$std_i] = '.';
			# we'll leave the original start/stop coordinates the same
		}
	}
}



sub print_final_stats {
	
	my ($data_ref, $mean_i) = @_;
	
	# collect all the determined sizes
	# walk through the nucleosomes
	my @sizes;
	for my $row (1 .. $data_ref->{'last_row'} ) {
		my $value = $data_ref->{'data_table'}->[$row][$mean_i];
		unless ($value eq '.' or $value == 0) {
			# no null or 0 values
			push @sizes, $value;
		}
	}
	
	# print stats
	printf " The mean nucleosome size is %.0f", mean(@sizes);
	printf " +/- %.0f bp\n", stddevp(@sizes);
	printf " The median is %.0f bp\n", median(@sizes);
	printf " The range is ", min(@sizes), "..", max(@sizes), " bp\n";
}




__END__

=head1 NAME

get_actual_nuc_sizes.pl

A script to pull out actual nucleosome fragments and enumerate their sizes.

=head1 SYNOPSIS

  get_actual_nuc_sizes.pl --in <file1.txt> --bam <file2.bam> [--options]
  
  Options:
  --in <filename>
  --bam <filename.bam>
  --min <integer>
  --max <integer>
  --win <integer>
  --at
  --gff
  --type <gff_type>
  --source <gff_source>
  --out <filename>
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a nucleosome data file generated by the script 
B<map_nucleosomes.pl>. Other data files will likely not work.

=item --bam <filename>

Specify the file name of a binary BAM file containing the original 
paired-end sequence alignment pairs representing nucleosome fragments. 
The file should be sorted and indexed.

=item --min <integer>

Optionally specify the minimum size of fragment to include when determing 
fragment lengths.

=item --max <integer>

Optionally specify the maximum size of fragment to include when determing 
fragment lengths.

=item --win <integer>

Optionally specify the window size when searching for corresponding 
sequence alignment pairs. The window is determined as the mapped nucleosome 
midpoint +/- the specified value. The default value is the 
calculated fuzziness value determined when mapping the nucleosome.

=item --at

Boolean option to indicate that only fragments whose paired sequence reads 
end in a [AT] nucleotide should be included in the output GFF file. 
Micrococcal nuclease (MNase) cuts (almost?) exclusively at AT dinucleotides; 
this option ensures that the fragment is more likely derived from a MNase 
cut. Default is false where all fragments are taken.

=item --gff

Indicate whether a GFF file should be written in addition to the standard 
text data file. The GFF file version is 3. Default is false (no GFF written).

=item --type <gff_type>

Provide the text to be used as the GFF type (or method) used in 
writing the GFF file. The default value is the Sequence Ontology term 
'histone_binding_site'.

=item --source <gff_source>

Provide the text to be used as the GFF source used in writing the 
GFF file. The default value is the name of this program.

=item --out <filename>

Provide a new output file name. By default it overwrites the input file.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will determine actual nucleosome fragment sizes based on the 
original paired-end sequence alignments. It searches a BAM file for all 
aligned read-pairs whose midpoints are within a specific window centered 
around the midpoint of a mapped nucleosome. The program accepts as input 
data the mapped nucleosomes identified using the script B<map_nucleosomes.pl>. 
The window size may be specified explicitly, or by default it uses the 
fuzziness value identified in the mapping program. Once the read-pairs are 
identified, then the mean length of all fragments is determined. The nucleosme 
start and stop coordinates are then updated to accurately reflect the real 
length. The midpoint coordinate is not updated, nor is it checked for 
accuracy; it is assumed to be accurately mapped.

In addition to updating the start and stop coordinates, three additional 
columns of data are appended to the data table. These include the count of 
sequence read-pair fragments, and the standard deviation of the fragment 
lengths. In addition to writing a new data file, it can optionally write 
a GFF3 file.

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
