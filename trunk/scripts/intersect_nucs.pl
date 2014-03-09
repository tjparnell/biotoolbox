#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Statistics::Lite qw(min max mean stddevp);
use Pod::Usage;
use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure
	index_data_table
	find_column_index
	format_with_commas
);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	write_tim_data_file
	convert_and_write_to_gff_file
);
my $VERSION = '1.15';

print "\n This script will intersect two lists of nucleosomes\n\n";



### Quick help
unless (@ARGV) { # when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Get command line options and initialize values

# Initialize values
my (
	$infile1,
	$infile2,
	$outfile,
	$set_strand,
	$gff,
	$type,
	$source,
	$gz,
	$help,
	$print_version,
); 

# Command line options
GetOptions( 
	'in1=s'      => \$infile1, # input file one
	'in2=s'      => \$infile2, # input file two
	'out=s'      => \$outfile, # output filename
	'force_strand|set_strand' => \$set_strand, # artificially enforce a strand for target
				# force_strand is preferred option, but respect the old option
	'gff!'       => \$gff, # output gff file
	'type=s'     => \$type, # the gff type
	'source=s'   => \$source, # the gff source
	'gz!'        => \$gz, # gzip status
	'help'       => \$help, # help
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
	print " Biotoolbox script intersect_nucs.pl, version $VERSION\n\n";
	exit;
}



# Check for required values
unless ($infile1 and $infile2) {
	$infile1 = shift @ARGV or die " no input file specified!\n";
	$infile2 = shift @ARGV or die " two input files must be specified!\n";
}
unless (defined $gz) {$gz = 0}




### Load input files
my $data1 = load_tim_data_file($infile1) 
	or die " no data loaded from file '$infile1'!\n";
print " Loaded $data1->{last_row} features from file '$infile1'\n"; 

my $data2 = load_tim_data_file($infile2) 
	or die " no data loaded from file '$infile2'!\n";
print " Loaded $data2->{last_row} features from file '$infile2'\n\n"; 




### Intersection

# we will assign the datasets as target and reference based on the number 
# of features, reference always has more

# intersect the two lists
my $output;
if ($data1->{last_row} <= $data2->{last_row}) {
	
	# data1 is target, data2 is reference
	print " Intersecting target list from '$data1->{basename}'\n" .
		"    with\n reference list from '$data2->{basename}'...\n";
	$output = intersect_nucs($data1, $data2);
	
	# generate report
	print_statistics($output, $data1, $data2);
}

else {
	
	# data2 is target, data1 is reference
	print " Intersecting target list from '$data2->{basename}'\n" .
		"    with\n reference list from '$data1->{basename}'...\n";
	$output = intersect_nucs($data2, $data1);
	
	# generate report
	print_statistics($output, $data2, $data1);
}





### Output

# write data file
unless ($outfile) {
	# make up a filename
	$outfile = 'intersection_' . $data1->{'basename'} . '_' . 
		$data2->{'basename'};
}
my $success = write_tim_data_file(
	'data'     => $output,
	'filename' => $outfile,
);
if ($success) {
	print " Wrote data file '$success'\n";
}
else {
	print " Unable to write data file!\n";
}

# write gff file
if ($gff) {
	# write a gff file if requested
	unless ($type) {
		$type = 'nucleosome_intersection';
	}
	unless ($source) {
		# set default source, the name of this program
		$source = 'intersect_nucs.pl';
	}
	$success = convert_and_write_to_gff_file(
		'data'     => $output,
		'filename' => $outfile,
		'version'  => 3,
		'score'    => 5,
		'strand'   => 3,
		'type'     => $type,
		'source'   => $source,
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




# Main nucleosome intersection subroutine
sub intersect_nucs {
	
	# get the source data structures
	my ($target_data, $reference_data) = @_;
	
	# initialize the output data structure
	my $out_data = generate_tim_data_structure(
		'nucleosome',
		qw(
			Chromosome
			Start
			Stop
			Direction
			Length
			Delta_Occupancy
			Target_ID
			Reference_ID
		)
	) or die " unable to generate tim data structure!\n";
	
	# add metadata to output data structure
	$out_data->{6}{'file'} = $target_data->{'filename'}; 
	$out_data->{7}{'file'} = $reference_data->{'filename'}; 
		
	
	# index the reference data
	unless ( index_data_table($reference_data, 3000) ) {
		# I'm hard encoding the index increment value to 3 kb
		# that should be about 20 nucleosomes
		die " unable to index reference data table for '" . 
			$reference_data->{filename} . "'\n";
	}
# 	{
# 		open FH, ">intersect_dumper.txt";
# 		print FH Dumper($reference_data);
# 		close FH;
# 	}
	
	# identify indices
	my $t_chromo_i  = find_column_index($target_data, '^chr|seq|refseq');
	my $t_start_i   = find_column_index($target_data, 'start');
	my $t_stop_i    = find_column_index($target_data, 'stop|end');
	my $t_name_i    = find_column_index($target_data, 'name|NucleosomeID');
	my $t_score_i   = find_column_index($target_data, 'score|occupancy');
	my $r_chromo_i  = find_column_index($reference_data, '^chr|seq|refseq');
	my $r_start_i   = find_column_index($reference_data, 'start');
	my $r_stop_i    = find_column_index($reference_data, 'stop|end');
	my $r_name_i    = find_column_index($reference_data, 'name|NucleosomeID');
	my $r_score_i   = find_column_index($reference_data, 'score|occupancy');
	unless (
		defined $t_chromo_i and 
		defined $t_start_i and 
		defined $t_stop_i and 
		defined $r_chromo_i and 
		defined $r_start_i and 
		defined $r_stop_i 
	) {
		die " unable to identify one or more required indexs in the source files!\n";
	}
	
	# set the strand for the target nucleosomes if requested
	my $t_strand_i;
	if ($set_strand) {
		$t_strand_i   = find_column_index($target_data, 'strand');
		
		# skip strand if we can't find it
		unless (defined $t_strand_i) {
			warn " unable to identify strand column for target nucleosome file"
				. " \n   '" . $target_data->{filename} . "'! ignoring strand\n";
			$set_strand = 0;
		}
	}
	
	
	# begin intersection
	my $t = $target_data->{'data_table'};
	my $r = $reference_data->{'data_table'};
	for (my $trow = 1; $trow <= $target_data->{'last_row'}; $trow++) {
		
		# calculate index lookup value
		my $lookup = int( $t->[$trow][$t_start_i] / 3000);
		my $t_chr = $t->[$trow][$t_chromo_i];
		
		# check whether there is an index value
		if (exists 
			$reference_data->{'index'}{$t_chr}{$lookup}
		) {
			# there exists an index value
			# that means a possibility for overlap
			
			# calculate the target midpoint
				# since we are working with nucleosomes here, we primarily 
				# use the midpoint in calculating positions
			my $t_mid = sprintf "%.0f", ($t->[$trow][$t_start_i] + 
				$t->[$trow][$t_stop_i] ) / 2;
				
			
			# find intersecting nucs
			for (
				# walk through the list of reference nucleosomes starting at 
				# the index position
				my $rrow = $reference_data->{'index'}{$t_chr}{$lookup}; 
				$rrow <= $reference_data->{'last_row'};
				$rrow++
			) {
				# this loop will proceed all the way to the end of the data 
				# table
				# there are conditionals below which break the loop once 
				# we've passed the region where we might find overlapping 
				# nucleosomes
				
				
				# check intersection
				if (
					$t_chr eq $r->[$rrow][$r_chromo_i] and
					$t_mid >= $r->[$rrow][$r_start_i] and 
					$t_mid <= $r->[$rrow][$r_stop_i]
				) {
					# we have intersection!
					
					# calculate the reference midpoint
					my $r_mid = sprintf "%.0f", 
						( $r->[$rrow][$r_start_i] + $r->[$rrow][$r_stop_i] )/2;
					
					# determine extent and coordinates of overlap
					# according to implied strand
					my ($length, $direction, $start, $stop);
					if ($set_strand) {
						# implied strand
						
						if ($t->[$trow][$t_strand_i] =~ /^1|\+|f|w/i) {
							# target is forward strand
							# no need to flip coordinates
							
							$length = $t_mid - $r_mid;
							if ($t_mid == $r_mid) {
								$direction = '.';
								$start = $t_mid;
								$stop = $t_mid;
							}
							elsif ($t_mid < $r_mid) {
								# leftward shift in position
								$direction = 'r';
								$start = $t_mid;
								$stop = $r_mid;
							}
							elsif ($t_mid > $r_mid) {
								# rightward shift in position
								$direction = 'f';
								$start = $r_mid;
								$stop = $t_mid;
							}
						}
						
						elsif ($t->[$trow][$t_strand_i] =~ /^\-|r|c/i) {
							# target is reverse strand
							# we will essentially flip the coordinates around 
							
							$length = -($t_mid - $r_mid);
							if ($t_mid == $r_mid) {
								$direction = '.';
								$start = $t_mid;
								$stop = $t_mid;
							}
							elsif ($t_mid < $r_mid) {
								# leftward becomes rightward shift in position
								$direction = 'f';
								$start = $r_mid;
								$stop = $t_mid;
							}
							elsif ($t_mid > $r_mid) {
								# rightward becomes leftward shift in position
								$direction = 'r';
								$start = $t_mid;
								$stop = $r_mid;
							}
						}
						
						else {
							# huh!!!????
							die " unrecgnizable strand symbol '" . 
								$t->[$trow][$t_strand_i] . "' at target data" .
								" row $trow!\n";
						}
					}
					
					else {
						# no strand implied
						
						$length = $t_mid - $r_mid;
						if ($t_mid == $r_mid) {
							$direction = '.';
							$start = $t_mid;
							$stop = $t_mid;
						}
						elsif ($t_mid < $r_mid) {
							# leftward shift in position
							$direction = 'r';
							$start = $t_mid;
							$stop = $r_mid;
						}
						elsif ($t_mid > $r_mid) {
							# rightward shift in position
							$direction = 'f';
							$start = $r_mid;
							$stop = $t_mid;
						}
					}
					
					# determine names if available
					my $t_name = get_nucleosome_name(
						$target_data, $trow, $t_name_i);
					my $r_name = get_nucleosome_name(
						$reference_data, $rrow, $r_name_i);
					
					# determine change in score
					my $score;
					if (defined $t_score_i and defined $r_score_i) {
						if (
							$t->[$trow][$t_score_i] ne '.' and 
							$r->[$rrow][$r_score_i] ne '.'
						) {
							# calculate difference in score value if values 
							# are not nulls
							$score = $t->[$trow][$t_score_i] - 
								$r->[$rrow][$r_score_i];
						}
						else {
							$score = '.';
						}
					}
					else {
						$score = '.';
					}
					
					# record the information
					push @{ $out_data->{'data_table'} }, [
						$t_chr,
						$start,
						$stop,
						$direction,
						$length,
						$score,
						$t_name,
						$r_name
					];
					
				}
				
				# check chromosome 
				elsif ($t_chr ne $r->[$rrow][$r_chromo_i]) {
					# appears we have moved off of the chromosome
					last;
				}
				
				# check position
				elsif ($r->[$rrow][$r_start_i] > $t->[$trow][$t_stop_i]) {
					# we have moved beyond the end of the target nucleosome
					last;
				}
			}
		
		}
		else {
			# print " no index for $t->[$trow][$t_chromo_i] -> $lookup\n";
		}
		
	} # finished intersecting
	
	# update data
	$out_data->{'last_row'} = scalar @{ $out_data->{'data_table'} } - 1;
	
	#done
	return $out_data;
}




# Identify the target or reference nucleosome name from the source data
sub get_nucleosome_name {
	my ($data_ref, $row, $name_i) = @_;
	my $name;
	
	# first check whether we have an index
	if (defined $name_i) {
		$name = $data_ref->{'data_table'}->[$row][$name_i];
	}
	
	# check whether the source file was GFF
	elsif ($data_ref->{gff}) {
		# we'll need to extract from the group column
		foreach (split /\s*;\s*/, $data_ref->{'data_table'}->[$row][8]) {
			if (/Name=(.+)/) {
				$name = $1;
				last;
			}
		}
		unless ($name) {
			# default in case there is no identifiable name field
			$name = '.';
		}
	}
	
	# else, no name
	else {
		$name = '.';
	}
	
	return $name;
}




# Print summary information regarding the intersection
sub print_statistics {
	my ($data_ref, $target_ref, $reference_ref) = @_;
	
	### initialize counts
	my %names;
	my $once = 0;
	my $twice = 0;
	my $thrice = 0;
	my $more = 0;
	my $left = 0;
	my $right = 0;
	my $noshift = 0;
	my @shifts;
	my @occupancies;
	my $number = $data_ref->{'last_row'};
	
	
	### collect the data
	my $table = $data_ref->{'data_table'};
	for my $row (1 .. $number) {
		
		# check the reference name
		if ($table->[$row][7] ne '.') {
			if (exists $names{ $table->[$row][7] } ) {
				# we've seen this name before
				# this reference nuc was intersected more than once
				$names{ $table->[$row][7] } += 1;
			}
			else {
				# this nuc has only been intersected once so far
				$names{ $table->[$row][7] } = 1;
			}
		}
		
		# check direction
		if ($table->[$row][3] eq 'f') {
			$right++;
		}
		elsif ($table->[$row][3] eq 'r') {
			$left++;
		}
		else {
			$noshift++;
		}
		
		# record change amounts
		push @shifts, $table->[$row][4];
		push @occupancies, $table->[$row][5] if $table->[$row][5] ne '.';
	}
	
	# count the number of hits
	foreach (keys %names) {
		if ($names{$_} == 1) {
			$once++;
		}
		elsif ($names{$_} == 2) {
			$twice++;
		}
		elsif ($names{$_} == 3) {
			$thrice++;
		}
		else {
			$more++;
		}
	}
		
	
	### print the results
	
	# total numbers
	print "\n Identifed " . format_with_commas($data_ref->{last_row}) . 
		" nucleosome intersections\n";
	printf "   %.1f%% of target nucleosomes\n", 
		($number / $target_ref->{'last_row'}) * 100;
	printf "   %.1f%% of reference nucleosomes\n", 
		($number / $reference_ref->{'last_row'}) * 100;
	
	if ($number) {
		
		# number of nucleosome hits
		print "  " . format_with_commas($once);
		printf " (%.1f%%) nucleosomes intersected unique reference nucleosomes\n", 
			( ($once/$number) * 100 );
		if ($twice) {
			print "  " . format_with_commas($twice);
			printf " (%.1f%%) intersected the same reference nucleosome twice\n", 
				( ($twice/$number) * 100 );
		}
		if ($thrice) {
			print "  " . format_with_commas($thrice);
			printf " (%.1f%%) intersected the same reference nucleosome three times\n", 
				( ($thrice/$number) * 100 );
		}
		if ($more) {
			print "  " . format_with_commas($more);
			printf " (%.1f%%) intersected the same reference nucleosome four or more times\n", 
				( ($more/$number) * 100 );
		}
		
		# direction
		print " " . format_with_commas($left);
		printf " nucleosomes (%.1f%%) had a leftward shift\n", 
			( ($left/$number) * 100 );
		print " " . format_with_commas($right);
		printf " nucleosomes (%.1f%%) had a rightward shift\n", 
			( ($right/$number) * 100 );
		print " " . format_with_commas($noshift);
		printf " nucleosomes (%.1f%%) had no shift\n", 
			( ($noshift/$number) * 100 );
		print "  The mean shift was ", sprintf( "%.0f", mean(@shifts) ), 
			" +/- ", sprintf( "%.0f", stddevp(@shifts) ), " bp, range ", 
			min(@shifts), "..", max(@shifts), " bp\n";
	}
	
	# occupancy
	if (@occupancies) {
		print " The mean nucleosome occupancy change was ", 
			sprintf( "%.0f", mean(@occupancies) ), " +/- ", 
			sprintf( "%.0f", stddevp(@occupancies) ), "\n";
	}
	print "\n"; # add an extra line of space
}





__END__

=head1 NAME

intersect_nucs.pl

A script to intersect two lists of nucleosomes.

=head1 SYNOPSIS

intersect_nucs.pl [--options...] <filename_1> <filename_2>
  
  Options:
  --in1 <filename1>
  --in2 <filename2>
  --out <filename>
  --force_strand
  --gff
  --type <gff_type>
  --source <gff_source>
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in1 <filename>

=item --in2 <filename>

Specify two files of nucleosome lists. The files must contain sorted
genomic position coordinates for each nucleosome. Supported file formats
include any text data file with chromosome, start, stop, and name. The
file with the least number of nucleosomes is automatically designated as
the target, while the file with the most is designated as the reference
list. When files with equivalent numbers are provided, the first file 
is target.

=item --out <filename>

Specify the output file name. The default is "intersection_" appended 
with both input names.

=item --force_strand

Force the target nucleosomes to be considered as stranded. This enforces 
an orientation and affects the direction of any reported nucleosome shift. 
A column with a label including 'strand' is required in the target file. 
The default is false.

=item --gff

Indicate whether a GFF file should be written in addition to the standard 
text data file. The GFF file version is 3. The default value is false.

=item --type <gff_type>

Provide the text to be used as the GFF type (or method) used in 
writing the GFF file. The default value is "nucleosome_intersection".

=item --source <gff_source>

Provide the text to be used as the GFF source used in writing the 
GFF file. The default value is the name of this program.

=item --gz

Specify whether (or not) the output files should be compressed 
with gzip. 

=item --version

Print the version number.

=item --help

Display the POD documentation of the script. 

=back

=head1 DESCRIPTION

This program will intersect two lists of nucleosomes. It will identify which 
nucleosomes overlap, the direction and extent of shift of their midpoints, 
and the delta change of the occupancy (or score value). The file with the 
least number of nucleosomes is automatically designated as the target list, 
and the file with the most number is designated as the reference list. The 
reference file, at least, must be sorted in genomic order; otherwise, 
intersection will yield undesirable results.

The program will output a tim data text file of the intersections, which
include the start and stop points that indicate the positions and extent of
the midpoint shift, the direction of shift, and the name of the
intersecting nucleosomes. Optionally a GFF file may also be written as well.

The target nucleosomes may have strand optionally imposed. This is useful 
when working with nucleosomes that are associated with stranded genomic 
features, for example, nucleosomes flanking a transcription start site.

A summary and statistics of the intersection are printed to standard output 
upon completion.

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
