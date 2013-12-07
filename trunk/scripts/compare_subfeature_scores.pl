#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::Lite qw(min max range);
use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure
	find_column_index
);
use Bio::ToolBox::file_helper qw(
	open_tim_data_file 
	write_tim_data_file 
);
my $VERSION = '1.14';

print "\n This program will compare scores from multiple subfeatures\n\n";

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
	$gene_i,
	$transcript_i,
	$score_i,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the input data file
	'out=s'     => \$outfile, # name of output file 
	'gene=i'    => \$gene_i, # gene column index
	'transcript=i' => \$transcript_i, # transcript column index
	'score=i'   => \$score_i, # score column index
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
	print " Biotoolbox script compare_subfeature_scores.pl, version $VERSION\n\n";
	exit;
}



### Check for requirements
unless ($infile) {
	$infile = shift @ARGV or
		die " no input file! use --help for more information\n";
}
unless ($outfile) {
	die " please define an output file name!\n";
}
unless (defined $gz) {
	$gz = 0;
}



### Open input file
my ($in_fh, $in_md) = open_tim_data_file($infile) or 
	die " unable to open input file!\n";

# Check for columns
identify_columns();

# Load the gene tree
my $tree = load_gene_tree();
$in_fh->close;



### Process the tree
my $outdata = process_gene_tree($tree);



### Finished
my $success = write_tim_data_file(
	'data'     => $outdata,
	'filename' => $outfile,
	'gz'       => $gz,
);
if ($success) {
	print " wrote file '$success'\n";
}
else {
	# failure! the subroutine will have printed error messages
	print " unable to write file!\n";
}




########################   Subroutines   ###################################


sub identify_columns {
	
	# we will ask user if they didn't provide us with the required indices
	unless (defined $gene_i and defined $transcript_i and defined $score_i) {
		
		# identify possibilites
		unless ($gene_i) {
			$gene_i = find_column_index($in_md, 'gene|parent');
		}
		unless ($transcript_i) {
			$transcript_i = find_column_index($in_md, 'transcript|name');
		}
		unless ($score_i) {
			$score_i = find_column_index($in_md, 'score');
		}
		
		# present the list to the user
		print " These are the columns in the input file:\n";
		map { print "  $_\t", $in_md->{$_}{'name'}, "\n"} 
			(0 .. ($in_md->{number_columns} - 1) );
		
		# process answers
		$gene_i = process_answer('parent feature name', $gene_i);
		$transcript_i = process_answer('subfeature name', $transcript_i);
		$score_i = process_answer('score', $score_i);
	}
}



sub process_answer {
	my ($name, $index) = @_;
	
	# request from user
	print " Enter the index for the $name [$index]   ";
	my $answer = <STDIN>;
	chomp $answer;
	if ($answer =~ /^\d+$/ and exists $in_md->{$answer}) {
		$index = $answer;
	}
	# else we use the default
	
	unless (defined $index) {
		die " no column index for $name defined!\n";
	}
	return $index;
}



sub load_gene_tree {
	
	# we will generate a hash of hashes
	# first key will be the gene name
	# second key will be the transcript name
	# the value will be the score
	my %tree;
	
	# process the input file and load the tree
	while (my $line = $in_fh->getline) {
		chomp $line;
		my @data = split /\t/, $line;
		
		# check the tree first and load
		if (exists $tree{ $data[$gene_i] }{ $data[$transcript_i] } ) {
			warn " gene $data[$gene_i] transcript $data[$transcript_i]" .
				" exists more than once! skipping\n";
			next;
		}
		else {
			$tree{ $data[$gene_i] }{ $data[$transcript_i] } = $data[$score_i];
		}
	}
	
	return \%tree;
}



sub process_gene_tree {
	
	my $tree = shift;
	
	# generate output data structure
	my $output = generate_tim_data_structure( 
		'Gene_transcripts',
		qw(
			Parent
			Number_subfeatures
			Min_name
			Min_score
			Max_name
			Max_score
			Range
	) );
	push @{ $output->{'other'} }, " # original_input_file $infile\n";
	
	# walk through each gene
	foreach my $gene (sort {$a cmp $b} keys %$tree) {
		
		# determine the number of transcripts
		my $number = scalar keys %{ $tree->{$gene} };
		
		# only 1 transcript found
		if ($number == 1) {
			
			# get the transcript id
			my ($transcript) = keys %{ $tree->{$gene} };
			
			# record it
			push @{ $output->{'data_table'} },
				[ 
					$gene,
					$number,
					$transcript,
					$tree->{$gene}{$transcript},
					$transcript,
					$tree->{$gene}{$transcript},
					0
				];
			$output->{'last_row'} += 1;
		}
		
		# more than 1 transcript
		else {
			# we need to identify which is the minimum and maximum
			
			# determine the scores
			my $min = min(values %{ $tree->{$gene} });
			my $max = max(values %{ $tree->{$gene} });
			
			# identify the min max ids
			my ($min_id, $max_id);
			foreach my $id (sort {$a cmp $b} keys %{ $tree->{$gene} } ) {
				unless ($min_id) {
					$min_id = $id if $tree->{$gene}{$id} == $min;
				}
				unless ($max_id) {
					$max_id = $id if $tree->{$gene}{$id} == $max;
				}
			}
			
			# record it
			push @{ $output->{'data_table'} },
				[
					$gene,
					$number,
					$min_id,
					$min,
					$max_id,
					$max,
					range(values %{ $tree->{$gene} }), # the range of the scores
				 ];
			$output->{'last_row'} += 1;
		}
	}
	
	# done
	return $output;
}





__END__

=head1 NAME

compare_subfeature_scores.pl

A script to compare the scores between one or more subfeatures.

=head1 SYNOPSIS

compare_subfeature_scores.pl --in <filename> --out <filename>
  
  Options:
  --in <filename>
  --out <filename> 
  --parent <index>
  --subfeature <index>
  --score <index>
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input file. It should be a tab-delimited text file with 
headers. At least three columns must be present: the name of the 
parent (gene) feature, the name of subfeature (transcript) and the 
score. The file may be compressed with gzip.

=item --out <filename>

Specify the output filename. 

=item --parent <index>

Optionally specify the index column for the parent or gene name. If 
not specified, the program will interactively present a list of columns 
to choose from.

=item --subfeature <index>

Optionally specify the index column for the subfeature or transcript 
name. If not specified, the program will interactively present a list 
of columns to choose from.

=item --score <index>

Optionally specify the index column for the score to compare. If 
not specified, the program will interactively present a list of columns 
to choose from.

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will compare the scores of all the subfeatures of a parent 
feature. For example, comparing RNA expression of all of the alternative 
transcripts from a gene. As input, it expects a tab-delimited text file 
with at least three columns: the name of the parent (e.g. gene) feature 
(not expected to be unique in the file), the name of the subfeature (e.g. 
transcript; must be unique with respect to each parent feature), and a 
score for each subfeature. Such a file may be generated using the 
biotoolbox script L<get_gene_regions.pl> followed by L<get_datasets.pl>.

The program will output a new file. Each line will represent one 
parent feature. The columns include the parent feature name, number of 
subfeatures, the minimum and maximum subfeature names and scores, and 
the range of scores.

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
