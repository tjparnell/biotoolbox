#!/usr/bin/perl

# This script will generate a tim data file of genomic bins
# This essentially duplicates the functionality of 
# tim_db_helper::get_new_genome_list
# but allows for massive files without memory constraints

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	generate_tim_data_structure
);
use tim_db_helper qw(
	open_db_connection
);
use tim_file_helper qw(
	write_tim_data_file
	open_to_write_fh
);
use tim_db_helper::config;
my $VERSION = '1.5.8';

print "\n This program will generate a file of genomic bins\n";

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
	$database,
	$win,
	$step,
	$keep_mito,
	$split,
	$max,
	$outfile,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'db=s'      => \$database, # the name of Bio::DB database
	'win=i'     => \$win, # the window or interval size
	'step=i'    => \$step, # the step size
	'mito!'     => \$keep_mito, # keep the mitochondrial chromosome
	'split!'    => \$split, # split the file into segments
	'max=i'     => \$max, # maximum number of lines per split file part
	'out=s'     => \$outfile, # name of output file 
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
	print " Biotoolbox script generate_genomic_bins.pl, version $VERSION\n\n";
	exit;
}



### Check for requirements
unless ($database) {
	die "  No database specified! see help\n";
}
unless ($outfile) {
	die "  No output file name specified! see help\n";
}
unless (defined $gz) {
	$gz = 0;
}
unless (defined $keep_mito) {
	$keep_mito = 0;
}


### Set up defaults
unless (defined $win) {
	$win = 
		$TIM_CONFIG->param("$database\.window") ||
		$TIM_CONFIG->param('default_db.window');
	print "  Using default window size of $win bp\n";
}
unless (defined $step) {
	$step = $win;
	print "  Using default step size of $step\n";
}



### Open database connection
my $db = open_db_connection($database) or 
	die " unable to open database connection!\n";



### Set up data structure and open file

# new data structure
my $data_ref = generate_tim_data_structure(
		'genome',
		'Chromosome',
		'Start',
		'Stop',
) or die " unable to generate tim data structure!\n";
$data_ref->{'db'}      = $database; # the db name
$data_ref->{1}{'win'}  = $win; 
$data_ref->{1}{'step'} = $step; 

# write a tim data file
my $new_output = write_tim_data_file( {
	'data'       => $data_ref,
	'filename'   => $outfile,
	'gz'         => $gz,
	'format'     => 'text',
} );
unless ($new_output) {
	die " unable to write output file!\n";
}

# reopen file for writing the data stream
my $out_fh = open_to_write_fh($new_output, $gz, 1) or 
	die " unable to re-open output file for writing data stream!\n";



### Collect the windows
collect_genomic_windows();

# finished
$out_fh->close;



### Split if requested
if ($split) {
	# to make this easy, we'll be simply calling split_data_file.pl
	print " splitting the file....\n";
	
	# build executable
	my $executable = $Bin . '/split_data_file.pl';
	$executable .= ' --index 0'; # splitting by chromosome
	if ($max) {
		# maximum file size is requested
		$executable .= " --max $max"; 
	}
	$executable .= " --in $new_output";
	
	# execute
	!system $executable or 
		die " unable to run split_data_file.pl!\n"; 
	
	# remove the original file
	unlink $new_output;
}


### Done
print " Finished\n";








########################   Subroutines   ###################################


sub collect_genomic_windows {
	# this is essentially copy and pasted from tim_db_helper::get_new_genome_list
	
	# Collect the chromosomes from the db
	my @chromosomes = $db->seq_ids; 
	unless (@chromosomes) {
		die " unable to retrieve chromosomes from the database!\n";
	}
	
	# Get the names of chromosomes to avoid
	my @excluded_chromosomes = 
			$TIM_CONFIG->param("$database\.chromosome_exclude") ||
			$TIM_CONFIG->param('default_db.chromosome_exclude');
	
	# Collect the genomic windows
	print "   Generating $win bp windows in $step bp increments\n";
	my $count = 0;
	foreach my $chr (@chromosomes) {
		
		# check for excluded chromosomes
		my $skip_chr = 0;
		foreach (@excluded_chromosomes) {
			if ($chr eq $_) {
				$skip_chr = 1;
				last;
			}
		}
		next if $skip_chr;
		
		my $chrobj = $db->segment($chr);
		my $length = $chrobj->length;
		for (my $start = 1; $start <= $length; $start += $step) {
			# set the end point
			my $end = $start + $win - 1; 
			
			if ($end > $length) {
				# fix end to the length of the chromosome
				$end = $length;
			} 
			
			# write to the output
			print {$out_fh} join("\t", $chr, $start, $end), "\n";
			$count++;
		}
	}
	print "  Collected $count genomic windows\n";
}








__END__

=head1 NAME

generate_genomic_bins.pl

=head1 SYNOPSIS

generate_genomic_bins.pl [--options...] --db <database> --out <file>
  
  Options:
  --db <database>
  --win <integer>
  --step <integer>
  --mito
  --split
  --max <integer>
  --out <filename> 
  --(no)gz
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <database>

Specify the name of the BioPerl gff database to use as source. This is required.

=item --win <integer>

Optionally specify the window size. The default size is defined in the
configuration file, biotoolbox.cfg.

=item --step <integer>

Optionally indicate the step size. The default is equal to the window size.

=item --mito

Specify whether the mitochondrial chromosome should be included. The 
default is false.

=item --split

Specify whether the output file should be split into individual chromosome 
files. 

=item --max <integer>

When splitting the output file, optionally specify the maximum number of 
lines in each file part.

=item --out <filename>

Specify the output filename. By default it uses 

=item --(no)gz

Specify whether (or not) the output file should be compressed with gzip.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will generate a tim data format file of genomic bins or 
intervals. Generating bins of small size (say, 1 or 10 bp in size) 
for an entire genome, particularly metazoan genomes, produces 
extremely large data files, which demands large memory resources to 
load. This simple script avoids the memory demands by writing 
directly to file as the bins are generated. It will optionally also 
split the file by chromosome (using the biotoolbox script 
split_data_file.pl). 



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


