#!/usr/bin/perl

# This script will run the k-means cluster analysis

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename qw(fileparse);
use Algorithm::Cluster::Record;
my $VERSION = '1.0.2';

print "\n A script to run the k-means cluster analysis\n\n";

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
	$number,
	$runs,
	$method,
	$distribution,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the solexa data file
	'out=s'     => \$outfile, # name of base output file 
	'num=i'     => \$number, # the number of clusters
	'run=i'     => \$runs, # the number of runs to perform the
	'method=s'  => \$method, # the method to perform
	'dist=s'    => \$distribution, # similiarity metric of measuring distance
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
	print " Biotoolbox script run_cluster.pl, version $VERSION\n\n";
	exit;
}




### Check for requirements
# check input file
unless ($infile) {
	$infile = shift @ARGV or
		die "  OOPS! No source data file specified! \n use --help\n";
}
my ($basename, $path, $extension) = fileparse($infile, qw(.txt .cdt));

unless ( 
	($basename =~ /_tview$/ and $extension eq '.txt') or $extension eq '.cdt'
) {
	die " The input file does not look like the appropriate file type\n" . 
		" Consider using 'manipulate_datasets.pl' with the treeview export function\n" . 
		" Use --help for more information\n";
}

# set defaults
unless ($outfile) {
	$outfile = $basename;
}
unless ($number) {
	# I like six clusters - an informative but not overwhelming number
	$number = 6;
}
unless ($runs) {
	# I like 200 runs for some reason - a lot but not too much
	$runs = 200;
}
unless ($method) {
	# default is k-means
	$method = 'a'; # arithmetic mean
}
unless ($distribution) {
	# default is Euclidean distance
	# it looks like this is the default for the module as well
	$distribution = 'e';
}


### Load the file
my $record = Algorithm::Cluster::Record->new() or 
	die " unable to intialize Cluster::Record object!\n";

open INPUT, $infile;
$record->read(*INPUT);
close INPUT;



### Run the cluster
print " running cluster analysis....\n";
my ($clusterid, $error, $nfound) = $record->kcluster(
	'nclusters'     => $number,
	'npass'         => $runs,
	'method'        => $method,
	'dist'          => $distribution,
);
print " An optimal solution was identified $nfound times\n";



### Output results
$record->save(
	'jobname'       => $outfile,
	'geneclusters'  => $clusterid,
);
print " Finished with job '$outfile'\n";




__END__

=head1 NAME

run_cluster.pl

=head1 SYNOPSIS

run_cluster.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --out <jobname> 
  --num <integer>
  --run <integer>
  --method [a|m]
  --dist [c|a|u|x|s|k|e|b]
  --version
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input file. The file should be a simple tab-delimited text 
file, genes (features) as rows, experimental data (microarray or sequencing 
data) should be columns. The first column contains unique gene identifiers. 
A column header row is expected. Standard tim data text files with 
metadata lines should be exported to a compatible format using the treeview 
function in the B<manipulate_datasets.pl> script. A .cdt file generated 
from this may also be used.

=item --out <jobname>

Specify the output jobname, which will be the basename of the output files. 
By default it uses the input base filename.

=item --num <integer>

Specify the number of clusters to identify. Default value is 6 – enough to be 
informative but not overwhelming.

=item --run <integer>

Enter the number of times to run the cluster algorithm to find a solution. 
The default value is 200.

=item --method [a|m]

Specify the method of finding the center of a cluster. Two values are 
allowed, arithmetic mean (a) and median (m). Default is mean.

=item --dist [c|a|u|x|s|k|e|b]

Specify the distance function to be used. Several options are available.
	
	c  correlation
	a  absolute value of the correlation
	u  uncentered correlation
	x  absolute uncentered correlation
	s  Spearman's rank correlation
	k  Kendall's tau
	e  Euclidean distance
	b  City-block distance
	
The default value is 'e', Euclidean distance.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program is a wrapper around the Cluster 3.0 C library, which identifies 
clusters between genes. Currently the program performs the k-means or 
k-medians functions, although other functions could be implemented if 
requested.



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


