#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::data_helper qw(find_column_index);
use Bio::ToolBox::file_helper qw(open_tim_data_file);
my $cluster_ok;
eval {
	require Algorithm::Cluster::Record;
	$cluster_ok = 1;
};

my $VERSION = '1.14';

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
unless ($cluster_ok) {
	die "Module Algorithm::Cluster must be installed to run this script.\n";
} 

# check input file
unless ($infile) {
	$infile = shift @ARGV or
		die "  OOPS! No source data file specified! \n use --help\n";
}


# set defaults
unless ($number) {
	# I like six clusters - an informative but not overwhelming number
	$number = 6;
}
unless ($runs) {
	$runs = 500;
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


### Check the input file format
# open to read metadata
# we're not actually loading the file
# we just want to verify it looks ok
my ($in_fh, $metadata) = open_tim_data_file($infile);
my $check = 1; # assume ok to begin with
my $error;

# check comment lines
if (scalar @{ $metadata->{'other'}} != 0) {
	$check = 0;
	$error .= "  file has extraneous comment lines\n";
}

# check first column
if ($metadata->{0}{'name'} !~ /name|id|gene|transcript/i) {
	# may not be lethal
	$error .= "  first column name is unusual\n";
}

# check for column data
for (my $i = 0; $i < $metadata->{'number_columns'}; $i++) {
	if (not exists $metadata->{$i}{'AUTO'}) {
		# no automatically generated column metadata 
		# suggests there was column metadata in the file
		$check = 0;
		$error .= "  file column $i has extra metadata\n";
	}
}

# check for extraneous data columns
foreach (qw(chr seq start stop end strand type class source phase)) {
	my $i = find_column_index($metadata, "^$_");
	if (defined $i) {
		$check = 0;
		$error .= "  file has extraneous column '$_' at position $i\n";
	}
}

# print errors
if ($check == 0 and $error) {
	print " input file did not pass validation for the following reasons\n";
	print $error;
	print " consider using the treeview export function in manipulate_datasets.pl\n";
	exit;
}
elsif ($check == 1 and $error) {
	print " input file may not be valid for the following reasons\n";
	print $error;
}
else {
	print " input file appears to be valid\n";
}
$in_fh->close;




### Load the file
my $start_time = time;
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
printf " An optimal solution was identified %s times in %.1f minutes\n",
	$nfound, (time - $start_time)/60;



### Output results
unless ($outfile) {
	$outfile = $metadata->{'path'} . $metadata->{'basename'};
}
$record->save(
	'jobname'       => $outfile,
	'geneclusters'  => $clusterid,
);
print " Finished with job '$outfile'\n";




__END__

=head1 NAME

run_cluster.pl

A script to run the k-means cluster analysis.

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

Specify the number of clusters to identify. Default value is 6.

=item --run <integer>

Enter the number of iterations to run the cluster algorithm to find an 
optimal solution. The default value is 500.

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

Please refer to the Cluster 3 documentation for more detailed information 
regarding the implementation and detailed methods. Documentation may be 
found at L<http://bonsai.hgc.jp/~mdehoon/software/cluster/>.

Select the desired number of clusters that are appropriate for your dataset 
and an appropriate number of iterations. The default values are fine to 
start with, but should be customized for your dataset. In general, empirically 
test a range of cluster numbers, e.g. 2 to 12, to find the optimal cluster 
number that is both informative and manageable. Increasing the number of 
iterations will increase confidence at the expense of compute time. The goal  
is to find an optimal solution more than once; the more times a solution 
has been found, the higher the confidence. Note that noisy or very large 
datasets may never yield more than 1 solution.
 
The resulting CDT files may be visualized using the Java Treeview program, 
found at L<http://jtreeview.sourceforge.net>. 

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
