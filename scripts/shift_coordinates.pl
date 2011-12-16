#!/usr/bin/perl 

# This script will shift the genomic coordinates by a specified value
# 

use strict;
use File::Copy;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper qw(
	open_tim_data_file
	write_tim_data_file
	open_to_write_fh
);
my $VERSION = '1.0.0';

print "\n A script to shift genomic coordinates by a specified value\n";

# useage
unless (@ARGV) {
	print " Usage: shift_coordinates.pl <filename> [number]\n";
	exit;
}

# run parameters
my $filename = $ARGV[0];
my $value = $ARGV[1];



# Ask for which conversion
unless ($value) {
	print " enter the value by which to shift the coordinates   ";
	$value = <STDIN>;
	chomp $value;
}



# Open the input file
my ($fh, $metadata_ref) = open_tim_data_file($filename);

# Identify the columns
my $pos_index; # first position, or start
my $pos2_index; # a second position, or end
if ($metadata_ref->{extension} =~ /sgr/i) { 
	# file is a *.sgr
	$pos_index = 1;
} 
elsif ($metadata_ref->{extension} =~ /gff/i) { 
	# file is a gff file
	$pos_index = 3;
	$pos2_index = 4;
} 
elsif ($metadata_ref->{extension} =~ /bed/i) { 
	# file is a bed file
	$pos_index = 1;
	$pos2_index = 2;
} 
else { 
	# a non-standard file
	# identify the indices to the genomic coordinates
	my ($chromo_index, $pos_index, $pos2_index);
	for (my $i = 0; $i < $metadata_ref->{'number_columns'}; $i++) {
		# check the names of each column
		# looking for chromo, start, stop
		if ($metadata_ref->{$i}{'name'} =~ /start|position/i) {
			$pos_index = $i;
		}
		elsif ($metadata_ref->{$i}{'name'} =~ /stop|end/i) {
			$pos2_index = $i;
		}
	}
	unless ( defined $pos_index ) {
		# check that we have these basic genomic coordinates
		die " unable to identify chromosome, start, and/or stop datasets!";
	}	
}




# Open the output file
my $gz;
if ($metadata_ref->{extension} =~ /gz/) {
	$gz = 1;
}
else {
	$gz = 0;
}
push @{ $metadata_ref->{other} }, 
	# record a miscelleneous header line to explain the conversion
	"# Shifted coordinates by $value";
$metadata_ref->{data_table} = []; # create an empty data table
push @{ $metadata_ref->{data_table} }, $metadata_ref->{'column_names'};

# write the new file
my $new_filename = $metadata_ref->{basename} . '.shifted' . 
	$metadata_ref->{extension};
unless ( write_tim_data_file( {
		'data'       => $metadata_ref,
		'filename'   => $new_filename,
		'gz'         => $gz,
	} )
) {
	die " unable to open file for writing!\n";
}
my $out_fh = open_to_write_fh($new_filename, $gz, 1) or 
	die " unable to re-open file for writing!\n";



# Do the conversion
my $totalcount = 0; # a running tally
while (my $line = $fh->getline) {
	
	# process each line
	chomp $line;
	my @data = split /\t/, $line;
	$data[$pos_index] = $data[$pos_index] + $value;
	if ($pos2_index) {
		$data[$pos2_index] = $data[$pos2_index] + $value;
	}
			
	# put back
	$totalcount++;
	print {$out_fh} join("\t", @data) . "\n";	
}
close INFILE;
print " converted $totalcount features to new coordinates\n";


# Finish
$fh->close;
$out_fh->close;
print " wrote converted file '$new_filename'\n";


