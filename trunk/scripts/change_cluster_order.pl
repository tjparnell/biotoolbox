#!/usr/bin/env perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::ToolBox::data_helper qw(parse_list);
use Bio::ToolBox::file_helper qw(
	open_to_read_fh
	open_to_write_fh
);
my $VERSION = 1.14;

print "\n A tool for switching the order of cluster files\n";



### Quick help
unless (@ARGV) { # when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Get command line options and initialize values
my ( # command line option variables
	$infile, 
	$outfile, 
	$given_order,
	$keep_group,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # name of input file
	'out=s'     => \$outfile, # name of new output file 
	'order=s'   => \$given_order, # provide a new order
	'keep!'     => \$keep_group, # keep the original group identifiers
	'help'      => \$help, # request help
	'version'   => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Get file name
unless ($infile) {
	$infile = shift @ARGV;
}


### Print help if requested
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script manipulate_datasets.pl, version $VERSION\n\n";
	exit;
}




### Check input files
my ($cdt_file, $kgg_file, $cdt_outfile, $kgg_outfile);
check_defaults();



### Load the input files
my $kgg_data = load_input_kgg_file();
my $cdt_data = load_input_cdt_file();



### Get the new order
my @order = get_new_order();



### Print the output files
print_new_kgg_file();
print_new_cdt_file();




########################   Subroutines   ###################################


sub check_defaults {
	
	# check input files
	if ($infile =~ /\.kgg$/i) {
		$kgg_file = $infile;
		$cdt_file = $infile;
		$cdt_file =~ s/kgg$/cdt/i;
	}
	elsif ($infile =~ /\.cdt$/i) {
		$cdt_file = $infile;
		$kgg_file = $infile;
		$kgg_file =~ s/cdt$/kgg/i;
	}
	unless (-e $kgg_file and -e $cdt_file) {
		die " Cannot find both KGG and CDT files! see help\n";
	}
	
	# check output files
	if ($outfile) {
		# user provided, strip the any extension
		$outfile =~ s/\.(?:kgg|cdt)$//i;
		$cdt_outfile = "$outfile\.cdt";
		$kgg_outfile = "$outfile\.kgg";
	}
	else {
		# overwrite the input files
		$cdt_outfile = $cdt_file;
		$kgg_outfile = $kgg_file;
	}
}



sub load_input_kgg_file {
	my $fh = open_to_read_fh($kgg_file) or 
		die " unable to open input KGG file '$kgg_file'!\n";
	
	# check the header line, which we can dispense
	my $header = $fh->getline; 
	unless ($header =~ /^id\sgroup$/i) {
		die " supposed KGG file does not have ID GROUP column headers!\n";
	}
	
	# load the KGG file
	my %data;
	$data{'original'} = []; # the original group order
	while (my $line = $fh->getline) {
		chomp $line;
		my ($id, $group) = split /\t/, $line;
		if (exists $data{$group}) {
			push @{ $data{$group} }, $id;
		}
		else {
			$data{$group} = [ $id ];
			push @{ $data{'original'} }, $group;
		}
	}
	$fh->close;
	
	return \%data;
}


sub load_input_cdt_file {
	my $fh = open_to_read_fh($cdt_file) or 
		die " unable to open input CDT file '$cdt_file'\n";
	
	# check the header line 
	my $header = $fh->getline;
	unless ($header =~ /^ID/) {
		die " input CDT file does not have a proper header line!?\n";
	}
	
	# initialize data structure
	my %data;
	$data{'header'} = $header;
	
	# load the file
	while (my $line = $fh->getline) {
		my ($id, $remainder) = split(/\t/, $line, 2);
		# remember that there could be an 'EWEIGHT' line
		$data{$id} = $line; 
	}
	$fh->close;
	return \%data;
}


sub get_new_order {
	
	# get the order
	my @order;
	if ($given_order) {
		@order = parse_list($given_order);
	}
	else {
		print " The current group order: ", join(" ", @{$kgg_data->{'original'}} ), "\n";
		print " Please enter a new order as comma-delimited or range\n  ";
		my $response = <STDIN>;
		chomp $response;
		@order = parse_list($response);
	}
	
	# check the order
	unless (@order) {
		die " No new order provided!\n";
	}
	foreach (@order) {
		unless (exists $kgg_data->{$_}) {
			die " non-existant group number provided! Cannot proceed!\n";
		}
	}
	
	return @order;
}


sub print_new_kgg_file {
	
	# open file handle
	my $fh = open_to_write_fh($kgg_outfile) or 
		die " unable to write new KGG file $kgg_outfile!\n";
	$fh->print("ID\tGROUP\n");
	
	# write the new file
	if ($keep_group) {
		# we need to keep the original group identifiers
		foreach my $group (@order) {
			$fh->print(
				join("\n", map {"$_\t$group"} @{ $kgg_data->{$group} } ), "\n"
			);
		}
	}
	else {
		# assign new group identifiers
		my $new_group = 0;
		foreach my $group (@order) {
			$fh->print(
				join("\n", map {"$_\t$new_group"} @{ $kgg_data->{$group} } ), "\n"
			);
			$new_group++;
		}
	}
	$fh->close;
	print " Wrote file $kgg_outfile\n";
}


sub print_new_cdt_file {
	
	# open file handle
	my $fh = open_to_write_fh($cdt_outfile) or 
		die " unable to write new CDT file $cdt_outfile!\n";
	
	# write headers
	$fh->print( $cdt_data->{'header'} );
	if (exists $cdt_data->{'EWEIGHT'}) {
		# this is an optional data line
		# it is written by the Cluster algorithm
		$fh->print( $cdt_data->{'EWEIGHT'} );
	}
	
	# write data
	foreach my $group (@order) {
		foreach my $id ( @{ $kgg_data->{$group} } ) {
			$fh->print( $cdt_data->{$id} );
		}
	}
	$fh->close;	
	print " Wrote file $cdt_outfile\n";
}


__END__

=head1 NAME

change_cluster_order.pl

A script to change the order of gene cluster groups in a file.

=head1 SYNOPSIS

change_cluster_order.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --out <basename>
  --order <numbers,range>
  --keep
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify either the input CDT or KGG file. The program assumes both 
files exist with the same basename and either a .kgg or .cdt file 
extension.

=item --out <basename>

Specify the output filename. By default it uses the base name of the 
input file. An appropriate .kgg and .cdt extension will be added.

=item ---order <numbers,range>

Optionally provide the new order of gene cluster groups. A comma 
delimited list and/or range may be provided, without spaces. For 
example, 1-3,0,4,5.

=item --keep

Optionally keep the same cluster group numbers in the KGG file as 
the original, just in the new order. The default is to renumber 
the group numbers.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will change the order of gene cluster groups in KGG and CDT files. 
These files are generated by the Cluster algorithm and represent k-means 
clusters of genes based on collected experimental data. The cluster groups 
are defined in the KGG file, which is a simple text format with two columns, 
the gene ID and the GROUP number. Sometimes, when comparing data between two  
or more separate cluster analyses, it is useful to re-order the clusters such 
that similar clusters are ranked in a similar order. This script will 
accomplish that goal.

A KGG or CDT file is provided (both are required, but only one needs to be 
provided as a command line argument). The current order of the clusters is 
presented, and a new order is then requested. The new groups are re-labeled 
with new identifiers, or the old group numbers may be retained if requested.

CDT files may be visualized using the Java Treeview program, found at 
L<http://jtreeview.sourceforge.net>. 

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
