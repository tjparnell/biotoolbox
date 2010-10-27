#!/usr/bin/perl

# a script to convert the UCSC chromInfo file into a GFF3 file

# a simple file: chr, size, path

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper qw(
	open_to_read_fh
);


print "\n A script to convert UCSC chromInfo tables to GFF3\n\n";



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
	$type,
	$source,
	$outfile,
#	$gz,
	$help, 
);
GetOptions( 
	'in=s'       => \$infile, # the input file
	'type=s'     => \$type, # the GFF type
	'source=s'   => \$source, # the GFF source
	'out=s'      => \$outfile, # output file name
#	'gz!'        => \$gz, # compress file
	'help'       => \$help, # request help
);



### Check requirements and defaults
unless ($infile) {
	$infile = shift @ARGV or 
		die " input file is required!\n";
}
unless ($type) {
	$type = 'chromosome';
}
unless ($source) {
	$source = 'UCSC';
}
if ($outfile) {
	# add extension as necessary
	unless ($outfile =~ m/\.gff3?$/) {
		$outfile .= '.gff';
	}
}
else {
	# assign default output file
	$outfile = $infile;
	$outfile =~ s/\.txt(?:\.gz)$/.gff/;
}



### Conversion
my $fh = open_to_read_fh($infile) or 
	die " unable to open input file '$infile'!\n";

my %name2chr;
while (my $line = $fh->getline) {
	next if ($line =~ /^#/);
	chomp $line;
	my ($chr, $end, $path) = split /\t/, $line;
	
	# generate chromosome object
	my $chromosome =  Bio::SeqFeature::Generic->new(
		-seq_id        => $chr,
		-source        => $source,
		-primary_tag   => $type,
		-start         => 1,
		-end           => $end,
		-strand        => 0,
		-frame         => '.',
		-display_name  => $chr,
	);
	
	# add ID
	$chromosome->add_tag_value('ID', $chr);

	# add Name
	$chromosome->add_tag_value('Name', $chr);
	
	# determine simple name for loading
	if ($chr =~ /^chr(\d+)$/) {
		# a simple chr number
		my $name = sprintf "%03d", $1;
		$name2chr{$name} = $chromosome;
	}
	elsif ($chr =~ /^chr(\d+)(_.+)$/) {
		# a partial simple chromosome, such as chr1_random
		my $name = sprintf "%03d", $1;
		$name .= $2;
		$name2chr{$name} = $chromosome;
	}
	else {
		# a complex name
		$name2chr{$chr} = $chromosome;
	}
}


### Write out GFF
my $gff =Bio::Tools::GFF->new(
	-file         => ">$outfile",
	-gff_version  => 3,
);

foreach (sort {$a cmp $b} keys %name2chr) {
	$gff->write_feature($name2chr{$_});
}

print " Done. Wrote file '$outfile'\n";













