#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;
use Bio::ToolBox::data_helper qw(find_column_index);
use Bio::ToolBox::db_helper qw(
	open_db_connection
);
use Bio::ToolBox::file_helper qw(
	open_tim_data_file
	open_to_write_fh
);

my $VERSION = '1.14';

print "\n This program will convert a data file to fasta\n\n";

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
	$database,
	$id_i,
	$desc_i,
	$seq_i,
	$chr_i,
	$start_i,
	$stop_i,
	$strand_i,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the solexa data file
	'out=s'     => \$outfile, # name of output file 
	'db|fasta=s'=> \$database, # database name or genomic fasta file
	'id=i'      => \$id_i, # id index
	'desc=i'    => \$desc_i, # description index
	'seq=i'     => \$seq_i, # sequence index
	'chr=i'     => \$chr_i, # chromosome index
	'start=i'   => \$start_i, # start index
	'stop|end=i' => \$stop_i, # stop index
	'strand=i'  => \$strand_i, # strand index
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
	print " Biotoolbox script data2fasta.pl, version $VERSION\n\n";
	exit;
}



### Check for requirements
unless ($infile) {
	$infile = shift @ARGV or
		die " no input file! use --help for more information\n";
}
unless (defined $gz) {
	$gz = 0;
}




### Open file ####
my ($in_fh, $metadata_ref) = open_tim_data_file($infile) or 
	die "unable to open input file!\n";
unless ($database) {
	$database = $metadata_ref->{'db'};
}



### Identify columns ####
unless (defined $id_i) {
	$id_i = find_column_index($metadata_ref, 'ID|Name');
}
unless (defined $seq_i) {
	$seq_i = find_column_index($metadata_ref, 'sequence');
}
unless (defined $desc_i) {
	$desc_i = find_column_index($metadata_ref, 'description|note');
}
unless (defined $chr_i) {
	$chr_i = find_column_index($metadata_ref, '^chr|seq|ref');
}
unless (defined $start_i) {
	$start_i = find_column_index($metadata_ref, '^start');
}
unless (defined $stop_i) {
	$stop_i = find_column_index($metadata_ref, '^stop|end');
}
unless (defined $strand_i) {
	$strand_i = find_column_index($metadata_ref, '^strand');
}



### Determine mode ###
if (defined $id_i and defined $seq_i) {
	# sequence is already in the source file
	write_direct_fasta();
}
elsif (defined $chr_i and $start_i and $stop_i) {
	# need to collect sequence
	fetch_seq_and_write_fasta();
}
else {
	die " unable to identify appropriate columns! see help\n";
}



### Finished
$in_fh->close;
print " wrote file '$outfile'\n";


########################   Subroutines   ###################################

sub write_direct_fasta {
	
	# open output file
	my $seq_io = open_output_fasta();
	
	while (my $line = $in_fh->getline) {
		chomp $line;
		my @data = split /\t/, $line;
		
		# create seq object
		my $seq = Bio::Seq->new(
			-id     => $data[$id_i],
			-seq    => $data[$seq_i],
		);
		
		if (defined $desc_i) {
			$seq->desc( $data[$desc_i] );
		}
		
		# print sequence
		$seq_io->write_seq($seq);
	}
}



sub fetch_seq_and_write_fasta {
	
	# Open fasta database
	unless ($database) {
		die " Must provide a database or genomic fasta file(s)!\n";
	}
	print " Indexing Fasta database for sequence collection....\n";
	my $db = open_db_connection($database) or 
		die " Unable to open database '$database'!\n";
	
	# open output file
	my $seq_io = open_output_fasta();
	
	# collect sequences and write
	while (my $line = $in_fh->getline) {
		chomp $line;
		my @data = split /\t/, $line;
		
		# coordinates
		my $seq_id = $data[$chr_i];
		my $start  = $data[$start_i];
		my $stop   = $data[$stop_i];
		
		# name
		my $name;
		if (defined $id_i) {
			$name = $data[$id_i];
		}
		else {
			$name = "$seq_id:$start..$stop";
		}
		
		# create seq object
		my $seq = Bio::Seq->new(
			-id     => $name,
			-seq    => $db->seq($seq_id, $start, $stop),
		);
		
		# description
		if (defined $desc_i) {
			$seq->desc( $data[$desc_i] );
		}
		
		# reverse if necessary
		if (defined $strand_i and $data[$strand_i] =~ /^[\-r]/) {
			$seq = $seq->revcom();
		}
		
		# write out
		$seq_io->write_seq($seq);
	}
}



sub open_output_fasta {
	
	# get filename
	unless ($outfile) {
		$outfile = $metadata_ref->{'basename'} . '.fa';
	}
	if ($gz and $outfile !~ /\.gz$/i) {
		$outfile .= '.gz';
	}
	
	# open for writing
	my $out_fh = open_to_write_fh($outfile, $gz) or 
		die "unable to open '$outfile' for writing!\n";
	
	# open SeqIO object
	my $seq_io = Bio::SeqIO->new(
		-fh     => $out_fh,
		-format => 'fasta',
	);
	
	return $seq_io;
}


__END__

=head1 NAME

data2fasta.pl

A script to retrieve sequences from a list of features

=head1 SYNOPSIS

data2fasta.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --db <name|file|directory>
  --id <index>
  --seq <index>
  --desc <index>
  --chr <index>
  --start <index>
  --stop <index>
  --strand <index>
  --out <filename> 
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input data file. The file should be a tab-delimited text file 
with columns representing the sequence id or name, sequence, description, 
chromosome, start, stop, and/or strand information. The file may be 
compressed with gzip.

=item --db <name|file|directory>

Provide the name of a Bio::DB::SeqFeature::Store database from which to 
collect the genomic sequence. A relational database name, SQLite file, or 
GFF3 file with sequence may be provided. Alternatively, provide the name 
of an uncompressed Fasta file (multi-fasta is ok) or directory containing 
multiple fasta files representing the genomic sequence. The directory 
must be writeable for a small index file to be written. The database 
name may be obtained from the input file metadata. Required.

=item --id <index>

Optionally specify the index for the name or ID column. It may be 
automatically determined from the column header.

=item --seq <index>

Optionally specify the index for the sequence column. It may be 
automatically determined from the column header.

=item --desc <index>

Optionally specify the index of the description column. It may be 
automatically determined from the column header.

=item --chr <index>

Optionally specify the index for the chromosome column. It may be 
automatically determined from the column header.

=item --start <index>

Optionally specify the index for the start position column. It may be 
automatically determined from the column header.

=item --stop <index>

Optionally specify the index for the stop position column. It may be 
automatically determined from the column header.

=item --strand <index>

Optionally specify the index for the strand column. It may be 
automatically determined from the column header.

=item --out <filename>

Specify the output filename. By default it uses the input file basename.

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will take a tab-delimited text file (tim data formated file, 
for example) and generate a multi-sequence fasta file containing the 
sequences of each feature defined in the input file. This program has 
two modes.

If the name and sequence is already present in the file, it will generate 
the fasta file directly from the file content.

If only genomic position information (chromosome, start, stop, and 
optionally strand) is present in the file, then the sequence will be 
retrieved from a database, either a Bio::DB::SeqFeature::Store database, 
a genomic sequence multi-fasta, or a directory of multiple fasta files. 
If strand information is provided, then the sequence reverse complement 
is returned for reverse strand coordinates.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
