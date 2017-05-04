#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;
use Bio::ToolBox::db_helper qw(open_db_connection);
use Bio::ToolBox::Data::Stream;
my $BAM_OK;
eval { 
	# we want access to Bio::DB::Sam::Fai
	require Bio::DB::Sam;
	$BAM_OK = 1;
};

my $VERSION =  '1.34';

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
	$interbase,
	$extend,
	$concatenate,
	$pad,
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
	'extend=i'  => \$extend, # extend sequence by given bp
	'cat!'      => \$concatenate, # concatenate sequences into one
	'pad=i'     => \$pad, # pad concatenate sequences with given N bp
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
	print " Biotoolbox script data2fasta.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
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
my $Input = Bio::ToolBox::Data::Stream->new(in => $infile) or 
	die "unable to open input file!\n";
unless ($database) {
	$database = $Input->database;
}



### Identify columns ####
unless (defined $id_i) {
	$id_i = $Input->name_column;
}
unless (defined $seq_i) {
	$seq_i = $Input->find_column('sequence');
}
unless (defined $desc_i) {
	$desc_i = $Input->find_column('description|note');
}
my $coords;
my $do_feature;
if (defined $chr_i and defined $start_i and defined $stop_i) {
	# user defined coordinates
	$coords = 1;
}
elsif ($Input->feature_type eq 'coordinate') {
	# Input has coordinate columns
	$coords = 1;
}
elsif ($Input->feature_type eq 'named') {
	# Input has named features that presumably have coordinates in a database
	$coords = 1;
	$do_feature = 1;
}
if (defined $start_i and substr($Input->name($start_i), -1) eq '0') {
	# name suggests $interbase
	$interbase = 1;
}
printf " Found ID column %s\n", defined $id_i ? $id_i : '-';
printf " Found Sequence column %s\n", defined $seq_i ? $seq_i : '-';
printf " Found Description column %s\n", defined $desc_i ? $desc_i : '-';


### Determine mode ###
if (defined $id_i and defined $seq_i and $concatenate) {
	# sequence is already in the source file
	print " writing a single concatenated fasta with the provided sequence\n";
	write_direct_single_fasta();
}
elsif (defined $id_i and defined $seq_i) {
	# sequence is already in the source file
	print " writing a multi-fasta with the provided sequence\n";
	write_direct_multi_fasta();
}
elsif ($coords and $concatenate) {
	# collect sequences and concatenate into single
	print " fetching sequence from $database and writing a concatenated fasta\n";
	fetch_seq_and_write_single_fasta();
}
elsif ($coords) {
	# need to collect sequence
	print " fetching sequence from $database and writing a multi-fasta\n";
	fetch_seq_and_write_multi_fasta();
}
else {
	die " unable to identify appropriate columns! see help\n";
}



### Finished
print " wrote file '$outfile'\n";


########################   Subroutines   ###################################

sub write_direct_single_fasta {
	# concatenate each of the provided sequences
	my $concat_seq;
	while (my $row = $Input->next_row) {
		$concat_seq .= $row->value($seq_i);
		$concat_seq .= 'N' x $pad if $pad;
	}
	
	# create final sequence object
	my $seq = Bio::Seq->new(
		-id     => $Input->basename,
		-desc   => "Concatenated sequences",
		-seq    => $concat_seq,
	);
	
	# write out
	my $seq_io = open_output_fasta();
	$seq_io->write_seq($seq);
}


sub write_direct_multi_fasta {
	# write multi-fasta with the provided sequences
	my $seq_io = open_output_fasta();
	while (my $row = $Input->next_row) {
		# create seq object
		my $seq = Bio::Seq->new(
			-id     => $row->value($id_i),
			-seq    => $row->value($seq_i),
		);
		if (defined $desc_i) {
			$seq->desc( $row->value($desc_i) );
		}
		$seq_io->write_seq($seq);
	}
}


sub fetch_seq_and_write_single_fasta {
	# fetch sequence from database and write concatenated fasta file
	
	# Open fasta database
	unless ($database) {
		die " Must provide a database or genomic fasta file(s)!\n";
	}
	my $db = open_sequence_db() or 
		die " Unable to open database '$database'!\n";
	
	# collect concatenated sequences and write
	my $concat_seq;
	while (my $row = $Input->next_row) {
		
		# collect sequence
		my ($sequence, $seq_id, $start, $stop) = fetch_sequence($row, $db);
		unless ($sequence) {
			printf "no sequence for $seq_id:$start..$stop! skipping\n";
			next;
		}
	
		# reverse if necessary
		if ($row->strand < 0) {
			# we will create a quick Bio::Seq object to do the reverse complementation
			my $seq = Bio::Seq->new(
				-id     => 'temp',
				-seq    => $sequence,
			);
			my $rev = $seq->revcom;
			$sequence = $rev->seq;
		}
		
		# concatenate the sequence
		$concat_seq .= $sequence;
		$concat_seq .= 'N' x $pad if $pad;
	}
	
	# create final sequence object
	my $seq = Bio::Seq->new(
		-id     => $Input->basename,
		-desc   => "Concatenated sequences",
		-seq    => $concat_seq,
	);
	
	# write out
	my $seq_io = open_output_fasta();
	$seq_io->write_seq($seq);
}



sub fetch_seq_and_write_multi_fasta {
	# fetch sequence from database and write multi-fasta file
	
	# Open fasta database
	unless ($database) {
		die " Must provide a database or genomic fasta file(s)!\n";
	}
	my $db = open_sequence_db() or 
		die " Unable to open database '$database'!\n";
	
	# open output file
	my $seq_io = open_output_fasta();
	
	# collect sequences and write
	while (my $row = $Input->next_row) {
		
		# collect sequence
		my ($sequence, $seq_id, $start, $stop) = fetch_sequence($row, $db);
		unless ($sequence) {
			print "no sequence for $seq_id:$start..$stop! skipping\n";
			next;
		}
	
		# create seq object
		my $seq = Bio::Seq->new(
			-id     => $row->name || "$seq_id:$start..$stop",
			-seq    => $sequence,
		);
		if (defined $desc_i) {
			$seq->desc( $row->value($desc_i) );
		}
		
		# reverse if necessary
		if ($row->strand < 0) {
			$seq = $seq->revcom();
		}
		
		# write out
		$seq_io->write_seq($seq);
	}
}


sub fetch_sequence {
	my ($row, $db) = @_;
	
	my $f = $row->feature if $do_feature;
	my $seq_id = defined $chr_i ? $row->value($chr_i) : $row->seq_id;
	my $start  = defined $start_i ? $row->value($start_i) : $row->start;
	my $stop   = defined $stop_i ? $row->value($stop_i) : $row->stop;
	$start += 1 if $interbase;
	if ($extend) {
		$start -= $extend;
		$start = 1 if $start < 0;
		$stop += $extend;
	}

	# collect sequence
	return ($db->seq($seq_id, $start, $stop), $seq_id, $start, $stop);
}


sub open_output_fasta {
	
	# get filename
	unless ($outfile) {
		$outfile = $Input->path . $Input->basename . '.fa';
	}
	unless ($outfile =~ /\.fa(?:sta)?(?:\.gz)?/i) {
		$outfile .= '.fasta';
	}
	if ($gz and $outfile !~ /\.gz$/i) {
		$outfile .= '.gz';
	}
	
	# open for writing
	my $out_fh = Bio::ToolBox::Data::Stream->open_to_write_fh($outfile, $gz) or 
		die "unable to open '$outfile' for writing!\n";
	
	# open SeqIO object
	my $seq_io = Bio::SeqIO->new(
		-fh     => $out_fh,
		-format => 'fasta',
	);
	return $seq_io;
}


sub open_sequence_db {
	my $db;
	if ($database =~ /\.fa(?:sta)?$/i and $BAM_OK) {
		# this is a limited but very fast sequence accessor
		$db = Bio::DB::Sam::Fai->open($database);
	}
	else {
		$db = open_db_connection($database);
	}
	return $db;
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
  --extend <integer>
  --zero
  --cat
  --pad <integer>
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

Provide the name of an uncompressed Fasta file (multi-fasta is ok) or 
directory containing multiple fasta files representing the genomic 
sequence. The directory must be writeable for a small index file to be 
written. If Bam support is available, the fasta file can be 
indexed with samtools, allowing for a faster fasta experience. 
Alternatively, the name of a Bio::DB::SeqFeature::Store 
annotation database that contains genomic sequence may be provided. 
The database name may be obtained from the input file metadata. 
Required only if collecting sequence from genomic coordinates.

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

=item --extend <integer>

Optionally provide the number of extra base pairs to extend the start 
and stop positions. This will then include the given number of base 
pairs of flanking sequence from the database. This only applies when 
sequence is obtained from the database.

=item --zero

Input file is in interbase or 0-based coordinates. This should be 
automatically detected for most known file formats, e.g. BED.

=item --cat

Optionally indicate that all of the sequences should be concatenated 
into a single Fasta sequence. The default is to write a multi-fasta 
file with separate sequences.

=item --pad <integer>

When concatenating sequences into a single Fasta sequence, optionally 
indicate the number of 'N' bases to insert between the individual 
sequences. The default is zero.

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

This program will take a tab-delimited text file (BED file, 
for example) and generate either a multi-sequence fasta file containing the 
sequences of each feature defined in the input file, or optionally a single 
concatenated fasta file. If concatenating, the individual sequences may be 
padded with the given number of 'N' bases. 

This program has two modes. If the name and sequence is already present in 
the file, it will generate the fasta file directly from the file content.

Alternatively, if only genomic position information (chromosome, start, 
stop, and optionally strand) is present in the file, then the sequence will 
be retrieved from a database, either a Bio::DB::SeqFeature::Store database, 
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
it under the terms of the Artistic License 2.0.  
