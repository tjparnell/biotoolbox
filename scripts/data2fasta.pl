#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use Bio::ToolBox::Data;
my $bio;
eval {
	require Bio::Seq;
	require Bio::SeqIO;
	$bio = 1;
};

my $VERSION = '1.63';

print "\n This program will convert a data file to fasta\n\n";

### Quick help
unless (@ARGV) {

	# when no command line options are present
	# print SYNOPSIS
	pod2usage(
		{
			'-verbose' => 0,
			'-exitval' => 1,
		}
	);
}

### Get command line options and initialize values
my (
	$infile, $outfile,  $database, $feature,     $subfeature,
	$id_i,   $desc_i,   $seq_i,    $chr_i,       $start_i,
	$stop_i, $strand_i, $extend,   $concatenate, $pad,
	$gz,     $help,     $print_version,
);

# Command line options
GetOptions(
	'i|in=s'          => \$infile,           # the solexa data file
	'o|out=s'         => \$outfile,          # name of output file
	'd|db|fasta=s'    => \$database,         # database name or genomic fasta file
	'f|feature=s'     => \$feature,          # feature from input file
	'u|subfeature=s'  => \$subfeature,       # subfeature to collect sequence
	'n|name|id=i'     => \$id_i,             # id index
	'desc=i'          => \$desc_i,           # description index
	's|seq=i'         => \$seq_i,            # sequence index
	'c|chr=i'         => \$chr_i,            # chromosome index
	'b|begin|start=i' => \$start_i,          # start index
	'e|stop|end=i'    => \$stop_i,           # stop index
	't|strand=i'      => \$strand_i,         # strand index
	'x|extend=i'      => \$extend,           # extend sequence by given bp
	'cat!'            => \$concatenate,      # concatenate sequences into one
	'pad=i'           => \$pad,              # pad concatenate sequences with given N bp
	'z|gz!'           => \$gz,               # compress output
	'h|help'          => \$help,             # request help
	'v|version'       => \$print_version,    # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {

	# print entire POD
	pod2usage(
		{
			'-verbose' => 2,
			'-exitval' => 1,
		}
	);
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
unless ($bio) {
	print <<END;
 This program requires Bio::Seq to ensure properly formatted fasta files.
 Please install the Bio::Perl package. You will need this to also obtain
 fasta indexing and retrieval modules, including Bio::DB::HTS, Bio::DB::Sam, 
 or Bio::DB::Fasta.
END
	exit;
}
unless ($infile) {
	$infile = shift @ARGV
		or die " no input file! use --help for more information\n";
}
unless ( defined $gz ) {
	$gz = 0;
}

### Open file ####
my $Data = Bio::ToolBox::Data->new(
	file       => $infile,
	parse      => 1,
	feature    => $feature,
	subfeature => $subfeature,
) or die " unable to load input file '$infile'\n";
unless ($database) {
	$database = $Data->database;
}

### Identify columns ####
unless ( defined $id_i ) {
	$id_i = $Data->name_column;
}
unless ( defined $seq_i ) {
	$seq_i = $Data->find_column('sequence');
}
unless ( defined $desc_i ) {
	$desc_i = $Data->find_column('description|note');
}
my $coords;
my $do_feature;
if ( defined $chr_i and defined $start_i and defined $stop_i ) {

	# user defined coordinates
	$coords = 1;
}
elsif ( $Data->feature_type eq 'coordinate' ) {

	# Input has coordinate columns
	$coords = 1;
}
elsif ( $Data->feature_type eq 'named' ) {

	# Input has named features that presumably have coordinates in a database
	$coords     = 1;
	$do_feature = 1;
}

### Determine mode ###
if ( defined $id_i and defined $seq_i and $concatenate ) {

	# sequence is already in the source file
	printf " Found Sequence column %s\n",    defined $seq_i  ? $seq_i  : '-';
	printf " Found Description column %s\n", defined $desc_i ? $desc_i : '-';
	print " writing a single concatenated fasta with the provided sequence\n";
	write_direct_single_fasta();
}
elsif ( defined $id_i and defined $seq_i ) {

	# sequence is already in the source file
	printf " Found Sequence column %s\n",    defined $seq_i  ? $seq_i  : '-';
	printf " Found Description column %s\n", defined $desc_i ? $desc_i : '-';
	print " writing a multi-fasta with the provided sequence\n";
	write_direct_multi_fasta();
}
elsif ( $coords and $concatenate ) {

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

	$Data->iterate(
		sub {
			my $row = shift;
			$concat_seq .= $row->value($seq_i);
			$concat_seq .= 'N' x $pad if $pad;
		}
	);

	# create final sequence object
	my $seq = Bio::Seq->new(
		-id   => $Data->basename,
		-desc => "Concatenated sequences",
		-seq  => $concat_seq,
	);

	# write out
	my $seq_io = open_output_fasta();
	$seq_io->write_seq($seq);
}

sub write_direct_multi_fasta {

	# write multi-fasta with the provided sequences
	my $seq_io = open_output_fasta();
	$Data->iterate(
		sub {
			my $row = shift;

			# create seq object
			my $seq = Bio::Seq->new(
				-id  => $row->value($id_i),
				-seq => $row->value($seq_i),
			);
			if ( defined $desc_i ) {
				$seq->desc( $row->value($desc_i) );
			}
			$seq_io->write_seq($seq);
		}
	);
}

sub fetch_seq_and_write_single_fasta {

	# fetch sequence from database and write concatenated fasta file

	# Open fasta database
	my $db;
	if ($database) {
		$db = $Data->open_new_database($database);
		unless ($db) {
			die " Could not open database '$database' to use!\n";
		}
	}
	elsif ( $Data->database ) {

		# cool, database defined in metadata, we'll use that
		# hope it works....
	}
	else {
		die " A sequence or fasta database must be provided to collect sequence!\n";
	}

	# collect concatenated sequences and write
	my $concat_seq;
	$Data->iterate(
		sub {
			my $row = shift;

			# make sure we parse and/or fetch the seqfeature if need be
			# this isn't necessarily automatic....
			my $f = $row->seqfeature if $do_feature;

			# collect provided arguments for generating sequence
			my @args;
			push @args, ( 'db',     $db )                    if $db;
			push @args, ( 'start',  $row->value($start_i) )  if defined $start_i;
			push @args, ( 'stop',   $row->value($stop_i) )   if defined $stop_i;
			push @args, ( 'seq_id', $row->value($chr_i) )    if defined $chr_i;
			push @args, ( 'strand', $row->value($strand_i) ) if defined $strand_i;
			push @args, ( 'extend', $extend )                if $extend;

			# collect sequence using provided arguments as necessary
			my $sequence = $row->get_sequence(@args);
			unless ($sequence) {
				printf "no sequence for line %d, %s", $row->line_number,
					$row->name || $row->coordinate;
				next;
			}

			# concatenate the sequence
			$concat_seq .= $sequence;
			$concat_seq .= 'N' x $pad if $pad;
		}
	);

	# create final sequence object
	my $seq = Bio::Seq->new(
		-id   => $Data->basename,
		-desc => "Concatenated sequences",
		-seq  => $concat_seq,
	);

	# write out
	my $seq_io = open_output_fasta();
	$seq_io->write_seq($seq);
}

sub fetch_seq_and_write_multi_fasta {

	# fetch sequence from database and write multi-fasta file

	# Open fasta database
	my $db;
	if ($database) {
		$db = $Data->open_new_database($database);
		unless ($db) {
			die " Could not open database '$database' to use!\n";
		}
	}
	elsif ( $Data->database ) {

		# cool, database defined in metadata, we'll use that
		# hope it works....
	}
	else {
		die " A sequence or fasta database must be provided to collect sequence!\n";
	}

	# open output file
	my $seq_io = open_output_fasta();

	# collect sequences and write
	$Data->iterate(
		sub {
			my $row = shift;

			# make sure we parse and/or fetch the seqfeature if need be
			# this isn't necessarily automatic....
			my $f = $row->seqfeature if $do_feature;

			# collect provided arguments for generating sequence
			my @args;
			push @args, ( 'subfeature', $subfeature )            if $subfeature;
			push @args, ( 'db',         $db )                    if $db;
			push @args, ( 'start',      $row->value($start_i) )  if defined $start_i;
			push @args, ( 'stop',       $row->value($stop_i) )   if defined $stop_i;
			push @args, ( 'seq_id',     $row->value($chr_i) )    if defined $chr_i;
			push @args, ( 'strand',     $row->value($strand_i) ) if defined $strand_i;
			push @args, ( 'extend',     $extend )                if $extend;

			# collect sequence based on values obtained above
			my $sequence = $row->get_sequence(@args);
			unless ($sequence) {
				printf "no sequence for line %d, %s", $row->line_number,
					$row->name || $row->coordinate;
				next;
			}

			# create seq object
			my $seq = Bio::Seq->new(
				-id  => $row->name || $row->coordinate,
				-seq => $sequence,
			);
			if ( defined $desc_i ) {
				$seq->desc( $row->value($desc_i) );
			}

			# write out
			$seq_io->write_seq($seq);
		}
	);
}

sub open_output_fasta {

	# get filename
	unless ($outfile) {
		$outfile = $Data->path . $Data->basename . '.fa';
	}
	unless ( $outfile =~ /\.fa(?:sta)?(?:\.gz)?/i ) {
		$outfile .= '.fasta';
	}
	if ( $gz and $outfile !~ /\.gz$/i ) {
		$outfile .= '.gz';
	}

	# open for writing
	my $out_fh = Bio::ToolBox::Data->open_to_write_fh( $outfile, $gz )
		or die "unable to open '$outfile' for writing!\n";

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

A program to retrieve sequences from a list of features

=head1 SYNOPSIS

data2fasta.pl [--options...] <filename>
  
  File Options:
  -i --in <filename>                input file: txt, gff, bed, ucsc, vcf, etc
  -o --out <filename>               output file name
  
  Database:
  -d --db <name|fasta>              annotation database with sequence or fasta
  
  Feature selection:
  -f --feature <text>               feature when parsing gff3, gtf, or ucsc input
  -u --subfeature [exon|cds|        collect over subfeatures 
        5p_utr|3p_utr] 
  
  Column indices:
  -n --name --id <index>            name or ID column
  -s --seq <index>                  column with sequence
  -c --chr <index>                  chromosome column
  -b --begin --start <index>        start coordinate column
  -e --end --stop <index>           stop coordinate column
  -t --strand <index>               strand column
  -x --extend <integer>             extend coordinates in both directions
  --desc <index>                    description column
  
  Fasta output options:
  --cat                             concatenate all sequences into one
  --pad <integer>                   pad concatenated sequences with Ns
  
  General options:
  -z --gz                           compress output fasta file
  -v --version                      print version and exit
  -h --help                         show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 File options

=over 4

=item --in E<lt>filenameE<gt>

Specify the input data file. The file may be a tab-delimited text file 
with coordinate columns for fetching genomic sequence. Alternatively it 
may be an annotation file such as GFF3, GTF, refFlat, genePred, etc, 
in which case sequence may be selected from certain features and 
subfeatures, for example mRNA and CDS. B<Note> that no collapsing of 
redundant or overlapping subfeatures is performed; see L<get_features.pl>. 
Finally, text files with sequence in a column, for example oligo sequences, 
may be used, skipping the need for database sequence retrieval.
The file may be compressed with gzip.

=item --out E<lt>filenameE<gt>

Specify the output filename. By default it uses the input file basename.

=back

=head2 Database

=over 4

=item --db E<lt>name|fastaE<gt>

Provide the name of an uncompressed Fasta file (multi-fasta is ok) or 
directory containing multiple fasta files representing the genomic 
sequence. If Bam file support is available, then the fasta will be 
indexed and searched with a fasta index .fai file. If not, then the 
fasta can by indexed by the older L<Bio::DB::Fasta> adapter, which 
also supports a directory of multiple fasta files. If the index is 
not present, then the parent directory must be writeable.
Alternatively, the name of a L<Bio::DB::SeqFeature::Store> 
annotation database that contains genomic sequence may also be provided. 
The database name may be obtained from the input file metadata. 
Required only if collecting sequence from genomic coordinates.

=back

=head2 Feature selection

=over 4

=item --feature E<lt>textE<gt>

When parsing a gene annotation file such as a GFF3, GTF, or UCSC format 
file, provide a feature type to select features if desired.

=item --subfeature [exon|cds|5p_utr|3p_utr]

When collecting from subfeatures, indicate the subfeature type from 
list available. No merging of overlapping or redundant subfeatures 
is performed here. See L<get_features.pl>.

=back

=head2 Column indices

=over 4

=item --name --id <column_index>

Optionally specify the index for the name or ID column. It may be 
automatically determined from the column header.

=item --seq E<lt>column_indexE<gt>

Optionally specify the index for the sequence column. It may be 
automatically determined from the column header.

=item --chr E<lt>column_indexE<gt>

Optionally specify the index for the chromosome column. It may be 
automatically determined from the column header.

=item --start E<lt>column_indexE<gt>

=item --begin E<lt>column_indexE<gt>

Optionally specify the index for the start position column. It may be 
automatically determined from the column header.

=item --stop E<lt>column_indexE<gt>

=item --end E<lt>column_indexE<gt>

Optionally specify the index for the stop position column. It may be 
automatically determined from the column header.

=item --strand E<lt>column_indexE<gt>

Optionally specify the index for the strand column. It may be 
automatically determined from the column header.

=item --extend E<lt>integerE<gt>

Optionally provide the number of extra base pairs to extend the start 
and stop positions. This will then include the given number of base 
pairs of flanking sequence from the database. This only applies when 
sequence is obtained from the database.

=item --desc E<lt>column_indexE<gt>

Optionally specify the index of the description column. It may be 
automatically determined from the column header.

=back

=head2 Fasta output options

=over 4

=item --cat

Optionally indicate that all of the sequences should be concatenated 
into a single Fasta sequence. The default is to write a multi-fasta 
file with separate sequences.

=item --pad E<lt>integerE<gt>

When concatenating sequences into a single Fasta sequence, optionally 
indicate the number of 'N' bases to insert between the individual 
sequences. The default is zero.

=back

=head2 General options

=over 4

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
be retrieved from a database. Multiple database adapters are supported for 
indexing genomic fastas, including the L<Bio::DB::HTS> package, the 
L<Bio::DB::Sam> package, or the BioPerl L<Bio::DB::Fasta> adapter. Annotation 
databases such as L<Bio::DB::SeqFeature::Store> are also supported.  
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
