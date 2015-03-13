#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::big_helper qw(bed_to_bigbed_conversion);
use Bio::ToolBox::data_helper qw(format_with_commas);
use Bio::ToolBox::file_helper qw(open_to_write_fh);
eval {
	# check for bam support
	require Bio::ToolBox::db_helper::bam;
	Bio::ToolBox::db_helper::bam->import;
};
my $VERSION = '1.15';


print "\n A script to convert Bam alignments to GFF or BED files\n\n";


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
	$outfile,
	$paired_end,
	$gff,
	$bed,
	$bigbed,
	$type, 
	$source,
	$random_strand,
	$bb_app_path,
	$gz,
	$help,
	$print_version,
);
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
	'pe'         => \$paired_end, # bam is paired end reads
	'gff'        => \$gff, # output a GFF v.3 file
	'bed'        => \$bed, # ouput a BED file
	'bigbed|bb'  => \$bigbed, # generate a binary bigbed file
	'type=s'     => \$type, # the GFF type or method
	'source=s'   => \$source, # GFF source field
	'randstr'    => \$random_strand, # assign strand randomly
	'bbapp=s'    => \$bb_app_path, # path to bedToBigBed utility
	'gz!'        => \$gz, # gzip the output file
	'help'       => \$help, # request help
	'version'    => \$print_version, # print the version
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
	print " Biotoolbox script bam2gff_bed.pl, version $VERSION\n\n";
	exit;
}



### Check for required values and set defaults
# input file
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		die "  An input BAM file must be specified!\n";
	}
}
unless ($infile =~ /\.bam$/) {
	die " The input file must be a .bam file!\n";
}

# output type
if (($gff and $bed) or ($gff and $bigbed)) {
	die " Can not, must not, specify both GFF and BED simultaneously!\n";
}
unless ($gff or $bed or $bigbed) {
	print " Generating default BED formatted file\n";
	$bed = 1;
}
if ($bigbed) {
	$bed = 1;
}

# set GFF options
if ($gff) {
	unless ($type) {
		# derive the type from the input file name
		$type = $infile;
		$type =~ s/\.bam$//;
		$type =~ s/\.sorted$//;
		if ($paired_end) {
			$type .= '_paired_reads';
		}
		else {
			$type .= '_reads';
		}
	}
	unless ($source) {
		# reasonable, but biased, default
		$source = 'Illumina';
	}
}

# generate default output file name
unless ($outfile) {
	if ($gff) {
		# gff file uses the type
		$outfile = $type;
	}
	else { 
		# bed file simply reuses the input file name
		$outfile = $infile;
		$outfile =~ s/\.bam$//;
		$outfile =~ s/\.sorted$//;
		if ($paired_end) {
			$outfile .= '_paired_reads';
		}
		else {
			$outfile .= '_reads';
		}
	}
}

# default is no compression
unless (defined $gz) {
	$gz = 1;
}
if ($bigbed and $gz) {
	# enforce no compression if wanting to write bigbed file
	warn " compression not allowed when converting to BigBed format\n";
	$gz = 0; 
}






### Load the BAM file
print " Opening bam file....\n";
unless (exists &open_bam_db) {
	die " unable to load Bam file support! Is Bio::DB::Sam installed?\n"; 
}
my $sam = open_bam_db($infile) or die " unable to open bam file '$infile'!\n";




### Open the output file
# check extension
if ($gff) {
	unless ($outfile =~ /\.gff3?$/) {
		$outfile .= '.gff';
	}
}
else {
	unless ($outfile =~ /\.bed$/) {
		$outfile .= '.bed';
	}
}
if ($gz) {
	$outfile .= '.gz';
}

# open handle
my $out = open_to_write_fh($outfile, $gz) or 
		die " unable to open output file '$outfile'!\n";

# write metadata
if ($gff) {
	# print headers for GFF file
	$out->print("##gff-version 3\n");
	$out->print("# Program $0\n");
	$out->print("# Converted from source file $infile\n");
}
elsif ($bed and !$bigbed) {
	# bed file
	$out->print("# Program $0\n");
	$out->print("# Converted from source file $infile\n");
}


### Initialize global variables
# determine the appropriate callback sub for processing bam alignments
my $callback;
if ($paired_end) {
	$callback = \&process_paired_end;
}
else {
	$callback = \&process_single_end;
}

# determine appropriate output sub for writing data
my $write_feature;
if ($gff) {
	$write_feature = \&write_gff_feature;
}
else {
	$write_feature = \&write_bed_feature;
}

# Initialize output variables
my $data_count = 0;
my @output_data;
my %gff_ids; # a hash to make unique identifiers




### Process the reads
# Loop through the chromosomes
for my $tid (0 .. $sam->n_targets - 1) {
	# each chromosome is internally represented in the bam file as a numeric
	# target identifier
	# we can easily convert this to an actual sequence name
	# we will force the conversion to go one chromosome at a time
	
	# sequence name
	my $seq_id = $sam->target_name($tid);
	print " Converting reads on $seq_id...\n";
	
	# process the reads
	$sam->fetch($seq_id, $callback);
}

# Final write of output
incremental_write_data();

# Finished with output file
print " Wrote " . format_with_commas($data_count) . 
	" features to file '$outfile'.\n";
$out->close;




### Convert to BigBed format
if ($bed and $bigbed) {
	# requested to continue and generate a binary bigbed file
	print " converting to bigbed file....\n";
	
			
	# perform the conversion
	my $bb_file = bed_to_bigbed_conversion(
			'bed'       => $outfile,
			'db'        => $sam,
			'bbapppath' => $bb_app_path,
	);

	
	# confirm
	if ($bb_file) {
		print " BigBed file '$bb_file' generated\n";
		unlink $outfile; # remove the bed file
	}
	else {
		warn " BigBed file not generated! see standard error\n";
	}
	
}

print " Finished!\n";









########################   Subroutines   ###################################


### Process single end reads
sub process_single_end {
	my $a = shift;
	
	# check alignment
	return if $a->unmapped;
	
	# collect alignment data
	my $seq    = $a->seq_id;
	my $start  = $a->start;
	my $end    = $a->end;
	my $name   = $a->display_name;
	my $strand = $a->strand == 1 ? '+' : '-';
	my $score  = $a->qual;
	
	# write out the feature
	&{$write_feature}(
		$seq,
		$start,
		$end,
		$name,
		$score,
		$strand
	);
}



### Process paired end reads
sub process_paired_end {
	my $a = shift;
	
	# check alignment
	return if $a->unmapped;
	return unless $a->proper_pair;
	
	# we only need to process one of the two pairs, 
	# so only take the left (forward strand) read
	return unless $a->strand == 1;
	
	# collect alignment data
	my $seq    = $a->seq_id;
	my $start  = $a->start;
	my $name   = $a->display_name;
	my $score  = $a->qual;
	
	# calculate end
		# I occasionally get errors if I call mate_end method
		# just use the reported insert size listed in the original bam file
	my $end = $start + $a->isize - 1;
	
	# identify strand if possible
	my $strand;
	if ($random_strand) {
		if ( rand(1) > 0.5 ) {
			$strand = '+';
		}
		else {
			$strand = '-';
		}
	}
	else {
		# this is for paired-end RNA-Seq data aligned by TopHat, which records
		# the original strand in a BAM record attribute under the flag XS
		# this is either + or -
		
		# the default value will be +, as the proper paired-end fragments 
		# don't have inherent strand
		$strand = $a->get_tag_values('XS') || '+';
	}
	
	# write out the feature
	&{$write_feature}(
		$seq,
		$start,
		$end,
		$name,
		$score, 
		$strand
	);
}



### Write GFF features
sub write_gff_feature {
	
	# passed arguments
		# $seq,
		# $start,
		# $end,
		# $name,
		# $score,
		# $strand
	
	# generate unique ID
	my $id = $_[3]; # id will start off same as the name
	if (exists $gff_ids{$id}) {
		# must be duplicate alignment for this sequence read
		# we will need to make it unique
		# increment the counter by one
		$gff_ids{$id}++;
		# and then append this digit to the name to make a unique id
		$id = $id . '.' . $gff_ids{$id};
	}
	else {
		# record this name to check in the future
		$gff_ids{$id} = 0;
	}
	
	# generate GFF fields
	push @output_data, [ (
		$_[0],
		$source,
		$type,
		$_[1],
		$_[2],
		$_[4],
		$_[5],
		'.', # phase
		"ID=$id; Name=$_[3]"
	) ];
	
	# increment counter
	$data_count += 1;
	
	# periodically write data to file every 5000 lines
	if (scalar(@output_data) == 5000) {
		incremental_write_data();
	}
}


### Write BED features
sub write_bed_feature {
	
	# passed arguments
		# $seq,
		# $start,
		# $end,
		# $name,
		# $score,
		# $strand
	
	# generate BED fields, they're already in bed order
	# also convert to interbase format, per the specification
	push @output_data, [ (
		$_[0],
		$_[1] - 1,
		$_[2],
		$_[3],
		$_[4],
		$_[5]
	) ];
	
	# increment counter
	$data_count += 1;
	
	# periodically write data to file every 5000 lines
	if (scalar(@output_data) == 5000) {
		incremental_write_data();
	}
}



### Periodically write to output file
sub incremental_write_data {
	
	# write the output data
	while (@output_data) {
		my $data_ref = shift @output_data;
		
		# print line
		my $string = join("\t", @{ $data_ref }) . "\n";
		$out->print($string);
	}
}






__END__

=head1 NAME

bam2gff_bed.pl

A script to convert bam paired_reads to a gff or bed file.

=head1 SYNOPSIS

bam2gff_bed.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --bed | --gff | --bigbed | --bb
  --pe
  --type <text>
  --source <text>
  --randstr
  --out <filename> 
  --bbapp </path/to/bedToBigBed>
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a .bam alignment file. The file should be indexed. 
If not, the program will attempt to automatically index it.

=item --bed

=item --gff

=item --bigbed

=item --bb

Specify the output file format. Either a BED file (6-columns) or GFF v.3 
file will be written. You can also specify a bigBed format (bigbed or bb), 
where a BED text file will first be generated and then automatically 
converted to a compressed, binary BigBed file. The default behavior is to 
write a BED file.

=item --pe

Indicate that the BAM data is paired end, and that the generated feature 
should represent the insert fragment.

=item --type <text>

For GFF output files, specify the GFF type. The default is to use the 
input file base name, appended with either '_reads' or '_paired_reads'.

=item --source <text>

For GFF output files, specify the source tag. The default is 'Illumina'.

=item --randstr

For paired-end Bam files, specify that the strand of the alignment be 
randomly assigned to either the forward or reverse strand.

=item --out

Specify the output file name. The default for GFF output files is to use 
the GFF type; for BED output files the input file base name is used. A 
'.bed' extension is added if necessary.

=item --bbapp </path/to/bedToBigBed>

Specify the path to the Jim Kent's bedToBigBed conversion utility. The 
default is to first check the biotoolbox.cfg configuration file for 
the application path. Failing that, it will search the default 
environment path for the utility. If found, it will automatically 
execute the utility to convert the bed file.

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will convert alignments in a BAM formatted file into either a
BED or GFF formatted feature file. For BED output, a 6-column BED file is
generated. For GFF output, a v.3 GFF file is generated. In both cases, the
alignment name, strand, and a score value is recorded.

Both single- and paired-end alignments may be converted. For paired-end 
alignments, the genomic coordinates of the entire insert fragment is recorded. 

The mapping quality score is recorded as the GFF or BED score. For paired-end 
alignments, only the left score is recorded for efficiency.

The strand information is retained from the original alignment, except for 
proper paired-end alignments. For paired-end alignments, the TopHat-specific 
attribute XS is searched, and, if found, is used for the strand. Otherwise, 
all features default to the forward strand. An option is available to 
randomly assign a strand to paired-end features.

For BED files, coordinates are adjusted to interbase format, according to 
the specification.

An option exists to further convert the BED file to an indexed, binary BigBed 
format. Jim Kent's bedToBigBed conversion utility must be available.

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
