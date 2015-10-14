#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(open_db_connection);
use Bio::ToolBox::utility;
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};
my $BAM_OK;
eval { 
	# we want access to Bio::DB::Sam::Fai
	require Bio::DB::Sam;
	$BAM_OK = 1;
};
my $VERSION = '1.33';

print "\n This program will calculate observed & expected CpGs\n\n";

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
	$database,
	$window,
	$outfile,
	$gz,
	$cpu,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the input data file
	'db|fasta=s'=> \$database, # a SeqFeature::Store database or fasta file
	'win=i'     => \$window, # window size to take
	'out=s'     => \$outfile, # name of output file 
	'gz!'       => \$gz, # compress output
	'cpu=i'     => \$cpu, # number of execution threads
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
	print " Biotoolbox script CpG_calculator.pl, version $VERSION\n\n";
	exit;
}



### Check for requirements
unless ($window) {
	# default window size
	# could get from biotoolbox.cfg but I'm lazy right now
	$window = 1000;
	unless ($infile) {
		# print the obligatory default statement if no input file 
		print " using default window size of $window bp\n";
	}
}
if (!$infile and !$outfile) {
	die " must define an output file name!\n";
}
unless ($outfile) {
	$outfile = $infile;
}
unless (defined $gz) {
	$gz = 0;
}

# check parallel support
if ($parallel) {
	# conservatively enable 2 cores
	$cpu ||= 2;
}
else {
	# disable cores
	print " disabling parallel CPU execution, no support present\n" if $cpu;
	$cpu = 0;
}

my $start_time = time;
	



### Prepare the database and main data structure
my $Data;
if ($infile) {
	# an input file of regions is provided
	$Data = Bio::ToolBox::Data->new(file => $infile) or 
		die " unable to open input file '$infile'!\n";
	printf " Loaded %s features from $infile.\n", format_with_commas( $Data->last_row );
	
	# check database
	unless ($database) {
		$database = $Data->database or 
			die " no database or fasta file given! use --help for more information\n";
	}
}
else {
	# make a new genome list based on the type of database we're using
	if ($database) {
		$Data = Bio::ToolBox::Data->new(
			'feature' => 'genome',
			'db'      => $database,
			'win'     => $window,
		) or die " unable to generate genome window list!\n";
	}
	else {
		# no database, cannot continue
		die " no database or fasta file given! use --help for more information\n";
	}
	$Data->program("$0, v $VERSION");
}

# open database in case we need to do something, like index a fasta
my $db = open_sequence_db(0) or 
	die " unable to open database connection to '$database'!\n";



# check whether it is worth doing parallel execution
if ($cpu > 1) {
	while ($cpu > 1 and ($Data->last_row / $cpu) < 500) {
		# We need at least 1000 lines in each fork split to make 
		# it worthwhile to do the split, otherwise, reduce the number of 
		# splits to something more worthwhile
		$cpu--;
	}
}



### Process regions
print " Calculating CpG statistics....\n";
if ($cpu > 1) {
	# parallel execution
	print " Forking into $cpu children for parallel execution\n";
	parallel_execution();
}

else {
	# single threaded execution
	single_execution();
}



### Finished
# write the data file
my $written_file = $Data->write_file(
	'filename' => $outfile,
	'gz'       => $gz,
);
if ($written_file) {
	print " Wrote data file '$written_file' ";
}
else {
	print " unable to write data file! ";
}
printf " in %.2f minutes\n", (time - $start_time) / 60;




########################   Subroutines   ###################################


sub parallel_execution {
	my $pm = Parallel::ForkManager->new($cpu);
	$pm->run_on_start( sub { sleep 1; }); 
		# give a chance for child to start up and open databases, files, etc 
		# without creating race conditions
	
	# generate base name for child processes
	my $child_base_name = $outfile . ".$$"; 

	# Split the input data into parts and execute in parallel in separate forks
	for my $i (1 .. $cpu) {
		$pm->start and next;
	
		#### In child ####
	
		# splice the data structure
		$Data->splice_data($i, $cpu);
		
		# re-open database objects to make them clone safe
		# pass true to avoid cached database objects
		$db = open_sequence_db(1);
		
		# Collect the data
		process_regions();
		
		# write out result
		my $success = $Data->write_file(
			'filename' => sprintf("$child_base_name.%03s",$i),
			'gz'       => 0, # faster to write without compression
		);
		if ($success) {
			printf " wrote child file $success\n";
		}
		
		# Finished
		$pm->finish;
	}
	$pm->wait_all_children;
	
	# reassemble children files into output file
	my @files = glob "$child_base_name.*";
	unless (@files) {
		die "unable to find children files!\n";
	}
	unless (scalar @files == $cpu) {
		die "only found " . scalar(@files) . " child files when there should be $cpu!\n";
	}
	my $count = $Data->reload_children(@files);
	printf " reloaded %s features from children\n", format_with_commas($count);
}


sub single_execution {
	# execute
	process_regions();
}


sub process_regions {
	
	# Add new columns
	# Fraction gc
	my $fgc_i = $Data->add_column('Fraction_GC');
	
	# number of CpG
	my $cg_i = $Data->add_column('Number_CpG');
	
	# expected number of CpG
	my $exp_i = $Data->add_column('Expected_CpG');
	
	# observed/expected ratio
	my $oe_i = $Data->add_column('Obs_Exp_Ratio');
	
	
	# Process the regions
	$Data->iterate( sub {
		
		my $row = shift;
		# get the region subsequence
		# we are using our own database and not $Data's internal database connection
		# since we might be using the faster Bio::DB::Sam::Fai module
		my $seq = $db->seq(
			$row->seq_id,
			$row->start,
			$row->end + 1,
			# we add 1 bp so that we can count CpG that cross a window border
		) || undef;
		unless ($seq) {
			# no sequence, possibly because of difference in chromosome name convention?
			my $chrom = $row->seq_id;
			if ($chrom =~ /^chr(.+)$/i) {
				$chrom = $1;
			}
			else {
				$chrom = "chr$chrom";
			}
			# try again with different chromosome name
			$seq = $db->seq(
				$chrom,
				$row->start,
				$row->end + 1,
				# we add 1 bp so that we can count CpG that cross a window border
			) || undef;
		}
		unless ($seq) {
			# this may happen if 0 or >1 chromosomes match the name
			# or possibly coordinates are off the end, although I thought this was 
			# checked by the db adaptor
			my $w = sprintf 
				"No sequence available for segment %s:%s..%s at row %s, skipping\n", 
				$row->seq_id, $row->start, $row->end, $row->row_index;
			warn $w;
			
			# fill out null data
			$row->value($fgc_i, '.');
			$row->value($cg_i, '.');
			$row->value($exp_i, '.');
			$row->value($oe_i, '.');
			next;
		}
		
		# count frequencies
		# we could use the transliterate tr function to count single 
		# nucleotides quite efficiently, but it will NOT count dinucleotides
		# therefore, we will use the slightly more intensive approach of using 
		# substr to march through the sequence and count
		my $numC  = 0;
		my $numG  = 0;
		my $numCG = 0;
		for (my $i = 0; $i < length($seq) - 1; $i++) {
			my $dinuc = substr($seq, $i, 2);
			$numC  += 1 if $dinuc =~ m/^c/i;
			$numG  += 1 if $dinuc =~ m/^g/i;
			$numCG += 1 if $dinuc =~ m/^cg/i;
		}
		
		# record the statistics
		# we subtract 1 from the length because we added 1 when we generated the seq
		if (length($seq) > 1) {
			# must have reasonable length to avoid div by 0 errors
			$row->value($fgc_i, 
				sprintf "%.3f", ($numC + $numG) / (length($seq) - 1) );
		
			$row->value($cg_i, $numCG);
		
			$row->value($exp_i,  
				sprintf "%.0f", ($numC * $numG) / (length($seq) - 1) );
		
			$row->value($oe_i, 
				$row->value($exp_i) ? # avoid div by 0
				sprintf("%.3f", $numCG / $row->value($exp_i) ) : 0);
		}
		else {
			# a sequence of 1 bp? odd, just record default values
			$row->value($fgc_i, 0);
			$row->value($cg_i, $numCG); 
			$row->value($exp_i, 0);
			$row->value($oe_i, 0);
		}
	} );
	
}


# special sub to open a database or sequence file
# with option to open a genomic fasta file using a fast sequence accessor
sub open_sequence_db {
	my $nocache = shift || 0;
	my $db;
	print "Bam status is $BAM_OK\n";
	if ($database =~ /\.fa(?:sta)?$/i and $BAM_OK) {
		# this is a limited but very fast sequence accessor
		# based on samtools fasta index
		print "opening Fai db\n";
		$db = Bio::DB::Sam::Fai->open($database);
	}
	else {
		# otherwise we use a standard database connection
		print "opening ordinary db\n";
		$db = open_db_connection($database, $nocache);
	}
	return $db;
}



__END__

=head1 NAME

CpG_calculator.pl

A script to calculate observed vs expected CpG dinucleotides

=head1 SYNOPSIS

CpG_calculator.pl --fasta <directory|filename> [--options...]

CpG_calculator.pl --db <text> [--options...]
  
  Options:
  --db <name|file|directory>
  --fasta <file|directory>
  --in <filename>
  --win <integer>
  --out <filename> 
  --gz
  --cpu <integer>
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <name|file|directory>

=item --fasta <file|directory>

Provide the name of a Bio::DB::SeqFeature::Store database from which to 
collect the genomic sequence. Alternatively, provide the name 
of an uncompressed Fasta file (multi-fasta is ok) or directory containing 
multiple fasta files representing the genomic sequence. The directory 
must be writeable for a small index file to be written. For more information 
about using databases, see 
L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. 
The database may be provided in the metadata of an input file.

=item --in <filename>

Optionally specify an input file containing either a list of database features 
or genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. The 
first row should be column headers. Text files generated by other 
B<BioToolBox> scripts are acceptable. Files may be gzipped compressed.

=item --win <integer>

Optionally provide the window size in bp with which to scan the genome. 
Option is ignored if an input file is provided. Default is 1000 bp.

=item --out <filename>

Specify the output filename. By default it uses the input file base 
name if provided. Required if no input file is provided.

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --cpu <integer>

Specify the number of CPU cores to execute in parallel. This requires 
the installation of Parallel::ForkManager. With support enabled, the 
default is 2. Disable multi-threaded execution by setting to 1. 

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will calculate percent GC composition, number of CpG 
dinucleotide pairs, number of expected CpG dinucleotide pairs based 
on GC content, and the ratio of observed / expected CpG pairs. 
Calculations are performed on either windows across the entire genome 
(default behavior using 1000 bp windows) or user-provided regions in an 
input file (BED, GFF, or custom text file are supported). 

Genomic sequence may be provided in two ways. First, a Fasta file or 
directory of Fasta files may be provided. A small index file will be 
written to assist in random access using the Bio::DB::Fasta module. 
Alternatively, a Bio::DB::SeqFeature::Store database with sequence may 
be provided. Depending on the database driver and implementation, the 
fasta option is usually faster.

The four additional columns of information are appended to the input 
or generated file.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
