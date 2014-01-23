#!/usr/bin/env perl

# documentation at end of file

use strict;
use Pod::Usage;
use Getopt::Long;
use IO::Dir;
use File::Basename qw(fileparse);
use Bio::ToolBox::db_helper::config qw($BTB_CONFIG add_program);
use Bio::ToolBox::file_helper qw(
	open_to_read_fh
	open_to_write_fh
);
my $parallel;
eval {
	# check for parallel support
	require Parallel::ForkManager;
	$parallel = 1;
};

my $VERSION = '1.14';


print "\n This script is a wrapper for the Novoalign program\n\n";

### Quick help
unless (@ARGV) { # when no command line options are present
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Get command line options and initialize values

# Initialize values
my @infiles;
my (
	$index,
	$paired,
	$repeat,
	$split,
	$cores,
	$options,
	$novo_path,
	$sam_path,
	$help,
	$print_version,
); # command line variables



# Command line options
GetOptions( 
	'in=s'        => \@infiles, # input file(s)
	'index=s'     => \$index, # the novo index to use
	'pe!'         => \$paired, # paired end alignments
	'repeat=s'    => \$repeat, # which repeats to include
	'split!'      => \$split, # split the input file to fake
	'cores=i'     => \$cores, # number of cores to use
	'opt=s'       => \$options, # additional novoalign options
	'novo=s'      => \$novo_path, # location of novoalign executable
	'sam=s'       => \$sam_path, # location of samtools executable
	'help'        => \$help, # print the help
	'version'     => \$print_version, # print the version
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
	print " Biotoolbox script novo_wrapper.pl, version $VERSION\n\n";
	exit;
}




### Check for general required values
check_defaults();





### Perform the alignments
# record start time
my $start_time = time;

# list of files to clean up
my @to_delete;

if ($paired and $split) {
	run_paired_split_alignments();
}
elsif ($paired and !$split) {
	run_paired_alignments();
}
elsif ($split) {
	run_single_split_alignments();
}
else {
	run_single_alignments();
}


# finished
foreach my $f (@to_delete) {
	unlink $f;
}
printf " Finished in %.1f min\n", (time - $start_time)/60;







########################   Subroutines   ###################################

sub check_defaults {

	# files
	if (@infiles) {
		if (scalar @infiles == 1) {
			# only one file provided, but may be comma delimited list
			@infiles = split /,/, shift @infiles;
		}
	}
	else {
		# file list was provided on the command line
		@infiles = @ARGV or
			die "  OOPS! No source data files specified! \n use $0 --help\n";
	}
	
	# novoalign index
	unless ($index and -e $index and -r $index) {
		die " Unable to verify novo index database!\n";
	}
	
	# cores
	if ($split) {
		$split = 0 unless $parallel;
	}
	if ($cores) {
		if ($split) {
			$cores = 1 unless $parallel;
		}
	}
	else {
		if ($split and $parallel) {
			# running multiple instances
			# pretend everyone has a quad-core processor
			$cores = 4;
		}
		else {
			# no multiple processor
			$cores = 1;
		}
	}
	
	# repeats
	unless ($repeat) {
		$repeat = 'None';
	}
	
	
	# novo executable
	unless ($novo_path) {
		$novo_path = $BTB_CONFIG->param('applications.novoalign') || undef;
		
		# try the environment path
		unless ($novo_path) {
			eval {
				use File::Which;
				$novo_path = which('novoalign');
			};
			add_program($novo_path) if $novo_path; # remember for next time
		}
		
		# fail
		unless ($novo_path) {
			die " Unable to determine path to Novocraft's novoalign application!\n";
		}
	}
	
	# samtools executable
	unless ($sam_path) {
		$sam_path = $BTB_CONFIG->param('applications.samtools') || undef;
		
		# try the environment path
		unless ($sam_path) {
			$sam_path = `which samtools` || undef;
			chomp $sam_path;
		}
		
		# fail
		unless ($sam_path) {
			die " Unable to identify the samtools path!\n";
		}
	}
	
}


sub run_paired_split_alignments {
	die " paired alignments not implemented yet! Complain to Tim....\n";
}


sub run_paired_alignments {
	die " paired alignments not implemented yet! Complain to Tim....\n";
}



sub run_single_alignments {
	
	# run through input files
	while (@infiles) {
		
		# files
		my $file = shift @infiles;
		warn "#### Processing '$file'... ####\n";
		my ($basename, $path, $extension) = fileparse($file, 
			qw(.txt .txt.gz .fq .fq.gz .fastq .fastq.gz) );
		unless ($extension) {
			warn "#### WARNING: input file has no known extension!!!???? may not work\n";
		}
		
		
		# check if we need to decompress first
		if ($cores == 1 and $extension =~ m/^(.+)\.gz$/) {
			my $newfile = $path . $basename . $1;
			unless (-s $newfile) {
				warn " de-compressing file...\n";
				system('gunzip', '-c', $file, '>', $newfile);
				unless (-s $newfile) {
					die " decompressing $file to $newfile failed!\n";
				}
				push @to_delete, $newfile;
			}
			$file = $newfile;
		}
		
		# prepare the novo command
		my $outbam = $path . $basename . '.unsorted.bam';
		my $novo_command = "$novo_path -d $index -o SAM -r $repeat -s 2 -k" .
			" $options -f $file";
			# output to SAM
			# read trimming of 2 nt for unaligned reads
			# use base quality calibration
		
		# cpu cores
		if ($cores > 1) {
			$novo_command .= " -c $cores";
		}
		
		# piping to samtools
		$novo_command .= " | $sam_path view -bS - > $outbam";
		
		# executing alignment
		warn " #### Executing alignment: \"$novo_command\"\n";
		system($novo_command) == 0 or
			die " unable to execute novoalign: $!\n";
		unless (-s $outbam) {
			die "alignment failed!\n";
		}
		
		# post-process alignment
		post_process($path . $basename, $outbam);
	}
	
}


sub run_single_split_alignments {
	# run through input files
	while (@infiles) {
		
		# files
		my $file = shift @infiles;
		warn "#### Processing '$file' ####\n";
		my ($basename, $path, $extension) = fileparse($file, 
			qw(.txt .txt.gz .fq .fq.gz .fastq .fastq.gz) );
		
		# Split the files
		my $newfile = $path . $basename . '_';
		
		# check if we need to decompress first
		if ($extension =~ m/\.gz$/) {
			# run split through gunzip
			my $split_command = "gunzip -c $file | split -l 2000000 -a 2" .
				" - $newfile";
				# decompress input file
				# split into 2 million lines per file
				# use two letters for unique file suffix
			warn "#### Splitting compressed file\n";
			system($split_command) == 0 or 
				die " unable to decompress and split input file: $!\n";
		}
		else {
			# run split 
			my $split_command = "split -l 2000000 -a 2 $file $newfile";
				# split into 2 million lines per file
				# use two letters for unique file suffix
			warn "#### Splitting file\n";
			system($split_command) == 0 or 
				die " unable to split input file: $!\n";
		}
			
		# check split
		my $dir = IO::Dir->new($path) or
			die " unable to read directory $path!\n";
		my @split_files;
		while (my $f = $dir->read) {
			if ($f =~ m/$basename\_\w{2} \Z/x) {
				push @split_files, $f;
				push @to_delete, $f;
			}
		}
		unless (@split_files) {
			die " unable to find split files for $basename!\n";
		}
		
		# Prepare Parallel::ForkManager
		print " Forking into $cores children for parallel alignment\n";
		my $pm = Parallel::ForkManager->new($cores);
		
		foreach my $split_file (@split_files) {
			# run each split file in a separate fork
			$pm->start and next;
		
			### in child ###
			# prepare novoalign execution  
			my $novo_command = "$novo_path -d $index -o SAM -r $repeat -s 2 $options " . 
				"-f $split_file | $sam_path view -bS - > $split_file.unsorted.bam";
				# execute novoalign with index
				# output to SAM
				# read trimming of 2 nt for unaligned reads
				# pipe to samtools
			
			# executing alignment
			warn "#### Executing alignment: \"$novo_command\"\n";
			system($novo_command) == 0 or
				warn " unable to execute novoalign: $!\n";
			$pm->finish;
		}
		$pm->wait_all_children;
		
		
		# check for alignment files
		my @alignments;
		$dir->rewind;
		while (my $f = $dir->read) {
			if ($f =~m/ $basename \_ \w{2} \. unsorted \. bam \Z/x) {
				push @alignments, $f;
			}
		}
		unless (@alignments) {
			die " no output alignment files found!?\n";
		}
		
		# post-process alignment
		post_process($path . $basename, @alignments);
	}
}


sub post_process {
	my $basename = shift;
	my @files = @_;
	
	my $unsorted_file;
	if (scalar @files > 1) {
		$unsorted_file = $basename . '.merged.bam';
		warn "#### Merging Bam files\n";
		
		# execute cat
		system(
			$sam_path,
			'cat',
			'-o',
			$unsorted_file,
			@files
		) == 0 or die " unable to concatenate Bam files: $!\n";
		
		# check
		if (-s $unsorted_file) {
			push @to_delete, @files;
		}
		else {
			die " unable to concatenate split Bam files!\n";
		}
	}
	else {
		# no need to merge
		$unsorted_file = shift @files;
	}

	# sort bam file
	my $sorted_file = $basename . '.sorted'; # actually sorted file basename
	warn "#### Sorting Bam file\n";
	system(
		$sam_path,
		'sort',
		'-m',
		'4000000000',
		$unsorted_file,
		$sorted_file
	) == 0 or die " unable to sort file '$unsorted_file': $!\n";
	
	# check sorted file
	$sorted_file .= '.bam'; # samtools adds extension automatically
	if (-s $sorted_file) {
		push @to_delete, $unsorted_file;
	}
	else {
		die " could not sort $basename Bam file!\n";
	}
	
	# fix the header
	warn "#### Fixing the sort flag in Bam header\n";
	my $header_file = $basename . 'header.sam';
	system(
		$sam_path,
		'view',
		'-H',
		'-o',
		$header_file,
		$sorted_file,
	) == 0 or die "unable to export header for '$sorted_file': $!\n";
	
	my $final_file;
	if (-s $header_file) {
		
		# open files
		my $in = open_to_read_fh($header_file);
		my $fixed_header_file = $basename . 'fixed.header.sam';
		my $out = open_to_write_fh($fixed_header_file);
		
		# fix the Sort header
		while (my $line = $in->getline) {
			$line =~ s/SO:unsorted/SO:Coordinate/;
			$out->print($line);
		}
		$in->close;
		$out->close;
		
		# rehead
		$final_file = $basename . '.bam';
		system("$sam_path reheader $fixed_header_file $sorted_file" .
			" > $final_file") == 0 or 
			die " can't rehead $basename Bam file: $!\n";
		unless (-s $final_file) {
			die " could not rehead $basename Bam with fixed header!\n";
		}
		
		# clean up
		push @to_delete, (
			$sorted_file,
			$header_file,
			$fixed_header_file
		);
	}
	else {
		warn " could not fix Sort flag in Bam header!\n";
		$final_file = $sorted_file;
	}
	
	# Index file
	warn "#### Indexing final Bam file\n";
	system(
		$sam_path,
		'index',
		$final_file
	) == 0 or die " unable to index file '$final_file': $!\n";
	
	warn "#### Finished with final file '$final_file'\n\n\n";
}



__END__

=head1 NAME

novo_wrapper.pl

A parallelized wrapper program for Novocraft's novoaligner.

=head1 SYNOPSIS
 
 novo_wrapper.pl --index <file> [--options] <seq_file1> ...
  
  Options:
  --in <seq_file>
  --index </path/to/novoalign_index>
  --repeat [None | All | Random | Exhaustive]
  --split
  --cores <integer>
  --opt <"text">
  --novo </path/to/novoalign>
  --sam </path/to/samtools>
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <seq_file>

Specify the name(s) of the input Fastq sequence files. Multiple files may be 
specified using multiple --in arguments, with a single --in argument and a 
comma delimited list, or as a list at the end of all options. The files may 
kept gzipped; they will automatically be uncompressed.

=item --index <novo_index>

Specify the path to genome database index file to which the sequences will 
be aligned. The file should be built using Novocraft's 'novoindex' program.

=item --repeat [None | All | Random | Exhaustive]

Specify the novoalign option of dealing with sequences which align to 
multiple locations in the genome (repeats). See the Novoalign 
documentation for more information. Default value is None.

=item --split

Indicate that the input file should be split into numerous files (2 
million lines each) and multiple instances of Novoalign should be 
executed. This requires the GNU 'parallel' utility to be installed. 

=item --cores <integer>

When used with the --split option, this will limit the number of 
Novoalign instances executed at once. Default is 4.

When using a paid, licensed version of Novoalign, skip the --split 
option and specify as many CPU cores you want Novoalign to use in 
multi-threaded execution. Default is 1 (assume free academic version 
of Novoalign).

=item --opt <"text">

Specify any additional options you want to passed to the Novoalign 
executable. See the documentation for Novoalign for a full list of 
available options. Must enclose in quotes.

=item --novo </path/to/novoalign>

Optionally specify the full path to Novocraft's 'novoalign' executable 
file. It may also be specified in the 'biotoolbox.cfg' file, or it may 
be automatically identified by searching the environment path. 

=item --sam </path/to/samtools>

Optionally specify the full path to the 'samtools' executable 
file. It may also be specified in the 'biotoolbox.cfg' file, or it may 
be automatically identified by searching the environment path. 

=item --version

Print the program version number.

=item --help

This help text.

=back

=head1 DESCRIPTION

This program is a wrapper for Novocraft's Novoalign alignment tool. It also 
performs a number of post-alignment functions, including converting the 
alignments into a sorted, indexed, Bam file using samtools. 

When using the unlicensed, free, academic version of Novoalign that is 
limited to single-thread execution, novo_wrapper.pl can split the input 
file into numerous smaller versions and execute multiple instances of 
Novoalign. This requires the Parallel::ForkManager module for managing 
multiple executions. Note that only a limited number of instances 
should be run simultaneously; too many and your system may come to a 
standstill.

Licensed versions of Novoalign are multi-threaded and do not need to be 
split.

Novoalign is executed with the "-s 2", and "-r None" as default 
options. Additional custom options may be passed on to Novoalign by 
using the --opt flag above.

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
