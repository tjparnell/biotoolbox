#!/usr/bin/perl

# a wrapper program for the Novo Aligner

use strict;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
use File::Copy;
use File::Basename qw(fileparse);


print "\n This script is a wrapper for the Novoaligner program\n\n";

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
	$novo_path,
	$sam_path,
	$repeat,
	$illumina,
	$tmpdir,
	$help,
); # command line variables



# Command line options
GetOptions( 
	'in=s'        => \@infiles, # input file(s)
	'index=s'     => \$index, # the novo index to use
	'novo=s'      => \$novo_path, # location of novoalign executable
	'sam=s'       => \$sam_path, # location of samtools executable
	'repeat=s'    => \$repeat, # which repeats to include
	'il!'         => \$illumina, # sequence is Illumina 1.3 format
	'temp=s'      => \$tmpdir, # the temp directory
	'help'        => \$help, # print the help
);


# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}



### Check for general required values

# files
if (@infiles) {
	if (scalar @infiles == 1) {
		# only one file provided, but may be comma delimited list
		my $file = shift @infiles;
		@infiles = split /,/, $file;
	}
}
else {
	# file list was provided on the command line
	@infiles = @ARGV or
		die "  OOPS! No source data files specified! \n use $0 --help\n";
}

# novo executable
unless ($novo_path) {
	$novo_path = `which novoalign` or 
		die " Unable to identify the novoalign executable!\n";
	chomp $novo_path;
}

# novoalign index
unless ($index and -e $index and -r $index) {
	die " Unable to verify novo index database!\n";
}

# repeats
unless ($repeat) {
	$repeat = 'None';
}

# samtools executable
unless ($sam_path) {
	$sam_path = `which samtools` or 
		die " Unable to identify the samtools executable!\n";
	chomp $sam_path;
}

# temporary directory
	# we'll be using a temporary directory to uncompress the source files and 
	# write the alignment files
	# this avoids automatic hourly backups and indexing of large uncompressed
	# sequence files
unless ($tmpdir) {
	$tmpdir = File::Spec->tmpdir();
}

# current directory
my $curdir = File::Spec->curdir();



### Perform the alignment(s)
foreach my $file (@infiles) {
	warn " processing '$file'...\n";
	
	# check the file name
	my ($basename, $path, $extension) = fileparse($file, 
		qw(.txt .txt.gz) );
	
	# move to the temp directory
	warn " ... copying sequence file to temp directory...\n";
	my $temp_seqfile = File::Spec->catfile($tmpdir, 
		$basename . $extension);
	copy($file, $temp_seqfile);
	
	# uncompress if necessary
	if ($extension eq '.txt.gz') {
		warn " ... uncompressing sequence file...\n";
		my $gunzip = `which gunzip` 
			or die " Unable to identify gunzip executable in path!\n";
		chomp $gunzip;
		system $gunzip, $temp_seqfile;
		$temp_seqfile =~ s/\.gz$//;
		unless (-e $temp_seqfile) {
			die " Unable to uncompress temp sequence file!\n" . 
				" command was '$gunzip $temp_seqfile.gz'\n";
		}
	}
	
	# perform the alignment
	my $temp_sam = File::Spec->catfile($tmpdir, $basename . '.sam');
	my $novo = "$novo_path -d $index -f $temp_seqfile -r $repeat -o SAM";
	if ($illumina) {
		$novo .= " -F ILMFQ";
	}
	$novo .= " >$temp_sam";
		# we need to put the entire novoalign command line into a single 
		# scalar value, this forces system to execute a sub shell first
		# the reason is because novoalign outputs to STDOUT, and we want 
		# to redirect it to a file
		# I guess we could capture it through perl and pass it to an 
		# open filehandle, but that seems messy.....
	warn " ... aligning...\n";
	if (system $novo) {;
		unlink $temp_seqfile;
		die " Unable to execute novoalign!\n command was '$novo'\n";
	}
	
	# convert to bam file
	my $temp_bam = File::Spec->catfile($tmpdir, $basename . '.bam');
	my @bam = ($sam_path, 'view', '-b', '-h', '-o', $temp_bam, '-S', $temp_sam);
	warn " ... converting to bam...\n";
	if (system @bam) {
		unlink $temp_seqfile;
		unlink $temp_sam;
		die " Unable to execute samtools!\n command was '" . 
			join(" ", @bam) . "'\n";
	}
	
	# sort the bam file
	warn " ... sorting the bam file...\n";
	my $temp_sbam = File::Spec->catfile($tmpdir, $basename . '.sorted');
	@bam = ($sam_path, 'sort', $temp_bam, $temp_sbam);
	system @bam;
	$temp_sbam .= '.bam'; # add the extension
	
	# index the file
	warn " ... indexing the bam file...\n";
	@bam = ($sam_path, 'index', $temp_sbam);
	my $temp_bai = $temp_sbam . '.bai'; # the index file name
	system @bam;
	
	# clean up
	move($temp_sbam, $curdir);
	move($temp_bai, $curdir);
	unlink $temp_seqfile;
	unlink $temp_sam;
	unlink $temp_bam;
	warn " ... finished.\n\n\n";
}





__END__

=head1 NAME

novo_wrapper.pl

=head1 SYNOPSIS
 
 novo_wrapper.pl --index <file> [--options] <seq_file1> ...
  
  Options:
  --in <seq_file>
  --index <novo_index>
  --novo </path/to/novoalign>
  --sam </path/to/samtools>
  --repeat [None|All|Random|Exhaustive]
  --(no)il
  --temp </path/to/temp/directory>
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <seq_file>

Specify the name(s) of the input sequence files. Multiple files may be 
specified using multiple --in arguments, with a single --in argument and a 
comma delimited list, or as a list at the end of all options. The files may 
kept gzipped; they will automatically be uncompressed.

=item --index <novo_index>

Specify the indexed genome database file to which the sequences will be 
aligned. The file should be built using Novocraft's 'novoindex' program.

=item --novo </path/to/novoalign>

Specify the full path to Novocraft's 'novoalign' executable file. This 
argument is not necessary if it is located in the path.

=item --sam </path/to/samtools>

Specify the full path to the 'samtools' executable file. This 
argument is not necessary if it is located in the path.

=item --repeat [None|All|Random|Exhaustive]

Specify the novoalign option of dealing with sequences which align to 
multiple locations in the genome (repeats). See the novoalign documenation 
for more information. Default value is None.

=item --(no)il

Indicate whether the input sequence file is (not) an Illumina sequence file 
from the Illumina Pipeline v. 1.3. See the novoalign documentation for more 
information. Default value is true.

=item --temp </path/to/temp/directory>

Specify a temp directory for writing uncompressed sequence files and 
intermediate files. By default it uses the temp directory specified in 
the user's environment.

=item --help

This help text.

=back

=head1 DESCRIPTION

This program is a simple wrapper for the novoalign and samtools programs, and 
performs the same functions as a custom shell script, including aligning the 
sequence reads to the genome using Novacraft's novoalign and then converting 
the alignments to an indexed, sorted, binary bam file using samtools, but 
perhaps in a slightly more friendly, if bloated, script.

The temp directory is used for temporary uncompressed sequence files and 
intermediate files to avoid unnecessary scheduled backups and indexing.


=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112
























