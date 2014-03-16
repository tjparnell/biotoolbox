#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
my $BAM_OK = 0;
eval { 
	require Bio::ToolBox::db_helper::bam;
	Bio::ToolBox::db_helper::bam->import;
	$BAM_OK = 1;
};
my $VERSION = '1.16';

print "\n A script to report the alignment sequence nucleotide frequencies\n\n";

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
	$help,
	$print_version,
);
my @size_list;
GetOptions( 
	'in=s'       => \$infile, # the input bam file path
	'out=s'      => \$outfile, # name of output file 
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
	print " Biotoolbox script get_bam_seq_stats.pl, version $VERSION\n\n";
	exit;
}



### Check for required values and set defaults
# input file
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		die " An input BAM file must be specified!\n";
	}
}
unless ($infile =~ /\.bam$/i) {
	die " Input file must be a .bam file!\n";
}
unless ($outfile) {
	$outfile = $infile;
	$outfile =~ s/\.bam$//;
	$outfile =~ s/\.sorted//;
	$outfile .= '.seq_stats.txt';
}
unless ($BAM_OK) {
	die "Unable to load Bam file support! Is Bio::DB::Sam installed?\n"; 
}


### Run the program
# open the bam object
my $sam = open_bam_db($infile) or die " unable to open input bam file '$infile'!\n";
	
my $seq_counter = get_bam_seq_stats::counter->new($sam);

$seq_counter->count_stats;

$seq_counter->print_stats($outfile);

$seq_counter->time_difference;

exit 0;






### Internal packages
package get_bam_seq_stats::counter;
use strict;
use IO::File;
1;

sub new {
	my $class = shift;
	my $sam  = shift;
	
	# return the object
	my $self = {
		'sam'      => $sam,
		'data'     => {},
		'no_align' => 0,
		'total'    => 0,
		'align'    => 0,
		'start_time'  => time,
	};
	return bless($self, $class);
}

sub count_stats {
	my $self = shift;
	my $sam = $self->{sam};
	
	# walk through each chromosome
	for my $tid (0 .. $sam->n_targets - 1) {
		print "  counting ", $sam->target_name($tid), "....\n";
		$sam->bam_index->fetch(
			$sam->bam, 
			$tid, 
			0, 
			$sam->target_len($tid), 
			\&callback,
			$self,
		);
	}
}


sub callback {
	# THIS IS NOT A CLASS METHOD!!!!!!
	# this is a callback subroutine for processing an alignment
	my ($a, $counter) = @_;
	
	# increase counters 
	$counter->{total}++;
	if ($a->unmapped) {
		$counter->{no_align}++;
		return;
	}
	$counter->{align}++;
	
	# record the stats on the alignment query sequence
	$counter->record_length($a->qseq);
}


sub record_length {
	my $self = shift;
	my $seq  = shift;
	$seq = lc $seq; # convert to lower case to make things easier
	my $len  = length $seq;
	$self->check_length($len);
	my @nucs = split //, $seq;
	for my $p (0 .. $#nucs) {
		if ($nucs[$p] eq 'a') {
			$self->{data}{$len}{$p}{a} += 1;
			next;
		}
		elsif ($nucs[$p] eq 'c') {
			$self->{data}{$len}{$p}{c} += 1;
			next;
		}
		elsif ($nucs[$p] eq 'g') {
			$self->{data}{$len}{$p}{g} += 1;
			next;
		}
		elsif ($nucs[$p] eq 't') {
			$self->{data}{$len}{$p}{t} += 1;
			next;
		}
		else {
			$self->{data}{$len}{$p}{n} += 1;
		}
	}
	$self->{data}{$len}{count} += 1;
}

sub check_length {
	my $self = shift;
	my $length = shift;
	return if exists $self->{data}{$length};
	$self->{data}{$length} = {
		'count' => 0,
	};
	for my $p (0 .. $length-1) {
		$self->{data}{$length}{$p} = {
			'a'     => 0,
			'c'     => 0,
			'g'     => 0,
			't'     => 0,
			'n'     => 0,
		};
	}
}

sub print_stats {
	my $self = shift;
	my $file = shift;
	
	my $fh = IO::File->new(">$file") or 
		die "unable to write output file $file!\n";
	
	# header
	# the bam file path is an undocumented component of the sam object
	$fh->print("# File ", $self->{sam}->{'bam_path'}, "\n");
	$fh->print("# ", $self->{'align'}, " total aligned reads\n");
	$fh->print("# ", $self->{'no_align'}, " total non-aligned reads\n\n");
	
	# data tables
	# there is a separate one for each size, one right after another, separated by space
	foreach my $size (sort {$a <=> $b} keys %{ $self->{data} } ) {
		my $count = $self->{data}{$size}{count};
		$fh->print("### Sequence length $size\n");
		$fh->print(
			"# $count (", 
			sprintf("%.2f%%", ($self->{data}{$size}{count} / $self->{align}) * 100 ), 
			") aligned reads\n"
		);
		$fh->print( join("\t", 'Nuc', (1 .. $size) ), "\n");
		foreach my $n (qw(a c g t n)) {
			$fh->print( 
				join("\t", $n, map { $self->{data}{$size}{$_}{$n} / $count } 
					(0 .. $size - 1) 
				), "\n"
			);
		}
		$fh->print("\n\n");
	}
}

sub time_difference {
	my $self = shift;
	printf " Completed in %.1f minutes\n", (time - $self->{start_time}) / 60;
}


__END__

=head1 NAME

get_bam_seq_stats.pl

A script to report the alignment sequence nucleotide frequencies.

=head1 SYNOPSIS

get_bam_seq_stats.pl <file.bam>
  
  Options:
  --in <file.bam>
  --out <filename>
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <file.bam>

Specify the file name of a binary Bam file of alignments as 
described for Samtools. It will be automatically indexed if 
necessary.

=item --out <filename>

Optionally specify the base name of the output file. The default is to use 
input base name. The output file names are appended with '.seq_stats.txt'. 

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will generate some statistics about the alignment 
sequences associated with a Bam file. This is using the the 
query sequence reported in the Bam file, not the genomic 
sequence or alignment. Only aligned sequences are analyzed.

The number and fraction of total for each length of the query 
sequences are reported. Additionally, the nucleotide composition 
for each position in the query sequences are also reported in 
a table, which should be suitable for generating a sequence logo, 
if desired.

The input file must be a BAM file as described by the Samtools 
project (http://samtools.sourceforge.net). 

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
