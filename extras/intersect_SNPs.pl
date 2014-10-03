#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename qw(fileparse);
use Bio::ToolBox::file_helper qw(
	open_tim_data_file 
	write_tim_data_file 
);
use Bio::ToolBox::data_helper qw(generate_tim_data_structure);
my $VERSION = '1.15';

print "\n A script to intersect and pull unique SNPs from multiple lists\n\n";

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
	$gz,
	$help,
	$print_version,
);
my @infiles;

# Command line options
GetOptions( 
	'in=s'      => \@infiles, # the input data files
	'gz!'       => \$gz,
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
	print " Biotoolbox script <NAME>.pl, version $VERSION\n\n";
	exit;
}





### Check input files

# get the list of files
unless (@infiles) {
	if (@ARGV) {
		@infiles = @ARGV;
	}
	else {
		die " must provide input files!\n";
	}
}

# check for more than one 
if (scalar @infiles == 1) {
	if ($infiles[0] =~ /,/) {
		@infiles = split /,/, shift @infiles;
	}
	else {
		die " must provide more than one input file!\n";
	}
}




### Collect the SNP data from all the files

# initialize data hashes
my %snps; # structure: chr => pos => snp => basename = line
my %samples; # structure: basename => input_file metadata

# load the SNPs into the data hashes
load_snp_files();




### Sort the SNP data

# Prepare the common output structure
my $common = prepare_common_data();

# Sort through all the SNPs
print " Sorting SNPs....\n";
sort_snps();



### Write the data files

# Output the separate SNPs
foreach my $name (sort {$a cmp $b} keys %samples) {
	
	# report count
	print " There were ", $samples{$name}{'last_row'}, 
		" unique SNPs for '$name'\n";
	
	# print
	my $success = write_tim_data_file(
		'data'     => $samples{$name},
		'filename' => $name . '_unique.vcf',
		'gz'       => $gz,
	);
	if ($success) {
		print " Wrote file $success\n";
	}
	else {
		print " Failed to write file!\n";
	}
}

# write the common SNPs
{
	# report count
	print " There were ", $common->{'last_row'}, " common SNPs\n";
	
	# print
	my $success = write_tim_data_file(
		'data'     => $common,
		'filename' => join('_', keys %samples) . '_common.vcf',
		'gz'       => $gz,
	);
	if ($success) {
		print " Wrote file $success\n";
	}
	else {
		print " Failed to write file!\n";
	}
}



########################   Subroutines   ###################################

sub load_snp_files {
	
	foreach my $file (@infiles) {
		
		print " Loading '$file'... ";
		
		# open file
		my ($fh, $metadata) = open_tim_data_file($file) or 
			die " could not open file '$file'!\n";
		
		# check for VCF header
		my $vcf_format = 0;
		foreach ( @{ $metadata->{'other'} } ) {
			if (/^##fileformat=VCF/) {
				$vcf_format = 1;
				last;
			}
		}
		unless ($vcf_format) {
			warn " '$file' is not a valid VCF format! skipping\n";
			next;
		}
		
		# prep metadata
		$metadata->{'program'} = undef;
		push @{ $metadata->{'data_table'} }, $metadata->{'column_names'};
		
		# store the metadata for use later
		$samples{ $metadata->{'basename'} } = $metadata;
		
		# walk through file
		my $count = 0;
		while (my $line = $fh->getline) {
			chomp $line;
			
			# collect relevant information
			my ($chr, $pos, $snp) = (split /\t/, $line)[0,1,4];
			
			# store the snps
			$snps{$chr}{$pos}{$snp}{ $metadata->{'basename'} } = $line;
			
			$count++;
		}
		
		# finish
		$fh->close;
		print " $count SNPs\n";
	}
}



sub prepare_common_data {
	# generate an empty VCF structure to load the common SNPs
	
	# pick one input at random to copy the metadata
	my $name = (keys %samples)[0];
	
	# generate clean structure
	my $common = generate_tim_data_structure( undef, qw(
		#CHROM
		POS
		ID
		REF
		ALT
		QUAL
		FILTER
		INFO
		samples
	) );
	
	# add other metadata
	foreach ( @{ $samples{$name}->{'other'} } ) {
		# all metadata lines except FORMAT for the genotype
		unless (/^##FORMAT/) {
			push @{ $common->{'other'} }, $_;
		}
	}
	for my $i (0..8) {
		# add the AUTO metadata to prevent column metadata from being written
		$common->{$i}{'AUTO'} = 3;
	}
	
	# remove program information
	$common->{'program'} = undef;
	
	# done
	return $common;
}



sub sort_snps {
	
	foreach my $chr (sort { $a cmp $b } keys %snps) {
		# sort chromosomes asciibetically
		
		foreach my $pos (sort {$a <=> $b} keys %{ $snps{$chr} } ) {
			# sort by increasing position
			
			foreach my $snp (sort {$a cmp $b} keys %{ $snps{$chr}{$pos} } ) {
				# sort by the type of snp
				# in all liklihood this will only be one, but just in case
				
				# check how many samples have this SNP
				my @names = keys %{ $snps{$chr}{$pos}{$snp} };
				if (scalar( @names ) == 1) {
					# this is an unique SNP
					
					# break the line out into data array
					my @data = split /\t/, $snps{$chr}{$pos}{$snp}{$names[0]};
					
					# store it in the appropriate data table
					push @{ $samples{ $names[0] }->{'data_table'} }, \@data;
					$samples{ $names[0] }->{'last_row'} += 1;
				}
				
				else {
					# this is a common SNP
					
					# find the SNP with the greatest read depth and use it
					# as the example
					my $rd = 0;
					my $example;
					foreach my $n (@names) {
						if ($snps{$chr}{$pos}{$snp}{$n} =~ /DP=(\d+)/) {
							if ($1 > $rd) {
								$rd = $1;
								$example = $n;
							}
						}
					}
					unless ($example) {
						# read depth wasn't recorded!?
						# take the first one then
						$example = $names[0];
					}
					
					# generate array for the example data
					my @data = split /\t/, $snps{$chr}{$pos}{$snp}{$example};
					
					# remove extraneous columns to the required 8
					splice(@data, $#data, (8 - $#data));
					
					# add sample names
					push @data, join(',', @names);
					
					# store in the common data structure
					push @{ $common->{'data_table'} }, \@data;
					$common->{'last_row'} += 1;
				}
			}
		}
	}
}


__END__

=head1 NAME

intersect_SNPs.pl

A script to identify unique and common SNPs from multiple sequence runs.

=head1 SYNOPSIS
 
intersect_SNPs.pl <file1.vcf> <file2.vcf> ...
  
  Options:
  --in <filename>
  --gz
  --version
  --help
  
=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input SNP files. The files should be in the .vcf 
format, and may be gzipped. Each SNP file is assumed to contain 
one sample or strain only. 

=item --gz

Specify whether (or not) the output files should be compressed 
with gzip. 

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This simple program will identify unique and common Single Nucleotide 
Polymorphisms (SNPs) between two or more sequenced strains. This is 
useful, for example, in isolating unique SNPs that may be responsible 
for a mutant phenotype from background polymorphisms common to the 
strains. 

Each strain should have a separate SNP file in the Variant Call Format 
(VCF) 4.0 or 4.1 format, a tab-delimited text file with metadata. 
See L<http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41> 
for more information about the format. The files may be gzipped.

Numerous SNP callers are capable of generating the VCF format from 
sequence (usually Bam) files. The Samtools program is one such program, 
using the "mpileup" function in conjunction with it's "bcftools" tool. 
See the Samtools site at L<http://samtools.sourceforge.net> for more 
information.

Note that this program currently loads all SNPs into memory, thus for 
large genomes extensive memory requirements may be required.

SNPs are identified as unique vs common based on the reported coordinate 
and the alternate sequence. Overlapping SNPs will likely be treated 
separately. The unique SNPs are written to a new file with the file's 
base name appended with "_unique". The VCF format and headers are 
maintained. 

The common SNPs are written to a separate VCF file, with the file name 
comprised of input file base names appended with "_common". Genotype 
data, if present, is stripped from the common SNP, and only one 
representative is recorded.

Note that certain sequence variants may be reported as unique when 
in fact identical alternate sequences are also present in other strains. 
Usually in these cases, One of the strains has an additional variant 
sequence not present in the other.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
