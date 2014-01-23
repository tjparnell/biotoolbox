#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Tools::CodonTable;
use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure
	find_column_index
);
use Bio::ToolBox::db_helper qw(
	open_db_connection
);
use Bio::ToolBox::file_helper qw(
	open_tim_data_file 
	write_tim_data_file
);
my $VERSION = '1.14';



print "\n A script to locate SNPs and identify codon changes\n\n";



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
my (
	$database,
	$featurelist,
	$help,
	$print_version,
); # command line variables
my @infiles; 


# Command line options
GetOptions( 
	'in=s'        => \@infiles, # input file(s)
	'db=s'        => \$database, # the name of the database to use
	'features=s'  => \$featurelist, # list of features to look for
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
	print " Biotoolbox script locate_SNPs.pl, version $VERSION\n\n";
	exit;
}



### Check for general required values

# files
unless (@infiles) {
	# file list was provided on the command line
	@infiles = @ARGV or
		die " No input files! use --help for more information\n";
}
if (scalar @infiles == 1) {
	# only one file provided, but may be comma delimited list
	@infiles = split /,/, shift @infiles;
}

# database
unless ($database) {
	die " Must define an annotation database! use --help for more information\n";
}

# search features
my @search_types;
if ($featurelist) {
	@search_types = split /,/, $featurelist;
}
else {
	@search_types = qw(gene);
	print " Using default search feature of 'gene'\n";
}



### Initializations

# annotation database
my $db = open_db_connection($database) or
	die " unable to connect to database!\n";

# open codontable connection
my $codontable = Bio::Tools::CodonTable->new() or 
	die " unable to initialize Bio::Tools::CodonTable\n";



# Process each of the requested SNP files
foreach my $infile (@infiles) {
	
	## Open and process the SNP file
	print " processing SNP file '$infile'....\n";
	my ($fh, $metadata) = open_tim_data_file($infile) 
		or die "unable to open file!\n";
	
	# check input file
	my ($vcf_format, $genotype_exists);
	foreach ( @{ $metadata->{'other'} } ) {
		if (/^##fileformat=VCF/) {
			$vcf_format = 1;
		}
		elsif (/^##FORMAT=<ID=GT/) {
			$genotype_exists = 1;
		}
	}
	unless ($vcf_format) {
		warn " '$infile' is not a valid VCF format! skipping\n";
		next;
	}
	my $info_index = find_column_index($metadata, 'INFO');
	my $format_index = find_column_index($metadata, 'FORMAT');
	
	
	## Initialize data structure
	my $output = generate_tim_data_structure(
		'SNPs',
		qw(
			Variation_Type
			Overlapping_Feature
			Subfeature
			Codon_Change
			Chromosome
			Position
			Reference_Seq
			Variant_Seq
			Number_Supporting_Reads
			Total_Number_Reads
			Percent_Supporting_Reads
			Genotype
			Feature_description
		)
	) or die " unable to generate tim data structure!\n";
	$output->{'db'} = $database;
	my $table = $output->{'data_table'};
	
	
	## Walk through the file
	while (my $line = $fh->getline) {
		
		# collect the basic data
		chomp $line;
		my @data = split /\t/, $line;
		my ($chr, $pos, $ref, $snp) = @data[0,1,3,4];
		
		# collect the read depth
		my $total_depth = '0'; 
		my $supporting_depth = '0';
		foreach my $info ( split /;/, $data[$info_index] ) {
			# to find read depth and supporting read depth
			# we need to parse the INFO attributes
			if ($info =~ /^DP=(\d+)$/) {
				$total_depth = $1;
			}
			elsif ($info =~ /^DP4=(.+)/) {
				# samtools specific count for reads
				# comma separated ref-f, ref-r, alt-f, alt-r counts
				# not sure if this exists for other SNP callers
				# this may need editing for other formats
				my (undef, undef, $af, $ar) = split /,/, $1;
				$supporting_depth = $af + $ar;
			}
		}
		
		# collect the genotype
		my $genotype = '.'; # default null value
		if ($genotype_exists and $format_index) {
			
			# determine genotype index
			my $gt_index;
			my $i = 0;
			foreach (split /:/, $data[$format_index]) {
				if (/GT/) {
					$gt_index = $i;
					last;
				}
				$i++;
			}
			
			# pull out the genotype
				# we are using the sample immediately following the FORMAT
				# column, but there may be others in the file
			$genotype = ( split /:/, $data[$format_index + 1] )[$gt_index];
		}
		else {
			# no genotype present
			$genotype = '.'; # internal null
		}
		
		# determine the type of SNP
		my $snp_type;
		my @snps = split /,/, $snp;
			# SNP type is based only on the first variant seq, not subsequent
		if (length $ref == length $snps[0]) {
			$snp_type = 'substitution';
		}
		elsif (length $ref > length $snps[0]) {
			$snp_type = 'deletion';
		}
		elsif (length $ref < length $snps[0]) {
			$snp_type = 'insertion';
		}
		else {
			$snp_type = 'unknown';
		}
		
		# determine overlapping features
		my ($snp_feature, $snp_subfeature, $codon_change, $description) = 
			find_overlapping_features($chr, $pos, $snp, $snp_type);
		
		# calculate percent supporting depth
		my $percent = '0';
		if ($total_depth != 0 and $supporting_depth != 0) {
			$percent = sprintf "%.0f", ($supporting_depth / $total_depth) * 100;
		}
		
		# record for output
		push @{ $table }, [ 
			(
				$snp_type, 
				$snp_feature,
				$snp_subfeature,
				$codon_change, 
				$chr, 
				$pos, 
				$ref,
				$snp, 
				$supporting_depth, 
				$total_depth, 
				$percent, 
				$genotype,
				$description
			) 
		];

	}
	
	# finished with the file
	$fh->close;
	
	# re-sort the data table in decreasing order of confidence and position
	my $header = shift @{ $table }; # move the header line
	my @sorted_table = sort {
		$b->[10] <=> $a->[10] or # decreasing percentage
		$b->[8] <=> $a->[8] or # decreasing supporting depth
		$a->[4] cmp $b->[4] or # chromosome
		$a->[5] <=> $b->[5] # position
	} @{$table};
	unshift @sorted_table, $header;
	$output->{'data_table'} = \@sorted_table; # replace the table reference
	
	# update the number of records
	$output->{'last_row'} = scalar(@sorted_table) - 1;
		
	# write the output
	my $outfile = $infile;
	$outfile = $metadata->{'basename'};
	$outfile .= '_summary.txt';
	my $written_file = write_tim_data_file(
		'data'     => $output,
		'filename' => $outfile,
		'format'   => 'simple',
	);
	if ($written_file) {
		print " wrote file '$written_file'\n";
	}
	else {
		print " unable to write file!\n";
	}
}

# The end




############# Subroutines #########################

# find the feature, if any, that the SNP affects
sub find_overlapping_features {
	my ($chr, $pos, $snp, $snp_type) = @_;
	
	# establish a database segment based on position and length
	my $segment = $db->segment($chr, $pos, ($pos + length $snp - 1));
	unless (defined $segment) {
		warn "  unable to find segment for '$chr:$pos'!\n";
		return ('.', '.');
	}
	
	# collected information for the found features
	my @snp_features;
	my @snp_subfeatures;
	my @codon_changes;
	my $description;
	
	# identify overlapping features
	my @features = $segment->features(
			-type => \@search_types,
	);
	
	# Process features
	if (@features) {
		
		# we may have more than one overlapping feature
		# unlikely, but possible
		
		# walk through the the features
		foreach my $feature (@features) {
			
			# get info
			my $name = $feature->display_name;
			my $type = $feature->primary_tag;
			my @aliases = $feature->get_tag_values('Alias') if 
				$feature->has_tag('Alias');
			my @qual = $feature->get_tag_values('orf_classification') if 
				$feature->has_tag('orf_classification');
				# primarily for SGD features
			
			# generate name for the overlapping feature
			my $snp_feature = "$type $name";
			$snp_feature .= " (" . join(",", @aliases) . ")" if @aliases;
			if (@qual and $qual[0] =~ /dubious/i) {
				$snp_feature = "Dubious $snp_feature";
			}
			push @snp_features, $snp_feature;
			
			# get the description
			# except for dubious genes, or if we already have a description
			unless (@qual and $qual[0] =~ /dubious/i) {
				unless ($description) {
					$description = $feature->desc || 'None';
				}
			}
			
			# identify subfeatures
			my @subfeatures = $feature->get_SeqFeatures;
			if (@subfeatures) {
				
				# walk through 1st level subfeatures
				foreach my $subfeat (@subfeatures) {
					# the gene feature object may be comprised of subfeatures
					# we're looking for the CDS subfeature to look for 
					# codon changes
					
					# CDS Subfeature
					if (
						$subfeat->primary_tag eq 'CDS' and
						$subfeat->overlaps($segment)
					) {
						# we have found the specific CDS that overlaps our SNP
						
						push @snp_subfeatures, 'CDS';
						push @codon_changes, determine_codon_change(
							$snp_type,
							$snp,
							$subfeat,
							$segment,
						);
					}
					
					# RNA subfeature
					elsif ($subfeat->primary_tag =~ /rna/i) {
						# an mRNA, ncRNA, etc subfeature
						
						# get 2nd level subfeatures
						my @subfeatures2 = $subfeat->get_SeqFeatures;
						if (@subfeatures2) {
							
							# walk through each of the 2nd level subfeatures
							my $cds_found = 0;
							foreach my $subfeat2 (@subfeatures2) {
								# we may be looking at a mRNA feature
								# so go down one more level to look for CDS
								
								if (
									$subfeat2->primary_tag eq 'CDS' and
									$subfeat2->overlaps($segment)
								) {
									# we have found the CDS
									push @snp_subfeatures, 'CDS';
									
									# determine codon changes
									push @codon_changes, determine_codon_change(
										$snp_type,
										$snp,
										$subfeat,
										$segment,
									);
									
									# no need to go on
									$cds_found = 1;
									last;
								}
							}
							
							# no CDS found?
							unless ($cds_found) {
								if (scalar @subfeatures2 == 1) {
									# one subfeature2 found
									# use this one
									# could be UTR, or something like that
									push @snp_subfeatures, 
										$subfeatures2[0]->primary_tag;
									push @codon_changes, 'NA';
								}
								
								else {
									# use the first level subfeature instead
									push @snp_subfeatures, $subfeat->primary_tag;
									push @codon_changes, 'NA';
								}
							}
						}
						
						# no 2nd level subfeatures? maybe Intron?
						else {
							push @snp_subfeatures, 'Intron?';
							push @codon_changes, 'NA';
						}
					}
					
					# something else: UTR, non-coding exon, etc
					else {
						push @snp_subfeatures, $subfeat->primary_tag;
						push @codon_changes, 'NA';
					}
				}
			}
			
			# No subfeatures, but current feature is CDS
			elsif ($snp_feature =~ /^cds$/i) {
				# we're already in a CDS feature
				
				# no subfeature
				push @snp_subfeatures, 'None';
				
				# determine codon changes
				push @codon_changes, determine_codon_change(
					$snp_type,
					$snp,
					$feature,
					$segment,
				);
			}
			
			# No subfeatures, but current feature is gene/rna
			elsif ($snp_feature =~ /^gene|rna$/i) {
				# we're in gene or RNA feature, but can't find subfeatures
				# perhaps we're in an intron
				push @snp_subfeatures, 'Intron?';
				push @codon_changes, 'NA';
			}
			
			# No subfeatures at all, intergenic maybe?
			else {
				push @snp_subfeatures, 'None';
				push @codon_changes, 'NA';
			}
		}
	}
	
	# No features found!
	else {
		push @snp_features, 'None';
		push @snp_subfeatures, 'None';
		push @codon_changes, 'NA';
		$description = 'None';
	}
	
	# finished
	$description ||= 'None';
	return (
		join(',', @snp_features), 
		join(',', @snp_subfeatures), 
		join(',', @codon_changes),
		$description,
	);
}


sub determine_codon_change {
	# determine codon changes in the CDS
	
	# get passed arguments
	my ($snp_type, $given_snp, $feature, $segment) = @_; 
	
	# the return values
	my @codon_changes;
	
	# check each given SNP, there may be more than one
	# we're assuming that all the SNPs are the same type
		# may be bad assumption, but generally true
	foreach my $snp (split /,/, $given_snp) {
		# determine the type of change
		if ($snp_type eq 'substitution') {
			# this may be a real codon change
			
			# get more feature info
			my $phase = $feature->phase;
			my $start = $feature->start;
			my $stop = $feature->stop;
			
			# determine the SNP position
			my $pos = $segment->start;
			my $chr = $segment->seq_id;
			
			# collect the original codon
			my $codon_segment;
			my $pos_phase;
			my $codon;
			if ($feature->strand >= 0) {
				# watson or forward strand
				$pos_phase = ($pos - $start - $phase) % 3;
				if ($pos_phase == 0) {
					$codon_segment = $db->segment($chr, $pos, $pos + 2);
				} 
				elsif ($pos_phase == 1) {
					$codon_segment = $db->segment($chr, $pos - 1, $pos + 1);
				} 
				elsif ($pos_phase == 2) {
					$codon_segment = $db->segment($chr, $pos - 2, $pos);
				}
				
				# obtain the codon sequence
				$codon = $codon_segment->seq->seq;
			}
			else {
				# crick or reverse strand
				$pos_phase = ($stop - $pos - $phase) % 3;
				if ($pos_phase == 0) {
					$codon_segment = $db->segment($chr, $pos - 2, $pos);
				} 
				elsif ($pos_phase == 1) {
					$codon_segment = $db->segment($chr, $pos - 1, $pos + 1);
				} 
				elsif ($pos_phase == 2) {
					$codon_segment = $db->segment($chr, $pos, $pos + 2);
				}
				
				# obtain the reverse complement codon sequence
				$codon = $codon_segment->seq->revcom->seq;
			}
			
			# generate the corresponding translated product
			my $aa = $codontable->translate($codon);
			
			# make the substitution into the mutant
			if ($feature->strand < 0) {
				# don't forget to take the complement of the SNP
				# the SNP will always be the top strand base
				$snp =~ tr/agctAGCT/tcgaTCGA/;
				
			}
			my $mutant_codon = $codon;
			substr($mutant_codon, $pos_phase, 1, $snp);
			my $mutant_aa = $codontable->translate($mutant_codon);
			
			# report the change
			if ($aa eq $mutant_aa) {
				# no change
				push @codon_changes, 'silent';
			}
			else {
				# real change!
				push @codon_changes, "$aa->$mutant_aa";
			}
		}
		
		else {
			# insertion/deletion
			push @codon_changes, 'frameshift';
		}
	}
	
	return join(',', @codon_changes);
}





__END__

=head1 NAME

locate_SNPs.pl

A script to locate the position of SNPs and identify codon changes.

=head1 SYNOPSIS
 
 locate_SNPs.pl --db <database> <snp_file1> ...
  
  Options:
  --in <snp_file>
  --db <database>
  --features type1,type2,...
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <snp_file>

Specify the name(s) of the input SNP file(s). SNP files should be in the 
Variation Call Format (VCF) format. They may be gzipped. Multiple files 
may be specified; they will be processed sequentially. Provide multiple 
--in arguments, a comma delimited list, or a free list at the end of 
the command. 

=item --db <database>

Specify the name of a Bio::DB::SeqFeature::Store database that contains 
the genome annotation and sequence. Alternatively, for small genomes, 
a single GFF3 genome annotation file may be provided for loading into memory.

=item --features type1,type2,...

Provide a comma-delimited list (no spaces) of the GFF3 feature types 
of the features to intersect with the SNPs. Complex gene structures 
(gene->mRNA->CDS) can be parsed to look for amino-acid changes. The 
default feature is "gene". 

=item --version

Print the version number.

=item --help

This help text.

=back

=head1 DESCRIPTION

This program will locate SNPs and other sequence variants and identify 
the feature that overlaps the variant. In most cases, the features are 
genes, in which case the corresponding coding sequence (CDS), if 
present, is evaluated for a change in coding potential due to the 
sequence variation. Codon changes, silent changes, and frame shifts 
may be reported.

The input SNP files must be in the Variant Call Format (VCF) 4.0 or 
4.1 format, a tab-delimited text file with metadata. 
See L<http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41> 
for more information about the format. The files may be gzipped.

Numerous SNP callers are capable of generating the VCF format from 
sequence (usually Bam) files. The Samtools program is one such program, 
using the "mpileup" function in conjunction with it's "bcftools" tool. 
See the Samtools site at L<http://samtools.sourceforge.net> for more 
information.

The output file is a simple tab-delimited text file with headers, 
suitable for opening in spread-sheet program. Each row represents a SNP. 
The following columns are include.
  
  Variation_Type (substitution, insertion, deletion)
  Overlapping_Feature (feature type, name, alias)
  Subfeature (exon, CDS, etc)
  Codon_Change (Reference AA -> mutant AA)
  Chromosome
  Position (start position, 1-base, forward strand)
  Reference_Seq
  Variant_Seq
  Number_Supporting_Reads (if reported)
  Total_Number_Reads (if reported)
  Percent_Supporting_Reads (if reported)
  Genotype (if reported)
  Feature_description (if present in the database)

If more than one variant sequence is reported at a position, then each 
is evaluated for effects on the coding potential. If more than one 
feature is found overlapping the SNP, then both features are reported.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
