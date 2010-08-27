#!/usr/bin/perl

# a script to identify and locate the position of SNPs

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Tools::CodonTable;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_db_helper qw(
	open_db_connection
);
use tim_file_helper;



print "A script to identify the location of SNPs\n";



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
	$help,
); # command line variables
my @infiles; 


# Command line options
GetOptions( 
	'in=s'        => \@infiles, # input file(s)
	'db=s'        => \$database, # the name of the database to use
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

# database
unless ($database) {
	die " Must define a database name or GFF3 annotation file!\n";
}




### Initializations

# annotation database
my $db = open_db_connection($database) or
	die " unable to connect to database!\n";

# open codontable connection
my $codontable = Bio::Tools::CodonTable->new() or 
	die " unable to load Codon Table! Requires Bio::Tools::CodonTable\n";



# Process each of the requested SNP files
foreach my $infile (@ARGV) {
	
	# initialize data structure
	my $output = initialize_output_data_structure($infile);
	
	my $table = $output->{'data_table'};
	
	# load and process the SNP file
	print " processing SNP file '$infile'....\n";
	my $fh = open_to_read_fh($infile) 
		or die "unable to open file!\n";
	while (my $line = $fh->getline) {
		
		# collect data
		chomp $line;
		my @data = split /\t/, $line;
		my ($chr, $pos, $ref, $snp, $cons_q, $snp_q, $number_reads) = 
			@data[0,1,2,3,4,5,7];
		
		# determine type
		my $snp_type;
		if ($snp =~ /\-/) {
			$snp_type = 'deletion';
		}
		elsif ($snp =~ /\+/) {
			$snp_type = 'insertion';
		}
		else {
			$snp_type = 'substitution';
		}
		
		# determine overlapping features
		my $snp_feature;
		my $codon_change;
		my $segment = $db->segment($chr, $pos, $pos + 1);
		my @search_types = qw(gene ncRNA snRNA snoRNA tRNA);
		my @features = $segment->features(
				-type => \@search_types,
		);
		if (@features) {
			my @collections;
			my @codon_changes;
			foreach my $feature (@features) {
				# get info
				my $name = $feature->name;
				my $type = $feature->primary_tag;
				my @aliases = $feature->get_tag_values('Alias') if 
					$feature->has_tag('Alias');
				my @qual = $feature->get_tag_values('orf_classification') if 
					$feature->has_tag('orf_classification');
				
				# generate name
				my $id;
				$id .= "$qual[0] " if (@qual);
				$id .= "$type $name";
				$id .= " (" . join(" ", @aliases) . ")" if @aliases;
				push @collections, $id;
				
				# identify subfeatures
				foreach my $subfeat ($feature->get_SeqFeatures) {
					# the gene feature object may be comprised of subfeatures
					# we're looking for the CDS subfeature to look for 
					# codon changes
					
					if (
						$subfeat->primary_tag eq 'CDS' and
						$subfeat->overlaps($segment)
					) {
						# we have found the specific CDS that overlaps our SNP
						
						my $codon_change = determine_codon_change(
							$snp_type,
							$snp,
							$subfeat,
							$segment,
						);
						if ($codon_change) {
							push @codon_changes, "$name: $codon_change";
						}
					}
					else {
						# we haven't found it yet
						
						# make sure this feature doesn't have more subfeatures
						foreach my $subfeat2 ($subfeat->get_SeqFeatures) {
							# we may be looking at a mRNA feature
							# so go down one more level to look for CDS
							
							if (
								$subfeat2->primary_tag eq 'CDS' and
								$subfeat2->overlaps($segment)
							) {
								# we have found the specific CDS that overlaps our SNP
								
								my $codon_change = determine_codon_change(
									$snp_type,
									$snp,
									$subfeat,
									$segment,
								);
								if ($codon_change) {
									push @codon_changes, "$name: $codon_change";
								}
							}
						}
					}
				}
				
				
			}
			$snp_feature = join ", ", @collections;
			$codon_change = join ", ", @codon_changes;
		}
		else {
			$snp_feature = 'no features';
		}
		
		# determine the number of supporting reads
		my $supporting_count = 0;
		if ($snp_type eq 'substitution') {
			if ($snp eq 'A') {
				$supporting_count = ($data[8] =~ tr/Aa//);
			}
			elsif ($snp eq 'T') {
				$supporting_count = ($data[8] =~ tr/Tt//);
			}
			elsif ($snp eq 'C') {
				$supporting_count = ($data[8] =~ tr/Cc//);
			}
			elsif ($snp eq 'G') {
				$supporting_count = ($data[8] =~ tr/Gg//);
			}
		}
		else { 
			# insertion/deletion
			$supporting_count = $data[10];
		}
		
		# record for output
		my $percent = sprintf "%.0f", ($supporting_count / $number_reads) * 100;
		push @{ $table }, [ 
			(
					$snp_type, 
					$snp_feature, 
					$codon_change, 
					$chr, 
					$pos, 
					$ref,
			 		$snp, 
			 		$supporting_count, 
			 		$number_reads, 
			 		$percent, 
			 		$cons_q, 
			 		$snp_q
			) 
		];

	}
	$fh->close;
	
	# re-sort the data table in decreasing order of confidence and position
	my $header = shift @{ $table }; # move the header line
	my @sorted_table = sort {
		$b->[7] <=> $a->[7] or # decreasing supporting count
		$b->[9] <=> $a->[9] or # decreasing percentage
		$a->[3] cmp $b->[3] or # chromosome
		$a->[4] <=> $b->[4] # position
	} @{$table};
	unshift @sorted_table, $header;
	$output->{'data_table'} = \@sorted_table; # replace the table reference
	
		
	# write the output
	my $outfile = $infile;
	$outfile =~ s/\.txt(?:\.gz)?$//; # strip extension
	$outfile .= '_summary.txt';
	my $written_file = write_tim_data_file( {
		'data'     => $output,
		'filename' => $outfile,
		'format'   => 'simple',
	} );
	if ($written_file) {
		print " wrote file '$written_file'\n";
	}
	else {
		print " unable to write file!\n";
	}
}



sub initialize_output_data_structure {
	my $filename = shift;
	
	# generate the data hash
	my %datahash;
	
	# populate the standard data hash keys
	$datahash{'program'}        = $0;
	$datahash{'feature'}        = 'SNPs';
	$datahash{'gff'}            = 0;
	$datahash{'number_columns'} = 12;
	
	# set column metadata
	$datahash{0} = {
		'name'     => 'Variation_Type',
		'index'    => 0,
	};
	$datahash{1} = {
		'name'     => 'Overlapping_Feature',
		'index'    => 1,
	};
	$datahash{2} = {
		'name'     => 'Codon_Change',
		'index'    => 2,
	};
	$datahash{3} = {
		'name'     => 'Chromosome',
		'index'    => 3,
	};
	$datahash{4} = {
		'name'      => 'Start',
		'index'     => 4,
	};
	$datahash{5} = {
		'name'     => 'Reference_Base',
		'index'    => 5,
	};
	$datahash{6} = {
		'name'     => 'Variation',
		'index'    => 6,
	};
	$datahash{7} = {
		'name'     => 'Number_Supporting_Reads',
		'index'    => 7,
	};
	$datahash{8} = {
		'name'     => 'Total_Number_Reads',
		'index'    => 8,
	};
	$datahash{9} = {
		'name'     => 'Percent_Supporting',
		'index'    => 9,
	};
	$datahash{10} = {
		'name'     => 'Consensus_Quality',
		'index'    => 10,
	};
	$datahash{11} = {
		'name'     => 'SNP_Quality',
		'index'    => 11,
	};
	
	# Set the data table
	my @data_table = ( [ qw(
		Variation_Type
		Overlapping_Feature
		Codon_Change
		Chromosome
		Start
		Reference_Base
		Variation
		Number_Supporting_Reads
		Total_Number_Reads
		Percent_Supporting
		Consensus_Quality
		SNP_Quality
	) ] );
	$datahash{'data_table'} = \@data_table;
	
	# return the reference to the generated data hash
	return \%datahash;
}



sub determine_codon_change {
	# determine codon changes in the CDS
	
	# get passed arguments
	my ($snp_type, $snp, $feature, $segment) = @_; 
	
	# the return value
	my $codon_change;
	
	# determine the type of change
	if ($snp_type eq 'substitution') {
		# this may be a real codon change
		
		# get more feature info
		my $strand = $feature->strand;
		my $phase = $feature->phase;
		my $start = $feature->start;
		my $stop = $feature->stop;
		
		# determine the SNP position
		my $pos = $segment->start;
		my $chr = $segment->seq_id;
		
		# collect the original codon
		my $codon_segment;
		my $pos_phase;
		if ($strand > 0) {
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
		}
		else {
			# crick or reverse strand
			$pos_phase = ($stop - $pos - $phase) % 3;
			if ($pos_phase == 0) {
				$codon_segment = $db->segment($chr, $pos + 2, $pos);
			} 
			elsif ($pos_phase == 1) {
				$codon_segment = $db->segment($chr, $pos + 1, $pos - 1);
			} 
			elsif ($pos_phase == 2) {
				$codon_segment = $db->segment($chr, $pos, $pos - 2);
			}
		}
		my $codon = $codon_segment->seq->seq;
		my $aa = $codontable->translate($codon);
		
		# make the substitution into the mutant
		my @triplet = split //, $codon;
		$triplet[$pos_phase] = $snp; # change the appropriate codon
		my $mutant_codon = join q(), @triplet;
		my $mutant_aa = $codontable->translate($mutant_codon);
		
		# report the change
		if ($aa eq $mutant_aa) {
			# no change
			$codon_change = 'silent';
		}
		else {
			# real change!
			$codon_change = "$aa->$mutant_aa";
		}
	}
	
	else {
		# insertion/deletion
		$codon_change = 'frameshift';
	}
	
	return $codon_change;
}





__END__

=head1 NAME

locate_SNPs.pl

=head1 SYNOPSIS
 
 locate_SNPs.pl --db <database> <snp_file1> ...
  
  Options:
  --in <snp_file>
  --db <database>
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <snp_file>

Specify the name(s) of the input SNP file. Multiple files may be 
specified using multiple --in arguments, with a single --in argument and a 
comma delimited list, or as a list at the end of all options. The files may 
kept gzipped; they will automatically be uncompressed.

=item --db

Specify the name of a single GFF3 genome annotation file. This is required.
Ideally, a database would be specified, but I am having difficulty accessing 
the child CDS features of genes to get at the coding frame to identify 
codon changes.

=item --help

This help text.

=back

=head1 DESCRIPTION

This program will locate SNPs and other sequence variants and identify 
which gene the variant is located in, the type of polymorphism (insertion, 
deletion, substitution), and, if the variant is located within a CDS, whether 
a codon change is generated.

The input files should be files generated using the varFilter function of
the samtools.pl script, included with the SamTools distribution
L<http://samtools.sourceforge.net>. That script will generate a list of the
sequence variations that differ from the reference genome. No verification
of the source file format is performed. The files may be gzipped.


=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112



















