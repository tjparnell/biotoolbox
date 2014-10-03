#!/usr/bin/perl

# documentation at end of file

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use Statistics::Lite qw(:all);
use File::Basename qw(fileparse);
use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure
	find_column_index
);
use Bio::ToolBox::db_helper qw(
	open_db_connection 
	verify_or_request_feature_types
	get_chromo_region_score
);
use Bio::ToolBox::file_helper qw(
	open_tim_data_file
	write_tim_data_file
	convert_and_write_to_gff_file
);
my $VERSION = '1.15';

print "\n This script will map transcription-enriched windows to gene transcripts\n\n";

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

# Initialize command line values
my (
	$database,
	$outfile,
	$fdata,
	$rdata,
	$ffile,
	$rfile,
	$win,
	$step,
	$threshold,
	$tolerance,
	$minfragment,
	$log,
	$gff,
	$source,
	$raw,
	$help,
	$print_version,
); 

# Command line options
GetOptions( 
	'db=s'        => \$database, # name of the database
	'out=s'       => \$outfile, # output file name
	'fdata=s'     => \$fdata, # name of the forward dataset in the db
	'rdata=s'     => \$rdata, # name of the reverse dataset in the db
	'ffile=s'     => \$ffile, # name of the foreward data file from 'find_enriched_regions.pl'
	'rfile=s'     => \$rfile, # name of the reverse data file from 'find_enriched_regions.pl'
	'win=i'       => \$win, # window size to find enriched regions
	'step=i'      => \$step, # step size to find enriched regions
	'threshold=f' => \$threshold, # threshold value to find enriched regions
	'tolerance=i' => \$tolerance, # tolerance value to call a transcript complete
	'min=i'       => \$minfragment, # minimum fragment to accept as transcript
	'log!'        => \$log, # log2 status of data
	'gff!'        => \$gff, # indicate a gff file should be written
	'source=s'    => \$source, # the gff source
	'raw'         => \$raw, # indicate a raw file should be written for debug purposes
	'help'        => \$help, # print help
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
	print " Biotoolbox script map_transcripts.pl, version $VERSION\n\n";
	exit;
}




### Check requirements
unless ($outfile) {
	die " Output file not defined!\n";
}
$outfile =~ s/\.txt$//; # remove pre-existing extension if exists

unless ($ffile and $rfile) {
	# we will need to generate window data if it's not specified
	# in which case we need parameters
	unless ($database and $win and defined $threshold and $fdata and $rdata) {
		die " One or more of the search parameters (db, win, step, thresh, " .
			"fdata, rdata) are missing!\n";
	}
	unless ($step) {
		# assign the default value of win
		$step = $win;
	}
}



### Set up data variables

# assign default values
unless ($tolerance) {
	# set default tolerance limit to 20
	$tolerance = 20;
} 
unless (defined $minfragment) {
	$minfragment = 100;  
	# the minimum size of a non-gene fragment to keep.
	# There are a lot of really short antisense and non-gene 
	# fragments (likely noise and non-specific signals) that 
	# I'm just trying to minimize, so simply imposing a hard
	# cutoff size
}
unless (defined $gff) {
	$gff = 1;
}
unless ($source) {
	$source = 'map_transcripts.pl';
}
unless (defined $log) {
	$log = 0;
}




### global values
my %transcripts; # a hash of all the transcription fragments (enriched windows)
my %sensetranscripts; # a hash of gene names with associated sense transfrags
my %antisensetranscripts; # a hash of gene names with associated antisense transfrags
my ($fdata_name, $rdata_name); # basename for the datasets

# open raw output file if asked
if ($raw) {
	open RAWFILE, ">$outfile\_raw.txt";
}




### Generate (or load) the transcript data
my ($fcount, $rcount) = prepare_regions();




### Open database connection
	# we do this here because the database could be specified in the 
	# provided $ffile metadata
my $db = open_db_connection($database) or 
	die " unable to open database connection!\n";





### Identify features associated with the transcription fragments 
print " Identifying associated features with the transcription fragments...\n";
identify_transcripts();

if ($raw) { # raw data dump of transcripts hash after identification 
	print RAWFILE "\n\n### Raw Transcript Hash dump ###\n";
	foreach my $cv (sort {$a <=> $b} keys %transcripts) {
		foreach my $sv (sort {$a <=> $b} keys %{ $transcripts{$cv} }) {
			my @data;
			push @data, "cv $cv";
			push @data, "sv $sv";
			foreach my $element (sort {$a cmp $b} keys %{ $transcripts{$cv}{$sv} }) {
				push @data, "$element $transcripts{$cv}{$sv}{$element}";
			}
			print RAWFILE join("; ", @data) . "\n";
		}
	}
}



### Collect transcripts from transcription fragments
print " Generating complete merged transcripts....\n";
if ($raw) {print RAWFILE "\n\n### Processing Transcripts ###\n"}
process_transcripts(\%sensetranscripts);
process_transcripts(\%antisensetranscripts);



### Determine UTRs and re-calculate extent of merged transcripts
print " Calculating UTRs....\n";
calculate_utr_and_extent();



### Prepare for output
# generate data structure
my $main_data_ref = convert_transcripts_to_tim_data_structure();
if (defined $main_data_ref) {
	# we no longer need the transcripts hash
	undef %transcripts;
}

# write main data file
{
	my $written_file = write_tim_data_file(
		'data'     => $main_data_ref,
		'filename' => $outfile,
	);
	if ($written_file) {
		print " Wrote data file '$written_file'\n";
	}
	else {
		print " unable to write data file!\n";
	}
}

# write gff file
if ($gff) {	
	my $gff_file = convert_and_write_to_gff_file(
		'data'     => $main_data_ref,
		'filename' => $outfile,
		'version'  => 3,
		'name'     => 10,
		'score'    => 11,
		'strand'   => 3,
		'source'   => $source,
		'type'     => 6,
		'tags'     => [ (5, 8, 12, 14) ] 
			# Tags: Transcript_Type, Extent, 5prime_UTR, 3prime_UTR
	);
	if ($gff_file) {
		print " Wrote GFF file '$gff_file'\n";
	}
	else {
		print " unable to write GFF file!\n";
	}
}

# Print the report to standard out
report_processing();


print "All done!\n";







######### Subroutines ########################################################


sub prepare_regions {
	# We will check the regions file
	# if none was provided, then execute find_enriched_regions.pl
	
	# Forward strand
	if ($ffile) {
		# previously run enriched regions file
		print " Using forward strand enriched regions data file '$ffile'\n";
	} 
	else {
		# generate enriched regions
		print " No forward strand data file specified. Generating new data by\n" .
			"  running 'find_enriched_regions.pl' script...\n";
		
		# need to generate the file name
		unless ($fdata_name) {
			($fdata, $fdata_name) = check_dataset($fdata);
		}
		$ffile = "$fdata_name\_w$win\_s$step\_t$threshold\_f.txt";
		
		# find_enriched_regions.pl should be found at the same location
		my $find_enriched_program = $Bin . '/find_enriched_regions.pl';
		$find_enriched_program .= " --db $database --out $ffile --win $win" . 
			" --step $step --nofeat --thresh $threshold --method median" . 
			" --trim --strand f";
		if ($fdata =~ /^file:(.+)$/) {
			$find_enriched_program .= " --data $1";
		}
		else {
			$find_enriched_program .= " --data $fdata";
		}
		$find_enriched_program .= ' --log' if $log;
			
		# execute
		!system "$find_enriched_program" or
				die "can't launch '$find_enriched_program'\n"; 
				# successful run of program will return 0, not 1, 
				# so negate system with !
	}
	# Reverse strand
	if ($rfile) {
		# previously run enriched regions file
		print "   using reverse strand enriched regions data file '$rfile'\n";
	} 
	else {
		# generate enriched regions
		print " No reverse strand data file specified. Generating new data by\n" .
			"  running find_enriched_regions.pl script...";
		
		# need to generate the file name
		unless ($rdata_name) {
			($rdata, $rdata_name) = check_dataset($rdata);
		}
		$rfile = "$rdata_name\_w$win\_s$step\_t$threshold\_r.txt";
		
		# find_enriched_regions.pl should be found at the same location
		my $find_enriched_program = $Bin . '/find_enriched_regions.pl';
		$find_enriched_program .= " --db $database --out $rfile --win $win" . 
			" --step $step --nofeat --thresh $threshold --method median" . 
			" --trim --strand r";
		if ($rdata =~ /^file:(.+)$/) {
			$find_enriched_program .= " --data $1";
		}
		else {
			$find_enriched_program .= " --data $rdata";
		}
		$find_enriched_program .= ' --log' if $log;
		
		# execute
		!system "$find_enriched_program " or
				die "can't launch '$find_enriched_program'\n";
				# successful run of program will return 0, not 1, 
				# so negate system with !
	}
	
	# Loading the enriched regions transcript data
	print " Loading enriched regions transcript data...\n";
	my %chromosomes; # for remembering chromosome names
	$chromosomes{'dummy_no_name_chromosome'} = 0;
	
	# load the f regions and return the count
	my $fcount = load_regions($ffile, 1, \%chromosomes); 
	print "  loaded $fcount forward transcripts\n";
	
	# load the r regions and return the count
	my $rcount = load_regions($rfile, -1, \%chromosomes); 
	print "  loaded $rcount forward transcripts\n";
	
	# finished loading
	return ($fcount, $rcount);
}

sub check_dataset {
	my $data = shift;
	
	# a remote file
	if ($data =~ /^(?:http|ftp):\/\/(.+)$/) {
		# no test for validity, we'll find out later when we try to use it
		# use the file's basename as the dataset name
		my ($basename, undef, undef) = fileparse($1, qw(\.bb \.bw \.bam));
		return ($data, $basename);
	}
	
	# local file
	elsif ($data =~ /^file:(.+)$/) {
		# a local file
		unless (-e $1) {
			die " requested dataset file $1 does not exist!\n";
		}
		my ($basename, undef, undef) = fileparse($data, 
			qw(\.bb \.bw \.bam));
		return ($data, $basename);
	}
	
	# local file or database feature
	elsif ($data =~ /\.(?:bw|bb|bam)$/i) {
		# check if we have a file 
		if (-e $data) {
			# file exists
			my ($basename, undef, undef) = fileparse($data, 
				qw(\.bb \.bw \.bam));
			return ("file:$data", $basename);
		}
		else {
			# maybe it's a funny named dataset?
			# verify with the database
			$data = verify_or_request_feature_types(
				'db'      => $db,
				'feature' => $data,
			);
			unless ($data) {
				# if it wasn't returned, it is not valid
				die " The requested file or dataset '$data' " . 
					"neither exists or is valid!\n";
			}
			return ($data, $data);
		}
	}
}


sub load_regions {
	my ($filename, $strand, $chromosomes) = @_;
	print " Loading regions from file '$filename'....\n";
	
	# First open the enriched windows file
	my $count = 0;
	my ($fh, $metadata_ref) = open_tim_data_file($filename) or 
		die " unable to open regions file '$filename!'\n";
	
	# check the feature of the loaded file
	unless ($metadata_ref->{'feature'} eq 'enriched_regions') {
		die " file does not have features of 'enriched_regions'!\n";
	}
	
	
	# Determine indices
	my $id_idx = find_column_index($metadata_ref, '^window|id') || 0;
	my $chr_idx = find_column_index($metadata_ref, '^chr|seq') || 1;
	my $start_idx = find_column_index($metadata_ref, '^start') || 2;
	my $stop_idx = find_column_index($metadata_ref, '^stop|end') || 3;
	my $size_idx = find_column_index($metadata_ref, '^size|length') || 4;
	my $score_idx = find_column_index($metadata_ref, '^score') or 
		die " unable to identify score column in '$filename'!\n";
	unless (defined $chr_idx and defined $start_idx and defined $stop_idx) {
		die " unable to identify coordinate columns in '$filename'!\n";
	}
	
	# Update missing information from the metadata
		# these references are hard coded because they should be coming from 
		# 'find_enriched_regions.pl'. Woe is me if I ever change that 
		# program!
	if (!$fdata and $strand == 1) {
		# forward dataset
		$fdata = $metadata_ref->{$score_idx}{dataset};
		unless ($fdata_name) {
			($fdata, $fdata_name) = check_dataset($fdata);
		}
	}
	if (!$rdata and $strand == -1) {
		# reverse dataset
		$rdata = $metadata_ref->{$score_idx}{dataset};
		unless ($rdata_name) {
			($rdata, $rdata_name) = check_dataset($rdata);
		}
	}
	unless (defined $database) {
		# define database
		$database = $metadata_ref->{db};
	}
	if (defined $threshold) {
		# check that thresholds match
		unless ($metadata_ref->{$score_idx}{threshold} == $threshold) {
			warn " whoa! the threshold values don't match between this program and file!\n";
		}
	}
	else {
		# define theshold
		$threshold = $metadata_ref->{$score_idx}{threshold};
	}
	
	
	
	# Load the regions data into the transcripts hash
	while (my $line = $fh->getline) {
		chomp $line;
		my @line_data = split /\t/, $line;
		
		# identify the chromosome and value
		# to simplify sorting later on, we will assign the chromosome 
		# name a simple number
		unless (exists $chromosomes->{ $line_data[$chr_idx] } ) {
			# generate a new number for this chromosome
			$chromosomes->{ $line_data[$chr_idx] } = 
				max( values %{$chromosomes} ) + 1;
		}
		my $cv = $chromosomes->{ $line_data[$chr_idx] };
		
		my $sv = $line_data[$start_idx]; # start value
		while (exists $transcripts{$cv}{$sv}) {
			# check that the start position is unique
			$sv += 0.01; # if not, make it unique by adding a small value
		}
		$transcripts{$cv}{$sv} = {
				# these are hardcoded because they're defined in the program
				# 'find_enriched_regions.pl'
				'transfrag' => $line_data[$id_idx],
				'chromo'    => $line_data[$chr_idx],
				'start'     => $line_data[$start_idx],
				'end'       => $line_data[$stop_idx],
				'size'      => $line_data[$size_idx],
				'score'     => $line_data[$score_idx],
				'strand'    => $strand,
		};
		$count++;
	}
	$fh->close;
	return $count;
}




sub identify_transcripts {
	# Next look for transcripts
	
	
	# chromosome value
	foreach my $cv (sort { $a <=> $b } keys %transcripts) { 
		# we will first sort by increasing chromosome number at the first 
		# level of hashes
		print "   ...on chromosome $cv\n";
		
		
		# start value
		foreach my $sv (sort { $a <=> $b } keys %{ $transcripts{$cv} } ) { 
			# next we will sort by increasing start position of the enriched windows
			my $strand = $transcripts{$cv}{$sv}{strand}; 
			
			# Collect known ORFs in the transcript fragment region
			my $region = $db->segment(
					-name   => $transcripts{$cv}{$sv}{chromo},
					-start  => $transcripts{$cv}{$sv}{start},
					-end    => $transcripts{$cv}{$sv}{end},
			);
			my @orflist;
			foreach my $f ( $region->features(-types => 'gene') ) {
				# collect all the genes in the region
				# check that they are not dubious in nature
				if ( $f->has_tag('orf_classification') ) {
					# feature has the orf classification tag, check it
					my @tags = $f->get_tag_values('orf_classification');
					unless ($tags[0] eq 'Dubious') {
						# only keep the non-dubious genes
						push @orflist, $f;
					}
				}
				else {
					# otherwise put them all in the array
					push @orflist, $f; 
				}
			}
			
						
			# Initialize variables
			my $transcript_type;
			my $name;
			my $type;
			my $sense;
			
			# report raw data
			if ($raw) {
				print RAWFILE "cv $cv\tsv $sv\ttransfrag $transcripts{$cv}{$sv}{transfrag}\tsegment $transcripts{$cv}{$sv}{chromo}:$transcripts{$cv}{$sv}{start}..$transcripts{$cv}{$sv}{end}";
				print RAWFILE "\torflist: ", join " ", @orflist;
			}
			
			
			# Process the transcript
			if (scalar @orflist == 0) {
				# there are no ORFs found
				# so then look for other non-orf genes that may overlap
				my @rnalist = $region->features(-types => [
					qw(ncRNA snRNA snoRNA tRNA rRNA)
				]);
				
				if (scalar @rnalist == 0) {
					
					# No ORFs or RNAs, so next check for repetitive elements
					my @repeatlist = $region->features(-types => [
						qw(
							long_terminal_repeat	
							retrotransposon 
							telomere
							repeat_region
							repeat
						)
					]);
					if (@repeatlist) {
						$transcript_type = 'repeat';
						my @namelist;
						my @typelist;
						foreach (@repeatlist) {
							push @namelist, $_->name;
							push @typelist, $_->type;
							if ($_->strand == $strand) {
								# check the strand of the feature; if more 
								# than one, unlikely to be opposite
								$sense = 'sense';
							} 
							else {
								$sense = 'anti-sense';
								$transcript_type = 'anti-sense repeat';
							}
						}
						$name = join ", ", @namelist;
						$type = join ", ", @typelist;
						
						# report raw data
						if ($raw) {print RAWFILE "\trepeatlist: ", join " ", @repeatlist}
					
						
					# no known feature is overlapping, true unidentified intergenic transcript
					} 
					else {
						if ($transcripts{$cv}{$sv}{size} < $minfragment) {
							 # skip all very small non-gene fragments
							delete $transcripts{$cv}{$sv};
							if ($raw) {print RAWFILE "\tnon-gene transfrag too small - skipping\n"}
							next;
						}
						$transcript_type = 'non-gene';
						$name = 'TransFrag' . ($transcripts{$cv}{$sv}{transfrag} =~ /(\d+)/)[0];
									# generate name using the original enriched window number 
						$type = '.';
						$sense = 'none';
					}
				
				} 
				
				elsif (scalar @rnalist == 1) {
					# one ncRNA, this is good and easy
					$transcript_type = 'single-ncRNA';
					my $rna = shift @rnalist;
					$name = $rna->name;
					$type = $rna->type;
					if ($rna->strand == $strand) {
						$sense = 'sense';
					} else {
						$sense = 'anti-sense';
						$transcript_type = 'anti-sense ncRNA';
					}
				
				} 
				
				else {
					# more than one ncRNA
					$transcript_type = 'multi-ncRNA';
					
					# check whether the genes are all on the same strand or not
					my $check;
					my $current = $rnalist[0]->strand;
					for my $i (1..$#rnalist) {
						if ($rnalist[$i]->strand == $current) {
							$check = 'same';
						} else {
							$check = 'different';
							last; # at least two are on different strands, 
								  # that's all we need to know
						}
					}
					
					# collect the names
					my @namelist; # an array for all the names
					my @typelist;
					if ($check eq 'same') { 
						# multiple genes all on same strand
						foreach (@rnalist) {
							push @namelist, $_->name;
							push @typelist, $_->type;
						}
						if ($current == $strand) { 
							# check whether sense or anti-sense
							$sense = 'sense';
						} 
						else {
							$sense = 'anti-sense';
							$transcript_type = 'anti-sense multi-ncRNA';
						}
						$name = join " ", @namelist;
						$type = join " ", @typelist;
					} elsif ($check eq 'different') { 
						# multiple genes not all on same strand
						# in this case only take the sense gene
						foreach (@rnalist) {
							if ($_->strand == $strand) {
								# only take the gene(s) on the sense strand
								push @namelist, $_->name;
								push @typelist, $_->type;
							}
						}
						$sense = 'sense';
						# re-evaluate the multi-strand status
						if (scalar @namelist == 1) { 
							# we now only have 1 gene
							$transcript_type = 'single-ncRNA'; # re-set type to single
							$name = $namelist[0];
							$type = $typelist[0];
						} 
						else { 
							# still have multiple genes
							$name = join " ", @namelist;
							$type = join " ", @typelist;
						}
					}
				}
				
				# report raw data
				if ($raw) {print RAWFILE "\trnalist: ", join " ", @rnalist}
				
				
			} 
			
			elsif (scalar @orflist == 1) { 
				# there was only one orf found - that is good and fine and easy
				$transcript_type = 'single-orf';
				my $orf = shift @orflist;
				$name = $orf->name;
				$type = $orf->type;
				if ($orf->strand == $strand) {
					$sense = 'sense';
				} 
				else {
					$sense = 'anti-sense';
					$transcript_type = 'anti-sense ORF';
				}
				
			} 
			
			else {
				# there are more than one orf - an ambiguous situation at best
				$transcript_type = 'multi-orf';
				$type = '.'; 
				
				# check whether the genes are all on the same strand or not
				my $check;
				my $first_strand = $orflist[0]->strand;
				for my $i (1..$#orflist) {
					if ($orflist[$i]->strand == $first_strand) {
						$check = 'same';
					} 
					else {
						$check = 'different';
						last; # at least two are on different strands, 
							  # that's all we need to know
					}
				}
				
				# collect the names
				my @namelist; # an array for all the names
				if ($check eq 'same') { 
					# multiple genes all on same strand
					foreach (@orflist) {
						push @namelist, $_->name;
					}
					if ($orflist[0]->strand == $strand) { 
						# check whether sense or anti-sense
						$sense = 'sense';
					} 
					else {
						$sense = 'anti-sense';
						$transcript_type = 'anti-sense ORF';
					}
					$name = join "_", @namelist;
				} 
				elsif ($check eq 'different') { 
					# multiple genes not all on same strand
					# in this case only take the sense gene
					foreach (@orflist) {
						if ($_->strand == $strand) {
							# only take the gene(s) on the sense strand
							push @namelist, $_->name;
						}
					}
					$sense = 'sense'; # reset because we're taking sense genes
					
					# re-evaluate the multi-strand status
					if (scalar @namelist == 1) { 
						# we now only have 1 gene
						$transcript_type = 'single-orf'; # reset type to single
						$name = $namelist[0];
						
						# now need to get the actual type
						foreach (@orflist) {
							# we have the gene objects in orflist array
							# but have to find the right one by name
							if ($_->name eq $name) {
								# this is it
								$type = $_->type;
								last;
							}
						}
					} 
					
					else { 
						# still have multiple genes
						$name = join "_", @namelist;
					}
				}
			}
			
			# store transcript associations in the hash
			$transcripts{$cv}{$sv}{transcript_type} = $transcript_type;
			$transcripts{$cv}{$sv}{type} = $type;
			$transcripts{$cv}{$sv}{name} = $name;
			$transcripts{$cv}{$sv}{sense} = $sense;
			
			
			# store the transcript fragment id in specific hashes
			# store these in an array in a hash keyed by the gene name
			my $fragmentID = "$cv\t$sv";
			if ($sense eq 'sense') {
				push @{ $sensetranscripts{$name} }, $fragmentID;
			} 
			elsif ($sense eq 'anti-sense') {
				push @{ $antisensetranscripts{$name} }, $fragmentID;
			}
			
			if ($raw) {
				print RAWFILE "\tfinal: $sense $transcript_type $type $name\n";
			}
		}
	}
}


sub process_transcripts {
	# Process the transcription fragments into full transcripts. Sometimes
	# there are multiple transcription fragments in the list for the same gene.
	# This is most likely due to either multiple exons, lack of consistent 
	# expression across the gene, or very low levels of expression. 
	# We will fix these by merging the fragments into the same transcript.
	# We will only merge transcription fragments if they are associated with the
	# same genomic feature. Otherwise, we have no basis that they should be merged.
	
	# Identify and process the multi gene transcripts
	my $transcriptref = shift;
	foreach my $gene (sort {$a cmp $b} keys %{ $transcriptref }) {
		my @fragments = @{ $transcriptref->{$gene} };
		if (scalar @fragments > 1) {
			# 
			# there are more than one transcription fragments associated with this gene
			# pass this multi-transcipt gene to another subroutine that will merge
			merge_multi_transfrags(@fragments);
		} 
		else {
			# There is only one transcript for this gene; yeah!
			# check the size of anti-sense transcripts
			my ($cv, $sv) = split /\t/, $fragments[0];
			my $size = $transcripts{$cv}{$sv}{size};
			if (
				$transcripts{$cv}{$sv}{sense} eq 'anti-sense' and 
				$size < $minfragment
			) {
				# antisense transcript is less than the allowed minimum fragment size
				# want to delete the very small anti-sense transcripts
				delete $transcripts{$cv}{$sv}; # delete this transcript 
				if ($raw) {
					print RAWFILE "Transcript for gene $gene is antisense and too small (size $size); Transcript deleted\n";
				}
			} 
			else {
				if ($raw) {
					print RAWFILE "Gene $gene was found in one transcript\n";
				}
			}
		}
	}
	
}



# Subroutine to do perform the merging of multiple transcription fragments
sub merge_multi_transfrags {
	# we will assume they are actually all part of one primary transcript
	# (no alternative splicing in yeast) so simply assign the final transcript
	# the smallest start value and the biggest stop value
	
	my @fragments = @_; # an array of the transcript fragments ids to merge
	my @allstarts;
	my @allstops;
	my ($pcv, $psv); # the primary chromo, start values for a representitive transcript
					 # that will become the primary transcript for this gene
	foreach ( @fragments ) {
		# get the transfrag id info
		my ($cv, $sv) = split /\t/;
		# collect the starts and stops for each transcription fragment
		push @allstarts, $transcripts{$cv}{$sv}{start};
		push @allstops, $transcripts{$cv}{$sv}{end};
		# define the primary representitive transcription fragment
		if ($pcv) { 
			# already defined
			$transcripts{$pcv}{$psv}{transfrag} .= ", $transcripts{$cv}{$sv}{transfrag}";
			# after we get the information we need from this transfrag, delete it
			delete $transcripts{$cv}{$sv};
		} 
		else {
			# define it
			$pcv = $cv;
			$psv = $sv;
		}
		
	}
	
	# determine new positions for the merged transcript
	my $start = min(@allstarts); # the transcript start will be the smallest transfrag start
	my $stop = max(@allstops); # the transcript stop will be the biggest transfrag stop
	my $size = $stop - $start + 1;
	
	# check the size of anti-sense transcripts
	if (
		$transcripts{$pcv}{$psv}{sense} eq 'anti-sense' and 
		$size < $minfragment
	) {
		# antisense transcript is less than the allowed minimum fragment size
		# want to delete the very small anti-sense transcripts
		delete $transcripts{$pcv}{$psv}; # delete this transcript 
		if ($raw) {
			print RAWFILE "Transcript for gene $transcripts{$pcv}{$psv}{name} is antisense and too small (size $size); Transcript deleted\n";
		}
		return;
	}
	
	# generate score for the reason
	my $score = sprintf "%.3f", get_chromo_region_score(
				'db'      => $db,
				'chromo'  => $transcripts{$pcv}{$psv}{chromo},
				'start'   => $start,
				'end'     => $stop,
				'method'  => 'median',
				'stranded' => 'all', # force it to look at only this dataset (inherently stranded)
				'dataset' => # assign the appropriate dataset based on strand
					$transcripts{$pcv}{$psv}{strand} == 1 ? $fdata : $rdata,
				'log'     => $log,
	); # calculate a new score for the merged transcript
	
	
	# update the information on the merged transcript
	if (defined $score and $score >= $threshold) { 
		# only if the merged fragment passes the threshold
		$transcripts{$pcv}{$psv}{start} = $start;
		$transcripts{$pcv}{$psv}{end}   = $stop;
		$transcripts{$pcv}{$psv}{size}  = $size;
		$transcripts{$pcv}{$psv}{score} = $score;
		if ($raw) {
			print RAWFILE "Transcripts for gene $transcripts{$pcv}{$psv}{name} was merged: now $start\..$stop\n";
		}
	} 
	else {
		# this merged transcript does not have the signal to pass the threshold
		if ($raw) {
			print RAWFILE "Transcripts for gene $transcripts{$pcv}{$psv}{name} was merged (now $start\..$stop) but signal ($score) was below threshold; Transcripts deleted\n";
		}
		delete $transcripts{$pcv}{$psv}; # delete this transcript 
	}
	
}



# Subroutine to calculate the UTRs and the extent of merged fragments
sub calculate_utr_and_extent {
	for my $cv (sort {$a <=> $b} keys %transcripts) {
		# for each chromosome value
		for my $sv (sort {$a <=> $b} keys %{ $transcripts{$cv} } ) {
			# for each start value
			
			# single sense strand genes
			if ( 
					(
					$transcripts{$cv}{$sv}{transcript_type} eq 'single-orf' or
					$transcripts{$cv}{$sv}{transcript_type} eq 'single-ncRNA'
					) 
				and
				$transcripts{$cv}{$sv}{sense} eq 'sense'
			) {
				# only check sense single genes, other features don't need checking
				
				my @genes = $db->features( 
						# set the gene object
						-name => $transcripts{$cv}{$sv}{name},
						-type => $transcripts{$cv}{$sv}{type},
				);
				unless (@genes) {
					warn " found no gene features for $transcripts{$cv}{$sv}{type} $transcripts{$cv}{$sv}{name}!\n";
					# put in null values for their UTRs
					$transcripts{$cv}{$sv}{extent} = ".";
					$transcripts{$cv}{$sv}{utr5} = "."; # no need for UTR
					$transcripts{$cv}{$sv}{utr3} = ".";
					$transcripts{$cv}{$sv}{utr5length} = '.';
					$transcripts{$cv}{$sv}{utr3length} ='.';
					next;
				}
					
				my $gene_seg = $genes[0]->segment; # convert to segment
				my $transcript = $db->segment( # set the transcript object 
						# extended by $tolerance value
						-name  => $transcripts{$cv}{$sv}{chromo},
						-start => $transcripts{$cv}{$sv}{start} - $tolerance,
						-stop  => $transcripts{$cv}{$sv}{end} + $tolerance
				);
				unless ($transcript) { 
					die " can't establish transcript object at chr$cv:$sv (hash value: $transcripts{$cv}{$sv}{name})!\n";
				}
				
				if ($transcript->contains($gene_seg) ) { # function to check for complete overlap
					# gene is completely within the  transcript
					$transcripts{$cv}{$sv}{extent} = 'complete'; # set this regardless of previous value
					
					# determine the 5' and 3' UTRs
					if ($transcripts{$cv}{$sv}{transcript_type} eq 'single-orf') {
						
						# forward strand
						if ($transcripts{$cv}{$sv}{strand} == 1) { 
							# remember, we're working with absolute coordinates 
							# here, but trying to get relative coordinates
							my $utr5start = 1;
							my $utr5stop = $gene_seg->start - 
								$transcripts{$cv}{$sv}{start} - 1;
							my $utr3start = $gene_seg->end - 
								$transcripts{$cv}{$sv}{start} + 1;
							my $utr3stop = $transcripts{$cv}{$sv}{end} - 
								$transcripts{$cv}{$sv}{start};
							$transcripts{$cv}{$sv}{utr5} = 
								"$utr5start\..$utr5stop";
							$transcripts{$cv}{$sv}{utr3} = 
								"$utr3start\..$utr3stop";
							$transcripts{$cv}{$sv}{utr5length} = 
								$utr5stop - $utr5start + 1;
							$transcripts{$cv}{$sv}{utr3length} = 
								$utr3stop - $utr3start + 1;
						} 
						
						# reverse strand
						elsif ($transcripts{$cv}{$sv}{strand} == -1) { 
							my $utr5start = 1;
							my $utr5stop = $transcripts{$cv}{$sv}{end} - 
								$gene_seg->end + 1;
							my $utr3start = $transcripts{$cv}{$sv}{end} - 
								$gene_seg->start - 1;
							my $utr3stop = $transcripts{$cv}{$sv}{end} - 
								$transcripts{$cv}{$sv}{start};
							$transcripts{$cv}{$sv}{utr5} = 
								"$utr5start\..$utr5stop";
							$transcripts{$cv}{$sv}{utr3} = 
								"$utr3start\..$utr3stop";
							$transcripts{$cv}{$sv}{utr5length} = 
								$utr5stop - $utr5start + 1;
							$transcripts{$cv}{$sv}{utr3length} = 
								$utr3stop - $utr3start + 1;
						}
					} 
					
					# no UTRs for ncRNAs
					else { 
						$transcripts{$cv}{$sv}{extent} = 'complete'; # set this regardless of previous value
						$transcripts{$cv}{$sv}{utr5} = '.';
						$transcripts{$cv}{$sv}{utr3} = '.';
						$transcripts{$cv}{$sv}{utr5length} = '.';
						$transcripts{$cv}{$sv}{utr3length} ='.';
					}
				} 
				
				# gene is not within the extended transcript
				else {
					$transcripts{$cv}{$sv}{extent} = 'overlap';
					$transcripts{$cv}{$sv}{utr5} = "."; # unable to determine UTR
					$transcripts{$cv}{$sv}{utr3} = ".";
					$transcripts{$cv}{$sv}{utr5length} = '.';
					$transcripts{$cv}{$sv}{utr3length} ='.';
				}
			
			} 
			
			# all other transcripts: multi-gene, anti-sense, repeat
			else { 
				# put in null values for their UTRs
				$transcripts{$cv}{$sv}{extent} = ".";
				$transcripts{$cv}{$sv}{utr5} = "."; # no need for UTR
				$transcripts{$cv}{$sv}{utr3} = ".";
				$transcripts{$cv}{$sv}{utr5length} = '.';
				$transcripts{$cv}{$sv}{utr3length} ='.';
			}
		}	
			
	}
}



# Convert the transcripts hash into a final tim data structure
sub convert_transcripts_to_tim_data_structure {
	
	# generate the data hash
	my $data = generate_tim_data_structure(
		'mapped_transcripts',
		qw(
			Chromosome
			Start
			End
			Strand
			Size
			Transcript_Type
			GFF_Type
			Sense
			Extent
			Gene_Type
			Gene_Name
			Expression_Score
			5prime_UTR
			5prime_UTR_Length
			3prime_UTR
			3prime_UTR_Length
			Transcription_Fragment
		)
	) or die " unable to generate tim data structure!\n";
	
	# set metadata
	$data->{'db'} = $database;
	$data->{10}{'tolerance'} = $tolerance;
	$data->{11}{'threshold'} = $threshold;
	$data->{11}{'f_dataset'} = $fdata;
	$data->{11}{'r_dataset'} = $rdata;
	$data->{11}{'log'}       = $log;
	
	# load the data table
	foreach my $cv (sort {$a <=> $b} keys %transcripts) { 
		# chromo value
		foreach my $sv (sort {$a <=> $b} keys %{ $transcripts{$cv} } ) { 
			# start value
			
			# determine GFF type
			my $gff_type;
			if ($transcripts{$cv}{$sv}{transcript_type} =~ /antisense/) {
				$gff_type = 'Transcription_fragment';
			} 
			elsif (
				$transcripts{$cv}{$sv}{transcript_type} eq 'single-orf' 
				or 
				$transcripts{$cv}{$sv}{transcript_type} eq 'single-ncRNA'
			) {
				$gff_type = 'Transcript';
			}
			elsif (
				$transcripts{$cv}{$sv}{transcript_type} eq 'multi-orf' 
				or 
				$transcripts{$cv}{$sv}{transcript_type} eq 'multi-ncRNA'
			) {
				$gff_type = 'Multi_transcript';
			}
			else {
				$gff_type = 'Transcription_fragment';
			}
			
			push @{$data->{'data_table'}}, [ (
				$transcripts{$cv}{$sv}{chromo},
				$transcripts{$cv}{$sv}{start},
				$transcripts{$cv}{$sv}{end},
				$transcripts{$cv}{$sv}{strand},
				$transcripts{$cv}{$sv}{size},
				$transcripts{$cv}{$sv}{transcript_type},
				$gff_type,
				$transcripts{$cv}{$sv}{sense},
				$transcripts{$cv}{$sv}{extent},
				$transcripts{$cv}{$sv}{type},
				$transcripts{$cv}{$sv}{name},
				$transcripts{$cv}{$sv}{score},
				$transcripts{$cv}{$sv}{utr5},
				$transcripts{$cv}{$sv}{utr5length},
				$transcripts{$cv}{$sv}{utr3},
				$transcripts{$cv}{$sv}{utr3length},
				$transcripts{$cv}{$sv}{transfrag},
			) ];
			
			$data->{'last_row'}++;
		}
	}
	
	return $data;
}





# Generate the final report
sub report_processing {
	# initialize count variables
	my $single_rna_within = 0;
	my $single_rna_overlap = 0;
	my $multi_rna = 0;
	my $single_orf_within = 0;
	my $single_orf_overlap = 0;
	my $multi_orf = 0;
	my $repeat = 0;
	my $antisense = 0;
	my $nongene = 0;
	my $total = 0;
	
	# count each type of transcript
	foreach my $row (1 .. $main_data_ref->{last_row}) {
		my $type = $main_data_ref->{data_table}->[$row][5];
		my $extent = $main_data_ref->{data_table}->[$row][8];
		if ($type eq 'non-gene') {
			$nongene++;
		} elsif ($type eq 'repeat') {
			$repeat++;
		} elsif ($type =~ /anti-sense/) {
			$antisense++;
		} elsif ($type eq 'single-ncRNA') {
			if ($extent eq 'complete') {
				$single_rna_within++;
			} elsif ($extent eq 'overlap') {
				$single_rna_overlap++;
			}
		} elsif ($type eq 'multi-ncRNA') {
			$multi_rna++;
		} elsif ($type eq 'single-orf') {
			if ($extent eq 'complete') {
				$single_orf_within++;
			} elsif ($extent eq 'overlap') {
				$single_orf_overlap++;
			}
		} elsif ($type eq 'multi-orf') {
			$multi_orf++;
		}
		$total++;
	}
	
	# generate statistics
	my $totalin = $fcount + $rcount;
	my $one_rna_total = $single_rna_within + $single_rna_overlap;
	my $rna_total = $one_rna_total + $multi_rna;
	my $single_orf_total = $single_orf_within + $single_orf_overlap;
	my $orf_total = $single_orf_total + $multi_orf;
	
	# print the report
	open REPORT, ">$outfile\_report.txt";
	print REPORT "Summary results from transcript mapping\n\n";
	print REPORT " Parameters for transcript mapping:\n";
	print REPORT " # Database $database\n";
	if ($fdata) {print REPORT " # Forward strand dataset $fdata\n"}
	if ($rdata) {print REPORT " # Reverse strand dataset $rdata\n"}
	print REPORT " # Threshold $threshold\n";
	if ($win) {
		print REPORT " # Window $win\n";
		print REPORT " # Step $step\n";
	}
	print REPORT " # Forward strand datafile $ffile\n";
	print REPORT " # Reverse strand datafile $rfile\n";
	print REPORT " # Within tolerance $tolerance\n";
	print REPORT "$totalin transcript fragments loaded\n";
	print REPORT " $fcount from forward strand\n";
	print REPORT " $rcount from reverse strand\n";
	print REPORT "$total transcripts identified after final merging\n";
	print REPORT " $orf_total transcripts mapped to ORFs\n";
	print REPORT "  $single_orf_total transcripts mapped to single ORF\n";
	print REPORT "   $single_orf_within has orf completely within\n";
	print REPORT "   $single_orf_overlap has orf overlap\n";
	print REPORT "  $multi_orf transcripts mapped to multiple ORFs\n";
	print REPORT " $rna_total transcripts mapped to ncRNAs\n";
	print REPORT "  $one_rna_total transcripts mapped to a single ncRNA\n";
	print REPORT "   $single_rna_within has RNA completely within\n";
	print REPORT "   $single_rna_overlap has RNA overlap\n";
	print REPORT "  $multi_rna transcripts mapped to multiple ncRNAs\n";
	print REPORT " $repeat transcripts mapped to repetitive elements\n";
	print REPORT " $antisense transcripts appear to be antisense products\n";
	print REPORT " $nongene transcripts mapped to no known gene features\n";
	close REPORT;
}



__END__

=head1 NAME

map_transcripts.pl

A script to associate enriched regions of transcription with gene annotation.

=head1 SYNOPSIS

map_transcripts.pl --ffile <filename> --rfile <filename> --out <filename>

map_transcripts.pl --db <name> --win <i> --step <i> --fdata <name> 
--rdata <name> --thresh <i> --out <filename>
  
  Options:
  --out <filename>
  --ffile <filename>
  --rfile <filename>
  --db <db_name>
  --win <integer>
  --step <integer>
  --thresh <number>
  --fdata <dataset | filename>
  --rdata <dataset | filename>
  --tol <integer>
  --min <integer>
  --(no)gff
  --source <text>
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --out <filename>

Provide the base output filename for the transcript data file, GFF file, 
and summary file.

=item --ffile <filename>

=item --rfile <filename>

Provide the names of eniched regions files representing transcription 
fragments. The files should be generated by the program 
C<find_enriched_regions.pl>. If these files are not provided, they 
will be automatically generated by running the find_enriched_regions.pl 
program using the parameters defined here. Two files should be provided: 
each representing transcription fragments on the forward and reverse 
strands. 

=item --db <database_name>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. For more 
information about using annotation databases, see 
L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. 

=item --win <integer>

Specify the window size in bp to scan the genome for enriched regions. 
When the find_enriched_regions.pl program is run, trimming is turned 
on so that the final regions are not limited to multiples of 
window size.

=item --step <integer>

Specify the step size in bp for advancing the window across the genome 
when scanning for enriched regions.

=item --thresh <number>

Specify the threshold value to identify regions that are enriched, 
i.e. expressed. A relatively low value should be set to ensure 
low expressed transcripts are also included. The value may be gleaned 
from the metadata in the provided enriched regions files.

=item --fdata <dataset>

=item --rdata <dataset>

Specify the name of the datasets in the database representing the 
forward and reverse transcription data. Alternatively, 
the paths of data files may be provided. Supported formats include 
Bam (.bam), BigBed (.bb), or BigWig (.bw). Files may be local or 
remote (http:// or ftp://).  

=item --tol <integer>

Specify the tolerance distance in bp with which the ends of a 
transcription fragment and gene need to be within to call the 
transcription fragment a complete transcript. If the ends 
exceed the tolerance value, then the transcription fragment is 
labeled as 'overlap'. The default value is 20 bp.

=item --min <integer>

Small, spurious transcription fragments may be identified, particularly 
if the threshold is set too low or the dataset is noisy. This 
parameter sets the minimum size of the transcription fragment to 
initiate a search for associated known genes. The default value is 
100 bp. 

=item --(no)gff

Indicate whether a GFF file should also be written. A version 3 file 
is generated. The default is true.

=item --source <text>

Specify the GFF source value when writing a GFF file. The default 
value is the name of this program.

=item --version

Print the version number.

=item --help

Display the POD documentation.

=back

=head1 DESCRIPTION

This program will identify transcription fragments that correspond to 
gene transcripts using available transcriptome microarray data. Specifically, 
it will identify the start and stop coordinates of transcribed regions 
that correspond or overlap an annotated gene or ORF. It does not identify 
exons or introns. 

The program was initially written to address the lack of officially mapped 
transcripts in the S. cerevisiae genome, which has few and small introns. 
It may work well other similar genomes, but probably not complex 
metazoan genomes with very large introns.

Transcription fragments are identified as windows of enrichment using the 
script 'find_enriched_regions.pl'. This script can either be automatically
executed using the specified parameters, or run separately. If run separately, 
the output files should be indicated (--ffile and --rfile).

Each enriched window, or transcription fragment, is checked for corresponding 
or overlapping genomic features. The name(s) of the overlapping gene(s) are 
reported, as well as a classification for the transcript. Transcripts that 
completely contain (within a default 20 bp tolerance) a single known gene 
(ORF or ncRNA) are labeled as 'complete'. Transcription fragments that 
simply overlap a known gene are labeled as 'overlap'. Transcription fragments 
that overlap more than one annotated gene are labeled as 'multi-orf'.

Transcription fragments that only overlap annotated genes on the opposite 
strand are labeled as 'anti-sense'. Transcription fragments overlapping 
repetitive elements or no known feature are also reported.

Some rudimentary calculations are performed to identify the length of 
the 5' and 3' UTRs. 

The program writes out a tab delimited text file. It will also write out
a gff file for the genome browser. It also writes out a summary report file.

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
