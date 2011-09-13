#!/usr/bin/perl

# convert ucsc refseq table file to gff3


use strict;
#use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqFeature::Lite;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	find_column_index
	format_with_commas
);
use tim_file_helper qw(
	open_tim_data_file
	open_to_read_fh
	open_to_write_fh
);
#use Data::Dumper;

print "\n A script to convert UCSC tables to GFF3 files\n\n";





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
	$genetablef,
	$refseqstatusf,
	$refseqsumf,
	$chromof,
	$source,
	$do_gene,
	$do_utr,
	$do_codon,
	$outfile,
	$gz,
	$help, 
);
GetOptions( 
	'table=s'    => \$genetablef, # the input refseq file
	'status=s'   => \$refseqstatusf, # the refseqstatus file
	'sum=s'      => \$refseqsumf, # the refseqsummary file
	'chromo=s'   => \$chromof, # a chromosome file
	'source=s'   => \$source, # the GFF source
	'gene!'      => \$do_gene, # include genes in output
	'utr!'       => \$do_utr, # include UTRs in output
	'codon!'     => \$do_codon, # include start & stop codons in output
	'out=s'      => \$outfile, # output file name
	'gz!'        => \$gz, # compress file
	'help'       => \$help, # request help
) or die " unrecognized options! use --help\n";

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}





### Check requirements and defaults
unless ($genetablef or $chromof) {
	die " Input table file not specified!\n";
}
unless ($source) {
	$source = 'UCSC';
}
unless (defined $do_gene) {
	$do_gene = 1;
}
unless (defined $do_utr) {
	$do_utr = 1;
}
my $start_time = time;




### Open files

# input files
my $refseqsum = load_refseq_summary_data($refseqsumf);

my $refseqstat = load_refseq_status_data($refseqstatusf);

my ($table_fh, $table_metadata) = open_tim_data_file($genetablef) or
	die " unable to open gene table file!\n";

# ouput file
if ($outfile) {
	# add extension as necessary
	unless ($outfile =~ m/\.gff3?$/) {
		$outfile .= '.gff';
	}
}
else {
	# assign default output file
	$outfile = $table_metadata->{'basename'} . '.gff';
}
my $gff_fh = open_to_write_fh($outfile, $gz) or
	die " unable to open file for writing!\n";

# print comments
$gff_fh->print( "##gff-version 3\n");
$gff_fh->print( "# Generated on " . localtime(time) . "\n");
$gff_fh->print( "# UCSC gene table file $genetablef\n");
$gff_fh->print( "# UCSC RefSeq Status file $refseqstatusf\n") if $refseqstatusf;
$gff_fh->print( "# UCSC RefSeq Summary file $refseqsumf\n") if $refseqsumf;
$gff_fh->print( "# UCSC chromosome file $chromof\n") if $chromof;



### Print chromosomes
if ($chromof) {
	print " Writing chromosome features....\n";
	print_chromosomes();
}



### Process the files

print " Converting gene table features....\n";
my $count = process_gene_table();

# report outcomes
print " Converted ", format_with_commas($count->{gene}), 
	" gene features\n" if $count->{gene} > 0;
print " Converted ", format_with_commas($count->{mrna}), 
	" mRNA transcripts\n" if $count->{mrna} > 0;
print " Converted ", format_with_commas($count->{pseudogene}), 
	" pseudogene transcripts\n" if $count->{pseudogene} > 0;
print " Converted ", format_with_commas($count->{ncrna}), 
	" ncRNA transcripts\n" if $count->{ncrna} > 0;
print " Converted ", format_with_commas($count->{mirna}), 
	" miRNA transcripts\n" if $count->{mirna} > 0;
print " Converted ", format_with_commas($count->{snrna}), 
	" snRNA transcripts\n" if $count->{snrna} > 0;
print " Converted ", format_with_commas($count->{snorna}), 
	" snoRNA transcripts\n" if $count->{snorna} > 0;
print " Converted ", format_with_commas($count->{trna}), 
	" tRNA transcripts\n" if $count->{trna} > 0;
print " Converted ", format_with_commas($count->{rrna}), 
	" rRNA transcripts\n" if $count->{rrna} > 0;
print " Converted ", format_with_commas($count->{other}), 
	" other transcripts\n" if $count->{other} > 0;





### Finish
print " Completed. wrote file '$outfile'\n";
printf " Finished in %.1f minutes\n", (time - $start_time)/60;

exit;





#########################  Subroutines  #######################################

sub load_refseq_summary_data {
	
	# initialize 
	my $file = shift;
	my %sumdata;
	
	# load file
	if ($file) {
		my ($fh, $metadata) = open_tim_data_file( $file ) or 
			die " unable to open summary file '$file'!\n";
		
		# check headers
		if (
			$metadata->{0}{'name'} eq '#mrnaAcc' and 
			$metadata->{1}{'name'} eq 'completeness' and 
			$metadata->{2}{'name'} eq 'summary'
		) {
			# file looks good, proceed
		}
		else {
			# headers don't match
			
			if ($metadata->{'number_columns'} == 3) {
				# file has the correct number of columns
				# but not the right headers
				# the file may have been downloaded directly from the FTP site
				# and not through the table browser
				# the table browser adds the column headers as a comment line
				# whereas the FTP files do not include it 
				
				# we'll blindly assume the file is ok and proceed
				
				# need to add back the the column headers
				unshift @{ $metadata->{'data_table'} }, 
					[ qw(#mrnaAcc completeness summary) ];
			}
			else {
				die " Either the summary file format has changed or you've loaded the\n".
					" wrong file. Please confirm this and take appropriate action\n";
			}
		}
		
		# load into hash
		while (my $line = $fh->getline) {
			chomp $line;
			my @a = split /\t/, $line;
			if (exists $sumdata{ $a[0] } ) {
				warn " summary line for refseq RNA $a[0] exists twice!\n";
			}
			else {
				$sumdata{ $a[0] } = [ ( $a[1], $a[2] ) ];
			}
		}
		
		# finish
		print " loaded ", format_with_commas( scalar(keys %sumdata) ), 
			" gene summaries from '$file'\n";
		$fh->close;
	}
	
	return \%sumdata;
}



sub load_refseq_status_data {
	
	# initialize 
	my $file = shift;
	my %statdata;
	
	# load file
	if ($file) {
		my ($fh, $metadata) = open_tim_data_file( $file ) or 
			die " unable to open status file '$file'!\n";
		
		# check headers
		if (
			$metadata->{0}{'name'} eq '#mrnaAcc' and 
			$metadata->{1}{'name'} eq 'status' and 
			$metadata->{2}{'name'} eq 'mol'
		) {
			# file looks good, proceed
		}
		else {
			# headers don't match
			
			if ($metadata->{'number_columns'} == 3) {
				# file has the correct number of columns
				# but not the right headers
				# the file may have been downloaded directly from the FTP site
				# and not through the table browser
				# the table browser adds the column headers as a comment line
				# whereas the FTP files do not include it 
				
				# we'll blindly assume the file is ok and proceed
				
				# need to add back the the column headers
				unshift @{ $metadata->{'data_table'} }, 
					[ qw(#mrnaAcc status mol) ];
			}
			else {
				die " Either the status file format has changed or you've loaded the\n".
					" wrong file. Please confirm this and take appropriate action\n";
			}
		}
		
		# load into hash
		while (my $line = $fh->getline) {
			chomp $line;
			my @a = split /\t/, $line;
			if (exists $statdata{ $a[0] } ) {
				warn " status line for refseq RNA $a[0] exists twice!\n";
			}
			else {
				$statdata{ $a[0] } = [ ( $a[1], $a[2] ) ];
			}
		}
		
		# finish
		print " loaded ", format_with_commas( scalar(keys %statdata) ), 
			" RNA status from '$file'\n";
		$fh->close;
	}
	
	return \%statdata;
}




sub identify_indices {
	
	my %column2index;
	
	foreach (
		# column name, default index number
		['name', 1],
		['chrom', 2], 
		['strand', 3],
		['txStart', 4],
		['txEnd', 5],
		['cdsStart', 6],
		['cdsEnd', 7],
		['exonCount', 8],
		['exonStarts', 9],
		['exonEnds', 10],
		['name2', 12],
		['exonFrames', 15],
	) {
		# we need to identify which columns are which,
		# but this is a huge headache when the column headers are not present
		# as when the file was downloaded from UCSC FTP site rather than going
		# through the table browser (which adds the column headers as a comment 
		# line)
		# therefore, we need to make assumptions, and hope that the assumption
		# is valid, otherwise bad things will happen
		
		my ($name, $default) = @{$_};
		my $index;
		
		# check the default location
		if ($table_metadata->{$default}{'name'} eq $name) {
			# it's a match, presumably because it was a standard gene table
			# downloaded from the table browser
			$index = $default;
		}
		
		# else attempt to identify it
		else {
			# maybe a non-standard table? or it has not column headers
			$index = find_column_index($table_metadata, $name);
			
			unless (defined $index) {
				warn " unable to identify column index for '$name'," .
					" assuming default index of $default\n";
				if ($default < $table_metadata->{'number_columns'}) {
					# make sure we have at least this many columns
					$index = $default;
				}
				else {
					die " cannot assign '$name' to an index!\n";
				}
			}
		}
		
		# assign 
		$column2index{$name} = $index;
	}
	
	# finished
	return \%column2index;
	
	# column positions for a standard RefSeq table
		# 0  bin
		# 1  name
		# 2  chrom
		# 3  strand
		# 4  txStart
		# 5  txEnd
		# 6  cdsStart
		# 7  cdsEnd
		# 8  exonCount
		# 9  exonStarts
		# 10 exonEnds
		# 11 score
		# 12 name2
		# 13 cdsStartStat
		# 14 cdsEndStat
		# 15 exonFrames
}




sub process_line_data {
	
	my ($c2i, $line) = @_;
	my %data;
	
	# load the relevant data from the table line into the hash
	# using the identified column indices
	chomp $line;
	my @linedata = split /\t/, $line;
	foreach my $key (keys %{$c2i}) {
		$data{$key} = $linedata[ $c2i->{$key} ];
	}
	
	return \%data;
}




sub process_gene_table {
	
	# initialize 
	my $current_chrom;
	my %gene2seqf; # hash to assemble genes and/or transcripts for this chromosome
	my %id2count; # hash to aid in generating unique primary IDs
	my %counts = (
		'gene'       => 0,
		'mrna'       => 0,
		'pseudogene' => 0,
		'ncrna'      => 0,
		'mirna'      => 0,
		'snrna'      => 0,
		'snorna'     => 0,
		'trna'       => 0,
		'rrna'       => 0,
		'other'      => 0,
	);
	
	# identify relevant column indices
	my $c2i = identify_indices();
	
	
	
	#### Main Loop
	while (my $line = $table_fh->getline) {
		
		## process the row from the gene table
		my $linedata = process_line_data($c2i, $line);
		
		
		## check the chromosome
		unless (defined $current_chrom) {
			$current_chrom = $linedata->{chrom};
		}
		if ($linedata->{chrom} ne $current_chrom) {
			# moved on to next chromosome
			# the table should be sorted by chromosome
			
			# print the current gene list and prepare for next chromosome
			print "  writing ", format_with_commas( scalar keys %gene2seqf ),
				$do_gene ? " genes" : " transcripts",
				" for chromosome $current_chrom\n";
			print_current_gene_list(\%gene2seqf);
			
			$current_chrom = $linedata->{chrom};
			%gene2seqf = ();
		}
		
		
		## generate the transcript
		my $transcript = generate_new_transcript($linedata, \%id2count);
		
		
		## count the transcript type
		my $type = $transcript->primary_tag;
		if ($type eq 'mRNA') {
			$counts{mrna}++;
		}
		elsif ($type eq 'pseudogene') {
			$counts{pseudogene}++;
		}
		elsif ($type eq 'ncRNA') {
			$counts{ncrna}++;
		}
		elsif ($type eq 'miRNA') {
			$counts{mirna}++;
		}
		elsif ($type eq 'snRNA') {
			$counts{snrna}++;
		}
		elsif ($type eq 'snoRNA') {
			$counts{snorna}++;
		}
		elsif ($type eq 'tRNA') {
			$counts{trna}++;
		}
		elsif ($type eq 'rRNA') {
			$counts{rrna}++;
		}
		else {
			$counts{other}++;
		}
		
		
		## add transcript to gene if requested
		if ($do_gene) {
			# assemble transcripts into genes
			# multiple transcripts may be associated with a single gene
			# genes are store in the gene2seqf hash
			# there may be more than one gene with each gene name, each with 
			# a non-overlapping transcript (how complicated!)
			
			my $gene;
			if (exists $gene2seqf{ lc $linedata->{name2} }) {
				# we already have a gene for this transcript
				
				# pull out the gene seqfeature(s)
				my $genes = $gene2seqf{ lc $linedata->{name2} };
				
				# check that the current transcript intersects with the gene
				# sometimes we can have two separate transcripts with the 
				# same gene name, but located on opposite ends of the chromosome
				# part of a gene family, but unlikely the same gene 200 Mb in 
				# length
				foreach my $g (@$genes) {
					if ( $g->overlaps($transcript, 'strong') ) {
						# gene and transcript overlap on the same strand
						# we found the intersecting gene
						$gene = $g;
						last;
					}
				}
				
				# we have a gene for our transcript
				if ($gene) {
					# update the gene coordinates if necessary
					if ( ($linedata->{txStart} + 1) < $gene->start) {
						# update the transcription start position
						$gene->start( $linedata->{txStart} + 1 );
					}
					if ($linedata->{txEnd} > $gene->end) {
						# update the transcription stop position
						$gene->end( $linedata->{txEnd} );
					}
				}
				
				# NONE of the genes and our transcript overlap
				else {
					# must make a new gene
					$gene = generate_new_gene($linedata, \%id2count);
					$counts{gene}++;
					
					# store the new gene oject into the gene hash
					push @{ $genes }, $gene;
				}
			}
			
			else {
				# generate new gene SeqFeature object
				$gene = generate_new_gene($linedata, \%id2count);
				$counts{gene}++;
				
				# store the gene oject into the gene hash
				$gene2seqf{ lc $linedata->{name2} } = [ $gene ];
			} 
			
			
			# associate our transcript with the gene
			$gene->add_SeqFeature($transcript);
			
			# Update gene note if necessary
			unless ($gene->has_tag('Note')) {
				# unless it somehow wasn't available for a previous transcript, 
				# but is now, we'll add it now
				# we won't check if transcripts for the same gene have the 
				# same note or not, why wouldn't they????
				if (
					exists $refseqsum->{ $linedata->{name} }
					and
					defined $refseqsum->{ $linedata->{name} }->[1]
				) {
					# check for a summary for this mRNA
					
					# make note if it exists
					my $text = $refseqsum->{ $linedata->{name} }->[1];
					$gene->add_tag_value('Note', $text);;
				}
			} 
			
		}
		else {
			# do not assemble transcripts into genes
			# we will still use the gene2seqf hash, just organized by 
			# transcript name
			# there may be more than one transcript with the same name
			
			if (exists $gene2seqf{ lc $linedata->{name} }) {
				push @{ $gene2seqf{ lc $linedata->{name} } }, $transcript;
			}
			else {
				$gene2seqf{ lc $linedata->{name} } = [ $transcript ];
			}
		}
		
		
	} # Finished working through the table
	
	
	
	#### Finished
	
	# print remaining current genes and transcripts
	print "  writing ", format_with_commas( scalar keys %gene2seqf ),
		$do_gene ? " genes" : " transcripts",
		" for chromosome $current_chrom\n";
	print_current_gene_list(\%gene2seqf);
	
	# return the counts
	return \%counts;
}



sub print_current_gene_list {
	my $gene2seqf = shift;
	
	# first store all seqfeatures in a new hash by their start position 
	my %start2seqf;
	foreach my $g (keys %{ $gene2seqf }) {
		# each gene id points to an array of gene/transcript seqfeatures
		
		foreach my $t (@{ $gene2seqf->{$g} }) {
			my $start = $t->start;
			
			while (exists $start2seqf{$start}) {
				# make sure unique start positions
				$start++;
			}
			
			# store the seqfeature
			$start2seqf{$start} = $t;
		}
	}
	
	# then print all the seqfeatures in ascending order by start position
	foreach my $s (sort {$a <=> $b} keys %start2seqf) {
		
		# set gff version
		$start2seqf{$s}->version(3); 
		
		# print the seqfeature recursively
		$gff_fh->print( $start2seqf{$s}->gff_string(1) . "\n");
			# the gff_string method is undocumented in the POD, but is a 
			# valid method. Passing 1 should force a recursive action to 
			# print both parent and children.
	}
	
	# print directive to close out all previous genes
	$gff_fh->print("###\n"); 
}



sub generate_new_gene {
	my ($linedata, $id2counts) = @_;
	
	# make sure we have a gene name
	# some genes, notably some ncRNA genes, have no gene or name2 entry
	unless ($linedata->{name2}) {
		# we'll fake it and assign the transcript name
		# change it in linedata hash to propagate it in downstream code
		$linedata->{name2} = $linedata->{name};
	}
	
	# Uniqueify the gene ID
	my $id;
	if (exists $id2counts->{ lc $linedata->{name2} }) {
		# we've encountered this transcript ID before
		
		# now need to make ID unique by appending a number
		$id = $linedata->{name2} . '.' . $id2counts->{ lc $linedata->{name2} };
		
		# remember this one
		$id2counts->{ lc $linedata->{name2} } += 1;
	}
	else {
		# this is the first transcript with this id
		$id = $linedata->{name2};
		$id2counts->{lc $id} = 1;
	}
	
	
	# generate the gene SeqFeature object
	my $gene = Bio::SeqFeature::Lite->new(
		-seq_id        => $linedata->{chrom},
		-source        => $source,
		-primary_tag   => 'gene',
		-start         => $linedata->{txStart} + 1,
		-end           => $linedata->{txEnd},
		-strand        => $linedata->{strand} eq '+' ? 1 : -1,
		-phase         => '.',
		-display_name  => $linedata->{name2},
		-primary_id    => $id,
	);
	
	# add status if possible
	if (exists $refseqstat->{  $linedata->{name} } ) {
		$gene->add_tag_value('status', 
			$refseqstat->{ $linedata->{name} }->[0] );
	}
	
	# add Note if possible
	if (
		exists $refseqsum->{ $linedata->{name} }
		and
		defined $refseqsum->{ $linedata->{name} }->[1]
	) {
		# add the the note
		my $text = $refseqsum->{ $linedata->{name} }->[1];
		$gene->add_tag_value('Note', $text);
	}
	
	# finished
	return $gene;
}



sub generate_new_transcript {
	my ($linedata, $id2counts) = @_;
	
	# Uniqueify the transcript ID
	my $id;
	if (exists $id2counts->{ lc $linedata->{name} } ) {
		# we've encountered this transcript ID before
		
		# now need to make ID unique by appending a number
		$id = $linedata->{name} . '.' . $id2counts->{ lc $linedata->{name} };
		
		# remember this one
		$id2counts->{ lc $linedata->{name} } += 1;
	}
	else {
		# this is the first transcript with this id
		$id = $linedata->{name};
		$id2counts->{lc $id} = 1;
	}
	
	# Generate the transcript SeqFeature object
	my $transcript = Bio::SeqFeature::Lite->new(
		-seq_id        => $linedata->{chrom},
		-source        => $source,
		-start         => $linedata->{txStart} + 1,
		-end           => $linedata->{txEnd},
		-strand        => $linedata->{strand} eq '+' ? 1 : -1,
		-phase         => '.',
		-display_name  => $linedata->{name},
		-primary_id    => $id,
	);
	
	# Attempt to identify the transcript type
	if (
		$linedata->{cdsStart} == $linedata->{txEnd} and 
		$linedata->{cdsEnd} == $linedata->{txEnd} 
	) {
		# there appears to be no coding potential when 
		# txEnd = cdsStart = cdsEnd
		# if you'll look, all of the exon phases should also be -1
		
		# if we are using a RefSeq gene table, we may be able to infer
		# some certain types from the gene name
		
		# in the mean time, we'll try to deduce from the gene name
		if ($linedata->{name2} =~ /^LOC\d+/) {
			# empirical testing seems to suggest that all the noncoding 
			# genes with a name like LOC123456 are pseudogenes
			# well, at least with hg18, it may not be true for others
			$transcript->primary_tag('pseudogene');
		}
		if ($linedata->{name2} =~ /^mir/i) {
			# a noncoding gene whose name begins with mir is likely a 
			# a micro RNA
			$transcript->primary_tag('miRNA');
		}
		elsif ($linedata->{name2} =~ /^snr/i) {
			# a noncoding gene whose name begins with snr is likely a 
			# a snRNA
			$transcript->primary_tag('snRNA');
		}
		elsif ($linedata->{name2} =~ /^sno/i) {
			# a noncoding gene whose name begins with sno is likely a 
			# a snoRNA
			$transcript->primary_tag('snoRNA');
		}
		else {
			# a generic ncRNA
			$transcript->primary_tag('ncRNA');
		}
	}
	else {
		# the transcript has an identifiable CDS
		$transcript->primary_tag('mRNA');
	}
	
	
	# add gene name as an alias
	if (defined $linedata->{name2}) {
		$transcript->add_tag_value('Alias', $linedata->{name2});
	}
	
	# add a status for the transcript
	if (exists $refseqstat->{ $linedata->{name} } ) {
		$transcript->add_tag_value('status', 
			$refseqstat->{ $linedata->{name} }->[0] );
	}
	
	# add the completeness value for the tag
	if (exists $refseqsum->{ $linedata->{name} } ) {
		$transcript->add_tag_value('completeness', 
			$refseqsum->{ $linedata->{name} }->[0] );
	}
	
	# Add the exons
	if ($linedata->{exonCount} > 1) { 
		
		# determine how to process the exons by the feature type 
		# and the strand orientation
		
		if (
			$transcript->primary_tag eq 'mRNA' and 
			$transcript->strand > 0
		) {
			# forward strand coding exon
			process_forward_exons($transcript, $linedata);
		}
		elsif (
			$transcript->primary_tag eq 'mRNA' and 
			$transcript->strand < 0
		) {
			# reverse strand coding exon
			process_reverse_exons($transcript, $linedata);
		}
		
		# otherwise assume exons are part of noncoding transcript
		else {
			process_noncoding_exons($transcript, $linedata);
		}
		
	}
	
	# transcript is complete
	return $transcript;
}



sub process_forward_exons {
	my ($transcript, $linedata) = @_;
	
	# strip trailing comma from exon lists
	$linedata->{exonStarts} =~ s/,$//; 
	$linedata->{exonEnds} =~ s/,$//; 
	$linedata->{exonFrames} =~ s/,$//;
	
	# get the exon coordinates
	my @starts = split /,/, $linedata->{exonStarts};
	my @ends   = split /,/, $linedata->{exonEnds};
	my @phases = split /,/, $linedata->{exonFrames};
	unless (scalar @starts == scalar @ends) {
		die " source table error! not equal number of exon starts" .
			" for transcript $linedata->{name}\n";
	}
	my $cdsStart = $linedata->{cdsStart} + 1;
	my $cdsStop = $linedata->{cdsEnd};
	
	# adjust for 0-based
	@starts = map {$_ + 1} @starts;
	
	# add the exons and start/stop codons
	add_exons_codons_subfeatures($transcript, \@starts, \@ends, $cdsStart, 
		$cdsStop);
	
	
	# process the exons for cds and utr
	for (my $i = 0; $i <= $#starts; $i++) {
		
		my $exon; # this will be the exon SeqFeature object
		
		# we need to determine whether the exon is UTR, CDS, or split both
		
		#### 5'UTR only ####
		if (
			$starts[$i] < $cdsStart # cdsStart
			and
			$ends[$i] < $cdsStart
		) {
			# the exon start/end is entirely before the cdsStart
			# we have a 5'UTR
			
			# build the utr object
			if ($do_utr) {
				$exon = Bio::SeqFeature::Lite->new(
					-seq_id        => $linedata->{chrom},
					-source        => $source,
					-start         => $starts[$i],
					-end           => $ends[$i],
					-strand        => 1,
					-phase         => '.',
					-primary_tag   => 'five_prime_UTR',
					-primary_id    => $transcript->primary_id . ".utr$i",
					# -display_name  => $transcript->display_name . ".utr$i",
				);
			};
		}
		
		
		#### Split 5'UTR and CDS ####
		elsif (
			$starts[$i] < $cdsStart 
			and
			$ends[$i] >= $cdsStart
		) {
			# the start codon is in this exon
			# we need to make two features, 
			# one for the utr half, other for the cds
			
			# build the utr half of the object
			if ($do_utr) {
				my $utr_exon = Bio::SeqFeature::Lite->new(
					-seq_id        => $linedata->{chrom},
					-source        => $source,
					-start         => $starts[$i],
					-end           => $cdsStart - 1,
					-strand        => 1,
					-phase         => '.',
					-primary_tag   => 'five_prime_UTR',
					-primary_id    => $transcript->primary_id . ".utr$i",
					# -display_name  => $transcript->display_name . ".utr$i",
				);
				# since we're actually building two objects here, we have to 
				# complete the process of building the utr object and associate 
				# it with the transcript object
				# the cds half of the exon will be finished below
				
				# associate add the utr half to the parent transcript
				$transcript->add_SeqFeature($utr_exon);
			}
			
			# now build the cds half of the object
			$exon = Bio::SeqFeature::Lite->new(
				-seq_id        => $linedata->{chrom},
				-source        => $source,
				-start         => $cdsStart,
				-end           => $ends[$i],
				-strand        => 1,
				-phase         => $phases[$i],
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds$i",
				# -display_name  => $transcript->display_name . ".cds$i",
			);
		}
		
		#### CDS only ####
		elsif (
			$starts[$i] >= $cdsStart 
			and
			$ends[$i] <= $cdsStop
		) {
			# we are in the CDS
			# this will also work with genes that have no UTRs, where the 
			# the cdsStart == txStart, and same with End
			
			$exon = Bio::SeqFeature::Lite->new(
				-seq_id        => $linedata->{chrom},
				-source        => $source,
				-start         => $starts[$i],
				-end           => $ends[$i],
				-strand        => 1,
				-phase         => $phases[$i],
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds$i",
				# -display_name  => $transcript->display_name . ".cds$i",
			);
			
		}
		
		
		#### Split CDS and 3'UTR ####
		elsif (
			$starts[$i] <= $cdsStop 
			and
			$ends[$i] > $cdsStop
		) {
			# the stop codon is in this exon
			# we need to make two features, 
			# one for the cds half, other for the utr
			
			# build the cds half of the object
			$exon = Bio::SeqFeature::Lite->new(
				-seq_id        => $linedata->{chrom},
				-source        => $source,
				-start         => $starts[$i],
				-end           => $cdsStop,
				-strand        => 1,
				-phase         => $phases[$i],
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds$i",
				# -display_name  => $transcript->display_name . ".cds$i",
			);
			
			# now build the utr half of the object
			if ($do_utr) {
				my $utr_exon = Bio::SeqFeature::Lite->new(
					-seq_id        => $linedata->{chrom},
					-source        => $source,
					-start         => $cdsStop + 1,
					-end           => $ends[$i],
					-strand        => 1,
					-phase         => '.',
					-primary_tag   => 'three_prime_UTR',
					-primary_id    => $transcript->primary_id . ".utr$i",
					# -display_name  => $transcript->display_name . ".utr$i",
				);
				
				# associate add the utr half to the parent transcript
				$transcript->add_SeqFeature($utr_exon);
			}
		}
		
		#### 3'UTR only ####
		elsif (
			$starts[$i] > $cdsStop 
			and
			$ends[$i] > $cdsStop
		) {
			# the exon start/end is entirely after the cdsStop
			# we have a 3'UTR
			
			# build the utr object
			if ($do_utr) {
				$exon = Bio::SeqFeature::Lite->new(
					-seq_id        => $linedata->{chrom},
					-source        => $source,
					-start         => $starts[$i],
					-end           => $ends[$i],
					-strand        => 1,
					-phase         => '.',
					-primary_tag   => 'three_prime_UTR',
					-primary_id    => $transcript->primary_id . ".utr$i",
					# -display_name  => $transcript->display_name . ".utr$i",
				);
			}
		}
		
		#### Something's wrong ####
		else {
			# just in case I goofed something up
			die " programming error! the exon coordinates don't match up with" .
				" CDS coordinates for transcript $linedata->{name}\n " . 
				" Forward strand exon coordinates $starts[$i] .. $ends[$i]\n" .
				" CDS coordinates $cdsStart .. $cdsStop\n";
		}
		
		
		# check that we generated an exon
		next unless $exon;
		
		# associate exon with the parent transcript
		$transcript->add_SeqFeature($exon);
		
	}
	
	# done with the exon list
}




sub process_reverse_exons {
	my ($transcript, $linedata) = @_;
		
	# strip trailing comma from exon lists
	$linedata->{exonStarts} =~ s/,$//; 
	$linedata->{exonEnds}   =~ s/,$//; 
	$linedata->{exonFrames} =~ s/,$//;
	
	# get the exon coordinates
	my @starts = split /,/, $linedata->{exonStarts};
	my @ends   = split /,/, $linedata->{exonEnds};
	my @phases = split /,/, $linedata->{exonFrames};
	unless (scalar @starts == scalar @ends) {
		die " source table error! not equal number of exon starts" .
			" and stops for transcript $linedata->{name}\n";
	}
	my $cdsStart = $linedata->{cdsStart} + 1;
	my $cdsStop = $linedata->{cdsEnd};
	
	# adjust for 0-based
	@starts = map {$_ + 1} @starts;
	
	# add the exons and start/stop codons
	add_exons_codons_subfeatures($transcript, \@starts, \@ends, $cdsStart, 
		$cdsStop);
	
	# process the exons for cds and utr
	for (my $i = 0; $i <= $#starts; $i++) {
		
		my $exon; # this will be the exon SeqFeature object
		
		# we need to determine whether the exon is UTR, CDS, or split both
		
		#### 3'UTR only ####
		if (
			$starts[$i] < $cdsStart # cdsStart
			and
			$ends[$i] < $cdsStart
		) {
			# the exon start/end is entirely before the cdsStart
			# we have a 3'UTR
			
			# build the utr object
			if ($do_utr) {
				$exon = Bio::SeqFeature::Lite->new(
					-seq_id        => $linedata->{chrom},
					-source        => $source,
					-start         => $starts[$i],
					-end           => $ends[$i],
					-strand        => -1,
					-phase         => '.',
					-primary_tag   => 'three_prime_UTR',
					-primary_id    => $transcript->primary_id . ".utr$i",
					# -display_name  => $transcript->display_name . ".utr$i",
				);
			}
		}
		
		
		#### Split 3'UTR and CDS ####
		elsif (
			$starts[$i] < $cdsStart 
			and
			$ends[$i] >= $cdsStart
		) {
			# the (stop) codon is in this exon
			# we need to make two features, 
			# one for the utr half, other for the cds
			
			# build the utr half of the object
			if ($do_utr) {
				my $utr_exon = Bio::SeqFeature::Lite->new(
					-seq_id        => $linedata->{chrom},
					-source        => $source,
					-start         => $starts[$i],
					-end           => $cdsStart - 1,
					-strand        => -1,
					-phase         => '.',
					-primary_tag   => 'three_prime_UTR',
					-primary_id    => $transcript->primary_id . ".utr$i",
					# -display_name  => $transcript->display_name . ".utr$i",
				);
				# since we're actually building two objects here, we have to 
				# complete the process of building the utr object and associate 
				# it with the transcript object
				# the cds half of the exon will be finished below
				
				# associate add the utr half to the parent transcript
				$transcript->add_SeqFeature($utr_exon);
			}
			
			# now build the cds half of the object
			$exon = Bio::SeqFeature::Lite->new(
				-seq_id        => $linedata->{chrom},
				-source        => $source,
				-start         => $cdsStart,
				-end           => $ends[$i],
				-strand        => -1,
				-phase         => $phases[$i],
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds$i",
				# -display_name  => $transcript->display_name . ".cds$i",
			);
		}
		
		#### CDS only ####
		elsif (
			$starts[$i] >= $cdsStart 
			and
			$ends[$i] <= $cdsStop
		) {
			# we are in the CDS
			# this will also work with genes that have no UTRs, where the 
			# the cdsStart == txStart, and same with End
			
			$exon = Bio::SeqFeature::Lite->new(
				-seq_id        => $linedata->{chrom},
				-source        => $source,
				-start         => $starts[$i],
				-end           => $ends[$i],
				-strand        => -1,
				-phase         => $phases[$i],
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds$i",
				# -display_name  => $transcript->display_name . ".cds$i",
			);
			
		}
		
		
		#### Split CDS and 5'UTR ####
		elsif (
			$starts[$i] <= $cdsStop 
			and
			$ends[$i] > $cdsStop
		) {
			# the (start) codon is in this exon
			# we need to make two features, 
			# one for the cds half, other for the utr
			
			# build the cds half of the object
			$exon = Bio::SeqFeature::Lite->new(
				-seq_id        => $linedata->{chrom},
				-source        => $source,
				-start         => $starts[$i],
				-end           => $cdsStop,
				-strand        => -1,
				-phase         => $phases[$i],
				-primary_tag   => 'CDS',
				-primary_id    => $transcript->primary_id . ".cds$i",
				# -display_name  => $transcript->display_name . ".cds$i",
			);
			
			# now build the utr half of the object
			if ($do_utr) {
				my $utr_exon = Bio::SeqFeature::Lite->new(
					-seq_id        => $linedata->{chrom},
					-source        => $source,
					-start         => $cdsStop + 1,
					-end           => $ends[$i],
					-strand        => -1,
					-phase         => '.',
					-primary_tag   => 'five_prime_UTR',
					-primary_id    => $transcript->primary_id . ".utr$i",
					# -display_name  => $transcript->display_name . ".utr$i",
				);
				
				# associate add the utr half to the parent transcript
				$transcript->add_SeqFeature($utr_exon);
			}
		}
		
		#### 5'UTR only ####
		elsif (
			$starts[$i] > $cdsStop 
			and
			$ends[$i] > $cdsStop
		) {
			# the exon start/end is entirely after the cdsStop
			# we have a 5'UTR
			
			# build the utr object
			if ($do_utr) {
				$exon = Bio::SeqFeature::Lite->new(
					-seq_id        => $linedata->{chrom},
					-source        => $source,
					-start         => $starts[$i],
					-end           => $ends[$i],
					-strand        => -1,
					-phase         => '.',
					-primary_tag   => 'five_prime_UTR',
					-primary_id    => $transcript->primary_id . ".utr$i",
					# -display_name  => $transcript->display_name . ".utr$i",
				);
			}
		}
		
		#### Something's wrong ####
		else {
			# just in case I goofed something up
			die " programming error! the exon coordinates don't match up with" .
				" CDS coordinates for transcript $linedata->{name}\n " . 
				" Reverse strand exon coordinates $starts[$i] .. $ends[$i]\n" .
				" CDS coordinates $cdsStart .. $cdsStop\n";
		}
		
		
		# check that we generated an exon
		next unless $exon;
		
		# associate exon with the parent transcript
		$transcript->add_SeqFeature($exon);
		
	}
	
	# done with the exon list
}




sub process_noncoding_exons {
	my ($transcript, $linedata) = @_;
	
	# strip trailing comma from exon lists
	$linedata->{exonStarts} =~ s/,$//; 
	$linedata->{exonEnds} =~ s/,$//; 
	
	# get the exon coordinates
	my @starts = split /,/, $linedata->{exonStarts};
	my @ends = split /,/, $linedata->{exonEnds};
	unless (scalar @starts == scalar @ends) {
		die " source table error! not equal number of exon starts" .
			" and stops for transcript $linedata->{name}\n";
	}
	
	# adjust for 0-based
	@starts = map {$_ + 1} @starts;
	
	# process the exons
	for (my $i = 0; $i <= $#starts; $i++) {
		
		# build the exon SeqFeature object
		my $exon = Bio::SeqFeature::Lite->new(
			-seq_id        => $linedata->{chrom},
			-source        => $source,
			-start         => $starts[$i],
			-end           => $ends[$i],
			-strand        => $transcript->strand,
			-phase         => '.',
			-primary_tag   => 'exon',
			-primary_id    => $transcript->primary_id . ".$i",
			# -display_name  => $transcript->display_name . ".$i",
		);
		
		# associate exon with the parent transcript
		$transcript->add_SeqFeature($exon);
	
	# done with the exon list
	}	
}


sub add_exons_codons_subfeatures {
	my ($transcript, $starts, $ends, $cdsStart, $cdsStop) = @_;
	
	# generate the start and stop codons
	if ($do_codon) {
		if ($transcript->strand == 1) {
			# forward strand
			
			# start codon
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
					-seq_id        => $transcript->seq_id,
					-source        => $transcript->source,
					-primary_tag   => 'start_codon',
					-start         => $cdsStart,
					-end           => $cdsStart + 2,
					-strand        => 1,
					-phase         => 0,
					-primary_id    => $transcript->primary_id . '.start_codon',
					# -display_name  => $transcript->display_name . '.start_codon',
				)
			);
			
			# stop codon
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
					-seq_id        => $transcript->seq_id,
					-source        => $transcript->source,
					-primary_tag   => 'stop_codon',
					-start         => $cdsStop - 2,
					-end           => $cdsStop,
					-strand        => 1,
					-phase         => 0,
					-primary_id    => $transcript->primary_id . '.stop_codon',
					# -display_name  => $transcript->display_name . '.stop_codon',
				)
			);
		}
		
		else {
			# reverse strand
			
			# stop codon
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
					-seq_id        => $transcript->seq_id,
					-source        => $transcript->source,
					-primary_tag   => 'stop_codon',
					-start         => $cdsStart,
					-end           => $cdsStart + 2,
					-strand        => -1,
					-phase         => 0,
					-primary_id    => $transcript->primary_id . '.stop_codon',
					# -display_name  => $transcript->display_name . '.stop_codon',
				)
			);
			
			# start codon
			$transcript->add_SeqFeature( 
				Bio::SeqFeature::Lite->new(
					-seq_id        => $transcript->seq_id,
					-source        => $transcript->source,
					-primary_tag   => 'start_codon',
					-start         => $cdsStop - 2,
					-end           => $cdsStop,
					-strand        => -1,
					-phase         => 0,
					-primary_id    => $transcript->primary_id . '.start_codon',
					# -display_name  => $transcript->display_name . '.start_codon',
				)
			);
		}
	}
	
	# add exons
	for (my $i = 0; $i < scalar(@$starts); $i++) {
		
		# add the exon as a subfeature
		$transcript->add_SeqFeature( 
			Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'exon',
				-start         => $starts->[$i],
				-end           => $ends->[$i],
				-strand        => $transcript->strand,
				-primary_id    => $transcript->primary_id . ".exon$i",
				# -display_name  => $transcript->display_name . ".$i",
			)
		);
		
	}
	
}



sub print_chromosomes {
	
	my $chromo_fh = open_to_read_fh($chromof) or die 
		"unable to open specified chromosome file '$chromof'!\n";
	
	# convert the chromosomes into GFF features
	while (my $line = $chromo_fh->getline) {
		next if ($line =~ /^#/);
		chomp $line;
		my ($chr, $end, $path) = split /\t/, $line;
		unless (defined $chr and $end =~ m/^\d+$/) {
			die " format of chromsome doesn't seem right! Are you sure?\n";
		}
		
		# generate seqfeature
		my $chrom = Bio::SeqFeature::Lite->new(
			-seq_id        => $chr,
			-source        => $source,
			-primary_tag   => $chr =~ m/scaffold/i ? 'scaffold' : 'chromosome',
			-start         => 1,
			-end           => $end,
			-primary_id    => $chr,
			-display_name  => $chr,
		);
		
		# print the gff
		$chrom->version(3);
		$gff_fh->print( $chrom->gff_string . "\n" );
	}
	
	# finished
	$gff_fh->print( "###\n" );
	$chromo_fh->close;
}




__END__

=head1 NAME ucsc_table2gff3.pl



=head1 SYNOPSIS

   ucsc_table2gff3.pl [--options] --table <file>
  
  Options:
  --table <filename>
  --status <filename>
  --sum <filename>
  --chromo <filename>
  --source <text>
  --(no)gene
  --(no)utr
  --codon
  --out <filename>
  --(no)gz
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --table <filename>

Provide the name of a UCSC gene table. It may be obtained preferably
through the UCSC Genome Table Browser, or downloaded from their FTP site.
The file should have specific columns, including txStart, txEnd, cdsStart,
cdsEnd, exonStarts, exonEnds, and others. The conversion may fail if
certain columns are not found. The refSeq and ensGene tables should work
great. The file may be gzipped.

=item --status <filename>

Optionally provide the name of the RefSeq Status file. It may be obtained
through the UCSC Genome Table Browser, or downloaded from their FTP site.
The file should have 3 columns, including 'mrnaAcc', 'status', and 'mol'.
It will only really work with RefSeq tables. The file may be gzipped.

=item --sum <filename>

Optionally provide the name of the RefSeq Summary file. It may be obtained
through the UCSC Genome Table Browser, or downloaded from their FTP site.
The file should have 3 columns, including 'mrnaAcc', 'completeness', and
'summary'. It will only really work with RefSeq tables. The file may be
gzipped.

=item --chromo <filename>

Optionally provide the name of the chromInfo text file. Chromosome 
and/or scaffold features will then be written at the beginning of the 
output GFF file. The chromInfo file may be obtained from the UCSC Downloads 
section or FTP site. The file may be gzipped.

=item --source <text>

Optionally provide the text to be used as the GFF source. The default is 
'UCSC'.

=item --(no)gene

Specify whether (or not) to assemble mRNA transcripts into genes. This 
will create the canonical gene->mRNA->(exon,CDS) heirarchical structure. 
Otherwise, mRNA transcripts are kept independent. The gene names, when 
available, are always associated with transcripts through the Alias tag. 
The default is true.

=item --(no)utr

Specify whether (or not) to include three_prime_utr and five_prime_utr 
features in the transcript heirarchy. If not defined, the GFF interpreter 
must infer the UTRs from the CDS and exon features. The default is true.

=item --codon

Specify whether (or not) to include start_codon and stop_codon features 
in the transcript heirarchy. The default is false.

=item --out <filename>

Optionally specify a new filename. By default it uses the basename of the 
table file.

=item --(no)gz

Specify whether the output file should be compressed with gzip.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will convert a UCSC gene table file into a GFF3 format file 
suitable for loading into GBrowse. Unlike other simple converters out there, 
it will attempt to build complex gene structures with gene, transcripts, 
UTRs, and exons. It will attempt to identify non-coding genes as to type 
using the gene name as inference. For RefSeq tables, it can also use 
additional data tables to identify the feature type, gene status, and include
gene descriptions.

Files should preferentially be downloaded through the UCSC table browser, if 
at all possible, for two good reasons. One, the table browser includes 
column headings, allowing for accurate identification of the columns and 
verification of the right file. Two, the genes are ordered, or at least 
grouped by chromosome, which is not the case with files downloaded from the 
FTP site. 

The program works best with the refGene (RefSeq Genes), ensGene 
(Ensembl Genes), and xenoRefGene tables. Other tables may work, if the 
appropriate columns can be identified in the column headings, but your results 
may vary and are not guaranteed. 

If provided, chromosome and/or scaffold features may also be written at the 
beginning of the GFF file. 

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
