#!/usr/bin/perl

# convert ucsc refseq table file to gff3
# why am I doing this anyway???? for fun, experience, completeness....
# the other options out there are extremely limited....

# the original version of this script used SeqFeature::Annotated and 
# SeqFeatureIO. However, it made the script pretty glacial with all the 
# excess object overhead, particularly the Ontology stuff

# this version uses much, much, simpler SeqFeature::Generic and Tools::GFF 
# and is over 2 orders of magnitude faster than the original

use strict;
#use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::Annotation::Comment;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	find_column_index
);
use tim_file_helper qw(
	open_tim_data_file
	open_to_write_fh
);
#use Data::Dumper;

print "\n A script to convert UCSC RefSeq gene tables to GFF3\n\n";





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
	$source,
	$link_parent,
	$valid_tags,
	$outfile,
#	$gz,
	$help, 
);
GetOptions( 
	'table=s'    => \$genetablef, # the input refseq file
	'status=s'   => \$refseqstatusf, # the refseqstatus file
	'sum=s'      => \$refseqsumf, # the refseqsummary file
	'source=s'   => \$source, # the GFF source
	'parent!'    => \$link_parent, # link transcript to parent
	'valid!'     => \$valid_tags, # only include valid tags
	'out=s'      => \$outfile, # output file name
#	'gz!'        => \$gz, # compress file
	'help'       => \$help, # request help
);

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}





### Check requirements and defaults
unless ($genetablef) {
	die " Gene table file not specified!\n";
}
unless ($source) {
	$source = 'UCSC';
}
unless (defined $link_parent) {
	# do not link the transcript to the parent gene
	$link_parent = 0;
}
unless (defined $valid_tags) {
	# include non valid tags too
	$valid_tags = 0;
}
my $start_time = time;




### Open files

# input files
my $refseqsum = load_refseq_summary_data($refseqsumf);

my $refseqstat = load_refseq_status_data($refseqstatusf);

my ($gene_fh, $gene_metadata) = open_tim_data_file($genetablef) or
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
	$outfile = $gene_metadata->{'basename'} . '.gff';
}
my $gff_fh = open_to_write_fh($outfile);
	# I'm not allowed to open a gzip compressed output file because 
	# Bio::Tools::GFF will not recognize the file object
	# I could filter through an external gzip but I need access to the raw
	# file handle to write comments and directives and Tools::GFF doesn't 
	# give me that access
my $gff =Bio::Tools::GFF->new(
	-fh           => $gff_fh,
	-gff_version  => 3,
);
unless ($gff) {
	die " unable to open output GFF3 file!\n";
}




### Identify indices
my (
	$chr_i, 
	$str_i,	
	$txStart_i, 
	$txEnd_i, 
	$cdStart_i, 
	$cdEnd_i, 
	$exCount_i,
	$exStarts_i,
	$exEnds_i,
	$exFrames_i,
	$trnscpt_i,
	$gene_i,
) = identify_indices($gene_metadata);





### Process the files
# print comments
print {$gff_fh} "##gff-version 3\n";
$gff->{'_first'} = 0; # pretend we've already written the first feature
					# this avoids writing two gff declaration lines
print {$gff_fh} "# UCSC gene table file $genetablef\n";
print {$gff_fh} "# UCSC RefSeq Status file $refseqstatusf\n" if $refseqstatusf;
print {$gff_fh} "# UCSC RefSeq Summary file $refseqsumf\n" if $refseqsumf;

print " Converting gene table features....\n";
my $count = process_gene_table();
print " Converted $count->[0] gene features\n" if $count->[0] > 0;
print " Converted $count->[1] mRNA transcripts\n" if $count->[1] > 0;
print " Converted $count->[2] pseudogene transcripts\n" if $count->[2] > 0;
print " Converted $count->[3] ncRNA transcripts\n" if $count->[3] > 0;
print " Converted $count->[4] miRNA transcripts\n" if $count->[4] > 0;
print " Converted $count->[5] snRNA transcripts\n" if $count->[5] > 0;
print " Converted $count->[6] snoRNA transcripts\n" if $count->[6] > 0;
print " Converted $count->[7] tRNA transcripts\n" if $count->[7] > 0;
print " Converted $count->[8] rRNA transcripts\n" if $count->[8] > 0;
print " Converted $count->[9] other transcripts\n" if $count->[9] > 0;





### Finish
print " completed. wrote file '$outfile'\n";
printf " finished in %.1f minutes\n", (time - $start_time)/60;

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
		print " loaded ", scalar(keys %sumdata), " gene summaries from '$file'\n";
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
		print " loaded ", scalar(keys %statdata), " RNA status from '$file'\n";
		$fh->close;
	}
	
	return \%statdata;
}




sub identify_indices {
	
	my $metadata = shift;
	my @indices;
	
	foreach (
		['chrom', 2],
		['strand', 3],
		['txStart', 4],
		['txEnd', 5],
		['cdsStart', 6],
		['cdsEnd', 7],
		['exonCount', 8],
		['exonStarts', 9],
		['exonEnds', 10],
		['exonFrames', 15],
		['name', 1]
	) {
		# we need to identify which columns are which,
		# but this is a huge headache when the column headers are not present
		# as when the file was downloaded from UCSC FTP site rather than going
		# through the table browser (which adds the column headers as a comment 
		# line)
		# therefore, we need to make assumptions, and hope that the assumption
		# is valid, otherwise bad things will happen
		
		# attempt to identify
		my $index = find_column_index($metadata, $_->[0]);
		
		# verify
		unless (defined $index) {
			warn " unable to identify column index for '$_->[0]'," .
				" assuming default index of $_->[1]\n";
			if ($_->[1] < $metadata->{'number_columns'}) {
				$index = $_->[1];
			}
			else {
				die " cannot assign '$_->[0]' to an index!\n";
			}
		}
		push @indices, $index;
	}
	
	# identify the gene index, if unable simply re-use name
	my $gene_index = find_column_index($metadata, 'name2');
	unless (defined $gene_index) {
		if ($metadata->{'number_columns'} == 16) {
			# table has the same column number as a standard gene table
			# make assumption
			warn " unable to identify gene name column index (name2), " . 
				"using default index 12\n";
			$gene_index = 12;
		}
		else {
			# non-standard table, can only hope this works
			warn " unable to identify gene name column index (name2), using name\n";
			$gene_index = $indices[9];
		}
	}
	push @indices, $gene_index;
	
	return @indices;
	
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
		# 11 id
		# 12 name2
		# 13 cdsStartStat
		# 14 cdsEndStat
		# 15 exonFrames
}





sub process_gene_table {
	
	# initialize 
	my $current_chrom;
	my $gene_count        = 0;
	my $mrna_count        = 0;
	my $pseudogene_count  = 0;
	my $ncrna_count       = 0;
	my $mirna_count       = 0;
	my $snrna_count       = 0;
	my $snorna_count      = 0;
	my $trna_count        = 0;
	my $rrna_count        = 0;
	my $other_count       = 0;
	my %gene2seqf;
	my %transcript_check;
	my @letters = qw(A B C D E F G H I J K L M N O P Q S T U V W X Y Z);
	
	#### Main Gene Loop
		# we're not using a defined for loop since we may have to read 
		# multiple lines for multi-mRNA genes
	while (my $line = $gene_fh->getline) {
		
		# process the row from the gene table
		chomp $line;
		my $linedata = [ split /\t/, $line ];
		
		# check the chromosome
		unless (defined $current_chrom) {
			$current_chrom = $linedata->[$chr_i];
		}
		if ($linedata->[$chr_i] ne $current_chrom) {
			# moved on to next chromosome
			
			# print all of the current genes and clear memory
			foreach my $id (
				map $_->[0],
				sort { $a->[1] <=> $b->[1] }
				map [$_, $gene2seqf{$_}->start], 
				keys %gene2seqf
			) {
				# first sort all of the genes by start position
				# count the genes
				$gene_count++;
				
				# count the transcript types
				foreach ( $gene2seqf{$id}->get_SeqFeatures ) {
					my $type = $_->primary_tag;
					if ($type eq 'mRNA') {
						$mrna_count++;
					}
					elsif ($type eq 'pseudogene') {
						$pseudogene_count++;
					}
					elsif ($type eq 'ncRNA') {
						$ncrna_count++;
					}
					elsif ($type eq 'miRNA') {
						$mirna_count++;
					}
					elsif ($type eq 'snRNA') {
						$snrna_count++;
					}
					elsif ($type eq 'snoRNA') {
						$snorna_count++;
					}
					elsif ($type eq 'tRNA') {
						$trna_count++;
					}
					elsif ($type eq 'rRNA') {
						$rrna_count++;
					}
					else {
						$other_count++;
					}
				}
				
				# then print the gene seqfeature object
				print_my_features_and_subfeatures( $gene2seqf{ $id } );
				
				# finally, delete the gene seqfeature object
				delete $gene2seqf{ $id };
				
			}
			
			print {$gff_fh} "###\n"; # directive to close out all previous genes
			$current_chrom = $linedata->[$chr_i];
		}
		
		# keep track of current gene name
		my $current_gene = $linedata->[$gene_i];
		
		# define the strand, must use 1 or -1
		my $strand;
		if ($linedata->[$str_i] eq '+') {
			# forward
			$strand = 1;
		}
		else {
			# reverse
			$strand = -1;
		}
		
		#### Generate the gene object if necessary
		unless (exists $gene2seqf{ $linedata->[$gene_i] }) {
		
			# generate the gene SeqFeature object
			my $gene = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $linedata->[$txStart_i] + 1,
				-end           => $linedata->[$txEnd_i],
				-strand        => $strand,
				-frame         => '.',
				-display_name  => $linedata->[$gene_i],
				-primary_tag   => 'gene',
			);
			
			
			# add ID
			$gene->add_tag_value('ID', 
				$linedata->[$chr_i] . '_' . $linedata->[$gene_i]);
				# we are prepending the chromosome name to the gene name to 
				# ensure this ID is unique in the database
				# occasionally sections of chromosomes are duplicated and 
				# treated as a separate chromosome, which could lead to 
				# gene duplications
				# for example, in hg18 there is chr6 and chr6_qbl_hap2
	
			# add Name, since display_name above doesn't seem to be doing it for me
			$gene->add_tag_value('Name', $linedata->[$gene_i]);
	
			
			# add status 
			if (exists $refseqstat->{  $linedata->[$trnscpt_i] } ) {
				
				# add the status tag unless we want only pure valid tags
				unless ($valid_tags) {
					$gene->add_tag_value('status', 
						$refseqstat->{ $linedata->[$trnscpt_i] }->[0] );
				}
				
			}
			
			# add Note if possible
			if (
				exists $refseqsum->{ $linedata->[$trnscpt_i] }
				and
				defined $refseqsum->{ $linedata->[$trnscpt_i] }->[1]
			) {
				# add the the note
				my $text = $refseqsum->{ $linedata->[$trnscpt_i] }->[1];
				my $note = Bio::Annotation::Comment->new( 
					-text => $text,
				);
				$gene->add_tag_value('Note',$note);
			}
			
			# Store the gene oject into the gene hash
			$gene2seqf{ $linedata->[$gene_i] } = $gene;
			
		} # finished building the gene object
		
		
		
		#### Process the Transcript
		
		# Pull out the gene SeqFeature object
		my $gene = $gene2seqf{ $linedata->[$gene_i] };
		
		# check the gene coordinates are current
		if ($linedata->[$txStart_i] < $gene->start) {
			# update the txstart position
			$gene->start( $linedata->[$txStart_i] );
		}
		if ($linedata->[$txEnd_i] > $gene->end) {
			# update the txstop position
			$gene->end( $linedata->[$txEnd_i] );
		}
		
		# Check for note
		unless ($gene->has_tag('Note')) {
			# unless it somehow wasn't available for a previous transcript, 
			# but is now, we'll add it now
			# we won't check if transcripts for the same gene have the 
			# same note or not, why wouldn't they????
			if (
				exists $refseqsum->{ $linedata->[$trnscpt_i] }
				and
				defined $refseqsum->{ $linedata->[$trnscpt_i] }->[1]
			) {
				# check for a summary for this mRNA
				
				# make note if it exists
				my $text = $refseqsum->{ $linedata->[$trnscpt_i] }->[1];
				my $note = Bio::Annotation::Comment->new( 
					-value => $text,
				);
				$gene->add_tag_value('Note',$note);;
			}
		} # finished with the note
		
		
		
		# build the transcript SeqFeature 
		my $transcript = Bio::SeqFeature::Generic->new(
			-seq_id        => $linedata->[$chr_i],
			-source        => $source,
			-start         => $linedata->[$txStart_i] + 1,
			-end           => $linedata->[$txEnd_i],
			-strand        => $strand,
			-frame         => '.',
			-display_name  => $linedata->[$trnscpt_i],
		);
		
		# Attempt to identify the transcript type
		if (
			$linedata->[$cdStart_i] == $linedata->[$txEnd_i] and 
			$linedata->[$cdEnd_i] == $linedata->[$txEnd_i] 
		) {
			# there appears to be no coding potential when 
			# txEnd = cdsStart = cdsEnd
			# if you'll look, all of the exon phases should also be -1
			
			# if we are using a RefSeq gene table, we may be able to infer
			# some certain types from the gene name
			
			# in the mean time, we'll try to deduce from the gene name
			if ($linedata->[$gene_i] =~ /^LOC\d+/) {
				# empirical testing seems to suggest that all the noncoding 
				# genes with a name like LOC123456 are pseudogenes
				# well, at least with hg18, it may not be true for others
				$transcript->primary_tag('pseudogene');
			}
			elsif ($linedata->[$gene_i] =~ /^mir/i) {
				# a noncoding gene whose name begins with mir is likely a 
				# a micro RNA
				$transcript->primary_tag('miRNA');
			}
			elsif ($linedata->[$gene_i] =~ /^snr/i) {
				# a noncoding gene whose name begins with snr is likely a 
				# a snRNA
				$transcript->primary_tag('snRNA');
			}
			elsif ($linedata->[$gene_i] =~ /^sno/i) {
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
		
		
		# add the parent if requested
		if ($link_parent) {
			# adding the parent will actually mess up display of the 
			# transcript using the process_transcript glyph in GBrowse
			# so it's best not to link it
			# It will however be linked to the gene for purposes of this 
			# script only
			my ($parent_tag) = $gene->get_tag_values('ID');
			$transcript->add_tag_value('Parent', $parent_tag);
		}
		
		# add the id
		{
			# To add insult to injury, there are sometimes
			# non-identical transcripts with the same ID!!!!!
			# so we really have to ensure we have unique transcript IDs
			
			# This also occurs sometimes when subsections of a chromosome are 
			# included as a separate chromosome, for example in hg18 there is
			# chr6 and chr6_cox_hap1 which have identical genes
			
			# generate unique transcript id
			my $transcript_id = $linedata->[$trnscpt_i];
			my $i = 0;
			while (exists $transcript_check{$transcript_id}) {
				# we'll append a letter to the id until we have a unique id
				$transcript_id = $linedata->[$trnscpt_i] . $letters[$i];
				$i++;
			}
			
			# the id is (finally) unique, add it
			$transcript_check{ $transcript_id } = 1;
			$transcript->add_tag_value('ID', $transcript_id);
		}
			
			
		
		# add the name, since display_name doesn't seem to be working for me
		$transcript->add_tag_value('Name', $linedata->[$trnscpt_i]);
		
		# add extra non-valid tag information
		unless ($valid_tags) {
			# add a status for the transcript
			if (exists $refseqstat->{ $linedata->[$trnscpt_i] } ) {
				$transcript->add_tag_value('status', 
					$refseqstat->{ $linedata->[$trnscpt_i] }->[0] );
			}
			# add the completeness value for the tag
			if (exists $refseqsum->{ $linedata->[$trnscpt_i] } ) {
				$transcript->add_tag_value('completeness', 
					$refseqsum->{ $linedata->[$trnscpt_i] }->[0] );
			}
		}
		
		# Add the exons
		if ($linedata->[$exCount_i] > 1) { 
			
			# determine how to process the exons by the feature type 
			# and the strand orientation
			
			if ($transcript->primary_tag eq 'mRNA' and $strand > 0) {
				# forward strand coding exon
				process_forward_exons($transcript, $linedata);
			}
			elsif ($transcript->primary_tag eq 'mRNA' and $strand < 0) {
				# reverse strand coding exon
				process_reverse_exons($transcript, $linedata);
			}
			
			# otherwise assume exons are part of noncoding transcript
			else {
				process_noncoding_exons($transcript, $linedata);
			}
			
		}
		
		# this transcript is now complete
		# associate with the gene
		$gene->add_SeqFeature($transcript);
		
	}
	
	
	
	#### Finished
	
	# print all of the current genes
	foreach my $id (
		map $_->[0],
		sort { $a->[1] <=> $b->[1] }
		map [$_, $gene2seqf{$_}->start], 
		keys %gene2seqf
	) {
		# first sort all of the genes by start position
		
		# count the genes
		$gene_count++;
		
		# count the transcript types
		foreach ( $gene2seqf{$id}->get_SeqFeatures ) {
			my $type = $_->primary_tag;
			if ($type eq 'mRNA') {
				$mrna_count++;
			}
			elsif ($type eq 'pseudogene') {
				$pseudogene_count++;
			}
			elsif ($type eq 'ncRNA') {
				$ncrna_count++;
			}
			elsif ($type eq 'miRNA') {
				$mirna_count++;
			}
			elsif ($type eq 'snRNA') {
				$snrna_count++;
			}
			elsif ($type eq 'snoRNA') {
				$snorna_count++;
			}
			elsif ($type eq 'tRNA') {
				$trna_count++;
			}
			elsif ($type eq 'rRNA') {
				$rrna_count++;
			}
			else {
				$other_count++;
			}
		}
		
		# then print the gene seqfeature object
		print_my_features_and_subfeatures( $gene2seqf{ $id } );
		
		# finally, delete the gene seqfeature object
		delete $gene2seqf{ $id };
	}
	print {$gff_fh} "###\n"; # directive to close out all previous genes
	return [
		$gene_count,
		$mrna_count,
		$pseudogene_count,
		$ncrna_count,
		$mirna_count,
		$snrna_count,
		$snorna_count,
		$trna_count,
		$rrna_count,
		$other_count
	];
;
}



sub process_forward_exons {
	my ($transcript, $linedata) = @_;
	
	# strip trailing comma from exon lists
	$linedata->[$exStarts_i] =~ s/,$//; 
	$linedata->[$exEnds_i] =~ s/,$//; 
	$linedata->[$exFrames_i] =~ s/,$//;
	
	# get the exon coordinates
	my @starts = split /,/, $linedata->[$exStarts_i];
	my @ends = split /,/, $linedata->[$exEnds_i];
	my @phases = split /,/, $linedata->[$exFrames_i];
	unless (scalar @starts == scalar @ends) {
		die " source table error! not equal number of exon starts" .
			" for transcript $linedata->[$trnscpt_i]\n";
	}
	my $cdsStart = $linedata->[$cdStart_i] + 1;
	my $cdsStop = $linedata->[$cdEnd_i];
	
	# adjust for 0-based
	@starts = map {$_ + 1} @starts;
	
	# get the parent id
	my ($parent_tag) = $transcript->get_tag_values('ID');
	
	# process the exons
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
			$exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $starts[$i],
				-end           => $ends[$i],
				-strand        => 1,
				-frame         => '.',
				-primary_tag   => 'five_prime_UTR',
			);
			
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
			my $utr_exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $starts[$i],
				-end           => $cdsStart - 1,
				-strand        => 1,
				-frame         => '.',
				-primary_tag   => 'five_prime_UTR',
			);
			# since we're actually building two objects here, we have to 
			# complete the process of building the utr object and associate 
			# it with the transcript object
			# the cds half of the exon will be finished below
			
			# add the utr half id
			$utr_exon->add_tag_value('ID', $parent_tag . ".$i" . 'u');
			
			# add the utr half parent
			$utr_exon->add_tag_value('Parent', $parent_tag);
	
			# associate add the utr half to the parent transcript
			$transcript->add_SeqFeature($utr_exon);
			
			# now build the cds half of the object
			$exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $cdsStart,
				-end           => $ends[$i],
				-strand        => 1,
				-frame         => $phases[$i],
				-primary_tag   => 'CDS',
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
			
			$exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $starts[$i],
				-end           => $ends[$i],
				-strand        => 1,
				-frame         => $phases[$i],
				-primary_tag   => 'CDS',
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
			$exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $starts[$i],
				-end           => $cdsStop,
				-strand        => 1,
				-frame         => $phases[$i],
				-primary_tag   => 'CDS',
			);
			
			# now build the utr half of the object
			my $utr_exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $cdsStop + 1,
				-end           => $ends[$i],
				-strand        => 1,
				-frame         => '.',
				-primary_tag   => 'three_prime_UTR',
			);
			
			# add the utr half id
			$utr_exon->add_tag_value('ID', $parent_tag . ".$i" . 'u');
			
			# add the utr half parent
			$utr_exon->add_tag_value('Parent', $parent_tag);
	
			# associate add the utr half to the parent transcript
			$transcript->add_SeqFeature($utr_exon);
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
			$exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $starts[$i],
				-end           => $ends[$i],
				-strand        => 1,
				-frame         => '.',
				-primary_tag   => 'three_prime_UTR',
			);
		}
		
		#### Something's wrong ####
		else {
			# just in case I goofed something up
			die " programming error! the exon coordinates don't match up with" .
				" CDS coordinates for transcript $linedata->[$trnscpt_i]\n " . 
				" Forward strand exon coordinates $starts[$i] .. $ends[$i]\n" .
				" CDS coordinates $cdsStart .. $cdsStop\n";
		}
		
		
		# add id
		$exon->add_tag_value('ID', $parent_tag . ".$i");
		
		# add parent
		$exon->add_tag_value('Parent', $parent_tag);

		# associate exon with the parent transcript
		$transcript->add_SeqFeature($exon);
		
	}
	
	# done with the exon list
}




sub process_reverse_exons {
	my ($transcript, $linedata) = @_;
		
	# strip trailing comma from exon lists
	$linedata->[$exStarts_i] =~ s/,$//; 
	$linedata->[$exEnds_i] =~ s/,$//; 
	$linedata->[$exFrames_i] =~ s/,$//;
	
	# get the exon coordinates
	my @starts = split /,/, $linedata->[$exStarts_i];
	my @ends = split /,/, $linedata->[$exEnds_i];
	my @phases = split /,/, $linedata->[$exFrames_i];
	unless (scalar @starts == scalar @ends) {
		die " source table error! not equal number of exon starts" .
			" and stops for transcript $linedata->[$trnscpt_i]\n";
	}
	my $cdsStart = $linedata->[$cdStart_i] + 1;
	my $cdsStop = $linedata->[$cdEnd_i];
	
	# adjust for 0-based
	@starts = map {$_ + 1} @starts;
	
	# get the parent id
	my ($parent_tag) = $transcript->get_tag_values('ID');
	
	# process the exons
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
			$exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $starts[$i],
				-end           => $ends[$i],
				-strand        => -1,
				-frame         => '.',
				-primary_tag   => 'three_prime_UTR',
			);
			
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
			my $utr_exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $starts[$i],
				-end           => $cdsStart - 1,
				-strand        => -1,
				-frame         => '.',
				-primary_tag   => 'three_prime_UTR',
			);
			# since we're actually building two objects here, we have to 
			# complete the process of building the utr object and associate 
			# it with the transcript object
			# the cds half of the exon will be finished below
			
			# add the utr half id
			$utr_exon->add_tag_value('ID', $parent_tag . ".$i" . 'u');
			
			# add the utr half parent
			$utr_exon->add_tag_value('Parent', $parent_tag);
	
			# associate add the utr half to the parent transcript
			$transcript->add_SeqFeature($utr_exon);
			
			# now build the cds half of the object
			$exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $cdsStart,
				-end           => $ends[$i],
				-strand        => -1,
				-frame         => $phases[$i],
				-primary_tag   => 'CDS',
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
			
			$exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $starts[$i],
				-end           => $ends[$i],
				-strand        => -1,
				-frame         => $phases[$i],
				-primary_tag   => 'CDS',
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
			$exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $starts[$i],
				-end           => $cdsStop,
				-strand        => -1,
				-frame         => $phases[$i],
				-primary_tag   => 'CDS',
			);
			
			# now build the utr half of the object
			my $utr_exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $cdsStop + 1,
				-end           => $ends[$i],
				-strand        => -1,
				-frame         => '.',
				-primary_tag   => 'five_prime_UTR',
			);
			
			# add the utr half id
			$utr_exon->add_tag_value('ID', $parent_tag . ".$i" . 'u');
			
			# add the utr half parent
			$utr_exon->add_tag_value('Parent', $parent_tag);
	
			# associate add the utr half to the parent transcript
			$transcript->add_SeqFeature($utr_exon);
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
			$exon = Bio::SeqFeature::Generic->new(
				-seq_id        => $linedata->[$chr_i],
				-source        => $source,
				-start         => $starts[$i],
				-end           => $ends[$i],
				-strand        => -1,
				-frame         => '.',
				-primary_tag   => 'five_prime_UTR',
			);
		}
		
		#### Something's wrong ####
		else {
			# just in case I goofed something up
			die " programming error! the exon coordinates don't match up with" .
				" CDS coordinates for transcript $linedata->[$trnscpt_i]\n " . 
				" Reverse strand exon coordinates $starts[$i] .. $ends[$i]\n" .
				" CDS coordinates $cdsStart .. $cdsStop\n";
		}
		
		
		# add id
		$exon->add_tag_value('ID', $parent_tag . ".$i");
		
		# add parent
		$exon->add_tag_value('Parent', $parent_tag);

		# associate exon with the parent transcript
		$transcript->add_SeqFeature($exon);
		
	}
	
	# done with the exon list
}




sub process_noncoding_exons {
	my ($transcript, $linedata) = @_;
	
	# strip trailing comma from exon lists
	$linedata->[$exStarts_i] =~ s/,$//; 
	$linedata->[$exEnds_i] =~ s/,$//; 
	
	# get the exon coordinates
	my @starts = split /,/, $linedata->[$exStarts_i];
	my @ends = split /,/, $linedata->[$exEnds_i];
	unless (scalar @starts == scalar @ends) {
		die " source table error! not equal number of exon starts" .
			" and stops for transcript $linedata->[$trnscpt_i]\n";
	}
	
	# adjust for 0-based
	@starts = map {$_ + 1} @starts;
	
	# get the parent id
	my ($parent_tag) = $transcript->get_tag_values('ID');
	
	# process the exons
	for (my $i = 0; $i <= $#starts; $i++) {
		
		# build the exon SeqFeature object
		my $exon = Bio::SeqFeature::Generic->new(
			-seq_id        => $linedata->[$chr_i],
			-source        => $source,
			-start         => $starts[$i],
			-end           => $ends[$i],
			-strand        => 1,
			-frame         => '.',
			-primary_tag   => 'exon',
		);
		
		# add id
		$exon->add_tag_value('ID', $parent_tag . ".$i");
		
		# add parent
		$exon->add_tag_value('Parent', $parent_tag);

		# associate exon with the parent transcript
		$transcript->add_SeqFeature($exon);
	
	# done with the exon list
	}	
}




sub print_my_features_and_subfeatures {
	# Bio::Tools::GFF will not recursively print all subfeatures
	# so we fake it here
	
	my $feat = shift;
	
	# print the passed SeqFeature object
	$gff->write_feature($feat);
	
	# print the subfeatures
	foreach ($feat->get_SeqFeatures) {
		print_my_features_and_subfeatures($_);
	}
}







__END__

=head1 NAME ucsc_refseq2gff3.pl



=head1 SYNOPSIS

   ucsc_table2gff3.pl [--options] --table <file>
  
  Options:
  --table <filename>
  --status <filename>
  --sum <filename>
  --source <text>
  --(no)parent
  --(no)valid
  --out <filename>
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --table <filename>

Provide the name of a UCSC gene table. It may be obtained preferably
through the UCSC Genome Table Browser, or downloaded from their FTP site.
The file should have specific columns, including txStart, txEnd, cdsStart,
cdsEnd, exonStarts, exonEnds, and others. The conversion may fail if
certain columns are not found. The RefSeq and ensGene tables should work
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

=item --source <text>

Optionally provide the text to be used as the GFF source. The default is 
'UCSC'.

=item --(no)parent

Indicate whether or not transcript features should be linked to a parent 
gene feature. While linking the transcript feature to a parent gene feature 
may generate a more complete gene object model, the processed_transcript 
glyph in GBrowse cannot handle this object. Without linking, two objects will 
be present in a SeqFeature::Store database, the gene feature (with no 
children) and the transcript feature (with exon and UTR children). The 
default behavior is false or no linking.

=item --(no)valid

Indicate whether or not valid GFF tags should only be included. When the 
refSeq Summary and/or refSeq Status files are included, two additional 
non-valid GFF tags may be included, 'completeness' as indicated in the 
refSeq Summary file, and 'status' in the refSeq Status file. The presence 
of these tags may technically invalidate proper GFF structure, but are 
tolerated by some programs (including Bio::SeqFeature::Store databases) 
and may prove useful in analysis. The default value is false or to include 
all available tags.

=item --out <filename>

Optionally specify a new filename. By default it uses the basename of the 
table file.

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


=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

























