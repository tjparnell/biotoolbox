#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::db_helper::gff3_parser;
use Bio::ToolBox::file_helper qw(
	open_to_read_fh
	open_to_write_fh
);

my $VERSION = '1.15';


print "\n This program will convert a GFF3 file to UCSC gene table\n";

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
	$outfile,
	$use_alias,
	$gz,
	$verbose,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the gff3 data file
	'out=s'     => \$outfile, # name of output file 
	'alias!'    => \$use_alias, # append aliases to geneName
	'gz!'       => \$gz, # compress output
	'verbose!'  => \$verbose, # verbose print statements
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
	print " Biotoolbox script gff3_to_ucsc_table.pl, version $VERSION\n\n";
	exit;
}



### Check for requirements
unless ($infile) {
	$infile = shift @ARGV or
		die " no input file! use --help for more information\n";
}
unless ($infile =~ /\.gff 3? (?: \.gz )? $/xi) {
	die " input file doesn't have a gff extension! Is this a GFF3 file?\n";
}

if ($outfile) {
	# be nice and add an extension for them if it's missing
	unless ($outfile =~ /\.txt$/) {
		$outfile .= '.txt';
	}
}
else {
	# define a new output file name for them
	$outfile = $infile;
	$outfile =~ s/\.gff 3? (?: \.gz )? $/.refFlat.txt/xi;
}
unless (defined $gz) {
	# mimic the input file as far as compression is concerned
	$gz = 1 if $infile =~ m/\.gz$/i;
}

### Open the input and ouput files
# Open the input GFF file
my $in_fh = open_to_read_fh($infile) or 
	die " unable to open input file '$infile'!\n";

# Open the output gene table file
my $outfh = open_to_write_fh($outfile, $gz) or 
	die " unable to open file handle for '$outfile'!\n";

# Print headers
$outfh->print( join("\t", qw(#geneName name chrom strand txStart txEnd cdsStart 
	cdsEnd exonCount exonStarts exonEnds) ) . "\n");



### Process the GFF3 table
# initialize global variables
my $count = 0; # a count of written features
my %unknowns; # a hash of unrecognized feature types to report at the end
process_gff_file_to_table();






### Finished

# Close output file
$outfh->close;

# print warnings about unknown feature types
if (%unknowns) {
	foreach (sort {$a cmp $b} keys %unknowns) {
		print " Found $unknowns{$_} unrecognized features of '$_'!\n";
	}
}

print " Finished! Wrote $count features to file '$outfile'\n";






########################   Subroutines   ###################################


sub process_gff_file_to_table {
	
	# Collect the top features for each sequence group.
	# Rather than going feature by feature through the gff,
	# we'll load the top features, which will collect all the features 
	# and assemble appropriate feature -> subfeatures according to the 
	# parent - child attributes.
	# This may (will?) be memory intensive. This can be limited by 
	# including '###' directives in the GFF3 file after each chromosome.
	# This directive tells the parser that all previously opened feature 
	# objects are finished and may be closed.
	# Without the directives, all feature objects loaded from the GFF3 file 
	# will be kept open until the end of the file is reached. 
	
	# open gff3 parser object
	my $parser = Bio::ToolBox::db_helper::gff3_parser->new($infile) or
		die " unable to open input file '$infile'!\n";
	
	
	# process the features
	while (my @top_features = $parser->top_features() ) {
		
		# processing statement
		if ($verbose) {
			print " Collected ", scalar(@top_features), " features from ";
			if ( $top_features[0]->seq_id eq $top_features[-1]->seq_id ) {
				# check whether all the features came from the same top sequence
				print "sequence ID '", $top_features[0]->seq_id, "'\n";
			}
			else {
				print "multiple sequences\n";
			}
		}
		
		# Process the top features
		while (@top_features) {
			my $feature = shift @top_features;
			my $type = lc $feature->primary_tag;
			
			# check for chromosome
			if ($type =~ /chromosome|contig|scaffold|sequence|region/) {
				next;
			}
			
			# process the features
			if ($type eq 'gene' or $type eq 'pseudogene') {
				# a gene object, we will need to process it's transcript subfeatures
				
				# process the gene object
				process_gene($feature);
				
			}
			
			elsif ($type eq 'mrna') {
				# a coding transcript
				
				# first need to initialize the ucsc table object
				my $ucsc_transcript = initialize_transcript($feature);
				
				# process the gene object
				process_transcript($ucsc_transcript, $feature);
			}
			
			elsif ($type =~ /rna|noncoding/) {
				# a non-coding RNA transcript
				
				# first need to initialize the ucsc table object
				my $ucsc_transcript = initialize_transcript($feature);
				
				# process the gene object
				process_nc_transcript($ucsc_transcript, $feature);
			}
			
			elsif ($type =~ /transcript/) {
				# we'll assume this is a protein_coding transcript?
				# hopefully noncoding transcripts will be caught earlier
				
				# first need to initialize the ucsc table object
				my $ucsc_transcript = initialize_transcript($feature);
				
				# process the gene object
				process_transcript($ucsc_transcript, $feature);
			}
			
			else {
				# catchall for unrecognized feature types
				# record and warn at end
				$unknowns{$type} += 1;
				next;
			}
			
		}
		
		# print running total
		print "  wrote $count features so far\n" if $verbose;
	}
}




### Initialize a data structure for the UCSC table line item
sub initialize_transcript {
	# a hash structure storing the fields to be used in an UCSC gene table
	my $feature = shift;
	
	# Fill the data hash
	my %transcript = (
		'name'         => undef,
		'name2'        => undef,
		'chr'          => $feature->seq_id,
		'strand'       => $feature->strand > 0 ? '+' : '-', # let's hope there is no 0
		'txStart'      => undef,
		'txEnd'        => undef,
		'cdsStart'     => undef,
		'cdsEnd'       => undef,
		'exonCount'    => 0,
		'exonStarts'   => [],
		'exonEnds'     => [],
	);
	
	# get score, rarely used, but just in case
	my $score = $feature->score;
	if ($score and $score ne '.') {
		$transcript{score} = $score;
	}
	
	return \%transcript;
}



### process a gene
sub process_gene {
	my $gene = shift;
	
	# we need to process the child subfeature transcripts
	# these should be the first level subfeatures
	# walk through each one
	my $transcript_success = 0;
	my $cds_count = 0;
	foreach my $subfeat ($gene->get_SeqFeatures) {
		
		# check the type, and process accordingly
		my $type = lc($subfeat->primary_tag);
		
		if ($type eq 'mrna') {
			# a protein coding transcript
			
			# initialize the gene table line structure
			my $ucsc = initialize_transcript($subfeat);
			
			# collect the names for the transcript and/or gene
			collect_names($ucsc, $gene, $subfeat);
			
			# process the transcript
			process_transcript($ucsc, $subfeat);
			
			$transcript_success++;
		}
		
		elsif ($type =~ /rna|pseudogene|transcript/) {
			# other types of non-coding RNA species
			# should cover misc_RNA snRNA snoRNA ncRNA rRNA tRNA miRNA etc.
			# also pseudogenes and processed_transcripts that do not have CDS
			
			# initialize the gene table line structure
			my $ucsc = initialize_transcript($subfeat);
			
			# collect the names for the transcript and/or gene
			collect_names($ucsc, $gene, $subfeat);
			
			# process the transcript
			process_nc_transcript($ucsc, $subfeat);
			
			$transcript_success++;
		}
		
		elsif ($subfeat->get_SeqFeatures) {
			# there are some oddball subfeatures, particularly from ensGene, 
			# that have their own defined exons - essentially their own transcripts
			# examples include retained_intron, sense_overlapping, etc
			# we should keep these
			
			my $exon_check = 0;
			my $cds_check  = 0;
			foreach my $sf ($subfeat->get_SeqFeatures) {
				$exon_check++ if $sf->primary_tag =~ /^exon$/i;
				$cds_check++  if $sf->primary_tag =~ /^cds$/i;
			}
			next unless $exon_check or $cds_check;
			
			# initialize the gene table line structure
			my $ucsc = initialize_transcript($subfeat);
			
			# collect the names for the transcript and/or gene
			collect_names($ucsc, $gene, $subfeat);
			
			# process the transcript
			if ($cds_check) {
				process_transcript($ucsc, $subfeat);
			}
			else {
				process_nc_transcript($ucsc, $subfeat);
			}
			
			$transcript_success++;
		}
		
		elsif ($type eq 'exon') {
			# an exon, skip for now
			next;
		}
		
		elsif ($type eq 'codon') {
			# start_codon or stop_codon, skip for now
			next;
		}
		
		elsif ($type eq 'cds') {
			# there is a CDS feature in this gene
			# no RNA transcript defined
			# presumed gene transcript is gene -> CDS
			$cds_count++;
		}
		
		else {
			# catchall for unrecognized feature types
			# record and warn at end
			$unknowns{$type} += 1;
			next;
		}
	}
	
	# check that we have made a RNA transcript
	if (!$transcript_success and $cds_count) {
		# if not, make a presumed transcript from any CDSs
		
		# initialize the gene table line structure
		my $ucsc = initialize_transcript($gene);
		
		# collect the names for the transcript and/or gene
		collect_names($ucsc, $gene, $gene);
		
		# assume the gene is a transcript and process as such
		process_transcript($ucsc, $gene);
	}
}



### Collect gene names
sub collect_names {
	my ($ucsc, $gene, $transcript) = @_;
	
	# these names may or may not be same
	# we assume we will always have a primary_id
	$ucsc->{'name'}  = $gene->display_name || $gene->primary_id;
	$ucsc->{'name2'} = $transcript->display_name || $transcript->primary_id;
	
	# append additional names if requested and if possible
	if ($use_alias) {
		if ($gene->primary_id ne $ucsc->{'name'}) {
			$ucsc->{'name'} .= "|" . $gene->primary_id;
		}
		if ($gene->has_tag('Alias')) {
			foreach my $a ($gene->get_tag_values('Alias')) {
				next if $a eq $gene->display_name;
				next if $a eq $gene->primary_id;
				$ucsc->{'name'} .= "|$a";
			}
		}
	}
	
	# clean up name
	$ucsc->{'name'} =~ s/ ?\(\d+ of \d+\)//; # remove (1 of 2)
	$ucsc->{'name'} =~ s/ /_/;
	$ucsc->{'name2'} =~ s/ ?\(\d+ of \d+\)//; # remove (1 of 2)
	$ucsc->{'name2'} =~ s/ /_/;
}



### process a transcript
sub process_transcript {
	my ($ucsc, $transcript) = @_;
	
	# Record coordinates for this transcript
	$ucsc->{'txStart'} = $transcript->start - 1;
	$ucsc->{'txEnd'}   = $transcript->end;
	
	# Get the subfeatures of the transcript
	# these should be utrs, cds, exons, introns
	my @exons;
	my @cds;
	my @utr;
	foreach my $subf ($transcript->get_SeqFeatures) {
		
		# check the type
		my $type = lc $subf->primary_tag;
		
		# collect accordingly
		if ($type eq 'cds') {
			push @cds, $subf;
		}
		elsif ($type eq 'exon') {
			push @exons, $subf;
		}
		elsif ($type =~ /utr|untranslated/) {
			push @utr, $subf;
		}
		elsif ($type =~ /codon|intron/) {
			# ignore these
		}
		else {
			# catchall for unrecognized feature types
			# record and warn at end
			$unknowns{$type} += 1;
			next;
		}
	}
	
	# process the subfeatures
	if (@exons and @cds) {
		# we prefer to use the exon and cds information 
		return unless process_exon_cds($ucsc, \@exons, \@cds);
	}
	elsif (@cds) {
		# we have at least CDSs and maybe UTRs
		# we can calculate both exons and txStart and txEnd from these
		# adjacent cds and utr in the same exon will be appropriately merged
		merge_cds_utr_exons($ucsc, @cds, @utr);
	}
	else {
		warn " transcript " . $ucsc->{'name'} . " does not have exon, cds, " .
			"and/or utr subfeatures to process. Skipping\n";
	}
	
	# the transcript should be finished now
	# print the transcript
	print_table_item($ucsc);
}



### process a noncoding transcript
sub process_nc_transcript {
	my ($ucsc, $transcript) = @_;
	
	# Record coordinates for this transcript
	$ucsc->{'txStart'}  = $transcript->start - 1;
	$ucsc->{'txEnd'}    = $transcript->end;
	$ucsc->{'cdsStart'} = $transcript->end; # no cds, so set to end
	$ucsc->{'cdsEnd'}   = $transcript->end;
	
	# Get the subfeatures of the transcript
	# these should be non-coding exons
	my @subfeatures = $transcript->get_SeqFeatures;
	
	# process subfeatures, if there are any
	if (@subfeatures) {
		foreach my $subf (@subfeatures) {
			
			# look for exons, could be exon or noncoding_exon
			if ($subf->primary_tag =~ /exon/i) {
				# record the coordinates and information for this exon
				push @{ $ucsc->{'exonStarts'} }, $subf->start - 1;
				push @{ $ucsc->{'exonEnds'} }, $subf->end;
				$ucsc->{'exonCount'} += 1;
			}
			else {
				# catchall for unrecognized feature types
				# record and warn at end
				$unknowns{$subf->primary_tag} += 1;
				next;
			}
		}
	}
	else {
		# no exon subfeatures? I guess it's just one exon
		push @{ $ucsc->{'exonStarts'} }, $transcript->start - 1;
		push @{ $ucsc->{'exonEnds'} }, $transcript->end;
		$ucsc->{'exonCount'} += 1;
	}
	
	# the transcript should be finished now
	# print the transcript
	print_table_item($ucsc);
}






### process both exons and cds together
sub process_exon_cds {
	my ($ucsc, $exon, $cds) = @_;
	return if (scalar @$cds > @$exon);
	
	# sort both exons and cds by start position
	my @exons = map {$_->[0]} sort {$a->[1] <=> $b->[1]} map {[$_, $_->start]} @$exon;
	my @cdss  = map {$_->[0]} sort {$a->[1] <=> $b->[1]} map {[$_, $_->start]} @$cds;
	
	while (@exons) {
		my $e = shift @exons;
		
		# check that we still have CDSs left
		unless (@cdss) {
			# no more CDSs left, all utr now
			
			# check that we had defined the cdsEnd
			unless (defined $ucsc->{'cdsEnd'}) {
				# the cdsEnd would be at the end of the previous exon
				$ucsc->{'cdsEnd'} = $ucsc->{'exonEnds'}->[-1];
			}
			
			# add the utr
			push @{ $ucsc->{'exonStarts'} }, $e->start - 1;
			push @{ $ucsc->{'exonEnds'} }, $e->end;
			$ucsc->{'exonCount'} += 1;
			next;
		}
		
		# check if this exon overlaps a cds
		if ($e->start < $cdss[0]->start and $e->end < $cdss[0]->start) {
			# utr only
			push @{ $ucsc->{'exonStarts'} }, $e->start - 1;
			push @{ $ucsc->{'exonEnds'} }, $e->end;
			$ucsc->{'exonCount'} += 1;
		}
		elsif ($e->start < $cdss[0]->start and $e->end > $cdss[0]->end) {
			# exon contains the entire CDS with UTR on either side
			push @{ $ucsc->{'exonStarts'} }, $e->start - 1;
			push @{ $ucsc->{'exonEnds'} }, $e->end;
			$ucsc->{'exonCount'} += 1;
			$ucsc->{'cdsStart'} = $cdss[0]->start - 1;
			$ucsc->{'cdsEnd'} = $cdss[0]->end;
			shift @cdss; # finished with this cds
		}
		elsif ($e->start < $cdss[0]->start and $e->end == $cdss[0]->end) {
			# exon with cdsStart and both cds and utr
			push @{ $ucsc->{'exonStarts'} }, $e->start - 1;
			push @{ $ucsc->{'exonEnds'} }, $e->end;
			$ucsc->{'exonCount'} += 1;
			$ucsc->{'cdsStart'} = $cdss[0]->start - 1;
			
			# check txEnd
			if ($e->end == $ucsc->{'txEnd'}) {
				# in case there is no UTR at the end
				$ucsc->{'cdsEnd'} = $e->end;
			}
			shift @cdss; # finished with this cds
		}
		elsif ($e->start == $cdss[0]->start and $e->end == $cdss[0]->end) {
			# a CDS exon
			push @{ $ucsc->{'exonStarts'} }, $e->start - 1;
			push @{ $ucsc->{'exonEnds'} }, $e->end;
			$ucsc->{'exonCount'} += 1;
			
			# check txStart and txEnd
			if ($e->start - 1 == $ucsc->{'txStart'}) {
				# in case there is no UTR at the beginning
				$ucsc->{'cdsStart'} = $e->start - 1;
			}
			if (not defined $ucsc->{'cdsStart'}) {
				# cds starts with this exon
				$ucsc->{'cdsStart'} = $e->start - 1;
			}
			if ($e->end == $ucsc->{'txEnd'}) {
				# in case there is no UTR at the end
				$ucsc->{'cdsEnd'} = $e->end;
			}
			
			shift @cdss; # finished with this cds
		}
		elsif ($e->start == $cdss[0]->start and $e->end > $cdss[0]->end) {
			# exon with cdsEnd and both cds and utr
			push @{ $ucsc->{'exonStarts'} }, $e->start - 1;
			push @{ $ucsc->{'exonEnds'} }, $e->end;
			$ucsc->{'exonCount'} += 1;
			$ucsc->{'cdsEnd'} = $cdss[0]->end;
			
			# check txStart
			if ($e->start - 1 == $ucsc->{'txStart'}) {
				# in case there is no UTR at the beginning
				$ucsc->{'cdsStart'} = $e->start - 1;
			}
			if (not defined $ucsc->{'cdsStart'}) {
				# cds starts with this exon
				$ucsc->{'cdsStart'} = $e->start - 1;
			}
			shift @cdss; # finished with this cds
		}
		else {
			warn "Programmer error! compare transcript " . $ucsc->{'name'} . ", exon " . 
				$e->start . '..' . $e->end . " with next cds " . $cdss[0]->start . 
				'..' . $cdss[0]->end . "\n";
			return;
		}
	}
	
	return 1;
}



sub merge_cds_utr_exons {
	my $ucsc = shift;
	
	# sort the mix of CDS and any UTR by increasing start position
	my @exons = map {$_->[1]} sort {$a->[0] <=> $b->[0]} map { [$_->start, $_] } @_;
	
	while (@exons) {
		my $e = shift @exons;
		
		# let's work backwards on this, always shift new exon and compare to previous
		# better able to handle utr-cds-utr in one exon without so many conditionals
		
		# check for previous exons
		if ( $ucsc->{'exonCount'} ) {
			my $previous_end = $ucsc->{'exonEnds'}->[-1];
			
			#check for adjoining exons
			if ($e->start - 1 == $previous_end) {
				# adjoining
				$ucsc->{'exonEnds'}->[-1] = $e->end;
			}
			else {
				# not adjoining
				push @{$ucsc->{'exonStarts'}}, $e->start - 1;
				push @{$ucsc->{'exonEnds'}}, $e->end;
				$ucsc->{'exonCount'} += 1;
			}
			
			# check for cdsEnd
			if (
				defined $ucsc->{'cdsStart'} and
				not defined $ucsc->{'cdsEnd'} and
				$e->primary_tag !~ /cds/i
			) {
				# looks like the beginning of the utr
				# set cdsEnd to previous end
				$ucsc->{'cdsEnd'} = $previous_end;
			}
		}
		else {
			# first exon
			push @{$ucsc->{'exonStarts'}}, $e->start - 1;
			push @{$ucsc->{'exonEnds'}}, $e->end;
			$ucsc->{'exonCount'} += 1;
		}
		
		# check for cdsStart
		if ($e->primary_tag =~ /cds/i and not defined $ucsc->{'cdsStart'}) {
			# previous was utr, and now we have CDS
			$ucsc->{'cdsStart'} = $e->start - 1;
		}
	}
	
	unless (defined $ucsc->{'cdsEnd'}) {
		# set to the last exon end
		$ucsc->{'cdsEnd'} = $ucsc->{'exonEnds'}->[-1];
	}
}


### print the completed UCSC table line item
sub print_table_item {
	my $ucsc = shift;
	
	# assemble the required components
	my @components = (
		# name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds
			$ucsc->{'name'},
			$ucsc->{'name2'},
			$ucsc->{'chr'},
			$ucsc->{'strand'},
			$ucsc->{'txStart'},
			$ucsc->{'txEnd'},
			$ucsc->{'cdsStart'},
			$ucsc->{'cdsEnd'},
			$ucsc->{'exonCount'},
			join(",", @{$ucsc->{'exonStarts'}} ) . ',', 
			join(",", @{$ucsc->{'exonEnds'}} ) . ',',
	);
	
	$outfh->print(join("\t", @components) . "\n");
	$count++; # increment global counter
}



__END__

=head1 NAME

gff3_to_ucsc_table.pl

A script to convert a GFF3 file to a UCSC style refFlat table

=head1 SYNOPSIS

gff3_to_ucsc_table.pl [--options...] <filename>
  
  Options:
  --in <filename>
  --out <filename> 
  --alias
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the input GFF3 file. The file may be compressed with gzip.

=item --out <filename>

Specify the output filename. By default it uses the input file base 
name appended with '.refFlat'.

=item --alias

Specify that any additional aliases, including the primary_ID, should 
be appended to the gene name. They are concatenated using the pipe "|"
symbol.

=item --gz

Specify whether (or not) the output file should be compressed with gzip. 
The default is to mimic the status of the input file

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will convert a GFF3 annotation file to a UCSC-style 
gene table, using the refFlat format. This includes transcription 
and translation start and stops, as well as exon start and stops, 
but does not include coding exon frames. See the documentation at 
L<http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat> 
for more information.

The program assumes the input GFF3 file includes standard 
parent-E<gt>child relationships using primary IDs and primary tags, 
including gene, mRNA, exon, CDS, and UTRs. Non-standard genes, 
including non-coding RNAs, will also be processed too. Chromosomes, 
contigs, and embedded sequence are ignored. Non-pertinent features are 
safely ignored but reported. Most pragmas are ignored, except for close 
feature pragmas (###), which may aid in processing very large files. 
See the documentation for the GFF3 file format at 
L<http://www.sequenceontology.org/resources/gff3.html> for more 
information. 

Previous versions of this script attempted to export in the UCSC 
genePredExt table format, often with inaccurate results. Users 
who need this format should investigate the C<gff3ToGenePred> 
program available at L<http://hgdownload.cse.ucsc.edu/admin/exe/>.

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
