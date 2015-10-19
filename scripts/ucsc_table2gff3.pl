#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Net::FTP;
use Bio::SeqFeature::Lite;
use Bio::ToolBox::utility qw(
	format_with_commas
	open_to_read_fh
	open_to_write_fh
);
use Bio::ToolBox::parser::ucsc;
my $VERSION = '1.33';

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
	$ftp_file,
	$database,
	$host,
	$do_chromo,
	$refseqstatusf,
	$refseqsumf,
	$ensemblnamef,
	$ensemblsourcef,
	$kgxreff,
	$chromof,
	$user_source,
	$do_gene,
	$do_cds,
	$do_utr,
	$do_codon,
	$share,
	$do_name,
	$gz,
	$help,
	$print_version,
);
my @genetables;
GetOptions( 
	'ftp=s'      => \$ftp_file, # which database table to retrieve
	'db=s'       => \$database, # which ucsc genome to use
	'host=s'     => \$host, # the ftp server to connect to
	'chr!'       => \$do_chromo, # include the chromosome file from ftp
	'table=s'    => \@genetables, # the input gene table files
	'status=s'   => \$refseqstatusf, # the refseqstatus file
	'sum=s'      => \$refseqsumf, # the refseqsummary file
	'kgxref=s'   => \$kgxreff, # the kgXref info file
	'ensname=s'  => \$ensemblnamef, # the ensemblToGeneName file
	'enssrc=s'   => \$ensemblsourcef, # the ensemblSource file
	'chromo=s'   => \$chromof, # a chromosome file
	'source=s'   => \$user_source, # user provided source
	'gene!'      => \$do_gene, # include genes in output
	'cds!'       => \$do_cds, # include CDS in output
	'utr!'       => \$do_utr, # include UTRs in output
	'codon!'     => \$do_codon, # include start & stop codons in output
	'share!'     => \$share, # share common exons and UTRs
	'name!'      => \$do_name, # assign names to CDSs, UTRs, and exons
	'gz!'        => \$gz, # compress file
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
	print " Biotoolbox script ucsc_table2gff3.pl, version $VERSION\n\n";
	exit;
}





### Check requirements and defaults
unless (@genetables or $ftp_file or $chromof) {
	die " Specify either an input table file, chromosome file, or a FTP table!\n";
}
if ($ftp_file) {
	unless ($ftp_file =~ m/^refgene|ensgene|xenorefgene|known|all$/i) {
		die " requested table '$ftp_file' by FTP not supported! see help\n";
	}
	unless (defined $database) {
		die " a UCSC genome database must be provided! see help\n";
	}
	unless (defined $do_chromo) {
		$do_chromo = 1;
	}
	unless (defined $host) {
		$host = 'hgdownload.cse.ucsc.edu';
	}
}
unless (defined $do_gene) {
	$do_gene = 1;
}
unless (defined $do_utr) {
	$do_utr = 0;
}
unless (defined $do_cds) {
	$do_cds = 1;
	unless (defined $do_codon) {
		$do_codon = 0;
	}
}
unless (defined $share) {
	$share = 1;
}
unless (defined $do_name) {
	$do_name = 0;
}
my $start_time = time;




### Fetch files if requested
if ($ftp_file) {
	
	# collect the requested files by ftp
	my @files = fetch_files_by_ftp();
	
	# push file names into appropriate variables
	foreach my $file (@files) {
		if ($file =~ /refgene|ensgene|knowngene/i) {
			push @genetables, $file;
		}
		elsif ($file =~ /summary/i) {
			$refseqsumf = $file;
		}
		elsif ($file =~ /status/i) {
			$refseqstatusf = $file;
		}
		elsif ($file =~ /ensembltogene/i) {
			$ensemblnamef = $file;
		}
		elsif ($file =~ /ensemblsource/i) {
			$ensemblsourcef = $file;
		}
		elsif ($file =~ /kgxref/i) {
			$kgxreff = $file;
		}
		elsif ($file =~ /chrom/i) {
			$chromof = $file;
		}
	}
}




### Initiate the parser
my $ucsc = Bio::ToolBox::parser::ucsc->new(
	do_gene     => $do_gene,
	do_cds      => $do_cds,
	do_utr      => $do_utr,
	do_codon    => $do_codon,
	do_name     => $do_name,
	share       => $share,
) or die "cannot initialize ucsc parser!";

# add options
if ($user_source) {
	$ucsc->source($user_source);
}
if ($refseqsumf) {
	my $c = $ucsc->load_extra_data($refseqsumf, 'refseqsummary');
	printf " Loaded %s transcripts from supplemental data file '$refseqsumf'\n", 
		format_with_commas($c);
}
if ($refseqstatusf) {
	my $c = $ucsc->load_extra_data($refseqstatusf, 'refseqstatus');
	printf " Loaded %s transcripts from supplemental data file '$refseqstatusf'\n", 
		format_with_commas($c);
}
if ($ensemblnamef) {
	my $c = $ucsc->load_extra_data($ensemblnamef, 'ensembltogene');
	printf " Loaded %s transcripts from supplemental data file '$ensemblnamef'\n", 
		format_with_commas($c);
}
if ($ensemblsourcef) {
	my $c = $ucsc->load_extra_data($ensemblsourcef, 'ensemblsource');
	printf " Loaded %s transcripts from supplemental data file '$ensemblsourcef'\n", 
		format_with_commas($c);
}
if ($kgxreff) {
	my $c = $ucsc->load_extra_data($kgxreff, 'kgxref');
	printf " Loaded %s transcripts from supplemental data file '$kgxreff'\n", 
		format_with_commas($c);
}




# initialize globals
my $chromosome_done = 0; # boolean indicating chromosomes are written

# walk through the input tables
foreach my $file (@genetables) {
	
	# open output file
	my ($outfile, $gff_fh) = open_output_gff($file);
	
	# process chromosome
	if ($chromof and !$chromosome_done) {
		# if there is only one genetable, we will prepend the chromosomes 
		# to that output file, otherwise we'll make a separate gff file
		# I'm making this assumption because the chromosomes only need to be 
		# defined once when loading Bio::DB::SeqFeature::Store database
		# If user is collecting multiple gene tables, then separate files 
		# are ok, probably preferable, than a gigantic one
		print " Writing chromosome features....\n";
		
		if (scalar @genetables > 1) {
			# let's write a separate chromosome gff file
			
			# open new filehandle
			my ($chromo_outfile, $chromo_gff_fh) = open_output_gff($chromof);
			
			# convert the chromosomes
			print_chromosomes($chromo_gff_fh);
			
			# done
			$chromo_gff_fh->close;
			print " Wrote chromosome GFF file '$chromo_outfile'\n"; 
			$chromosome_done = 1;
		}
		else {
			# let's write to one gff file
			print_chromosomes($gff_fh);
			$chromosome_done = 1;
		}
	}	
	
	# convert the table
	print " Converting gene table '$file' features....\n";
	unless ($ucsc->open_file($file)) {
		warn " Unable to open file!\n";
		next;
	}
	my $tops = $ucsc->top_features;
	print_current_gene_list($gff_fh, $tops);
	
	
	# report outcomes
	my $count = $ucsc->counts;
	print "  converted ", format_with_commas($count->{gene}), 
		" gene features\n" if $count->{gene} > 0;
	print "  converted ", format_with_commas($count->{mrna}), 
		" mRNA transcripts\n" if $count->{mrna} > 0;
	print "  converted ", format_with_commas($count->{pseudogene}), 
		" pseudogene transcripts\n" if $count->{pseudogene} > 0;
	print "  converted ", format_with_commas($count->{ncrna}), 
		" ncRNA transcripts\n" if $count->{ncrna} > 0;
	print "  converted ", format_with_commas($count->{mirna}), 
		" miRNA transcripts\n" if $count->{mirna} > 0;
	print "  converted ", format_with_commas($count->{snrna}), 
		" snRNA transcripts\n" if $count->{snrna} > 0;
	print "  converted ", format_with_commas($count->{snorna}), 
		" snoRNA transcripts\n" if $count->{snorna} > 0;
	print "  converted ", format_with_commas($count->{trna}), 
		" tRNA transcripts\n" if $count->{trna} > 0;
	print "  converted ", format_with_commas($count->{rrna}), 
		" rRNA transcripts\n" if $count->{rrna} > 0;
	print "  converted ", format_with_commas($count->{other}), 
		" other transcripts\n" if $count->{other} > 0;
	
	# Finished
	printf "  wrote file '$outfile' in %.1f minutes\n", 
		(time - $start_time)/60;
	
}



### Finish
exit;





#########################  Subroutines  #######################################

sub fetch_files_by_ftp {
	
	
	# generate ftp request list
	my @ftp_files;
	if ($ftp_file eq 'all') {
		@ftp_files = qw(
			refgene
			ensgene
			xenorefgene
			known
		);
	}
	elsif ($ftp_file =~ /,/) {
		@ftp_files = split /,/, $ftp_file;
	}
	else {
		push @ftp_files, $ftp_file;
	}
	
	# generate list of files
	my @files;
	foreach my $item (@ftp_files) {
		if ($item =~ m/^xeno/i) {
			push @files, qw(
				xenoRefGene.txt.gz 
				refSeqStatus.txt.gz 
				refSeqSummary.txt.gz
			);
		}
		elsif ($item =~ m/refgene/i) {
			push @files, qw(
				refGene.txt.gz 
				refSeqStatus.txt.gz 
				refSeqSummary.txt.gz
			);
		}
		elsif ($item =~ m/ensgene/i) {
			push @files, qw(
				ensGene.txt.gz 
				ensemblToGeneName.txt.gz
				ensemblSource.txt.gz
			);
		}
		elsif ($item =~ m/known/i) {
			push @files, qw(
				knownGene.txt.gz 
				kgXref.txt.gz 
			);
		}
	}
	# this might seem convulated....
	# but we're putting all the file names in a single array
	# instead of specific global variables
	# to make retrieving through FTP a little easier
	# plus, not all files may be available for each species, e.g. knownGene
	# we also rename the files after downloading them
	
	# we will sort out the list of downloaded files later and assign them 
	# to specific global filename variables
	
	# add chromosome file if requested
	if ($do_chromo) {
		push @files, 'chromInfo.txt.gz';
	}
	
	# set the path based on user provided database
	my $path = 'goldenPath/' . $database . '/database/';
	
	# initiate connection
	print " Connecting to $host....\n";
	my $ftp = Net::FTP->new($host) or die "Cannot connect! $@";
	$ftp->login or die "Cannot login! " . $ftp->message;
	
	# prepare for download
	$ftp->cwd($path) or 
		die "Cannot change working directory to '$path'! " . $ftp->message;
	$ftp->binary;
	
	# download requested files
	my @fetched_files;
	foreach my $file (@files) {
		print "  fetching $file....\n";
		# prepend the local file name with the database
		my $new_file = $database . '_' . $file;
		
		# fetch
		if ($ftp->get($file, $new_file) ) { 
			push @fetched_files, $new_file;
		}
		else {	
			my $message = $ftp->message;
			if ($message =~ /no such file/i) {
				print "   file unavailable\n";
			}
			else {
				warn $message;
			}
		}
	}
	$ftp->quit;
	
	print " Finished\n";
	return @fetched_files;
}




sub open_output_gff {
	
	# prepare output file name
	my $file = shift;
	my $outfile = $file;
	$outfile =~ s/\.txt(?:\.gz)?$//i; # remove the extension
	$outfile .= '.gff3';
	if ($gz) {
		$outfile .= '.gz';
	}
	
	# open file handle
	my $fh = open_to_write_fh($outfile, $gz) or
		die " unable to open file '$outfile' for writing!\n";
	
	# print comments
	$fh->print( "##gff-version 3\n");
	$fh->print( "##genome-build UCSC $database\n") if $database;
	$fh->print( "# UCSC table file $file\n");
	
	# finish
	return ($outfile, $fh);
}


sub print_current_gene_list {
	my ($gff_fh, $top_features) = @_;
	
	# we need to sort the genes in genomic order before writing the GFF
	printf "  Sorting %s top features....\n", format_with_commas(scalar(@$top_features));
	my %pos2seqf;
	foreach my $gene (@$top_features) {
		# get coordinates
		my $start = $gene->start;
		my $chr;
		my $key;
		
		# identify which key to put under
		if ($gene->seq_id =~ /^chr(\d+)$/i) {
			$chr = $1;
			$key = 'numeric_chr';
		}
		elsif ($gene->seq_id =~ /^chr(\w+)$/i) {
			$chr = $1;
			$key = 'other_chr';
		}
		elsif ($gene->seq_id =~ /(\d+)$/) {
			$chr = $1;
			$key = 'other_numeric';
		}
		else {
			$chr = $gene->seq_id;
			$key = 'other';
		}
		
		# make sure start positions are unique, just in case
		# these modifications won't make it into seqfeature object
		while (exists $pos2seqf{$key}{$chr}{$start}) {
			$start++;
		}
		
		# store the seqfeature
		$pos2seqf{$key}{$chr}{$start} = $gene;
	}
	
	# print in genomic order
	# the gff_string method is undocumented in the POD, but is a 
	# valid method. Passing 1 should force a recursive action to 
	# print both parent and children.
	print "  Writing features to GFF....\n";
	foreach my $chr (sort {$a <=> $b} keys %{$pos2seqf{'numeric_chr'}} ) {
		foreach my $start (sort {$a <=> $b} keys %{ $pos2seqf{'numeric_chr'}{$chr} }) {
			# print the seqfeature recursively
			$gff_fh->print( $pos2seqf{'numeric_chr'}{$chr}{$start}->gff3_string(1));
			
			# print directive to close out all previous features
			$gff_fh->print("\n###\n"); 
		}
	}
	foreach my $chr (sort {$a cmp $b} keys %{$pos2seqf{'other_chr'}} ) {
		foreach my $start (sort {$a <=> $b} keys %{ $pos2seqf{'other_chr'}{$chr} }) {
			$gff_fh->print( $pos2seqf{'other_chr'}{$chr}{$start}->gff3_string(1));
			$gff_fh->print("\n###\n"); 
		}
	}
	foreach my $chr (sort {$a <=> $b} keys %{$pos2seqf{'other_numeric'}} ) {
		foreach my $start (sort {$a <=> $b} keys %{ $pos2seqf{'other_numeric'}{$chr} }) {
			$gff_fh->print( $pos2seqf{'other_numeric'}{$chr}{$start}->gff3_string(1));
			$gff_fh->print("\n###\n"); 
		}
	}
	foreach my $chr (sort {$a cmp $b} keys %{$pos2seqf{'other'}} ) {
		foreach my $start (sort {$a <=> $b} keys %{ $pos2seqf{'other'}{$chr} }) {
			$gff_fh->print( $pos2seqf{'other'}{$chr}{$start}->gff3_string(1));
			$gff_fh->print("\n###\n"); 
		}
	}
}



sub print_chromosomes {
	
	my $out_fh = shift;
	
	# open the chromosome file
	my $chromo_fh = open_to_read_fh($chromof) or die 
		"unable to open specified chromosome file '$chromof'!\n";
	
	# convert the chromosomes into GFF features
	# UCSC orders their chromosomes by chromosome length
	# I would prefer to order by numeric ID if possible
	my %chromosomes;
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
			-source        => 'UCSC', # using a generic source here
			-primary_tag   => $chr =~ m/^chr/i ? 'chromosome' : 'scaffold',
			-start         => 1,
			-end           => $end,
			-primary_id    => $chr,
			-display_name  => $chr,
		);
		
		# store the chromosome according to name
		if ($chr =~ /^chr(\d+)$/i) {
			$chromosomes{'numeric_chr'}{$1} = $chrom;
		}
		elsif ($chr =~ /^chr(\w+)$/i) {
			$chromosomes{'other_chr'}{$1} = $chrom;
		}
		elsif ($chr =~ /(\d+)$/) {
			$chromosomes{'other_numeric'}{$1} = $chrom;
		}
		else {
			$chromosomes{'other'}{$chr} = $chrom;
		}
	}
	$chromo_fh->close;
	
	# print the chromosomes
	foreach my $key (sort {$a <=> $b} keys %{ $chromosomes{'numeric_chr'} }) {
		# numeric chromosomes
		$chromosomes{'numeric_chr'}{$key}->version(3);
		$out_fh->print( $chromosomes{'numeric_chr'}{$key}->gff_string . "\n" );
	}
	foreach my $key (sort {$a cmp $b} keys %{ $chromosomes{'other_chr'} }) {
		# other chromosomes
		$chromosomes{'other_chr'}{$key}->version(3);
		$out_fh->print( $chromosomes{'other_chr'}{$key}->gff_string . "\n" );
	}
	foreach my $key (sort {$a <=> $b} keys %{ $chromosomes{'other_numeric'} }) {
		# numbered contigs, etc
		$chromosomes{'other_numeric'}{$key}->version(3);
		$out_fh->print( $chromosomes{'other_numeric'}{$key}->gff_string . "\n" );
	}
	foreach my $key (sort {$a cmp $b} keys %{ $chromosomes{'other'} }) {
		# contigs, etc
		$chromosomes{'other'}{$key}->version(3);
		$out_fh->print( $chromosomes{'other'}{$key}->gff_string . "\n" );
	}
	
	# finished
	$out_fh->print( "###\n" );
}




__END__

=head1 NAME 

ucsc_table2gff3.pl

A script to convert UCSC gene tables to GFF3 annotation.

=head1 SYNOPSIS

   ucsc_table2gff3.pl --ftp <text> --db <text>
   
   ucsc_table2gff3.pl [--options] --table <filename>
  
  Options:
  --ftp [refgene|ensgene|xenorefgene|known|all]
  --db <text>
  --host <text>
  --table <filename>
  --status <filename>
  --sum <filename>
  --ensname <filename>
  --enssrc <filename>
  --kgxref <filename>
  --chromo <filename>
  --source <text>
  --(no)chr             (true)
  --(no)gene            (true)
  --(no)cds             (true)
  --(no)utr             (false)
  --(no)codon           (false)
  --(no)share           (true)
  --(no)name            (false)
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --ftp [refgene|ensgene|xenorefgene|known|all]

Request that the current indicated tables and supporting files be 
downloaded from UCSC via FTP. Four different tables may be downloaded, 
including I<refGene>, I<ensGene>, I<xenoRefGene> mRNA gene prediction 
tables, and the UCSC I<knownGene> table (if available). Specify all to 
download all four tables. A comma delimited list may also be provided.

=item --db <text>

Specify the genome version database from which to download the requested 
table files. See L<http://genome.ucsc.edu/FAQ/FAQreleases.html> for a 
current list of available UCSC genomes. Examples included hg19, mm9, and 
danRer7.

=item --host <text>

Optionally provide the host FTP address for downloading the current 
gene table files. The default is 'hgdownload.cse.ucsc.edu'.

=item --table <filename>

Provide the name of a UCSC gene or gene prediction table. Tables known 
to work include the I<refGene>, I<ensGene>, I<xenoRefGene>, and UCSC 
I<knownGene> tables. Both simple and extended gene prediction tables, as 
well as refFlat tables are supported. The file may be gzipped. When 
converting multiple tables, use this option repeatedly for each table. 
The C<--ftp> option is recommended over using this one.

=item --status <filename>

Optionally provide the name of the I<refSeqStatus> table file. This file 
provides additional information for the I<refSeq>-based gene prediction 
tables, including I<refGene>, I<xenoRefGene>, and I<knownGene> tables. 
The file may be gzipped. The C<--ftp> option is recommended over using this.

=item --sum <filename>

Optionally provide the name of the I<refSeqSummary> file. This file 
provides additional information for the I<refSeq>-based gene prediction 
tables, including I<refGene>, I<xenoRefGene>, and I<knownGene> tables. The 
file may be gzipped. The C<--ftp> option is recommended over using this.

=item --ensname <filename>

Optionally provide the name of the I<ensemblToGeneName> file. This file 
provides a key to translate the Ensembl unique gene identifier to the 
common gene name. The file may be gzipped. The C<--ftp> option is 
recommended over using this.

=item --enssrc <filename>

Optionally provide the name of the I<ensemblSource> file. This file 
provides a key to translate the Ensembl unique gene identifier to the 
type of transcript, provided by Ensembl as the source. The file may be 
gzipped. The C<--ftp> option is recommended over using this.

=item --kgxref <filename>

Optionally provide the name of the I<kgXref> file. This file 
provides additional information for the UCSC I<knownGene> gene table.
The file may be gzipped.

=item --chromo <filename>

Optionally provide the name of the chromInfo text file. Chromosome 
and/or scaffold features will then be written at the beginning of the 
output GFF file (when processing a single table) or written as a 
separate file (when processing multiple tables). The file may be gzipped.

=item --source <text>

Optionally provide the text to be used as the GFF source. The default is 
automatically derived from the source table file name, if recognized, or 
'UCSC' if not recognized.

=item --(no)chr

When downloading the current gene tables from UCSC using the C<--ftp> 
option, indicate whether (or not) to include the I<chromInfo> table. 
The default is true. 

=item --(no)gene

Specify whether (or not) to assemble mRNA transcripts into genes. This 
will create the canonical gene-E<gt>mRNA-E<gt>(exon,CDS) heirarchical 
structure. Otherwise, mRNA transcripts are kept independent. The gene name, 
when available, are always associated with transcripts through the Alias 
tag. The default is true.

=item --(no)cds

Specify whether (or not) to include CDS features in the output GFF file. 
The default is true.

=item --(no)utr

Specify whether (or not) to include three_prime_utr and five_prime_utr 
features in the transcript heirarchy. If not defined, the GFF interpreter 
must infer the UTRs from the CDS and exon features. The default is false.

=item --(no)codon

Specify whether (or not) to include start_codon and stop_codon features 
in the transcript heirarchy. The default is false.

=item --(no)share

Specify whether exons, UTRs, and codons that are common between multiple 
transcripts of the same gene may be shared in the GFF3. Otherwise, each 
subfeature will be represented individually. This will reduce the size of 
the GFF3 file at the expense of increased complexity. If your parser 
cannot handle multiple parents, set this to --noshare. Due to the 
possibility of multiple translation start sites, CDS features are never 
shared. The default is true.

=item --(no)name

Specify whether you want subfeatures, including exons, CDSs, UTRs, and 
start and stop codons to have display names. In most cases, this 
information is not necessary. The default is false.

=item --gz

Specify whether the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will convert a UCSC gene or gene prediction table file into a
GFF3 format file. It will build canonical gene-E<gt>transcript-E<gt>[exon, 
CDS, UTR] heirarchical structures. It will attempt to identify non-coding genes
as to type using the gene name as inference. Various additional
informational attributes may also be included with the gene and transcript
features, which are derived from supporting table files.

Four table files are supported. Gene prediction tables, including I<refGene>, 
I<xenoRefGene>, and I<ensGene>, are supported. The UCSC I<knownGene> gene 
table, if available, is also supported. Supporting tables include I<refSeqStatus>, 
I<refSeqSummary>, I<ensemblToGeneName>, I<ensemblSource>, and I<kgXref>. 

Tables obtained from UCSC are typically in the extended GenePrediction 
format, although simple genePrediction and refFlat formats are also 
supported. See L<http://genome.ucsc.edu/FAQ/FAQformat.html#format9> regarding
UCSC gene prediction table formats. 

The latest table files may be automatically downloaded using FTP from 
UCSC or other host. Since these files are periodically updated, this may 
be the best option. Alternatively, individual files may be specified 
through command line options. Files may be obtained manually through FTP, 
HTTP, or the UCSC Table Browser. However, it is B<highly recommended> to 
let the program obtain the necessary files using the C<--ftp> option, as 
using the wrong file format or manipulating the tables may prevent the 
program from working properly.

If provided, chromosome and/or scaffold features may also be written to a 
GFF file. If only one table is being converted, then the chromosome features 
are prepended to the GFF file; otherwise, a separate chromosome GFF file is 
written.

If you need to set up a database using UCSC annotation, you should first 
take a look at the BioToolBox script B<db_setup.pl>, which provides a 
convenient automated database setup based on UCSC annotation. You can also 
find more information about loading a database in a How To document at 
L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. 

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
