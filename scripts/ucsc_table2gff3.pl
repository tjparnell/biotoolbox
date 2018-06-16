#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use Net::FTP;
use Bio::ToolBox::utility qw(
	format_with_commas
	open_to_read_fh
	open_to_write_fh
);
use Bio::ToolBox::parser::ucsc;
use Bio::ToolBox::GeneTools qw(gtf_string);
my $VERSION = '1.60';

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
	$do_gtf,
	$gz,
	$help,
	$print_version,
);
my @genetables;
GetOptions( 
	'f|ftp=s'      => \$ftp_file, # which database table to retrieve
	'd|db=s'       => \$database, # which ucsc genome to use
	'h|host=s'     => \$host, # the ftp server to connect to
	'chr!'         => \$do_chromo, # include the chromosome file from ftp
	't|table=s'    => \@genetables, # the input gene table files
	'a|status=s'   => \$refseqstatusf, # the refseqstatus file
	's|sum=s'      => \$refseqsumf, # the refseqsummary file
	'k|kgxref=s'   => \$kgxreff, # the kgXref info file
	'n|ensname=s'  => \$ensemblnamef, # the ensemblToGeneName file
	'r|enssrc=s'   => \$ensemblsourcef, # the ensemblSource file
	'c|chromo=s'   => \$chromof, # a chromosome file
	'source=s'     => \$user_source, # user provided source
	'gene!'        => \$do_gene, # include genes in output
	'cds!'         => \$do_cds, # include CDS in output
	'utr!'         => \$do_utr, # include UTRs in output
	'codon!'       => \$do_codon, # include start & stop codons in output
	'share!'       => \$share, # share common exons and UTRs
	'name!'        => \$do_name, # assign names to CDSs, UTRs, and exons
	'g|gtf!'       => \$do_gtf, # write a gtf file instead
	'z|gz!'        => \$gz, # compress file
	'h|help'       => \$help, # request help
	'v|version'    => \$print_version, # print the version
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
	print " Biotoolbox script ucsc_table2gff3.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
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
	do_exon     => 1, # always
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




# walk through the input tables
foreach my $file (@genetables) {
	
	# open output file
	my ($outfile, $gff_fh) = open_output_gff($file);
	
	# process chromosome
	if ($chromof) {
		# we will write sequence-region pragmas for every gff file automatically
		# the current gff parser and Bio::DB::SeqFeature will correctly parse them
		print_chromosomes($gff_fh);
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
	$outfile .= $do_gtf ? '.gtf' : '.gff3';
	if ($gz) {
		$outfile .= '.gz';
	}
	
	# open file handle
	my $fh = open_to_write_fh($outfile, $gz) or
		die " unable to open file '$outfile' for writing!\n";
	
	# print comments
	$fh->printf("##gff-version %s\n", $do_gtf ? '2.5' : '3');
	$fh->print("##genome-build UCSC $database\n") if $database;
	$fh->print("# UCSC table file $file\n");
	
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
			if ($do_gtf) {
				$gff_fh->print(gtf_string($pos2seqf{'numeric_chr'}{$chr}{$start}));
			}
			else {
				$gff_fh->print( $pos2seqf{'numeric_chr'}{$chr}{$start}->gff3_string(1));
				$gff_fh->print("###\n");
			} 
		}
	}
	foreach my $chr (sort {$a cmp $b} keys %{$pos2seqf{'other_chr'}} ) {
		foreach my $start (sort {$a <=> $b} keys %{ $pos2seqf{'other_chr'}{$chr} }) {
			if ($do_gtf) {
				$gff_fh->print(gtf_string($pos2seqf{'other_chr'}{$chr}{$start}));
			}
			else {
				$gff_fh->print( $pos2seqf{'other_chr'}{$chr}{$start}->gff3_string(1));
				$gff_fh->print("###\n");
			} 
		}
	}
	foreach my $chr (sort {$a <=> $b} keys %{$pos2seqf{'other_numeric'}} ) {
		foreach my $start (sort {$a <=> $b} keys %{ $pos2seqf{'other_numeric'}{$chr} }) {
			if ($do_gtf) {
				$gff_fh->print(gtf_string($pos2seqf{'other_numeric'}{$chr}{$start}));
			}
			else {
				$gff_fh->print( $pos2seqf{'other_numeric'}{$chr}{$start}->gff3_string(1));
				$gff_fh->print("###\n");
			} 
		}
	}
	foreach my $chr (sort {$a cmp $b} keys %{$pos2seqf{'other'}} ) {
		foreach my $start (sort {$a <=> $b} keys %{ $pos2seqf{'other'}{$chr} }) {
			if ($do_gtf) {
				$gff_fh->print(gtf_string($pos2seqf{'other'}{$chr}{$start}));
			}
			else {
				$gff_fh->print( $pos2seqf{'other'}{$chr}{$start}->gff3_string(1));
				$gff_fh->print("###\n");
			} 
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
		
		# store the chromosome according to name
		if ($chr =~ /^chr(\d+)$/i) {
			$chromosomes{'numeric_chr'}{$1} = "$chr 1 $end";
		}
		elsif ($chr =~ /^chr(\w+)$/i) {
			$chromosomes{'other_chr'}{$1} = "$chr 1 end";
		}
		elsif ($chr =~ /(\d+)$/) {
			$chromosomes{'other_numeric'}{$1} = "$chr 1 end";
		}
		else {
			$chromosomes{'other'}{$chr} = "$chr 1 end";
		}
	}
	$chromo_fh->close;
	
	# print the chromosomes
	foreach my $key (sort {$a <=> $b} keys %{ $chromosomes{'numeric_chr'} }) {
		# numeric chromosomes
		$out_fh->printf("##sequence-region  %s\n", $chromosomes{'numeric_chr'}{$key});
	}
	foreach my $key (sort {$a cmp $b} keys %{ $chromosomes{'other_chr'} }) {
		# other chromosomes
		$out_fh->printf("##sequence-region  %s\n", $chromosomes{'other_chr'}{$key});
	}
	foreach my $key (sort {$a <=> $b} keys %{ $chromosomes{'other_numeric'} }) {
		# numbered contigs, etc
		$out_fh->printf("##sequence-region  %s\n", $chromosomes{'other_numeric'}{$key});
	}
	foreach my $key (sort {$a cmp $b} keys %{ $chromosomes{'other'} }) {
		# contigs, etc
		$out_fh->printf("##sequence-region  %s\n", $chromosomes{'other'}{$key});
	}
}




__END__

=head1 NAME 

ucsc_table2gff3.pl

A program to convert UCSC gene tables to GFF3 or GTF annotation.

=head1 SYNOPSIS

   ucsc_table2gff3.pl --ftp <text> --db <text>
   
   ucsc_table2gff3.pl [--options] --table <filename>
  
  UCSC database options:
  -f --ftp [refgene|ensgene|            specify what tables to retrieve from UCSC
            xenorefgene|known|all]
  -d --db <text>                        UCSC database name: hg19,hg38,danRer7, etc
  -h --host <text>                      specify UCSC hostname
  
  Input file options:
  -t --table <filename>                 name of table, repeat or comma list
  -a --status <filename>                refSeqStatus file
  -s --sum <filename>                   refSeqSummary file
  -n --ensname <filename>               ensemblToGeneName file
  -r --enssrc <filename>                ensemblSource file
  -k --kgxref <filename>                kgXref file
  -c --chromo <filename>                chromosome file
  
  Conversion options:
  --source <text>                       source text, default UCSC
  --chr   | --nochr         (true)      include chromosomes in output
  --gene  | --nogene        (true)      assemble into genes
  --cds   | --nocds         (true)      include CDS subfeatures
  --utr   | --noutr         (false)     include UTR subfeatures
  --codon | --nocodon       (false)     include start and stop codons
  --share | --noshare       (true)      share subfeatures
  --name  | --noname        (false)     include name
  -g --gtf                              convert to GTF instead of GFF3
  
  General options:
  -z --gz                               compress output
  -v --version                          print version and exit
  -h --help                             show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 UCSC database options

=over 4

=item --ftp [refgene|ensgene|xenorefgene|known|all]

Request that the current indicated tables and supporting files be 
downloaded from UCSC via FTP. Four different tables may be downloaded, 
including I<refGene>, I<ensGene>, I<xenoRefGene> mRNA gene prediction 
tables, and the UCSC I<knownGene> table (if available). Specify all to 
download all four tables. A comma delimited list may also be provided.

=item --db E<lt>textE<gt>

Specify the genome version database from which to download the requested 
table files. See L<http://genome.ucsc.edu/FAQ/FAQreleases.html> for a 
current list of available UCSC genomes. Examples included hg19, mm9, and 
danRer7.

=item --host E<lt>textE<gt>

Optionally provide the host FTP address for downloading the current 
gene table files. The default is 'hgdownload.cse.ucsc.edu'.

=back

=head2 Input file options

=over 4

=item --table E<lt>filenameE<gt>

Provide the name of a UCSC gene or gene prediction table. Tables known 
to work include the I<refGene>, I<ensGene>, I<xenoRefGene>, and UCSC 
I<knownGene> tables. Both simple and extended gene prediction tables, as 
well as refFlat tables are supported. The file may be gzipped. When 
converting multiple tables, use this option repeatedly for each table. 
The C<--ftp> option is recommended over using this one.

=item --status E<lt>filenameE<gt>

Optionally provide the name of the I<refSeqStatus> table file. This file 
provides additional information for the I<refSeq>-based gene prediction 
tables, including I<refGene>, I<xenoRefGene>, and I<knownGene> tables. 
The file may be gzipped. The C<--ftp> option is recommended over using this.

=item --sum E<lt>filenameE<gt>

Optionally provide the name of the I<refSeqSummary> file. This file 
provides additional information for the I<refSeq>-based gene prediction 
tables, including I<refGene>, I<xenoRefGene>, and I<knownGene> tables. The 
file may be gzipped. The C<--ftp> option is recommended over using this.

=item --ensname E<lt>filenameE<gt>

Optionally provide the name of the I<ensemblToGeneName> file. This file 
provides a key to translate the Ensembl unique gene identifier to the 
common gene name. The file may be gzipped. The C<--ftp> option is 
recommended over using this.

=item --enssrc E<lt>filenameE<gt>

Optionally provide the name of the I<ensemblSource> file. This file 
provides a key to translate the Ensembl unique gene identifier to the 
type of transcript, provided by Ensembl as the source. The file may be 
gzipped. The C<--ftp> option is recommended over using this.

=item --kgxref E<lt>filenameE<gt>

Optionally provide the name of the I<kgXref> file. This file 
provides additional information for the UCSC I<knownGene> gene table.
The file may be gzipped.

=item --chromo E<lt>filenameE<gt>

Optionally provide the name of the chromInfo text file. Chromosome 
and/or scaffold features will then be written at the beginning of the 
output GFF file (when processing a single table) or written as a 
separate file (when processing multiple tables). The file may be gzipped.

=back

=head2 Conversion options

=over 4

=item --source E<lt>textE<gt>

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
shared. This will have no effect with GTF output. The default is true. 

=item --(no)name

Specify whether you want subfeatures, including exons, CDSs, UTRs, and 
start and stop codons to have display names. In most cases, this 
information is not necessary. This will have no effect with GTF output. 
The default is false.

=item --gtf

Specify that a GTF (version 2.5) format file should be written instead of 
GFF3. Yes, the name of the program says GFF3, but now we can output GTF 
too, and changing the name of the program is too late now.

=back

=head2 General options

=over 4

=item --gz

Specify whether the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will convert a UCSC gene or gene prediction table file into a
GFF3 (or optionally GTF) format file. It will build canonical 
gene-E<gt>transcript-E<gt>[exon, CDS, UTR] heirarchical structures. It will 
attempt to identify non-coding genesas to type using the gene name as inference. 
Various additional informational attributes may also be included with the gene 
and transcriptfeatures, which are derived from supporting table files.

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

If provided, chromosome and/or scaffold features will be written as GFF3-style 
sequence-region pragmas (even for GTF files, just in case).

If you need to set up a database using UCSC annotation, you should first 
take a look at the BioToolBox script L<db_setup.pl>, which provides a 
convenient automated database setup based on UCSC annotation.  

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
