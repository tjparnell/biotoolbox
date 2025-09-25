#!/usr/bin/perl

# documentation at end of file

use warnings;
use strict;
use English      qw(-no_match_vars);
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use IO::Prompt::Tiny qw(prompt);
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
);
use Bio::ToolBox::Parser;
use Bio::ToolBox::GeneTools qw(:all);
use Bio::ToolBox::utility   qw(parse_list format_with_commas);

our $VERSION = '2.03';

print "\n This program will get specific regions from features\n\n";

### Quick help
unless (@ARGV) {

	# when no command line options are present
	# print SYNOPSIS
	pod2usage(
		{
			'-verbose' => 0,
			'-exitval' => 1,
		}
	);
}

### Get command line options and initialize values
my (
	$infile,             $outfile,  $database, $request,
	$transcript_type,    $tsl,      $gencode,  $tbiotype,
	$start_adj,          $stop_adj, $unique,   $slop,
	$chromosome_exclude, $bed,      $gz,       $bed_name,
	$help,               $print_version,
);
my @feature_types;

# Command line options
GetOptions(
	'i|in=s'          => \$infile,              # the input data file
	'o|out=s'         => \$outfile,             # name of output file
	'd|db=s'          => \$database,            # source annotation database
	'f|feature=s'     => \@feature_types,       # the features requested from database
	'r|region=s'      => \$request,             # the region requested
	't|transcript=s'  => \$transcript_type,     # which transcripts to take
	'tsl=s'           => \$tsl,                 # filter on transcript support level
	'gencode!'        => \$gencode,             # filter on gencode basic tag
	'biotype=s'       => \$tbiotype,            # filter on transcript biotype
	'b|begin|start=i' => \$start_adj,           # start coordinate adjustment
	'e|end|stop=i'    => \$stop_adj,            # stop coordinate adjustment
	'u|unique!'       => \$unique,              # boolean to ensure uniqueness
	'l|slop=i'        => \$slop,                # slop factor in bp to identify uniqueness
	'K|chrskip=s'     => \$chromosome_exclude,  # skip chromosomes
	'bed!'            => \$bed,                 # convert the output to bed format
	'bedname=s'       => \$bed_name,            # what to use for bed name
	'z|gz!'           => \$gz,                  # compress output
	'h|help'          => \$help,                # request help
	'v|version'       => \$print_version,       # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {

	# print entire POD
	pod2usage(
		{
			'-verbose' => 2,
			'-exitval' => 1,
		}
	);
}

# Print version
if ($print_version) {
	print " Biotoolbox script get_gene_regions.pl, version $VERSION\n\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

### Check for requirements and set defaults
unless ( $infile or $database ) {
	print STDERR
" FATAL: must define a database or input gene table file! use --help for more information\n";
	exit 1;
}
if (    $database
	and $database =~
	/\. (?: gtf | gff3? | txt | refflat | genepred | ucsc ) (?:\.gz)? $/xi )
{

	# looks like a gene table file was specified as the database
	# intercept and assign to input file name
	$infile   = $database;
	$database = undef;
}

unless ( defined $outfile ) {
	print STDERR
		" FATAL: must define an output file name! use --help for more information\n";
	exit 1;
}

unless ( defined $slop ) {
	$slop = 0;
}

unless ( defined $gz ) {
	$gz = 0;
}

# one or more feature types may have been provided
# check if it is a comma delimited list
if ( scalar @feature_types == 1 and $feature_types[0] =~ /,/ ) {
	@feature_types = split /,/, shift @feature_types;
}

# boolean values for different transcript types to take
my (
	$do_mrna, $do_mirna, $do_ncrna,   $do_snrna,   $do_snorna,
	$do_trna, $do_rrna,  $do_miscrna, $do_lincrna, $do_all_rna
);

### Determine methods and transcript types
# Determine request_method and put into global variable
my $method;
determine_method();

### Collect feature regions
# collect the regions based on provided database or input file
# using the method determined above
printf " Collecting %s$request regions...\n", $unique ? 'unique ' : q();
my $outdata;
if ($database) {
	$outdata = collect_from_database();
}
elsif ($infile) {
	$outdata = collect_from_file();
}

printf " Collected %s regions\n", format_with_commas( $outdata->last_row );
print " Sorting...\n";
$outdata->gsort_data;

### Write the file
my $success;
if ($bed) {

	# output stream
	my $Stream = Bio::ToolBox::Data->new(
		stream => 1,
		out    => $outfile,
		bed    => 6,
	) or die "unable to write output file '$outfile'!";

	# determine what to use for bed name column
	my $name_i;
	if ($bed_name) {
		if ( $bed_name =~ /geneid/i ) {
			$name_i = 1;
		}
		elsif ( $bed_name =~ /genename/i ) {
			$name_i = 2;
		}
		elsif ( $bed_name =~ /transcriptid/i ) {
			$name_i = 3;
		}
		elsif ( $bed_name =~ /transcriptname/i ) {
			$name_i = 4;
		}
		elsif ( $bed_name =~ /featurename/i ) {
			$name_i = 5;
		}
		else {
			printf " ERROR: unrecognized bed name '%s'! Using default feature name\n";
			$name_i = 5;
		}
	}
	else {
		$name_i = 5;
	}

	# generate bed
	$outdata->iterate(
		sub {
			my $row = shift;
			my $bedstring =
				$row->bed_string( bed => 6, name => $row->value($name_i) );
			$Stream->write_row($bedstring);
		}
	);
	$success = $Stream->filename;
}
else {
	$success = $outdata->write_file(
		'filename' => $outfile,
		'gz'       => $gz,
	);
}
if ($success) {
	print " wrote file '$success'\n";
}
else {
	# failure! the subroutine will have printed error messages
	print " unable to write file!\n";
}

########################   Subroutines   ###################################

sub determine_method {

	# determine the region request from user if necessary
	unless ($request) {
		$request = collect_method_from_user();
	}
	$request = lc $request;

	# determine the method
	# also change the name of the request from short to long form
	if ( $request eq 'tss' ) {
		$request = 'transcription start site';
		$method  = \&collect_tss;
	}
	elsif ( $request eq 'transcription start site' ) {
		$method = \&collect_tss;
	}
	elsif ( $request eq 'tts' ) {
		$request = 'transcription stop site';
		$method  = \&collect_tts;
	}
	elsif ( $request eq 'transcription stop site' ) {
		$method = \&collect_tts;
	}
	elsif ( $request =~ /^splices?$/x ) {
		$request = 'splice sites';
		$method  = \&collect_splice_sites;
	}
	elsif ( $request =~ /^exons?$/x ) {
		$request = 'exon';
		$method  = \&collect_exons;
	}
	elsif ( $request =~ /^collapsed .? exon s?$/x ) {
		$request = 'collapsed exon';
		$method  = \&collect_collapsed_exons;
	}
	elsif ( $request =~ /^first.?exon$/x ) {
		$request = 'first exon';
		$method  = \&collect_first_exon;
	}
	elsif ( $request =~ /^last.?exon$/x ) {
		$request = 'last exon';
		$method  = \&collect_last_exon;
	}
	elsif ( $request =~ /^alt .* exons?$/x ) {
		$request = 'alternate exon';
		$method  = \&collect_alt_exons;
	}
	elsif ( $request =~ /^common .? exons?$/x ) {
		$request = 'common exon';
		$method  = \&collect_common_exons;
	}
	elsif ( $request =~ /^uncommon .? exons?$/x ) {
		$request = 'uncommon exon';
		$method  = \&collect_uncommon_exons;
	}
	elsif ( $request =~ /^introns?$/ ) {
		$method = \&collect_introns;
	}
	elsif ( $request =~ /^collapsed .? introns?$/x ) {
		$request = 'collapsed intron';
		$method  = \&collect_collapsed_introns;
	}
	elsif ( $request =~ /^first .? intron$/x ) {
		$request = 'first intron';
		$method  = \&collect_first_intron;
	}
	elsif ( $request =~ /^last .? intron$/x ) {
		$request = 'last intron';
		$method  = \&collect_last_intron;
	}
	elsif ( $request =~ /^alt .* introns?$/x ) {
		$request = 'alternate intron';
		$method  = \&collect_alt_introns;
	}
	elsif ( $request =~ /^common .? intron s?$/x ) {
		$request = 'common intron';
		$method  = \&collect_common_introns;
	}
	elsif ( $request =~ /^uncommon .? intron s?$/x ) {
		$request = 'uncommon intron';
		$method  = \&collect_uncommon_introns;
	}
	elsif ( $request eq 'utr' ) {
		$request = 'UTRs';
		$method  = \&collect_utrs;
	}
	elsif ( $request =~ /^cds .? start$/x ) {
		$request = 'CDS start';
		$method  = \&collect_cds_start;
	}
	elsif ( $request =~ /^cds .? stop$/x ) {
		$request = 'CDS stop';
		$method  = \&collect_cds_stop;
	}
	else {
		print STDERR " FATAL: unknown region request!\n";
		exit 1;
	}

	return $method;
}

sub collect_method_from_user {

	my %list = (
		1  => 'transcription start site',
		2  => 'transcription stop site',
		3  => 'exons',
		4  => 'collapsed exons',
		5  => 'first exon',
		6  => 'last exon',
		7  => 'alternate exons',
		8  => 'uncommon exons',
		9  => 'common exons',
		10 => 'introns',
		11 => 'collapsed introns',
		12 => 'first intron',
		13 => 'last intron',
		14 => 'alternate introns',
		15 => 'uncommon introns',
		16 => 'common introns',
		17 => 'splice sites',
		18 => 'UTRs',
		19 => 'CDS start',
		20 => 'CDS stop',
	);

	# request feature from the user
	print " These are the available feature types in the database:\n";
	foreach my $i ( sort { $a <=> $b } keys %list ) {
		print "   $i\t$list{$i}\n";
	}
	my $p      = " Enter the type of region to collect:  ";
	my $answer = prompt($p);

	# verify and return answer
	if ( exists $list{$answer} ) {
		return $list{$answer};
	}
	else {
		print STDERR " unknown request!\n";
		exit 1;
	}
}

sub determine_transcript_types {

	# if we are collecting from a database, the user may have already selected
	# an RNA feature type, which makes this selection redundant.
	my @features = @_;

	# collect all the transcript types requested
	my @types;
	if ($transcript_type) {

		# provided by the user from the command line
		if ( $transcript_type =~ /,/ ) {
			@types = split /,/, $transcript_type;
		}
		else {
			push @types, $transcript_type;
		}
	}
	elsif (@features) {

		# user selected types from a database
		foreach (@features) {
			my ( $p, $s ) = split /:/;    # take only the primary tag if both present
			push @types, $p if $p =~ /rna | transcript/xi;
		}
	}

	# request from the user if not otherwise provided
	unless (@types) {

		print " Genes may generate different types of RNA transcripts.\n";
		my $i = 1;
		my %i2tag;
		foreach (qw(all mRNA ncRNA snRNA snoRNA tRNA rRNA miRNA lincRNA misc_RNA)) {
			print "   $i\t$_\n";
			$i2tag{$i} = $_;
			$i++;
		}
		my $p        = ' Select one or more RNA types to include:  ';
		my $response = prompt($p);
		@types = map { $i2tag{$_} || undef } parse_list($response);
	}

	# string for visual output
	my $string = " Collecting transcript types:";

	foreach (@types) {
		if (m/^all$/i) {
			$do_all_rna = 1;
			$string .= ' all RNAs';
			last;
		}
		if (m/^mRNA$/i) {
			$do_mrna = 1;
			$string .= ' mRNA';
		}
		if (m/^miRNA$/i) {
			$do_mirna = 1;
			$string .= ' miRNA';
		}
		if (m/^ncRNA$/i) {
			$do_ncrna = 1;
			$string .= ' ncRNA';
		}
		if (m/^snRNA$/i) {
			$do_snrna = 1;
			$string .= ' snRNA';
		}
		if (m/^snoRNA$/i) {
			$do_snorna = 1;
			$string .= ' snoRNA';
		}
		if (m/^tRNA$/i) {
			$do_trna = 1;
			$string .= ' tRNA';
		}
		if (m/^rRNA$/i) {
			$do_rrna = 1;
			$string .= ' rRNA';
		}
		if (m/^misc_RNA$/i) {
			$do_miscrna = 1;
			$string .= ' misc_RNA';
		}
		if (m/^lincRNA$/i) {
			$do_lincrna = 1;
			$string .= ' lincRNA';
		}
	}
	print "$string\n";
}

sub collect_from_database {

	# open database connection
	my $db = open_db_connection($database)
		or die " unable to open database connection!\n";

	# get feature type if necessary
	my $prompt = <<PROMPT;
 Select one or more database features (typically genes) from which to collect regions. 
PROMPT
	my @features = verify_or_request_feature_types(
		'db'      => $db,
		'feature' => \@feature_types,
		'prompt'  => $prompt,
		'single'  => 0,
		'limit'   => 'gene|rna',
	) or die "No valid gene feature type was provided! see help\n";

	# get transcript_type
	determine_transcript_types(@features);

	# generate output data
	my $Data = generate_output_structure();
	$Data->database($database);

	# generate a seqfeature stream
	my $iterator = $db->features(
		-type     => \@features,
		-iterator => 1,
	);

	# process the features
	while ( my $seqfeat = $iterator->next_seq ) {

		# skip unwanted chromosomes
		if ( $chromosome_exclude and $seqfeat->seq_id =~ $chromosome_exclude ) {
			next;
		}

		# collect the regions based on the primary tag and the method re
		if ( $seqfeat->primary_tag eq 'gene' ) {

			# gene
			my @regions = process_gene($seqfeat);
			foreach (@regions) {

				# each element is an anon array of found feature info
				$Data->add_row($_);
			}
		}
		elsif ( $seqfeat->primary_tag =~ /rna | transcript/xi ) {

			# transcript
			my @regions = process_transcript($seqfeat);

			# remove duplicates if requested
			if ($unique) {
				remove_duplicates( \@regions );
			}

			# add regions to collection
			foreach (@regions) {
				$Data->add_row($_);
			}
		}
	}

	# finished
	return $Data;
}

sub collect_from_file {

	# get transcript_type
	unless ( defined $transcript_type ) {

		# user providing a file, so we'll just take everything in here
		$transcript_type = 'all';
	}
	determine_transcript_types();

	# Collect the top features for each sequence group.
	# Rather than going feature by feature through the gff,
	# we'll load the top features, which will collect all the features
	# and assemble appropriate feature -> subfeatures according to the
	# parent - child attributes.

	# generate output data
	my $Data = generate_output_structure();
	$Data->add_comment("Source data file $infile");

	# open appropriate parser object
	my $Parser = Bio::ToolBox::Parser->new( file => $infile, );
	unless ($Parser) {
		die " unable to open and parse '$infile'!\n";
	}

	# set parser attributes - most of these are probably redundant or not necessary
	my $type_string = $Parser->typelist;
	$Parser->do_gene(1);    # assume we always want genes?
	$Parser->do_exon(1);
	$Parser->do_cds(1);
	if ( $request =~ /utr/i ) {
		if ( $type_string =~ /utr/i ) {
			$Parser->do_utr(1);
		}
	}
	elsif ( $request =~ /cds st/i ) {
		if ( $type_string =~ /codon/ ) {
			$Parser->do_codon(1);
		}
	}
	if ( $tsl or $type_string !~ /rna/i ) {

		# we need the extra attributes for transcript_support_level and biotype
		$Parser->simplify(0);
	}
	else {
		$Parser->simplify(1);
	}

	# process the features
	my @bad_features;
	while ( my $seqfeat = $Parser->next_top_feature ) {

		# skip unwanted chromosomes
		if ( $chromosome_exclude and $seqfeat->seq_id =~ $chromosome_exclude ) {
			next;
		}

		# collect the regions based on the primary tag
		my $type = $seqfeat->primary_tag;
		if ( $type =~ /gene$/i ) {

			# gene, including things like gene, miRNA_gene, pseudogene, etc
			my @regions = process_gene($seqfeat);
			foreach (@regions) {
				$Data->add_row($_);
			}
		}
		elsif ( $type =~ /rna | transcript/xi ) {

			# any sort of RNA transcript
			my @regions = process_transcript($seqfeat);

			# remove duplicates if requested
			if ($unique) {
				remove_duplicates( \@regions );
			}
			foreach (@regions) {
				$Data->add_row($_);
			}
		}
		elsif ( $type =~ /chromosome | contig | scaffold | sequence/xi ) {
			next;    # safely ignore these
		}
		else {
			push @bad_features, $seqfeat;
		}
	}

	# finished
	if (@bad_features) {
		my %bad_types;
		foreach (@bad_features) {
			$bad_types{ $_->primary_tag } += 1;
		}
		printf " skipped %d unrecognized top feature types:\n%s\n",
			scalar(@bad_features),
			join( "\n", map {"  $bad_types{$_} $_"} sort { $a cmp $b } keys %bad_types );
	}
	return $Data;
}

sub generate_output_structure {
	my $Data = Bio::ToolBox::Data->new(
		feature => "region",
		columns => [
			qw(GeneID GeneName TranscriptID TranscriptName Name Chromosome Start
				Stop Strand)
		],
	);
	$Data->program("$PROGRAM_NAME, v $VERSION");
	my $r = $request;
	$r =~ s/\s/_/g;    # remove spaces
	$Data->metadata( 3, 'type', $transcript_type );
	$Data->metadata( 5, 'type', $r );
	if ($start_adj) {
		$Data->metadata( 7, 'start_adjusted', $start_adj );
	}
	if ($stop_adj) {
		$Data->metadata( 8, 'stop_adjusted', $stop_adj );
	}
	if ($unique) {
		$Data->metadata( 5, 'unique', 1 );
		$Data->metadata( 5, 'slop',   $slop );
	}

	return $Data;
}

sub process_gene {

	# passed gene object
	my $gene = shift;
	my @regions;

	# need to pull out the appropriate transcript types from the gene
	my @transcripts;
	foreach my $t ( get_transcripts($gene) ) {
		next unless acceptable_transcript($t);
		push @transcripts, $t;
	}
	return unless @transcripts;

	# filter for gencode transcripts if requested
	if ($gencode) {
		my $new_transcripts = filter_transcript_gencode_basic( \@transcripts );
		return unless scalar @{$new_transcripts};
		@transcripts = @{$new_transcripts};
	}

	# filter for transcript support level if requested
	if ($tsl) {
		my $new_transcripts = filter_transcript_support_level( \@transcripts, $tsl );
		return unless scalar @{$new_transcripts};
		@transcripts = @{$new_transcripts};
	}

	# filter for biotype if requested
	if ($tbiotype) {
		my $new_transcripts = filter_transcript_biotype( \@transcripts, $tbiotype );
		return unless scalar @{$new_transcripts};
		@transcripts = @{$new_transcripts};
	}

	# alternate or common exons require working with multiple transcripts
	if ( $request =~ /alternate | common | collapsed/xi ) {

		# pass all the transcripts together
		@regions = &{$method}(@transcripts);
	}
	else {
		# do each transcript one at a time
		foreach my $t (@transcripts) {
			push @regions, &{$method}($t);
		}
	}
	return unless @regions;

	# finally add gene name
	foreach my $region (@regions) {
		unshift @{$region}, $gene->primary_id, $gene->display_name;
	}

	# remove duplicates if requested
	if ($unique) {
		remove_duplicates( \@regions );
	}

	# return the regions
	return @regions;
}

sub process_transcript {

	# passed objects
	my $transcript = shift;

	# filter transcripts
	return unless acceptable_transcript($transcript);

	# filter for gencode transcripts if requested
	if ($gencode) {
		$transcript = filter_transcript_gencode_basic($transcript);
		return unless scalar @{$transcript};
	}

	# filter for transcript support level if requested
	if ($tsl) {
		$transcript = filter_transcript_support_level( $transcript, $tsl );
		return unless scalar @{$transcript};
	}

	# filter for biotype if requested
	if ($tbiotype) {
		$transcript = filter_transcript_biotype( $transcript, $tbiotype );
		return unless scalar @{$transcript};
	}

	# call appropriate method
	my @regions = &{$method}($transcript);

	# add non-existent gene id and name
	foreach my $region (@regions) {
		unshift @{$region}, '.', '.';
	}
	return @regions;
}

sub collect_tss {

	# get seqfeature objects
	my $transcript = shift;

	# get coordinates
	my $chromo = $transcript->seq_id;
	my ( $start, $stop, $strand );
	if ( $transcript->strand == 1 ) {

		# forward strand

		$strand = 1;
		$start  = $transcript->start;
		$stop   = $transcript->start;
	}
	elsif ( $transcript->strand == -1 ) {

		# reverse strand

		$strand = -1;
		$start  = $transcript->end;
		$stop   = $transcript->end;
	}
	else {
		die " poorly formatted transcript seqfeature object with strand 0!\n";
	}

	# get name
	my $name = $transcript->display_name . '_TSS';

	return _adjust_positions(
		[
			$transcript->primary_id, $transcript->display_name, $name, $chromo,
			$start, $stop, $strand
		]
	);
}

sub collect_tts {

	# get seqfeature objects
	my $transcript = shift;

	# get coordinates
	my $chromo = $transcript->seq_id;
	my ( $start, $stop, $strand );
	if ( $transcript->strand == 1 ) {

		# forward strand

		$strand = 1;
		$start  = $transcript->end;
		$stop   = $transcript->end;
	}
	elsif ( $transcript->strand == -1 ) {

		# reverse strand

		$strand = -1;
		$start  = $transcript->start;
		$stop   = $transcript->start;
	}
	else {
		die " poorly formatted transcript seqfeature object with strand 0!\n";
	}

	# get name
	my $name = $transcript->display_name . '_TTS';

	return _adjust_positions(
		[
			$transcript->primary_id, $transcript->display_name, $name, $chromo,
			$start, $stop, $strand
		]
	);
}

sub collect_exons {
	my $transcript = shift;

	# process and adjust the exons
	my @exons;
	foreach my $e ( get_exons($transcript) ) {
		push @exons,
			_adjust_positions(
				[
					$transcript->primary_id, $transcript->display_name,
					$e->display_name,        $e->seq_id,
					$e->start,               $e->end,
					$e->strand
				]
			);
	}

	return @exons;
}

sub collect_collapsed_exons {
	my $transcript = collapse_transcripts(@_);
	return collect_exons($transcript);
}

sub collect_first_exon {
	my $transcript = shift;

	# find the exons and/or CDSs
	my @list = get_exons($transcript);
	return unless @list;

	# the first exon
	my $first = $transcript->strand >= 0 ? shift @list : pop @list;

	# identify the exon name if it has one
	my $name = $first->display_name
		|| $transcript->display_name . "_firstExon";

	# finished
	return _adjust_positions(
		[
			$transcript->primary_id, $transcript->display_name,
			$name,                   $first->seq_id,
			$first->start,           $first->end,
			$first->strand
		]
	);
}

sub collect_last_exon {
	my $transcript = shift;

	# find the exons and/or CDSs
	my @list = get_exons($transcript);
	return unless @list;

	# the last exon
	my $last = $transcript->strand >= 0 ? pop @list : shift @list;

	# identify the exon name if it has one
	my $name = $last->display_name
		|| $transcript->display_name . "_lastExon";

	# finished
	return _adjust_positions(
		[
			$transcript->primary_id, $transcript->display_name,
			$name,                   $last->seq_id,
			$last->start,            $last->end,
			$last->strand
		]
	);
}

sub collect_alt_exons {
	my $ac_exons = get_alt_common_exons(@_);

	# we need the transcript name, so can't use the simpler get_alt_exons()
	my @exons;
	foreach my $transcript ( keys %{$ac_exons} ) {
		next if $transcript eq 'common';
		next if $transcript eq 'uncommon';
		foreach my $e ( @{ $ac_exons->{$transcript} } ) {
			push @exons,
				_adjust_positions(
					[
						$transcript->primary_id, $transcript, $e->display_name,
						$e->seq_id,              $e->start,   $e->end,
						$e->strand
					]
				);
		}
	}
	return @exons;
}

sub collect_uncommon_exons {
	my @exons;
	foreach my $e ( get_uncommon_exons(@_) ) {
		push @exons, _adjust_positions(
			[
				'uncommon',    # more than one transcript, so put generic identifier
				'uncommon',
				$e->display_name,
				$e->seq_id,
				$e->start,
				$e->end,
				$e->strand,
			]
		);
	}
	return @exons;
}

sub collect_common_exons {
	my @exons;
	foreach my $e ( get_common_exons(@_) ) {
		push @exons, _adjust_positions(
			[
				'common',    # more than one transcript, so put generic identifier
				'common',
				$e->display_name,
				$e->seq_id,
				$e->start,
				$e->end,
				$e->strand,
			]
		);
	}
	return @exons;
}

sub collect_splice_sites {

	# seqfeature object
	my $transcript = shift;

	# find the exons and/or CDSs
	my $list = get_exons($transcript);
	return unless $list;
	return if ( scalar @{$list} == 1 );

	# identify the last exon index position
	my $last = scalar @{$list} - 1;

	# collect the splice sites
	my @splices;

	# forward strand
	if ( $transcript->strand == 1 ) {

		# walk through each exon
		for ( my $i = 0; $i <= $last; $i++ ) {

			# get the exon name
			my $exon = $list->[$i];
			my $name = $exon->display_name
				|| $transcript->display_name . ".exon$i";

			# first exon
			if ( $i == 0 ) {
				push @splices,
					_adjust_positions(
						[
							$transcript->primary_id, $transcript->display_name,
							$name . '_3p',           $exon->seq_id,
							$exon->end + 1,          $exon->end + 1
						]
					);
			}

			# last exon
			elsif ( $i == $last ) {
				push @splices,
					_adjust_positions(
						[
							$transcript->primary_id, $transcript->display_name,
							$name . '_5p',           $exon->seq_id,
							$exon->start - 1,        $exon->start - 1,
							$exon->strand
						]
					);

			}

			# middle exons
			else {

				# 5' splice
				push @splices,
					_adjust_positions(
						[
							$transcript->primary_id, $transcript->display_name,
							$name . '_5p',           $exon->seq_id,
							$exon->start - 1,        $exon->start - 1,
							$exon->strand
						]
					);

				# 3' splice
				push @splices,
					_adjust_positions(
						[
							$transcript->primary_id, $transcript->display_name,
							$name . '_3p',           $exon->seq_id,
							$exon->end + 1,          $exon->end + 1,
							$exon->strand
						]
					);
			}
		}
	}

	# reverse strand
	else {

		# walk through each exon
		for ( my $i = 0; $i <= $last; $i++ ) {

			# get the exon name
			my $exon = $list->[$i];
			my $name = $exon->display_name
				|| $transcript->display_name . ".exon$i";

			# first exon
			if ( $i == 0 ) {
				push @splices,
					_adjust_positions(
						[
							$transcript->primary_id, $transcript->display_name,
							$name . '_3p',           $exon->seq_id,
							$exon->start - 1,        $exon->start - 1,
							$exon->strand
						]
					);
			}

			# last exon
			elsif ( $i == $last ) {
				push @splices,
					_adjust_positions(
						[
							$transcript->primary_id, $transcript->display_name,
							$name . '_5p',           $exon->seq_id,
							$exon->end + 1,          $exon->end + 1,
							$exon->strand
						]
					);

			}

			# middle exons
			else {

				# 5' splice
				push @splices,
					_adjust_positions(
						[
							$transcript->primary_id, $transcript->display_name,
							$name . '_5p',           $exon->seq_id,
							$exon->end + 1,          $exon->end + 1,
							$exon->strand
						]
					);

				# 3' splice
				push @splices,
					_adjust_positions(
						[
							$transcript->primary_id, $transcript->display_name,
							$name . '_3p',           $exon->seq_id,
							$exon->start - 1,        $exon->start - 1,
							$exon->strand
						]
					);
			}
		}
	}

	# finished
	return @splices;
}

sub collect_introns {
	my $transcript = shift;

	# collect the introns
	my @introns;
	foreach my $int ( get_introns($transcript) ) {
		push @introns,
			_adjust_positions(
				[
					$transcript->primary_id, $transcript->display_name,
					$int->display_name,      $int->seq_id,
					$int->start,             $int->end,
					$int->strand
				]
			);
	}

	# finished
	return @introns;
}

sub collect_collapsed_introns {
	my $transcript = collapse_transcripts(@_);
	return collect_introns($transcript);
}

sub collect_first_intron {
	my $transcript = shift;

	# find the introns
	my @list = get_introns($transcript);
	return unless @list;

	# the first intron
	my $first = $transcript->strand >= 0 ? shift @list : pop @list;
	return _adjust_positions(
		[
			$transcript->primary_id, $transcript->display_name,
			$first->display_name,    $first->seq_id,
			$first->start,           $first->end,
			$first->strand
		]
	);
}

sub collect_last_intron {
	my $transcript = shift;

	# get the introns
	my @list = get_introns($transcript);
	return unless @list;

	# the last intron
	my $last = $transcript->strand >= 0 ? pop @list : shift @list;
	return _adjust_positions(
		[
			$transcript->primary_id, $transcript->display_name,
			$last->display_name,     $last->seq_id,
			$last->start,            $last->end,
			$last->strand
		]
	);
}

sub collect_alt_introns {
	my $ac_introns = get_alt_common_introns(@_);
	my @introns;
	foreach my $transcript ( keys %{$ac_introns} ) {
		next if $transcript eq 'common';
		next if $transcript eq 'uncommon';
		foreach my $i ( @{ $ac_introns->{$transcript} } ) {
			push @introns,
				_adjust_positions(
					[
						$transcript, $transcript, $i->display_name, $i->seq_id,
						$i->start,   $i->end,     $i->strand
					]
				);
		}
	}
	return @introns;
}

sub collect_uncommon_introns {
	my @introns;
	foreach my $i ( get_uncommon_introns(@_) ) {
		push @introns, _adjust_positions(
			[
				'uncommon',    # more than one transcript, so put generic identifier
				'uncommon',
				$i->display_name,
				$i->seq_id,
				$i->start,
				$i->end,
				$i->strand,
			]
		);
	}
	return @introns;
}

sub collect_common_introns {
	my @introns;
	foreach my $i ( get_common_introns(@_) ) {
		push @introns, _adjust_positions(
			[
				'common',    # more than one transcript, so put generic identifier
				'common',
				$i->display_name,
				$i->seq_id,
				$i->start,
				$i->end,
				$i->strand,
			]
		);
	}
	return @introns;
}

sub collect_utrs {
	my $transcript = shift;

	# process and adjust the UTRs
	my @utrs;
	foreach my $u ( get_utrs($transcript) ) {
		push @utrs,
			_adjust_positions(
				[
					$transcript->primary_id, $transcript->display_name,
					$u->display_name,        $u->seq_id,
					$u->start,               $u->end,
					$u->strand
				]
			);
	}

	return @utrs;
}

sub collect_cds_start {

	# get seqfeature objects
	my $transcript = shift;

	# get the cds start
	my $pos = get_cdsStart($transcript);
	return unless $pos;

	# get other things
	my $chromo = $transcript->seq_id;
	my $strand = $transcript->strand;
	my $name   = $transcript->display_name . '_cdsStart';

	return _adjust_positions(
		[
			$transcript->primary_id, $transcript->display_name, $name, $chromo, $pos,
			$pos, $strand
		]
	);
}

sub collect_cds_stop {

	# get seqfeature objects
	my $transcript = shift;

	# get the cds start
	my $pos = get_cdsEnd($transcript);
	return unless $pos;

	# get other things
	my $chromo = $transcript->seq_id;
	my $strand = $transcript->strand;
	my $name   = $transcript->display_name . '_cdsStart';

	return _adjust_positions(
		[
			$transcript->primary_id, $transcript->display_name, $name, $chromo, $pos,
			$pos, $strand
		]
	);
}

sub acceptable_transcript {
	my $transcript = shift;
	my $t          = $transcript->primary_tag;
	return 1
		if (    $t =~ m/rna | transcript | retained_intron | antisense | nonsense/xi
			and $do_all_rna );
	return 1 if ( is_coding($transcript) and $do_mrna );
	return 1 if ( $t =~ /mirna/i  and $do_mirna );
	return 1 if ( $t =~ /ncrna/i  and $do_ncrna );
	return 1 if ( $t =~ /snrna/i  and $do_snrna );
	return 1 if ( $t =~ /snorna/i and $do_snorna );
	return 1 if ( $t =~ /trna/i   and $do_rrna );
	return 1 if ( $t =~ /rrna/i   and $do_rrna );
	return 1
		if (    $t =~ m/misc_rna | transcript | retained_intron | antisense | nonsense/xi
			and $do_miscrna );
	return 1 if ( $t =~ /lincrna/i and $do_lincrna );
	return 0;
}

sub _adjust_positions {

	my $region = shift;

	# region is an anonymous array of 7 elements
	# [$transcript_id, $transcript_name, $name, $chromo, $start, $stop, $strand]

	# adjust the start and end positions according to strand
	if ( $region->[6] == 1 ) {

		# forward strand

		if ($start_adj) {
			$region->[4] += $start_adj;
		}
		if ($stop_adj) {
			$region->[5] += $stop_adj;
		}
	}
	elsif ( $region->[6] == -1 ) {

		# reverse strand

		if ($start_adj) {
			$region->[5] -= $start_adj;
		}

		# stop
		if ($stop_adj) {
			$region->[4] -= $stop_adj;
		}
	}

	# return adjusted region coordinates
	return $region;
}

sub remove_duplicates {

	my $regions = shift;

	# regions are an anonymous arrays of 9 elements
	# start/stop are indexes 6,7

	# look for duplicates using a quick hash of seen positions
	my %seenit;
	my @to_remove;

	if ($slop) {

		# to simplify, we only use the start position when using slop, since
		# the end position is going to be too variable and flapping in the
		# breeze, so to speak
		# not entirely accurate, but that's slop for you
		for my $i ( 0 .. $#{$regions} ) {

			# we will be using the start position as the unique identifier
			# to account for the slop factor, we'll be adding/subtracting the
			# slop value to/from the start position
			# if this position matches anything else, we'll assume it's a duplicate

			foreach my $pos (

				# generate an array of possible start positions
				# with a default slop of 0, this will only be 1 position
				( $regions->[$i]->[6] - $slop ) .. ( $regions->[$i]->[6] + $slop )
				)
			{
				if ( exists $seenit{$pos} ) {
					push @to_remove, $i;
				}
				else {
					$seenit{$pos} = 1;
				}
			}
		}
	}
	else {
		# without slop, we look for precise matches based on both start and end
		for my $i ( 0 .. $#{$regions} ) {
			my ( $s, $e ) = ( $regions->[$i]->[6], $regions->[$i]->[7] );
			if ( exists $seenit{$s}{$e} ) {
				push @to_remove, $i;
			}
			else {
				$seenit{$s}{$e} = 1;
			}
		}
	}

	# remove the duplicates
	while (@to_remove) {
		my $i = pop @to_remove;

		# take from end to avoid shifting regions array
		splice( @{$regions}, $i, 1 );
	}
}

__END__

=head1 NAME

get_gene_regions.pl

A program to collect specific, often un-annotated, regions from genes.

=head1 SYNOPSIS

get_gene_regions.pl [--options...] --in E<lt>filenameE<gt> --out E<lt>filenameE<gt>
  
get_gene_regions.pl [--options...] --db <text> --out E<lt>filenameE<gt>

  Source data:
  -i --in <filename>            input annotation: GFF3, GTF, genePred, etc
  -d --db <name | filename>     database: name, file.db, or file.sqlite
  
  Feature selection:
  -f --feature <type>           optionally specify gene type or type:source
  -t --transcript               specify the transcript type
       [all|mRNA|ncRNA|snRNA|
       snoRNA|tRNA|rRNA|miRNA|
       lincRNA|misc_RNA]
  -r --region                   specify the gene region to collect
       [tss|tts|cdsStart|cdsStop|
       splice|UTR|exon|
       collapsedExon|altExon|
       uncommonExon|commonExon|
       firstExon|lastExon|intron|
       collapsedIntron|altIntron|
       uncommonIntron|commonIntron|
       firstIntron|lastIntron]
  --gencode                     include only GENCODE tagged genes
  --biotype <regex>             include only specific biotype
  --tsl                         select transcript support level
       [best|best1|best2|best3|
       best4|best5|1|2|3|4|5|NA]
  -u --unique                   select only unique regions
  -l --slop <integer>           duplicate region if within X bp
  -K --chrskip <regex>          skip features from certain chromosomes
  
  Adjustments:
  -b --begin --start integer    specify adjustment to start coordinate
  -e --end --stop integer       specify adjustment to stop coordinate
  
  Output options:
  -o --out <filename>           specify output name
  --bed                         output as a bed6 format
  --bedname                     specify what to use for bed name column
       [genename|geneid|            default is 'featurename'
       transcriptname|transcriptid
       featurename]
  -z --gz                       compress output

  General options:
  -v --version                  print version and exit
  -h --help

=head1 OPTIONS

The command line flags and descriptions:

=head2 Source data

=over 4

=item --in E<lt>filenameE<gt>

Provide a gene table or annotation file, including GTF, GFF, GFF3, UCSC 
refFlat, UCSC genePred or genePredExt, or UCSC knownGene table. Files 
may be gzipped.

=item --db E<lt>textE<gt>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be obtained. Only required if 
an input gene table is not provided.

=back

=head2 Feature selection

=over 4

=item --feature E<lt>typeE<gt>

Specify the parental gene feature type (C<primary_tag>) or C<type:source> when
using a database. If not specified, a list of available types will be
presented interactively to the user for selection. This is not relevant for
GFF3 source files (all gene or transcript features are considered). This is 
helpful when gene annotation from multiple sources are present in the same 
database, e.g. refSeq and ensembl sources. More than one feature may be 
included, either as a comma-delimited list or multiple options.

=item --transcript E<lt>typeE<gt>

Specify the transcript type (usually a gene subfeature) from which to  
collect the regions. Multiple types may be specified as a comma-delimited 
list, or 'all' may be specified. If not specified, an interactive list 
will be presented from which the user may select. Available options include:

  all
  mRNA
  ncRNA
  snRNA
  snoRNA
  tRNA
  rRNA
  miRNA
  lincRNA
  misc_RNA
 
=item --region E<lt>regionE<gt>

Specify the type of region to retrieve. If not specified on the command 
line, the list is presented interactively to the user for selection. The 
possibilities are listed below.
     
  tss           The first base of transcription
  tts           The last base of transcription
  exon          The exons of each transcript
  collapsedExon The exons after collapsing all gene transcripts
  firstExon     The first exon of each transcript
  lastExon      The last exon of each transcript
  altExon       Exons unique to one of several transcripts from a gene
  uncommonExon  Exons shared by 2 or more but not all transcripts
  commonExon    Exons shared by all transcripts from a gene
  intron        Each intron (usually not defined in the GFF3)
  collapsedIntron Introns after collapsing all gene transcripts
  firstIntron   The first intron of each transcript
  lastIntron    The last intron of each transcript
  altIntron     Introns unique to one of several transcripts from a gene
  uncommonIntron Introns shared by 2 or more but not all transcripts
  commonIntron  Introns shared by all transcripts of a gene
  splice        The first and last base of each intron
  UTR           The untranslated regions of each coding transcript
  cdsStart      The first base of the CDS
  cdsStop       The last base of the CDS

=item --gencode

Boolean option to filter transcripts as part of the GENCODE specification. 
These are marked in Ensembl GTF/GFF3 annotation files as the C<tag> attribute 
with value "basic". Typically, at least one transcript for every gene is 
marked as part of the GENCODE set. Transcripts not marked as such usually 
lack sufficient experimental evidence.

=item --biotype E<lt>regex<gt> 

Filter transcripts using the C<transcript_biotype> or C<biotype> 
GTF/GFF3 attribute, typically found in Ensembl annotation files. Provide 
a regex compatible string which must match the biotype value to keep the 
transcripts. For example, to keep specify "miRNA" to keep all micro-RNA 
transcripts. This works on a subfeature level as well, so that C<gene> 
may be specified as the feature to collect, and only the gene transcripts 
belonging to the indicating biotype are retained.

=item --tsl E<lt>levelE<gt>

Filter transcripts on the Ensembl GTF/GFF3 attribute 'transcript_support_level', 
which is described at L<Ensembl TSL glossary entry|http://uswest.ensembl.org/info/website/glossary.html>.
Provide a level of support to filter. Values include: 
    
    1       All splice junctions supported by evidence
    2       Transcript flagged as suspect or only support from multiple ESTs
    3       Only support from single EST
    4       Best supporting EST is suspect
    5       No support
    best    Transcripts at the best (lowest) available level are taken
    best1   The word followed by a digit 1-5, indicating any transcript 
            at or better (lower) than the indicated level
    NA      Only transcripts without a level (NA) are retained.

=item --unique

Compare start and stop coordinates of each collected region from 
each feature and remove duplicate regions. When the --slop option 
is provided, only the start coordinate plus/minus the slop factor 
is checked. 

=item --slop E<lt>integerE<gt>

When identifying unique regions, specify the number of bp to 
add and subtract to the start position (the slop or fudge factor) 
of the regions when considering duplicates. Any other region 
within this window will be considered a duplicate. Useful, for 
example, when start sites of transcription are not precisely mapped, 
but not useful with defined introns and exons. This does not take 
into consideration transcripts from other genes, only the current 
gene. The default is 0 (no sloppiness).

=item --chrskip E<lt>regexE<gt>

Exclude features from the output whose sequence ID or chromosome matches 
the provided regex-compatible string. Expressions should be quoted or 
properly escaped on the command line. Examples might be 
    
    'chrM'
    'scaffold.+'
    'chr.+alt|chrUn.+|chr.+_random'

=back

=head2 Adjustments

=over 4

=item --start E<lt>integerE<gt>

=item --begin E<lt>integerE<gt>

=item --stop E<lt>integerE<gt>

=item --end E<lt>integerE<gt>

Optionally specify adjustment values to adjust the reported start and 
end coordinates of the collected regions. A negative value is shifted 
upstream (5' direction), and a positive value is shifted downstream.
Adjustments are made relative to the feature's strand, such that 
a start adjustment will always modify the feature's 5'end, either 
the feature startpoint or endpoint, depending on its orientation. 

=back

=head2 Output options

=over 4

=item --out E<lt>filenameE<gt>

Specify the output filename.

=item --bed

Automatically convert the output file to a BED file.

=item --bedname E<lt>name<gt>

Specify what to use for the Name column in the output BED file.
Several options are available, including:

    geneid          - The Primary ID of the parent Gene feature
    genename        - The Display Name of the parent Gene feature
    transcriptid    - The Primary ID of the parent Transcript feature
    transcriptname  - The Display Name of the parent Transcript feature
    featurename     - The generated name of the feature (default)

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=back

=head2 General options

=over 4

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will collect specific regions from annotated genes and/or 
transcripts. Often these regions are not explicitly defined in the 
source GFF3 annotation, necessitating a script to pull them out. These 
regions include the start and stop sites of transcription, introns, 
the splice sites (both 5' and 3'), exons, the first (5') or last (3') 
exons, or all alternate or common exons of genes with multiple 
transcripts. Importantly, unique regions may only be reported, 
especially important when a single gene may have multiple alternative 
transcripts. A slop factor is included for imprecise annotation.

The program will report the chromosome, start and stop coordinates, 
strand, name, and parent and transcript names for each region 
identified. The reported start and stop sites may be adjusted with 
modifiers. A standard biotoolbox data formatted text file is generated. 
This may be converted into a standard BED or GFF file using the 
appropriate biotoolbox scripts. The file may also be used directly in 
data collection. 

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
