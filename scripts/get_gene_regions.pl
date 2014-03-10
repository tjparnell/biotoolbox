#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use Bio::ToolBox::data_helper qw(
	generate_tim_data_structure
	format_with_commas
);
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
);
use Bio::ToolBox::db_helper::gff3_parser;
use Bio::ToolBox::file_helper qw(
	open_to_read_fh
	write_tim_data_file
);
my $VERSION = '1.15';

print "\n This program will get specific regions from features\n\n";

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
	$database,
	$feature,
	$request,
	$transcript_type,
	$start_adj,
	$stop_adj,
	$unique,
	$slop,
	$bed,
	$gz,
	$help,
	$print_version,
);

# Command line options
GetOptions( 
	'in=s'      => \$infile, # the input data file
	'out=s'     => \$outfile, # name of output file 
	'db=s'      => \$database, # source annotation database
	'feature=s' => \$feature, # the gene feature from the database
	'region=s'  => \$request, # the region requested
	'transcript=s' => \$transcript_type, # which transcripts to take
	'start=i'   => \$start_adj, # start coordinate adjustment
	'stop=i'    => \$stop_adj, # stop coordinate adjustment
	'unique!'   => \$unique, # boolean to ensure uniqueness
	'slop=i'    => \$slop, # slop factor in bp to identify uniqueness
	'bed!'      => \$bed, # convert the output to bed format
	'gz!'       => \$gz, # compress output
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
	print " Biotoolbox script get_gene_regions.pl, version $VERSION\n\n";
	exit;
}



### Check for requirements and set defaults
unless ($infile or $database) {
	die " must define a database or input GFF3 file! use --help for more information\n";
}
if ($database =~ /\.gff3?(?:\.gz)?$/) {
	# a gff3 was specified as the database
	# intercept and assign to input file name
	# faster than trying to load the gff3 file into a memory db
	$infile = $database;
	$database = undef;
}

unless ($outfile) {
	die " must define an output file name! use --help for more information\n";
}

unless ($transcript_type) {
	$transcript_type = 'mRNA';
}

unless (defined $slop) {
	$slop = 0;
}

unless (defined $gz) {
	$gz = 0;
}



### Determine methods and transcript types
# Determine request_method
my $method = determine_method();

# detemine transcript type to take
my ($do_mrna, $do_mirna, $do_ncrna, $do_snrna, $do_snorna, $do_trna, $do_rrna);
determine_transcript_types();


### Collect feature regions
# collection
print " Collecting ";
print $unique ? "unique " : "";
print "$request regions...\n";

my $outdata;
if ($database) {
	$outdata = collect_from_database($method);
}
elsif ($infile) {
	$outdata = collect_from_file($method);
}
print "  collected ", format_with_commas($outdata->{'last_row'}), " regions\n";





### Finished
my $success = write_tim_data_file(
	'data'     => $outdata,
	'filename' => $outfile,
	'gz'       => $gz,
);
if ($success) {
	print " wrote file '$success'\n";
}
else {
	# failure! the subroutine will have printed error messages
	print " unable to write file!\n";
}


### Convert to bed format if requested
# rather than taking the time to modify the data structures and all the 
# the data collection subroutines to a BED format, we'll just simply 
# take advantage of the data2bed.pl program as a convenient cop-out
if ($bed and $success) {
	system(
		"$Bin/data2bed.pl",
		"--chr",
		3,
		"--start",
		4,
		"--stop",
		5,
		"--strand",
		6,
		"--name",
		2,
		"--in",
		$success,
		$gz ? "--gz" : "",
	) == 0 or warn " unable to execute data2bed.pl for converting to bed!\n";
}



########################   Subroutines   ###################################

sub determine_method {
	
	# determine the region request from user if necessary
	unless ($request) {
		$request = collect_method_from_user();
	}
	
	# determine the method
	# also change the name of the request from short to long form
	my $method;
	if ($request =~ /^first ?exon$/i) {
		$request = 'first exon';
		$method = \&collect_first_exon;
	}
	elsif ($request =~ /^second ?exon$/i) {
		$request = 'last exon';
		$method = \&collect_second_exon;
	}
	elsif ($request =~ /tss/i) {
		$request = 'transcription start site';
		$method = \&collect_tss;
	}
	elsif ($request =~ /start site/i) {
		$method = \&collect_tss;
	}
	elsif ($request =~ /tts/i) {
		$request = 'transcription stop site';
		$method = \&collect_tts;
	}
	elsif ($request =~ /stop site/i) {
		$method = \&collect_tts;
	}
	elsif ($request =~ /^splices?/i) {
		$request = 'splice sites';
		$method = \&collect_splice_sites;
	}
	elsif ($request =~ /^introns?$/i) {
		$method = \&collect_introns;
	}
	elsif ($request =~ /^first ?intron/i) {
		$request = 'first intron';
		$method = \&collect_first_intron;
	}
	elsif ($request =~ /^last ?intron/i) {
		$request = 'last intron';
		$method = \&collect_last_intron;
	}
	else {
		die " unknown region request!\n";
	}
	
	return $method;
}



sub collect_method_from_user {
	
	my %list = (
		1	=> 'first exon',
		2	=> 'last exon',
		3	=> 'transcription start site',
		4	=> 'transcription stop site',
		5	=> 'splice sites',
		6	=> 'introns',
		7   => 'first intron',
		8   => 'last intron',
	);
	
	# request feature from the user
	print " These are the available feature types in the database:\n";
	foreach my $i (sort {$a <=> $b} keys %list ) {
		print "   $i\t$list{$i}\n";
	}
	print " Enter the type of region to collect   ";
	my $answer = <STDIN>;
	chomp $answer;
	
	# verify and return answer
	if (exists $list{$answer}) {
		return $list{$answer};
	}
	else {
		die " unknown request!\n";
	}
}



sub determine_transcript_types {
	
	# string for visual output
	my $string = " Collecting transcript types:";
	
	# collect all the transcript types requested
	my @types;
	if ($transcript_type =~ /,/) {
		@types = split ",", $transcript_type;
	}
	else {
		push @types, $transcript_type;
	}
	
	foreach (@types) {
		if (m/^mRNA|all$/i) {
			$do_mrna   = 1;
			$string .= ' mRNA';
		}
		if (m/^miRNA|all$/i) {
			$do_mirna  = 1;
			$string .= ' miRNA';
		}
		if (m/^ncRNA|all$/i) {
			$do_ncrna  = 1;
			$string .= ' ncRNA';
		}
		if (m/^snRNA|all$/i) {
			$do_snrna  = 1;
			$string .= ' snRNA';
		}
		if (m/^snoRNA|all$/i) {
			$do_snorna = 1;
			$string .= ' snoRNA';
		}
		if (m/^tRNA|all$/i) {
			$do_trna   = 1;
			$string .= ' tRNA';
		}
		if (m/^rRNA|all$/i) {
			$do_rrna   = 1;
			$string .= ' rRNA';
		}
	}
	print "$string\n";
}



sub collect_from_database {
	
	# collection method
	my $method = shift;
	
	# open database connection
	my $db = open_db_connection($database) or 
		die " unable to open database connection!\n";
	
	# get feature type if necessary
	$feature = verify_or_request_feature_types(
		'db'      => $db,
		'feature' => $feature,
		'prompt'  => 'Enter the gene feature from which to collect regions   ',
		'single'  => 1,
	) or die "No valid gene feature type was provided! see help\n";
	
	# generate output data
	my $output = generate_output_structure();
	$output->{'db'} = $database;
	$output->{0}{'feature'} = $feature; # parent column
	
	# generate a seqfeature stream
	my $iterator = $db->features(
		-type     => $feature,
		-iterator => 1,
	);
	
	# process the features
	while (my $seqfeat = $iterator->next_seq) {
		# collect the regions based on the primary tag and the method re
		if ($seqfeat->primary_tag eq 'gene') {
			# gene
			push @{ $output->{'data_table'} },
				process_gene($seqfeat, $method);
		}
		elsif ($seqfeat->primary_tag =~ /rna/i) {
			# transcript
			my @regions = process_transcript($seqfeat, $method);
			
			# add the parent name
			map { unshift @$_, $seqfeat->display_name } @regions;
			
			push @{ $output->{'data_table'} }, @regions;
		}
	}
	$output->{'last_row'} = scalar( @{ $output->{'data_table'} }) - 1;
	
	# finished
	return $output;
}



sub collect_from_file {

	# collection method
	my $method = shift;
	
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
	
	# generate output data
	my $output = generate_output_structure();
	push @{ $output->{'other'} }, "Source data file $infile\n";
	
	# open gff3 parser object
	my $parser = Bio::ToolBox::db_helper::gff3_parser->new($infile) or
		die " unable to open input file '$infile'!\n";
	
	
	# process the features
	while (my @top_features = $parser->top_features() ) {
		
		# Process the top features
		while (@top_features) {
			my $seqfeat = shift @top_features;
		
			# collect the regions based on the primary tag and the method re
			if ($seqfeat->primary_tag =~ /^gene$/i) {
				# gene
				push @{ $output->{'data_table'} },
					process_gene($seqfeat, $method);
			}
			elsif ($seqfeat->primary_tag =~ /rna/i) {
				# transcript
				my @regions = process_transcript($seqfeat, $method);
				
				# add the parent name
				map { unshift @$_, $seqfeat->display_name } @regions;
				
				push @{ $output->{'data_table'} }, @regions;
			}
		}
	}
	$output->{'last_row'} = scalar( @{ $output->{'data_table'} }) - 1;
	
	# finished
	return $output;
}



sub generate_output_structure {
	
	# generate
	my $data = generate_tim_data_structure(
		"region",
		qw(
			Parent
			Transcript
			Name
			Chromosome
			Start
			Stop
			Strand
		)
	);
	
	# add metadata
	$data->{1}{'type'} = $request;
	$data->{1}{'type'} =~ s/\s/_/; # replace spaces
	$data->{2}{'type'} = $transcript_type;
	if ($start_adj) {
		$data->{4}{'start_adjusted'} = $start_adj;
	}
	if ($stop_adj) {
		$data->{5}{'stop_adjusted'} = $stop_adj;
	}
	if ($unique) {
		$data->{2}{'unique'} = 1; # Name
		$data->{2}{'slop'} = $slop;
	}
	
	return $data;
}



sub process_gene {
	
	# passed objects
	my ($gene, $method) = @_;
	
	# look for transcripts for this gene
	my @regions;
	foreach my $subfeat ($gene->get_SeqFeatures) {
		if ($subfeat->primary_tag =~ /rna/i) {
			my @r = process_transcript($subfeat, $method);
			if ($r[0]) {
				push @regions, @r;
			}
		}
	}
	return unless @regions;
# 	warn "  found " . scalar(@regions) . " regions for " . $gene->display_name . "\n";
# 	print Dumper(\@regions);
	
	# remove duplicates if requested
	if ($unique) {
		remove_duplicates(\@regions);
	}
	
	# add the parent name
	for my $i (0 .. $#regions) {
		unshift @{ $regions[$i] }, $gene->display_name;
	}
	
	# return the regions
	return @regions;
}



sub process_transcript {
	
	# passed objects
	my ($transcript, $method) = @_;
	
	if (
		($transcript->primary_tag =~ /mrna/i and $do_mrna) or
		($transcript->primary_tag =~ /mirna/i and $do_mirna) or
		($transcript->primary_tag =~ /ncrna/i and $do_ncrna) or
		($transcript->primary_tag =~ /snrna/i and $do_snrna) or
		($transcript->primary_tag =~ /snorna/i and $do_snorna) or
		($transcript->primary_tag =~ /trna/i and $do_rrna) or
		($transcript->primary_tag =~ /rrna/i and $do_rrna)
	) {
		return &{$method}($transcript);
	}
}



sub collect_tss {
	
	# get seqfeature objects
	my $transcript = shift;
	
	# get coordinates
	my $chromo = $transcript->seq_id;
	my ($start, $stop, $strand);
	if ($transcript->strand == 1) {
		# forward strand
		
		$strand = 1;
		$start = $transcript->start;
		$stop = $transcript->start;
	}
	elsif ($transcript->strand == -1) {
		# reverse strand
		
		$strand = -1;
		$start = $transcript->end;
		$stop = $transcript->end;
	}
	else {
		die " poorly formatted transcript seqfeature object with strand 0!\n";
	}
	
	# get name
	my $name = $transcript->display_name . '_TSS';
	
	return _adjust_positions( 
		[$transcript->display_name, $name, $chromo, $start, $stop, $strand] 
	);
}



sub collect_tts {
	
	# get seqfeature objects
	my $transcript = shift;
	
	# get coordinates
	my $chromo = $transcript->seq_id;
	my ($start, $stop, $strand);
	if ($transcript->strand == 1) {
		# forward strand
		
		$strand = 1;
		$start = $transcript->end;
		$stop = $transcript->end;
	}
	elsif ($transcript->strand == -1) {
		# reverse strand
		
		$strand = -1;
		$start = $transcript->start;
		$stop = $transcript->start;
	}
	else {
		die " poorly formatted transcript seqfeature object with strand 0!\n";
	}
	
	# get name
	my $name = $transcript->display_name . '_TTS';
	
	return _adjust_positions( 
		[$transcript->display_name, $name, $chromo, $start, $stop, $strand] 
	);
}



sub collect_first_exon {
	
	my $transcript = shift;
	
	# find the exons and/or CDSs
	my $list = _collect_exons($transcript);
	return unless $list;
	
	# the first exon
	my $first = shift @{ $list };
	
	# identify the exon name if it has one
	my $name = $first->display_name || 
		$transcript->display_name . "_lastExon";
	
	# finished
	return _adjust_positions( [ 
		$transcript->display_name,
		$name, 
		$first->seq_id, 
		$first->start, 
		$first->end,
		$first->strand,
	] );
}



sub collect_last_exon {
	
	my $transcript = shift;
	
	# find the exons and/or CDSs
	my $list = _collect_exons($transcript);
	return unless $list;
	
	# the last exon
	my $last = pop @{ $list };
	
	# identify the exon name if it has one
	my $name = $last->display_name || 
		$transcript->display_name . "_lastExon";
	
	# finished
	return _adjust_positions( [ 
		$transcript->display_name,
		$name, 
		$last->seq_id, 
		$last->start, 
		$last->end,
		$last->strand,
	] );
}



sub collect_splice_sites {
	
	# seqfeature object
	my $transcript = shift;
	
	# find the exons and/or CDSs
	my $list = _collect_exons($transcript);
	return unless $list;
	return if (scalar(@$list) == 1);
	
	# identify the last exon index position
	my $last = scalar(@$list) - 1;
	
	# collect the splice sites
	my @splices;
	
	# forward strand
	if ($transcript->strand == 1) {
		
		# walk through each exon
		for (my $i = 0; $i <= $last; $i++) {
			
			# get the exon name
			my $exon = $list->[$i];
			my $name = $exon->display_name || 
				$transcript->display_name . ".exon$i";
			
			# first exon
			if ($i == 0) {
				push @splices, _adjust_positions( [ 
					$transcript->display_name,
					$name . '_3\'', 
					$exon->seq_id, 
					$exon->end + 1, 
					$exon->end + 1,
					$exon->strand,
				] );
			}
			
			# last exon
			elsif ($i == $last) {
				push @splices, _adjust_positions( [ 
					$transcript->display_name,
					$name . '_5\'', 
					$exon->seq_id, 
					$exon->start - 1, 
					$exon->start - 1,
					$exon->strand,
				] );
			
			}
			
			# middle exons
			else {
				
				# 5' splice
				push @splices, _adjust_positions( [ 
					$transcript->display_name,
					$name . '_5\'', 
					$exon->seq_id, 
					$exon->start - 1, 
					$exon->start - 1,
					$exon->strand,
				] );
				
				# 3' splice
				push @splices, _adjust_positions( [ 
					$transcript->display_name,
					$name . '_3\'', 
					$exon->seq_id, 
					$exon->end + 1, 
					$exon->end + 1,
					$exon->strand,
				] );
			}
		}
	}
	
	# reverse strand
	else {
		
		# walk through each exon
		for (my $i = 0; $i <= $last; $i++) {
			
			# get the exon name
			my $exon = $list->[$i];
			my $name = $exon->display_name || 
				$transcript->display_name . ".exon$i";
			
			# first exon
			if ($i == 0) {
				push @splices, _adjust_positions( [ 
					$transcript->display_name,
					$name . '_3\'', 
					$exon->seq_id, 
					$exon->start - 1, 
					$exon->start - 1,
					$exon->strand,
				] );
			}
			
			# last exon
			elsif ($i == $last) {
				push @splices, _adjust_positions( [ 
					$transcript->display_name,
					$name . '_5\'', 
					$exon->seq_id, 
					$exon->end + 1, 
					$exon->end + 1,
					$exon->strand,
				] );
			
			}
			
			# middle exons
			else {
				
				# 5' splice
				push @splices, _adjust_positions( [ 
					$transcript->display_name,
					$name . '_5\'', 
					$exon->seq_id, 
					$exon->end + 1, 
					$exon->end + 1,
					$exon->strand,
				] );
				
				# 3' splice
				push @splices, _adjust_positions( [ 
					$transcript->display_name,
					$name . '_3\'', 
					$exon->seq_id, 
					$exon->start - 1, 
					$exon->start - 1,
					$exon->strand,
				] );
			}
		}
	}
	
	# finished
	return @splices;
}



sub collect_introns {
	
	# seqfeature object
	my $transcript = shift;
	
	# find the exons and/or CDSs
	my $exons = _collect_exons($transcript);
	return unless $exons;
	return if (scalar(@$exons) == 1);
	
	# identify the last exon index position
	my $last = scalar(@$exons) - 1;
	
	# collect the introns
	my @introns;
	
	# forward strand
	if ($transcript->strand == 1) {
		
		# walk through each exon
		for (my $i = 0; $i < $last; $i++) {
			push @introns, _adjust_positions( [ 
				$transcript->display_name,
				$transcript->display_name . ".intron$i", 
				$transcript->seq_id, 
				$exons->[$i]->end + 1, 
				$exons->[$i + 1]->start - 1,
				$transcript->strand,
			] );
		}
	}
	
	# reverse strand
	else {
	
		# walk through each exon
		for (my $i = 0; $i < $last; $i++) {
			push @introns, _adjust_positions( [ 
				$transcript->display_name,
				$transcript->display_name . ".intron$i", 
				$transcript->seq_id, 
				$exons->[$i + 1]->end + 1,
				$exons->[$i]->start - 1, 
				$transcript->strand,
			] );
		}
	}
	
	# finished
	return @introns;
}


sub collect_first_intron {
	
	# seqfeature object
	my $transcript = shift;
	
	# collect all of the introns
	my @introns = collect_introns($transcript);
	
	# return the first one
	return shift @introns;
}


sub collect_last_intron {
	
	# seqfeature object
	my $transcript = shift;
	
	# collect all of the introns
	my @introns = collect_introns($transcript);
	
	# return the first one
	return pop @introns;
}


sub _collect_exons {
	
	# initialize
	my $transcript = shift;
	my @exons;
	my @cdss;
	
	# go through the subfeatures
	foreach my $subfeat ($transcript->get_SeqFeatures) {
		if ($subfeat->primary_tag eq 'exon') {
			push @exons, $subfeat;
		}
		elsif ($subfeat->primary_tag =~ /cds|utr|untranslated/i) {
			push @cdss, $subfeat;
		}
	}
	
	# check which array we'll use
	# prefer to use actual exon subfeatures, but those may not be defined
	my $list;
	if (@exons) {
		$list = \@exons;
	}
	elsif (@cdss) {
		$list = \@cdss;
	}
	else {
		# nothing found!
		return;
	}
	
	# sort the list using a Schwartzian transformation by stranded start position
	my @sorted;
	if ($transcript->strand == 1) {
		# forward strand, sort by increasing start positions
		@sorted = 
			map { $_->[0] }
			sort { $a->[1] <=> $b->[1] }
			map { [$_, $_->start] } 
			@{ $list };
	}
	else {
		# reverse strand, sort by decreasing end positions
		@sorted = 
			map { $_->[0] }
			sort { $b->[1] <=> $a->[1] }
			map { [$_, $_->end] }
			@{ $list };
	}
	
	return \@sorted;
}



sub _adjust_positions {
	
	my $region = shift;
	# region is an anonymous array of 5 elements
	# [$transcript_name, $name, $chromo, $start, $stop, $strand]
	
	# adjust the start and end positions according to strand
	if ($region->[5] == 1) {
		# forward strand
		
		if ($start_adj) {
			$region->[3] += $start_adj;
		}
		if ($stop_adj) {
			$region->[4] += $stop_adj;
		}
	}
	elsif ($region->[5] == -1) {
		# reverse strand
		
		if ($start_adj) {
			$region->[4] -= $start_adj;
		}
		
		# stop
		if ($stop_adj) {
			$region->[3] -= $stop_adj;
		}
	}
	
	# return adjusted region coordinates
	return $region;
}



sub remove_duplicates {
	
	my $regions = shift;
	
	# look for duplicates using a quick hash of seen positions
	my %seenit;
	my @to_remove;
	for my $i (0 .. $#{ $regions } ) {
		# we will be using the start position as a unique identifier
		# to account for the slop factor,
		# we'll be adding/subtracting the slop value to/from the start position
		# if this position matches anything else, we'll assume it's a duplicate
		
		foreach my $pos ( 
			# generate an array of possible start positions
			# with a default slop of 0, this will only be 1 position
			($regions->[$i]->[3] - $slop) .. ($regions->[$i]->[3] + $slop)
		) {
			if (exists $seenit{ $pos }) {
				push @to_remove, $i;
			}
			else {
				$seenit{ $pos } = 1;
			}
		}
	}
	
	# remove the duplicates
	while (@to_remove) {
		my $i = pop @to_remove; 
			# take from end to avoid shifting regions array
		splice( @{$regions}, $i, 1);
	}
}


__END__

=head1 NAME

get_gene_regions.pl

A script to collect specific, often un-annotated regions from genes.

=head1 SYNOPSIS

get_gene_regions.pl [--options...] --db <text> --out <filename>

get_gene_regions.pl [--options...] --in <filename> --out <filename>
  
  Options:
  --db <text>
  --in <filename>
  --out <filename> 
  --feature <type | type:source>
  --transcript [all|mRNA|miRNA|ncRNA|snRNA|snoRNA|tRNA|rRNA]
  --region [tss|tts|firstExon|lastExon|splice|intron|firstIntron|lastIntron]
  --start=<integer>
  --stop=<integer>
  --unique
  --slop <integer>
  --bed
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <text>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A database is 
required for generating new data files with features. For more information 
about using annotation databases, 
see L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. 
Also see C<--in> as an alternative.

=item --in <filename>

Alternative to a database, a GFF3 annotation file may be provided. 
For best results, the database or file should include hierarchical 
parent-child annotation in the form of gene -> mRNA -> [exon or CDS]. 
The GFF3 file may be gzipped.

=item --out <filename>

Specify the output filename.

=item --feature <type | type:source>

Specify the parental gene feature type (primary_tag) or type:source when
using a database. If not specified, a list of available types will be
presented interactively to the user for selection. This is not relevant for
GFF3 source files (all gene or transcript features are considered). Helpful
when gene annotation from multiple sources are listed in the same database,
e.g. refSeq and ensembl sources.

=item --transcript [all|mRNA|miRNA|ncRNA|snRNA|snoRNA|tRNA|rRNA]

Specify the transcript feature or gene subfeature type from which to  
collect the regions. Multiple types may be specified as a comma-delimited 
list, or 'all' may be specified. The default value is mRNA.

=item --region [tss|tts|firstExon|lastExon|splice|intron|firstIntron|lastIntron]

Specify the type of region to retrieve. If not specified on the command 
line, the list is presented interactively to the user for selection. Six 
possibilities are possible.
     
     tss         The first base of transcription
     tts         The last base of transcription
     firstExon   The first exon of each transcript
     lastExon    The last exon of each transcript
     splice      The first and last base of each intron
     intron      Each intron (usually not defined in the GFF3)
     firstIntron The first intron of each transcript
     lastIntron  The last intron of each transcript

=item --start=<integer>

=item --stop=<integer>

Optionally specify adjustment values to adjust the reported start and 
end coordinates of the collected regions. A negative value is shifted 
upstream (5' direction), and a positive value is shifted downstream.
Adjustments are made relative to the feature's strand, such that 
a start adjustment will always modify the feature's 5'end, either 
the feature startpoint or endpoint, depending on its orientation. 

=item --unique

For gene features only, take only the unique regions. Useful when 
multiple alternative transcripts are defined for a single gene.

=item --slop <integer>

When identifying unique regions, specify the number of bp to 
add and subtract to the start position (the slop or fudge factor) 
of the regions when considering duplicates. Any other region 
within this window will be considered a duplicate. Useful, for 
example, when start sites of transcription are not precisely mapped, 
but not useful with defined introns and exons. This does not take 
into consideration transcripts from other genes, only the current 
gene. The default is 0 (no sloppiness).

=item --bed

Automatically convert the output file to a BED file.

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

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
the splice sites (both 5' and 3'), and the first and last exons. 
Importantly, unique regions may only be reported, especially important 
when a single gene may have multiple alternative transcripts. A 
slop factor is included for imprecise annotation.

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
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
