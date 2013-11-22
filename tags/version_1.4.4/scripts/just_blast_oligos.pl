#!/usr/bin/perl -w
# This script will use BLASTn to align a list of oligos back to the genome
# and simply report back the hits
# Customized for any microarray

# warning: this script is based on some very old code and has some less than
# desirable coding styles, and is thus subsequently patched and updated many times.
# it's really due for a complete overhaul and rewrite, but I'm too lazy

use Statistics::Lite qw(:all);
use Getopt::Long;
use DateTime;

print "\n  This program will align microarray oligo probes to a genome using BLAST.\n";

# UPDATED!!!! This script is no longer limited to these specific columns. The only requirement
# is that the first line is column headings, and that two columns must be Probe and Sequence. 
# You can have as many or few columns in whatever ordeer as you desire! Enjoy!!!

# Get command line options --New!
my $filename;
my $database;
my $gap = 'T';
my $gapext = 2;
my $filter = 'F';
my $expect = 0.000001;
my $length =  45;
my $number_analyze = 100;
my $raw;
my $log;
my $help;

GetOptions('file=s' => \$filename,
			'database=s' => \$database,
			'number=i' => \$number_analyze,
			'gap=s' => \$gap,
			'gapext=s' => \$gapext,
			'filter=s' => \$filter,
			'expect=s' => \$expect,
			'length=i' => \$length,
			'raw' => \$raw,
			'log' => \$log,
			'help' => \$help,
			);

if ($help) { # print the online help documentation and exit the program
	print "\n Command line options for just_blast_oligos.pl\n";
	print "  --file      the filename for a tabbed delimited text file containing the oligo\n";
	print "              list. The first line must be a header line, and at least 2 columns\n";
	print "              must be present: Probe and Sequence (in any order)\n";
	print "  --database  name of the blast database, which must be found by BLAST\n";
	print "  --gap       T or F allow gaps in the search. default is T\n";
	print "  --gapext    Gap extension penalty. default is 2\n";
	print "  --filter    T or F filter the oligo query for poor sequences. default is F\n";
	print "  --expect    The Expect score cutoff. default is 1e-6\n";
	print "  --length    The length of homology expected. default is 45\n";
	print "  --number    The number of oligos to feed BLAST at a time. default is 100\n";
	print "  --raw       save the raw blast results too\n";
	print "  --log       print out a log file during the run\n";
	print "  --help      This help text\n";
	exit;
}
unless ($filename) {die "No file name specified! run again with --help flag\n"};
unless ($database) {die "No database specified! run again with --help flag\n"};
$filename =~ s/\.txt$//;

$start_time = time; #record start time
if ($log) {open LOG, ">$filename.blast_log.txt"};

# Delete any possible output files from  previous runs
if (-e "$filename.rawblast.txt") {unlink "$filename.rawblast.txt"};

# Input the oligo data and generate hash
$oligo_filename = "$filename.txt";
open OLIGO_FILE, $oligo_filename or die "Can't open $oligo_filename!\n";
@oligo_input = <OLIGO_FILE>;
close OLIGO_FILE;

$header = shift @oligo_input; # grab the column headers
chomp $header;
@header_names = split /\t/, $header;
$head_number = scalar @header_names;
my ($probe_index, $seq_index, $desc_index) = (9999, 9999, 9999);
$i = 0;
foreach (@header_names) {
	if (/probe/i) { $probe_index = $i } # look for which column header has the probe name or id
	if (/seq/i) { $seq_index = $i } # look for which column header has the sequence
	if ( (/desc/i) or (/coord/i) ) { $desc_index = $i } # look for which column header has the description
	$i++;
}
if ($probe_index == 9999) {die "\n Can't find which column has the probe name!\n"}
if ($seq_index == 9999) {die "\n Can't find which column has the sequence!\n"}
# we'll be adding on two more columns of data to the end: Number of hits and list of hits
my $numhit_index = $head_number; # the index number for the Number_Hits, which is $head_number + 1 - 1
my $hits_index = $head_number + 1; # the index number for the Number_Hits, which is $head_number + 2 - 1
foreach (@oligo_input) {
	chomp;
	@oligo_data = split /\t/;
	# adding default data elements: Number_Hits, Hits
	push @oligo_data, qw(0 .);
	# Each line now has (for example): Probe_ID, Description, Chromosome, Start, End, Sequence, Number_Hits, Hits
	$oligos{$oligo_data[$probe_index]} = join "\t", @oligo_data; # put it into a hash
}
@oligo_input = (); # clear memory
@oligo_keys = keys %oligos;
$number_keys = @oligo_keys;
if ($log) {print LOG "Loaded $number_keys oligos from file $oligo_filename.\n"};

$total_blast_number = 0;
$extra_hits = 0;

# Work through the oligos and blast them
while (@oligo_keys) {
	if (scalar(@oligo_keys) < $number_analyze) { # in case there isn't 100 left
		$number_analyze = scalar(@oligo_keys);
	}
	# get current batch of oligos to blast
	for ($count = 1; $count <= $number_analyze; $count++) {
		$current_oligo = shift @oligo_keys; 
		@oligo_data = split /\t/, $oligos{$current_oligo};
		push @temp_fasta, ">$oligo_data[0]\n$oligo_data[$seq_index]\n"; #fasta format
	}
	# write out temp fasta file
	open TEMP_FASTA_FILE, ">temporary.fasta";
	print TEMP_FASTA_FILE @temp_fasta;
	close TEMP_FASTA_FILE;
	
	# run the blastn program against pombe database using the fasta file as input
	my $gapparameter;
	if ($gap eq 'T') {
		$gapparameter = "-g T -E $gapext ";
	} elsif ($gap eq 'F') {
		$gapparameter = "-g F ";
	}
	system "/usr/local/blast/bin/blastall -p blastn -d $database -i temporary.fasta -a 2 $gapparameter -F $filter -e $expect -m 9 -o temporary.blast";
	open BLAST_RESULTS, "temporary.blast";
	@blast_input = <BLAST_RESULTS>;
	close BLAST_RESULTS;
	if ($raw) {system "cat temporary.blast >> $filename.rawblast.txt"};
	unlink "temporary.blast";
	unlink "temporary.fasta";
	$blast_number = scalar(@blast_input);
	$total_blast_number += $blast_number; # keep a running total of the number of hits
	$remainder = scalar @oligo_keys;
	if ($log) {print LOG " Found $blast_number hits for $number_analyze oligos. $remainder oligos remaining.\n"};
	
	# Parse the output and look for extra genomic hits
	foreach (@blast_input) {
		chomp;
		if (/^#/) {next}; # skip comment line
		@blast_data = split /\t/;
		# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
		if ($blast_data[3] < $length) {next} # ignore those partial matches of less than X bp
		@oligo_data = split /\t/, $oligos{$blast_data[0]}; # retrieve the oligo data for this BLAST hit from the oligo hash
		$blast_chr = $blast_data[1];
		$blast_chr =~ s/\D//g; # isolate the chromosome number from subject name
		if ($blast_data[8] < $blast_data[9]) { # determine strand
			$strand = '+';
			$blast_start = $blast_data[8];
			$blast_end = $blast_data[9];
		} else { 
			$strand = '-';
			$blast_start = $blast_data[9];
			$blast_end = $blast_data[8];
		}
		$identity = $blast_data[2]; # set percent identity
		$align_length = $blast_data[3]; # set the alignment length
		$mismatch = $blast_data[4]; # set the number of mismatches
		$gaps = $blast_data[5]; # set the number of gaps
		$oligo_data[$numhit_index]++; #  increment oligo hit counter
		# generate name of extra hit chromosome:start-end (alignment length, % identity, strand)
		if ($oligo_data[$hits_index] eq '.') { # no extragenomic hits were previously recorded
			$oligo_data[$hits_index] = "chr$blast_chr:$blast_start-$blast_end ($align_length, $identity, $strand)";
		} else {
			$oligo_data[$hits_index] .= "; chr$blast_chr:$blast_start-$blast_end ($align_length, $identity, $strand)";
		}
		
		# output new oligo data from blast hit
		# hit_output is chromosome, start, end, strand, probe_id, description, identity, alignment length, mismatch, gaps
		delete $oligos{$blast_data[0]};
		$oligos{$blast_data[0]} = (join "\t", @oligo_data); # put the edited oligo data back into the hash
		$blast_start = sprintf "%07d", $blast_start; # pad with numbers to make sorting easier
		$blast_end = sprintf "%07d", $blast_end;
		$blast_chr = sprintf "%02d", $blast_chr;
		my $description;
		if ($desc_index == 9999) {
			$description = '.';
		} else {
			$description = $oligo_data[$desc_index];
		}
		push @hit_output, "$blast_chr\t$blast_start\t$blast_end\t$strand\t$oligo_data[0]\t$description\t$identity\t$align_length\t$mismatch\t$gaps";
		$extra_hits++; # increment total accepted hits counter	
	} # end the processing of the blast results
	@temp_fasta = ();
}
	
# Output the oligo results
@oligo_output = values %oligos;
open OLIGO_FILE, ">$filename\_blast.txt";
print OLIGO_FILE "$header\tNumber_Hits\tHits\n";
foreach (sort @oligo_output) {
	print OLIGO_FILE "$_\n";
}
close OLIGO_FILE;

open HIT_FILE, ">$filename\_sorted_hits.txt";
print HIT_FILE "Chromo\tStart\tEnd\tStrand\tProbe_ID\tDescription\tIdentity\tAlignment_Length\tMismatch\tGaps\tNumber_Hits\n";
foreach (sort @hit_output) {
	@data = split /\t/; # get the current probe_id
	my $number_hit = (split /\t/, $oligos{$data[4]})[$numhit_index]; # get the number of blast hits for the current probe_id
	$data[0] =~ s/^0+//;
	$data[1] =~ s/^0+//;
	$data[2] =~ s/^0+//;
	print HIT_FILE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\t$data[8]\t$data[9]\t$number_hit\n";
}
close HIT_FILE;


open RESULTS, ">$filename\_results.txt";
my $dt = DateTime->now(time_zone => 'America/Denver');
$dt =~ s/T/ /; # replace T with space
print RESULTS "Running just_blast_oligos on $dt\n";
print RESULTS "Using oligo data file $filename.\n";
print RESULTS "BLASTed $number_keys oligos against the $database genome.\n";
print RESULTS "  Using following parameters:\n\t gapped search is $gap\n\t gap extension penalty of $gapext\n\t filter is $filter\n\t expect cutoff of $expect\n\t discarding hits < $length bp\n";
print RESULTS "Found $total_blast_number blast hits for $number_keys oligos.\n";
print RESULTS "Kept $extra_hits out of $total_blast_number blast hits.\n\n";


### Perform autopsy on results ###

# set values
$no_one = 0; # matches no one at all
$unique = 0; # matches only one spot
$multicopy = 0; # matches multiple spots
$total_oligos = 0;
$complete = 0; # the number of perfect matches
$partial = 0; # the number of oligos that don't have 100 % identity
$mismatches = 0; # the number of oligos with mismatches
$gaps = 0; # the number of oligos with gaps
$misgaps = 0; # the number of oligos with both gaps and mismatches

foreach (@oligo_output) {
	chomp;
	@data = split /\t/; # the data array should consist of the following
	
	if ($data[$numhit_index] == 0) {
		$no_one++;
		push @losers, $oligos{$data[0]};
	} elsif ($data[$numhit_index] == 1) {
		$unique++;
	} elsif ($data[$numhit_index] > 1) {
		$multicopy++;
		push @numbers_of_hits, $data[$numhit_index];
	}
	$total_oligos++;
}
foreach (@hit_output) {
	@data = split /\t/;
	if ($data[6] == 100.00) {$complete++} # count the number of perfect hits (100$ identity)
	else {$partial++} # count those with partial matches (not 100% identity)
	if ($data[8] != 0) {$mismatches++} # count those with mismatches
	if ($data[9] != 0) {$gaps++} # count those with gaps
	if ( ($data[8] != 0) and ($data[9] != 0) ) {$misgaps++}
	$total_hits++;
}

$mean_hit = sprintf "%.1f", &mean(@numbers_of_hits);
$min_hit = &min(@numbers_of_hits);
$max_hit = &max(@numbers_of_hits);
$median_hit = sprintf "%.1f", &median(@numbers_of_hits);
$mode_hit = &mode(@numbers_of_hits);

$percent_no_one = sprintf "%.2f", ($no_one/$total_oligos) * 100;
$percent_unique = sprintf "%.2f", ($unique/$total_oligos) * 100;
$percent_multicopy = sprintf "%.2f", ($multicopy/$total_oligos) * 100;

$percent_complete = sprintf "%.2f", ($complete/$total_hits) * 100;
$percent_partial = sprintf "%.2f", ($partial/$total_hits) * 100;
$percent_mismatches = sprintf "%.2f", ($mismatches/$total_hits) * 100;
$percent_gaps = sprintf "%.2f", ($gaps/$total_hits) * 100;
$percent_misgaps = sprintf "%.2f", ($misgaps/$total_hits) * 100;


print RESULTS "\nHere are the summaries for the $filename file:\n";
print RESULTS " Total number of oligos: $total_oligos\n";
print RESULTS " Number of oligos with no blast hit: $no_one ($percent_no_one %)\n";
print RESULTS " Number of unique oligos: $unique ($percent_unique %)\n";
print RESULTS " Number of multicopy oligos: $multicopy ($percent_multicopy %)\n";
print RESULTS "   Range of numbers of genomic hits: $min_hit to $max_hit\n";
print RESULTS "   Mean number of genomic hits: $mean_hit\n";
print RESULTS "   Median number of genomic hits: $median_hit\n";
print RESULTS "   Mode number of genomic hits: $mode_hit\n";
print RESULTS "\nThere were $total_hits blast hits of oligos to the genome\n";
print RESULTS " Number of hits with 100 % identity: $complete ($percent_complete %)\n";
print RESULTS " Number of hits with partial identity: $partial ($percent_partial %)\n";
print RESULTS " Number of hits with mismatches: $mismatches ($percent_mismatches %)\n";
print RESULTS " Number of hits with gaps: $gaps ($percent_gaps %)\n";
print RESULTS " Number of hits with both mismatches and gaps: $misgaps ($percent_misgaps %)\n";
print RESULTS "\nThat's all.\n\n";

$elapsed = (time - $start_time) / 60;
$elapsed = sprintf "%.2f", $elapsed;

print "\n I'm all done in $elapsed minutes!\n";
print RESULTS "\n Performed in $elapsed minutes.\n";
if ($log) {close LOG};
close RESULTS;

# report the losers...er, I mean, ones without any matches
if (@losers) {
	open LOSER_FILE, ">$filename\_no_hits.txt";
	print LOSER_FILE "Probe_ID\tDescription\tChromosome\tStart\tEnd\tSequence\tNumber_Hits\tHits\n";
	foreach (sort @losers) {
		print LOSER_FILE "$_\n";
	}
	close LOSER_FILE;
}
# THE END

# arrays indexing

#  index	oligo_data   -> Be aware that this script is now column number agnostic
#  0	Probe_ID  -> must be the first column now
#  1	Description -> looked up and referred to as $desc_index
#  2	Chromosome
#  3	Start
#  4	End
#  5	Sequence  -> looked up and referred to as $seq_index
#  6	Number_Hits -> referred to as $numhit_index
#  7	Hits -> referred to as $hits_index
#  
#  index	blast_data
#  0	Query id
#  1	Subject id
#  2	% identity
#  3	alignment length
#  4	mismatches
#  5	gap openings
#  6	q. start
#  7	q. end
#  8	s. start
#  9	s. end
#  10	e-value
#  11	bit score
#  
#  index	hit_output
#  0	Chromo
#  1	Start
#  2	End
#  3	Strand
#  4	Probe_ID
#  5	Description
#  6	Identity
#  7	alignment length
#  8	mismatch
#  9	gaps
