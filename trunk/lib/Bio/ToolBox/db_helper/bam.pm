package Bio::ToolBox::db_helper::bam;

# modules
require Exporter;
use strict;
use Carp;
use File::Copy;
use Statistics::Lite qw(mean);
use Bio::DB::Sam;
our $parallel;
eval {
	# check for parallel support, when counting bam alignments
	require Parallel::ForkManager;
	$parallel = 1;
};
our $VERSION = '1.15';

# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	open_bam_db
	check_bam_index
	write_new_bam_file
	collect_bam_scores
	collect_bam_position_scores
	sum_total_bam_alignments
);

# Hash of Bigfile chromosomes
our %BAM_CHROMOS;
	# sometimes user may request a chromosome that's not in the bigfile
	# that could lead to an exception
	# we will record the chromosomes list in this hash
	# $BAM_CHROMOS{bigfile}{chromos}

# Opened Bam db objects
our %OPENED_BAM;
	# a cache for opened Bam files, primarily for collecting scores

# The true statement
1; 



### Modules ###



### Open a bigWig database connection
sub open_bam_db {
	my $bamfile = shift;
	
	# check the path
	my $path = $bamfile;
	$path =~ s/^file://; # strip the file prefix if present
	
	# check for bam index
	check_bam_index($path);
	
	# open the bam database object
	my $sam;
	eval {
		$sam = Bio::DB::Sam->new(-bam => $path);
	};
	return unless $sam;
	
	# done
	return $sam;
}



sub check_bam_index {
	# I find that relying on -autoindex yields a flaky Bio::DB::Sam object that 
	# doesn't always work as expected. Best to create the index BEFORE opening
	
	my $bamfile = shift;
	return if ($bamfile =~ /^(?:http|ftp)/i); # I can't do much with remote files
	
	# we will check the modification time to make sure index is newer
	my $bam_mtime = (stat($bamfile))[9];
	
	# optional index names
	my $bam_index = "$bamfile.bai"; # .bam.bai
	my $alt_index = $bamfile;
	$alt_index =~ s/bam$/bai/i; # picard uses .bai instead of .bam.bai as samtools does
	
	# check for existing index
	if (-e $bam_index) {
		if ( (stat($bam_index))[9] < $bam_mtime) {
			# index is older than bam file !? re-index
			Bio::DB::Bam->index_build($bamfile);
		}
	}
	elsif (-e $alt_index) {
		if ( (stat($alt_index))[9] < $bam_mtime) {
			# index is older than bam file !? re-index
			Bio::DB::Bam->index_build($bamfile);
		}
		else {
			# reuse this index
			copy($alt_index, $bam_index)
		}
	}
	else {
		# make a new index
		Bio::DB::Bam->index_build($bamfile);
	}
}


### Write a new bam file
sub write_new_bam_file {
	my $file = shift;
	$file .= '.bam' unless $file =~ /\.bam$/i;
	my $bam = Bio::DB::Bam->open($file, 'w');
	carp "unable to open bam file $file!\n" unless $bam;
	return $bam;
}



### Collect Bam scores only
sub collect_bam_scores {
	
	# set the do_index boolean to false
	# return the scores
	return _collect_bam_data(0, @_);
}




### Collect positioned Bam scores
sub collect_bam_position_scores {
	
	# collect the raw data
	# set the do_index boolean to true
	my %bam_data = _collect_bam_data(1, @_);
	
	# grab the method from the passed arguments
	my $method = $_[3];
	
	# combine multiple datapoints at the same position
	if ($method eq 'length') {
		# each value is an array of one or more datapoints
		# we will take the simple mean
		foreach my $position (keys %bam_data) {
			$bam_data{$position} = mean( @{$bam_data{$position}} );
		}
	}
	
	# return collected data
	return %bam_data;
}




### Actual collection of scores
sub _collect_bam_data {
	
	# pass the required information
	unless (scalar @_ >= 8) {
		confess " At least eight arguments must be passed to collect Bam data!\n";
	}
	my (
		$do_index, 
		$chromo,
		$start,
		$stop,
		$strand, 
		$stranded, 
		$value_type, 
		@bam_features
	) = @_;
		# method can be score, count, or length
	
	# initialize score structures
	# which one is used depends on the $do_index boolean variable
	my %pos2data; # either position => count or position => [scores]
	my @scores; # just scores
	
	# look at each bamfile
	# usually there is only one, but there may be more than one
	foreach my $bamfile (@bam_features) {
	
		## Open the Bam File
		my $bam;
		if (exists $OPENED_BAM{$bamfile}) {
			# re-use the cached bam object
			$bam = $OPENED_BAM{$bamfile};
		}
		else {
			# open and cache the bam file
			$bam = open_bam_db($bamfile);
			$OPENED_BAM{$bamfile} = $bam;
			%{ $BAM_CHROMOS{$bamfile} } = map { $_ => 1 } $bam->seq_ids;
		}
			
		# first check that the chromosome is present
		unless (exists $BAM_CHROMOS{$bamfile}{$chromo}) {
			next;
		}
		
		# convert coordinates into low level coordinates
		# consumed by the low level Bam API
		my ($tid, $zstart, $end) = 
			$bam->header->parse_region("$chromo:$start\-$stop");
	
		
		## Collect the data according to the requested value type
		# we will either use simple coverage (score method) or
		# process the actual alignments (count or length)
		
		## Coverage
		if ($value_type eq 'score') {
			# collecting scores, or in this case, basepair coverage of 
			# alignments over the requested region
			
			# generate the coverage, this will ignore strand
			my $coverage = $bam->bam_index->coverage(
				$bam->bam,
				$tid,
				$zstart, # 0-based coordinates
				$end,
			);
			
			# convert the coverage data
			# by default, this should return the coverage at 1 bp resolution
			if (scalar @$coverage) {
				
				# check whether we need to index the scores
				if ($do_index) {
					for (my $i = $start; $i <= $stop; $i++) {
						# move the scores into the position score hash
						$pos2data{$i} += $coverage->[ $i - $start ];
					}
				}
				else {
					@scores = @$coverage;
				}
			}
		}
		
		
		## Alignments
		else {
			# either collecting counts or length
			# working with actual alignments
			
			## Set the callback and a callback data structure
			my $callback = _assign_callback($stranded, $strand, $value_type, $do_index);
			my %data = (
				'scores' => \@scores,
				'index'  => \%pos2data,
				'start'  => $start,
				'stop'   => $stop,
			);
			
			# get the alignments
			# we are using the low level API to eke out performance
			$bam->bam_index->fetch($bam->bam, $tid, $zstart, $end, $callback, \%data);
			
		}
	}

	
	## Return collected data
	if ($do_index) {
		return %pos2data;
	}
	else {
		return @scores;
	}
}


### Determine total number of alignments in a bam file
sub sum_total_bam_alignments {
	
	# Passed arguments;
	my $sam_file = shift;
	my $min_mapq = shift || 0; # by default we take all alignments
	my $paired   = shift || 0; # by default we assume all alignments are single-end
	my $cpu      = shift || 2; # number of forks to execute in parallel
	$cpu = 1 unless ($parallel);
	unless ($sam_file) {
		carp " no Bam file or bam db object passed!\n";
		return;
	}
	
	
	# Open Bam file if necessary
	my $sam;
	my $sam_ref = ref $sam_file;
	if ($sam_ref =~ /Bio::DB::Sam/) {
		# we have an opened sam db object
		$sam = $sam_file;
	}
	else {
		# we have a name of a sam file
		$sam = open_bam_db($sam_file);
		return unless ($sam);
	}
	
	# prepare the counting subroutine
	my $counter = sub {
		my $tid = shift;
		my $number = 0;
		
		# process the reads according to single or paired-end
		# paired end alignments
		if ($paired) {
			$sam->bam_index->fetch(
				$sam->bam, 
				$tid, 
				0, 
				$sam->target_len($tid), 
				sub {
					my ($a, $number) = @_;
					
					# check paired alignment
					return unless $a->proper_pair;
					return if $a->reversed; # only count left alignments
					return if $a->qual < $min_mapq;
					
					# count this fragment
					$$number++;
				}, 
				\$number
			);
		}
		
		# single end alignments
		else {
			$sam->bam_index->fetch(
				$sam->bam, 
				$tid, 
				0, 
				$sam->target_len($tid), 
				sub {
					my ($a, $number) = @_;
					
					# check alignment
					return if $a->unmapped;
					return if $a->qual < $min_mapq;
					
					# count this fragment
					$$number++;
				}, 
				\$number
			);
		}
		return $number;
	};
	
	# Count the alignments on each chromosome
	my $total_read_number = 0;
	if ($cpu > 1) {
		# count each chromosome in multiple parallel threads to speed things up 
		
		# generate relatively equal lists of chromosome for each process based on length
		my @chromosomes = map { $_->[0] }
			sort { $b->[1] <=> $a->[1] }
			map { [$_, $sam->target_len($_)] } 
			(0 .. $sam->n_targets - 1);
		my @list; # array of arrays, [process][chromosome id]
		my $i = 1;
		while (@chromosomes) {
			push @{ $list[$i] }, shift @chromosomes;
			$i++;
			$i = 1 if $i > $cpu;
		}
		
		# we will use Parallel ForkManager for convenience
		my $pm = Parallel::ForkManager->new($cpu);
		$pm->run_on_finish( sub {
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $count) = @_;
			$total_read_number += $$count;
		});
		
		# Count the chromosomes in parallel processess
		foreach my $n (1 .. $cpu) {
			$pm->start and next;
	
			### In child
			$sam->clone; # to make it fork safe
			my $count = 0;
			foreach ( @{$list[$n]} ) {
				# count each chromosome in this process list
				$count += &{$counter}($_);
			}
			$pm->finish(0, \$count); 
		}
		$pm->wait_all_children;
	}
	
	else {
		# loop through all the chromosomes in one execution thread
		for my $tid (0 .. $sam->n_targets - 1) {
			# each chromosome is internally represented in the bam file as 
			# a numeric target identifier
			$total_read_number += &{$counter}($tid);
		}
	}
	
	# done
	return $total_read_number;
}


### Generate callback subroutine for walking through Bam alignments
sub _assign_callback {
	# generate the callback code depending on whether we want to look at 
	# stranded data, collecting counts or length, or whether indexed data
	# is wanted.
	
	# we perform a check of whether the alignment midpoint is within the 
	# search region for indexed data only
	
	# these subroutines are designed to work with the low level fetch API
	
	# there are so many different subroutines because I want to increase 
	# efficiency by limiting the number of conditional tests in one generic subroutine
	
	my ($stranded, $strand, $value_type, $do_index) = @_;
	
	# all alignments
	if (
		$stranded eq 'all' and 
		$value_type eq 'count' and 
		$do_index
	) {
		return \&_all_count_indexed;
	}
	elsif (
		$stranded eq 'all' and 
		$value_type eq 'count' and 
		!$do_index
	) {
		return \&_all_count_array;
	}
	elsif (
		$stranded eq 'all' and 
		$value_type eq 'length' and 
		$do_index
	) {
		return \&_all_length_indexed;
	}
	elsif (
		$stranded eq 'all' and 
		$value_type eq 'length' and 
		!$do_index
	) {
		return \&_all_length_array;
	}
	
	
	# sense, forward strand 
	elsif (
		$stranded eq 'sense' and 
		$strand == 1 and 
		$value_type eq 'count' and 
		$do_index
	) {
		return \&_sense_forward_count_indexed;
	}
	elsif (
		$stranded eq 'sense' and 
		$strand == 1 and 
		$value_type eq 'count' and 
		!$do_index
	) {
		return \&_sense_forward_count_array;
	}
	elsif (
		$stranded eq 'sense' and 
		$strand == 1 and 
		$value_type eq 'length' and 
		$do_index
	) {
		return \&_sense_forward_length_indexed;
	}
	elsif (
		$stranded eq 'sense' and 
		$strand == 1 and 
		$value_type eq 'length' and 
		!$do_index
	) {
		return \&_sense_forward_length_array;
	}
	
	
	# sense, reverse strand
	elsif (
		$stranded eq 'sense' and 
		$strand == -1 and 
		$value_type eq 'count' and 
		$do_index
	) {
		return \&_sense_reverse_count_indexed;
	}
	elsif (
		$stranded eq 'sense' and 
		$strand == -1 and 
		$value_type eq 'count' and 
		!$do_index
	) {
		return \&_sense_reverse_count_array;
	}
	elsif (
		$stranded eq 'sense' and 
		$strand == -1 and 
		$value_type eq 'length' and 
		$do_index
	) {
		return \&_sense_reverse_length_indexed;
	}
	elsif (
		$stranded eq 'sense' and 
		$strand == -1 and 
		$value_type eq 'length' and 
		!$do_index
	) {
		return \&_sense_reverse_length_array;
	}
	
	
	# anti-sense, forward strand 
	if (
		$stranded eq 'antisense' and 
		$strand == 1 and 
		$value_type eq 'count' and 
		$do_index
	) {
		return \&_antisense_forward_count_indexed;
	}
	elsif (
		$stranded eq 'antisense' and 
		$strand == 1 and 
		$value_type eq 'count' and 
		!$do_index
	) {
		return \&_antisense_forward_count_array;
	}
	elsif (
		$stranded eq 'antisense' and 
		$strand == 1 and 
		$value_type eq 'length' and 
		$do_index
	) {
		return \&_antisense_forward_length_indexed;
	}
	elsif (
		$stranded eq 'antisense' and 
		$strand == 1 and 
		$value_type eq 'length' and 
		!$do_index
	) {
		return \&_antisense_forward_length_array;
	}
	
	
	# anti-sense, reverse strand
	elsif (
		$stranded eq 'antisense' and 
		$strand == -1 and 
		$value_type eq 'count' and 
		$do_index
	) {
		return \&_antisense_reverse_count_indexed;
	}
	elsif (
		$stranded eq 'antisense' and 
		$strand == -1 and 
		$value_type eq 'count' and 
		!$do_index
	) {
		return \&_antisense_reverse_count_array;
	}
	elsif (
		$stranded eq 'antisense' and 
		$strand == -1 and 
		$value_type eq 'length' and 
		$do_index
	) {
		return \&_antisense_reverse_length_indexed;
	}
	elsif (
		$stranded eq 'antisense' and 
		$strand == -1 and 
		$value_type eq 'length' and 
		!$do_index
	) {
		return \&_antisense_reverse_length_array ;
	}
	
	
	# I goofed
	else {
		confess("Programmer error: stranded $stranded, strand $strand, value_type ". 
				"$value_type, index $do_index\n");
	}
}


#### Callback subroutines 
# the following are all of the callback subroutines 

sub _all_count_indexed {
	my ($a, $data) = @_;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++ if 
		( $pos >= $data->{'start'} and $pos <= $data->{'stop'} );
}

sub _all_count_array {
	my ($a, $data) = @_;
	push @{ $data->{'scores'} }, 1;
}

sub _all_length_indexed {
	my ($a, $data) = @_;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	if ( $pos >= $data->{'start'} and $pos <= $data->{'stop'} ) {
		push @{ $data->{'index'}{$pos} }, ($a->calend - $a->pos);
	}
}

sub _all_length_array {
	my ($a, $data) = @_;
	push @{ $data->{'scores'} }, ($a->calend - $a->pos);
}

sub _sense_forward_count_indexed {
	my ($a, $data) = @_;
	return if $a->reversed;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++ if 
		( $pos >= $data->{'start'} and $pos <= $data->{'stop'} );
}

sub _sense_forward_count_array {
	my ($a, $data) = @_;
	return if $a->reversed;
	push @{ $data->{'scores'} }, 1;
}

sub _sense_forward_length_indexed {
	my ($a, $data) = @_;
	return if $a->reversed;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	if ( $pos >= $data->{'start'} and $pos <= $data->{'stop'} ) {
		push @{ $data->{'index'}{$pos} }, ($a->calend - $a->pos);
	}
}

sub _sense_forward_length_array {
	my ($a, $data) = @_;
	return if $a->reversed;
	push @{ $data->{'scores'} }, ($a->calend - $a->pos);
}

sub _sense_reverse_count_indexed {
	my ($a, $data) = @_;
	return unless $a->reversed;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++ if 
		( $pos >= $data->{'start'} and $pos <= $data->{'stop'} );
}

sub _sense_reverse_count_array {
	my ($a, $data) = @_;
	return unless $a->reversed;
	push @{ $data->{'scores'} }, 1;
}

sub _sense_reverse_length_indexed {
	my ($a, $data) = @_;
	return unless $a->reversed;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	if ( $pos >= $data->{'start'} and $pos <= $data->{'stop'} ) {
		push @{ $data->{'index'}{$pos} }, ($a->calend - $a->pos);
	}
}

sub _sense_reverse_length_array {
	my ($a, $data) = @_;
	return unless $a->reversed;
	push @{ $data->{'scores'} }, ($a->calend - $a->pos);
}

sub _antisense_forward_count_indexed {
	my ($a, $data) = @_;
	return unless $a->reversed;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++ if 
		( $pos >= $data->{'start'} and $pos <= $data->{'stop'} );
}

sub _antisense_forward_count_array {
	my ($a, $data) = @_;
	return unless $a->reversed;
	push @{ $data->{'scores'} }, 1;
}

sub _antisense_forward_length_indexed {
	my ($a, $data) = @_;
	return unless $a->reversed;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	if ( $pos >= $data->{'start'} and $pos <= $data->{'stop'} ) {
		push @{ $data->{'index'}{$pos} }, ($a->calend - $a->pos);
	}
}

sub _antisense_forward_length_array {
	my ($a, $data) = @_;
	return unless $a->reversed;
	push @{ $data->{'scores'} }, ($a->calend - $a->pos);
}

sub _antisense_reverse_count_indexed {
	my ($a, $data) = @_;
	return if $a->reversed;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++ if 
		( $pos >= $data->{'start'} and $pos <= $data->{'stop'} );
}

sub _antisense_reverse_count_array {
	my ($a, $data) = @_;
	return if $a->reversed;
	push @{ $data->{'scores'} }, 1;
}

sub _antisense_reverse_length_indexed {
	my ($a, $data) = @_;
	return if $a->reversed;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	if ( $pos >= $data->{'start'} and $pos <= $data->{'stop'} ) {
		push @{ $data->{'index'}{$pos} }, ($a->calend - $a->pos);
	}
}

sub _antisense_reverse_length_array {
	my ($a, $data) = @_;
	return if $a->reversed;
	push @{ $data->{'scores'} }, ($a->calend - $a->pos);
}



__END__

=head1 NAME

Bio::ToolBox::db_helper::bam

=head1 DESCRIPTION

This module is used to collect the dataset scores from a binary 
bam file (.bam) of alignments. Bam files may be local or remote, 
and are usually prefixed with 'file:', 'http://', of 'ftp://'.

Collected data values may be restricted to strand by specifying the desired 
strandedness (sense, antisense, or all), 
depending on the method of data collection. Collecting scores, or basepair 
coverage of alignments over the region of interest, does not currently support 
stranded data collection (as of this writing). However, enumerating 
alignments (count method) and collecting alignment lengths do support 
stranded data collection. Alignments are checked to see whether their midpoint 
is within the search interval before counting or length collected. 

Currently, paired-end bam files are treated as single-end files. There are 
some limitations regarding working with paired-end alignments that don't 
work well (search, strand, length, etc). If paired-end alignments are to 
be analyzed, they should be processed into another format (BigWig or BigBed). 
See the biotoolbox scripts 'bam2gff_bed.pl' or 'bam2wig.pl' for solutions.

=head1 USAGE

The module requires Lincoln Stein's Bio::DB::Sam to be installed. 

Load the module at the beginning of your program.

	use Bio::ToolBox::db_helper::bam;

It will automatically export the name of the subroutines. 

=over

=item open_bam_db()

This subroutine will open a Bam database connection. Pass either the 
local path to a Bam file (.bam extension) or the URL of a remote Bam 
file. A remote bam file must be indexed. A local bam file may be 
automatically indexed upon opening if the user has write permissions 
in the parent directory. 

It will return the opened database object.

=item check_bam_index()

This subroutine will check whether a bam index file is present and, 
if not, generate one. The Bio::DB::Sam module uses the samtools 
style index extension, F<.bam.bai>, as opposed to the picard style 
extension, F<.bai>. If a F<.bai> index is present, it will copy the 
file as F<.bam.bai> index. Unfortunately, a F<.bai> index cannot be 
used directly.

This method is called automatically prior to opening a bam database. 

=item write_new_bam_file()

This subroutine will open a new empty Bam file. Pass the name of the 
new file as the argument. It will return a Bio::DB::Bam object to 
which you can write a header followed by alignments. Be sure you know 
what to do before using this method! 

=item collect_bam_scores

This subroutine will collect only the data values from a binary bam file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed seven or more arguments in the following order:
    
    1) The chromosome or seq_id
    2) The start position of the segment to collect 
    3) The stop or end position of the segment to collect 
    4) The strand of the original feature (or region), -1, 0, or 1.
    5) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       or "all". Only those scores which match the indicated 
       strandedness are collected.
    6) The type of data collected. 
       Acceptable values include 'score' (returns the basepair coverage
       of alignments over the region of interest), 'count' (returns the 
       number of alignments found at each base position in the region, 
       recorded at the alignment's midpoint), or 'length' (returns the 
       mean lengths of the alignments found at each base position in 
       the region, recorded at the alignment's midpoint). 
    7) The paths, either local or remote, to one or more Bam files.

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_bam_position_scores

This subroutine will collect the score values from a binary bam file 
for the specified database region keyed by position. 

The subroutine is passed the same arguments as collect_bam_scores().

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. The feature midpoint is used 
as the key position. When multiple features are found at the same 
position, a simple mean (for length data methods) or sum 
(for count methods) is returned.

=item sum_total_bam_alignments()

This subroutine will sum the total number of properly mapped alignments 
in a bam file. Pass the subroutine one to four arguments. 
    
    1) The name of the Bam file which should be counted. Alternatively,  
       an opened Bio::DB::Sam object may also be given. Required.
    2) Optionally pass the minimum mapping quality of the reads to be 
       counted. The default is 0, where all alignments are counted.
    3) Optionally pass a boolean value (1 or 0) indicating whether 
       the Bam file represents paired-end alignments. Only proper 
       alignment pairs are counted. The default is to treat all 
       alignments as single-end.
    4) Optionally pass the number of parallel processes to execute 
       when counting alignments. Walking through a Bam file is 
       time consuming but can be easily parallelized. The module 
       Parallel::ForkManager is required, and the default is a 
       conservative two processes when it is installed.
       
The subroutine will return the number of alignments.

=back

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

