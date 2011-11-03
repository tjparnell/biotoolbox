package tim_db_helper::bam;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(mean);
use Bio::DB::Sam;


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_bam_scores
	collect_bam_position_scores
	open_bam_db
	sum_total_bam_alignments
);

# Hashes of opened file objects
our %OPENED_BAMFILES; # opened bam file objects
	# in empirical testing, this doesn't really seem to speed things up
	# like I thought it would
	# oh well, keep it anyway????


# The true statement
1; 



### Modules ###



### Collect Bam scores only
sub collect_bam_scores {
	
	# we will collect positioned values but
	# only return the values
	
	# grab the method from the passed arguments
	my $method = $_[3];
	
	# collect the raw data
	my %bam_data = _collect_bam_data(@_);
	
	# combine multiple datapoints at the same position
	my @values;
	if ($method eq 'length') {
		# each hash value is an array of one or more datapoints
		# dump them all into the final values array
		foreach my $position (keys %bam_data) {
			push @values, @{ $bam_data{$position} };
		}
	}
	else {
		# score (coverage) or count
		# each value is a count
		@values = values %bam_data;
	}
	
	# return collected data
	return @values;
	
}




### Collect positioned Bam scores
sub collect_bam_position_scores {
	
	# grab the method from the passed arguments
	my $method = $_[3];
	
	# collect the raw data
	my %bam_data = _collect_bam_data(@_);
	
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
	unless (scalar @_ >= 5) {
		croak " At least five arguments must be passed to collect Bam data!\n";
	}
	my ($region, $region_strand, $stranded, $method, @bam_features) = @_;
		# method can be score, count, or length
	
	# set up hash, either position => count or position => [scores]
	my %bam_data;
	
	# look at each bamfile
	# usually there is only one, but there may be more than one
	foreach my $feature (@bam_features) {
	
		## Get the name of the bam file
		my $bamfile;
		
		if ($feature =~ /^file:(.+)$/) {
			# the passed feature appears to specify a file
			$bamfile = $1;
			
			# check file
			unless (-e $bamfile) {
				croak " Bam file '$bamfile' does not exist!\n";
				return;
			}
		}
		elsif ($feature =~ /^http|ftp/i) {
			# a remote file
			
			# this should be supported by Bio::DB::Sam
			$bamfile = $feature;
		}
		else {
			# otherwise we assume the passed feature is a database object
			
			# get bedfile name
			($bamfile) = $feature->get_tag_values('bamfile');
		}
		croak " no bamfile specified!\n" unless $bamfile;
		
		
		## Open the Bam File
		my $bam;
		if (exists $OPENED_BAMFILES{$bamfile} ) {
			# this file is already opened, use it
			$bam = $OPENED_BAMFILES{$bamfile};
		}
		else {
			# this file has not been opened yet, open it
			$bam = open_bam_db($bamfile) or
				croak " unable to open Bam file '$bamfile'";
			
			# store the opened object for later use
			$OPENED_BAMFILES{$bamfile} = $bam;
		}
			
		
		# Set the code to filter alignments based on strand 
		my $filter;
		if ($stranded eq 'sense' and $region_strand == 1) {
			$filter = sub {
				my $a = shift;
				return $a->strand == 1 ? 1 : 0;
			};
		}
		elsif ($stranded eq 'sense' and $region_strand == -1) {
			$filter = sub {
				my $a = shift;
				return $a->strand == -1 ? 1 : 0;
			};
		}
		elsif ($stranded eq 'antisense' and $region_strand == 1) {
			$filter = sub {
				my $a = shift;
				return $a->strand == -1 ? 1 : 0;
			};
		}
		elsif ($stranded eq 'antisense' and $region_strand == -1) {
			$filter = sub {
				my $a = shift;
				return $a->strand == 1 ? 1 : 0;
			};
		}
		else {
			# no strand requested, take all
			$filter = sub {
				return 1;
			};
		}
		
		
		## Collect the data according to the requested method
		if ($method eq 'score') {
			# collecting scores, or in this case, basepair coverage of 
			# alignments over the requested region
			
			my $coverage;
			if ($stranded eq 'sense' or $stranded eq 'antisense') {
				# Cannot currently collect stranded data with the coverage 
				# method. I will keep the filter in here anyway to future-proof 
				# in case Lincoln ever adds this support (don't hold your 
				# breath!)
				
				($coverage) = $bam->features(
					-type     => 'coverage',
					-seq_id   => $region->seq_id,
					-start    => $region->start,
					-end      => $region->end,
					-filter   => $filter,
				);
			}
			else {
				# no stranded data wanted
				($coverage) = $bam->features(
					-type     => 'coverage',
					-seq_id   => $region->seq_id,
					-start    => $region->start,
					-end      => $region->end,
				);
			}
			
			# convert the coverage data
			# by default, this should return the coverage at 1 bp resolution
			if ($coverage) {
				my @scores = $coverage->coverage;
				for (my $i = $region->start; $i <= $region->end; $i++) {
					$bam_data{$i} += shift @scores;
				}
			}
		}
		
		else {
			# either collecting counts or length
			# working with actual alignments
			
			my @alignments;
			if ($stranded eq 'sense' or $stranded eq 'antisense') {
				
				@alignments = $bam->features(
					-type     => 'match',
					-seq_id   => $region->seq_id,
					-start    => $region->start,
					-end      => $region->end,
					-filter   => $filter,
				);
			}
			else {
				# no stranded data wanted
				@alignments = $bam->features(
					-type     => 'match',
					-seq_id   => $region->seq_id,
					-start    => $region->start,
					-end      => $region->end,
				);
			}
			
			# process the alignments
				# want to avoid those whose midpoint are not technically 
				# within the region of interest
			if ($method eq 'count') {
				foreach my $a (@alignments) {
					# enumerate at the alignment's midpoint
					my $position = int( ( ($a->start + $a->end) / 2) + 0.5);
					
					# check that the midpoint is within the requested region
					if (
						$position >= $region->start and 
						$position <= $region->end
					) {
						$bam_data{$position} += 1;
					}
				}
			}
			elsif ($method eq 'length') {
				foreach my $a (@alignments) {
					# record length at the alignment's midpoint
					my $position = int( ( ($a->start + $a->end) / 2) + 0.5);
					
					# check that the midpoint is within the requested region
					if (
						$position >= $region->start and 
						$position <= $region->end
					) {
						push @{ $bam_data{$position} }, $a->length;
					}
				}
			}
		
		}
	}

	
	# return collected data
	return %bam_data;
}


### Open a bigWig database connection
sub open_bam_db {
	
	my $path = shift;
	$path =~ s/^file://; # clean up file prefix if present
	
	# open the database connection 
	my $db;
	eval {
		$db = Bio::DB::Sam->new(
				-bam         => $path,
				-autoindex   => 1,
		);
	};
	
	if ($db) {
		return $db;
	}
	else {
		carp " ERROR: can't open BAM file '$path'!\n";
		return;
	}
}



### Determine total number of alignments in a bam file
sub sum_total_bam_alignments {
	
	# Passed arguments;
	my $sam_file = shift;
	my $min_mapq = shift || 0; # by default we take all alignments
	my $paired   = shift || 0; # by default we assume all alignments are single-end
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
	
	# Count the number of alignments
	my $total_read_number = 0;
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as 
		# a numeric target identifier
		# we can easily convert this to an actual sequence name
		# we will force the conversion to go one chromosome at a time
		
		# sequence name
		my $seq_id = $sam->target_name($tid);
		
		# process the reads according to single or paired-end
		# paired end alignments
		if ($paired) {
			$sam->fetch($seq_id, 
				sub {
					my $a = shift;
					
					# check paired alignment
					return if $a->unmapped;
					return unless $a->proper_pair;
					return if $a->qual < $min_mapq;
					
					# we're only counting forward reads of a pair 
					return if $a->strand != 1; 
					
					# count this fragment
					$total_read_number++;
				}
			);
		}
		
		# single end alignments
		else {
			$sam->fetch($seq_id, 
				sub {
					my $a = shift;
					
					# check paired alignment
					return if $a->unmapped;
					return if $a->qual < $min_mapq;
					
					# count this fragment
					$total_read_number++;
				}
			);
		}
	}
	
	# done
	return $total_read_number;
}





__END__




=head1 NAME

tim_db_helper::bam

=head1 DESCRIPTION

This module is used to collect the dataset scores from a binary 
bam file (.bam) of alignments. The bam file may be identified in one of 
multiple ways. First, a local file may be specified directly by prefixing 
the file name with "file:", for example "file:/my/path/to/file.bam". 
Second, a remote file may be specifie with a URL, for example 
"http://my.server.com/path/file.bam". Third, the bam file may be referenced 
in the database. Typically, a single feature representing the dataset is 
present across each chromosome. The 
feature should contain an attribute ('bamfile') that references the 
location of the binary file representing the alignments. 
In either case, the file is read using the Bio::DB::Sam module, and 
the values extracted from the region of interest. 

Collected data values may be restricted to strand by specifying the desired 
strandedness, 
depending on the method of data collection. Collecting scores, or basepair 
coverage of alignments over the region of interest, does not currently support 
stranded data collection (as of this writing). However, enumerating 
alignments (count method) and collecting alignment lengths do support 
stranded data collection.

If stranded coverage is desired, the best solution is to split the bam file 
into two files according to alignment strand using the biotoolbox script 
'split_bam_by_strand.pl'. 

Currently, paired-end bam files are treated as single-end files. There are 
some limitations regarding working with paired-end alignments that don't 
work well (search, strand, length, etc). If paired-end alignments are to 
be analyzed, they should be processed into another format (BigWig or BigBed). 
See the biotoolbox scripts 'bam2gff_bed.pl' or 'bam2wig.pl' for solutions.

To speed up the program and avoid repetitive opening and 
closing of the files, the opened bam file object is stored in a global 
hash in case it is needed again.

=head1 USAGE

The module requires Lincoln Stein's Bio::DB::Sam to be installed. 

Load the module at the beginning of your program.

	use tim_db_helper::bam;

It will automatically export the name of the subroutines. 

=over

=item collect_bam_scores

This subroutine will collect only the data values from a binary bam file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed five or more arguments in the following order:
    
    1) The database object representing the genomic region of interest. 
       This should be a Bio::DB::SeqFeature object that supports the 
       start, end, and strand methods. Alternatively, a bam file may 
       also be directly specified, prefixed with "file:", "http://", or 
       "ftp://".
    2) The strand of the original feature (or region), -1, 0, or 1.
    3) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       "none" or "no". Only those scores which match the indicated 
       strandedness are collected.
    4) The method or type of data collected. 
       Acceptable values include 'score' (returns the basepair coverage
       of alignments over the region of interest), 'count' (returns the 
       number of alignments found at each base position in the region, 
       recorded at the alignment's midpoint), or 'length' (returns the 
       mean lengths of the alignments found at each base position in 
       the region, recorded at the alignment's midpoint). 
    5) One or more database feature objects that contain the reference 
       to the .bam file. They should contain the attribute 'bamfile' 
       which has the path to the Bam file. Alternatively, pass one 
       or more filenames of .bam files. Each filename should be 
       prefixed with 'file:' to indicate that it is a direct file 
       reference, and not a database object.

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_bam_position_scores

This subroutine will collect the score values from a binary bam file 
for the specified database region keyed by position. 

The subroutine is passed five or more arguments in the following order:
    
    1) The database object representing the genomic region of interest. 
       This should be a Bio::DB::SeqFeature object that supports the 
       start, end, and strand methods. Alternatively, a bam file may 
       also be directly specified, prefixed with "file:", "http://", or 
       "ftp://".
    2) The strand of the original feature (or region), -1, 0, or 1.
    3) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       "none" or "no". Only those scores which match the indicated 
       strandedness are collected.
    4) The method or type of data collected. 
       Acceptable values include 'score' (returns the basepair coverage
       of alignments over the region of interest), 'count' (returns the 
       number of alignments found at each base position in the region, 
       recorded at the alignment's midpoint), or 'length' (returns the 
       mean lengths of the alignments found at each base position in 
       the region, recorded at the alignment's midpoint). 
    5) One or more database feature objects that contain the reference 
       to the .bam file. They should contain the attribute 'bamfile' 
       which has the path to the Bam file. Alternatively, pass one 
       or more filenames of .bam files. Each filename should be 
       prefixed with 'file:' to indicate that it is a direct file 
       reference, and not a database object.

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. The feature midpoint is used 
as the key position. When multiple features are found at the same 
position, a simple mean (for length data methods) or sum 
(for count methods) is returned.

=item open_bam_db()

This subroutine will open a Bam database connection. Pass either the 
local path to a Bam file (.bam extension) or the URL of a remote Bam 
file. It will return the opened database object.

=item sum_total_bam_alignments()

This subroutine will sum the total number of properly mapped alignments 
in a bam file. Pass the subroutine one to three arguments. 
    
    1) The name of the Bam file which should be counted. Alternatively,  
       an opened Bio::DB::Sam object may also be given. Required.
    2) Optionally pass the minimum mapping quality of the reads to be 
       counted. The default is 0, or all alignments are counted.
    3) Optionally pass a boolean value (1 or 0) indicating whether 
       the Bam file represents paired-end alignments. Only proper 
       alignment pairs are counted. The default is to treat all 
       alignments as single-end.
       
The subrouting will return the number of alignments.

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

