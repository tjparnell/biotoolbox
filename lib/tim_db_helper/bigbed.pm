package tim_db_helper::bigbed;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(mean);
use Bio::DB::BigBed;


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_bigbed_scores
	collect_bigbed_position_scores
	bed_to_bigbed_conversion
);

# Hashes of opened file objects
our %OPENED_BEDFILES; # opened bigbed file objects
	# in empirical testing, this doesn't really seem to speed things up
	# like I thought it would
	# oh well, keep it anyway????


# The true statement
1; 



### Modules ###



### Collect BigBed scores only
sub collect_bigbed_scores {
	
	# we will collect positioned values but
	# only return the values
	
	# grab the method from the passed arguments
	my $method = $_[3];
	
	# collect the raw data
	my %bed_data = _collect_bigbed_data(@_);
	
	# combine multiple datapoints at the same position
	my @values;
	if ($method eq 'score' or $method eq 'length') {
		# each value is an array of one or more datapoints
		foreach my $position (keys %bed_data) {
			push @values, mean( @{ $bed_data{$position} } );
		}
	}
	elsif ($method eq 'count') {
		# each value is a count
		@values = values %bed_data;
	}
	
	# return collected data
	return @values;
	
}




### Collect positioned BigBed scores
sub collect_bigbed_position_scores {
	
	# grab the method from the passed arguments
	my $method = $_[3];
	
	# collect the raw data
	my %bed_data = _collect_bigbed_data(@_);
	
	# combine multiple datapoints at the same position
	if ($method eq 'score' or $method eq 'length') {
		# each value is an array of one or more datapoints
		# we will take the simple mean
		foreach my $position (keys %bed_data) {
			$bed_data{$position} = mean( @{$bed_data{$position}} );
		}
	}
	
	# return collected data
	return %bed_data;
}



### Actual collection of scores
sub _collect_bigbed_data {
	
	# pass the required information
	unless (scalar @_ >= 5) {
		croak " At least five arguments must be passed to collect BigBed data!\n";
	}
	my ($region, $region_strand, $stranded, $method, @bed_features) = @_;
		# method can be score, count, or length
	
	# set up hash, either position => count or position => [scores]
	my %bed_data;
	
	# look at each bedfile
	# usually there is only one, but for stranded data there may be 
	# two bedfiles (+ and -), so we'll check each wig file for strand info
	foreach my $feature (@bed_features) {
	
		# Get the name of the bigbed file
		my $bedfile;
		
		if ($feature =~ /^file:(.+)$/) {
			# the passed feature appears to specify a file
			$bedfile = $1;
			
			# check the file
			unless (-e $bedfile) {
				croak " BigBed file '$bedfile' does not exist!\n";
				return;
			}
		}
		elsif ($feature =~ /^http|ftp/i) {
			# a remote file
			
			# this should be supported by Bio::DB::BigBed
			$bedfile = $feature;
		}
		else {
			# otherwise we assume the passed feature is a database object
			
			# get bedfile name
			my ($bedfile) = $feature->get_tag_values('bigbedfile');
		}
		
		# check for opened bedfile
		my $bb;
		if (exists $OPENED_BEDFILES{$bedfile} ) {
			# this file is already opened, use it
			$bb = $OPENED_BEDFILES{$bedfile};
		}
		else {
			# this file has not been opened yet, open it
			$bb = Bio::DB::BigBed->new($bedfile);
			unless ($bb) {
				croak " unable to open data BigBed file '$bedfile'";
			}
			
			# store the opened object for later use
			$OPENED_BEDFILES{$bedfile} = $bb;
		}
			
		# collect the features overlapping the region
		my $bb_stream = $bb->features(
			-seq_id   => $region->seq_id, 
			-start    => $region->start, 
			-end      => $region->end,
			-iterator => 1,
		);
		
		# process each feature
		while (my $bed = $bb_stream->next_seq) {
			
			# First check whether the strand is acceptable
			if (
				$stranded eq 'all' # all data is requested
				or $bed->strand == 0 # unstranded data
				or ( 
					# sense data
					$region_strand == $bed->strand 
					and $stranded eq 'sense'
				) 
				or (
					# antisense data
					$region_strand != $bed->strand  
					and $stranded eq 'antisense'
				)
			) {
				# we have acceptable data to collect
			
				# determine position to record
				my $position;
				if ($bed->start == $bed->end) {
					# just one position recorded
					$position = $bed->start;
				}
				else {
					# calculate the midpoint
					$position = int( 
						( ($bed->start + $bed->end) / 2) + 0.5
					);
				}
				
				# store the appropriate datapoint
				# for score and length, we're putting these into an array
				if ($method eq 'score') {
					# perform addition to force the score to be a scalar value
					push @{ $bed_data{$position} }, $bed->score + 0;
				}
				elsif ($method eq 'count') {
					$bed_data{$position} += 1;
				}
				elsif ($method eq 'length') {
					# I hope that length is supported, but not sure
					# may have to calculate myself
					push @{ $bed_data{$position} }, $bed->length;
				}
			}
		}
	}

	
	# return collected data
	return %bed_data;
}



### Bed to BigBed file conversion
sub bed_to_bigbed_conversion {
	
	# Collect passed arguments
	my $argument_ref = shift;
	unless ($argument_ref) {
		carp "no arguments passed!";
		return;
	}
	
	# wigfile
	my $bedfile = $argument_ref->{'bed'} || undef;
	unless ($bedfile) {
		carp "no bed file passed!";
		return;
	}
	
	# identify bigbed conversion utility
	my $bb_app_path = $argument_ref->{'bbapppath'} || undef;
	unless ($bb_app_path) {
		carp " bedToBigBed converter utility not specified; unable to proceed!\n";
		return;
	}
	
	# Generate list of chromosome sizes if necessary
	my $chromo_file = $argument_ref->{'chromo'} || undef;
	unless ($chromo_file) {
		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		print " generating chromosome file....\n";
		
		# check that we have tim_db_helper loaded
		my $db = $argument_ref->{'db'} || undef;
		unless ($db) {
			carp " database or chromosome file not specified! " . 
				"Unable to convert!\n";
			return;
		};
		
		# determine reference sequence type
		my $ref_seq_type = $argument_ref->{'seq_type'} || 'chromosome';
			# the annotation gff may have the reference sequences labeled
			# as various types, such as chromosome, sequence, 
			# contig, scaffold, etc
			# this is set in the configuration file
			# this could pose problems if more than one is present
		
		# generate chromosome lengths file
		my @chromos = $db->features(-type => $ref_seq_type);
		unless (@chromos) {
			die " no '$ref_seq_type' features identified in database!\n";
		}
		open CHR_FILE, ">tim_helper_chr_lengths.txt";
		foreach (@chromos) {
			print CHR_FILE $_->name, "\t", $_->length, "\n";
		}
		close CHR_FILE;
		$chromo_file = "tim_helper_chr_lengths.txt";
	}
	
	# generate the bw file name
	my $bb_file = $bedfile;
	$bb_file =~ s/\.bed$/.bb/;
	
	# generate the bigbed file using Jim Kent's utility
	# Bio::DB::BigFile does not support BigBed conversion, unlike BigWig files 
	# execute
	print " converting $bedfile to BigBed....\n";
	system $bb_app_path, $bedfile, $chromo_file, $bb_file;
	
	# check the result
	if (-e $bb_file) {
		# conversion successful
		if ($chromo_file eq 'tim_helper_chr_lengths.txt') {
			# we no longer need our temp chromosome file
			unlink $chromo_file;
		}
		return $bb_file;
	}
	else {
		print " Conversion failed. You should try manually and watch for errors\n";
		if ($chromo_file eq 'tim_helper_chr_lengths.txt') {
			# leave the temp chromosome file as a courtesy
			print " Leaving temporary chromosome file '$chromo_file'\n";
		}
		return;
	}
}




__END__




=head1 NAME

tim_db_helper::bigbed

=head1 DESCRIPTION

This module supports the use of bigbed file in the biotoolbox scripts, both 
in the collection of data from a bigbed file, as well as the generation of 
bigbed files.

=head2 Data collection

This module is used to collect the dataset scores from a binary 
bigbed file (.bb). The file may be identified in one of two ways. First,
it may be referenced in the database. Typically, a single 
feature representing the dataset is present across each chromosome. The 
feature should contain an attribute ('bigbedfile') that references the 
location of the binary file representing the dataset scores. Second, 
the local location of the file may be directly passed to the subroutine. 

In either case, the file is read using the Bio::DB::BigBed module, and 
the values extracted from the region of interest. 

Scores may be restricted to strand by specifying the desired strandedness. 
For example, to collect transcription data over a gene, pass the strandedness 
value 'sense'. If the strand of the region database object (representing the 
gene) matches the strand of the bed feature, then the data for that bed 
feature is collected.  

For loading bigbed files into a Bio::DB database, see the biotoolbox perl 
script 'big_filegff3.pl'.

To speed up the program and avoid repetitive opening and 
closing of the files, the opened bigbed file object is stored in a global 
hash in case it is needed again.

=head2 File generation

This module also supports the generation of bigbed files. This is dependent 
on Jim Kent's UCSC commandline utility. It automates the collection of 
chromosome information in prepration of conversion, if necessary.

=head1 USAGE

The module requires Lincoln Stein's Bio::DB::BigBed to be installed. 

Load the module at the beginning of your program.

	use tim_db_helper::bigbed;

It will automatically export the name of the subroutines. 

=over

=item collect_bigbed_scores

This subroutine will collect only the data values from a binary bigbed file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed five or more arguments in the following order:
    
    1) The database object representing the genomic region of interest. 
       This should be a Bio::DB::SeqFeature object that supports the 
       start, end, and strand methods.
    2) The strand of the original feature (or region), -1, 0, or 1.
    3) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       "none" or "no". Only those scores which match the indicated 
       strandedness are collected.
    4) The method or type of data collected. 
       Acceptable values include 'score' (returns the bed feature 
       score), 'count' (returns the number of bed features found), or 
       'length' (returns the length of the bed features found). 
    5) One or more database feature objects that contain the reference 
       to the .bb file. They should contain the attribute 'bigbedfile' 
       which has the path to the BigBed file. Alternatively, pass one 
       or more filenames of .bb files. Each filename should be 
       prefixed with 'file:' to indicate that it is a direct file 
       reference, and not a database object.

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_bigbed_position_scores

This subroutine will collect the score values from a binary wig file 
for the specified database region keyed by position. 

The subroutine is passed five or more arguments in the following order:
    
    1) The database object representing the genomic region of interest. 
       This should be a Bio::DB::SeqFeature object that supports the 
       start, end, and strand methods.
    2) The strand of the original feature (or region), -1, 0, or 1.
    3) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       "none" or "no". Only those scores which match the indicated 
       strandedness are collected.
    4) The method or type of data collected. 
       Acceptable values include 'score' (returns the bed feature 
       score), 'count' (returns the number of bed features found), or 
       'length' (returns the length of the bed features found). 
    5) One or more database feature objects that contain the reference 
       to the .bb file. They should contain the attribute 'bigbedfile' 
       which has the path to the BigBed file. Alternatively, pass one 
       or more filenames of .bb files. Each filename should be 
       prefixed with 'file:' to indicate that it is a direct file 
       reference, and not a database object.

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. The feature midpoint is used 
as the key position. When multiple features are found at the same 
position, a simple mean (for score or length data methods) or sum 
(for count methods) is returned.

=item bed_to_bigbed_conversion

This subroutine will convert a bed file to a bigBed file. See the UCSC 
documentation regarding bed (http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED)
and bigBed (http://genome.ucsc.edu/goldenPath/help/bigBed.html) file formats. 
It uses Jim Kent's bedToBigBed utility to perform the conversion. This 
must be present on the system for the conversion to succeed. 

The conversion requires a list of chromosome name and sizes in a simple text 
file, where each line is comprised of two columns, "<chromosome name> 
<size in bases>". This file may be specified, or automatically generated if 
given a Bio::DB database name (preferred to ensure genome version 
compatibility).

The function returns the name of the bigBed file, which will be the 
input bed file basename with the extension ".bb". Note that the it does 
not check for success of writing the bigbed file. Check STDERR for errors 
in bigbed file generation.

Pass the function an anonymous hash of arguments, including the following:

  Required:
  bed         => The name of the bed source file. 
  db          => Provide an opened database object from which to generate 
                 the chromosome sizes information.
  Optional: 
  seq_type    => The GFF type of the reference sequence in the database. 
                 This is typically "chromosome", but could also be 
                 "sequence", "scaffold", "contig", etc. The default 
                 value is "chromosome". This value may be provided by 
                 checking the entry in tim_db_helper.cfg configuration 
                 file.
  chromo      => The name of the chromosome sizes text file, described 
                 above, as an alternative to providing the database name.
  bbapppath   => Provide the full path to Jim Kent's bedToBigBed  
                 utility. This parameter may instead be defined in the 
                 configuration file "tim_db_helper.cfg". 

Example

	my $wig_file = 'example_wig';
	my $bw_file = wig_to_bigwig_conversion( {
			'wig'   => $wig_file,
			'db'    => $database,
	} );
	if (-e $bw_file) {
		print " success! wrote bigwig file $bw_file\n";
		unlink $wig_file; # no longer necessary
	}
	else {
		print " failure! see STDERR for errors\n";
	};

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



