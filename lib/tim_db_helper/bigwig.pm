package tim_db_helper::bigwig;

# modules
require Exporter;
use strict;
use Carp;
use Bio::DB::BigWig;
use Bio::DB::BigFile;


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_bigwig_scores
	collect_bigwig_position_scores
	wig_to_bigwig_conversion
);


# Hashes of opened file objects
our %OPENED_BIGFILES; # opened bigwig file objects
	# in empirical testing, this doesn't really seem to speed things up
	# like I thought it would
	# oh well, keep it anyway????


# The true statement
1; 



### Modules ###



### Collect BigWig scores only
sub collect_bigwig_scores {
	
	# we will actually call collect_bigwig_position_scores()
	# but only return the values
	
	my %wig_data = collect_bigwig_position_scores(@_);
	
	# return the values
	return values %wig_data;
}




### Collect positioned BigWig scores
sub collect_bigwig_position_scores {
	
	# pass the required information
	my ($region, $region_strand, $stranded, @wig_features) = @_;
	
	# set up hash, position => score
	my %wig_data;
	
	# look at each wigfile
	# usually there is only one, but for stranded data there may be 
	# two wigfiles (+ and -), so we'll check each wig file for strand info
	foreach my $feature (@wig_features) {
	
		# Check which data to take based on strand
		if (
			$stranded eq 'all' # stranded data not requested
			or $feature->strand == 0 # unstranded data
			or ( 
				# sense data
				$region_strand == $feature->strand 
				and $stranded eq 'sense'
			) 
			or (
				# antisense data
				$region_strand != $feature->strand  
				and $stranded eq 'antisense'
			)
		) {
			# we have acceptable data to collect
			
			# collect from wigfile if present
			if ($feature->has_tag('bigwigfile') ) {
				
				# get wigfile name
				my @wigfiles = $feature->get_tag_values('bigwigfile');
				my $wigfile = shift @wigfiles;
				
				# check for opened wigfile
				my $bw;
				if (exists $OPENED_BIGFILES{$wigfile} ) {
					# this file is already opened, use it
					$bw = $OPENED_BIGFILES{$wigfile};
				}
				else {
					# this file has not been opened yet, open it
					unless (-e $wigfile) {
						croak " BigWig file '$wigfile' does not exist!\n";
						return;
					}
					$bw = Bio::DB::BigWig->new($wigfile);
					unless ($bw) {
						croak " unable to open data BigWig file '$wigfile'";
					}
					
					# store the opened object for later use
					$OPENED_BIGFILES{$wigfile} = $bw;
				}
				
				# We're not adjusting start and end points as with wig data
				# because the bigwig file is by default set up to cover the 
				# entire chromosome (chromosome information is required for 
				# generating bigwig files)
				
				# collect the features and values
				my @features = $bw->get_features_by_location(
					$region->seq_id, 
					$region->start, 
					$region->end
				);
				foreach (@features) {
					my $pos = $_->start;
					$wig_data{$pos} = $_->score;
				}
			}
		}
	}
	
	# return collected data
	return %wig_data;
}




### Wig to BigWig file conversion
sub wig_to_bigwig_conversion {
	
	# Collect passed arguments
	my $argument_ref = shift;
	unless ($argument_ref) {
		carp "no arguments passed!";
		return;
	}
	
	# wigfile
	my $wigfile = $argument_ref->{'wig'} || undef;
	unless ($wigfile) {
		carp "no wig file passed!";
		return;
	}
	
	# identify bigwig conversion utility
	my $bw_app_path = $argument_ref->{'bwapppath'} || undef;
	unless ($bw_app_path) {
		print " wigToBigWig utility not specified; using Bio::DB::BigFile\n";
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
	my $bw_file = $wigfile;
	$bw_file =~ s/\.wig$/.bw/;
	
	# generate the bigwig file 
	if ($bw_app_path) {
		# we found Kent's utility
		# this is arguably the best method for converting
		# execute
		print " converting $wigfile to bigWig....\n";
		system $bw_app_path, '-clip', $wigfile, $chromo_file, $bw_file;
	}
	else {
		# we are using the Bio::DB::BigFile module to generate the 
		# bigwig file
		# however, Lincoln notes that this method may be deprecated
		# in future versions
		# for the time being we will use this method as it avoids
		# having to hunt down Jim Kent's utility in the path
		
		# we'll use Lincoln's default values, which are slightly
		# different from Kent's default values in his utility
		# but I'm not sure the reasoning behind the differences
		Bio::DB::BigFile->createBigWig(
			$wigfile, 
			$chromo_file,
			$bw_file
		);
	}
	
	# check the result
	if (-e $bw_file) {
		# conversion successful
		if ($chromo_file eq 'tim_helper_chr_lengths.txt') {
			# we no longer need our temp chromosome file
			unlink $chromo_file;
		}
		return $bw_file;
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

tim_db_helper::bigwig

=head1 DESCRIPTION

This module supports the use of bigwig file in the biotoolbox scripts, both 
in the collection of data from a bigwig file, as well as the generation of 
bigwig files.

=head2 Data collection

This module is used to collect the dataset scores from a binary 
bigwig file (.bw) that is referenced in the database. Typically, a single 
feature representing the dataset is present across each chromosome. The 
feature should contain an attribute ('bigwigfile') that references the 
location of the binary file representing the dataset scores. The file is 
read using the Bio::DB::BigWig module, and the values extracted from the 
region of interest. 

Scores may be restricted to strand by specifying the desired strandedness. 
For example, to collect transcription data over a gene, pass the strandedness 
value 'sense'. If the strand of the region database object (representing the 
gene) matches the strand of the wig file data feature, then the data is 
collected.

For loading bigwig files into a Bio::DB database, see the biotoolbox perl 
script 'bw2gff3.pl'.

To speed up the program and avoid repetitive opening and 
closing of the files, the opened bigwig file object is stored in a global 
hash in case it is needed again.

=head2 File generation

This module also supports the generation of bigwig files. This is dependent 
on either Jim Kent's UCSC commandline utility, or Lincoln Stein's 
Bio::DB::BigFile support. It automates the collection of chromosome 
information in prepration of conversion, if necessary.

=head1 USAGE

The module requires Lincoln Stein's Bio::DB::BigWig to be installed. 

Load the module at the beginning of your program.

	use tim_db_helper::bigwig;

It will automatically export the name of the subroutines. 

=over

=item collect_bigwig_scores

This subroutine will collect only the score values from a binary bigwig file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed three or more arguments in the following order:
    
    1) The database object representing the genomic region of interest. 
       This should be a Bio::DB::SeqFeature object that supports the 
       start, end, and strand methods.
    2) The strand of the original feature (or region), -1, 0, or 1.
    3) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       "none" or "no". Only those scores which match the indicated 
       strandedness are collected.
    4) One or more database feature objects that contain the reference 
       to the wib file. They should contain the attribute 'wigfile'.

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_bigwig_position_scores

This subroutine will collect the score values from a binary wig file 
for the specified database region keyed by position. 

The subroutine is passed three or more arguments in the following order:
    
    1) The database object representing the genomic region of interest. 
       This should be a Bio::DB::SeqFeature object that supports the 
       start, end, and strand methods.
    2) The strand of the original feature (or region), -1, 0, or 1.
    3) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       "none" or "no". Only those scores which match the indicated 
       strandedness are collected.
    4) One or more database feature objects that contain the reference 
       to the .bw file. They should contain the attribute 'bigwigfile'.

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. Note that only one value is 
returned per position, regardless of the number of dataset features 
passed. Usually this isn't a problem as only one dataset is examined at a 
time.

=item wig_to_bigwig_conversion()

This subroutine will convert a wig file to a bigWig file. See the UCSC 
documentation regarding wig (http://genome.ucsc.edu/goldenPath/help/wiggle.html)
and bigWig (http://genome.ucsc.edu/goldenPath/help/bigWig.html) file formats. 
It preferentially uses Jim Kent's wigToBigWig utility to perform the 
conversion, although Lincoln Stein's Bio::DB::BigFile module may alternatively 
be used. One of these must be present on the system for the conversion to 
succeed. 

The conversion requires a list of chromosome name and sizes in a simple text 
file, where each line is comprised of two columns, "<chromosome name> 
<size in bases>". This file may be specified, or automatically generated if 
given a Bio::DB database name (preferred to ensure genome version 
compatibility).

The function returns the name of the bigWig file, which will be the 
input wig file basename with the extension ".bw". Note that the it does 
not check for success of writing the bigwig file. Check STDERR for errors 
in bigwig file generation.

Pass the function an anonymous hash of arguments, including the following:

  Required:
  wig         => The name of the wig source file. 
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
  bwapppath   => Provide the full path to Jim Kent's wigToBigWig 
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



