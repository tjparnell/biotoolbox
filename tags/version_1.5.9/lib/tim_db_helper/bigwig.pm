package tim_db_helper::bigwig;

# modules
require Exporter;
use strict;
use Carp;
use Bio::DB::BigWig;
use Bio::DB::BigFile;
use Bio::DB::BigWigSet;


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_bigwig_scores
	collect_bigwig_position_scores
	wig_to_bigwig_conversion
	open_bigwig_db
	open_bigwigset_db
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
	unless (scalar @_ >= 4) {
		confess " At least four arguments must be passed to collect BigWig data!\n";
	}
	my ($region, $region_strand, $stranded, @wig_features) = @_;
	
	# initialize 
	my %wig_data; # hash of position => score
	
	
	# Open a bigwig database object to collect the data
	# the region object passed to us might be from a BigWig database 
	# or it could be from a Bio::DB::SeqFeature::Store database
	# we can directly use the former, but if the latter then we need to open
	# a bigwig db connection
	# also, if there are more than one bigwig features to work with, then we 
	# will need to open them individually
	
	# region object is from BigWig database
	if (
		scalar @wig_features == 1 and
		ref $region eq 'Bio::DB::BigFile::Segment'
	) {
		# since bigwig files do not store strand information
		# we cannot check for strand information
		
		# collect the features from this region
		my @features = $region->features('region');
		
		# convert list of values to position => value
		foreach (@features) {
			my $pos = $_->start;
			$wig_data{$pos} = $_->score;
		}
	}
	
	# region is not from BigWig database or there are multiple bigwig features
	else {
		
		# there is usually only one bigwig feature sent, but there may be 
		# more to combine data
		foreach my $feature (@wig_features) {
		
			# Get the name of the bigwig file and check the strand
			my $wigfile;
			my $strand_check;
			
			if ($feature =~ /^file:(.+)$/) {
				# the passed feature appears to specify a file
				$wigfile = $1;
				
				# check the file
				unless (-e $wigfile) {
					confess " BigWig file '$wigfile' does not exist!\n";
					return;
				}
				
				# since it's a passed data file name, we can't check for strand
				# the bigwig file does not inherently have strand
				# a BigWigSet may support strand as feature attribute, but 
				# that's for another day
				# assume the strand is good, onus on the user
				$strand_check = 1;
			}
			elsif ($feature =~ /^http|ftp/i) {
				# a remote file
				# this should be supported by Bio::DB::BigWig
				$wigfile = $feature;
				
				# again, like the local file, can't perform strand check
				$strand_check = 1;
			}
			else {
				# passed feature could be a SeqFeature object
				
				# get wigfile name
				($wigfile) = $feature->get_tag_values('bigwigfile') or
					confess " passed feature '$feature' is not a SeqFeature object!\n";
				
				# check strand
					# database features will support the strand method
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
					$strand_check = 1;
				}
			}
			
			# confirm that we have acceptable data to collect
			next unless $strand_check == 1;
			confess " no wigfile passed!\n" unless $wigfile;
			
			# Open the BigWig file
			my $bw;
			if (exists $OPENED_BIGFILES{$wigfile} ) {
				# this file is already opened, use it
				$bw = $OPENED_BIGFILES{$wigfile};
			}
			else {
				# this file has not been opened yet, open it
				$bw = open_bigwig_db($wigfile) or 
					confess " unable to open data BigWig file '$wigfile'";
				
				# store the opened object for later use
				$OPENED_BIGFILES{$wigfile} = $bw;
			}
			
			# Collect from bigwig file
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
			
			# convert list of values to position => value
			foreach (@features) {
				my $pos = $_->start;
				$wig_data{$pos} = $_->score;
			}
		
		} # end of foreach loop
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
		
		# check that we a specified database
		my $db = $argument_ref->{'db'} || undef;
		unless ($db) {
			carp " database or chromosome file not specified! " . 
				"Unable to convert!\n";
			return;
		};
		
		# generate chromosome lengths file
		my @chromos = $db->seq_ids;
		unless (@chromos) {
			carp " no chromosome sequences identified in database!\n";
			return;
		}
		open CHR_FILE, ">tim_helper_chr_lengths.txt";
		foreach my $chr (@chromos) {
			my $segment = $db->segment($chr);
			print CHR_FILE "$chr\t", $segment->length, "\n";
		}
		close CHR_FILE;
		$chromo_file = "tim_helper_chr_lengths.txt";
	}
	
	# generate the bw file name
	# we can substitute one of three possible names for bw
	my $bw_file = $wigfile;
	$bw_file =~ s/\.(?:bed|bedgraph|wig)$/.bw/;
	
	# generate the bigwig file 
	if ($bw_app_path) {
		# we found Kent's utility
		# this is arguably the best method for converting
		# execute
		print " converting $wigfile to bigWig....\n";
		if ($bw_app_path =~ /wigToBigWig$/) {
			# include the -clip option in case there are any positions 
			# out of bounds of the chromosome
			# it will just warn instead of fail
			system $bw_app_path, '-clip', $wigfile, $chromo_file, $bw_file;
		}
		elsif ($bw_app_path =~ /bedGraphToBigWig$/) {
			# this doesn't have the -clip option, too bad
			system $bw_app_path, $wigfile, $chromo_file, $bw_file;
		}
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
	if (-e $bw_file and -s $bw_file) {
		# conversion successful
		if ($chromo_file eq 'tim_helper_chr_lengths.txt') {
			# we no longer need our temp chromosome file
			unlink $chromo_file;
		}
		return $bw_file;
	}
	else {
		warn " Conversion failed. You should try manually and watch for errors\n";
		if (-e $bw_file) {
			# 0-byte file was created
			unlink $bw_file;
		}
		if ($chromo_file eq 'tim_helper_chr_lengths.txt') {
			# leave the temp chromosome file as a courtesy
			warn " Leaving temporary chromosome file '$chromo_file'\n";
		}
		return;
	}
}



### Open a bigWig database connection
sub open_bigwig_db {
	
	my $path = shift;
	$path =~ s/^file://; # clean up file prefix if present
	
	# open the database connection 
	my $db;
	eval {
		$db = Bio::DB::BigWig->new( -bigwig => $path);
	};
	
	if ($db) {
		return $db;
	}
	else {
		carp " ERROR: can't open BigWig file '$path'!\n";
		return;
	}
}




### Open a bigWigSet database connection
sub open_bigwigset_db {
	
	my $directory = shift;
	
	# open the database connection 
	# we're using the region feature type because that's what the rest of 
	# tim_db_helper modules expect and work with
	my $bws = Bio::DB::BigWigSet->new(
				-dir            => $directory,
				-feature_type   => 'region',
	);
	
	# check that we haven't just opened a new empty bigwigset object
	my @paths = $bws->bigwigs;
	
	if (@paths) {
		# we have bigwig files, must be a valid bigwigset directory
		return $bws;
	}
	else {
		warn " ERROR: can't open BigWigSet directory '$directory'!\n";
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
bigwig file (.bw). The file may be identified in one of two ways. First,
it may be referenced in the database. Typically, a single 
feature representing the dataset is present across each chromosome. The 
feature should contain an attribute ('bigwigfile') that references the 
location of the binary file representing the dataset scores. Second, 
the local location of the file may be directly passed to the subroutine. 

In either case, the file is read using the Bio::DB::BigWig module, and 
the values extracted from the region of interest. 

Scores may be restricted to strand by specifying the desired strandedness. 
For example, to collect transcription data over a gene, pass the strandedness 
value 'sense'. If the strand of the region database object (representing the 
gene) matches the strand of the referencing wig file feature, then the data is 
collected. NOTE that this only works when the bigwig file is referenced through 
a database; directly referenced bigwig files do not have inherent strand 
information - all data will be collected!

For loading bigwig files into a Bio::DB database, see the biotoolbox perl 
script 'big_file2gff3.pl'.

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

This subroutine will collect only the score values from a binary BigWig file 
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
       to the BigWig file. They should contain the attribute 'bigwigfile'
       which has the path to the BigWig file. Alternatively, pass one 
       or more filenames of .bw files. Each filename should be 
       prefixed with 'file:' to indicate that it is a direct file 
       reference, and not a database object.

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_bigwig_position_scores

This subroutine will collect the score values from a binary BigWig file 
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
       to the BigWig file. They should contain the attribute 'bigwigfile'
       which has the path to the BigWig file. Alternatively, pass one 
       or more filenames of .bw files. Each filename should be 
       prefixed with 'file:' to indicate that it is a direct file 
       reference, and not a database object.

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
input wig file basename with the BigWig ".bw". Note that the it does 
not check for success of writing the bigwig file. Check STDERR for errors 
in bigwig file generation.

Pass the function an anonymous hash of arguments, including the following:

  Required:
  wig         => The name of the wig source file. 
  db          => Provide an opened database object from which to generate 
                 the chromosome sizes information.
  Optional: 
  chromo      => The name of the chromosome sizes text file, described 
                 above, as an alternative to providing the database name.
  bwapppath   => Provide the full path to Jim Kent's wigToBigWig 
                 utility. This parameter may instead be defined in the 
                 configuration file "biotoolbox.cfg". 

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

=item open_bigwig_db()

This subroutine will open a BigWig database connection. Pass either the 
local path to a bigWig file (.bw extension) or the URL of a remote bigWig 
file. It will return the opened database object.

=item open_bigwigset_db()

This subroutine will open a BigWigSet database connection using a directory 
of BigWig files and one metadata index file, as described in 
Bio::DB::BigWigSet. Essentially, this treats a directory of BigWig files as 
a single database with each BigWig file representing a different feature 
with unique attributes (type, source, strand, etc). 

Pass the subroutine a scalar value representing the local path to the 
directory. It presumes a feature_type of 'region', as expected by the other 
tim_db_helper subroutines and modules. It will return the opened database 
object.


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



