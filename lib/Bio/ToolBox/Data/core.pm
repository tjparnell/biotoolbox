package Bio::ToolBox::Data::core;

use warnings;
use strict;
use Carp qw(carp cluck croak confess);
use base 'Bio::ToolBox::Data::file';
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	use_bam_adapter
	use_big_adapter
);
use Module::Load;

our $VERSION = '1.70';

#### Initialization and verification ###############################################

sub new {
	my $class = shift;

	# in case someone calls this from an established object
	if ( ref $class ) {
		$class = ref $class;
	}

	# Initialize the hash structure
	my %data = (
		'program'           => undef,
		'feature'           => undef,
		'feature_type'      => undef,
		'db'                => undef,
		'format'            => q(),
		'gff'               => 0,
		'bed'               => 0,
		'ucsc'              => 0,
		'vcf'               => 0,
		'number_columns'    => 0,
		'last_row'          => 0,
		'headers'           => 1,
		'filename'          => q(),
		'basename'          => q(),
		'extension'         => q(),
		'path'              => q(),
		'comments'          => [],
		'data_table'        => [],
		'header_line_count' => 0,
		'zerostart'         => 0,
	);

	# Finished
	return bless \%data, $class;
}

sub verify {

	# this function does not rely on any self functions for two reasons
	# this is a low level integrity checker
	# this is very old code from before the days of an OO API of Bio::ToolBox
	# although a lot of things have been added and changed since then....
	my $self    = shift;
	my $silence = shift || 0;    # default is to yell

	# check for data table
	unless ( defined $self->{'data_table'}
		and ref( $self->{'data_table'} ) eq 'ARRAY' )
	{
		carp sprintf "\n DATA INTEGRITY ERROR: No data table in %s object!", ref $self
			unless $silence;
		return;
	}

	# check for last row index
	if ( defined $self->{'last_row'} ) {
		my $number = scalar( @{ $self->{'data_table'} } ) - 1;
		if ( $self->{'last_row'} != $number ) {
			$self->{'last_row'} = $number;
		}
	}
	else {
		# define it for them
		$self->{'last_row'} = scalar( @{ $self->{'data_table'} } ) - 1;
	}

	# check for consistent number of columns
	if ( defined $self->{'number_columns'} ) {
		my $number = $self->{'number_columns'} + 1;   # account for base 1 column indexing
		my @problems;
		my $too_low  = 0;
		my $too_high = 0;
		for my $row ( 0 .. $self->{'last_row'} ) {
			my $count = scalar @{ $self->{'data_table'}->[$row] };
			if ( $count != $number ) {
				push @problems, $row;
				$too_low++  if $count < $number;
				$too_high++ if $count > $number;

				# add null values
				if ( $count < $number ) {
					for ( $count .. $number ) {
						push @{ $self->{'data_table'}->[$row] }, '.';
					}
				}
			}
		}

		# we found errors
		if (@problems) {

			unless ($silence) {

				# collapse problem list into compact string
				# from http://www.perlmonks.org/?node_id=87538
				my $problem = join( ',', @problems );
				$problem =~ s/(?<!\d) (\d+) (?: ,((??{$++1})) (?!\d) )+/$1-$+/xg;

				if ($too_low) {
					print "
 COLUMN INCONSISTENCY ERRORS: $too_low rows had fewer than expected columns!
 padded rows $problem with null values\n";
				}
				if ($too_high) {
					print "
 COLUMN INCONSISTENCY ERRORS: $too_high rows had more columns than 
 expected!\n  Problem rows: $problem\n";
				}
			}
		}
	}
	else {
		# this wasn't set???? then set it
		$self->{'number_columns'} =
			scalar @{ $self->{'data_table'}->[0] } - 1;
	}

	# check metadata and table names
	my $mdcheck = 0;
	for my $i ( 1 .. $self->{'number_columns'} ) {
		if ( not exists $self->{$i} ) {
			unless ($silence) {
				print " METADATA ERROR: No metadata hash for column $i!";
			}
			$mdcheck++;
			next;    # the other tests will fail anyway
		}
		if ( $self->{$i}{'name'} ne $self->{'data_table'}->[0][$i] ) {
			unless ($silence) {
				printf "
TABLE/METADATA MISMATCH ERROR: Column header names don't match metadata 
name values for index $i! compare '%s' with '%s'\n",
					$self->{'data_table'}->[0][$i], $self->{$i}{'name'};
			}
			$mdcheck++;
		}
		if ( $self->{$i}{'index'} != $i ) {
			unless ($silence) {
				printf "
 METADATA INDEX ERROR: index $i metadata doesn't match its index value %s\n",
					$self->{$i}{'index'};
				$mdcheck++;
			}
		}
	}
	return if $mdcheck;

	### Defined file format structure integrity
	my $error;

	# check for proper gff structure
	if ( $self->{'gff'} ) {

		# if any of these checks fail, we will reset the gff version to
		# the default of 0, or no gff
		my $gff_check = 1;    # start with assumption it is true

		# check number of columns
		if ( $self->{'number_columns'} != 9 ) {
			$gff_check = 0;
			$error .= ' Number of columns not 9.';
		}

		# check column indices
		if (
			# column 1 should look like chromosome
			exists $self->{1}
			and $self->{1}{'name'} !~
			m/^\#? (?: chr | chromo | chromosome | seq | seq.id )/xi
			)
		{
			$gff_check = 0;
			$error .= ' Column 1 name not chromosome-like.';
		}
		if (
			# column 4 should look like start
			exists $self->{4} and $self->{4}{'name'} !~ m/start | pos | position/xi
			)
		{
			$gff_check = 0;
			$error .= ' Column 4 name not start-like.';
		}
		if (
			# column 5 should look like end
			exists $self->{5}
			and $self->{5}{'name'} !~ m/stop | end | pos | position/xi
			)
		{
			$gff_check = 0;
			$error .= ' Column 5 name not stop-like.';
		}
		if (
			# column 7 should look like strand
			exists $self->{7} and $self->{7}{'name'} !~ m/strand/i
			)
		{
			$gff_check = 0;
			$error .= ' Column 7 name not strand-like.';
		}

		# check column data
		unless ( $self->_column_is_integers( 4, 5 ) ) {
			$gff_check = 0;
			$error .= ' Columns 4 and 5 not integers.';
		}
		unless ( $self->_column_is_numeric(6) ) {
			$gff_check = 0;
			$error .= ' Column 6 not numeric.';
		}
		unless ( $self->_column_is_stranded(7) ) {
			$gff_check = 0;
			$error .= ' Column 7 not strand values.';
		}

		# update metadata as necessary
		if ( $gff_check == 1 ) {
			$self->{'bed'}     = 0;
			$self->{'ucsc'}    = 0;
			$self->{'vcf'}     = 0;
			$self->{'headers'} = 0;
			if ( $self->{'format'} and $self->{'format'} !~ m/g[tvf]f/i ) {
				$self->{'format'} = 'gff';    # generic
			}
			if ( $self->{'extension'} and $self->{'extension'} !~ m/\.g[tvf]f/i ) {
				$self->{'extension'} = '.gff';    # generic????
			}
		}
		else {
			$self->{'gff'}       = 0;
			$self->{'headers'}   = 1;
			$self->{'extension'} = '.txt';
			$self->{'format'}    = q();

			# remove the AUTO key from the metadata
			for my $i ( 1 .. $self->{'number_columns'} ) {
				if ( exists $self->{$i}{'AUTO'} ) {
					delete $self->{$i}{'AUTO'};
				}
			}
			print "\n GFF FILE FORMAT ERROR: $error\n" unless $silence;
		}
	}

	# check for proper BED structure
	if ( $self->{'bed'} ) {

		# if any of these checks fail, we will reset the bed flag to 0
		# to make it not a bed file format
		my $bed_check = 1;    # start with assumption it is correct

		# check number of columns
		if ( $self->{'number_columns'} < 3 ) {
			$bed_check = 0;
			$error .= ' Number of columns not at least 3.';
		}

		# check column index names
		if ( exists $self->{1}
			and $self->{1}{'name'} !~
			m/^\#? (?: chr | chromo | chromosome | seq | seq.id )/xi )
		{
			$bed_check = 0;
			$error .= ' Column 1 name not chromosome-like.';
		}
		if ( exists $self->{2}
			and $self->{2}{'name'} !~ m/start | pos | position/xi )
		{
			$bed_check = 0;
			$error .= ' Column 2 name not start-like.';
		}
		if ( exists $self->{3}
			and $self->{3}{'name'} !~ m/stop | end | pos | position/xi )
		{
			$bed_check = 0;
			$error .= ' Column 3 name not stop-like.';
		}
		if ( exists $self->{6}
			and $self->{6}{'name'} !~ m/strand/i )
		{
			$bed_check = 0;
			$error .= ' Column 6 name not strand-like.';
		}
		if (    exists $self->{7}
			and $self->{'format'}  !~ /narrow | broad/xi
			and $self->{7}{'name'} !~ m/start | thick | cds/xi )
		{
			$bed_check = 0;
			$error .= ' Column 7 name not thickStart-like.';
		}
		if (    exists $self->{8}
			and $self->{'format'}  !~ /narrow | broad/xi
			and $self->{8}{'name'} !~ m/end | stop | thick | cds/xi )
		{
			$bed_check = 0;
			$error .= ' Column 8 name not thickEnd-like.';
		}
		if (    exists $self->{9}
			and $self->{'format'}  !~ /narrow | broad/xi
			and $self->{9}{'name'} !~ m/item | rgb | color/xi )
		{
			$bed_check = 0;
			$error .= ' Column 9 name not itemRGB-like.';
		}
		if (    exists $self->{10}
			and $self->{'format'}   !~ /narrow | broad/xi
			and $self->{10}{'name'} !~ m/count | number | block | exon/xi )
		{
			$bed_check = 0;
			$error .= ' Column 10 name not blockCount-like.';
		}
		if ( exists $self->{11}
			and $self->{11}{'name'} !~ m/size | length | block | exon/xi )
		{
			$bed_check = 0;
			$error .= ' Column 11 name not blockSizes-like.';
		}
		if ( exists $self->{12}
			and $self->{12}{'name'} !~ m/start | block | exon/xi )
		{
			$bed_check = 0;
			$error .= ' Column 12 name not blockStarts-like.';
		}

		# check column data
		unless ( $self->_column_is_integers( 2, 3 ) ) {
			$bed_check = 0;
			$error .= ' Columns 2 and 3 not integers.';
		}
		if ( $self->{'number_columns'} >= 5 ) {

			# only check if it is actually present, since could be optional
			unless ( $self->_column_is_numeric(5) ) {
				$bed_check = 0;
				$error .= ' Column 5 not numeric.';
			}
		}
		if ( $self->{'number_columns'} >= 6 ) {

			# only check if it is actually present, since could be optional
			unless ( $self->_column_is_stranded(6) ) {
				$bed_check = 0;
				$error .= ' Column 6 not strand values.';
			}
		}
		if (    $self->{'format'}
			and $self->{'format'} =~ /narrow | broad/xi )
		{
			unless ( $self->_column_is_numeric( 7, 8, 9 ) ) {
				$bed_check = 0;
				$error .= ' Columns 7,8,9 not numeric.';
			}
		}
		if (    $self->{'format'}
			and $self->{'format'} =~ /narrow/i )
		{
			unless ( $self->_column_is_integers(10) ) {
				$bed_check = 0;
				$error .= ' Column 10 is not integers.';
			}
		}
		if ( $self->{'number_columns'} == 12 ) {

			# bed12 has extra special limitations
			unless ( $self->_column_is_integers( 7, 8, 10 ) ) {
				$bed_check = 0;
				$error .= ' Column 7,8,10 not integers.';
			}
			unless ( $self->_column_is_comma_integers( 11, 12 ) ) {
				$bed_check = 0;
				$error .= ' Column 11,12 not comma-delimited integers.';
			}
		}
		if (    $self->{'number_columns'} == 15
			and $self->{'format'} =~ /gapped/i )
		{
			# gappedPeak has extra special limitations
			unless ( $self->_column_is_integers( 7, 8, 10 ) ) {
				$bed_check = 0;
				$error .= ' Column 7,8,10 not integers.';
			}
			unless ( $self->_column_is_comma_integers( 11, 12 ) ) {
				$bed_check = 0;
				$error .= ' Column 11,12 not comma-delimited integers.';
			}
			unless ( $self->_column_is_numeric( 13, 14, 15 ) ) {
				$bed_check = 0;
				$error .= ' Columns 13,14,15 not numeric.';
			}
		}

		# peak file format
		if (    $self->{'format'}
			and $self->{'format'} =~ /narrowpeak/i
			and $self->{'number_columns'} != 10 )
		{
			$bed_check = 0;
			$error .= ' NarrowPeak has 10 columns.';
		}
		if (    $self->{'format'}
			and $self->{'format'} =~ /broadpeak/i
			and $self->{'number_columns'} != 9 )
		{
			$bed_check = 0;
			$error .= ' BroadPeak has 9 columns.';
		}
		if (    $self->{'format'}
			and $self->{'format'} =~ /gappedpeak/i
			and $self->{'number_columns'} != 15 )
		{
			$bed_check = 0;
			$error .= ' GappeddPeak has 15 columns.';
		}

		# update metadata as necessary
		if ( $bed_check == 1 ) {
			$self->{'bed'}     = $self->{'number_columns'};
			$self->{'gff'}     = 0;
			$self->{'ucsc'}    = 0;
			$self->{'vcf'}     = 0;
			$self->{'headers'} = 0;
			if (    $self->{'format'}
				and $self->{'format'} !~ m/ (?: bed | peak ) /xi )
			{
				$self->{'format'} = sprintf "bed%d", $self->{'number_columns'};
			}
			if (    $self->{'extension'}
				and $self->{'extension'} !~ m/ (?: bed | peak ) /xi )
			{
				$self->{'extension'} = sprintf ".%s",
					$self->{'format'} =~ /peak/ ? lc $self->{'format'} : 'bed';
			}
		}
		else {
			$self->{'bed'}       = 0;
			$self->{'headers'}   = 1;
			$self->{'extension'} = '.txt';
			$self->{'format'}    = q();

			# remove the AUTO key from the metadata
			for my $i ( 1 .. $self->{'number_columns'} ) {
				if ( exists $self->{$i}{'AUTO'} ) {
					delete $self->{$i}{'AUTO'};
				}
			}
			print "\n BED FILE FORMAT ERROR: $error\n" unless $silence;
		}
	}

	# check refFlat or genePred gene structure
	if ( $self->{'ucsc'} ) {

		# if any of these checks fail, we will reset the extension
		my $ucsc_check = 1;    # start with assumption it is correct
		my $ucsc_type;

		# check number of columns
		my $colnumber = $self->{number_columns};
		if ( $colnumber == 16 ) {
			$ucsc_type = 'genePredExtBin';

			# genePredExt with bin
			# bin name chrom strand txStart txEnd cdsStart cdsEnd
			# exonCount exonStarts exonEnds score name2 cdsStartSt
			# cdsEndStat exonFrames
			unless ( $self->{3}{name} =~
				m/^\#? (?: chr | chromo | chromosome | seq | seq.id )/xi )
			{
				$ucsc_check = 0;
				$error .= ' Column 3 name not chromosome-like.';
			}
			unless ( $self->{4}{name} =~ m/strand/i ) {
				$ucsc_check = 0;
				$error .= ' Column 4 name not strand-like.';
			}
			unless ( $self->{5}{name} =~ m/start | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 5 name not start-like.';
			}
			unless ( $self->{6}{name} =~ m/stop | end | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 6 name not stop-like.';
			}
			unless ( $self->{7}{name} =~ m/start | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 7 name not start-like.';
			}
			unless ( $self->{8}{name} =~ m/stop | end | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 8 name not stop-like.';
			}
			unless ( $self->_column_is_integers( 5, 6, 7, 8, 9 ) ) {
				$ucsc_check = 0;
				$error .= ' Columns 5,6,7,8,9 not integers.';
			}
			unless ( $self->_column_is_comma_integers( 10, 11 ) ) {
				$ucsc_check = 0;
				$error .= ' Columns 10,11 not comma-delimited integers.';
			}
			unless ( $self->_column_is_stranded(4) ) {
				$ucsc_check = 0;
				$error .= ' Column 4 not strand values.';
			}
		}
		elsif ( $colnumber == 15 or $colnumber == 12 ) {
			$ucsc_type = $colnumber == 15 ? 'genePredExt' : 'knownGene';

			# GenePredExt
			# name chrom strand txStart txEnd cdsStart cdsEnd
			# exonCount exonStarts exonEnds score name2 cdsStartSt
			# cdsEndStat exonFrames
			# or knownGene
			# name chrom strand txStart txEnd cdsStart cdsEnd
			# exonCount exonStarts exonEnds proteinID alignID
			unless ( $self->{2}{name} =~
				m/^\#? (?: chr | chromo | chromosome | seq | seq.id )/xi )
			{
				$ucsc_check = 0;
				$error .= ' Column 2 name not chromosome-like.';
			}
			unless ( $self->{3}{name} =~ m/strand/i ) {
				$ucsc_check = 0;
				$error .= ' Column 3 name not strand-like.';
			}
			unless ( $self->{4}{name} =~ m/start | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 4 name not start-like.';
			}
			unless ( $self->{5}{name} =~ m/stop | end | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 5 name not stop-like.';
			}
			unless ( $self->{6}{name} =~ m/start | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 6 name not start-like.';
			}
			unless ( $self->{7}{name} =~ m/stop | end | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 7 name not stop-like.';
			}
			unless ( $self->_column_is_integers( 4, 5, 6, 7, 8 ) ) {
				$ucsc_check = 0;
				$error .= ' Columns 4,5,6,7,8 not integers.';
			}
			unless ( $self->_column_is_comma_integers( 9, 10 ) ) {
				$ucsc_check = 0;
				$error .= ' Columns 9,10 not comma-delimited integers.';
			}
			unless ( $self->_column_is_stranded(3) ) {
				$ucsc_check = 0;
				$error .= ' Column 3 not strand values.';
			}
		}
		elsif ( $colnumber == 11 ) {
			$ucsc_type = 'refFlat';

			# geneName transcriptName chrom strand txStart txEnd
			# cdsStart cdsEnd exonCount exonStarts exonEnds
			unless ( $self->{3}{name} =~
				m/^\#? (?: chr | chromo | chromosome | seq | seq.id )/xi )
			{
				$ucsc_check = 0;
				$error .= ' Column 3 name not chromosome-like.';
			}
			unless ( $self->{4}{name} =~ m/strand/i ) {
				$ucsc_check = 0;
				$error .= ' Column 4 name not strand-like.';
			}
			unless ( $self->{5}{name} =~ m/start | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 5 name not start-like.';
			}
			unless ( $self->{6}{name} =~ m/stop | end | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 6 name not stop-like.';
			}
			unless ( $self->{7}{name} =~ m/start | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 7 name not start-like.';
			}
			unless ( $self->{8}{name} =~ m/stop | end | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 8 name not stop-like.';
			}
			unless ( $self->_column_is_integers( 5, 6, 7, 8, 9 ) ) {
				$ucsc_check = 0;
				$error .= ' Columns 5,6,7,8,9 not integers.';
			}
			unless ( $self->_column_is_comma_integers( 10, 11 ) ) {
				$ucsc_check = 0;
				$error .= ' Columns 10,11 not comma-delimited integers.';
			}
			unless ( $self->_column_is_stranded(4) ) {
				$ucsc_check = 0;
				$error .= ' Column 4 not strand values.';
			}
		}
		elsif ( $colnumber == 10 ) {
			$ucsc_type = 'genePred';

			# name chrom strand txStart txEnd cdsStart cdsEnd
			# exonCount exonStarts exonEnds
			unless ( $self->{2}{name} =~
				m/^\#? (?: chr | chromo | chromosome | seq | seq.id )/xi )
			{
				$ucsc_check = 0;
				$error .= ' Column 2 name not chromosome-like.';
			}
			unless ( $self->{3}{name} =~ m/strand/i ) {
				$ucsc_check = 0;
				$error .= ' Column 3 name not strand-like.';
			}
			unless ( $self->{4}{name} =~ m/start | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 4 name not start-like.';
			}
			unless ( $self->{5}{name} =~ m/stop | end | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 5 name not stop-like.';
			}
			unless ( $self->{6}{name} =~ m/start | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 6 name not start-like.';
			}
			unless ( $self->{7}{name} =~ m/stop | end | pos | position/xi ) {
				$ucsc_check = 0;
				$error .= ' Column 7 name not stop-like.';
			}
			unless ( $self->_column_is_integers( 4, 5, 6, 7, 8 ) ) {
				$ucsc_check = 0;
				$error .= ' Columns 4,5,6,7,8 not integers.';
			}
			unless ( $self->_column_is_comma_integers( 9, 10 ) ) {
				$ucsc_check = 0;
				$error .= ' Columns 9,10 not comma-delimited integers.';
			}
			unless ( $self->_column_is_stranded(3) ) {
				$ucsc_check = 0;
				$error .= ' Column 3 not strand values.';
			}
		}
		else {
			$ucsc_type  = 'UCSC';
			$ucsc_check = 0;
			$error .= ' Wrong # of columns.';
		}

		# update metadata as necessary
		if ( $ucsc_check == 1 ) {
			$self->{'bed'}    = 0;
			$self->{'gff'}    = 0;
			$self->{'vcf'}    = 0;
			$self->{'format'} = $ucsc_type;
			if (    $self->{'extension'}
				and $self->{'extension'} !~ m/(?: ucsc | txt | $ucsc_type )/xi )
			{
				$self->{'extension'} = ".ucsc";    # generic extension
			}
		}
		else {
			$self->{'ucsc'}      = 0;
			$self->{'extension'} = '.txt';
			$self->{'format'}    = q();

			# remove the AUTO key
			for my $i ( 1 .. $self->{'number_columns'} ) {
				if ( exists $self->{$i}{'AUTO'} ) {
					delete $self->{$i}{'AUTO'};
				}
			}
			print "\n $ucsc_type FILE FORMAT ERROR: $error\n" unless $silence;
		}
	}

	# check VCF format
	if ( $self->{vcf} ) {

		# if any of these checks fail, we will reset the vcf flag to 0
		# to make it not a vcf file format
		my $vcf_check = 1;    # start with assumption it is correct

		# check number of columns
		if ( $self->{'number_columns'} < 8 ) {
			$vcf_check = 0;
			$error .= ' Number of Columns is too few.';
		}

		# check column index names
		if ( $self->{1}{'name'} !~ m/chrom/i ) {
			$vcf_check = 0;
			$error .= ' Column 1 name not chromosome.';
		}
		if ( $self->{2}{'name'} !~ m/^pos|start/i ) {
			$vcf_check = 0;
			$error .= ' Column 2 name not position.';
		}

		# check column data
		unless ( $self->_column_is_integers(2) ) {
			$vcf_check = 0;
			$error .= ' Columns 2 not integers.';
		}

		# reset the vcf tag value as appropriate
		# update metadata as necessary
		if ( $vcf_check == 1 ) {
			$self->{'bed'}  = 0;
			$self->{'gff'}  = 0;
			$self->{'ucsc'} = 0;
			if ( $self->{'format'} and $self->{'format'} !~ m/vcf/i ) {
				$self->{'format'} = 'vcf';    # generic
			}
			if ( $self->{'extension'} and $self->{ 'extension' !~ /\.vcf/ } ) {
				$self->{'extension'} = '.vcf';
			}
		}
		else {
			# reset metadata
			$self->{'vcf'}       = 0;
			$self->{'headers'}   = 1;
			$self->{'format'}    = q();
			$self->{'extension'} = '.txt';

			# remove the AUTO key from the metadata
			for my $i ( 1 .. $self->{'number_columns'} ) {
				if ( exists $self->{$i}{'AUTO'} ) {
					delete $self->{$i}{'AUTO'};
				}
			}
			print "\n VCF FILE FORMAT ERROR: $error\n" unless $silence;
		}
	}

	# check proper SGR file structure
	if (    exists $self->{'extension'}
		and defined $self->{'extension'}
		and $self->{'extension'} =~ /sgr/i )
	{
		# there is no sgr field in the data structure
		# so we're just checking for the extension
		# we will change the extension as necessary if it doesn't conform
		my $sgr_check = 1;
		if ( $self->{'number_columns'} != 3 ) {
			$sgr_check = 0;
			$error .= ' Column number is not 3.';
		}
		if ( $self->{1}{'name'} !~
			m/^\#? (?: chr | chromo | chromosome | seq | seq.id )/xi )
		{
			$sgr_check = 0;
			$error .= ' Column 1 name not chromosome-like.';
		}
		if ( $self->{2}{'name'} !~ m/start | pos | position/xi ) {
			$sgr_check = 0;
			$error .= ' Column 2 name not start-like.';
		}
		unless ( $self->_column_is_integers(2) ) {
			$sgr_check = 0;
			$error .= ' Columns 2 not integers.';
		}

		# update metadata as necessary
		if ( $sgr_check == 1 ) {
			$self->{'bed'}     = 0;
			$self->{'gff'}     = 0;
			$self->{'ucsc'}    = 0;
			$self->{'vcf'}     = 0;
			$self->{'headers'} = 0;
			if ( $self->{'format'} and $self->{'format'} !~ m/sgr/i ) {
				$self->{'format'} = 'sgr';
			}
			if ( $self->{'extension'} and $self->{ 'extension' !~ /\.sgr/ } ) {
				$self->{'extension'} = '.sgr';
			}
		}
		else {

			# doesn't smell like a SGR file
			# change the extension so the write subroutine won't think it is
			# make it a text file
			$self->{'extension'} =~ s/sgr/txt/i;
			$self->{'headers'} = 1;
			$self->{'format'}  = q();

			# remove the AUTO key from the metadata
			for my $i ( 1 .. $self->{'number_columns'} ) {
				if ( exists $self->{$i}{'AUTO'} ) {
					delete $self->{$i}{'AUTO'};
				}
			}
			print "\n SGR FILE FORMAT ERROR: $error\n" unless $silence;
		}
	}

	# if we have made it here, then there were no major structural problems
	# file is verified, any minor issues should have been fixed
	return 1;
}

# internal method to check if a column is nothing but integers, i.e. start, stop
sub _column_is_integers {
	my $self  = shift;
	my @index = @_;
	return 1 if ( $self->{last_row} == 0 );    # can't check if table is empty
	foreach (@index) {
		return 0 unless exists $self->{$_};
	}
	for my $row ( 1 .. $self->{last_row} ) {
		for my $i (@index) {
			return 0 unless ( $self->{data_table}->[$row][$i] =~ /^\d+$/ );
		}
	}
	return 1;
}

# internal method to check if a column appears numeric, i.e. scores
sub _column_is_numeric {
	my $self  = shift;
	my @index = @_;
	return 1 if ( $self->{last_row} == 0 );    # can't check if table is empty
	foreach (@index) {
		return 0 unless exists $self->{$_};
	}
	for my $row ( 1 .. $self->{last_row} ) {
		for my $i (@index) {

			# we have a very loose definition of numeric: exponents, signs, commas
			return 0 unless ( $self->{data_table}->[$row][$i] =~ m/^ [\d\-\+\.,eE]+ $/x );
		}
	}
	return 1;
}

# internal method to check if a column is nothing but comma delimited integers
sub _column_is_comma_integers {
	my $self  = shift;
	my @index = @_;
	return 1 if ( $self->{last_row} == 0 );    # can't check if table is empty
	foreach (@index) {
		return 0 unless exists $self->{$_};
	}
	for my $row ( 1 .. $self->{last_row} ) {
		for my $i (@index) {
			return 0 unless ( $self->{data_table}->[$row][$i] =~ /^[\d,]+$/ );
		}
	}
	return 1;
}

# internal method to check if a column looks like a strand column
sub _column_is_stranded {
	my ( $self, $index ) = @_;
	return unless exists $self->{$index};
	for my $row ( 1 .. $self->{last_row} ) {
		return 0
			if ( $self->{data_table}->[$row][$index] !~
				m/^(?: \-1 | 0 | 1 | \+ | \- | \. )$/x );
	}
	return 1;
}

#### Database methods ##############################################################

sub open_database {
	my $self = shift;
	if ( not defined $_[0] ) {
		return $self->open_meta_database;
	}
	elsif ( $_[0] eq '0' or $_[0] eq '1' ) {

		# likely a boolean value to indicate force
		return $self->open_meta_database( $_[0] );
	}
	elsif ( $_[0] =~ /[a-zA-Z]+/ ) {

		# likely the name of a database
		return $self->open_new_database(@_);
	}
	else {
		# original default
		return $self->open_meta_database(@_);
	}
}

sub open_meta_database {
	my $self  = shift;
	my $force = shift || 0;
	return unless $self->{db};
	return if $self->{db} =~ /^Parsed:/;    # we don't open parsed annotation files
	if ( exists $self->{db_connection} ) {
		return $self->{db_connection} unless $force;
	}
	my $db = open_db_connection( $self->{db}, $force );
	return unless $db;
	$self->{db_connection} = $db;
	return $db;
}

sub open_new_database {
	my $self     = shift;
	my $database = shift;
	my $force    = shift || 0;
	return open_db_connection( $database, $force );
}

sub verify_dataset {
	my ( $self, $dataset, $database ) = @_;
	return unless $dataset;
	if ( exists $self->{verfied_dataset}{$dataset} ) {
		return $self->{verfied_dataset}{$dataset};
	}
	else {
		if ( $dataset =~ /^(?: file | http | ftp )/x ) {

			# local or remote file already verified?
			$self->{verfied_dataset}{$dataset} = $dataset;
			return $dataset;
		}
		$database ||= $self->open_meta_database;
		my ($verified) = verify_or_request_feature_types(

			# normally returns an array of verified features, we're only checking one
			db      => $database,
			feature => $dataset,
		);
		if ($verified) {
			$self->{verfied_dataset}{$dataset} = $verified;
			return $verified;
		}
	}
	return;
}

#### Column Manipulation ####

sub delete_column {
	my $self = shift;

	# check for Stream
	if ( ref($self) eq 'Bio::ToolBox::Data::Stream' ) {
		unless ( $self->mode ) {
			cluck 'ERROR: We have a read-only Stream object, cannot add columns';
			return;
		}
		if ( defined $self->{fh} ) {

			# Stream file handle is opened
			cluck 'ERROR: Cannot modify columns when a Stream file handle is opened!';
			return;
		}
	}
	unless (@_) {
		carp 'ERROR: must provide a list';
		return;
	}

	# identify columns to retain
	my @deletion_list = @_;
	if ( scalar @_ > 1 ) {
		@deletion_list = sort { $a <=> $b } @_;
	}
	my @retain_list;
	for my $i ( 1 .. $self->{'number_columns'} ) {

		# compare each current index with the first one in the list of
		# deleted indices. if it matches, delete. if not, keep
		if ( defined $deletion_list[0] and $i == $deletion_list[0] ) {
			shift @deletion_list;    # this particular index should be deleted
		}
		else {
			push @retain_list, $i;    # this particular index should be kept
		}
	}
	return $self->reorder_column(@retain_list);
}

sub reorder_column {
	my $self = shift;

	# check for Stream
	if ( ref($self) eq 'Bio::ToolBox::Data::Stream' ) {
		unless ( $self->mode ) {
			cluck 'ERROR: We have a read-only Stream object, cannot add columns';
			return;
		}
		if ( defined $self->{fh} ) {

			# Stream file handle is opened
			cluck 'ERROR: Cannot modify columns when a Stream file handle is opened!';
			return;
		}
	}

	# reorder data table
	unless (@_) {
		carp 'ERROR: must provide a list';
		return;
	}
	my @order = @_;
	unshift @order, 0;    # to maintain 1-based indexing
	for my $r ( 0 .. $self->last_row ) {
		my @old = $self->row_values($r);    # this is an array copy
		unshift @old, q();                  # to maintain 1-based indexing
		my @new = map { $old[$_] } @order;
		splice( @{ $self->{data_table} }, $r, 1, \@new );
	}
	$self->{data_table}->[0][0] = 'BLANK';    # not necessary but consistent

	# reorder metadata
	# copy the current metadata info hash into a temporary hash
	# then copy back from the old_metadata into the main data hash
	# using the new index number in the @order array
	# must regenerate the hash, not just link to the old anonymous hash, in
	# case we're duplicating columns
	my %old_metadata;
	for my $i ( 1 .. $self->{'number_columns'} ) {
		$old_metadata{$i} = $self->{$i};
		delete $self->{$i};
	}
	for my $i ( 1 .. ( scalar(@order) - 1 ) ) {

		#specifically skipping 0 because it doesn't have metadata
		$self->{$i} = {};
		foreach my $k ( keys %{ $old_metadata{ $order[$i] } } ) {
			$self->{$i}{$k} = $old_metadata{ $order[$i] }{$k};
		}
		$self->{$i}{'index'} = $i;    # reset index
	}

	# additional metadata
	$self->{'number_columns'} = scalar(@order) - 1;    # ignore first index
	delete $self->{column_indices} if exists $self->{column_indices};

	return 1;
}

#### General Metadata ####

sub feature {
	my $self = shift;
	if (@_) {
		$self->{feature} = shift;
	}
	return $self->{feature};
}

sub feature_type {
	my $self = shift;
	carp 'ERROR: feature_type is a read only method' if @_;
	if ( defined $self->{feature_type} ) {
		return $self->{feature_type};
	}
	my $feature_type;
	if ( $self->chromo_column and $self->start_column ) {
		$feature_type = 'coordinate';
	}
	elsif ($self->id_column
		or ( $self->type_column and $self->name_column )
		or ( $self->feature     and $self->name_column ) )
	{
		$feature_type = 'named';
	}
	else {
		$feature_type = 'unknown';
	}
	$self->{feature_type} = $feature_type;
	return $feature_type;
}

sub program {
	my $self = shift;
	if (@_) {
		$self->{program} = shift;
	}
	return $self->{program};
}

sub database {
	my $self = shift;
	if (@_) {
		$self->{db} = shift;
		if ( exists $self->{db_connection} and $self->{db_connection} ) {
			$self->open_meta_database(1);
		}
	}
	return $self->{db};
}

sub bam_adapter {
	my $self = shift;
	return use_bam_adapter(@_);
}

sub big_adapter {
	my $self = shift;
	return use_big_adapter(@_);
}

sub format {
	my $self = shift;
	if ( defined $_[0] ) {
		$self->{format} = $_[0];
	}
	return $self->{format};
}

sub gff {
	my $self = shift;
	if ( defined $_[0] and $_[0] =~ m/^(?: 0 | 1 | 2 | 2\.[2|5] | 3 )$/x ) {
		$self->{gff} = $_[0];
	}
	return $self->{gff};
}

sub bed {
	my $self = shift;
	if ( defined $_[0] and $_[0] =~ /^\d+$/ ) {
		$self->{bed} = $_[0];
	}
	return $self->{bed};
}

sub ucsc {
	my $self = shift;
	if ( defined $_[0] and $_[0] =~ /^\d+$/ ) {
		$self->{ucsc} = $_[0];
	}
	return $self->{ucsc};
}

sub vcf {
	my $self = shift;
	if ( defined $_[0] and $_[0] =~ /^[\d\.]+$/ ) {
		$self->{vcf} = $_[0];
	}
	return $self->{vcf};
}

sub number_columns {
	my $self = shift;
	carp 'ERROR: number_columns is a read only method' if @_;
	return $self->{'number_columns'};
}

sub number_rows {
	my $self = shift;
	carp 'ERROR: number_rows is a read only method' if @_;
	return $self->{last_row};
}

sub last_column {
	my $self = shift;
	carp 'ERROR: last_column is a read only method' if @_;
	return $self->{'number_columns'};
}

sub last_row {
	my $self = shift;
	carp 'ERROR: last_row is a read only method' if @_;
	return $self->{last_row};
}

sub filename {
	my $self = shift;
	carp 'ERROR: filename is a read only method. Use add_file_metadata().' if @_;
	return $self->{filename};
}

sub basename {
	my $self = shift;
	carp 'ERROR: basename is a read only method. Use add_file_metadata().' if @_;
	return $self->{basename};
}

sub path {
	my $self = shift;
	carp 'ERROR: path is a read only method. Use add_file_metadata().' if @_;
	return $self->{path};
}

sub extension {
	my $self = shift;
	carp 'ERROR: extension() is a read only method. Use add_file_metadata().' if @_;
	return $self->{extension};
}

#### General Comments ####

sub comments {
	my $self     = shift;
	my @comments = @{ $self->{comments} };
	foreach (@comments) {s/[\r\n]+//g}

	# comments are not chomped when loading
	# side effect of dealing with rare commented header lines with null values at end
	return @comments;
}

sub add_comment {
	my $self    = shift;
	my $comment = shift or return;

	# comment is not required to be prefixed with "# ", it will be added when saving
	push @{ $self->{comments} }, $comment;
	return 1;
}

sub delete_comment {
	my $self  = shift;
	my $index = shift;
	if ( defined $index ) {
		eval { splice @{ $self->{comments} }, $index, 1 };
	}
	else {
		$self->{comments} = [];
	}
}

sub vcf_headers {
	my $self = shift;
	return unless $self->vcf;
	return $self->{vcf_headers} if exists $self->{vcf_headers};
	my $headers = {};
	foreach my $comment ( $self->comments ) {
		my ( $key, $value );
		if ( $comment =~ m/^ \#\# ( [\w\-\.]+ ) = (.+) $/x ) {
			$key   = $1;
			$value = $2;
		}
		else {
			# invalid vcf header format!?
			next;
		}
		if ( $value !~ /^<.+>$/ ) {

			# simple value
			$headers->{$key} = $value;
		}
		else {
			# process complex values
			# extract ID with regex which should have
			my $id = ( $value =~ m/ID = ( [\w\-\.:]+ )/x )[0];
			$headers->{$key}{$id} = $value;
		}
	}

	# store and return
	$self->{vcf_headers} = $headers;
	return $headers;
}

sub rewrite_vcf_headers {
	my $self = shift;
	return unless $self->vcf;
	return unless exists $self->{vcf_headers};
	my @newComments;

	# file format always comes first
	push @newComments, sprintf( "##fileformat=%s\n", $self->{vcf_headers}{fileformat} );

	# common attributes
	foreach my $key ( sort { $a cmp $b } keys %{ $self->{vcf_headers} } ) {
		next if $key eq 'fileformat';
		if ( ref( $self->{vcf_headers}{$key} ) eq 'HASH' ) {

			# we have a complex VCF header with multiple keys
			# we will rewrite for each ID
			foreach my $id (
				sort { $a cmp $b }
				keys %{ $self->{vcf_headers}{$key} }
				)
			{
				# to avoid complexity of writing correct formatting
				push @newComments,
					sprintf( "##%s=%s\n", $key, $self->{vcf_headers}{$key}{$id} );
			}
		}
		else {
			# a simple value
			push @newComments, sprintf( "##%s=%s\n", $key, $self->{vcf_headers}{$key} );
		}
	}

	# replace the headers
	$self->{comments} = \@newComments;
}

#### Column Metadata ####

sub list_columns {
	my $self = shift;
	carp 'ERROR: list_columns is a read only method' if @_;
	my @list = map { $self->name($_) } ( 1 .. $self->{'number_columns'} );
	return wantarray ? @list : \@list;
}

sub name {
	my $self = shift;
	my ( $index, $new_name ) = @_;
	return unless defined $index;
	return unless exists $self->{$index};
	if ( defined $new_name ) {
		$self->{$index}{name} = $new_name;
		$self->{data_table}->[0][$index] = $new_name;
	}
	return $self->{$index}{'name'};
}

sub metadata {
	my $self = shift;
	my ( $index, $key, $value ) = @_;
	return unless defined $index;
	return unless exists $self->{$index};
	if ( $key and $key eq 'name' ) {
		return $self->name( $index, $value );
	}
	if ( $key and defined $value ) {

		# we are setting a new value
		$self->{$index}{$key} = $value;
		return $value;
	}
	elsif ( $key and not defined $value ) {
		if ( exists $self->{$index}{$key} ) {

			# retrieve a value
			return $self->{$index}{$key};
		}
		else {
			# key does not exist
			return;
		}
	}
	else {
		my %hash = %{ $self->{$index} };
		return wantarray ? %hash : \%hash;
	}
}

sub delete_metadata {
	my $self = shift;
	my ( $index, $key ) = @_;
	return unless defined $index;
	if ( defined $key ) {
		if ( exists $self->{$index}{$key} ) {
			return delete $self->{$index}{$key};
		}
	}
	else {
		# user wants to delete the metadata
		# but we need to keep the basics name and index
		foreach my $key ( keys %{ $self->{$index} } ) {
			next if $key eq 'name';
			next if $key eq 'index';
			delete $self->{$index}{$key};
		}
	}
}

sub copy_metadata {
	my ( $self, $source, $target ) = @_;
	return unless ( exists $self->{$source}{name} and exists $self->{$target}{name} );
	my $md = $self->metadata($source);
	delete $md->{name};
	delete $md->{'index'};
	delete $md->{'AUTO'} if exists $md->{'AUTO'};   # presume this is no longer auto index
	foreach ( keys %{$md} ) {
		$self->{$target}{$_} = $md->{$_};
	}
	return 1;
}

#### Column Indices ####

sub find_column {
	my ( $self, $name ) = @_;
	return unless $name;

	# the $name variable will be used as a regex in identifying the name
	# fix it so that it will possible accept a # character at the beginning
	# without a following space, in case the first column has a # prefix
	# also place the remainder of the text in a non-capturing parentheses for
	# grouping purposes while maintaining the anchors
	$name =~ s/ \A (\^?) (.+) (\$?)\Z /$1#?(?:$2)$3/x;

	# walk through each column index and check the name
	my $index;
	for my $i ( 1 .. $self->{'number_columns'} ) {
		if ( $self->{$i}{'name'} =~ /$name/i ) {
			$index = $i;
			last;
		}
	}
	return $index;
}

sub _find_column_indices {
	my $self = shift;

	# these are hard coded index name regex to accomodate different possibilities
	# these do not include parentheses for grouping
	# non-capturing parentheses will be added later in the sub for proper
	# anchoring and grouping - long story why, don't ask
	# some of these have progressively more relaxed criteria to avoid false positives
	my $name = $self->find_column('^name$');
	unless ( defined $name ) {

		# try more agressive search for names
		$name = $self->find_column('^gene.?name|transcript.?name$');
	}
	unless ( defined $name ) {

		# take any column with "name" in it
		$name = $self->find_column('name');
	}
	unless ( defined $name ) {

		# try even more agressive for annotation identifiers
		$name = $self->find_column('^gene|gene.?id|transcript.?id$');
	}
	unless ( defined $name ) {

		# use an id
		$name = $self->find_column('^id$');
	}
	my $type   = $self->find_column('^type|class|primary_tag|biotype');
	my $id     = $self->find_column('^primary_id$');    # database specific primary id
	my $chromo = $self->find_column('^chromosome|chromo|chrom|chr|seq.*id$');
	unless ( defined $chromo ) {
		$chromo = $self->find_column('chr');
	}
	my $start = $self->find_column('^start$');
	unless ( defined $start ) {
		$start = $self->find_column('start');
	}
	unless ( defined $start ) {
		$start = $self->find_column('^position|pos');
	}
	my $stop = $self->find_column('^stop|end$');
	unless ( defined $stop ) {
		$stop = $self->find_column('stop|end');
	}
	my $strand = $self->find_column('^strand');
	my $score  = $self->find_column('^score$');

	# check for coordinate string
	my $coord = $self->find_column('^coordinate$');
	if ( not $coord and not $chromo and not $start ) {

		# check if ID or name id looks like a coordinate string
		if ( $id and $self->{data_table}->[1] ) {
			if ( $self->{data_table}->[1][$id] =~
				m/^ [\w\-\.]+: \d+ (?: [\-\.]{1,2} \d+ )?$/x )
			{

				# first value looks like a coordinate string
				$coord = $id;
			}
		}
		if (    not $coord
			and $name
			and defined $self->{data_table}->[1] )
		{
			if ( $self->{data_table}->[1][$name] =~
				m/^ [\w\-\.]+: \d+ (?: [\-\.]{1,2} \d+ )?$/x )
			{
				# first value looks like a coordinate string
				$coord = $name;
			}
		}
	}

	# determine zero-based start coordinates
	if (    $self->{zerostart} == 0
		and $start
		and substr( $self->name($start), -1 ) eq '0' )
	{
		$self->{zerostart} = 1;
	}

	# cache the indexes
	$self->{column_indices} = {
		'name'   => $name,
		'type'   => $type,
		'id'     => $id,
		'seq_id' => $chromo,
		'chromo' => $chromo,
		'start'  => $start,
		'stop'   => $stop,
		'end'    => $stop,
		'strand' => $strand,
		'score'  => $score,
		'coord'  => $coord,
	};
	return 1;
}

sub chromo_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	if ( defined $_[0] ) {
		if ( $_[0] =~ /^\d+$/ and $_[0] <= $self->{'number_columns'} ) {
			$self->{column_indices}{chromo} = $_[0];
		}
		else {
			carp sprintf( "chromo_column %s is not a valid index!", $_[0] );
		}
	}
	return $self->{column_indices}{chromo};
}

sub start_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	if ( defined $_[0] ) {
		if ( $_[0] =~ /^\d+$/ and $_[0] <= $self->{'number_columns'} ) {
			$self->{column_indices}{start} = $_[0];
		}
		else {
			carp sprintf( "ERROR: start_column %s is not a valid index!", $_[0] );
		}
	}
	return $self->{column_indices}{start};
}

sub stop_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	if ( defined $_[0] ) {
		if ( $_[0] =~ /^\d+$/ and $_[0] <= $self->{'number_columns'} ) {
			$self->{column_indices}{stop} = $_[0];
		}
		else {
			carp sprintf( "ERROR: stop_column %s is not a valid index!", $_[0] );
		}
	}
	return $self->{column_indices}{stop};
}

*end_column = \&stop_column;
*end_column if 0;    # avoid once warning

sub coord_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	if ( defined $_[0] ) {
		if ( $_[0] =~ /^\d+$/ and $_[0] <= $self->{'number_columns'} ) {
			$self->{column_indices}{coord} = $_[0];
		}
		else {
			carp sprintf( "ERROR: coord_column %s is not a valid index!", $_[0] );
		}
	}
	return $self->{column_indices}{coord};
}

sub strand_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	if ( defined $_[0] ) {
		if ( $_[0] =~ /^\d+$/ and $_[0] <= $self->{'number_columns'} ) {
			$self->{column_indices}{strand} = $_[0];
		}
		else {
			carp sprintf( "ERROR: strand_column %s is not a valid index!", $_[0] );
		}
	}
	return $self->{column_indices}{strand};
}

sub name_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	if ( defined $_[0] ) {
		if ( $_[0] =~ /^\d+$/ and $_[0] <= $self->{'number_columns'} ) {
			$self->{column_indices}{name} = $_[0];
		}
		else {
			carp sprintf( "ERROR: name_column %s is not a valid index!", $_[0] );
		}
	}
	return $self->{column_indices}{name};
}

sub type_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	if ( defined $_[0] ) {
		if ( $_[0] =~ /^\d+$/ and $_[0] <= $self->{'number_columns'} ) {
			$self->{column_indices}{type} = $_[0];
		}
		else {
			carp sprintf( "ERROR: type_column %s is not a valid index!", $_[0] );
		}
	}
	return $self->{column_indices}{type};
}

sub id_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	if ( defined $_[0] ) {
		if ( $_[0] =~ /^\d+$/ and $_[0] <= $self->{'number_columns'} ) {
			$self->{column_indices}{id} = $_[0];
		}
		else {
			carp sprintf( "ERROR: id_column %s is not a valid index!", $_[0] );
		}
	}
	return $self->{column_indices}{id};
}

sub score_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	if ( defined $_[0] ) {
		if ( $_[0] =~ /^\d+$/ and $_[0] <= $self->{'number_columns'} ) {
			$self->{column_indices}{score} = $_[0];
		}
		else {
			carp sprintf( "ERROR: score_column %s is not a valid index!", $_[0] );
		}
	}
	return $self->{column_indices}{score};
}

*zero_start = \&interbase;
*zero_start if 0;    # avoid once warning

sub interbase {
	my $self = shift;
	if (@_) {
		my $i = $self->start_column;
		my $n = $self->name($i);
		if ( $_[0] eq '1' and lc $n eq 'start' ) {
			$self->{zerostart} = 1;
			$self->name( $i, 'Start0' );
		}
		elsif ( $_[0] eq '0' and $n eq 'start0' ) {
			$self->{zerostart} = 0;
			$self->name( $i, 'Start' );
		}
	}
	return $self->{zerostart};
}

#### Special Row methods ####

# Why is this in core and not in Data? I keep asking myself.
# Because this can get called from a Data::Feature object, which might be
# associated with a Data::Stream object. No, Stream objects don't have stored
# SeqFeatures, but I don't want the entire program to crash because of an
# undefined method because some doofus forgot. Since both Data and Stream
# objects inherit from Data::core, this is in here.
sub get_seqfeature {
	my ( $self, $row ) = @_;
	return unless ( ref($self) eq 'Bio::ToolBox::Data' );
	return unless ( $row and $row <= $self->{last_row} );
	return unless exists $self->{SeqFeatureObjects};
	return $self->{SeqFeatureObjects}->[$row] || undef;
}

1;

__END__

=head1 NAME

Bio::ToolBox::Data::core - Common functions to Bio:ToolBox::Data family

=head1 DESCRIPTION

Common methods for metadata and manipulation in a L<Bio::ToolBox::Data> 
data table and L<Bio::ToolBox::Data::Stream> file stream. This module 
should not be used directly. See the respective modules for more information.

=head1 METHODS

For quick reference only. Please see L<Bio::ToolBox::Data> for implementation.

=over 4

=item new

Generate new object. Used as a common base for L<Bio::ToolBox::Data> 
and L<Bio::ToolBox::Data::Stream>.

=item verify

Verify the integrity of the Data object. Checks multiple things, 
including metadata, table integrity (consistent number of rows and 
columns), and special file format structure.

=item open_database

This is wrapper method that tries to do the right thing and passes 
on to either L</open_meta_database> or L</open_new_database> methods. 
Basically a legacy method for L</open_meta_database>.

=item open_meta_database

Open the database that is listed in the metadata. Returns the 
database connection. Pass a true value to force a new database 
connection to be opened, rather than returning a cached connection 
object (useful when forking).

=item open_new_database

Convenience method for opening a second or new database that is 
not specified in the metadata, useful for data collection. This 
is a shortcut to L<Bio::ToolBox::db_helper/open_db_connection>.
Pass the database name.

=item verify_dataset

Verifies the existence of a dataset or data file before collecting 
data from it. Multiple datasets may be verified. This is a convenience 
method to L<Bio::ToolBox::db_helper/verify_or_request_feature_types>.
Pass the name of the dataset to verify.

=item delete_column

Delete one or more columns in a data table. Pass a list of the 
indices to delete.

=item reorder_column

Reorder the columns in a data table. Allows for skipping (deleting) and 
duplicating columns. Pass a list of the new index order.

=item feature

Returns or sets the string of the feature name listed in the metadata. 

=item feature_type

Returns "named", "coordinate", or "unknown" based on what kind of feature 
is present in the data table.

=item program

Returns or sets the program string in the metadata.

=item database

Returns or sets the name of the database in the metadata.

=item bam_adapter

Returns or sets the short name of bam adapter being used: "sam" or "hts".

=item big_adapter

Returns or sets the short name of the bigWig and bigBed adapter being 
used: "ucsc" or "big".

=item format

Returns a text string describing the format of the file contents, such 
as C<gff3>, C<gtf>, C<bed>, C<genePred>, C<narrowPeak>, etc.

=item gff

Returns or sets the GFF version value in the metadata.

=item bed

Returns or sets the number of BED columns in the metadata. 

=item ucsc

Returns or sets the number of columns in a UCSC-type file 
format, including genePred and refFlat.

=item vcf

Returns or sets the VCF version value in the metadata.

=item number_columns

Returns the number of columns in the data table.

=item number_rows

Returns the number of rows in the data table.

=item last_column

Returns the array index of the last column in the data table.

=item last_row

Returns the array index of the last row in the data table.

=item filename

Returns the complete filename listed in the metadata.

=item basename

Returns the base name of the filename listed in the metadata.

=item path

Returns the path portion of the filename listed in the metadata.

=item extension

Returns the recognized extension of the filename listed in the metadata.

=item comments

Returns an array of comment lines present in the metadata.

=item add_comment

Adds a string to the list of comments to be included in the metadata.

=item delete_comment

Deletes the indicated array index from the metadata comments array.

=item vcf_headers

Partially parses VCF metadata header lines into a hash structure.

=item rewrite_vcf_headers

Rewrites the vcf headers back into the metadata comments array.

=item list_columns

Returns an array of the column names

=item name

Returns or sets the name of the column. Pass the index, and optionally new name.

=item metadata($index, $key)

Returns or sets the metadata keyE<sol>value pair for a specific column.
Pass the index, key, and optionally new value.

=item delete_metadata

Deletes the metadata key for a column. Pass the index and key.

=item copy_metadata

Copies the metadata values from one column to another column.
Pass the source and target indices.

=item find_column

Returns the column index for the column with the specified name. Name 
searches are case insensitive and can tolerate a # prefix character. 
The first match is returned. Pass the name to search.

=item chromo_column

Returns the index of the column that best represents the chromosome column.

=item start_column

Returns the index of the column that best represents the start, position, 
or transcription start column.

=item stop_column

=item end_column

Returns the index of the column that best represents the stop or end 
column. 

=item strand_column

Returns the index of the column that best represents the strand.

=item name_column

Returns the index of the column that best represents the name.

=item type_column

Returns the index of the column that best represents the type.

=item id_column

Returns the index of the column that represents the Primary_ID 
column used in databases.

=item score_column

Returns the index of the column that represents the Score 
column in certain formats, such as GFF, BED, bedGraph, etc.

=item zero_start

=item interbase

Returns true (1) or false (0) if the coordinate system appears to 
be an interbase, half-open, or zero-based coordinate system. This is 
based on file type, e.g. F<.bed>, or if the start coordinate column 
name is C<start0>. The coordinate system can also be explicitly changed 
by passing an appropriate value; note that this will also change the 
start coordinate column name as appropriate. 

=item get_seqfeature

Returns the stored SeqFeature object for a given row.

=back

=head1 SEE ALSO

L<Bio::ToolBox::Data>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
