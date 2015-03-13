package Bio::ToolBox::big_helper;

### modules
require Exporter;
use strict;
use Carp qw(carp cluck);
use File::Temp;
use Bio::ToolBox::db_helper qw(get_chromosome_list);
use Bio::ToolBox::db_helper::config qw($BTB_CONFIG add_program);



### Export
our @ISA = qw(Exporter);
our @EXPORT = qw(
);
our @EXPORT_OK = qw(
	wig_to_bigwig_conversion
	bed_to_bigbed_conversion
	generate_chromosome_file
);


our $VERSION = '1.14';

1;


### Wig to BigWig file conversion
sub wig_to_bigwig_conversion {
	
	# Collect passed arguments
	my %args = @_; 
	unless (%args) {
		cluck "no arguments passed!";
		return;
	}
	
	# wigfile
	$args{'wig'} ||= undef;
	unless ($args{'wig'}) {
		cluck "no wig file passed!";
		return;
	}
	
	
	# Identify bigwig conversion utility
	my $utility = $args{'wig'} =~ /\.(?:bdg|bedgraph|bed)$/i ? 'bedGraphToBigWig' :
				'wigToBigWig';
	$args{'bwapppath'} ||= undef;
	unless ($args{'bwapppath'}) {
		# check for an entry in the configuration file
		$args{'bwapppath'} = $BTB_CONFIG->param("applications.$utility") || 
				undef;
	}
	unless ($args{'bwapppath'}) {
		# try checking the system path as a final resort
		eval {
			require File::Which;
			File::Which->import;
			$args{'bwapppath'} = which($utility);
		};
		add_program($args{'bwapppath'}) if $args{'bwapppath'};
	}
	unless ($args{'bwapppath'}) {
		carp " Utility '$utility' not specified and can not be found!" . 
			" Conversion failed!\n";
		return;
	}
	
	# verify utility
	unless ($args{'bwapppath'} =~ /$utility\Z/) {
		carp " Wrong utility for the type of wig file! Unable to convert!\n";
		return;
	}
	
	# Generate list of chromosome sizes if necessary
	$args{'chromo'} ||= undef;
	unless ($args{'chromo'}) {
		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{'db'} ||= undef;
		unless ($args{'db'}) {
			carp " No requisite database or chromosome info file provided!" .
				" Conversion failed\n";
			return;
		}
		$args{'chromo'} = generate_chromosome_file($args{'db'});
		unless ($args{'chromo'}) {
			carp " Cannot generate chromosome info file! Conversion failed\n";
			return;
		}
	}
	
	
	# Generate the bw file name
	# we can substitute one of three possible names for bw
	my $bw_file = $args{'wig'};
	$bw_file =~ s/\.(?:bed|bdg|bedgraph|wig)$/.bw/;
	
	
	# Generate the bigwig file 
	print " converting $args{'wig'} to bigWig....\n";
	if ($args{'bwapppath'} =~ /wigToBigWig$/) {
		# include the -clip option in case there are any positions 
		# out of bounds of the chromosome
		# it will just warn instead of fail
		system $args{'bwapppath'}, '-clip', $args{'wig'}, $args{'chromo'}, $bw_file;
	}
	elsif ($args{'bwapppath'} =~ /bedGraphToBigWig$/) {
		# this doesn't have the -clip option, too bad
		system $args{'bwapppath'}, $args{'wig'}, $args{'chromo'}, $bw_file;
	}
	
	# check the result
	if (-e $bw_file and -s $bw_file) {
		# conversion successful
		if ($args{'chromo'} =~ /^chr_sizes_\w{5}/) {
			# we no longer need our temp chromosome file
			unlink $args{'chromo'};
		}
		return $bw_file;
	}
	else {
		warn " Conversion failed. You should try manually and watch for errors\n";
		if (-e $bw_file) {
			# 0-byte file was created
			unlink $bw_file;
		}
		if ($args{'chromo'} =~ /^chr_sizes_\w{5}/) {
			# leave the temp chromosome file as a courtesy
			warn " Leaving temporary chromosome file '$args{'chromo'}'\n";
		}
		return;
	}
}


### Bed to BigBed file conversion
sub bed_to_bigbed_conversion {
	
	# Collect passed arguments
	my %args = @_; 
	unless (%args) {
		cluck "no arguments passed!";
		return;
	}
	
	# bedfile
	$args{'bed'} ||= undef;
	unless ($args{'bed'}) {
		carp "no bed file passed!";
		return;
	}
	
	
	# identify bigbed conversion utility
	$args{'bbapppath'} ||= undef;
	unless ($args{'bbapppath'}) {
		# check for an entry in the configuration file
		$args{'bbapppath'} = $BTB_CONFIG->param('applications.bedToBigBed') || 
			undef;
	}
	unless ($args{'bbapppath'}) {
		# try checking the system path as a final resort
		eval {
			require File::Which;
			File::Which->import;
			$args{'bbapppath'} = which('bedToBigBed');
		};
		add_program($args{'bbapppath'}) if $args{'bbapppath'};
	}
	unless ($args{'bbapppath'}) {
		carp " Utility 'bedToBigBed' not specified and can not be found!" . 
			" Conversion failed!\n";
		return;
	}
	
	
	# Generate list of chromosome sizes if necessary
	$args{'chromo'} ||= undef;
	unless ($args{'chromo'}) {
		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{'db'} ||= undef;
		unless ($args{'db'}) {
			carp " No requisite database or chromosome info file provided!" .
				" Conversion failed\n";
			return;
		}
		$args{'chromo'} = generate_chromosome_file($args{'db'});
		unless ($args{'chromo'}) {
			carp " Cannot generate chromosome info file! Conversion failed\n";
			return;
		}
	}
	
	
	# Generate the bb file name
	my $bb_file = $args{'bed'};
	$bb_file =~ s/\.bed$/.bb/;
	
	
	# Generate the bigBed file using Jim Kent's utility
	print " converting $args{'bed'} to BigBed....\n";
	system $args{'bbapppath'}, $args{'bed'}, $args{'chromo'}, $bb_file;
	
	
	# Check the result
	if (-e $bb_file and -s $bb_file) {
		# conversion successful
		if ($args{'chromo'} =~ /^chr_sizes_\w{5}/) {
			# we no longer need our temp chromosome file
			unlink $args{'chromo'};
		}
		return $bb_file;
	}
	else {
		warn " Conversion failed. You should try manually and watch for errors\n";
		if (-e $bb_file) {
			# 0-byte file was created
			unlink $bb_file;
		}
		if ($args{'chromo'} =~ /^chr_sizes_\w{5}/) {
			# leave the temp chromosome file as a courtesy
			warn " Leaving temporary chromosome file '$args{'chromo'}'\n";
		}
		return;
	}
}



sub generate_chromosome_file {
	
	my $database = shift;
	print " generating chromosome file....\n";
	
	# generate chromosome lengths file
	my @chromosomes = get_chromosome_list($database);
	unless (@chromosomes) {
		carp " no chromosome sequences identified in database!\n";
		return;
	}
	
	# prepare temp file
	my $chr_fh = new File::Temp(
		'UNLINK'   => 0,
		'TEMPLATE' => 'chr_sizes_XXXXX',
	);
	my $chromo_file = $chr_fh->filename;

	# write out
	foreach my $chr (@chromosomes) {
		# chromosome name and size
		$chr_fh->print( $chr->[0] . "\t" . $chr->[1] . "\n");
	}
	$chr_fh->close;
	
	return $chromo_file;
}




__END__

=head1 NAME

Bio::ToolBox::big_helper

=head1 DESCRIPTION

This module helps in the conversion of wig and bed files to bigWig and 
bigBed files, respectively. It uses external applications to 
accomplish this, taking care of generating a chromosome file from a 
database if necessary. 

Two exported subroutines are available for wig and bed conversions. 

=head1 USAGE

Load the module at the beginning of your program and include the name or 
names of the subroutines to export. None are automatically exported.

	use Bio::ToolBox::big_helper qw(wig_to_bigwig_conversion);


=over

=item wig_to_bigwig_conversion()

This subroutine will convert a wig file to a bigWig file. See the UCSC 
documentation regarding wig (http://genome.ucsc.edu/goldenPath/help/wiggle.html)
and bigWig (http://genome.ucsc.edu/goldenPath/help/bigWig.html) file formats. 
It uses Jim Kent's wigToBigWig or bedGraphToBigWig utility to perform the 
conversion, depending on the format of the wig file. The utility must be 
available on the system for the conversion to succeed. 

The conversion requires a list of chromosome name and sizes in a simple text 
file, where each line is comprised of two columns, "<chromosome name> 
<size in bases>". This file may be specified, or automatically generated if 
given a Bio::DB database name (preferred to ensure genome version 
compatibility).

After running the utility, the existence of a non-zero byte bigWig file 
is checked. If it does, then the name of the file is returned. If not, 
an error is printed and nothing is returned. 

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

After running the utility, the existence of a non-zero byte bigBed file 
is checked. If it does, then the name of the file is returned. If not, 
an error is printed and nothing is returned. 

Pass the function an anonymous hash of arguments, including the following:

  Required:
  bed         => The name of the bed source file. 
  db          => Provide an opened database object from which to generate 
                 the chromosome sizes information.
  Optional: 
  chromo      => The name of the chromosome sizes text file, described 
                 above, as an alternative to providing the database name.
  bbapppath   => Provide the full path to Jim Kent's bedToBigBed  
                 utility. This parameter may instead be defined in the 
                 configuration file "biotoolbox.cfg". 

Example

	my $bed_file = 'example.bed';
	my $bb_file = bed_to_bigbed_conversion( {
			'bed'   => $bed_file,
			'db'    => $database,
	} );
	if ($bb_file) {
		print " success! wrote bigBed file $bb_file\n";
	}
	else {
		print " failure! see STDERR for errors\n";
	};

=item generate_chromosome_file

This subroutine will generate a chromosome sizes files appropriate for 
the big file conversion utilities from an available database. It is a 
two column text file, the first column is the chromosome name, and the 
second column is the length in bp. The file is written in the 
current directory with a name of "chr_sizesXXXXX", where X are random 
characters as defined by File::Temp. 

The chromosome names and lengths are obtained from a Bio::DB 
database using the C<Bio::ToolBox::db_helper::get_chromosome_list()> 
subroutine.

Pass the subroutine a database name, path to a supported database file, 
or opened Bio::DB object.

The file will be written, closed, and the filename returned.

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  


