package Bio::ToolBox::big_helper;
our $VERSION = '1.54';

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

=item wig_to_bigwig_conversion

This subroutine will convert a wig file to a bigWig file. See the UCSC 
documentation regarding L<wig|http://genome.ucsc.edu/goldenPath/help/wiggle.html>
and L<bigWig|http://genome.ucsc.edu/goldenPath/help/bigWig.html> file formats. 
It uses the UCSC I<wigToBigWig> utility to perform the conversion. The utility 
must be available on the system for the conversion to succeed. 

For bedGraph format wig files, the utility I<bedGraphToBigWig> may be substituted 
if desired, but I<wigToBigWig> can sufficiently handle all wig formats. When 
no utility is available but L<Bio::DB::BigFile> is installed, then the module 
may be used for generating the bigWig file.

The conversion requires a list of chromosome name and sizes in a simple text 
file, where each line is comprised of two columns, "C<chromosome_name> 
<size_in_bases>". This file may be specified, or automatically generated if 
given a C<Bio::DB> database name (preferred to ensure genome version 
compatibility).

After running the utility, the existence of a non-zero byte bigWig file 
is checked. If it does, then the name of the file is returned. If not, 
an error is printed and nothing is returned. 

Pass the function an array of key =E<gt> value arguments, including the 
following:

  Required:
  wig         => The name of the wig source file. 
  db          => Provide an opened database object from which to generate 
                 the chromosome sizes information.
  Optional: 
  chromo      => The name of the chromosome sizes text file, described 
                 above, as an alternative to providing the database name.
  bwapppath   => Provide the full path to Jim Kent's wigToBigWig 
                 utility. This parameter may instead be defined in the 
                 configuration file C<biotoolbox.cfg>. 

Example

	my $wig_file = 'example_wig';
	my $bw_file = wig_to_bigwig_conversion(
		'wig'   => $wig_file,
		'db'    => $database,
	);
	if (-e $bw_file) {
		print " success! wrote bigwig file $bw_file\n";
		unlink $wig_file; # no longer necessary
	}
	else {
		print " failure! see STDERR for errors\n";
	};

=item open_wig_to_bigwig_fh

This subroutine will open a forked process to the UCSC I<wigToBigWig> utility 
as a file handle, allowing wig lines to be "printed" to the utility for 
conversion. This is useful for writing directly to a bigWig file without 
having to write a temporary wig file first. This is also useful when you 
have multiple wig files, for example individual wig files from separate 
forked processes, that need to be combined into a bigWig file. 

Note that the I<wigToBigWig> utility does not handle errors gracefully 
and will immediately fail upon encountering errors, usually also bringing 
the main Perl process with it. Make sure the chromosome file is accurate 
and the wig lines are properly formatted and in order! 

Pass the function an array of key =E<gt> value arguments. An L<IO::File> 
object will be returned. Upon the closing the file handle, the 
I<wigToBigWig> utility will generate the bigWig file.

  Required:
  bw          => The output file name for the bigWig file.
                 Also accepts the keys file, wig, and out. 
  chromo      => The name of the chromosome sizes text file, described 
                 in wig_to_bigwig_conversion()
  Optional: 
  db          => Alternatively, provide an opened database object from which 
                 to generate a temporary chromosome sizes file. It is up to the 
                 user to delete this file.
  bwapppath   => Provide the full path to the UCSC I<wigToBigWig>utility. 
                 The path may be obtained from the configuration file 
                 F<.biotoolbox.cfg>. 

  Example:
	my $bw_file = 'example.bw';
	my $chromo_file = generate_chromosome_file($db);
	my $bwfh = open_wig_to_bigwig_fh(
		file    => $bw_file,
		chromo  => $chromo_file,
	);
	foreach (@wig_lines) {
		$bwfh->print("$_\n");
	}
	$bwfh->close;
		# this signals the forked wigToBigWig process to write 
		# the bigWig file, which may take a few seconds to minutes
	unlink $chromo_file;

=item bed_to_bigbed_conversion

This subroutine will convert a bed file to a bigBed file. See the UCSC 
documentation regarding L<bed|http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED>
and L<bigBed|http://genome.ucsc.edu/goldenPath/help/bigBed.html> file formats. 
It uses the UCSC I<bedToBigBed> utility to perform the conversion. This 
must be present on the system for the conversion to succeed. 

The conversion requires a list of chromosome name and sizes in a simple text 
file, where each line is comprised of two columns, "C<chromosome_name> 
C<size_in_bases>". This file may be specified, or automatically generated if 
given a C<Bio::DB> database name (preferred to ensure genome version 
compatibility).

After running the utility, the existence of a non-zero byte bigBed file 
is checked. If it does, then the name of the file is returned. If not, 
an error is printed and nothing is returned. 

Pass the function an array of key =E<gt> value arguments, including the 
following:

  Required:
  bed         => The name of the bed source file. 
  db          => Provide an opened database object from which to generate 
                 the chromosome sizes information.
  Optional: 
  chromo      => The name of the chromosome sizes text file, described 
                 above, as an alternative to providing the database name.
  bbapppath   => Provide the full path to the UCSC bedToBigBed  
                 utility. This parameter may instead be defined in the 
                 configuration file "biotoolbox.cfg". 

Example

	my $bed_file = 'example.bed';
	my $bb_file = bed_to_bigbed_conversion(
		'bed'   => $bed_file,
		'db'    => $database,
	);
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
current directory with a name of F<chr_sizesXXXXX>, where X are random 
characters as defined by L<File::Temp>. 

The chromosome names and lengths are obtained from a C<Bio::DB> 
database using the L<Bio::ToolBox::db_helper/get_chromosome_list> 
subroutine.

Pass the subroutine a database name, path to a supported database file, 
or opened C<Bio::DB> object. 

Optionally pass a second value, a regular expression compatible string 
or C<qr> for skipping specific chromosomes or chromosome classes, such as 
mitochondrial or unmapped contigs. The default is to return all chromosomes. 

The file will be written, closed, and the filename returned.

=back

=cut


### modules
require Exporter;
use strict;
use Carp qw(carp cluck);
use File::Temp;
use IO::File;
use Bio::ToolBox::db_helper qw(get_chromosome_list);
use Bio::ToolBox::db_helper::config qw($BTB_CONFIG add_program);



### Export
our @ISA = qw(Exporter);
our @EXPORT = qw(
);
our @EXPORT_OK = qw(
	wig_to_bigwig_conversion
	open_wig_to_bigwig_fh
	bed_to_bigbed_conversion
	generate_chromosome_file
);


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
	$args{wig} ||= undef;
	unless ($args{wig}) {
		cluck "no wig file passed!";
		return;
	}
	
	
	# Identify bigwig conversion utility
	$args{bwapppath} ||= undef;
	unless ($args{bwapppath}) {
		# check for an entry in the configuration file
		$args{bwapppath} = $BTB_CONFIG->param("applications.wigToBigWig") || 
				undef;
	}
	unless ($args{bwapppath}) {
		# try checking the system path as a final resort
		eval {
			require File::Which;
			File::Which->import;
			$args{bwapppath} = which('wigToBigWig');
		};
		add_program($args{bwapppath}) if $args{bwapppath};
	}
	unless ($args{bwapppath}) {
		# last attempt to use Bio::DB::BigFile
		# this is not a recommended method, and produces different sized 
		# bw files, I think because of different internal settings
		# this may be deprecated
		eval {
			require Bio::DB::BigFile;
			if (Bio::DB::BigFile->can('createBigWig')) {
				$args{bwapppath} = 'BioDBBigFile';
			}
		};
	}
	unless ($args{bwapppath}) {
		warn " Utility 'wigToBigWig' not specified and can not be found!" . 
			" Conversion failed!\n";
		return;
	}
	
	
	# Generate list of chromosome sizes if necessary
	$args{chromo} ||= undef;
	unless ($args{chromo}) {
		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{db} ||= undef;
		unless ($args{db}) {
			carp " No requisite database or chromosome info file provided!" .
				" Conversion failed\n";
			return;
		}
		$args{chromo} = generate_chromosome_file($args{db});
		unless ($args{chromo}) {
			carp " Cannot generate chromosome info file! Conversion failed\n";
			return;
		}
	}
	
	
	# Generate the bw file name
	# we can substitute one of three possible names for bw
	my $bw_file = $args{wig};
	$bw_file =~ s/\.(?:bed|bdg|bedgraph|wig)$/.bw/;
	
	
	# Generate the bigwig file 
	printf " converting %s to bigWig....\n", $args{wig};
	if ($args{bwapppath} =~ /wigToBigWig$/) {
		# include the -clip option in case there are any positions 
		# out of bounds of the chromosome
		# it will just warn instead of fail
		system $args{bwapppath}, '-clip', $args{wig}, $args{chromo}, $bw_file;
	}
	elsif ($args{bwapppath} =~ /bedGraphToBigWig$/) {
		# this doesn't have the -clip option, too bad
		system $args{bwapppath}, $args{wig}, $args{chromo}, $bw_file;
	}
	elsif ($args{bwapppath} eq 'BioDBBigFile') {
		Bio::DB::BigFile->createBigWig($args{wig}, $args{chromo}, $bw_file);
	}
	
	# check the result
	if (-e $bw_file and -s $bw_file) {
		# conversion successful
		if ($args{chromo} =~ /^chr_sizes_\w{5}/) {
			# we no longer need our temp chromosome file
			unlink $args{chromo};
		}
		return $bw_file;
	}
	else {
		warn " Conversion failed. You should try manually and watch for errors\n";
		if (-e $bw_file) {
			# 0-byte file was created
			unlink $bw_file;
		}
		if ($args{chromo} =~ /^chr_sizes_\w{5}/) {
			# leave the temp chromosome file as a courtesy
			warn " Leaving temporary chromosome file '$args{chromo}'\n";
		}
		return;
	}
}


### Open a file handle to wigToBigWig
sub open_wig_to_bigwig_fh {
	
	# Collect passed arguments
	my %args = @_; 
	unless (%args) {
		cluck "no arguments passed!";
		return;
	}
	
	# bigWig output file
	$args{bw} ||= $args{wig} || $args{out} || $args{file} || undef;
	unless ($args{bw}) {
		cluck "no output bw file name passed!";
		return;
	}
	unless ($args{bw} =~ /\.bw$/i) {
		$args{bw} .= '.bw';
	}
	
	# Identify bigwig conversion utility
	$args{bwapppath} ||= undef;
	unless ($args{bwapppath}) {
		# check for an entry in the configuration file
		$args{bwapppath} = $BTB_CONFIG->param("applications.wigToBigWig") || undef;
	}
	unless ($args{bwapppath}) {
		# try checking the system path as a final resort
		eval {
			require File::Which;
			File::Which->import;
			$args{bwapppath} = which('wigToBigWig');
		};
		add_program($args{bwapppath}) if $args{bwapppath};
	}
	unless ($args{bwapppath} =~ /ToBigWig$/) {
		warn " Utility 'wigToBigWig' not specified and can not be found!" . 
			" Conversion failed!\n";
		return;
	}

	# Generate list of chromosome sizes if necessary
	$args{chromo} ||= undef;
	unless ($args{chromo}) {
		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{db} ||= undef;
		unless ($args{db}) {
			carp " No requisite database or chromosome info file provided!" .
				" Conversion failed\n";
			return;
		}
		$args{chromo} = generate_chromosome_file($args{db});
		unless ($args{chromo}) {
			carp " Cannot generate chromosome info file! Conversion failed\n";
			return;
		}
	}
	
	# open the filehandle
	my $command = sprintf "%s stdin %s %s", $args{bwapppath}, $args{chromo}, 
		$args{bw};
	my $bwfh = IO::File->new("| $command") or 
		die sprintf("cannot open %s!\n", $args{bwapppath});
		# wigToBigWig will always die anyway if something is wrong
		# cannot trap it with an eval, since it doesn't just error out
		# but actually exits, dragging the whole Perl process with it
	return $bwfh;
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
	$args{bed} ||= undef;
	unless ($args{bed}) {
		carp "no bed file passed!";
		return;
	}
	
	
	# identify bigbed conversion utility
	$args{bbapppath} ||= undef;
	unless ($args{bbapppath}) {
		# check for an entry in the configuration file
		$args{bbapppath} = $BTB_CONFIG->param('applications.bedToBigBed') || 
			undef;
	}
	unless ($args{bbapppath}) {
		# try checking the system path as a final resort
		eval {
			require File::Which;
			File::Which->import;
			$args{bbapppath} = which('bedToBigBed');
		};
		add_program($args{bbapppath}) if $args{bbapppath};
	}
	unless ($args{bbapppath}) {
		carp " Utility 'bedToBigBed' not specified and can not be found!" . 
			" Conversion failed!\n";
		return;
	}
	
	
	# Generate list of chromosome sizes if necessary
	$args{chromo} ||= undef;
	unless ($args{chromo}) {
		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{db} ||= undef;
		unless ($args{db}) {
			carp " No requisite database or chromosome info file provided!" .
				" Conversion failed\n";
			return;
		}
		$args{chromo} = generate_chromosome_file($args{db});
		unless ($args{chromo}) {
			carp " Cannot generate chromosome info file! Conversion failed\n";
			return;
		}
	}
	
	
	# Generate the bb file name
	my $bb_file = $args{bed};
	$bb_file =~ s/\.bed$/.bb/;
	
	
	# Generate the bigBed file using Jim Kent's utility
	print " converting $args{bed} to BigBed....\n";
	system $args{bbapppath}, $args{bed}, $args{chromo}, $bb_file;
	
	
	# Check the result
	if (-e $bb_file and -s $bb_file) {
		# conversion successful
		if ($args{chromo} =~ /^chr_sizes_\w{5}/) {
			# we no longer need our temp chromosome file
			unlink $args{chromo};
		}
		return $bb_file;
	}
	else {
		warn " Conversion failed. You should try manually and watch for errors\n";
		if (-e $bb_file) {
			# 0-byte file was created
			unlink $bb_file;
		}
		if ($args{chromo} =~ /^chr_sizes_\w{5}/) {
			# leave the temp chromosome file as a courtesy
			warn " Leaving temporary chromosome file '$args{chromo}'\n";
		}
		return;
	}
}



sub generate_chromosome_file {
	
	my $database = shift;
	my $chr_exclude = shift || undef;
	print " generating chromosome file....\n";
	
	# generate chromosome lengths file
	my @chromosomes = get_chromosome_list($database, $chr_exclude);
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

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  


