package Bio::ToolBox::big_helper;

use warnings;
use strict;
use Carp qw(carp cluck confess);
use File::Temp;
use File::Which;
use IO::File;
use Bio::ToolBox::db_helper qw(get_chromosome_list);
use Bio::ToolBox::db_helper::config qw($BTB_CONFIG);
require Exporter;

our $VERSION = '1.68';

### Export
our @ISA    = qw(Exporter);
our @EXPORT_OK = qw(
	wig_to_bigwig_conversion
	open_wig_to_bigwig_fh
	open_bigwig_to_wig_fh
	bed_to_bigbed_conversion
	generate_chromosome_file
);

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
	unless ( $args{wig} ) {
		cluck "no wig file passed!";
		return;
	}

	# Identify bigwig conversion utility
	$args{bwapppath} ||= undef;
	unless ( $args{bwapppath} ) {

		# check for an entry in the configuration file
		$args{bwapppath} = $BTB_CONFIG->param("applications.wigToBigWig")
			|| undef;
	}
	unless ( $args{bwapppath} ) {

		# try checking the system path
		$args{bwapppath} = which('wigToBigWig');
	}
	unless ( $args{bwapppath} ) {

		# last attempt to use Bio::DB::BigFile
		# this is not a recommended method, and produces different sized
		# bw files, I think because of different internal settings
		# this may be deprecated
		eval {
			require Bio::DB::BigFile;
			if ( Bio::DB::BigFile->can('createBigWig') ) {
				$args{bwapppath} = 'BioDBBigFile';
			}
		};
	}
	unless ( $args{bwapppath} ) {
		warn " Utility 'wigToBigWig' not specified and can not be found!"
			. " Conversion failed!\n";
		return;
	}

	# Generate list of chromosome sizes if necessary
	$args{chromo} ||= undef;
	unless ( $args{chromo} ) {

		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{db} ||= undef;
		unless ( $args{db} ) {
			carp " No requisite database or chromosome info file provided!"
				. " Conversion failed\n";
			return;
		}
		$args{chromo} = generate_chromosome_file( $args{db} );
		unless ( $args{chromo} ) {
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
	if ( $args{bwapppath} =~ /wigToBigWig$/ ) {

		# include the -clip option in case there are any positions
		# out of bounds of the chromosome
		# it will just warn instead of fail
		system $args{bwapppath}, '-clip', $args{wig}, $args{chromo}, $bw_file;
	}
	elsif ( $args{bwapppath} =~ /bedGraphToBigWig$/ ) {

		# this doesn't have the -clip option, too bad
		system $args{bwapppath}, $args{wig}, $args{chromo}, $bw_file;
	}
	elsif ( $args{bwapppath} eq 'BioDBBigFile' ) {
		Bio::DB::BigFile->createBigWig( $args{wig}, $args{chromo}, $bw_file );
	}

	# check the result
	if ( -e $bw_file and -s $bw_file ) {

		# conversion successful
		if ( $args{chromo} =~ /^chr_sizes_\w{5}/ ) {

			# we no longer need our temp chromosome file
			unlink $args{chromo};
		}
		return $bw_file;
	}
	else {
		warn " Conversion failed. You should try manually and watch for errors\n";
		if ( -e $bw_file ) {

			# 0-byte file was created
			unlink $bw_file;
		}
		if ( $args{chromo} =~ /^chr_sizes_\w{5}/ ) {

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
	unless ( $args{bw} ) {
		cluck "no output bw file name passed!";
		return;
	}
	unless ( $args{bw} =~ /\.bw$/i ) {
		$args{bw} .= '.bw';
	}

	# Identify bigwig conversion utility
	$args{bwapppath} ||= undef;
	unless ( $args{bwapppath} ) {

		# check for an entry in the configuration file
		$args{bwapppath} = $BTB_CONFIG->param("applications.wigToBigWig") || undef;
	}
	unless ( $args{bwapppath} ) {

		# try checking the system path
		$args{bwapppath} = which('wigToBigWig');
	}
	unless ( $args{bwapppath} =~ /ToBigWig$/ ) {
		carp " Utility 'wigToBigWig' not specified and can not be found!\n";
		return;
	}

	# Generate list of chromosome sizes if necessary
	$args{chromo} ||= undef;
	unless ( $args{chromo} ) {

		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{db} ||= undef;
		unless ( $args{db} ) {
			cluck " No requisite database or chromosome info file provided!"
				. " Conversion failed\n";
			return;
		}
		$args{chrskip} ||= undef;
		$args{chromo} = generate_chromosome_file( $args{db}, $args{chrskip} );
		unless ( $args{chromo} ) {
			cluck " Cannot generate chromosome info file! Conversion failed\n";
			return;
		}
	}

	# open the filehandle
	my $command = sprintf "%s stdin %s %s", $args{bwapppath}, $args{chromo}, $args{bw};
	my $bwfh    = IO::File->new("| $command")
		or confess sprintf( "cannot open %s!\n", $args{bwapppath} );

	# wigToBigWig will always die anyway if something is wrong
	# cannot trap it with an eval, since it doesn't just error out
	# but actually exits, dragging the whole Perl process with it
	# printf "we have a filehandle %s\n", ref $bwfh;
	confess "unable to execute command '$command'" unless ref $bwfh;

	# we will still get an IO::File handle back even with a failed convertor - sigh
	return $bwfh;
}

### Open a file handle from bigWigToWig
sub open_bigwig_to_wig_fh {

	# Collect passed arguments
	my %args = @_;
	unless (%args) {
		cluck "no arguments passed!";
		return;
	}

	# bigWig output file
	$args{bw} ||= $args{wig} || $args{file} || undef;
	unless ( $args{bw} ) {
		cluck "no input bw file name passed!";
		return;
	}
	unless ( $args{bw} =~ /\.bw$/i ) {
		$args{bw} .= '.bw';
	}

	# Identify bigwig conversion utility
	$args{bwapppath} ||= undef;
	unless ( $args{bwapppath} ) {

		# check for an entry in the configuration file
		$args{bwapppath} = $BTB_CONFIG->param("applications.bigWigToWig") || undef;
	}
	unless ( $args{bwapppath} ) {

		# try checking the system path
		$args{bwapppath} = which('bigWigToWig') || which('bigWigToBedGraph');
	}
	unless ( $args{bwapppath} =~ /bigWigTo(?:Wig|BedGraph)$/ ) {
		carp " Utility 'bigWigToWig' not specified and can not be found!\n";
		return;
	}

	# open the filehandle
	my $command = sprintf "%s %s stdout", $args{bwapppath}, $args{bw};
	my $bwfh    = IO::File->new("$command |")
		or confess sprintf( "cannot open %s!\n", $args{bwapppath} );

	# bigWigToWig will always die anyway if something is wrong
	# cannot trap it with an eval, since it doesn't just error out
	# but actually exits, dragging the whole Perl process with it
	# printf "we have a filehandle %s\n", ref $bwfh;
	confess "unable to execute command '$command'" unless ref $bwfh;

	# we will still get an IO::File handle back even with a failed convertor - sigh
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
	unless ( $args{bed} ) {
		cluck "no bed file passed!";
		return;
	}

	# identify bigbed conversion utility
	$args{bbapppath} ||= undef;
	unless ( $args{bbapppath} ) {

		# check for an entry in the configuration file
		$args{bbapppath} = $BTB_CONFIG->param('applications.bedToBigBed')
			|| undef;
	}
	unless ( $args{bbapppath} ) {

		# try checking the system path
		$args{bbapppath} = which('bedToBigBed');
	}
	unless ( $args{bbapppath} ) {
		carp " Utility 'bedToBigBed' not specified and can not be found!\n";
		return;
	}

	# Generate list of chromosome sizes if necessary
	$args{chromo} ||= undef;
	unless ( $args{chromo} ) {

		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{db} ||= undef;
		unless ( $args{db} ) {
			cluck " No requisite database or chromosome info file provided!\n";
			return;
		}
		$args{chromo} = generate_chromosome_file( $args{db} );
		unless ( $args{chromo} ) {
			carp " Cannot generate chromosome info file!\n";
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
	if ( -e $bb_file and -s $bb_file ) {

		# conversion successful
		if ( $args{chromo} =~ /^chr_sizes_\w{5}/ ) {

			# we no longer need our temp chromosome file
			unlink $args{chromo};
		}
		return $bb_file;
	}
	else {
		carp " Conversion failed. You should try manually and watch for errors\n";
		if ( -e $bb_file ) {

			# 0-byte file was created
			unlink $bb_file;
		}
		if ( $args{chromo} =~ /^chr_sizes_\w{5}/ ) {

			# leave the temp chromosome file as a courtesy
			print STDERR " Leaving temporary chromosome file '$args{chromo}'\n";
		}
		return;
	}
}

sub generate_chromosome_file {

	my $database    = shift;
	my $chr_exclude = shift || undef;
	print " generating chromosome file....\n";

	# generate chromosome lengths file
	my @chromosomes = get_chromosome_list( $database, $chr_exclude );
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
		$chr_fh->printf( "%s\t%d\n", $chr->[0], $chr->[1] );
	}
	$chr_fh->close;

	return $chromo_file;
}

1;

__END__

=head1 NAME

Bio::ToolBox::big_helper

=head1 DESCRIPTION

This module helps in the conversion of wig and bed files to bigWig and 
bigBed files, respectively. It uses external applications to 
accomplish this, taking care of generating a chromosome file from a 
database if necessary. 

For wig to bigWig conversion, see the UCSC 
documentation regarding L<wig|http://genome.ucsc.edu/goldenPath/help/wiggle.html>
and L<bigWig|http://genome.ucsc.edu/goldenPath/help/bigWig.html> file formats. 
It uses the UCSC I<wigToBigWig> utility to perform the conversion. The utility 
must be available on the system for the conversion to succeed. It may be downloaded 
from the L<UCSC Genome utilities|http://hgdownload.soe.ucsc.edu/admin/exe/> page.

For bed to bigBed conversion, See the UCSC 
documentation regarding L<bed|http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED>
and L<bigBed|http://genome.ucsc.edu/goldenPath/help/bigBed.html> file formats. 
It uses the UCSC I<bedToBigBed> utility to perform the conversion. This 
must be present on the system for the conversion to succeed. It may be downloaded 
from the L<UCSC Genome utilities|http://hgdownload.soe.ucsc.edu/admin/exe/> page.

In both cases, the conversion requires a list of chromosome name and sizes in a 
simple text file, where each line is comprised of two columns, "C<chromosome_name> 
<size_in_bases>". This file may be specified, or automatically generated if 
given a C<Bio::DB> database name (preferred to ensure genome version 
compatibility).

=head1 USAGE

Load the module at the beginning of your program and include the name or 
names of the subroutines to export. None are automatically exported.

	use Bio::ToolBox::big_helper qw(wig_to_bigwig_conversion);

There are are five available exported subroutines.

=over 4

=item wig_to_bigwig_conversion

This subroutine will convert a wig file to a bigWig file. 

For bedGraph format wig files, the utility I<bedGraphToBigWig> may be substituted 
if desired, but I<wigToBigWig> can sufficiently handle all wig formats. When 
no utility is available but L<Bio::DB::BigFile> is installed, then the module 
may be used for generating the bigWig file.

After running the utility, the existence of a non-zero byte bigWig file 
is checked. If it does, then the name of the file is returned. If not, 
an error is printed and nothing is returned. 

Pass the function an array of key =E<gt> value arguments, including the 
following:

=over 4

=item wig

Pass the name of the wig source file. This is required.

=item chromo

Pass the path to a chromosome sizes text file, described as above. This is 
required. Alternatively, a database object may provided.

=item db

Provide an opened database object from which to generate the chromosome 
sizes information. This will be passed to L<generate_chromosome_file> to 
generate a chromosome file. This is a convenience option and alternative to 
providing an existing chromosome file.

=item bwapppath

Provide the full path to the UCSC F<wigToBigWig> utility. If not provided, the 
default C<PATH> will be searched for the utility. The path may also 
be defined in the configuration file F<.biotoolbox.cfg>. 

=back

Example:

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

This subroutine will open a forked process to the UCSC F<wigToBigWig> utility 
as a file handle, allowing wig lines to be "printed" to the utility for 
conversion. This is useful for writing directly to a bigWig file without 
having to write a temporary wig file first. This is also useful when you 
have multiple wig files, for example individual wig files from separate 
forked processes, that need to be combined into a bigWig file. 

Note that the F<wigToBigWig> utility does not handle errors gracefully 
and will immediately fail upon encountering errors, usually also bringing 
the main Perl process with it. Make sure the chromosome file is accurate 
and the wig lines are properly formatted and in order! 

Pass the function an array of key =E<gt> value arguments. An L<IO::File> 
object will be returned. Upon the closing the file handle, the 
F<wigToBigWig> utility will generate the bigWig file.

=over 4

=item bw

The output file name for the bigWig file. Also accepts the keys C<file>, 
C<wig>, and C<out>. This is required.
 
=item chromo

Pass the path to a chromosome sizes text file, described as above. This is 
required. Alternatively, a database object may provided.

=item db

Provide an opened database object from which to generate the chromosome 
sizes information. This will be passed to L<generate_chromosome_file> to 
generate a chromosome file. This is a convenience option and alternative to 
providing an existing chromosome file. B<Note> that you will need to clean 
up the file yourself; this function will not do it for you!

=item chrskip

Provide a regular-expression compatible string for any chromosomes that should 
be skipped when generating a chromosome sizes file from a provided database.

=item bwapppath

Provide the full path to the UCSC F<wigToBigWig> utility. If not provided, the 
default C<PATH> will be searched for the utility. The path may also 
be defined in the configuration file F<.biotoolbox.cfg>. 

=back

Example:

	my $bw_file = 'example.bw';
	my $chromo_file = generate_chromosome_file($db);
	my $bwfh = open_wig_to_bigwig_fh(
		bw      => $bw_file,
		chromo  => $chromo_file,
	);
	foreach (@wig_lines) {
		$bwfh->print("$_\n");
	}
	$bwfh->close;
		# this signals the forked wigToBigWig process to write 
		# the bigWig file, which may take a few seconds to minutes
	unlink $chromo_file;

=item open_bigwig_to_wig_fh

This subroutine will open a forked process from the UCSC F<bigWigToWig> utility 
as a file handle, allowing a bigWig file to be converted to standard text 
wig format and processed as an input stream. Note that the entire file will 
be converted in this manner, not specific locations. This is intended for 
working with the wig file as a whole. 

Note that F<bigWigToWig> will output a wig file in whatever format as the 
bigWig was originally generated with, i.e. fixedStep, varStep, or bedGraph. To 
export explicitly as a bedGraph, which may be useful in certain circumstances, 
the UCSC F<bigWigToBedGraph> utility is also supported.

Note that the UCSC utilities do not always handle errors gracefully 
and will immediately fail upon encountering errors, usually also bringing 
the main Perl process with it. 

Pass the function an array of key =E<gt> value arguments. An L<IO::File> 
object will be returned.

=over 4

=item bw

The output file name for the bigWig file. Also accepts the keys C<file> 
and C<wig>. This is required.
 
=item bwapppath

Provide the full path to the UCSC F<bigWigToWig> or F<bigWigToBedGraph> utility. 
If not provided, the default C<PATH> will be searched for the utility. The path may 
also be defined in the configuration file F<.biotoolbox.cfg>. 

=back

Example:

	my $bw_file = 'example.bw';
	my $bwfh = open_bigwig_to_wig_fh(
		bw    => $bw_file,
	);
	while (my $line = $bwfh->getline) {
		# do something with wig line
	}
	$bwfh->close;

=item bed_to_bigbed_conversion

This subroutine will convert a bed file to a bigBed file. 

After running the utility, the existence of a non-zero byte bigBed file 
is checked. If it does, then the name of the file is returned. If not, 
an error is printed and nothing is returned. 

Pass the function an array of key =E<gt> value arguments, including the 
following:

=over 4

=item bed

The path and name for the bed file. Only standard Bed files with 3-12 
columns are supported. Additional columns, e.g. C<bed6+4> formats, are 
not supported. This value is required.
 
=item chromo

Pass the path to a chromosome sizes text file, described as above. This is 
required. Alternatively, a database object may provided.

=item db

Provide an opened database object from which to generate the chromosome 
sizes information. This will be passed to L<generate_chromosome_file> to 
generate a chromosome file. This is a convenience option and alternative to 
providing an existing chromosome file.

=item bbapppath

Provide the full path to the UCSC F<bedToBigBed> utility. If not provided, the 
default C<PATH> will be searched for the utility. The path may also 
be defined in the configuration file F<.biotoolbox.cfg>. 

=back

Example:

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

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  


