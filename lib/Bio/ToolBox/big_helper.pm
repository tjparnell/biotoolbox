package Bio::ToolBox::big_helper;

use warnings;
use strict;
use English qw(-no_match_vars);
use Carp    qw(carp cluck croak);
use File::Temp;
use File::Which;
use IO::File;
use Bio::ToolBox::db_helper         qw(get_chromosome_list);
use Bio::ToolBox::db_helper::config qw($BTB_CONFIG);
require Exporter;

our $VERSION = '2.02';

### Export
our @ISA       = qw(Exporter);
our @EXPORT_OK = qw(
	get_bed_to_bigbed_app
	get_bigwig_to_bdg_app
	get_bigwig_to_wig_app
	get_wig_to_bigwig_app
	check_wigToBigWig_version
	open_wig_to_bigwig_fh
	open_bigwig_to_wig_fh
	bed_to_bigbed_conversion
	generate_chromosome_file
	wig_to_bigwig_conversion
);

### Find Converter external applications
sub get_wig_to_bigwig_app {

	# we prefer the wigToBigWig convertor

	# first check for an entry in the configuration file then environment path
	my $app =
		   $BTB_CONFIG->param('applications.wigToBigWig')
		|| which('wigToBigWig')
		|| undef;

	unless ($app) {

		# last attempt is to use the older Bio::DB::BigFile module
		# this may not be installed and will be slower since it runs in same process
		eval {
			require Bio::DB::BigFile;
			if ( Bio::DB::BigFile->can('createBigWig') ) {
				$app = 'BioDBBigFile';
			}
		};
	}
	return $app;
}

sub get_bigwig_to_wig_app {

	# first check for an entry in the configuration file then environment path
	return
		   $BTB_CONFIG->param('applications.bigWigToWig')
		|| which('bigWigToWig')
		|| which('bigWigToBedGraph')
		|| undef;
}

sub get_bigwig_to_bdg_app {

	# first check for an entry in the configuration file then environment path
	return
		   $BTB_CONFIG->param('applications.bigWigToBedGraph')
		|| which('bigWigToBedGraph')
		|| undef;
}

sub get_bed_to_bigbed_app {

	# first check for an entry in the configuration file then environment path
	return
		   $BTB_CONFIG->param('applications.bedToBigBed')
		|| which('bedToBigBed')
		|| undef;
}

### Check the version of wigToBigWig
# older versions support reading from stdin while newer versions do not
sub check_wigToBigWig_version {
	my $app = shift;
	return unless $app;
	return 0 if $app eq 'BioDBBigFile';
	if ( not -e $app ) {
		carp " '$app' does not exist!";
		exit 1;
	}
	if ( not -x _ ) {
		carp " '$app' is not executable!";
		exit 1;
	}
	if ( $app =~ /bedGraphToBigWig$/x ) {

		# this app is not always reliable
		# generally does not accept stdin pipes, but even if it's old enough to do so
		# it may still complain about line seeking and/or chromosome sort order
		# best just not support it
		return 0;
	}
	my $result = qx($app 2>&1);
	if ( $result =~ /wigToBigWig \s v \s ([\d\.]+)/x ) {
		return 1 if $1 eq '2.8';    # older acceptable version
		return 0 if $1 eq '2.9';    # newer unsupported version
		return 1 if $1 eq '4';      # very old, this is actually the bbi version I think
		return 0;                   # do not trust anything else unfortunately
	}
	return 0;
}

### Wig to BigWig file conversion
sub wig_to_bigwig_conversion {

	# Collect passed arguments
	my %args = @_;
	unless (%args) {
		carp 'no arguments passed!';
		return;
	}

	# wigfile
	$args{wig} ||= undef;
	unless ( $args{wig} ) {
		carp 'no wig file passed!';
		return;
	}

	# Identify bigwig conversion utility
	my $bwapp = $args{bwapppath} || undef;
	unless ($bwapp) {
		$bwapp = get_wig_to_bigwig_app();
	}
	unless ($bwapp) {
		carp
q(Utility 'wigToBigWig' not specified and can not be found! Conversion failed!);
		return;
	}

	# Generate list of chromosome sizes if necessary
	my $chromfile = $args{chromo} ||= undef;
	unless ($chromfile) {

		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{db} ||= undef;
		unless ( $args{db} ) {
			carp
'No requisite database or chromosome info file provided! Conversion failed';
			return;
		}
		$chromfile = generate_chromosome_file( $args{db} );
		unless ($chromfile) {
			carp 'Cannot generate chromosome info file! Conversion failed';
			return;
		}
	}

	# Generate the bw file name
	# we can substitute one of three possible names for bw
	my $bw_file = $args{wig};
	$bw_file =~ s/\.(?: bed | bdg | bedgraph | wig) $/.bw/x;

	# Generate the bigwig file
	printf " Converting %s to %s...\n", $args{wig}, $bw_file;
	if ( $bwapp =~ /wigToBigWig$/x ) {

		# include the -clip option in case there are any positions
		# out of bounds of the chromosome it will just warn instead of fail
		if ( system( $bwapp, '-clip', $args{wig}, $chromfile, $bw_file ) ) {
			my $c = join q( ), $bwapp, '-clip', $args{wig}, $chromfile, $bw_file;
			print STDERR "conversion command '$c' failed - check error";
			return;
		}

	}
	elsif ( $bwapp =~ /bedGraphToBigWig$/x ) {

		# this doesn't have the -clip option
		if ( system( $bwapp, $args{wig}, $chromfile, $bw_file ) ) {
			my $c = join q( ), $bwapp, $args{wig}, $chromfile, $bw_file;
			print STDERR "conversion command '$c' failed - check error";
			return;
		}
	}
	elsif ( $bwapp eq 'BioDBBigFile' ) {

		# this cannot be caught if it fails, it will terminate the perl process
		print " Converting with Bio::DB::BigFile - this will terminate if errors\n";
		Bio::DB::BigFile->createBigWig( $args{wig}, $chromfile, $bw_file );
	}

	# check the result
	if ( -e $bw_file and -s $bw_file ) {

		# conversion assumed successful
		if ( not $args{chromo} ) {

			# we no longer need our temp chromosome file
			unlink $chromfile;
		}
		return $bw_file;
	}
	else {
		if ( -e $bw_file ) {
			unlink $bw_file;    # remove any partial file
		}
		if ( not $args{chromo} ) {

			# leave the temp chromosome file as a courtesy
			carp
" Conversion failed. You should try manually and watch for errors\n Leaving temporary chromosome file '$chromfile'.";
		}
		else {
			carp 'Conversion failed. You should try manually and watch for errors.';
		}
		return;
	}
}

### Open a file handle to wigToBigWig
sub open_wig_to_bigwig_fh {

	# Collect passed arguments
	my %args = @_;
	unless (%args) {
		carp 'no arguments passed!';
		return;
	}

	# bigWig output file
	$args{bw} ||= $args{wig} || $args{out} || $args{file} || undef;
	unless ( $args{bw} ) {
		carp 'no output bw file name passed!';
		return;
	}
	unless ( $args{bw} =~ /\.bw$/i ) {
		$args{bw} .= '.bw';
	}

	# Identify bigwig conversion utility
	my $bwapp = $args{bwapppath} || get_wig_to_bigwig_app() || undef;
	if ($bwapp) {

		# check the version of the utility
		if ( not check_wigToBigWig_version($bwapp) ) {
			carp
" $bwapp does not support stdin file handles. Either downgrade or use wig_to_bigwig_conversion()\n";
			return;
		}
	}
	else {
		print STDERR q( Utility 'wigToBigWig' not specified and can not be found!);
		return;
	}

	# Generate list of chromosome sizes if necessary
	$args{chromo} ||= undef;
	unless ( $args{chromo} ) {

		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{db} ||= undef;
		unless ( $args{db} ) {
			carp
'No requisite database or chromosome info file provided! Conversion failed';
			return;
		}
		$args{chrskip} ||= undef;
		$args{chromo} = generate_chromosome_file( $args{db}, $args{chrskip} );
		unless ( $args{chromo} ) {
			carp ' Cannot generate chromosome info file! Conversion failed';
			return;
		}
	}

	# open the filehandle
	my $command = sprintf "%s stdin %s %s", $bwapp, $args{chromo}, $args{bw};
	my $bwfh    = IO::File->new("| $command")
		or croak sprintf( "cannot open %s! $OS_ERROR", $args{bwapppath} );

	# wigToBigWig will always die anyway if something is wrong
	# cannot trap it with an eval, since it doesn't just error out
	# but actually exits, dragging the whole Perl process with it
	# printf "we have a filehandle %s\n", ref $bwfh;
	croak "unable to execute command '$command'" unless ref $bwfh;

	# we will still get an IO::File handle back even with a failed convertor - sigh
	return $bwfh;
}

### Open a file handle from bigWigToWig
sub open_bigwig_to_wig_fh {

	# Collect passed arguments
	my %args = @_;
	unless (%args) {
		cluck 'no arguments passed!';
		return;
	}

	# bigWig output file
	$args{bw} ||= $args{wig} || $args{file} || undef;
	unless ( $args{bw} ) {
		carp 'no input bw file name passed!';
		return;
	}

	# Identify bigwig conversion utility
	my $bwapp = $args{bwapppath} || get_bigwig_to_wig_app() || undef;
	unless ($bwapp) {

		# check for an entry in the configuration file
		$args{bwapppath} = $BTB_CONFIG->param('applications.bigWigToWig') || undef;
	}
	unless ( $args{bwapppath} ) {

		# try checking the system path
		$args{bwapppath} = which('bigWigToWig') || which('bigWigToBedGraph');
	}
	unless ( $args{bwapppath} =~ /bigWigTo(?: Wig | BedGraph)$/x ) {
		carp q(Utility 'bigWigToWig' not specified and can not be found!);
		return;
	}

	# open the filehandle
	my $command = sprintf "%s %s stdout", $args{bwapppath}, $args{bw};
	my $bwfh    = IO::File->new("$command |")
		or croak sprintf( "cannot open %s! $OS_ERROR", $args{bwapppath} );

	# bigWigToWig will always die anyway if something is wrong
	# cannot trap it with an eval, since it doesn't just error out
	# but actually exits, dragging the whole Perl process with it
	# printf "we have a filehandle %s\n", ref $bwfh;
	croak "unable to execute command '$command'" unless ref $bwfh;

	# we will still get an IO::File handle back even with a failed convertor - sigh
	return $bwfh;
}

### Bed to BigBed file conversion
sub bed_to_bigbed_conversion {

	# Collect passed arguments
	my %args = @_;
	unless (%args) {
		carp 'no arguments passed!';
		return;
	}

	# bedfile
	$args{bed} ||= undef;
	unless ( $args{bed} ) {
		carp 'no bed file passed!';
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
		carp q(Utility 'bedToBigBed' not specified and can not be found!);
		return;
	}

	# Generate list of chromosome sizes if necessary
	$args{chromo} ||= undef;
	unless ( $args{chromo} ) {

		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		$args{db} ||= undef;
		unless ( $args{db} ) {
			carp 'No requisite database or chromosome info file provided!';
			return;
		}
		$args{chromo} = generate_chromosome_file( $args{db} );
		unless ( $args{chromo} ) {
			carp 'Cannot generate chromosome info file!';
			return;
		}
	}

	# Generate the bb file name
	my $bb_file = $args{bed};
	$bb_file =~ s/\.bed$/.bb/;

	# Generate the bigBed file using Jim Kent's utility
	printf " converting %s to BigBed....\n", $args{bed};
	system( $args{bbapppath}, $args{bed}, $args{chromo}, $bb_file );

	# Check the result
	if ( -e $bb_file and -s $bb_file ) {

		# conversion successful
		if ( $args{chromo} =~ /^chr_sizes_\w{5}/x ) {

			# we no longer need our temp chromosome file
			unlink $args{chromo};
		}
		return $bb_file;
	}
	else {
		if ( -e $bb_file ) {
			unlink $bb_file;    # remove any partial file
		}
		if ( $args{chromo} =~ /^chr_sizes_\w{5}/x ) {

			# leave the temp chromosome file as a courtesy
			carp
" Conversion failed. You should try manually and watch for errors\n Leaving temporary chromosome file '$args{chromo}'.";
		}
		else {
			carp 'Conversion failed. You should try manually and watch for errors.';
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
		carp 'no chromosome sequences identified in database!';
		return;
	}

	# prepare temp file
	my $chr_fh = new File::Temp(
		'UNLINK'   => 0,
		'TEMPLATE' => 'chr_sizes_XXXXX',
	) or croak "unable to open chromosome temp file! $OS_ERROR";
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

=head3 Note regarding compatibility

Some of these methods open a pipe to these external utilities, either
for reading (L</open_wig_to_bigwig_fh>) or writing (L<open_bigwig_to_wig_fh>). 
Newer (current) versions of these utilities, for
example C<wigToBigWig> or C<bedGraphToBigWig>, no longer support reading from 
standard input, resulting in failure with L<open_bigwig_to_wig_fh>. This change
occurred with version release 439 of the UCSC UserApps, released circa 2022-11-14.

Being able to write to the utility directly has its advantages, namely avoiding 
writing an intermediate text file and speed. 

To check whether your version of C<wigToBigWig> supports reading from C<stdin>,
use the L</check_wigToBigWig_version> method. If it does not, simply write to a
text wig file and use the L</wig_to_bigwig_conversion> function. 

Otherwise, you can always try compiling an older version.

The C<bedToBigBed> and C<bedGraphtoBigWig> utilities do not support reading from
C<stdin>, primarily because they perform read seek on the input file.


=head1 USAGE

Load the module at the beginning of your program and include the name or 
names of the subroutines to export. None are automatically exported.

	use Bio::ToolBox::big_helper qw(wig_to_bigwig_conversion);

There are are nine available exported subroutines.

=head2 Find and check applications

These functions look for UCSC utility applications in your environment C<PATH>
or defined in your F<.biotoolbox.cfg> config file (if present) and returns
the path or undefined.

=over 4

=item get_bed_to_bigbed_app

Looks for the C<bedToBigBed> application for converting text bed files to
binary bigBed.

=item get_bigwig_to_bdg_app

Looks for the C<bigWigToBedGraph> application for converting binary bigWig
to text bedGraph.

=item get_bigwig_to_wig_app

Looks for the C<bigWigToWig> application for converting binary bigWig to
the original text wig format (fixedStep, varStep, or bedGraph).

=item get_wig_to_bigwig_app

Looks for the C<wigToBigWig> application for converting text wig (fixedStep,
varStep, or bedGraph) to binary bigWig.

=item check_wigToBigWig_version

This function checks the version of the C<wigToBigWig> application and
whether it supports a direct pipe to C<stdin>. Pass the path of the utility.
Returns a boolean 1 or 0.

=back


=head2 Open file handles

These functions will open a either a read or write filehandle to a C<bigWig>
file through an appropriate external utility. The respective input or output
should be properly formatted formatted C<wig> text, including C<fixedStep>,
C<variableStep>, or C<bedGraph>.

=over 4

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

B<NOTE>: Recent versions of F<wigToBigWig> no longer support C<stdin> as
a file handle and cannot be used here.

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

=back

=head2 Convert files

These functions easily handle the conversion of a text file to a binary 
bigWig or bigBed file.

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

=item bed_to_bigbed_conversion

This subroutine will convert a bed file to a bigBed file. 

After running the utility, the existence of a non-zero byte bigBed file 
is checked. If it does, then the name of the file is returned. If not, 
an error is printed and nothing is returned. 

Note that this utility requires chromosomes to be sorted in ASCIbetical
order, and not necessarily numeric order.

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


