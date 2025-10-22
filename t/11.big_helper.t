#!/usr/bin/env perl

# Test script for Bio::ToolBox::big_helper

use Test2::V0 -no_srand => 1;
plan(40);
use English qw(-no_match_vars);
use File::Spec;
use FindBin '$Bin';

BEGIN {
	## no critic
	$ENV{'BIOTOOLBOX'} = File::Spec->catfile( $Bin, 'Data', 'biotoolbox.cfg' );
	## use critic
}

use IO::File;
use Bio::ToolBox::big_helper qw(
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

# sample files
my $bamfile  = File::Spec->catfile( $Bin, 'Data', 'sample1.bam' );
my $bwfile   = File::Spec->catfile( $Bin, 'Data', 'sample2.bw' );
my $bbfile   = File::Spec->catfile( $Bin, 'Data', 'sample1.bb' );
my $wigfile  = File::Spec->catfile( $Bin, 'Data', 'example.wig' );
my $chrfile  = File::Spec->catfile( $Bin, 'Data', 'chroms.txt' );
my $outfile1 = File::Spec->catfile( $Bin, 'Data', 'example.bw' );
my $bedfile  = File::Spec->catfile( $Bin, 'Data', 'example.bed' );
my $outfile2 = File::Spec->catfile( $Bin, 'Data', 'example.bb' );

## generate chromosome file using variety of database adapters
SKIP: {
	eval { require Bio::DB::HTS };
	skip( "Bio::DB::HTS not installed", 3 ) if $EVAL_ERROR;
	my $file = generate_chromosome_file($bamfile);
	like( $file, qr/chr_sizes_ \w{5} $/x, 'HTS BAM chromosome file' );
	my $fh   = IO::File->new( $file, 'r' );
	my $line = $fh->getline;
	chomp $line;
	my @bits = split /\t/, $line;
	is( $bits[0], 'chrI', 'chromosome correct' );
	is( $bits[1], 230208, 'size chromosome' );
	$fh->close;
	unlink $file;
}
SKIP: {
	eval { require Bio::DB::Sam };
	skip( "Bio::DB::Sam not installed", 3 ) if $EVAL_ERROR;
	my $file = generate_chromosome_file($bamfile);
	like( $file, qr/chr_sizes_ \w{5} $/x, 'Sam BAM chromosome file' );
	my $fh   = IO::File->new( $file, 'r' );
	my $line = $fh->getline;
	chomp $line;
	my @bits = split /\t/, $line;
	is( $bits[0], 'chrI', 'chromosome correct' );
	is( $bits[1], 230208, 'size chromosome' );
	$fh->close;
	unlink $file;
}
SKIP: {
	eval { require Bio::DB::BigWig };
	skip( "Bio::DB::BigWig not installed", 3 ) if $EVAL_ERROR;
	my $file = generate_chromosome_file($bwfile);
	like( $file, qr/chr_sizes_ \w{5} $/x, 'UCSC BigWig chromosome file' );
	my $fh   = IO::File->new( $file, 'r' );
	my $line = $fh->getline;
	chomp $line;
	my @bits = split /\t/, $line;
	is( $bits[0], 'chrI', 'chromosome correct' );
	is( $bits[1], 230208, 'size chromosome' );
	$fh->close;
	unlink $file;
}
SKIP: {
	eval { require Bio::DB::BigBed };
	skip( "Bio::DB::BigBed not installed", 3 ) if $EVAL_ERROR;
	my $file = generate_chromosome_file($bbfile);
	like( $file, qr/chr_sizes_ \w{5} $/x, 'UCSC BigBed chromosome file' );
	my $fh   = IO::File->new( $file, 'r' );
	my $line = $fh->getline;
	chomp $line;
	my @bits = split /\t/, $line;
	is( $bits[0], 'chrI', 'chromosome correct' );
	is( $bits[1], 230208, 'size chromosome' );
	$fh->close;
	unlink $file;
}
SKIP: {
	eval { require Bio::DB::Big };
	skip( "Bio::DB::Big not installed", 6 ) if $EVAL_ERROR;

	# bigWig
	my $file = generate_chromosome_file($bwfile);
	like( $file, qr/chr_sizes_ \w{5} $/x, 'Big BigWig chromosome file' );
	my $fh   = IO::File->new( $file, 'r' );
	my $line = $fh->getline;
	chomp $line;
	my @bits = split /\t/, $line;
	is( $bits[0], 'chrI', 'chromosome correct' );
	is( $bits[1], 230208, 'size chromosome' );
	$fh->close;

	# bigbed
	unlink $file;
	undef $file;
	undef $fh;
	undef $line;
	@bits = ();
	$file = generate_chromosome_file($bbfile);
	like( $file, qr/chr_sizes_ \w{5} $/x, 'Big BigBed chromosome file' );
	$fh   = IO::File->new( $file, 'r' );
	$line = $fh->getline;
	chomp $line;
	@bits = split /\t/, $line;
	is( $bits[0], 'chrI', 'chromosome correct' );
	is( $bits[1], 230208, 'size chromosome' );
	$fh->close;
	unlink $file;
}

## obtain paths to helper applications and open filehandles using them

# bigWigToWig filehandle pipe
# this method could actually return bigWigToBedGraph if available and bigWigToWig is not
SKIP: {
	my $bw2w_app = get_bigwig_to_wig_app();

	skip( 'bigWigToWig not available', 5 ) if not defined $bw2w_app;
	like( $bw2w_app, qr/bigWigTo (?:Wig | BedGraph) $/x, 'bigWigToWig application path' );
	my $fh = open_bigwig_to_wig_fh(
		bw        => $bwfile,
		bwapppath => $bw2w_app
	);
	isa_ok( $fh, ['IO::File'], 'opened bigWigToWig file handle pipe object' );
	is( $fh->opened, 1, 'file handle is opened' );
	my $line   = $fh->getline;
	my $expect = <<END;
variableStep chrom=chrI span=1
END
	is( $line, $expect, 'varStep header line' );
	$line   = $fh->getline;
	$expect = <<END;
55043	-0.551
END
	is( $line, $expect, 'varStep data line' );
	$fh->close;
}

# bigWigToBedGraph filehandle pipe
# this can only return bigWigToBedGraph path
SKIP: {
	my $bw2w_app = get_bigwig_to_bdg_app();
	skip( 'bigWigToBedGraph not available', 4 ) if not defined $bw2w_app;
	like( $bw2w_app, qr/bigWigToBedGraph $/x, 'bigWigToBedGraph application path' );
	my $fh = open_bigwig_to_wig_fh(
		bw        => $bwfile,
		bwapppath => $bw2w_app
	);
	isa_ok( $fh, ['IO::File'], 'opened bigWigToBedGraph file handle pipe object' );
	is( $fh->opened, 1, 'file handle is opened' );
	my $line   = $fh->getline;
	my $expect = <<END;
chrI	55042	55043	-0.551
END
	is( $line, $expect, 'bedGraph data line' );
	$fh->close;
}

## Convert files to big files using external utility

# wigToBigWig tests
SKIP: {
	my $w2bw_app = get_wig_to_bigwig_app();
	skip( 'wigToBigWig not available', 9 ) if not defined $w2bw_app;
	like(
		$w2bw_app, qr/ (?: wigToBigWig | BioDBBigFile ) $/x,
		'wigToBigWig application path'
	);

	# write test wig and chromosome files
	my $fh = IO::File->new( $wigfile, 'w' );
	$fh->print( var_data() );
	$fh->close;
	$fh = IO::File->new( $chrfile, 'w' );
	$fh->print("chrI\t230208\n");
	$fh->close;

	# convert externally
	my $success = wig_to_bigwig_conversion(
		wig       => $wigfile,
		bwapppath => $w2bw_app,
		chromo    => $chrfile,
	);
	is( $success, $outfile1, 'wig to bigWig external conversion success' );
	ok( -e $outfile1,         'output bigWig file1 exists' );
	ok( -s $outfile1 > 25000, 'output bigWig file1 is not empty' );
	if ( -e $outfile1 ) {
		unlink $outfile1;
	}
	unlink $wigfile;

	# convert directly through file handle
SKIP: {
		my $version;
		eval { $version = check_wigToBigWig_version($w2bw_app) };
		skip( 'wigToBigWig does not support stdin', 5 ) if not $version;
		my $bwfh = open_wig_to_bigwig_fh(
			bwapppath => $w2bw_app,
			chromo    => $chrfile,
			bw        => $outfile1
		);
		isa_ok( $bwfh, ['IO::File'], 'opened wigToBigWig pipe object' );
		is( $fh->opened, q(), 'pipe has no valid file descriptor' );
		$bwfh->print( var_data() );
		my $complete = $bwfh->close;
		ok( $complete,            'closed bigWig pipe handle ok' );
		ok( -e $outfile1,         'output bigWig file2 exists' );
		ok( -s $outfile1 > 25000, 'output bigWig file2 is not empty' );
		if ( -e $outfile1 ) {
			unlink $outfile1;
		}
	}
	unlink $chrfile;
}

# bedToBigBed conversion
SKIP: {
	my $b2bb_app = get_bed_to_bigbed_app();
	skip( 'bedToBigBed not available', 4 ) unless defined $b2bb_app;
	like( $b2bb_app, qr/bedToBigBed $/x, 'bedToBigBed application' );

	# write test bed and chromosome files
	my $fh = IO::File->new( $bedfile, 'w' );
	$fh->print( bed_data() );
	$fh->close;
	$fh = IO::File->new( $chrfile, 'w' );
	$fh->print("chrI\t230208\n");
	$fh->close;

	# convert externally
	my $success = bed_to_bigbed_conversion(
		bed       => $bedfile,
		bbapppath => $b2bb_app,
		chromo    => $chrfile,
	);
	is( $success, $outfile2, 'bed to bigBed external conversion success' );
	ok( -e $outfile2,         'output bigBed file3 exists' );
	ok( -s $outfile2 > 18000, 'output bigBed file3 is not empty' );
	if ( -e $outfile2 ) {
		unlink $outfile2;
	}
	unlink $bedfile;
	unlink $chrfile;
}

sub var_data {
	return <<END;
variableStep chrom=chrI span=1
55043	-0.551
55116	-0.55
55138	-0.618
55177	-0.453
55229	-0.241
55268	-0.288
55294	-0.217
55333	-0.166
55385	-0.322
55435	-0.285
55474	-0.318
55502	-0.256
55541	-0.252
55583	-0.239
55612	-0.177
55647	-0.165
55684	-0.177
55709	-0.142
55735	-0.139
55760	-0.109
55782	-0.162
55843	-0.242
55881	-0.265
55922	-0.223
55957	-0.307
55980	-0.266
56009	-0.3
56048	-0.203
56072	-0.109
56098	-0.146
56132	-0.158
56183	-0.135
56208	-0.115
56258	-0.041
56293	0.053
56324	0.103
56371	0.062
56401	-0.035
56444	-0.081
56502	-0.025
56538	0.119
56565	-0.037
56608	0.028
56666	0.149
56730	0.228
56761	0.116
56800	0.2
56822	0.25
56855	0.302
56878	0.294
56912	0.465
56969	0.341
57013	0.432
57056	0.321
57086	0.52
57118	0.514
57165	0.471
57243	0.58
57289	0.482
57377	0.495
57401	0.48
57436	0.537
57469	0.524
57499	0.41
57560	0.392
57583	0.339
57614	0.412
57638	0.532
57691	0.599
57730	0.724
57763	0.75
57807	0.785
57843	0.649
57868	0.585
57924	0.642
57960	0.696
57995	0.601
58025	0.577
58054	0.639
58101	0.686
58149	0.706
58178	0.496
58212	0.52
58244	0.62
58267	0.584
58306	0.51
58338	0.445
58380	0.47
58464	0.487
58504	0.403
58532	0.34
58593	0.09
58626	0.198
58683	0.188
58710	0.037
58739	0.262
58762	0.253
58800	0.105
58822	0.287
58876	0.244
58929	0.307
58974	0.301
59010	0.239
59046	0.195
59099	0.26
59153	0.103
59181	0.192
59205	0.234
59281	0.238
59304	0.199
59336	0.127
59363	0.145
59402	0.094
59451	0.125
59476	0.045
59501	0.062
59548	0.085
59572	0.069
59594	-0.066
59626	0.026
59654	0.083
59684	0.034
59715	0.001
59757	0.021
59779	-0.023
59831	-0.225
59861	-0.03
59945	-0.031
59975	-0.081
60003	0.011
60027	-0.015
60091	0.05
60169	0.137
60212	0.16
60244	0.122
60287	0.224
60315	0.2
60342	0.103
60366	0.26
60393	0.237
60444	0.318
60468	0.289
60530	0.317
60552	0.35
60587	0.437
60613	0.39
60651	0.374
60700	0.441
60731	0.461
60761	0.443
60834	0.405
60894	0.341
60940	0.418
60994	0.407
END
}

sub bed_data {
	return <<END;
chrI	55927	56000	align1	150	-
chrI	55928	56001	align2	150	+
chrI	55928	56001	align3	150	+
chrI	55935	56008	align4	150	-
chrI	55936	56009	align5	150	+
chrI	55939	56012	align6	150	+
chrI	55939	56012	align7	150	+
chrI	55941	56014	align8	150	+
chrI	55944	56017	align9	150	+
chrI	55948	56021	align10	150	+
chrI	55956	56029	align11	150	+
chrI	55967	56040	align12	150	+
chrI	55968	56041	align13	150	+
chrI	55973	56046	align14	150	+
chrI	55975	56048	align15	150	-
chrI	55976	56049	align16	150	+
chrI	55982	56055	align17	150	+
chrI	55982	56038	align18	150	-
chrI	55983	56028	align19	60	-
chrI	55985	56058	align20	150	+
chrI	55988	56061	align21	150	+
chrI	55989	56062	align22	150	+
chrI	55993	56066	align23	150	-
chrI	55999	56048	align24	103	-
chrI	56001	56074	align25	150	+
chrI	56003	56050	align26	77	-
chrI	56004	56077	align27	150	-
chrI	56004	56077	align28	150	+
chrI	56008	56051	align29	36	-
chrI	56010	56083	align30	150	+
chrI	56010	56056	align31	71	-
chrI	56012	56056	align32	36	+
chrI	56014	56087	align33	150	-
chrI	56014	56087	align34	150	-
chrI	56015	56088	align35	150	+
chrI	56015	56088	align36	150	+
chrI	56026	56099	align37	150	-
chrI	56027	56100	align38	150	+
chrI	56027	56100	align39	150	-
chrI	56031	56104	align40	150	-
chrI	56034	56107	align41	150	-
chrI	56035	56108	align42	150	+
chrI	56038	56111	align43	150	+
chrI	56038	56111	align44	150	-
chrI	56042	56115	align45	150	-
chrI	56042	56115	align46	150	-
chrI	56043	56116	align47	150	+
chrI	56047	56120	align48	150	-
chrI	56048	56121	align49	150	-
chrI	56048	56121	align50	150	-
END
}

