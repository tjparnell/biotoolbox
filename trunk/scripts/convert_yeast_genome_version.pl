#!/usr/bin/env perl

# documentation at end of file 

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::data_helper qw(find_column_index);
use Bio::ToolBox::file_helper qw(
	open_tim_data_file
	write_tim_data_file
	open_to_write_fh
);
use Bio::ToolBox::big_helper qw(wig_to_bigwig_conversion);
use Bio::ToolBox::db_helper::config qw($BTB_CONFIG add_program);

my $VERSION = '1.14';

print "\n This script will convert yeast genomic coordinates\n";


### Quick help
unless (@ARGV) { # when no command line options are present
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Get command line options and initialize values

# Initialize values
my (
	$version,
	$convert_to_roman,
	$database,
	$help,
	$print_version,
); # command line variables
my @infiles; 


# Command line options
GetOptions( 
	'in=s'        => \@infiles, # input file(s)
	'convert=i'   => \$version, # the number of the conversion subroutine
	'roman!'      => \$convert_to_roman, # convert chromosome names
	'db=s'        => \$database, # for bigwig conversions
	'help'        => \$help, # print the help
	'version'     => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";


# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script convert_yeast_genome_version.pl, version $VERSION\n\n";
	exit;
}



### Check for general required values

# files
unless (@infiles) {
	# file list was provided on the command line
	@infiles = @ARGV or
		die " No input files! use --help for more information\n";
}
if (scalar @infiles == 1) {
	# only one file provided, but may be comma delimited list
	@infiles = split /,/, shift @infiles;
}

# Prepare for converting to roman numerals
my %old2new = (
	'chr1'	=>	'chrI',
	'chr2'	=>	'chrII',
	'chr3'	=>	'chrIII',
	'chr4'	=>	'chrIV',
	'chr5'	=>	'chrV',
	'chr6'	=>	'chrVI',
	'chr7'	=>	'chrVII',
	'chr8'	=>	'chrVIII',
	'chr9'	=>	'chrIX',
	'chr10'	=>	'chrX',
	'chr11'	=>	'chrXI',
	'chr12'	=>	'chrXII',
	'chr13'	=>	'chrXIII',
	'chr14'	=>	'chrXIV',
	'chr15'	=>	'chrXV',
	'chr16'	=>	'chrXVI',
	'chrMito'	=>	'chrMito', # no change for these
	'2-micron' =>	'2-micron',
);

# find bigWigToBedGraph utility
	# check for an entry in the configuration file
my $bdg_app_path = $BTB_CONFIG->param('applications.bigWigToBedGraph') || undef;
unless ($bdg_app_path) {
	# next check the system path
	eval {
		require File::Which;
		File::Which->import;
		$bdg_app_path = which('bigWigToBedGraph');
	};
	add_program($bdg_app_path) if $bdg_app_path; # remember for next time
}

# find bedGraphToBigWig utility
	# check for an entry in the configuration file
my $bw_app_path = $BTB_CONFIG->param('applications.bedGraphToBigWig') || undef;
unless ($bw_app_path) {
	# next check the system path
	eval {
		require File::Which;
		File::Which->import;
		$bw_app_path = which('bigWigToBedGraph');
	};
	add_program($bw_app_path) if $bw_app_path; # remember for next time
}
		



### Determine conversion method

# Ask for which conversion
unless ($version) {
	&print_version_choices;
	print " enter the number of the conversion you want   ";
	$version = <STDIN>;
	chomp $version;
}

# determine method
my $goconvert; # this will be a reference to the appropriate version subroutine
my $old_version_number; # the date of the originating version
my $new_version_number; # the date of the new version
determine_method();





### Convert files
foreach my $f (@infiles) {
	process_file($f);
}

print "Done!\n";

############# Subroutines ######################################################


sub determine_method {
	if ($version == 1) {
		$goconvert = \&version_20050806_20070113;
		$old_version_number = '20050806';
		$new_version_number = '20070113';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	elsif ($version == 2) {
		$goconvert = \&version_20051106_20070113;
		$old_version_number = '20051106';
		$new_version_number = '20070113';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	elsif ($version == 3) {
		$goconvert = \&version_20060204_20070113;
		$old_version_number = '20060204';
		$new_version_number = '20070113';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	elsif ($version == 4) {
		$goconvert = \&version_20031001_20070113;
		$old_version_number = '20031001';
		$new_version_number = '20070113';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	elsif ($version == 5) {
		$goconvert = \&version_200601xx_20070113;
		$old_version_number = '200601xx';
		$new_version_number = '20070113';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	elsif ($version == 6) {
		$goconvert = \&version_20060916_20070113;
		$old_version_number = '20060916';
		$new_version_number = '20070113';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	elsif ($version == 7) {
		$goconvert = \&version_20070113_20100109;
		$old_version_number = '20070113';
		$new_version_number = '20100109';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	elsif ($version == 8) {
		$goconvert = \&version_20050806_20100109;
		$old_version_number = '20050806';
		$new_version_number = '20100109';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	elsif ($version == 9) {
		$goconvert = \&version_20070113_20100109;
		# the 20070113 and 20070901 versions are identical so reuse the sub
		$old_version_number = '20070901';
		$new_version_number = '20100109';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	elsif ($version == 10) {
		$goconvert = \&version_20100109_20110203;
		$old_version_number = '20100109';
		$new_version_number = '20110203';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	elsif ($version == 11) {
		$goconvert = \&version_20080606_20110203;
		$old_version_number = '20080606';
		$new_version_number = '20110203';
		print " converting from version $old_version_number to " . 
			"$new_version_number...\n";
	} 
	else {
		die " unknown version conversion!\n";
	}
}


sub print_version_choices {
print " Available conversions:
   1) 20050806 (SGD R43) to 20070113 (SGD R55)
   2) 20051106 (SGD R46) to 20070113 (SGD R55)
   3) 20060204 (SGD R52) to 20070113 (SGD R55)
   4) 20031001 (SGD R27, UCSC SacCer1) to 20070113 (SGD R55)
   5) 200601xx (SGD R52?) to 20070113 (SGD R55)
   6) 20060916 (Pugh's GeneTrack, SGD R53) to 20070113 (SGD R55)
   7) 20070113 (SGD R55) to 20100109 (SGD R63)
   8) 20050806 (perocchi..steinmetz 2007, SGD R43?) to 20100109 (SGD R63)
   9) 20070901 (xu..steinmetz 2009, SGD R56) to 20100109 (SGD R63)
   10) 20100109 (SGD R63) to 20110203 (SGD R64, UCSC SacCer3)
   11) 20080606 (SGD R61, UCSC SacCer2) to 20110203 (SGD R64, UCSC SacCer3)
";
}


sub process_file {
	my $filename = shift;
	print " converting file '$filename'....\n";
	
	# pre-processing bigWig files
	my $bigwig = 0;
	if ($filename =~ /\.bw$/) {
		print "  first converting bigWig to bedGraph....\n";
		
		# check for utility
		unless ($bdg_app_path) {
			die " no bedWigToBedGraph utility found! skipping\n";
		}
		
		# convert to bedgraph
		my $bedgraph = $filename;
		$bedgraph =~ s/\.bw$/.bdg/;
		system($bdg_app_path, $filename, $bedgraph) == 0 or 
			die " can't execute bigWigToBedGraph!\n";
		if (-s $bedgraph) {
			# success
			$filename = $bedgraph;
			$bigwig = 1;
		}
		else {
			die " unable to convert bigWig file!\n";
		}
	}
	
	
	# Open the input file
	my ($fh, $metadata_ref) = open_tim_data_file($filename);
	
	# Identify the columns
	my $chromo_index;
	my $pos_index; # first position, or start
	my $pos2_index; # a second position, or end
	if ($metadata_ref->{extension} =~ /sgr/i) { 
		# file is a *.sgr
		$chromo_index = 0;
		$pos_index = 1;
	} 
	elsif ($metadata_ref->{gff}) { 
		# file is a gff file
		$chromo_index = 0;
		$pos_index = 3;
		$pos2_index = 4;
	} 
	elsif ($metadata_ref->{bed}) { 
		# file is a bed file
		$chromo_index = 0;
		$pos_index = 1;
		$pos2_index = 2;
	} 
	else { 
		# a non-standard file
		# identify the indices to the genomic coordinates
		$chromo_index = find_column_index($metadata_ref, '^chr|seq|ref');
		$pos_index    = find_column_index($metadata_ref, '^start|pos');
		$pos2_index   = find_column_index($metadata_ref, '^stop|end$');
		unless ( 
			# check that we have these basic genomic coordinates
			defined $chromo_index and
			defined $pos_index
		) {
			warn " unable to identify chromosome, start, and/or stop datasets! skipping!";
			return;
		}
		print " Using column $chromo_index for chromosome, $pos_index for start";
		if (defined $pos2_index) {
			print ", $pos2_index for stop\n";
		}
		else {
			print "\n";
		}
	}
	
	
	
	# Open the output file
	my $gz;
	if ($metadata_ref->{extension} =~ /gz/) {
		$gz = 1;
	}
	else {
		$gz = 0;
	}
	push @{ $metadata_ref->{other} }, 
		# record a miscelleneous header line to explain the conversion
		"# Converted genome coordinates from $old_version_number to $new_version_number";
	$metadata_ref->{data_table} = []; # create an empty data table
	push @{ $metadata_ref->{data_table} }, $metadata_ref->{'column_names'};
	my $new_filename = $metadata_ref->{path} . $metadata_ref->{basename} . '.' . 
		$new_version_number . $metadata_ref->{extension};
	unless ( write_tim_data_file(
			'data'       => $metadata_ref,
			'filename'   => $new_filename,
			'gz'         => $gz,
		)
	) {
		die " unable to open file for writing!\n";
	}
	my $out_fh = open_to_write_fh($new_filename, $gz, 1) or 
		die " unable to re-open file for writing!\n";
	
	
	
	# Do the conversion
	my $totalcount = 0; # a running tally
	my $changecount = 0; # a running tally
	while (my $line = $fh->getline) {
		my $changed = 0;
		
		# process each line
		chomp $line;
		my @data = split /\t/, $line;
		my $chromo = $data[$chromo_index];
		my $position = $data[$pos_index];
		
		# convert the first position
		my $newposition = &{$goconvert}(convert_chr($chromo), $position);
		if ($newposition != $position) {
			$data[$pos_index] = $newposition;
			$changed = 1;
		}
		
		# convert the second position if necessary
		if (defined $pos2_index) {
			my $position2 = $data[$pos2_index];
			if ($position == $position2) {
				# they're the same position, often the case in my gff files
				# then match the new position if necessary
				if ($newposition != $position) {
					$data[$pos2_index] = $newposition;
				}
			} 
			else {
				# using different start and end positions
				my $newposition2 = 
					&{$goconvert}(convert_chr($chromo), $position2); 
				if ($newposition2 != $position2) {
					$data[$pos2_index] = $newposition2;
					$changed = 1;
				}
			}
		}
		
		# convert chromosome name
		if ($convert_to_roman) {
			$data[$chromo_index] = convert_chr($chromo);
		}
		
		# put back
		$changecount++ if $changed;
		$totalcount++;
		print {$out_fh} join("\t", @data) . "\n";	
	}
	close INFILE;
	print " converted $changecount out of $totalcount features to new coordinates\n";
	
	
	# Finish
	$fh->close;
	$out_fh->close;
	
	# post process
	if ($bigwig) {
		print "  post-converting bedGraph back to bigWig....\n";
		
		# perform the conversion
		my $bw_file = wig_to_bigwig_conversion(
				'wig'       => $new_filename,
				'db'        => $database,
				'bwapppath' => $bw_app_path,
		);
	
		
		# confirm
		if ($bw_file) {
			print " bigwig file '$bw_file' generated\n";
			unlink $new_filename; # remove the converted bedgraph file
			unlink $filename; # actually the original bedgraph file
			$new_filename = $bw_file;
		}
		else {
			die " bigwig file not generated! see standard error, no files removed\n";
		}
	}
	
	print " wrote converted file '$new_filename'\n";
}



sub convert_chr {
	my $chromo = shift;
	
	# return if looks ok
	return $chromo if $chromo =~ /^chr[IVX]+$/;
	
	# add chr prefix if missing (numerals only)
	if ($chromo =~ m/^[\dIXV]+$/) {
		$chromo = 'chr' . $chromo;
	}
	
	# normalize prefix case
	$chromo =~ s/^chr/chr/i;
		
	# remove preceding 0 if present
	$chromo =~ s/0(\d)$/$1/;
	
	if (exists $old2new{$chromo}) {
		return $old2new{$chromo};
	}
	else {
		warn " do not recognize '$chromo' to convert to Roman numeral!\n";
		return $chromo;
	}
}




# the conversion subroutines
	# the adjustment numbers below were obtained by running the EMBOSS program 
	# 'diffseq' on the two sequences. I used -wordsize of 8 or 10 and -globaldifferences Y
	# I wrote a shell script find_diffs.sh to expediate the search
	# The diffseq program will report the location and sequence changes between the two 
	# sequence versions. I am using the highest position of the 2003 sequence as 
	# my trigger below, and using the difference in length between the two mismatched 
	# regions. If 2007 is longer, then +; if 2003 is longer, then -.
	
	# new version subroutines may be added by duplicating the subroutine,
	# renaming it, and making the necessary changes. Don't forget to add it to the
	# choice menu below
	
	# here is the list of methods I've been using lately after running diffseq
	# grep ^(chr\w+) ([\d\-]+) Length: (\d+) ?\(?(\w*)\)?\rSequence: \w*\rSequence: \w*\rchr\w+ ([\d\-]+) Length: (\d+)[ \w\(\)]*$
	# replace with \1\t\2\t\3\t\5\t\6\t\4
	# toss empty and comment lines
	# extract indels with perl -w -n -e 'chomp; next if /^#|C/; my @a = split /\t/; next if $a[2] == $a[4]; $pos = $1 if $a[1] =~ /^(\d+)/; if ($a[5]) {print join("\t", $a[0], $pos, $a[4]) . "\n";}  elsif ($a[2] == 2 and $a[4] == 1) { print join("\t", $a[0], $pos, 1) . "\n";} else {print join("\t", $a[0], $pos, $a[4] - $a[2]) . "\n";}' file.table.txt | pbcopy
	# new document from clipboard
	# manually add new chromosome name by itself at beginning of each chromo block
	# grep "^chr(\w+)\t(\d+)\t(\d+)$"
	# replace with "		$new += \3 if ($position >= \2);"
	# grep "^chr(\w+)\t(\d+)\t\-(\d+)$"
	# replace with "		$new -= \3 if ($position >= \2);"
	# grep "^chr(\w+)$"
	# replace with "\t}\r\t\r	elsif ($chromo eq '\1') {"
	

sub version_20050806_20070113 {
	my ($chromo, $position) = @_;
	my $new = $position;
	if ($chromo eq 'chrI') { # Chromosome 1
		# no changes
	} elsif ($chromo eq 'chrII') { # Chromosome 2
		# no changes
	} elsif ($chromo eq 'chrIII') { # Chromosome 3
		if ($position >= 76147) {$new += 1}
	} elsif ($chromo eq 'chrIV') { # Chromosome 4
		if ($position >= 478026) {$new += 1}
		if ($position >= 613205) {$new += 1}
	} elsif ($chromo eq 'chrV') { # Chromosome 5
		# no changes
	} elsif ($chromo eq 'chrVI') { # Chromosome 6
		# no changes
	} elsif ($chromo eq 'chrVII') { # Chromosome 7
		if ($position >= 946899) {$new += 1}
		if ($position >= 962548) {$new -= 1}
	} elsif ($chromo eq 'chrVIII') { # Chromosome 8
		if ($position >= 54132) {$new += 1}
	} elsif ($chromo eq 'chrIX') { # Chromosome 9
		# no changes
	} elsif ($chromo eq 'chrX') { # Chromosome 10
		if ($position >= 101726) {$new += 1}
		if ($position >= 120804) {$new += 6}
		if ($position >= 120890) {$new += 72}
	} elsif ($chromo eq 'chrXI') { # Chromosome 11
		if ($position >= 185987) {$new -= 2}
		if ($position >= 255409) {$new += 1}
		if ($position >= 363801) {$new += 1}
		if ($position >= 552426) {$new += 1}
		if ($position >= 552516) {$new -= 1}
	} elsif ($chromo eq 'chrXII') { # Chromosome 12
		if ($position >= 922648) {$new += 1}
	} elsif ($chromo eq 'chrXIII') { # Chromosome 13
		# no changes
	} elsif ($chromo eq 'chrXIV') { # Chromosome 14
		if ($position >= 60402) {$new += 1}
		if ($position >= 68346) {$new += 1} # this one is not present in 20051106
	} elsif ($chromo eq 'chrXV') { # Chromosome 15
		if ($position >= 42534) {$new += 1}
		if ($position >= 803745) {$new += 1}
	} elsif ($chromo eq 'chrXVI') { # Chromosome 16
		# no changes
	}
	# no changes for Chromosome MT
	return $new;
}



sub version_20051106_20070113 {
	my ($chromo, $position) = @_;
	my $new = $position;
	if ($chromo eq 'chrI') { # Chromosome 1
		# no changes
	} elsif ($chromo eq 'chrII') { # Chromosome 2
		# no changes
	} elsif ($chromo eq 'chrIII') { # Chromosome 3
		if ($position >= 76147) {$new += 1}
	} elsif ($chromo eq 'chrIV') { # Chromosome 4
		if ($position >= 478026) {$new += 1}
		if ($position >= 613205) {$new += 1}
	} elsif ($chromo eq 'chrV') { # Chromosome 5
		# no changes
	} elsif ($chromo eq 'chrVI') { # Chromosome 6
		# no changes
	} elsif ($chromo eq 'chrVII') { # Chromosome 7
		if ($position >= 946899) {$new += 1}
		if ($position >= 962548) {$new -= 1}
	} elsif ($chromo eq 'chrVIII') { # Chromosome 8
		if ($position >= 54132) {$new += 1}
	} elsif ($chromo eq 'chrIX') { # Chromosome 9
		# no changes
	} elsif ($chromo eq 'chrX') { # Chromosome 10
		if ($position >= 101726) {$new += 1}
		if ($position >= 120804) {$new += 6}
		if ($position >= 120890) {$new += 72}
	} elsif ($chromo eq 'chrXI') { # Chromosome 11
		if ($position >= 185987) {$new -= 2}
		if ($position >= 255409) {$new += 1}
		if ($position >= 363801) {$new += 1}
		if ($position >= 552426) {$new += 1}
		if ($position >= 552516) {$new -= 1}
	} elsif ($chromo eq 'chrXII') { # Chromosome 12
		if ($position >= 922648) {$new += 1}
	} elsif ($chromo eq 'chrXIII') { # Chromosome 13
		# no changes
	} elsif ($chromo eq 'chrXIV') { # Chromosome 14
		if ($position >= 60402) {$new += 1}
	} elsif ($chromo eq 'chrXV') { # Chromosome 15
		if ($position >= 42534) {$new += 1}
		if ($position >= 803745) {$new += 1}
	} elsif ($chromo eq 'chrXVI') { # Chromosome 16
		# no changes
	}
	# no changes for Chromosome MT
	return $new;
}



sub version_20060204_20070113 {
	my ($chromo, $position) = @_;
	my $new = $position;
	if ($chromo eq 'chrI') { # Chromosome 1
		# no changes
	} elsif ($chromo eq 'chrII') { # Chromosome 2
		# no changes
	} elsif ($chromo eq 'chrIII') { # Chromosome 3
		# no changes
	} elsif ($chromo eq 'chrIV') { # Chromosome 4
		#if ($position >= 478023) {$new += 1}
		# comment this out for Segal's 2006 version only
	} elsif ($chromo eq 'chrV') { # Chromosome 5
		# no changes
	} elsif ($chromo eq 'chrVI') { # Chromosome 6
		# no changes
	} elsif ($chromo eq 'chrVII') { # Chromosome 7
		# no changes
	} elsif ($chromo eq 'chrVIII') { # Chromosome 8
		# no changes
	} elsif ($chromo eq 'chrIX') { # Chromosome 9
		# no changes
	} elsif ($chromo eq 'chrX') { # Chromosome 10
		if ($position >= 120805) {$new += 6}
		if ($position >= 120891) {$new += 72}
	} elsif ($chromo eq 'chrXI') { # Chromosome 11
		# no changes
	} elsif ($chromo eq 'chrXII') { # Chromosome 12
		# no changes
	} elsif ($chromo eq 'chrXIII') { # Chromosome 13
		# no changes
	} elsif ($chromo eq 'chrXIV') { # Chromosome 14
		# no changes
	} elsif ($chromo eq 'chrXV') { # Chromosome 15
		# no changes
	} elsif ($chromo eq 'chrXVI') { # Chromosome 16
		# no changes
	}
	# no changes for Chromosome MT
	return $new;
}


sub version_200601xx_20070113 {
	my ($chromo, $position) = @_;
	my $new = $position;
	if ($chromo eq 'chrI') { # Chromosome 1
		# no changes
	} elsif ($chromo eq 'chrII') { # Chromosome 2
		# no changes
	} elsif ($chromo eq 'chrIII') { # Chromosome 3
		# no changes
	} elsif ($chromo eq 'chrIV') { # Chromosome 4
		if ($position >= 478023) {$new += 1}
	} elsif ($chromo eq 'chrV') { # Chromosome 5
		# no changes
	} elsif ($chromo eq 'chrVI') { # Chromosome 6
		# no changes
	} elsif ($chromo eq 'chrVII') { # Chromosome 7
		# no changes
	} elsif ($chromo eq 'chrVIII') { # Chromosome 8
		# no changes
	} elsif ($chromo eq 'chrIX') { # Chromosome 9
		# no changes
	} elsif ($chromo eq 'chrX') { # Chromosome 10
		if ($position >= 120805) {$new += 6}
		if ($position >= 120891) {$new += 72}
	} elsif ($chromo eq 'chrXI') { # Chromosome 11
		# no changes
	} elsif ($chromo eq 'chrXII') { # Chromosome 12
		# no changes
	} elsif ($chromo eq 'chrXIII') { # Chromosome 13
		# no changes
	} elsif ($chromo eq 'chrXIV') { # Chromosome 14
		# no changes
	} elsif ($chromo eq 'chrXV') { # Chromosome 15
		# no changes
	} elsif ($chromo eq 'chrXVI') { # Chromosome 16
		# no changes
	}
	# no changes for Chromosome MT
	return $new;
}

sub version_20060916_20070113 {
	my ($chromo, $position) = @_;
	my $new = $position;
	if ($chromo eq 'chrI') { # Chromosome 1
		# no changes
	} elsif ($chromo eq 'chrII') { # Chromosome 2
		# no changes
	} elsif ($chromo eq 'chrIII') { # Chromosome 3
		# no changes
	} elsif ($chromo eq 'chrIV') { # Chromosome 4
		# no changes
	} elsif ($chromo eq 'chrV') { # Chromosome 5
		# no changes
	} elsif ($chromo eq 'chrVI') { # Chromosome 6
		# no changes
	} elsif ($chromo eq 'chrVII') { # Chromosome 7
		# no changes
	} elsif ($chromo eq 'chrVIII') { # Chromosome 8
		# no changes
	} elsif ($chromo eq 'chrIX') { # Chromosome 9
		# no changes
	} elsif ($chromo eq 'chrX') { # Chromosome 10
		if ($position >= 120805) {$new += 6}
		if ($position >= 120891) {$new += 72}
	} elsif ($chromo eq 'chrXI') { # Chromosome 11
		# no changes
	} elsif ($chromo eq 'chrXII') { # Chromosome 12
		# no changes
	} elsif ($chromo eq 'chrXIII') { # Chromosome 13
		# no changes
	} elsif ($chromo eq 'chrXIV') { # Chromosome 14
		# no changes
	} elsif ($chromo eq 'chrXV') { # Chromosome 15
		# no changes
	} elsif ($chromo eq 'chrXVI') { # Chromosome 16
		# no changes
	}
	# no changes for Chromosome MT
	return $new;
}


sub version_20031001_20070113 {
	my ($chromo, $position) = @_;
	my $new = $position;
	if ($chromo eq 'chrI') { # Chromosome 1
		if ($position >= 41778) {$new += 1}
		if ($position >= 51691) {$new -= 2}
		if ($position >= 130247) {$new += 1}
		if ($position >= 143846) {$new += 1}
		if ($position >= 166768) {$new -= 1}
	} elsif ($chromo eq 'chrII') { # Chromosome 2
		if ($position >= 18364) {$new += 1}
		if ($position >= 21346) {$new += 4}
		if ($position >= 40893) {$new -= 2}
		if ($position >= 93896) {$new += 1}
		if ($position >= 96957) {$new += 37}
		if ($position >= 203175) {$new += 1}
		if ($position >= 216779) {$new -= 1}
		if ($position >= 320082) {$new -= 2}
		if ($position >= 366381) {$new -= 1}
		if ($position >= 457296) {$new += 1}
		if ($position >= 553996) {$new += 1}
		if ($position >= 742173) {$new += 2}
	} elsif ($chromo eq 'chrIII') { # Chromosome 3
		if ($position >= 76147) {$new += 1}
		if ($position >= 101692) {$new += 4}
		if ($position >= 105969) {$new -= 1}
	} elsif ($chromo eq 'chrIV') { # Chromosome 4
		if ($position >= 235713) {$new += 1}
		if ($position >= 478022) {$new += 1}
		if ($position >= 503839) {$new += 1}
		if ($position >= 613203) {$new += 1}
		if ($position >= 819477) {$new += 1}
		if ($position >= 1192051) {$new -= 1}
	# no changes for Chromosome 5
	} elsif ($chromo eq 'chrVI') { # Chromosome 6
		if ($position >= 231670) {$new -= 1}
		if ($position >= 242179) {$new += 1}
	} elsif ($chromo eq 'chrVII') { # Chromosome 7
		if ($position >= 53811) {$new += 1}
		if ($position >= 93036) {$new += 1}
		if ($position >= 93253) {$new += 1}
		if ($position >= 93423) {$new -= 1}
		if ($position >= 124278) {$new += 1}
		if ($position >= 124330) {$new -= 1}
		if ($position >= 130105) {$new += 1}
		if ($position >= 130297) {$new -= 1}
		if ($position >= 130414) {$new -= 1}
		if ($position >= 131043) {$new += 1}
		if ($position >= 506671) {$new += 1}
		if ($position >= 946896) {$new += 1}
		if ($position >= 952554) {$new -= 1}
		if ($position >= 962546) {$new -= 1}
	} elsif ($chromo eq 'chrVIII') { # Chromosome 8
		if ($position >= 54135) {$new += 1}
		if ($position >= 98460) {$new += 1}
		if ($position >= 217748) {$new -= 1}
		if ($position >= 367891) {$new += 1}
		if ($position >= 455307) {$new += 1}
		if ($position >= 455502) {$new += 1}
	# no changes for Chromosome 9
	} elsif ($chromo eq 'chrX') { # Chromosome 10
		if ($position >= 89893) {$new += 1}
		if ($position >= 101724) {$new += 1}
		if ($position >= 118302) {$new -= 1}
		if ($position >= 120804) {$new += 6}
		if ($position >= 120890) {$new += 72}
		if ($position >= 121257) {$new += 220}
		if ($position >= 411267) {$new += 1}
		if ($position >= 460191) {$new -= 1}
	} elsif ($chromo eq 'chrXI') { # Chromosome 11
		if ($position >= 48828) {$new += 1}
		if ($position >= 58939) {$new += 3}
		if ($position >= 185983) {$new -= 2}
		if ($position >= 255405) {$new += 1}
		if ($position >= 363798) {$new += 1}
		if ($position >= 374307) {$new += 1}
		if ($position >= 550828) {$new += 1}
		if ($position >= 552420) {$new += 1}
		if ($position >= 552510) {$new -= 1}
		if ($position >= 603770) {$new += 1}
		if ($position >= 611272) {$new += 1}
		if ($position >= 638899) {$new += 1}
	} elsif ($chromo eq 'chrXII') { # Chromosome 12
		if ($position >= 553548) {$new += 1}
		if ($position >= 553629) {$new += 1}
		if ($position >= 760764) {$new -= 2}
		if ($position >= 899751) {$new += 1}
		if ($position >= 922644) {$new += 1}
	} elsif ($chromo eq 'chrXIII') { # Chromosome 13
		if ($position >= 804684) {$new -= 1}
	} elsif ($chromo eq 'chrXIV') { # Chromosome 14
		if ($position >= 60402) {$new += 1}
		if ($position >= 68346) {$new += 1}
		if ($position >= 254767) {$new += 1}
		if ($position >= 254925) {$new += 1}
		if ($position >= 271648) {$new += 1}
		if ($position >= 472581) {$new -= 1}
		if ($position >= 505511) {$new += 1}
	} elsif ($chromo eq 'chrXV') { # Chromosome 15
		if ($position >= 42534) {$new += 1}
		if ($position >= 803745) {$new += 1}
		if ($position >= 877622) {$new += 2}
	} elsif ($chromo eq 'chrXVI') { # Chromosome 16
		if ($position >= 347285) {$new += 1}
		if ($position >= 347373) {$new += 1}
	}
	# no changes for Chromosome MT
	return $new;
}



sub version_20070113_20100109 {
	my ($chromo, $position) = @_;
	my $new = $position;
	if ($chromo eq 'chrI') { # Chromosome 1
		# no changes
	} elsif ($chromo eq 'chrII') { # Chromosome 2
		# no changes
	} elsif ($chromo eq 'chrIII') { # Chromosome 3
		# no changes
	} elsif ($chromo eq 'chrIV') { # Chromosome 4
		if ($position >= 382379) {$new += 1}
	} elsif ($chromo eq 'chrV') { # Chromosome 5
		# no changes
	} elsif ($chromo eq 'chrVI') { # Chromosome 6
		# no changes
	} elsif ($chromo eq 'chrVII') { # Chromosome 7
		if ($position >= 1038054) {$new += 1}
	} elsif ($chromo eq 'chrVIII') { # Chromosome 8
		# no changes
	} elsif ($chromo eq 'chrIX') { # Chromosome 9
		# no changes
	} elsif ($chromo eq 'chrX') { # Chromosome 10
		if ($position >= 126896) {$new -= 1}
		if ($position >= 599597) {$new -= 1}
		if ($position >= 599640) {$new -= 1}
		if ($position >= 599715) {$new -= 1}
	} elsif ($chromo eq 'chrXI') { # Chromosome 11
		# no changes
	} elsif ($chromo eq 'chrXII') { # Chromosome 12
		# no changes
	} elsif ($chromo eq 'chrXIII') { # Chromosome 13
		# no changes
	} elsif ($chromo eq 'chrXIV') { # Chromosome 14
		if ($position >= 472585) {$new += 1}
	} elsif ($chromo eq 'chrXV') { # Chromosome 15
		# no changes
	} elsif ($chromo eq 'chrXVI') { # Chromosome 16
		# no changes
	}
	# no changes for Chromosome MT
	return $new;
}



sub version_20050806_20100109 {
	my ($chromo, $position) = @_;
	my $new = $position;
	if ($chromo eq 'chrI') { # Chromosome 1
		# no changes
	} elsif ($chromo eq 'chrII') { # Chromosome 2
		# no changes
	} elsif ($chromo eq 'chrIII') { # Chromosome 3
		if ($position >= 76147) {$new += 1}
	} elsif ($chromo eq 'chrIV') { # Chromosome 4
		if ($position >= 382381) {$new += 1}
		if ($position >= 478023) {$new += 1}
		if ($position >= 613205) {$new += 1}
	} elsif ($chromo eq 'chrV') { # Chromosome 5
		# no changes
	} elsif ($chromo eq 'chrVI') { # Chromosome 6
		# no changes
	} elsif ($chromo eq 'chrVII') { # Chromosome 7
		if ($position >= 946899) {$new += 1}
		if ($position >= 962548) {$new -= 1}
		if ($position >= 1038054) {$new += 1}
	} elsif ($chromo eq 'chrVIII') { # Chromosome 8
		if ($position >= 54132) {$new += 1}
	} elsif ($chromo eq 'chrIX') { # Chromosome 9
		# no changes
	} elsif ($chromo eq 'chrX') { # Chromosome 10
		if ($position >= 101726) {$new += 1}
		if ($position >= 120804) {$new += 6}
		if ($position >= 120890) {$new += 72}
		if ($position >= 126817) {$new -= 1}
		if ($position >= 599518) {$new -= 1}
		if ($position >= 599561) {$new -= 1}
		if ($position >= 599636) {$new -= 1}
	} elsif ($chromo eq 'chrXI') { # Chromosome 11
		if ($position >= 185979) {$new -= 1}
		if ($position >= 185987) {$new -= 1}
		if ($position >= 255409) {$new += 1}
		if ($position >= 363801) {$new += 1}
		if ($position >= 552426) {$new += 1}
		if ($position >= 552516) {$new -= 1}
	} elsif ($chromo eq 'chrXII') { # Chromosome 12
		if ($position >= 922648) {$new += 1}
	} elsif ($chromo eq 'chrXIII') { # Chromosome 13
		# no changes
	} elsif ($chromo eq 'chrXIV') { # Chromosome 14
		if ($position >= 60402) {$new += 1}
		if ($position >= 68346) {$new += 1}
		if ($position >= 472579) {$new += 1}
	} elsif ($chromo eq 'chrXV') { # Chromosome 15
		if ($position >= 42534) {$new += 1}
		if ($position >= 803745) {$new += 1}
	} elsif ($chromo eq 'chrXVI') { # Chromosome 16
		# no changes
	}
	# no changes for Chromosome MT
	return $new;
}


sub version_20100109_20110203 {
	my ($chromo, $position) = @_;
	my $new = $position;
	if ($chromo eq 'chrI') {
		$new += -1 if ($position >= 3836);
		$new += 1 if ($position >= 5926); 
		$new += -1 if ($position >= 6454);
		$new += -1 if ($position >= 16468);
		$new += -1 if ($position >= 16498);
		$new += 1 if ($position >= 21530);
		$new += 1 if ($position >= 23240);
		$new += -1 if ($position >= 28069);
		$new += 1 if ($position >= 28948);
		$new += -2 if ($position >= 130488);
		$new += 1 if ($position >= 134124);
		$new += 1 if ($position >= 158933);
		$new += 1 if ($position >= 167086);
		$new += 1 if ($position >= 167113);
		$new += 4 if ($position >= 167135);
		$new += 1 if ($position >= 171882);
		$new += -5 if ($position >= 171969);
		$new += -2 if ($position >= 172018);
		$new += 2 if ($position >= 172042);
		$new += 1 if ($position >= 172169);
		$new += 1 if ($position >= 175305);
		$new += -1 if ($position >= 175372);
		$new += 1 if ($position >= 177136);
		$new += -1 if ($position >= 178602);
		$new += 1 if ($position >= 178617);
		$new += -1 if ($position >= 178639);
		$new += 4 if ($position >= 180960); # manual adjustment
		$new += 1 if ($position >= 198909);
		$new += 1 if ($position >= 202933);
		$new += 1 if ($position >= 203326);
		$new += 1 if ($position >= 229570);
	}
	elsif ($chromo eq 'chrII') {
		$new += 1 if ($position >= 21375);
		$new += -2 if ($position >= 23855);
		$new += 1 if ($position >= 23945);
		$new += 2 if ($position >= 28886); # manual insertion 
		$new += 1 if ($position >= 28906);
		$new += 1 if ($position >= 59405);
		$new += -1 if ($position >= 61051);
		$new += -1 if ($position >= 80743);
		$new += -3 if ($position >= 92916);
		$new += -1 if ($position >= 93620);
		$new += -1 if ($position >= 113437);
		$new += 1 if ($position >= 237020);
		$new += 1 if ($position >= 248622);
		$new += 1 if ($position >= 254543);
		$new += 2 if ($position >= 254720);
		$new += -1 if ($position >= 265902);
		$new += -1 if ($position >= 265920);
		$new += 1 if ($position >= 310456);
		$new += -1 if ($position >= 315071);
		$new += 1 if ($position >= 324361);
		$new += 1 if ($position >= 328378);
		$new += 1 if ($position >= 371419);
		$new += 1 if ($position >= 382728);
		$new += 1 if ($position >= 382740);
		$new += 1 if ($position >= 390032);
		$new += -1 if ($position >= 390442);
		$new += 1 if ($position >= 391272);
		$new += -2 if ($position >= 392341);
		$new += 1 if ($position >= 392818);
		$new += 1 if ($position >= 394725);
		$new += -1 if ($position >= 561435);
		$new += -1 if ($position >= 623445);
		$new += 1 if ($position >= 625434);
		$new += -1 if ($position >= 757291);
		$new += 1 if ($position >= 786760);
		$new += 1 if ($position >= 792275);
	}
	elsif ($chromo eq 'chrIII') {
		$new += 1 if ($position >= 110880);
		$new += 1 if ($position >= 111717);
		$new += 2 if ($position >= 155960); # manual insertion
		$new += -1 if ($position >= 242460);
	}
	elsif ($chromo eq 'chrIV') {
		$new += -1 if ($position >= 34006);
		$new += 1 if ($position >= 215995);
		$new += -1 if ($position >= 330573);
		$new += 1 if ($position >= 332985);
		$new += 2 if ($position >= 369164);
		$new += 1 if ($position >= 371107);
		$new += -1 if ($position >= 381865);
		$new += 1 if ($position >= 396409);
		$new += -1 if ($position >= 396437);
		$new += -1 if ($position >= 402651);
		$new += 1 if ($position >= 402658);
		$new += -1 if ($position >= 542033);
		$new += 1 if ($position >= 569993);
		$new += -1 if ($position >= 578227);
		$new += 1 if ($position >= 592031);
		$new += 1 if ($position >= 894307);
		$new += -1 if ($position >= 1084519);
		$new += 2 if ($position >= 1117109); # manual insertion
		$new += 1 if ($position >= 1140693);
		$new += 1 if ($position >= 1176394);
		$new += 1 if ($position >= 1194953);
		$new += 1 if ($position >= 1502733);
		$new += 1 if ($position >= 1509885);
		$new += 4 if ($position >= 1516809);
		$new += 2 if ($position >= 1516838);
		$new += -1 if ($position >= 1522501);
	}
	elsif ($chromo eq 'chrV') {
		$new += 1 if ($position >= 141111);
		$new += 1 if ($position >= 268856);
		$new += 1 if ($position >= 303531);
		$new += 1 if ($position >= 305709);
		$new += 1 if ($position >= 449958);
	}
	elsif ($chromo eq 'chrVI') {
		$new += 1 if ($position >= 54506);
		$new += 1 if ($position >= 65120);
		$new += 1 if ($position >= 97513);
		$new += 3 if ($position >= 98795); # manual adjustment
		$new += -1 if ($position >= 118606);
		$new += 1 if ($position >= 155994);
		$new += 1 if ($position >= 164823);
		$new += 1 if ($position >= 168868);
		$new += 4 if ($position >= 181038); # manual adjustment
		$new += -1 if ($position >= 192393);
		$new += 1 if ($position >= 194784);
		$new += 1 if ($position >= 226683);
	}
	elsif ($chromo eq 'chrVII') {
		$new += 1 if ($position >= 72711);
		$new += -4 if ($position >= 89587);
		$new += 1 if ($position >= 89600);
		$new += -1 if ($position >= 89610);
		$new += -2 if ($position >= 89749);
		$new += 1 if ($position >= 90017);
		$new += -1 if ($position >= 93660);
		$new += 1 if ($position >= 94594);
		$new += -2 if ($position >= 95082);
		$new += -1 if ($position >= 95444);
		$new += 1 if ($position >= 95470);
		$new += 2 if ($position >= 95478);
		$new += -1 if ($position >= 122568);
		$new += -1 if ($position >= 131261);
		$new += -1 if ($position >= 152961);
		$new += 1 if ($position >= 153008);
		$new += 1 if ($position >= 206465);
		$new += -1 if ($position >= 269097);
		$new += 1 if ($position >= 286354);
		$new += -1 if ($position >= 286727);
		$new += 2 if ($position >= 289526);
		$new += 1 if ($position >= 289566);
		$new += -1 if ($position >= 389965);
		$new += -1 if ($position >= 390007);
		$new += -1 if ($position >= 393811);
		$new += -1 if ($position >= 397868);
		$new += 1 if ($position >= 398013);
		$new += 2 if ($position >= 398027);
		$new += 1 if ($position >= 398065);
		$new += 1 if ($position >= 412392);
		$new += -1 if ($position >= 413089);
		$new += 1 if ($position >= 413409);
		$new += -1 if ($position >= 413784);
		$new += -1 if ($position >= 413969);
		$new += -1 if ($position >= 418633);
		$new += -1 if ($position >= 419054);
		$new += 1 if ($position >= 428248);
		$new += 1 if ($position >= 485888);
		$new += -1 if ($position >= 727145);
		$new += 2 if ($position >= 924486);
		$new += -1 if ($position >= 938725);
		$new += 1 if ($position >= 938791);
		$new += -1 if ($position >= 958924);
		$new += -1 if ($position >= 962981);
		$new += 1 if ($position >= 999335);
		$new += -1 if ($position >= 999676);
		$new += -1 if ($position >= 1004036);
		$new += -1 if ($position >= 1004345);
		$new += 1 if ($position >= 1004545);
		$new += -1 if ($position >= 1039741);
	}
	elsif ($chromo eq 'chrVIII') {
		$new += 1 if ($position >= 2249);
		$new += 1 if ($position >= 8357);
		$new += 1 if ($position >= 18566);
		$new += -1 if ($position >= 25811);
		$new += 1 if ($position >= 64382);
		$new += 1 if ($position >= 78732);
		$new += 1 if ($position >= 102563);
		$new += 1 if ($position >= 126325);
		$new += 1 if ($position >= 131450);
		$new += -1 if ($position >= 207353);
		$new += -8 if ($position >= 212258);
		$new += -2 if ($position >= 287116);
		$new += 1 if ($position >= 287175);
		$new += 1 if ($position >= 293373);
		$new += -1 if ($position >= 294429);
		$new += 1 if ($position >= 369888);
		$new += -1 if ($position >= 369985);
		$new += 1 if ($position >= 423722);
		$new += 1 if ($position >= 441868);
		$new += 1 if ($position >= 445613);
	}
	elsif ($chromo eq 'chrIX') {
		$new += 3 if ($position >= 128404);
		$new += -1 if ($position >= 249637);
		$new += -1 if ($position >= 257479);
		$new += 1 if ($position >= 300238);
		$new += 1 if ($position >= 333320);
		$new += 1 if ($position >= 427191);
		$new += -1 if ($position >= 435031);
	}
	elsif ($chromo eq 'chrX') {
		$new += 1 if ($position >= 59470);
		$new += 1 if ($position >= 79349);
		$new += 2 if ($position >= 97251);
		$new += 1 if ($position >= 97582);
		$new += 1 if ($position >= 97647);
		$new += 4 if ($position >= 99526);
		$new += -2 if ($position >= 99569);
		$new += -2 if ($position >= 99634);
		$new += -1 if ($position >= 110898);
		$new += 1 if ($position >= 110918);
		$new += -1 if ($position >= 111288);
		$new += -1 if ($position >= 111645);
		$new += -1 if ($position >= 180995);
		$new += -1 if ($position >= 187024);
		$new += 1 if ($position >= 195501);
		$new += -1 if ($position >= 198588);
		$new += -1 if ($position >= 198629);
		$new += 1 if ($position >= 204518);
		$new += 3 if ($position >= 204533);
		$new += 1 if ($position >= 247043);
		$new += 1 if ($position >= 407436);
		$new += -1 if ($position >= 422266);
		$new += 1 if ($position >= 422448);
		$new += 1 if ($position >= 438629);
		$new += 2 if ($position >= 613905);
		$new += 1 if ($position >= 613917);
		$new += 1 if ($position >= 622414);
		$new += -1 if ($position >= 627319);
		$new += -1 if ($position >= 628638);
	}
	elsif ($chromo eq 'chrXI') {
		$new += -1 if ($position >= 394);
		$new += -1 if ($position >= 12941);
		$new += 2 if ($position >= 14461); # manual insertion
		$new += -1 if ($position >= 24799);
		$new += -1 if ($position >= 67956);
		$new += -1 if ($position >= 68314);
		$new += -1 if ($position >= 70487);
		$new += -1 if ($position >= 70683);
		$new += -1 if ($position >= 158182);
		$new += 1 if ($position >= 158689);
		$new += -1 if ($position >= 189342);
		$new += 1 if ($position >= 189366);
		$new += 6 if ($position >= 190399);
		$new += 2 if ($position >= 199357);
		$new += 1 if ($position >= 199367);
		$new += 352 if ($position >= 200376); # manual adjustment
		$new += 1 if ($position >= 359603);
		$new += 1 if ($position >= 496922);
		$new += 2 if ($position >= 630164);
		$new += 2 if ($position >= 630275);
	}
	elsif ($chromo eq 'chrXII') {
		$new += 1 if ($position >= 32901);
		$new += -1 if ($position >= 187394);
		$new += -1 if ($position >= 187419);
		$new += -1 if ($position >= 491156);
		$new += 2 if ($position >= 760763); # manual insertion
		$new += 1 if ($position >= 822958);
		$new += -1 if ($position >= 878172);
		$new += 1 if ($position >= 924691);
		$new += 1 if ($position >= 933638);
		$new += 1 if ($position >= 949229);
		$new += 1 if ($position >= 1028607);
		$new += -1 if ($position >= 1032533);
		$new += -1 if ($position >= 1038769);
	}
	elsif ($chromo eq 'chrXIII') {
		$new += -1 if ($position >= 9529);
		$new += 1 if ($position >= 26591);
		$new += 1 if ($position >= 283465);
		$new += 1 if ($position >= 904640);
	}
	elsif ($chromo eq 'chrXIV') {
		$new += 1 if ($position >= 23470);
		$new += -1 if ($position >= 56218);
		$new += -1 if ($position >= 129279);
		$new += -1 if ($position >= 287983);
		$new += -1 if ($position >= 618198);
		$new += 1 if ($position >= 721803);
		$new += 1 if ($position >= 731362);
	}
	elsif ($chromo eq 'chrXV') {
		$new += 1 if ($position >= 8474);
		$new += -1 if ($position >= 58558);
		$new += 1 if ($position >= 218999);
		$new += -1 if ($position >= 220195);
		$new += 1 if ($position >= 220203);
		$new += -2 if ($position >= 259230);
		$new += 1 if ($position >= 343135);
		$new += -1 if ($position >= 414326);
		$new += -1 if ($position >= 423750);
		$new += 1 if ($position >= 492978);
		$new += 1 if ($position >= 521070);
		$new += -1 if ($position >= 521098);
		$new += 1 if ($position >= 588317);
		$new += -1 if ($position >= 588344);
		$new += 1 if ($position >= 808201);
		$new += 1 if ($position >= 816958);
		$new += 1 if ($position >= 877794);
		$new += 1 if ($position >= 892172);
		$new += -1 if ($position >= 1004738);
	}
	elsif ($chromo eq 'chrXVI') {
		$new += 1 if ($position >= 126767);
		$new += 1 if ($position >= 347527);
		$new += 1 if ($position >= 347758);
		$new += -1 if ($position >= 482784);
		$new += 1 if ($position >= 520377);
		$new += 1 if ($position >= 520646);
		$new += -1 if ($position >= 520818);
		$new += -1 if ($position >= 698541);
		$new += 1 if ($position >= 778302);
		$new += 1 if ($position >= 778375);
	}
	
	return $new;
}




sub version_20080606_20110203 {
	my ($chromo, $position) = @_;
	my $new = $position;
	
	if ($chromo eq 'chrI') {
		$new -= 1 if ($position >= 3836);
		$new += 1 if ($position >= 5926);
		$new -= 1 if ($position >= 6454);
		$new -= 1 if ($position >= 16468);
		$new -= 1 if ($position >= 16498);
		$new += 1 if ($position >= 21530);
		$new += 1 if ($position >= 23240);
		$new -= 1 if ($position >= 28069);
		$new += 1 if ($position >= 28948);
		$new -= 2 if ($position >= 130488);
		$new += 1 if ($position >= 134124);
		$new += 1 if ($position >= 158933);
		$new += 1 if ($position >= 167086);
		$new += 1 if ($position >= 167113);
		$new += 4 if ($position >= 167135);
		$new += 1 if ($position >= 171882);
		$new -= 5 if ($position >= 171969);
		$new -= 2 if ($position >= 172018);
		$new += 2 if ($position >= 172042);
		$new += 1 if ($position >= 172169);
		$new += 1 if ($position >= 175305);
		$new -= 1 if ($position >= 175372);
		$new += 1 if ($position >= 177136);
		$new -= 1 if ($position >= 178602);
		$new += 1 if ($position >= 178617);
		$new -= 1 if ($position >= 178639);
		$new += 4 if ($position >= 180960);
		$new += 1 if ($position >= 198909);
		$new += 1 if ($position >= 202933);
		$new += 1 if ($position >= 203326);
		$new += 1 if ($position >= 229570);
	}
	
	elsif ($chromo eq 'chrII') {
		$new += 1 if ($position >= 21375);
		$new -= 2 if ($position >= 23855);
		$new += 1 if ($position >= 23945);
		$new += 2 if ($position >= 28886); # manual insertion 
		$new += 1 if ($position >= 28906);
		$new += 1 if ($position >= 59405);
		$new -= 1 if ($position >= 61051);
		$new -= 1 if ($position >= 80743);
		$new -= 3 if ($position >= 92916);
		$new -= 1 if ($position >= 93620);
		$new -= 1 if ($position >= 113437);
		$new += 1 if ($position >= 237020);
		$new += 1 if ($position >= 248622);
		$new += 1 if ($position >= 254543);
		$new += 2 if ($position >= 254720);
		$new -= 1 if ($position >= 265902);
		$new -= 1 if ($position >= 265920);
		$new += 1 if ($position >= 310456);
		$new -= 1 if ($position >= 315071);
		$new += 1 if ($position >= 324361);
		$new += 1 if ($position >= 328378);
		$new += 1 if ($position >= 371419);
		$new += 1 if ($position >= 382728);
		$new += 1 if ($position >= 382740);
		$new += 1 if ($position >= 390032);
		$new -= 1 if ($position >= 390442);
		$new += 1 if ($position >= 391272);
		$new -= 2 if ($position >= 392341);
		$new += 1 if ($position >= 392818);
		$new += 1 if ($position >= 394725);
		$new -= 1 if ($position >= 561435);
		$new -= 1 if ($position >= 623445);
		$new += 1 if ($position >= 625434);
		$new -= 1 if ($position >= 757291);
		$new += 1 if ($position >= 786760);
		$new += 1 if ($position >= 792275);
	}
	
	elsif ($chromo eq 'chrIII') {
		$new += 1 if ($position >= 110880);
		$new += 2 if ($position >= 155960); # manual insertion
		$new += 1 if ($position >= 111717);
		$new -= 1 if ($position >= 242460);
	}
	
	elsif ($chromo eq 'chrIV') {
		$new -= 1 if ($position >= 34006);
		$new += 1 if ($position >= 215995);
		$new -= 1 if ($position >= 330573);
		$new += 1 if ($position >= 332985);
		$new += 2 if ($position >= 369164);
		$new += 1 if ($position >= 371107);
		$new -= 1 if ($position >= 381865);
		$new += 1 if ($position >= 396409);
		$new -= 1 if ($position >= 396437);
		$new -= 1 if ($position >= 542033);
		$new += 1 if ($position >= 569993);
		$new -= 1 if ($position >= 578227);
		$new += 1 if ($position >= 592031);
		$new += 1 if ($position >= 894307);
		$new -= 1 if ($position >= 1084519);
		$new += 2 if ($position >= 1117109); # manual insertion
		$new += 1 if ($position >= 1140693);
		$new += 1 if ($position >= 1176394);
		$new += 1 if ($position >= 1194953);
		$new += 1 if ($position >= 1502733);
		$new += 1 if ($position >= 1509885);
		$new += 4 if ($position >= 1516809);
		$new += 2 if ($position >= 1516838);
		$new -= 1 if ($position >= 1522501);
	}
	
	elsif ($chromo eq 'chrV') {
		$new += 1 if ($position >= 141111);
		$new += 1 if ($position >= 268856);
		$new += 1 if ($position >= 303531);
		$new += 1 if ($position >= 305709);
		$new += 1 if ($position >= 449958);
	}
	
	elsif ($chromo eq 'chrVI') {
		$new += 1 if ($position >= 54506);
		$new += 1 if ($position >= 65120);
		$new += 1 if ($position >= 97513);
		$new += 3 if ($position >= 98795);
		$new -= 1 if ($position >= 118606);
		$new += 1 if ($position >= 155994);
		$new += 1 if ($position >= 164823);
		$new += 1 if ($position >= 168868);
		$new += 4 if ($position >= 181038);
		$new -= 1 if ($position >= 192393);
		$new += 1 if ($position >= 194784);
		$new += 1 if ($position >= 226683);
	}
	
	elsif ($chromo eq 'chrVII') {
		$new += 1 if ($position >= 72711);
		$new -= 4 if ($position >= 89587);
		$new -= 2 if ($position >= 89749);
		$new += 1 if ($position >= 90017);
		$new -= 1 if ($position >= 93660);
		$new += 1 if ($position >= 94594);
		$new -= 2 if ($position >= 95082);
		$new -= 1 if ($position >= 95444);
		$new += 3 if ($position >= 95472);
		$new -= 1 if ($position >= 122568);
		$new -= 1 if ($position >= 131261);
		$new -= 1 if ($position >= 152961);
		$new += 1 if ($position >= 153008);
		$new += 1 if ($position >= 206465);
		$new -= 1 if ($position >= 269097);
		$new += 1 if ($position >= 286354);
		$new -= 1 if ($position >= 286727);
		$new += 2 if ($position >= 289526);
		$new += 1 if ($position >= 289566);
		$new -= 1 if ($position >= 389965);
		$new -= 1 if ($position >= 390007);
		$new -= 1 if ($position >= 393811);
		$new -= 1 if ($position >= 397868);
		$new += 1 if ($position >= 398013);
		$new += 2 if ($position >= 398027);
		$new += 1 if ($position >= 398065);
		$new += 1 if ($position >= 412392);
		$new -= 1 if ($position >= 413089);
		$new += 1 if ($position >= 413409);
		$new -= 1 if ($position >= 413784);
		$new -= 1 if ($position >= 413969);
		$new -= 1 if ($position >= 418633);
		$new -= 1 if ($position >= 419054);
		$new += 1 if ($position >= 428248);
		$new += 1 if ($position >= 485888);
		$new -= 1 if ($position >= 727145);
		$new += 2 if ($position >= 924486);
		$new -= 1 if ($position >= 938725);
		$new += 1 if ($position >= 938791);
		$new -= 1 if ($position >= 958924);
		$new -= 1 if ($position >= 962981);
		$new += 1 if ($position >= 999335);
		$new -= 1 if ($position >= 999676);
		$new -= 1 if ($position >= 1004036);
		$new -= 1 if ($position >= 1004345);
		$new += 1 if ($position >= 1004545);
		$new -= 1 if ($position >= 1039741);
	}
	
	elsif ($chromo eq 'chrVIII') {
		$new += 1 if ($position >= 2249);
		$new += 1 if ($position >= 8357);
		$new += 1 if ($position >= 18566);
		$new -= 1 if ($position >= 25811);
		$new += 1 if ($position >= 64382);
		$new += 1 if ($position >= 78732);
		$new += 1 if ($position >= 102563);
		$new += 1 if ($position >= 126325);
		$new += 1 if ($position >= 131450);
		$new -= 1 if ($position >= 207353);
		$new -= 8 if ($position >= 212258);
		$new -= 2 if ($position >= 287116);
		$new += 1 if ($position >= 287175);
		$new += 1 if ($position >= 293373);
		$new -= 1 if ($position >= 294429);
		$new += 1 if ($position >= 369888);
		$new -= 1 if ($position >= 369985);
		$new += 1 if ($position >= 423722);
		$new += 1 if ($position >= 441868);
		$new += 1 if ($position >= 445613);
	}
	
	elsif ($chromo eq 'chrIX') {
		$new += 3 if ($position >= 128404);
		$new -= 1 if ($position >= 249637);
		$new -= 1 if ($position >= 257479);
		$new += 1 if ($position >= 300238);
		$new += 1 if ($position >= 333320);
		$new += 1 if ($position >= 427191);
		$new -= 1 if ($position >= 435031);
	}
	
	elsif ($chromo eq 'chrX') {
		$new += 1 if ($position >= 59470);
		$new += 1 if ($position >= 79349);
		$new += 2 if ($position >= 97251);
		$new += 1 if ($position >= 97582);
		$new += 1 if ($position >= 97647);
		$new += 4 if ($position >= 99526);
		$new -= 2 if ($position >= 99569);
		$new -= 2 if ($position >= 99634);
		$new -= 1 if ($position >= 110898);
		$new += 1 if ($position >= 110918);
		$new -= 1 if ($position >= 111288);
		$new -= 1 if ($position >= 111645);
		$new -= 1 if ($position >= 126896);
		$new -= 1 if ($position >= 180996);
		$new -= 1 if ($position >= 187025);
		$new += 1 if ($position >= 195502);
		$new -= 1 if ($position >= 198589);
		$new -= 1 if ($position >= 198630);
		$new += 1 if ($position >= 204519);
		$new += 3 if ($position >= 204534);
		$new += 1 if ($position >= 247044);
		$new += 1 if ($position >= 407437);
		$new -= 1 if ($position >= 422267);
		$new += 1 if ($position >= 422449);
		$new += 1 if ($position >= 438630);
		$new += 2 if ($position >= 613906);
		$new += 1 if ($position >= 613918);
		$new += 1 if ($position >= 622415);
		$new -= 1 if ($position >= 627320);
		$new -= 1 if ($position >= 628639);
	}
	
	elsif ($chromo eq 'chrXI') {
		$new -= 1 if ($position >= 394);
		$new -= 1 if ($position >= 12941);
		$new += 2 if ($position >= 14461); # manual insertion
		$new -= 1 if ($position >= 24799);
		$new -= 1 if ($position >= 67956);
		$new -= 1 if ($position >= 68314);
		$new -= 1 if ($position >= 70487);
		$new -= 1 if ($position >= 70683);
		$new -= 1 if ($position >= 158182);
		$new += 1 if ($position >= 158689);
		$new -= 1 if ($position >= 189342);
		$new += 1 if ($position >= 189366);
		$new += 6 if ($position >= 190399);
		$new += 3 if ($position >= 199357);
		$new += 352 if ($position >= 200376);
		$new += 1 if ($position >= 359603);
		$new += 1 if ($position >= 496922);
		$new += 2 if ($position >= 630164);
		$new += 2 if ($position >= 630275);
	}
	
	elsif ($chromo eq 'chrXII') {
		$new += 1 if ($position >= 32901);
		$new -= 1 if ($position >= 187394);
		$new -= 1 if ($position >= 187419);
		$new -= 1 if ($position >= 491156);
		$new += 2 if ($position >= 760763); # manual insertion
		$new += 1 if ($position >= 822958);
		$new -= 1 if ($position >= 878172);
		$new += 1 if ($position >= 924691);
		$new += 1 if ($position >= 933638);
		$new += 1 if ($position >= 949229);
		$new += 1 if ($position >= 1028607);
		$new -= 1 if ($position >= 1032533);
		$new -= 1 if ($position >= 1038769);
	}
	
	elsif ($chromo eq 'chrXIII') {
		$new -= 1 if ($position >= 9529);
		$new += 1 if ($position >= 26591);
		$new += 1 if ($position >= 283465);
		$new += 1 if ($position >= 904640);
	}
	
	elsif ($chromo eq 'chrXIV') {
		$new += 1 if ($position >= 23470);
		$new -= 1 if ($position >= 56218);
		$new -= 1 if ($position >= 129279);
		$new -= 1 if ($position >= 287983);
		$new += 1 if ($position >= 472581);
		$new -= 1 if ($position >= 618197);
		$new += 1 if ($position >= 721802);
		$new += 1 if ($position >= 731361);
	}
	
	elsif ($chromo eq 'chrXV') {
		$new += 1 if ($position >= 8474);
		$new -= 1 if ($position >= 58558);
		$new += 1 if ($position >= 218999);
		$new -= 2 if ($position >= 259230);
		$new += 1 if ($position >= 343135);
		$new -= 1 if ($position >= 414326);
		$new -= 1 if ($position >= 423750);
		$new += 1 if ($position >= 492978);
		$new += 1 if ($position >= 521070);
		$new -= 1 if ($position >= 521098);
		$new += 1 if ($position >= 588317);
		$new -= 1 if ($position >= 588344);
		$new += 1 if ($position >= 808201);
		$new += 1 if ($position >= 816958);
		$new += 1 if ($position >= 877794);
		$new += 1 if ($position >= 892172);
		$new -= 1 if ($position >= 1004738);
	}
	
	elsif ($chromo eq 'chrXVI') {
		$new += 1 if ($position >= 126767);
		$new += 1 if ($position >= 347527);
		$new += 1 if ($position >= 347758);
		$new -= 1 if ($position >= 482784);
		$new += 1 if ($position >= 520377);
		$new += 1 if ($position >= 520646);
		$new -= 1 if ($position >= 520818);
		$new -= 1 if ($position >= 698541);
		$new += 1 if ($position >= 778302);
		$new += 1 if ($position >= 778375);
	}
	
	return $new;
}



__END__

=head1 NAME

convert_yeast_genome_version.pl

A script to convert genomic data between different S cerevisiae versions

=head1 SYNOPSIS

convert_yeast_genome_version.pl [options] <file1> <file2> ...
  
  Options:
  --in <file>
  --convert <integer>
  --roman
  --db <text>
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <file>

Provide the input file to be converted. Supported files include GFF, 
BED, BedGraph, SGR, and custom text files with chromosome, start, stop, 
and/or position columns. The columns in custom text files must be 
labeled with a header line. BigWig files are also supported by going 
through a conversion step.

More than one file may be provided. Files may be gzipped.

=item --convert <integer>

If you happen to know the custom index number for the conversion method, 
then Great! you can put it here to save time. Otherwise, the program 
will interactively present a list of available conversions from which 
to select.

=item --roman

Boolean flag to indicate whether or not the chromosome name should be 
converted from Arabic (normal) numerals to Roman numerals (used by SGD). 
The default is false.

=item --db <text>

If you're converting a bigWig file, provide a database name to collect 
chromosome names and sizes.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This script will shift coordinates due to indels from one genome release to 
another. It will accept GFF, BED, and bigWig files. It was primarily designed 
to work with genomic data (microarray values, etc) and simple segments, not 
complex features. It will convert both start and stop coordinates, ignoring 
strand. It does not take into account changes in the sequence or coding 
potential.

The chromosome name may also be converted to standard Roman numerals as 
used by SGD. Multiple chromosome name styles are supported: chr1, chr01, 
chrI, 1, or I. 

A new converted file will be written, maintaining the original format, 
with the file name appended with the new format date.

Only a limited number of version conversions are available. 
New conversions may be generated by comparing genomic sequence using the 
EMBOSS utility diffseq, with a word size of 8 or 10 bp, then doing some 
magic manipulations to identify the indels and calculate shifts, often with 
manual inspection of alignments. This is not a trivial process.

Check SGD (http://www.yeastgenome.org) to identify and obtain the available 
yeast genome version releases. In some cases, multiple conversions 
may be necessary to reach the desired version.

=head1 AVAILABLE CONVERSIONS

1) 20050806 (SGD R43) to 20070113 (SGD R55)

2) 20051106 (SGD R46) to 20070113 (SGD R55)

3) 20060204 (SGD R52) to 20070113 (SGD R55)

4) 20031001 (SGD R27, UCSC SacCer1) to 20070113 (SGD R55)

5) 200601xx (SGD R52?) to 20070113 (SGD R55)

6) 20060916 (Pugh's GeneTrack, SGD R53) to 20070113 (SGD R55)

7) 20070113 (SGD R55) to 20100109 (SGD R63)

8) 20050806 (perocchi..steinmetz 2007, SGD R43?) to 20100109 (SGD R63)

9) 20070901 (xu..steinmetz 2009, SGD R56) to 20100109 (SGD R63)

10) 20100109 (SGD R63) to 20110203 (SGD R64, UCSC SacCer3)

11) 20080606 (SGD R61, UCSC SacCer2) to 20110203 (SGD R64, UCSC SacCer3)

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
