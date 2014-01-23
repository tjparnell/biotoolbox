#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Copy;
use File::Spec;
use File::Path 'make_path';
use File::Basename qw(fileparse);
use Bio::ToolBox::data_helper qw(generate_tim_data_structure);
use Bio::ToolBox::file_helper qw(
	write_tim_data_file
	open_to_write_fh
	convert_genome_data_2_gff_data
);
eval {
	# check for bigWig support
	require Bio::ToolBox::db_helper::bigwig;
	Bio::ToolBox::db_helper::bigwig->import;
	
	# check for bigBed support
	require Bio::ToolBox::db_helper::bigbed;
	Bio::ToolBox::db_helper::bigbed->import;
};
eval {
	# check for bam support
	require Bio::ToolBox::db_helper::bam;
	Bio::ToolBox::db_helper::bam->import;
};
my $VERSION = '1.14';

print "\n This script will generate a GFF3 file for BigBed, BigWig or Bam files\n";


### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Get command line options and initialize values
my (
	$path,
	$source,
	$rename,
	$include_chromo,
	$write_metadata,
	$set_name,
	$write_conf,
	$help,
	$print_version,
);
my @infiles;
my @types;
my @names;
my @strands;


# Command line options
GetOptions( 
	'in=s'      => \@infiles, # name of input files
	'path=s'    => \$path, # path to move the bigwig file
	'source=s'  => \$source, # the gff source
	'type=s'    => \@types, # the gff type
	'name=s'    => \@names, # the gff name 
	'strand=s'  => \@strands, # indicate the strand for the feature
	'rename'    => \$rename, # rename the file
	'chromo'    => \$include_chromo, # include chromosomes in GFF
	'set!'      => \$write_metadata, # write a metadata index file for BigWigSet
	'setname=s' => \$set_name, # name for the bigwigset
	'conf!'     => \$write_conf, # write GBrowse conf stanzas
	'help'      => \$help, # request help
	'version'   => \$print_version, # print the version
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
	print " Biotoolbox script big_file2gff3.pl, version $VERSION\n\n";
	exit;
}



### Check for general required values
# files
if (@infiles) {
	if (scalar @infiles == 1 and $infiles[0] =~ /,/) {
		# a comma delimited list is provided
		my $file = shift @infiles;
		@infiles = split /,/, $file;
	}
}
else {
	# file list was provided on the command line
	@infiles = @ARGV or
		die "  OOPS! No source data files specified! \n use --help\n";
}

# types
if (scalar @types == 1 and $types[0] =~ /,/) {
	# a comma delimited list is provided
	my $type = shift @types;
	@types = split /,/, $type;
	if (scalar @types != scalar @infiles) {
		die " unequal number of types (" . scalar(@types) . ") and files (" . 
			scalar(@infiles) . ") provided!\n";
	}
}

# names
if (scalar @names == 1 and $names[0] =~ /,/) {
	# a comma delimited list is provided
	my $name = shift @names;
	@names = split /,/, $name;
	if (scalar @names != scalar @infiles) {
		die " unequal number of names (" . scalar(@names) . ") and files (" . 
			scalar(@infiles) . ") provided!\n";
	}
}

# strands
if (scalar @strands == 1 and $strands[0] =~ /,/) {
	# a separate strand value given for each file as comma delimited list
	my $strand = shift @strands;
	@strands = split /,/, $strand;
	if (scalar @strands != scalar @infiles) {
		die " unequal number of strands (" . scalar(@strands) . ") and files (" . 
			scalar(@infiles) . ") provided!\n";
	}
}
elsif (scalar @strands == 1) {
	# a single strand value
	# assign to each input file name
	my $strand = shift @strands;
	foreach (@infiles) {
		push @strands, $strand;
	}
}
if (@strands) {
	# check the strand values
	foreach my $strand (@strands) {
		unless ($strand =~ /^[frwc\.\-\+]$/) {
			die " unknown strand value '$strand'!\n use --help\n";
		}
	}
}


# target directory
if (defined $path) {
	
	# clean up path as necessary
	$path = File::Spec->rel2abs($path);
	$path = File::Spec->canonpath($path);
	
	# add the set name to the path to make a subdirectory
	if ($set_name) {
		unless ($path =~ m/$set_name\Z/) {
			$path = File::Spec->catdir($path, $set_name);
		}
	}
	
	unless (-e $path) {
		make_path($path) or die "unable to generate target directory: '$path'";
	}
	unless (-w $path) {
		die " target '$path' doesn't seem to be writeable!\n";
	}
}
else {
	# default is to use the current working directory
	$path = File::Spec->curdir();
	$path = File::Spec->rel2abs($path);
}

# my bigwigset name
unless ($set_name) {
	# we'll assume that the last directory in the defined path is the 
	# name we'll use for the set
	my @dirs = File::Spec->splitdir($path);
	$set_name = $dirs[-1];
}

# default source
unless ($source) {
	$source = 'data';
}


# prepare metadata and conf output arrays
my @metadata;
	# we'll be dumping the metadata in here for writing later after going 
	# through all the input files
	# it will be a simple array of lines to be written (or appended) to a 
	# metadata index file written to the output file path
my @confdata;
	# for the conf array, we'll either be writing individual GBrowse 
	# conf stanzas for each bigwig and bam file or one stanza for the bigwigset
my @chromodata;
	# if user wants chromosome features included in the GFF3 file,
	# this information will be added in here
my @subtracks;
	# for bigwigset databases, we'll defer writing the GBrowse conf stanzas 
	# until we've processed all of the bigwig files, and then we'll write a 
	# subtrack table in one bigwigset stanza
my $chromo_check = 0;
	# remember whether we've written the chromosome GFF3 lines as requested



### Processing each of the input files
while (@infiles) {
	
	# get infile name
	my $infile = shift @infiles;
	print "\n processing '$infile'...\n";
	
	
	### Determine file specific GFF variables
	# determine input file components
	my ($infile_basename, $infile_path, $infile_ext) = 
		fileparse($infile, '.bb', '.bw', '.bam');
		# these correspond to bigbed, bigwig, and bam
	unless ($infile_ext) {
		warn "   unrecognized file type!\n";
		next;
	}
	
	# determine gff name
	my $name;
	if (@names) {
		# hopefully the same number of names was provided as files!
		$name = shift @names;
	}
	else {
		# default is to use the base name
		$name = $infile_basename;
	}
	$name =~ s/[_\.\-]?sort(?:ed)?//i; # remove the sorted name if present
	
	# generate a display name without underscores, periods
	my $display_name = $name;
	$display_name =~ s/_|\./ /g; # substitute any underscores or periods with spaces
	
	
	# determine gff type
	my $gfftype;
	if (@types) {
		$gfftype = shift @types;
	}
	else {
		# default to use the name
		$gfftype = $infile_basename;
	}
	$gfftype =~ s/[\-\.\s]+/_/g; # substitute dashes, periods, spaces with underscores
	
	# determine strand
	my $strand;
	if (@strands) {
		$strand =  shift @strands;
	}
	else {
		# default is no strand
		$strand = '.';
	}
	
	### Generate the GFF data structure
	# Initialize the data structure
	my $main_data_ref = generate_tim_data_structure(
		'bigwig_features',
		'Chromosome',
		'Start',
		'Stop',
		'Strand',
		'File'
	) or die " unable to generate tim data structure!\n";
	
	# Determine the target file name
	my $target_basename;
	if ($rename) {
		$target_basename = "$source.$gfftype";
	}
	else {
		$target_basename = $infile_basename;
	}
	my $target_file = 
		File::Spec->catfile($path, "$target_basename$infile_ext");
	
	# Load the chromosome data
	if ($infile_ext eq '.bw') {
		# source data is a bigwig file
		collect_chromosomes_from_bigwig(
			$infile, $target_file, $main_data_ref, $strand);
	}
	elsif ($infile_ext eq '.bb') {
		# source data is a bigbed file
		collect_chromosomes_from_bigbed(
			$infile, $target_file, $main_data_ref, $strand);
	}
	elsif ($infile_ext eq '.bam') {
		# source data is a bam file
		collect_chromosomes_from_bam(
			$infile, $target_file, $main_data_ref, $strand);
	}
	
	# check that we have stored the chromosome information
	if ($include_chromo) {
		unless ($chromo_check) {
			# chromosome information should automatically have been loaded
			# from the collect_chromosomes_from_xxx subs
			$chromo_check = 1;
		}
	}
	
	
	### Move the Input files
	if (-e $target_file) {
		print "  target file already present\n";
	}
	else {
		# copy the file if isn't there already
		copy($infile, $target_file);
		
		# check
		if (-e $target_file) {
			print "  Copied target file '$target_file'\n";
		}
		else {
			warn "  attempted to copy target file but failed!?\n";
		}
	}
	
	
	### Write the GFF3 file
	
	# set the output gff file name
	my $gff_file;
	if ($write_metadata) {
		$gff_file = "$set_name\.gff3";
	} else {
		$gff_file = "$gfftype\.gff3";
	}
	
	# convert the data structure to GFF for writing
	convert_genome_data_2_gff_data(
		'data'       => $main_data_ref,
		'version'    => 3,
		'source'     => $source,
		'type'       => $gfftype,
		'name'       => $name,
		'strand'     => 3,
		'tags'       => [4],
	) or die " unable to convert data to GFF format!\n";
	
	# write new or append existing GFF file
	if (-e $gff_file) {
		# file exists, append to it
		my $gff_fh = open_to_write_fh($gff_file, 0, 1) or 
			die " can't append to GFF file!\n";
		
		# append the table contents to the file
		for my $row (1 .. $main_data_ref->{'last_row'}) {
			print {$gff_fh} join("\t", 
				@{ $main_data_ref->{'data_table'}->[$row] }), "\n";
		}
		$gff_fh->close;
		print "  appended to GFF3 file '$gff_file'\n";
	}
	else {
		# file doesn't exist, write new one
		
		# prepend the chromosome information if requested
		if ($include_chromo) {
			# insert the chromosome gff information into the now-converted
			# gff data table
			splice(
				@{ $main_data_ref->{'data_table'} },
				1,
				0,
				@chromodata
			);
			$main_data_ref->{'last_row'} += scalar(@chromodata);
		}
		
		my $success = write_tim_data_file(
			'data'       => $main_data_ref,
			'filename'   => $gff_file,
		);
		if ($success) {
			print "  wrote GFF3 file '$success'\n";
		}
		else {
			print "  unable to write output file!\n";
		}
	}
	
	
	### Add metadata to index file
	if ($write_metadata and $infile_ext eq '.bw') {
		# we'll be adding all of the input bigwig files to the metadata index
		
		# add metadata header block
			# this is the file name in square brackets
		push @metadata, "[$target_basename$infile_ext]\n";
		
		# add metadata
		push @metadata, "type         = $gfftype\n";
		push @metadata, "source       = $source\n";
		push @metadata, "display_name = $display_name\n";
		if ($strand =~ /^f|w|\+|1/) {
			push @metadata, "strand       = 1\n";
		}
		elsif ($strand =~ /^r|c|\-/) {
			push @metadata, "strand       = -1\n";
		}
		else {
			# in case user wants to change it later
			push @metadata, "strand       = 0\n";
		}
		push @metadata, "\n"; # empty line to make things purdy
	}
	
	### Write conf stanza data
	if ($write_conf) {
		
		# here we write individual stanzas for each big file
		
		# stanzas for a single bigwig files
		if ($infile_ext eq '.bw' and !$write_metadata) {
			# we are using the write metadata boolean variable as an indicator 
			# of whether to write an individual conf stanza for each bigfile
			# file or write a single stanza for the bigwig set
			
			# add the database stanza
			push @confdata, "[$source\_$target_basename\_db:database]\n";
			push @confdata, "db_adaptor   = Bio::DB::BigWig\n";
			push @confdata, "db_args      = -bigwig $target_file\n";
			
			# add the basic track stanza
			push @confdata, "[$source\_$name]\n";
			push @confdata, "database     = $source\_$target_basename\_db\n";
			push @confdata, "feature      = summary\n";
			push @confdata, "glyph        = wiggle_whiskers\n";
			push @confdata, "graph_type   = boxes\n";
			push @confdata, "autoscale    = global_clipped\n";
			push @confdata, "# min_score  = 0\n";
			push @confdata, "# max_score  = 50\n";
			push @confdata, "height       = 50\n";
			push @confdata, "key          = $display_name\n";
			push @confdata, "category     = $set_name\n";
			push @confdata, "citation     = Data file $infile_basename$infile_ext\n";
			push @confdata, "\n\n";
		}
		
		# stanzas for a bigwig set
		elsif ($infile_ext eq '.bw' and $write_metadata) {
			# store information for this bigwig to generate subtrack
			# tables later
			push @subtracks, [$gfftype, $name];
		}
		
		# stanzas for a bigbed file
		elsif ($infile_ext eq '.bb') {
			# add the database stanza
			push @confdata, "[$source\_$target_basename\_db:database]\n";
			push @confdata, "db_adaptor   = Bio::DB::BigBed\n";
			push @confdata, "db_args      = -bigbed $target_file\n";
			
			# add the basic track stanza
				# this is assuming you want the individual features and not 
				# binned summary of scores
				# if you really want that, why don't you go with bigwig?
			push @confdata, "[$source\_$name]\n";
			push @confdata, "database     = $source\_$target_basename\_db\n";
			push @confdata, "feature      = region\n";
			push @confdata, "glyph        = segments\n";
			push @confdata, "stranded     = 1\n";
			push @confdata, "label        = 1\n";
			push @confdata, "key          = $display_name\n";
			push @confdata, "category     = $set_name\n";
			push @confdata, "citation     = Data file $infile_basename$infile_ext\n";
			
			push @confdata, "\n\n";
		
		}
		elsif ($infile_ext eq '.bam') {
			# add the database stanza
			push @confdata, "[$source\_$target_basename\_db:database]\n";
			push @confdata, "db_adaptor     = Bio::DB::Sam\n";
			push @confdata, "db_args        = -bam $target_file\n";
			push @confdata, "#                -fasta /path/to/your/fasta_file.fa\n";
			
			# add the segments track stanza
				# this is assuming single-end alignments
			push @confdata, "[$source\_$name]\n";
			push @confdata, "database       = $source\_$target_basename\_db\n";
			push @confdata, "feature        = match\n";
			push @confdata, "glyph          = segments\n";
			push @confdata, "stranded       = 1\n";
			push @confdata, "label          = 1\n";
			push @confdata, "draw_target    = 1\n";
			push @confdata, "show_mismatch  = 1\n";
			push @confdata, "mismatch_color = red\n";
			push @confdata, "bgcolor        = blue\n";
			push @confdata, "fgcolor        = white\n";
			push @confdata, "key            = $display_name\n";
			push @confdata, "category       = $set_name\n";
			push @confdata, "citation       = Data file $infile_basename$infile_ext\n";
			push @confdata, "\n";
			
			# add the coverage track stanza
				# semantic zooming for above 1000 bp
			push @confdata, "[$source\_$name:1001]\n";
			push @confdata, "feature        = coverage\n";
			push @confdata, "glyph          = wiggle_xyplot\n";
			push @confdata, "graph_type     = boxes\n";
			push @confdata, "autoscale      = local\n";
			push @confdata, "# min_score      = 0\n";
			push @confdata, "# max_score      = 50\n";
			push @confdata, "height         = 50\n";
			push @confdata, "\n\n";
		
		}
	}
}


### Write metadata file
if ($write_metadata) {
	my $md_file = File::Spec->catfile($path, 'metadata.index');
	my $md_fh;
	if (-e $md_file) {
		# file already exists!?
		# append to the file, no compression
		$md_fh = open_to_write_fh($md_file, 0, 1) or 
			die " can't append to metadata index file!\n";
	}
	else {
		# write a new file
		$md_fh = open_to_write_fh($md_file) or 
			die " can't write metadata index file!\n";
	}
	foreach (@metadata) {
		print {$md_fh} $_;
	}
	$md_fh->close;
	print " wrote metadata index file '$md_file'\n";
}



### write the starter conf data
if ($write_conf) {
	
	
	# Generate the stanzas for a bigwigset if necessary
	if ($write_metadata) {
		# write a conf stanza for the bigwig set
		
		# generate stanza name
		my $stanza_name;
		if ($set_name eq $source) {
			# avoid duplication 
			$stanza_name = $source;
		}
		else {
			$stanza_name = "$source\_$set_name";
		}
		
		# add the database stanza, in reverse order
		push @confdata, "[$stanza_name\_db:database]\n";
		push @confdata, "db_adaptor   = Bio::DB::BigWigSet\n";
		push @confdata, "db_args      = -dir $path\n";
		push @confdata, "               -feature_type summary\n\n";
		
		# generate the conf stanzas 
		push @confdata, "[$stanza_name]\n";
		push @confdata, "database     = $stanza_name\_db\n";
		push @confdata, "feature      = " . 
			# using the gfftype
			join(" ", map {$_->[0]} @subtracks) . "\n";
		push @confdata, "subtrack select = Feature type\n";
		for my $i (0 .. $#subtracks) {
			# create subtrack table 
			# each item has gfftype and name in anon array
			# use the name as the label, it could have spaces
			# the gfftype will have no whitespace
			if ($i == 0) {
				push @confdata, "subtrack table = " . 
					":\"$subtracks[$i][1]\" $subtracks[$i][0] " . 
					"=$subtracks[$i][0] ;\n";
			}
			else {
				push @confdata, "               " . 
					":\"$subtracks[$i][1]\" $subtracks[$i][0] " . 
					"=$subtracks[$i][0] ;\n";
			}
		}
		for my $i (0 .. $#subtracks) {
			# subtrack select labels
			# just to try and prettify the names a little, not really necessary
			if ($i == 0) {
				push @confdata, "subtrack select labels = " . 
					"$subtracks[$i][0] \"$subtracks[$i][1]\" ;\n";
			}
			else {
				push @confdata, "               " . 
					"$subtracks[$i][0] \"$subtracks[$i][1]\" ;\n";
			}
		}
		push @confdata, "glyph        = wiggle_whiskers\n";
		push @confdata, "graph_type   = boxes\n";
		push @confdata, "height       = 50\n";
		push @confdata, "scale        = right\n";
		push @confdata, "# min_score  = 0\n";
		push @confdata, "# max_score  = 50\n";
		push @confdata, "autoscale    = global_clipped\n";
		push @confdata, "label        = 0\n";
		push @confdata, "group_label  = 1\n";
		push @confdata, "group_label_position = top\n";
		my $key = $set_name;
		$key =~ s/_/ /g;
		push @confdata, "key          = $key\n";
		push @confdata, "category     = $set_name\n";
		push @confdata, "citation     = Data file collection from $set_name\n";
		push @confdata, "\n\n";
	}
	
	# write the conf stanza file
	my $conf_file = 'conf_stanzas.txt';
	my $conf_fh;
	if (-e $conf_file) {
		# file already exists!?
		# append to it
		$conf_fh = open_to_write_fh($conf_file, 0, 1) or 
			die " can't append to GBrowse configuration stanza file!\n";
	}
	else {
		# write a new file
		$conf_fh = open_to_write_fh($conf_file) or 
			die " can't write GBrowse configuration stanza file!\n";
	}
	foreach (@confdata) {
		print {$conf_fh} $_;
	}
	$conf_fh->close;
	print " wrote sample GBrowse configuration stanza file '$conf_file'\n";
}
print " Finished\n";






#############################  Subroutines  ################################


sub collect_chromosomes_from_bigwig {
	my ($infile, $target_file, $data_ref, $strand) = @_;
	
	# adjust the File column name
	$data_ref->{4}{'name'} = 'bigwigfile';
	$data_ref->{'data_table'}->[0][4] = 'bigwigfile';
	
	# open the bigwig file
	unless (exists &open_bigwig_db) {
		die " unable to load BigWig file support! Is Bio::DB::BigWig installed?\n"; 
	}
	my $wig = open_bigwig_db($infile) or 
		die " unable to open bigwig file '$infile'!\n";
	
	
	# collect the chromosomes
	my @seq_ids = $wig->seq_ids;
	unless (@seq_ids) {
		die " no chromosomes in the bigwig file!?\n";
	}
	
	# fill out the data table
	foreach my $seq_id (@seq_ids) {
		push @{ $data_ref->{'data_table'} }, [ (
			$seq_id,
			1,
			$wig->length($seq_id),
			$strand,
			$target_file,
		) ];
	}
	
	# update
	$data_ref->{'last_row'} = scalar(@seq_ids);
	
	# record chromosome information
	if ($include_chromo and !$chromo_check) {
		foreach my $seq_id (@seq_ids) {
			push @chromodata, [ 
				(
					$seq_id,
					'.',
					'chromosome',
					1,
					$wig->length($seq_id),
					'.',
					'.',
					'.',
					"Name=$seq_id;ID=$seq_id",
				)
			];
		}
		$chromo_check = 1;
	}
}



sub collect_chromosomes_from_bigbed {
	my ($infile, $target_file, $data_ref, $strand) = @_;
	
	# adjust the File column name
	$data_ref->{4}{'name'} = 'bigbedfile';
	$data_ref->{'data_table'}->[0][4] = 'bigbedfile';
	
	# open the bigbed file
	unless (exists &open_bigbed_db) {
		die " unable to load BigBed file support! Is Bio::DB::BigBed installed?\n"; 
	}
	my $bb = open_bigbed_db($infile) or 
		die " unable to open bigbed file '$infile'!\n";
	
	# collect the chromosomes
		# BigBed.pm doesn't provide handy seq_ids() and length() methods 
		# like BigWig.pm does....
		# so we have to fake it using underlying BigFile methods
	my $chromList = $bb->bb->chromList or 
		die " unable to collect chromosomes from BigBed BigFile!\n";
	# this returns an object representing the chromosomes in the BigFile
	
	# fill out the data table
	# this follows an example from BigFile.pm pod
	for (my $c = $chromList->head; $c; $c = $c->next) {
		push @{ $data_ref->{'data_table'} }, [ (
			$c->name,
			1,
			$c->size,
			$strand,
			$target_file,
		) ];
	}
	
	# update number
	$data_ref->{'last_row'} = scalar(@{ $data_ref->{'data_table'} }) - 1;
	
	# record chromosome information
	if ($include_chromo and !$chromo_check) {
		for (my $c = $chromList->head; $c; $c = $c->next) {
			my $seq_id = $c->name;
			push @chromodata, [
				(
					$seq_id,
					'.',
					'chromosome',
					1,
					$c->size,
					'.',
					'.',
					'.',
					"Name=$seq_id;ID=$seq_id",
				)
			];
		}
		$chromo_check = 1;
	}
}



sub collect_chromosomes_from_bam {
	my ($infile, $target_file, $data_ref, $strand) = @_;
	
	# adjust the File column name
	$data_ref->{4}{'name'} = 'bamfile';
	$data_ref->{'data_table'}->[0][4] = 'bamfile';
	
	# open the bam file
	unless (exists &open_bam_db) {
		die " unable to load Bam file support! Is Bio::DB::Sam installed?\n"; 
	}
	my $sam = open_bam_db($infile) or 
		die " unable to open bam file '$infile'!\n";
	
	# collect the chromosomes
	my $seq_num = $sam->n_targets;
	unless ($seq_num) {
		die " no chromosomes in the bam file!?\n";
	}
	
	# fill out the data table
	foreach my $tid (0 .. ($seq_num - 1) ) {
		push @{ $data_ref->{'data_table'} }, [ (
			$sam->target_name($tid),
			1,
			$sam->target_len($tid),
			$strand,
			$target_file,
		) ];
	}
	
	# update
	$data_ref->{'last_row'} = $seq_num;
	
	# record chromosome information
	if ($include_chromo and !$chromo_check) {
		foreach my $tid (0 .. ($seq_num - 1) ) {
			my $seq_id = $sam->target_name($tid);
			push @chromodata, [
				(
					$seq_id,
					'.',
					'chromosome',
					1,
					$sam->target_len($tid),
					'.',
					'.',
					'.',
					"Name=$seq_id;ID=$seq_id",
				)
			];
		}
		$chromo_check = 1;
	}
}




__END__

=head1 NAME

big_file2gff3.pl

A script to generate GFF3 files for bigwig, bigbed, and bam files.

=head1 SYNOPSIS

big_file2gff3.pl [--options...] <filename1.bw> <filename2.bb> ...
  
  Options:
  --in <file> or <file1,file2,...>
  --path </destination/path/for/bigfiles/>
  --source <text>
  --name <text> or <text1,text2,...>
  --type <text> or <text1,text2,...>
  --strand [f|r|w|c|+|-|1|0|-1],...
  --chromo 
  --rename
  --set
  --setname <text>
  --conf
  --version
  --help
  
=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename.bw | filename.bb | filename.bam>

Provide the name(s) of the input BigFiles. Three types of files are 
supported: BigWig (.bw), BigBed (.bb), or BAM (.bam) files. 
The list may be specified three ways: by reiterating the 
--in option, providing a single comma-delimited list, or simply listing 
all files after the options (e.g. "data_*.bw"). 


=item --path </destination/path/for/bigfiles/>

Provide the destination directory name for the BigFile files. If the
destination does not exist, then it will created. This directory should be
writeable by the user and readable by all (or at least the Apache and MySQL
users). If the input files are not currently located here, they will be
copied to the directory for you. Note that when generating a BigWigSet, a
subdirectory with the set name (option --setname) will be made for you. The
default path is the current path for the input file.

=item --source <text>

Provide the name for the GFF feature's source for all of the input value. 
The default value is "data". Unique values for each input file is not 
supported.

=item --name <text>

Provide the name(s) for the GFF feature(s). This will be used as the GFF 
feature's name. A unique name should be provided for each file, and may be 
specified as a single comma-delimited list or by reiterating the --name 
option. The default value is to use the input file basename.

=item --type <text>

Provide the GFF type for the GFF features for all input files. A unique 
type should be provided for each file, and may be specified as a single 
comma-delimited list or by reiterating the --type option. By default, 
it re-uses the GFF name. 

=item --strand [f|r|w|c|+|-|1|0|-1],...

Indicate which strand the feature will be located. Acceptable values include 
(f)orward, (r)everse, (w)atson, (c)rick, +, -, 1, or -1. By default no strand 
is used. For mulitple input files with different strands, use the option 
repeatedly or provide a comma-delimited list. If only one value is provided, 
it is used for all input files.

=item --chromo

Include chromosome features in the output GFF3 file. These are necessary 
to make a complete GFF3 file compatible for loading into a Bio::DB database 
when chromosome information is not provided elsewhere.

=item --rename

Rename the input source file basenames as "$source.$name" when moving to the 
target directory. This may help in organization and clarity of file listings. 
Default is false.

=item --set

Indicate that all of the input BigWig files should be part of a BigWigSet,
which treats all of the BigWig files in the target directory as a single
database. A text file is written in the target directory with metadata for 
each BigWig file (feature, source, strand, name) as described in the 
Bio::DB::BigWigSet documentation. Additional metadata may be manually 
added if desired. The default is false.

=item --setname <text>

Optionally specify the name for the BigWigSet track when writing the 
GBrowse configuration stanza. It is also used as the basename for the 
GFF3 file, as well as the name of the new subdirectory in the target path 
for use as the BigWigSet directory. The default is to use the name of 
the last directory in the target path.

=item --conf

Write sample GBrowse database and track configuration stanzas. Each BigFile 
file will get individual stanzas, unless the --set option is enabled, where 
a single stanza with subtracks for the BigWigSet is generated. This is 
helpful when setting up GBrowse database and configurations. Default is false.

=item --version

Print the version number.

=item --help

Display this POD help.

=back

=head1 DESCRIPTION

This program will generate a GFF3 file with features representing a 
BigFile file. Two types of BigFiles are supported: BigWig and BigBed. BAM 
files are also supported. The generated features will encompass the entire 
length of each chromosome represented in the data file. The name of the data 
file and its absolute path is stored as a tag value in each feature. This 
tag value can be used by C<Bio::ToolBox::db_helper> to collect data from the file 
with respect to various locations and features in the database. 

The source data file is copied to the destination directory. The file may be 
renamed as "source.name" so as to avoid confusion when lots of files are 
dumped into the same directory. 

The BigWig files may also be designated as a BigWigSet, with unique metadata 
assigned to each file (source, type, name, strand). The BigWigSet may be 
treated as a single database with multiple bigwig data sources by GBrowse. A 
metadata index file is written in the target directory as described in 
Bio::DB::BigWigSet.

Multiple source files may be designated, and each may have its own name. 
This facilitates multiple file processing.

The generated GFF3 file is written to the current directory. One GFF file is 
written for each input file, or one GFF file for a BigWigSet. It uses the 
provided GFF name as the basename for the file.

Optionally, sample database and track GBrowse configuration stanzas may also be 
written to the current directory to facilitate setting up GBrowse. If a 
BigWigSet database is requested, then the track stanza will be set up with 
subtrack tables representing each BigWig file. 

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
