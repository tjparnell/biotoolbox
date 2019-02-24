#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(verify_or_request_feature_types);
use Bio::ToolBox::GeneTools qw(
	:export 
	:filter
	:transcript
);
use Bio::ToolBox::utility;
my $VERSION = '1.65';

print "\n This program will collect features from annotation sources\n\n";

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
	$input,
	$database,
	$get_subfeatures,
	$include_coordinates,
	$start_adj,
	$stop_adj,
	$position,
	$tsl,
	$gencode,
	$tbiotype,
	$collapse,
	$chromosome_exclude,
	$convert_to_bed,
	$convert_to_gff,
	$convert_to_gtf,
	$convert_to_refflat,
	$outfile,
	$sort_data,
	$gz,
	$bgz,
	$help,
	$print_version,
);
my @features;
my @include_tags;
my @exclude_tags;
my %exclude_tag2value;

# Command line options
GetOptions( 
	'i|in=s'      => \$input, # input table
	'd|db=s'      => \$database, # source annotation database
	'f|feature=s' => \@features, # the features to collect from the database
	'u|sub!'      => \$get_subfeatures, # collect subfeatures
	'coord!'    => \$include_coordinates, # collect coordinates
	'b|start=i'   => \$start_adj, # start coordinate adjustment
	'e|stop=i'    => \$stop_adj, # stop coordinate adjustment
	'p|pos=s'     => \$position, # relative position to adjust coordinates
	't|tag=s'     => \@include_tags, # attributes to include
	'x|exclude=s' => \@exclude_tags, # attribute and keys to exclude
	'tsl=s'     => \$tsl, # filter on transcript support level
	'gencode!'  => \$gencode, # filter on gencode basic tag
	'biotype=s' => \$tbiotype, # filter on transcript biotype
	'collapse!' => \$collapse, # collapse multi-transcript genes
	'K|chrskip=s' => \$chromosome_exclude, # skip chromosomes
	'B|bed!'      => \$convert_to_bed, # convert to bed format
	'G|gff|gff3!' => \$convert_to_gff, # convert to GFF3 format
	'g|gtf!'      => \$convert_to_gtf, # convert to gtf format
	'r|refflat!'  => \$convert_to_refflat, # convert to refFlat format
	'o|out=s'     => \$outfile, # name of output file 
	'sort!'       => \$sort_data, # sort the output file
	'z|gz!'       => \$gz, # compress output
	'Z|bgz!'      => \$bgz, # compress with bgzip
	'h|help'      => \$help, # request help
	'v|version'   => \$print_version, # print the version
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
	print " Biotoolbox script get_features.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}



### Check for requirements
check_requirements();


### Fill Data object
my $Data;
if ($input) {
	$Data = load_from_infile();
} else {
	$Data = load_from_database();
}


### Filter
filter_features();


### Write out file
if ($start_adj or $stop_adj) {
	printf " Adjusting start by %s bp and stop by %s bp relative to position %s\n",
		$start_adj, $stop_adj, $position eq '5' ? 'start' : $position eq '3' ? 'end' : 
		$position eq '4' or $position eq 'm' ? 'middle' : $position eq '53' ? 'both ends' :
		'neither';
}
if ($convert_to_bed) {
	print " Writing to bed file...\n";
	export_to_bed();
} elsif ($convert_to_gff) {
	print " Writing to GFF3 file...\n";
	export_to_gff();
} elsif ($convert_to_gtf) {
	print " Writing to GTF file...\n";
	export_to_gtf();
} elsif ($convert_to_refflat) {
	print " Writing to refFlat file...\n";
	export_to_ucsc();
} else {
	print " Writing to text file...\n";
	export_to_txt();
}





########################   Subroutines   ###################################

sub check_requirements {
	
	# check input
	unless ($database or $input) {
		die " Must provide an input file or database name! use --help for more information\n";
	}
	if ($input =~ /\.(?:sqlite|db)$/i) {
		# whoops! specifiying a database file as input
		$database = $input;
		undef $input;
	}
	
	# check if feature is a comma delimited list
	if (scalar @features == 1 and $features[0] =~ /,/) {
		@features = split ',', shift @features;
	}
	if (scalar(@features) > 1 and $input) {
		warn sprintf(" Only one feature allowed when parsing an input file! Using %s",
			$features[0]);
	}
	if (not @features and $input) {
		print " using default feature of 'gene'\n";
		$features[0] = 'gene';
	}
	
	# check conversions
	my $conversions = $convert_to_bed + $convert_to_gff + $convert_to_gtf + $convert_to_refflat;
	if ($conversions > 1) {
		die " Too many bed/gff/gtf/refFlat conversions specified!\n";
	}
	if ($convert_to_gff or $convert_to_gtf or $convert_to_refflat) {
		$get_subfeatures = 1 if ($input and not defined $get_subfeatures);
	}
	
	# check collapse
	if ($collapse) {
		unless ($get_subfeatures) {
			die " Cannot collapse transcript unless subfeatures are turned on!\n";
		}
		unless ($convert_to_gff or $convert_to_gtf or $convert_to_refflat or 
				$convert_to_bed
		) {
			die " Cannot collapse transcripts unless writing to BED12/GFF/GTF/refFlat!\n";
		}
	}
	
	# check adjustments
	if ($start_adj or $stop_adj) {
		# automatically include coordinates if we're adjusting them
		if ($get_subfeatures) {
			die " Cannot adjust coordinates when including subfeatures!\n";
		}
		$include_coordinates = 1;
	}
	if ($position) {
		unless ($position =~ /[543m]{1,2}/) {
			die " unrecognized position value '$position'! see help\n";
		}
	}
	else {
		# default is from both ends
		$position = '53';
	}
	
	# check include tags
	if (@include_tags and scalar(@include_tags) == 1 and $include_tags[0] =~ /,/) {
		@include_tags = split /,/, shift @include_tags;
	}
	
	# exclude tags
	if (@exclude_tags) {
		foreach (@exclude_tags) {
			my ($k, $v) = split /[,=]/, $_;
			$exclude_tag2value{$k} = $v;
		}
	}
	
	# check output
	unless ($outfile) {
		die " Must provide an output file!\n";
	}
	if ($bgz) {
		$gz = 2;
		$sort_data = 1;
	}
}


sub load_from_database {
	# validate and/or request features
	@features = verify_or_request_feature_types(
		'db'      => $database,
		'feature' => \@features,
		'prompt'  => " Enter the feature(s) to collect." . 
				" A comma de-limited list or range may be given\n",
	) or die " no valid features were provided! see help\n";
	
	# generate a list from the database
	my $D = Bio::ToolBox::Data->new(
		db         => $database,
		feature    => join(',', @features),
		chrskip    => $chromosome_exclude,
	) or die " unable to generate new feature list\n";
	
	if ($D->last_row) {
		printf " Loaded %s features from %s.\n", format_with_commas($D->last_row), 
			$input ? $input : $database;
	}
	else {
		die " No features loaded!\n";
	}
	return $D;
}


sub load_from_infile {
	# parse file
	my $D = Bio::ToolBox::Data->new(
		file       => $input, 
		parse      => 1,
		simplify   => 0, # we want everything!
		feature    => $features[0],
		subfeature => $get_subfeatures ? 'exon,cds,utr,codon' : '',
		chrskip    => $chromosome_exclude,
	) or die " unable to load input file '$input'\n";
	
	if ($D->last_row) {
		printf " Loaded %s features from $input.\n", format_with_commas($D->last_row);
	}
	else {
		die " No features loaded! Re-check your feature type. If you are attempting to \n" . 
			" parse subfeatures like exon or CDS, try the program get_gene_regions instead.\n";
	}
	return $D;
}


sub filter_features {
	# filter on specific tags
	if (%exclude_tag2value) {
		foreach my $k (keys %exclude_tag2value) {
			my $check = $exclude_tag2value{$k};
			print " Filtering tag $k for $check...\n";
			my @unwanted;
			$Data->iterate( sub {
				my $row = shift;
				my ($v) = $row->seqfeature(1)->get_tag_values($k);
				if ($v =~ /$check/i) {
					# the tag matches, so discard
					push @unwanted, $row->row_index;
				}
			});
			if (@unwanted) {
				$Data->delete_row(@unwanted);
			}
		}
		printf "  Kept %s features.\n", format_with_commas($Data->last_row);
	}
	
	# filter on gencode
	if ($tbiotype) {
		print " Filtering transcript biotype for $tbiotype...\n";
		my @unwanted;
		$Data->iterate( sub {
			my $row = shift;
			my $good = filter_transcript_biotype($row->seqfeature(1), $tbiotype);
			unless (defined $good) {
				push @unwanted, $row->row_index;
				next;
			}
			my @t = get_transcripts($good); # verify we have transcripts
			if (scalar @t) {
				$Data->store_seqfeature($row->row_index, $good);
			}
			else {
				push @unwanted, $row->row_index;
			}
		});
		if (@unwanted) {
			$Data->delete_row(@unwanted);
		}
		printf "  Kept %s features.\n", format_with_commas($Data->last_row);
	}
	
	# filter on tsl
	if ($tsl) {
		print " Filtering for transcript support level of $tsl...\n";
		my @unwanted;
		$Data->iterate( sub {
			my $row = shift;
			my $good = filter_transcript_support_level($row->seqfeature(1), $tsl);
			unless (defined $good) {
				push @unwanted, $row->row_index;
				next;
			}
			my @t = get_transcripts($good); # verify we have transcripts
			if (scalar @t) {
				$Data->store_seqfeature($row->row_index, $good);
			}
			else {
				push @unwanted, $row->row_index;
			}
		});
		if (@unwanted) {
			$Data->delete_row(@unwanted);
		}
		printf "  Kept %s features.\n", format_with_commas($Data->last_row);
	}
	
	# filter on gencode
	if ($gencode) {
		print " Filtering for GENCODE transcripts...\n";
		my @unwanted;
		$Data->iterate( sub {
			my $row = shift;
			my $good = filter_transcript_gencode_basic($row->seqfeature(1));
			unless (defined $good) {
				push @unwanted, $row->row_index;
				next;
			}
			my @t = get_transcripts($good); # verify we have transcripts
			if (scalar @t) {
				$Data->store_seqfeature($row->row_index, $good);
			}
			else {
				push @unwanted, $row->row_index;
			}
		});
		if (@unwanted) {
			$Data->delete_row(@unwanted);
		}
		printf "  Kept %s features.\n", format_with_commas($Data->last_row);
	}
	
	# collapse transcripts
	if ($collapse) {
		print " Collapsing alternate transcripts...\n";
		$Data->iterate( sub {
			my $row = shift;
			my $gene = collapse_transcripts($row->seqfeature(1));
			if ($gene) {
				$Data->store_seqfeature($row->row_index, $gene);
			}
		});
	}
}


sub export_to_bed {
	# prepare output
	my $outData = Bio::ToolBox::Data->new(
		bed => $get_subfeatures ? 12 : 6,
	) or die "unable to create output Data structure!\n";
	
	# Write method based on subfeatures or coordinate adjustment
	if ($start_adj or $stop_adj) {
		# adjust coordinates as necessary and write a BED6 file
		$Data->iterate( sub {
			my $row = shift;
			my $f = $row->seqfeature(1); # make sure we get the seqfeature
			my ($start, $stop) = adjust_coordinates($f); 
			my $string = $row->bed_string(
				start => $start,
				end   => $stop,
			);
			$outData->add_row($string);
			$Data->delete_seqfeature($row->row_index);
		});
	}
	elsif ($get_subfeatures) {
		# write transcripts as BED12
		if ($features[0] =~ /gene/i) {
			print " NOTE: gene information is discarded when writing BED12\n";
		}
		$Data->iterate( sub {
			my $row = shift;
			my $f = $row->seqfeature(1); # make sure we get the seqfeature
			my $string = bed_string($f);
			foreach (split /\n/, $string) {
				$outData->add_row($_);
			}
			$Data->delete_seqfeature($row->row_index);
		});
	}
	else {
		# write ordinary BED6
		$Data->iterate( sub {
			my $row = shift;
			my $f = $row->seqfeature(1); # make sure we get the seqfeature
			my $string = $row->bed_string;
			$outData->add_row($string);
			$Data->delete_seqfeature($row->row_index);
		});
	}
	print " done\n";
	
	# sort as necessary
	if ($sort_data) {
		print " Sorting data...\n";
		$outData->gsort_data;
	}
	
	# write
	unless ($outfile =~ /\.bed(?:\.gz)?$/i) {
		$outfile .= '.bed';
	}
	$outfile = $outData->write_file(
		filename => $outfile, 
		gz       => $gz,
	);
	print " wrote file $outfile\n";
}


sub export_to_gff {
	# check output filename
	unless ($outfile =~ /\.gff3?(?:\.gz)?$/i) {
		$outfile .= '.gff3';
	}
	
	# how we write the output depends on whether we need to sort the file or not
	if ($sort_data) {
		# we need to keep all gff lines in memory so that we can sort the lines
		
		my $outData = Bio::ToolBox::Data->new(
			gff => 3,
		) or die "unable to create output Data structure!\n";
		$outData->add_comment( sprintf("exported from %s\n", 
			$database ? $database : $input) );
		
		# iterate
		$Data->iterate( sub {
			my $row = shift;
			my $string = $row->seqfeature(1)->gff_string(1);
				# force retrieving the seqfeature, and recurse through subfeature
			foreach (split /\n/, $string) {
				$outData->add_row($_);
			}
			$Data->delete_seqfeature($row->row_index);
		});
		print " Sorting data...\n";
		$outData->gsort_data;
		$outfile = $outData->write_file(
			filename => $outfile, 
			gz       => $gz,
		);
	}
	else {
		# we can simply write out gff directly
		$outfile .= '.gz' if ($gz and $outfile !~ /\.gz$/i);
		my $fh = Bio::ToolBox::Data->open_to_write_fh($outfile, $gz) or 
			die "unable to open '$outfile' for writing! $!\n";
		$fh->print("##gff-version 3\n");
		$fh->printf("# exported from %s\n", $database ? $database : $input);
	
		# write to GFF
		$Data->iterate( sub {
			my $row = shift;
			my $string = $row->seqfeature(1)->gff_string(1);
				# force retrieving the seqfeature, and recurse through subfeature
			$fh->print( $string . "###\n"); # include pragma close lines
		});
		$fh->close;
	}
	printf " wrote file $outfile\n";
}


sub export_to_gtf {
	# check output filename
	unless ($outfile =~ /\.gtf?(?:\.gz)?$/i) {
		$outfile .= '.gtf';
	}
	
	# how we write the output depends on whether we need to sort the file or not
	if ($sort_data) {
		# we need to keep all gff lines in memory so that we can sort the lines
		
		my $outData = Bio::ToolBox::Data->new(
			gff => 2.5,
		) or die "unable to create output Data structure!\n";
		$outData->add_comment( sprintf("exported from %s\n", 
			$database ? $database : $input) );
		
		# iterate
		$Data->iterate( sub {
			my $row = shift;
			my $string = gtf_string( $row->seqfeature(1) );
			foreach (split /\n/, $string) {
				$outData->add_row($_);
			}
			$Data->delete_seqfeature($row->row_index);
		});
		print " Sorting data...\n";
		$outData->gsort_data;
		$outfile = $outData->write_file(
			filename => $outfile, 
			gz       => $gz,
		);
	}
	else {
		# we can simply write out gff directly
		$outfile .= '.gz' if ($gz and $outfile !~ /\.gz$/i);
		my $fh = Bio::ToolBox::Data->open_to_write_fh($outfile, $gz) or 
			die "unable to open '$outfile' for writing! $!\n";
		$fh->print("##gff-version 2.5\n");
		$fh->printf("# exported from %s\n", $database ? $database : $input);
	
		# write to GTF
		$Data->iterate( sub {
			my $row = shift;
			my $string = gtf_string( $row->seqfeature(1) );
			$fh->print($string);
		});
		$fh->close;
	}
	printf " wrote file $outfile\n";
}


sub export_to_txt {
	# adjust coordinates as necessary
	if ($start_adj or $stop_adj) {
		
		# make sure we establish the feature type first, before we add 
		# coordinate columns, otherwise features might not be returned properly
		my $ftype = $Data->feature_type;
		
		# add coordinate columns
		my $seq_i    = $Data->add_column('Chromosome');
		my $start_i  = $Data->add_column('Start');
		my $stop_i   = $Data->add_column('Stop');
		my $strand_i = $Data->add_column('Strand');
		
		# iterate
		$Data->iterate( sub {
			my $row = shift;
			my $f = $row->seqfeature(1);
			my ($start, $stop) = adjust_coordinates($f); 
			$row->value($seq_i, $f->seq_id);
			$row->value($start_i, $start);
			$row->value($stop_i, $stop);
			$row->value($strand_i, $f->strand);
		});
	}
	elsif ($include_coordinates) {
		# just include original coordinates
		
		# make sure we establish the feature type first, before we add 
		# coordinate columns, otherwise features might not be returned properly
		my $ftype = $Data->feature_type;
		
		# add coordinate columns
		my $seq_i    = $Data->add_column('Chromosome');
		my $start_i  = $Data->add_column('Start');
		my $stop_i   = $Data->add_column('Stop');
		my $strand_i = $Data->add_column('Strand');
		
		# iterate
		$Data->iterate( sub {
			my $row = shift;
			my $f = $row->seqfeature(1);
			$row->value($seq_i, $f->seq_id);
			$row->value($start_i, $f->start);
			$row->value($stop_i, $f->stop);
			$row->value($strand_i, $f->strand);
		});
	}
	
	# collect attribute tags, this does not recurse
	if (@include_tags) {
		foreach my $t (@include_tags) {
			my $i = $Data->add_column($t);
			$Data->iterate( sub {
				my $row = shift;
				# get the tag value from the feature and record it, null if not present
				my @v = $row->seqfeature->get_tag_values($t); # could be more than 1
				$row->value($i, @v ? join(',', @v) : '.');
			});
		}
	}
	
	# sort if requested
	if ($sort_data and $include_coordinates) {
		print " Sorting data...\n";
		$Data->gsort_data;
	}
	
	# write the output
	my $success = $Data->write_file(
		filename => $outfile,
		gz       => $gz,
	);
	printf " wrote file $success\n" if $success;
}


sub export_to_ucsc {
	# check output filename
	unless ($outfile =~ /\.(?:refflat|ucsc)(?:\.gz)?$/i) {
		$outfile .= '.refFlat';
	}
	
	# how we write the output depends on whether we need to sort the file or not
	if ($sort_data) {
		# we need to keep all gff lines in memory so that we can sort the lines
		
		my $outData = Bio::ToolBox::Data->new(
			ucsc => 11,
		) or die "unable to create output Data structure!\n";
		$outData->add_comment( sprintf("exported from %s\n", 
			$database ? $database : $input) );
		
		# iterate
		$Data->iterate( sub {
			my $row = shift;
			my $string = ucsc_string( $row->seqfeature(1) );
			foreach (split /\n/, $string) {
				$outData->add_row($_);
			}
			$Data->delete_seqfeature($row->row_index);
		});
		print " Sorting data...\n";
		$outData->gsort_data;
		$outfile = $outData->write_file(
			filename => $outfile, 
			gz       => $gz,
		);
	}
	else {
		# we can simply write out gff directly
		$outfile .= '.gz' if ($gz and $outfile !~ /\.gz$/i);
		my $fh = Bio::ToolBox::Data->open_to_write_fh($outfile, $gz) or 
			die "unable to open '$outfile' for writing! $!\n";
		$fh->printf("# exported from %s\n", $database ? $database : $input);
	
		# write to UCSC file
		$Data->iterate( sub {
			my $row = shift;
			my $string = ucsc_string( $row->seqfeature(1) );
			$fh->print($string);
		});
		$fh->close;
	}
	printf " wrote file $outfile\n";
}

sub adjust_coordinates {
	my $feature = shift;
	
	# we will always adjust relative coordinates based on strand
	# and not absolute coordinates
	my ($start, $end);
	
	# get the original coordinates
	my $fstart = $feature->start;
	my $fend   = $feature->end;
	
	# adjust from 5' end
	if ($position eq '5') {
	
		if ($feature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$start = $fstart + $start_adj;
			}
			else {
				$start = $fstart;
			}
			if ($stop_adj) {
				$end = $fstart + $stop_adj;
			}
			else {
				$end = $fstart;
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$end = $fend - $start_adj;
			}
			else {
				$end = $fend;
			}
			if ($stop_adj) {
				$start = $fend - $stop_adj;
			}
			else {
				$start = $fend;
			}
		}
	}
	
	# adjust from 3' end
	elsif ($position eq '3') {
	
		if ($feature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$end = $fend + $start_adj;
			}
			else {
				$end = $fend;
			}
			if ($stop_adj) {
				$start = $fend + $stop_adj;
			}
			else {
				$start = $fend;
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$start = $fstart - $start_adj;
			}
			else {
				$start = $fstart;
			}
			if ($stop_adj) {
				$end = $fstart - $stop_adj;
			}
			else {
				$end = $fstart;
			}
		}
	}
	
	# adjust from middle position
	elsif ($position eq 'm' or $position eq '4') {
		
		my $midpoint = int( ( ($fstart + $fend) / 2) + 0.5);
		if ($feature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$start = $midpoint + $start_adj;
			}
			else {
				$start = $midpoint;
			}
			if ($stop_adj) {
				$end = $midpoint + $stop_adj;
			}
			else {
				$end = $midpoint;
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$end = $midpoint - $start_adj;
			}
			else {
				$end = $midpoint;
			}
			if ($stop_adj) {
				$start = $midpoint - $stop_adj;
			}
			else {
				$start = $midpoint;
			}
		}
	}
	
	# adjust from both ends
	elsif ($position eq '53') {
	
		if ($feature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$start = $fstart + $start_adj;
			}
			else {
				$start = $fstart;
			}
			if ($stop_adj) {
				$end = $fend + $stop_adj;
			}
			else {
				$end = $fend;
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$end = $fend - $start_adj;
			}
			else {
				$end = $fend;
			}
			if ($stop_adj) {
				$start = $fstart - $stop_adj;
			}
			else {
				$start = $fstart;
			}
		}
	}
	
	# flip coordinates to make start and stop consistent with strand
	# sometimes when only one coordinate is changed, it flips the orientation
	# start must always be less than the stop coordinate
	# but always respect the given strand
	if ($start > $end) {
		my $newstart = $end;
		my $newend   = $start;
		$start = $newstart;
		$end = $newend;
	}
	
	# make sure no negative coordinates sneak through
	$start = 1 if $start < 1;
	
	# return
	return ($start, $end);
}


__END__

=head1 NAME

get_features.pl

A program to collect and filter annotated features from source files.

=head1 SYNOPSIS

get_features.pl --in E<lt>filenameE<gt> --out E<lt>filenameE<gt>

get_features.pl --db E<lt>nameE<gt> --out E<lt>filenameE<gt>
  
  Source data:
  -d --db <name | filename>     database: name, file.db, or file.sqlite
  -i --in <filename>            input annotation: GFF3, GTF, genePred, etc
  
  Selection:
  -f --feature <type>           feature: gene, mRNA, transcript, etc
  -u --sub                      include subfeatures (true if gff, gtf, refFlat)
  
  Filter features:
  -K --chrskip <regex>          skip features from certain chromosomes
  -x --exclude <tag=value>      exclude features with specific attribute value
  --biotype <regex>             include only specific biotype
  --tsl [best|best1|best2|      specify minimum transcript support level 
         best3|best4|best5|       primarily Ensembl annotation 
         1|2|3|4|5|NA]  
  --gencode                     include only GENCODE tagged genes
  
  Adjustments:
  -b --start=<integer>          modify start positions
  -e --stop=<integer>           modify stop positions
  -p --pos [ 5 | m | 3 | 53 ]   relative position from which to modify
  --collapse                    collapse subfeatures from alt transcripts
  
  Report format options:
  -B --bed                      write BED6 (no --sub) or BED12 (--sub) format
  -G --gff                      write GFF3 format
  -g --gtf                      write GTF format
  -r --refflat                  write UCSC refFlat format
  -t --tag <text>               include specific GFF attributes in text output
  --coord                       include coordinates in text output
  
  General options:
  -o --out <filename>           output file name
  --sort                        sort output by genomic coordinates
  -z --gz                       compress output
  -Z --bgz                      bgzip compress output
  -v --version                  print version and exit
  -h --help                     show full documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 Source data

=over 4

=item --db E<lt>textE<gt>

Specify the name of a L<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A SQLite file 
or a named relational database may be provided. Used as an alternate 
to an input file.

=item --in E<lt>filenameE<gt>

Specify the filename of a gene annotation file, including GTF, GFF3, 
or a UCSC-formatted file including, refFlat, genePred, or knownGene.
The file may be gzip compressed. Used as an alternate to a database.

=back

=head2 Selection

=over 4

=item --feature E<lt>typeE<gt>

Provide a feature type to collect. Typically, this corresponds to the 
GFF type, or 3rd column in a GFF3 or GTF file. Examples include 
C<gene>, C<mRNA>, or C<transcript>. The default value for input files 
is 'C<gene>'. For databases, an interactive list will be presented 
from which one or more may be chosen.

=item --sub

Optionally include all child subfeatures in the output. For example, 
transcript, CDS, and/or exon subfeatures of a gene. This option is 
automatically enabled with GFF, GTF, or refFlat output; it may be 
turned off with C<--nosub>. With BED output, it will force a BED12 
file to be written. It has no effect with standard text. 

=back

=head2 Filter features

=over 4

=item --chrskip E<lt>regexE<gt>

Exclude features from the output whose sequence ID or chromosome matches 
the provided regex-compatible string. Expressions should be quoted or 
properly escaped on the command line. Examples might be 
    
    'chrM'
    'scaffold.+'
    'chr.+alt|chrUn.+|chr.+_random'

=item --exclude E<lt>tag=valueE<gt>

Provide a GFF/GTF attribute tag on which to filter out the features 
matching the indicated value. For example, to filter on protein 
coding genes using C<gene_biotype>, specify "gene_biotype=protein_coding". 
This filter does not descend into subfeatures. More than one exclusion 
tag may be specified with multiple options or a comma-delimited list. 

=item --biotype E<lt>regex<gt> 

Filter transcripts using the C<transcript_biotype> or C<biotype> 
GTF/GFF3 attribute, typically found in Ensembl annotation files. Provide 
a regex compatible string which must match the biotype value to keep the 
transcripts. For example, to keep specify "miRNA" to keep all micro-RNA 
transcripts. This works on a subfeature level as well, so that C<gene> 
may be specified as the feature to collect, and only the gene transcripts 
belonging to the indicating biotype are retained.

=item --tsl E<lt>levelE<gt>

Filter transcripts on the Ensembl GTF/GFF3 attribute C<transcript_support_level>, 
which is described at L<Ensembl TSL glossary entry|http://uswest.ensembl.org/info/website/glossary.html>.
Provide a level of support to filter. Values include: 
    
    1       All splice junctions supported by evidence
    2       Transcript flagged as suspect or only support from multiple ESTs
    3       Only support from single EST
    4       Best supporting EST is suspect
    5       No support
    best    Transcripts at the best (lowest) available level are taken
    best1   The word followed by a digit 1-5, indicating any transcript 
            at or better (lower) than the indicated level
    NA      Only transcripts without a level (NA) are retained.

=item --gencode

Boolean option to filter transcripts as part of the GENCODE specification. 
These are marked in Ensembl GTF/GFF3 annotation files as the C<tag> attribute 
with value "basic". Typically, at least one transcript for every gene is 
marked as part of the GENCODE set. Transcripts not marked as such usually 
lack sufficient experimental evidence.

=back

=head2 Adjustments

=over 4

=item --start=<integer>

=item --stop=<integer>

Optionally specify adjustment values to adjust the reported start and 
end coordinates of the collected regions. A negative value is shifted 
upstream (towards the 5 prime end), and a positive value is shifted 
downstream (towards the 3 prime end). Adjustments are made relative 
to the indicated position (--pos option, below) based on the feature 
strand. Adjustments are only allowed when writing standard BED6 or 
standard text files.

=item --pos [ 5 | m | 3 | 53 ]

Indicate the relative position from which both coordinate adjustments 
are made. Both start and stop adjustments may be made from the respective 
5 prime, 3 prime, or middle position as dictated by the feature's strand 
value. Alternatively, specify '53' to indicate that the start adjustment 
adjusts the 5 prime end and the stop adjustment adjusts the 3 prime end. 
The default is '53'.

=item --collapse

Boolean option to collapse multiple alternate transcripts of a gene into 
a single artificial transcript, merging overlapping exons and minimizing 
introns, where appropriate. Genes without alternate transcripts are not 
collapsed.

=back

=head2 Report format options

=over 4

=item --bed

With subfeatures enabled, write a BED12 (12-column BED) file. 
Otherwise, write a standard 6-column BED format file. 

=item --gff

Write a GFF version 3 (GFF3) format output file. Subfeatures are 
automatically included and coordinate adjustments ignored.

=item --gtf

Write a GTF format (GFF version 2.2 or 2.5) output file. Subfeatures are 
automatically included and coordinate adjustments ignored.

=item --refflat 

=item --ucsc

Write a UCSC-style refFlat format (10 columns) gene annotation table. 
Subfeatures are automatically included and coordinate adjustments ignored.

=item --tag E<lt>textE<gt>

When writing a standard text file, optionally include additional 
GFF/GTF attribute tags. Specify as a comma-delimited list or as 
separate options.

=item --coord

When writing a standard text file, optionally include the chromosome, 
start, stop, and strand coordinates. These are automatically included 
in other formats. This is automatically included when adjusting 
coordinate positions.

=back

=head2 General options

=over 4

=item --out E<lt>filenameE<gt>

Specify the output file name. Default is the joined list of features. 

=item --sort

Sort the output file by genomic coordinates. Automatically enabled 
when compressing with bgzip. This may require more memory.

=item --gz

Specify whether the output file should be compressed with gzip.

=item --bgz

Specify whether the output file should be compressed with block gzip 
(bgzip) for tabix compatibility.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will extract a list of features from a database or input 
annotation file and write them out to a file. Features may be selected 
using their feature type (the 3rd column in a GFF or GTF file). When 
selecting features from a database, types may be selected interactively 
from a presented list. Features may be filtered based on various 
GFF attributes typically found in Ensembl annotation, including 
C<transcript_biotype>, C<transcript_support_level>, and GENCODE basic 
tags. Features may also be filtered by chromosome. 

Collected features may be written to a variety of formats, including 
GFF3, GTF, refFlat, simple 6-column BED, or a simple text format. With 
GFF, GTF, and refFlat formats, subfeatures such as exons are automatically 
included (which may also be turned off). With a simple text format, 
the source database or parsed input file are recorded in the header 
metadata for use in subsequent programs. Coordinates may be optionally 
included in the text file, which preempts using parsed features in other 
tools. 

=head2 Coordinate adjustments

Coordinates of the features may be adjusted as desired when writing to 
text or BED file formats. Adjustments may be made relative to either 
the 5 prime, 3 prime, both ends, or the feature midpoint. Positions 
are based on the feature strand. Use the following examples as a guide. 

=over 4

=item upstream 500 bp only

  get_features.pl --start=-500 --stop=-1 --pos 5

=item 1 kb total around 5 prime end

  get_features.pl --start=-500 --stop=500 --pos 5

=item last 500 bp of feature

  get_features.pl --start=-500 --pos 3

=item middle 500 bp of feature

  get_features.pl --start=-250 --stop=250 --pos m

=item entire feature plus 1 kb of flanking

  get_features.pl --start=-1000 --stop=1000 --pos 53

=back

Note that positions are always in base coordinates, and the resulting regions 
may be 1 bp longer depending on whether the reference base was included or not.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
