#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(verify_or_request_feature_types);
use Bio::ToolBox::GeneTools qw(
	:export 
	collapse_transcripts
	filter_transcript_support_level
	filter_transcript_gencode_basic
	filter_transcript_biotype
);
use Bio::ToolBox::utility;
my $VERSION = '1.54';

print "\n This program will collect features from a database\n\n";

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
	$gz,
	$help,
	$print_version,
);
my @features;
my @include_tags;
my @exclude_tags;
my %exclude_tag2value;

# Command line options
GetOptions( 
	'in=s'      => \$input, # input table
	'db=s'      => \$database, # source annotation database
	'feature=s' => \@features, # the features to collect from the database
	'sub!'      => \$get_subfeatures, # collect subfeatures
	'coord!'    => \$include_coordinates, # collect coordinates
	'start=i'   => \$start_adj, # start coordinate adjustment
	'stop=i'    => \$stop_adj, # stop coordinate adjustment
	'pos=s'     => \$position, # relative position to adjust coordinates
	'tag=s'     => \@include_tags, # attributes to include
	'exclude=s' => \@exclude_tags, # attribute and keys to exclude
	'tsl=s'     => \$tsl, # filter on transcript support level
	'gencode!'  => \$gencode, # filter on gencode basic tag
	'biotype=s' => \$tbiotype, # filter on transcript biotype
	'collapse!' => \$collapse, # collapse multi-transcript genes
	'chrskip=s' => \$chromosome_exclude, # skip chromosomes
	'bed!'      => \$convert_to_bed, # convert to bed format
	'gff|gff3!' => \$convert_to_gff, # convert to GFF3 format
	'gtf!'      => \$convert_to_gtf, # convert to gtf format
	'ucsc|refFlat!' => \$convert_to_refflat, # convert to refFlat format
	'out=s'     => \$outfile, # name of output file 
	'gz!'       => \$gz, # compress output
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
	print " Biotoolbox script <NAME>.pl, version $VERSION\n\n";
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
		unless ($convert_to_gff or $convert_to_gtf or $convert_to_refflat) {
			die " Cannot collapse transcripts unless writing to GFF/GTF/refFlat!\n";
		}
	}
	
	# check adjustments
	if ($start_adj or $stop_adj) {
		# automatically include coordinates if we're adjusting them
		if ($convert_to_gff or $convert_to_gtf or $convert_to_refflat) {
			die " Cannot adjust coordinates when converting to GFF/GTF/refFlat!\n";
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
}


sub load_from_database {
	# validate and/or request features
	@features = verify_or_request_feature_types(
		'db'      => $database,
		'feature' => [ @features ],
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
		printf " Loaded %s features from $input.\n", format_with_commas($D->last_row);
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
		subfeature => $get_subfeatures ? 'exon,cds,utr' : '',
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
				my ($v) = $row->seqfeature->get_tag_values($k);
				if ($v =~ /$check/i) {
					# the tag matches, so discard
					push @unwanted, $row->row_index;
				}
			});
			if (@unwanted) {
				$Data->delete_row(@unwanted);
			}
		}
	}
	
	# filter on gencode
	if ($tbiotype) {
		print " Filtering transcript biotype for $tbiotype...\n";
		my @unwanted;
		$Data->iterate( sub {
			my $row = shift;
			my $good = filter_transcript_biotype($row->seqfeature, $tbiotype);
			if ($good) {
				$Data->store_seqfeature($row->row_index, $good);
			}
			else {
				push @unwanted, $row->row_index;
			}
		});
		if (@unwanted) {
			$Data->delete_row(@unwanted);
		}
	}
	
	# filter on tsl
	if ($tsl) {
		print " Filtering for transcript support level of $tsl...\n";
		my @unwanted;
		$Data->iterate( sub {
			my $row = shift;
			my $good = filter_transcript_support_level($row->seqfeature, $tsl);
			if ($good) {
				$Data->store_seqfeature($row->row_index, $good);
			}
			else {
				push @unwanted, $row->row_index;
			}
		});
		if (@unwanted) {
			$Data->delete_row(@unwanted);
		}
	}
	
	# filter on gencode
	if ($gencode) {
		print " Filtering for GENCODE transcripts...\n";
		my @unwanted;
		$Data->iterate( sub {
			my $row = shift;
			my $good = filter_transcript_support_level($row->seqfeature, $tsl);
			if ($good) {
				$Data->store_seqfeature($row->row_index, $good);
			}
			else {
				push @unwanted, $row->row_index;
			}
		});
		if (@unwanted) {
			$Data->delete_row(@unwanted);
		}
	}
	
	# collapse transcripts
	if ($collapse) {
		print " Collapsing alternate transcripts...\n";
		$Data->iterate( sub {
			my $row = shift;
			my $gene = collapse_transcripts($row->seqfeature);
			if ($gene) {
				$Data->store_seqfeature($row->row_index, $gene);
			}
		});
	}
}


sub export_to_bed {
	# prepare output
	unless ($outfile =~ /\.bed(?:\.gz)?$/i) {
		$outfile .= '.bed';
	}
	$outfile .= '.gz' if ($gz and $outfile !~ /\.gz$/i);
	my $fh = Bio::ToolBox::Data->open_to_write_fh($outfile, $gz) or 
		die "unable to open '$outfile' for writing! $!\n";
	
	# adjust coordinates as necessary and write
	if ($start_adj or $stop_adj) {
		$Data->iterate( sub {
			my $row = shift;
			my ($start, $stop) = adjust_coordinates($row); 
			my $string = $row->bed_string(
				start => $start,
				end   => $stop,
			);
			$fh->print( $string . "\n" );
		});
	}
	else {
		$Data->iterate( sub {
			my $row = shift;
			$fh->print( $row->bed_string . "\n" );
		});
	}
	$fh->close;
	printf " wrote file $outfile\n";
}


sub export_to_gff {
	# prepare output
	unless ($outfile =~ /\.gff3?(?:\.gz)?$/i) {
		$outfile .= '.gff3';
	}
	$outfile .= '.gz' if ($gz and $outfile !~ /\.gz$/i);
	my $fh = Bio::ToolBox::Data->open_to_write_fh($outfile, $gz) or 
		die "unable to open '$outfile' for writing! $!\n";
	$fh->print("##gff-version 3\n");
	$fh->printf("# exported from %s\n", $database ? $database : $input);
	
	# write to gff3
	$Data->iterate( sub {
		my $row = shift;
		my $string = $row->seqfeature->gff_string(1);
		$fh->print( $string . "###\n"); # should already have a newline
	});
	
	$fh->close;
	printf " wrote file $outfile\n";
}


sub export_to_gtf {
	# prepare output
	unless ($outfile =~ /\.gtf(?:\.gz)?$/i) {
		$outfile .= '.gtf';
	}
	$outfile .= '.gz' if ($gz and $outfile !~ /\.gz$/i);
	my $fh = Bio::ToolBox::Data->open_to_write_fh($outfile, $gz) or 
		die "unable to open '$outfile' for writing! $!\n";
	$fh->print("##gff-version 2.5\n");
	$fh->printf("# exported from %s\n", $database ? $database : $input);
	
	# write to GTF
	$Data->iterate( sub {
		my $row = shift;
		my $string = gtf_string( $row->seqfeature );
		$fh->print( $string);
	});
	
	$fh->close;
	printf " wrote file $outfile\n";
}


sub export_to_txt {
	# adjust coordinates as necessary
	if ($start_adj or $stop_adj or $include_coordinates) {
		
		# add coordinate columns
		my $seq_i    = $Data->add_column('Chromosome');
		my $start_i  = $Data->add_column('Start');
		my $stop_i   = $Data->add_column('Stop');
		my $strand_i = $Data->add_column('Strand');
		
		# iterate
		$Data->iterate( sub {
			my $row = shift;
			my ($start, $stop) = adjust_coordinates($row); 
			$row->value($seq_i, $row->seq_id);
			$row->value($start_i, $start);
			$row->value($stop_i, $stop);
			$row->value($strand_i, $row->strand);
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
	
	# write the output
	my $success = $Data->write_file(
		filename => $outfile,
		gz       => $gz,
	);
	printf " wrote file $success\n" if $success;
}


sub export_to_ucsc {
	# prepare output
	unless ($outfile =~ /\.(?:refflat|ucsc)(?:\.gz)?$/i) {
		$outfile .= '.refFlat';
	}
	$outfile .= '.gz' if ($gz and $outfile !~ /\.gz$/i);
	my $fh = Bio::ToolBox::Data->open_to_write_fh($outfile, $gz) or 
		die "unable to open '$outfile' for writing! $!\n";
	$fh->printf("# exported from %s\n", $database ? $database : $input);
	
	# write to gff3
	$Data->iterate( sub {
		my $row = shift;
		my $string = ucsc_string( $row->seqfeature );
		$fh->print($string);
	});
	
	$fh->close;
	printf " wrote file $outfile\n";
}

sub adjust_coordinates {
	my $feature = shift;
	
	# we will always adjust relative coordinates based on strand
	# and not absolute coordinates
	
	# get the original coordinates
	my $start = $feature->start;
	my $end   = $feature->end;
	
	# adjust from 5' end
	if ($position eq '5') {
	
		if ($feature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$start = $start + $start_adj;
			}
			if ($stop_adj) {
				$end = $start + $stop_adj;
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$end = $end - $start_adj;
			}
			if ($stop_adj) {
				$start = $end - $stop_adj;
			}
		}
	}
	
	# adjust from 3' end
	elsif ($position eq '3') {
	
		if ($feature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$end = $end + $start_adj;
			}
			if ($stop_adj) {
				$start = $end + $stop_adj;
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$start = $start - $start_adj;
			}
			if ($stop_adj) {
				$end = $start - $stop_adj;
			}
		}
	}
	
	# adjust from middle position
	elsif ($position eq 'm' or $position eq '4') {
		
		my $midpoint = int( ( ($start + $end) / 2) + 0.5);
		if ($feature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$start = $midpoint + $start_adj;
			}
			if ($stop_adj) {
				$end = $midpoint + $stop_adj;
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$end = $midpoint - $start_adj;
			}
			if ($stop_adj) {
				$start = $midpoint - $stop_adj;
			}
		}
	}
	
	# adjust from both ends
	elsif ($position eq '53') {
	
		if ($feature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$start = $start + $start_adj;
			}
			if ($stop_adj) {
				$end = $end + $stop_adj;
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$end = $end - $start_adj;
			}
			if ($stop_adj) {
				$start = $start - $stop_adj;
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

A script to collect and filter annotated features from source files.

=head1 SYNOPSIS

get_features.pl --in E<lt>filenameE<gt> --out E<lt>filenameE<gt>

get_features.pl --db E<lt>nameE<gt> --out E<lt>filenameE<gt>
  
  Source data:
  --db <name | filename>    (name, file.db, or file.sqlite)
  --in <filename>           (gff, gtf, genePred, etc)
  
  Selection:
  --feature <type>          (gene, mRNA, transcript, etc)
  --sub
  
  Filter features:
  --exclude <tag=value>
  --biotype <regex>
  --tsl [best|best1|best2|best3|best4|best5|1|2|3|4|5|NA]
  --gencode
  --chrskip <regex>
  
  Adjustments:
  --collapse
  --start=<integer>
  --stop=<integer>
  --pos [ 5 | m | 3 | 53 ]
  
  Report format options:
  --coord
  --tag <text>
  --bed
  --gff or --gff3
  --gtf
  --refflat or --ucsc
  
  General options:
  --out <filename>
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db E<lt>textE<gt>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A SQLite file 
or a named relational database may be provided. Used as an alternate 
to an input file.

=item --in E<lt>filenameE<gt>

Specify the filename of a gene annotation file, including GTF, GFF3, 
or a UCSC-formatted file including, refFlat, genePred, or knownGene.
The file may be gzip compressed. Used as an alternate to a database.

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
turned off with "--nosub". It has no effect with standard text or BED output.

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

=item --chrskip E<lt>regexE<gt>

Exclude features from the output whose sequence ID or chromosome matches 
the provided regex-compatible string. Expressions should be quoted or 
properly escaped on the command line. Examples might be 
    
    'chrM'
    'scaffold.+'
    'chr.+alt|chrUn.+|chr.+_random'

=item --collapse

Boolean option to collapse multiple alternate transcripts of a gene into 
a single artificial transcript, merging overlapping exons and minimizing 
introns, where appropriate. Genes without alternate transcripts are not 
collapsed.

=item --start=<integer>

=item --stop=<integer>

Optionally specify adjustment values to adjust the reported start and 
end coordinates of the collected regions. A negative value is shifted 
upstream (towards the 5 prime end), and a positive value is shifted 
downstream (towards the 3 prime end). Adjustments are made relative 
to the indicated position (--pos option, below) based on the feature 
strand. Adjustments are only allowed when writing BED or standard 
text files.

=item --pos [ 5 | m | 3 | 53 ]

Indicate the relative position from which both coordinate adjustments 
are made. Both start and stop adjustments may be made from the respective 
5 prime, 3 prime, or middle position as dictated by the feature's strand 
value. Alternatively, specify '53' to indicate that the start adjustment 
adjusts the 5 prime end and the stop adjustment adjusts the 3 prime end. 
The default is '53'.

=item --coord

When writing a standard text file, optionally include the chromosome, 
start, stop, and strand coordinates. These are automatically included 
in other formats. This is automatically included when adjusting 
coordinate positions.

=item --tag E<lt>textE<gt>

When writing a standard text file, optionally include additional 
GFF/GTF attribute tags. Specify as a comma-delimited list or as 
separate options.

=item --bed

Write a standard 6-column BED format file. Subfeatures are not included.

=item --gff

Write a GFF version 3 (GFF3) format output file. Subfeatures are 
automatically included and coordinate adjustments ignored.

=item --gtf

Write a GTF format (GFF version 2.2 or 2.5) output file. Subfeatures are 
automatically included and coordinate adjustments ignored.

=item --refflat or --ucsc

Write a UCSC-style refFlat format (10 columns) gene annotation table. 
Subfeatures are automatically included and coordinate adjustments ignored.

=item --out <filename>

Specify the output file name. Default is the joined list of features. 

=item --gz

Specify whether the output file should be compressed with gzip.

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

=head1 COORDINATE ADJUSTMENTS

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
