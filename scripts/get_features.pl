#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	validate_included_feature
);
use Bio::ToolBox::file_helper qw(
	open_to_write_fh
);
use Bio::ToolBox::db_helper::config;
use Bio::ToolBox::utility;
my $VERSION = '1.20';

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
	$database,
	$all,
	$get_subfeatures,
	$include_coordinates,
	$start_adj,
	$stop_adj,
	$position,
	$convert_to_bed,
	$convert_to_gff,
	$outfile,
	$gz,
	$help,
	$print_version,
);
my @features;

# Command line options
GetOptions( 
	'db=s'      => \$database, # source annotation database
	'feature=s' => \@features, # the features to collect from the database
	'all!'      => \$all, # all features must be collected
	'sub!'      => \$get_subfeatures, # collect subfeatures
	'coord!'    => \$include_coordinates, # collect coordinates
	'start=i'   => \$start_adj, # start coordinate adjustment
	'stop=i'    => \$stop_adj, # stop coordinate adjustment
	'pos=s'     => \$position, # relative position to adjust coordinates
	'bed!'      => \$convert_to_bed, # convert to bed format
	'gff!'      => \$convert_to_gff, # convert to GFF3 format
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
unless ($database) {
	die " must provide a database name! use --help for more information\n";
}
unless (defined $gz) {
	$gz = 0;
}
unless (defined $convert_to_bed) {
	$convert_to_bed = 0;
}
unless (defined $convert_to_gff) {
	$convert_to_gff = 0;
}
unless (defined $get_subfeatures) {
	$get_subfeatures = 0;
}
if ($convert_to_bed and $convert_to_gff) {
	die " cannot convert to both GFF and BED formats!\n";
}
if ($start_adj or $stop_adj) {
	# automatically include coordinates if we're adjusting them
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


### Open database
my $db = open_db_connection($database) or 
	die " unable to open database connection to '$database'!\n";
# check database
{
	my $db_ref = ref $db;
	unless ($db_ref =~ m/Bio::DB::SeqFeature::Store/) {
		die " $db_ref databases are not supported! sorry\n";
	}
}


### Get features
# check if it is a comma delimited list
if (scalar @features == 1 and $features[0] =~ /,/) {
	@features = split /,/, shift @features;
}

# validate and/or request features
@features = verify_or_request_feature_types(
	'db'      => $db,
	'feature' => [ @features ],
	'prompt'  => " Enter the feature(s) to collect." . 
			" A comma de-limited list or range may be given\n",
) or die " no valid features were provided! see help\n";



### Prepare data structure and/or output file
my ($Data, $out_fh) = prepare_data_structure_or_output();



### Collect features
my $count;
if ($convert_to_gff) {
	print " Collecting features as GFF...\n";
	$count = collect_features(\&record_gff_feature);
}
elsif ($convert_to_bed) {
	print " Collecting features as BED...\n";
	$count = collect_features(\&record_bed_feature);
}
else {
	print " Collecting features...\n";
	if ($include_coordinates) {
		$count = collect_features(\&record_standard_coordinate_feature);
	}
	else {
		$count = collect_features(\&record_standard_feature);
	}
}
print "  Collected " . format_with_commas($count) . " features\n";



### Print output
if ($convert_to_gff) {
	$out_fh->close;
	print " Wrote file '$outfile'\n";
}
else {
	unless ($outfile) {
		$outfile = generate_file_name();
	}
	my $success = $Data->write_file(
		'filename' => $outfile,
		'gz'       => $gz,
	);
	if ($success) {
		print " Wrote file '$success'\n";
	}
	else {
		print " Failed to write file!\n";
	}
}




########################   Subroutines   ###################################

sub prepare_data_structure_or_output {
	# how we prepare the structure is dependent on the output format
	my ($Data, $fh);
	
	if ($convert_to_gff) {
		# no need for a data structure
		# we will be printing directly to outfile
		
		# check filename
		unless ($outfile) {
			$outfile = generate_file_name();
		}
		unless ($outfile =~ /\.gff3?(?:\.gz)?$/i) {
			$outfile .= '.gff3';
		}
		
		# open file handle
		$fh = open_to_write_fh($outfile, $gz) or 
			die " unable to open output file for writing!\n";
		
		# print GFF headers
		$fh->print("##gff-version 3\n");
		$fh->print("# Features collected from database $database\n");
	}
	else {
		# generate a tim data structure
		
		# structure dependent on output format
		if ($convert_to_bed) {
			# bed structure
			$Data = Bio::ToolBox::Data->new(
				feature   => 'region', 
				datasets  => [qw(Chromosome Start End Name Score Strand)],
			);
			$Data->bed = 6; # set the bed parameter
			
			# position adjustments
			if ($start_adj) {
				$Data->metadata(1, 'start_adjustment', $start_adj);
				$Data->metadata(1, 'position', $position);
			}
			if ($stop_adj) {
				$Data->metadata(2, 'start_adjustment', $start_adj);
				$Data->metadata(2, 'position', $position);
			}
		}
		
		elsif ($include_coordinates) {
			# generic tim data structure with coordinates
			
			$Data = Bio::ToolBox::Data->new(
				feature   => 'region', 
				datasets  => [qw(Name Type Chromosome Start End Strand)],
			);
		
			# position adjustments
			if ($start_adj) {
				$Data->metadata(4, 'start_adjustment', $start_adj);
				$Data->metadata(4, 'position', $position);
			}
			if ($stop_adj) {
				$Data->metadata(5, 'start_adjustment', $start_adj);
				$Data->metadata(5, 'position', $position);
			}
		}
			
		else {
			# just name and type, thank you very much
			$Data = Bio::ToolBox::Data->new(
				feature   => join(',', @features), 
				datasets  => [qw(Primary_ID Name Type)],
			);
		}
		
		# add database
		$Data->database($database);
	}
	return ($Data, $fh);
}



sub collect_features {
	
	# collection method
	my $method = shift;
	
	# Get the names of chromosomes to avoid
	my @excluded_chromosomes = 
		$BTB_CONFIG->param("$database\.chromosome_exclude");
	unless (@excluded_chromosomes) {
		@excluded_chromosomes = 
			$BTB_CONFIG->param('default_db.chromosome_exclude');
	}
	my %excluded_chr_lookup = map {$_ => 1} @excluded_chromosomes;
	
	# generate a seqfeature stream
	my $iterator = $db->features(
		-type     => \@features,
		-iterator => 1,
	);
	
	# process the features
	my $count = 0;
	while (my $seqfeat = $iterator->next_seq) {
		
		# skip it unless it is validated
		unless ($all) {
			# pass the validation attributes
			next unless validate_included_feature($seqfeat);
			
			# not on excluded chromosome
			next if (exists $excluded_chr_lookup{ $seqfeat->seq_id });
		}
		
		# record each feature according to the method
		&{$method}($seqfeat);
		$count++;
	}
	
	return $count;	
}



sub record_gff_feature {
	my $seqfeature = shift;
	
	# directly print the feature
	# do or do not include subfeatures
	$out_fh->print( $seqfeature->gff3_string($get_subfeatures) . "\n" );
}



sub record_bed_feature {
	my $seqfeature = shift;
	
	# adjust coordinates if requested
	if ($start_adj or $stop_adj) {
		adjust_coordinates($seqfeature);
	}
	
	# record the feature
	my @fdata = (
		$seqfeature->seq_id,
		$seqfeature->start - 1,
		$seqfeature->end,
		$seqfeature->display_name || $seqfeature->primary_id || 'region',
		$seqfeature->score || 0,
		$seqfeature->strand, # this will be converted to + or - upon writing
	);
	$Data->add_row(\@fdata);
	
	# record subfeatures if requested
	if ($get_subfeatures) {
		foreach my $subfeat ($seqfeature->get_SeqFeatures) {
			record_bed_feature($subfeat);
		}
	}
}



sub record_standard_coordinate_feature {
	my $seqfeature = shift;
	
	# adjust coordinates if requested
	if ($start_adj or $stop_adj) {
		adjust_coordinates($seqfeature);
	}
		
	# record the feature
	my @fdata = (
		$seqfeature->display_name,
		$seqfeature->type,
		$seqfeature->seq_id,
		$seqfeature->start,
		$seqfeature->end,
		$seqfeature->strand,
	);
	$Data->add_row(\@fdata);
	
	# record subfeatures if requested
	if ($get_subfeatures) {
		foreach my $subfeat ($seqfeature->get_SeqFeatures) {
			record_standard_coordinate_feature($subfeat);
		}
	}
}



sub record_standard_feature {
	my $seqfeature = shift;
	
	# record the feature
	my @fdata = (
		$seqfeature->primary_id,
		$seqfeature->display_name,
		$seqfeature->type,
	);
	$Data->add_row(\@fdata);
	
	# record subfeatures if requested
	if ($get_subfeatures) {
		foreach my $subfeat ($seqfeature->get_SeqFeatures) {
			record_standard_feature($subfeat);
		}
	}
}



sub adjust_coordinates {
	my $seqfeature = shift;
	
	# we will always adjust relative coordinates based on strand
	# and not absolute coordinates
	
	# get the original coordinates
	my $start = $seqfeature->start;
	my $end   = $seqfeature->end;
	
	# adjust from 5' end
	if ($position eq '5') {
	
		if ($seqfeature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$seqfeature->start($start + $start_adj);
			}
			if ($stop_adj) {
				$seqfeature->end($start + $stop_adj);
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$seqfeature->end($end - $start_adj);
			}
			if ($stop_adj) {
				$seqfeature->start($end - $stop_adj);
			}
		}
	}
	
	# adjust from 3' end
	elsif ($position eq '3') {
	
		if ($seqfeature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$seqfeature->end($end + $start_adj);
			}
			if ($stop_adj) {
				$seqfeature->start($end + $stop_adj);
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$seqfeature->start($start - $start_adj);
			}
			if ($stop_adj) {
				$seqfeature->end($start - $stop_adj);
			}
		}
	}
	
	# adjust from middle position
	elsif ($position eq 'm' or $position eq '4') {
		
		my $midpoint = int( ( ($start + $end) / 2) + 0.5);
		if ($seqfeature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$seqfeature->start($midpoint + $start_adj);
			}
			if ($stop_adj) {
				$seqfeature->end($midpoint + $stop_adj);
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$seqfeature->end($midpoint - $start_adj);
			}
			if ($stop_adj) {
				$seqfeature->start($midpoint - $stop_adj);
			}
		}
	}
	
	# adjust from both ends
	elsif ($position eq '53') {
	
		if ($seqfeature->strand >= 0) {
			# forward strand
			if ($start_adj) {
				$seqfeature->start($start + $start_adj);
			}
			if ($stop_adj) {
				$seqfeature->end($end + $stop_adj);
			}
		}
		else {
			# reverse strand
			if ($start_adj) {
				$seqfeature->end($end - $start_adj);
			}
			if ($stop_adj) {
				$seqfeature->start($start - $stop_adj)
			}
		}
	}
	
	# something else?
	else {
		die "unrecognized position value '$position'! see help\n";
	}
	
	# flip coordinates to make start and stop consistent with strand
	# sometimes when only one coordinate is changed, it flips the orientation
	# start must always be less than the stop coordinate
	# but always respect the given strand
	if ($seqfeature->start > $seqfeature->end) {
		my $start = $seqfeature->end;
		my $end   = $seqfeature->start;
		$seqfeature->start($start);
		$seqfeature->end($end);
	}
}


sub generate_file_name {
	my $filename = join(',', @features);
	$filename =~ s/:/_/g; # 
	if ($convert_to_gff) {
		$filename .= '.gff';
	}
	elsif ($convert_to_bed) {
		$filename .= '.bed';
	}
	else {
		$filename .= '.txt';
	}
	return $filename;
}



__END__

=head1 NAME

get_features.pl

A script to collect features from a BioPerl SeqFeature::Store database.

=head1 SYNOPSIS

get_features.pl --db <text> [--options...]
  
  Options:
  --db <text>
  --feature <type | type:source>
  --all
  --sub
  --coord
  --start=<integer>
  --stop=<integer>
  --pos [ 5 | m | 3 | 53 ]
  --out <filename>
  --bed
  --gff 
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <text>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A database is 
required for generating new data files with features. For more information 
about using annotation databases, 
see L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. 

=item --feature <type | type:source>

Optionally specify the type of the feature(s) to collect. The GFF 
primary_tag or primary_tag:source_tag should be specified. More than 
one feature may be specified at a time, either as a comma de-limited 
list or separate options. If not specified, then an interactive list 
will be presented to the user for selection.

=item --all

Optionally indicate that all features present in the database must 
be included. By default, certain features may be excluded based on 
parameters defined in the BioToolBox configuration file. See below 
for details.

=item --sub

Optionally include all subfeatures in the output. For example, 
transcript, CDS, and/or exon subfeatures of a gene.

=item --coord

When writing a standard format file, optionally include the chromosome, 
start, stop, and strand coordinates. These are automatically included 
when writing a BED or GFF format.

=item --start=<integer>

=item --stop=<integer>

Optionally specify adjustment values to adjust the reported start and 
end coordinates of the collected regions. A negative value is shifted 
upstream (towards the 5 prime end), and a positive value is shifted 
downstream (towards the 3 prime end). Adjustments are made relative 
to the indicated position (--pos option, below) based on the feature 
strand. Adjustments are ignored if a GFF file is written.

=item --pos [ 5 | m | 3 | 53 ]

Indicate the relative position from which both coordinate adjustments 
are made. Both start and stop adjustments may be made from the respective 
5 prime, 3 prime, or middle position as dictated by the feature's strand 
value. Alternatively, specify '53' to indicate that the start adjustment 
adjusts the 5 prime end and the stop adjustment adjusts the 3 prime end. 
The default is '53'.

=item --out <filename>

Specify the output file name. Default is the joined list of features. 

=item --bed

Optionally indicate that a 6-column BED format file should be 
written. Currently, 12-column BED formats with exon information is 
not supported (yet).

=item --gff

Optionally indicate that a GFF3 format file should be written. This 
option enables the database features to be written completely with 
all attributes. Coordinate adjustments are ignored.

=item --gz

Specify whether (or not) the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display this POD documentation.

=back

=head1 DESCRIPTION

This program will extract a list of features from a database and 
write them out to a file. Specifically, the requested features in 
a Bio::DB::SeqFeature::Store schema database are pulled and written 
as either a list of named features with or without coordinate information, 
a BED-style formatted file, or a GFF3-formatted file. The GFF option 
is essentially a database dump, as it enables a low-level option to 
write the features as original GFF3 lines complete with all attributes.

Features may be specified through their GFF type or primary_tag. They
may be specified as a command-line option or selected interactively
from a presented list. They may be restricted through two options
defined in the BioToolBox configuration file, biotoolbox.cfg. These
include a database-specific or default database option,
"chromosome_exclude", which excludes features located on the listed
chromosomes (such as the mitochondrial chromosome), and the
"exclude_tags", which are attribute keys and values to be avoided.
More information may be found in the configuration file itself.

=head1 COORDINATE ADJUSTMENTS

Coordinates of the features may be adjusted as desired. Adjustments 
may be made relative to either the 5 prime, 3 prime, both ends, or the 
feature midpoint. Positions are based on the feature strand. Use the 
following examples as a guide. 

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
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
