#!/usr/bin/perl

# This script will collect features from a database

# Inspired by Magda and all those confused by GFF3 annotation databases

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	format_with_commas
	generate_tim_data_structure
);
use tim_db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	validate_included_feature
);
use tim_file_helper qw(
	open_to_write_fh
	write_tim_data_file
);
use tim_db_helper::config;
my $VERSION = '1.9.1';

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
@features = verify_or_request_feature_types( {
	'db'      => $db,
	'feature' => [ @features ],
	'prompt'  => " Enter the feature(s) to collect." . 
			" A comma de-limited list or range may be given\n",
} ) or die " no valid features were provided! see help\n";



### Prepare data structure and/or output file
my ($data, $out_fh);
prepare_data_structure_or_output();



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
	my $success = write_tim_data_file( {
		'data'     => $data,
		'filename' => $outfile,
		'gz'       => $gz,
	} );
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
		$out_fh = open_to_write_fh($outfile, $gz) or 
			die " unable to open output file for writing!\n";
		
		# print GFF headers
		$out_fh->print("##gff-version 3\n");
		$out_fh->print("# Features collected from database $database\n");
	}
	else {
		# generate a tim data structure
		
		# structure dependent on output format
		if ($convert_to_bed) {
			# bed structure
			$data = generate_tim_data_structure(
				'region', qw(Chromosome Start End Name Score Strand) );
			$data->{'bed'} = 6; # set the bed parameter
			
			# add extra comments
			# normally added to main data hash, but we want these 
			# written to the bed file where they normally are not
			push @{ $data->{'other'} }, 
				"# Collected " . join(',', @features) . 
				" features from database $database\n";
			
			# position adjustments
			if ($start_adj) {
				$data->{1}{'start_adjustment'} = $start_adj;
			}
			if ($stop_adj) {
				$data->{2}{'stop_adjustment'} = $stop_adj;
			}
		}
		else {
			# generic structure
			
			# set the main feature string to be used in the output file
			my $feature_string;
			if ($include_coordinates) {
				# collecting coordinates, which takes precedence 
				# over the named features, especially if adjusting coordinates
				$feature_string = 'region';
			}
			else {
				# named features
				$feature_string = join(',', @features);
			}
			
			# depends on whether we want coordinates or not
			if ($include_coordinates) {
				# don't forget the coordinates
				$data = generate_tim_data_structure(
					$feature_string, 
					qw(Name Type Chromosome Start End Strand)
				);
			
				# position adjustments
				if ($start_adj) {
					$data->{3}{'start_adjustment'} = $start_adj;
				}
				if ($stop_adj) {
					$data->{4}{'stop_adjustment'} = $stop_adj;
				}
			}
			
			else {
				# just name and type, thank you very much
				$data = generate_tim_data_structure(
					$feature_string, qw(Name Type) );
			}
			
			# add database
			$data->{'db'} = $database;
		}
	}
}



sub collect_features {
	
	# collection method
	my $method = shift;
	
	# Get the names of chromosomes to avoid
	my @excluded_chromosomes = 
		$TIM_CONFIG->param("$database\.chromosome_exclude");
	unless (@excluded_chromosomes) {
		@excluded_chromosomes = 
			$TIM_CONFIG->param('default_db.chromosome_exclude');
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
	push @{ $data->{'data_table'} }, [
		$seqfeature->seq_id,
		$seqfeature->start - 1,
		$seqfeature->end,
		$seqfeature->display_name || $seqfeature->primary_id || 'region',
		$seqfeature->score || 0,
		$seqfeature->strand >= 0 ? '+' : '-',
	];
	$data->{'last_row'} += 1;
	
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
	push @{ $data->{'data_table'} }, [
		$seqfeature->display_name || $seqfeature->primary_id,
		$seqfeature->type,
		$seqfeature->seq_id,
		$seqfeature->start,
		$seqfeature->end,
		$seqfeature->strand,
	];
	$data->{'last_row'} += 1;
	
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
	push @{ $data->{'data_table'} }, [
		$seqfeature->display_name || $seqfeature->primary_id,
		$seqfeature->type,
	];
	$data->{'last_row'} += 1;
	
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
	if ($seqfeature->strand >= 0) {
		# forward strand
		if ($start_adj) {
			my $start = $seqfeature->start;
			$start += $start_adj;
			$seqfeature->start($start);
		}
		if ($stop_adj) {
			my $stop = $seqfeature->end;
			$stop += $stop_adj;
			$seqfeature->end($stop);
		}
	}
	else {
		# reverse strand
		if ($start_adj) {
			my $start = $seqfeature->end;
			$start -= $start_adj;
			$seqfeature->start($start);
		}
		if ($stop_adj) {
			my $stop = $seqfeature->start;
			$stop -= $stop_adj;
			$seqfeature->start($stop);
		}
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
  --out <filename>
  --bed
  --gff 
  --(no)gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --db <text>

Specify the name or SQLite file of the Bio::DB:SeqFeature::Store 
database from which to collect the features. Other database schemas 
are not currently supported. Required.

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
upstream (5' direction), and a positive value is shifted downstream.
Adjustments are made relative to the feature's strand. Adjustments 
are ignored if a GFF file is written.

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

=item --(no)gz

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

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
