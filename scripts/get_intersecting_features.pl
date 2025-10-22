#!/usr/bin/perl

# documentation at end of file

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use Bio::ToolBox::Data;
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
	get_chromosome_list
);
use Bio::ToolBox::utility qw(format_with_commas);

our $VERSION = '2.04';

print "\n A script to pull out overlapping features\n\n";

### Quick help
unless (@ARGV) {

	# when no command line options are present
	# print SYNOPSIS
	pod2usage(
		{
			'-verbose' => 0,
			'-exitval' => 1,
		}
	);
}

### Command line options
my ( $infile, $database, $start, $stop, $extend, $adjustment_position,
	$reference_position, $outfile, $gz, $help, $print_version, );
my @search_features;
GetOptions(
	'i|in=s'          => \$infile,             # the input nucleosome data file
	'd|db=s'          => \$database,           # the name of the database
	'f|feature=s'     => \@search_features,    # the feature(s) to look for
	'b|begin|start=i' => \$start,              # the relative start position
	'e|end|stop=i'    => \$stop,               # the relative stop position
	'p|pos=s'         =>
		\$adjustment_position,    # the coordinate to make start and stop adjustments
	'x|extend=i' => \$extend,                # extend lookup-feature by this amount
	'r|ref=s'    => \$reference_position,    # relative position when calc distance
	'o|out=s'    => \$outfile,               # output file name
	'z|gz!'      => \$gz,                    # compress file
	'h|help'     => \$help,                  # request help
	'v|version'  => \$print_version,         # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {

	# print entire POD
	pod2usage(
		{
			'-verbose' => 2,
			'-exitval' => 1,
		}
	);
}

# Print version
if ($print_version) {
	print " Biotoolbox script get_intersecting_features.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}

### Check for required values and assign defaults
unless ($infile) {
	if (@ARGV) {
		$infile = shift @ARGV;
	}
	else {
		print STDERR " FATAL: no input file! use --help for more information\n";
		exit 1;
	}
}
unless ($reference_position) {
	$reference_position = 'start';
}
unless ($adjustment_position) {
	$adjustment_position = '5';
}
unless ( defined $gz ) {
	$gz = 0;
}
if ( scalar @search_features == 1 and $search_features[0] =~ /,/ ) {

	# a comma-delimited list of features
	@search_features = split /,/, shift @search_features;
}

### Load the input file
my $Data = Bio::ToolBox::Data->new( file => $infile )
	or die "  Unable to load data file!\n";
printf " Loaded %s features from $infile.\n", format_with_commas( $Data->last_row );

### Open database connection
if ($database) {
	if ( $Data->database and $Data->database ne $database ) {
		print " WARNING: provided database '$database' does not match file metadata!\n"
			. " overriding metadata and using '$database'\n";
	}
	$Data->database($database);
}
elsif ( not $Data->database ) {
	die "No database defined! See help\n";
}
my $db = $Data->open_database;

### Identify the Features to Search
@search_features = verify_or_request_feature_types(
	'db'      => $db,
	'feature' => [@search_features],
	'prompt'  => "Enter the number(s) to the intersecting feature(s) to"
		. " search.\n Enter as comma delimited list and/or range   ",
);

### Collect the features
# search
print " Searching for intersecting " . join( ", ", @search_features ) . " features....\n";
find_overlapping_features();

### Write output
# write data file
unless ($outfile) {

	# overwrite the input file
	$outfile = $infile;
}
my $success = $Data->write_file(
	'filename' => $outfile,
	'gz'       => $gz,
);
if ($success) {
	print " Wrote data file '$success'\n";
}
else {
	print " Unable to write data file!\n";
}

### The End

############################# Subroutines #####################################

### Main starting point to find overlapping features
sub find_overlapping_features {

	if ( $Data->feature_type eq 'coordinate' ) {

		# we're working with genomic coordinates here
		intersect_genome_features();
	}
	elsif ( $Data->feature_type eq 'named' ) {

		# we're working with named features
		intersect_named_features();
	}
	else {
		# unable to identify
		print STDERR
			" FATAL: unable to identify feature information columns in source file "
			. "'$infile'\n No chromosome, start, stop, name, ID,  and/or type columns\n";
		exit 1;
	}
}

### Working with named features
sub intersect_named_features {

	# Named features

	# prepare new metadata columns
	my ( $number_i, $name_i, $type_i, $strand_i, $distance_i, $overlap_i ) =
		generate_new_metadata();

	# iterate through table
	my $stream = $Data->row_stream;
	while ( my $row = $stream->next_row ) {

		# identify feature first
		my $feature = $row->feature;
		unless ($feature) {
			process_no_feature( $row, $number_i, $name_i, $type_i, $strand_i,
				$distance_i, $overlap_i );
			next;
		}

		# Establish the region based on the found feature
		my $region;

		# extend the region
		if ($extend) {

			# we're adding an extension on either side of the feature

			# establish region
			$region = $db->segment(
				$feature->seq_id,
				$feature->start - $extend,
				$feature->end + $extend
			);

			# this segment will not have a strand, and I can not set it
		}

		# specific relative start, stop from indicated position
		elsif ( defined $start and defined $stop ) {

			# we'll adjust the coordinates specifically
			# this is relative to the indicated position

			# establish region based on the feature's orientation
			if ( $adjustment_position eq '5' and $feature->strand >= 0 ) {

				# Watson strand or unstranded
				$region = $db->segment(
					$feature->seq_id,
					$feature->start + $start,
					$feature->start + $stop
				);
			}
			elsif ( $adjustment_position eq '5' and $feature->strand < 0 ) {

				# Crick strand
				$region = $db->segment(
					$feature->seq_id,
					$feature->end - $stop,
					$feature->end - $start
				);
			}
			if ( $adjustment_position eq '3' and $feature->strand >= 0 ) {

				# Watson strand or unstranded
				$region = $db->segment(
					$feature->seq_id,
					$feature->end + $start,
					$feature->end + $stop
				);
			}
			elsif ( $adjustment_position eq '3' and $feature->strand < 0 ) {

				# Crick strand
				$region = $db->segment(
					$feature->seq_id,
					$feature->start - $stop,
					$feature->start - $start
				);
			}
			elsif ( $adjustment_position eq 'm' and $feature->strand >= 0 ) {

				# midpoint
				my $mid = int( ( $feature->start + $feature->end ) / 2 );
				$region = $db->segment( $feature->seq_id, $mid + $start, $mid + $stop );
			}
			elsif ( $adjustment_position eq 'm' and $feature->strand >= 0 ) {

				# midpoint
				my $mid = int( ( $feature->start + $feature->end ) / 2 );
				$region = $db->segment( $feature->seq_id, $mid - $start, $mid - $stop );
			}
		}

		# default is entire region
		else {
			# establish region as is
			$region = $feature->segment();

			# this segment will have strand
		}

		# check region
		if ($region) {

			# succesfully established a region, find features
			$region->start(1) if $region->start <= 0;    # sanity check
			process_region(
				$region,
				$feature->strand,    # the region may not support strand
									 # so need to pass this separately
				$row,
				$number_i,
				$name_i,
				$type_i,
				$strand_i,
				$distance_i,
				$overlap_i
			);
		}
		else {
			# no region defined
			printf " WARNING: unable to establish region for %s %s\n",
				$row->type, $row->name;

			# fill in table anyway
			process_no_feature( $row, $number_i, $name_i, $type_i, $strand_i,
				$distance_i, $overlap_i );
		}

	}

	# summarize the findings
	summarize_found_features($number_i);
}

### Working with genomic features
sub intersect_genome_features {

	# prepare new metadata columns
	my ( $number_i, $name_i, $type_i, $strand_i, $distance_i, $overlap_i ) =
		generate_new_metadata();

	# get chromosome list and their lengths
	# this is to ensure that regions don't go over length
	my %chrom2length;
	foreach ( get_chromosome_list($db) ) {

		# each element is [name, length]
		$chrom2length{ $_->[0] } = $_->[1];
	}

	# iterate through table
	my $stream = $Data->row_stream;
	while ( my $row = $stream->next_row ) {

		# adjust positions as necessary
		my ( $new_start, $new_stop );
		if ($extend) {

			# we're adding an extension on either side of the region
			$new_start = $row->start - $extend;
			$new_stop  = $row->end + $extend;
		}
		elsif ( defined $start and defined $stop ) {

			# we'll adjust the coordinates specifically
			# this is relative to the start position
			$new_start = $row->start + $start;
			$new_stop  = $row->start + $stop;
		}
		else {
			$new_start = $row->start;
			$new_stop  = $row->end;
		}

		# check new positions
		if ( $new_start < 1 ) {

			# limit to actual start
			$new_start = 1;
		}
		if ( $chrom2length{ $row->seq_id } and $new_stop > $chrom2length{ $row->seq_id } )
		{
			# limit to actual length
			$new_stop = $chrom2length{ $row->seq_id };
		}

		# establish region
		my $region = $db->segment( $row->seq_id, $new_start, $new_stop );

		# check region
		if ($region) {

			# succesfully established a region, find features
			process_region(
				$region,
				$row->strand,    # default is 0 if not defined
				$row,
				$number_i,
				$name_i,
				$type_i,
				$strand_i,
				$distance_i,
				$overlap_i
			);
		}
		else {
			# no region defined
			printf " WARNING: unable to establish region for %s:%s..%s\n",
				$row->seq_id, $row->start, $row->end;

			# fill in table anyway
			process_no_feature( $row, $number_i, $name_i, $type_i, $strand_i,
				$distance_i, $overlap_i );
		}
	}

	# summarize the findings
	summarize_found_features($number_i);
}

### Prepare new columns
sub generate_new_metadata {

	# count of features column
	my $number_i = $Data->add_column('Number_features');

	# Name column
	my $name_i = $Data->add_column('Target_Name');
	if ( defined $start and defined $stop ) {
		$Data->metadata( $name_i, 'Start', $start );
		$Data->metadata( $name_i, 'Stop',  $stop );
	}
	if ( defined $extend ) {
		$Data->metadata( $name_i, 'Extend', $extend );
	}

	# Type column
	my $type_i = $Data->add_column('Target_Type');

	# Strand column
	my $strand_i = $Data->add_column('Target_Strand');

	# Distance column
	my $distance_i = $Data->add_column('Target_Distance');
	$Data->metadata( $distance_i, 'reference', $reference_position );

	# Overlap column
	my $overlap_i = $Data->add_column('Target_Overlap');
	$Data->metadata( $overlap_i, 'reference', $reference_position );

	return ( $number_i, $name_i, $type_i, $strand_i, $distance_i, $overlap_i );
}

### Prepare new columns
sub process_region {

	my (
		$region, $region_strand, $row,        $number_i, $name_i,
		$type_i, $strand_i,      $distance_i, $overlap_i
	) = @_;

	# look for the requested features
	my @features = $region->features( -type => [@search_features], );

	#print "   found ", scalar(@features), "\n";

	# process depending on the number of features found
	if ( scalar @features == 0 ) {

		# no features found

		# put in null data
		process_no_feature( $row, $number_i, $name_i, $type_i, $strand_i,
			$distance_i, $overlap_i );
	}

	elsif ( scalar @features == 1 ) {

		# only one feature found, that's perfect
		my $f = shift @features;

		# record information
		$row->value( $number_i,   1 );
		$row->value( $name_i,     $f->display_name );
		$row->value( $type_i,     $f->type );
		$row->value( $strand_i,   $f->strand );
		$row->value( $distance_i, determine_distance( $region, $region_strand, $f ) );
		$row->value( $overlap_i,  determine_overlap( $region, $region_strand, $f ) );
	}

	elsif ( scalar @features > 1 ) {

		# more than one feature
		# need to identify the most appropriate one
		my $f;
		if ( scalar @features > 1 ) {
			my %overlap2f;
			foreach (@features) {
				my $overlap = determine_overlap( $region, $region_strand, $_ );
				$overlap2f{$overlap} = $_;

				# this may overwrite if two or more features have identical
				# amounts of overlap, but we'll simply use that as a means
				# of "randomly" selecting one of them ;)
				# otherwise how else to select one?
			}

			# now take the feature with the greatest overlap
			my $max = ( sort { $b <=> $a } keys %overlap2f )[0];
			$f = $overlap2f{$max};
		}
		else {
			# now there is only one, one of them must've been tossed because
			# it was dubious
			$f = $features[0];
		}

		# record the information
		$row->value( $number_i,   scalar(@features) );
		$row->value( $name_i,     $f->display_name );
		$row->value( $type_i,     $f->type );
		$row->value( $strand_i,   $f->strand );
		$row->value( $distance_i, determine_distance( $region, $region_strand, $f ) );
		$row->value( $overlap_i,  determine_overlap( $region, $region_strand, $f ) );
	}

	return;
}

### Fill in data table with null data
sub process_no_feature {
	my ( $row, $number_i, $name_i, $type_i, $strand_i, $distance_i, $overlap_i ) = @_;

	$row->value( $number_i,   0 );
	$row->value( $name_i,     '.' );
	$row->value( $type_i,     '.' );
	$row->value( $strand_i,   0 );
	$row->value( $distance_i, '.' );
	$row->value( $overlap_i,  '.' );
}

### Calculate the distance between target and reference features
sub determine_distance {

	my ( $reference, $ref_strand, $target ) = @_;

	# determine distance
	my $distance;

	# Calculating the distance between the 5' ends
	if ( $reference_position eq 'start' ) {

		# Calculation dependent on identifying the strand
		# basically calculating distance from target 5' end to reference 5' end
		if ( $ref_strand >= 0 and $target->strand >= 0 ) {

			# Both are Watson strand or unstranded
			$distance = $target->start - $reference->start;
		}
		elsif ( $ref_strand >= 0 and $target->strand < 0 ) {

			# Reference is Watson, target is Crick
			$distance = $target->end - $reference->start;
		}
		elsif ( $ref_strand < 0 and $target->strand >= 0 ) {

			# Reference is Crick, target is Watson
			$distance = $target->start - $reference->end;
		}
		elsif ( $ref_strand < 0 and $target->strand < 0 ) {

			# Both are Crick
			$distance = $target->end - $reference->end;
		}
	}

	# Calculating the distance between the midpoints
	# since midpoint is equidistant from 5' and 3' ends, strand doesn't matter
	elsif ( $reference_position eq 'mid' ) {
		my $reference_mid = ( $reference->end + $reference->start ) / 2;
		my $target_mid    = ( $target->end + $target->start ) / 2;
		$distance = int( $target_mid - $reference_mid + 0.5 );
	}

	return $distance;
}

### Calculate the overlap between target and reference features
sub determine_overlap {
	my ( $reference, $ref_strand, $target ) = @_;

	# apparently, using Bio::RangeI methods on a region doesn't
	# work (what!!!!????), intersection and overlap_extent fail to
	# give proper end values, just returns "-end"

	# the workaround is to just calculate the intersection ourself
	my ( $istart, $istop );
	if ( $reference->start <= $target->start ) {
		$istart = $reference->start;
	}
	else {
		$istart = $target->start;
	}
	if ( $reference->end >= $target->end ) {
		$istop = $reference->end;
	}
	else {
		$istop = $target->end;
	}
	return 0 if $istart > $istop;
	return ($istop - $istart + 1);
}

sub summarize_found_features {
	my $number_i = shift;

	# intialize counts
	my $none     = 0;
	my $one      = 0;
	my $multiple = 0;

	# count up
	$Data->iterate(
		sub {
			my $row = shift;
			if ( $row->value($number_i) == 0 ) {
				$none++;
			}
			elsif ( $row->value($number_i) == 1 ) {
				$one++;
			}
			elsif ( $row->value($number_i) > 1 ) {
				$multiple++;
			}
		}
	);

	# print summary
	printf " $one (%.1f%%) reference features intersected with unique target features\n",
		( ( $one / $Data->{'last_row'} ) * 100 )
		if $one;
	printf " $none (%.1f%%) reference features intersected with zero target features\n",
		( ( $none / $Data->{'last_row'} ) * 100 )
		if $none;
	printf
" $multiple (%.1f%%) reference features intersected with multiple target features\n",
		( ( $multiple / $Data->{'last_row'} ) * 100 )
		if $multiple;
}

__END__

=head1 NAME 

get_intersecting_features.pl

A program to pull out overlapping features from the database.

=head1 SYNOPSIS

get_intersecting_features.pl [--options] <filename>
  
  File options:
  -i --in <filename>                input file
  -o --out <filename>               optionally output file
  
  Database options:
  -d --db <database>                database to search: name or sqlite
  -f --feature <text>               db feature to search
  
  Modify search range:
  -b --begin --start <integer>      adjust relative search start coordinate
  -e --end --stop <integer>         adjust relative search stop coordinate
  -p --pos [5 | m | 3]              relative position of search coordinate
  -x --extend <integer>             extend search in both directions
  -r --ref [start | mid]            measure distance from which coordinate
  
  General options:
  -z --gz                           compress output
  -v --version                      print version and exit
  -h --help                         show extended documentation

=head1 OPTIONS

The command line flags and descriptions:

=head2 File options:

=over 4

=item --in E<lt>filenameE<gt>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. The 
first row should be column headers. Bed files are acceptable, as are 
text files generated by other B<BioToolBox> scripts. Files may be 
gzipped compressed.

=item --out E<lt>filenameE<gt>

Optionally specify a new filename. A standard tim data text file is written. 
The default is to rewrite the input file.

=back

=head2 Database options

=over 4

=item --db E<lt>databaseE<gt>

Specify the name of a L<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A database is 
required for generating new data files with features. This option may 
skipped when using coordinate information from an input file (e.g. BED 
file), or when using an existing input file with the database indicated 
in the metadata.  

=item --feature E<lt>textE<gt>

Specify the name of the target features to search for in the database that 
intersect with the list of reference features. The type may be a either a 
GFF "type" or a "type:method" string. If not specifed, then the database 
will be queried for potential GFF types and a list presented to the user to 
select one.

=back

=head2 Modify search range

=over 4

=item --start E<lt>integerE<gt>

=item --stop E<lt>integerE<gt>

=item --begin E<lt>integerE<gt>

=item --end E<lt>integerE<gt>

Optionally specify the relative start and stop positions from the 5' end 
(default) or the end specified by the "--pos" option with which to restrict 
the search region for target features. For example, specify "--start=-200 
--stop=0" to restrict to the promoter region of genes. Both positions must 
be specified. Default is to take the entire region of the reference feature.

=item --pos [ 5 | m | 3 ]

Indicate the relative position from which to make the adjustments to the 
search window. Both start and stop adjustments may be made from the 
respective 5 prime, 3 prime, or middle position as dictated by the feature's 
strand value. 

=item --extend E<lt>integerE<gt>

Optionally specify the number of bp to extend the reference feature's region 
on each side. Useful when you have small reference regions and you want to 
include a larger search region.

=item --ref [start | mid]

Indicate the reference point from which to calculate the distance between the 
reference and target features. The same reference point is used for both 
features. Valid options include "start" (or 5' end for stranded features) and 
"mid" (for midpoint). Default is "start".

=back

=head2 General options

=over 4

=item --gz

Specify whether the output file should (not) be compressed with gzip.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will take a list of reference features and identify 
target features which intersect them. The reference features 
may be either named features (name and type) or genomic regions (chromosome, 
start, stop). By default, the search region for each reference feature is the 
entire feature, but may be restricted or expanded in size with appropriate 
modifiers (--start, --stop, --extend). The target features are specifed as 
specific types. 

Several attributes of the found features are appended to the original input 
file data. First, the number of 
target features are reported. If more than one are found, the feature with 
the most overlap with the reference feature is preferentially listed. The name, 
type, and strand of the selected target feature is reported. Finally, the 
distance from the reference feature to the target feature is reported. The 
reference points for measuring the distance is by default the start or 5' end 
of the features, or optionally the midpoints. Note that the distance 
measurement is relative to the coordinates after adjustment with the --start, 
--stop, and --extend options.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
