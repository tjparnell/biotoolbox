#!/usr/bin/perl

# documentation at end of file

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::Range;
use Bio::ToolBox::data_helper qw(
	parse_list
	find_column_index
);
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types 
	get_chromosome_list
	validate_included_feature
	get_feature
);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	write_tim_data_file
);
my $VERSION = '1.15';


print "\n A script to pull out overlapping features\n\n";

### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}


### Command line options
my (
	$infile, 
	$database,
	$search_feature,
	$start,
	$stop,
	$extend,
	$reference_position,
	$outfile,
	$gz,
	$help,
	$print_version,
);
my @search_features;
GetOptions( 
	'in=s'       => \$infile, # the input nucleosome data file
	'db=s'       => \$database, # the name of the database
	'feature=s'  => \@search_features, # the feature(s) to look for
	'start=i'    => \$start, # the relative start position
	'stop=i'     => \$stop, # the relative stop position
	'extend=i'   => \$extend, # extend lookup-feature by this amount
	'ref=s'      => \$reference_position, # relative position when calc distance
	'out=s'      => \$outfile, # output file name
	'gz!'        => \$gz, # compress file
	'help'       => \$help, # request help
	'version'    => \$print_version, # print the version
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
	print " Biotoolbox script get_intersecting_features.pl, version $VERSION\n\n";
	exit;
}



### Check for required values and assign defaults
unless ($infile) {
	$infile = shift @ARGV || die " no input file specified\n";
}
unless ($reference_position) {
	$reference_position = 'start';
}
unless (defined $gz) {
	$gz = 0;
}
if (scalar @search_features == 1 and $search_features[0] =~ /,/) {
	# a comma-delimited list of features
	@search_features = split /,/, shift @search_features;
}



### Load the input file
my $main_data_ref = load_tim_data_file($infile) || 
	die "  Unable to load data file!\n";
print " loaded file '$main_data_ref->{filename}' with $main_data_ref->{last_row}" . 
	" features\n";




### Open database connection

# check database
unless ($database) {
	
	# check for database referenced in the input file metadata
	if (exists $main_data_ref->{'db'}) {
		$database = $main_data_ref->{'db'};
		print " using database '$database'\n";
	}
	else {
		# no database defined!
		die " no database defined!\n";
	}
}

# open connection
my $db = open_db_connection($database) || 
	die " unable to establish database connection!\n";





### Identify the Features to Search
@search_features = verify_or_request_feature_types(
	'db'      => $db,
	'feature' => [ @search_features ],
	'prompt'  => "Enter the number(s) to the intersecting feature(s) to" . 
				" search.\n Enter as comma delimited list and/or range   ",
);





### Collect the features
# search
print " Searching for intersecting " . join(", ", @search_features) . 
	" features....\n";
find_overlapping_features();





### Write output
# write data file
unless ($outfile) {
	# overwrite the input file
	$outfile = $infile;
}
my $success = write_tim_data_file(
	'data'      => $main_data_ref,
	'filename'  => $outfile,
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
	
	# identify the required indices of columns 
	# either genome coordinates or named features
	my $chrom_i  = find_column_index($main_data_ref, '^chr|seq');
	my $start_i  = find_column_index($main_data_ref, '^start');
	my $stop_i   = find_column_index($main_data_ref, '^stop|end');
	my $strand_i = find_column_index($main_data_ref, '^strand');
	my $name_i   = find_column_index($main_data_ref, '^name');
	my $type_i   = find_column_index($main_data_ref, '^type');
	my $id_i     = find_column_index($main_data_ref, '^primary_id');
	
	# genomic coordinates
	if (
		defined $chrom_i and 
		defined $start_i and 
		defined $stop_i
	) {
		# we're working with genomic coordinates here
		print " Input file '$infile' has genomic interval features\n";
		intersect_genome_features($chrom_i, $start_i, $stop_i, $strand_i);
	}
	
	# named database features
	elsif (
		defined $id_i or 
		(defined $name_i and defined $type_i)
	) {
		# we're working with named features
		print " Input file '$infile' has named features\n";
		intersect_named_features($id_i, $name_i, $type_i);
	}
	else {
		# unable to identify
		die " unable to identify feature information columns in source file " .
			"'$infile'\n No chromosome, start, stop, name, ID,  and/or type columns\n";
	}
}
	
	

### Working with named features 
sub intersect_named_features {
	# Named features 
	
	# search feature indices
	my ($search_id_i, $search_name_i, $search_type_i) = @_;
	
	# shortcut reference
	my $table = $main_data_ref->{'data_table'};

	# prepare new metadata columns 
	my ($number_i, $name_i, $type_i, $strand_i, $distance_i, $overlap_i) = 
		generate_new_metadata();
	
	
	# loop
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		
		# identify feature first
		my $feature = get_feature(
			'db'    => $db,
			'id'    => defined $search_id_i   ? $table->[$row][$search_id_i]   : undef,
			'name'  => defined $search_name_i ? $table->[$row][$search_name_i] : undef,
			'type'  => defined $search_type_i ? $table->[$row][$search_type_i] : undef,
		);
		unless ($feature) {
			process_no_feature(
				$row, 
				$number_i, 
				$name_i, 
				$type_i, 
				$strand_i, 
				$distance_i,
				$overlap_i,
			);
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
		
		# specific relative start, stop from feature 5' position
		elsif (defined $start and defined $stop) {
			# we'll adjust the coordinates specifically
			# this is relative to the start position
			
			# establish region based on the feature's orientation
			if ($feature->strand >= 0) {
				# Watson strand or unstranded
				$region = $db->segment(
					$feature->seq_id,
					$feature->start + $start, 
					$feature->start + $stop
				);
			}
			elsif ($feature->strand < 0) {
				# Crick strand
				$region = $db->segment(
					$feature->seq_id,
					$feature->end - $stop, 
					$feature->end - $start
				);
			}	
			# this segment will not have a strand, and I can not set it 
				
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
			process_region(
				$region,
				$feature->strand, 	# the region may not support strand
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
			warn " unable to establish region for $table->[$row][$search_type_i] " . 
				"$table->[$row][$search_name_i]\n";
			
			# fill in table anyway
			process_no_feature(
				$row, 
				$number_i, 
				$name_i, 
				$type_i, 
				$strand_i, 
				$distance_i,
				$overlap_i,
			);
		}
	
	}
	
	# summarize the findings
	summarize_found_features($number_i);
}



### Working with genomic features 
sub intersect_genome_features {
	
	# search feature indices
	my ($search_chrom_i, $search_start_i, $search_stop_i, $search_strand_i) = @_;
	
	# shortcut
	my $table = $main_data_ref->{'data_table'};

	# prepare new metadata columns 
	my ($number_i, $name_i, $type_i, $strand_i, $distance_i, $overlap_i) = 
		generate_new_metadata();
	
	# get chromosome list and their lengths
	# this is to ensure that regions don't go over length
	my %chrom2length;
	foreach (get_chromosome_list($db)) {
		# each element is [name, length]
		$chrom2length{ $_->[0] } = $_->[1];
	}
	
	# loop
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		my $region;
		
		# extend the region
		if ($extend) {
			# we're adding an extension on either side of the region
			
			# new start
			my $new_start = $table->[$row][$search_start_i] - $extend;
			if ($new_start < 1) {
				# limit to actual start
				$new_start = 1;
			}
			
			# new stop
			my $new_stop = $table->[$row][$search_stop_i] + $extend;
			if ($new_stop > $chrom2length{ $table->[$row][$search_chrom_i] }) {
				# limit to actual length
				$new_stop = $chrom2length{ $table->[$row][$search_chrom_i] };
			}
			
			# establish region
			$region = $db->segment(
				$table->[$row][$search_chrom_i], # chromosome
				$new_start,        # start
				$new_stop          # stop
			);
		}
		
		# specific relative start, stop
		if (defined $start and defined $stop) {
			# we'll adjust the coordinates specifically
			# this is relative to the start position
			
			# new start
			my $new_start = $table->[$row][$search_start_i] + $start;
			if ($new_start < 1) {
				# limit to actual start
				$new_start = 1;
			}
			
			# new stop
			my $new_stop = $table->[$row][$search_start_i] + $stop;
			
			# establish region
			$region = $db->segment(
				$table->[$row][$search_chrom_i], # chromosome
				$new_start,        # start
				$new_stop          # stop
			);
		}
		
		# default is entire region
		else {
			
			# establish region as is
			$region = $db->segment(
				$table->[$row][$search_chrom_i], # chromosome
				$table->[$row][$search_start_i], # start
				$table->[$row][$search_stop_i]   # stop
			);
		}
		
		# check region
		if ($region) {
			# succesfully established a region, find features
			process_region(
				$region,
				defined $search_strand_i ? $table->[$row][$search_strand_i] : 0, 
					# use strand if available in source data file, otherwise
					# it is non-stranded region
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
			warn " unable to establish region for " . 
				$table->[$row][$search_chrom_i] . ":" . 
				$table->[$row][$search_start_i] . ".." . 
				$table->[$row][$search_stop_i]. "\n";
			
			# fill in table anyway
			process_no_feature(
				$row, 
				$number_i, 
				$name_i, 
				$type_i, 
				$strand_i, 
				$distance_i, 
				$overlap_i
			);
		}
			
	}
	
	# summarize the findings
	summarize_found_features($number_i);
}






### Prepare new columns
sub generate_new_metadata {
	
	# count of features column
	my $number_i = $main_data_ref->{'number_columns'};
	$main_data_ref->{$number_i} = {
		'index'        => $number_i,
		'name'         => 'Number_features',
		'intersection' => $main_data_ref->{'feature'},
	};
	$main_data_ref->{'data_table'}->[0][$number_i] = 'Number_features';
	$main_data_ref->{'number_columns'} += 1;
	
	# Name column
	my $name_i = $main_data_ref->{'number_columns'};
	$main_data_ref->{$name_i} = {
		'index'        => $name_i,
		'name'         => 'Name',
		'intersection' => $main_data_ref->{'feature'},
	};
	$main_data_ref->{'data_table'}->[0][$name_i] = 'Name';
	$main_data_ref->{'number_columns'} += 1;
	if (defined $start and defined $stop) {
		$main_data_ref->{$name_i}{'Start'} = $start;
		$main_data_ref->{$name_i}{'Stop'} = $stop;
	}
	if (defined $extend) {
		$main_data_ref->{$name_i}{'Extend'} = $extend;
	}
	
	# Type column
	my $type_i = $main_data_ref->{'number_columns'};
	$main_data_ref->{$type_i} = {
		'index'        => $type_i,
		'name'         => 'Type',
		'intersection' => $main_data_ref->{'feature'},
	};
	$main_data_ref->{'data_table'}->[0][$type_i] = 'Type';
	$main_data_ref->{'number_columns'} += 1;
	
	# Strand column
	my $strand_i = $main_data_ref->{'number_columns'};
	$main_data_ref->{$strand_i} = {
		'index'        => $strand_i,
		'name'         => 'Strand',
		'intersection' => $main_data_ref->{'feature'},
	};
	$main_data_ref->{'data_table'}->[0][$strand_i] = 'Strand';
	$main_data_ref->{'number_columns'} += 1;
	
	# Distance column
	my $distance_i = $main_data_ref->{'number_columns'};
	$main_data_ref->{$distance_i} = {
		'index'        => $distance_i,
		'name'         => 'Distance',
		'intersection' => $main_data_ref->{'feature'},
		'reference'    => $reference_position,
	};
	$main_data_ref->{'data_table'}->[0][$distance_i] = 'Distance';
	$main_data_ref->{'number_columns'} += 1;
			
	
	# Overlap column
	my $overlap_i = $main_data_ref->{'number_columns'};
	$main_data_ref->{$overlap_i} = {
		'index'        => $overlap_i,
		'name'         => 'Overlap',
		'intersection' => $main_data_ref->{'feature'},
		'reference'    => $reference_position,
	};
	$main_data_ref->{'data_table'}->[0][$overlap_i] = 'Overlap';
	$main_data_ref->{'number_columns'} += 1;
	
	
	return ($number_i, $name_i, $type_i, $strand_i, $distance_i, $overlap_i);
}



### Prepare new columns
sub process_region {
	
	my ($region, $region_strand, $row, $number_i, $name_i, 
		$type_i, $strand_i, $distance_i, $overlap_i) = @_;
	
	# look for the requested features
	my @features = $region->features(
		-type  => [ @search_features ],
	);
	#print "   found ", scalar(@features), "\n";
	
	# process depending on the number of features found
	if (scalar @features == 0) {
		# no features found
		
		# put in null data
		process_no_feature(
			$row, 
			$number_i, 
			$name_i, 
			$type_i, 
			$strand_i, 
			$distance_i,
			$overlap_i,
		);
	}
	
	elsif (scalar @features == 1) {
		# only one feature found, that's perfect
		my $f = shift @features;
		
		# record information
		if ( validate_included_feature($f) ) {
			# the feature is ok to use (doesn't have tag to exclude it)
			$main_data_ref->{'data_table'}->[$row][$number_i]   = 1;
			$main_data_ref->{'data_table'}->[$row][$name_i]     = $f->display_name;
			$main_data_ref->{'data_table'}->[$row][$type_i]     = $f->type;
			$main_data_ref->{'data_table'}->[$row][$strand_i]   = $f->strand;
			$main_data_ref->{'data_table'}->[$row][$distance_i] = 
				determine_distance($region, $region_strand, $f);
			$main_data_ref->{'data_table'}->[$row][$overlap_i] = 
				determine_overlap($region, $region_strand, $f);
		}
		else {
			# the feature should be excluded
			process_no_feature(
				$row, 
				$number_i, 
				$name_i, 
				$type_i, 
				$strand_i, 
				$distance_i,
				$overlap_i,
			);
		}
	}
	
	elsif (scalar @features > 1) {
		# more than one feature
		# need to identify the most appropriate one
		
		# first check for excluded feature tags
		for (my $i = $#features - 1; $i >= 0; $i -= 1) {
			# walk through the list and delete any that should be excluded
			# we're working backwards to avoid indexing problems with splice
			unless ( validate_included_feature($features[$i]) ) {
				splice(@features, $i, 1);
			}
		}
		
		
		# next we'll take the one with the most overlap first
		my $f;
		if (scalar @features > 1) {
			my %overlap2f;
			foreach (@features) {
				my $overlap = determine_overlap($region, $region_strand, $_);
				$overlap2f{$overlap} = $_;
				# this may overwrite if two or more features have identical 
				# amounts of overlap, but we'll simply use that as a means 
				# of "randomly" selecting one of them ;)
				# otherwise how else to select one?
			}
			
			# now take the feature with the greatest overlap
			my $max = (sort {$b <=> $a} keys %overlap2f)[0];
			$f = $overlap2f{$max};
		}
		else {
			# now there is only one, one of them must've been tossed because
			# it was dubious
			$f = $features[0];
		}
		
		# record the information
		$main_data_ref->{'data_table'}->[$row][$number_i]   = scalar(@features);
		$main_data_ref->{'data_table'}->[$row][$name_i]     = $f->display_name;
		$main_data_ref->{'data_table'}->[$row][$type_i]     = $f->type;
		$main_data_ref->{'data_table'}->[$row][$strand_i]   = $f->strand;
		$main_data_ref->{'data_table'}->[$row][$distance_i] = 
			determine_distance($region, $region_strand, $f);
		$main_data_ref->{'data_table'}->[$row][$overlap_i] = 
			determine_overlap($region, $region_strand, $f);
		
	}
	
	return;
}






### Fill in data table with null data
sub process_no_feature {
	
	my ($row, $number_i, $name_i, $type_i, $strand_i, 
		$distance_i, $overlap_i) = @_;

	$main_data_ref->{'data_table'}->[$row][$number_i]   = 0;
	$main_data_ref->{'data_table'}->[$row][$name_i]     = '.';
	$main_data_ref->{'data_table'}->[$row][$type_i]     = '.';
	$main_data_ref->{'data_table'}->[$row][$strand_i]   = 0;
	$main_data_ref->{'data_table'}->[$row][$distance_i] = '.';
	$main_data_ref->{'data_table'}->[$row][$overlap_i]  = '.';

}





### Calculate the distance between target and reference features
sub determine_distance {
	
	my ($reference, $ref_strand, $target) = @_;
	
	# determine distance
	my $distance;
	
	# Calculating the distance between the 5' ends
	if ($reference_position eq 'start') {
		
		# Calculation dependent on identifying the strand
		# basically calculating distance from target 5' end to reference 5' end
		if ($ref_strand >= 0 and $target->strand >= 0) {
			# Both are Watson strand or unstranded
			$distance = $target->start - $reference->start;
		}
		elsif ($ref_strand >= 0 and $target->strand < 0) {
			# Reference is Watson, target is Crick
			$distance = $target->end - $reference->start;
		}
		elsif ($ref_strand < 0 and $target->strand >= 0) {
			# Reference is Crick, target is Watson
			$distance = $target->start - $reference->end;
		}
		elsif ($ref_strand < 0 and $target->strand < 0) {
			# Both are Crick
			$distance = $target->end - $reference->end;
		}
	}
	
	# Calculating the distance between the midpoints
	# since midpoint is equidistant from 5' and 3' ends, strand doesn't matter
	elsif ($reference_position eq 'mid') {
		my $reference_mid = ($reference->end + $reference->start) / 2;
		my $target_mid = ($target->end + $target->start) / 2;
		$distance = int($target_mid - $reference_mid + 0.5);
	}
	
	return $distance;
}



### Calculate the overlap between target and reference features
sub determine_overlap {
	my ($reference, $ref_strand, $target) = @_;
	
	# apparently, using Bio::RangeI methods on a region doesn't 
	# work (what!!!!????), intersection and overlap_extent fail to 
	# give proper end values, just returns "-end"
	
	# the workaround is to create a new simple Bio::Range object 
	# using the coordinates from the reference region, and then determine the 
	# overlap between it and the target feature
	my $a = Bio::Range->new(
		-start  => $reference->start,
		-end    => $reference->end,
		-strand => $ref_strand,
	);
	
	# find the overlap
	my $int = $a->intersection($target);
	return $int->length;
}



sub summarize_found_features {
	my $number_i = shift;
	
	# intialize counts
	my $none     = 0;
	my $one      = 0;
	my $multiple = 0;
	
	# count up
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		if ($main_data_ref->{'data_table'}->[$row][$number_i] == 0) {
			$none++;
		}
		elsif ($main_data_ref->{'data_table'}->[$row][$number_i] == 1) {
			$one++;
		}
		elsif ($main_data_ref->{'data_table'}->[$row][$number_i] > 1) {
			$multiple++;
		}
	}
	
	# print summary
	printf " $one (%.1f%%) reference features intersected with unique target features\n",
		(($one / $main_data_ref->{'last_row'}) * 100) if $one;
	printf " $none (%.1f%%) reference features intersected with zero target features\n",
		(($none / $main_data_ref->{'last_row'}) * 100) if $none;
	printf " $multiple (%.1f%%) reference features intersected with multiple target features\n",
		(($multiple / $main_data_ref->{'last_row'}) * 100) if $multiple;
}





__END__

=head1 NAME 

get_intersecting_features.pl

A script to pull out overlapping features from the database.

=head1 SYNOPSIS

get_intersecting_features.pl [--options] <filename>
  
  Options:
  --in <filename>
  --db <database>
  --feature <text>
  --start <integer>
  --stop <integer>
  --extend <integer>
  --ref [start | mid]
  --out <filename>
  --gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. The 
first row should be column headers. Bed files are acceptable, as are 
text files generated by other B<BioToolBox> scripts. Files may be 
gzipped compressed.

=item --db <database>

Specify the name of a C<Bio::DB::SeqFeature::Store> annotation database 
from which gene or feature annotation may be derived. A database is 
required for generating new data files with features. This option may 
skipped when using coordinate information from an input file (e.g. BED 
file), or when using an existing input file with the database indicated 
in the metadata. For more information about using annotation databases, 
see L<https://code.google.com/p/biotoolbox/wiki/WorkingWithDatabases>. 

=item --feature <text>

Specify the name of the target features to search for in the database that 
intersect with the list of reference features. The type may be a either a 
GFF "type" or a "type:method" string. If not specifed, then the database 
will be queried for potential GFF types and a list presented to the user to 
select one.

=item --start <integer>, --stop <integer>

Optionally specify the relative start and stop positions from the 5' end (or 
start coordinate for non-stranded features) with which to restrict the region 
when searching for target features. For example, specify "--start=-200 
--stop=0" to restrict to the promoter region of genes. Both positions must 
be specified. Default is to take the entire region of the reference feature.

=item --extend <integer>

Optionally specify the number of bp to extend the reference feature's region 
on each side. Useful when you have small reference regions and you want to 
include a larger search region.

=item --ref [start | mid]

Indicate the reference point from which to calculate the distance between the 
reference and target features. The same reference point is used for both 
features. Valid options include "start" (or 5' end for stranded features) and 
"mid" (for midpoint). Default is "start".

=item --out <filename>

Optionally specify a new filename. A standard tim data text file is written. 
The default is to rewrite the input file.

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

A standard tim data text file is written.

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
