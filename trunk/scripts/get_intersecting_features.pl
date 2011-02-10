#!/usr/bin/perl

# a script to pull out overlapping features from the database

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::Range;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_db_helper qw(
	open_db_connection
	validate_dataset_list 
	get_dataset_list 
	validate_included_feature
);
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
);


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
);
GetOptions( 
	'in=s'       => \$infile, # the input nucleosome data file
	'db=s'       => \$database, # the name of the database
	'feature=s'  => \$search_feature, # the feature(s) to look for
	'start=i'    => \$start, # the relative start position
	'stop=i'     => \$stop, # the relative stop position
	'extend=i'   => \$extend, # extend lookup-feature by this amount
	'ref=s'      => \$reference_position, # relative position when calc distance
	'out=s'      => \$outfile, # output file name
	'gz!'        => \$gz, # compress file
	'help'       => \$help, # request help
);

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
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
if ($search_feature) {
	# check the requested feature type
	if ( validate_dataset_list($db, $search_feature) ) {
		die " Requested feature '$search_feature' does not appear to be valid!\n";
	}
}
else {
	# ask the user to respond from a list
	$search_feature = request_features_from_user();
}





### Collect the features
# search
print " searching for intersecting features....\n";
find_overlapping_features($main_data_ref, $search_feature);





### Write output
# write data file
unless ($outfile) {
	# overwrite the input file
	$outfile = $infile;
}
my $success = write_tim_data_file( {
	'data'      => $main_data_ref,
	'filename'  => $outfile,
} );
if ($success) {
	print " Wrote data file '$success'\n";
}
else {
	print " Unable to write data file!\n";
}


### The End





############################# Subroutines #####################################

### Collect features from the database for the user
sub request_features_from_user {
	
	# first get a list of all features in the database
	my %features = get_dataset_list($db, 'all');
	unless (%features) {
		die " no features in the database!?\n";
	}
	
	# present list to user
	print " These are the available features in the database '$database'\n";
	foreach (sort {$a <=> $b} keys %features) {
		print "   $_\t$features{$_}\n";
	}
	
	# collect and verify the response
	print " Enter the number of the feature to search    ";
	my $answer = <STDIN>;
	chomp $answer;
	if (exists $features{$answer}) {
		return $features{$answer};
	}
	else {
		die " unknown response\n";
	}
}







### 
sub find_overlapping_features {
	
	my ($data_ref, $search_feature) = @_;
	
	
	# identify the type of features in the list we're looking up
	# either genome coordinates or named features
	if (
		$data_ref->{0}{'name'} =~ /^chr|seq/i and 
		$data_ref->{1}{'name'} =~ /start/i and 
		$data_ref->{2}{'name'} =~ /stop|end/i
	) {
		# we're working with genomic coordinates here
		print " Input file '$infile' has genomic interval features\n";
		intersect_genome_features($data_ref, $search_feature);
	}
	elsif (
		$data_ref->{0}{'name'} =~ /name/i and
		$data_ref->{1}{'name'} =~ /type/i
	) {
		# we're working with named features
		print " Input file '$infile' has named features\n";
		intersect_named_features($data_ref, $search_feature);
	}
	else {
		# unable to identify
		die " unable to identify features in source file '$infile'\n Beginning". 
			" columns are not labeled chr|seq, start, stop|end, or name, type\n";
	}
}
	
	

### Working with named features 
sub intersect_named_features {
	# Named features 
	
	my ($data_ref, $search_feature) = @_;
	my $table = $data_ref->{'data_table'};

	# prepare new metadata columns 
	my ($number_i, $name_i, $type_i, $strand_i, $distance_i) = 
		generate_new_metadata($data_ref);
	
	
	# loop
	for (my $row = 1; $row <= $data_ref->{'last_row'}; $row++) {
		
		# identify feature first
		my @features = $db->features(
			-name   => $table->[$row][0],
			-type   => $table->[$row][1],
		);
		my $feature;
		if (scalar @features == 0) {
			warn " no features found for $table->[$row][1] " . 
				"$table->[$row][0]\n";
			process_no_feature(
				$data_ref, 
				$row, 
				$number_i, 
				$name_i, 
				$type_i, 
				$strand_i, 
				$distance_i
			);
			next;
		}
		elsif (scalar @features > 1) {
			warn " more than one feature found for $table->[$row][1] " .
				"$table->[$row][0]; using first one\n";
			$feature = shift @features;
		}
		else {
			# only one feature - as it should be
			$feature = shift @features;
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
					$feature->end + $stop
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
				$search_feature,
				$data_ref, 
				$row, 
				$number_i, 
				$name_i, 
				$type_i, 
				$strand_i, 
				$distance_i
			);
		}
		else {
			# no region defined
			warn " unable to establish region for $table->[$row][1] " . 
				"$table->[$row][0]\n";
			
			# fill in table anyway
			process_no_feature(
				$data_ref, 
				$row, 
				$number_i, 
				$name_i, 
				$type_i, 
				$strand_i, 
				$distance_i
			);
		}
	
	}
	
	# summarize the findings
	summarize_found_features($data_ref, $number_i);
}



### Working with genomic features 
sub intersect_genome_features {
	
	my ($data_ref, $search_feature) = @_;
	my $table = $data_ref->{'data_table'};

	# prepare new metadata columns 
	my ($number_i, $name_i, $type_i, $strand_i, $distance_i) = 
		generate_new_metadata($data_ref);
	
	
	# loop
	for (my $row = 1; $row <= $data_ref->{'last_row'}; $row++) {
		my $region;
		
		# extend the region
		if ($extend) {
			# we're adding an extension on either side of the region
			
			# new start
			my $new_start = $table->[$row][1] - $extend;
			if ($new_start < 1) {
				# limit to actual start
				$new_start = 1;
			}
			
			# new stop
			my $new_stop = $table->[$row][2] + $extend;
			my ($chr) = $db->get_features_by_name( $table->[$row][0] );
			if ($new_stop > $chr->length) {
				# limit to actual length
				$new_stop = $chr->length;
			}
			
			# establish region
			$region = $db->segment(
				$table->[$row][0], # chromosome
				$new_start,        # start
				$new_stop          # stop
			);
		}
		
		# specific relative start, stop
		if (defined $start and defined $stop) {
			# we'll adjust the coordinates specifically
			# this is relative to the start position
			
			# new start
			my $new_start = $table->[$row][1] + $start;
			if ($new_start < 1) {
				# limit to actual start
				$new_start = 1;
			}
			
			# new stop
			my $new_stop = $table->[$row][1] + $stop;
			
			# establish region
			$region = $db->segment(
				$table->[$row][0], # chromosome
				$new_start,        # start
				$new_stop          # stop
			);
		}
		
		# default is entire region
		else {
			
			# establish region as is
			$region = $db->segment(
				$table->[$row][0], # chromosome
				$table->[$row][1], # start
				$table->[$row][2]   # stop
			);
		}
		
		# check region
		if ($region) {
			# succesfully established a region, find features
			process_region(
				$region,
				0, # region inherently has no strand
				$search_feature,
				$data_ref, 
				$row, 
				$number_i, 
				$name_i, 
				$type_i, 
				$strand_i, 
				$distance_i
			);
		}
		else {
			# no region defined
			warn " unable to establish region for $table->[$row][0]:" . 
				"$table->[$row][1]..$table->[$row][2]\n";
			
			# fill in table anyway
			process_no_feature(
				$data_ref, 
				$row, 
				$number_i, 
				$name_i, 
				$type_i, 
				$strand_i, 
				$distance_i
			);
		}
			
	}
	
	# summarize the findings
	summarize_found_features($data_ref, $number_i);
}






### Prepare new columns
sub generate_new_metadata {
	my $data_ref = shift;
	
	# count of features column
	my $number_i = $data_ref->{'number_columns'};
	$data_ref->{$number_i} = {
		'index'        => $number_i,
		'name'         => 'Number_features',
		'intersection' => $data_ref->{'feature'},
	};
	$data_ref->{'data_table'}->[0][$number_i] = 'Number_features';
	$data_ref->{'number_columns'} += 1;
	
	# Name column
	my $name_i = $data_ref->{'number_columns'};
	$data_ref->{$name_i} = {
		'index'        => $name_i,
		'name'         => 'Name',
		'intersection' => $data_ref->{'feature'},
	};
	$data_ref->{'data_table'}->[0][$name_i] = 'Name';
	$data_ref->{'number_columns'} += 1;
	if (defined $start and defined $stop) {
		$data_ref->{$name_i}{'Start'} = $start;
		$data_ref->{$name_i}{'Stop'} = $stop;
	}
	if (defined $extend) {
		$data_ref->{$name_i}{'Extend'} = $extend;
	}
	
	# Type column
	my $type_i = $data_ref->{'number_columns'};
	$data_ref->{$type_i} = {
		'index'        => $type_i,
		'name'         => 'Type',
		'intersection' => $data_ref->{'feature'},
	};
	$data_ref->{'data_table'}->[0][$type_i] = 'Type';
	$data_ref->{'number_columns'} += 1;
	
	# Strand column
	my $strand_i = $data_ref->{'number_columns'};
	$data_ref->{$strand_i} = {
		'index'        => $strand_i,
		'name'         => 'Strand',
		'intersection' => $data_ref->{'feature'},
	};
	$data_ref->{'data_table'}->[0][$strand_i] = 'Strand';
	$data_ref->{'number_columns'} += 1;
	
	# Distance column
	my $distance_i = $data_ref->{'number_columns'};
	$data_ref->{$distance_i} = {
		'index'        => $distance_i,
		'name'         => 'Distance',
		'intersection' => $data_ref->{'feature'},
		'reference'    => $reference_position,
	};
	$data_ref->{'data_table'}->[0][$distance_i] = 'Distance';
	$data_ref->{'number_columns'} += 1;
			
	
	# add extra metadata
	
	
	return ($number_i, $name_i, $type_i, $strand_i, $distance_i);
}



### Prepare new columns
sub process_region {
	
	my ($region, $region_strand, $search_feature, $data_ref, $row, 
		$number_i, $name_i, $type_i, $strand_i, $distance_i) = @_;
	
	# look for the requested features
	my @features = $region->features(
		-type  => $search_feature,
	);
	#print "   found ", scalar(@features), "\n";
	
	# process depending on the number of features found
	if (scalar @features == 0) {
		# no features found
		
		# put in null data
		process_no_feature(
			$data_ref, 
			$row, 
			$number_i, 
			$name_i, 
			$type_i, 
			$strand_i, 
			$distance_i
		);
	}
	
	elsif (scalar @features == 1) {
		# only one feature found, that's perfect
		my $f = shift @features;
		
		# record information
		if ( validate_included_feature($f) ) {
			# the feature is ok to use (doesn't have tag to exclude it)
			$data_ref->{'data_table'}->[$row][$number_i]   = 1;
			$data_ref->{'data_table'}->[$row][$name_i]     = $f->display_name;
			$data_ref->{'data_table'}->[$row][$type_i]     = $f->type;
			$data_ref->{'data_table'}->[$row][$strand_i]   = $f->strand;
			$data_ref->{'data_table'}->[$row][$distance_i] = 
				determine_distance($region, $region_strand, $f);
		}
		else {
			# the feature should be excluded
			process_no_feature(
				$data_ref, 
				$row, 
				$number_i, 
				$name_i, 
				$type_i, 
				$strand_i, 
				$distance_i
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
				
				# apparently, using Bio::RangeI methods on $region doesn't 
				# work, intersection and overlap_extent fail to give proper 
				# end values, just returns "-end"
				
				# the workaround is to create a new simple Bio::Range object 
				# using the coordinates from $region, and then determine the 
				# overlap between it and the list of found target features in $_
				
				my $a = Bio::Range->new(
					-start  => $region->start,
					-end    => $region->end,
					-strand => $region_strand,
				);
				
				
				my $int = $a->intersection($_);
				my $overlap = $int->length;
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
		$data_ref->{'data_table'}->[$row][$number_i]   = scalar(@features);
		$data_ref->{'data_table'}->[$row][$name_i]     = $f->display_name;
		$data_ref->{'data_table'}->[$row][$type_i]     = $f->type;
		$data_ref->{'data_table'}->[$row][$strand_i]   = $f->strand;
		$data_ref->{'data_table'}->[$row][$distance_i] = 
			determine_distance($region, $region_strand, $f);
		
	}
	
	return;
}






### Fill in data table with null data
sub process_no_feature {
	
	my ($data_ref, $row, $number_i, $name_i, $type_i, $strand_i, $distance_i) =
		@_;

	$data_ref->{'data_table'}->[$row][$number_i]   = 0;
	$data_ref->{'data_table'}->[$row][$name_i]     = '.';
	$data_ref->{'data_table'}->[$row][$type_i]     = '.';
	$data_ref->{'data_table'}->[$row][$strand_i]   = 0;
	$data_ref->{'data_table'}->[$row][$distance_i] = '.';

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





sub summarize_found_features {
	my ($data_ref, $number_i) = @_;
	
	# intialize counts
	my $none     = 0;
	my $one      = 0;
	my $multiple = 0;
	
	# count up
	for (my $row = 1; $row <= $main_data_ref->{'last_row'}; $row++) {
		if ($data_ref->{'data_table'}->[$row][$number_i] == 0) {
			$none++;
		}
		elsif ($data_ref->{'data_table'}->[$row][$number_i] == 1) {
			$one++;
		}
		elsif ($data_ref->{'data_table'}->[$row][$number_i] > 1) {
			$multiple++;
		}
	}
	
	# print summary
	printf " $one (%.1f%%) reference features intersected with unique target features\n",
		(($one / $data_ref->{'last_row'}) * 100) if $one;
	printf " $none (%.1f%%) reference features intersected with zero target features\n",
		(($none / $data_ref->{'last_row'}) * 100) if $none;
	printf " $multiple (%.1f%%) reference features intersected with multiple target features\n",
		(($multiple / $data_ref->{'last_row'}) * 100) if $multiple;
}





__END__

=head1 NAME get_intersecting_features.pl



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
  --(no)gz
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a list of reference features to find overlapping 
target features. These reference features may be genomic coordinates (chromo,
start, stop) or named features (name, type). A tim data formatted file is 
best used but other tab delimited text formats may be used.

=item --db <database>

Provide the name of a Bio::DB::SeqFeature store database to use when 
finding features. If not specified, the database specified in the input file 
metadata will be used.

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

=item --(no)gz

Specify whether the output file should (not) be compressed with gzip.

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
of the features, or optionally the midpoints.

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











