#!/usr/bin/env perl

# documentation at end of file

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::ToolBox::data_helper qw(
	find_column_index
	parse_list
);
use Bio::ToolBox::db_helper qw(
	open_db_connection
	get_feature
);
use Bio::ToolBox::file_helper qw(
	load_tim_data_file
	write_tim_data_file
);
my $VERSION = '1.14';

print "\n This script will collect information for a list of features\n\n";



### Quick help
unless (@ARGV) { # when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Get command line options and initialize values

# Initialize values
my (
	$infile,
	$outfile,
	$database,
	$attrib_request,
	$use_type,
	$gz,
	$help,
	$print_version,
); 

# Command line options
GetOptions( 
	'in=s'     => \$infile, # input file
	'attrib=s' => \$attrib_request, # attribute
	'out=s'    => \$outfile, # output filename
	'db=s'     => \$database, # database name
	'type=s'   => \$use_type, # force a type
	'gz!'      => \$gz, # gzip status
	'help'     => \$help, # help
	'version'  => \$print_version, # print the version
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
	print " Biotoolbox script get_feature_info.pl, version $VERSION\n\n";
	exit;
}



# Check for required values
unless ($infile) {
	$infile = shift @ARGV;
}
unless (defined $gz) {$gz = 0}




### Load the feature list

# load file
print " Loading feature list from '$infile'....\n";
my $main_data_ref = load_tim_data_file($infile);
unless ($main_data_ref) {
	die " No file data loaded!\n";
}

# identify indices
my $name_index = find_column_index($main_data_ref, '^name');
my $type_index = find_column_index($main_data_ref, '^type');
my $id_index   = find_column_index($main_data_ref, '^primary_id');



### Establish database connection
unless ($database) {
	# define database if it wasn't on the command line
	$database = $main_data_ref->{'db'};
}
my $db = open_db_connection($database);





### Determine the attribute to collect

# get the attribute list from the user
my @attribute_list = get_attribute_list_from_user();

# process the requests
collect_attributes_for_list(@attribute_list);







### Output the data
# assign file name if necessary
unless ($outfile) {
	$outfile = $infile;
}

# write the file
my $file_success = write_tim_data_file(
	'data'      => $main_data_ref,
	'filename'  => $outfile,
);
if ($file_success) {
	# success
	print " Wrote file '$file_success'\n";
}
else {
	# failure
	print " Unable to write output file!\n";
}
print " That's it!\n";







### Subroutines


sub get_attribute_list_from_user {
	# get the list from the user
	# either as a command line option or interactively
	
	my @list;
	
	# provided by command line argument
	if ($attrib_request) {
		@list = split /,/, $attrib_request;
	}
	
	# request interactively from user
	else {
		
		# get the list of features to check for examples
		
		# get the attributes for a sample of features
		# store the tag keys in an example hash
		my %tagexamples;
		for (my $i = 1; $i < 50; $i++) {
			last if $i == $main_data_ref->{'last_row'};
			
			my @examples = $db->features(
				-name     => $main_data_ref->{'data_table'}->[$i][$name_index],
				-type     => $main_data_ref->{'data_table'}->[$i][$type_index]
			);
			unless (@examples) {
				next;
			}
			foreach my $example (@examples) {
				my %taghash = $example->attributes();
				foreach (keys %taghash) {
					$tagexamples{$_} += 1;
				}
			}
		}
		
		# present list to user
		print " These are the attributes which may be collected:\n";
		my $i = 1;
		my %index2att;
		# standard attributes for any user
		foreach ( 
			qw(Chromosome Start Stop Strand Score Name Alias Note Type Primary_tag Source 
				Length Midpoint Phase RNA_count Exon_count Transcript_length Parent 
				Primary_ID 
			) 
		) {
			print "   $i\t$_\n";
			$index2att{$i} = $_;
			$i++;
		}
		# specific attributes for these features
		foreach (sort {$a cmp $b} keys %tagexamples) {
			# all other attributes
			print "   $i\t$_\n";
			$index2att{$i} = $_;
			$i++;
		}
		
		# collect the user response
		print " Enter the attribute number(s) to collect, comma-delimited or range  ";
		my $answer = <STDIN>;
		chomp $answer;
		my @answers = parse_list($answer);
		
		# check the answers
		foreach (@answers) {
			if (exists $index2att{$_}) {
				push @list, $index2att{$_};
			}
			else {
				warn " unknown response '$_'!\n";
			}
		}
	}
	
	return @list;
}



sub get_attribute_method {
	# a subroutine to get the appropriate subroutine method
	
	my $attrib = shift;
	
	# set the appropriate attribute collection subroutine
	my $method;
	if ($attrib =~ /^chromo/i) {
		$method = \&get_chromo;
	} 
	elsif ($attrib =~ /^start$/i) {
		$method = \&get_start;
	} 
	elsif ($attrib =~ /^stop$/i) {
		$method = \&get_stop;
	} 
	elsif ($attrib =~ /^midpoint$/i) {
		$method = \&get_midpoint;
	} 
	elsif ($attrib =~ /^length$/i) {
		$method = \&get_length;
	} 
	elsif ($attrib =~ /^transcript_length$/i) {
		$method = \&get_transcript_length;
	} 
	elsif ($attrib =~ /^strand$/i) {
		$method = \&get_strand;
	} 
	elsif ($attrib =~ /^phase$/i) {
		$method = \&get_phase;
	} 
	elsif ($attrib =~ /^score$/i) {
		$method = \&get_score;
	} 
	elsif ($attrib =~ /^rna_count$/i) {
		$method = \&get_rna_number;
	} 
	elsif ($attrib =~ /^exon_count$/i) {
		$method = \&get_exon_number;
	} 
	elsif ($attrib =~ /^parent$/i) {
		$method = \&get_parent;
	} 
	elsif ($attrib =~ /^primary_id$/i) {
		$method = \&get_primary_id;
	} 
	elsif ($attrib =~ /^name$/i) {
		$method = \&get_display_name;
	} 
	elsif ($attrib =~ /^source$/i) {
		$method = \&get_source;
	} 
	elsif ($attrib =~ /^type$/i) {
		$method = \&get_type;
	} 
	elsif ($attrib =~ /^primary_tag$/i) {
		$method = \&get_primary_tag;
	} 
	else {
		# unrecognized, must be tag key
		$method = \&get_tag_value;
	}
		
	return $method;
}


sub collect_attributes_for_list {
	
	my @list = @_;
	
	# get the attribute method(s)
	my @methods; # the methods are references to appropriate subroutines
	foreach (@list) {
		push @methods, get_attribute_method($_);
	}
	
	print " Retrieving ", join(", ", @list), "\n";
	my $table = $main_data_ref->{'data_table'}; # shortcut reference
	for my $row (1..$main_data_ref->{'last_row'}) {
		
		# get the name of the feature
		my $name = $table->[$row][$name_index];
		$name = (split(';', $name))[0] if $name =~ /;/; # take the first name only
		
		# pull the feature(s) from the database
		my $feature = get_feature(
			'db'    => $db,
			'name'  => defined $name_index ? 
				$main_data_ref->{'data_table'}->[$row][$name_index] : undef,
			'type'  => defined $type_index ? 
				$main_data_ref->{'data_table'}->[$row][$type_index] : 
				defined $use_type ? $use_type : undef,
			'id'    => defined $id_index ? 
				$main_data_ref->{'data_table'}->[$row][$id_index] : undef,
		);
		
		# get the attribute(s)
		if ($feature) {
			for (my $i = 0; $i < scalar @list; $i++) {
				# for each request in the list, we will collect the attribute
				# we'll use the method sub defined in the methods
				# pass both the feature and the name of the attribute
				# only the tag value actually needs the name of the attribute
				push @{ $table->[$row] }, 
					&{ $methods[$i] }($feature, $list[$i]);
			}
		}
		else {
			# no feature found, cannot collect attributes
			# record nulls
			for (my $i = 0; $i < scalar @list; $i++) {
				push @{ $table->[$row] }, '.';
			}
		}
	}
	
	# record the metadata
	foreach (@list) {
		record_metadata($_);
	}
}




sub get_chromo {
	my $feature = shift;
	return '.' unless $feature;
	return $feature->seq_id || '.';
}


sub get_start {
	my $feature = shift;
	return '.' unless $feature;
	return $feature->start || '.';
}



sub get_stop {
	my $feature = shift;
	return '.' unless $feature;
	return $feature->end || '.';
}


sub get_length {
	my $feature = shift;
	return 0 unless $feature;
	return $feature->length || 0;
}


sub get_transcript_length {
	my $feature = shift;
	return 0 unless $feature;
	my $exon_total = 0;
	my $cds_total  = 0;
	foreach my $subf ( $feature->get_SeqFeatures() ) {
		# feature may consist of multiple subfeature types
		# we're only interested in the exon subfeatures, or if those don't
		# exist, then the CDS subfeatures
		# ignore all other subfeatures
		
		# exon subfeature
		if ($subf->primary_tag eq 'exon') {
			$exon_total += $subf->length;
		}
		# cds subfeature
		elsif ($subf->primary_tag eq 'CDS') {
			$cds_total += $subf->length;
		}
	}
	# return most appropriate number
	return $exon_total > 0 ? $exon_total : $cds_total;
}


sub get_midpoint {
	my $feature = shift;
	return '.' unless $feature;
	return ($feature->start + int( $feature->length / 2) ) || '.';
}


sub get_strand {
	my $feature = shift;
	return '.' unless $feature;
	return $feature->strand || 0;
}


sub get_phase {
	my $feature = shift;
	return '.' unless $feature;
	return $feature->phase || '.';
}


sub get_score {
	my $feature = shift;
	return '.' unless $feature;
	return $feature->score || '.';
}


sub get_rna_number {
	my $feature = shift;
	return 0 unless $feature;
	my $rna_count = 0;
	foreach my $f ($feature->get_SeqFeatures) {
		if ($f->primary_tag =~ /rna/i) {
			# an RNA transcript
			$rna_count++;
		}
	}
	return $rna_count;
}


sub get_exon_number {
	my $feature = shift;
	return 0 unless $feature;
	my $exon_count = 0;
	my $cds_count = 0;
	foreach my $f ($feature->get_SeqFeatures) {
		# count both exons and CDSs
		if ($f->primary_tag =~ /exon/i) {
			$exon_count++;
		}
		elsif ($f->primary_tag =~ /cds/i) {
			$cds_count++;
		}
		elsif ($f->primary_tag =~ /rna/i) {
			# an RNA transcript, go one more level
			foreach my $sf ($f->get_SeqFeatures) {
				if ($sf->primary_tag =~ /exon/i) {
					$exon_count++;
				}
				elsif ($sf->primary_tag =~ /cds/i) {
					$cds_count++;
				}
			}
		}
	}
	# return exon_count if non-zero, else return cds_count, zero or non-zero
	return $exon_count ? $exon_count : $cds_count;
}


sub get_parent {
	my $feature = shift;
	return '.' unless $feature;
	if ($feature->has_tag('parent_id')) {
		# feature has a parent
		my ($parent_id) = $feature->get_tag_values('parent_id');
		# unfortunately we have to do a slow lookup through the database 
		# using the attribute load_id
		# this seems to be the only way to get to the parent name
		my @parents = $db->get_features_by_attribute('load_id' => $parent_id);
		if (@parents) {
			my $parent = shift @parents;
			if (@parents) {
				warn " more than feature found with load_id $parent_id!\n";
			}
			return $parent->display_name;
		}
		else {
			warn " can't find a parent with load_id '$parent_id' for " . 
				$feature->display_name . "!\n";
			return '.';
		}
	}
	else {
		return '.';
	}
}


sub get_primary_id {
	my $feature = shift;
	return $feature->primary_id || '.';
}


sub get_display_name {
	my $feature = shift;
	return $feature->display_name || '.';
}


sub get_type {
	my $feature = shift;
	return $feature->type || $feature->primary_tag || '.';
}


sub get_primary_tag {
	my $feature = shift;
	return $feature->primary_tag || $feature->type || '.';
}


sub get_source {
	my $feature = shift;
	return $feature->source_tag || '.';
}

sub get_tag_value {
	my $feature = shift;
	my $attrib = shift;
	return '.' unless $feature;
	if ($feature->has_tag($attrib)) {
		return join(';', ($feature->get_tag_values($attrib)) );
	}
	else {
		return '.';
	}
}


# subroutine to record the metadata for each new attribute dataset
sub record_metadata {
	my $attrib = shift;
	
	# determine new index
	my $new_index = $main_data_ref->{'number_columns'};
	# remember that the index counting is 0-based, so the new index is 
	# essentially number_columns - 1 + 1
	$main_data_ref->{'number_columns'} += 1; # update
	
	# generate new metadata hash for this column
	my %metadata = (
		'name'     => $attrib,
		'index'    => $new_index,
	);
	
	# add database name if different
	if ($database ne $main_data_ref->{'db'}) {
		$metadata{'db'} = $database;
	}
	
	# generate column name
	$main_data_ref->{'data_table'}->[0][$new_index] = $attrib;
	
	# place metadata hash into main data structure
	$main_data_ref->{$new_index} = \%metadata;
	
}





__END__

=head1 NAME

get_feature_info.pl

A script to collect feature information from a BioPerl SeqFeature::Store db.

=head1 SYNOPSIS

get_feature_info.pl <filename> 

  Options:
  --in <filename> 
  --db <name>
  --attrib <attribute1,attribute2,...>
  --type <primary_tag>
  --out <filename>
  --gz
  --version
  --help
  
  Attributes include:
   Chromosome
   Start
   Stop
   Strand
   Score
   Name
   Alias
   Note
   Type
   Primary_tag
   Source
   Length
   Midpoint
   Phase
   RNA_count
   Exon_count
   Transcript_length (sum of exon lengths)
   Parent (name)
   Primary_ID
   <tag>

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a previously generated feature dataset.
It should be in the tim data format that is generated by this
program and others, although other tab-delimited text data
formats may be usable. See the file description in C<Bio::ToolBox::file_helper>.

=item --db <name>

Specify the name of the BioPerl gff database to use as source. This is required 
for new feature data files. For pre-existing input data files, this argument 
is optional, but if given it overrides the database listed in the file; this 
is useful for collecting data from multiple databases.

=item --attrib <attribute>

Specify the attribute to collect for each feature. Standard GFF attributes 
may be collected, as well as values from specific group tags. These tags 
are found in the group (ninth) column of the source GFF file. Standard 
attributes include the following
   
   - Chromosome
   - Start
   - Stop
   - Strand
   - Score
   - Name
   - Alias
   - Note
   - Type
   - Primary_tag
   - Source
   - Length
   - Midpoint
   - Phase
   - RNA_count (number of RNA subfeatures)
   - Exon_count (number of exons, or CDS, subfeatures)
   - Transcript_length
   - Parent (name)
   - Primary_ID
   - <tag>

If attrib is not specified on the command line, then an interactive list 
will be presented to the user for selection. Especially useful when you 
can't remember the feature's tag keys in the database.
   
=item --type <primary_tag>

When the input file does not have a type column, a type or primary_tag 
may be provided. This is especially useful to restrict the database 
search when there are multiple features with the same name.

=item --out <filename>

Optionally specify an alternate output file name. The default is to 
overwrite the input file.

=item --gz

Indicate whether the output file should (not) be compressed by gzip. 
If compressed, the extension '.gz' is appended to the filename. If a compressed 
file is opened, the compression status is preserved unless specified otherwise.

=item --version

Print the version number.

=item --help

Display this help. 

=back

=head1 DESCRIPTION

This program will collect attributes for a list of features from the database. 
The attributes may be general attributes, such as chromsome, start, stop, 
strand, etc., or feature specific attributes stored in the original group 
field of the original source GFF file.

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