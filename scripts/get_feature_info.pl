#!/usr/bin/perl

# documentation at end of file

use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use Bio::ToolBox::db_helper qw(
	open_db_connection
);
use Bio::ToolBox::Data;
use Bio::ToolBox::utility;
my $VERSION = '1.67';

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
	'i|in=s'     => \$infile, # input file
	'a|attrib=s' => \$attrib_request, # attribute
	'o|out=s'    => \$outfile, # output filename
	'd|db=s'     => \$database, # database name
	't|type=s'   => \$use_type, # force a type
	'z|gz!'      => \$gz, # gzip status
	'h|help'     => \$help, # help
	'v|version'  => \$print_version, # print the version
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
	print " Biotoolbox script get_feature_info.pl, version $VERSION\n";
	eval {
		require Bio::ToolBox;
		my $v = Bio::ToolBox->VERSION;
		print " Biotoolbox package version $v\n";
	};
	exit;
}



# Check for required values
unless ($infile) {
	$infile = shift @ARGV;
}
unless (defined $gz) {$gz = 0}




### Load the feature list

# load file
my $Data = Bio::ToolBox::Data->new(file => $infile) or
	die " Unable to load data file!\n";
printf " Loaded %s features from $infile.\n", format_with_commas( $Data->last_row );

if ($use_type) {
	# set the general feature type, which will be used as a proxy for individual 
	# feature types in the Data method that looks for the database features.
	my $f = $Data->feature;
	if (defined $f) {
		print " Resetting feature type '$f' to requested '$use_type'\n";
	}
	$Data->feature($use_type);
}

# check columns
{
	my $name_i = $Data->name_column;
	my $id_i   = $Data->id_column;
	my $type_i = $Data->type_column;
	unless (defined $name_i or defined $id_i) {
		die " unable to identify a name or ID column. Database lookup unlikely to work\n";
	}
	unless (defined $type_i or $use_type) {
		warn " No type column present or type specified. Database lookup may have problems\n";
	}
}


### Establish database connection
my $db;
if ($database) {
	if ($Data->database and $Data->database ne $database) {
		print " Using provided database '$database' instead of metadata-specified database\n";
	}
	
	# check if we have a real file
	if (-e $database) {
		if ($database =~ /(?:db|database|sqlite)$/i) {
			# looks like maybe a SQLite file
			$Data->database($database);
			$db = $Data->open_meta_database or die "unable to open database!\n";
		}
		elsif ($Data->taste_file($database)) {
			# looks like a parsable annotation file
			# parse the file and associate with the table
			$Data->parse_table( {
				file     => $database,
				feature  => $Data->feature, # may be explicitly set by use_type above
				simplify => 0,
				subfeature => 'exon,cds,utr,codon',
			} ) or die " unable to parse annotation file '$database'!\n";
			# this may fail if not everything can be loaded
		}
		else {
			die " unrecognized database file!\n";
		}
	}
	else {
		# maybe the name of a database
		# attempt to open
		$Data->database($database);
		$db = $Data->open_meta_database or die "unable to open database!\n";
	}
}
elsif ($Data->database) {
	# use metadata-specified database
	if ($Data->database =~ /^Parsed:(.+)$/) {
		# there was a parsed annotation file, try and re-parse it
		$Data->parse_table( {
			file     => $1,
			feature  => $Data->feature, # may be explicitly set by use_type above
			simplify => 0,
			subfeature => 'exon,cds,utr,codon',
		} ) or die " unable to parse annotation file '$1' in metadata!\n try specifying a new database or annoation file\n";
		# this may fail if not everything can be loaded
	}
	else {
		$db = $Data->open_meta_database or die "unable to open database!\n";
	}
}
else {
	die "No database defined! See help\n";
}





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
my $file_success = $Data->write_file(
	'filename'  => $outfile,
	'gz'        => $gz,
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
		print " Collecting sample features to generate list of attributes....\n";
		
		# get the attributes for a sample of features
		# store the tag keys in an example hash
		my %tagexamples;
		my $stream = $Data->row_stream;
		for (1 .. 50) {
			my $row = $stream->next_row;
			last unless $row;
			my $f = $row->feature;
			next unless $f;
			my %taghash = $f->attributes();
			foreach (keys %taghash) {
				$tagexamples{$_} += 1;
			}
		}
		
		# present list to user
		print " These are the attributes which may be collected:\n";
		my $i = 1;
		my %index2att;
		# standard attributes for any user
		foreach ( 
			qw(Chromosome Start Stop Strand Score Name Alias Note Type Primary_tag Source 
				Length Midpoint Phase RNA_count Exon_count Gene_length Transcript_length 
				Parent Primary_ID 
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
	elsif ($attrib =~ /^gene.?length$/i) {
		$method = \&get_gene_length;
	} 
	elsif ($attrib =~ /^transcript.?length$/i) {
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
	elsif ($attrib =~ /^rna.?count$/i) {
		$method = \&get_rna_number;
	} 
	elsif ($attrib =~ /^exon.?count$/i) {
		$method = \&get_exon_number;
	} 
	elsif ($attrib =~ /^parent$/i) {
		$method = \&get_parent;
	} 
	elsif ($attrib =~ /^primary.?id$/i) {
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
	elsif ($attrib =~ /^primary.?tag$/i) {
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
	
	# prepare columns
	my @indices;
	foreach (@list) {
		my $i = $Data->add_column($_);
		push @indices, $i;
	}
	
	print " Retrieving ", join(", ", @list), "\n";
	my $stream = $Data->row_stream;
	while (my $row = $stream->next_row) {
		
		my $feature = $row->feature(1);
			# pass a true value to force the seqfeature lookup
		if ($feature) {
			# get the attribute(s)
			for (my $i = 0; $i < scalar @list; $i++) {
				# for each request in the list, we will collect the attribute
				# we'll use the method sub defined in the methods
				# pass both the feature and the name of the attribute
				# only the tag value actually needs the name of the attribute
				my $v = &{ $methods[$i] }($feature, $list[$i]);
				$row->value($indices[$i], $v);
			}
		}
		else {
			# no feature found, cannot collect attributes
			# record nulls
			foreach (@indices) { $row->value($_, '.') }
		}
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


sub get_gene_length {
	my $feature = shift;
	return 0 unless $feature;
	
	# collect all exons or CDSs for every transcript
	my @exons;
	my @cdss;
	foreach my $subfeat ( $feature->get_SeqFeatures() ) {
		# feature may consist of multiple transcripts
		if ($subfeat->primary_tag =~ /exon/i) {
			push @exons, $subfeat;
		}
		elsif ($subfeat->primary_tag =~ /utr|untranslated/i) {
			push @cdss, $subfeat;
		}
		elsif ($subfeat->primary_tag =~ /cds/i) {
			push @cdss, $subfeat;
		}
		elsif ($subfeat->primary_tag =~ /rna|transcript/i) {
			# an RNA subfeature, keep going down another level
			foreach my $f ($subfeat->get_SeqFeatures) {
				if ($f->primary_tag =~ /exon/i) {
					push @exons, $f;
				}
				elsif ($f->primary_tag =~ /utr|untranslated/i) {
					push @cdss, $f;
				}
				elsif ($f->primary_tag =~ /cds/i) {
					push @cdss, $f;
				}
			}
		}
	}
	
	# Determine which subfeatures to collect
	# we prefer to use exons because they are easier
	# if exons are not defined then we'll infer them from CDSs and UTRs
	my @features_to_check;
	if (@exons) {
		@features_to_check = @exons;
	}
	elsif (@cdss) {
		@features_to_check = @cdss;
	}
	else {
		# found neither exons, CDSs, or UTRs
		# possibly because there were no subfeatures
		# in this case we just take the whole thing
		push @features_to_check, $feature;
	}
	
	
	# Order and collapse the subfeatures
	# we don't want overlapping features
	my @sorted_features =   map {$_->[1]} 
							sort {$a->[0] <=> $b->[0]} 
							map { [$_->start, $_] } 
							@features_to_check;
	my @merged;
	my $current = shift @sorted_features;
	while ($current) {
		if (@sorted_features) {
			my $next = shift @sorted_features;
			if ($current->overlaps($next) ) {
				# overlapping features, so reassign the end point to that of the next
				$current->end($next->end);
			}
			else {
				# no overlap, keep current, next becomes current
				push @merged, $current;
				$current = $next;
			}
		}
		else {
			# no more features
			push @merged, $current;
			undef $current;
		}
	}
	
	# calculate the length
	my $gene_length = 0;
	foreach my $f (@merged) {
		$gene_length += $f->length;
	}
	
	# return most appropriate number
	return $gene_length;
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
		if ($f->primary_tag =~ /rna|transcript/i) {
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
		elsif ($f->primary_tag =~ /rna|transcript/i) {
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
	# this is tricky, because it's not always recorded in the SeqFeature object itself
	# and we can't easily recurse in memory an object and what it is linked to 
	# 
	if ($feature->has_tag('parent_id')) {
		# feature has a parent in a database
		my @parent_ids = $feature->get_tag_values('parent_id');
		if ($db) {
			# looks like we have SeqFeature database
			# unfortunately we have to do a slow lookup through the database 
			# using the attribute load_id
			# this seems to be the only way to get to the parent name
			my @parents;
			foreach my $id (@parent_ids) {
				foreach my $p ($db->get_features_by_attribute('load_id' => $id)) {
					push @parents, $p->display_name;
				}
			}
			return join(',', @parents);
		}
	}
	if ($feature->has_tag('Parent')) {
		my @parents = $feature->get_tag_values('Parent');
		return join(',', @parents);
	}
	return '.';
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





__END__

=head1 NAME

get_feature_info.pl

A program to collect feature information from a BioPerl SeqFeature::Store db.

=head1 SYNOPSIS

get_feature_info.pl <filename> 

  File options:
  -i --in <filename>                        input file of list db features
  -o --out <filename>                       optional output file
  
  Database options:
  -d --db <name>                            annotation database: mysql sqlite
                                              or annotation file: gtf gff ucsc
  -a --attrib <attribute1,attribute2,...>   list of attributes to collect
  -t --type <primary_tag>                   specify a feature type as needed
  
  General options:
  -z --gz                                   compress output file
  -v --version                              print version and exit
  -h --help                                 show extended documentation
  
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
   RNA_count (number of transcript subfeatures)
   Exon_count (number of exon subfeatures)
   Gene_length (sum of all merged, collapsed, transcript exon lengths)
   Transcript_length (sum of exon lengths)
   Parent (name)
   Primary_ID
   <tag>

=head1 OPTIONS

The command line flags and descriptions:

=head2 File options

=over 4

=item --in E<lt>filenameE<gt>

Specify an input file containing either a list of database features or 
genomic coordinates for which to collect data. The file should be a 
tab-delimited text file, one row per feature, with columns representing 
feature identifiers, attributes, coordinates, and/or data values. The 
first row should be column headers. Text files generated by other 
B<BioToolBox> scripts are acceptable. Files may be gzipped compressed.

=item --out E<lt>filenameE<gt>

Optionally specify an alternate output file name. The default is to 
overwrite the input file.

=back

=head2 Database options

=over 4

=item --db E<lt>nameE<gt>

Specify the name or SQLite file of a L<Bio::DB::SeqFeature::Store> 
annotation database from which the information may be derived. This may 
be stored in the metadata comments of the input file, or an alternative 
file may be provided. 

Alternatively, specify an annotation file, e.g. GTF, GFF3, or UCSC gene 
table, that may be parsed into memory. 

The input file should include Name or ID columns that match features in the 
provided database. If a Type column is not present, then a type should be 
provided with the C<--type> option. Note that mixing and matching files and  
databases may not always work as well as intended.

=item --attrib E<lt>attributeE<gt>

Specify the attribute to collect for each feature. Standard GFF attributes 
may be collected, as well as values from specific group tags. These tags 
are found in the group (ninth) column of the source GFF file. Standard 
attributes include the following

=over 4
 
=item * Chromosome

=item * Start

=item * Stop

=item * Strand

=item * Score

=item * Name

=item * Alias

=item * Note

=item * Type

=item * Primary_tag

=item * Source

=item * Length

=item * Midpoint

=item * Phase

=item * RNA_count (number of transcript subfeatures)

=item * Exon_count (number of exon subfeatures)

=item * Gene_length (sum of all merged, collapsed, transcript exon lengths)

=item * Transcript_length (sum of exon lengths)

=item * Parent (name)

=item * Primary_ID

=item * <tag>

=back

If attrib is not specified on the command line, then an interactive list 
will be presented to the user for selection. Especially useful when you 
can't remember the feature's tag keys in the database.
   
=item --type E<lt>primary_tagE<gt>

When the input file does not have a type column, a type or primary_tag 
may be provided. This is especially useful to restrict the database 
search when there are multiple features with the same name.

=back

=head2 General options

=over 4

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
it under the terms of the Artistic License 2.0.  
