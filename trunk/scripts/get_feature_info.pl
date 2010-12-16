#!/usr/bin/perl

# A script to pull out feature attributes from a database

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	find_column_index
);
use tim_db_helper qw(
	open_db_connection
);
use tim_file_helper qw(
	load_tim_data_file
	write_tim_data_file
);

print "\n This script will get additional information about features\n\n";



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
	$gz,
	$help
); 

# Command line options
GetOptions( 
	'in=s'     => \$infile, # input file
	'attrib=s' => \$attrib_request, # attribute
	'out=s'    => \$outfile, # output filename
	'db=s'     => \$database, # database name
	'gz!'      => \$gz, # gzip status
	'help'     => \$help, # help
);


# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
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
my $name_index = find_column_index($main_data_ref, 'name');
my $type_index = find_column_index($main_data_ref, 'type');
unless (
	defined $name_index and
	defined $type_index 
) {
	die 'unable to identify Name and/or Type columns in data table';
}




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
my $file_success = write_tim_data_file( {
	'data'      => $main_data_ref,
	'filename'  => $outfile,
} );
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
	
	my @list;
	
	# provided by command line argument
	if ($attrib_request) {
		
		# check if there are multiple items in the list
		if ($attrib_request =~ /,/) {
			@list = split /,/, $attrib_request;
		}
		else {
			push @list, $attrib_request;
		}
	}
	
	# request interactively from user
	else {
		
		# get the list of features to check for examples
		
		# get the attributes for a sample of features
		# store the tag keys in an example hash
		my %tagexamples;
		for (my $i = 1; $i < 50; $i++) {
			my @examples = $db->features(
				-name     => $main_data_ref->{'data_table'}->[$i][$name_index],
				-type     => $main_data_ref->{'data_table'}->[$i][$type_index]
			);
			unless (@examples) {
				next;
			}
			my %taghash = $examples[0]->attributes();
			foreach (keys %taghash) {
				$tagexamples{$_} += 1;
			}
		}
		
		# present list to user
		print " These are the attributes which may be collected:\n";
		my $i = 1;
		my %index2att;
		# standard attributes for any user
		foreach ( qw(chromo start stop length midpoint strand phase score) ) {
			print "   $i\t$_\n";
			$index2att{$i} = $_;
			$i++;
		}
		# specific attributes for these features
		foreach (sort {$a cmp $b} keys %tagexamples) {
			print "   $i\t$_\n";
			$index2att{$i} = $_;
			$i++;
		}
		
		# collect the user response
		print " Enter the attribute number(s) to collect, comma delimited   ";
		my $answer = <STDIN>;
		chomp $answer;
		$answer =~ s/\s//g;
		foreach (split /,/, $answer) {
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
	if ($attrib eq 'chromo') {
		$method = \&get_chromo;
	} 
	elsif ($attrib eq 'start') {
		$method = \&get_start;
	} 
	elsif ($attrib eq 'stop') {
		$method = \&get_stop;
	} 
	elsif ($attrib eq 'length') {
		$method = \&get_length;
	} 
	elsif ($attrib eq 'midpoint') {
		$method = \&get_midpoint;
	} 
	elsif ($attrib eq 'strand') {
		$method = \&get_strand;
	} 
	elsif ($attrib eq 'phase') {
		$method = \&get_phase;
	} 
	elsif ($attrib eq 'score') {
		$method = \&get_score;
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
		my @features = $db->features( 
			# define the region
			-name   => $table->[$row][$name_index],
			-type  => $table->[$row][$type_index],
		);
		if (scalar @features == 0) {
			warn " no features found for '$table->[$row][$name_index]'\n";
			
			# record null value(s)
			foreach (@list) {
				push @{ $table->[$row] }, '.'; 
			}
			next;
		}
		elsif (scalar @features > 1) {
			warn " multiple features found for '$table->[$row][$name_index]'" .
				". Using first one\n";
		}
		
		# get the attribute(s)
		for (my $i = 0; $i < scalar @list; $i++) {
 			# for each request in the list, we will collect the attribute
 			# we'll use the method sub defined in the methods
 			# pass both the feature and the name of the attribute
 			# only the tag value actually needs the name of the attribute
 			push @{ $table->[$row] }, 
 				&{ $methods[$i] }($features[0], $list[$i]);
 		}
	}
	
	# record the metadata
	foreach (@list) {
		record_metadata($_);
	}
}




sub get_chromo {
	my $feature = shift;
	return $feature->seq_id || '.';
}


sub get_start {
	my $feature = shift;
	return $feature->start || '.';
}



sub get_stop {
	my $feature = shift;
	return $feature->end || '.';
}


sub get_length {
	my $feature = shift;
	return $feature->length || '.';
}


sub get_midpoint {
	my $feature = shift;
	return ($feature->start + int( $feature->length / 2) ) || '.';
}


sub get_strand {
	my $feature = shift;
	return $feature->strand || '.';
}


sub get_phase {
	my $feature = shift;
	return $feature->phase || '.';
}


sub get_score {
	my $feature = shift;
	return $feature->score || '.';
}


sub get_tag_value {
	my $feature = shift;
	my $attrib = shift;
	my @values = $feature->get_tag_values($attrib);
	return $values[0] || '.';
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

A script to feature information from a Bioperl SeqFeature::Store db.

=head1 SYNOPSIS

get_feature_info.pl <filename> 

  --in <filename> 
  --attrib <attribute1,attribute2,...>
  --db <name>
  --out filename
  --(no)gz
  --help

Attributes include:
   chromo
   start
   stop
   length
   midpoint
   strand
   phase
   score
   <tag>


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --in <filename>

Specify the file name of a previously generated feature dataset.
It should be in the tim data format that is generated by this
program and others, although other tab-delimited text data
formats may be usable. See the file description in C<tim_db_helper.pm>.

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
   
   -chromo
   -start
   -stop
   -length
   -midpoint
   -strand
   -phase
   -score

If attrib is not specified on the command line, then an interactive list 
will be presented to the user for selection. Especially useful when you 
can't remember the feature's tag keys in the database.
   
=item --out <filename>

Optionally specify an alternate output file name. The default is to 
overwrite the input file.

=item --(no)gz

Indicate whether the output file should (not) be compressed by gzip. 
If compressed, the extension '.gz' is appended to the filename. If a compressed 
file is opened, the compression status is preserved unless specified otherwise.

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


=head1 TODO

Finish the coding!

