#!/usr/bin/perl

# a script to convert my data GFF v.2 files to GFF3 files

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper qw(
	open_tim_data_file
	open_to_write_fh
);
my $VERSION = '1.5.7';


print "\n This script will convert my data GFF v.2 files to GFF3\n";


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
	$type,
	$source,
	$mito,
	$help,
	$print_version,
);


# Command line options
GetOptions( 
	'type=s'    => \$type, # the new gff type
	'source=s'  => \$source, # the new source
	'mt=s'      => \$mito, # new chr17 name
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
	print " Biotoolbox script my_gff2gff3.pl, version $VERSION\n\n";
	exit;
}


### Check for required values
unless (@ARGV) {
	die " No input files defined!\n";
}

### Convert the files
foreach my $infile (@ARGV) {
	
	print " Converting file '$infile'... ";
	
	# Open Input file
	my ($in_fh, $metadata_ref) = open_tim_data_file($infile) or 
		warn "\n unable to open input file!\n";
	unless ($in_fh) {next;}
	unless ($metadata_ref->{gff}) {
		warn "\n input file must be a gff file!\n";
		next;
	}
	if ($metadata_ref->{gff} == 3) {
		warn "\n the input file is already version 3!\n";
		next;
	}
	
	
	# Prepare Output file
	my $outfile = $infile;
	$outfile =~ s/\.gff/.gff3/;
	my $out_fh = open_to_write_fh($outfile) or 
		warn "\n unable to open output file!\n";
	unless ($out_fh) {next;}
	
	# Write output headers
	{
		print {$out_fh} "##gff-version 3\n";
		if ($metadata_ref->{'program'}) {
			# write program header if present
			print {$out_fh} '# Program ' . $metadata_ref->{'program'} . "\n";
		}
		if ($metadata_ref->{'db'}) {
			# write database header if present
			print {$out_fh} '# Database ' . $metadata_ref->{'db'} . "\n";
		}
		if ($metadata_ref->{'feature'}) {
			# write feature header if present
			print {$out_fh} '# Feature ' . $metadata_ref->{'feature'} . "\n";
		}
		
		# Write the miscellaneous headers
		foreach ( @{ $metadata_ref->{'other'} } ) {
			# write remaining miscellaneous header lines if present
			print {$out_fh} $_;
		}
		
		# Write the column metadata headers
		for (my $i = 0; $i < $metadata_ref->{'number_columns'}; $i++) {
			if (scalar( keys %{ $metadata_ref->{$i} } ) > 2) {
				# more than two keys are present
				# we will put each key=value pair into @pairs, listed asciibetically
				my @pairs; # an array of the key value pairs from the metadata hash
				foreach (sort {$a cmp $b} keys %{ $metadata_ref->{$i} } ) {
					push @pairs,  $_ . '=' . $metadata_ref->{$i}{$_};
				}
				print {$out_fh} "# Column_", $i+1, " ", join(";", @pairs), "\n";
			}
		}
	}
	
	
	# Do the conversion
	my %id_lookup;
	my $id_count = 0;
	while (my $line = $in_fh->getline) {
		chomp $line;
		my @data = split /\t/, $line;
		
		
		# update group
		my @groupdata = split /\s*;\s*/, $data[8];
		my ($class, $id) = split / /, shift @groupdata, 2;
			# the GFF2 groupd field usually has a "class id" as the first 
			# element in the group field
			# we'll ignore the class, since it's usually redundant with the
			# method
			# but we want the id
		$id =~ s/"//g; # strip the quotation marks
		$id =~ s/ /_/g; # replace any spaces with underscores, why are these 
						# these here anyway!!!???
		$id =~ s/,/./g; # replace any commas with periods, why are these 
						# these here anyway!!!???
		my $newid = $id;
		if ($class eq 'Experiment') {
			# this is a generic class that I used for microarray experiments
			# I'm using this if test as a shortcut to speed up processing
			# It is highly likely that all the features have the identical
			# id, and so it is very expensive to check for id uniqueness
			# this shortcut simply appends a unique value on the end of the 
			# of the id
			$id_count++;
			$newid .= '.' . $id_count;
			if (exists $id_lookup{$newid}) {
				warn " uh oh, we're in trouble, $newid exists!\n";
				next;
			}
		}
		else {
			# expensive search looking for a unique id
			my $i = 0;
			while (exists $id_lookup{$newid}) {
				# look up based first on current type, then id
				# goal is to make the id unique
				# this will become very expensive with a big data file with all 
				# unique ids, as each one will have to loop through looking for 
				# the unique one
				$i++;
				$newid = $id . '.' . $i;
			}
		}
		$id_lookup{$newid} = $id; # remember this new id
		my $group = "ID=$newid;Name=$id";
			# the ID must be a unique value, the name does not have to be
		foreach (@groupdata) {
			# there may be additional tags remaining in the group field
			# simply convert these into key=value tags
			if (/^(\w+) (.+)$/) {
				my $key = $1;
				my $value = $2;
				$value =~ s/"//g; # strip quotation marks
				$group .= ";$key=$value";
			}
		}
			
		
		
		# update
		if ($type) {
			# only if a new type was defined
			$data[2] = $type;
		}
		if ($source) {
			# only if a new source was defined
			$data[1] = $source;
		}
		$data[8] = $group;
		
		
		# fix chr17
		if ($data[0] eq 'chr17' and $mito) {
			$data[0] = $mito;
		}
		
		
		# print
		print {$out_fh} join("\t", @data) . "\n";
	}
	
	
	# Finish
	$in_fh->close;
	$out_fh->close;
	print "Done.\n";
}	



__END__

=head1 NAME

my_gff2gff3.pl

A script to convert my data GFF v.2 files to GFF3

=head1 SYNOPSIS

my_gff2gff3.pl [--options ...] <file1.gff> <file2.gff> ...
  
  --type <text>
  --source <text>
  --mt <text>
  --help


=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --type <text>

Provide a text string to bue used in the method/class/type
column. Since we are generating a v.3 GFF file, the type 
value should map to a valid Sequence Ontology Term. Common 
data terms include 'microarray_oligo', 'tag', or 'STS'. See
http://www.sequenceontology.org for more information. 
Default is to use the original.

=item --source <text>

Optionally provide a new text string to replace the 'source' 
field. Default is to use the original.

=item --mt <text>

I used to rename the mitochrondrial chromosome in cerevisiae to 
'chr17' to help facilitate analysis. This was not very 
professional. A new name may be provided for renaming all 'chr17' 
features (back) to the new (original) name.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will convert a simple data version 2 GFF formatted file into a GFF3 
file. Typical GFF data files that I had generated are in GFF v.2 format. 
Additionally, I used the 'type' or 'method' column (the third column) as 
a means of identifying and grouping values into a single data set. These 
were further often identified by the source tag of 'data'. Finally, the 
group field (ninth column) was simply in the format 'Experiment $type'.

With bioperl, GBrowse, and other programs moving towards format v.3 
(GFF3), these practices are no longer acceptable. This program attempts 
to convert my GFF v.2 data files into valid GFF3 files. First, the type 
column must be a valid Sequence Ontology term. For microarray data, this 
can be 'microarray_oligo'. For ChIP-seq data, the best terms may either be 
'tag' or 'STS' (for Sequence Tag Site). Next, the group field should 
contain both Name and ID tags. The old 'type' field will become the new 
Name tag and will be used for grouping. The ID tag must be made unique 
(accomplished by appending a unique number to the Name). Finally, the 
source field could be better utilized with more unique information.

The program will accept and process one or more gff input files. It will 
write a new file, renaming the extension to '.gff3' for each input file. 

NOTE that this program will not deal with parent-child relationships. It 
is strictly designed to work with very simple one-line features, such as 
representing microarray data.



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







