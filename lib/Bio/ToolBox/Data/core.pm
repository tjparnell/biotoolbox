package Bio::ToolBox::Data::core;
our $VERSION = 1.27;

=head1 NAME

Bio::ToolBox::Data::core - Common functions to Bio:ToolBox::Data family

=head1 DESCRIPTION

Common methods for metadata and manipulation in a L<Bio::ToolBox::Data> 
data table and L<Bio::ToolBox::Data::Stream> file stream. This module 
should not be used directly. See the respective modules for more information.

=cut

use strict;
use Carp qw(carp cluck croak confess);
use Bio::ToolBox;
use base 'Bio::ToolBox::Data::file';
use Bio::ToolBox::db_helper qw(
	open_db_connection
	verify_or_request_feature_types
);

1;

#### Initialization and verification ###############################################

sub new {
	my $class = shift;
	
	# Initialize the hash structure
	my $ToolBoxVersion = Bio::ToolBox->VERSION;
	my %data = (
		'program'        => "$0, lib v$ToolBoxVersion",
		'feature'        => undef,
		'feature_type'   => undef,
		'db'             => undef,
		'gff'            => 0,
		'bed'            => 0,
		'ucsc'           => 0,
		'number_columns' => 0,
		'last_row'       => 0,
		'headers'        => 1,
		'column_names'   => [],
		'filename'       => undef,
		'basename'       => undef,
		'extension'      => undef,
		'path'           => undef,
		'comments'       => [],
		'data_table'     => [],
	);
	
	# Finished
	return bless \%data, $class;
}


sub verify {
	# this function does not rely on any self functions for two reasons
	# this is a low level integrity checker
	# this is very old code from before the days of an OO API of Bio::ToolBox
	my $self = shift;
	
	# check for data table
	unless (
		defined $self->{'data_table'} and 
		ref $self->{'data_table'} eq 'ARRAY'
	) {
		cluck " No data table in passed data structure!";
		return;
	}
	
	# check for last row index
	if (defined $self->{'last_row'}) {
		my $number = scalar( @{ $self->{'data_table'} } ) - 1;
		if ($self->{'last_row'} != $number) {
			cluck " data table last_row index [$number] doesn't match " . 
				"metadata value [" . $self->{'last_row'} . "]!\n";
			# fix it for them
			$self->{'last_row'} = $number;
		}
	}
	else {
		# define it for them
		$self->{'last_row'} = 
			scalar( @{ $self->{'data_table'} } ) - 1;
	}
	
	# check for consistent number of columns
	if (defined $self->{'number_columns'}) {
		my $number = $self->{'number_columns'};
		my @problems;
		my $too_low = 0;
		my $too_high = 0;
		for (my $row = 0; $row <= $self->{'last_row'}; $row++) {
			my $count = scalar @{ $self->{'data_table'}->[$row] };
			if ($count != $number) {
				push @problems, $row;
				$too_low++ if $count < $number;
				$too_high++ if $count > $number;
				while ($count < $number) {
					# we can sort-of-fix this problem
					$self->{'data_table'}->[$row][$count] = '.';
					$count++;
				}
			}
		}
		if ($too_low) {
			cluck " $too_low rows in data table had fewer than expected columns!\n" . 
				 "  padded rows " . join(',', @problems) . " with null values\n";
		}
		if ($too_high) {
			cluck " $too_high rows in data table had more columns than expected!\n" . 
				" rows " . join(',', @problems) . "\n";
			return;
		}
	}
	else {
		$self->{'number_columns'} = 
			scalar @{ $self->{'data_table'}->[0] };
	}
	
	# check metadata
	for (my $i = 0; $i < $self->{'number_columns'}; $i++) {
		unless (
			$self->{$i}{'name'} eq 
			$self->{'data_table'}->[0][$i]
		) {
			cluck " incorrect or missing metadata!  Column header names don't" .
				" match metadata name values for index $i!" . 
				" compare '" . $self->{$i}{'name'} . "' with '" .
				$self->{'data_table'}->[0][$i] . "'\n";
			return;
		}
	}
	
	# check for proper gff structure
	if ($self->{'gff'}) {
		# if any of these checks fail, we will reset the gff version to 
		# the default of 0, or no gff
		my $gff_check = 1; # start with assumption it is true
		
		# check number of columns
		if ($self->{'number_columns'} != 9) {
			$gff_check = 0;
		}
		
		# check column indices
		if (
			exists $self->{0} and
			$self->{0}{'name'} !~ 
			m/^#?(?:chr|chromo|seq|refseq|ref_seq|seq|seq_id)/i
		) {
			$gff_check = 0;
		}
		if (
			exists $self->{1} and
			$self->{1}{'name'} !~ m/^source/i
		) {
			$gff_check = 0;
		}
		if (
			exists $self->{2} and
			$self->{2}{'name'} !~ m/^type|method/i
		) {
			$gff_check = 0;
		}
		if (
			exists $self->{3} and
			$self->{3}{'name'} !~ m/^start/i
		) {
			$gff_check = 0;
		}
		if (
			exists $self->{4} and
			$self->{4}{'name'} !~ m/^stop|end/i
		) {
			$gff_check = 0;
		}
		if (
			exists $self->{5} and
			$self->{5}{'name'} !~ m/^score|value/i
		) {
			$gff_check = 0;
		}
		if (
			exists $self->{6} and
			$self->{6}{'name'} !~ m/^strand/i
		) {
			$gff_check = 0;
		}
		if (
			exists $self->{7} and
			$self->{7}{'name'} !~ m/^phase/i
		) {
			$gff_check = 0;
		}
		if (
			exists $self->{8} and
			$self->{8}{'name'} !~ m/^group|attribute/i
		) {
			$gff_check = 0;
		}
		
		# update gff value as necessary
		if ($gff_check == 0) {
			# reset metadata
			$self->{'gff'} = 0;
			$self->{'headers'} = 1;
			
			# remove the AUTO key from the metadata
			for (my $i = 0; $i < $self->{'number_columns'}; $i++) {
				if (exists $self->{$i}{'AUTO'}) {
					delete $self->{$i}{'AUTO'};
				}
			}
		}
	}
	
	# check for proper BED structure
	if ($self->{'bed'}) {
		# if any of these checks fail, we will reset the bed flag to 0
		# to make it not a bed file format
		my $bed_check = 1; # start with assumption it is correct
		
		# check number of columns
		if (
			$self->{'number_columns'} < 3 and 
			$self->{'number_columns'} > 12 
		) {
			$bed_check = 0;
		}
		
		# check column index names
		if (
			exists $self->{0} and
			$self->{0}{'name'} !~ 
			m/^#?(?:chr|chromo|seq|refseq|ref_seq|seq|seq_id)/i
		) {
			$bed_check = 0;
		}
		if (
			exists $self->{1} and
			$self->{1}{'name'} !~ m/^start/i
		) {
			$bed_check = 0;
		}
		if (
			exists $self->{2} and
			$self->{2}{'name'} !~ m/^stop|end/i
		) {
			$bed_check = 0;
		}
		
		# the remaining columns are tricky, as they may or may not be 
		# named as I expect, especially if it was generated de novo
		# so only check these if the original file extension was bed
		if (
			exists $self->{'extension'} and 
			$self->{'extension'} =~ /bed|bdg/i
		) {
			if (
				exists $self->{3} and
				$self->{3}{'name'} !~ m/^name|id|score/i
				# for bed this should be name or ID
				# for bedgraph this should be score
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{4} and
				$self->{4}{'name'} !~ m/^score|value/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{5} and
				$self->{5}{'name'} !~ m/^strand/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{6} and
				$self->{6}{'name'} !~ m/^thickstart/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{7} and
				$self->{7}{'name'} !~ m/^thickend/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{8} and
				$self->{8}{'name'} !~ m/^itemrgb/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{9} and
				$self->{9}{'name'} !~ m/^blockcount/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{10} and
				$self->{10}{'name'} !~ m/^blocksizes/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{11} and
				$self->{11}{'name'} !~ m/^blockstarts/i
			) {
				$bed_check = 0;
			}
		} 
		elsif (
			exists $self->{'extension'} and 
			$self->{'extension'} =~ /peak/i
		) {
			# some sort of peak file
			# narrowpeak: signalValue pValue qValue peak
			# broadpeak: signalValue pValue qValue
			if (
				exists $self->{3} and
				$self->{3}{'name'} !~ m/^name|id|score/i
				# for bed this should be name or ID
				# for bedgraph this should be score
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{4} and
				$self->{4}{'name'} !~ m/^score|value/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{5} and
				$self->{5}{'name'} !~ m/^strand/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{6} and
				$self->{6}{'name'} !~ m/^signalvalue/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{7} and
				$self->{7}{'name'} !~ m/^pvalue/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{8} and
				$self->{8}{'name'} !~ m/^qvalue/i
			) {
				$bed_check = 0;
			}
			if (
				exists $self->{9} and
				$self->{9}{'name'} !~ m/^peak/i
			) {
				$bed_check = 0;
			}
		}
		
		# reset the BED tag value as appropriate
		if ($bed_check) {
			$self->{'bed'} = $self->{'number_columns'};
		}
		else {
			# reset metadata
			$self->{'bed'} = 0;
			$self->{'headers'} = 1;
			
			# remove the AUTO key from the metadata
			for (my $i = 0; $i < $self->{'number_columns'}; $i++) {
				if (exists $self->{$i}{'AUTO'}) {
					delete $self->{$i}{'AUTO'};
				}
			}
		}
	}
	
	# check refFlat or genePred gene structure
	if ($self->{'ucsc'}) {
		# if any of these checks fail, we will reset the extension
		my $ucsc_check = 1; # start with assumption it is correct
		
		if ($self->{'number_columns'} == 16) {
			my @names = qw(bin name chrom strand txStart txEnd cdsStart cdsEnd 
				exonCount exonStarts exonEnds score name2 cdsStartSt 
				cdsEndStat exonFrames);
			for my $i (0 .. 15) {
				unless ($self->{$i}{'name'} =~ /$names[$i]/i) {
					$ucsc_check = 0;
					last;
				}
			}
		}		
		elsif ($self->{'number_columns'} == 15) {
			my @names = qw(name chrom strand txStart txEnd cdsStart cdsEnd 
				exonCount exonStarts exonEnds score name2 cdsStartSt 
				cdsEndStat exonFrames);
			for my $i (0 .. 14) {
				unless ($self->{$i}{'name'} =~ /$names[$i]/i) {
					$ucsc_check = 0;
					last;
				}
			}
		}		
		elsif ($self->{'number_columns'} == 12) {
			my @names = qw(name chrom strand txStart txEnd cdsStart cdsEnd 
						exonCount exonStarts exonEnds proteinID alignID);
			for my $i (0 .. 11) {
				unless ($self->{$i}{'name'} =~ /$names[$i]/i) {
					$ucsc_check = 0;
					last;
				}
			}
		}		
		elsif ($self->{'number_columns'} == 11) {
			my @names = qw(geneName transcriptName chrom strand txStart txEnd 
						cdsStart cdsEnd exonCount exonStarts exonEnds);
			for my $i (0 .. 10) {
				unless ($self->{$i}{'name'} =~ /$names[$i]/i) {
					$ucsc_check = 0;
					last;
				}
			}
		}		
		elsif ($self->{'number_columns'} == 10) {
			my @names = qw(name chrom strand txStart txEnd cdsStart cdsEnd 
						exonCount exonStarts exonEnds);
			for my $i (0 .. 9) {
				unless ($self->{$i}{'name'} =~ /$names[$i]/i) {
					$ucsc_check = 0;
					last;
				}
			}
		}
		
		if ($ucsc_check == 0) {
			# failed the check
			my $ext = $self->{'extension'};
			$self->{'filename'} =~ s/$ext/.txt/;
			$self->{'extension'} = '.txt';
			$self->{'ucsc'} = 0;
			
			# remove the AUTO key
			for (my $i = 0; $i < $self->{'number_columns'}; $i++) {
				if (exists $self->{$i}{'AUTO'}) {
					delete $self->{$i}{'AUTO'};
				}
			}
		}	
	}
	
	# check proper SGR file structure
	if (defined $self->{'extension'} and 
		$self->{'extension'} =~ /sgr/i
	) {
		# there is no sgr field in the data structure
		# so we're just checking for the extension
		# we will change the extension as necessary if it doesn't conform
		if (
			$self->{'number_columns'} != 3 or
			$self->{0}{'name'} !~ /^chr|seq|ref/i or
			$self->{1}{'name'} !~ /^start|position/i
		) {
			# doesn't smell like a SGR file
			# change the extension so the write subroutine won't think it is
			# make it a text file
			$self->{'extension'} =~ s/sgr/txt/i;
			$self->{'filename'}  =~ s/sgr/txt/i;
			$self->{'headers'} = 1;
			
			# remove the AUTO key from the metadata
			for (my $i = 0; $i < $self->{'number_columns'}; $i++) {
				if (exists $self->{$i}{'AUTO'}) {
					delete $self->{$i}{'AUTO'};
				}
			}
		}
	}
	
	# if we haven't made it here yet, then there was a problem
	return 1;
}



#### Database methods ##############################################################

sub open_database {
	my $self = shift;
	my $force = shift || 0;
	return unless $self->{db};
	if (exists $self->{db_connection}) {
		return $self->{db_connection} unless $force;
	}
	my $db = open_db_connection($self->{db}, $force);
	return unless $db;
	$self->{db_connection} = $db;
	return $db;
}

sub verify_dataset {
	my $self = shift;
	my $dataset = shift;
	my $database = shift; # name or object?
	return unless $dataset;
	if (exists $self->{verfied_dataset}{$dataset}) {
		return $self->{verfied_dataset}{$dataset};
	}
	else {
		if ($dataset =~ /^(?:file|http|ftp)/) {
			# local or remote file already verified?
			$self->{verfied_dataset}{$dataset} = $dataset;
			return $dataset;
		}
		$database ||= $self->open_database;
		my ($verified) = verify_or_request_feature_types(
			# normally returns an array of verified features, we're only checking one
			db      => $database,
			feature => $dataset,
		);
		if ($verified) {
			$self->{verfied_dataset}{$dataset} = $verified;
			return $verified;
		}
	}
	return;
}



#### Column Manipulation ####

sub delete_column {
	my $self = shift;
	if (defined $self->{fh}) {
		# Stream file handle is opened
		cluck "Cannot modify columns when a Stream file handle is opened!";
		return;
	}
	
	my @deletion_list = sort {$a <=> $b} @_;
	my @retain_list; 
	for (my $i = 0; $i < $self->number_columns; $i++) {
		# compare each current index with the first one in the list of 
		# deleted indices. if it matches, delete. if not, keep
		if ( $i == $deletion_list[0] ) {
			# this particular index should be deleted
			shift @deletion_list;
		}
		else {
			# this particular index should be kept
			push @retain_list, $i;
		}
	}
	return $self->reorder_column(@retain_list);
}

sub reorder_column {
	my $self = shift;
	if (defined $self->{fh}) {
		# Stream file handle is opened
		cluck "Cannot modify columns when a Stream file handle is opened!";
		return;
	}
	
	# reorder data table
	my @order = @_;
	for (my $row = 0; $row <= $self->last_row; $row++) {
		my @old = $self->row_values($row);
		my @new = map { $old[$_] } @order;
		splice( @{ $self->{data_table} }, $row, 1, \@new);
	}
	
	# reorder metadata
	my %old_metadata;
	for (my $i = 0; $i < $self->number_columns; $i++) {
		# copy the metadata info hash into a temporary hash
		$old_metadata{$i} = $self->{$i};
		delete $self->{$i}; # delete original
	}
	for (my $i = 0; $i < scalar(@order); $i++) {
		# now copy back from the old_metadata into the main data hash
		# using the new index number in the @order array
		$self->{$i} = { %{ $old_metadata{ $order[$i] } } };
		# assign new index number
		$self->{$i}{'index'} = $i;
	}
	$self->{'number_columns'} = scalar @order;
	delete $self->{column_indices} if exists $self->{column_indices};
	return 1;
}



#### General Metadata ####

sub feature {
	my $self = shift;
	if (@_) {
		$self->{feature} = shift;
	}
	return $self->{feature};
}

sub feature_type {
	my $self = shift;
	if (defined $self->{feature_type}) {
		return $self->{feature_type};
	}
	my $feature_type;
	if (defined $self->chromo_column and defined $self->start_column) {
		$feature_type = 'coordinate';
	}
	elsif (defined $self->id_column or 
		( defined $self->type_column and defined $self->name_column ) or 
		( defined $self->feature and defined $self->name_column )
	) {
		$feature_type = 'named';
	}
	else {
		$feature_type = 'unknown';
	}
	$self->{feature_type} = $feature_type;
	return $feature_type;
}

sub program {
	my $self = shift;
	if (@_) {
		$self->{program} = shift;
	}
	return $self->{program};
}

sub database {
	my $self = shift;
	if (@_) {
		$self->{db} = shift;
		if (exists $self->{db_connection}) {
			my $db = open_db_connection($self->{db});
			$self->{db_connection} = $db if $db;
		}
	}
	return $self->{db};
}

sub gff {
	my $self = shift;
	if ($_[0] and $_[0] =~ /^[123]/) {
		$self->{gff} = $_[0];
	}
	return $self->{gff};
}

sub bed {
	my $self = shift;
	if ($_[0] and $_[0] =~ /^\d+/) {
		$self->{bed} = $_[0];
	}
	return $self->{bed};
}

sub number_columns {
	my $self = shift;
	if (@_) {
		carp "number_columns() is a read only method";
	}
	return $self->{number_columns};
}

sub last_row {
	my $self = shift;
	if (@_) {
		carp "last_row() is a read only method";
	}
	return $self->{last_row};
}

sub filename {
	my $self = shift;
	if (@_) {
		carp "filename() is a read only method. Use add_file_metadata().";
	}
	return $self->{filename};
}

sub basename {
	my $self = shift;
	if (@_) {
		carp "basename() is a read only method. Use add_file_metadata().";
	}
	return $self->{basename};
}

sub path {
	my $self = shift;
	if (@_) {
		carp "path() is a read only method. Use add_file_metadata().";
	}
	return $self->{path};
}

sub extension {
	my $self = shift;
	if (@_) {
		carp "extension() is a read only method. Use add_file_metadata().";
	}
	return $self->{extension};
}



#### General Comments ####

sub comments {
	my $self = shift;
	my @comments = @{ $self->{comments} };
	foreach (@comments) {s/[\r\n]+//g}
	# comments are not chomped when loading
	# side effect of dealing with rare commented header lines with null values at end
	return @comments;
}

sub add_comment {
	my $self = shift;
	my $comment = shift or return;
	# comment is not required to be prefixed with "# ", it will be added when saving
	push @{ $self->{comments} }, $comment;
	return 1;
}

sub delete_comment {
	my $self = shift;
	my $index = shift;
	if (defined $index) {
		eval {splice @{$self->{comments}}, $index, 1};
	}
	else {
		$self->{comments} = [];
	}
}



#### Column Metadata ####

sub list_columns {
	my $self = shift;
	my @list;
	for (my $i = 0; $i < $self->number_columns; $i++) {
		push @list, $self->{$i}{'name'};
	}
	return wantarray ? @list : \@list;
}

sub name {
	my $self = shift;
	my ($index, $new_name) = @_;
	return unless defined $index;
	return unless exists $self->{$index}{name};
	if (defined $new_name) {
		$self->{$index}{name} = $new_name;
		if (exists $self->{data_table}) {
			$self->{data_table}->[0][$index] = $new_name;
		}
		elsif (exists $self->{column_names}) {
			$self->{column_names}->[$index] = $new_name;
		}
	}
	return $self->{$index}{name};
}

sub metadata {
	my $self = shift;
	my ($index, $key, $value) = @_;
	return unless defined $index;
	return unless exists $self->{$index};
	if ($key and $key eq 'name') {
		return $self->name($index, $value);
	}
	if ($key and defined $value) { 
		# we are setting a new value
		$self->{$index}{$key} = $value;
		return $value;
	}
	elsif ($key and not defined $value) {
		if (exists $self->{$index}{$key}) {
			# retrieve a value
			return $self->{$index}{$key};
		}
		else {
			# set a new empty key
			$self->{$index}{$key} = q();
			return 1;
		}
	}
	else {
		my %hash = %{ $self->{$index} };
		return wantarray ? %hash : \%hash;
	}
}

sub delete_metadata {
	my $self = shift;
	my ($index, $key) = @_;
	return unless defined $index;
	if (defined $key) {
		if (exists $self->{$index}{$key}) {
			return delete $self->{$index}{$key};
		}
	}
	else {
		# user wants to delete the metadata
		# but we need to keep the basics name and index
		foreach my $key (keys %{ $self->{$index} }) {
			next if $key eq 'name';
			next if $key eq 'index';
			delete $self->{$index}{$key};
		}
	}
}

sub copy_metadata {
	my ($self, $source, $target) = @_;
	return unless (exists $self->{$source}{name} and exists $self->{$target}{name});
	my $md = $self->metadata($source);
	delete $md->{name};
	delete $md->{'index'};
	delete $md->{'AUTO'} if exists $md->{'AUTO'}; # presume this is no longer auto index
	foreach (keys %$md) {
		$self->{$target}{$_} = $md->{$_};
	}
	return 1;
}



#### Column Indices ####

sub find_column {
	my ($self, $name) = @_;
	return unless $name;
	
	# the $name variable will be used as a regex in identifying the name
	# fix it so that it will possible accept a # character at the beginning
	# without a following space, in case the first column has a # prefix
	# also place the remainder of the text in a non-capturing parentheses for 
	# grouping purposes while maintaining the anchors
	$name =~ s/ \A (\^?) (.+) (\$?)\Z /$1#?(?:$2)$3/x;
	
	# walk through each column index
	my $index;
	for (my $i = 0; $i < $self->{'number_columns'}; $i++) {
		# check the names of each column
		if ($self->{$i}{'name'} =~ /$name/i) {
			$index = $i;
			last;
		}
	}
	return $index;
}

sub _find_column_indices {
	my $self = shift;
	# these are hard coded index name regex to accomodate different possibilities
	# these do not include parentheses for grouping
	# non-capturing parentheses will be added later in the sub for proper 
	# anchoring and grouping - long story why, don't ask
	my $name   = $self->find_column('^name|geneName|transcriptName|geneid|id|alias');
	my $type   = $self->find_column('^type|class|primary_tag');
	my $id     = $self->find_column('^primary_id');
	my $chromo = $self->find_column('^chr|seq|ref|ref.?seq');
	my $start  = $self->find_column('^start|position|pos|txStart$');
	my $stop   = $self->find_column('^stop|end|txEnd');
	my $strand = $self->find_column('^strand');
	$self->{column_indices} = {
		'name'      => $name,
		'type'      => $type,
		'id'        => $id,
		'seq_id'    => $chromo,
		'chromo'    => $chromo,
		'start'     => $start,
		'stop'      => $stop,
		'end'       => $stop,
		'strand'    => $strand,
	};
	return 1;
}

sub chromo_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{chromo};
}

sub start_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{start};
}

sub stop_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{stop};
}

sub end_column {
	return shift->stop_column;
}

sub strand_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{strand};
}

sub name_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{name};
}

sub type_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{type};
}

sub id_column {
	my $self = shift;
	$self->_find_column_indices unless exists $self->{column_indices};
	return $self->{column_indices}{id};
}



__END__

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  
