package Bio::ToolBox::db_helper::alignment_callbacks;

# modules
require Exporter;
use strict;
use Carp;
use Bio::ToolBox::db_helper::constants;
our $VERSION = '1.51';

# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	assign_callback
	_all_count_indexed
	_all_precise_count_indexed
	_all_name_indexed
	_all_count_array
	_all_precise_count_array
	_all_name_array
	_forward_count_indexed
	_forward_precise_count_indexed
	_forward_name_indexed
	_forward_count_array
	_forward_precise_count_array
	_forward_name_array
	_reverse_count_indexed
	_reverse_precise_count_indexed
	_reverse_name_indexed
	_reverse_count_array
	_reverse_precise_count_array
	_reverse_name_array
);

# Lookup hash for caching callback methods
our %CALLBACKS;

# The true statement
1; 



### Generate callback subroutine for walking through Bam alignments
sub assign_callback {
	# generate the callback code depending on whether we want to look at 
	# stranded data, collecting counts or length, or whether indexed data
	# is wanted.
	
	# we perform a check of whether the alignment midpoint is within the 
	# search region for indexed data only
	
	# these subroutines are designed to work with the low level fetch API
	
	# there are so many different subroutines because I want to increase 
	# efficiency by limiting the number of conditional tests in one generic subroutine
	
	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, value, db, dataset
	my $param = shift;
	
	
	# check the current list of calculated callbacks
		# cache the calculated callback method to speed up subsequent data 
		# collections it's likely only one method is ever employed in an 
		# execution, but just in case we will cache all that we calculate
	my $string = sprintf "%s_%s_%s_%d", $param->[STND], $param->[STR], 
		$param->[METH], $param->[RETT];
	return $CALLBACKS{$string} if exists $CALLBACKS{$string};
	
	# determine the callback method based on requested criteria
	my $callback;
	
	# all reads, either strand
	if (
		$param->[STND] eq 'all' and 
		$param->[METH] eq 'count' and 
		$param->[RETT] == 2
	) {
		$callback = \&_all_count_indexed;
	}
	elsif (
		$param->[STND] eq 'all' and 
		$param->[METH] eq 'pcount' and 
		$param->[RETT] == 2
	) {
		$callback = \&_all_precise_count_indexed;
	}
	elsif (
		$param->[STND] eq 'all' and 
		$param->[METH] eq 'count' and 
		$param->[RETT] != 2
	) {
		$callback = \&_all_count_array;
	}
	elsif (
		$param->[STND] eq 'all' and 
		$param->[METH] eq 'pcount' and 
		$param->[RETT] != 2
	) {
		$callback = \&_all_precise_count_array;
	}
	elsif (
		$param->[STND] eq 'all' and 
		$param->[METH] eq 'ncount' and 
		$param->[RETT] != 2
	) {
		$callback = \&_all_name_array;
	}
	elsif (
		$param->[STND] eq 'all' and 
		$param->[METH] eq 'ncount' and 
		$param->[RETT] == 2
	) {
		$callback = \&_all_name_indexed;
	}
	
	
	# sense, forward strand 
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'count' and 
		$param->[RETT] == 2
	) {
		$callback = \&_forward_count_indexed;
	}
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'pcount' and 
		$param->[RETT] == 2
	) {
		$callback = \&_forward_precise_count_indexed;
	}
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'count' and 
		$param->[RETT] != 2
	) {
		$callback = \&_forward_count_array;
	}
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'pcount' and 
		$param->[RETT] != 2
	) {
		$callback = \&_forward_precise_count_array;
	}
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'ncount' and 
		$param->[RETT] != 2
	) {
		$callback = \&_forward_name_array;
	}
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'ncount' and 
		$param->[RETT] == 2
	) {
		$callback = \&_forward_name_indexed;
	}
	
	
	# sense, reverse strand
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'count' and 
		$param->[RETT] == 2
	) {
		$callback = \&_reverse_count_indexed;
	}
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'pcount' and 
		$param->[RETT] == 2
	) {
		$callback = \&_reverse_precise_count_indexed;
	}
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'ncount' and 
		$param->[RETT] == 2
	) {
		$callback = \&_reverse_name_indexed;
	}
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'count' and 
		$param->[RETT] != 2
	) {
		$callback = \&_reverse_count_array;
	}
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'pcount' and 
		$param->[RETT] != 2
	) {
		$callback = \&_reverse_precise_count_array;
	}
	elsif (
		$param->[STND] eq 'sense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'ncount' and 
		$param->[RETT] != 2
	) {
		$callback = \&_reverse_name_array;
	}
	
	
	# anti-sense, forward strand 
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'count' and 
		$param->[RETT] == 2
	) {
		$callback = \&_reverse_count_indexed;
	}
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'pcount' and 
		$param->[RETT] == 2
	) {
		$callback = \&_reverse_precise_count_indexed;
	}
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'ncount' and 
		$param->[RETT] == 2
	) {
		$callback = \&_reverse_name_indexed;
	}
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'count' and 
		$param->[RETT] != 2
	) {
		$callback = \&_reverse_count_array;
	}
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'pcount' and 
		$param->[RETT] != 2
	) {
		$callback = \&_reverse_precise_count_array;
	}
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] >= 0 and 
		$param->[METH] eq 'ncount' and 
		$param->[RETT] != 2
	) {
		$callback = \&_reverse_name_array;
	}
	
	
	# anti-sense, reverse strand
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'count' and 
		$param->[RETT] == 2
	) {
		$callback = \&_forward_count_indexed;
	}
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'pcount' and 
		$param->[RETT] == 2
	) {
		$callback = \&_forward_precise_count_indexed;
	}
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'ncount' and 
		$param->[RETT] == 2
	) {
		$callback = \&_forward_name_indexed;
	}
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'count' and 
		$param->[RETT] != 2
	) {
		$callback = \&_forward_count_array;
	}
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'pcount' and 
		$param->[RETT] != 2
	) {
		$callback = \&_forward_precise_count_array;
	}
	elsif (
		$param->[STND] eq 'antisense' and 
		$param->[STR] == -1 and 
		$param->[METH] eq 'ncount' and 
		$param->[RETT] != 2
	) {
		$callback = \&_forward_name_array;
	}
	
	# I goofed
	else {
		confess sprintf "Programmer error: stranded %s, strand %s, method %s, do_index %d", 
			$param->[STND], $param->[STR], $param->[METH], $param->[RETT];
	}
	
	# remember next time 
	$CALLBACKS{$string} = $callback;
	
	return $callback;
}


#### Callback subroutines 
# the following are all of the callback subroutines 

# we explicitly check that alignment start and/or stop overlap the search coordinates
# to avoid e.g. weird alignments with massive splice junctions that would nevertheless 
# be pulled out by the bam index search

sub _all_count_indexed {
	my ($a, $data) = @_;
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	if ($a->reversed) {
		$data->{'index'}{$e}++;
	}
	else {
		$data->{'index'}{$s}++;
	}
}

sub _all_precise_count_indexed {
	my ($a, $data) = @_;
	my $s = $a->pos + 1;
	return unless ($s >= $data->{start} and $a->calend <= $data->{stop} );
	if ($a->reversed) {
		$data->{'index'}{$s}++;
	}
	else {
		$data->{'index'}{$a->calend}++;
	}
}

sub _all_name_indexed {
	my ($a, $data) = @_;
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	if ($a->reversed) {
		# since we're working with names, only record the 5' end of fragment
		return if ($a->paired and $a->proper_pair); # not the end we're looking for
		$data->{'index'}{$e} ||= [];
		push @{ $data->{'index'}{$e} }, $a->qname;
	}
	else {
		$data->{'index'}{$s} ||= [];
		push @{ $data->{'index'}{$s} }, $a->qname;
# 		$data->{'index'}{$s} .= $a->qname . "\t";
# 		$data->{'index'}{$a->qname} = $s;
	}
}

sub _all_count_array {
	my ($a, $data) = @_;
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	push @{$data->{scores}}, 1;
}

sub _all_precise_count_array {
	my ($a, $data) = @_;
	return unless ($a->pos + 1 >= $data->{start} and $a->calend <= $data->{stop} );
	push @{$data->{scores}}, 1;
}

sub _all_name_array {
	my ($a, $data) = @_;
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	push @{ $data->{scores} }, $a->qname;
}

sub _forward_count_indexed {
	my ($a, $data) = @_;
	my $r = $a->reversed;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and $r);
		return if (not $first and not $r);
	}
	else {
		return if $r;
	}
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	if ($r) {
		$data->{'index'}{$e}++;
	}
	else {
		$data->{'index'}{$s}++;
	}
}

sub _forward_precise_count_indexed {
	my ($a, $data) = @_;
	my $r = $a->reversed;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and $r);
		return if (not $first and not $r);
	}
	else {
		return if $r;
	}
	my $s = $a->pos + 1;
	return unless ($s >= $data->{start} and $a->calend <= $data->{stop} );
	if ($r) {
		$data->{'index'}{$s}++;
	}
	else {
		$data->{'index'}{$a->calend}++;
	}
}

sub _forward_name_indexed {
	my ($a, $data) = @_;
	my $r = $a->reversed;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and $r);
		return if (not $first and not $r);
	}
	else {
		return if $r;
	}
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	if ($r) {
		$data->{'index'}{$e} ||= [];
		push @{ $data->{'index'}{$e} }, $a->qname;
	}
	else {
		$data->{'index'}{$s} ||= [];
		push @{ $data->{'index'}{$s} }, $a->qname;
	}
}

sub _forward_count_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and $a->reversed);
		return if (not $first and not $a->reversed);
	}
	else {
		return if $a->reversed;
	}
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	push @{$data->{scores}}, 1;
}

sub _forward_precise_count_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and $a->reversed);
		return if (not $first and not $a->reversed);
	}
	else {
		return if $a->reversed;
	}
	return unless ($a->pos + 1 >= $data->{start} and $a->calend <= $data->{stop} );
	push @{$data->{scores}}, 1;
}

sub _forward_name_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and $a->reversed);
		return if (not $first and not $a->reversed);
	}
	else {
		return if $a->reversed;
	}
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	push @{ $data->{scores} }, $a->qname;
}

sub _reverse_count_indexed {
	my ($a, $data) = @_;
	my $r = $a->reversed;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and not $r);
		return if (not $first and $r);
	}
	else {
		return unless $r;
	}
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	if ($r) {
		$data->{'index'}{$e}++;
	}
	else {
		$data->{'index'}{$s}++;
	}
}

sub _reverse_precise_count_indexed {
	my ($a, $data) = @_;
	my $r = $a->reversed;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and not $r);
		return if (not $first and $r);
	}
	else {
		return unless $r;
	}
	my $s = $a->pos + 1;
	return unless ($s >= $data->{start} and $a->calend <= $data->{stop} );
	if ($r) {
		$data->{'index'}{$s}++;
	}
	else {
		$data->{'index'}{$a->calend}++;
	}
}

sub _reverse_name_indexed {
	my ($a, $data) = @_;
	my $r = $a->reversed;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and not $r);
		return if (not $first and $r);
	}
	else {
		return unless $r;
	}
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	if ($r) {
		$data->{'index'}{$e} ||= [];
		push @{ $data->{'index'}{$e} }, $a->qname;
	}
	else {
		$data->{'index'}{$s} ||= [];
		push @{ $data->{'index'}{$s} }, $a->qname;
	}
}

sub _reverse_count_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and not $a->reversed);
		return if (not $first and $a->reversed);
	}
	else {
		return unless $a->reversed;
	}
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	push @{$data->{scores}}, 1;
}

sub _reverse_precise_count_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and not $a->reversed);
		return if (not $first and $a->reversed);
	}
	else {
		return unless $a->reversed;
	}
	return unless ($a->pos + 1 >= $data->{start} and $a->calend <= $data->{stop} );
	push @{$data->{scores}}, 1;
}

sub _reverse_name_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->flag & 0x40; # true if FIRST_MATE
		return if ($first and not $a->reversed);
		return if (not $first and $a->reversed);
	}
	else {
		return unless $a->reversed;
	}
	my $s = $a->pos + 1;
	my $e = $a->calend;
	return unless ( ($s >= $data->{start} and $s <= $data->{stop}) or 
					($e >= $data->{start} and $e <= $data->{stop}) );
	push @{ $data->{scores} }, $a->qname;
}

__END__

=head1 NAME

Bio::ToolBox::db_helper::alignment_callbacks

=head1 DESCRIPTION

This module provides common callback subroutines for working with bam alignments. 
It is generalized and may be used with either L<Bio::DB::Sam> or L<Bio::DB::HTS> 
alignments. 

It's not meant to be used by individuals. 
Please see L<Bio::ToolBox::db_helper::bam> or L<Bio::ToolBox::db_helper::hts>.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  



