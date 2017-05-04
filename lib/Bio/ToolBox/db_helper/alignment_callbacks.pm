package Bio::ToolBox::db_helper::alignment_callbacks;

# modules
require Exporter;
use strict;
use Carp;
use constant {
	CHR  => 0,  # chromosome
	STRT => 1,  # start
	STOP => 2,  # stop
	STR  => 3,  # strand
	STND => 4,  # strandedness
	METH => 5,  # method
	RETT => 6,  # return type
	DB   => 7,  # database object
	DATA => 8,  # first dataset, additional may be present
};
our $VERSION = '1.50';

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

sub _all_count_indexed {
	my ($a, $data) = @_;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++ if 
		( $pos >= $data->{start} and $pos <= $data->{stop} );
}

sub _all_precise_count_indexed {
	my ($a, $data) = @_;
	return unless ($a->pos >= $data->{start} and $a->calend <= $data->{stop} );
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++;
}

sub _all_name_indexed {
	my ($a, $data) = @_;
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos} ||= [];
	push @{ $data->{'index'}{$pos} }, $a->qname;
}

sub _all_count_array {
	my ($a, $data) = @_;
	push @{$data->{scores}}, 1;
}

sub _all_precise_count_array {
	my ($a, $data) = @_;
	return unless ($a->pos >= $data->{start} and $a->calend <= $data->{stop} );
	push @{$data->{scores}}, 1;
}

sub _all_name_array {
	my ($a, $data) = @_;
	push @{ $data->{scores} }, $a->qname; # query or read name
}

sub _forward_count_indexed {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and $a->reversed);
		return if (not $first and not $a->reversed);
	}
	else {
		return if $a->reversed;
	}
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++ if 
		( $pos >= $data->{start} and $pos <= $data->{stop} );
}

sub _forward_precise_count_indexed {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and $a->reversed);
		return if (not $first and not $a->reversed);
	}
	else {
		return if $a->reversed;
	}
	return unless ($a->pos >= $data->{start} and $a->calend <= $data->{stop} );
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++;
}

sub _forward_name_indexed {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and $a->reversed);
		return if (not $first and not $a->reversed);
	}
	else {
		return if $a->reversed;
	}
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos} ||= [];
	push @{ $data->{'index'}{$pos} }, $a->qname;
}

sub _forward_count_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and $a->reversed);
		return if (not $first and not $a->reversed);
	}
	else {
		return if $a->reversed;
	}
	push @{$data->{scores}}, 1;
}

sub _forward_precise_count_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and $a->reversed);
		return if (not $first and not $a->reversed);
	}
	else {
		return if $a->reversed;
	}
	return unless ($a->pos >= $data->{start} and $a->calend <= $data->{stop} );
	push @{$data->{scores}}, 1;
}

sub _forward_name_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and $a->reversed);
		return if (not $first and not $a->reversed);
	}
	else {
		return if $a->reversed;
	}
	push @{ $data->{scores} }, $a->qname;
}

sub _reverse_count_indexed {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and not $a->reversed);
		return if (not $first and $a->reversed);
	}
	else {
		return unless $a->reversed;
	}
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++ if 
		( $pos >= $data->{start} and $pos <= $data->{stop} );
}

sub _reverse_precise_count_indexed {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and not $a->reversed);
		return if (not $first and $a->reversed);
	}
	else {
		return unless $a->reversed;
	}
	return unless ($a->pos >= $data->{start} and $a->calend <= $data->{stop} );
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos}++;
}

sub _reverse_name_indexed {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and not $a->reversed);
		return if (not $first and $a->reversed);
	}
	else {
		return unless $a->reversed;
	}
	my $pos = int( ($a->pos + 1 + $a->calend) / 2);
	$data->{'index'}{$pos} ||= [];
	push @{ $data->{'index'}{$pos} }, $a->qname;
}

sub _reverse_count_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and not $a->reversed);
		return if (not $first and $a->reversed);
	}
	else {
		return unless $a->reversed;
	}
	push @{$data->{scores}}, 1;
}

sub _reverse_precise_count_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and not $a->reversed);
		return if (not $first and $a->reversed);
	}
	else {
		return unless $a->reversed;
	}
	return unless ($a->pos >= $data->{start} and $a->calend <= $data->{stop} );
	push @{$data->{scores}}, 1;
}

sub _reverse_name_array {
	my ($a, $data) = @_;
	if ($a->paired) {
		my $first = $a->get_tag_values('FIRST_MATE');
		return if ($first and not $a->reversed);
		return if (not $first and $a->reversed);
	}
	else {
		return unless $a->reversed;
	}
	push @{ $data->{scores} }, $a->qname;
}


