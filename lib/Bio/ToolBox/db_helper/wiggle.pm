package Bio::ToolBox::db_helper::wiggle;

use warnings;
use strict;
use Carp;
use Bio::ToolBox::db_helper::constants;
use Bio::Graphics::Wiggle;
require Exporter;

our $VERSION = '1.51';

# Exported names
our @ISA    = qw(Exporter);

## no critic
## this is never intended to be used directly by end users
## and exporting everything is required
our @EXPORT = qw(
	collect_wig_scores
	collect_wig_position_scores
);
## use critic

# Hashes of opened file objects
my %OPENED_WIGFILES;  # opened wigfile objects
                       # in empirical testing, this doesn't really seem to speed things up
                       # like I thought it would
                       # oh well, keep it anyway????
    # I think this is safe to keep opened wigfiles cached, even across forks,
    # since they are being opened only during data collection, which should
    # only occur within child processes, and there is no explicit db open

sub collect_wig_scores {

	# we will actually call collect_wig_position_scores()
	# but only return the values
	my $wig_data = collect_wig_position_scores(shift);
	return unless $wig_data;

	# return the values
	my @values = values %{ $wig_data };
	return wantarray ? @values : \@values;
}

sub collect_wig_position_scores {

	# passed parameters as array ref
	# chromosome, start, stop, strand, strandedness, method, db, dataset
	my $param = shift;

	# look at each wigfile
	# usually there is only one, but for stranded data there may be
	# two wigfiles (+ and -), so we'll check each wig file for strand info
	my %pos2score;    # position => score
	for ( my $d = DATA; $d < scalar @{ $param }; $d++ ) {

		my $feature = $param->[$d];
		confess 'dataset is not a seqfeature object!'
			unless ref($feature) =~ /seqfeature/i;

		# Check which data to take based on strand
		if (
			$param->[STND] eq 'all'     # all data is requested
			or $feature->strand == 0    # unstranded data
			or (
				# sense data
				$param->[STR] == $feature->strand and $param->[STND] eq 'sense'
			)
			or (
				# antisense data
				$param->[STR] != $feature->strand and $param->[STND] eq 'antisense'
			)
			)
		{
			# we have acceptable data to collect

			# collect from wigfile if present
			if ( $feature->has_tag('wigfile') ) {

				# get wigfile name
				my @wigfiles = $feature->get_tag_values('wigfile');
				my $wigfile  = shift @wigfiles;    # there should only be one wigfile
				confess ' no wigfile passed!' unless $wigfile;

				# check for opened wigfile
				my $wig;
				if ( exists $OPENED_WIGFILES{$wigfile} ) {

					# this file is already opened, use it
					$wig = $OPENED_WIGFILES{$wigfile};
				}
				else {
					# this file has not been opened yet, open it
					unless ( -e $wigfile ) {
						confess " Binary wiggle file '$wigfile' does not exist!\n";
					}
					$wig = Bio::Graphics::Wiggle->new( $wigfile, 0 );
					unless ($wig) {
						confess " unable to open data wigfile '$wigfile'";
					}

					# store the opened object for later use
					$OPENED_WIGFILES{$wigfile} = $wig;
				}

				# adjust as necessary to avoid wig errors
				if ( $param->[STRT] < $wig->start ) {

					# adjust the start position
					$param->[STRT] = $wig->start;
				}
				elsif ( $param->[STRT] > $wig->end ) {

					# nothing we can do here, no values
					return;
				}
				if ( $param->[STOP] > $wig->end ) {

					# adjust the end position
					$param->[STOP] = $wig->end;
				}
				elsif ( $param->[STOP] < $wig->start ) {

					# nothing we can do here, no values
					return;
				}

				# collect the wig values
				my $scores_ref = $wig->values( $param->[STRT] => $param->[STOP] );

				# re-associate position with the scores
				my $step =
					$wig->step || 1;    # step should always be defined but just in case
				my $pos = $param->[STRT];
				foreach my $s ( @{$scores_ref} ) {

					#print Dumper($s);
					if ( defined $s ) {

						# the binary wig file (.wib) is usually set up with
						# a step of 1 bp, even if the original wig file was not
						# this can result in lots of undefined values at the
						# positions where there was no original data
						# hence the defined check here
						# store a real value in the hash keyed under the position
						if ( $param->[METH] eq 'count' ) {
							$pos2score{$pos} = 1;
						}
						else {
							$pos2score{$pos} = $s;
						}
					}

					# adjust position by the step size,
					$pos += $step;
				}
			}
		}
	}

	# return the wig data hash
	return wantarray ? %pos2score : \%pos2score;
}

1;

__END__


=head1 NAME

Bio::ToolBox::db_helper::wiggle

=head1 DESCRIPTION

This module provides support for legacy binary wiggle (F<.wib>) files. These are 
not BigWig files, which are supported through the L<Bio::ToolBox::Data::db_helper::bigwig> 
module. Rather, these are condensed binary representations of text wiggle files 
developed for use with the L<GBrowse browser|http://gmod.org/gbrowse>. Use the 
modern BigWig adapter for improved and more efficient access.

=head1 USAGE

The module requires L<Bio::Perl> and L<Bio::Graphics> to be installed.

In general, this module should not be used directly. Use the methods 
available in L<Bio::ToolBox::db_helper> or <Bio::ToolBox::Data>.  

All subroutines are exported by default.

=over

=item collect_wig_scores

This subroutine will collect dataset scores from a binary wig file (.wib).

The subroutine is passed a parameter array reference. See below for details.

The subroutine returns an array or array reference of the requested dataset 
values found within the region of interest. 

=item collect_wig_position_scores

This subroutine will collect the score values from a binary wig file 
for the specified database region keyed by position. 

The subroutine is passed a parameter array reference. See 
L</"Data Collection Parameters Reference"> below for details.

The subroutine returns a hash or hash reference of the requested dataset values 
found within the region of interest keyed by position. Note that only one 
value is returned per position, regardless of the number of dataset features 
passed.

=back

=head2 Binary wiggle files

Binary wiggle files (F<.wib>) files are referenced via a SeqFeature object. 
These features are typically stored from a L<Bio::DB::SeqFeature::Store> 
database. A single feature representing the dataset is present across 
each chromosome. The feature should contain an attribute C<wigfile> that 
references the location of the binary file representing the dataset scores. 
The file is opened and the values extracted from the region of interest. 

For loading wig files into a L<Bio::DB::SeqFeature::Store> database, see 
the perl script 'wiggle2gff3.pl' included with the L<Bio::Graphics> 
distribution, as well as L<Bio::Graphics::Wiggle::Loader>.

To speed up the program and avoid repetitive opening and 
closing of the files, the opened wig file object is stored in a global 
hash in case it is needed again.

=head2 Data Collection Parameters Reference

The data collection subroutines are passed an array reference of parameters. 
The recommended  method for data collection is to use the 
L<Bio::ToolBox::db_helper/get_segment_score> method. 

The parameters array reference includes these items:

=over 4

=item 1. chromosome

=item 2. start coordinate

=item 3. stop coordinate 

Coordinates are in BioPerl-style 1-base system.

=item 4. strand

Should be standard BioPerl representation: -1, 0, or 1.

=item 5. strandedness

A scalar value representing the desired strandedness of the data 
to be collected. Acceptable values include "sense", "antisense", 
or "all". Only those scores which match the indicated 
strandedness are collected.

=item 6. score method

Acceptable values include score and count.

   * score returns the basepair coverage of alignments over the 
   region of interest
   
   * count returns the number of alignments that overlap the 
   search region. 

=item 7. A database object.

Not used here.

=item 8. Database features.

These are the SeqFeature objects that contain the file path to the 
binary wig files. Provide the C<type> or C<primary_tag>.

=back

=head1 SEE ALSO

L<Bio::ToolBox::Data::Feature>, L<Bio::ToolBox::db_helper>, 
L<Bio::Graphics>, L<Bio::DB::SeqFeature::Store>

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  




