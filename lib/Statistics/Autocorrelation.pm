package Statistics::Autocorrelation;

use warnings;
use strict;
use Carp 'croak';
use Statistics::Lite qw(mean variance);

=head1 NAME

Statistics-Autocorrelation - Coefficients for any lag

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Statistics-Autocorrelation - Coefficients for any lag

    use Statistics::Autocorrelation;
    $acorr = Statistics::Autocorrelation->new();
    $coeff = $acorr->coefficient(data => \@data, lag => integer (from 1 to N-1), exact => 0, simple => 1);

=head1 DESCRIPTION

Calculates autocorrelation coefficients for a single series of numerical data, for any valid length of I<lag>.

=head1 SUBROUTINES/METHODS

=head2 new

 $acorr = Statistics::Autocorrelation->new();

Return a new module object for accessing its methods.

=cut

sub new {
	my $class = shift;
	my $self = {};
	bless $self, $class;
    #$self->init(@_);
    return $self;
}

=head2 coefficient

 $autocorr->coefficient(data => \@data, lag => integer (from 1 to N-1), exact => 0|1, simple => 1|0);

Corporate stats programs (e.g., SPSS/PASW), and examples of autocorrelation on the web (e.g., L<http://www.itl.nist.gov/div898/handbook/eda/section3/eda35c.htm|itl.nist.gov>), often yield/demonstrate a value of the autocorrelation coefficient that assumes that the underlying series is stationary (has no linear or curvilinear trend, no periodicity), and so can be calculated as the population autocorrelation coefficient, the ratio of the autocovariance to the overall variance. To be a valid estimate of the sample autocorrelation coefficient, it is assumed that the number of observations, I<N>, in the sample is "reasonably large" - so that all the observations in each summation in the numerator are taken relative to the mean of the whole series, rather than the exact value at lag I<k>, and also that the variance used in the denominator can be that of the whole series - instead of using completely pairwise products, i.e., making the coefficient dependent on the values of I<u1>, ..., I<un-k>, and I<ui+k>, ..., I<un> at each summation. Additionally, these sources also commonly drop the factor I<I>/(I<N> - 1), assuming that the value is close enough to 1 for large I<N>.

By default, then, this method returns an estimate of the population autocorrelation coefficient, and the divisor factors are dropped; i.e., the default values of the options C<exact> = 0, and C<simple> = 1. The default value returned from this method, then, is equivalent to those returned from such corporate software, and demonstrated via such URLs, as cited; and this is also the default form of calculating the coefficient as used in texts such as Chatfield (1975). If you want, however, to keep these divisors, you need to specify C<simple> = 0; and if you want, furthermore, the exact sample autocorrelation coefficient, then specify C<exact> = 1; then you get the coefficient as calculated by Kendall's (1973) Eq. 3.35.

A croak will be heard if no value is given for C<data> (expects an array reference), or if the array is empty. A value for C<data> is the only one for which there is no default value and must be given.

If a value is not given for C<lag>, or it equals zero, it is set to the value of 1 by default. 

=cut

sub coefficient {
	my $self = shift;
    my (%args) = @_;
	my $dat = ref $args{'data'} ? $args{'data'} : croak 'No value for data for calculating coefficient';
	my $n = scalar(@{$dat}) || croak 'No data are available for calculating coefficient';
	croak 'Can\'t autocorrelate a dataset of only one element' if $n < 2;
	my $lag = $args{'lag'} || 1;
	#croak "Can\'t autocorrelate with a lag of < $lag > for only < $n > data" if $n - $lag == 1 && !$args{'circular'};
    $args{'exact'} ||= 0;

	my $rk;
	if ($args{'exact'} == 1) {
		if ($args{'circular'}) {
			$rk = _calc_exact_circ($dat, $n, $lag);
		}
		else {
			$rk = _calc_exact($dat, $n, $lag);
		}
	}
	else {
		if ($args{'circular'}) {
			$rk = _calc_approx_circ($dat, $n, $lag, $args{'simple'});
		}
		else {
			$rk = _calc_approx($dat, $n, $lag, $args{'simple'});
		}
	}
	return $rk;
}

# Kendall's (1973) Eq. 3.35., p. 40

sub _calc_exact {
	my ($dat, $n, $k) = @_;

	my $c0 = 1 / ($n - $k);
	my ($sum_ui, $sum_uik, $sum_prodsum, $sum_sq_ui, $sum_sq_uik, $i) = (0, 0, 0, 0);

	for ($i = 0; $i < $n - $k; $i++) {
		$sum_ui += $dat->[$i];
	}
	for ($i = 0; $i < $n - $k; $i++) {
		$sum_uik += $dat->[$i + $k];
	}

	my $c1 = $c0 * $sum_ui;
	my $c2 = $c0 * $sum_uik;

	# Numerator:
	for ($i = 0; $i < $n - $k; $i++) {
		$sum_prodsum += ($dat->[$i] - $c1) * ($dat->[$i + $k] - $c2);
	}
	my $numerator = $c0 * $sum_prodsum;

	# Denominator:
	for ($i = 0; $i < $n - $k; $i++) {
		$sum_sq_ui += ($dat->[$i] - $c1)**2;
	}
	for ($i = 0; $i < $n - $k; $i++) {
		$sum_sq_uik += ($dat->[$i + $k] - $c2)**2;
	}

	my $denominator = sqrt($c0 * $sum_sq_ui * $c0 * $sum_sq_uik);

	my $rk = $numerator / $denominator;

	return $rk;

}

# Kendall Eq. 3.36

sub _calc_approx {
	my ($dat, $n, $k) = @_;
	my $simple = shift;
	my $mean = mean(@{$dat});
	my ($rk, $sum_resid, $sum_sq, $i) = ();

	for ($i = 0; $i < $n - $k; $i++) {
		$sum_resid += ($dat->[$i] - $mean)*($dat->[$i + $k] - $mean);
	}

	$sum_sq = _sum_sq($dat, $n, $mean);
	
	if (defined $simple && !$simple) {
		$rk = ($sum_resid / ($n - $k)) / ($sum_sq/$n);
	}
	else {
		$rk = $sum_resid / $sum_sq;
	}
	
	return $rk;
}

sub _calc_approx_circ {
	my ($dat, $n, $k) = @_;
	my $simple = shift || 0;
	my $mean = mean(@{$dat});
	my ($rk, $sum_resid, $sum_sq, $i) = ();

	for ($i = 0; $i < $n; $i++) { # go right up to the last value in the data
		my $i2 = $i + $k >= $n ? $k : $i + $k;
		$sum_resid += ($dat->[$i] - $mean)*($dat->[$i2] - $mean);
	}
	
	$sum_sq = _sum_sq($dat, $n, $mean);
	
	if (!$simple) {
		$rk = ($sum_resid / ($n - $k)) / ($sum_sq / $n);
	}
	else {
		$rk = $sum_resid / $sum_sq;
	}
	
	return $rk;
}

sub _sum_sq {
	my ($dat, $n, $mean) = @_;
	my ($i, $sum) = ();
	for ($i = 0; $i < $n; $i++) {  
		$sum += ($dat->[$i] - $mean)**2;
	}
	return $sum;
}

=head1 REFERENCES

Chatfield, C. (1975). I<The analysis of time series: Theory and practice>. London, UK: Chapman and Hall.

Kendall, M. G. (1973). I<Time-series>. London, UK: Griffin.

=head1 AUTHOR

Roderick Garton, C<< <rgarton at cpan.org> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-statistics-autocorrelation-0.01 at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Statistics-Autocorrelation-0.01>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Statistics::Autocorrelation


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Statistics-Autocorrelation-0.01>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Statistics-Autocorrelation-0.01>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Statistics-Autocorrelation-0.01>

=item * Search CPAN

L<http://search.cpan.org/dist/Statistics-Autocorrelation-0.01/>

=back

=head1 LICENSE AND COPYRIGHT

Copyright 2011 Roderick Garton.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.

=cut

1; # End of Statistics::Autocorrelation
