use strict;
use warnings;
use Test::More tests => 4;
use constant EPS => 1e-9;

BEGIN { use_ok('Statistics::Autocorrelation') };

my $acorr = Statistics::Autocorrelation->new();
isa_ok($acorr, 'Statistics::Autocorrelation');

my %refvals = (
    coeff_1 => -.5,
);

my @data = (1, 2);

my $coeff;

$coeff = $acorr->coefficient(data => \@data);
ok( about_equal($coeff, $refvals{'coeff_1'}), "Coefficient lag 1: $coeff = $refvals{'coeff_1'}");

# Check alias:
$coeff = $acorr->coeff(data => \@data);
ok( about_equal($coeff, $refvals{'coeff_1'}), "Coefficient lag 1 (by alias): $coeff = $refvals{'coeff_1'}");

sub about_equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}



