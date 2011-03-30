use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Statistics::Autocorrelation',
    license             => 'perl',
    dist_author         => q{Roderick Garton <rgarton@cpan.org>},
    dist_version_from   => 'lib/Statistics/Autocorrelation.pm',
    build_requires => {
        'Test::More' => 0,
        'Statistics::Lite' => '3.0',
    },
    add_to_cleanup      => [ 'Statistics-Autocorrelation-0.02-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
