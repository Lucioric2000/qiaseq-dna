#!/bin/bash
## Perl
#cpan install CPAN
#cpan reload cpan
#conda_home=$1
#conda_preffix=$2
#conda_env=$3
#sudo ${conda_home} install -n ${conda_env} -c bioconda perl-archive-zip
#sudo ${conda_home} install -n ${conda_env} -c bioconda perl-datetime perl-dbi perl-dbd-sqlite perl-env-path perl-file-chdir perl-getopt-long-descriptive
#sudo ${conda_home} install -n ${conda_env} -c bioconda perl-sort-naturally perl-config-inifiles perl-data-dump-color perl-data-table-excel
#sudo ${conda_home} install -n ${conda_env} -c bioconda perl-hash-merge perl-file-slurp
#sudo ${conda_home} install -n ${conda_env} -c bioconda perl-data-dumper perl-scalar-util-numeric perl-data-dump
#sudo ${conda_home} install -n ${conda_env} -c dan_blanchard perl-config-inifiles
#sudo ${conda_preffix}/bin/cpanm install -n ${conda_env} Data::Table::Excel Data::Dump::Color Config::IniFiles

cpancmd='/usr/bin/cpanm install --sudo'
${cpancmd} HTTP::Tiny Archive::Zip Test::MockModule DBI Module::Runtime HTTP::Date Test::Pod Sub::Identify Data::Dump::Color Hash::Merge
${cpancmd} XML::Twig XML::XPath XML::Parser HTML::TreeBuilder YAML HTML::Entities HTML::HeadParser HTML::Element HTML::Tree Data::Table::Excel
${cpancmd} CPAN::Meta CPAN::Meta::YAML ExtUtils::CBuilder Module::Metadata Parse::CPAN::Meta Perl::OSType TAP::Harness JSON::PP File::Slurp
${cpancmd} IO::Socket::SSL DateTime DBI DBD::SQLite Env::Path File::chdir Getopt::Long::Descriptive Sort::Naturally Config::IniFiles
