#!/bin/bash
## Perl
#cpan install CPAN
#cpan reload cpan
sudo /srv/conda/bin/cpanm install HTTP::Tiny Archive::Zip Test::MockModule
sudo /srv/conda/bin/cpanm install XML XML::Twig XML::XPath XML::Parser HTML::TreeBuilder YAML HTML::Entities HTML::HeadParser HTML::Element HTML::Tree
sudo /srv/conda/bin/cpanm install CPAN::Meta CPAN::Meta::YAML ExtUtils::CBuilder Module::Metadata Parse::CPAN::Meta Perl::OSType TAP::Harness JSON::PP

sudo /srv/conda/bin/cpanm install Module::Runtime HTTP::Date Test::Pod Sub::Identify
sudo /srv/conda/bin/cpanm install IO::Socket::SSL DateTime DBI DBD::SQLite Env::Path File::chdir Getopt::Long::Descriptive Sort:Naturally Config::IniFiles Data::Dump::Color Data::Table::Excel Hash::Merge File::Slurp
