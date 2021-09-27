#!/usr/bin/env perl
use warnings;
use strict;
use ExtUtils::Installed;

my $inst= ExtUtils::Installed->new();

my @modules = $inst->modules();

foreach(@modules)

{

        my $ver = $inst->version($_) || "???";

        printf("%-12s --  %s\n", $_, $ver);   

}
