#!/usr/bin/perl -w
use strict;
use File::Basename;

my $dirname = dirname(__FILE__);

print "$dirname\n";

my ($name, $path, $suffix)=fileparse($dirname);

print "$name\n";
print "$path\n";
print "$suffix\n";
