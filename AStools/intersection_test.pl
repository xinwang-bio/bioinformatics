#!/usr/bin/perl -w
use strict;
use Set::IntSpan;
my $set = new Set::IntSpan('1-10, 30-40, 100-200');

#$set->add( '1-10, 30-40, 100-200' );

my $a = new Set::IntSpan('5-8, 41-51, 92-150');
print"$set\n$a\n";
#$a-> add ('5-8, 41-51, 92-150' );

 my $in = $set->I($a);
 print "$in\n";

 $in = $set->X($a);
 print "$in\n";
