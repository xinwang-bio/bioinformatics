#!/usr/bin/env perl
use strict;
use warnings;

if(@ARGV < 1){
	&usage;
};

for my $i(@ARGV){
	my $cmd = "ps aux \| grep $i \| awk \'{print \$2}\' \| xargs kill \-9";
	print "$cmd\n";
	system $cmd;
}



sub usage {
	print "Usage: $0 -OPTIONS VALUES\n\n";
    print "batch_kill.pl keyword1 keyword2 keyword3 ......\n";
    exit;
}