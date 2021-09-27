#!/usr/bin/perl -w
use strict;

my ($in,$out) = @ARGV;

open IN,"<",$in;
open OUT,">",$out;

while(<IN>){
	chomp;
	my @line = split/\t/,$_;
	if(($line[1] != 0) && ($line[1] != 16) && ($line[1] != 4)){
		print OUT"$_\n";
	}
}

close IN;
close OUT;