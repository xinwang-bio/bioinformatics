#!/usr/bin/perl -w
use strict;

my ($in,$tpm,$out) = @ARGV;

my %ha;
open IN,"<",$in;

while(<IN>){
	chomp;
	if(/>/){
		
	}
	else{
		my @line = split;
		$ha{$line[0]} = 0;
	}
}
close IN;

open OUT,">",$out;
open IN,"<",$tpm;

while(<IN>){
	chomp;
	if(/\#/){
		print OUT"$_\n";
	}
	else{
		my @line = split/\t/,$_;
		if(exists $ha{$line[1]}){
			print OUT"$_\n";
		}
	}
}

close IN;
close OUT;
