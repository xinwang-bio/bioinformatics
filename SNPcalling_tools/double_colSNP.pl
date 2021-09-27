#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 LA483.10_3_0p3.comm.snp\n"; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1) {
		$ta[0] eq 'chr' or die "$_\n"; 
		print STDOUT "$_\n"; 
		next; 
	}
	for (my $i=3; $i<@ta; $i++) {
		if ($ta[$i] eq 'N') {
			$ta[$i] = './.'; 
		} elsif ($ta[$i] =~ m!^[ATGC*]$!) {
			$ta[$i] = "$ta[$i]/$ta[$i]"; 
		} elsif ($ta[$i] =~ s!^([ATGC*])([ATGC*])$!$1/$2!) {
			; 
		} elsif ($ta[$i] =~ m!^([ATGC])\+([ATGC]+)$!) {
			$ta[$i] = "$1$2"; 
		} elsif ($ta[$i] =~ m!^[ATGC]{3,}$!) {
			$ta[$i] = './.'; 
		} else {
			die "Bad [$ta[$i]] : $_\n"; 
		}
	}
	print join("\t", @ta)."\n"; 
}

