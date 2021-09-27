#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 CheesRIL_p204_LA483_diffHomo_wi35offs.tab.indv_cnt > CheesRIL_p204_LA483_diffHomo_wi35offs.tab.indv_cntR\n"; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1) {
		$ta[0] eq 'chr' or die "$_\n"; 
		for my $tb (@ta[4..($#ta-1)]) {
			$tb =~ s!_N$!_R!i; 
		}
		print STDOUT join("\t", @ta)."\n"; 
		next; 
	}
	$ta[8] > 0 or next; 
	for my $tb (@ta[4 .. ($#ta-1)]) {
		$tb = sprintf("%02.2f", $tb/$ta[8]*100); 
	}
	print STDOUT join("\t", @ta)."\n"; 
}
