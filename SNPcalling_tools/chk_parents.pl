#!/usr/bin/perl
use strict;
use warnings; 

-t and !@ARGV and die "perl $0 CheesRIL_p204_LA483_35offs_rmfilt00V_use_snp.tab\n"; 

my $p1_c = 3; 
my $p2_c = 4; 
my @hh; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1) {
		$ta[0] eq 'chr' or die "$_\n"; 
		print STDOUT join("\t", qw/chr pos/, $ta[$p1_c], $ta[$p2_c], qw/GoodCnt BadCnt MissCnt Hete% BadIDs/)."\n"; 
		@hh = @ta; 
		next; 
	}
	my %p_al; 
	my %off_cnt; 
	my @out_acc; 
	for my $tb (@ta[$p1_c,$p2_c]) {
		for my $tc (split(/\//, $tb)) {
			$p_al{$tc} ++; 
		}
	}
	for (my $i=3; $i<@ta; $i++) {
		$i == $p1_c and next; 
		$i == $p2_c and next; 
		my $is_good = 1; 
		if ($ta[$i] =~ m!^([ATGC])/([ATGC])$!) {
			$1 ne $2 and $off_cnt{'heteN'} ++; 
		}
		for my $tc (split(/\//, $ta[$i])) {
			$tc eq '.' and do { $is_good = 2; last; }; 
			defined $p_al{$tc} or do { $is_good = 0; last; }; 
		}
		$is_good != 2 and $off_cnt{'typeN'} ++; 
		if ($is_good == 1) {
			$off_cnt{'good'} ++; 
		} elsif ($is_good == 0) {
			$off_cnt{'out'} ++; 
			push(@out_acc, $hh[$i]); 
		} elsif ($is_good == 2) {
			$off_cnt{'miss'} ++; 
		} else {
			die "is_good=$is_good\n"; 
		}
	}
	for my $td (qw/good out miss heteN typeN/) {
		$off_cnt{$td} //= 0; 
	}
	$off_cnt{'heteR'} = ( $off_cnt{'typeN'} == 0 ) ? -1 : sprintf("%0.2f",$off_cnt{'heteN'}/$off_cnt{'typeN'}*100) ; 
	print STDOUT join("\t", @ta[0,1], $ta[$p1_c], $ta[$p2_c], @off_cnt{qw/good out miss heteR/}, join(";;",@out_acc))."\n"; 
}

