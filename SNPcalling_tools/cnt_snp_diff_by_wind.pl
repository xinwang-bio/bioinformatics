#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 bin_size p1_col p2_col offsprint_col in_comb.snp\n"; 

my $binL = shift; 
my $c1 = shift; 
my $c2 = shift; 
my $c3 = shift; 

my %miss = qw(
N    1
./.  1
.    1
); 

my ($id1, $id2, $id3); 
my %cnt; # P1/P2/Hete/Miss
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1) {
		( $ta[0] eq 'chr' and $ta[1] eq 'pos' ) or die "$_\n"; 
		$id1 = $ta[$c1]; 
		$id2 = $ta[$c2]; 
		$id3 = $ta[$c3]; 
		next; 
	}
	my $binID = int(($ta[1]-1)/$binL); 
	my ($a1, $a2, $a3) = @ta[$c1, $c2, $c3]; 
	$a1 =~ s!^([ATGC*])/\1$!$1! or next; 
	$a2 =~ s!^([ATGC*])/\1$!$1! or next; 
	if (defined $miss{$a3}) {
		$cnt{$ta[0]}{$binID}{'Miss'} ++; 
	} elsif ($a3 eq "$a1/$a1") {
		$cnt{$ta[0]}{$binID}{'P1'} ++; 
	} elsif ( $a3 eq "$a2/$a2" ) {
		$cnt{$ta[0]}{$binID}{'P2'} ++; 
	} elsif ( $a3 eq "$a1/$a2" or $a3 eq "$a2/$a1" ) {
		$cnt{$ta[0]}{$binID}{'Hete'} ++; 
	} else {
		&tsmsg("[Wrn] Skip bad genotype [$ta[0] $ta[1] $a3] : $_\n"); 
		next; 
	}
}

print STDOUT join("\t", qw/chr binID wind_S wind_E P1_N P2_N Hete_N Miss_N sumN/)."\n"; 
for my $chrID (sort keys %cnt) {
	for my $binID (sort {$a<=>$b} keys %{$cnt{$chrID}}) {
		my %v = %{ $cnt{$chrID}{$binID} }; 
		my $cutP1 = $binID*$binL+1; 
		my $cutP2 = ($binID+1) * $binL; 
		my $sum = 0; 
		for my $tt (qw/P1 P2 Hete Miss/) {
			$v{$tt} //= 0; 
			$sum += $v{$tt}; 
		}
		$sum - $v{'Miss'} > 0 or next; 
		print STDOUT join("\t", $chrID, $binID, $cutP1, $cutP2, @v{qw/P1 P2 Hete Miss/}, $sum)."\n"; 
	}
}
