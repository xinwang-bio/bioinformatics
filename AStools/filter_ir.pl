#!/usr/bin/perl -w
use strict;

my ($in,$out) = @ARGV;

my %doner;
my %acceptor;

open IN,"<",$in;
while(<IN>){
	chomp;
	my @line = split/\t/,$_;
#	my @infor = split/\_/,$line[0];
	my $length = @line;
	my $total =0;
	for my $t(3..($length -1)){
		$line[$t] =~s/spanned\://;
		$line[$t] =~s/within\://;
		$total += $line[$t];
	}
	#print "$infor[0]\t$infor[1]\t$infor[2]\t$total\n";
	$doner{$line[0]}{$line[1]}{$line[2]} = $total;
	$acceptor{$line[0]}{$line[2]}{$line[1]} = $total;
}

close IN;

open OUT,">",$out;

for my $i(keys %doner){
	for my $m(keys %{$doner{$i}}){
		my @c = keys %{$doner{$i}{$m}};
		my $c = keys %{$doner{$i}{$m}};
		my @a = keys %{$acceptor{$i}{$c[0]}};
		my $a = keys %{$acceptor{$i}{$c[0]}};
		#print "@c\t$c\t$a\n";
		if(($a ==1) && ($c ==1) && ($doner{$i}{$m}{$c[0]} >=10)){
			print OUT"$i\t$m\t$c[0]\t\t$doner{$i}{$m}{$c[0]}\n";
		}	
	}
}
close OUT;
