#!/usr/bin/perl -w
use strict;



my ($yi,$xin) = @ARGV;

my %total;
my %yi;

open IN,"<",$yi;
while(<IN>){
	chomp;
	my @line = split;
	my @pos = split/\#/,$line[0];
	$yi{$pos[0]}{$pos[1]}{$pos[2]} = $_;
	$total{$pos[0]}{$pos[1]}{$pos[2]} =0;
}
close IN;

my %xin;

open IN,"<",$xin;
while(<IN>){
	chomp;
	my @line = split;
	next if(/\#/);
	my @pos = split/\_/,$line[0];
	$xin{$pos[0]}{$pos[1]}{$pos[2]} = $_;
	$total{$pos[0]}{$pos[1]}{$pos[2]} =0;
}
close IN;

open XIN,">","xin_spec.txt";

open YI,">","yi_spec.txt";

open OUT,">","xin_and_yi.txt";

for my $i(sort{$a cmp $b} keys %total){
	for my $m(sort{$a <=> $b} keys %{$total{$i}}){
		for my $n(sort{$a <=> $b} keys %{$total{$i}{$m}}){
			if((exists $xin{$i}{$m}{$n}) && (exists $yi{$i}{$m}{$n})){
				print OUT"$i\t$m\t$n\t$xin{$i}{$m}{$n}\t\#\t$yi{$i}{$m}{$n}\n";
			}
			elsif((exists $xin{$i}{$m}{$n}) && (!exists $yi{$i}{$m}{$n})){
				print XIN"$i\t$m\t$n\t$xin{$i}{$m}{$n}\n";
			}
			elsif((!exists $xin{$i}{$m}{$n}) && (exists $yi{$i}{$m}{$n})){
				print YI"$i\t$m\t$n\t$yi{$i}{$m}{$n}\n";
			}
		}
	}
}
close OUT;
close YI;
close XIN;









