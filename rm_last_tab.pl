#!/usr/bin/perl -w
use strict;

my ($in) = @ARGV;


open IN,"<",$in;

my $maxLenth=16;
my @a = (0..9,'a'..'z','A'..'Z');
my $out = join '', map { $a[int rand @a] } 0..($maxLenth-1);


open OUT,">",$out;
while(<IN>){
	chomp;
	my @line = split/\t/,$_;
	my $length = @line;
	for my $i(0..($length -2)){
		print OUT"$line[$i]\t";
	}
	print OUT"$line[$length -1]\n";
}

close OUT;
system "mv $out $in";