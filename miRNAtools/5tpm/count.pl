#!usr/bin/perl -w
use strict;

my($in)=@ARGV;
my %ha;

open IN,"<",$in;
while(<IN>){
	chomp;
	my @line=split/\t/,$_;
	if(exists $ha{$line[1]}){
		$ha{$line[1]}++;
	}
	else{
		$ha{$line[1]} = 1;
	}
}
close IN;


for my $i(keys %ha){
	print "$i\t$ha{$i}\n";
}
