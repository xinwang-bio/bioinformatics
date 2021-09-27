#!/use/bin/perl -w
use strict;

my ($yi,$wx) = @ARGV;



open YI,">","yi.fa";
open WX,">","wx.fa";


open IN,"<",$yi;
while(<IN>){
	chomp;
	if(/>/){
		print YI"$_\n";
	}
	else{
		$_ =~ s/n/N/g;
		print YI"$_\n";
	}
}
close IN;
close YI;


open IN,"<",$wx;
while(<IN>){
	chomp;
	if(/>/){
		print WX"$_\n";
	}
	else{
		$_ =~ s/n/N/g;
		print WX"$_\n";
	}
}
close IN;
close WX;

