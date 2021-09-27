#!/usr/bin/env perl
use warnings;
use strict;

#open OUT,">$ARGV[1]";

my %hash =(
        "0,1-2^" => "ES",
        "0,1^2-" => "IR",
        "1-,2-" => "AA",
        "1^,2^" => "AD",
);
#open OUT,">","test.txt";
for my $i(@ARGV){
	my ($es,$ir,$aa,$ad,$oth) = (0,0,0,0,0);
	open IN,"$i";
	while(<IN>){
       my ($key)=$_=~ /structure "(\S+)";/;
       if(exists $hash{$key}){
       	if($hash{$key} eq "ES"){
       		$es++;
       	}elsif($hash{$key} eq "IR"){
       		$ir++;
       	}elsif($hash{$key} eq "AA"){
       		$aa++;
       	}elsif($hash{$key} eq "AD"){
       		$ad++;
       	}
       }
       else{
       	$oth++;
       }
	}
	close IN;
	print"#####################$i########################\n";
	print "ES\t$es\nIR\t$ir\nAA\t$aa\nAD\t$ad\nOther\t$oth\n";
}


=p
my ($es,$ir,$aa,$ad,$oth) = (0,0,0,0,0);
while(<IN>){
                my ($key)=$_=~ /structure "(\S+)";/;
                if(exists $hash{$key}){
                        if($hash{$key} eq "ES"){
                                $es++;
                        }elsif($hash{$key} eq "IR"){
                                $ir++;
                        }elsif($hash{$key} eq "AA"){
                                $aa++;
                        }elsif($hash{$key} eq "AD"){
                                $ad++;
                        }
                }else{
                                $oth++;
			#	print OUT"$key\n";
                }
                                
}
print "ES\t$es\nIR\t$ir\nAA\t$aa\nAD\t$ad\nOther\t$oth\n";
