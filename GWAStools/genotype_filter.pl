#!usr/bin/perl

use strict;
use File::Basename;

#get the reference and remove missing data > 80 % and only two allele postion is left

#Edit by WangXin,9/11/2013

my ($in,$out)=@ARGV;

# $long how many rows in one line(individual numbers)
#open TE,">","D:\\test\\test.txt";
my $long;

open IN,"<",$in;
open OUT,">",$out;

my $label = basename($out);
my $path = fileparse($out);
my $dirname = dirname($out);

#generate a information file named **.info
open INFO,">","$dirname/$label.info";

#calculate each postion frequence and missing data
while(<IN>){
	chomp;
	my @line =split;
	my %snp;
	my %allele;
	my $ch=$line[0];
	my $pos=$line[1];
	#print"$seq\n";
	$long=@line;
	$allele{$ch}{$pos}{"-"}=0;
	for my $k(2..($long-1)){
		if($line[$k] eq "M"){
			$allele{$ch}{$pos}{"A"}++;
			$allele{$ch}{$pos}{"C"}++;	
		}
		elsif($line[$k] eq "R"){
			$allele{$ch}{$pos}{"A"}++;
			$allele{$ch}{$pos}{"G"}++;	
		}
		elsif($line[$k] eq "W"){
			$allele{$ch}{$pos}{"A"}++;
			$allele{$ch}{$pos}{"T"}++;
		}
		elsif($line[$k] eq "S"){
			$allele{$ch}{$pos}{"C"}++;
			$allele{$ch}{$pos}{"G"}++;						
		}
		elsif($line[$k] eq "Y"){	
			$allele{$ch}{$pos}{"T"}++;
			$allele{$ch}{$pos}{"C"}++;
		}
		elsif($line[$k] eq "K"){
			$allele{$ch}{$pos}{"T"}++;
			$allele{$ch}{$pos}{"G"}++;	
		}
		elsif($line[$k] eq "-"){
			$allele{$ch}{$pos}{"-"}=$allele{$ch}{$pos}{"-"}+2;
		}	
		else{
			$allele{$ch}{$pos}{$line[$k]}=$allele{$ch}{$pos}{$line[$k]}+2;			
		}
	}	
#output the information of the postion	

	print INFO"$ch\t$pos\t";
	for my $m (sort {$a cmp $b} keys %{$allele{$ch}{$pos}}){
		my $frequ= $allele{$ch}{$pos}{$m}/(($long-2)*2);
		print INFO"$m\t$frequ\t";
		$snp{$ch}{$pos}{$m}=$frequ;
	}
	print INFO"\n";		
	
#sort the frenqu in this postion

	
	
#remove missing data > 80% ,remove maf of allele < 0.05 ,only two allele in one postion is left
	my @count;
	for my $n (keys %{$snp{$ch}{$pos}}){
		if($n ne "-"){
			push (@count,$snp{$ch}{$pos}{$n});	
		}
	}
	@count=sort{$a <=> $b} @count;
	my $allele_length = @count;
	my $maf=shift @count;
	#print TE"@count\t$maf\t$allele_length\n";
	if($allele{$ch}{$pos}{"-"} < ($long-2)*2*0.6 && $allele_length ==2 && $maf >= 0.05){
		@line=join"\t",@line;
		print OUT"@line\n";	
	}		
}	
print"$label filter is finished\n";
close IN;	
close OUT;
close INFO;


