#!usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
#SL2.40ch12:18056736-18063004

my ($infile,$chr,$start,$end,$strand) = @ARGV;
my $out = $chr."_".$start."_".$end."_".$strand;

my $result ="$out.fa";
open OUT,">",$result;
my $db = Bio::DB::Fasta -> new ($infile);
	my $seq;
	my $na = $chr;#14026758-14026957,14028250-14028650 regin:14026000 14038000
	my $lian = $strand;
	my $s = $start;#14026758-2000;  #43349267;43349396 15797714-15797042
	my $e = $end;
	if ($lian eq "+"){
		$seq =$db -> subseq("$na",$s,$e);
		print OUT">";
		print OUT"$na\t$s-$e$lian\n";
		print OUT"$seq\n";
	}
	if ($lian eq "-"){
		$seq =$db -> subseq("$na",$s,$e);
		$seq = &revcom($seq);
		print OUT">";
		print OUT"$na\t$s-$e$lian\n";
		print OUT"$seq\n";
	}
	
	
sub revcom {
	my $seq = shift;
	my $rev = reverse $seq;
	(my $revcom = $rev) =~ tr/ACGTacgt/TGCAtgca/;
	return $revcom;
}
	