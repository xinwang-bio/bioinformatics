#!/usr/local/bin/perl -w

use strict;
use Bio::SeqIO;
use Bio::DB::Fasta;
use 5.010;

my ($fasta, $out) = @ARGV;
my $tmp = "$out.tmp";
my $file = Bio::SeqIO -> new(-format => "fasta",-file => $fasta);
my $outfile = Bio::SeqIO -> new (-format => "fasta", -file => ">$tmp");

while (my $inseq = $file->next_seq) {
    $outfile->write_seq($inseq);
}




my $db = Bio::DB::Fasta -> new($tmp);
open OUT, '>', $out;

# print head line;
my @ids = $db -> ids;
my $length = $db -> length($ids[0]);
say OUT " ".($#ids+1)." $length";

for my $l (0..int($length/50)) {
	for my $name (@ids) {
		# print first sequence line
		if ($l == 0) {
			printf OUT "%-12s", $name;
		}
		else {
			print OUT ' ' x 12;
		}
		if ($l < int($length/50)) {
			for my $i (1..5) {
				print OUT ' '.($db->subseq($name, $l*50+($i-1)*10+1, $l*50+$i*10));
			}
		}
		else {
			my $last = $length % 50;
			if ($last != 0) {
				for my $j (1..int($last/10)) {
					print OUT ' '.($db->subseq($name, $l*50+($j-1)*10+1, $l*50+$j*10));
				}
				print OUT ' '.($db->subseq($name, 10*int($length/10)+1, $length)) if $length % 10 != 0;
			}
		}
		print OUT "\n";
	}
	print OUT "\n";
}
close OUT;
system "rm $tmp"
