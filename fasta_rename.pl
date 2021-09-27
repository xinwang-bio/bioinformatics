#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;


my ($fasta,$out,$prefix)=@ARGV;

my $l =0;

my $seq_in = Bio::SeqIO->new(
                             -file   => "<$fasta",
                             -format => "fasta",
                             );
open OUT,">","$out.tmp";
open NUM,">","$out.id";

while (my $inseq = $seq_in->next_seq) {
        my $id= $inseq ->id;
        my $seq = $inseq -> seq;
        print OUT">$prefix\_$l\n$seq\n";
        print NUM"$prefix\_$l\t$id\n";
        $l++;
}
close OUT;
close NUM;

my $file = Bio::SeqIO -> new(-format => "fasta",-file => "$out.tmp");
my $outfile = Bio::SeqIO -> new (-format => "fasta", -file => ">$out");

while (my $inseq = $file->next_seq) {
    $outfile->write_seq($inseq);
}
system "rm $out.tmp";


