#!/usr/bin/env perl
use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;
my ($input,$help);
GetOptions(
	'help|h!' => \$help,
	'input|i=s' => \$input,
);

my $usage = "
length_fasta.pl:
This program can stat the fasta sequecne length.
Author: Xin Wang  Email: wangxinbio\@gmail.com

Usage:
	perl length_fasta.pl input.fasta
Options:
	-h|--help         print this help message.
	-i|--input        inputfile (required, format should be fasta).
  
";
unless ($input) {
	print "input file is not exist!!\n";
	print "$usage";
	exit;
}
$input=$input;
my $total =0 ;
my $total_non_N = 0;

my $seq_in =Bio::SeqIO->new(-format =>"fasta",-file =>$input);
print "\#chr\tATCGN\tATCG\n";								
while( my $seq =$seq_in ->next_seq()){
	my $chrom = $seq -> id;
	my $length = $seq -> length;
	my $seq_nt = $seq -> seq;
	$total = $total + $length;
	$seq_nt =~s/N//ig;
	my $atcg = length($seq_nt);
	$total_non_N = $total_non_N + $atcg;
	print"$chrom\t$length\t$atcg\n";
}
#print "Total\t$total\t$total_non_N\n";
