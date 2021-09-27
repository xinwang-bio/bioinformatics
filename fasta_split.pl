#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Bio::SeqIO;
my ($help,$input,$output,$force,$overwrite);
my $reads_end_turn;
GetOptions(
	'help|h!' => \$help,
	'input|i=s' => \$input,
	'output|o=s' => \$output,
	'force|f!' => \$force,
	'overwrite|r!' => \$overwrite,
);

my $usage = "
fasta_split.pl:
This program can split the fasta sequecne into files that contain one fasta.
Author: Xin Wang  Email: wangxinbio\@gmail.com

Usage:
	perl fasta_split.pl [options] input.fasta
Options:
	-h|--help         print this help message.
	-i|--input        inputfile for split (required).
	-o|--output       output document.		[default: output] 
        -f|--force        ignore if the output document is existed.
        -r|--overwrite    overwrite the output document. WARNING: will delete the files already exist in it.
  
";
unless ($input) {
	print "$usage";
	exit;
}
$input=$input;
if($output){
	$output=$output;
}
else{
	$output="output";
}


if(-e $output){
	if($force){
		
	}
	elsif($overwrite){
		system"rm -rf $output";
		system "mkdir $output";
	}
	else{
		print "WARNING: The $output document is existed !!!!!\n";
		exit();
	}
		
}
else{
	system"mkdir $output";
}

my $in=Bio::SeqIO -> new (-file => "<$input",-format =>"fasta");

while(my $inseq = $in -> next_seq){
	my $id = $inseq ->id;
	my $seq_out = Bio::SeqIO->new(-file   => ">$output/$id.fa",-format => "fasta");	
	$seq_out->write_seq($inseq);
}

