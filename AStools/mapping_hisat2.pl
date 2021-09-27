#!/usr/bin/perl -w
use warnings;
use strict;

my $in=shift;

open IN,"<",$in;
system"mkdir 01_mapping_hisat2";
while(<IN>){
	chomp;
	my $name=$_;
	my $input=$name;
	$name=~s/.merged.fastq.gz//;
	#system "gzip ../02_clean/03_rm_virus/clean/$name.clean.fastq";
	system "mkdir 01_mapping_hisat2/$name";
	system "hisat2 -p 40 --phred33 --dta --rna-strandness R -x /data/xinwang/genome/watermelon/watermelon_v1 -U ../02_clean/03_rm_virus/clean/$name.clean.fastq.gz --min-intronlen 50 --max-intronlen 20000 --ignore-quals --mp 2,2 --rdg 2,1 --rfg 2,1 --score-min C,-4,0 -S 01_mapping_hisat2/$name/$name.sam 2> 01_mapping_hisat2/$name/$name.out ";
	system "samtools view -Sb 01_mapping_hisat2/$name/$name.sam >01_mapping_hisat2/$name/$name.bam";
	system "samtools sort -@ 8 01_mapping_hisat2/$name/$name.bam 01_mapping_hisat2/$name/$name.sorted";
	system "rm 01_mapping_hisat2/$name/$name.sam 01_mapping_hisat2/$name/$name.bam";
}
close IN;
