#!/usr/bin/perl -w

use strict;
use Bio::DB::Fasta;


my($siR_seq, $ref, $out) = @ARGV;

my $cutoff = 20;

system "mkdir $out";

system "bowtie -v 0 -k 21 -S -f $ref $siR_seq $out/sRNA.map.genome.sam";

&filter_match("$out/sRNA.map.genome.sam","$out/sRNA.filter.sam",$cutoff);

&get_genome_regions("200", "$out/sRNA.filter.sam", "$ref", "$out/sR_precursor_candidate.fa");




sub filter_match
{
    my ($genome_match_sam, $filter_match_sam, $cutoff) = @_;

    my %num_hit;
    open IN,"<",$genome_match_sam || die $!;
    while(<IN>){
    	next if $_ =~ m/^@/;
    	my @a = split(/\t/, $_);
    	next if $a[1] & 0x4;
    	if(exists $num_hit{$a[0]}){
    		$num_hit{$a[0]}++
    	}
    	else{
    		$num_hit{$a[0]} = 1;
    	}
    }
    close IN;

	# filter the result and output
	open OUT,">",$filter_match_sam || die $!;
	open IN,"<",$genome_match_sam  || die $!;
    while(<IN>) {
		print OUT "$_" and next if $_ =~ m/^@/;
		my @a = split(/\t/, $_);
		next if $a[1] & 0x4;
		die "[ERR]undef sRNA $a[0]\n" unless defined $num_hit{$a[0]};
		print OUT $_ if $num_hit{$a[0]} < $cutoff;
	}
	close IN;
	close OUT;
}



sub get_genome_regions{
	
	my ($flanking_len, $srna_map_sam, $genome_seq, $precursor_candidate) = @_;	
	my $genome = Bio::DB::Fasta -> new ($genome_seq);
	
	open OUT,">",$precursor_candidate || die $!;
	open IN,"<",$srna_map_sam || die $!;
	while(<IN>) {
		chomp;
		next if $_ =~ m/^@/;
		my @line = split(/\t/, $_);
		next if $line[1] & 0x4;
		my ($sid, $seq, $chr, $flag, $start, $end, $hit_seq) = ($line[0], $line[9], $line[2], $line[1], $line[3], $line[3]+length($line[9])-1, $line[9]);
		my $sub_start = $start-$flanking_len;
		my $sub_end   = $end+$flanking_len;
		my $sub_seq = $genome -> subseq($chr, $sub_start, $sub_end);
		my $strand;
		my $new_id;
		if ($flag == 0) { 
			$strand = '+'; 
			$new_id = $sid."_".$chr."_".$start.$strand."\t".$chr."\t".$sub_start."\t".$sub_end."\t".$seq."\t$flanking_len\t".($flanking_len + length($line[9]) - 1);
		}
		elsif ($flag == 16) { 
			$strand = '-'; 
			$seq = revcom($seq);
			$sub_seq = revcom($sub_seq);
			$new_id = $sid."_".$chr."_".$end.$strand."\t".$chr."\t".$sub_end."\t".$sub_start."\t".$seq."\t$flanking_len\t".($flanking_len + length($line[9]) - 1);	
		}
		#$seq =~s/n/N/g;
		print OUT ">".$new_id."\n".$sub_seq."\n";
	}
	close IN;
	close OUT;
}










	
sub revcom {
	my $seq = shift;
	my $rev = reverse $seq;
	(my $revcom = $rev) =~ tr/ACGTacgt/TGCAtgca/;
	return $revcom;
}
