#!usr/bin/perl
use strict;
use warnings;
use Bio::DB::Fasta;
use Getopt::Long;
###Author : Xin Wang  03/31/2017


my ($input, $output, $anchor_match_length, $reference, $min_maq, $intron_reads_number, $junction_reads, $keep_tmp, $help);

GetOptions(
        "i|input=s" => \$input,
        "o|output=s" => \$output,
        "a|min_anchor=i" => \$anchor_match_length,  # minimum anchor length for junctions
        "q|maq=i" => \$min_maq,
        "r|reference=s" => \$reference,
        "n|intron_reads=i" => \$intron_reads_number,   ###########
        "j|junction_reads=i" => \$junction_reads,
        "t|tmp!" => \$keep_tmp,
        "h|?|help"		=> \$help,
);



unless ( $input && $output && $reference) {
	&usage;
	exit(1);
}

$anchor_match_length = $anchor_match_length || 10;
$min_maq = $min_maq || 30;
$intron_reads_number = $intron_reads_number || 3;
$junction_reads = $junction_reads || 3;

system "mkdir $output";
system "mkdir $output/tmp";

#my ($in,$ref,$out)=@ARGV;

##### sam2junction ######
open TEMP,">","$output/tmp/check_junction.tmp";
open TEMP1,">","$output/tmp/check_non_junction.tmp";
my %junction_site;
my %non_junction;
my $db = Bio::DB::Fasta -> new ($reference);
open IN,"<",$input;

while(<IN>){
	chomp;
	next if $_ =~ m/^@/;
	my @line = split(/\t/, $_);
	next if $line[1] eq '4';    ##remove the unmapped reads
	next if ($line[4] < $min_maq);
	my $align_pos = $line[3];
	my $chr = $line[2];
	my $cigar = $line[5];
	next if $cigar =~ m/S|H/;	# filter out soft clipping or hard clipping
	my $ed = 0;
	if ($_ =~ m/NM:i:(\d+)/i) { $ed = $1; } else { die "[ERR]in NM $_\n"; }
	my $strand;
	if(/XS:A:\-/){
		$strand = "-";
	}
	elsif(/XS:A:\+/){
		$strand = "+";
	}
	unless($cigar =~ m/N/){
		my @cigar_number = split/\D+/,$cigar;
		my @cigar_type = split/\d+/,$cigar;
		my $cigar_length = @cigar_type;
		my $insertion;
		my $site =0 ;
		my $temp = @cigar_number;
	#print "$cigar_type[0]\n";
		for my $i(1..($cigar_length-1)){
			if($cigar_type[$i] eq "M"){
				#print"$i\t$cigar_number[$i]\n";
				$site = $site + $cigar_number[$i-1];
			}
			elsif($cigar_type[$i] eq "D"){
				$site = $site + $cigar_number[$i-1];
			}
		}
		my $match_start = $align_pos;
		my $match_end = $align_pos + $site -1;
		#print "$match_start\t$match_end\n";
		my $name = $chr."_".$match_start."_".$match_end;
		print TEMP1"$line[0]\t$chr\t$strand\t$match_start\t$match_end\t$cigar\n";
		if(exists $non_junction{$name}){
			$non_junction{$name}{$strand}++;
		}
		else{
			$non_junction{$name}{$strand} = 1;
		}
	}
	##########parse cigar######################
	my @cigar_number = split/\D+/,$cigar;
	my @cigar_type = split/\d+/,$cigar;
	my $cigar_length = @cigar_type;
	my $insertion;
	my $site =0 ;
	my @junction_start;
	my @junction_end;
	my $temp = @cigar_number;
	#print "$cigar_type[0]\n";
	for my $i(1..($cigar_length-1)){
		if($cigar_type[$i] eq "M"){
			#print"$i\t$cigar_number[$i]\n";
			$site = $site + $cigar_number[$i-1];
		}
		elsif($cigar_type[$i] eq "D"){
			$site = $site + $cigar_number[$i-1];
		}
		elsif($cigar_type[$i] eq "N"){
			my $start = $site;
			my $end = $site + $cigar_number[$i-1] -1;
			push(@junction_start,$start);
			push(@junction_end,$end);
			$site = $site + $cigar_number[$i-1];
		}
		
	}
	#print OUT"$line[0]\t$cigar\t@junction_start\t#######@junction_end\n";
	
	my $j_count = @junction_start;
	for my $s(0..($j_count-1)){
		my $j_start = $align_pos + $junction_start[$s];
		my $j_end = $align_pos + $junction_end[$s];
		my $j_start_nt = $db -> subseq($chr,$j_start,$j_start + 1);
		my $j_end_nt = $db -> subseq($chr,$j_end - 1,$j_end);
		if($strand eq "\-"){
			my $temp = $j_start_nt;
			$j_start_nt = &revcom($j_end_nt);
			$j_end_nt = &revcom($temp);
			
		}
		my $name = $chr."_".$j_start."_".$j_end;
		if(exists $junction_site{$name}{$strand}){
			$junction_site{$name}{$strand}++;	
		}
		else{
			$junction_site{$name}{$strand} = 1;
			$junction_site{$name}{$strand} = 1;
					
		}
		print TEMP"$line[0]\t$chr\t$strand\t$j_start\t$j_end\t$j_start_nt\t$j_end_nt\n";
	}
}

close IN;
close TEMP;

open OUT,">","$output/junction.txt";
print OUT"\#name\tstrand\tplus\tminus\ttotal\n";
for my $k(sort {$a cmp $b} keys %junction_site){
	my $total;
	my $sense;
	my $antisense;
	my @strand = keys $junction_site{$k};
	#next if(($junction_site{$k}{"\+"} < $junction_reads) && ($junction_site{$k}{"\-"} < $junction_reads));
	if(exists $junction_site{$k}{"\+"}){
		$sense = $junction_site{$k}{'+'};
	}
	else{
		$sense = 0;
	}
	if(exists $junction_site{$k}{"-"}){
		$antisense = $junction_site{$k}{'-'};
	}
	else{
		$antisense = 0;
	}
	next if(($sense < $junction_reads) && ($antisense < $junction_reads));
	
	$total = $sense + $antisense;
	print OUT"$k\t$strand[0]\t$sense\t$antisense\t$total\n";
}


close OUT;


open OUT,">","$output/non_junction.txt";
print OUT"\#name\tplus\tminus\ttotal\n";
for my $non(sort{$a cmp $b} keys %non_junction){
	my $total;
	my $sense;
	my $antisense;
	if(exists $non_junction{$non}{"\+"}){
		$sense = $non_junction{$non}{'+'};
	}
	else{
		$sense = 0;
	}
	if(exists $non_junction{$non}{"-"}){
		$antisense = $non_junction{$non}{'-'};
	}
	else{
		$antisense = 0;
	}
	$total = $sense + $antisense;
	print OUT"$non\t$sense\t$antisense\t$total\n";
}

close OUT;







my $ir_return = &AS_event_IR(\%junction_site, \%non_junction);

open TEST,">","ir_result.txt";
my %ir_return = %$ir_return;

for my $t(sort {$a cmp $b} keys %ir_return){
	print TEST"$t\t";
	for my $p(sort {$a <=> $b} keys %{$ir_return{$t}}){
		print TEST"$p\:$ir_return{$t}{$p}\t";
	}
	print TEST"\n";
}
close TEST;






sub AS_event_IR{
	my ($junction,$non_junction) = @_;
	###### input files start #####
	my %non_junction = %$non_junction;
	my %junction = %$junction;
	###### input files end #####
	##### returned results hash start #####
	my %ir_return;
	##### returned results hash end #####
	
	##### AS analysis start #####
	for my $j_n(sort {$a cmp $b} keys %junction){
		my @j_infor = split/\_/,$j_n;
		my $j_chr = $j_infor[0];
		my $j_start = $j_infor[1];
		my $j_end = $j_infor[2];
		my @j_strand = keys $junction{$j_n};
		my $j_strand = $j_strand[0];
		for my $nj_n(sort {$a cmp $b} keys %non_junction){
			my @nj_infor = split/\_/,$nj_n;
			my $nj_chr = $nj_infor[0];
			my $nj_start = $nj_infor[1];
			my $nj_end = $nj_infor[2];
			my @nj_strand = keys $non_junction{$nj_n};
			my $nj_strand = $nj_strand[0];
			if(($j_chr eq $nj_chr) && ($j_strand eq $nj_strand)){
				if(($nj_start >= $j_start) && ($nj_start <= $j_end)){
					if(($nj_end >= $j_start) && ($nj_end <= $j_end)){
						print "$j_n\twithin\t$nj_n\n";
						if(exists $ir_return{$j_n}{"within"}){
							$ir_return{$j_n}{"within"} = $ir_return{$j_n}{"within"} + $non_junction{$nj_n}{$nj_strand};	
						}
						else{
							$ir_return{$j_n}{"within"} = $non_junction{$nj_n}{$nj_strand};
						}
					}
					else{
						print "$j_n\tspanned\t$nj_n\n";
						if(exists $ir_return{$j_n}{"spanned"}){
							$ir_return{$j_n}{"spanned"} = $ir_return{$j_n}{"spanned"} + $non_junction{$nj_n}{$nj_strand};
						}
						else{
							$ir_return{$j_n}{"spanned"} = $non_junction{$nj_n}{$nj_strand};
						}
					}
				}
				elsif(($nj_end >= $j_start) && ($nj_end <= $j_end)){
					if(($nj_start >= $j_start) && ($nj_start <= $j_end)){
						print "$j_n\twithin\t$nj_n\n";
						if(exists $ir_return{$j_n}{"within"}){
							$ir_return{$j_n}{"within"} = $ir_return{$j_n}{"within"} + $non_junction{$nj_n}{$nj_strand};	
						}
						else{
							$ir_return{$j_n}{"within"} = $non_junction{$nj_n}{$nj_strand};
						}
					}
					else{
						print "$j_n\tspanned\t$nj_n\n";
						if(exists $ir_return{$j_n}{"spanned"}){
							$ir_return{$j_n}{"spanned"} = $ir_return{$j_n}{"spanned"} + $non_junction{$nj_n}{$nj_strand};
						}
						else{
							$ir_return{$j_n}{"spanned"} = $non_junction{$nj_n}{$nj_strand};
						}
					}
				}
			}
		}
	}
	return \%ir_return;
}








sub revcom {
	my $seq = shift;
	my $rev = reverse $seq;
	(my $revcom = $rev) =~ tr/ACGTacgt/TGCAtgca/;
	return $revcom;
}


sub usage {
print"Usage: $0 -OPTIONS VALUES
This script is used to caculate the SNPindex and plot it on every chromosome.
Options:
     -i  YES  snp file produced by callSNPfromBam script
     -o  YES  prefix of output file name
     -p  YES  prefix of the chromosome name
     -d  NO   minimum read depth required for calling SNP, default = 10
     -w  NO   windowsize to caculate the average SNPindex(megabase unit), default = 3Mb
     -s  NO   step size, should be greater than 0kb and less than windowsize(kilobase unit), default = 500kb 
     \n"
}









