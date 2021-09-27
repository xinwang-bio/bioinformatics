#!/usr/bin/env perl
use strict;
use warnings;
use Bio::DB::Fasta;
use Getopt::Long;
###Author : Xin Wang  03/31/2017
open TEST,">","log.txt";

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
$intron_reads_number = $intron_reads_number || 10;
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
system "samtools index $input";
system "samtools view $input > $output/tmp/$output.sam";
open IN,"<","$output/tmp/$output.sam";

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
	#if ($_ =~ m/NM:i:(\d+)/i) { $ed = $1; } else { die "[ERR]in NM $_\n"; }
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
			if(($cigar_number[$i] >= $anchor_match_length) && ($cigar_number[$i-2] >= $anchor_match_length)){  ##### anchor matche number filter ($cigar_number[$i] downstream $cigar_number[$i-2] upstream )######
				push(@junction_start,$start);
				push(@junction_end,$end);	
			}
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


##### output the junction site and non-junction reads ######

open OUT,">","$output/junction.txt";
print OUT"\#name\tstrand\tplus\tminus\ttotal\n";
for my $k(sort {$a cmp $b} keys %junction_site){
	my $total;
	my $sense;
	my $antisense;
	my $strand;
	#my @strand = keys $junction_site{$k};
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
	next if(($sense < $junction_reads) && ($antisense < $junction_reads));  ##### junction reads filter ######
	
	$total = $sense + $antisense;
	if(($sense > $antisense) && ($antisense/$total <= 0.1)){
		print OUT"$k\t\+\t$sense\t$antisense\t$total\n";
	}
	elsif(($sense < $antisense) && ($sense/$total <= 0.1)){
		print OUT"$k\t\-\t$sense\t$antisense\t$total\n";
	}
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

##### junction site to hash

my %plus_acceptor;
my %plus_doner;

my %minus_acceptor;
my %minus_doner;

my %graph_postion;
my %graph_type; 

open IN,"<","$output/junction.txt";
while(<IN>){
	chomp;
	next if(/\#/);
	my @line = split/\t/,$_;
	my @infor = split/\_/,$line[0];
	
	if($line[1] eq "+"){
		$plus_doner{$infor[0]}{$infor[1]}{$infor[2]} = $line[4]; 
		$plus_acceptor{$infor[0]}{$infor[2]}{$infor[1]} = $line[4];
		$graph_postion{$infor[0]}{$line[1]}{$infor[1]} = "doner";
		$graph_postion{$infor[0]}{$line[1]}{$infor[2]} = "acceptor";		
	}
	else{
		$minus_doner{$infor[0]}{$infor[2]}{$infor[1]} = $line[4];
		$minus_acceptor{$infor[0]}{$infor[1]}{$infor[2]} = $line[4];
		$graph_postion{$infor[0]}{$line[1]}{$infor[1]} = "acceptor";
		$graph_postion{$infor[0]}{$line[1]}{$infor[2]} = "doner";    ##########$graph_postion{chr}{strand}{postion} = 
	}

}
close IN;


########generate the grafh for each chromosome

my %index;
my %junction_postion;

open GRA,">","$output/$output.graph";
for my $c(keys %graph_postion){
	for my $s(sort {$a cmp $b} keys %{$graph_postion{$c}}){ 	
		print GRA"$c\t$s\t";
		my $index =0;
		my @pos;
		if($s eq "+"){
			@pos = sort {$a <=> $b} keys %{$graph_postion{$c}{$s}};
			$junction_postion{$c}{$s} = [@pos];
		}
		else{
			@pos = sort {$b <=> $a} keys %{$graph_postion{$c}{$s}};
			$junction_postion{$c}{$s} = [@pos];
		}
		for my $p(@pos){
			$index{$c}{$s}{$p} = $index;
			$index++;
			#push(@[$junction_postion{$c}{$s}],$p);
			print GRA"$p\t";
		}
		print GRA"\n";
		print GRA"$c\t$s\t";
		for my $t(@pos){
			print GRA"$graph_postion{$c}{$s}{$t}\t";
		}
		print GRA"\n";
		print GRA"$c\t$s\t";
		for my $i(@pos){
			print GRA"$index{$c}{$s}{$i}\t";
		}
		print GRA"\n";	
	}
}

close GRA;
#####################


my $sort_bam = $input;
my $strand = "\+";

my ($ir_return_plus,$read_infor_return_plus) = &AS_event_IR(\%plus_doner, \%plus_acceptor, $sort_bam, $strand, $min_maq);

my ($aa_return_plus,$aa_multi_return_plus) = &AS_event_AA(\%plus_doner, \%plus_acceptor,\%index);

my ($ad_return_plus,$ad_multi_return_plus) = &AS_event_AD(\%plus_doner, \%plus_acceptor,\%index);

my ($es_return_plus,$es_multi_return_plus) = &AS_event_ES(\%plus_doner, \%plus_acceptor, \%index, \%junction_postion,\%graph_postion,$strand);;



$strand = "\-";

my ($ir_return_minus, $read_infor_return_minus)= &AS_event_IR(\%minus_doner, \%minus_acceptor, $sort_bam, $strand, $min_maq);

my ($aa_return_minus,$aa_multi_return_minus) = &AS_event_AA(\%minus_doner, \%minus_acceptor,\%index);

my ($ad_return_minus,$ad_multi_return_minus) = &AS_event_AD(\%minus_doner, \%minus_acceptor,\%index);

my ($es_return_minus,$es_multi_return_minus) = &AS_event_ES(\%minus_doner, \%minus_acceptor, \%index, \%junction_postion,\%graph_postion,$strand);;


###### output results ######
my ($t,$m,$p,$l,$c,$cc,$n);

###### output IR results start ######
open OUT,">","$output/$output.IR.txt";

my %ir_return; 
my %read_infor;

%ir_return = %$ir_return_plus;
%read_infor = %$read_infor_return_plus;

for $t(sort {$a cmp $b} keys %ir_return){
	for $p(sort {$a <=> $b} keys %{$ir_return{$t}}){
		for $m(sort {$a <=> $b} keys %{$ir_return{$t}{$p}}){
			my $total = 0;
			for $l(sort {$a cmp $b} keys %{$ir_return{$t}{$p}{$m}}){	
				$total = $total + $ir_return{$t}{$p}{$m}{$l};
			}
			if($total >= $intron_reads_number){
				print OUT"$t\t\+\t$p\t$m\t$total\t";
				print TEST"$t\t\+\t$p\t$m\t$total\t";
				for $l(sort {$a cmp $b} keys %{$ir_return{$t}{$p}{$m}}){
					print OUT"$l\:$ir_return{$t}{$p}{$m}{$l}\t";
				}
				my $total_cov = 0;
				my $uniq_cov = 0;
				for $c(sort {$a cmp $b} keys %{$ir_return{$t}{$p}{$m}}){
					my $uniq_num = 0;
					#my $uniq_num = 0;
					if(exists $read_infor{$t}{$p}{$m}{$c}){
						$uniq_num = keys $read_infor{$t}{$p}{$m}{$c};
						print OUT"$c\:$uniq_num\t";
						print TEST"$c\:$uniq_num\t"; 
						for $cc(keys $read_infor{$t}{$p}{$m}{$c}){
							print TEST"$cc\t$read_infor{$t}{$p}{$m}{$c}{$cc}[0]\t$read_infor{$t}{$p}{$m}{$c}{$cc}[1]\n";
							$total_cov = $total_cov + ($read_infor{$t}{$p}{$m}{$c}{$cc}[0]*$read_infor{$t}{$p}{$m}{$c}{$cc}[1]);
							$uniq_cov = $uniq_cov + $read_infor{$t}{$p}{$m}{$c}{$cc}[1];
						}
					}
					else{
						print OUT"$c\:0\t";
					}
				}
				print OUT"$total_cov\t$uniq_cov\t";
				print OUT"\n";	
			}

		}
	}

}

%ir_return = %$ir_return_minus;

for $t(sort {$a cmp $b} keys %ir_return){
	for $p(sort {$a <=> $b} keys %{$ir_return{$t}}){
		for $m(sort {$a <=> $b} keys %{$ir_return{$t}{$p}}){
			my $total = 0;
			for $l(sort {$a cmp $b} keys %{$ir_return{$t}{$p}{$m}}){	
				$total = $total + $ir_return{$t}{$p}{$m}{$l};
			}
			if($total >= $intron_reads_number){
				print OUT"$t\t\-\t$p\t$m\t$total\t";
				for $l(sort {$a cmp $b} keys %{$ir_return{$t}{$p}{$m}}){
					print OUT"$l\:$ir_return{$t}{$p}{$m}{$l}\t";
				}
				print OUT"\n";	
			}

		}
	}

}
close OUT;

###### output IR results end ######


###### output AA results start ######
my %aa_return;
my %aa_multi_return;

%aa_return = %$aa_return_plus;
%aa_multi_return = %$aa_multi_return_plus;
open OUT,">","$output/$output.AA.txt";
open OUT1,">","$output/$output.AA.multi.txt";

for $t(sort {$a cmp $b} keys %aa_return){
	for $p(sort {$a <=> $b} keys %{$aa_return{$t}}){
		print OUT"$t\t\+\t$p\t$aa_return{$t}{$p}[0]\t$aa_return{$t}{$p}[1]\n";
	}
}

for $c(sort {$a cmp $b} keys %aa_multi_return){
	for $cc(sort {$a <=> $b} keys %{$aa_multi_return{$c}}){
		print OUT1"$c\t\+\t$cc\t$aa_multi_return{$c}{$cc}[0]\t$aa_multi_return{$c}{$cc}[1]\t$aa_multi_return{$c}{$cc}[2]\n";
	}
}


%aa_return = %$aa_return_minus;
%aa_multi_return = %$aa_multi_return_minus;

for $t(sort {$a cmp $b} keys %aa_return){
	for $p(sort {$a <=> $b} keys %{$aa_return{$t}}){
		print OUT"$t\t\-\t$p\t$aa_return{$t}{$p}[0]\t$aa_return{$t}{$p}[1]\n";
	}
}

for $c(sort {$a cmp $b} keys %aa_multi_return){
	for $cc(sort {$a <=> $b} keys %{$aa_multi_return{$c}}){
		print OUT1"$c\t\-\t$cc\t$aa_multi_return{$c}{$cc}[0]\t$aa_multi_return{$c}{$cc}[1]\t$aa_multi_return{$c}{$cc}[2]\n";
	}
}

close OUT;
close OUT1;

###### output AA results end ######


###### output AD results start ######
my %ad_return;
my %ad_multi_return;

%ad_return = %$ad_return_plus;
%ad_multi_return = %$ad_multi_return_plus;
open OUT,">","$output/$output.AD.txt";
open OUT1,">","$output/$output.AD.multi.txt";

for $t(sort {$a cmp $b} keys %ad_return){
	for $p(sort {$a <=> $b} keys %{$ad_return{$t}}){
		print OUT"$t\t\+\t$ad_return{$t}{$p}[0]\t$ad_return{$t}{$p}[1]\t$ad_return{$t}{$p}[2]\n";
	}
}

for $c(sort {$a cmp $b} keys %ad_multi_return){
	for $cc(sort {$a <=> $b} keys %{$ad_multi_return{$c}}){
		print OUT1"$c\t\+\t$ad_multi_return{$c}{$cc}[0]\t$cc\t$ad_multi_return{$c}{$cc}[1]\t$ad_multi_return{$c}{$cc}[2]\n";
	}
}


%ad_return = %$ad_return_minus;
%ad_multi_return = %$ad_multi_return_minus;

for $t(sort {$a cmp $b} keys %ad_return){
	for $p(sort {$a <=> $b} keys %{$ad_return{$t}}){
		print OUT"$t\t\-\t$ad_return{$t}{$p}[0]\t$ad_return{$t}{$p}[1]\t$ad_return{$t}{$p}[2]\n";
	}
}

for $c(sort {$a cmp $b} keys %ad_multi_return){
	for $cc(sort {$a <=> $b} keys %{$ad_multi_return{$c}}){
		print OUT1"$c\t\-\t$ad_multi_return{$c}{$cc}[0]\t$cc\t$ad_multi_return{$c}{$cc}[1]\t$ad_multi_return{$c}{$cc}[2]\n";
	}
}


close OUT;
close OUT1;

###### output AD results end ######


###### output ES results start ######
my %es_return;
my %es_multi_return;
open OUT,">","$output/$output.ES.txt";
open OUT1,">","$output/$output.ES.multi.txt";

%es_return = %$es_return_plus;
%es_multi_return = %$es_multi_return_plus;

for $t(sort {$a cmp $b} keys %es_return){
	for $p(sort {$a <=> $b} keys %{$es_return{$t}}){
		print OUT"$t\t\+\t$p\t$es_return{$t}{$p}[0]\t$es_return{$t}{$p}[1]\t$es_return{$t}{$p}[2]\t$es_return{$t}{$p}[3]\t$es_return{$t}{$p}[4]\n";
	}
}

for $c(sort {$a cmp $b} keys %es_multi_return){
	for $m(sort {$a <=> $b} keys %{$es_multi_return{$c}}){
		for $n(sort {$a <=> $b} keys %{$es_multi_return{$c}{$m}}){
			print OUT1"$c\t\+\t$m\t$n\t$es_multi_return{$c}{$m}{$n}\n";
		} 
	}
}


%es_return = %$es_return_minus;
%es_multi_return = %$es_multi_return_minus;

for $t(sort {$a cmp $b} keys %es_return){
	for $p(sort {$a <=> $b} keys %{$es_return{$t}}){
		print OUT"$t\t\-\t$p\t$es_return{$t}{$p}[0]\t$es_return{$t}{$p}[1]\t$es_return{$t}{$p}[2]\t$es_return{$t}{$p}[3]\t$es_return{$t}{$p}[4]\n";
	}
}

for $c(sort {$a cmp $b} keys %es_multi_return){
	for $m(sort {$a <=> $b} keys %{$es_multi_return{$c}}){
		for $n(sort {$a <=> $b} keys %{$es_multi_return{$c}{$m}}){
			print OUT1"$c\t\-\t$m\t$n\t$es_multi_return{$c}{$m}{$n}\n";
		} 
	}
}

close OUT;
close OUT1;


###### output ES results end ######




close TEST;



#########################sub scripts##############################



sub AS_event_IR{
	my ($doner,$acceptor,$sort_bam,$strand,$min_maq) = @_;
	###### input files start #####
	my %doner_acceptor = %$doner;
	my %acceptor_doner = %$acceptor;
	###### input files end #####
	##### returned results hash start #####
	my %ir_return;
	my %read_infor;
	##### returned results hash end #####
	
	##### AS analysis start #####
	for my $j_chr(keys %doner_acceptor){																			#######chromosome
		for my $doner(keys %{$doner_acceptor{$j_chr}}){   														#######doner site	
			my @acceptor = sort {$a <=> $b} keys %{$doner_acceptor{$j_chr}{$doner}};
			my $acceptor_number = @acceptor;
			if($acceptor_number == 1){
				my @doner = keys %{$acceptor_doner{$j_chr}{$acceptor[0]}};
				my $doner_number = keys %{$acceptor_doner{$j_chr}{$acceptor[0]}};
				if(($doner_number ==1)){
					my $region;
					my $j_start;
					my $j_end;
					my $j_strand = $strand;
					if($j_strand eq "+"){
						$j_start = $doner;
						$j_end = $acceptor[0];	
					}
					else{
						$j_start = $acceptor[0];
						$j_end = $doner;
					}
					$ir_return{$j_chr}{$doner}{$acceptor[0]}{"within"} = 0;
					$ir_return{$j_chr}{$doner}{$acceptor[0]}{"spanned"} = 0;
					$region = $j_chr."\:".$j_start."-".$j_end;
					my $samtools_out = `samtools view -q $min_maq $sort_bam $region`;
					my @list = split(/\n/, $samtools_out);
					for my $line(@list){
						chomp $line;
						next if $line =~ m/^@/;
						my @sam = split(/\t/, $line);
						next if $sam[1] eq '4';    ##remove the unmapped reads
						#next if ($line[4] < $min_maq);
						my $nj_n = $sam[0];
						my $align_pos = $sam[3];
						my $nj_chr = $sam[2];
						my $cigar = $sam[5];
						next if $cigar =~ m/S|H/;	# filter out soft clipping or hard clipping
						my $nj_strand;
						if($line =~ /XS:A:\-/){
							$nj_strand = "-";
						}
						elsif($line =~ /XS:A:\+/){
							$nj_strand = "+";
						}
						my $ed =0;
						if ($line =~ m/NM:i:(\d+)/i) { $ed = $1; } else { die "[ERR]in NM $_\n"; }
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
							my $nj_start = $align_pos;
							my $nj_end = $align_pos + $site -1;
							my $nj_name = $nj_start."_".$nj_end."_".$ed;

							if(($j_chr eq $nj_chr) && ($j_strand eq $nj_strand)){
								if(($nj_start >= $j_start) && ($nj_start <= $j_end) && ($nj_end >= $j_start) && ($nj_end <= $j_end)){
									#print "$j_chr\t$j_start\t$j_end\t$j_strand\twithin\t$nj_n\t$nj_strand\t$nj_start\t$nj_end\t$cigar\n";
									if(exists $ir_return{$j_chr}{$doner}{$acceptor[0]}{"within"}){
										$ir_return{$j_chr}{$doner}{$acceptor[0]}{"within"} = $ir_return{$j_chr}{$doner}{$acceptor[0]}{"within"} + 1;	
									}
									else{
										$ir_return{$j_chr}{$doner}{$acceptor[0]}{"within"} = 1;
									}
									if(exists $read_infor{$j_chr}{$doner}{$acceptor[0]}{"within"}{$nj_name}){
										$read_infor{$j_chr}{$doner}{$acceptor[0]}{"within"}{$nj_name}[0] = $read_infor{$j_chr}{$doner}{$acceptor[0]}{"within"}{$nj_name}[0] + 1;
										$read_infor{$j_chr}{$doner}{$acceptor[0]}{"within"}{$nj_name}[1] = $site;
										#print "$nj_name\t$read_infor{$j_chr}{$doner}{$acceptor[0]}{within}{$nj_name}[0]\n";
									}
									else{
										$read_infor{$j_chr}{$doner}{$acceptor[0]}{"within"}{$nj_name}[0] = 1; ###duplicate reads count
										$read_infor{$j_chr}{$doner}{$acceptor[0]}{"within"}{$nj_name}[1] = $site;	####duplicate reads length
									}
								}
								else{
									#print "$j_chr\t$j_start\t$j_end\t$j_strand\tspanned\t$nj_n\t$nj_n\t$nj_strand\t$nj_start\t$nj_end\t$cigar\n";
									if(exists $ir_return{$j_chr}{$doner}{$acceptor[0]}{"spanned"}){
										$ir_return{$j_chr}{$doner}{$acceptor[0]}{"spanned"} = $ir_return{$j_chr}{$doner}{$acceptor[0]}{"spanned"} + 1;
									}
									else{
										$ir_return{$j_chr}{$doner}{$acceptor[0]}{"spanned"} = 1;
									}
									if(exists $read_infor{$j_chr}{$doner}{$acceptor[0]}{"spanned"}{$nj_name}){
										$read_infor{$j_chr}{$doner}{$acceptor[0]}{"spanned"}{$nj_name}[0] = $read_infor{$j_chr}{$doner}{$acceptor[0]}{"spanned"}{$nj_name}[0] + 1;
										$read_infor{$j_chr}{$doner}{$acceptor[0]}{"spanned"}{$nj_name}[1] = $site;
									}
									else{
										$read_infor{$j_chr}{$doner}{$acceptor[0]}{"spanned"}{$nj_name}[0] = 1; ###duplicate reads count
										$read_infor{$j_chr}{$doner}{$acceptor[0]}{"spanned"}{$nj_name}[1] = $site;	####duplicate reads length
									}
								}
							}
						}
					}
				}
			}
		}
	}			#my $name = $i."_".$m;			
	return (\%ir_return,\%read_infor);
}





sub AS_event_AA{																			####### Alternative acceptor scan
	my ($doner,$acceptor,$index,$strand) = @_;
	
	#####input files ----> sub files start##### 
	my %doner_acceptor = %$doner;
	my %acceptor_doner = %$acceptor;
	my %index = %$index;
	#####input files ----> sub files end##### 
	
	#####returned result hash start#####
	my %aa_return;
	my %aa_multi_return;
	#####returned result hash end#####
	
	#####AS analysis start#####
	for my $chr(keys %doner_acceptor){																			#######chromosome
		for my $doner(keys %{$doner_acceptor{$chr}}){   														#######doner site	
			my @acceptor = sort {$a <=> $b} keys %{$doner_acceptor{$chr}{$doner}};
			my $acceptor_number = @acceptor;
	##### AA analysis #####		
			if($acceptor_number == 2){
				#print "$i\t$m\t@acceptor\n";
				my @doner0 = keys %{$acceptor_doner{$chr}{$acceptor[0]}};
				my @doner1 = keys %{$acceptor_doner{$chr}{$acceptor[1]}};
				my $doner0_number = keys %{$acceptor_doner{$chr}{$acceptor[0]}};
				my $doner1_number = keys %{$acceptor_doner{$chr}{$acceptor[1]}};
				if(($doner0_number ==1) && ($doner1_number ==1)){
					#my $name = $i."_".$m;
					my $acceptor_return = join"\,",@acceptor;
					my $reads_count = $acceptor_doner{$chr}{$acceptor[0]}{$doner}."\,".$acceptor_doner{$chr}{$acceptor[1]}{$doner};
					$aa_return{$chr}{$doner}[0] = $acceptor_return;
					$aa_return{$chr}{$doner}[1] = $reads_count;
					#print "$i\t$m\t@acceptor\n";	
				}
			}
	##### multiple AA analysis #####		
			elsif($acceptor_number > 2){
				my $check =0;
				my $multi_acceptor_return =  join"\,",@acceptor;
				my @multi_acceptor_reads;
				for my $a(@acceptor){
					my $doner_number = keys %{$acceptor_doner{$chr}{$a}};
					push(@multi_acceptor_reads,$acceptor_doner{$chr}{$a}{$doner});
					if($doner_number > 1){
						$check =1;
					}
				}
				if($check ==0){
					my $multi_acceptor_reads = join",",@multi_acceptor_reads;
					$aa_multi_return{$chr}{$doner}[0] = $multi_acceptor_return;
					$aa_multi_return{$chr}{$doner}[1] = "Mutil_AA";
					$aa_multi_return{$chr}{$doner}[2] = $multi_acceptor_reads;
				}
			}
		}
	}
	return (\%aa_return, \%aa_multi_return);
}


sub AS_event_AD{
	my ($doner,$acceptor,$index,$strand) = @_;
	
	#####input files ----> sub files start##### 
	my %doner_acceptor = %$doner;
	my %acceptor_doner = %$acceptor;
	my %index = %$index;
	#####input files ----> sub files end##### 
	
	#####returned result hash start#####
	my %ad_return;
	my %ad_multi_return;
	#####returned result hash end#####
	
	#####AS analysis start#####
	for my $chr(keys %acceptor_doner){																				#######chromosome
		for my $acceptor (keys %{$acceptor_doner{$chr}}){   														#######doner site	
			my @doner = sort {$a <=> $b} keys %{$acceptor_doner{$chr}{$acceptor}};
			my $doner_number = @doner;
	##### AD analysis ######		
			if($doner_number ==2){
				my @acceptor0 = keys %{$doner_acceptor{$chr}{$doner[0]}};
				my @acceptor1 = keys %{$doner_acceptor{$chr}{$doner[1]}};
				my $acceptor0_number = keys %{$doner_acceptor{$chr}{$doner[0]}};
				my $acceptor1_number = keys %{$doner_acceptor{$chr}{$doner[1]}};
				if(($acceptor0_number ==1) && ($acceptor1_number ==1)){
					#my $name = $i."_".$m;
					my $doner_return = join"\,",@doner;
					my $reads_count = $doner_acceptor{$chr}{$doner[0]}{$acceptor}."\,".$doner_acceptor{$chr}{$doner[1]}{$acceptor};
					$ad_return{$chr}{$acceptor}[0] = $doner_return;
					$ad_return{$chr}{$acceptor}[1] = $acceptor;
					$ad_return{$chr}{$acceptor}[2] = $reads_count;
					#print "$i\t$m\t@acceptor\n";	
				}
			}
	##### multiple AD analysis #####		
			elsif($doner_number > 2){
				my $check =0;
				my $multi_doner_return =  join"\,",@doner;
				my @multi_doner_reads;
				for my $d(@doner){
					my $acceptor_number = keys %{$doner_acceptor{$chr}{$d}};
					push(@multi_doner_reads,$doner_acceptor{$chr}{$d}{$acceptor});
					if($acceptor_number > 1){
						$check =1;
					}
				}
				if($check ==0){
					my $multi_acceptor_reads = join",",@multi_doner_reads;
					$ad_multi_return{$chr}{$acceptor}[0] = $multi_doner_return;
					$ad_multi_return{$chr}{$acceptor}[1] = "Mutil_DD";
					$ad_multi_return{$chr}{$acceptor}[2] = $multi_acceptor_reads;
				}
			}
		}
	}
	return (\%ad_return,\%ad_multi_return);
}




sub AS_event_ES{																			#######only for one exon skipping
	my ($doner,$acceptor,$index,$junction_postion,$graph_postion,$strand) = @_;
	
	##### input files --------> sub files start ######
	my %doner_acceptor = %$doner;
	my %acceptor_doner = %$acceptor;
	my %junction_postion = %$junction_postion;
	my %graph_postion = %$graph_postion;
	my %index = %$index;
	##### input files --------> sub files end ######
	
	##### returned results hash start #####
	my %es_return;
	my %es_multi_return;
	##### returned results hash end #####
	
	##### AS analysis start ######
	for my $chr(keys %doner_acceptor){																			#######chromosome
		for my $doner (keys %{$doner_acceptor{$chr}}){   														#######doner site	
			my @acceptor = sort {$a <=> $b} keys %{$doner_acceptor{$chr}{$doner}};
			my $acceptor_number = @acceptor;
			if($acceptor_number ==2){
				my $ds1 = $doner;
				my $as4;
				my $as2;
				my $ds3;
	##### search the last acceptor ######				
				my $delta1 = abs($acceptor[0] - $doner);
				my $delta2 = abs($acceptor[1] - $doner);
				if($delta1 > $delta2){
					$as2 = $acceptor[1];
					$as4 = $acceptor[0];
				}
				else{
					$as2 = $acceptor[0];
					$as4 = $acceptor[1];
				}
		##### last acceptor decided #####	
				my $delta_index = abs($index{$chr}{$strand}{$as4} - $index{$chr}{$strand}{$ds1}); ##### how many junction site between ds1 and as4. junction site number + 1 = delta_index 
				my @as4_doner = keys %{$acceptor_doner{$chr}{$as4}};
				my $as4_doner = keys %{$acceptor_doner{$chr}{$as4}};
				if($as4_doner ==2){
					for my $as_d(@as4_doner){
						if($as_d != $ds1){
							$ds3 = $as_d;
							my $as2_doner = keys %{$acceptor_doner{$chr}{$as2}};
							my $ds3_acceptor = keys %{$doner_acceptor{$chr}{$ds3}};
		##### only one exon skipping analysis #####
							if(($delta_index ==3) && ($as2_doner ==1) && ($ds3_acceptor ==1)){
								#print "$count_acceptor\t$as2_doner\t$ds3_acceptor\t$as4_doner\n";
								#print "$i\t$strand\t$ds1\t$as2\t$ds3\t$as4\t@acceptor\t@as4_doner\n";	
								my $ds1_acceptor = $as2."\,".$as4;
								my $ds1_reads = $doner_acceptor{$chr}{$ds1}{$as2}."\,".$doner_acceptor{$chr}{$ds1}{$as4};
								my $as4_doner = $ds1."\,".$ds3;
								my $as4_reads = $acceptor_doner{$chr}{$as4}{$ds1}."\,".$acceptor_doner{$chr}{$as4}{$ds3};
								$es_return{$chr}{$ds1}[0] = $ds1_acceptor; 
								$es_return{$chr}{$ds1}[1] = $ds1_reads;
								$es_return{$chr}{$ds1}[2] = $as4;
								$es_return{$chr}{$ds1}[3] = $as4_doner; 
								$es_return{$chr}{$ds1}[4] = $as4_reads;
							}
			##### for multi exon skipping ######				
							elsif(($delta_index > 3) && ($delta_index %2 ==1)){
								my $total_score =0;
								my $score =0;
								#print "$chr\t$strand\t$ds1\t$as4\n";
								$es_multi_return{$chr}{$ds1}{$as4} = "Multi_ES";
								#for my $index_check(0..($delta_index-2)){
									#$total_score++;
									#print "$total_score\t";
									#if (($index_check %2 ==0) && ($graph_postion{$chr}{$strand}{$junction_postion{$chr}{$strand}[$index{$chr}{$strand}{$doner} +1 +$index_check]} eq "acceptor")){            ######### acceptor
									#	$score++;
									#}
									#elsif(($index_check %2 ==0) && ($graph_postion{$chr}{$strand}{$junction_postion{$chr}{$strand}[$index{$chr}{$strand}{$doner} +1 +$index_check]} eq "doner")){								 ######### doner
									#	$score++;
									#}
									#print "$score\n";
								#}
								#if($total_score == $score){	
								#}
							}
						}
					}
				}	
			}
		}
	}
	return (\%es_return,\%es_multi_return);
}


sub revcom {
	my $seq = shift;
	my $rev = reverse $seq;
	(my $revcom = $rev) =~ tr/ACGTacgt/TGCAtgca/;
	return $revcom;
}


sub usage {
print"Usage: $0 -OPTIONS VALUES
This script is used to scan the alternative splicing (AS) event from specific strand RNA-seq.

Version: 0.1beta.

Author : Xin Wang. 
Email: wangxinbio\@gmail.com || xw96\@cornell.edu

Options:
     -i|input  YES  The sorted bam file.
     -o  YES  output folder
     -r  YES  genome file (fasta format)
     -q  NO   minumum mapping quality for scaning AS
     -a  NO   minimum anchor for junction reads for junction site scan, default = 10 bp
     -n  NO   minimum reads number cover the intron region, default = 10
     -j  NO   minimum junction coverage, default = 3
     \n"
}









