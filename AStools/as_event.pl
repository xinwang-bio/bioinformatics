#!/usr/bin/perl -w
use strict;

my ($in,$out) = @ARGV;


my %plus_acceptor;
my %plus_doner;

my %minus_acceptor;
my %minus_doner;

my %graph_postion;
my %graph_type; 
open IN,"<",$in;
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

open GRA,">","graph_result.txt";
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



my $sort_bam = "pf10-1.tophat.bam";

###### filter junction site, only one doner ------> one acceptor (one to one) were used to identified for intron retention (IR)

my $strand = "\+";
my $ir_return = &AS_event_IR(\%plus_doner, \%plus_acceptor, $sort_bam, $strand);

open TEST,">","ir_result0419.txt";
my %ir_return = %$ir_return;

for my $t(sort {$a cmp $b} keys %ir_return){
	for my $p(sort {$a <=> $b} keys %{$ir_return{$t}}){
		for my $m(sort {$a <=> $b} keys %{$ir_return{$t}{$p}}){
			print TEST"$t\t$p\t$m\t";
			for my $l(sort {$a cmp $b} keys %{$ir_return{$t}{$p}{$m}}){	
				print TEST"$l\:$ir_return{$t}{$p}{$m}{$l}\t";
			}	
		}
		print TEST"\n";
	}

}
close TEST;




sub AS_event_IR{
	my ($doner,$acceptor,$sort_bam,$strand) = @_;
	###### input files start #####
	my %doner_acceptor = %$doner;
	my %acceptor_doner = %$acceptor;
	###### input files end #####
	##### returned results hash start #####
	my %ir_return;
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
					if($j_strand eq "\+"){
						$j_start = $doner;
						$j_end = $acceptor[0];	
					}
					else{
						$j_start = $acceptor[0];
						$j_end = $doner;
					}
					$region = $j_chr."\:".$j_start."-".$j_end;
					my $samtools_out = `samtools view $sort_bam $region`;
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
							if(($j_chr eq $nj_chr) && ($j_strand eq $nj_strand)){
								if(($nj_start >= $j_start) && ($nj_start <= $j_end) && ($nj_end >= $j_start) && ($nj_end <= $j_end)){
									print "$j_chr\t$j_start\t$j_end\t$j_strand\twithin\t$nj_n\t$nj_strand\t$nj_start\t$nj_end\t$cigar\n";
									if(exists $ir_return{$j_chr}{$doner}{$acceptor[0]}{"within"}){
										$ir_return{$j_chr}{$doner}{$acceptor[0]}{"within"} = $ir_return{$j_chr}{$doner}{$acceptor[0]}{"within"} + 1;	
									}
									else{
										$ir_return{$j_chr}{$doner}{$acceptor[0]}{"within"} = 1;
									}
								}
								else{
									print "$j_chr\t$j_start\t$j_end\t$j_strand\tspanned\t$nj_n\t$nj_n\t$nj_strand\t$nj_start\t$nj_end\t$cigar\n";
									if(exists $ir_return{$j_chr}{$doner}{$acceptor[0]}{"spanned"}){
										$ir_return{$j_chr}{$doner}{$acceptor[0]}{"spanned"} = $ir_return{$j_chr}{$doner}{$acceptor[0]}{"spanned"} + 1;
									}
									else{
										$ir_return{$j_chr}{$doner}{$acceptor[0]}{"spanned"} = 1;
									}
								}
							}
						}
					}
				}
			}
		}
	}			#my $name = $i."_".$m;			
	return \%ir_return;
}



=p
open OUT,">",$out;

#my %as_doner = &AS_acceptor(\%plus_doner, \%plus_acceptor);
#my %as_doner = &AS_acceptor(\%minus_doner, \%minus_acceptor);

#my %as_doner = &AS_doner(\%plus_doner, \%plus_acceptor);
#my %as_doner = &AS_doner(\%minus_doner, \%minus_acceptor);

for my $i(sort {$a cmp $b} keys %as_doner){
	for my $m(sort {$a <=> $b} keys %{$as_doner{$i}}){
		print OUT"$i\t$m\t$as_doner{$i}{$m}[0]\t\[$as_doner{$i}{$m}[1]\]\t\[$as_doner{$i}{$m}[2]\]\n";	
	}
}


close IN;
close OUT;
=cut
#my ($aa_return,$aa_multi_return) = &AS_event_AD(\%plus_doner, \%plus_acceptor,\%index);


=p
my ($aa_return,$aa_multi_return) = &AS_event_AD(\%minus_doner, \%minus_acceptor,\%index);


open TEST,">","test_result.txt";
open MULTI,">","test_multi.txt";
my %return = %$aa_return;
my %multi_return = %$aa_multi_return;

for my $t(sort {$a cmp $b} keys %return){
	for my $p(sort {$a <=> $b} keys %{$return{$t}}){
		print TEST"$t\t$p\t$return{$t}{$p}[0]\t$return{$t}{$p}[1]\t$return{$t}{$p}[2]\n";
	}
}
close TEST;

for my $c(sort {$a cmp $b} keys %multi_return){
	for my $b(sort {$a <=> $b} keys %{$multi_return{$c}}){
		print MULTI"$c\t$b\t$multi_return{$c}{$b}[0]\t$multi_return{$c}{$b}[1]\t$multi_return{$c}{$b}[2]\t$multi_return{$c}{$b}[3]\n";
	}
}
close MULTI;

=cut





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
					$aa_return{$chr}{$doner}[0] = $chr;
					$aa_return{$chr}{$doner}[1] = $acceptor_return;
					$aa_return{$chr}{$doner}[2] = $reads_count;
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
					$aa_multi_return{$chr}{$doner}[0] = $chr;
					$aa_multi_return{$chr}{$doner}[1] = $multi_acceptor_return;
					$aa_multi_return{$chr}{$doner}[2] = "Mutil_AA";
					$aa_multi_return{$chr}{$doner}[3] = $multi_acceptor_reads;
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
					$ad_multi_return{$chr}{$acceptor}[0] = $chr;
					$ad_multi_return{$chr}{$acceptor}[1] = $multi_doner_return;
					$ad_multi_return{$chr}{$acceptor}[2] = "Mutil_DD";
					$ad_multi_return{$chr}{$acceptor}[3] = $multi_acceptor_reads;
				}
			}
		}
	}
	return (\%ad_return,\%ad_multi_return);
}

=p
my $strand = "\+";


my ($aa_return,$aa_multi_return) = &AS_event_ES(\%plus_doner, \%plus_acceptor, \%index, \%junction_postion,\%graph_postion,$strand);;


open TEST,">","test_result.txt";
open MULTI,">","test_multi.txt";
my %return = %$aa_return;
my %multi_return = %$aa_multi_return;

for my $t(sort {$a cmp $b} keys %return){
	for my $p(sort {$a <=> $b} keys %{$return{$t}}){
		print TEST"$t\t$p\t$return{$t}{$p}[0]\t$return{$t}{$p}[1]\t$return{$t}{$p}[2]\t$return{$t}{$p}[3]\n";
	}
}
close TEST;

for my $c(sort {$a cmp $b} keys %multi_return){
	for my $m(sort {$a <=> $b} keys %{$multi_return{$c}}){
		for my $n(sort {$a <=> $b} keys %{$multi_return{$c}{$m}}){
			print MULTI"$c\t$m\t$n\t$multi_return{$c}{$m}{$n}\n";
		} 
	}
}
close MULTI;

#&AS_event_ES(\%plus_doner, \%plus_acceptor, \%index, \%junction_postion,\%graph_postion,$strand);
#my $strand = "\-";
#&AS_event_ES(\%minus_doner, \%minus_acceptor, \%index, \%junction_postion, $strand);

=cut

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
								$es_return{$chr}{$ds1}[2] = $as4_doner; 
								$es_return{$chr}{$ds1}[3] = $as4_reads;
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














