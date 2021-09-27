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
		my $index =1;
		my @pos;
		if($s eq "+"){
			@pos = sort {$a <=> $b} keys %{$graph_postion{$c}{$s}};
		}
		else{
			@pos = sort {$b <=> $a} keys %{$graph_postion{$c}{$s}};
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

open OUT,">",$out;

#my %as_doner = &AS_acceptor(\%plus_doner, \%plus_acceptor);
#my %as_doner = &AS_acceptor(\%minus_doner, \%minus_acceptor);

my %as_doner = &AS_doner(\%plus_doner, \%plus_acceptor);
#my %as_doner = &AS_doner(\%minus_doner, \%minus_acceptor);

for my $i(sort {$a cmp $b} keys %as_doner){
	for my $m(sort {$a <=> $b} keys %{$as_doner{$i}}){
		print OUT"$i\t$m\t$as_doner{$i}{$m}[0]\t\[$as_doner{$i}{$m}[1]\]\t\[$as_doner{$i}{$m}[2]\]\n";	
	}
}


close IN;
close OUT;





sub AS_acceptor{
	my ($doner,$acceptor) = @_;
	my %doner = %$doner;
	my %acceptor = %$acceptor;
	my %as_acceptor;
	for my $i(keys %doner){																		#######chromosome
		for my $m (keys %{$doner{$i}}){   														#######doner site	
			my @acceptor = sort {$a <=> $b} keys %{$doner{$i}{$m}};
			my $count_acceptor = @acceptor;
			if($count_acceptor == 2){
				#print "$i\t$m\t@acceptor\n";
				my @doner0 = keys %{$acceptor{$i}{$acceptor[0]}};
				my @doner1 = keys %{$acceptor{$i}{$acceptor[1]}};
				my $doner0 = keys %{$acceptor{$i}{$acceptor[0]}};
				my $doner1 = keys %{$acceptor{$i}{$acceptor[1]}};
				#print "$i\t$m\t$doner0\t$doner1\n";
				#print "$i\t$m\t@acceptor\t$acceptor{$i}{$acceptor[0]}{$m}\t$acceptor{$i}{$acceptor[1]}{$m}\n";
				if(($doner0 ==1) && ($doner1 ==1)){
					#my $name = $i."_".$m;
					my $acceptor = join"\,",@acceptor;
					my $reads_count = $acceptor{$i}{$acceptor[0]}{$m}."\,".$acceptor{$i}{$acceptor[1]}{$m};
					$as_doner{$i}{$m}[0] = $m;
					$as_doner{$i}{$m}[1] = $acceptor;
					$as_doner{$i}{$m}[2] = $reads_count;
					#print "$i\t$m\t@acceptor\n";	
				}
			}
			elsif($count_acceptor > 2){
				
			}
		}
	}
	return %as_acceptor;
}





sub AS_doner{
	my ($doner,$acceptor) = @_;
	my %doner = %$doner;
	my %acceptor = %$acceptor;
	my %as_doner;
	for my $i(keys %acceptor){																		#######chromosome
		for my $m (keys %{$acceptor{$i}}){   														#######doner site	
			my @doner = sort {$a <=> $b} keys %{$acceptor{$i}{$m}};
			my $count_doner = @doner;
			if($count_doner == 2){
				#print "$i\t$m\t@acceptor\n";
				my @acceptor0 = keys %{$doner{$i}{$doner[0]}};
				my @acceptor1 = keys %{$doner{$i}{$doner[1]}};
				my $acceptor0 = keys %{$doner{$i}{$doner[0]}};
				my $acceptor1 = keys %{$doner{$i}{$doner[1]}};
				#print "$i\t$m\t$doner0\t$doner1\n";
				#print "$i\t$m\t@acceptor\t$acceptor{$i}{$acceptor[0]}{$m}\t$acceptor{$i}{$acceptor[1]}{$m}\n";
				if(($acceptor0 ==1) && ($acceptor1 ==1)){
					#my $name = $i."_".$m;
					my $doner = join"\,",@doner;
					my $reads_count = $doner{$i}{$doner[0]}{$m}."\,".$doner{$i}{$doner[1]}{$m};
					$as_doner{$i}{$m}[0] = $doner;
					$as_doner{$i}{$m}[1] = $m;
					$as_doner{$i}{$m}[2] = $reads_count;
					#print "$i\t$m\t@acceptor\n";	
				}
			}
		}
	}
	return %as_doner;
}


#my $strand = "\+";
#&exon_skipping(\%plus_doner, \%plus_acceptor, \%index, $strand);
my $strand = "\-";
&exon_skipping(\%minus_doner, \%minus_acceptor, \%index, $strand);


sub exon_skipping{																			#######only for one exon skipping
	my ($doner,$acceptor,$index,$strand) = @_;
	my %doner = %$doner;
	my %acceptor = %$acceptor;
	my %index = %$index;
	my %as_acceptor;
	for my $i(keys %doner){																		#######chromosome
		for my $m (keys %{$doner{$i}}){   														#######doner site	
			my @acceptor = sort {$a <=> $b} keys %{$doner{$i}{$m}};
			my $count_acceptor = @acceptor;
			if($count_acceptor == 2){
				my $as4;
				my $as2;
				my $ds1 = $m;
				my $ds3;
				my $delta1 = abs($acceptor[0] - $m);
				my $delta2 = abs($acceptor[1] - $m);
				#print "$delta1\t$delta2\n";
				if($delta1 > $delta2){
					$as2 = $acceptor[1];
					$as4 = $acceptor[0];
				}
				else{
					$as2 = $acceptor[0];
					$as4 = $acceptor[1];
				}
				####delta index
				#print "$m\t";
			#	print "$ds1\t$as2\t$as4\n";
				my $delta_index = abs($index{$i}{$strand}{$as4} - $index{$i}{$strand}{$ds1});
				#print "$delta_index\n";
				my @as4_doner = keys %{$acceptor{$i}{$as4}};
				my $as4_doner = keys %{$acceptor{$i}{$as4}};
				if($as4_doner ==2){
					for my $as_d(@as4_doner){
						if($as_d != $ds1){
							$ds3 = $as_d;
							my $as2_doner = keys %{$acceptor{$i}{$as2}};
							my $as4_doner = keys %{$acceptor{$i}{$as4}};
							my $ds3_acceptor = keys %{$doner{$i}{$ds3}};
							if($delta_index ==3){
								print "$count_acceptor\t$as2_doner\t$ds3_acceptor\t$as4_doner\n";
								print "$i\t$strand\t$ds1\t$as2\t$ds3\t$as4\t@acceptor\t@as4_doner\n";	
							}
						}
					}
				}	
			}
		}
	}
	#return %as_acceptor;
}





=p
open OUT,">",$out;
open NONAS,">","non_as.txt";
open MIL,">","multi.txt";
for my $i(keys %plus_doner){
	for my $m (keys %{$plus_doner{$i}}){
		my @acceptor = sort {$a <=> $b} keys %{$plus_doner{$i}{$m}};
		my $count_acceptor = @acceptor;
		if($count_acceptor == 2){
			my @doner = keys %{$plus_acceptor{$i}{$acceptor[1]}};
			my $count_doner = keys %{$plus_acceptor{$i}{$acceptor[1]}};
			if($count_doner ==2){
				print OUT"$i\t$m\t$acceptor[0]\,$acceptor[1]\t$acceptor[1]\t$doner[0]\,$doner[1]\n";	
			}
		}
		elsif($count_acceptor == 1){
			print NONAS"$i\t$m\t$acceptor[0]\n"
		}
		elsif($count_acceptor >2){
			print MIL"$i\t$m\t@acceptor\n";
		}
	}
}

close OUT;
close NONAS;
close MIL;












