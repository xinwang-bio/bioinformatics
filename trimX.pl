#!/usr/bin/perl -w
use strict;

my ($in,$out,$window,$mismatch) = @ARGV;


my $min_length = -10;
#my $window = 15;
#my $mismatch = 1;
my $extend = 3;

my $seq1 = "ATTTGCTTTAAAAAATTGTTTATCGACTAAAAAAAAAAAAAAAAAAAAAACGGGATATGTTTTTAAAAATTTCAAAACTTTTATTGCGCCAAAAAATTTTT";
my $seq2 = "TGCTTCACTGCATCTAATACGACTCACTATAGGGAGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTAATGATCAATGACAGCCCCAGATTTTTTTGTGCCC";
#my $time1 = time();
#for my $file(@ARGV){
	open IN,"<",$in;
	open OUT,">",$out;
	while(<IN>){
		chomp;
		my $name = $_;
		my $seq = <IN>;
		chomp($seq);
		my $third = <IN>;
		chomp($third);
		my $quality = <IN>;
		chomp($quality);
		my $length = length($seq);
		#print "$name\t$length\n$seq\n";
		my ($check_polyA_result,$check_polyT_result) = &check_polyX($seq,$window,$mismatch);
		#print"$check_polyA_result\t$check_polyT_result\n$quality\n";
		my ($trim_seq,$trim_quality) =("N","N");
		if(($check_polyA_result==1) && ($check_polyT_result==0)){
			($trim_seq,$trim_quality) = &trim_polyX($seq,$quality,$window,$mismatch,$extend);
			#print "$trim_seq\n$trim_quality\n\n";
		}
		elsif(($check_polyA_result==0) && ($check_polyT_result==1)){
			my $rc_seq = &revcom($seq);
			my $rc_quality = reverse($quality);
			($trim_seq,$trim_quality) = &trim_polyX($rc_seq,$rc_quality,$window,$mismatch,$extend);
			$trim_seq = &revcom($trim_seq);
			$trim_quality = reverse($trim_quality);
			#print "$trim_seq\n$trim_quality\n\n";
		}
		elsif(($check_polyA_result==0) && ($check_polyT_result==0)){
			$trim_seq = $seq;
			$trim_quality = $quality;
		}
		my $trim_length = length($trim_seq);
		if($trim_length>=$min_length){
			print OUT"$name\n$trim_seq\n$third\n$trim_quality\n";
		}
	}
	close IN;
#}
#my $time2 = time();
#my $time  = sprintf "%.4f", ( $time2 - $time1 ) / 60;
#print "$0 finished! $time mins elapsed.\n";
	

sub trim_polyX{
	my ($seq,$quality,$window,$mismatch,$extend) = @_;
	my $length = length($seq);
	my $start_postion = $length-1;
	my $trim_seq ="N";
	my $trim_quality = "N";
	for my $i(0..($length -$window)){		######check the polyA/T postion
		my $start_nt = substr($seq, $i, 1);
		$start_nt =~ s/a/A/;
		my $sub_seq = substr($seq, $i, $window);
		my $extend_seq = substr($seq, $i, $extend);
		#print "$sub_seq\t";
		my $count_a = 0;
		my $count_t = 0;
		if($sub_seq =~ /A|a/){
			$count_a = ($sub_seq =~ s/A/A/ig);
		}
		if(($window - $count_a <= $mismatch) && ($start_nt eq "A") && ($extend_seq eq "AAA")){
			$start_postion = $i;
			last;
		}
	}
	#print "$start_postion\n";
	if($start_postion !=0){
		$trim_seq = substr($seq,0,$start_postion);
		$trim_quality = substr($quality,0,$start_postion);
	}
	return ($trim_seq,$trim_quality); 
}

sub check_polyX{
	my ($seq,$window,$mismatch) = @_; 
	my $check_polyA = 0;
	my $check_polyT = 0;
	my $length =length($seq);
	for my $i(0..($length -$window)){		######check the polyA/T

		my $start_nt = substr($seq, $i, 1);
		my $sub_seq = substr($seq, $i, $window);
		#print "$sub_seq\t";
		#print "$start_nt\t";
		my $count_a = 0;
		my $count_t = 0;
		if($sub_seq =~ /A|a/){
			$count_a = ($sub_seq =~ s/A/A/ig);
		}
		if($sub_seq =~ /T|t/){
			$count_t = ($sub_seq =~ s/T/T/ig);
		}
		#print "$count_a\t$count_t\t";
		if($window - $count_a <= $mismatch){
			$check_polyA = 1;
		}
		if($window - $count_t <= $mismatch){
			$check_polyT = 1;
		}
	}
	return($check_polyA,$check_polyT);
	#print "\n";
}
sub revcom {
	my $seq = shift;
	my $rev = reverse $seq;
	(my $revcom = $rev) =~ tr/ACGTacgt/TGCAtgca/;
	return $revcom;
}



