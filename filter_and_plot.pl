#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;

# zhanglei

our ($opt_i, $opt_p, $opt_s, $opt_d, $opt_w, $opt_o);
getopt("s:i:p:d:w:o:");

unless ( $opt_i && $opt_o && $opt_p) {
	&usage;
	exit(1);
}

my $time1 = time();
print "Starting $0...\n";

# default parameters
$opt_w = $opt_w || 3;
$opt_s = $opt_s || 500;
$opt_d = $opt_d || 10;

# output file
my $filter_snp = $opt_o . ".filtered_SNP";

# minimum number of SNP
my $min_SNP = int($opt_w * 30 / 2 + 1);

# filter the SNP
open IN,"<$opt_i"; open OUT,">$filter_snp";
print OUT "#chr\tpos\trefBase\tmut_freq\twild_freq\tdeltaSNP\tmut_cov\tmut_A\tmut_C\tmut_G\tmut_T\twild_cov\twild_A\twild_C\twild_G\twild_T\n";
print "Step1: filtering the SNP in two pools.\n";
my $total_number = my $true_number = 0;
my %hash = ();
while (<IN>) {
	chomp;
	next if /#/;
	$total_number++;
	my ($chr, $pos, $refBase, $mut_cov, %mut, %wild, $wild_cov) = ();
	($chr, $pos, $refBase, $mut_cov, $mut{'A'}, $mut{'C'}, $mut{'G'}, $mut{'T'}, $wild_cov, $wild{'A'}, $wild{'C'}, $wild{'G'}, $wild{'T'}) = (split /\s+/);
	next unless $chr =~ /^$opt_p/;
	$chr = $1 if $chr =~ /$opt_p(\d+)/;
	# minimum read depth
	next if (($mut_cov < $opt_d) || ($wild_cov < $opt_d));	
	my @mut_alleles = sort { $mut{$b} <=> $mut{$a} } keys %mut;
	my @wild_alleles = sort { $wild{$b} <=> $wild{$a} } keys %wild;
	my $mut_freq = sprintf "%.4f", $mut{ $mut_alleles[0] } / $mut_cov;
	my $wild_freq  = sprintf "%.4f", $wild{ $mut_alleles[0] } / $wild_cov;
	my $deltaSNP = sprintf "%.4f", $mut_freq - $wild_freq;
	# SNPs with SNP-index < 0.3 in both pools are filtered out because they may be spurious SNPs
	next if ($mut_freq < 0.3 && $wild_freq < 0.3);
	# SNPs with SNP-index > 0.7 in both pools are filtered out because they may be spurious SNPs
	next if ($mut_freq > 0.7 && $wild_freq > 0.7);
		
	$true_number++;
	$hash{$chr}{$pos} = $deltaSNP;
	print OUT "$chr\t$pos\t$refBase\t$mut_freq\t$wild_freq\t$deltaSNP\t$mut_cov\t$mut{'A'}\t$mut{'C'}\t$mut{'G'}\t$mut{'T'}\t$wild_cov\t$wild{'A'}\t$wild{'C'}\t$wild{'G'}\t$wild{'T'}\n";
}
close OUT;
close IN;
print "\tinput: $total_number SNPs, after filtered: $true_number SNPs.\n";

# using a sliding-window method to split genome
my @chr_array = keys %hash;
print "Step2: using a sliding-window method to split genome, " . "chromosome number: " . @chr_array . ", sliding window: $opt_w mb, step size: $opt_s kb.\n";
my %window_hash = ();
$opt_w = $opt_w * 1000000;
$opt_s = $opt_s * 1000;
foreach my $chr (sort @chr_array) {
	my $chr_len = (sort {$a <=> $b} keys %{$hash{$chr}})[-1];
#	print "\tprocessing $opt_p$chr, length: $chr_len, ";
	# get the window number, window start, window end
	my $window_num = 0;
	if ($chr_len <= $opt_w) {               #windowsize greater than chro_length
		$window_num = 1;
		($window_hash{$chr}{$window_num}{start}, $window_hash{$chr}{$window_num}{end}) = (1, $chr_len);
	}else{                                  #windowsize less than chro_length
		$window_num = &split_genome($chr, $chr_len);
	}
#	print "window_number: $window_num.\n";
}

# caculate the average SNPindex of each window located in each chromosome
print "Step3: caculate the average SNPindex of each window located in each chromosome. The minimum numbers of SNP used to caculate SNPindex in a window is $min_SNP\n";
open ALL, ">$opt_o.ALL.window_SNPindex" or die $!;
print ALL "chr\tstart\tend\taverage_SNP_index\n";
foreach my $chr (sort @chr_array) {
	open OUT,">$opt_o.$opt_p$chr.window_SNPindex";
	print OUT "chr\twindow_number\tstart\tend\taverage_SNP_index\tSNP_number\n";
	my @array = sort {$a <=> $b} keys %{$hash{$chr}};
	my ($new_window_number, $start_t, $end_t, $SNP_number_t, $average_SNPindex_t) = (0, 0, 0, 0, 0);
	foreach my $window_number (sort {$a <=> $b} keys %{$window_hash{$chr}}) {
		my $average_SNPindex = 0;
		my $SNP_number = 0;
		my ($start, $end) = ($window_hash{$chr}{$window_number}{start}, $window_hash{$chr}{$window_number}{end});
		foreach (@array) {
			next if $_ < $start;
			last if $_ >= $end;
			$SNP_number++;
			$average_SNPindex += $hash{$chr}{$_};
		}
		# if the number of SNPs within the $opt_w Mb window was < $min_SNP, then merge the interval
		if ($SNP_number < $min_SNP && $SNP_number_t < $min_SNP) {
			$start_t = $start if $start_t == 0;
			$average_SNPindex_t += $average_SNPindex;
			$SNP_number_t += $SNP_number;
		}
		else {
#			print "SNP_number: $SNP_number;  SNP_number_t: $SNP_number_t\n";
			$start_t = $start if $start_t == 0;
			$end_t = $end;
			$SNP_number_t += $SNP_number;
			$average_SNPindex_t += $average_SNPindex;
			$new_window_number++;
			$average_SNPindex_t = sprintf "%.4f", $average_SNPindex_t/$SNP_number_t;
			print OUT "$opt_p$chr\t$new_window_number\t$start_t\t$end_t\t$average_SNPindex_t\t$SNP_number_t\n";
			$average_SNPindex_t = abs($average_SNPindex_t);
			print ALL "$chr\t$start_t\t$end_t\t$average_SNPindex_t\n";
			($start_t, $end_t, $SNP_number_t, $average_SNPindex_t) = (0, 0, 0, 0);
		}
	}
	close OUT;
}
close ALL;

# plot the average SNP-index along the coordinate
print "Step4: plot the average SNP-index.\n\tprocessing ";
foreach my $chr (sort @chr_array) {
	my $file = "$opt_o.$opt_p$chr.window_SNPindex";
	print "$opt_p$chr...";
	&plot_R($file, $chr);
}
print "done\n";

# plot all the chr together
print "Step5: plot the combined results.\n";
&plot_combined("$opt_o.ALL.window_SNPindex");


my $time2 = time();
my $time  = sprintf "%.4f", ( $time2 - $time1 ) / 60;
print "$0 finished! $time mins elapsed.\n";

sub split_genome {
	my ($chro_id, $chro_length) = @_;
	my $window_number = 0;
	if((($chro_length - $opt_w)/$opt_s) =~ /^\d+$/){
		$window_number = int(($chro_length - $opt_w)/$opt_s) + 1;
	}else{
		$window_number = int(($chro_length - $opt_w)/$opt_s) + 2;
	}
	my $window_start = $window_hash{$chro_id}{'1'}{start} = 1;                    #the first window
	my $window_end   = $window_hash{$chro_id}{'1'}{end}   = $opt_w;
	for (my $i = 2; $i < $window_number; $i++) {                                  #the middle window
		$window_start += $opt_s;
		$window_end += $opt_s;
		$window_hash{$chro_id}{$i}{start} = $window_start;
		$window_hash{$chro_id}{$i}{end} = $window_end;
	}
	$window_start += $opt_s;      $window_end += $opt_s;                          #the last window
	$window_end = $chro_length if $window_end > $chro_length;
	my $last_window_size = $window_end - $window_hash{$chro_id}{$window_number - 1}{start} + 1;   #check the last window size, if less than 500kb, then merge it to the previous window
	if ($last_window_size <= 500000) {
		$window_number = $window_number - 1;
		$window_hash{$chro_id}{$window_number}{end} = $window_end;
	}else{
		($window_hash{$chro_id}{$window_number}{start}, $window_hash{$chro_id}{$window_number}{end}) = ($window_start, $window_end);
	}
	return $window_number;
}

sub plot_combined {
	my $file = shift @_;
	my $R_script = <<END;
	data <- read.table("$file",header=TRUE)
	pdf(w=16, h=7,"$opt_o\_combined_chr.pdf")
	plot(1:nrow(data), data\$average_SNP_index, pch=16, cex=.8, col=2-data[,1]%%2, ylim=c(0, 1), ylab="delta (SNPindex)", xlab="Chromosome", xaxt='n', xaxs='i')
	breaks <- NULL
	for (row in 1:(nrow(data)-1)){
		if (data\$chr[row] != data\$chr[row+1]) {
			breaks <- append(breaks, row+1)
		}
	}
	breaks <- append(breaks, nrow(data))
	abline(v=(breaks[1:length(breaks)-1]+2), col="grey")
	
	lastbreak <- 1
	labelpos <- NULL
	for( b in breaks){
		labelpos <- append(labelpos, ((b+1)-lastbreak)/2+lastbreak)
		lastbreak <- b
	}
	mtext(unique(gsub("\\\\D", "", data\$chr[!is.na(data\$chr)])), at = labelpos, side=1, cex=.8)
	dev.off()
END
    open OUT, ">plot_combine.R" or die $!;
    print OUT "$R_script";
    close OUT;
    `R CMD BATCH ./plot_combine.R`;
   # `rm plot.R`;
   # `rm plot.Rout`; 
   # `rm $file`;
}


sub plot_R {
	my ($file, $chr) = @_;
	my $R_script = <<END;
	data <- read.table("$file",header=TRUE)
	thr <- 0.4
	start <- data\$start
	end <- data\$end
	pos <- (end - start) / 2 + start
	index <- data\$average_SNP_index
	smaller <- which(index < thr)
	larger <- which(index >= thr)
	
	pdf(w=16, h=9,"$chr.pdf")
	y <- seq(-10, 10, 2)
	y_y <- y/10
	x <- round(pos/1000000, 0)
	x_max <- max(x)
	x <- seq(0, x_max, 10)
	plot(pos/1000000, index, ylab="delta (SNPindex)", xlab="$opt_p$chr (MB)", xlim=c(0, x_max), ylim=c(-1, 1), type="o", lwd=3, pch=16 , cex=0.8, xaxt="n", yaxt="n")
	axis(side=2, y_y, tcl=-0.2)
	axis(side=1, x, tcl=-0.2)
	points(pos[smaller]/1000000, index[smaller], col="black", pch=16, cex=0.8)
	points(pos[larger]/1000000, index[larger], col="red", pch=16, cex=0.8)
#	abline(h=c(0.4, -0.4), lty = 5, col = "red")
	abline(h=c(0.0), lty = 5, col = "black")
	dev.off()
END
    open OUT, ">plot.R" or die $!;
    print OUT "$R_script";
    close OUT;
    `R CMD BATCH ./plot.R`;
  #  `rm plot.R`;
  #  `rm plot.Rout`;
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
