#!/usr/bin/env perl 

# Author: Xin Wang

use warnings;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Getopt::Std;

&main;
exit;

sub main {
  &usage if (@ARGV < 1);
  my $command = shift(@ARGV);
  my %func = (SNPfromPool=>\&SNPfromPool, VCFfilter=>\&VCFfilter, Cutoff=>\&cutoff, filter_and_plot=>\&filter_and_plot, plotCutoff=>\&plotCutoff,
  			   hapmap2vcf=>\&hapmap2vcf,
			  );
  die("Unknown command \"$command\".\n") if (!defined($func{$command}));
  &{$func{$command}};
}


sub SNPfromPool{
		
	our %opt;
	
	getopt('r:w:m:q:Q:p:d:o:', \%opt);
	
	unless ( $opt{r} && $opt{w} && $opt{m} && $opt{o} ) {
	print"Usage: $0 -OPTIONS VALUES\n
	This script is used to call SNP between wild pool and mutant pool.
	Options:
	     -r  YES  reference fasta file
	     -w  YES  bamfile of wild(dominant trait) pool
	     -m  YES  bamfile of mutant(recessive trait) pool
	     -o  YES  output snp file name
	     -d  NO   minimum read depth required for calling SNP, default = 10
	     -q  NO   minimum mapping quality required for calling SNP, default = 30
	     -Q  NO   minimum base quality required for calling SNP, default = 20
	     \n";
		exit(1);
	}
	
	# default parameters
	$opt{q} = $opt{q} || 30;
	$opt{q} = $opt{q} || 20;
	$opt{d} = $opt{d} || 10;
	
	open OUT,">$opt{o}";
	print OUT "#chr\tpos\trefbase\tmut_cov\tmut_A\tmut_C\tmut_G\tmut_T\twt_cov\twt_A\twt_C\twt_G\twt_T\n";
	
	## main program
	print "Running samtools mpileup command below:\n";
	print "samtools mpileup -q $opt{q} -Q $opt{q} -f $opt{r} $opt{m} $opt{w}\n";
	print "Parsing the output in real time...\n";
	open IN,"samtools mpileup -q $opt{q} -Q $opt{q} -f $opt{r} $opt{m} $opt{w} |" or die $!;
	while (my $line = <IN>) {
		chomp $line;
		next if $line =~ /^\[/;
		next if $line =~ /^</;
		my @Contents = (split /\s+/,$line);
		
		# Omit any pos where only one sample mapped
		next if @Contents < 9;
		my ( $Chr, $Pos, $refBase ) = @Contents[ 0, 1, 2 ];
		
		# For convenience, Sample1 for mutant pool, Sample2 for wild pool
		my ( $mutcov0, $mutbases, $wtcov0, $wtbases ) = @Contents[ 3, 4, 6, 7 ];
		
		# Omit low-coverage position
		next if ( $mutcov0 < $opt{d} || $wtcov0 < $opt{d});
		
		# Calculate reads coverage for the position
		my @mutCounts = &base_counter( $mutbases, $refBase );
		my @wtCounts  = &base_counter( $wtbases,  $refBase );
		my %mut       = ('A' => $mutCounts[0],'C' => $mutCounts[1],'G' => $mutCounts[2],'T' => $mutCounts[3]);
		my %wt        = ('A' => $wtCounts[0],'C' => $wtCounts[1],'G' => $wtCounts[2],'T' => $wtCounts[3]);
		
		my $mutcov1 = $mut{'A'} + $mut{'C'} + $mut{'G'} + $mut{'T'};
		my $wtcov1  = $wt{'A'} + $wt{'C'} + $wt{'G'} + $wt{'T'};
		
		# the true read depth
		next if ( $mutcov1 < $opt{d} || $wtcov1 < $opt{d} );
		
		print OUT  "$Chr\t$Pos\t$refBase\t$mutcov1\t".join("\t",@mutCounts)."\t$wtcov1\t".join("\t",@wtCounts)."\n";
	}
	close OUT;
	print "\tdone\n";
	
	
	## subroutine
	# base_counter: calculate base counts for each base in (A,C,G,T) order
	sub base_counter {
		my ( $sample_bases, $refbase ) = @_;
		
		# Convert all dot and comma symbol to ref base
		$sample_bases =~ s/\.|,/$refbase/gi;
		
		# Remove patterns that represents INDELs
		while ( $sample_bases =~ /(.*?)[+-](\d+)[ATCG.,]+/ig ) {
			$sample_bases =~ s/(.*?)[+-](\d+)[ATCGNatcgn]{$2}(.*)/$1$3/i;
		}
	
		# count Aa/Cc/Gg/Tt
		my $baseA = ($sample_bases =~ tr/Aa//);
		my $baseC = ($sample_bases =~ tr/Cc//);
		my $baseG = ($sample_bases =~ tr/Gg//);
		my $baseT = ($sample_bases =~ tr/Tt//);
		
		return ( $baseA, $baseC, $baseG, $baseT );
	}
	
=p	
	sub howto {
	print"Usage: $0 -OPTIONS VALUES\n
	This script is used to call SNP between wild pool and mutant pool.
	Options:
	     -r  YES  reference fasta file
	     -w  YES  bamfile of wild(dominant trait) pool
	     -m  YES  bamfile of mutant(recessive trait) pool
	     -o  YES  output snp file name
	     -d  NO   minimum read depth required for calling SNP, default = 10
	     -q  NO   minimum mapping quality required for calling SNP, default = 30
	     -Q  NO   minimum base quality required for calling SNP, default = 20
	     \n"
	}
=cut
}


sub VCFfilter{
	my $file = $ARGV[0];
	my $out  = $ARGV[1];
die "Usage: $0 file.vcf outfile\n\nThis script is used to filter the VCF files
	Note: Before run this script, VCF file must be called and filterd by samtools, the commands are as follows:
	Step1: samtools mpileup <-r region[chr:start-end]> -C 50 -q 30 -Q 20 -u -f ref.fa file1.bam file2.bam | bcftools view -Nbvcg -> file.raw.bcf
	Step2: bcftools view file.raw.bcf | vcfutils.pl varFilter -Q 20 -d 10 -D 100 > var.flt.vcf\n" if @ARGV < 2;
	
	### main program ###
	open IN, "<$file" or die $!;
	open OUT, ">$out" or die $!;
	my ($Sample1, $Sample2) = ();
	while (<IN>) {
		chomp;
		next if /^##/;
		($Sample1, $Sample2) = (split)[9, 10] and print OUT "#chr\tpos\ttype\treference\t$Sample1\t$Sample2\n" and next if /^#CHROM\s/;
		my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $genotype1, $genotype2) = (split);
		
		my ($sense_variant, $anti_variant) = $INFO =~ /;DP4=\d+,\d+,(\d+),(\d+);/;
		my $total_variant = $sense_variant + $anti_variant;
		next if $total_variant < 7;
		
		my $PL1 = (split /:/, $genotype1)[1];
		my $PL2 = (split /:/, $genotype2)[1];
		my @PL1s = (split /\,/, $PL1);
		my @PL2s = (split /\,/, $PL2);
		my $value1 = checkPL(@PL1s);
		next if $value1;
		my $value2 = checkPL(@PL2s);
		next if $value2;
		
		my @genotypes;
		my @alts = (split /\,/, $ALT);
		push @genotypes, $REF;
		push @genotypes, @alts;
		my $geno_number = scalar @genotypes;
		my $index1 = assign_genotype($geno_number, @PL1s);
		my $index2 = assign_genotype($geno_number, @PL2s);
		next if $index1 eq "FALSE" || $index2 eq "FALSE";
		next if $index1 eq $index2;
		$genotype1 = $genotypes[$index1];
		$genotype2 = $genotypes[$index2];
		
		next if $REF =~ /N/i;
		next if $genotype1 =~ /N/i;
		next if $genotype2 =~ /N/i;
		
		next if $REF =~ /A{4,}|T{4,}|C{4,}|G{4,}/i;
		next if $genotype1 =~ /A{4,}|T{4,}|C{4,}|G{4,}/i;
		next if $genotype2 =~ /A{4,}|T{4,}|C{4,}|G{4,}/i;
		
		my $ref_len = length($REF);
		my $geno1_len = length($genotype1);
		my $geno2_len = length($genotype2);
		next if (($ref_len > 7) || ($geno1_len > 7) || ($geno2_len > 7));
		
		unless ($INFO =~ /^INDEL/) {
			print OUT "$CHROM\t$POS\tSNP\t$REF\t$genotype1\t$genotype2\n";
		}
		else {
			if (($ref_len == $geno1_len) || ($ref_len == $geno2_len)) {
				my $diff_len = $geno1_len == $ref_len ? $geno2_len : $geno1_len;
				if ($ref_len - $diff_len < 0) {                                   # INSERTION
					print OUT "$CHROM\t$POS\tINS\t$REF\t$genotype1\t$genotype2\n";
				}
				elsif ($ref_len - $diff_len > 0) {                                # DELETION
					print OUT "$CHROM\t$POS\tDEL\t$REF\t$genotype1\t$genotype2\n";
				}
			}
			else {                                                                # INSERTION
				print OUT "$CHROM\t$POS\tINS\t$REF\t$genotype1\t$genotype2\n";
			}
		}
	}
	close IN;
	close OUT;
	
	### subroutines ###
	sub checkPL {
		my @array = @_;
		my $number = 0;
		foreach my $value (@array) {
			$number++ if $value == 0;
		}
		if ($number == 1) {
			return "";
		}
		else {
			return "1";
		}
	}
	
	sub assign_genotype {
		my ($number, @PL1s) = @_;
		$number = $number - 1;
		foreach my $geno_index (0..$number) {
			my $PL_index = ($geno_index * ($geno_index + 1) / 2) + $geno_index;
			return "$geno_index" if $PL1s[$PL_index] == 0;
		}
		return "FALSE";
	}
}


sub filter_and_plot{
		
	our ($opt_i, $opt_p, $opt_s, $opt_d, $opt_w, $opt_o);
	getopt("s:i:p:d:w:o:");
	
	unless ( $opt_i && $opt_o && $opt_p) {
		print"Usage: $0 -OPTIONS VALUES
	This script is used to caculate the SNPindex and plot it on every chromosome.
	Options:
	     -i  YES  snp file produced by callSNPfromBam script
	     -o  YES  prefix of output file name
	     -p  YES  prefix of the chromosome name
	     -d  NO   minimum read depth required for calling SNP, default = 10
	     -w  NO   windowsize to caculate the average SNPindex(megabase unit), default = 3Mb
	     -s  NO   step size, should be greater than 0kb and less than windowsize(kilobase unit), default = 500kb 
	     \n";
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
	my $min_SNP = int($opt_w * 30 / 2);
	
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
	our %window_hash = ();	#####change the my to our to avoid the Variable "%window_hash" will not stay shared 
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
		&plot_R($file, "$opt_o.$opt_p$chr");
	}
	print "done\n";
	
	# plot all the chr together
	print "Step5: plot the combined results.\n";
	&plot_combined("$opt_o.ALL.window_SNPindex");
	
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
	
		
	
	my $time2 = time();
	my $time  = sprintf "%.4f", ( $time2 - $time1 ) / 60;
	print "$0 finished! $time mins elapsed.\n";
=p
	sub howto {
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
=cut		
}

sub cutoff{
	our ($opt_i, $opt_o, $opt_p, $opt_m, $opt_w, $opt_r, $opt_s, $opt_f, $opt_d);
	getopt("i:o:p:m:w:r:p:f:d");
	
	unless ( $opt_i && $opt_o && $opt_m && $opt_w) {
		print"Usage: $0 -OPTIONS VALUES
	This script is used to add the simulate value to the SNPindex result.
	Options:
	     -i  YES  inputfile (filtered snp index result) 
	     -o  YES  output file name
	     -p  YES  prefix of the chromosome name(such as SL2.40ch), defaule = SL2.40ch
	     -m  YES  mutation number
	     -w  YES  wild number
	     -r  YES  replicate number, defaule = 1000
	     -s  YES  population structure(F2 or RIL), defaule = F2
	     -f  YES  SNP index filter value, defaule = 0.3
	     -d  YES  max depth for each site, defaule = 300
	     \n";
		exit(1);
	}
=p
	his script is used to add the simulate value to the SNPindex result.
	Options:
	     -i  YES  inputfile (filtered snp index result) 
	     -o  YES  output file name
	     -p  NO  prefix of the chromosome name(such as SL2.40ch), defaule = SL2.40ch
	     -m  YES  mutation number
	     -w  YES  wild number
	     -r  NO  replicate number, defaule = 1000
	     -p  NO  population structure(F2 or RIL), defaule = F2
	     -f  NO  SNP index filter value, defaule = 0.3
	     -d  NO  max depth for each site, defaule = 300
	
=cut
	# default parameters
	$opt_p = $opt_p || "SL2.40ch";
	$opt_r = $opt_r || 1000;
	$opt_s = $opt_s || "F2";
	$opt_f = $opt_f || 0.3;
	$opt_d = $opt_d || 300;
	
	my $input=$opt_i;
	my $out =$opt_o;
	my $mut_num=$opt_m;
	my $wild_num=$opt_w;
	my $replicate=$opt_r;
	my $population_structure=$opt_s;
	my $filter_value=$opt_f;
	
	
	my $time1=time();
	system "mkdir tmp";
	system "mkdir tmp/a";
	system "mkdir tmp/b";
	
	#open TE1,">","test1.txt";
	#open TE2,">","test2.txt";
	#open TE,">","test.txt";
	
	
	my $a_population="a";
	my $b_population="b";
	&simulate_test($a_population);
	system"Rscript tmp/$a_population/simulate_test.r $mut_num $replicate $population_structure $opt_d";
	&simulate_test($b_population);
	system"Rscript tmp/$b_population/simulate_test.r $wild_num $replicate $population_structure $opt_d";
	#print "$replicate\t$population_structure\n";
	my %ha;
	for my $m(1..$opt_d){   #### 300 is depth
		my $file="$mut_num"."_"."$m"."individuals.txt";
		#print"$file\n";
		open IN,"<","tmp/a/$file";
		while(<IN>){
			chomp;
			if(/DEPTH/){
				
			}
			else{
				my @line=split/\t/,$_;
				#print"$line[1]\n";
				push (@{$ha{$a_population}{$m}},$line[1]);
			}
		}
		close IN;
	}
	for my $n(1..$opt_d){   #### 300 is depth
		my $file="$wild_num"."_"."$n"."individuals.txt";
		#print"$file\n";
		open IN,"<","tmp/b/$file";
		while(<IN>){
			chomp;
			if(/DEPTH/){
				
			}
			else{
				my @line=split/\t/,$_;
				#print"$line[1]\n";
				push (@{$ha{$b_population}{$n}},$line[1]);
			}
		}	
		close IN;
	}
	
	
	
	open IN,"<",$input;
	open OUT,">",$out;
	my $head;
	while(<IN>){
		chomp;
		if(/#/){
			$head=$_;
			print OUT "$head\tpl_99\tpl_95\tph_95\tph_99\n";
		}
		else{
			my @line=split/\t/,$_;
			my $mut_depth=$line[6];
			my $wild_depth=$line[11];
			next if(($mut_depth > $opt_d)||($wild_depth > $opt_d));
			my @delta_snp_index;
			for my $i(0..($replicate-1)){
				if($ha{$a_population}{$mut_depth}[$i] >= $filter_value || $ha{$b_population}{$wild_depth}[$i] >= $filter_value){
					my $snp_index_delta=$ha{$a_population}{$mut_depth}[$i]-$ha{$b_population}{$wild_depth}[$i];
					push(@delta_snp_index,$snp_index_delta);
				}
			}
			my $delta_snp_length=@delta_snp_index;
			@delta_snp_index=sort{$a <=> $b} @delta_snp_index;
			#print TE"@delta_snp_index\n";
			my $array_index;
			
			####do not -1#### $array_index =(int($delta_snp_length*0.005) == ($delta_snp_length*0.005) ? (int($delta_snp_length*0.005)) : (int($delta_snp_length*0.005) - 1));
			$array_index =(int($delta_snp_length*0.005));
			#print"$array_index\t";
			my $left_1=$delta_snp_index[($array_index-1)];
			
		####do not -1#### 	$array_index =(int($delta_snp_length*0.025) == ($delta_snp_length*0.025) ? (int($delta_snp_length*0.025)) : (int($delta_snp_length*0.025) - 1));
			$array_index =(int($delta_snp_length*0.025)); 
			my $left_5=$delta_snp_index[($array_index-1)];
			#print"$array_index\t";
			$array_index =(int($delta_snp_length*0.975) == ($delta_snp_length*0.975) ? (int($delta_snp_length*0.975)) : (int($delta_snp_length*0.975) + 1));
			my $right_1=$delta_snp_index[($array_index-1)];
			#print"$array_index\t";
			$array_index =(int($delta_snp_length*0.995) == ($delta_snp_length*0.995) ? (int($delta_snp_length*0.995)) : (int($delta_snp_length*0.995) + 1));
			my $right_5=$delta_snp_index[($array_index-1)];
			#print"$array_index\t";
			my $chrom=shift(@line);
			@line=join"\t",@line;
			print OUT"$opt_p$chrom\t@line\t$left_1\t$left_5\t$right_1\t$right_5\n";
		}
	}
	
	close IN;
	close OUT;
	#close TE;
	
	
	
	
	my $time2=time();
	
	my $delat_time=($time2-$time1);
	
	print"Total time is $delat_time s \n";
	
	sub simulate_test {
		my $document = shift @_;
		my $R_script = <<END;
		setwd("tmp/$document")
	genotype<-function(){
		count<-0
		
		
		if (population_structure=="RIL"){
			x<-runif(1) 
			if (x<=0.5){
				count<-1
			}else{
				count<-0
			}	
			
		}else{
			for(i in 1:2){
				x<-runif(1) 
				if (x<=0.5){
					number<-0.5
				}else{
					number<-0
				}	
				if(number == 0.5){
					count<- count+0.5
				}
			}
		}
		return(count)
	}
	individuals_genotype<-function(number_of_total_individuals){
		
		ratio_of_genotype<-c()
		for(i in 1:number_of_total_individuals){
			ratio_of_genotype<-c(ratio_of_genotype,genotype())
			
		}
		return(mean(ratio_of_genotype))	
	}
	
	
	snp_index<-function(read_depth,ratio_of_genotype_in_the_population_in_A){
		x1<-rbinom(1,read_depth,ratio_of_genotype_in_the_population_in_A)
		return(x1/read_depth)	
	}
	
	Arg<-commandArgs(TRUE)
	###########input###############################################
	
	individual_analysis<- c(Arg[1])
	reprication<-c(Arg[2])
	population_structure<-c(Arg[3]) #RIL or F2
	depth_analysis<-c(Arg[4])
	#depth_analysis<-c(1:300)
	
	for (key_individual in individual_analysis){
		individual_number<-key_individual 
		depth_data<-c()
		for (key_depth in 1:depth_analysis){
			depth_data<-c(depth_data,key_depth)
			depth<-key_depth
			data_of_delta_snp_index<-c()  
			Snp_index_of_A_array<-c()
			Snp_index_of_B_array<-c()
			for(i in 1:reprication){    
	        ratio_of_genotype_in_the_population_in_A<-individuals_genotype(key_individual)
	        Snp_index_of_A<-snp_index(key_depth,ratio_of_genotype_in_the_population_in_A)
	        Snp_index_of_A_array<-c(Snp_index_of_A_array,Snp_index_of_A)					
	        ratio_of_genotype_in_the_population_in_B<-individuals_genotype(key_individual)
	        Snp_index_of_B<-snp_index(key_depth,ratio_of_genotype_in_the_population_in_B)			
	        Snp_index_of_B_array<-c(Snp_index_of_B_array,Snp_index_of_B)
			}
			FINAL_DATA<-data.frame(DEPTH=depth,Snp_index_of_B_array)
	    table_name<-paste("./",individual_number,sep="")
	    table_name<-paste(table_name,"_",sep="")
	    table_name<-paste(table_name,key_depth,sep="")
	    table_name<-paste(table_name,"individuals.txt",sep="")
	    write.table(FINAL_DATA,table_name,sep="\t", quote=F, append=F,row.name=F)
	  } 
	}  
END
	    open OUT, ">tmp/$document/simulate_test.r" or die $!;
	    print OUT "$R_script";
	    close OUT;
	   # `R CMD BATCH ./simulate_test.R $mut_num $replicate $population_structure`;
	    #`rm simulate_test.R`;
	    #`rm simulate_test.Rout`; 
	}
	
	#$replicate,$population_structure,$filter_value,$out
=p
	sub usage {
	print"Usage: $0 -OPTIONS VALUES
	This script is used to add the simulate value to the SNPindex result.
	Options:
	     -i  YES  inputfile (filtered snp index result) 
	     -o  YES  output file name
	     -p  YES  prefix of the chromosome name(such as SL2.40ch), defaule = SL2.40ch
	     -m  YES  mutation number
	     -w  YES  wild number
	     -r  YES  replicate number, defaule = 1000
	     -s  YES  population structure(F2 or RIL), defaule = F2
	     -f  YES  SNP index filter value, defaule = 0.3
	     -d  YES  max depth for each site, defaule = 300
	     \n"
	}
=cut		
}

sub plotCutoff{
		
	our ($opt_i, $opt_p, $opt_s, $opt_d, $opt_w, $opt_o);
	getopt("s:i:p:d:w:o:");
	
	unless ( $opt_i && $opt_o && $opt_p) {
		print"Usage: $0 -OPTIONS VALUES
	This script is used to caculate the SNPindex and plot it on every chromosome.
	Options:
	     -i  YES  snp file produced by callSNPfromBam script
	     -o  YES  prefix of output file name
	     -p  YES  prefix of the chromosome name
	     -d  NO   minimum read depth required for calling SNP, default = 10
	     -w  NO   windowsize to caculate the average SNPindex(megabase unit), default = 3Mb
	     -s  NO   step size, should be greater than 0kb and less than windowsize(kilobase unit), default = 500kb 
	     \n";
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
	my $min_SNP = int($opt_w * 30 / 2);
	
	# filter the SNP
	open IN,"<$opt_i"; open OUT,">$filter_snp";
	print OUT "#chr\tpos\trefBase\tmut_freq\twild_freq\tdeltaSNP\tmut_cov\tmut_A\tmut_C\tmut_G\tmut_T\twild_cov\twild_A\twild_C\twild_G\twild_T\n";
	print "Step1: filtering the SNP in two pools.\n";
	my $total_number = my $true_number = 0;
	my %hash = ();
	my %pl_99 = ();
	my %pl_95 = ();
	my %ph_99 = ();
	my %ph_95 = ();
	
	while (<IN>) {
		chomp;
		next if /#/;
		$total_number++;
		my ($chr, $pos, $refBase, $mut_freq, $wild_freq, $deltaSNP, $mut_cov, %mut, %wild, $wild_cov, $pl_99,$pl_95, $ph_95, $ph_99) = ();
		($chr, $pos, $refBase, $mut_freq, $wild_freq, $deltaSNP, $mut_cov, $mut{'A'}, $mut{'C'}, $mut{'G'}, $mut{'T'}, $wild_cov, $wild{'A'}, $wild{'C'}, $wild{'G'}, $wild{'T'}, $pl_99, $pl_95, $ph_95, $ph_99) = (split /\s+/);
		next unless $chr =~ /^$opt_p/;
		$chr = $1 if $chr =~ /$opt_p(\d+)/;
		# minimum read depth
		next if (($mut_cov < $opt_d) || ($wild_cov < $opt_d));	
		my @mut_alleles = sort { $mut{$b} <=> $mut{$a} } keys %mut;
		my @wild_alleles = sort { $wild{$b} <=> $wild{$a} } keys %wild;
		#my $mut_freq = sprintf "%.4f", $mut{ $mut_alleles[0] } / $mut_cov;
		#my $wild_freq  = sprintf "%.4f", $wild{ $mut_alleles[0] } / $wild_cov;
		#my $deltaSNP = sprintf "%.4f", $mut_freq - $wild_freq;
		# SNPs with SNP-index < 0.3 in both pools are filtered out because they may be spurious SNPs
		next if ($mut_freq < 0.3 && $wild_freq < 0.3);
		# SNPs with SNP-index > 0.7 in both pools are filtered out because they may be spurious SNPs
		next if ($mut_freq > 0.7 && $wild_freq > 0.7);
			
		$true_number++;
		$hash{$chr}{$pos} = $deltaSNP;
		$pl_99{$chr}{$pos} = $pl_99;
		$pl_95{$chr}{$pos} = $pl_95;
		$ph_95{$chr}{$pos} = $ph_95;
		$ph_99{$chr}{$pos} = $ph_99;
		print OUT "$chr\t$pos\t$refBase\t$mut_freq\t$wild_freq\t$deltaSNP\t$mut_cov\t$mut{'A'}\t$mut{'C'}\t$mut{'G'}\t$mut{'T'}\t$wild_cov\t$wild{'A'}\t$wild{'C'}\t$wild{'G'}\t$wild{'T'}\n";
	}
	close OUT;
	close IN;
	print "\tinput: $total_number SNPs, after filtered: $true_number SNPs.\n";
	
	# using a sliding-window method to split genome
	my @chr_array = keys %hash;
	print "Step2: using a sliding-window method to split genome, " . "chromosome number: " . @chr_array . ", sliding window: $opt_w mb, step size: $opt_s kb.\n";
	our %window_hash = ();
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
			$window_num = &split_genome_cutoff($chr, $chr_len);
		}
	#	print "window_number: $window_num.\n";
	}
	
	# caculate the average SNPindex of each window located in each chromosome
	print "Step3: caculate the average SNPindex of each window located in each chromosome. The minimum numbers of SNP used to caculate SNPindex in a window is $min_SNP\n";
	open ALL, ">$opt_o.ALL.window_SNPindex" or die $!;
	print ALL "chr\tstart\tend\taverage_SNP_index\tpl_99\tpl_95\tph_95\tph_99\n";
	foreach my $chr (sort @chr_array) {
		open OUT,">$opt_o.$opt_p$chr.window_SNPindex";
		print OUT "chr\twindow_number\tstart\tend\taverage_SNP_index\tSNP_number\tpl_99\tpl_95\tph_95\tph_99\n";
		my @array = sort {$a <=> $b} keys %{$hash{$chr}};
		my ($new_window_number, $start_t, $end_t, $SNP_number_t, $average_SNPindex_t) = (0, 0, 0, 0, 0);
		my $pl_99_average_t=0;
		my $pl_95_average_t=0;
		my $ph_99_average_t=0;
		my $ph_95_average_t=0;
		
		foreach my $window_number (sort {$a <=> $b} keys %{$window_hash{$chr}}) {
			my $average_SNPindex = 0;
			my $pl_99_average=0;
			my $pl_95_average=0;
			my $ph_99_average=0;
			my $ph_95_average=0;
			
			my $SNP_number = 0;
			my ($start, $end) = ($window_hash{$chr}{$window_number}{start}, $window_hash{$chr}{$window_number}{end});
			foreach (@array) {
				next if $_ < $start;
				last if $_ >= $end;
				$SNP_number++;
				$average_SNPindex += $hash{$chr}{$_};
				$pl_99_average +=$pl_99{$chr}{$_};
				$pl_95_average +=$pl_95{$chr}{$_};
				$ph_99_average +=$ph_99{$chr}{$_};
				$ph_95_average +=$ph_95{$chr}{$_};	
			}
			# if the number of SNPs within the $opt_w Mb window was < $min_SNP, then merge the interval
			if ($SNP_number < $min_SNP && $SNP_number_t < $min_SNP) {
				$start_t = $start if $start_t == 0;
				$average_SNPindex_t += $average_SNPindex;
				$pl_99_average_t +=$pl_99_average;
				$pl_95_average_t +=$pl_95_average;
				$ph_99_average_t +=$ph_99_average;
				$ph_95_average_t +=$ph_95_average;
				$SNP_number_t += $SNP_number;
			}
			else {
	#			print "SNP_number: $SNP_number;  SNP_number_t: $SNP_number_t\n";
				$start_t = $start if $start_t == 0;
				$end_t = $end;
				$SNP_number_t += $SNP_number;
				$average_SNPindex_t += $average_SNPindex;
				$pl_99_average_t +=$pl_99_average;
				$pl_95_average_t +=$pl_95_average;
				$ph_99_average_t +=$ph_99_average;
				$ph_95_average_t +=$ph_95_average;
				$new_window_number++;
				$average_SNPindex_t = sprintf "%.4f", $average_SNPindex_t/$SNP_number_t;
				$pl_99_average_t =$pl_99_average_t/$SNP_number_t;
				$pl_95_average_t =$pl_95_average_t/$SNP_number_t;
				$ph_99_average_t =$ph_99_average_t/$SNP_number_t;
				$ph_95_average_t =$ph_95_average_t/$SNP_number_t;
				print OUT "$opt_p$chr\t$new_window_number\t$start_t\t$end_t\t$average_SNPindex_t\t$SNP_number_t\t$pl_99_average_t\t$pl_95_average_t\t$ph_95_average_t\t$ph_99_average_t\n";
				#$average_SNPindex_t = abs($average_SNPindex_t);
				print ALL "$chr\t$start_t\t$end_t\t$average_SNPindex_t\t$pl_99_average_t\t$pl_95_average_t\t$ph_95_average_t\t$ph_99_average_t\n";
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
		&plot_R_cutoff($file, "$opt_o.$opt_p$chr");
	}
	print "done\n";
	
	# plot all the chr together
	print "Step5: plot the combined results.\n";
	&plot_combined_cutoff("$opt_o.ALL.window_SNPindex");
	
	
	my $time2 = time();
	my $time  = sprintf "%.4f", ( $time2 - $time1 ) / 60;
	print "$0 finished! $time mins elapsed.\n";
	
	sub split_genome_cutoff {
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
	
	sub plot_combined_cutoff {
		my $file = shift @_;
		my $R_script = <<END;
		data <- read.table("$file",header=TRUE)
		pdf(w=16, h=7,"$opt_o\_combined_chr.pdf")
		plot(1:nrow(data), data\$average_SNP_index, pch=16, cex=.8, col=2-data[,1]%%2, ylim=c(-1, 1), ylab="delta (SNPindex)", xlab="Chromosome", xaxt='n', xaxs='i')
		points(1:nrow(data), data\$pl_99, pch=16, cex=0.01, col="blue",xaxt='n', xaxs='i',type="o", lwd=0.3)
		points(1:nrow(data), data\$pl_95, pch=16, cex=0.01, col="green",xaxt='n', xaxs='i',type="o", lwd=0.3)
		points(1:nrow(data), data\$ph_99, pch=16, cex=0.01, col="blue",xaxt='n', xaxs='i',type="o", lwd=0.3)
		points(1:nrow(data), data\$ph_95, pch=16, cex=0.01, col="green",xaxt='n', xaxs='i',type="o", lwd=0.3)
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
	
	
	sub plot_R_cutoff {
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
		plot(pos/1000000, index, ylab="delta (SNPindex)", xlab="$opt_p$chr (MB)", xlim=c(0, x_max), ylim=c(-1, 1), type="o", lwd=3, pch=16 , cex=0.8, xaxt="n", yaxt="n", )
		points(pos/1000000, data\$pl_99, pch=16, cex=0.01, col="blue",xaxt='n', xaxs='i',type="o", lwd=0.3)
		points(pos/1000000, data\$pl_95, pch=16, cex=0.01, col="green",xaxt='n', xaxs='i',type="o", lwd=0.3)
		points(pos/1000000, data\$ph_99, pch=16, cex=0.01, col="blue",xaxt='n', xaxs='i',type="o", lwd=0.3)
		points(pos/1000000, data\$ph_95, pch=16, cex=0.01, col="green",xaxt='n', xaxs='i',type="o", lwd=0.3)
		axis(side=2, y_y, tcl=-0.2)
		axis(side=1, x, tcl=-0.2)
		for (row in 1:nrow(data)){
			if (data\$average_SNP_index[row] >= data\$ph_99[row]) {
				points(((data\$end[row]- data\$start[row])/2+data\$start[row])/1000000, data\$average_SNP_index[row], col="red", pch=16, cex=0.8)
			}
		}	
		#points(pos[smaller]/1000000, index[smaller], col="black", pch=16, cex=0.8)
		#points(pos[larger]/1000000, index[larger], col="red", pch=16, cex=0.8)
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
	
=p
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
=cut
}








sub usage {
  die(qq/
Usage:   BSAtools.pl <command> [<arguments>]\n
Command: SNPfromPool      boxcox transformation for phenotype (finished)
         VCFfilter        convert the genotype format to hapmap (finished);
         Cutoff           convert the vcf format to genotype file ;
         filter_and_plot  convert the hapmap format to ped file (finished);
         plotCutoff       convert the hapmap format to map file;

Notes: Commands with description endting with (*) may need bcftools
       specific annotations.
\n/);
}



	







