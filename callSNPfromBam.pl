#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;

our %opt;

getopt('r:w:m:q:Q:p:d:o:', \%opt);

unless ( $opt{r} && $opt{w} && $opt{m} && $opt{o} ) {
	&usage;
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
print "samtools-0.19 mpileup -q $opt{q} -Q $opt{q} -f $opt{r} $opt{m} $opt{w}\n";
print "Parsing the output in real time...\n";
open IN,"samtools-0.19 mpileup -q $opt{q} -Q $opt{q} -f $opt{r} $opt{m} $opt{w} |" or die $!;
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


sub usage {
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
