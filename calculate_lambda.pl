#!/usr/bin/env perl
use warnings;
use strict;
use Statistics::R;

my $association = $ARGV[0];             # EMMAX ps file, full path
die "Usage: $0 EMMAX.ps\n" if @ARGV < 1;

print "processing $association ...";
# transform file
my $tmp = "$association.tmp";
open AS, "<$association";
open OUT, ">$tmp";
while (<AS>) {
	chomp;
	my @array = (split /\s+/, $_);
	# remove NA
	next unless $array[3] =~ /\d/;
	my ($chr, $pos) = (split /\_/, shift @array);
	print OUT "$chr\t$pos\t" . join("\t", @array) . "\n" if $pos;
	print OUT "$chr\t" . join("\t", @array) . "\n" unless $pos;
}
close AS;
close OUT;
`mv $tmp $association`;

my $lambda = "$association.lambda";
# calculate lambda
#$tmp = int(rand(10000));
my $R = Statistics::R->new();
$R->startR;
$R->send(qq`library(GenABEL)`);
$R->send(qq`data <- read.table("$association", header = FALSE)`);
$R->send(qq`write.table(estlambda(data[,5], filter = F, plot = T, method = "median")\$estimate, file = "$lambda", row.names = FALSE, col.names = FALSE)`);
$R->stopR;

print " Done\n";


=pod
input header:
SNP ID
Beta (1 is effect allele)
SE(beta)
p-value

output header:
Chr
pos
Beta (1 is effect allele)
SE(beta)
p-value
