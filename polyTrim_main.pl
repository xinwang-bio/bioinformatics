#!/usr/bin/perl -w
use strict;
use File::Basename;
use File::Spec;
use Getopt::Long;

my $window;
my $mismatch;
my $threads;



GetOptions(
"threads|t:i"=>\$threads,
"window|w:i"=>\$window,
"mismatch|m:i"=> \$mismatch,
);

die(qq/Usage: poly_main.pl [option] 1.fastq 2.fastq.gz 3.fastq ... \n


Options: -w|threads 	number of threads [40]\n
         -w|window      window used for PolyA\/T detecting [15]\n
         -m|mismatch    mismatch for the window [1]\n
/) unless ( $ARGV[0]);



$window = $window || 15;
$mismatch =$mismatch || 1;
$threads = $threads || 40;

	

my $path_curf = File::Spec->rel2abs(__FILE__);
#print "C PATH = ",$path_curf,"\n";
my ($vol, $dirs, $file) = File::Spec->splitpath($path_curf);
#print "C Dir = ", $dirs,"\n";


#my ($file1,$out) = @ARGV;
my $time1 = time();

#elsif($type eq "PE"){
#	push @file,$file1;
#	push @file,$file2;
#}

##split file into tmp file
if(-e $ARGV[0]){
	#die "Input files $ARGV[0] not found\n";
}
else{
	die "Input files $ARGV[0] not found\n";
}
for my $file(@ARGV){
	my $time3 = time();
	my $reads_number;
	if ($file =~ m/\.gz$/) {
		open (IN,'-|', "gzip -cd $file") or die;
		$reads_number = `gzip -cd $file | wc -l`;
	}
	else{
		open (IN, '<', "$file") or die; 
		$reads_number = `wc -l < $file`;
	}
	#open IN,"<",$file;
	 
	$reads_number = $reads_number/4;
	print "$file total Reads Number\:$reads_number\n";
	my $reads_split_number = $reads_number;
	my $test =int($reads_number/($threads-1));
	#print "$test\n";
	if((int($reads_number/($threads-1))) >=1){
		$reads_split_number = int($reads_number/($threads-1));
		
	}
	else{
		$reads_split_number = $reads_number;
	
	}
	print "$file Reads number in each thread\:$reads_split_number\n";
	my $file_index =1;
	my $reads_count =0;
	open OUT,">","$file.tmp$file_index";
	while(<IN>){
		chomp;
		my $name = $_;
		my $seq = <IN>;
		chomp($seq);
		my $third = <IN>;
		chomp($third);
		my $quality = <IN>;
		chomp($quality);
		$reads_count++;
		if($reads_count >$reads_split_number){
			$reads_count =0;
			print OUT"$name\n$seq\n$third\n$quality\n";	
			close OUT;
			$file_index++;
			open OUT,">","$file.tmp$file_index";
		}
		else{
			print OUT"$name\n$seq\n$third\n$quality\n";	
		}
	}
	close IN;
	close OUT;
	open SPLIT,">","$file\_split.sh";
	for my $run(1..$file_index){
		#print "$file_index\n";
		print SPLIT "perl $dirs/trimX.pl $file.tmp$run $file.tmp$run.trim $window $mismatch \& \n";
	}
	print SPLIT"wait\n";
	system "sh $file\_split.sh";
	my $merge = "cat";
	for my $run(1..$file_index){
		$merge = $merge." $file.tmp$run.trim";
	}
	my $output = $file;
	$output =~s/\.fastq.gz//;
	$output =~s/\.fq.gz//;
	$output =~s/\.fastq//;
	$merge= $merge." >$output.rmPolyAT.fastq";
	system "$merge";
	system "rm $file.tmp\*";
	system "rm $file\_split.sh";
	my $time4 = time();
	my $file_time  = sprintf "%.4f", ( $time4 - $time3 ) / 60;
	print "$file finished! $file_time mins elapsed.\n";
}





my $time2 = time();
my $time  = sprintf "%.4f", ( $time2 - $time1 ) / 60;
print "$0 finished! $time mins elapsed.\n";





