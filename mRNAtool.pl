#!/usr/bin/perl

=head1
 mRNAtool.pl -- tools for mRNA analysis
=cut
use strict;
use warnings;
use IO::File;
use FindBin;
use Statistics::Basic qw(:all nofill);
use Getopt::Std;
use Bio::SeqIO;
use Bio::SearchIO;

my $version = 0.1;
my $debug = 0;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);

unless (defined $options{'t'} ) { usage($version); }

# checking parameters
my %ss = ('SE'=>1,'SS'=>1,'PE'=>1,'PS'=>1);
$options{'s'} = "SS" unless defined $options{'s'};
die "[ERR]sequencing-method: $options{'s'}\n" unless defined $ss{$options{'s'}};

if (defined $options{'l'} && $options{'l'} ne "fr-firststrand" && $options{'l'} ne "fr-secondstrand" && $options{'l'} ne "fr-unstranded" ) {
	die "[ERR]library-type: $options{'l'}\n";
} else {
	$options{'l'} = "fr-firststrand";
}

if	($options{'t'} eq 'align')	{ rnaseq_align(\%options, \@ARGV); }	# parse multi dataset
elsif	($options{'t'} eq 'tport')	{ rnaseq_tport(\%options, \@ARGV); }	# parse multi dataset
elsif	($options{'t'} eq 'count')	{ rnaseq_count(\%options, \@ARGV); }	# parse multi dataset
elsif	($options{'t'} eq 'norm')	{ rnaseq_norm(\%options, \@ARGV);  }	# parse single dataset
elsif	($options{'t'} eq 'corre')	{ rnaseq_corre(\%options, \@ARGV); }	# parse single dataset
elsif	($options{'t'} eq 'unmap')	{ rnaseq_unmap(\%options, \@ARGV); }	# parse single dataset
elsif   ($options{'t'} eq 'isize')      { rnaseq_isize(\%options, \@ARGV); }    # parse single dataset, get insert size for paired end dataset
elsif   ($options{'t'} eq 'Dassembly')	{ rnaseq_Dassem(\%options, \@ARGV);}	# denovo assemble RNASeq reads
elsif   ($options{'t'} eq 'Rassembly')  { rnaseq_Rassem(\%options, \@ARGV);}    # denovo assemble RNASeq reads
elsif   ($options{'t'} eq 'mapping')	{ rnaseq_map(\%options, \@ARGV);   }	# align read to cDNA reference
elsif   ($options{'t'} eq 'ctgFeature')	{ rnaseq_ctgFeature(\%options, \@ARGV);}# generate feature bed for reads aligned to cDNA reference		
elsif	($options{'t'} eq 'blastn')	{ rnaseq_map(\%options, \@ARGV);   }	# blast assembled contigs to nt to remove contanmiantion
elsif   ($options{'t'} eq 'novoclean')	{ rnaseq_novoclean(\%options, \@ARGV);}  # clean assembled contigs using seqclean
elsif	($options{'t'} eq 'translate')	{ rnaseq_translate(\%options, \@ARGV);} # translate EST to protein
elsif   ($options{'t'} eq 'annotate')   { rnaseq_annotate(\%options, \@ARGV);}  # generate annotation command by AHRD
elsif   ($options{'t'} eq 'fixanno')  	{ rnaseq_fixanno(\%options, \@ARGV);} 	# fix annotation command by AHRD
elsif	($options{'t'} eq 'pipeline')	{ pipeline(); }				# print pipelines
else	{ usage($version); } 

#################################################################
# kentnf: subroutine for denovo pipeline			#
#################################################################

=head2
 rnaseq_novoclean: clean contamination after de novo assembly 
=cut
sub rnaseq_novoclean
{
	my ($options, $files) = @_;
	
	my $usage = qq'
USAGE: $0 -t novoclean Trinity_kn.fasta
* this tool is not a script like others, but a pipeline description.  


';
	print $usage and exit unless (@$files == 2);
	my $seqclean_bin = $FindBin::RealBin."/bin/seqclean/seqclean";
	die "[ERR]seqclean not exist\n" unless -s $seqclean_bin;

	print qq'
$usage

# A. remove contamination by blastn
\$blastTool.pl 
  [integrade honghe blast method to this method]

# B. remove reads without complete degradation for strand specific samples. 

# C. remove rRNA for de novo asembled transcripts:
$seqclean_bin $$files[1] -s $$files[0] -c 15 -l 200 -o seqclean_output_cycle1
	-s rRNA sequence file (rRNA_silva111)
	-l min seq length (seq will be discard short than -l)
	-c cpu number

# D. remove duplication
	\$ iAssembler.pl

';
	exit;
}
=head2
 rnaseq_Dassem: denovo assembly using trinity
=cut
sub rnaseq_Dassem
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 Dassembly [options] input1.fq input2.fq ...... | input1_r1_fq,input1_r2.fq input2_r1_fq,input2_r2.fq ...

example: Trinity --seqType fq --JM 50G --left EB_l1_1.clean.PE.fq,FB_l1_1.clean.PE.fq --right EB_l1_2.clean.PE.fq,FB_l1_2.clean.PE.fq --SS_lib_type RF --output assembly_k3   --min_kmer_cov 3 --CPU 8 --bflyCPU 8 --inchworm_cpu 8 --path_reinforcement_distance 25 > assembly_SS_k3.log

';
}

=head2
 rnaseq_Rassem: reference based assemly using cufflinks
=cut
sub rnaseq_Rassem
{

}

=head2
 rnaseq_mapping: map rnaseq reads to cDNA reference
=cut
sub rnaseq_mapping
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t mapping [options] reference.fasta [input1.fq input2.fq ... | input1_r1_fq,input1_r2.fq input2_r1_fq,input2_r2.fq ...]

	-p cpu
	-a BWA|bowtie
	-n edit distance (default:0)
	-m mode (1: keep unmapped reads, 2: generate aligned bam file)
	-l fr-firststrand, fr-secondstrand, fr-unstranded
	
';

	# check input files
	my @files = @$files;
	print $usage and exit if (scalar(@files) < 2);
	my $reference = shift @files;
	print "[ERR]no reference\n $usage" and exit unless -s $reference;
	my $err = 0;
	foreach my $f (@files) {
		my @a = split(/,/, $f);
		if (scalar(@a) == 1)	{ print "[ERR]no input read $f\n" unless -s $f; $err = 1; }
		elsif (scalar(@a) == 2) { print "[ERR]no input read1 $a[0]\n" unless -s $a[0]; print "[ERR]no input read2 $a[1]\n" unless -s $a[1]; $err = 1; }
		else { print "[ERR]input read $f\n"; $err = 1; }
		exit if $err;
	}

	# check input parameters
	my $distance = 0; $distance = $$options{'n'} if (defined $$options{'n'} && $distance > 0 && $distance < 4);
	my $cpu = 32; $cpu = $$options{'p'} if (defined $$options{'p'} && $$options{'p'} > 0);
	my $program = 'bowtie'; $program = $$options{'a'} if (defined $$options{'a'} && $$options{'a'} eq 'bwa');
	my $maxins = 100000;
	

	my $cmd;
	foreach my $f (@files) {
		my @a = split(/,/, $f);	
		if (scalar(@a) == 1) {
			my $unmap = $a[0].".unmap";
			my $out_sam = $a[0].".sam";
			$cmd = "bowtie -v $distance -k 1 -p $cpu --un $unmap -S $reference $a[0] $out_sam";
		} else {
			my $unmap = $a[0].".unmap";
			my $out_sam = $a[0].".sam";
			$cmd = "bowtie -v $distance -X $maxins -k 1 -p $cpu --un $unmap -S $reference -m1 $a[0] -m2 $a[1] $out_sam";
		}
	}
}

=head2
 rnaseq_ctgFeature: convert contigs to feature
=cut
sub rnaseq_ctgFeature
{
	my ($options, $files) = @_;	
	my $usage = qq'
USAGE: $0 -t ctgFeature contig.fasta > feature.bed

* after generate the feature, you could count the reads using count tool 

';
	print $usage and exit unless defined $$files[0];
	die "[ERR]file not exist $$files[0]\n" unless -s $$files[0];
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$$files[0]);
	while(my $inseq = $in->next_seq) {
		print $inseq->id,"\t0\t",$inseq->length,"\t", $inseq->id,"\t",$inseq->length,"\t+\n",
	}
}

=head2
 rnaseq_translate: translate transcript to protein
=cut
sub rnaseq_translate
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t translate [options] transcript.fasta > protein.fasta

	-m method used for translate 
		6frame: 6frame translate (default)
		3frame: 3frame translate (for strand-specific)
		estscan: translate using estscan (for denovo assembled transcripts)
	-n min length of protein (default: 100)
	-a [0 or 1]require start codon (M) for output protein
	-e [0 or 1]require stop codon (*) for output protein
    -b output best (longest) protein
	-t matrix_file for estscan

	* the output translated proteins need to be select for best result

';
	# check input files and parameters
	print $usage and exit unless defined $$files[0];
	print "[ERR]no input transcripts\n" and exit unless -s $$files[0];
	my $method = '6frame';
	if (defined $$options{'m'} && ($$options{'m'} eq '3frame' || $$options{'m'} eq 'estscan' || $$options{'m'} eq '6frame')) {
		$method = $$options{'m'};
	} else {
		print "[ERR]method $$options{'m'}\n" and exit;
	}

	my $min_len = 100;
	$min_len = $$options{'n'} if (defined $$options{'n'} && $$options{'n'} > 5);

	# array for order of 6 frames
	my @frames = ('0F','1F', '2F','0R','1R','2R');
	
	my ($require_start, $require_stop, $output_best) = (0, 0, 0);
	$require_start = 1 if (defined $$options{'a'} && $$options{'a'} eq '1');
	$require_stop = 1 if (defined $$options{'e'} && $$options{'e'} eq '1');
	$output_best = 1 if (defined $$options{'b'} && $$options{'b'} eq '1');

	if ($method eq '6frame' || $method eq '3frame')
	{
		# put the translated prot to this hash
		# key1: tid
		# key2: frame
		# array: [a_start, a_end, aseq, tstart, tend, tseq]
		my %trans_prot;

		my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$$files[0]);
		while(my $inseq = $in->next_seq)
		{
			my $sid = $inseq->id;
			my $sequence = $inseq->seq;
			# my $revcom_seq = reverse($sequence);
			# $revcom_seq =~ tr/ATCGNatcgn/TAGCNtagcn/;
			if ( $inseq->alphabet eq "rna" || $inseq->alphabet eq "dna" )
			{
				# translate to proteins
				my @prots;
				if ($method eq '6frame') {
					@prots = Bio::SeqUtils->translate_6frames($inseq);
				} elsif ($method eq '3frame') {
					@prots = Bio::SeqUtils->translate_3frames($inseq);
				} else {
					print $usage and exit;
				}

				# filter translated proteins
				for (my $i = 0; $i < @prots; $i++)
				{
					my $frame = $frames[$i];
					my $frame_num = $frame; $frame_num =~ s/F//;
					my $tseq = $prots[$i]->seq;		 # translated seq
					my $tid = $inseq->id."_".$frame; # translated ID -- seqID_frame

					#print ">$tid\n\t$frame_num$tseq\n"; exit;
					my @t = split(//, $tseq);
					my ($n_start, $n_end, $n_len, $a_len, $a_seq, $a_switch);
					$a_seq = ''; $a_switch = 0;

					for( my $j=0; $j<@t; $j++)
					{
						if ($t[$j] ne "*" && $a_switch == 0) {
							$a_switch = 1;
							$n_start = $j * 3 + 1 + $frame_num;
						}

						$a_seq.=$t[$j] if $a_switch == 1;

						if (($t[$j] eq "*" && $a_switch == 1) || 
							($j == scalar(@t)-1 && $a_switch == 1))
						{
							$a_switch = 0;
							
							if ($t[$j] eq "*") { $n_end = $j * 3 + 3 + $frame_num; }
							else			   { $n_end = $j * 3 + $frame_num; }
							
							# output protein seq
							# print "[ERR]no protein start $a_seq\n" and exit unless $a_seq =~ m/^M/;
							# print "[ERR]no protein end $a_seq\n" and exit unless $a_seq =~ m/\*$/;
                        	                        
							$a_len = length($a_seq) - 1;
							$n_len = $n_end - $n_start + 1;
                                        	        
							#if ($a_len >= $min_len) {
							#	print ">$tid $n_start-$n_end:$n_len translated to $a_len\n$a_seq\n";
							#}
							# push it to hash
							push(@{$trans_prot{$sid}{$frame_num}}, [$n_start, $n_end, $n_len, $a_len, $a_seq]);

							$a_seq = '';
							$n_start = '';
							$n_end = '';
						}
					}
				}
			}
			else
			{
				die "[ERR]input is not nucleotied $inseq->id\n";
			}
		}

		# output result
		
		foreach my $sid (sort keys %trans_prot) 
		{
			my $best_len = 0;
			my $best_m;
			my $best_frame;

			foreach my $frame_num (sort keys %{$trans_prot{$sid}}) 
			{
				my $frame = $frame_num."F";
				my @m = @{$trans_prot{$sid}{$frame_num}};
				#print "$sid $frame_num\n";
				
				foreach my $m ( @m ) 
				{
					my $n_start = $m->[0];
					my $n_end = $m->[1];
					my $n_len = $m->[2];
					my $a_len = $m->[3];
					my $a_seq = $m->[4];

					# skip 
					next if ($a_len < $min_len);
					next if ($require_start eq '1' && $a_seq !~ m/^M/);
					next if ($require_stop  eq '1' && $a_seq !~ m/\*$/);

					if ($a_len > $best_len) {
						$best_len = $a_len;
						$best_frame = $frame;
						$best_m = $m;
					}

					unless ($output_best eq '1') {
						print ">$sid $frame $n_start-$n_end:$n_len translated to $a_len\n$a_seq\n";
					}
				}
			}

			if ($output_best eq '1' && defined $best_m) {
				print ">$sid $best_frame $best_m->[0]-$best_m->[1]:$best_m->[2] translated to $best_m->[3]\n$best_m->[4]\n";
			}
		}

	} # end of 6frame and 3 frame method
	else
	{
		print qq'estscan \
  -M /home/kentnf/software/ESTScan2-2.1/Arabidopsis_thaliana.smat \
  -l 30 \
  -t mango_unigene_proteinA.fasta \
  -o mango_unigene_transA.fasta \
  input.fasta';
	}
}
=head2
 rnaseq_annotate: generate command for rnaseq contig annotation 
=cut
sub rnaseq_annotate
{
	my ($options, $files) = @_;
	my $usage = qq'
>> pipeline for annotation
1. AHRD, input file is input_ctg.fa
    \$ blastall -p blastx -i input_ctg.fa -o input_ctg_tr.pairwise -d uniprot_trembl_plant -e 0.0001 -v 200 -b 200 -m 0 -a 64
    \$ blastall -p blastx -i input_ctg.fa -o input_ctg_sp.pairwise -d uniprot_sprot.fasta -e 0.0001 -v 200 -b 200 -m 0 -a 64
    \$ blastall -p blastx -i input_ctg.fa -o input_ctg_at.pairwise -d TAIR10_pep_20110103_representative_gene_model_parse_desc -e 0.0001 -v 200 -b 200 -m 0 -a 64
    \$ mRNAtool.pl -t annotate -n 1 input_ctg.fa input_ctg_sp.pairwise input_ctg_at.pairwise input_ctg_tr.pairwise
    * the output file is input_ctg.ahrd.csv
    * the blast file should be input by order (SP, AT, and TR)
    * please use "-n" if the input seq is nucleotie and blastx 

2. GO annotation (using blast result from AHRD)
    \$ parse_blast.pl input_ctg_tr.pairwise 5 > input_ctg_tr.blast.table
    \$ parse_blast.pl input_ctg_sp.pairwise 5 > input_ctg_sp.blast.table
    \$ cat input_ctg_tr.blast.table input_ctg_sp.blast.table > input_ctg_uniport.blast.table
    \$ go_link_gene.pl input_ctg_uniport.blast.table parsed_go_mapping_file[option] > gene_GO
    \$ go_generate_associate.pl gene_GO organism > associate_file
    \$ go_enrichment.pl -i list_gene_id -a associate_file
    \$ go_slim.pl -a associate_file 

3. Pathway annotation
    * parse the AHRD csv file then fed to PathwayTools to perform analysis
 
';
	my @files = @$files;
	print $usage and exit unless scalar (@files) == 4;
	foreach my $f (@files) {
		print "[ERR]file not exist $f\n" and exit unless -s $f;
		# parse blast file if from code quest
	}

	my $output_file = $files[0];
	$output_file =~ s/\.fa$//;
	$output_file =~ s/\.fasta$//;
	$output_file .= ".ahrd.csv";

	my ($p1, $p2, $p3) = (0.5, 0.3, 0.2);
	if (defined $$options{'n'}) { ($p1, $p2, $p3) = (0.6, 0.4, 0.0); }

	my $ahrd_fdr = ${FindBin::RealBin}."/bin/AHRD";

	# generate temp xml file for AHRD
	my $temp_ahrd_yml = 'temp_ahrd.yml';
	my $fh = IO::File->new(">".$temp_ahrd_yml) || die $!;
	print $fh qq'
proteins_fasta: $files[0]
blast_dbs:
  swissprot:
    weight: 100
    file: $files[1]
    blacklist: $ahrd_fdr/blacklist_descline.txt
    filter: $ahrd_fdr/filter_descline_sprot.txt
    token_blacklist: $ahrd_fdr/blacklist_token.txt
    description_score_bit_score_weight: 0.2

  tair:
    weight: 50
    file: $files[2]
    blacklist: $ahrd_fdr/blacklist_descline.txt
    filter: $ahrd_fdr/filter_descline_tair.txt
    token_blacklist: $ahrd_fdr/blacklist_token.txt
    description_score_bit_score_weight: 0.4

  trembl:
    weight: 10
    file: $files[3]
    blacklist: $ahrd_fdr/blacklist_descline.txt
    filter: $ahrd_fdr/filter_descline_trembl.txt
    token_blacklist: $ahrd_fdr/blacklist_token.txt
    description_score_bit_score_weight: 0.4

#interpro_database: $ahrd_fdr/interpro_31.xml
#interpro_result: $ahrd_fdr/interpro_result.raw
#gene_ontology_result: $ahrd_fdr/go_results.csv
token_score_bit_score_weight: $p1
token_score_database_score_weight: $p2
token_score_overlap_score_weight: $p3
description_score_relative_description_frequency_weight: 0.6
output: $output_file
';

	# run AHRD
	run_cmd("java -Xmx20g -jar $ahrd_fdr/ahrd.jar $temp_ahrd_yml") unless -s $output_file;
}

=head2
 rnaseq_fixanno: fix ahrd annotation result
=cut
sub rnaseq_fixanno
{
	my ($options, $files) = @_;
        my $usage = qq'
USAGE: $0 -t fixanno [options] input_blast1 input_blast2 ... input_blastN
	-a ahrd_file (output of ahrd annotation)
	-e evalue (report blast hit report) [default: 1e-5]

* the output is ahrd_file.fix

';

	print $usage and exit unless $$files[0];
	foreach my $blast ( @$files ) {
		die "[ERR]blast file not exist: $blast\n" unless -s $blast;
	}


	my $evalue = 1e-5;
	$evalue = $$options{'e'} if (defined $$options{'e'} && $$options{'e'} > 0);

	my $ahrd_result;
	$ahrd_result = $$options{'a'} if defined $$options{'a'};
	die "[ERR]ahrd file not exist: $ahrd_result\n" unless -s $ahrd_result;
	$evalue = 'NA' if $ahrd_result;	# disable evalue cutoff

	# main: put blast hit to hash
	# key: protein/transcript ID
	# value: hit
	my %id_hit;
	my $id_hit = \%id_hit;
	foreach my $blast ( @$files )
	{
        	print "Parsing blast file: $blast ....\n";
	        $id_hit = parse_blast($blast, $id_hit, $evalue);
	}

	# parse ahrd result
	if (defined $ahrd_result) {

		my $out = IO::File->new(">".$ahrd_result.".fixed") || die $!;
		my $fh = IO::File->new($ahrd_result) || die $!;
		for (1 .. 3) {
	        	my $line = <$fh>;
		        print $out $line;
		}

		while(<$fh>)
		{
		        chomp;
		        # Protein-Accession     Blast-Hit-Accession     AHRD-Quality-Code       Human-Readable-Description      Interpro-ID (Description)       Gene-Ontology-ID (Name)
	        	my @a = split(/\t/, $_);
		        my $id = $a[0];
			$a[3] = "No hit" unless defined $$id_hit{$id};
		        print $out join("\t", @a)."\n";
		}
		$fh->close;
		$out->close;
	}
	print "done\n";
}

# parse blast reulst 
sub parse_blast
{
        my ($blast, $id_hit, $evalue) = @_;

        my ($hit_num, $query_name);
	my $query_hit_num = 0;
        my $in = Bio::SearchIO->new(-format=>'blast', -file=>$blast);
        while(my $result = $in->next_result)
        {
                $hit_num = 0;
                $query_name = $result->query_name;
                while(my $hit = $result->next_hit)
                {
			if ($evalue eq 'NA') {
				if ($query_name ne $hit->name) {
					$hit_num++;
				}
			}
			else
			{
				if (($query_name ne $hit->name) && ($hit->significance < $evalue)) {
                                	$hit_num++;
				}
                        }
                }

                if ($hit_num > 0)
                {
                        $$id_hit{$query_name} = 1;
			$query_hit_num++;
                }
        }

	print "No. of hit for $blast: $query_hit_num\n";
        return $id_hit;
}

#################################################################
# kentnf: subroutine for reference pipeline			#
#################################################################
=head2
 rnaseq_align: align rnaseq_read to reference genome, generate report
=cut
sub rnaseq_align
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 [options] input1.fq input2.fq ...... | input1_r1_fq,input1_r2.fq input2_r1_fq,input2_r2.fq ...
	
 -a  [String]	aligner name, tophat or hisat (default: tophat)
 -s  [String]   Sequencing method for read files: (required)
        PE (paired-end);
        SE (single-end);
        PS (paired strand-specific);
        SS (single strand-specific);
 -d  [String]   Indexed database of genome sequences (required)
 -l  [String]   Library type, fr-unstranded, fr-firststrand, fr-secondstrand (required)
 -p  [Integer]  number of CPUs used for megablast clustering (default = 1)
 -r  [Integer]  mate-inner-dist (default = 100)
 -v  [Integer]  mate-std-dev (default = 20)
 -n  [Integer]  read mismatches (default = 0)
 -g  [Integer]  max-multihits (default : 20)
 -f  [String]	gene feature file (optional)
 -k  [1 or 0]   skip align just print command (1 skip, 0 normal)

* input file could be fq, fa, or gzip file

';
	# check & correct parameters
	die "[ERR]no database\n$usage" unless defined $$options{'d'};
	my @index_files = ( 
		$$options{'d'}.".1.bt2", $$options{'d'}.".2.bt2", 
		$$options{'d'}.".3.bt2", $$options{'d'}.".4.bt2", 
		$$options{'d'}.".5.bt2", $$options{'d'}.".6.bt2",
		$$options{'d'}.".rev.1.bt2", $$options{'d'}.".rev.2.bt2",
		$$options{'d'}.".rev.5.bt2", $$options{'d'}.".rev.6.bt2" 
	);
	my $cpu = 20; $cpu = $$options{'p'} if defined $$options{'p'};
	my $mismatch = 1; $mismatch = $$options{'n'} if defined $$options{'n'};

	foreach my $index_file (@index_files) {
        	unless(-s $index_file) { die "[ERR] $index_file not exist\n"; }
	}
	
	# set the bin folder to ENV path
	$ENV{'PATH'} = ${FindBin::RealBin}."/bin".":".$ENV{'PATH'};

	my ($tophat_bin, $hisat_bin) = ('tophat2', 'hisat');
	if (-s "${FindBin::RealBin}/bin/tophat-2.0.11.Linux_x86_64/tophat2" ) {
        	$tophat_bin = "${FindBin::RealBin}/bin/tophat-2.0.11.Linux_x86_64/tophat2";
	} 
	#print "Topaht2:\n".$tophat_bin."\n";
	my $aligner = 'hisat';
	$aligner = 'tophat' if (defined $$options{'a'} && $$options{'a'} eq 'tophat');
	
	# check input files
	# key: output_suffix; value: input_file
	my @input_files;
	foreach my $f ( @$files ) 
	{ 
		if ($f =~ m/,/) 
		{
			my @a = split(/,/, $f);
			die "[ERR]file $a[0] or $a[1] not exist\n" unless (-s $a[0] && -s $a[1]);
			my $fm1 = check_seq_format($a[0]);
			my $fm2 = check_seq_format($a[1]);
			die "[ERR]file $a[0] and $a[1] have diff format\n" if $fm1 ne $fm2;

			push(@input_files, [$a[0],$a[1],$fm1]);
		}
		else
		{
			die "[ERR]file $f not exist\n" unless -s $f;
			my $fm1 = check_seq_format($f);
			push(@input_files, [$f,'',$fm1]);
		}
	}
	die "[ERR]no input files\n" if (scalar(@input_files) == 0);
	
	# main 
	# 1. create header of report file
	my $report_stat = "";

	if ($$options{'s'} =~ m/^S/) {
		$report_stat.="#sample\ttotal\tmapped\t%mapped\tmhit\t%mhit\n";
	} else {
		$report_stat.="#sample\ttotal\tmapped\t%mapped\tmhit\t%mhit";
		$report_stat.="\tLeftMap\t%LeftMap\tLeftMhit\t%LeftMhit\tRightMap\t%RightMap\tRightMhit\t%RightMhit\n";
	}

	# 2. tophat mapping
	foreach my $f (@input_files) 
	{
		my ($r1, $r2, $seq_type) = @$f;
		my ($input, $output, $format, $library, $report_file);

		# run tophat
		my $cmd_align = ''; 

		if ($aligner eq 'tophat') 
		{
			if ($r2) { $input = "$r1 $r2"; } else { $input = $r1; }
			$output = $r1."_tophat";
			$report_file = $output."/align_summary.txt";

			$cmd_align = "$tophat_bin -o $output --library-type $$options{'l'} -p $cpu ";
			$cmd_align.= "--segment-mismatches $mismatch --read-mismatches $mismatch --read-edit-dist $mismatch --read-gap-length $mismatch ";
			$cmd_align.= "--max-multihits 20 --segment-length 25 $$options{'d'} $input";
		} 
		else 
		{
			if ($r2) { $input = "-1 $r1 -2 $r2"; } else { $input = "-U $r1"; }
			$output = $r1."_hisat";

			if 	($seq_type eq 'fasta') { $format = '-f'; } 
			elsif	($seq_type eq 'fastq') { $format = '-q'; }
			else	{ die "[ERR]seq froamt $r1, $r2, $seq_type\n"; }

			if ($$options{'l'} eq 'fr-firststrand') {
				if ($r2) { $library = 'RF'; } else { $library = "R"; }
			} elsif ($$options{'l'} eq 'fr-secondstrand') {
				if ($r2) { $library = 'FR'; } else { $library = "F"; }
			} else {
				$library = '';		
			}

			$report_file = "$output.report.txt";

			my $penalty = 1;
			my $score_min = $penalty * $mismatch;

			if ($library) { $library = "--rna-strandness $library"; }
			$cmd_align = "$hisat_bin --time -p $cpu --no-unal $library $format ";
			$cmd_align.= "-k 20 --mp $penalty,$penalty --rdg 0,$penalty --rfg 0,$penalty --np $penalty --score-min C,-$score_min,0 --ignore-quals ";
			$cmd_align.= "-x $$options{'d'} $input -S $output.sam 1>$report_file 2>&1";
		}

		if (defined $$options{'k'} && $$options{'k'} == 1) {
			print $cmd_align."\n"; next;
		}
		run_cmd($cmd_align);
		
		# parse report file
		my $report_info = parse_align_report_file($report_file, $r1);
		$report_stat.= $report_info;
	}

	# output report information to file
	my $stat_file = "report_ralign_stat.txt";
	my $out = IO::File->new(">".$stat_file) || die $!;
	print $out $report_stat;
	$out->close;
}

=head2
 parse_align_report_file: child of rnaseq_align
=cut
sub parse_align_report_file 
{
	my ($report_file, $f) = @_;

	my $rinfo = `cat $report_file`;
	chomp($rinfo); my @r = split(/\n/, $rinfo);
	my $report_info = '';

	my ($total, $mapped, $mhit, $lefttotal, $leftmap, $leftmhit, $righttotal, $rightmap, $rightmhit, $unmap, $single_hit);
	if ($report_file =~ m/align_summary\.txt/)
	{
		if (scalar(@r) > 10)
                {
                        if ( $r[1] =~ m/Input\s+:\s+(\d+)/ )    { $lefttotal    = $1; }
                        if ( $r[2] =~ m/Mapped\s+:\s+(\d+)/ )   { $leftmap      = $1; }
                        if ( $r[3] =~ m/of these:\s+(\d+)/ )    { $leftmhit     = $1; }
                        if ( $r[5] =~ m/Input\s+:\s+(\d+)/ )    { $righttotal   = $1; }
                        if ( $r[6] =~ m/Mapped\s+:\s+(\d+)/ )   { $rightmap     = $1; }
                        if ( $r[7] =~ m/of these:\s+(\d+)/ )    { $rightmhit    = $1; }
                        if ( $r[10] =~ m/Aligned pairs:\s+(\d+)/ ) { $mapped = $1; }
                        if ( $r[11] =~ m/of these:\s+(\d+)/ )   { $mhit         = $1; }
                        print "[WARN]Left ($lefttotal) is not same as right ($righttotal)\n" if $lefttotal ne $righttotal;
                        $total = $lefttotal;
                        $report_info.="$f\t$total\t$mapped\t".sprintf("%.2f", ($mapped/$total)*100);
                        $report_info.="\t$mhit\t".sprintf("%.2f", ($mhit/$mapped)*100);
                        $report_info.="\t$leftmap\t".sprintf("%.2f", ($leftmap/$lefttotal)*100);
                        $report_info.="\t$leftmhit\t".sprintf("%.2f", ($leftmhit/$leftmap)*100);
                        $report_info.="\t$rightmap\t".sprintf("%.2f", ($rightmap/$righttotal)*100);
                        $report_info.="\t$rightmhit\t".sprintf("%.2f", ($rightmhit/$rightmap)*100)."\n";
                }
                else
                {
                        if ( $r[1] =~ m/Input\s+:\s+(\d+)/ )    { $total = $1; }
                        if ( $r[2] =~ m/Mapped\s+:\s+(\d+)/ )   { $mapped = $1;}
                        if ( $r[3] =~ m/of these:\s+(\d+)/ )    { $mhit = $1; }
                        $report_info.="$f\t$total\t$mapped\t".sprintf("%.2f", ($mapped/$total)*100);
                        $report_info.="\t$mhit\t".sprintf("%.2f", ($mhit/$mapped)*100)."\n";
                }
	} 
	else 
	{
	    if ($rinfo =~ m/concordantly/) 
	    {
		foreach my $r (@r) 
		{
			if ($r =~ m/(\d+) reads; of these:/) { $total = $1; }
			if ($r =~ m/(\d+) \(\S+\) aligned concordantly 0 times/) { $unmap = $1; }
			if ($r =~ m/(\d+) \(\S+\) aligned concordantly exactly 1 time/) { $single_hit = $1; }
			if ($r =~ m/(\d+) \(\S+\) aligned concordantly >1 times/) { $mhit = $1; }
		}
	    }
	    else
	    {
		foreach my $r (@r)
                {
                        if ($r =~ m/(\d+) reads; of these:/) { $total = $1; }
                        if ($r =~ m/(\d+) \(\S+\) aligned 0 times/) { $unmap = $1; }
                        if ($r =~ m/(\d+) \(\S+\) aligned exactly 1 time/) { $single_hit = $1; }
                        if ($r =~ m/(\d+) \(\S+\) aligned >1 times/) { $mhit = $1; }
                }
	    }
			
		$mapped = $single_hit + $mhit;
		$report_info.="$f\t$total\t$mapped\t".sprintf("%.2f", ($mapped/$total)*100);
		$report_info.="\t$mhit\t".sprintf("%.2f", ($mhit/$mapped)*100)."\n";
	}

	return $report_info;
}

=head2
 rnaseq_tport: parse tophat report file to generate result
=cut
sub rnaseq_tport
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 align_summary1.txt align_summary2.txt ... align_summary3.txt > report_tophat.txt

';
	print $usage and exit if scalar(@$files) == 0;

	# checking input files	
	foreach my $f (@$files) {
		print "[ERR]no file $f\n" and exit unless -s $f;
	}

	# parse report file
	my $report_tophat;
	my $report_line_num = 0;
	foreach my $f (@$files) {
		my $rinfo = `cat $f`;
		chomp($rinfo);
		my @r = split(/\n/, $rinfo);
		my ($total, $mapped, $mhit, $lefttotal, $leftmap, $leftmhit, $righttotal, $rightmap, $rightmhit);

		if (scalar @r == 14 || scalar @r == 13)
		{
			if ( $r[1] =~ m/Input\s+:\s+(\d+)/ )    { $lefttotal    = $1; }
			if ( $r[2] =~ m/Mapped\s+:\s+(\d+)/ )   { $leftmap      = $1; }
			if ( $r[3] =~ m/of these:\s+(\d+)/ )    { $leftmhit     = $1; }
			if ( $r[5] =~ m/Input\s+:\s+(\d+)/ )    { $righttotal   = $1; }
			if ( $r[6] =~ m/Mapped\s+:\s+(\d+)/ )   { $rightmap     = $1; }
			if ( $r[7] =~ m/of these:\s+(\d+)/ )    { $rightmhit    = $1; }
			if ( $r[10] =~ m/Aligned pairs:\s+(\d+)/ ) { $mapped = $1; }
			if ( $r[11] =~ m/of these:\s+(\d+)/ )   { $mhit         = $1; }
			print "[WARN]Left ($lefttotal) is not same as right ($righttotal)\n" if $lefttotal ne $righttotal;
			$total = $lefttotal;
			$report_tophat.="$f\t$total\t$mapped\t".sprintf("%.2f", ($mapped/$total)*100);
			$report_tophat.="\t$mhit\t".sprintf("%.2f", ($mhit/$mapped)*100);
			$report_tophat.="\t$leftmap\t".sprintf("%.2f", ($leftmap/$lefttotal)*100);
			$report_tophat.="\t$leftmhit\t".sprintf("%.2f", ($leftmhit/$leftmap)*100);
			$report_tophat.="\t$rightmap\t".sprintf("%.2f", ($rightmap/$righttotal)*100);
			$report_tophat.="\t$rightmhit\t".sprintf("%.2f", ($rightmhit/$rightmap)*100)."\n";
			$report_line_num = 14 if scalar(@r) > $report_line_num;
		}
		elsif (scalar @r == 5)
		{
			if ( $r[1] =~ m/Input\s+:\s+(\d+)/ )    { $total = $1; }
			if ( $r[2] =~ m/Mapped\s+:\s+(\d+)/ )   { $mapped = $1;}
			if ( $r[3] =~ m/of these:\s+(\d+)/ )    { $mhit = $1; }
			$report_tophat.="$f\t$total\t$mapped\t".sprintf("%.2f", ($mapped/$total)*100);
			$report_tophat.="\t$mhit\t".sprintf("%.2f", ($mhit/$mapped)*100)."\n";
			$report_line_num = 5 if scalar(@r) > $report_line_num;
                }
		else
		{
			print "[WARN]erro report format $f\n";
		}
	}

	my $report_head = '';
	if ($report_line_num == 14) {
		$report_head.="#sample\ttotal\tmapped\t%mapped\tmhit\t%mhit";
                $report_head.="\tLeftMap\t%LeftMap\tLeftMhit\t%LeftMhit\tRightMap\t%RightMap\tRightMhit\t%RightMhit\n";
	} else {
		$report_head.="#sample\ttotal\tmapped\t%mapped\tmhit\t%mhit\n";
	}

	$report_tophat = $report_head.$report_tophat;

	my $tophat_report_file = "report_tophat.txt";
        my $out = IO::File->new(">".$tophat_report_file) || die $!;
        print $out $report_tophat;
        $out->close;
}

=head2
 rnaseq_count: count and normalization for rnaseq analysis
 generate raw count and rpkm dataset for statistics analysis
=cut
sub rnaseq_count
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 [options] input1.bam input2.bam ......
        
 -s  [String]   Sequencing method for read files: (required)
        PE (paired-end);
        SE (single-end);
        PS (paired strand-specific);
        SS (single strand-specific);
 -f  [String]	gene feature file
 -a  [Num]      extend 5\' length of feature position
 -b  [Num]	extend 3\' length of feature position
 -p  [0 or 1]	polarity of extend
 -l  [String]   Library type, fr-unstranded, fr-firststrand, fr-secondstrand (required)
 -o  [String]	prefix of output file (default : exp)

';
	# check parameters
	die "[ERR]no input file\n$usage" if (scalar(@$files) == 0);
	die "[ERR]no feature file\n$usage" unless defined $$options{'f'} && -s $$options{'f'};
	my $feature_bed = $$options{'f'};
	
	die "[ERR]sequencing method\n$usage" unless defined $$options{'s'};
	my $sequencing = $$options{'s'};
	die "[ERR]library type\n$usage" unless defined $$options{'l'};
	my $library = $$options{'l'}; 

	my $output_prefix = "exp";
	$output_prefix = $$options{'o'} if defined $$options{'o'};
	my $output_exp_raw  = $output_prefix."_raw_count.txt";
	die "[ERR]raw count exist\n" if -s $output_exp_raw;

	my ($p5extend, $p3extend) = (0,0);
	$p5extend = $$options{'a'} if defined $$options{'a'} && $$options{'a'} > 0;
	$p3extend = $$options{'b'} if defined $$options{'b'} && $$options{'b'} > 0;

	my $polarity = '0';
	$polarity = 1 if defined $$options{'p'} && $$options{'p'} == 1;
	if ($polarity == 1) {
		die "[ERR]polarity only works for extend 5p or 3p\n" if ($p5extend > 0 && $p3extend > 0);
	}

	# check alignment file
	foreach my $f ( @$files ) { 
		die "[ERR]no alignment file $f\n" unless -s $f;
		die "[ERR]no bam file $f\n" unless $f =~ m/\.(bam|sam)$/;
	}

	# load feature into to hash
	# set bin size
	# change it base on genome size and density of feature
	# set it small will use more memory, set it big will take more time to locate read to feature
	my $bin_size = 1000;
	my %feature_struc = load_feature_bed($feature_bed, $bin_size, $p5extend, $p3extend, $polarity);

	# count expression for alignment
	# key: file_name, feature_id, strand (sense or antisense)
	# value: raw count
	my %exp;

	my %flag_plus = ('147' => 1, '99' => 1, '145' => 1, '97' => 1,
			'355' => 1, '403'=> 1, '353' => 1, '401'=> 1);
	my %flag_minus= ('163' => 1, '83' => 1, '161' => 1, '81' => 1,
			'339' => 1, '419'=> 1, '337' => 1, '417'=> 1);

	foreach my $f ( @$files ) {
		my $sam = $f; 
		$sam =~ s/\.bam/\.sam/ if $f =~ m/\.bam$/;
		run_cmd("samtools view -@ 20 -h -o $sam $f") unless -s $sam;

		# uniq reads for Paired reads
		my %PE_uniq = ();

		# code for count expression
		my $fh = IO::File->new($sam) || die $!;
		while(<$fh>)
		{
			chomp;
			next if $_ =~ m/^@/;
			my @a = split(/\t/, $_);
			my $chr = $a[2];
			my $read_start = $a[3];				# get read start for single end
			my $mapped_len = parse_cigar($a[5]);
			my $read_end = $a[3] + $mapped_len - 1;		# get read end for single end
			my $flag = $a[1];

			# uniq the Paired reads, count only once for each fragment
			if ($sequencing =~ m/^P/) {
				my $key = '';
				if ($a[6] eq '=') # paired reads aligned to same reference
				{ 
					$a[6] = $a[2];
					if ($a[3] < $a[7]) { 
						$key = $a[0]."\t".$a[2]."\t".$a[3]."\t".$a[6]."\t".$a[7];
						$read_start = $a[3];
						$read_end   = $a[3] + $a[8] - 1;
					} else { 
						$key = $a[0]."\t".$a[6]."\t".$a[7]."\t".$a[2]."\t".$a[3];
						$read_start = $a[7];
						$read_end   = $a[7] + abs($a[8]) - 1;
					}
				} 
				else  # paired reads aligned to diff reference
				{
					$key = $a[0]."\t".$a[2]."\t".$a[3]."\t".$a[6]."\t".$a[7];
				}

				if ( defined $PE_uniq{$key} ) { next; }
				else { $PE_uniq{$key} = 1; }
			}

			my $bin_start_pos = int($read_start/$bin_size);
			my $bin_end_pos   = int($read_end/$bin_size);
			my %feature_pos = (
				$bin_start_pos-1 => 1, 
				$bin_start_pos =>1, 
				$bin_start_pos+1 => 1,
				$bin_end_pos-1 => 1,
				$bin_end_pos => 1,
				$bin_end_pos+1 => 1
			);

			my $feature = '';
			foreach my $fpos (sort {$a<=>$b} keys %feature_pos ) {
				$feature .= $feature_struc{$chr."\t".$fpos} if defined $feature_struc{$chr."\t".$fpos};
			}
			next unless $feature; # the read do not overlapped with feature

			my %uniq_feature;
			my $uniq_feature = '';
			chomp ($feature);
			my @fea = split(/\n/, $feature);
			foreach my $fea (@fea) {
				unless ( defined $uniq_feature{$fea} ) {
					$uniq_feature{$fea} = 1;
					$uniq_feature .= $fea."\n";
				}
			}

			# main: count reads number with strand-specific method
			# just count reads at this step
			my $fstat = "";
			if ($sequencing =~ m/^S/) {
				$fstat = "P" if ($flag eq '0'  || $flag eq '256');
				$fstat = "M" if ($flag eq '16' || $flag eq '272');
			} else {
				$fstat = "P" if defined $flag_plus{$flag};
				$fstat = "M" if defined $flag_minus{$flag};
			}
			next unless $fstat; # skip if do n

			chomp($uniq_feature); @fea = split(/\n/, $uniq_feature);
			foreach my $fea (@fea) {
				my @b = split(/\t/, $fea);
				if (($read_start >= $b[1] && $read_start <= $b[2]) || ($read_end >= $b[1] && $read_end <= $b[2]))
				{
					if ( defined $exp{$f}{$b[3]}{$fstat} ) { $exp{$f}{$b[3]}{$fstat}++; } 
					else { $exp{$f}{$b[3]}{$fstat}=1; }
				}
			}
		}
		$fh->close;
		unlink($sam) if $f =~ m/\.bam$/; # remove sam for saving space
	}

	# output raw count
	my $sense_raw_table = "FeatureID";
	my $antisense_raw_table = "FeatureID";
	foreach my $f (@$files) {
		$sense_raw_table.="\t".$f;
		$antisense_raw_table.="\t".$f;
	}
	$sense_raw_table.="\n";
	$antisense_raw_table.="\n";

	my $fh = IO::File->new($feature_bed) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		my ($fid, $length, $strand) = ($a[3], $a[4], $a[5]);
		$sense_raw_table.=$fid;
		$antisense_raw_table.=$fid;

		foreach my $f (@$files) 
		{
			my ($ct_plus, $ct_minus) = (0,0);
			$ct_plus  = $exp{$f}{$fid}{'P'} if defined $exp{$f}{$fid}{'P'};
			$ct_minus = $exp{$f}{$fid}{'M'} if defined $exp{$f}{$fid}{'M'};
			my $ct = $ct_plus + $ct_minus;

			if ($sequencing =~ m/E$/) 
			{
				$sense_raw_table.= "\t".$ct;
			} 
			else # count raw count for strand specific read
			{
				if ( $library eq 'fr-firststrand' ) 
				{
					if ($strand eq "+") {
						$sense_raw_table.= "\t".$ct_minus;
						$antisense_raw_table.= "\t".$ct_plus;
					} else {
						$sense_raw_table.= "\t".$ct_plus;
						$antisense_raw_table.= "\t".$ct_minus;
					}
				} 
				else 
				{
					if ($strand eq "+") {
						$sense_raw_table.= "\t".$ct_plus;
						$antisense_raw_table.= "\t".$ct_minus;
					} else {
						$sense_raw_table.= "\t".$ct_minus;
						$antisense_raw_table.= "\t".$ct_plus;
					}
				}
			}
		}
		$sense_raw_table.="\n";
		$antisense_raw_table.="\n";
	}
	$fh->close;

	my $out1 = IO::File->new(">".$output_exp_raw) || die $!;
	print $out1 $sense_raw_table;
	$out1->close;

	if ($sequencing =~ m/S$/) {
		my $output_exp_raw_antisense = $output_prefix."_raw_count_anti.txt";
		my $out2 = IO::File->new(">".$output_exp_raw_antisense) || die $!;
		print $out2 $antisense_raw_table;
		$out2->close;
	}
}

=head2
 rnaseq_norm: normalization of expression file
=cut
sub rnaseq_norm
{
	my ($options, $file) = @_;
	my $usage = qq'
USAGE: $0 -t norm -f feature.bed -m method(default:1) -u report_tophat.txt -d decimal(default: 2) project_raw_count.txt > project_rpkm.txt
	normalization method: 
	1 -- RPKM
	2 -- RPM 

';
	# check input file
	print $usage and exit unless defined $$file[0];
	my $exp_file = $$file[0];
	print "[ERR]no exp file $exp_file\n" and exit unless -s $exp_file;
	print "[ERR]no feature file\n$usage" and exit unless defined $$options{'f'};
	print "[ERR]no tophat report file\n$usage" and exit unless defined $$options{'u'};

	my $decimal = 2;
	$decimal = $$options{'d'} if defined $$options{'d'} && $$options{'d'} > 0;

	my $method = '1';
	$method = $$options{'m'} if defined $$options{'m'} && $$options{'m'} == 2;

	# check libsize and sample name
	my %sample_libsize;
	my $fh1 = IO::File->new($$options{'u'}) || die $!;
	<$fh1>;	# skip the title
	while(<$fh1>)
	{
		chomp;
		my @a = split(/\t/, $_);
		$sample_libsize{$a[0]} = $a[2];
	}
	$fh1->close;

	my $sname = `head -n 1 $exp_file`;
	chomp($sname);
	my @s = split(/\t/, $sname); shift @s;
	foreach my $s (@s) {
		print "[ERR]no libsize for sample: $s\n" and exit unless defined $sample_libsize{$s};
	}

	# load feature length from bed file
	my %feature_length;
	my $fh2 = IO::File->new($$options{'f'}) || die $!;
	while(<$fh2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		# SL2.40ch00	16437	18189	Solyc00g005000.2.1	1693	+
		my @a = split(/\t/, $_);
		$feature_length{$a[3]} = $a[4];
	}
	$fh2->close;

	# normlization using RPKM/FPKM method
	my $fh3 = IO::File->new($exp_file) || die $!;
	my $title = <$fh3>; chomp($title);
	print $title."\n";
	my @title = split(/\t/, $title);
	while(<$fh3>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		print "[ERR]no feature length $a[0]\n" and exit unless defined $feature_length{$a[0]};
		my $length = $feature_length{$a[0]};
		my $output = $a[0];

		for(my $i=1; $i<@a; $i++) 
		{
			print "[ERR]no libsize for sample $title[$i]\n" and exit unless defined $sample_libsize{$title[$i]};
			my $lib_size = $sample_libsize{$title[$i]};
			my $norm;
			if ($method == 1) { 
				$norm = ($a[$i] * 1000 * 1000000) / ($length * $lib_size);
			} elsif ($method == 2) {
				$norm = ($a[$i] * 1000000) / $lib_size;
			}
			$norm = sprintf("%.$decimal"."f" , $norm) if $norm > 0;
			$output.="\t".$norm;
		}

		print $output."\n";
	}
}

=head2
 rnaseq_corre: correlation analysis for rnaseq dataset
=cut
sub rnaseq_corre
{
	my ($options, $file) = @_;

	my $usage = qq'
USAGE: $0 -t corre [options] input_exp.txt > output_corre.txt

';

	print $usage and exit unless defined $$file[0];
	my $exp_rpkm_file = $$file[0];
	print "[ERR]no exp file $exp_rpkm_file\n" and exit unless -s $exp_rpkm_file;

	if (scalar(@$file) == 2) 
	{

	}
	elsif (scalar(@$file) == 1) 
	{
		my %sample_exp;
		my $fh = IO::File->new($exp_rpkm_file) || die $exp_rpkm_file;
		# parse first title line
		my $title = <$fh>; chomp($title);
		my @title = split(/\t/, $title);
		while(<$fh>)
		{
        		chomp;
			my @a = split(/\t/, $_);
			for(my $i=1; $i<@a; $i++)
			{
				my $tt = $title[$i];
				if (defined $sample_exp{$tt})   { $sample_exp{$tt}.="\t".$a[$i]; }
				else                            { $sample_exp{$tt} = $a[$i]; }
	       	 	}
		}
		$fh->close;

		#################################################################
		# produce correlation for every line data                       #
		#################################################################

		print "Sample";
		shift @title;
		foreach my $tt ( @title ) { print "\t".$tt; }
		print "\n";

		foreach my $tt_i ( @title )
		{
	        	print $tt_i;
			my @exp_i = split(/\t/, $sample_exp{$tt_i});

			foreach my $tt_j ( @title )
			{
	                	my @exp_j = split(/\t/, $sample_exp{$tt_j});
				my $correlation = corr([@exp_i], [@exp_j]);
				print "\t".$correlation;
			}
			print "\n";
		}

	} 
	else 
	{

	}
}

=head2
 load_feature_bed : load feature informat with bed format
 data structure: 
 key: chr position(kb) (chr:1)
 value: structure of feature in the key position
=cut
sub load_feature_bed
{
	my ($feature_bed, $bin_size, $p5extend, $p3extend, $polarity) = @_;

	# set bin size
	# change it base on genome size and density of feature
	# set it small will use more memory, set it big will take more time to locate read to feature
	# my $bin_size = 1000;

	# put feature to feature_struc hash
	my %feature_struc;
	my ($ref, $start, $end, $strand, $feature_id, $length, $kpos);
	my $fh = IO::File->new($feature_bed) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		($ref, $start, $end, $feature_id, $length, $strand) = ($a[0], $a[1]+1, $a[2], $a[3], $a[4], $a[5]);

		# works for count p5extend + genebody + p3extend
		if ($strand eq '-') {
			# works for polarity
			if ($polarity == 1) {
				if ($p3extend > 0) {
					$start = $start - $p3extend;
					$end = $a[1];
                } else {
					$start = $a[2]+1;
					$end = $end + $p5extend;
				}							
			} else {
				$start = $start - $p3extend;
				$end = $end + $p5extend;
			}
		} 
		else 
		{
			# works for polarity
			if ($polarity == 1) {
				if ($p5extend > 0) {
					$start = $start - $p5extend;
					$end = $a[1];
				} else {
					$start = $a[2]+1;
					$end = $end + $p3extend;
				}
			}
			else 
			{
				$start = $start - $p5extend;
				$end = $end + $p3extend;
			}
		}

		my $new_line = "$ref\t$start\t$end\t$feature_id\t$length\t$strand";

		for(my $p = $start; $p < $end; $p = $p + $bin_size)
		{
			$kpos = int( $p / $bin_size );
			if ( defined $feature_struc{$ref."\t".$kpos} ) {
				$feature_struc{$ref."\t".$kpos}.= "$new_line\n";
			} else {
				$feature_struc{$ref."\t".$kpos} = $new_line."\n";
			}
		}
	}
	$fh->close;

	return %feature_struc;
}

#################################################################
# toolbox for rnaseq analysis					#
#################################################################
=head2
 rnaseq_unmap: get unmapped read for checking 
=cut
sub rnaseq_unmap
{
	my ($options, $file) = @_;
	
	my $usage = qq'
USAGE: $0 -t unmap [options] bam,read ... | bam,read_1,read2 ... 
	
	-u  [1-2] level for paired end mapped
	    1. only paired aligned
	    2. paired align + single align(paired)

	* the input could be sam or bam
';
	if (scalar(@$file) == 0) {
		print $usage; exit;
	}

	my $u = 2;
	$u = $$options{'u'} if (defined $$options{'u'} && $$options{'u'} > 0 && $$options{'u'} < 3);

	# check input files
	my ($format1, $format2);
	foreach my $set (@$file) {
		my @a = split(/,/, $set);
		if (($a[0] =~ m/\.bam$/ || $a[0] =~ m/\.sam$/) && -s $a[0]) {}
		else {
			print "[ERR]no bam/sam file: $set\n";
		}

		if (scalar(@a) == 2) {
			print "[ERR]no read1 file: $set\n" and exit unless -s $a[1];
			$format1 = check_seq_format($a[1]);
			print "[ERR]seq format: $format1 for $a[1]\n" and exit if ($format1 ne 'fastq' && $format1 ne 'fasta');
		} elsif (scalar(@a) == 3) {
			print "[ERR]no read1 file: $set\n" and exit unless -s $a[1];
			print "[ERR]no read2 file: $set\n" and exit unless -s $a[2];
			$format1 = check_seq_format($a[1]);
			$format2 = check_seq_format($a[2]);
			print "[ERR]seq format: $format1 for $a[1]\n" and exit if ($format1 ne 'fastq' && $format2 ne 'fasta');
			print "[ERR]seq format: $format2 for $a[2]\n" and exit if ($format2 ne 'fastq' && $format2 ne 'fasta');
		} else {
			print "[ERR]set $set\n" and exit;
		}
	}

	# main: get unmapped reads
	foreach my $set (@$file) 
	{
		my @a = split(/,/, $set);
		my $bam = $a[0];
		my $sam = $bam;
		if ($bam =~ m/\.bam$/) {
			$sam =~ s/\.bam$/\.sam/;
			run_cmd("samtools view $bam > $sam");
		}

		# load reads flag to hash
		my %mapped_read;
		my $fh1 = IO::File->new($sam) || die $!;
		while(<$fh1>)
		{
			chomp;
			next if $_ =~ m/^@/;
			my @a = split(/\t/, $_);
			if ($a[1] & 0x1) # paired-end
			{
				if ($u == 1) 
				{
					if (($a[1] & 0x4) || ($a[1] & 0x8)) {}
					else { $mapped_read{$a[0]} = 1; }
				}
				elsif ($u == 2)
				{
					$mapped_read{$a[0]} = 1 unless ($a[1] & 0x4);
				}
			}
			else		 # single-end
			{
				$mapped_read{$a[0]} = 1 unless ($a[1] & 0x4);
			}
		}
		$fh1->close;

		print "No. of mapped reads: ".scalar(keys(%mapped_read))."\n";

		if (scalar(@a) == 2)
		{
			my $unmap_file = $a[1].".unmap";
			my $out = IO::File->new(">".$unmap_file) || die $!;
			my $fh2;
			if ($a[1] =~ m/\.gz$/) {
				$fh2 = IO::File->new("gunzip -c $a[1] | ") || die $!;
			} else {
				$fh2 = IO::File->new($a[1]) || die $!;
			}
			while(<$fh2>)
			{
				my $id = $_; chomp($id); my $format;
				if 	($id =~ m/^>/) { $format = 'fasta'; $id =~ s/^>//; }
				elsif	($id =~ m/^@/) { $format = 'fastq'; $id =~ s/^@//; }
				else	{ die "[ERR]seq format $id\n"; }
				my $tid = $id; $tid =~ s/ .*//;
				$tid =~ s/\/\d$// if $tid =~ m/\/\d$/;
				my $seq = <$fh2>; chomp($seq);
				if ($format eq 'fastq') { <$fh2>; <$fh2>; }
				print $out ">$id\n$seq\n" unless defined $mapped_read{$tid};
			}
			$fh2->close;
			$out->close;
		}
		elsif (scalar(@a) == 3)
		{
			my $unmap_file1 = $a[1].".unmap";
			my $unmap_file2 = $a[2].".unmap";
                        my $out1 = IO::File->new(">".$unmap_file1) || die $!;
			my $out2 = IO::File->new(">".$unmap_file2) || die $!;

			my $fh2;
			if ($a[1] =~ m/\.gz$/) {
				$fh2 = IO::File->new("gunzip -c $a[1] | ") || die $!;
			} else {
				$fh2 = IO::File->new($a[1]) || die $!;
			}
			while(<$fh2>)
			{
				my $id = $_; chomp($id); my $format;
                                if      ($id =~ m/^>/) { $format = 'fasta'; $id =~ s/^>//; }
                                elsif   ($id =~ m/^@/) { $format = 'fastq'; $id =~ s/^@//; }
                                else    { die "[ERR]seq format $id\n"; }
                                my $tid = $id; $tid =~ s/ .*//;
				$tid =~ s/\/\d$// if $tid =~ m/\/\d$/;
                                my $seq = <$fh2>; chomp($seq);
                                if ($format eq 'fastq') { <$fh2>; <$fh2>; }
				print $out1 ">$id\n$seq\n" unless defined $mapped_read{$tid};
			}
			$fh2->close;

			my $fh3;
			if ($a[2] =~ m/\.gz$/) {
				$fh3 = IO::File->new("gunzip -c $a[2] | ") || die $!;
			} else {
				$fh3 = IO::File->new($a[2]) || die $!;
			}
			while(<$fh3>)
			{
				my $id = $_; chomp($id); my $format;
                                if      ($id =~ m/^>/) { $format = 'fasta'; $id =~ s/^>//; }
                                elsif   ($id =~ m/^@/) { $format = 'fastq'; $id =~ s/^@//; }
                                else    { die "[ERR]seq format $id\n"; }
                                my $tid = $id; $tid =~ s/ .*//;
				$tid =~ s/\/\d$// if $tid =~ m/\/\d$/;
                                my $seq = <$fh3>; chomp($seq);
                                if ($format eq 'fastq') { <$fh3>; <$fh3>; }
				print $out2 ">$id\n$seq\n" unless defined $mapped_read{$tid};
			}
			$fh3->close;
			
			$out1->close;
			$out2->close;
		}
		else
		{
			print "[WARN]set has error: $set\n";
		}

		## unlink($sam); # do not delete sam files
	}
}

=head2
 check_seq_format: check seq format 
=cut
sub check_seq_format
{
	my $seq_file = shift;
	my $format;
	my $seq_top = `less $seq_file | head -n 8`;	# could read gzip file 
	chomp($seq_top); 
	my @a = split(/\n/, $seq_top);	
	if ($a[0] =~ m/^>/ && $a[2] =~ m/^>/ && $a[4] =~ m/^>/ && $a[6] =~ m/^>/) {
		$format = 'fasta';
	} elsif ( $a[0] =~ m/^@/ && $a[2] =~ m/^\+/ && $a[4] =~ m/^@/ && $a[6] =~ m/^\+/) {
		$format = 'fastq';
	}
	return $format;
}

=head2
 parse_flag: parse sam flag
=cut
sub parse_flag
{
	my $flag_num = shift;
	my @flag_stat = '';
	if ($flag_num & 0x1)    { push(@flag_stat, "paired"); }
	if ($flag_num & 0x2)    { push(@flag_stat, "GoodAligned"); } 	
	if ($flag_num & 0x4)    { push(@flag_stat, "unmap1stSeg"); }
	if ($flag_num & 0x8)    { push(@flag_stat, "unmap2ndSeg"); }
	if ($flag_num & 0x10)   { push(@flag_stat, "1stSEQ=RC"); }
	if ($flag_num & 0x20)   { push(@flag_stat, "2ndSEQ=RC"); }
	if ($flag_num & 0x40)   { push(@flag_stat, "Left"); }
	if ($flag_num & 0x80)   { push(@flag_stat, "Right"); }
	if ($flag_num & 0x100)  { push(@flag_stat, "Mhit"); }
	if ($flag_num & 0x200)  { push(@flag_stat, "lowQC"); }
	if ($flag_num & 0x400)  { push(@flag_stat, "Dup"); }
	if ($flag_num & 0x800)  { push(@flag_stat, "supplementary alignment"); } 
	return @flag_stat;
}

=head2
 parse_cigar: parse cigar to get mapped length 
=cut
sub parse_cigar
{
	my $cigar = shift;

	my $str_len = length($cigar);

	my $num = "";; my $mapped_length = 0;

	for(my $i=0; $i<$str_len; $i++)
	{
		my $str = substr($cigar, $i, 1);

		if ($str =~ m/\d+/)
		{
			$num = $num.$str;
		}
		elsif ($str eq "M" || $str eq "N" || $str eq "I")
		{
			$mapped_length = $mapped_length + $num;
			$num = "";
			
		}
		elsif ($str eq "D") 
		{
			$num = "";
		}
	}

	return $mapped_length;
}

=head2
 run_cmd: run other program using tophat
=cut
sub run_cmd
{
	my $cmd = shift;
	print $cmd."\n" and return(1) if $debug;
	system($cmd) && die "[ERR]cmd: $cmd\n";
}

=head2
 usage: print usage info
=cut
sub usage
{
	my $version = shift;
	print qq'
USAGE: $0 -t [tool] [options] input file

	align		align read to reference using tophat
	tport		summary tophat report info to table
	count		count aligned reads for each feature
	norm		normalization of expression file
	corre		generate correlation tables
	unmap		extract unampped reads (tophat could do it)
	denovo		denovo assembly using Trinity
	ctgFeature	generate feature file from contigs
	annotate	function annotation by AHRD (blast is not inlude)	
	fixanno		fix ahrd annotation or report blast hits
	novoclean	remove contamination after denovo assembly
	translate	translate EST to protein
	

'; 

	exit;
}

=head2
 pipeline: print pipeline
=cut
sub pipeline
{
	print qq'
>> RNASeq reference based pipeline:

1. align RNASeq read to reference (Tophat2)
   \$mRNAtool.pl -t align -s PS -d /home/database/tomato_genome -l fr-firststrand -p 24 -n 2 inputA_1.fastq,inputA_2.fastq ... inputN_1.fastq,inputN_2.fastq
   \$mRNAtool.pl -t tport *_tophat/*.txt > report_tophat.txt

2. count aligned reads for each feature, and normalization
   \$mRNAtool.pl -t count -s PS -l fr-firststrand -f feature.bed -o project_exp inputA.bam ... inputN.bam
   \$mRNAtool.pl -t norm  -f feature.bed -u report_tophat.txt project_raw_count.txt > project_rpkm.txt

3. correlation
   \$mRNAtool.pl -t corre project_rpkm.txt > project_corr.txt

4. statistics analysis

>> RNASeq denovo pipeline

1. assembly RNASeq read to contigs (Trinity)

2. align RNASeq read to reference (bowtie2/bwa)

3. count aligned reads for each contigs

4. 5. normalization and correlation like reference based 3. 4.

>> RNASeq debug utils

A. get unmapped reads
   \$mRNAtool.pl -t unmap -u 2 input_bam,input_read_r1,input_read_r2 ...


';
	exit;
}

