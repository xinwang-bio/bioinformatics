#!/usr/bin/perl

=head

 AStool.pl -- generate AS event result and stat the result

 20150719: ass AS differential analysis 
 20150419: add AS_splice_site_stat.pl to AStool
 20150418: change and combine script to AStool
 20140710: init
=cut

use strict;
use warnings;
use FindBin;
use IO::File;
use Getopt::Std;
use Bio::SeqIO;
use Text::NSP::Measures::2D::Fisher::twotailed;

my $version = 0.1;
my $debug = 0;
my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }

if		($options{'t'} eq 'analyze')	{ as_event_analyze(\%options, \@ARGV); }	# as event analyze
elsif	($options{'t'} eq 'spsites')	{ as_splice_site(\%options, \@ARGV); }		# as splice site
elsif	($options{'t'} eq 'stat')		{ as_event_stat(@ARGV); }		# count the number of AS events
elsif	($options{'t'} eq 'specific')	{ as_event_specific(\%options, \@ARGV); }	# find specific AS event to sample, not done
elsif   ($options{'t'} eq 'isodiff')	{ as_iso_diff(\%options, \@ARGV); }   		# find changed AS
else	{ usage($version); }

=head2
 as_iso_diff: find differentially expressed isoforms 
=cut
sub as_iso_diff
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t iso_ASI iso_diff iso_exp 

* iso_ASI is the output of AStalavista
* iso_diff is the output of cuffdiff
* iso_exp is the output of cuffnorm
';
	print $usage unless @$files == 3;
	my ($gtf_file, $diff_file,  $exp_file) = @$files;

	# load paired transcripts with AS events to hash
	my %pair_tid;
	
	my $f0 = IO::File->new($gtf_file) || die $!;
	while(<$f0>) {
		chomp;
		my @a = split(/\t/, $_);
		if ($a[8] =~ m/transcript_id "(\S+)";/) {
			my $t = $1;
			my @b = split(/,/, $t);
			die "[ERR]TID\n" unless (scalar(@b) == 2);
			my @c1 = split(/\//, $b[0]);
			my @c2 = split(/\//, $b[1]);
			foreach my $c1 (@c1) {
				foreach my $c2 (@c2) {
					$pair_tid{$c1."#".$c2} = 1;
					$pair_tid{$c2."#".$c1} = 1;
				}
			}
		}
	}
	$f0->close;

	my $pair_num = scalar(keys(%pair_tid))/2;
	print "No. of paired isoforms with AS events $pair_num\n";

	# load iso diff table to hash, the table is generated from cuffdiff
	my %sample_cmb;	# combination of sample
	my %gid_tid;	# key: gid, value: array of tids
	my %tid_sig;	# key1: tid, key2: sample1 # sample2, value: sig (yes or no)
	my $f1 = IO::File->new($diff_file) || die $!;
	<$f1>;
	while(<$f1>) {
		chomp;
		my @a = split(/\t/, $_);
		my ($tid, $gid, $sample1, $sample2, $sig) = ($a[0], $a[1], $a[4], $a[5], $a[13]);
		push(@{$gid_tid{$gid}}, $tid);
		$tid_sig{$tid}{$sample1."#".$sample2} = $sig;
		$sample_cmb{$sample1."#".$sample2} = 1 unless defined $sample_cmb{$sample1."#".$sample2};
	}
	$f1->close;

	#print scalar(keys(%gid_tid))."\n";
	#print scalar(keys(%tid_sig))."\n";
	print "No. of sample comparison ".scalar(keys(%sample_cmb))."\n";

	# load expression to hash
	# key1: iso ID, key2: sample, value expression
	my %iso_exp; 
	my $f2 = IO::File->new($exp_file) || die $!;
	my $title = <$f2>; chomp($title); 
	my @t = split(/\t/, $title);
	while(<$f2>) {
		chomp;
		my @a = split(/\t/, $_);
		my $tid = $a[0];
		for(my $i=1; $i<@a; $i++) {
			my $exp = $a[$i];
			$exp =~ s/\..*//;
			$exp++;

			my $sample = $t[$i];
			$iso_exp{$sample}{$tid} = $exp;
		}
	}
	$f2->close;

	# main
	foreach my $gid (sort keys %gid_tid) 
	{
		my @t_member = @{$gid_tid{$gid}}; # this array is not uniq
		my %t_mem; foreach my $t (@t_member) { $t_mem{$t} = 1; } 
		next if (scalar(keys(%t_mem)) == 1);

		# print $gid."\t";
		# print join("\t", keys(%t_mem)),"\n";
		my %checked;
		foreach my $t1 (sort keys %t_mem) 
		{
			foreach my $t2 (sort keys %t_mem) 
			{
				next if $t1 eq $t2;
				next unless defined $pair_tid{$t1."#".$t2};
				next if defined $checked{$t1."#".$t2};
				$checked{$t1."#".$t2} = 1;
				$checked{$t2."#".$t1} = 1;

				# check if the t1 t2 show sig change among different conditions or tissues
				my $sig = 0;
				my %ss; # hash of comparison of two conditions/tissues
				foreach my $s (sort keys %sample_cmb) {
					$sig++ and $ss{$s} = 1 if (defined $tid_sig{$t1}{$s} and  $tid_sig{$t1}{$s} eq 'yes');
					$sig++ and $ss{$s} = 1 if (defined $tid_sig{$t2}{$s} and  $tid_sig{$t2}{$s} eq 'yes');
				}
				# next if $sig == 0;

				# check each sample combination
				foreach my $s (sort keys %sample_cmb) {

				#foreach my $s (sort keys %ss) {
					my @sc = split(/#/, $s);
					my ($s1, $s2) = ($sc[0], $sc[1]);

					print "[ERR]undef exp $s1 $t1\n" and exit unless defined $iso_exp{$s1}{$t1};
					print "[ERR]undef exp $s1 $t2\n" and exit unless defined $iso_exp{$s1}{$t2};
					print "[ERR]undef exp $s2 $t1\n" and exit unless defined $iso_exp{$s2}{$t1};
					print "[ERR]undef exp $s2 $t2\n" and exit unless defined $iso_exp{$s2}{$t2};


					# s means different condition, tissues
					# t means different isoforms
					my $s1t1 = $iso_exp{$s1}{$t1};
					my $s1t2 = $iso_exp{$s1}{$t2};
					my $s2t1 = $iso_exp{$s2}{$t1};
					my $s2t2 = $iso_exp{$s2}{$t2};

					# Fisher Exact Test, two side
					my $p_FET = 1;	# set default P value

					# FET
					#         word2   ~word2
					#  word1    n11      n12 | n1p
					# ~word1    n21      n22 | n2p
					#           --------------
					#           np1      np2   npp
					# my $p1_ct = 3;
					# my $p2_ct = 158;
					# my $f1a_p1_ct = 612;
					# my $f1a_p2_ct = 21;

					# code for FET (left and right is same)
					#my $npp = $p1_ct + $p2_ct + $f1a_p1_ct + $f1a_p2_ct;
					#my $n1p = $p1_ct + $f1a_p1_ct;
					#my $np1 = $p1_ct + $p2_ct;
					#my $n11 = $p1_ct;
					my $npp = $s1t1 + $s1t2 + $s2t1 + $s2t2; 
					my $n1p = $s1t1 + $s2t1; 
					my $np1 = $s1t1 + $s1t2; 
					my $n11 = $s1t1;
					$p_FET = calculateStatistic( n11=>$n11, n1p=>$n1p, np1=>$np1, npp=>$npp) if $npp > 0;

					# output result	
					print "$gid\t$s1\t$s2\t$t1\t$t2\t$s1t1\t$s1t2\t$s2t1\t$s2t2\t$p_FET\n";
				}
			}
		}
	}
}

=head2
 as_event_specific: find specific as event for any condition/tissue
=cut
sub as_event_specific
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t specific [options] AS_gtf sample1_gtf sample2_gtf ... sampleN_gtf

* the AS_gtf is generated by AStalavista
* the sample gtf is assembled by cufflinks

';
	print $usage and exit unless @$files < 2;
	my $as_gtf = shift @$files;
	die "[ERR]file not exist $as_gtf\n";
	foreach my $f (@$files) { die "[ERR]file not exist $f\n"; }

	# put AS event to hash
	# key: splice_chain, AS type, 
	# value: 
	my $fh1 = IO::File->new($as_gtf) || die $!;
	while(<$fh1>)
	{

	}
	$fh1->close;

	# check if the splicing site is specific (need fast method)
	foreach my $f (@$files ) {

	}
	
}

=head2
 as_event_analyze: identify and classify as event 
=cut
sub as_event_analyze
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t analyze [options] input_gtf > output_report

';
	print $usage and exit unless defined $$files[0];
	my $input_gtf = $$files[0];
	die "[ERR]file not exist $input_gtf\n" unless -s $input_gtf;
	die "[ERR]file suffix $input_gtf\n" unless $input_gtf =~ m/\.gtf/;

	my $cutoff = 50;
	$cutoff = $$options{'c'} if (defined $$options{'c'} && $$options{'c'} > 0);

	# sort and asta
	my $sort_gtf = $input_gtf; $sort_gtf =~ s/\.gtf$/_sort\.gtf/;
	my $asi_gtf  = $input_gtf; $asi_gtf  =~ s/\.gtf$/_ASI\.gtf/;
	my $asi_gtf_gz = $asi_gtf.".gz";

	my $cmd1 = "astalavista -t sortGTF -i $input_gtf -o $sort_gtf";
	my $cmd2 = "astalavista -t asta -i $sort_gtf -o $asi_gtf_gz";
	my $cmd3 = "gunzip $asi_gtf_gz";

	run_cmd($cmd1, 1);
	run_cmd($cmd2, 1);
	run_cmd($cmd3, 1);
	#unlink($sort_gtf);

	# get stat of AS event
	#my $as_stat_perl = $FindBin::RealBin."/AS_event_stat.pl";
	#my $cmd4 = "perl $as_stat_perl $asi_gtf $cutoff";
	#run_cmd($cmd4, 1);
	as_event_stat($asi_gtf, 1);		
}

=head2
 as_event_stat: generate statistics info for as result generated by astalavista
 input: input_gtf, cutoff
 ouput: no. of 4 main AS events, and other AS events
=cut
sub as_event_stat
{
	my ($input_gtf, $cutoff) = @_;

	my $usage = qq'
USAGE: perl $0 -t stat input_gtf(output of ASTA) cutoff_display_AS_event

* output four main type number, and other number
  the other number should be total number - 4 main type number

';
	print $usage and exit unless -s $input_gtf;

	$cutoff = 50 unless $cutoff;

	# init number 
	my ($AA, $AD, $ES, $IR, $other, $total) = (0,0,0,0,0,0);

	$total = `wc -l $input_gtf`; chomp($total); $total =~ s/^\s+//; $total =~ s/\s+.*//ig;
	die "[ERR]input gtf file\n" unless $total > 0;

	# init hash
	my %type_stat;
	my %trans_type_stat;
	my %struc_type_count;	# key: structure, AS type; value: count

	# main
	my $fh = IO::File->new($input_gtf) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		die "Error in line $_\n" if scalar @a < 9;
		if ($a[8] =~ m/structure "(\S+)";/)
		{
			my $struc = $1;
			my @type = as_type($struc);
			#print join(", ", @type),"\t$struc\n";
		
			my $type_all = join(",", @type);
			if ( defined $struc_type_count{$struc}{$type_all} )
			{
				$struc_type_count{$struc}{$type_all}++;
			}
			else
			{
				$struc_type_count{$struc}{$type_all} = 1
			}
	
			my %type;
			my @uniq_type;
			foreach my $t (@type) {
				if ( defined $type_stat{$t} ) { $type_stat{$t}++; } 
				else { $type_stat{$t} = 1; }

				unless ( defined $type{$t} ) {
					$type{$t} = 1;
					push(@uniq_type, $t);
				}
			}

			@uniq_type = sort @uniq_type;
			my $uniq_type = join(", ", @uniq_type);

			if ( defined $trans_type_stat{$uniq_type} )
			{
				$trans_type_stat{$uniq_type}++;
			}
			else
			{
				$trans_type_stat{$uniq_type} = 1;
			}
		}
		else
		{
			die "Error in line $_\n";
		}
	}
	$fh->close;

	foreach my $t (sort keys %type_stat) 
	{
		#print $t."\t".$type_stat{$t}."\n";
	}

	#print "-----------------\n";

	#foreach my $t (sort keys %trans_type_stat)
	#{
	#	print $t."\t".$trans_type_stat{$t}."\n";
	#}

	# the others will be split to different events
	my ($AAx, $ADx, $ESx, $IRx) = (0,0,0,0);

	foreach my $stru (sort keys %struc_type_count)
	{
		foreach my $t ( sort keys %{$struc_type_count{$stru}})
		{
			if ($t eq 'AA' && $stru eq "1-,2-" )	{ $AA = $struc_type_count{$stru}{$t}; }
			if ($t eq 'AD' && $stru eq "1^,2^" )	{ $AD = $struc_type_count{$stru}{$t}; }
			if ($t eq 'ES' && $stru eq "0,1-2^" )	{ $ES = $struc_type_count{$stru}{$t}; }
			if ($t eq 'IR' && $stru eq "0,1^2-" )	{ $IR = $struc_type_count{$stru}{$t}; }

			my @m = split(/,/, $t);
			foreach my $m (@m) {
				if      ($m eq 'AA') { $AAx+= $struc_type_count{$stru}{$t}; }
				elsif   ($m eq 'AD') { $ADx+= $struc_type_count{$stru}{$t}; }
				elsif   ($m eq 'ES') { $ESx+= $struc_type_count{$stru}{$t}; }
				elsif   ($m eq 'IR') { $IRx+= $struc_type_count{$stru}{$t}; }
				else    { die "[ERR]Event $t -> $m\n"; }
			}

			if ($struc_type_count{$stru}{$t} >= $cutoff)
			{
				print $stru,"=>",$t,"=>",$struc_type_count{$stru}{$t},"\n";
			}
		}
	} 
	$other = $total - $AA - $AD - $ES - $IR;
	my $totalx = $AAx + $ADx + $ESx + $IRx;
	print "$input_gtf\tAA:$AA\tAD:$AD\tES:$ES\tIR:$IR\tOther:$other\tTotal:$total\n";
	print "$input_gtf\tAAx:$AAx\tADx:$ADx\tESx:$ESx\tIRx:$IRx\tTotal:$totalx\n";
}

# kentnf: subroutine
sub as_type
{
	my $struc = shift;

	my @type; # array for type, the members should be AD, AA, ES, IR
	# construct the hash for each point (splicing site):
	# key: number of splicing site
	# value of transcript {ts} : 1(left), 2(right)
	# value of as code    {as} : as code (^ ~)
	# could get stat_change_num, source_stat, sink_stat from the hash
	# the number of changed site, and the as code (^ -) could be traced 
	my %site;
	my @a = split(/,/, $struc);
	die "[Error]Struc: $struc\n" unless scalar(@a) == 2;	
	my @s1 = $a[0] =~ m/(\d+)/g; 
	my @s2 = $a[1] =~ m/(\d+)/g;
	my @as1 = $a[0] =~ m/(\^|-)/g;
	my @as2 = $a[1] =~ m/(\^|-)/g;
	if ($a[0] eq '0' ) { @s1 = (); @as1 = (); }
	die "[Error]left AS code do not match site number: $a[0]\n" unless scalar(@s1) == scalar(@as1);
	die "[Error]left AS code do not match site number: $a[1]\n" unless scalar(@s2) == scalar(@as2);
	for (my $i=0; $i<@s1; $i++) {
		$site{$s1[$i]}{'ts'} = 1;
		$site{$s1[$i]}{'as'} = $as1[$i];
	}
	for (my $i=0; $i<@s2; $i++) {
		$site{$s2[$i]}{'ts'} = 2;
		$site{$s2[$i]}{'as'} = $as2[$i];
	}
	my $total_site_num = scalar(@s1) + scalar(@s2);
	die "[Error]num of splicing site is not even number: $total_site_num\n$struc\n" unless ($total_site_num % 2 == 0);
	die "[Error]splicing site not start from 1: $struc\n" unless defined $site{"1"}{'ts'};

	for (my $i=1; $i<=$total_site_num; $i = $i+2) 
	{
		my $j = $i+1;	# i,j should be a pair of splicing events
		die "[Error]undef splicing site num: $i\n$struc\n" unless defined $site{$i}{'ts'};
		die "[Error]undef splicing site code:$i\n$struc\n" unless defined $site{$i}{'as'};
		die "[Error]undef splicing site num: $j\n$struc\n" unless defined $site{$j}{'ts'};
		die "[Error]undef splicing site code:$j\n$struc\n" unless defined $site{$j}{'as'};
		my $i_tsstat = $site{$i}{'ts'};
		my $i_ascode = $site{$i}{'as'};
		my $j_tsstat = $site{$j}{'ts'};
		my $j_ascode = $site{$j}{'as'};

		if ( $i_tsstat ne $j_tsstat ) 	# for event AD or AA
		{
			die "[Error]Diff AS code: $i $j\t$i_ascode $j_ascode\n$struc\n" if $i_ascode ne $j_ascode;
			#if 	($i_ascode eq '^') { push(@type, 'AD') if scalar @type == 0 || $type[scalar(@type)-1] ne 'AD'; }
			#elsif 	($i_ascode eq '-') { push(@type, 'AA') if scalar @type == 0 || $type[scalar(@type)-1] ne 'AA'; }
			if     ($i_ascode eq '^') { push(@type, 'AD'); }
			elsif  ($i_ascode eq '-') { push(@type, 'AA'); }

		}
		else				# for event ES and IR
		{
			die "[Error]Same AS code: $i $j\t$i_ascode $j_ascode\n$struc\n" if $i_ascode eq $j_ascode;

			# add code for detect another type of AS mu
			#if	($i_ascode eq '^') { push(@type, 'IR') if scalar @type == 0 || $type[scalar(@type)-1] ne 'IR'; }
			#elsif   ($i_ascode eq '-') { push(@type, 'ES') if scalar @type == 0 || $type[scalar(@type)-1] ne 'ES'; }
			if	($i_ascode eq '^') { push(@type, 'IR'); }
			elsif   ($i_ascode eq '-') { push(@type, 'ES'); }
		}
	}
	return @type;
}

=head1
 AS_splice_site_stat.pl -- get number of different splice site
 Yi Zheng
 201406:init 
=cut
sub as_splice_site
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: perl $0 -t spliceSite -g genome_sequence input_gtf 
	
';
	print $usage and exit unless defined $$files[0];
	my $input_gtf = $$files[0];
	die "[ERR]file not exist\n" unless -s $input_gtf;
	
	print $usage and exit unless defined $$options{'g'};
	my $genome = $$options{'g'};
	die "[ERR]file not exist\n" unless -s $genome;

	# put the transcript exon info to hash
	# key: tid; value: exon1_start exon1_end exon2_start exon2_end ....
	my %tid_info;
	my $fh = IO::File->new($input_gtf) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		# SL2.40ch00	Cufflinks	exon	3300	3326	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1"; oId "CUFF.1.1"; class_code "u"; tss_id "TSS1"; lincRNA "1";
		die "[ERR]col num $_\n" if scalar @a < 9;

		if ($a[2] eq 'exon')
		{
			if ($a[8] =~ m/transcript_id "(\S+)"/) 
			{
				my $tid = $1;
				if (defined $tid_info{$tid}{'exon'})
				{
					$tid_info{$tid}{exon}.="\t".$a[3]."\t".$a[4];
				}
				else
				{
					$tid_info{$tid}{exon} = $a[3]."\t".$a[4];
				}

				if (defined $tid_info{$tid}{'chr'}) 
				{
					die "Error, chr for $tid\n" if $tid_info{$tid}{'chr'} ne $a[0];
				}
				else
				{
					$tid_info{$tid}{'chr'} = $a[0];
				}

				if (defined $tid_info{$tid}{'strand'})
				{
					die "Error, strand for $tid\n" if $tid_info{$tid}{'strand'} ne $a[6];
				}
				else
				{
					$tid_info{$tid}{'strand'} = $a[6];
				}
			}
			else
			{
				die "Error in line $_\n";
			}
		}
	}
	$fh->close;

	# convert tid exon hash to split site hash, 
	# key: Chr exon_end 
	# value: exon start1, ...... 
	# key: chr exon_end, exon_start
	# value: strand
	my %site;
	my %site_strand;

	foreach my $tid (sort keys %tid_info)
	{
		my $chr = $tid_info{$tid}{'chr'};
		my $strand = $tid_info{$tid}{'strand'};
		my @exon = split(/\t/, $tid_info{$tid}{'exon'});
		@exon = sort {$a<=>$b} @exon;
		next if scalar @exon == 2;
		die "Error in exon num for $tid\n$tid_info{$tid}{'exon'}\n" if scalar @exon % 2 == 1;

		# check the exon order
		my $pre_exon = $exon[0];
		for(my $i=1; $i<@exon; $i++) {
			if ( $exon[$i] < $pre_exon ) {
				die "Error in exon order for $tid\n$tid_info{$tid}{'exon'}\n";
			} {
				$pre_exon = $exon[$i];
				next;
			}
		}

		shift @exon;
		pop @exon;

		# coreate hash for site
		for(my $i=0; $i<@exon; $i=$i+2) {
			my $start = $exon[$i];
			my $end = $exon[$i+1];
		
			if (defined $site_strand{$chr."\t".$start."\t".$end} ) 
			{
				if ( $site_strand{$chr."\t".$start."\t".$end} eq $strand ) {
	
				} else {
					#$site{$chr."\t".$start."\t".$end} = ".";
				}
				next;
			}
			else
			{
				$site_strand{$chr."\t".$start."\t".$end} = $strand;
			}

			if ( defined $site{$chr."\t".$start} )
			{
				$site{$chr."\t".$start}.= "\t".$end;
			}
			else
			{
				$site{$chr."\t".$start}	= $end;
			}
		}
	}

	# get splicing site, then put it to junc_stat hash
	# key: GT-AG
	# value: number of this site
	my %junc_stat;

	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$genome);
	while(my $inseq = $in->next_seq)
	{
		my $chr = $inseq->id;
		my $seq = $inseq->seq;
		my $len = $inseq->length;

		# below start is the exon end
		# below end is the exon start
		for(1 .. $len) 
		{
			my $start = $_;

			if ( defined $site{$chr."\t".$start} )
			{
				my @end = split(/\t/, $site{$chr."\t".$start});
				foreach my $end (@end) 
				{
					my $strand = $site_strand{$chr."\t".$start."\t".$end};

					my $start_base = uc(substr($seq, $start, 2));
					my $end_base = uc(substr($seq, $end-3, 2));

					my $junction;
					if ($strand eq "+") 
					{
						$junction = $start_base."-".$end_base;
					}
					elsif ($strand eq "-") 
					{
						$junction = $start_base."-".$end_base;
						my $rev_junc = reverse($junction);
						$rev_junc =~ tr/ACGTacgt-/TGCAtgca-/;
						$junction = $rev_junc;
					}
					else
					{
						next;
					}

					if (defined $junc_stat{$junction})
					{
						$junc_stat{$junction}++;
					}
					else
					{
						$junc_stat{$junction} = 1;
					}
				}
			}
		}
	}

	# report junc start
	foreach my $j (sort keys %junc_stat)
	{
		print $j."\t".$junc_stat{$j}."\n";
	}
}

=head2
 run_cmd : run command
=cut
sub run_cmd
{
        my ($cmd, $debug) = @_;
        print "[DEBUG]".$cmd."\n" if $debug;
        system($cmd) && die "[Error][CMD]$cmd\n";
}

=head2
 usage: print usage information
=cut
sub usage
{
        my $usage = qq'
USAGE: $0 -t tool 

	analyze		analyze AS events number
	stat		check the type of AS event
	spsites		count the number of splice site
	isodiff		statistics analysis for AS event to find change isoform within gene

';
        print $usage; exit;
}

