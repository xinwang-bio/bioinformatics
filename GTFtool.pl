#!/usr/bin/perl

=head1
 GTFtool -- parse GTF file
=cut
use strict;
use warnings;
use IO::File;
use Getopt::Std;
use FindBin;

my $version = 0.1;
my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }

if      ($options{'t'} eq 'stats')		{ gtf_stats(\%options, \@ARGV); }
elsif   ($options{'t'} eq 'convert')	{ gtf_convert(\%options, \@ARGV); }
elsif	($options{'t'} eq 'extract')	{ gtf_extract(\%options, \@ARGV); }
elsif	($options{'t'} eq 'gffread')    { gtf_gffread(); }
elsif	($options{'t'} eq 'sortgff')    { gtf_sortgff(\%options, \@ARGV); }
elsif   ($options{'t'} eq 'togpd')    	{ gtf_to_gpd(\%options, \@ARGV); }
elsif   ($options{'t'} eq 'tojunc')		{ gtf_to_junc(\%options, \@ARGV); }
elsif	($options{'t'} eq 'compare')	{ gtf_compare(\%options, \@ARGV); }
elsif	($options{'t'} eq 'findptc')	{ gtf_findptc(\%options, \@ARGV); }
else    { usage($version); }

#############################################
# kentnf: subroutine						#
#############################################

=head2
 gtf_sortgff -- sort GFF by order
=cut
sub gtf_sortgff
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGe: $0 -t sortgff  input_gff  order_gID > output_gff

this program only sort GFF with correct order of single gene structure
eg. the gene are disordered in chr level, but the mRNA, utr, exon, cds 
are well ordered in gene level

';
	print $usage and exit unless @$files == 2;
	my ($input_gff, $order_gid) = @$files;	
	
	# load gene to hash
	my %gid_gff; # key: gid, value: gff of the gid
	my $gid;
	open(FH, $input_gff) || die $!;
	while(<FH>)	
	{
		chomp;
		my @a = split(/\t/, $_);
	
		if ($a[2] eq 'gene') {
			my @b = split(/;/, $a[8]);
			my $ingid;
			foreach my $b (@b) {
				if ($b =~ m/ID=(\S+)/) {
					$ingid = $1;	
				}
			}
			
			if (defined $ingid) {
				$gid = $ingid;
			} else {
				die "[ERR]no gid $_\n";
			}
			$gid_gff{$gid} = $_."\n";
		} else {
			$gid_gff{$gid}.= $_."\n";
		}
	}
	close(FH);

	# output GFF with correct order
	open(IN, $order_gid) || die $!;
	while(<IN>) {
		chomp;
		my $gid = $_;
		die "[ERR]undef gene id: $gid\n" unless defined $gid_gff{$gid};
		print $gid_gff{$gid};
	}
	close(IN);
}

=head2
 gtf_to_bed: convert GTF to bed12 format
=cut
sub gtf_to_bed {
	
}

=head2
 gtf_findptc -- find ptc code
 1. cluster input GTF and annotation_GTF
 2. remove clusters with >=2 Solyc ID
 3. locate the start and end CDS position of Solyc in each cluster
 4. translate isoform of each cluster
 5. compare the region of best translation with region of Solyc
 Out: 
 TID, START, END, CDS-START, CDS-END, PTC?
 This is init program, will improve it later
=cut 
sub gtf_findptc
{
	my ($options, $files) = @_;	
	my $usage = qq'
USAGE: $0 -t findptc reference annotation_GTF input_GTF

* the annotation GTF file must contain CDS feature;

';		
	print $usage and exit if @$files != 3;
	my ($refseq, $anno_gtf, $input_gtf) = @$files;	


	# locate the start and end of Solyc CDS
	# key: mRNAID, value: [start, end]
	my %cds_region;		
}

=head2
 gtf_compare -- compare two GTF files
=cut

sub gtf_compare
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t compare input_GTF1 input_GTF2 output_prefix

* the output includes one shared GTF, two specific GTF 

';
	print $usage and exit unless @$files >= 3;
	my ($gtf1, $gtf2, $output_prefix) = @$files;

	# set bin folder	
	my $bedtools = $FindBin::RealBin."/bin/bedtools";
	my $gtftools = $FindBin::RealBin."/GTFtool.pl";

	# convert GTF to bed;
	# combine 2 bed files, sort and cluster 
	my $cmd_A2bed = "perl $gtftools -t convert -o A -f bed $gtf1";
	my $cmd_B2bed = "perl $gtftools -t convert -o B -f bed $gtf2";
	run_cmd($cmd_A2bed);
	run_cmd($cmd_B2bed);
	run_cmd("cat A.bed B.bed | $bedtools sort -i stdin | $bedtools cluster -i stdin > C.bed");

	# load cluster to hash
	my %cls; # key: CLS ID, value: array of tid
	my $f1 = IO::File->new("C.bed") || die $!;
	while(<$f1>) {
		chomp;
		my @a = split(/\t/, $_);
		push(@{$cls{$a[6]}}, $a[3]);
	}
	$f1->close;

	# load tid information
	run_cmd("cat $gtf1 $gtf2 > temp.gtf");
	my %trans_info = parse_gtf("temp.gtf");
	
	# put compare result to hash
	my %class; #key tid; value, 0, share 1 specific 
	my %single_exon_iso; # key: ID

	# remove duplicate transcript according to intron
	foreach my $cid (sort {$a<=>$b} keys %cls) {
		my @tids = @{$cls{$cid}};
		next if scalar(@tids) == 1; # loci without AS		

		# definition of exon boundaries for all transcripts in cluster
		my %p;
		foreach my $tid (@tids) {
			my @exon = split("\t", $trans_info{$tid}{'exon'});
			@exon = sort {$a<=>$b} @exon;
			if (@exon == 2) {
				$single_exon_iso{$tid} = 1;
				next;
			}
			die "[ERR]exon num $tid\n" unless (scalar(@exon) % 2 == 0);
			shift @exon; pop @exon;
			foreach my $e (@exon) {
				$p{$e} = 1;
			}
		}

		# give the order of p
		my $n = 0;
		foreach my $p (sort {$a<=>$b} keys %p) {
			$n++;
			$p{$p} = "p$n";
		}

		# put exon boundaries to each transcripts
		my %tid_p;
		foreach my $tid (@tids) {
			my @exon = split("\t", $trans_info{$tid}{'exon'});
			@exon = sort {$a<=>$b} @exon;
			next if @exon == 2;
			die "[ERR]exon num $tid\n" unless (scalar(@exon) % 2 == 0);
			# construct position str for each transcripts
			my $pstr = '';
			for (my $i=0; $i<@exon; $i=$i+2) 
			{
				my $es = $exon[$i];
				my $en = $exon[$i+1];
				foreach my $p (sort {$a<=>$b} keys %p) {
					next if ($p < $es || $p > $en);
					my $pname = $p{$p};
					if ($p == $es) {
						$pstr.= "$pname+,";
					} elsif ($p == $en) {
						$pstr.= "$pname-,";
					} else {
						$pstr.= "$pname-,$pname+,";
					}

				}
			}
			$pstr =~ s/,$//;
			$tid_p{$tid} = $pstr;
		}

		# sort hash %tid_p by value length
		my %pstrlen_tid;
		foreach my $tid (sort keys %tid_p) {
			my $pstr = $tid_p{$tid};
			my $len = length($pstr);
			push(@{$pstrlen_tid{$len}}, $tid);
		}

		my @tid_order;
		foreach my $len (sort {$b <=> $a} keys %pstrlen_tid) {
			my @tids = @{$pstrlen_tid{$len}};
			foreach my $t (@tids) { push(@tid_order, $t); }	
		}

		# remove duplicate according to pstr
		my %uniq; # key: tid, value: pstr
		foreach my $tid ( @tid_order )
		{
			my $pstr = $tid_p{$tid};
			
			if (scalar(keys(%uniq)) == 0) 
			{
				$uniq{$tid} = $pstr;
			}
			else
			{
				# check if there is any duplicate GTF in file
				my $share = 0;
				foreach my $u_tid (sort keys %uniq) {
					my $u_pstr = $uniq{$u_tid};
					#print "$u_tid, $tid\n";
					my ($type, $cmb_pstr) = compare_pstr($u_pstr, $pstr, 0);
					if ($type == 0) {
						my $c_tid = $u_tid."#".$tid;
						$uniq{$c_tid} = $cmb_pstr;	# use combined one
						delete $uniq{$tid};			# delete 
						delete $uniq{$u_tid};		# delete
						$share = 1;					# determine the share stat
						last;						# exit if find the shared seq
					}
				}

				# put the tid to uniq if it not shared with the GTF in uniq set
				$uniq{$tid} = $pstr if $share == 0; 
			}	

		}

		# convert uniq to shared or specific
		my $exit = 0;
		my %sub_class;
		foreach my $t (sort keys %uniq) {
			#$exit = 1 if $t =~ m/HF_MG\.152\.1/;
			if ($t =~ m/#/) {
				my @m = split(/#/, $t);
				foreach my $m (@m) { $class{$m} = 0; $sub_class{$m} = 0; }
			} else {
				$class{$t} = '1';
				$sub_class{$t} = 1;
			}
		}
		#if ($exit == 1 ) { foreach my $m (sort keys %sub_class) { print "$m\t$class{$m}\n"; } exit; }
	}

	# output files;
	my $out_share = IO::File->new(">".$output_prefix.".share.gtf") || die $!;
	my $out_spec1 = IO::File->new(">".$output_prefix.".spec1.gtf") || die $!;
	my $out_spec2 = IO::File->new(">".$output_prefix.".spec2.gtf") || die $!;

	my $fh1 = IO::File->new($gtf1) || die $!;
    while(<$fh1>)
    {
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
			
		# analyzing the attribute 
		my %attr_value;
		my @b = split(/; /, $a[8]);
		foreach my $b (@b) {
			$b =~ s/;$//;
			my @c = split(/ /, $b, 2);
			die "[ERR]attr $b in $a[8]\n" unless @c == 2;
			$c[1] =~ s/"//ig;
			$attr_value{$c[0]} = $c[1];
		}

		die "[ERR]No transcript id $_\n" unless defined $attr_value{'transcript_id'};
		die "[ERR]No gene id $_\n" unless defined $attr_value{'gene_id'};
		my ($tid, $gid) = ($attr_value{'transcript_id'}, $attr_value{'gene_id'});

		if (defined $class{$tid} ) {
			print $out_spec1 $_."\n" if $class{$tid} == 1;
			print $out_share $_."\n" if $class{$tid} == 0;
		}
	}
	$fh1->close;

	my $fh2 = IO::File->new($gtf2) || die $!;	
	while(<$fh2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);

		# analyzing the attribute
		my %attr_value;
		my @b = split(/; /, $a[8]);
		foreach my $b (@b) {
			$b =~ s/;$//;
			my @c = split(/ /, $b, 2);
			die "[ERR]attr $b in $a[8]\n" unless @c == 2;
			$c[1] =~ s/"//ig;
			$attr_value{$c[0]} = $c[1];
		}

		die "[ERR]No transcript id\n" unless defined $attr_value{'transcript_id'};
		die "[ERR]No gene id\n" unless defined $attr_value{'gene_id'};
		my ($tid, $gid) = ($attr_value{'transcript_id'}, $attr_value{'gene_id'});

		if (defined $class{$tid} ) {
			print $out_spec2 $_."\n" if $class{$tid} == 1;
			print $out_share $_."\n" if $class{$tid} == 0;
		}
	}
	$fh2->close;

	$out_share->close;
	$out_spec1->close;
	$out_spec2->close;

	################################
	# compare single-exon isoforms #
	################################
	
	# 1. extract single exon iso from bad
	my @bed_files = qw/A.bed B.bed/;
	foreach my $bed (@bed_files) {
		my $output = $bed;
		$output =~ s/\.bed$/_singleISO\.bed/;
		my $out = IO::File->new(">".$output) || die $!;
		my $in = IO::File->new($bed) || die $!;
		while(<$in>) {
			chomp;
			my @a = split(/\t/, $_);
			print $out $_."\n" if defined $single_exon_iso{$a[3]};
		}
		$in->close;
		$out->close;
	}

	# compare them using bed tools
	my $cmd_S1 = "$bedtools intersect -v  -a A_singleISO.bed -b B_singleISO.bed > A_specSingleISO.bed";
	my $cmd_S2 = "$bedtools intersect -wa -a A_singleISO.bed -b B_singleISO.bed > A_shareSingleISO.bed";
	my $cmd_S3 = "$bedtools intersect -v  -a B_singleISO.bed -b A_singleISO.bed > B_specSingleISO.bed"; 
	my $cmd_S4 = "$bedtools intersect -wb -a A_singleISO.bed -b B_singleISO.bed > B_shareSingleISO.bed";
	run_cmd($cmd_S1);
	run_cmd($cmd_S2);
	run_cmd($cmd_S3);
	run_cmd($cmd_S4);
}

sub compare_pstr 
{
	my ($pstr1, $pstr2, $debug) = @_;
	print "$pstr1, $pstr2\n" if $debug == 1;

	my @p1 = split(/,/, $pstr1);
	my @p2 = split(/,/, $pstr2);
	#print "$p1[0]\t$p2[0]\n";

	# check which p[0] is left
	my $left_p1 = $p1[0]; $left_p1 =~ s/p//; $left_p1 =~ s/\+//; $left_p1 =~ s/-//;
	my $left_p2 = $p2[0]; $left_p2 =~ s/p//; $left_p2 =~ s/\+//; $left_p2 =~ s/-//;
	print "$left_p1\t$left_p2\n" if $debug == 1;

	# determine start from which points
	my $left = $left_p2;
	my @px = @p2;
	my @py = @p1;
	if ($left_p1 <= $left_p2) { $left = $left_p1; @px = @p1; @py = @p2; }
	print "$px[0]\t$py[0]\n" if $debug == 1;

	# start from px, then find the 1st match of py;
	my @cmb_p = ();
	my $type = 0;
	my $py_start_count = 0; # switch for count py
	my $match = 0;
	while(@px > 0) {
		my $x = shift @px;
		#unless (defined $py[0]) {  print $x."\n"; foreach my $py (@py) { print $py."\n";} exit;}
		if ($py_start_count == 0) {
			if ($x eq $py[0]) { $py_start_count = 1; $match = 1;} # 1st match point
		}
		if ($py_start_count == 1 && scalar @py > 0) {
			my $y = shift @py;
			if ( $x ne $y )	{
				$type = 1;
				last;
			}
		} 
		push(@cmb_p, $x);
	}
	$type = 1 if $match == 0;

	# put remaining py to combine p
	while(@py > 0){
		my $y = shift @py;
		push(@cmb_p, $y);
	}
	my $cmb_pstr = join(",", @cmb_p);

	print $type."\n" if $debug == 1;
	return ($type , $cmb_pstr);
}

=head2
 gtf_to_junc -- extract junction from GTF
=cut
sub gtf_to_junc
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t tojunc input_GTF > ouput_file

# format of output
Chr#juncStart#juncEnd

* using the ouput information could get the stat if intron length

';
	print $usage and exit unless defined $$files[0];
        
	# load junctions from gene annotation file
	my $gtf_file = $$files[0];
	my %gtf_junction;
	my %trans_info = parse_gtf($gtf_file);
                
	foreach my $tid (sort keys %trans_info) {
		my $chr = $trans_info{$tid}{'chr'};
		my @exon = split("\t", $trans_info{$tid}{'exon'});
		@exon = sort {$a<=>$b} @exon;
		next if @exon == 2;
		die "[ERR]exon num $tid\n" unless (scalar(@exon) % 2 == 0);
		shift @exon; pop @exon;
                        
		for(my $i=0; $i<@exon; $i=$i+2) {
			my $jstart = $exon[$i]+1;
			my $jend = $exon[$i+1]-1;
			$gtf_junction{$chr."#".$jstart."#".$jend} = 1;
		}
	}
	#print scalar(keys(%gtf_junction)), " of junctions has been loaded from $$options{'g'}\n"; #exit;
        
	foreach my $jk (sort keys %gtf_junction) {
		print $jk."\n";
	}
}

=head2
  gtf_to_gpd -- convert gtf to gpd
=cut
sub gtf_to_gpd
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t togpd -i [gpd/gtf] input_file > ouput_file

example1 $0 -t togpd -i gpd gpd > gtf
example1 $0 -t togpd -i gtf gtf > gpd

';
	print $usage and exit unless defined $$files[0];
	my $input_file = $$files[0];
	die "[ERR]file not exist\n" unless -s $input_file;

	print $usage and exit unless defined $$options{'i'};
	
	if ($$options{'i'} eq 'gtf') 
	{
		my %trans_info = parse_gtf($input_file);
		foreach my $tid (sort keys %trans_info)
		{
			my $chr		= $trans_info{$tid}{'chr'};
			my $gid		= $trans_info{$tid}{'gid'};
			my $strand 	= $trans_info{$tid}{'strand'};
			my @exon	= split("\t",$trans_info{$tid}{'exon'});
			@exon = sort {$a<=>$b} @exon;
			die "[ERR]exon num\n" unless ((scalar @exon) % 2 == 0);
			my $exon_num = (scalar @exon) / 2;

			my @cds;
			if (defined $trans_info{$tid}{'cds'}) {
				@cds = split("\t",$trans_info{$tid}{'cds'});
			} else {
				@cds = @exon; 
			}

			@cds = sort {$a<=>$b} @cds;
			die "[ERR]CDS num\n" unless ((scalar @cds) % 2 == 0);

			my $start_e = $exon[0];
			my $end_e = $exon[scalar(@exon)-1];
			my $start_c = $cds[0];
			my $end_c = $cds[scalar(@cds)-1];
			my $exon_start = '';
			my $exon_end = '';

			for(my $i=0; $i<@exon; $i=$i+2)
			{
				$exon_start.=$exon[$i].",";
				$exon_end.=$exon[$i+1].",";
			}

			print "$gid\t$tid\t$chr\t$strand\t$start_e\t$end_e\t$start_c\t$end_c\t$exon_num\t$exon_start\t$exon_end\n";
		}
	} 
	elsif ($$options{'i'} eq 'gpd')
	{
		my $fh = IO::File->new($input_file) || die $!;
		while(<$fh>)
		{
			chomp;
			next if $_ =~ m/^#/;
			my @a = split(/\t/, $_);
			die "[ERR]col num: $_\n" if (scalar @a < 11);
			my ($gid, $tid, $chr, $strand, $e_start, $e_end) = ($a[0], $a[1], $a[2], $a[3], $a[9], $a[10]);
			my @m = split(/,/, $e_start);
			my @n = split(/,/, $e_end);
			die "[ERR]exon num\n" unless (scalar(@m) == scalar(@n));
			print "$chr\tGPD\ttranscript\t$a[4]\t$a[5]\t.\t$strand\t.\tgene_id \"$gid\"; transcript_id \"$tid\";\n";
			for(my $i=0; $i<@m; $i++)
			{
				print "$chr\tGPD\texon\t$m[$i]\t$n[$i]\t.\t$strand\t.\tgene_id \"$gid\"; transcript_id \"$tid\";\n";
			}	
		}
		$fh->close;
	}
	else
	{
		die "[ERR]parameter i $$options{'i'}\n";
	}
}

=head2
 gtf_extract -- extract gtf information by ID;
=cut
sub gtf_extract
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t extract [options] GTF list/BED

* listID could be gene ID, transcript ID;
* key for gene ID: gene_id "";
* key for transcript ID: transcript_id "";

';
	print $usage and exit unless scalar @$files == 2;
	my ($input_gtf, $list_file) = @$files;
	die "[ERR] file not exit $input_gtf" unless -s $input_gtf;
	die "[ERR] file not exit $list_file" unless -s $list_file;

	my $column = 1;
	$column = 4 if ($list_file =~ m/\.bed$/);

	# get listID from file
	my %list_id;
	my $fh1 = IO::File->new($list_file) || die $!;
	while(<$fh1>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		$list_id{$a[$column-1]} = 1;
	}
	$fh1->close;	

	my ($gid, $tid);

	my $fh2 = IO::File->new($input_gtf) || die $!;
	while(<$fh2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		$gid = ''; $tid = '';
		if ($_ =~ m/gene_id "(\S+)"; /) { $gid = $1; }
		if ($_ =~ m/transcript_id "(\S+)"; /) { $tid = $1; }
		else { die "Error, can not find gid or tid $_\n"; }
		die "[ERR]Dup ID $_\n" if $gid eq $tid;
		print $_."\n" if $gid && defined $list_id{$gid};
		print $_."\n" if $tid && defined $list_id{$tid};
	}
	$fh2->close;
}

=head2
 stat -- generate statistics information
=cut
sub gtf_stats
{
	my ($options, $files) = @_;
	my $subUsage = qq'
USAGE: $0 stats [options]
	-i	input file 
	-x	report exon/exon length statistics info
		(0 or 1, default: 0)
	-r	remove feature by type, and output new GTF file
		s (remove gene with inconsistency strand)		

';
	print $subUsage and exit unless $$options{'i'};
	my ($inFile, $out_prefix);
	$inFile = $$options{'i'};

	my $remove_type;
	$remove_type = $$options{'r'} if defined $$options{'r'};

	my $report_length_stat = 0;
	$report_length_stat = 1 if (defined $$options{'x'} && $$options{'x'} == 1);

	# save input GTF file to hash
	my %trans_info = parse_gtf($inFile);

	my %gene;
	my %trans;
	my %exon;
	my ($chr, $strand);
	my %gene_strand;

	# save exon/intron to hash, the exon/intron is only for multi-exon gene
	my %exon_list;
	my %intron_list;

	my ($trans_num_single_exon, $trans_num_multi_exon) = (0,0);

	foreach my $tid (sort keys %trans_info)
	{
		$trans{$tid} = 1;
		my $gid = $trans_info{$tid}{'gid'};
		defined $gene{$gid} and $gene{$gid}++ or $gene{$gid} = 1;
		$chr = $trans_info{$tid}{'chr'};
		$strand = $trans_info{$tid}{'strand'};
		
		if (defined $gene_strand{$gid}) {
			print "[WARN]inconsistency strand: $gid\n" and $gene_strand{$gid} = "+-" if ($strand ne $gene_strand{$gid});
		} else {
			$gene_strand{$gid} = $strand;
		}

		my @exon = split(/\t/, $trans_info{$tid}{'exon'});
		@exon = sort {$a<=>$b} @exon;
		my $num_exon_p = scalar(@exon);
		die "[ERR]exon num $tid\n" unless $num_exon_p % 2 == 0;
		
		if ($num_exon_p / 2 == 1) {
			$trans_num_single_exon++;
		} else {
			$trans_num_multi_exon++;
		}

		for(my $i=0; $i<scalar(@exon)-1; $i=$i+2)
		{
			my $start = $exon[$i];
			my $end = $exon[$i+1];
			my $key = "$chr\t$start\t$end";
			$exon{$key} = 1;
		}

		# put exon/intron to hash;
		if ($num_exon_p/2 > 1) {
			for(my $i=0; $i<scalar(@exon)-1; $i=$i+2) {
				my $start = $exon[$i]-1; # 0base
				my $end = $exon[$i+1];
				my $key = "$chr\t$start\t$end";
				# output bedformat
				$exon_list{$key} = $end - $start;
			}

			shift @exon; pop @exon;
			for(my $i=0; $i<scalar(@exon)-1; $i=$i+2) {
				my $start = $exon[$i]; 	# exon-end, is 0base of intron start
				my $end = $exon[$i+1]-1;# exon-start, minus 1 is the intron end
				my $key = "$chr\t$start\t$end";
				# output bedformat
				$intron_list{$key} = $end - $start;
			}
		}
	}

	print $inFile,"\tGene:",scalar(keys(%gene)),"\tTrans:",scalar(keys(%trans)),"\tExon:",scalar(keys(%exon)),"\n";
	print "Multiple-exon Trans: $trans_num_multi_exon\tSingle-exon Trans: $trans_num_single_exon\n";

	# report iso num for each gene
	my %gene_iso_stat;
	foreach my $gid (sort keys %gene) {
		my $num = $gene{$gid};
		$gene_iso_stat{$num} and $gene_iso_stat{$num}++ or $gene_iso_stat{$num} = 1;
	}
	print "LocuIsoNum:";
	foreach my $n (sort {$a<=>$b} sort keys %gene_iso_stat) {
		print " *include single exon*" if $n == 1;
		print " $n($gene_iso_stat{$n});";
	}
	print "\n";

	# output exon and intron length statistics info
	if ($report_length_stat == 1) {
	
		my @plist = (0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999);

		print "[Note]Below Percentile only count multi-exon transcripts\n";
		my ($max_exon, $min_exon) = (1,1000);
		my $o1 = IO::File->new(">".$inFile.".exon") || die $!;
		my @array_exon;
		foreach my $exon (sort keys %exon_list) {
			print $o1 $exon."\t".$exon_list{$exon}."\n";
			push(@array_exon, $exon_list{$exon});
			$max_exon = $exon_list{$exon} if $exon_list{$exon} > $max_exon;
			$min_exon = $exon_list{$exon} if $exon_list{$exon} < $min_exon; 
		}
		$o1->close;
		print "Max Exon:$max_exon\tMin Exon:$min_exon\n";
	
		foreach my $p (@plist) {
			my $p100 = $p * 100;
			my $percentile = percentile(\@array_exon, $p);
			print "[PERCENTILE]exon\t$p100\t$percentile\n";
		}

		my ($max_intron, $min_intron) = (1,1000);
		my @array_intron;
		my $o2 = IO::File->new(">".$inFile.".intron") || die $!;
		foreach my $intron (sort keys %intron_list) {
			print $o2 $intron."\t".$intron_list{$intron}."\n";
			push(@array_intron, $intron_list{$intron});
			$max_intron = $intron_list{$intron} if $intron_list{$intron} > $max_intron;
			$min_intron = $intron_list{$intron} if $intron_list{$intron} < $min_intron;
		}
		$o2->close;

		print "Max Intron:$max_intron\tMin Intron:$min_intron\n";
	
		foreach my $p (@plist) {
			my $p100 = $p * 100;
			my $percentile = percentile(\@array_intron, $p);
			print "[PERCENTILE]intron\t$p100\t$percentile\n";
		}
	}

	# output file if remove type works 
	if (defined $remove_type) {

		my $output = $inFile;
		$output =~ s/\.gtf$//ig;
		$output = $output.".filter.gtf";
		my $out = IO::File->new(">".$output) || die $!;
		my $in = IO::File->new($inFile) || die $!;
		while(<$in>)
		{
			chomp;
			next if $_ =~ m/^#/;
			my @a = split(/\t/, $_);
			# next if $a[2] ne 'exon';
				
			# analyzing the attribute 
			my %attr_value;
			my @b = split(/; /, $a[8]);
			foreach my $b (@b) {
				$b =~ s/;$//;
				my @c = split(/ /, $b, 2);
				die "[ERR]attr $b in $a[8]\n" unless @c == 2;
				$c[1] =~ s/"//ig;
				$attr_value{$c[0]} = $c[1];
			}

			my ($tid, $gid) = ($attr_value{'transcript_id'}, $attr_value{'gene_id'});
			
			if ($remove_type =~ m/s/) {
				next if $gene_strand{$gid} eq '+-';
			}
			print $out $_."\n";
		}
		$in->close;
		$out->close;
	}
}

=head2
 convert -- convert GTF to GFF,BED format
=cut
sub gtf_convert
{
	my ($options, $files) = @_;
	my $subUsage = qq'
USAGE $0 -t convert [options] input_gtf

	-f	output format (must be bed, or gff format, default: input.bed)
	-o	output file name (prefix)
	-b	feature type for bed: transcript_id, gene_id, exon (default: transcript_id)
	-a	add chromosome/scaffold seq feature to gtf

';
	print $subUsage and exit unless defined $$files[0] ;
	my $input_gtf = $$files[0];
	die "[ERR]file not exist: $input_gtf\n" unless -s $input_gtf;
	die "[ERR]file format: $input_gtf\n" unless $input_gtf =~ m/\.gtf$/i;
	
	my $output_file = $input_gtf; $output_file =~ s/\.gtf$/\.bed/i;
	$output_file = $$options{'o'} if defined $$options{'o'};

	my $output_format = 'bed';
	$output_format = $$options{'f'} if defined $$options{'f'};
	die "[ERR]file format: $output_format\n" unless ($output_format eq "bed" || $output_format eq "gff");

	$output_file.=".".$output_format;
	die "[ERR]file exist: $output_file\n" if -s $output_file;

	my $out_format = $output_format; # will delete this var later

	my $bed_type = 'transcript_id';
	$bed_type = $$options{'b'} if defined $$options{'b'};
	die "[ERR]bed feature type $bed_type\n" unless $bed_type =~ m/(transcript_id|gene_id|exon)/i;

	my %trans_info = parse_gtf($input_gtf);
	my %gene_info;
	# constract gene info according to trans info
	my ($chr, $gid, $start, $end, $strand);
	foreach my $tid (sort keys %trans_info)
	{
		$chr = $trans_info{$tid}{'chr'};
		$strand = $trans_info{$tid}{'strand'};
		$gid = $trans_info{$tid}{'gid'};
		my @exon = split(/\t/, $trans_info{$tid}{'exon'});
		@exon = sort {$a<=>$b} @exon;
		$start = $exon[0];
		$end = $exon[scalar(@exon)-1];

		# generate gene info
		if (defined $gene_info{$gid}{'start'}) {
			$gene_info{$gid}{'start'} = $start if $start < $gene_info{$gid}{'start'};
		} else {
			$gene_info{$gid}{'start'} = $start;
		}
		
		if (defined $gene_info{$gid}{'end'}) {
			$gene_info{$gid}{'end'} = $end if $end > $gene_info{$gid}{'end'};
		} else {
			$gene_info{$gid}{'end'} = $end;
		}	

		if (defined $gene_info{$gid}{'chr'}) {
			die "[ERR]inconsistency chr for $gid $tid\n" if $chr ne $gene_info{$gid}{'chr'};
		} else {
			$gene_info{$gid}{'chr'} = $chr;
		}
		
		if (defined $gene_info{$gid}{'strand'}) {
			die "[ERR]inconsistency strand for $gid $tid\n" if $strand ne $gene_info{$gid}{'strand'};
		} else {
			$gene_info{$gid}{'strand'} = $strand;
		}

		if (defined $gene_info{$gid}{'tid'}) {
			$gene_info{$gid}{'tid'}.="\t".$tid;
		} else {
			$gene_info{$gid}{'tid'} = $tid;
		}
	}

	# convert to bed/tab/gff format
	my $out = IO::File->new(">".$output_file) || die $!;
	my %gid_uniq;

	foreach my $tid (sort keys %trans_info)
        {
                $chr = $trans_info{$tid}{'chr'};
                $strand = $trans_info{$tid}{'strand'};
                $gid = $trans_info{$tid}{'gid'};
                my @exon = split(/\t/, $trans_info{$tid}{'exon'});
                @exon = sort {$a<=>$b} @exon;
                $start = $exon[0];
                $end = $exon[scalar(@exon)-1];	

		if ($out_format eq 'bed') 
		{
			if ($bed_type eq 'transcript_id') {
				print $out $chr,"\t",$start-1,"\t",$end,"\t",$tid,"\t.\t",$strand,"\n";	
			} elsif ( $bed_type eq 'gene_id') {
				unless (defined $gid_uniq{$gid}) {
					my $gstart = $gene_info{$gid}{'start'};
					my $gend = $gene_info{$gid}{'end'};
					$gid_uniq{$gid} = 1;
					print $out $chr,"\t",$gstart-1,"\t",$gend,"\t",$gid,"\t.\t",$strand,"\n";
				}			
			} else {
				my $exon_num = 0;
				for(my $i=0; $i<@exon; $i=$i+2) {
					$exon_num++;
					print $out $chr,"\t",$exon[$i]-1,"\t",$exon[$i+1],"\t",$tid.".exon".$exon_num,"\t.\t",$strand,"\n";
				}
			}
		} 
		elsif ($out_format eq 'gff')
		{
			unless (defined $gid_uniq{$gid}) {
				my $gstart = $gene_info{$gid}{'start'};
				my $gend = $gene_info{$gid}{'end'};
				$gid_uniq{$gid} = 1;
				print $out "$chr\tGTFtool\tgene\t$gstart\t$gend\t.\t$strand\t.\tID=$gid;Name=$gid;\n";
			}

			print $out "$chr\tGTFtool\tmRNA\t$start\t$end\t.\t$strand\t.\tID=$tid;Name=$tid;Parent=$gid;\n";
			for(my $i=0; $i<scalar(@exon)-1; $i=$i+2)
			{
				my $e_start = $exon[$i];
				my $e_end = $exon[$i+1];
				my $e_num = ($i/2)+1;
				print $out "$chr\tGTFtool\tCDS\t$e_start\t$e_end\t.\t$strand\t.\tID=$tid-exon$e_num;Name=$tid-exon$e_num;Parent=$tid;\n";
			}
		}
		else 
		{
			die "[ERR]output format\n";
		}
	}
	$out->close;
}

=head2
 parse_gtf -- parse gtf file, return gtf information
=cut
sub parse_gtf
{
	my $input_file = shift;

	my %trans_info; # key: tid, chr, exon, gene, strand
	
	my $fh = IO::File->new($input_file) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);

		if ($a[2] eq 'exon') 
		{
			# analyzing the attribute 
			my %attr_value;
			my @b = split(/; /, $a[8]); 
			foreach my $b (@b) {
				$b =~ s/;$//;
				my @c = split(/ /, $b, 2);
				die "[ERR]attr $b in $a[8]\n" unless @c == 2;
				$c[1] =~ s/"//ig;
				$attr_value{$c[0]} = $c[1];
			}

			die "[ERR]No transcript id\n" unless defined $attr_value{'transcript_id'};
			die "[ERR]No gene id $_\n" unless defined $attr_value{'gene_id'};
			my ($tid, $gid) = ($attr_value{'transcript_id'}, $attr_value{'gene_id'});

			if ( defined $trans_info{$tid}{'chr'} ) {
				die "[ERR]inconsistency chr for $tid\n"	if $trans_info{$tid}{'chr'} ne $a[0];
			} else {
				$trans_info{$tid}{'chr'} = $a[0];
			}

			if ( defined $trans_info{$tid}{'gid'} ) {
				die "[ERR]inconsistency gid for $tid\n" if $trans_info{$tid}{'gid'} ne $gid;
			} else {
				$trans_info{$tid}{'gid'} = $gid;
			}

			if ( defined $trans_info{$tid}{'strand'} ) {
				die "[ERR]inconsistency strand for $tid\n" if $trans_info{$tid}{'strand'} ne $a[6];
			} else {
				$trans_info{$tid}{'strand'} = $a[6];
			}

			if ( defined $trans_info{$tid}{'exon'}) {
				$trans_info{$tid}{'exon'}.="\t".$a[3]."\t".$a[4];
			} else {
				$trans_info{$tid}{'exon'} = $a[3]."\t".$a[4];
			}
			
			if ( $a[3] > $a[4] ) {
				print "[WARN]exon swap $_\n";
			}
		}

		if ($a[2] eq 'CDS') {
			# analyzing the attribute
			my %attr_value;
			my @b = split(/; /, $a[8]);
			foreach my $b (@b) {
				$b =~ s/;$//;  
				my @c = split(/ /, $b, 2);
				die "[ERR]attr $b in $a[8]\n" unless @c == 2;
				$c[1] =~ s/"//ig;
				$attr_value{$c[0]} = $c[1];
            }

            die "[ERR]No transcript id\n" unless defined $attr_value{'transcript_id'};
            die "[ERR]No gene id\n" unless defined $attr_value{'gene_id'};
            my ($tid, $gid) = ($attr_value{'transcript_id'}, $attr_value{'gene_id'});
			
			if ( defined $trans_info{$tid}{'cds'}) {
				$trans_info{$tid}{'cds'}.="\t".$a[3]."\t".$a[4];
			} else {
				$trans_info{$tid}{'cds'} = $a[3]."\t".$a[4];
			}
		}
	}
	$fh->close;

	return %trans_info;
}

=head2
 gffread -- print function for gffread
=cut
sub gtf_gffread
{
	print qq'
USAGE of gffread, which included in tophat package

1. convert GFF to GTF 
	gtffread -T -o output.gtf input.gff

2. extract transcript from GTF
	gffread -w transcript.fasta -g genome.fasta input.gtf

';

	exit;
}

=head2
 percentile
=cut
sub percentile {
	my ($array, $p) = @_;
	my $num = scalar(@$array);
	my $r = $p * ($num+1);
	my @a = split(/\./, $r);
	my ($ir, $fr) = ($a[0], 0);
	$fr = $a[1] if (defined $a[1] && $a[1] > 0);
	$fr = "0.".$fr;
	my @sort_array = sort {$a<=>$b} @$array;
	my $precentile = $fr * ($sort_array[$ir-2] - $sort_array[$ir-1]) + $sort_array[$ir-1];
	return $precentile;
}

=head2
 run_cmd -- run command
=cut
sub run_cmd
{
	my ($cmd, $debug) = @_;
	print $cmd."\n";
	exit if $debug;
	system($cmd) && die "Error in CMD: $cmd\n";
}

=head2
 usage -- print usage information
=cut
sub usage
{
	print qq'

Program: GTFtools (Tools for GTF file)
Version: $version

USAGE: $0 <command> [options] 
Command:
	stats       statistics for GTF file 
	convert     convert GTF to BED/GFF format
	extract     extract GTF by list
	gffread     usage of gffread
	sortgff		sort GFF file by input gene ID
	togpd       convert GTF to GPD format
	tojunc      get junction from GTF
	compare     compare two GTF file, output shared and specific GTF

* the gtf file must have exon feature, transcript_id, and gene_id attributes

';
	exit;
}

