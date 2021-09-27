#!/usr/bin/perl

=head2
 junctool --- generate and compare junction in GTF and sam file
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Std;

my $version = 0.1;
my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }

if		($options{'t'} eq 'sam2junc')		{ jc_sam2junc(\%options, \@ARGV); }
elsif	($options{'t'} eq 'juncComp')		{ jc_juncComp(\%options, \@ARGV); }
elsif	($options{'t'} eq 'shareDep')		{ jc_shareDep(\%options, \@ARGV); }
elsif	($options{'t'} eq 'juncGTF')		{ jc_juncGTF(\%options, \@ARGV); }
elsif	($options{'t'} eq 'hisat2junc')		{ jc_hisat2junc(\%options, \@ARGV); }
else	{ usage($version); }

sub usage {
	my $usage = qq'
USAGE: $0 -t [tool]
	sam2junc -- extract junc from sam
	juncComp -- compare junction
	shareDep -- get the depth of shared junctions
	hisat2junc -- convert hisat junc format to my format
	juncGTF  -- filter GTF using junction 
				(the transcript supported by junc will be kept)
';
	print $usage;
	exit;
}

#########################
# ====== kentnf ======= #
#########################
=head2
 jc_hisat2junc
=cut
sub jc_hisat2junc
{
	my ($options, $files) = @_;
    my $usage = qq'
USAGE: $0 -t hisat2junc [options] hisat_junc.txt > output_junc.txt

* convert hisat junc format to my format

';
	print $usage and exit unless defined $$files[0];

	open(FH, $$files[0]) || die $!;
	while(<FH>) {
		chomp;
		my @a = split(/\t/, $_);
		$a[1] = $a[1] + 2;
		print $a[0]."#".$a[1]."#".$a[2]."\n";
	}
	close(FH);
}

=head2
 jc_juncGTF
=cut
sub jc_juncGTF
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t juncGTF [options] junc1 input.gtf > output.gtf

	-a	type: one/all [default: one]
	-e	isform expression [for filter single exon transcript]
	-f	cutoff for FPKM [default: 0.1]
	-r	cutoff for replicate num [default: 2]
	-x	cutoff for multi-exon expression [default: 0]

one: the transcript in GTF will be selected if one junction include in input junction file
all: the transcript in GTF will be selected if all junctions include in input junction file
** the 2 replicates of all sample expressed > 0.1 will be selected as expressed single exon transcript
';
	# load junction to hash
	print $usage and exit unless scalar @$files == 2;
	my ($junc_file, $input_gtf) = @$files;

	my $type = 'one';
	$type = 'all' if (defined $$options{'a'} && $$options{'a'} eq 'all');

	# load expression to hash
	my $exp_cutoff = 0.1;
	$exp_cutoff = $$options{'f'} if (defined $$options{'f'} && $$options{'f'} > 0);

	my $rep_cutoff = 2;
	$rep_cutoff = $$options{'r'} if (defined $$options{'r'} && $$options{'r'} >=2);

	my $iso_exp_cutoff = 0;
	$iso_exp_cutoff = $$options{'x'} if (defined $$options{'x'} && $$options{'x'} > 0);

	# key: sample_name; value: 1
	my %sample;
	my $tid_num = 0;

	# key1: tid;
	# key2: sample;
	# value: number of replicate meet the exp_cutoff
	my %tid_exp;
	
	if (defined $$options{'e'} && -s $$options{'e'}) 
	{
		my $efh = IO::File->new($$options{'e'}) || die $!;
		my $title = <$efh>;	chomp($title);		
		my @t = split(/\t/, $title);
		while(<$efh>) {
			chomp;
			my @a = split(/\t/, $_);
			my $tid = $a[0];
			$tid_num++;
			for(my $i=1; $i<@a; $i++) {
				my $sample = $t[$i];
				$sample =~ s/_\d+$//;
				$sample{$sample} = 1 unless defined $sample{$sample};
				my $expression = $a[$i];
				next if $expression < $exp_cutoff;

				if (defined $tid_exp{$tid}{$sample}) {
					$tid_exp{$tid}{$sample}++;
				} else {
					$tid_exp{$tid}{$sample}=1;
				}
			}
		}
		$efh->close;
	}

	# filter the tid by expression
	# key: tid; value: 1
	my %tid_exp_select;
	foreach my $tid (sort keys %tid_exp) {
		my $nn = 0;
		foreach my $s (sort keys %sample) {
			$nn++ if (defined $tid_exp{$tid}{$s} && $tid_exp{$tid}{$s} >= $rep_cutoff);
		}
		$tid_exp_select{$tid} = 1 if $nn == scalar(keys(%sample));
	}

	print STDERR "sample: ";
	foreach my $s (sort keys %sample) { print STDERR " $s"; }
	print STDERR "\n";
	print STDERR "tid: $tid_num\n";
	print STDERR "tid (exp > $exp_cutoff):".scalar(keys(%tid_exp_select))."\n";

	# load junction to hash
	my %junc;
	my $fh = IO::File->new($junc_file) || die $!;
	while(<$fh>) {
		chomp;
		my @a = split(/\t/, $_);
		$junc{$a[0]} = 1;
	}
	$fh->close;

	# parse trans info
	my %trans_supp;		# key: tid: value: 1	-- splicing junction supported trancript
	my %final_tid;
	my %trans_supp_se;	# key: tid: value: 1	-- expression supported single exon
	my %trans_info = parse_gtf($input_gtf);
	foreach my $tid (sort keys %trans_info)	{
		my $chr = $trans_info{$tid}{'chr'};
		my $fpkm = $trans_info{$tid}{'fpkm'};
		my @exon = split("\t", $trans_info{$tid}{'exon'});
		@exon = sort {$a<=>$b} @exon;

		# parse single exon transcript
		if (@exon == 2) {
			if (defined $$options{'e'} && defined $tid_exp_select{$tid} && $tid_exp_select{$tid} == 1) {
				$trans_supp_se{$tid} = 1;
			}
			next;
		}

		die "[ERR]exon num $tid\n" unless (scalar(@exon) % 2 == 0);
		shift @exon;
		pop @exon;	

		# check if support
		my $support;

		if ($type eq 'one') {
			$support = 0;
			for(my $i=0; $i<@exon; $i=$i+2) {
				my $start = $exon[$i] + 1;
				my $end = $exon[$i+1] - 1;
				$support = 1 if defined $junc{$chr."#".$start."#".$end};
			}
		}
		else {
			$support = 1;
			for(my $i=0; $i<@exon; $i=$i+2) {
				my $start = $exon[$i] + 1;
				my $end = $exon[$i+1] - 1;
				$support = 0 unless defined $junc{$chr."#".$start."#".$end};
			}
		}

		if ($support == 1) {
			$trans_supp{$tid} = 1;
			if ($fpkm >= $iso_exp_cutoff) {
				$final_tid{$tid} = 1;
			}
		}
	}

	print STDERR "Supported Trans Num:".scalar(keys(%trans_supp))."\n";
	print STDERR "Supported Trans Num (Single Exon):".scalar(keys(%trans_supp_se))."\n";
	print STDERR "Final Trans Num:".scalar(keys(%final_tid))."\n";

	# output select transcript
    my $in = IO::File->new($input_gtf) || die $!;
	while(<$in>)
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
		die "[ERR]No gene id $_\n" unless defined $attr_value{'gene_id'};
		my ($tid, $gid) = ($attr_value{'transcript_id'}, $attr_value{'gene_id'});

		if (defined $final_tid{$tid} || defined $trans_supp_se{$tid}) {
			print $_."\n";
			$final_tid{$tid} = 1;
		}
	}
	$in->close;
}

=head2
 jc_shareDep
=cut
sub jc_shareDep
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t shareExp [options] junc1 junc2 juncN ...

	-n	number of supported samples as share [default: input junc1]
	-o  output prefix [default: juncShareDep]

';
	# check input file
	print $usage and exit if scalar @$files < 2;
	foreach my $f (@$files) {
		print "[ERR]file not exist $f\n" unless -s $f;
	}

	my $number = scalar @$files;
	$number = $$options{'n'} if (defined $$options{'n'} && $$options{'n'} > 1 && $$options{'n'} <= (scalar @$files));

	my $prefix = "juncShareDep";
	$prefix = $$options{'o'} if defined $$options{'o'};

	# main
	# save the depth to hash.
	# key1: juncID; key2: juncFileName; value: depth 	
    my %depth;
    foreach my $m (@$files) {
		my $in = IO::File->new($m) || die $!;
		while(<$in>) {
			chomp;
			my @a = split(/\t/, $_);
			$depth{$a[0]}{$m} = $a[1];
		}
		$in->close;
    }

    # find share, and output result 
	my $out1 = IO::File->new(">".$prefix."_share_dep.txt") || die $!;
	my $out2 = IO::File->new(">".$prefix."_spec_dep.txt") || die $!;

    print $out1 "juncID";
    print $out2 "juncID\tsuppNum";
	foreach my $m (@$files) {
        print $out1 "\t$m";
        print $out2 "\t$m";
    }
    print $out1 "\n";
    print $out2 "\n";

	foreach my $j (sort keys %depth) {
        my $c = 0;
        my $e = '';
        foreach my $m (@$files) {
            if (defined $depth{$j}{$m}) {
                $c++;
                $e.="\t$depth{$j}{$m}";
            } else {
				$e.="\t0";
			}
        }

		if ($c >= $number) {
			print $out1 $j.$e."\n";
		}
		else {
			print $out2 $j."\t$c"."$e\n";
		}
	}

	$out1->close;
	$out2->close;
}

=head2
 jc_juncComp
=cut
sub jc_juncComp
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t juncComp [options] junc1 junc2 junc3

	-c cutoff of junction count
	-s print select junction with select flag 

';
	print $usage and exit unless defined $$files[1];

	my $cutoff = 0;
	$cutoff = $$options{'c'} if defined $$options{'c'} and $$options{'c'} > 0;
	my $select_flag;
	$select_flag = $$options{'s'} if defined $$options{'s'} and $$options{'s'} > 0;
	# put junction to hash
	# key: junction
	# value: flag in samples
	# 	1: only in junc1 file
	# 	2: only in junc2 file
	# 	4: only in junc3 file
	# 	3: in junc1 and junc2
	# 	5: in junc1 and junc3
	# 	6: in junc2 and junc3
	# 	7: in all junc file
	my %jc_hash;

	# defined the flag description
	my %flag_desc;
	$flag_desc{'1'}	= $$files[0];
	$flag_desc{'2'} = $$files[1];
	$flag_desc{'3'} = "$$files[0] + $$files[1]";

	if (defined $$files[2]) {
		$flag_desc{'4'} = $$files[2];
		$flag_desc{'5'} = "$$files[0] + $$files[2]";
		$flag_desc{'6'} = "$$files[1] + $$files[2]";
		$flag_desc{'7'} = "$$files[0] + $$files[1] + $$files[2]" ;
	}

	my $flag = 0;
	foreach my $f (@$files)
	{
		if ($flag == 0) { $flag = 1; }
		else { $flag = $flag + $flag;}

		my $fh = IO::File->new($f) || die $!;
		while(<$fh>) {
			chomp;
			my @a = split(/\t/, $_);
			if ($a[1] > $cutoff) {
				if (defined $jc_hash{$a[0]}) {
					$jc_hash{$a[0]}+=$flag;
				} else {
					$jc_hash{$a[0]}=$flag;
				}
			}
		}
		$fh->close;
	} 

	# stat the flag we have
	my %flag_stat;
	foreach my $jc (sort keys %jc_hash) {
		my $flag = $jc_hash{$jc};

		if (defined $select_flag && $flag == $select_flag) {
			print $jc."\n";
		}

		defined $flag_stat{$flag} and $flag_stat{$flag}++ or $flag_stat{$flag}=1;
	} 

	# report result
	foreach my $flag (sort {$a<=>$b} keys %flag_stat) {
		print $flag."\t".$flag_stat{$flag}."\t".$flag_desc{$flag}."\n";
	}
} 

=head2
 sam2junc
=cut
sub jc_sam2junc
{
	my ($options, $files) = @_;
    my $usage = qq'
USAGE: $0 -t sam2junc [options] input.sam > output.txt

';
    print $usage and exit unless defined $$files[0];

	my $sr_sam = $$files[0];
	#####################################
	# load junctions					#
	#####################################
	# load illumina junctions by reads or by alignment software
	my %sr_junction;
	if (-s $sr_sam) {
		my $fh1 = IO::File->new($sr_sam) || die $!;
		while(<$fh1>)
		{
			chomp;
			next if $_ =~ m/^@/;
			my @a = split(/\t/, $_);
			next if $a[1] eq '4';
			my $apos = $a[3];
			my $cigar = $a[5];
			next if $cigar =~ m/S/;	# filter out soft clipping
			my $ed = 0;
			if ($_ =~ m/NM:i:(\d+)/i) { $ed = $1; } else { die "[ERR]in NM $_\n"; }
			#next if $ed > 0;	# filter out reads with mismatch
			my ($junction, $insertion) = parse_cigar($cigar);
			next if scalar @$junction == 0;
			foreach my $j (@$junction) {
				my $jstart = $apos + $j->[1];
				my $jend   = $apos + $j->[2];
				if (defined $sr_junction{$a[2]."#".$jstart."#".$jend} ) {
					$sr_junction{$a[2]."#".$jstart."#".$jend}++;
				} else {
					$sr_junction{$a[2]."#".$jstart."#".$jend} = 1;
				}
			}
		}
		$fh1->close;
	}

	foreach my $jk (sort keys %sr_junction) {
		print $jk."\t".$sr_junction{$jk}."\n";
	}
}

# child of pb_correctF
# parse cigar to get junction cite
#            |     | <--- will get this pos 
# +++++++++++-------+++++++
sub parse_cigar
{
	my $cigar = shift;

	my @junction = ();
	my %insertion = ();
	my $num = ""; 
	my $point0 = 0; # point for junction site in transcriptome
	my $point1 = 0;	# point for junction site in genome
	my $point2 = 0; # point for insertion in transcriptome

	for(my $i=0; $i<length($cigar); $i++)
	{
		my $str = substr($cigar, $i, 1);
		if ($str =~ m/\d+/) {
			$num = $num.$str;
		}
		elsif ($str eq "M")
		{
			$point0 = $point0 + $num;
			$point1 = $point1 + $num;
			$point2 = $point2 + $num;
			$num = "";
		}
		elsif ($str eq "N")
		{
			$point1 = $point1 + $num;
			push(@junction, [$point0, $point1-$num, $point1-1]);
			$num = "";
		}
		elsif ($str eq "D")
		{
			$point0 = $point0 + $num;
			$point1 = $point1 + $num;
			$num = "";
		}
		elsif ($str eq "I")
		{
			$insertion{$point2+1} = $num;
			$point2 = $point2 + $num;
			$num = "";
		}
		else 
		{
			$num = 0;
		}
	}
	return(\@junction, \%insertion);
}

# parse the GTF
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

		if ($a[2] eq 'transcript')
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
			die "[ERR]No FPKM $_\n" unless defined $attr_value{'FPKM'};

			my ($tid, $gid, $fpkm) = ($attr_value{'transcript_id'}, $attr_value{'gene_id'}, $attr_value{'FPKM'});
			$trans_info{$tid}{'fpkm'} = $fpkm;
		}
        
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
				die "[ERR]inconsistency chr for $tid\n" if $trans_info{$tid}{'chr'} ne $a[0];
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
