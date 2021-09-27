#!/usr/bin/env perl

=head1
# ==========================================================================
# smRNAtarget_Pred.pl: small RNA target prediction tool.
#
# This program returns potential small RNAs for an input 
# gene as well as potential sites for an input small RNA.  
#
# Usage: smRNAtarget_Pred.pl INPUT_FILE DB_FILE SOURCE_SEQ_FILE OUTPUT_FILE
# Example: smRNAtarget_Pred.pl query.fas tomato tomato.seq output.tab
#
# The Author should be Fei or Joung, modified by Yi
# 20151003 check_blastall(), do not need put blastll in bin  
# 20151003 make blast_formatted_db source_db(fasta) as one parameters 
# 20150423 reconstruct the program  
# 20150330 make it could analyze multiple files simultaneously
# ==========================================================================
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Std;

###########################################################
# Configuration section

# The path of BLAST tool

# The path of database files 
# 	- Save BLAST formatted database files of mRNAs and small RNAs 
# 	- In case of small RNAs, please make format database files 
#     with their reverse complement

###########################################################

my @sel_list_disc = ();
my @sel_list_pos_mRNA_st = ();
my @sel_list_pos_sRNA_st = ();
my @sel_list_pos_sRNA_end = ();
my @sel_list_direction = ();
my @description_list = ();

# $mRNA_fas = "";	#  Tomato unigene file (not use)
# $smRNA_fas = "";	#  Tomato smRNA	(not use)

=head2 insert_seq
 put fasta sequence to hash. 
 key: seqID ; value: sequence
=cut
sub insert_seq ($) {
	my $fn = $_[0];
	my %seq_hash;
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$fn);
	while(my $inseq = $in->next_seq) {
		my $sid = $inseq->id;
		my $seq = $inseq->seq;
		$seq_hash{$sid} = $seq;
	}
	return %seq_hash;
}

=head2 query_seq_from_file
 put fasta sequence to hash, 
 key: seqID ; value: sequence 
=cut
sub query_seq_from_file ($) {
	my $fn = $_[0];

	my %qseq_hash;
    my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$fn);
    while(my $inseq = $in->next_seq) {
        my $sid = $inseq->id;
        my $seq = $inseq->seq;
        $seq = uc($seq);
        $qseq_hash{$sid} = $seq;
        push (@description_list, $sid);
    }
    return %qseq_hash;
}

=head2 reverse_complement
 function: get the reverse complement sequence
=cut
sub reverse_complement {
	my $seq = shift;
	$seq =~ s/U/T/g;
	$seq =~ s/u/T/g;
	$seq = reverse($seq);
	$seq =~ tr/atcgATCG/tagcTAGC/;
	return $seq;
}

=head2 lsa (Local Sequence Alignment)
 perform local sequence alignment base on sw-algorithm
=cut
sub lsa {
	my ($adisc, $bdisc, $mRseq, $smRseq, $pos, $direction, $sc_cutoff, $max_score_cutoff) = @_;

	my $match = 1;
	my $wobble = 0.5; 
	my $mismatch = -1;
	my $indel = -2;

	my $first_i = 0;
	my $first_j = 0;
	my $last_i = 0;
	my $last_j = 0;

	my $offset = 3;
	my $mRNA_st = $sel_list_pos_mRNA_st[$pos]+0;
	my $sRNA_st = $sel_list_pos_sRNA_st[$pos]+0;
	my $sRNA_end = $sel_list_pos_sRNA_end[$pos]+0;
	my $start = $mRNA_st - ($sRNA_st + $offset);
	if ($direction eq "-") {
		$start = $mRNA_st - ((length($smRseq)-$sRNA_end) + $offset);
	}
	if ($start < 0) { $start = 0; }
	my $seqlen = length($smRseq)+$offset;
	my $end = $start+$seqlen;
	if ($end >= length($mRseq)) { $seqlen = length($mRseq)-$start+1; }

	my $a = substr($mRseq, $start, $seqlen);
	my $b = $smRseq;

	my $S = [];   # An array of scores
	my $R = [];   # An array of backtracking arrows
	my $n = length($a);
	my $m = length($b);

	# We need to work in letters, not in strings.  This is a simple way
	# to turn a string of letters into an array of letters.
	my @a = split // => $a;
	my @b = split // => $b;

	# These are "constants" which indicate a direction in the backtracking array.
	my $UP_AND_LEFT = "\\";
	my $UP          = "|";
	my $LEFT        = "-";
	my $NEITHER     = "o";

	my $max_score = 0;
	my $sc = 0.0;
	my $mc = 0;
	my $wc = 0;
	my $ic = 0;

	my $answer = 0;
	my $bind = "";
	
	# Initialization
	my $score = 0;
	for(my $j = 0; $j <= $m; $j++) { 
		$S->[0][$j] = 0;
		$R->[0][$j] = $NEITHER;
	}

	# This is the main dynamic programming loop that computes the score.
	for(my $i = 1; $i <= $n; $i++) { 
		$S->[$i][0] = 0;
		$R->[$i][0] = $NEITHER;
		for(my $j = 1; $j <= $m; $j++) { 
			$S->[$i][$j] = 0;
			$R->[$i][$j] = $NEITHER;

			if ( $a[$i-1] eq $b[$j-1] ) { $score = $match; }
			elsif ( is_wobble($a[$i-1], $b[$j-1], $direction) eq 1 ) { $score = $wobble; }
			else { $score = $mismatch; }

			if (($S->[$i-1][$j-1]+$score) > $S->[$i][$j]) {
				$S->[$i][$j] = $S->[$i-1][$j-1]+$score;
				$R->[$i][$j] =$UP_AND_LEFT;
			}
			if (($S->[$i][$j-1]+$indel) > $S->[$i][$j]) {
				$S->[$i][$j] = $S->[$i][$j-1]+$indel;
				$R->[$i][$j] =$LEFT;
			}
			if (($S->[$i-1][$j]+$indel) > $S->[$i][$j]) {
				$S->[$i][$j] = $S->[$i-1][$j]+$indel;
				$R->[$i][$j] =$UP;
			}
			if ($S->[$i][$j] > $max_score) {
				$max_score = $S->[$i][$j];
				$last_i = $i;
				$last_j = $j;
			}
		}
	}

	my $seqa = "";
	my $seqb = "";

	if ( $max_score > $max_score_cutoff ) {

		my $i = $last_i; 
		my $j = $last_j;

		# Trace the backtracking matrix.
		while( $i > 0 and $j > 0 ) {

			if( ($R->[$i][$j] eq $UP_AND_LEFT) and (($a[$i-1] eq $b[$j-1]) or (is_wobble($a[$i-1], $b[$j-1], $direction ) eq 1 )) ) {

				$seqa = $a[$i-1].$seqa;
				$seqb = $b[$j-1].$seqb;
				$i--; $j--;
			}
			elsif( $R->[$i][$j] eq $LEFT ) {
				$seqa = "-".$seqa;
				$seqb = $b[$j-1].$seqb;
				$j--;
			}
			elsif( $R->[$i][$j] eq $UP ) {
				$seqa = $a[$i-1].$seqa;
				$seqb = "-".$seqb;
				$i--;
			}
			else {
				$seqa = $a[$i-1].$seqa;
				$seqb = $b[$j-1].$seqb;
				$i--; $j--;
			}
		}

		$first_i = $i;
		$first_j = $j;

		# Insert postfixes of seqb 
		for(my $j = 1; $j <= length($smRseq)-$last_j; $j++) { 
			$seqb = $seqb.$b[$last_j+$j-1];
		}

		my $lp = get_last_pos($first_i, $seqa);

		# Insert prefixes 
		my $k = 1;
		my $ind = 0;
		for(my $j = $first_j-1; $j >= 0; $j--) { 
				if (($first_i-$k) >= 0) {
					$seqa = $a[$first_i-$k-1].$seqa;
				}	
				else {
					$seqa = "-".$seqa;
					$ind = 1;
				}
				$seqb = $b[$j].$seqb;
				$k++;
		}
		if ($ind eq 1) { $first_i = 0; }
		else { $first_i = $first_i-$k+1; }
		$first_j = $first_j-$k-1;

		# Delete postfix indels
		for(my $i = length($seqa)-1; $i >= 0; $i--) { 
			if ( substr($seqa,$i,1) eq "-" ) { 
				substr($seqa,$i,1) = "";
		 	}
			else { last; }
		}
		chomp($seqa);

		# Insert postfixes
		$k = 1;
		for(my $i = length($seqa); $i < length($seqb); $i++) { 
				if (($lp+$k) < length($a) )    {
					$seqa = $seqa.$a[$lp+$k];
				}
				else {
					$seqa = $seqa."-";
				}
				$k++;
		}

		# Scoring of binding pair
		my $ca = "";
		my $cb = "";
		for(my $i = 0; $i < length($seqb); $i++) { 
			$ca = substr($seqa,$i,1);
			$cb = substr($seqb,$i,1);
			if ( is_match($ca, $cb) eq 1 ) { $bind = $bind."|"; }
			elsif ( is_wobble($ca, $cb, $direction) eq 1 ) { $wc++; $sc=$sc+0.5; $bind = $bind."o"; }
			elsif ( is_indel($ca, $cb) eq 1 ) { $ic++; $sc=$sc+2.0; $bind = $bind."b"}
			else  { $sc=$sc+1.0;  $mc++; $bind = $bind."b"; }
		}
		if ( $sc <= $sc_cutoff ) { $answer = 1; }
	}

	my @seqinfo = ($answer, $sc, $mc, $start+$first_i+1, $first_j+1, $seqa, $seqb, $bind, $adisc, $bdisc, $direction, $wc, $ic);
	return @seqinfo;
}

sub get_last_pos (@) {
        my ($start, $seq) = @_;

        my $barc = 0;
	for(my $i = 0; $i < length($seq); $i++) { 
		if (substr($seq,$i,1) eq "-") { $barc++; }
	}
        my $last_pos = $start+length($seq)-$barc-1;
        return $last_pos;
}

sub is_match (@) {
	my ($ca, $cb) = @_;
	my $answer = 0;
	if ( $ca eq $cb ) { $answer = 1; }
	return $answer;
}

sub is_wobble (@) {
	my ($ca, $cb, $strand) = @_;
	my $answer = 0;

	#if ( $strand eq $pos_direction ) {
	if ( $strand eq '+' ) {		# change it by kentnf
		if ( (($ca eq "T") and ($cb eq "C")) or (($ca eq "G") and ($cb eq "A" )) ) {
			$answer = 1;
		}
	}
	else {
		if ( (($ca eq "C") and ($cb eq "T" )) or (($ca eq "A") and ($cb eq "G")) ) {
			$answer = 1;
		}
	}
	return $answer;
}

sub is_indel (@) {
	my ($ca, $cb) = @_;
	my $answer = 0;
	if ( (($ca eq "-") and ($cb ne "-")) or (($ca ne "-") and ($cb eq "-" )) ) {
		$answer = 1;
	}
	return $answer;
}

sub is_mismatch (@) {
	my ($ca, $cb) = @_;
	my $answer = 0;
	if ( $ca ne $cb ) { $answer = 1; }
	return $answer;
}

sub change_sequences (@) {
	my ($sad, $sbd, $sa, $sb, $is_mRNA) = @_;
	return 	($sad, $sbd, $sa, $sb);
}

sub split_discription ($) {
	my $sd = $_[0];
	my @list = split(" ", $sd);
	my $title = $list[0];
	return $title;
}

sub blast_input_save (@) {
	my ($disc, $qseq, $is_mRNA, $tempb) = @_;

	my $answer = 1;
        
	if ( $is_mRNA eq 0 ) {
		$qseq = reverse_complement($qseq);
	}

	chomp($disc);   				
	chomp($qseq);   				
	open (OUTFILE, ">$tempb") || die $!;
	print OUTFILE ">$disc\n$qseq\n";
	close OUTFILE;

	return $answer;
}

sub blast_run (@) {
	my ($blastdb, $queryfile, $direction) = @_;	
	my $blast_out = `blastall -p blastn -d $blastdb -i $queryfile -W 7 -q -1 -S $direction -e 100 -m 8`;
	#my $blast_out = `blastall -p blastn -d $blastdbname -i $queryfile -W 7 -q -1 -S $direction -e 200 -m 8`;
	return $blast_out;
}

sub blast (@) {
	my ($query_disc, $query_seq, $queryf, $is_mRNA, $search_direction, $blastdb_smRNA, $blastdb_mRNA) = @_;
	my $i = 0;
	my $db = "";
	my $forward = 1;
	my $reverse = 2;
	my $pos_direction = "+"; 
	my $neg_direction = "-"; 

	if ( $is_mRNA eq 1 ) { 
		$db = $blastdb_smRNA;
		my $rv = blast_input_save($query_disc, $query_seq, $is_mRNA, $queryf);

		if ( ($search_direction eq "forward") or ($search_direction eq "both") ) { 
			my $result = blast_run($db, $queryf, $forward);
			my @list = split(/\n/, $result);
			foreach my $line (@list) {
				my @entry = split(/\t/,$line);
				my $name = $entry[1];
				push @sel_list_disc, $name;
				push @sel_list_pos_mRNA_st, $entry[6];
				push @sel_list_pos_sRNA_st, $entry[8];
				push @sel_list_pos_sRNA_end, $entry[9];
				push @sel_list_direction, $pos_direction;
				$i++;
			}
		}
		if ( ($search_direction eq "reverse") or ($search_direction eq "both") ) { 

			my $result = blast_run($db, $queryf, $reverse);
			my @list = split(/\n/, $result);
			foreach my $line (@list) {
				my @entry = split(/\t/,$line);
				my $name = $entry[1];
				push @sel_list_disc, $name;
				push @sel_list_pos_mRNA_st, $entry[6];
				#push @f, $entry[9];
				push @sel_list_pos_sRNA_st, $entry[9];	# change from push @f, $entry[9];
				push @sel_list_pos_sRNA_end, $entry[8];
				push @sel_list_direction, $neg_direction;
				$i++;
			}
		}
	}
	else { 
		$db = $blastdb_mRNA;
		my $rv = blast_input_save($query_disc, $query_seq, $is_mRNA, $queryf);
		if ( ($search_direction eq "forward") or ($search_direction eq "both") ) { 

			my $result = blast_run($db, $queryf, $forward);
			my @list = split(/\n/, $result);
			foreach my $line (@list) {
				my @entry = split(/\t/,$line);
				my $name = $entry[1];
				push @sel_list_disc, $name;
				push @sel_list_pos_mRNA_st, $entry[8];
				push @sel_list_pos_sRNA_st, $entry[6];
				push @sel_list_pos_sRNA_end, $entry[7];
				push @sel_list_direction, $pos_direction;
				$i++;
			}		
		}

		if ( ($search_direction eq "reverse") or ($search_direction eq "both") ) { 

			my $result = blast_run($db, $queryf, $reverse);
			my @list = split(/\n/, $result);
			foreach my $line (@list) {
				my @entry = split(/\t/,$line);
				my $name = $entry[1];
				push @sel_list_disc, $name;
				push @sel_list_pos_mRNA_st, $entry[9];
				push @sel_list_pos_sRNA_st, $entry[6];
				push @sel_list_pos_sRNA_end, $entry[7];
				push @sel_list_direction, $neg_direction;
				$i++;
			}
		}
	}

	return $i;
}

sub print_tab {
	my ($schash, $resulthash) = @_;

	# unpack the hash
	my %schash = %$schash;
	my %resulthash = %$resulthash;

	my $summary = "";
	my $i = 1;
	foreach my $value (sort {$schash{$a} <=> $schash{$b}} keys %schash) {
		my ($istarget, $sc, $mc, $sa, $sb, $seqa, $seqb, $bind, $adisc, $bdisc, $direction, $wc, $ic) = @{$resulthash{$value}};
		my  $lp = get_last_pos($sa, $seqa);

		$bdisc =~ s/^\>//;
		$adisc =~ s/^\>//;
		my  $sRNA_ID = $bdisc;
		my  $unigene_ID = $adisc;
		
		if ( $direction eq "-") {
			$seqa = reverse_complement($seqa);
			$seqb = reverse($seqb);
        	$seqa =~ tr/tT/uU/;
        	$seqb =~ tr/tT/uU/;
			$bind = reverse($bind);
			my $tmp_pos = $sa;
			$sa =$lp;
			$lp = $tmp_pos;
		}
		else {
			$seqa =~ tr/tT/uU/;
			$seqb =~ tr/atcgATCG/uagcUAGC/;
		}
		my  $target_start = $sa;
		my  $align_string = $seqb;
		my  $align_score = $sc;
		my  $align_query = $seqa;
		my  $align_hit = $bind;
		my  $seq_strand = $direction;
		
		$summary.= "$sRNA_ID\t$unigene_ID\t$target_start\t$lp\t$align_score\t$align_query\t$align_string\t$align_hit\t$seq_strand\n";		
		$i++;
	}
	return $summary;
}

=head1 filter_output_result
=cut
sub filter_output_result
{
    my ($result_tab, $max_mismatch_cutoff, $not_mismatch_position) = @_;

    my @result_tabs = split(/\n/, $result_tab);
    my $result_tab_filtered = "";
    my @position = split(/\t/, $not_mismatch_position);

    foreach my $line (@result_tabs)
    {   
        my @a = split(/\t/, $line);
        my $align = reverse $a[7];

        my $mismatch = 0;
        for (my $i = 0; $i <= 8; $i++)
        {   
            if (substr($align, $i, 1) eq "b") { $mismatch++; }
        }

        my $mismatch_meet_position = 0;
        foreach my $pos (@position)
        {
            if ( substr($align,$pos,1) eq "b" ) { $mismatch_meet_position++; }
        }

        if ($mismatch < $max_mismatch_cutoff && $mismatch_meet_position == 0)
        {
            $result_tab_filtered.= $line."\n";
        }
    }

    return $result_tab_filtered;
}

=head2 check_blastall
=cut
sub check_blastall
{
    my $blast_ver = `blastall`;
    die "[ERR]need install blastall\n" unless ($blast_ver =~ m/blastall\s+\S+\s+arguments:/);
}

=head2 run_target_score
	main for miRNA target 
=cut
sub run_target_score 
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 sRNA_seq  reference_seq  output(options)

* sRNA_seq -- fasta file of sRNA/miRNA
* reference_seq -- fasta file of transcriptome, and blast db index
* output -- output_prefix (default: sRNA_seq.target.txt)

';
	print $usage and exit if @$files < 2;

	my $query_file  = $$files[0];
	my $target_file = $$files[1];
	my $out_file    = $query_file.".target.txt";
    $out_file = $$files[2] if defined $$files[2];
    warn "[WARN]output file exist\n" if -s $out_file;

	my $sc_cutoff = 4.5;				# cutoff of target score (0-10)
	my $max_score_cutoff = 8;			# cutoff value of alignment by dynamic programming
	my $seq_type = 0;					# type of query sequence (smRNA: 0, mRNA: 1)
	my $choice_direction = "forward";	# Direction of mRNA sequence (forward, reverse, both)
	my $tempb = "$query_file.blast.tmp";# Temorary file for blast input
	
	$sc_cutoff = $$options{'s'} if (defined $$options{'s'} && $$options{'s'} >= 0 && $$options{'s'} <= 10);	
	$max_score_cutoff = $$options{'m'} if (defined $$options{'m'} && $$options{'m'} >= 0);
	$seq_type = $$options{'p'} if (defined $$options{'p'} && $$options{'p'} == 1);
	$choice_direction = $$options{'d'} if defined $$options{'d'};

	# set blast smRNA and blast mRNA base on seq type
	my ($blastdb_smRNA, $blastdb_mRNA);
	if ($seq_type eq 1) { 
		$blastdb_smRNA = $target_file;	# smRNA as database
		$blastdb_mRNA  = 'NA';
	} else { 
		$blastdb_mRNA  = $target_file;	# reference(mRNA) as database
		$blastdb_smRNA = 'NA';
	}

	# load sRNA and mRNA sequence to hash
	my %seq_querys = query_seq_from_file($query_file);	# query seq to hash
	my %seqsa = insert_seq($target_file);				# target seq to hash

	# set the search direction
	my $direction = "+";
	my $search_direction = "both";
	
	if ($choice_direction eq "forward") { $search_direction = "forward"; }
	elsif ($choice_direction eq "reverse") { $search_direction = "reverse"; }	

	# main	
	my $seq_num = 1;
	
	my %resulthash = ();	# create hash for result
	my %schash = ();		# create hash for score cutoff
	my @index_list = ();

	my $output_content = "";
	foreach my $description (@description_list) {
		
		my $query_n = $description;
		print $seq_num."-th input query: ".$query_n;

		my $result_tab = "";
		
		my @seqb = ($description, $seq_querys{$description});
	
		if ( (length($seqb[1]) <= 30) and ($seq_type eq 1) ) { print "\nPlease check the query sequence type!";  exit; }
		if ( (length($seqb[1]) > 30) and ($seq_type eq 0) ) { print "\nPlease check the query sequence type!";  exit; }

		my $smRNA = "";
		my $smRNArc = "";
		my $mRNA = "";
		my $mRNAdisc = "";
		my $smRNAdisc = "";
		
		@sel_list_disc=();
		@sel_list_pos_mRNA_st=();
		@sel_list_pos_sRNA_st=();
		@sel_list_pos_sRNA_end=();
		@sel_list_direction=();

		my $count = blast($seqb[0], $seqb[1], $tempb, $seq_type, $search_direction, $blastdb_smRNA, $blastdb_mRNA);
		
		my $i = 0;

		foreach my $disc (@sel_list_disc) {

			$direction = $sel_list_direction[$i];

			$mRNAdisc = split_discription($mRNAdisc);
			if ( $seq_type eq 1 ) {
				($mRNAdisc, $smRNAdisc, $mRNA, $smRNA) = change_sequences($seqb[0], $disc, $seqb[1], $seqsa{$disc}, $seq_type);
			}
			else {
				($mRNAdisc, $smRNAdisc, $mRNA, $smRNA) = change_sequences($disc, $seqb[0], $seqsa{$disc}, $seqb[1], $seq_type);
			}

			if ( $direction eq "+" ) { $smRNArc = reverse_complement($smRNA); }
			else { $smRNArc = $smRNA; }

			chomp($smRNArc);
			chomp($mRNA);

			# ($answer, $sc, $mc, $start+$first_i+1, $first_j+1, $seqa, $seqb, $bind, $adisc, $bdisc, $direction, $wc, $ic);
			#  1        3.5  3    113                -1                               Csa3M778220.1 M00001 + 1, 0
			my @sinfo = lsa($mRNAdisc, $smRNAdisc, $mRNA, $smRNArc, $i, $direction, $sc_cutoff, $max_score_cutoff);
			if ( $sinfo[0] eq 1 ) {
				my $start_pos=sprintf("%d", $sinfo[3]);
				my $index = join "_", $mRNAdisc, $smRNAdisc, $start_pos;
				$resulthash{$index} = \@sinfo;
				$schash{$index} = $sinfo[1];
				push (@index_list, $index);
			}
			$i++;
		}

		$result_tab = print_tab(\%schash, \%resulthash);
		# filter the miR target result tab info
		# allowing max two mismatches for the alignment
		# the 9-th, 10-th position do not have mismatch
		if ($result_tab) {
			$result_tab = filter_output_result($result_tab, 2, "9\t10");
			$output_content.=$result_tab;
		}
		
		$seq_num = $seq_num + 1;
		print "...Done\n";
		
		# delete the array and hash
		foreach my $index (@index_list) {
			delete $resulthash{$index}; 
			delete $schash{$index};
		}
		@index_list = ();
	}

	# output file
	open (OUTFILE, ">$out_file") || die $!;
	print OUTFILE $output_content;
	close OUTFILE;

	# remove temp file
	unlink($tempb);
}

my %options;
getopts('a:b:c:d:e:g:i:j:k:l:m:n:o:p:q:r:s:t:v:w:x:y:z:fuh', \%options);
run_target_score(\%options, \@ARGV);	# add input files 	 
