#!/usr/bin/env perl
use strict;

use Getopt::Long;


my ( $rr, $ff, $RR, $ii, $II, $dd, $vv);
 
GetOptions(
        'f!'  => \$ff,
        'r!'     => \$rr,
        'R!'   => \$RR,
        'i!'     => \$ii,
        'I!'    => \$II,
        'd!'    => \$dd,
        'v!' => \$vv,
);
 



for my $i(@ARGV){
	my @line = split/\//,$i;
	my $file_name = $line[-1];
	if(-e "/data/xinwang/.Trash/$file_name"){
		for my $c(1..99999999){
			my $target_name = $file_name."(".$c.")";
			my $target = $file_name."\\(".$c."\\)";
			#print "$target_name\n";
			if(-e "/data/xinwang/.Trash/$target_name"){
				next;
			}
			else{
				system "mv $i /data/xinwang/.Trash/$target";
				print "$i\n";
				last;
			}
			
		}
	}
	else{
		system "mv $i /data/xinwang/.Trash/$file_name";
		print "$i\n";
	}
}