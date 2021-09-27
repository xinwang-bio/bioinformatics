#!/usr/bin/env perl
use warnings;
use strict;
use Excel::Writer::XLSX;


for my $i(@ARGV){
	my $input = $i;
	my $output = $input;
	if($output =~/\.xlsx$/){
		
	}
	else{
		$output =~s/\.txt$//i;
		$output =$output."."."xlsx";
		my $workbook = Excel::Writer::XLSX->new($output);
		my $worksheet = $workbook->add_worksheet();
		#my $format = $workbook->add_format();
		#$format->set_bold();
		#$fo0rmat->set_color( 'red' );
		open IN,"<",$input;
		my $row =0;
		while(<IN>){
			chomp;
			my @line = split/\t/,$_;
			my $col =0;
			for my $m(@line){
				$worksheet->write( $row, $col, $m);
				$col++;
			}
			$row++;
		}
		close IN;			
	}
}


