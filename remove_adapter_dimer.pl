#!/usr/bin/perl -w
use strict;

my $adapter_seq = "AGATCGGAAGAG";

my @files = <*.fastq_trimed_filtered>;
die "No file matched in current directory\n" unless @files;
foreach my $f(@files)
{
	print STDERR "Processing file $f ...\n";
	my $out_f = $f . '_AdaptorRemoved';
	open (IN, "$f") or die $!;
	open (OUT, ">$out_f") or die $!;
	my $count = 0;
	my @temp;
	while(<IN>)
	{
		next if /^\s+$/;
		chomp;
		push @temp, $_;
		if (@temp == 4)
		{
			my $seq = $temp[1];
			my $score = $temp[3];
			my @ad_pos; 
			while($seq=~/$adapter_seq/g)
			{
				my $pos = pos $seq;
				push @ad_pos, $pos;
			}
			@ad_pos = sort{$a<=>$b} @ad_pos;
			$temp[1] = defined $ad_pos[0]?substr($seq, 0, ($ad_pos[0]-(length $adapter_seq))):$seq;
			$temp[3] = substr($score, 0, length $temp[1]);
			print OUT join("\n", @temp),"\n" if length $temp[1] >= 50;
			@temp = ();
		}
	}
	close IN;
	close OUT;
}