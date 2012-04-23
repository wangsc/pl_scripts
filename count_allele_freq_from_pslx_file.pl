#!/usr/bin/perl -w
use strict;
use File::Basename;

my $sUsage = qq(

psLayout version 3

match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts  tStarts
        match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
---------------------------------------------------------------------------------------------------------------------------------------------------------------
43      2       0       0       0       0       0       0       +       MERCURE_0135:1:1101:1228:1961#CGATGT/1  101     56      101     tplb0024d04:398:957_G   61      0       45      1       45,     56,   0,       aaggccccctacgtggccaaggccaacaagctcaagggcgaggac,  aaggccccctacgtggccaaggccaacaaggtcaagggcgagtac,
44      1       0       0       0       0       0       0       +       MERCURE_0135:1:1101:1228:1961#CGATGT/1  101     56      101     tplb0024d04:398:957_C   61      0       45      1       45,     56,   0,       aaggccccctacgtggccaaggccaacaagctcaagggcgaggac,  aaggccccctacgtggccaaggccaacaagctcaagggcgagtac,
57      1       0       1       1       1       1       2       +       MERCURE_0135:1:1101:2741:1972#CGATGT/1  101     27      87      CAP12_mira1_rep_c7107:284:471_C 61      0       61      2       17,42,27,45,   0,19,   gttgttgacggggtcgg,gaggtggtcggcgaggttctcgagggggcccttgccggtgac,   gttgttgacngggtcgg,gaggtggtcggcgaggttctcgagtgggcccttgccggtgac,
....

Usage:
perl $0 <pslx file(s)>
);

die $sUsage unless @ARGV > 0; 

my @files = @ARGV;

foreach my $f (@files)
{
	my %allele_freq = count_allele_freq($f);
	foreach  my $id (keys %allele_freq)
	{
		my @tmp = map{$_."_".$allele_freq{$id}{$_}}keys %{$allele_freq{$id}};
		print join("\t", ($id, @tmp)), "\n";
	}
}

sub count_allele_freq
{
	my $file = shift;
	open (IN, $file) or die "can't open file $file \n";
	my %return;
	my $flag = 0;
	while(<IN>)
	{
		$flag = 1 if /^----/;
		next unless $flag;
		my @t = split /\s+/,$_; 
		my $target_id = $t[13];
		my ($target_start, $target_end) = @t[15, 16];
		next unless $target_start < 30 and $target_end > 30;
	}
	
}
