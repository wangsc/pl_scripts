#!/usr/bin/perl -w
use strict;
use Statistics::ChisqIndep;

my $sUsage = "perl $0 <input>\n";
die $sUsage unless @ARGV;

my $input = shift;
open (IN, $input) or die "can't open file \n";

my $chisq_test = new Statistics::ChisqIndep;
while (<IN>)
{
	chomp; 
	my @t = split /\s+/, $_; 
	#cluster_id      pvalue_A        pvalue_B        pvalue_D        count_A count_B count_D test_AB_D       A_effect        B_effect        D_effect
	#RAC875_mira1_rep_c113520        NA      0.17676 -9.2508e-10     NA      290/504 0/95    -1.8469e-15     A_NA    B_BAL   D_SILENCE
	
	my ($id, $count_A, $count_B) = @t[0, 4, 5];
	if(/cluster/)
	{
		print join("\t", ($id, $count_A, $count_B, "A/B","pvalue", "A_effect", "B_effect")), "\n";
		next;
	}
	my @arr = (0, 0);
	unless($t[4] =~ /NA/)
	{
		my @tmp = split /\//,$t[4];
		$arr[0] += $tmp[0];
		$arr[1] += $tmp[1];
	}
	unless($t[5] =~ /NA/)
	{
		my @tmp = split /\//,$t[5];
		$arr[0] += $tmp[1];
		$arr[1] += $tmp[0];
	}
	
	my $total = $arr[0] + $arr[1];
	next if $total == 0;
	$chisq_test->load_data([[@arr], [int($total/2), int($total/2)]]);
	my $pvalue = $chisq_test->p_value;
	my @effect = qw(A_BAL B_BAL);
	if($pvalue < 0.05)
	{
		if($arr[0]>$arr[1])
		{
			@effect = qw(A_UP B_DOWN);
			$effect[1] = "B_SILENCE" if $arr[1] < 0.1*$arr[0];
		}
		else
		{
			@effect = qw(A_DOWN B_UP);
			$effect[0] = "A_SILENCE" if $arr[0] < 0.1*$arr[1];			
		}
	}
	print join("\t", ($id, $count_A, $count_B, $arr[0]."/".$arr[1], $pvalue, @effect)), "\n";
}