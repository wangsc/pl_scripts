#!/usr/bin/perl -w
use strict;
use Statistics::ChisqIndep;

my $sUsage = qq(

perl $0 
<homeolog_variation_sites.out>
<vcf file>
<minimum coverage depth>
<output file>
);

die $sUsage unless @ARGV >= 4;
my ($variation_file, $vcf_file, $mid_dep, $output) = @ARGV;

my $chisq_test = new Statistics::ChisqIndep;
my %snps = read_homeolog_variations($variation_file);
open (OUT, ">$output") or die;
print OUT join("\t", qw(SNP Allele_A Allele_B Depth_A Depth_B Pvalue1 Pvalue2)), "\n";

open (V, "$vcf_file") or die;

while (<V>)
{
	# BobWhite_mira1_c10      610     .       G       A       81      .       DP=54;VDB=0.0181;AF1=1;AC1=2;DP4=4,1,36,11;MQ=20;FQ=-90;PV4=1,5.1e-16,1,1       GT:PL:GQ        1/1:114,63,0:99
	next if /^\#/;
	my @t = split /\s+/,$_;
	next if length $t[4] > 1;
	my $id = join(":", @t[0,1]);
	next unless exists $snps{$id};
	my ($rf, $rb, $af, $ab) = $_=~/DP4=(\d+)\,(\d+)\,(\d+)\,(\d+)/;
	my $total = sum($rf, $rb, $af, $ab);
	next if $total < $mid_dep;
	my $ref = $rf + $rb;
	my $alt = $af + $ab;
	$chisq_test->load_data([[sort{$a<=>$b}($ref, $alt)], [int($total/3), int($total*2/3)]]);
	my $pvalue_1 = $chisq_test->p_value;
	my $pvalue_2;
	if($snps{$id} eq $t[3])
	{
		$chisq_test->load_data([[$ref, $alt], [int($total/3), int($total*2/3)]]);
		$pvalue_2 =  $chisq_test->p_value;
	}
	else
	{
		$chisq_test->load_data([[$alt, $ref], [int($total/3), int($total*2/3)]]);
		$pvalue_2 =  $chisq_test->p_value;		
	}
	print OUT join("\t", (join(":", @t[0,1]), @t[3,4], $ref, $alt, $pvalue_1, $pvalue_2)), "\n";
}

close OUT;
close V;


sub read_homeolog_variations
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while (<IN>)
	{
		# Kukri_mira1_c31978:394:809      C       G
		# RAC875_mira1_c35029:1161:1390   G       C,G
		my $id = $1 if(/^(\S+?:\d+)/);
		my @t = split /\s+/,$_; 
		if(length $t[1] == length $t[2])
		{
			$return{$id} = $t[1]
		}
		else
		{
			my @arr = split /,/, $t[2];
			foreach (@arr)
			{
				$return{$id} = $_ unless $t[1] eq $_;
			}
		}
	}
	close IN;
	return %return;
}

sub sum
{
	my $return = 0;
	map{$return += $_} @_;
	return $return;
}