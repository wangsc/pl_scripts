#!/usr/bin/perl -w
use strict;
my $sUsage = qq(
perl $0
<ABD_allele.out>
<snps_fitBinomTest.list>
);

die $sUsage unless @ARGV == 2;
my ($allele_file, $snp_file) = @ARGV;

my %snps = read_snp_list($snp_file);

open (IN, $allele_file) or die "can't open file $allele_file\n";
while (<IN>)
{
	chomp; 
	# Excalibur_mira1_c3994   1849    GG      C
	# Excalibur_mira1_c93480  108     -       A
	my @t = split /\s+/, $_;
	my $id = join("\t", @t[0, 1]);
	my ($allele_m, $allele_l) = @{$snps{$id}};
	if (/-/)
	{
		if($t[2] eq '-')
		{
			if($allele_l eq $t[3])
			{
				$t[2] = $allele_m.$allele_m;
			}
			else
			{
				$t[2] = $allele_m.$allele_l;
			}
		}
		elsif($t[3] eq '-')
		{
			my @arr = split //, $t[2];
			if($arr[0] eq $arr[1])
			{
				$t[3] = $allele_l if $arr[0] eq $allele_m;
			}
			else
			{
				$t[3] = $arr[0] eq $allele_m?$arr[0]:$arr[1]
			}
		}
		print join("\t", @t), "\n" unless /-/;
	}
	else
	{
		my %h;
		map{$h{$_}++}split //, $t[2].$t[3];
		next if (keys %h) == 1;
		my @arr = sort{$h{$a}<=>$h{$b}} keys %h;		
		print $_, "\n" if $arr[0] eq $allele_l and $arr[1] eq $allele_m;
	}	
}
close IN;



sub read_snp_list
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while (<IN>)
	{
		chomp;
		#BobWhite_mira1_c1       39      T       A       49      22      72      0.319444444444444
		my @t = split /\s+/, $_; 
		$return{join("\t", @t[0,1])} = $t[4]>$t[5]?[@t[2,3]]:[@t[3,2]];
	}
	close IN;
	return %return;
}






