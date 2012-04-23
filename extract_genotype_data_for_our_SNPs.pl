#!/usr/bin/perl -w
use strict;
# extract genotype data from the transformed imputation gbs data for SNPs selected from our GBS tags;
my $sUsage = "perl $0 <our_snps_id file> <genotype files>\n";
die $sUsage unless @ARGV;

my ($id_file, @genotype_files) = @ARGV;

my %ids;
open (ID, $id_file) or die;
while(<ID>)
{
	chomp;
	my @t = split /\t/,$_;
	$ids{$t[1]}=1
}
close ID;

my %genotypes;
foreach my $file (@genotype_files)
{
	my $lines = 0;
	my %snp_index;
	open (IN, $file) or die;
	while (<IN>)
	{
		$lines++;
		chomp;
		my @t = split /,/, $_;
		if($lines == 1)
		{
			foreach my $index (0..$#t)
			{
				$snp_index{$index} = $t[$index] if exists $ids{$t[$index]};
			}
			next;
		}
		next if $lines == 2;
		foreach my $index ( 0.. $#t)
		{
			next unless exists $snp_index{$index};
			push @{$genotypes{$snp_index{$index}}}, $t[$index];
		}
	}
	close IN;
}

# Output
foreach my $snp (keys %genotypes)
{
	print join(",", ($snp, @{$genotypes{$snp}})), "\n";
}




