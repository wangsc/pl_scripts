#!/usr/bin/perl -w
use strict;

# extract genotype data for each chromosome for MACH

my $sUsage = "perl $0 <genotype_gbs_snp_combined.out> <1A_genotyped_snps> <output file>\n";
die $sUsage unless @ARGV >= 3;
my ($genotype_file, $snp_file, $output) = @ARGV;

my %genotyped_snps;
my %geno = read_genotype_file($genotype_file);
my @snps = read_snp_file($snp_file);

# output in pedigree format
open (OUT, ">$output") or die $!;
foreach my $acc (keys %geno)
{
	my @genotypes = map{exists $geno{$acc}{$_}?$geno{$acc}{$_}:'0/0'}@snps;
	print OUT join(" ", ($acc, $acc, (0,0,2), @genotypes)), "\n";
}
close OUT;


# Subroutines
sub read_snp_file
{
	my $file = shift;
	open (IN, $file) or die;
	my @return;
	while(<IN>)
	{
		chomp; 
		push @return, $_; 
	}	
	close IN;
	return @return;
}


sub read_genotype_file
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	my %acc_index;
	while (<IN>)
	{
		chomp; 
		my @t = split /\t/,$_;
		if(/^SNP/)
		{
			%acc_index = map{$_, $t[$_] }(1..$#t);
			next
		}
		if(/^wsnp/)
		{
			map{$t[$_]="NA" if !($t[$_]=~/NA/) and $t[$_]==0}1..$#t
		}
		my @alleles = unique(@t[1..$#t]);
		@alleles = ("NA", @alleles);
		my %allele_value = map{$alleles[$_], $_}0..$#alleles;
		foreach (1..$#t)
		{
			$return{$acc_index{$_}}{$t[0]} = join("/", ($allele_value{$t[$_]}, $allele_value{$t[$_]}));
		}
	}
	close IN;
	return %return;
}

sub unique
{
	my %t;
	foreach (@_)
	{
		next if /NA/;
		$t{$_}=1
	};
	return keys %t;
}