#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<allele frequency of 9k polysites in Excalibur>
<allele frequency of 9k polysites in Chara>
);
die $sUsage unless @ARGV ;
my ($ex_freq_file, $cha_freq_file) = @ARGV;
my %ex_freq = read_blat_freq_file($ex_freq_file);
my %cha_freq = read_blat_freq_file($cha_freq_file);

foreach my $cutoff (5..20)
{
	my $file = "cutoff_".$cutoff."_genotypes.out";
	open (OUT, ">$file") or die;
	print OUT join("\t",qw(ID Excalibur Chara) ),"\n";
	my %ex_genotypes = call_genotype_from_blat_results(\%ex_freq, $cutoff);
	my %cha_genotypes = call_genotype_from_blat_results(\%cha_freq, $cutoff);
	my ($total, $called) = (0, 0);
	foreach my $id (keys %ex_genotypes)
	{
		next unless exists $cha_genotypes{$id};
		print OUT join("\t", ($id, $ex_genotypes{$id}, $cha_genotypes{$id})),"\n";
		$total++; 
		$called++ unless $cha_genotypes{$id} eq $ex_genotypes{$id}
	}
	close OUT;
	my $call_rate  = sprintf("%.2f", $called/$total);
	print join("\t", ($cutoff, $total, $called, $call_rate)),"\n";
}


sub read_blat_freq_file
{
	my $file = shift;
	my $acc_name = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		# Excalibur_mira1_c57011:579	A	G	20	30
		# next unless /^$acc_name/;
		chomp;
		my @t = split /\t/, $_;
		my $id = $t[0];
		if(@t == 5)
		{			
			my %allele_count = @t[1,3,2,4];
			$return{$id} = {%allele_count};			
		}
		elsif(@t == 3)
		{
			$return{$id}{$t[1]} = $t[2];
		}
	}
	close IN;
	return %return;
}


sub call_genotype_from_blat_results
{
	my $snp_freq = shift;
	my $total_depth_cutoff = shift ;
	my %return;
	foreach my $id (keys %{$snp_freq})
	{
		my $genotype;
		my %allele_count = %{$snp_freq->{$id}};
		my @alleles = sort {$allele_count{$a} <=> $allele_count{$b}} keys %allele_count;
		if(@alleles == 1)
		{
			$genotype = $alleles[0] if $allele_count{$alleles[0]} >= $total_depth_cutoff;
		}
		else
		{
			my $total = $allele_count{$alleles[0]} + $allele_count{$alleles[1]};
			next unless $total >= $total_depth_cutoff;
			next if $total > 200;
			if($allele_count{$alleles[0]} >= 5 and $allele_count{$alleles[1]} >= 5) {$genotype = join("", sort{$a cmp $b}@alleles)}
			#if(abs($allele_count{$alleles[0]} -$allele_count{$alleles[1]})<=0.6*$total ) {$genotype = join("", sort{$a cmp $b}@alleles)}
			elsif($allele_count{$alleles[0]} == 0 and $allele_count{$alleles[1]} >= $total_depth_cutoff) {$genotype = $alleles[1]}
		}
		$return{$id} = $genotype if defined $genotype;
	}
	return %return;
}