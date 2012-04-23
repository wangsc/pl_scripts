#!/usr/bin/perl -w
use strict;

# Given genotype for SNPs were called from reuslts of short flanking sequence of SNPs BLAT against Genomic data
# Comparing the previous genotype calling results based on: call AB if A>=3 and B>=3; call A if A>=N and B=0
# SW Nov.07.2011

my $sUsage = qq(
perl $0
<combined_all_snps_unique_freq>
<accession name>
<coverage depth N for call monomorphism>
<Allele frequency called from BLAT results>

);
die $sUsage unless @ARGV >= 4;
my ($all_snp_freq_file, $acc_name, $coverage_depth, $blat_snp_freq_file) = @ARGV;

my %all_snp_freq = read_freq_file($all_snp_freq_file, $acc_name);
my %blat_snp_freq = read_blat_freq_file($blat_snp_freq_file);
my %blat_snp_genotype = call_genotype_from_blat_results(\%blat_snp_freq);

foreach my $dp ($coverage_depth .. 20)
{
	my %snp_genotype = call_genotype(\%all_snp_freq, $dp);
	my $output = $acc_name . "_depth_". $dp . "_compare_to_blat_genotype.out";
	open (OUT, ">$output") or die;
	my ($both_called, $both_called_same) = (0, 0);
	foreach my $id (keys %snp_genotype)
	{
		next unless exists $blat_snp_genotype{$id};
		$both_called++;
		$both_called_same++ if $snp_genotype{$id} eq $blat_snp_genotype{$id};
		print OUT join("\t", ($id, $snp_genotype{$id}, $blat_snp_genotype{$id})),"\n";
	}
	close OUT;
	print join(" ", ($dp, $both_called, $both_called_same)),"\t", sprintf("%.2f", $both_called_same/$both_called),"\n";
}

sub read_freq_file
{
	my $file = shift;
	my $acc_name = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		# Excalibur	Excalibur_mira1_c57011	579	A	G	20	30	NA
		next unless /^$acc_name/;
		chomp;
		my @t = split /\t/, $_;
		my $id = join(":", @t[1,2]);
		my %allele_count = @t[3,5,4,6];
		$return{$id} = {%allele_count};
	}
	close IN;
	return %return;
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
	my $total_depth_cutoff = 8;
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
			#next unless $total >= $total_depth_cutoff;
			#if($allele_count{$alleles[0]} >= 2 and $allele_count{$alleles[1]} >= 2) {$genotype = join("", sort{$a cmp $b}@alleles)}
			if(abs($allele_count{$alleles[0]} -$allele_count{$alleles[1]})<=0.6*$total ) {$genotype = join("", sort{$a cmp $b}@alleles)}
			elsif($allele_count{$alleles[0]} == 0 and $allele_count{$alleles[1]} >= $total_depth_cutoff) {$genotype = $alleles[1]}
		}
		$return{$id} = $genotype if defined $genotype;
	}
	return %return;
}

sub call_genotype
{
	my ($snp_freq, $dp) = @_;
	my %return;
	foreach my $id (keys %{$snp_freq})
	{
		my $genotype;
		my %allele_count = %{$snp_freq->{$id}};
		my @alleles = sort {$allele_count{$a} <=> $allele_count{$b}} keys %allele_count;
		if(@alleles == 1)
		{
			$genotype = $alleles[0] if $allele_count{$alleles[0]} >= $dp;
		}
		else
		{
			if($allele_count{$alleles[0]} >= 3 and $allele_count{$alleles[1]} >= 3) {$genotype = join("", sort{$a cmp $b}@alleles)}
			elsif($allele_count{$alleles[0]} == 0 and $allele_count{$alleles[1]} >= $dp) {$genotype = $alleles[1]}
		}
		$return{$id} = $genotype if defined $genotype;
	}
	return %return;
}





