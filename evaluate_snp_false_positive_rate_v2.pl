#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<combined_all_snps_unique_freq_labeled>
<accession 1 name>
<accession 2 name>
<accession 1 allele frequency from BLAT>
<accession 2 allele frequency from BLAT>
);

die $sUsage unless @ARGV >= 5;
my ($all_freq_file, $acc1, $acc2, $acc1_freq_file, $acc2_freq_file) = @ARGV;

my %acc1_freq = read_all_freq_file($all_freq_file, $acc1);
my %acc2_freq = read_all_freq_file($all_freq_file, $acc2);

my %acc1_genotype_blat = call_genotype_from_blat_freq_file($acc1_freq_file);
my %acc2_genotype_blat = call_genotype_from_blat_freq_file($acc2_freq_file);

print join("\t", qw(Depth DC SC DR Total)),"\n";
foreach my $d (8..20)
{
	my $output = join("_", ($acc1, $acc2,"genotypes_")) .$d . "_comparison.out";
	open (OUT, ">$output") or die;
	my %acc1_genotype = call_genotype(\%acc1_freq, $d);
	my %acc2_genotype = call_genotype(\%acc2_freq, $d);
	my($dc, $sc, $dr) = (0,0,0,0);
	foreach (keys %acc1_genotype)
	{
		next unless exists $acc1_genotype{$_} and exists $acc2_genotype{$_} and 
								exists $acc1_genotype_blat{$_} and exists $acc2_genotype_blat{$_};
		my $our_snp = $acc1_genotype{$_} eq $acc2_genotype{$_}?0:1;
		my $blat_snp = $acc1_genotype_blat{$_} eq $acc2_genotype_blat{$_}?0:1;
		if($acc1_genotype{$_} eq $acc1_genotype_blat{$_})
		{
			$dc++ if $acc2_genotype{$_} eq $acc2_genotype_blat{$_};
			$sc++ if $acc2_genotype{$_} ne $acc2_genotype_blat{$_}
		}
		else
		{
			$sc++ if $acc2_genotype{$_} eq $acc2_genotype_blat{$_};
			$dr++ if $acc2_genotype{$_} ne $acc2_genotype_blat{$_}
		}

		print OUT 
		join("\t",($_, $acc1_genotype{$_}, $acc2_genotype{$_}, $acc1_genotype_blat{$_}, $acc2_genotype_blat{$_})), "\n";
	}
	close OUT;
	print join("\t", ($d, $dc, $sc, $dr,sum($dc, $sc, $dr))),"\n";

}

sub sum
{
	my $sum = 0;
	map{$sum += $_}@_;
	return $sum;
}

sub read_all_freq_file
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

sub call_genotype_from_blat_freq_file
{
	my $file = shift;
	my %blat_freq = read_blat_freq_file($file);
	my %genotypes = call_genotype_from_blat_results(\%blat_freq);
	return %genotypes;
}

sub read_blat_freq_file
{
	my $file = shift;
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
			next unless $total >= $total_depth_cutoff;
			next if $total >= 200;
			#if($allele_count{$alleles[0]} >= 2 and $allele_count{$alleles[1]} >= 2) {$genotype = join("", sort{$a cmp $b}@alleles)}
			if(abs($allele_count{$alleles[0]} -$allele_count{$alleles[1]})<=0.6*$total ) {$genotype = join("", sort{$a cmp $b}@alleles)}
			elsif($allele_count{$alleles[0]} == 0 and $allele_count{$alleles[1]} >= $total_depth_cutoff) {$genotype = $alleles[1]}
		}
		$return{$id} = $genotype if defined $genotype;
	}
	return %return;
}