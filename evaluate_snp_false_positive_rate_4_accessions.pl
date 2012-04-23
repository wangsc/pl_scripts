#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<combined_all_snps_unique_freq_labeled>
<accession names, seperated by ,>
<allele frequency files from BLAT fore each accession>
);

die $sUsage unless @ARGV >= 3;
my ($all_freq_file, $acc_names, $accs_freq_files) = @ARGV;
my ($acc1, $acc2, $acc3, $acc4) = split /,/, $acc_names;
my ($acc1_freq_file, $acc2_freq_file, $acc3_freq_file, $acc4_freq_file) = split /,/, $accs_freq_files;
my %acc1_freq = read_all_freq_file($all_freq_file, $acc1);
my %acc2_freq = read_all_freq_file($all_freq_file, $acc2);
my %acc3_freq = read_all_freq_file($all_freq_file, $acc3);
my %acc4_freq = read_all_freq_file($all_freq_file, $acc4);

my %acc1_genotype_blat = call_genotype_from_blat_freq_file($acc1_freq_file);
my %acc2_genotype_blat = call_genotype_from_blat_freq_file($acc2_freq_file);
my %acc3_genotype_blat = call_genotype_from_blat_freq_file($acc3_freq_file);
my %acc4_genotype_blat = call_genotype_from_blat_freq_file($acc4_freq_file);

print join("\t", qw(Depth TN TP FN FP Total Correct_rate FP_rate)),"\n";
foreach my $d (8..20)
{
	my $output = join("_", ($acc1, $acc2, $acc3, $acc4,"genotypes_")) .$d . "_comparison.out";
	open (OUT, ">$output") or die;
	my %acc1_genotype = call_genotype(\%acc1_freq, $d);
	my %acc2_genotype = call_genotype(\%acc2_freq, $d);
	my %acc3_genotype = call_genotype(\%acc3_freq, $d);
	my %acc4_genotype = call_genotype(\%acc4_freq, $d);
	
	my($tn, $tp, $fn, $fp) = (0,0,0,0);
	foreach my $id (keys %acc1_genotype)
	{
		next unless genotyped_all_accessions($id, 4, \%acc1_genotype, \%acc2_genotype, \%acc3_genotype, \%acc4_genotype) and 
								genotyped_all_accessions($id, 2, \%acc1_genotype_blat, \%acc2_genotype_blat, \%acc3_genotype_blat, \%acc4_genotype_blat);
		my $our_snp = check_if_is_snp($id, \%acc1_genotype, \%acc2_genotype, \%acc3_genotype, \%acc4_genotype);
		my $blat_snp = check_if_is_blat_snp($id, \%acc1_genotype_blat, \%acc2_genotype_blat, \%acc3_genotype_blat, \%acc4_genotype_blat);
		
		if($our_snp == 1)
		{
			$tp++ if $blat_snp == 1;
			$fp++ if $blat_snp == 0;
		}
		else
		{
			$tn++ if $blat_snp == 0;
			$fn++ if $blat_snp == 1;
		}
		my @temp = map {exists $_->{$id}?$_->{$id}:"NA"} (\%acc1_genotype, \%acc2_genotype, \%acc3_genotype, \%acc4_genotype, \%acc1_genotype_blat, \%acc2_genotype_blat, \%acc3_genotype_blat, \%acc4_genotype_blat);
		print OUT join("\t",($id, @temp)), "\n";
	}
	close OUT;

	print join("\t", ($d, $tn, $tp, $fn, $fp, sum($tn, $tp, $fn, $fp), 
						sprintf("%.2f", sum($tn+$tp)/sum($tn, $tp, $fn, $fp)), 
						sprintf("%.2f", $fp/sum($tp, $fp)) ) ),"\n";
}

sub check_if_is_snp
{
	my $id = shift;
	my %genotypes;
	map{$genotypes{$_->{$id}}++} @_;
	my @types = keys %genotypes;
	return 1 if @types == 2 and $genotypes{$types[0]} == 2;
	return 0;
}

sub check_if_is_blat_snp
{
	my $id = shift;
	my %genotypes;
	map{$genotypes{$_->{$id}}++ if exists $_->{$id}} @_;
	my @types = keys %genotypes;
	return 1 if @types == 2;
	return 0;
}


sub genotyped_all_accessions
{
	my $id = shift;
	my $min = shift;
	my $count = 0; 
	foreach (@_)
	{
		$count++ if exists $_->{$id}
	}
	return $count>=$min?1:0;
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
		next unless /^$acc_name/i;
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
			my $total = $allele_count{$alleles[0]} + $allele_count{$alleles[1]};
			next unless $total >= $dp;
			#if($allele_count{$alleles[0]} >= 1 and $allele_count{$alleles[1]} >= 1) {$genotype = join("", sort{$a cmp $b}@alleles)}
			if(abs($allele_count{$alleles[0]} -$allele_count{$alleles[1]})<=0.6*$total ) {$genotype = join("", sort{$a cmp $b}@alleles)}
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