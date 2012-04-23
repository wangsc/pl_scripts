#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
This script will randomly select N high quality SNPs which show different genotypes in two accessions for validation. 

perl $0
<high qulaity snps>
<allele freq file>
<SNP exon-intron junction indication file>
<accession A>
<accession B>
<number of SNPs to be selected>
<output file>
);
die $sUsage unless @ARGV >= 7;

my($hq_snp_file, $allele_freq_file, $juc_file, $acc_a, $acc_b, $num_snp, $output_file) = @ARGV;

my %hq_snps = read_snp_file($hq_snp_file);
my %snp_junc = read_juc_file($juc_file); 
my %snps_between_ab = get_snps_between_ab($allele_freq_file, $acc_a, $acc_b, \%snp_junc);

my @selected_snps = random_select_snps(\%hq_snps, \%snps_between_ab, $num_snp);
map{ $selected_snps[$_] =~ s/:/\t/ }0..$#selected_snps;
# Output
open (OUT, ">$output_file") or die;
print OUT join("\n", @selected_snps), "\n";
close OUT;

sub read_juc_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		chomp; 
		my @t= split /\s+/,$_;
		my $id = $1 if $t[0] =~ /(\S+:\S+):/;
		$return{$id} = $t[1];
	}
	close IN;
	return %return;
}


sub random_select_snps
{
	my ($hq_snp, $ab_snp, $num) = @_;
	my @hq_ab = grep{exists $ab_snp->{$_}}keys %$hq_snp;
	if(@hq_ab < $num){die "!! Not enough SNPs for selection \n", scalar @hq_ab, "\n"}
	my @random_index = map{int(rand(scalar @hq_ab))}1..$num;
	return @hq_ab[@random_index];
}

sub read_snp_file
{
	my $file  = shift;
	my %return;
	open (IN, $file) or die;
	while (<IN>)
	{
		chomp; 
		next if /^\s+$/;
		$return{join(":", (split /\s+/,$_))} = 1
	}
	close IN;
	return %return;
}

sub get_snps_between_ab
{
	my ($file, $acca, $accb, $snp_junc) = @_;
	my %return;
	my %record;
	open (IN, $file) or die;
	while(<IN>)
	{
		next unless /^$acca|$accb/;
		chomp;
		my @t = split /\s+/,$_;
		if (@t == 6)
		{
			if($t[4]>=3 and $t[5]>=3)
			{
				$record{$t[1]}{$t[0]} = join("", sort{$a cmp $b}@t[2,3]);
			}
		}
		else
		{
			$record{$t[1]}{$t[0]} = $t[2] if $t[3] >= 8;
		}
	}
	close IN;
	
	foreach my $snp (keys %record)
	{
		next unless exists $snp_junc->{$snp};
		next unless $snp_junc->{$snp} == 2;
		if(exists $record{$snp}{$acca} and exists $record{$snp}{$accb})
		{
			$return{$snp} = 1 unless  $record{$snp}{$acca} eq  $record{$snp}{$accb} 
		}
	}
	
	return %return;
}



