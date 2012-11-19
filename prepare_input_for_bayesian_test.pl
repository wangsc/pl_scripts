#!/usr/bin/perl -w
use strict;
use Statistics::ChisqIndep;

my $sUsage = qq(
perl $0
<A_B_D_allele.out>
<43k_clusters.list>
<Grandin_sorted_allele_freq.out>
);
die $sUsage unless @ARGV >= 3;

my($abd_allele_file, $cluster_file, $freq_file, $dp_cutoff, $pvalue_outfile) = @ARGV;
#open (PO, ">$pvalue_outfile") or die "can't open file $pvalue_outfile\n";
#open (CO, ">$count_outfile") or die "can't open file $count_outfile\n";

my $chisq_test = new Statistics::ChisqIndep;
my ($genome_spec_ids, $id_alleles) = read_allele_file($abd_allele_file);
my ($cluster_ids, $id_to_cluster) = read_cluster_file($cluster_file);
my %count_allele_each_cluster = count_allele_freq($freq_file, $genome_spec_ids, $id_alleles, $cluster_ids, $id_to_cluster, $dp_cutoff);

# subroutines

sub sum
{
	my $return = 0;
	map{$return += $_}@_;
	return $return;
}

sub get_file_handle
{
	my ($file, @arr) = @_;
	my %return_hash;
	foreach (0..$#arr)
	{
		local *FH;
		my $out = $file . "_". $arr[$_];
		open (FH, ">$out") or die "can't open file $out\n";
		print FH join("\t", qw(geneNameColumn snpIndexColumn alleleOneCount alleleTwoCount)), "\n";
		$return_hash{$_} = *FH{IO};
	}
	return %return_hash;
}

sub count_allele_freq
{
	my ($file, $genome_spec_ids, $id_alleles, $cluster_ids, $id_to_cluster, $dp_cutoff) = @_;
	my %count_allele_each_cluster; # $count_allele_each_cluster{$cluster_id}{$genome} = [#A, #B];
	my @subgenomes = qw(A B D);
	my %out_fh = get_file_handle($file, @subgenomes);
	open (F, "$file") or die "can't open file $file\n";
	while(<F>)
	{
		# Kukri_mira1_c5005       806     G_104   C_53 N_1
		my @t = split /\s+/, $_;
		next unless exists $id_to_cluster->{$t[0]};
		my %cont = map{split /_/, $_} @t[2..$#t];
		my $cluster = $id_to_cluster->{$t[0]};
		my $id = join(":", @t[0,1]);
		next unless exists $id_alleles->{$id};
		my $sub_genome = $genome_spec_ids->{$id}; # 0->A; 1->B; 2->D
		my @alleles = @{$id_alleles->{$id}};
		my $total = 0;
		map{$total += exists $cont{$_}?$cont{$_}:0}@alleles;
		my @num_alleles = map{exists $cont{$_}?$cont{$_}:0}@alleles;
		
		print {$out_fh{$sub_genome}} join("\t", ($cluster, $id, @num_alleles)), "\n"; 
	}
	close F;
	return %count_allele_each_cluster;
}


sub read_allele_file
{
	my $file = shift;
	my %genome_specific_ids;
	my %id_alleles;
	open (IN, $file) or die "can't open file $file\n";
	while(<IN>)
	{
		# Excalibur_mira1_c62758  258     A       T       T
		chomp; 
		my @t = split /\s+/, $_;
		my $id = join(":", @t[0, 1]);
		my @arr = @t[2..4];
		my %allele_cnt;
		map{$allele_cnt{$_} = 0 unless exists $allele_cnt{$_}; $allele_cnt{$_}++}@arr;
		#print STDERR $_, "\n" if keys %allele_cnt > 2;
		next if keys %allele_cnt > 2;
		foreach my $index (0..$#arr)
		{
			if ($allele_cnt{$arr[$index]} == 1)
			{
				$genome_specific_ids{$id} = $index;
				last;
			}			
		}
		
		$id_alleles{$id} = [sort{$allele_cnt{$a}<=>$allele_cnt{$b}}keys %allele_cnt];
	}
	close IN;
	return (\%genome_specific_ids, \%id_alleles)
}

sub read_cluster_file
{
	my $file = shift;
	my (%cluster_ids, %id_cluster);
	open (IN, $file) or die "can't open file $file\n";
	while(<IN>)
	{
		# RAC875_mira1_rep_c85909 RAC875_mira1_c18853:RAC875_mira1_rep_c85909
		# Kukri_mira1_c35572      Kukri_mira1_c35572
		chomp; 
		my @t = split /\s+/, $_;
		my @ids = split /:/, $t[1];
		$cluster_ids{$t[0]} = $t[1];
		map{$id_cluster{$_} = $t[0]}@ids;
	}
	close IN;
	return (\%cluster_ids, \%id_cluster);
}



