#!/usr/bin/perl -w
use strict;
use Statistics::ChisqIndep;

my $sUsage = qq(
perl $0
<A_B_D_allele.out>
<43k_clusters.list>
<Grandin_sorted_allele_freq.out>
<min_coverage depth for chiseq test>
<pvalue output file>
);
die $sUsage unless @ARGV >= 4;

my($abd_allele_file, $cluster_file, $freq_file, $dp_cutoff, $pvalue_outfile) = @ARGV;
open (PO, ">$pvalue_outfile") or die "can't open file $pvalue_outfile\n";
#open (CO, ">$count_outfile") or die "can't open file $count_outfile\n";

my $chisq_test = new Statistics::ChisqIndep;
my ($genome_spec_ids, $id_alleles) = read_allele_file($abd_allele_file);
my ($cluster_ids, $id_to_cluster) = read_cluster_file($cluster_file);
my %count_allele_each_cluster = count_allele_freq($freq_file, $genome_spec_ids, $id_alleles, $cluster_ids, $id_to_cluster, $dp_cutoff);

&run_chisq_test(\%count_allele_each_cluster, $chisq_test);

# subroutines
sub run_chisq_test
{
	my $hashref = shift;
	my $test_handle = shift;
	
	foreach my $cluster (keys %$hashref)
	{
		my @pvalues = qw(NA NA NA);
		my @cnt = qw(NA NA NA);
		my @sub_genome = sort{$a<=>$b} keys %{$hashref->{$cluster}};
		foreach my $sub (@sub_genome)
		{
			my @allele_cnt = @{$hashref->{$cluster}{$sub}}[0, 1];
			map{$allele_cnt[$_] = int($allele_cnt[$_])}0..$#allele_cnt;
			#print STDERR  $cluster, "\n" if @allele_cnt > 2;
			$cnt[$sub] = join("/", @allele_cnt);
			my $total = sum(@allele_cnt);
			$chisq_test->load_data([[@allele_cnt], [int($total/3), int($total*2/3)]]);# 1:2
			my $pvalue = $chisq_test->p_value;
			$pvalue *= -1 if $allele_cnt[0]<($total/3) and $pvalue < 0.05;
			$pvalues[$sub] = $pvalue;
		}
		my @effect = check_effect(@pvalues, @cnt);
		print PO join("\t", ($cluster, @pvalues, @cnt, @effect)), "\n";
		#print CO join("\t", ($cluster, @cnt)), "\n";
	}	
}
close PO;
#close CO;

# subroutines
sub check_effect
{
	my @arr = @_;
	my @gn = qw(A B D);
	my @return;
	foreach my $ind (0..2)
	{
		my $p = $arr[$ind];
		if($p eq "NA"){$return[$ind] = $gn[$ind] . "_NA"; next}
		if ($p >0.05)
		{
			$return[$ind] = $gn[$ind] . "_BAL";
		}
		else
		{
			if($p > 0)
			{
				$return[$ind] = $gn[$ind] . "_UP";
			}
			else
			{
				my ($ca, $cb) = split /\//, $arr[$ind+3];
				if($ca < $cb*0.1)
				{
					$return[$ind] = $gn[$ind] . "_SILENCE";
				}
				else
				{
					$return[$ind] = $gn[$ind] . "_DOWN";
				}
			}
		}		
	}
	return @return;
}


sub sum
{
	my $return = 0;
	map{$return += $_}@_;
	return $return;
}


sub count_allele_freq
{
	my ($file, $genome_spec_ids, $id_alleles, $cluster_ids, $id_to_cluster, $dp_cutoff) = @_;
	my %count_allele_each_cluster; # $count_allele_each_cluster{$cluster_id}{$genome} = [#A, #B];
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
		next if $total < $dp_cutoff;
		$count_allele_each_cluster{$cluster}{$sub_genome} = [0, 0, 0] unless exists $count_allele_each_cluster{$cluster}{$sub_genome};
		#map{$count_allele_each_cluster{$cluster}{$sub_genome}[$_] += exists $cont{$alleles[$_]}?$cont{$alleles[$_]}:0} 0 .. $#alleles;
		$count_allele_each_cluster{$cluster}{$sub_genome}[2]++;
		foreach (0 .. $#alleles)
		{
			my $add = exists $cont{$alleles[$_]}?$cont{$alleles[$_]}:0;
			$count_allele_each_cluster{$cluster}{$sub_genome}[$_] = 
			($count_allele_each_cluster{$cluster}{$sub_genome}[$_] * ($count_allele_each_cluster{$cluster}{$sub_genome}[2] - 1) + $add)/$count_allele_each_cluster{$cluster}{$sub_genome}[2]
		}
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



