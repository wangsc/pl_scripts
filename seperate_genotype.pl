#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<SNPmap.2110_final_MAGIC.csv>
<combinedPM_mod_final_genotype_mod4.csv>
);

die $sUsage unless @ARGV;

my ($map_file, $genotype_file) = @ARGV;

my %map_ids = read_map_file($map_file);
my ($snp_index_hashref, $genotype_arrayref, $groups) = read_genotype_file($genotype_file);

# OUTPUT
my @chr_ids = map{$_."A", $_."B", $_."D"}(1..7);

my %dist_files;
foreach my $g (@$groups)
{	
	foreach my $chr (@chr_ids)
	{
		my $outfile = $g . "_" . $chr . "_genotype.out";
		open (OUT, ">$outfile") or die "can't open file $outfile\n";
		my @snpids = @{$map_ids{$chr}->[0]};
		my @index;
		my @snpid_index;
		foreach (0.. $#snpids)
		{
			my $snpid = $snpids[$_];
			next if(not exists $snp_index_hashref->{$snpid});
			push @index, $snp_index_hashref->{$snpid};
			push @snpid_index, $_;
		}
		print OUT join("\t", @snpids[@snpid_index]),"\n";
		foreach (@$genotype_arrayref)
		{
			my @data = @$_;
			if($g=~/landrace/i)
			{
				next unless $data[2] =~ /landrace/i
			}
			else
			{
				next unless $g eq join("_", @data[1,2])
			}			
			print OUT join("\t", @data[@index]),"\n"
		}
		close OUT;
		
		unless (exists $dist_files{$chr})
		{
			$dist_files{$chr} = 1;
			my $dist = "sub_" . $chr ."_dist.out";
			open (OUT, ">$dist") or die;
			my %remain_snps = map{$_,1} @snpids[@snpid_index];
			foreach (@{$map_ids{$chr}->[1]})
			{
				my ($id, $d) = split/\t/,$_;
				print OUT $_,"\n" if exists $remain_snps{$id};
			}
		}
		close OUT;	
	}
}

# Subroutines

sub read_map_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		# "1A35","wsnp_Ex_c19353_28290651","1A",38.1228643055790
		chomp;
		s/\"//g;
		my @t = split /,/, $_;
		push @{$return{$t[2]}->[0]}, $t[1];
		push @{$return{$t[2]}->[1]}, join("\t", ($t[1], $t[3]));
	}
	close IN;
	return %return;
}


sub read_genotype_file
{
	my $file = shift;
	open (IN, $file) or die;
	my $count = 0;
	my (%index_hash, @genotypes, %group_count);
	while(<IN>)
	{
		$count++;
		chomp;
		my @t = split /,/, $_;
		if($count==1)
		{
			foreach (3..$#t)
			{
				$index_hash{$t[$_]} = $_;
			}
			next;
		}
		if($t[2]=~/landrace/i)
		{
			$group_count{$t[2]}++
		}
		else
		{
			$group_count{join("_", @t[1,2])}++
		}
		push @genotypes, [@t];		
	}	
	close IN;
	return (\%index_hash, \@genotypes, [keys %group_count])
}



















