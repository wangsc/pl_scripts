#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<SNPmap.2110_final_MAGIC.csv>
<combinedPM_mod_final_genotype_mod4.csv>
<gruping file>
);

die $sUsage unless @ARGV;

my ($map_file, $genotype_file, $grp_file) = @ARGV;

my %grps = read_grp_file($grp_file) if $grp_file;
my %map_ids = read_map_file($map_file);
my ($snp_index_hashref, $genotype_arrayref, $groups) = read_genotype_file($genotype_file);

# OUTPUT
#my @chr_ids = map{$_."A", $_."B", $_."D"}(1..7);
my @chr_ids = unique(keys %map_ids);
my %dist_files;
my @grps = qw(cultivar landraces);
foreach my $g (@grps)
{	
	foreach my $chr (@chr_ids)
	{
		my $outfile = $g . "_" . $chr . "_genotype.out";
		open (OUT, ">$outfile") or die "can't open file $outfile\n";
		next unless exists $map_ids{$chr};
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
			next unless exists $grps{$g}{$data[0]};	
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
sub read_grp_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		# acc,Index_MH,lines,status,country,state,groups_4,growth_habit,groups_5
		# Atlas66,1,G1,cultivar,USA,KS,MIDWEST_winter,winter,MIDWEST_winter
		next if /^acc/;
		chomp; 
		my @t = split /,/,$_;
		$return{lc($t[3])}{$t[2]}=1;
	}
	close IN;
	return %return;
}


sub read_map_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		# wsnp_Ku_c183_358844,17.31,1A,1,A,11.58846353,1
		chomp;
		s/\"//g;
		next if /^Marker/i;
		my @t = split /,/, $_;
		push @{$return{$t[2]}->[0]}, $t[0];
		push @{$return{$t[2]}->[1]}, join("\t", ($t[0], $t[5]));
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
		# lines,groups_4,growth_habit,wsnp_AJ612027A_Ta_2_1,wsnp_AJ612027A_Ta_2_5
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
			$group_count{$t[1]}++
		}
		push @genotypes, [@t];		
	}	
	close IN;
	return (\%index_hash, \@genotypes, [keys %group_count])
}

sub unique
{
	my %hash = map{$_, 1} @_;
	return keys %hash;
}















