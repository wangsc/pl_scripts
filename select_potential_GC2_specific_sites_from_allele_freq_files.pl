#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <CS_allele_freq.out> <GC2_allele_freq.out> <coverage_depth_cutoff, default 15>\n";
die $sUsage unless @ARGV >= 3;

my ($cs_freq_file, $gc2_freq_file, $coverage_depth_cutoff) = @ARGV;
$coverage_depth_cutoff = 15 unless defined $coverage_depth_cutoff; # Min depth to call A or B genotype

my %cs_genotype = read_freq_file($cs_freq_file);
my %gc2_genotype = read_freq_file($gc2_freq_file);

print join("\t",qq(Contig SNP_pos GC2_genotype CS_genotype)),"\n";
foreach (keys %cs_genotype)
{
	next unless exists $gc2_genotype{$_};
	#print join("\t",($_, $gc2_genotype{$_}, $cs_genotype{$_})),"\n" unless $gc2_genotype{$_} eq $cs_genotype{$_};
	print join("\t",($_, $gc2_genotype{$_}, $cs_genotype{$_})),"\n" unless $gc2_genotype{$_} eq $cs_genotype{$_};
}


sub read_freq_file
{
	my $file = shift;
	my %return;
	open (IN, "$file") or die;
	while(<IN>)
	{
		chomp;
		my @t = split /\t/, $_;
		next if @t >= 5;
		#comp11721_c1_seq1	294	T_4	C_4
		#comp11721_c1_seq1	349	T_2
		#comp11721_c1_seq1	2246	G_9	A_1
		my $genotype = call_genotype(@t[2..$#t]);
		$return{join("\t", @t[0,1])} = $genotype unless $genotype eq "NA";		
	}
	close IN;
	return %return;
}


sub call_genotype
{
	my @arr = @_;
	my (@alleles, @count);
	foreach (@arr)
	{
		my @tmp = split /_/, $_;
		push @alleles, $tmp[0];
		push @count, $tmp[1];
	}
	if(@count==1)
	{
		return $alleles[0] if $count[0] >= $coverage_depth_cutoff
	}
	else
	{
		if($count[0]>=3 and $count[1]>=3){return join("", sort{$a cmp $b}@alleles)}
		else{return $alleles[1] if $count[1] >= $coverage_depth_cutoff and $count[0] == 0}
	}
	return "NA";
}