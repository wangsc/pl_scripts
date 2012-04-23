#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <position file> <genotype file>\n";
die $sUsage unless @ARGV >=2;
my ($pos_file, $genotype_file) = @ARGV;

my (%chr_start_end, %chr_markers);
open(P, $pos_file) or die "can't open file $pos_file\n";
my $count = -1;
while(<P>)
{
	# SNP6458	67.63008433	1A	wsnp_Ku_c12441_20092792	0.052863436
	chomp;
	next if /^\s+$/;
	my @t = split /\t/,$_;
	$count++;
	push @{$chr_start_end{$t[2]}}, $count;
	push @{$chr_markers{$t[2]}}, [@t[0,1]];
}
close P;

open (G, $genotype_file) or die "can't open $genotype_file\n";
my @genotype;
while(<G>)
{
	#print $_;
	chomp;
	my @t=split/\t/,$_;
	#print scalar @t,"\n";
	@t= @t[6..$#t];
	#print scalar @t,"\n";
	#print join("**", @t),"\n";
	foreach (0..$#t)
	{
		$t[$_] =~ s/1/A/g;
		$t[$_] =~ s/2/T/g;
		$t[$_] =~ s/ /\//;
		$t[$_] = "NA" if $t[$_]=~/0/;
	}
	#print join("\t", @t),"\n";
	push @genotype, [@t];
}
close G;

# OUTPUT

foreach my $chr (keys %chr_start_end)
{
	my $dist_file = $chr . '_dist.out';
	my $genotype_out = $chr . "_genotype.out";
	open (D, ">$dist_file") or die;
	open (G, ">$genotype_out") or die;
	my @markers;
	foreach (@{$chr_markers{$chr}})
	{
		print D join("\t", @$_),"\n";
		push @markers, $_->[0];
	}
	print G join("\t", @markers),"\n";
	foreach (@genotype)
	{
		my @array = @$_;
		print G join("\t", @array[@{$chr_start_end{$chr}}]),"\n";
	}
	close D;
	close G;	
}
