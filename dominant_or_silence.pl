#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 *_rmdup_allele_freq.out_pvalue2_adjusted_addeffect\n";
die $sUsage unless @ARGV;

my $infile = shift;
open (IN, $infile) or die "can't open file $infile\n";

while(<IN>)
{
	chomp;
	next if /^cluster/;
	my @t = split /,/, $_;
	my $id = $t[0];
	my @abd_eff = @t[-3, -2, -1];
	my @silence;
	my @dominant;
	my @balance;
	map{push @silence, $1 if /(\S+)_SILENCE/} @abd_eff;
	map{push @dominant, $1 if /(\S+)_UP/} @abd_eff;
	map{push @dominant, $1 if /(\S+)_BAL/} @abd_eff;
	if(@silence == 2)
	{
		$dominant[0] = get_the_else([qw(A B D)], [@silence])
	}
	
	if(@dominant == 3)
	{
		$dominant[0] = "NA";
	}
	
	my $s = @silence>=1?join("", sort{$a cmp $b}@silence):"NA";
	my $d = (@dominant>=1 and $dominant[0] ne "NA")?join("", sort{$a cmp $b}@dominant):"NA";
	
	print join("\t", ($id, $s, $d)), "\n";
}

sub get_the_else
{
	my $arr1 = shift;
	my $arr2 = shift;
	my %hash = map {$_, 1} @$arr2;
	foreach (@$arr1)
	{
		return $_ unless exists $hash{$_}
	}
}