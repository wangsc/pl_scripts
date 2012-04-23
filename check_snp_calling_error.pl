#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <Excalibur_Chara_genotypes_8_comparison.out>\n";
my $file = shift or die $sUsage;
my ($at_t, $t_at, $others) = (0, 0);
my $total = 0;
open (IN, $file) or die;
while(<IN>)
{
	chomp;
	my @t=split/\t/,$_;
	# Excalibur_mira1_c9752:289	T	C	T	CT
	# CAP11_mira1_c852:366	A	G	A	AG
	next unless $t[1] ne $t[2] and $t[3] eq $t[4];
	$total++;
	if($t[1] eq $t[3])
	{
		$at_t++ if length $t[2] == 2 and length $t[4] == 1;
		$t_at++ if length $t[4] == 2 and length $t[2] == 1;
	}
	elsif($t[1] eq $t[4])
	{
		$at_t++ if length $t[2] == 2 and length $t[3] == 1;
		$t_at++ if length $t[3] == 2 and length $t[2] == 1;
	}
	elsif($t[2] eq $t[4])
	{
		$at_t++ if length $t[1] == 2 and length $t[3] == 1;
		$t_at++ if length $t[3] == 2 and length $t[1] == 1;		
	}
	elsif($t[2] eq $t[3])
	{
		$at_t++ if length $t[1] == 2 and length $t[4] == 1;
		$t_at++ if length $t[4] == 2 and length $t[1] == 1;		
	}
	elsif(length $t[4] == 2 and length $t[1] == 1 and length $t[3] == 2 and length $t[2] == 1)
	{
		$t_at++;
	}
	else
	{
		$others++;
		print $_,"\n";
	}
}

print join("\t", qw(at_t t_at others total)),"\n";
print join("\t", ($at_t, $t_at, $others, $total)),"\n"