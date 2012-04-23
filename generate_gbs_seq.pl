#!/usr/bin/perl -w
use strict;
my $file = shift or die "perl $0 <gbs tag file>\n";
open (IN, $file) or die;
my $n++;
while (<IN>)
{
	# J_tag1,TGCAG[A/C]TTATCAAGGGCTTATTTCCACCGCAGCTTATGTAGAAGCCATGGCTCAAAAGAACTGG,1A,0
	chomp; 
	my @t= split /,/,$_;
	my $seq = $t[1];
	my @comb = generate_seq($seq);
	foreach (0..$#comb)
	{
		print ">",$t[0],"_", $_+1, "\n";
		print $comb[$_],"\n";
	}
}



sub generate_seq
{
	# TGCAGAAAA[A/T]AAAAATGAAACGAAAGGGACAGTCTACCTCTCCCTTCCCTT[T/A]AAAAACAGTTTG
	my $seq = shift;
	my @return;
	my @pos;
	while ($seq =~ /\[/g)
	{
		push @pos, pos($seq)
	}
	
	my @snps;
	foreach (@pos)
	{
		my $str = substr($seq, $_-1, 5);
		my @alleles = $str=~/\[(\S)\/(\S)\]/;
		#print scalar @alleles, "\n";
		push @snps, [@alleles];
	}
	#print scalar @snps, "\n";
	$seq =~ s/\[\S{3}\]/X/g;
	@pos = ();
	while ($seq =~ /X/g)
	{
		push @pos, pos($seq)
	}
	
	my @comb = all_combination(@snps);
	
	foreach (@comb)
	{
		my @tmp = split //,$_;
		foreach (0..$#tmp)
		{
			substr($seq, $pos[$_]-1, 1) = $tmp[$_]
		}
		push @return, $seq;
	}
	return @return;
}

sub all_combination
{
	my @arr = @_;
	#print scalar @arr, "\n";
	if(@arr == 1)
	{
		return @{$arr[0]};
	}
	else
	{
		my $crr = pop @arr;
		my @pre_comb = all_combination(@arr) if @arr > 0;
		my @return;
		foreach my $pre(@pre_comb)
		{
			foreach my $now (@$crr)
			{
				push @return, $pre.$now;
			}
		}
		return @return;
	}
}