#!/usr/bin/perl -w
#use strict;

open (IN, $key_file) or die;
my %hash;
while (<IN>)
{
	chomp; 
	# op1
	$hash{$_} = '';
	foreach (1..140)
	{
		vec($hash{$id}, $_, 6) = 0
	}
	
}

foreach my $f(@files)
{
	open (IN, $f) or die;
	while (<IN>)
	{
		# opata_1 A
		# synta_90 T
		
	}
}

foreach my $acc_id (keys %hash)
{
	my $string = $hash{$acc_id};
	foreach $ind (0..$max)
	{
		$bit_value = vec($string, $ind, 6);
		$allele = $code{$bit_value}
	}
}