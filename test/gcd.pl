#!/usr/bin/perl -w
use strict;

# Given two numbers, find out the greatest common divisor

die "perl $0 number_a number_b\n" unless @ARGV == 2;
my ($num_a, $num_b) = @ARGV;
if($num_a =~ /\D/ or $num_b =~ /\D/)
{
	print "Inputs: $num_a $num_b \n";
	print "ERROR: input should be digits only! \n";
	die "perl $0 number_a number_b\n"
}

print join("\t", ($num_a, $num_b)), ":\n";
my $gcd = gcd($num_a, $num_b);
print "GCD: ",$gcd , "\n";
my $mcm = $num_a * $num_b /$gcd;
print "MCM: ", $mcm, "\n";

sub gcd
{
	my ($na, $nb) = @_;
	return $na if $nb == 0;
	return $nb if $na == 0;
	return 1 if $na == 1 or $nb == 1;
	if($na==$nb)
	{
		return $nb;
	}
	else
	{
		my $res = abs($nb - $na);
		gcd($res, $nb>$na?$na:$nb);
	}
}
