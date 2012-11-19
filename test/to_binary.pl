#!/usr/bin/perl -w
use strict;
my $str;
while (1)
{
	$str = "";
	print "Please input an integer(or q to quit): ";
	my $num = <STDIN>;
	chomp($num);
	exit if $num =~/^q/i;
	unless($num =~/\D/)
	{
		my $binary = to_binary($num);
		print "Binary number is: ", $binary, "\n";
	}	
}

sub to_binary
{
	my $num = shift;
	my $r = $num%8;
	
	if($num>=8)
	{
		to_binary(int($num/8))
	}
	$str .= $r;
}