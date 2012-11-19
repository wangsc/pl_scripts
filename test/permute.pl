#!/usr/bin/perl -w
use strict;
my $string = shift;
my @arr = split //, $string;
die "Usage: perl $0 <string>\n" unless @arr;
my @index;
map{$index[$_]=0}0..$#arr;
my @str;
permute(0, (scalar @arr));
sub permute
{
	my $n = shift;
	my $len = shift;
	foreach my $ind (0..($len-1))
	{
		if ($index[$ind] == 0)
		{
			$index[$ind] = 1;
			$str[$n] = $arr[$ind];
			if($n<($len-1))
			{
				permute($n+1, $len);				
			}
			else
			{
				print @str, "\n";
			}
			$index[$ind] = 0;			
		}
	}
	
}