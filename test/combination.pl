#!/usr/bin/perl -w
use strict;

my @arr;
func(1);

sub func
{
	my $n = shift;
	foreach my $i (1..3)
	{
		$arr[$n] = $i;
		if($n<3)
		{
			func($n+1);
		}
		else
		{
			print join(" ", @arr[1..3]), "\n";
		}
		print "\n";
	}
}

=head
void Func(n) 
{ 
    for i = 1 to 3 
    { 
       A[n] = i 
       if (n<3) 
            Func(n+1) 
	else
	    print A[1],A[2],A[3]
    } 
} 
=cut