#!/usr/bin/perl -w
use strict;

my $file = shift or die "perl $0 [file] \n";
my @columns = 1..11;
my @plot = (0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2);
open (IN, "$file") or die "$! $file\n";
my %record;
my @types;
my $totalgenes = 0;
while (<IN>)
{
	next if /^\s+$/;
	if(/ID/)
	{
		chomp;
		@types = split /\t/, $_;
		next;
	}
	$totalgenes++;
	chomp;
	my @data = split /\t/, $_;
	foreach my $index (@columns)
	{
		die '>1: ',$data[$index],"\n" if $data[$index]>1;
		foreach (0..($#plot -1))
		{
			$record{$index}->[$_] = 0 unless defined $record{$index}->[$_];
			if($data[$index] >= $plot[$_] and $data[$index] < $plot[$_+1])
			{				
				$record{$index}->[$_]++;
			}
		}
	}
}
close IN;
foreach my $col (sort{$a<=>$b}keys %record)
{
	print "\t",$types[$col], "\n";
	my @counts = @{$record{$col}};
	foreach my $i (0..$#counts)
	{
		print $i*0.2, '-',($i+1)*0.2,"\t", $counts[$i],"\n";
	}
	print 'Total',"\t", $totalgenes,"\n\n\n";
}