#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
< ld.out >
< sub_dist.out >
< window_size, default 5 >
< step_size, default 1 >
);
die $sUsage unless @ARGV;
my ($ld_file, $dist_file, $window_size, $step_size) = @ARGV;
$window_size = 5 unless defined $window_size;
$step_size = 1 unless defined $step_size;

my %ld = read_ld_file($ld_file);

# read dist file
my (@snps, @pos);
open (D, $dist_file) or die;
while(<D>)
{
	chomp; 
	my @t = split /\s+/, $_;
	push @snps, $t[0];
	push @pos, $t[1];
}
close D;

# output
my $index = 0;
while($index < $#snps)
{
	last if ($index + $window_size - 1) > $#snps;
	my @pairs = pairs(@snps[$index..($index + $window_size - 1)]);
	my $average = average(@ld{@pairs});
	print $snps[$index], "\t", $pos[$index], "\t", $average, "\n";
	$index += $step_size;
}

# subroutines
sub read_ld_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		chomp;
		my @t=split /\s+/, $_;
		$return{join("\t", sort{$a cmp $b}@t[0,1])} = $t[2];
	}
	close IN;
	return %return
}


sub average
{
	my @arr = @_;
	my $sum = 0;
	map{$sum+=$_}@arr;
	return $sum/(scalar @arr);
}

sub pairs
{
	my @arr = @_;
	my @return;
	foreach my $i (0..$#arr-1)
	{
		foreach my $j ($i+1 .. $#arr)
		{
			push @return, join("\t", sort{$a cmp $b} @arr[$i, $j]);
		}
	}
	return @return;
}