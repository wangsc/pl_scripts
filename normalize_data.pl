#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;

my $file = shift or die "perl $0 input\n";
my %mean_var = calculate_mean_var($file);

open (IN, $file) or die;
while(<IN>)
{
	if(/^Name/i){print $_; next}
	chomp;
	my @t=split/,/,$_;
	foreach my $index (1..$#t)
	{
		my $norm = ($t[$index]-$mean_var{$index}->[0])/sqrt($mean_var{$index}->[1]);
		$t[$index] = $norm unless $t[$index] == -10;
	}
	print join(",", @t),"\n";
}
close IN;

sub calculate_mean_var
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while(<IN>)
	{
		next if /^name/i;
		chomp; 
		my @t=split/,/, $_;
		foreach my $ind (1..$#t)
		{
			push @{$return{$ind}}, $t[$ind] unless $t[$ind] == -10;
		}
	}
	close IN;
	foreach (keys %return)
	{
		my @data = @{$return{$_}};
		my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@data); 
    my $mean = $stat->mean();
    my $var  = $stat->variance();
    $return{$_} = [$mean, $var];
	}
	return %return;
}




