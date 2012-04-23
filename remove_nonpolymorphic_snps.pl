#!/usr/bin/perl -w
use strict;
use Tie::File;
my $file = shift or die "perl <file>\n";

tie my @file_array, "Tie::File", $file or die $!;
my @count;
my %record_polymorph;
foreach my $line (1..$#file_array)
{
	my $data = $file_array[$line];
	$data =~ s/\n//;
	my @t = split /\t/, $data;
	foreach (1..$#t)
	{
		$record_polymorph{$_}{$t[$_]}++;;
	}
}

my ($un_poly, $not_work) = (0, 0);
my %removed_index;
foreach my $index (keys %record_polymorph)
{
	my @g = keys %{$record_polymorph{$index}};
	next if @g > 2;
	if (@g == 1)
	{
		$removed_index{$index} = 0;
		($g[0]=~/\?/)?$not_work++:$un_poly++;
	}
	else
	{
		if (($g[0] =~ /\?/) or ($g[1] =~ /\?/)){$un_poly++; $removed_index{$index} = 0;}
	}
}

print STDERR 'Unpolymorphic SNPs: ', $un_poly,"\n";
print STDERR 'Not working: ', $not_work,"\n";

foreach my $line (0..$#file_array)
{
	my $data = $file_array[$line];
	$data =~ s/\n//;
	my @t = split /\t/, $data;
	print $t[0];
	foreach (1..$#t)
	{
		print "\t", $t[$_] unless exists $removed_index{$_};
	}
	print "\n";
}









