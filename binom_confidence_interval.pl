#!/usr/bin/perl -w
use strict;

sub calcualte_interval_std
{
	my ($p, $n) = @_;
	my $z = 1.96; # 95% interval
	my $residue_value = $z * sqrt($p*(1-$p)/$n);
	return ($p - $residue_value, $p + $residue_value);	
}

sub calcualte_interval_wilson
{
	my ($p, $n) = @_;
	my $z = 1.96; # 95% interval
	my $residue = $z*sqrt($p*(1-$p)/$n + $z*$z/(4*$n*$n));
	return (($p - $residue)/(1+$z*$z/$n), ($p + $residue)/(1+$z*$z/$n))
}

#my ($p, $n) = @ARGV;
#print join("\t", calcualte_interval_std($p, $n)), "\n";
#print join("\t", calcualte_interval_wilson($p, $n)), "\n";

my $file = shift or die "perl $0 <truman_homeolog_sites_allele_freq.out_cleaned_covDep_ratio> <Min coverage>\n";
my $min_dep = shift;
$min_dep = 30 unless defined $min_dep;
open (IN, $file) or die "can't open file $file\n";
while (<IN>)
{
	chomp;
	my @t = split /\s+/,$_;
	my ($ratio, $count) = @t[-1, -2];
	next unless $count >= $min_dep;
	my @interval = calcualte_interval_std(1/3, $count);
	print $_, "\n" if $ratio >= $interval[0] and $ratio <= $interval[1];
}
close IN;
