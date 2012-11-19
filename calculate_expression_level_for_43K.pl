#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<43k cluster list>
<gene expression file>
);
die $sUsage unless @ARGV;
my ($cluster_file, $exp_file) = @ARGV;
my %expression_values = read_expression_file($exp_file);

open (IN, $cluster_file) or die "can't open file $cluster_file \n";

while (<IN>)
{
	chomp;
	my @t = split /\s+/, $_;
	my @ids = split /:/, $t[1];
	my $cnt = 0;
	my $exp = 0;
	foreach my $id (@ids)
	{
		next unless exists $expression_values{$id};
		$cnt++;
		$exp += $expression_values{$id};
	}
	print $t[0], "\t", sprintf("%.1f", $exp/$cnt), "\n";
}
close IN;

# Subroutine
sub read_expression_file
{
	my $file = shift;
	open (IN, $file) or die "can't open file $file\n";
	my %return;
	while (<IN>)
	{
		my @t = split /\s+/, $_; 
		$return{$t[0]} = $t[1]
	}
	close IN;
	return %return;
}