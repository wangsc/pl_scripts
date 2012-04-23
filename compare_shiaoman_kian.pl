#!/usr/bin/perl -w
use strict;
my $sUsage = "perl $0 <9000 SNP file> <shiaoman's list> <kian's cleaned file>\n";
die $sUsage unless @ARGV >= 3;
my ($big_file, $shiao_file, $kian_file) = @ARGV;

my %big_index = read_9k_file($big_file, 5);
my %shiao_index = read_9k_file($shiao_file, 0);
my %kian_index = read_9k_file($kian_file, 2);

# 
my %shiao_not_called;
my %kian_not_called;
my %shiao_not_kian_yes;
my %shiao_yes_kian_no;

foreach my $ind (keys %big_index)
{
	$kian_not_called{$ind} = 1 unless exists $kian_index{$ind};
	$shiao_not_called{$ind} = 1 unless exists $shiao_index{$ind}; 
}

foreach my $ind (keys %shiao_not_called)
{
	$shiao_not_kian_yes{$ind} = 1 if exists $kian_index{$ind};
}

foreach my $ind (%kian_not_called)
{
	$shiao_yes_kian_no{$ind} = 1 if $shiao_index{$ind};
}

my @shiao_not = sort{$a<=>$b} keys %shiao_not_called;
my @kian_not = sort {$a<=>$b} keys %kian_not_called;
my @shiao_not_k_yes = sort {$a<=>$b} keys %shiao_not_kian_yes;
my @shiao_yes_k_no = sort {$a<=>$b} keys %shiao_yes_kian_no;
my $max = max(scalar @shiao_not, scalar @kian_not, scalar @shiao_not_k_yes, scalar @shiao_yes_k_no);

print join("\t",qw(Shiao_not_called Kian_not_called Shiao_not_Kian_called Shiao_called_Kian_not)),"\n";
foreach (0..($max-1))
{
	defined $shiao_not[$_]?print $shiao_not[$_]:print ' ';
	defined $kian_not[$_]?print "\t", $kian_not[$_]:print "\t,";
	defined $shiao_not_k_yes[$_]?print "\t",$shiao_not_k_yes[$_]:print "\t";
	defined $shiao_yes_k_no[$_]?print "\t", $shiao_yes_k_no[$_]:print "\t";
	print "\n";
}




sub max
{
	my $max = shift;
	foreach (@_)
	{
		$max = $_ if $_>$max
	}
	return $max;
}

sub read_9k_file
{
	my $file = shift;
	my $index = shift;
	my %return;
	open (IN, $file) or die "$file $!";
	while (<IN>)
	{
		next unless /wsnp_/;
		my @data;
		if (/,/)
		{
			@data = split /,/,$_;
		}else
		{
			@data = split /\s+/,$_;
		}
		$return{$data[$index]} = 1;
	}
	close IN;
	return %return;
}
