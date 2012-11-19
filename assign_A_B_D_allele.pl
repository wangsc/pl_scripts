#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<A_B_allele.out>
<ABD_allele_checked.out>
);
die $sUsage unless @ARGV >=2;
my ($ab_file, $abd_file) = @ARGV;
my %ab_alleles = read_ab_file($ab_file);
open (IN, $abd_file) or die;
while (<IN>)
{
	# Kukri_mira1_rep_c99352  140     GT      G
	# Kukri_mira1_rep_c68280  1029    GG      A
	chomp; 
	my @t = split /\s+/, $_;
	my $id = join("\t", @t[0, 1]);
	my @arr = split //, $t[2];
	if($arr[0] eq $arr[1])
	{
		print join("\t", (@t[0,1], @arr, $t[3])), "\n";
		next;
	}
	next unless exists $ab_alleles{$id};
	my($a_allele, $b_allele) = @{$ab_alleles{$id}};
	if($arr[0] eq $a_allele)
	{
		print join("\t", (@t[0,1], @arr, $t[3])), "\n";
	}
	if ($arr[1] eq $a_allele)
	{
		print join("\t", (@t[0,1], @arr[1, 0], $t[3])), "\n";
	}
		
}
close IN;

sub read_ab_file
{
	my $file = shift;
	open (IN, $file) or die "can't open file $file \n";
	my %return;
	while(<IN>)
	{
		chomp;
		my @t = split /\s+/, $_;
		my $id = join("\t", @t[0, 1]);
		$return{$id} = [@t[2,3]];	
	}
	close IN;
	return %return;	
}


