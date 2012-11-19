#!/usr/bin/perl -w
use strict;
use Statistics::ChisqIndep;

my $sUsage = qq(

perl $0 
<allele_freq file>
<output file>
);

die $sUsage unless @ARGV >= 2;
my ($allele_freq_file, $output) = @ARGV;

my $chisq_test = new Statistics::ChisqIndep;
my %snps = read_freq_file($allele_freq_file);
open (OUT, ">$output") or die;

sub read_freq_file
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while (<IN>)
	{
		# Kukri_mira1_c25815      354     A_71    T_30
		# Kukri_mira1_c25815      381     C_37    T_55
		chomp;
		my @t = split /\s+/, $_;
		
	}	
}


sub sum
{
	my $return = 0;
	map{$return += $_} @_;
	return $return;
}