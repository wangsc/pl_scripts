#!/usr/bin/perl -w
use strict;

# compare the output of calculate_PHS.pl for spring and winter wheat

my $sUsage = "perl $0 <spring phs file> <winter phs file>\n";
die $sUsage unless @ARGV >= 2;
my ($sping_phs_file, $winter_phs_file) = @ARGV;

my %sphs = read_phs_file($sping_phs_file);
my %wphs = read_phs_file($winter_phs_file);

sub read_phs_file
{
	my @files = @_;
	my %return;
	foreach my $file (@files)
	{
		open (IN, $file) or die $!;
		while (<IN>)
		{
			# Chr	SNP_index	SNP_name	Genetic_dist	PHS_A	PHS_B	Freq_A	Freq_B	Block_A	Block_B
			next if /^chr/i;
			chomp;
			my @data = split /\s+/,$_;
			push @{$return{$data[2]}}, [ [@data[4, 6, 8]], [@data[5, 7, 9]] ];
		}
		close IN;
	}
	return %return;
}