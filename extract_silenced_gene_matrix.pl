#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
A/B/D
);
die $sUsage unless @ARGV;
my $subgenome = shift;

my @accs = qw(Grandin ID0444 Jaypee Jupateco Steele Stephens Truman);
my %record_hash;

foreach my $acc (@accs)
{
	my $notNA_file = $acc . "_" . $subgenome . "_notNA.txt";
	my $silen_file = $acc . "_" . $subgenome . "_SILENCE.txt";
	my %notNA_list = read_file($notNA_file);
	my %silen_list = read_file($silen_file);
	map{$record_hash{$acc}{$_} = exists $silen_list{$_}?0:1} keys %notNA_list;	
}

# output
my %total_genes;
foreach (keys %record_hash)
{
	map{ $total_genes{$_}=1 } keys %{$record_hash{$_}};
}

print "Gene\t", join("\t", @accs), "\n";
foreach my $g (keys %total_genes)
{
	my @arr;
	foreach (@accs)
	{
		push @arr, exists $record_hash{$_}{$g}?$record_hash{$_}{$g}:"NA";
	}
	print join("\t", ($g, @arr)), "\n";
}


# Subroutines
sub read_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die "can't open file $file \n";
	while(<IN>)
	{
		chomp;
		$return{$_}=1;
	}
	close IN;
	return %return;
}