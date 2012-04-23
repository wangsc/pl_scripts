#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 *.csv\n";
die $sUsage unless @ARGV >= 1;
my @csv_files = @ARGV;

my %recorder;
my @accessions;
# read all csv files
#
foreach my $file (@csv_files)
{
	open (IN, $file) or die;
	my $counter = 0;
	my %index_accession;
	while(<IN>)
	{
		$counter++; 
		if($counter == 1)
		{
			chomp;
			my @t = split/,/, $_;
			foreach (2..$#t)
			{
				next unless $t[$_] =~ /\S/;
				push @accessions, $t[$_];
				$index_accession{$_} = $t[$_];
			}
			next;
		}
		
		chomp;
		my @t = split/,/, $_;
		my $id = join(",", @t[0,1]);
		foreach (2..$#t)
		{
			my $acc = $index_accession{$_};
			$recorder{$id}{$acc} = $t[$_];
		}					 
	}
	close IN;
}

# OUTPUT
my %unique_acc = map{$_,1} @accessions;
@accessions = keys %unique_acc;

print join(",", (0,0,@accessions)),"\n";
foreach my $id (keys %recorder)
{
	my @genotypes;
	foreach (@accessions)
	{
		my $gn = exists $recorder{$id}{$_}?$recorder{$id}{$_}:"--";
		push @genotypes, $gn;
	} 
	print join(",", ($id, @genotypes)),"\n";
}














