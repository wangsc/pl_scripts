#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 *_rmdup_allele_freq.out.pvalue2_adjusted_addeffect_balanced";
die $sUsage unless @ARGV >= 2;

my @files = @ARGV;


my %recorder;
foreach my $f (@files)
{
	open (IN, $f) or die "can't open file $f \n";
	while (<IN>)
	{
		next unless /SILEN/i;
		chomp;
		my @t = split /,/, $_;
		my $gene = $t[0];
		my @effects = @t[-3, -2, -1];
		
		foreach my $m (0..$#effects)
		{
			next if $effects[$m] =~ /NA/;
			my @arr = map{$_ unless $_ == $m}0..$#effects;
			my $eff = join("\t", @effects[@arr]);
			if($effects[$m] =~ /SILEN/)
			{
				if($eff=~/BAL|UP/)
				{
					$recorder{$gene}[$m][0] = 0 unless defined $recorder{$gene}[$m][0];
					$recorder{$gene}[$m][0]++
				}
			}
			else
			{
				if($eff=~/SILEN/)
				{
					$recorder{$gene}[$m][1] = 0 unless defined $recorder{$gene}[$m][1];
					$recorder{$gene}[$m][1]++
				}				
			}
		}
	}
	close IN;
}

my @arr = qw(A B D);

foreach my $gene (keys %recorder)
{
	print "! $gene\n";
	my @tmp = @{$recorder{$gene}};
	my @out;
	foreach (0..$#tmp)
	{
		next unless defined $tmp[$_];
		next unless @{$tmp[$_]} > 0;
		my ($sil, $nsil) = @{$tmp[$_]};
		next unless defined $sil and defined $nsil;
		if($sil > 0 and $nsil > 0)
		{
			push @out, join("|", ($arr[$_], $sil, $nsil));
		}
	}
	print join("\t", ($gene, @out)), "\n" if @out;
}






