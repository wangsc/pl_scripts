#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <blat output files>\n";
die $sUsage unless @ARGV;
my @files = @ARGV;

my $mapped_length_cutoff = 35;
my %blat_result;
foreach my $file (@files)
{
	open (IN,"<$file") or die "can't open file $file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my @line_data = split /\t/, $_;
		my ($oligo_id, $contig, $mapped_length, $mismatch, $gap) = @line_data[0, 1, 3..5];
		next if $mapped_length < $mapped_length_cutoff;
		$blat_result{$oligo_id} = [0, 0, 0, 0] unless exists $blat_result{$oligo_id} ;
		if ($mismatch == 0 and $gap == 0)
		{
			$blat_result{$oligo_id}->[0]++;
		}
		elsif ($mismatch != 0 and $gap == 0)
		{
			$blat_result{$oligo_id}->[1]++;
		}
		elsif ($mismatch == 0 and $gap != 0)
		{
			$blat_result{$oligo_id}->[2]++;
		}
		else
		{
			$blat_result{$oligo_id}->[3]++;
		}
	}
	close IN;	
}


foreach (keys %blat_result)
{
	print $_,"\t";
	print join("\t", @{$blat_result{$_}}),"\n";
}