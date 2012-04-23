#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <All_chr_concatenated.fasta>\n";
my $input_file = shift or die $sUsage;

my %file_handle;
foreach my $chr (1..7)
{
	my $outfile = 'Chr_' . $chr.".fasta";
	local *FH;
	open (FH, ">$outfile") or die;
	$file_handle{$chr} = *FH;
}

open (IN, $input_file) or die;
my $handle;
while(<IN>)
{
	next if /^\s+$/;
	if(/>(\d)/)
	{
		$handle = $file_handle{$1};
	}
	print $handle $_;
}


