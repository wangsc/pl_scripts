#!/usr/bin/perl -w
use strict;

# generate files for each fasta ID.

my $fasta_in = shift or die "perl $0 <fasta file>\n";

open (IN, $fasta_in) or die;

while (<IN>)
{
	if(/>(\S+)/)
	{
		my $out = $1 . ".fasta";
		open (OUT, ">$out") or die;
		print OUT $_;
	}
	else
	{
		print OUT $_
	}
}