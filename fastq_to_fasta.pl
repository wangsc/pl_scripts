#!/usr/bin/perl -w
use strict;
# For quick convert fastq to fasta, without considering the quality

my $sUsage = "perl $0 fastq_file fasta_file\n";
die $sUsage unless @ARGV >= 2;
my($fastq_file, $fasta_file) = @ARGV;

open (IN, "$fastq_file") or die;
open (OUT, ">$fasta_file") or die;
my $count=0;
my @temp;
while(<IN>)
{
	next if /^\s+$/;
	$count++;
	if($count==1){s/^\@/\>/; print OUT $_}
	elsif($count==2){print OUT $_}
	elsif($count==4){$count=0;next}
}
close IN;
close OUT;
