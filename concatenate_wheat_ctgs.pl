#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <input fasta> <output fasta>\n";
die $sUsage unless @ARGV >= 2;
my ($input, $output) = @ARGV;

open (IN, $input) or die;
open (OUT, ">$output") or die;

my $joint = 'X'x50;
my %input_seq;
my $id;
my $flag;
while(<IN>)
{
	chomp;
	next if /^\s+$/;
	if(/>(\S{2})/)
	{
		$id = $1;
		$flag = 1;
		next;
	}
	my $seq = $_;
	if($flag){$seq = $joint . $seq if(exists $input_seq{$id}); $flag = 0}
	$input_seq{$id} .= $seq;
}
close IN;

foreach (sort{$a cmp $b} keys %input_seq)
{
	print OUT ">", $_,"\n", $input_seq{$_},"\n";
}

close OUT;