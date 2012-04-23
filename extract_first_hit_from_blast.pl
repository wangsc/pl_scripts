#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $sUsage = qq(
# Read the blast result file and get the best hit.
perl $0 
<blast_result_file>
<output file>
<evalue cutoff, like 1e-05, optional>
<minimum length of hit, optional>
);
die $sUsage unless @ARGV >= 2;

my ($blast_file, $out_file, $evalue, $minlen) = @ARGV;
$evalue = 1e-05 unless defined $evalue;
$minlen = 100 unless defined $minlen;
my $hsp_len_cutoff = 1000;
open (OUT, ">$out_file") or die $!;

my $blast = Bio::SearchIO->new(-foramt=>'blast', -file => $blast_file);
while (my $result = $blast->next_result)
{
	print STDERR 'Query: ', $result->query_name, "\n";
	my $count_hit = 0;
	while (my $hit = $result->next_hit)
	{
		last if $hit->significance > $evalue;
		last unless $hit->length >= $minlen;
		print OUT $result->query_name, "\t", $hit->name,"\n";
		$count_hit++;
		while(my $hsp = $hit->next_hsp)
		{
			last if $hsp->hsp_length < $hsp_len_cutoff;
			print OUT "\t\t", $hsp->length('hit'),"\t", $hsp->num_identical,"\n";
		}
		last if $count_hit;
	}
}