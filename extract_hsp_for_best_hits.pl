#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $sUsage = qq(
# Read the blast result file and get the best hit.
perl $0 
<blast_result_file>
<blast result format: must be blastxml, ...>
<output file>
<evalue cutoff, like 1e-05, optional>
<minimum length of hit, optional>
);
die $sUsage unless @ARGV >= 3;

my ($blast_file, $blast_format, $out_file, $evalue, $minlen) = @ARGV;
$evalue = 1e-05 unless defined $evalue;
$minlen = 30 unless defined $minlen;
open (OUT, ">$out_file") or die $!;

open (IN, $blast_file) or die;
my $flag = 0;
while (<IN>)
{
	next if /^\s+$/;
	if(/\<\?xml/)
	{
		open (TMP, ">tmp.out") or die;
	}
	print TMP $_;
	if(/\<\/BlastOutput\>/)
	{
		close TMP;
		parser();
	}	
}


sub parser
{
	my $blast = Bio::SearchIO->new(-format=>$blast_format, -file => 'tmp.out');
	while (my $result = $blast->next_result)
	{
		#print STDERR 'Query: ', $result->query_description, "\n";
		my $count_hit = 0;
		while (my $hit = $result->next_hit)
		{
			last if $hit->significance > $evalue;
			last unless $hit->length >= $minlen;
			print OUT $result->query_description, "\t", $hit->name,"\t", $hit->description, "\t", $hit->significance;
			$count_hit++;
			while(my $hsp = $hit->next_hsp)
			{
				print OUT "\t", $hsp->start('query'), "_", $hsp->end('query'),"_", ($hsp->query->frame + 1) * $hsp->query->strand;
			}
			print OUT "\n";
			last if $count_hit;
		}
		print $result->query_description, "\tNo_hit\n" unless $count_hit;
	}	
}
