#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $sUsage = qq(
# Read the blast result file and get the best hit.
perl $0 
<blast_result_file>
<output file>
<evalue cutoff, like 1e-10, optional>
);
die $sUsage unless @ARGV >= 2;

my ($blast_file, $out_file, $evalue) = @ARGV;
$evalue = 1e-20 unless defined $evalue;
my $percent_aligned_cutoff = 0.6; # 50% aligned;
open (OUT, ">$out_file") or die $!;
my @ids;
my %redudant_seq;
my $blast = Bio::SearchIO->new(-foramt=>'blast', -file => $blast_file);
my $debug = 1;
while (my $result = $blast->next_result)
{
	my $query_name = $result->query_name;
	print '$query_name ',$query_name,"\n" if $debug;
	push @ids, $query_name;
	next if exists $redudant_seq{$query_name};
	my $query_length = $result -> query_length;
	while (my $hit = $result->next_hit)
	{
		next if $hit->name eq $query_name;
		last unless defined $hit;
		print '$hit_name ', $hit->name,"\n" if $debug; 
		print '$hit->length ', $hit->length,"\n" if $debug; $debug=0;
		die "$query_name \t ", $hit->name,"\n" if $hit->length == 0;
		last if $hit->significance > $evalue;
		next unless $hit->length <= $query_length;
		my $total_hsp_length = 0;
		while(my $hsp = $hit->next_hsp)
		{
			last unless defined $hsp;
			$total_hsp_length += $hsp->length('hit');
		}
		if(($total_hsp_length/$hit->length) >= $percent_aligned_cutoff)
		{
			$redudant_seq{$hit->name} = 1;
		}
	}
}

foreach my $id (@ids)
{
	print OUT $id,"\n" unless exists $redudant_seq{$id}
}
close OUT;