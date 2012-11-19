#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $sUsage = qq(
perl $0 <blast output>
);
my $blast_file = shift or die $sUsage;
my %allele_count;
my $searchio = Bio::SearchIO->new(-format => 'blast', file => "$blast_file" );
while (my $result = $searchio->next_result())
{
	last unless defined $result;
	my $query_name = $result->query_name;
	while(my $hit = $result->next_hit())
	{
		next unless defined $hit;
		while(my $hsp = $hit->next_hsp())
		{
			next unless $hsp->length('query') > 50;
			my ($start, $end) = $hsp->range('query');
			print STDERR "!!$start\t$end\n" if $start > $end;
			my $seq = $hsp->hit_string;
			my $allele = substr($seq, (51-$start), 1);
			#print  $query_name, "\t", $allele, "\n";
			$allele_count{$query_name}{$allele}++;
		}		
	}
}

foreach my $id (keys %allele_count)
{
	my @arr = sort {$allele_count{$id}{$b} <=> $allele_count{$id}{$a}} keys %{$allele_count{$id}};
	print join ("\t", ($id, @arr, @{$allele_count{$id}}{@arr})), "\n";
}