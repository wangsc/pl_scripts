#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

parse_blast_results(shift);

sub parse_blast_results
{
	my $file = shift or die;
	my %return_hash;
	my $searchio = Bio::SearchIO->new(-format => 'blast', file => "$file" );
	while (my $result = $searchio->next_result())
	{
		last unless defined $result;
		my $query_name = $result->query_name;
		print ">", $query_name,"\n";
		my $flag = 0;
		while(my $hit = $result->next_hit())
		{
			my $hit_name = $hit->name;
			while(my $hsp = $hit->next_hsp())
			{
				my $len = $hsp->length('hit');
				if($len >= 100)
				{
					print $hit_name,"\n"; last;
				}
			}
		}
		print "\n";
	}
	return %return_hash;
}