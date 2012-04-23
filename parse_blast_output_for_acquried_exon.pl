#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

die "perl $0 <blast_file>\n" unless @ARGV;

parse_blastx_results(shift);

sub parse_blastx_results
{
	my $file = shift;
	my %return_hash;
	my $searchio = Bio::SearchIO->new(-format => 'blast', file => "$file" );
	while (my $result = $searchio->next_result())
	{
		last unless defined $result;
		my $query_name = $result->query_name;
		print $query_name;
		my $flag = 0;
		while(my $hit = $result->next_hit())
		{
			$flag = 1;
			print "\t",$hit->name,"\t",$hit->description,"\t", $hit->significance,"\n";
		}
		print "\n";
	}
	return %return_hash;
}