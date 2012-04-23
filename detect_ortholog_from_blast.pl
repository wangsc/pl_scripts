#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $blast_file = shift or die "perl $0 blast_result\n";

my %orthologs = parse_blast_results($blast_file);

foreach my $id (keys %orthologs)
{
	print join("\t", ($id, @{$orthologs{$id}})), "\n";
}


sub parse_blast_results
{
	my $file = shift;
	my %return_hash;
	my $searchio = Bio::SearchIO->new(-format => 'blast', file => "$file" );
	while (my $result = $searchio->next_result())
	{
		last unless defined $result;
		my $query_name = $result->query_name;
		my $query_length = $result->query_length;
		next unless $query_length >= 300;
		my $hit = $result->next_hit();		
		next unless defined $hit;
		my $hit_name = $hit->name;
		my $hit_length = $hit->length;
		my $cip; 
		my $calp;
		my $id_len; 
		my $aligned_len;
		while(my $hsp = $hit->next_hsp)
		{
			last unless defined $hsp;
			$id_len += $hsp->num_identical;
			$aligned_len += $hsp->length('query');
		}
		next unless defined $aligned_len;
		$cip = 3*$id_len/$aligned_len;
		$calp = $aligned_len/($query_length);
		#print STDERR '$cip ', $cip,"\n", '$calp ', $calp,"\n";
		push @{$return_hash{$query_name}}, $hit_name if $cip >=0.6 and $calp>=0.7; # $cip >=0.6 and $calp>=0.7
		print STDERR join("\t", ($query_name, $query_length, $hit_length, $aligned_len, $id_len, $cip, $calp, $hit_name)), "\n";
		
	}
	return %return_hash;
}