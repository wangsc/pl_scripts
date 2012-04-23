#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <pasa_genes_with_acquired_exons.blastx.out> <pasa_genes_with_acquired_exons.fasta>\n";
die $sUsage unless @ARGV >= 2;
my($blast_result, $fasta) = @ARGV;
my %reading_frame = get_reading_frame();




sub get_reading_frame
{
	my $file = shift;
	my %return_hash;
	my $searchio = Bio::SearchIO->new(-format => 'blast', file => "$file" );
	while (my $result = $searchio->next_result())
	{
		last unless defined $result;
		my $query_name = $result->query_name;
	#	print $query_name, "\t";
		my $hit = $result->next_hit();
		next unless defined $hit;
		my $hsp = $hit->next_hsp();
		next unless defined $hsp;
		my $blastframe = ($hsp->query->frame + 1) * $hsp->query->strand;
	#	print $hit->name, "\t", $blastframe,"\n";
	#	$return_hash{$query_name} = [] unless exists $return_hash{$query_name};
		$return_hash{$query_name} = [$result->query_length, $blastframe];
	}
	return %return_hash;
}

sub check_stop_codon
{
	my $seq = shift;
	my $s = 0;
	while ($s+3<=(length $seq))
	{
		my $aa = codon_table(substr($seq, $s, 3));
		return 1 if $aa eq '_';
		$s += 3;
	}
	return 0;
}