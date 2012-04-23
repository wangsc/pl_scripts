#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $sUsage = qq(
# Read the blast result file and get the best hit.
perl $0 
<blast_result_file>
<output file>
<pasa_chr3A_clean_masked.pasa_assemblies.denovo_transcript_isoforms.gff3>
<evalue cutoff, like 1e-05, optional>
<minimum length of amino acid for best HSP, optional>
);
die $sUsage unless @ARGV >= 3;

my ($blast_file, $out_file, $iso_gff, $evalue, $minlen) = @ARGV;
$evalue = 1e-05 unless defined $evalue;
$minlen = 60 unless defined $minlen;
my %iso_gene = read_isoform_gff($iso_gff);

open (OUT, ">$out_file") or die $!;

my $blast = Bio::SearchIO->new(-foramt=>'blast', -file => $blast_file);
while (my $result = $blast->next_result)
{
	print STDERR 'Query: ', $result->query_name, "\n";
	while (my $hit = $result->next_hit)
	{
		last if $hit->significance > $evalue;
		next unless $hit->name =~ /Os01g|bradi2g/i; # rice
		my $hsp = $hit->next_hsp;
#		print '$hsp->length(hit): ', $hsp->length('hit'),"\n";
		next unless $hsp->length('hit') >= $minlen;
		print OUT $result->query_name, "\t", $hit->name, "\t", $iso_gene{$result->query_name},"\n";
		last;
	}
}

sub read_isoform_gff
{
	my $file = shift;
	my %iso_gene;
	my %iso_contig;
	open (IN, "$file") or die "$! $file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my ($gene, $iso_id) = $_=~/ID=(.*?)\-(.*?)\;/;
		$iso_gene{$iso_id} = $gene;
		my @data = split /\t/, $_;
		$iso_contig{$gene} = shift @data;
	}
	close IN;
	return %iso_gene;
}