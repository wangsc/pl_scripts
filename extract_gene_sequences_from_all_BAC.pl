#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage = qq(
perl $0
<gene structure file in gff>
<genomic fasta file>
<cDNA output file>
);
die $sUsage unless @ARGV >= 3;
my ($gff_file, $genome_fasta, $output) = @ARGV;

my $gnobj = Bio::DB::Fasta->new($genome_fasta);
my %gene_structure = read_gff($gff_file);

open (OUT, ">$output") or die;
foreach my $id (keys %gene_structure)
{
	my $bac_id = $1 if $id =~ /^(\S+)_\d+$/;
	my $seq = "";
	foreach (@{$gene_structure{$id}})
	{
		my ($start, $end) = @$_;
		$seq .= $gnobj->seq($bac_id, $start=>$end);
	}
	print OUT ">$id\n", $seq, "\n";
}
close OUT;

# 
sub read_gff
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while (<IN>)
	{
		chomp;
		my @t = split /\s+/, $_;
		my $id = $1 if /ID=(\S+?)\;/ ;
		push @{$return{$id}}, [@t[3,4]];
	}
	close IN;
	return %return;
}