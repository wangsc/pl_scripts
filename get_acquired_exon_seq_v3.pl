#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage = qq(
perl $0
<acquired exon list>
<3A_contigs_scaffold_masked.fa>
<pasa_chr3A_clean_masked.pasa_assemblies.denovo_transcript_isoforms.gff3>
);
die $sUsage unless @ARGV >= 3;

my ($list_file, $ctg_fasta, $gff_file) = @ARGV;

my %exons = read_list_file($list_file);
my $ctg = Bio::DB::Fasta->new($ctg_fasta);
my %gene_ctg = read_gff($gff_file);

foreach my $gene (keys %exons)
{
	my $ctg_id = $gene_ctg{$gene};
	foreach (@{$exons{$gene}})
	{
		my ($start, $end) = sort{$a<=>$b}@$_;
		my $seq = $ctg->seq("$ctg_id", $start => $end);
		print ">", join("_",($gene, $start, $end)),"\n", $seq,"\n";
	}	
}

sub read_list_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		next if /^\s+$/;
		s/>//;
		chomp;
		my @t = split /_/, $_;
		push @{$return{$t[0]}}, [@t[1,2]]
	}
	close IN;
	return %return;
}

sub read_gff
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		next if /^\s+$/;
		my $gene = $1 if /ID=(S\d+)\-/;
		my $ctg_id = (split/\s+/,$_)[0];
		$return{$gene} = $ctg_id;
	}
	close IN;
	return %return;
}











