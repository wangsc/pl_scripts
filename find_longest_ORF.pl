#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;

my $fasta_file = shift or die "Need a fasta file as input!\n";
my $seqin = Bio::SeqIO->new(-format => 'fasta', -file =>$fasta_file);
my @frames = (1, 2, 3, -1, -2, -3);

while (my $seqobj = $seqin->next_seq())
{
	print "\n", $seqobj->display_id(),"\n";
	foreach my $frame (@frames)
	{
		print 'Reading frame: ', $frame,"\n";
		print 'Peptide:       ', translate($seqobj->seq(), $frame),"\n";
	}
}


sub translate
{
	my ($seq, $frame) = @_;
	if ($frame < 0){$seq = rev_comp($seq)}
	my $start = abs($frame) - 1;
	my $pro = '';
	while ($start < ((length $seq)-1) )
	{
		my $codon = substr($seq, $start, 3);
		last if (length $codon < 3);
		my $aa = codon_table($codon);
		$pro .= $aa;
		$start += 3;
	}
	return $pro;
}

sub rev_comp
{
	my $seq = shift;
	$seq = reverse $seq;
	$seq =~ tr/[ATGCatgc]/[TACGtacg]/;
	return $seq;
}

sub codon_table
{
	my %g=(
		'TCA'=>'S',#Serine
		'TCC'=>'S',#Serine
		'TCG'=>'S',#Serine
		'TCT'=>'S',#Serine
		'TTC'=>'F',#Phenylalanine
		'TTT'=>'F',#Phenylalanine
		'TTA'=>'L',#Leucine
		'TTG'=>'L',#Leucine
		'TAC'=>'Y',#Tyrosine
		'TAT'=>'Y',#Tyrosine
		'TAA'=>'_',#Stop
		'TAG'=>'_',#Stop
		'TGC'=>'C',#Cysteine
		'TGT'=>'C',#Cysteine
		'TGA'=>'_',#Stop
		'TGG'=>'W',#Tryptophan
		'CTA'=>'L',#Leucine
		'CTC'=>'L',#Leucine
		'CTG'=>'L',#Leucine
		'CTT'=>'L',#Leucine
		'CCA'=>'P',#Proline
		'CAT'=>'H',#Histidine
		'CAA'=>'Q',#Glutamine
		'CAG'=>'Q',#Glutamine
		'CGA'=>'R',#Arginine
		'CGC'=>'R',#Arginine
		'CGG'=>'R',#Arginine
		'CGT'=>'R',#Arginine
		'ATA'=>'T',#Isoleucine
		'ATC'=>'T',#Isoleucine
		'ATT'=>'T',#Isoleucine
		'ATG'=>'M',#Methionine
		'ACA'=>'T',#Threonine
		'ACC'=>'T',#Threonine
		'ACG'=>'T',#Threonine
		'ACT'=>'T',#Threonine
		'AAC'=>'N',#Asparagine
		'AAT'=>'N',#Asparagine
		'AAA'=>'K',#Lysine
		'AAG'=>'K',#Lysine
		'AGC'=>'S',#Serine#Valine
		'AGT'=>'S',#Serine
		'AGA'=>'R',#Arginine
		'AGG'=>'R',#Arginine
		'CCC'=>'P',#Proline
		'CCG'=>'P',#Proline
		'CCT'=>'P',#Proline
		'CAC'=>'H',#Histidine
		'GTA'=>'V',#Valine
		'GTC'=>'V',#Valine
		'GTG'=>'V',#Valine
		'GTT'=>'V',#Valine
		'GCA'=>'A',#Alanine
		'GCC'=>'A',#Alanine
		'GCG'=>'A',#Alanine
		'GCT'=>'A',#Alanine
		'GAC'=>'D',#AsparticAcid
		'GAT'=>'D',#AsparticAcid
		'GAA'=>'E',#GlutamicAcid
		'GAG'=>'E',#GlutamicAcid
		'GGA'=>'G',#Glycine
		'GGC'=>'G',#Glycine
		'GGG'=>'G',#Glycine
		'GGT'=>'G',#Glycine
	);
	my $codon = uc shift;
	my $aminoacid = (exists $g{$codon})?$g{$codon}:'?';
	return $aminoacid;
}
