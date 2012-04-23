#!/usr/bin/perl -w
use strict;

my $file = shift or die "perl $0  <ac_exon fasta file>\n";
my %exon_fasta = read_fasta($file);

foreach my $id (keys %exon_fasta)
{
	my @count_stop;
	my $seq = $exon_fasta{$id};
	my $rev_seq = reverse_complemantary($seq);
	my @seq_trans = (translate($seq), translate($rev_seq));
	foreach (@seq_trans)
	{
		my $n = s/_/_/g;
		$n = 0 unless $n;
		push @count_stop, $n;
		my $seq = $_;
		$seq =~ s/_//g;
	}
	print join("\t", ($id, @count_stop)),"\n";
}

sub translate
{
	my $seq = shift;
	my @array;
	foreach my $frame (0..2)
	{
		my $pos = $frame;
		my $aa_str = '';
		while( ($pos+3)<=(length $seq) )
		{
			my $code = substr($seq, $pos, 3);
			my $aa = codon_table($code);
			$aa_str .= $aa;
			$pos += 3;
		}
		push @array, $aa_str;
	}
	return @array;
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

sub read_fasta
{
	my $file = shift;
	open (IN, "$file") or die "$! $file\n";
	my %return_hash;
	my $id;
	my $debug = 1;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		if (/^>(\S+)/)
		{
			$id = $1;
			print STDERR 'fasta id: ', $id ,"\n" if $debug; $debug = 0;
			$return_hash{$id} = '';
			next;
		}
		$return_hash{$id} .= $_;
	}
	close IN;
	return %return_hash;
}

sub reverse_complemantary
{
	my $seq = shift;
	$seq = reverse $seq;
	$seq =~ tr/[ATGC]/[TACG]/;
	return $seq;
}
