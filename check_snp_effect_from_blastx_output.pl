#!/usr/bin/perl -w 
use strict;
use Bio::SearchIO;

# Aim: parse the blastx output to determine the reading frame of the sequences which have SNPs located in;
#      then try to see the amino acid alteration caussed by SNPs.
# Input files: 1) the blastx output file; 2) the fasta file used to do blast;
# Output: The effect of SNPs (synonymous or nonsynonymous %->%);
# SW 04.29.2011

my $sUsage = qq(
perl $0 
<blastx result filename>
<fasta file used to do blst>
<output file for snp effect>
);
die $sUsage unless @ARGV >= 3;
my ($blast_result_file, $fasta_file, $output_file) = @ARGV;

# Main
print STDERR "parsing blastx results...\n";
my %blastx_results_hash = parse_blastx_results($blast_result_file);
print STDERR "check SNP effect ...\n";
my %snp_effect_hash = determine_snp_effect(\%blastx_results_hash, $fasta_file);
print STDERR "Output results ...\n";
output_snp_effect(\%snp_effect_hash, $output_file);

# Subroutines

sub output_snp_effect
{
	my ($snpeffect_hashref, $outfile) = @_;
	open (OUT, ">$outfile") or die "can't open file $outfile\n";
	foreach my $snpid (keys %$snpeffect_hashref)
	{
		print OUT $snpid, "\t";
		if ($snpeffect_hashref->{$snpid} eq 'no_hit')
		{
			print OUT $snpeffect_hashref->{$snpid},"\n";
			next;
		}
		my @data = split /:/, $snpeffect_hashref->{$snpid};
		print OUT join("\t", @data),"\n";
	}
	close OUT;
}

sub parse_blastx_results
{
	my $file = shift;
	my %return_hash;
	my $searchio = Bio::SearchIO->new(-format => 'blast', file => "$file" );
	while (my $result = $searchio->next_result())
	{
		last unless defined $result;
		my $query_name = $result->query_name;
		#print $query_name, "\t";
		my $hit = $result->next_hit();
		next unless defined $hit;
		my $hsp = $hit->next_hsp();
		next unless defined $hsp;
		my $blastframe = ($hsp->query->frame + 1) * $hsp->query->strand;
	#	print $hit->name, "\t", $blastframe,"\n";
	#	$return_hash{$query_name} = [] unless exists $return_hash{$query_name};
		$return_hash{$query_name} = [$hit->name, $blastframe];
	}
	return %return_hash;
}

sub determine_snp_effect
{
	my ($blastx_results_hashref, $fastafile) = @_;
	my %return_hash;
	my %fastafile_hash = read_fasta_file($fastafile);
	foreach my $snpid (keys %fastafile_hash)
	{
		unless (exists $blastx_results_hashref->{$snpid}) # and (scalar @{$blastx_results_hashref->{$snpid}}) == 2)
		{
			$return_hash{$snpid} = 'no_hit';
			next;
		}
		# snpid format: Ex_mira1_c37_80868_83059_81536_81536_[A/G]
		my ($refbase, $altbase) = $snpid=~/\[(\w)\/(\w)\]/;
		print STDERR '$refbase, $altbase: ', $refbase, "\t", $altbase, "\n" if $snpid eq 'Ex_mira1_rep_c115305_95195302_95196125_95195380_95195380_[A/G]';
		my ($length_left_flanking, $length_right_flanking);
		my @iddata = split /_/, $snpid;
		if ($snpid =~ /^wsnp/) #wsnp_BQ159467B_Ta_2_1_61_[T/C]
		{
			my $seq_length = length $fastafile_hash{$snpid};
			my $snp_pos = $iddata[-2];
			$length_left_flanking = $snp_pos-1;
			$length_right_flanking = $seq_length - $snp_pos;
		}
		else
		{
			$length_left_flanking = ($iddata[-3]-$iddata[-5])>=100?100:$iddata[-3]-$iddata[-5];
			$length_right_flanking = ($iddata[-4]-$iddata[-3])>=100?100:$iddata[-4]-$iddata[-3];			
		}

		my $hitname = $blastx_results_hashref->{$snpid}->[0];
		my $frame = $blastx_results_hashref->{$snpid}->[1];
#		print join("****", @{$blastx_results_hashref->{$snpid}}),"\n" if defined $frame;
		($length_left_flanking, $length_right_flanking) = ($length_right_flanking, $length_left_flanking) if $frame<0;
		my $seq = $frame<0?rev_comp($fastafile_hash{$snpid}):$fastafile_hash{$snpid};
		$altbase =~ tr/[ATGC]/[TACG]/ if $frame<0;
# 		print $snpid, "\n",$seq,"\n", $frame,"\n", translate($seq, $frame),"\n";
		my $ref_codon = get_codon($seq, $frame, $length_left_flanking);
		substr($seq, $length_left_flanking, 1) = $altbase;
		my $alt_codon = get_codon($seq, $frame, $length_left_flanking);
		print STDERR $snpid, ' $ref_codon: ', $ref_codon, "\t", '$alt_codon: ', $alt_codon, "\n" if $snpid eq 'Ex_mira1_rep_c115305_95195302_95196125_95195380_95195380_[A/G]';
		$return_hash{$snpid} = $hitname . ':'. codon_table($ref_codon) . '->' . codon_table($alt_codon);
													 
	}
	return %return_hash;
}

sub get_codon
{
	my ($seq, $frame, $left_length) = @_;
	my $codon_start = $left_length - ($left_length - abs($frame) + 1)%3;
	my $codon = substr $seq, $codon_start, 3;
#	if ((length $codon) < 3){print 'seq length: ', length $seq, "\n", 'left length: ', $left_length,"\n"}
	return $codon;
}


sub rev_comp
{
	my $seq = shift;
	$seq = reverse $seq;
	$seq =~ tr/[ATGCatgc]/[TACGtacg]/;
	return $seq;
}

sub read_fasta_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file\n";
	my $id;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		if (/^>(\S+)/)
		{
			$id = $1;
			$return_hash{$id} = '';
			next;
		}
		else
		{
			$return_hash{$id} .= $_;
		}
	}
	close IN;
	return %return_hash; 
}

sub translate
{
	my ($seq, $frame) = @_;
	#if ($frame < 0){$seq = rev_comp($seq)}
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
		'ATA'=>'I',#Isoleucine
		'ATC'=>'I',#Isoleucine
		'ATT'=>'I',#Isoleucine
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

