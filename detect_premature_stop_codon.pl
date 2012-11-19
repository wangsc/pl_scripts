#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SearchIO;

my $sUsage = qq(
perl $0
<cdna sequence fasta>
<blastx output>
<output file>
);
die $sUsage unless @ARGV >= 3;
my ($cdna_file, $blastx_out, $out_file) = @ARGV;

my $gnobj = Bio::DB::Fasta->new($cdna_file);
my $blastio = Bio::SearchIO->new("-file"=>$blastx_out, "-format" => "blast");
my %parse_blast = parse_blastx_results($blastx_out);
#exit; # debug
open (OUT, ">$out_file") or die;
foreach my $seq_id (keys %parse_blast)
{
	my ($pep_start_end, $frame, $query_pos) = @{$parse_blast{$seq_id}};
	
	my $seq = $gnobj->get_Seq_by_id($seq_id)->seq;
	
	if($pep_start_end->[0] == -1 and $pep_start_end->[0] == -1)
	{
		next;
	}
	elsif($pep_start_end->[0] != -1 and $pep_start_end->[0] != -1)
	{
		my ($start, $end) = sort{$a<=>$b} @{$pep_start_end};
		$seq = $gnobj->seq($seq_id, $start => $end);
	}
	elsif($pep_start_end->[0] == -1) # no start
	{
		if($frame < 0)
		{
			$seq = substr($seq, $pep_start_end->[1]-1);
		}
		else
		{
			$seq = substr($seq, 0, $pep_start_end->[1]);
		}		
	}
	else # no end
	{
		if($frame < 0)
		{
			$seq = substr($seq, 0, $pep_start_end->[0]);
			
		}
		else
		{
			$seq = substr($seq, $pep_start_end->[0]-1);
		}				
	}
	
	my $start = abs($frame) - 1;
	$start = 0 if $pep_start_end->[0] != -1;
	my $peptide = translate($seq, $frame, $start);
	#if ($seq_id =~ /S6297/)
	#{
	#	print STDERR "S6297 ", join("\t", (@{$pep_start_end}, $frame)), "\n";
	#	print STDERR $peptide, "\n";
	#}
#	print STDERR ">$seq_id", "\n", $peptide, "\n"; last;
	my $pep_len = length $peptide;
	my @stop_pos;
	while($peptide =~ /_/g)
	{
		my $pos = pos($peptide);
		push @stop_pos, $pos unless $pos >= $pep_len - 3;
	}
	if(@stop_pos)
	{
		print OUT join("\t", ($seq_id, @stop_pos)), "\n";
	}
	else
	{
		print OUT $seq_id, "\t", "No Stop_codon\n";
	}
	print STDERR $seq_id, "\t", $frame, "\n";
}
close OUT;


# Subroutines

sub parse_blastx_results
{
	my $file = shift;
	my %return_hash;
	my $searchio = Bio::SearchIO->new(-format => 'blast', file => "$file" );
	while (my $result = $searchio->next_result())
	{
		last unless defined $result;
		my $query_name = $result->query_name;
		my $query_length = $result->query_length;
		#next unless $query_length >= 300;
		my $hit = $result->next_hit();		
		next unless defined $hit;
		my $hit_name = $hit->name;
		my $hit_length = $hit->length;
		my $cip; 
		my $calp;
		my $id_len; 
		my $aligned_len;
		my $blastframe;
		my $num_hsp = 0;
		my @query_pos;
		my @subject_pos;
		while(my $hsp = $hit->next_hsp)
		{
			last unless defined $hsp;
			$num_hsp++;
			$blastframe = ($hsp->query->frame + 1) * $hsp->query->strand unless defined $blastframe;
			$id_len += $hsp->num_identical;
			$aligned_len += $hsp->length('query');
			
			push @query_pos, [$hsp->start('query'), $hsp->end('query')];
			push @subject_pos, [$hsp->start('subject'), $hsp->end('subject')];
		}
		#print STDERR join("\t", ($query_name, $hit_name, @{$subject_pos[0]})), "\n";
		next unless defined $aligned_len;
		$cip = 3*$id_len/$aligned_len;
		$calp = $aligned_len/($query_length);
		#print STDERR '$cip ', $cip,"\n", '$calp ', $calp,"\n";
		#push @{$return_hash{$query_name}}, $hit_name if $cip >=0.6 and $calp>=0.7; # $cip >=0.6 and $calp>=0.7
		next unless $cip >=0.6; # and $calp>=0.7;
		#print STDERR join("\t", ($query_name, $query_length, $hit_length, $aligned_len, $id_len, $cip, $calp, $hit_name)), "\n";
		
		my ($query_pep_start, $query_pep_end) = (-1, -1);
		foreach my $index (0..$#subject_pos)
		{
			if($subject_pos[$index][0] <= 10)
			{
				$query_pep_start = $query_pos[$index][0]
			}
			if($subject_pos[$index][1] >= $hit_length-10 and $subject_pos[$index][1] <= $hit_length+10)
			{
				$query_pep_end = $query_pos[$index][1];
			}
		} 
		
		$return_hash{$query_name} = [[$query_pep_start, $query_pep_end], $blastframe, [@query_pos] ];
		
	}
	return %return_hash;
}

sub rev_comp
{
	my $seq = shift;
	$seq = reverse $seq;
	$seq =~ tr/[ATGCatgc]/[TACGtacg]/;
	return $seq;
}

sub translate
{
	my ($seq, $frame, $start) = @_;
	if ($frame < 0){$seq = rev_comp($seq)}

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