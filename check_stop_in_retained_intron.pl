#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $sUsage = qq(
perl $0
<pasa_db_asmbl_splice_db_dump.csv>
<pasa_assemblies.denovo_transcript_isoforms.gff3>
<pasa_assemblies transcript isoform fasta file>
<blastx result file>
);
die $sUsage unless @ARGV >= 4;
my ($db_dump, $iso_gff, $asmbl_fasta_file, $blast_result_file) = @ARGV;
my %asmbl_fasta = read_fasta($asmbl_fasta_file);
my %id_reading_frame = get_reading_frame($blast_result_file);
my %asmbl_ids = read_db_dump($db_dump);

print join("\t",qw(Type Num_asmbl Num_asmbl_with_stop Num_intron Num_intron_with_stop_codon)), "\n";
foreach my $feature (keys %asmbl_ids)
{
	print STDERR '$feature ',$feature,"\n";
	my $count_id = 0;
	my $count_id_stop = 0;
	my $count_frag = 0;
	my $count_stop = 0;
	my $ori_feature_pos = $asmbl_ids{$feature};
	my %feature_positions = transform_cordinate($ori_feature_pos, $iso_gff);
	foreach my $id (keys %feature_positions)
	{
		next unless exists $id_reading_frame{$id};
		my ($query_length, $frame) = @{$id_reading_frame{$id}};
		$count_id++;
		my $seq = $asmbl_fasta{$id};
		my $flag = 0 ;
		foreach (@{$feature_positions{$id}})
		{
			$count_frag++;
			my ($start, $end) = @$_;
			print STDERR join("\t",($start, $end) ),"\n" if $id eq 'asmbl_1514';
			if($frame >0)
			{
				my $remain = ($start - 1 -$frame)%3;
				my $sub_start = $start - $remain -1;
				my $sub_seq = substr($seq, $sub_start, ($end-$start+1+$remain));
				if(($sub_start+$end-$start+1)>length $seq){print '$id ', $id,"\n";die}
				my $have_stop = check_stop_codon($sub_seq);
				$count_stop++ if $have_stop;
				$flag = 1 if $have_stop;
			}
			else
			{
				my $remain = ($query_length - $end - abs($frame))%3;
				my $sub_seq = substr($seq, $start, ($end-$start+1+$remain));
				$sub_seq = reverse_complemantary($sub_seq);
				my $have_stop = check_stop_codon($sub_seq);
				$count_stop++ if $have_stop;
				$flag = 1 if $have_stop;
			}
		}
		$count_id_stop++ if $flag;
	}
	print join("\t", ($feature, $count_id, $count_id_stop, $count_frag, $count_stop)),"\n";
}


# Subroutines
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

sub read_db_dump
{
	my ($db_file) = @_;
	my $feature = 'retained_intron';
	my %return;
	open (IN, $db_file) or die $!;
	while(<IN>)
	{
		s/\"//g;
		next unless /$feature/;
		my @t = split /,/, $_;
		my $feature = $t[5];
		push @{$return{$feature}{$t[1]}}, [@t[2, 3]];
	}
	close IN;
	return %return;
}

sub transform_cordinate
{
	my ($asmbl_ref, $gff) = @_;
	my %return;
	my ($asmbl_pos, $strand_ref) = read_iso_gff($gff);
	my $debug = 1;
	foreach my $id (keys %$asmbl_ref)
	{
		foreach (@{$asmbl_ref->{$id}})
		{
			my ($start, $end) = @$_;
			print STDERR join("\t",($start, $end) ),"\n" if $id eq 'asmbl_1514';
			foreach (@{$asmbl_pos->{$id}})
			{
				next unless $start >= $_->[0] and $_->[1]>=$end;
				my $ind = ($strand_ref->{$id} eq '-')?1:0;
				my $offset = $_->[$ind];
				$start = abs($offset-$start) + $_->[2];
				$end = abs($offset-$end) + $_->[2];
				last;
			}
			($start, $end) = ($end, $start) if $start > $end;
			push @{$return{$id}}, [$start, $end];
		}
	}
	return %return;	
}

sub read_iso_gff
{
	my $file = shift;
	open (IN, $file) or die $!;
	my %return;
	my %strand;
	while (<IN>)
	{
		my $id = $1 if /Target=(\S+)/;
		my ($start, $end) = ($1, $2) if /Target=\S+\s(\d+)\s+(\d+)/;
		next unless defined $id;
		my @t = split /\t/, $_;		
		push @{$return{$id}}, [@t[3, 4], $start, $end];
		$strand{$id} = $t[6];
	}
	close IN;
	return(\%return, \%strand);

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
		if (/^>.*(asmbl_\d+)/)
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
