#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
Usage:
perl $0
<ref sequence fasta file>
<replacement information file, each line like: 9065 T C>
<SNP position in fasta file (starting from 1), like 21>;
<output filename>

);
die $sUsage unless @ARGV >= 4;
my($fasta_file, $replace_infor, $snp_pos, $out_file) = @ARGV;

my %fasta_seqs = read_fasta_file($fasta_file);
my %replace_info = read_replace_infor_file($replace_infor);
replace_and_output(\%fasta_seqs, \%replace_info, $snp_pos, $out_file);

# Subroutines
sub read_fasta_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file\n";
	my $id;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		if (/^>/)
		{
			$id = $_;
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

sub read_replace_infor_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my ($id, $ref, $alt) = split /\s+/,$_;
		$return_hash{$id} = [$ref, $alt];
	}
	close IN;
	return %return_hash;
}

sub replace_and_output
{
	my ($fasta_hashref, $replaceinfo_hashref, $snp_pos ,$out_file) = @_;
	open (OUT, ">$out_file") or die "can't open file $out_file\n";
	foreach my $id (keys %$replaceinfo_hashref)
	{
		my $seq = $fasta_hashref->{'>'. $id};
		my $ref = substr($seq, $snp_pos-1, 1);
		map{substr($seq, $snp_pos-1, 1) = $_ unless $ref eq $_} @{$replaceinfo_hashref->{$id}};
		print OUT '>' . $id,"\n", $seq,"\n";
	}
	close OUT;
}





