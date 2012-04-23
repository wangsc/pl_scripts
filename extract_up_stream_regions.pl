#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage = qq(
# This script will extract up_stream regions of genes and store the sequences in a fasta file
Usage: 
	perl $0
	<gff3 file>
	<genome fasta file>
	<length of up_stream region, 5000 or ... >
	<fasta file stor storing up_stream sequences>
);
die $sUsage unless @ARGV >= 4;
my($gff_file, $genome_file, $length_up, $out_fasta) = @ARGV;

# Main
my $genome_obj = Bio::DB::Fasta->new($genome_file);

my %up_stream_regions = parse_gff_file($gff_file, $length_up, $genome_obj);

output(\%up_stream_regions, $genome_obj, $out_fasta);

# Subroutines
sub parse_gff_file
{
	my ($gff, $len_up, $genome) = @_;
	my %return_hash;
	open (IN ,"$gff") or die "can't open file $gff\n";
	while(<IN>)
	{
		next if /^#/;
		chomp;
		my @line_data = split /\t/, $_;
		my ($chr_id, $type, $gen_start, $gen_end, $strand, $description) = @line_data[0, 2, 3, 4, 6, 8];
		my $gene_id = $1 if $description =~ /ID=(.*?)\;/;
		next unless $type eq 'gene';
		my ($up_start, $up_end);
		if ($strand eq '+')
		{
			$up_start = $gen_start - $len_up + 1;
			$up_start = $up_start>0?$up_start:0;
			$up_end = $gen_start;
			push @{$return_hash{$chr_id}}, [$gene_id, $up_start, $up_end];
		}
		else
		{
			$up_start = $gen_end + $len_up - 1;
			$up_start = $up_start>$genome->length($chr_id)?$genome->length($chr_id):$up_start;
			$up_end = $gen_end;
			push @{$return_hash{$chr_id}}, [$gene_id, $up_start, $up_end];
		}
	}
	close IN;
	return %return_hash;
}

sub output
{
	my ($up_str_ref, $genome, $out_file) = @_;
	open (OUT, ">$out_file") or die "can't open file $out_file \n";
	foreach my $chr_id (keys %$up_str_ref)
	{
		foreach (@{$up_str_ref->{$chr_id}})
		{
			my ($gene_id, $start, $end) = @$_;
			print OUT ">", join("_",($chr_id, $gene_id)),"\n",
							$genome->seq($chr_id, $start => $end),"\n";
		}
	}
	close OUT;
}











