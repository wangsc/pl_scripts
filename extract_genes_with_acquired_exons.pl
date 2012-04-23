#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage = qq(
perl $0
<acquired exon combined fasta>
<isoform gff3>
<contig fasta>
);
die $sUsage unless @ARGV >= 3;

my ($ac_exon_fasta, $iso_gff, $contig_fasta) = @ARGV;

my %ctg_seq = read_fasta($contig_fasta);
my %gene_exon = read_exon_fasta($ac_exon_fasta);
my ($gene_region, $gene_ctg) =  read_iso_gff($iso_gff, \%gene_exon);
foreach my $gene (keys %gene_exon)
{
	my $contig = $gene_ctg->{$gene};
	my ($start, $end) = @{$gene_region->{$gene}};
	my @new_pos;
	foreach(@{$gene_exon{$gene}})
	{
		my @t=split /_/,$_;
		@t=map{$_-$start}@t;
		push @new_pos, join('_', @t);
	}
	print '>', join('_', ($gene, @new_pos) ),"\n";
	print substr($ctg_seq{$contig}, $start-1, ($end-$start+1)),"\n";
}

sub read_exon_fasta
{
	my $file = shift;
	open (IN, $file) or die $!;
	my %return;
	while (<IN>)
	{
		chomp;
		next unless /^>/;
		if(/>(\S+)/)
		{
			my $full_name = $1;
			my ($gene, $region) = $full_name =~ /(.*?)_(\S+)/;
			push @{$return{$gene}}, $region;
		}
	}
	close IN;
	return %return;
}

sub read_iso_gff
{
	my $file = shift;
	my $exon_ref = shift;
	open (IN, $file) or die $!;
	my %return;
	my %gene_ctg;
	while (<IN>)
	{
		next if /^\s+$/;
		my ($contig, $gene_id) = /^(\S+).*ID=(.*?)\-/;
		$gene_ctg{$gene_id} = $contig;
		my @t = split /\t/, $_;
		push @{$return{$gene_id}}, @t[3,4];
	}
	close IN;
	foreach (keys %return)
	{
		my @array = sort{$a<=>$b}@{$return{$_}};
		$return{$_} = [@array[0, -1]]
	}
	return(\%return, \%gene_ctg);
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
			print 'fasta id: ', $id ,"\n" if $debug; $debug = 0;
			$return_hash{$id} = '';
			next;
		}
		$return_hash{$id} .= $_;
	}
	close IN;
	return %return_hash;
}







