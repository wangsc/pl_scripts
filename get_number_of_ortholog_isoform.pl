#!/usr/bin/perl -w
use strict;
my $sUsage = qq(
perl $0
<ortholog.list for rice>
<pasa_chr3A_cdna.pasa_assemblies.denovo_transcript_isoforms.gff3>
<rice.gff3>
);
die $sUsage unless @ARGV >= 3;

my ($ortho_file, $isoform_gff, $rice_gff) = @ARGV;

my ($iso_arrayref, $rice_hashref) = read_ortho_file($ortho_file, $rice_gff);
my %iso_gene = read_gff_file($isoform_gff);
my $num_iso_wheat = scalar (@$iso_arrayref);
my %unique_gene = map{$iso_gene{$_}, 1} @$iso_arrayref;
my $num_gene_wheat = scalar (keys %unique_gene);
my $num_gene_rice = scalar (keys %$rice_hashref);
my $num_iso_rice = 0;
map{$num_iso_rice += $_} values %{$rice_hashref};
print '$num_gene_wheat ', $num_gene_wheat,"\n";
print '$num_iso_wheat ', $num_iso_wheat,"\n";
print '$num_gene_rice ', $num_gene_rice,"\n";
print '$num_iso_rice ', $num_iso_rice,"\n";

##########
sub read_ortho_file
{
	my $file = shift;
	my $rice_gff = shift;
	my @return_array;
	my %record;
	my %iso;
	open (IN, "$file") or die "$! $file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my @data  = split /\t/, $_;
		push @return_array, $data[0];
		my $id = $data[1];
		if ($id =~ /\|/)
		{
			$id = $1 if $data[1] =~ /^(.*?)\|/; 
			$iso{$id} = 1; 
			$id = $1 if $id=~/(.*)\.\d/;
			$record{$id}=0
		}
		else
		{
			$iso{$id} = 1; 
			$id = $1 if $id=~/(.*)\.\d/;
			#print '** ', $id,"*\n";
			$record{$id}=0
		}
	}
	close IN;
	open (RICE, "$rice_gff") or die $!;
	while (<RICE>)
	{
		next if /^\s+$/ or /^\#/;	
		my @t = split /\t/, $_;
		next unless $t[2] =~ /mRNA/i;
		my $line = $_;
		if ($rice_gff =~ /rice/i)
		{
			next unless /Alias/;
			my $id = $1 if $line =~ /Alias=(.*?)\.\d/;
			$record{$id}++ if exists $record{$id};
		}
		else
		{

			my $id = $1 if $line =~ /Parent=(\S+)/;
			$record{$id}++ if exists $record{$id};
		}		
	}
	close RICE;
	return (\@return_array, \%record);
}

sub read_gff_file
{
	my $file = shift;
	my %iso_gene;
	my %iso_contig;
	open (IN, "$file") or die "$! $file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my ($gene, $iso_id) = $_=~/ID=(.*?)\-(.*?)\;/;
		$iso_gene{$iso_id} = $gene;
		my @data = split /\t/, $_;
		$iso_contig{$gene} = shift @data;
	}
	close IN;
	return %iso_gene;
}



