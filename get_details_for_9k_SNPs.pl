#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <wheat csv file> <integrated SNP effect> <output file>\n";
die $sUsage unless @ARGV == 3;
my ($csv_file, $snp_effect_file, $out_file) = @ARGV;

my %csv = read_csv_file($csv_file);
my %snp_effect = read_snp_effect_file($snp_effect_file);
output(\%csv, \%snp_effect, $out_file);

sub read_csv_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file\n";
	while (<IN>)
	{
		next unless /^wsnp_/;
		chomp;
		my @line_data = split /,/, $_;
		my ($snp_id, $flanking_seq, $assay_design_id, $illum_id) = @line_data[0, 1, 13, 14];
		$return_hash{$snp_id} = [$flanking_seq, $assay_design_id, $illum_id];
	}
	close IN;
	return %return_hash;
}

sub read_snp_effect_file
{
	my $file = shift;
	open (IN, "$file") or die "can't open file $file\n";
	my %return_hash;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my @line_data = split /\t/, $_;
		my $id = shift @line_data;
		$return_hash{$id} = [@line_data];		
	}
	close IN;
	return %return_hash;
}

sub output
{
	my ($csv_ref, $snpeffect_ref, $out_file) = @_;
	open (OUT, ">$out_file") or die "can't open file $out_file\n";
	print OUT join("\t", qw(Locus_Name Flanking_Sequence Assay_Design_ID Ilmn_ID Blast_Hit AA_alter)),"\n";
	foreach my $snp_id (keys %{$csv_ref})
	{
	#	print $snp_id,"\n";
		print OUT join("\t", ($snp_id, @{$csv_ref->{$snp_id}}, @{$snpeffect_ref->{$snp_id}})),"\n";
	}
	close OUT;
}