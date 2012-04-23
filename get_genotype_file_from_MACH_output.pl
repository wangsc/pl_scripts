#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0 
<snp file, ie. mach1.out.erate>
<imputated genotype, ie. mach1.out.mlgeno>
<imputation quality, ie. mach1.out.mlqc>
<genotype output>
);

die $sUsage unless @ARGV >= 4;
my ($snp_list, $geno_file, $qc_file, $out_file) = @ARGV;

my $qc_cutoff = 0; # cutoff for quality of imputation

my @snps = read_snp_list($snp_list);
my %impu_geno = read_geno_file($geno_file);
my %impu_qc = read_qc_file($qc_file);

open (OUT, ">$out_file") or die $!;
print OUT join("\t", ("<Marker>", @snps)), "\n";
foreach my $acc (keys %impu_geno)
{
	my @quality = @{$impu_qc{$acc}};
	my @genos = @{$impu_geno{$acc}};
	map{$genos[$_]="?" if $quality[$_]<$qc_cutoff }0..$#genos;
	map{$genos[$_]=~s/1/A/g;  $genos[$_]=~s/2/B/g; $genos[$_]=~s/\//:/;}0..$#genos;
	print OUT join("\t", ($acc, @genos)), "\n"
}
close OUT;

sub read_snp_list
{
	my $file = shift;
	open (IN, $file) or die;
	my @return;
	while(<IN>)
	{
		chomp;
		next if /^Marker/i;
		my @t=split/\s+/, $_;
		push @return, $t[0];
	}
	close IN;
	return @return;
}

sub read_geno_file
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while(<IN>)
	{
		# Above->Above ML_GENO 2/2 1/1 1/1 1/1 1/1 1/1 1/2 2/2 2/2 1/1 2/2 1/1 2/1
		next if /^\s+$/;
		chomp;
		my @t = split /\s+/,$_; 
		my $acc = $1 if /^(\S+)\-\>/;
		$return{$acc} = [@t[2..$#t]];
	}	
	close IN;
	return %return;
}

sub read_qc_file
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while(<IN>)
	{
		# Above->Above ML_GENO 0.5 0.9 1.00 0.34
		next if /^\s+$/;
		chomp;
		my @t = split /\s+/,$_; 
		my $acc = $1 if /^(\S+)\-\>/;
		$return{$acc} = [@t[2..$#t]];
	}	
	close IN;
	return %return;
}


