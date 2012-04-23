#!/usr/bin/perl -w
use strict;

my $file = shift or die "perl $0 genotype_file\n";

open (IN, $file) or die;
while(<IN>)
{
	if(/^SNP/){ print $_;	next}
	chomp;
	my @transformed_data = transform($_);
	print join("\t", @transformed_data),"\n" if @transformed_data;
}
close IN;

sub transform
{
	my $line = shift;
	#print STDERR $line, "\n";
	my @data = split /\t/, $line;
	my $snp_id = shift @data;
	my %geno_count;
	map{$geno_count{$_}++}@data;
	my @genotypes;
	map{push @genotypes, $_ unless /NA/}keys %geno_count;
	#print join("*\n", @genotypes),"**\n"; exit 0;
	return () unless @genotypes;
	my %alleles;
	foreach (@genotypes)
	{
		next if /NA/;
		map{$alleles{$_}=1}(split//,$_);
	}
	my @alls = sort{$a cmp $b}keys %alleles;
	unless(@alls){print $line,"\n"; die "no alleles\n"}
	my %value;
	$value{$alls[0]} = 1;
	$value{$alls[1]} = -1 if defined $alls[1];
	if(@genotypes == 3)
	{		
		foreach (@genotypes)
		{
			next unless length $_ == 2;
			my $gn = join("", sort{$a cmp $b}(split//,$_));
			next unless $gn eq join("", @alls);
			$value{$_} = $geno_count{$alls[0]}>$geno_count{$alls[1]}?$value{$alls[1]}:$value{$alls[0]};
		}		
	}
	elsif (@genotypes == 2)
	{
		foreach (@genotypes)
		{
			next unless length $_ == 2;
			my $gn = join("", sort{$a cmp $b}(split//,$_));
			next unless $gn eq join("", @alls);
			$value{$_} = exists $geno_count{$alls[0]}?$value{$alls[1]}:$value{$alls[0]};
		}		
	}
	else
	{		
	}
	
	my @transformed_data = map{exists $value{$_}?$value{$_}:"NA" }@data;
	return ($snp_id, @transformed_data);	
}