#!/usr/bin/perl -w
use strict;
use File::Basename;

# $mapping_out_dir $all_vcf $geno_out
my ($mapping_out_dir, $output) = @ARGV;

my $min_dep = 1;
my @allele_freq_files = <$mapping_out_dir/*allele_freq.out>;

my @acc_names = get_acc_list(@allele_freq_files);
my %snp_genotype = read_allele_freq_files($min_dep, @allele_freq_files);

# Output
open (OUT, ">$output") or die $!;
print OUT "SNP\t", join("\t", @acc_names),"\n";
foreach my $snp (keys %snp_genotype)
{
	my @tmp = ();
	push @tmp, $snp;
	foreach my $acc (@acc_names)
	{
		my $genotype = exists $snp_genotype{$snp}{$acc}?$snp_genotype{$snp}{$acc}:"NA";
		push @tmp, $genotype
	}
	print OUT join("\t", @tmp), "\n";
}
close OUT;

# Sub
sub get_acc_list
{
	my @fs = @_;
	my @return;
	foreach (@fs)
	{
		my $acc = $1 if basename($_) =~ /(\S+)_allele_freq/;
		push @return, $acc
	}
	return @return;
}

sub read_allele_freq_files
{
	my ($mindep, @files) = @_;
	my %return;
	foreach my $file (@files)
	{
		my $acc = $1 if basename($file) =~ /(\S+)_allele_freq/;
		open (IN, $file) or die;
		while (<IN>)
		{
			chomp; 
			my @t = split /\t/,$_; 
			#Cheyenne_287213	35	C_4	T_2
			my $snp = join(":", @t[0, 1]);
			my @alleles;
			foreach (2..$#t)
			{
				my ($g, $count) = split /_/, $t[$_];
				push @alleles, $g if $count >= $mindep;
			}
			my $genotype = join("", sort{$a cmp $b}@alleles);
			$return{$snp}{$acc} = $genotype;
		}
		close IN;
	}
	return %return;
}