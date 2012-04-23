#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<minimum depth>
<accessions list>
<allele frequency files>
);
die $sUsage unless @ARGV >= 3;
my ($min_dep, $acc_list_file, @allele_freq_files) = @ARGV;

my @acc_names = read_acc_list_file($acc_list_file);
my %snp_genotype = read_allele_freq_files($min_dep, @allele_freq_files);

# Output
print "SNP\t", join("\t", @acc_names),"\n";
foreach my $snp (keys %snp_genotype)
{
	my @tmp = ();
	push @tmp, $snp;
	foreach my $acc (@acc_names)
	{
		my $genotype = exists $snp_genotype{$snp}{$acc}?$snp_genotype{$snp}{$acc}:"NA";
		push @tmp, $genotype
	}
	print join("\t", @tmp);
}

# Sub
sub read_acc_list_file
{
	my $file = shift;
	my @return;
	open (IN, $file) or die;
	while(<IN>)
	{
		chomp;
		push @return, $_;
	}
	close IN;
	return @return;
}

sub read_allele_freq_files
{
	my ($mindep, @files) = @_;
	my %return;
	foreach my $file (@files)
	{
		open (IN, $file) or die;
		while (<IN>)
		{
			chomp; 
			my @t = split /\t/,$_; 
			# TS0531	Cheyenne_287213	35	C_4	T_2
			my $acc = $t[0];
			my $snp = join(":", @t[1,2]);
			my @alleles;
			foreach (3..$#t)
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








