#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<homeolog_D_AB_allele file>
<homeolog_A_B_allele file>
<output>
);
die $sUsage unless @ARGV >= 3;

my ($dab_file, $ab_file, $output) = @ARGV;
open (OUT, ">$output") or die "can't open file $output \n";
my %dab = read_dab_file($dab_file);
my %ab = read_ab_file($ab_file);

my %fab;
foreach my $id (keys %dab)
{
	my $ab_combo = $dab{$id}->[1];
	$ab_combo = join("",  sort{$a cmp $b}(split /,/, $ab_combo));	
	if(length $ab_combo == 1)
	{
		$fab{$id} = [$ab_combo, $ab_combo];
		print OUT join("\t", ($id, $dab{$id}->[0], $ab_combo, $ab_combo)), "\n";
		next
	}
	next unless exists $ab{$id};
	my $ab_allele = join("", sort{$a cmp $b}@{$ab{$id}});
	if($ab_combo eq $ab_allele)
	{
		$fab{$id} = $ab{$id};
		print OUT join("\t", ($id, $dab{$id}->[0], @{$ab{$id}})), "\n";
	}
	elsif($ab_allele =~ /\*/)
	{
		my @tmp = split //, $ab_combo;
		my $ind = $ab{$id}[0] eq '*'?1:0;
		@tmp = ($tmp[$ind] eq $ab_allele)?@tmp:reverse @tmp;
		$fab{$id} = [@tmp];
		print OUT join("\t", ($id, $dab{$id}->[0], @tmp)), "\n";
	}	
}

close OUT;

sub read_dab_file
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while (<IN>)
	{
# Kukri_mira1_c3611       717     1556    T       C
# Excalibur_mira1_c547    1030    1458    T       G
# Kukri_mira1_c48906      138     1478    G       T,G
# Kukri_mira1_c6545       1182    1747    A       T
		chomp;
		my @t= split /\s+/,$_; 
		$return{join(":", @t[0,1])}=[@t[3,4]];
	}
	close IN;
	return %return;
}

sub read_ab_file
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while (<IN>)
	{
# Excalibur_mira1_c11667:2396     T       T
# Kukri_mira1_c18219:812          T       T
# Excalibur_mira1_c2714:2283      A       A
# Kukri_mira1_rep_c102148:444     *       C
		
		chomp;
		my @t= split /\s+/,$_; 
		$return{$t[0]}=[@t[1,2]];
	}
	close IN;
	return %return;	
}