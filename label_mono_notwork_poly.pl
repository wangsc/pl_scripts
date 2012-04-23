#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<mono snps in 9k>
<not worked snps in 9k>
<polymorphic snps in 9k>
<PrivKSU_WheatCons_9k_11497518_A_formatted_blat_all_snps_50bp_oligo.out>
<acc_snp_freq>
);

die $sUsage unless @ARGV >= 5;
my ($mono_file, $notwork_file, $poly_file, $blat_file, $snp_freq_file, $out_file) = @ARGV;

my %lables = read_mono_notwork_poly($notwork_file, $mono_file, $poly_file);

my %blat= read_blat($blat_file);

open (IN, "$snp_freq_file") or die;
while(<IN>)
{
	chomp; next if /^\s+$/;
	my @t = split/\t/,$_;
	my $id = join(":", @t[1,2]);
	if(exists $blat{$id})
	{
		print $_,"\t", $lables{$blat{$id}}, "\n" if exists $lables{$blat{$id}};
	}
	else
	{
		print $_,"\t", "NA", "\n";
	}
}

sub read_mono_notwork_poly
{
	my @files = @_;
	my %return;
	foreach my $index (0..$#files)
	{
		open(IN, "$files[$index]") or die;
		while(<IN>)
		{
			chomp;
			next if /^\s+$/;
			chomp; 
			my @t=split/\t/,$_;
			$return{$t[0]} = $index;
		}
		close IN;
	}
	return %return;
}

sub read_blat
{
	my $file = shift;
	my %return;
	open(IN, $file) or die;
	while(<IN>)
	{
		my @t = split/\t/,$_;
		my $id = $1 if $t[1]=~/(\S+:\S+)\:/;
		$return{$id} = $t[0] unless exists $return{$id};
	}
	close IN;
	return %return;
}





