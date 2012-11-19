#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<All_Asso_MAF0.04_MapPosition_4-10-2012.csv>
<pairwise_ld.out>
<additional list, optional>
);

die $sUsage unless @ARGV >= 2;
my ($asso_file, $ld_file, $additional_file) = @ARGV;
my $ld_cutoff = 0.6;

my ($mapped_snps_ref, $unmapped_snps_ref) = read_asso_file($asso_file);
get_additional($additional_file) if defined $additional_file;
map{delete($mapped_snps_ref->{$_}) if exists $unmapped_snps_ref->{$_}}keys %$mapped_snps_ref;

my %related;
open (IN, $ld_file) or die "can't open file $ld_file\n";
while (<IN>)
{
	chomp;
	next unless /\S/;
	# wsnp_Ex_c3478_6369892   wsnp_Ra_c17221_26052231 0.00333653107913019
	my @t = split /\s+/, $_;
	next if exists $unmapped_snps_ref->{$t[0]} and exists $unmapped_snps_ref->{$t[1]};
	next if exists $mapped_snps_ref->{$t[0]} and exists $mapped_snps_ref->{$t[1]};
	if(exists $unmapped_snps_ref->{$t[0]})
	{
		push @{$related{$t[0]}}, [@t[1,2]] if $t[2] >= $ld_cutoff;
	}
	else
	{
		push @{$related{$t[1]}}, [@t[0,2]] if $t[2] >= $ld_cutoff;
	}
}
close IN;

# output
foreach my $snp (keys %$unmapped_snps_ref)
{
	unless (exists $related{$snp})
	{
		print $snp, "\t No_linked\n";
		next;
	}
	my @linked_snps = sort{$a->[1] <=> $b->[1]}@{$related{$snp}};
	print join("\t", ($snp, @{$linked_snps[0]}, $mapped_snps_ref->{$linked_snps[0][0]})), "\n";	
}

#
sub read_asso_file
{
	my $file = shift;
	my %mapped;
	my %unmapped;
	open (IN, $file) or die;
	while (<IN>)
	{
		chomp;
		my @t = split /,/,$_;
		my ($snp, $chr, $pos, $id) = @t[2, 4, 7, 10];
		$snp = $id unless $id =~ /SNP\d+/;
		if($chr =~ /UM/i or (not $id =~ /SNP\d+/))
		{
			$unmapped{$snp} = 1;
		}
		else
		{
			$mapped{$snp} = join("\t", ($chr, $pos));
		}
	}
	close IN;
	return(\%mapped, \%unmapped)
}

sub get_additional
{
	my $file = shift;
	open (IN, $file) or die;
	while (<IN>)
	{
		chomp;
		$unmapped_snps_ref->{$_}=1
	}
}
