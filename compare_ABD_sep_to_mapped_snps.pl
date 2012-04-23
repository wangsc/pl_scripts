#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<ABD_coverage_percent.out>
<mapped_9k_snps>
);
die $sUsage unless @ARGV >= 2;
my($coverage_file, $mapped_chr_file) = @ARGV;

my %mapped_snps = get_mapped_chr($mapped_chr_file);
my %coverage = read_coverage_file($coverage_file);

my ($min, $max) = (0.5, 0);

while ($min < 1)
{
	$max = 0;
	while ($max < 1)
	{
		my %assigned_chr = assign_chr(\%coverage, $min, $max);
		my @comp_results = compare(\%mapped_snps, \%assigned_chr);
		print join("\t", ($min, $max, @comp_results)), "\n";
		$max += 0.1;
	}
	print "-"x30, "\n";
	$min += 0.1;
}

# Subroutines
sub read_coverage_file
{
	my $file = shift;
	my %return;
	open (IN, $file ) or die $!;
	while (<IN>)
	{
		chomp; 
		my @t = split /\s+/, $_;
		$return{$t[0]} = [@t[1,2]];
	}
	close IN;
	return %return;
}

sub compare
{
	my ($snp_hashref, $abd_hashref) = @_;
	my ($total, $correct, $wrong) = (0, 0, 0);
	my $debug =1 ;
	foreach my $ctg (keys %$abd_hashref)
	{
		next unless exists $snp_hashref->{$ctg};
		$total++;
		if ($snp_hashref->{$ctg} eq $abd_hashref->{$ctg})
		{
			$correct++ ;
			print STDERR join("\t", ($ctg, $snp_hashref->{$ctg}, $abd_hashref->{$ctg})), "\n";
		}
		else
		{
			$wrong++
		}
	}
	return ($total, $correct, $wrong);
}


sub assign_chr
{
	my ($hashref, $min, $max) = @_;
	my %return;
	foreach my $id (keys %$hashref)
	{
		my ($ab_cov, $d_cov) = @{$hashref->{$id}};
		if($ab_cov > $d_cov)
		{
			$return{$id} = "AB" if $ab_cov >= $min and $d_cov <= $max;
			#print STDERR $id, "\tAB", "\n" if $min == 0.7;
		}
		elsif ($d_cov > $ab_cov)
		{
			$return{$id} = "D" if $ab_cov <= $max and $d_cov >= $min;
			#print STDERR $id, "\tD", "\n" if $min == 0.7;
		}
	}
	return %return;
}


sub get_mapped_chr
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	my $debug = 1;
	while (<IN>)
	{
		#1A      Ex_c10657_17376448
		chomp;
		my @t = split /\s+/, $_;
		my $ctg_id = $1 if $t[1] =~ /^(\S+)_\d+$/;
		my $chr = $t[0];
		my $transformed = transfrom_name($ctg_id);
		next if $transformed eq 'none';
		#print STDERR 'transform: ', $transformed, "\n" if $debug; $debug=0;
		$return{$transformed} = ($chr=~ /D/i)?"D":"AB";
	}
	close IN;
	return %return;
}

sub transfrom_name
{
	
	my $name = shift; # Ex_c10657 to Excalibur_mira1_c10657
	my %alias = (
			"Ra","RAC875",
			"Ku","Kukri",
			"JD","BobWhite",
			"Ex","Excalibur",
			"CAP12","CAP12",
			"JG","Jagger",
			"CAP7","CAP7",
			"CAP11","CAP11",
			"CAP8","CAP8"							
	);
	
	my @t = split /_/, $name;
	return 'none' unless exists $alias{$t[0]};
	return join("_", ($alias{$t[0]}, "mira1", @t[1..$#t]));		
}
