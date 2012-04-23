#!/usr/bin/perl -w
use strict;

my $sUsage = "perl <final report> <deletion line file>\n";
die $sUsage unless @ARGV >= 2;
my($report_file, $del_list) = @ARGV;

my %snp_index;
my %report = read_report_file($report_file, \%snp_index);
my %deletion = read_del_file($del_list);
my %detect_del;
foreach my $snpid (keys %report )
{
	$detect_del{$snpid} = [] unless exists $detect_del{$snpid};
	foreach my $sample (keys %deletion)
	{
		push @{$detect_del{$snpid}}, $deletion{$sample} unless exists $report{$snpid}{$sample};
	}
	push @{$detect_del{$snpid}}, $snp_index{$snpid};
}
map{print $_, "\t", join("\t", @{$detect_del{$_}}),"\n"} keys %detect_del;



sub read_report_file
{
	my $file = shift;
	my $snpindex = shift;
	my %return;
	open (IN, "$file") or die "can't open file $file\n";
	while(<IN>)
	{
		next if /^\s+$/;
		next unless /^wsnp/;
		chomp;
		# wsnp_AJ612027A_Ta_2_1,5904359121_R01C01,A,A,0.9289,1,A,A,A,A,A,A,-,-,1.0000
		my @line_data = split /,/, $_;
		next if $line_data[2] eq '-';
		$return{$line_data[0]}{$line_data[1]} = 1;
		$snpindex->{$line_data[0]} = $line_data[5];
	}
	close IN;
	return %return;
}


sub read_del_file
{
	my $file = shift;
	my %return;
	open (IN, "$file") or die "can't open file $file\n";
	while(<IN>)
	{
		# Dt1AS	5904367078_R11C02
		next if /^\s+$/;
		chomp;
		my ($chr, $sample) = split /\s+/, $_;
		$return{$sample} = $chr;
	}
	close IN;
	return %return;
}