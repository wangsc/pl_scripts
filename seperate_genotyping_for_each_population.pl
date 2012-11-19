#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<grouping_5.csv>
<genotyping file>
);
die $sUsage unless @ARGV >= 2;

my ($grp_file, $genotyping_file) = @ARGV;

my %grps = read_grp_file($grp_file);
my %file_handles = get_file_handles (values %grps);

open (IN, $genotyping_file) or die;

while (my $line = <IN>)
{
	if($line =~ /^lines/)
	{
		foreach my $fh(values %file_handles)
		{
			print {$fh} $line;
		}
		next;
	}
	
	my $id = $1 if $line=~/^(\S+)/;
	my $fh = $file_handles{$grps{$id}};
	print {$fh} $line;
}

sub read_grp_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while (<IN>)
	{
		chomp;
		next if /^acc/;
		my @t = split /,/,$_; 
		$return{$t[2]} = $t[-1];
	}
	close IN;
	return %return;
}

sub get_file_handles
{
	my @files = @_;
	@files = unique(@files);
	my %return_hash;
	foreach (@files)
	{
		local *FH;
		my $f = $_. "_genotyping.out";
		open (FH, ">$f") or die "can't open file $f\n";
		$return_hash{$_} = *FH{IO};
	}
	return %return_hash;
}

sub unique
{
	my %h = map{$_, 1}@_;
	return keys %h;
}

