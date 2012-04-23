#!/usr/bin/perl -w
use strict;
my $sUsage = qq(
perl $0
<wr_acquired_exons.list>
<wr_acquired_exons.list>
<output file>
);

die $sUsage unless @ARGV >= 3;

my ($wr_file, $wb_file, $output) = @ARGV;
open (OUT, ">$output") or die;
my %wr = read_file_and_construct_vec($wr_file);
my %wb = read_file_and_construct_vec($wb_file);

foreach my $gene (keys %wr)
{
	next unless exists $wb{$gene};
	my ($wb_vec, $wb_start, $wb_end) = @{$wb{$gene}};
	my ($wr_vec, $wr_start, $wr_end) = @{$wr{$gene}};
	my ($start, $end) = (sort{$a<=>$b} ($wb_start, $wb_end, $wr_start, $wr_end))[0,-1];
	foreach ($start.. $end)
	{
		my $vb = vec($wb_vec, $_,1);
		vec($wb_vec, $_,1) = $vb==1?$vb:0;
		my $vr = vec($wr_vec, $_,1);
		vec($wr_vec, $_,1) = $vr==1?$vr:0;
	}
	my $comb_vec = $wb_vec & $wr_vec;
	my ($comb_start, $comb_end);
	foreach ($start..$end)
	{
		if(vec($comb_vec, $_,1)==1)
		{
			$comb_start = $_ unless defined $comb_start;
		}
		else
		{
			if (defined $comb_start)
			{
				$comb_end = $_ - 1;
				print OUT join("_", ($gene, $comb_start, $comb_end)),"\n" if ($comb_end-$comb_start+1) >= 30;
				undef $comb_start;				
			}
		}
	}
}


sub read_file_and_construct_vec
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while (<IN>)
	{
		chomp; 
		my @t=split /_/, $_;
		push @{$return{$t[0]}}, [@t[1,2]];
	}
	close IN;
	map{$return{$_} = construct_vec($return{$_})} keys %return;
	return %return;
}


sub construct_vec
{
	my $arrayref = shift;
	my $vec = '';
	my $max;
	my $min;
	my $debug =1 ;
	foreach (@$arrayref)
	{
		my @d = sort{$a<=>$b}@$_;
#		print '@d: ', join("\t", @d),"\n" if $debug; $debug=0;
		foreach ($d[0]..$d[1])
		{
			$max = $_ unless defined $max;
			$max = $_ if $_ > $max;
			$min = $_ unless defined $min;
			$min = $_ if $_ < $min;
			vec($vec,$_,1) = 0b1;
		}
	}
	foreach ($min..$max)
	{
		my $value = vec($vec, $_, 1);
		vec($vec, $_, 1) = $value==1?$value:0;
	}
	return [$vec, $min, $max];
}