#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <acquired_exons_gene_asmbl_region.out> <genewise_output.gff>\n";
die $sUsage unless @ARGV >= 2;

my ($ac_region_file, $genewise_gff) = @ARGV;

my %ac_regions = read_ac_region_file($ac_region_file);
my %genewise = read_genewise_gff($genewise_gff);

my %results = map{$_, 0 }keys %ac_regions;
foreach my $exon_id (keys %ac_regions)
{
	my ($start, $end) = (split/_/, $exon_id)[1,2];
	my @asmbles = keys %{$ac_regions{$exon_id}};
	my $overlap = 0;
	foreach my $asmbl (@asmbles)
	{
		next unless (exists $genewise{$asmbl});
		$overlap = check_overlap($start, $end, $genewise{$asmbl});
		last if $overlap == 1;
	}
	$results{$exon_id} = $overlap;
}

foreach (keys %results)
{
	print $_, "\t", $results{$_}, "\n";
}



#
sub check_overlap
{
	my ($start, $end, $arr_ref) = @_;
	my $flag = 0;
	foreach (@$arr_ref)
	{
		my ($arr_start, $arr_end) = @$_;
		if ( ($start>=$arr_start and $start<=$arr_end) or ($arr_start>=$start and $arr_start<=$end) )
		{
			$flag=1;
			last;
		}
	}
	return $flag;
}


sub read_ac_region_file
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while(<IN>)
	{
		chomp; 
		my @t = split/\s+/,$_;
		my $asmbl = $1 if /(asmbl_\d+)/;
		$return{$t[0]}{$asmbl} = 1;
	}
	close IN;
	return %return;
}

sub read_genewise_gff
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while(<IN>)
	{
		chomp;
		next unless /match/;
		my @t = split /\t/, $_;
		push @{$return{$t[0]}}, [sort {$a<=>$b} @t[3,4]];
	}
	close IN;
	return %return;
}