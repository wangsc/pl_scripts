#!/usr/bin/perl -w
use strict;
use Statistics::PointEstimation;
use Statistics::TTest;


my $sUsage = qq(
perl $0
<All_expresisons.out>
<output file>
);

die $sUsage unless @ARGV == 2;
my ($expr_file, $out_file) = @ARGV;
open (OUT, ">$out_file") or die "can't open file $out_file\n";
my %expr = get_expr($expr_file);
my $ttest = new Statistics::TTest;
$ttest->set_significance(95);


my @files = <*_rmdup_allele_freq.out.pvalue2_adjusted_addeffect_balanced>;
my %record;

foreach my $f (@files)
{
	my $id = $1 if $f =~ /(\S+)_rmdup_/;
	open (IN, $f) or die "can't open file $f\n";
	while (<IN>)
	{
		# cluster_id,pvalue_A,pvalue_B,pvalue_D,count_A,count_B,count_D,test_AB_D,A_effect,B_effect,D_effect
		# BobWhite_mira1_s66937,NA,NA,0.224317991371872,NA,NA,43/53,0.46988,A_NA,B_NA,D_BAL
		chomp;
		next if /^cluster_id/;
		my @t = split /,/, $_; 
		next if /A_NA/ and /B_NA/ and /D_NA/;
		if(/SILENCE/)
		{
			push @{$record{$t[0]}{SIL}}, $id
		}
		else
		{
			push @{$record{$t[0]}{BAL}}, $id
		}
	}
	close IN;
}

foreach my $ctg (keys %record)
{
	my @bal_accs = exists $record{$ctg}{BAL}?@{$record{$ctg}{BAL}}:();
	my @sil_accs = exists $record{$ctg}{SIL}?@{$record{$ctg}{SIL}}:();
	next unless @bal_accs > 1 and @sil_accs > 1;
	my @bal_exp = map{$expr{$ctg}{$_} if exists $expr{$ctg}{$_}}@bal_accs;
	my @sil_exp = map{$expr{$ctg}{$_} if exists $expr{$ctg}{$_}}@sil_accs;
	next unless @bal_exp > 1 and @sil_exp > 1;
	#print "Got one!\n";
	#print join("\t", @bal_exp), "\n", join("\t", @sil_exp), "\n"; exit;
	$ttest->load_data(\@bal_exp, \@sil_exp);
	my $alpha = $ttest->alpha * 2;
	print OUT $ctg, "\t", join(",", @bal_exp), "\t", join(",", @sil_exp), "\t", $alpha, "\t (have the same means is)", $ttest->null_hypothesis(), "\n";
}
close OUT;

sub get_expr
{
	my $file = shift;
	open (IN, $file) or die "can't open file $file \n";
	my %return;
	my $n = 0;
	my @accs;
	while (<IN>)
	{
		chomp;
		my @t = split /\t/, $_;
		if($n==0)
		{
			@accs = @t;
			$n++; 
			next
		}
		map{$return{$t[0]}{$accs[$_]}=$t[$_]}1..$#t;
	}
	close IN;
	return %return;
}


