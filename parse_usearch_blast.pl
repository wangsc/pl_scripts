#!/usr/bin/perl -w
use strict;

# parse blast-like file generated by "usearch -blastout";

my $sUsage = qq(
perl $0 
<usearch blast file>
<minimum cumulative aligned length proprotion, 0.6>
<minimum cumulative identity proprotion, 0.8>
);
die $sUsage unless @ARGV >= 3;
my ($file, $min_CALP, $min_CIP) = @ARGV;
parse_usearch_blast($file, $min_CALP, $min_CIP);

sub parse_usearch_blast
{
	my ($file, $min_calp, $min_cip) = @_;
	my $target_id;
	my $target_length;
	my $query_id;
	my $query_length;
	my $aln_length;
	my $idn_length;
	open (IN, $file) or die;
	while (<IN>)
	{
		next if /^\#/;
		if(/^Query\s+(\d+)nt\s+>(\S+)/)
		{
			if(defined $query_id)
			{
				my $calp = $aln_length/($target_length>$query_length?$query_length:$target_length);
				my $cip = $idn_length/$aln_length;
				print join("\t", ($query_id, $target_id, $query_length, $target_length, $aln_length, $idn_length)), "\n" if $calp>=$min_calp and $cip>=$min_cip;
				$query_length = $1;
				$query_id = $2;
				$aln_length = 0;
				$idn_length = 0;
 			}
 			else
 			{
				$query_length = $1;
				$query_id = $2; 
				$aln_length = 0;
				$idn_length = 0;				
 			}
 			next;
		}
		
		if(/^Target\s+(\d+)nt\s+>(\S+)/)
		{
			$target_length = $1;
			$target_id = $2;
			next;
		}
		
		if(/^(\d+)\s+cols\,\s+(\d+)\s+ids/) # 669 cols, 580 ids
		{
			$aln_length += $1;
			$idn_length += $2;
		}
	}
	close IN;
	my $calp = $aln_length/($target_length>$query_length?$query_length:$target_length);
	my $cip = $idn_length/$aln_length;
	print join("\t", ($query_id, $target_id, $query_length, $target_length, $aln_length, $idn_length)), "\n" if $calp>=$min_calp and $cip>=$min_cip;
}