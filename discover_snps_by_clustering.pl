#!/usr/bin/perl -w
use strict;

print <<END unless @ARGV;
**********************************************************************************************
Though lack of reference genome, SNP discovery in the GBS data could be done by clustering.
This script will call cd-hit-454 to cluster the GBS reads at first,
then call cdhit-cluster-consensus to report the alignment of reads for each cluster.
After these two stes, SNPs in each cluster will be detected and reported.

Usage:
perl $0 
<Input file in fasta format>
<cutoff of similarity of reads (0.75-1), for cd-hit-454 -c >
<SNP output file>

Notes:
1. The input file should be quality filtered and adapter-dimer removed;
   If there have multiple accessions in the file, the format should be like:
   >accession-ID_read-id (ie, >acc-1_R6587)
2. Make sure cd-hit-454 and cdhit-cluster-consensus are in the PATH;

**********************************************************************************************
END

my ($input_fasta, $sim_cutoff, $snp_outfile) = @ARGV;

my $align_file = run_cluster($input_fasta, $sim_cutoff);
report_snps($input_fasta, $align_file, $snp_outfile);

# Subroutines

sub run_cluster
{
	my ($infile, $cutoff) = @_;
	my $cmd = "cd-hit-454 -i $infile -o $infile" . "_". $cutoff. " -c $cutoff ". "-d 0 -M 0 -T 0 -g 1";
	die "$cmd failed: $!\n" if system($cmd);
	my $in_for_align = $infile . "_".$cutoff . '.clstr';
	my $build_aln_cmd = "cdhit-cluster-consensus -clustfile=".$in_for_align. " -fastafile=" . $infile . "-aln_output=" . 
											 $in_for_align."_align -maxlen=1";
	die "$build_aln_cmd failed: $!\n" if system($build_aln_cmd);
	return $in_for_align."_align";
}

sub report_snps
{
	my ($infasta, $align_file, $output) = @_;
	
}





