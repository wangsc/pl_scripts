#!/usr/bin/perl -w
use strict;
use Cwd;
use File::Basename;

my $sUsage = q(
# Suppose we have *_Merged_MappedReads.bam in the current folder
perl $0 

);

my $samtools_bin = "/homes/bioinfo/bioinfo_software/samtools/samtools ";

my $cwd = getcwd();
my @files = <*Merged_MappedReads.bam>;

my @sh_files;
foreach my $f (@files)
{
	my $acc = $1 if $f=~/(\S+)_Merged_MappedReads.bam/;
	my $basename = basename($f, ".bam");
	my $sort_output = $basename . "_sorted";
	my $rmdup_output = $sort_output . "_rmdup.bam";
	my $sort_cmd = $samtools_bin . "sort $f $sort_output";
	my $rmdup_cmd = $samtools_bin . "rmdup $sort_output" . ".bam " . "$rmdup_output";
	
	my $out_sh = $acc . ".sh";
	$out_sh = "q" . $out_sh if $out_sh =~/^\d/;
	push @sh_files, $out_sh;
	
	open (OUT, ">$out_sh") or die;
	print OUT "#!/bin/bash\n", "cd $cwd\n", $sort_cmd, "\n", $rmdup_cmd, "\n";
	close OUT;
}

map{print "qsub -cwd -l mem=2G,h_rt=10:0:0 -pe mpi-1 1 $_\n"}@sh_files;