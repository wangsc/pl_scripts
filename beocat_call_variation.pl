#!/usr/bin/perl -w
use strict;
use Cwd;
use File::Basename;

my $cwd = getcwd();

my @bam_files = <*bam>;
my $reference_directory = "/homes/wangsc/Seq_Capture/wheat_reference/chr1_2_sep";

my @ref_files = <$reference_directory/*fasta>;

my $samtools_bin = "/homes/bioinfo/bioinfo_software/samtools/samtools ";
my $bcftools_bin = "/homes/bioinfo/bioinfo_software/samtools/bcftools/bcftools ";
my $vcfutils_bin = "/homes/bioinfo/bioinfo_software/samtools/bcftools/vcfutils.pl ";

my @sh_files;
foreach my $ref (@ref_files)
{
	my $chr = basename($ref, ".fasta");
	my $out_sh = $chr . ".sh"; 
	$out_sh = "q".$out_sh if $out_sh =~ /^\d/;
	open (OUT, ">$out_sh") or die;
	push @sh_files, $out_sh;
	# $cmd=" samtools mpileup -uf ../all9lines_and_flcDNA.fasta2 $f" . "|bcftools view -bvcgIN - >$u"."_raw_var.bcf;". " bcftools view 
  # $u"."_raw_var.bcf". "|vcfutils.pl varFilter -d 4  >$u"."_var.filtered.vcf";
	my $cmd_1 =  $samtools_bin . "mpileup  -uDf $ref " . join(" ", @bam_files) . "|$bcftools_bin view -bvcgIN - >$chr" . "_raw_var.bcf";
	my $cmd_2 = $bcftools_bin . " view $chr"."_raw_var.bcf"."|".$vcfutils_bin." varFilter -d 4 >$chr"."_var.filtered.vcf";
	
	#print $cmd_1, "\n$cmd_2\n";
	print OUT "#!/bin/bash\n\n", "cd $cwd\n\n", $cmd_1, "\n\n", $cmd_2, "\n\n";
	close OUT;
}

map{print "qsub -cwd -l mem=10G,h_rt=10:0:0 -pe mpi-1 1 $_\n"}@sh_files;