#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
< files of read_1: AC_Barrie_1.fastq,Alsen_1.fastq,... >
< files of read_2: AC_Barrie_2.fastq,Alsen_2.fastq,... >
);

die $sUsage unless @ARGV;
my($files_left, $files_right) = @ARGV;

my $tophat_bin_dir = "/usr/local/bin/";
my $cufflinks_bin_dir = "/usr/local/bin/";
my $tophat_bin = $tophat_bin_dir . "tophat";
my $cufflinks_bin = $cufflinks_bin_dir . "cufflinks";

my $tophat_out_dir = "./tophat_out/";
my $cufflinks_out_dir = "./cufflinks_out/";

my $bowtie_index_dir = "/home/DNA/ctgs/bowtie_index_chr";
#  wheat contigs were spliced into 10 files and build index seperately
my $num_index = 7;
foreach my $ind (1 .. $num_index)
{
	print_time_comment("******* Processing Index $ind ...");
	my $tophat_out = $tophat_out_dir . $ind;
	my $cufflinks_out = $cufflinks_out_dir . $ind;
	my $index_name = "Chr_" . $ind;
	$index_name = $bowtie_index_dir . "/" . $index_name;
	my $tophat_parameters = "-p 20 --no-convert-bam -o $tophat_out $index_name $files_left $files_right";
	my $tophat_cmd = join(" ", ($tophat_bin, $tophat_parameters));
	die if system($tophat_cmd);
	my $sam_file = "accepted_hits.sam";
	my $sorted_sam = "accepted_hits_sorted.sam";
	my $sort_cmd = "sort -k3,3 -k4,4n $tophat_out". "/". $sam_file . ">" .$tophat_out ."/".$sorted_sam;
	die if system($sort_cmd);
	my $cufflinks_cmd = $cufflinks_bin . " -p 20 -o $cufflinks_out ". $tophat_out ."/".$sorted_sam;
	die if system($cufflinks_cmd);
}



sub print_time_comment
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}