#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
-v or -n
Reference_base
R1 fastq file
R2 fastq file
);
die $sUsage unless @ARGV;
my ($vn_mode, $ref, $r1_fastq, $r2_fastq) = @ARGV;
my @vns=0..3;
my @ls=map{28+$_}0..20; 
my @es=map{50+$_*10}0..15;

foreach my $vn (@vns)
{
	foreach my $l (@ls)
	{
		foreach my $e (@es)
		{
			my $param = "P" . join("-", ($vn_mode.$vn, "l".$l, "e".$e));
			my $log_dir = "logs";
			mkdir($log_dir) unless -d $log_dir;
			my $log_file = $log_dir . "/" . $param . ".log";
			my $cmd="bowtie -a -m 1 $vn_mode $vn -l $l -e $e $ref -1 $r1_fastq -2 $r2_fastq -p 20 --sam tmp.sam 2>$log_file";
			print STDERR $cmd, "\n";
			die if system($cmd);
			open (IN, $log_file) or die;
			my ($unique_aling, $supressed) = (0, 0);
			while(<IN>)
			{
				# reads with at least one reported alignment: 169615 (0.40%)
				if (/reads with at least one reported alignment:.*\((\S+)\)/)
				{
					$unique_aling = $1;
				}
				if(/reads with alignments suppressed due to -m:.*\((\S+)\)/)
				{
					$supressed  = $1;
				}				
			}
			close IN;
			print join("\t", ($param, $unique_aling, $supressed)), "\n";
		}
	}
}