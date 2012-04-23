#!/usr/bin/perl -w
use strict;
use File::Basename;
use Parallel::ForkManager;

my ($fastq_dir, $output_dir, $min_qual, $min_length) = @ARGV;
my @fastq_files  = <$fastq_dir/*.fastq>;

# check if fastx tool kits installed
my $t = `which fastq_quality_trimmer`;
unless ($t)
{
	die <<END;
	Error: cannot call fastq_qulaity_trimmer!
	fastx tool kits not installed or not in the PATH!
END
}
my @cmds;

foreach my $file (@fastq_files)
{
	my $basename = basename($file);
	my $out = $output_dir . "/" . $basename . "_trimmed_filtered";
	my $cmd = "fastq_quality_trimmer -Q33 -i $file -t $min_qual -l $min_length|".
            "fastq_quality_filter -Q33 -q $min_qual -p 80 -o $out";
  push @cmds, $cmd;
}

my $MAX_PROCESSES = 5;
my $pm = new Parallel::ForkManager($MAX_PROCESSES);

foreach my $cmd (@cmds)
{
	my $pid = $pm->start and next;
	if(system($cmd))
	{
		print STDERR "$cmd failed!!\n";
	}
	$pm->finish;
}