#!/usr/bin/perl -w
use strict;
use File::Basename;

my $samtools_bin = "/homes/bioinfo/bioinfo_software/samtools/samtools";
my $picard_dir = "/homes/wangsc/Tools/picard/";
my $picard_bin = "java -jar $picard_dir/MergeSamFiles.jar MSD=true";
my $NCPU = 8;
my $log_dir = "cmd_logs.$$";

my @bam_files = @ARGV;

# sort bam files
my @sort_cmds;
my @sorted_bam_files;
foreach my $bam (@bam_files)
{
	my $base = basename($bam, ".bam");
	my $out_prefix = $base . "_sorted";
	push @sorted_bam_files, $out_prefix . ".bam";
	my $cmd = "$samtools_bin sort $bam $out_prefix";
	print $cmd, "\n";
	push @sort_cmds, $cmd;
}
process_cmds(\@sort_cmds, $NCPU);


# rmdup for bam files
my @rmdup_cmds;
foreach my $bam (@sorted_bam_files)
{
	my $base = basename($bam, ".bam");
	my $out = $base . "_rmdup.bam";
	my $cmd = "$samtools_bin rmdup $bam $out";
	print $cmd, "\n";
	push @rmdup_cmds, $cmd;
}
process_cmds(\@rmdup_cmds, $NCPU);

# Subroutines
sub process_cmds
{
	my ($cmds_ref, $ncpu) = @_;
	my %job_tracker;
	my @failed_jobs;
	
	my $num_running = 0;
	my $cmd_counter = 0;
	foreach my $cmd (@$cmds_ref)
	{
		$cmd_counter++;
		$num_running++;
		
		my $child = fork();
		if($child)
		{
			# parent
			$job_tracker{$cmd_counter} = $cmd;
		}
		else
		{
			#child
			my $ret = run_cmd($cmd_counter, $cmd);
			exit ($ret);
		}
		
		if($num_running >= $ncpu)
		{
			wait();
			my $num_finished = collect_jobs(\%job_tracker, \@failed_jobs);
			$num_running -= $num_finished;
		}
	}
	
	while ((wailt() != -1)) {} # wait till all child processes finished
	
	collect_jobs(\%job_tracker, \@failed_jobs);
	
	`rm -rf $log_dir`;
   my $num_failed_jobs = scalar @failed_jobs;
   if ($num_failed_jobs == 0) {
                print "\n\nAll $cmd_counter jobs completed successfully! :) \n\n";
                exit(0);
   }
   else {
                # write all failed commands to a file.
                my $failed_cmd_file = "failed_cmds.$$.txt";
                open (my $ofh, ">$failed_cmd_file") or die "Error, cannot write to $failed_cmd_file";
                print $ofh join("\n", @failed_jobs) . "\n";
                close $ofh;

                print "\n\nSorry, $num_failed_jobs of $cmd_counter jobs failed.\n\n"
                        . "Failed commands written to file: $failed_cmd_file\n\n";
                exit(1);
   }

}

sub run_cmd {
        my ($index, $cmd) = @_;

        print "\nRUNNING: $cmd\n";

        my $ret = system($cmd);

        if ($ret) {
                print STDERR "Error, command: $cmd died with ret $ret";
        }

        open (my $log_fh, ">$log_dir/$index.ret") or die "Error, cannot write to log file for $index.ret";
        print $log_fh $ret;
        close $log_fh;
        return($ret);
}


sub collect_jobs {
	my ($job_tracker_href, $failed_jobs_aref) = @_;
	
	my @job_indices = keys %$job_tracker_href;
	
	my $num_finished = 0;
	
	foreach my $index (@job_indices) {
	
		my $log_file = "$log_dir/$index.ret";
	
	  if (-s $log_file) {
	  	my $ret_val = `cat $log_file`;
	    chomp $ret_val;
	    my $job = $job_tracker_href->{$index};
	    if ($ret_val == 0) {
	     # hurray, job succeded.
	    print "SUCCESS[$index]: $job\n";
	
	    }
	   else {  
	   # job failed.
	   	print "FAILED[$index]: $job\n";
	    push (@$failed_jobs_aref, $job_tracker_href->{$index});
	   }
	   unlink $log_file;
	   $num_finished++;
	  }
	}
	
	return($num_finished);
}





