#!/usr/bin/perl -w
use strict;
use File::Basename;

my $samtools_bin = "/homes/bioinfo/bioinfo_software/samtools/samtools";
my $picard_dir = "/homes/wangsc/Tools/picard/";
my $picard_bin = "java -jar $picard_dir/MergeSamFiles.jar MSD=true";
my $NCPU = 16;

my @sam_files = @ARGV;

my @bam_files;
my @sam2bam_cmds;
foreach my $sam (@sam_files)
{
        my $out = basename($sam , ".sam") . "_MappedOnly.bam";
        my $sam_to_bam_cmd = "$samtools_bin view -F 4 -Sbh $sam > $out";
        push @sam2bam_cmds, $sam_to_bam_cmd;
        #print STDERR "\t", $sam_to_bam_cmd, "\n";
        #die $! if system($sam_to_bam_cmd);
        push @bam_files, $out;
        #print STDERR "\t", "Remove sam file: $sam", "\n";
        #die if system("rm $sam");
}
my $log_dir = "cmds_log.$$";
mkdir($log_dir) unless -d $log_dir;

process_cmds(\@sam2bam_cmds, $NCPU);

my %prefix_files;
map{my $basename = basename($_); my $prefix = $1 if $basename=~/(\S+)_aligned/; push @{$prefix_files{$prefix}}, $_}@bam_files;

# Combine bam files;
my @combine_cmds;
foreach my $prefix (keys %prefix_files)
{
	print_time_stamp("\t Merge BAM files for $prefix ... \n");
	my @bam_files = @{$prefix_files{$prefix}};
	my @for_picard = map {"I=".$_;} @bam_files;
	my $out = "O=" . $prefix . "_Merged_MappedReads.bam";
	my $picard_cmd = join(" ", ($picard_bin, @for_picard, $out));
	push @combine_cmds, $picard_cmd;
	#print STDERR "\t", $picard_cmd, "\n";
	#die $! if system($picard_cmd);	
}

process_cmds(\@combine_cmds, $NCPU);

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





