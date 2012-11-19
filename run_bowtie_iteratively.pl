#!/usr/bin/perl -w
use strict;
use File::Basename;

my $sUsage = qq(
perl $0
<R1.fastq>
<R2.fastq>
<reference file (if multiple reference, separate with comma: ref1,ref2,ref3)>
<prefix for the final bam file>
);
die $sUsage unless @ARGV;
my ($r1_fastsq, $r2_fastq, $ref_str, $prefix) = @ARGV;

my $max_rd = 3;

my @references = split /,/, $ref_str;

# bowtie Parameters
# n, l ,e
my @params = (
[0, 48, 50],
[3, 28, 50],
[3, 15, 70]
);

# Iterative bowtie
my @sam_files;
foreach my $ref_index (0 .. $#references)
{
	
	my $round = 0;
	while ($round < $max_rd)
	{
		my $input_r1;
		my $input_r2;
		my $out_sam = $prefix . "_aligned_ref" . $ref_index . "_round". $round . ".sam";
		my $unalign = $prefix . "_unaligned_ref" . $ref_index . "_round". $round . ".fastq";
		#my $unalign_r2 = "unaligned_ref" . $ref_index . "_round". $round . ".fastq";		
		
		# determine input fastq to bowtie
		if($round == 0)
		{
			if($ref_index == 0)
			{
				($input_r1, $input_r2) = ($r1_fastsq, $r2_fastq); # just started				
			}
			else
			{
				my $tmp = $ref_index-1;
				$input_r1 = $prefix . "_unaligned_ref" . $tmp . "_round". ($max_rd-1) . "_1.fastq";
				$input_r2 = $prefix . "_unaligned_ref" . $tmp . "_round". ($max_rd-1) . "_2.fastq";
			}
		}
		else
		{
			my $tmp = $round-1;
			$input_r1 = $prefix . "_unaligned_ref" . $ref_index . "_round". $tmp . "_1.fastq";
			$input_r2 = $prefix . "_unaligned_ref" . $ref_index . "_round". $tmp . "_2.fastq";			
		}
		
		my ($n, $l, $e) = @{$params[$round]};
		
		my $log = $prefix . "_ref" . $ref_index . "_round". $round . ".log";
		# my $cmd="bowtie -a -m 1 $vn_mode $vn -l $l -e $e $ref -1 $r1_fastq -2 $r2_fastq -p 20 --sam tmp.sam 2>$log_file";
		my $ref = $references[$ref_index];
		my $bowtie_cmd = "bowtie -a -m 1 -n $n -l $l -e $e -X 500 $ref -1 $input_r1 -2 $input_r2 -p 20  --un $unalign --sam $out_sam 2>$log";
		print_time_stamp("Running comand:");
		print STDERR "CMD: ", $bowtie_cmd, "\n";
		die $! if system($bowtie_cmd);
		push @sam_files, $out_sam;
		my ($unique_aling, $supressed) = read_log($log);
		print STDERR "ref ", $ref_index, " round ", $round, "\n";
		print STDERR "\tParameters(n, l, e) : ", join(" ", (($n, $l, $e) )), "\n";
		print STDERR "\tunique aligned: $unique_aling\tsupressed: $supressed\n";
		
		print_time_stamp("Finished alignment for Ref ". $ref_index. " and round " . $round . "!\n");
		
		$round++;
		
	}	
}

# post_precess
print_time_stamp("*** Processing SAM files ... \n");
# remove unmapped reads
print_time_stamp("\t Remove unmapped reads ... \n");
my @bam_files;
foreach my $sam (@sam_files)
{
	my $out = basename($sam , ".sam") . "_MappedOnly.bam";
	my $sam_to_bam_cmd = "samtools view -F 4 -Sbh $sam > $out";
	print STDERR "\t", $sam_to_bam_cmd, "\n";
	die $! if system($sam_to_bam_cmd);
	push @bam_files, $out;
	print STDERR "\t", "Remove sam file: $sam", "\n";
	die if system("rm $sam");
}
# Combine bam files;
print_time_stamp("\t Merge BAM files ... \n");
my $picard_bin = "java -jar /home/DNA/Tools/picard/MergeSamFiles.jar MSD=true";
my @for_picard = map {"I=".$_;} @bam_files;
my $out = "O=" . $prefix . "_Merged_MappedReads.bam";
my $picard_cmd = join(" ", ($picard_bin, @for_picard, $out));
print STDERR "\t", $picard_cmd, "\n";
die $! if system($picard_cmd);

print_time_stamp("****** Finished... \n"); 


# subroutines

sub print_time_stamp
{
	my $comments = shift;
	my $time = localtime();
	print STDERR $time, "\t", $comments, "\n";
}

sub read_log
{
	my $log_file = shift;
	open (IN, $log_file) or die;
	my ($unique_aling, $supressed) = (0, 0);
	while(<IN>)
	{
		# reads with at least one reported alignment: 169615 (0.40%)
		if (/reads with at least one reported alignment:\s+(\d+)\s+\((\S+)\)/)
		{
			$unique_aling = join(" ", ($1,$2));
		}
		if(/reads with alignments suppressed due to -m:\s+(\d+)\s+\((\S+)\)/)
		{
			$supressed  = join(" ", ($1,$2));
		}				
	}
	close IN;
	return ($unique_aling, $supressed);	
}

