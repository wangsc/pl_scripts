#!/usr/bin/perl -w
use strict;
use File::Basename;

my $sUsage = qq(
perl $0
<accession ID>
<bam files: 1.bam,2.bam,2.bam,...>
<R1 fastq files: R1_1.fastq,R1_2.fastq,R1_3.fastq,...>
<R2 fastq files: R2_1.fastq,R2_2.fastq,R2_3.fastq,...>
);

die $sUsage unless @ARGV == 4;

my $samtools_bin = "/homes/bioinfo/bioinfo_software/samtools/samtools ";

my ($acc_id, $bam_file_str, $R1_fastq_file_str, $R2_fastq_file_str) = @ARGV;
my @bam_files = split /,/, $bam_file_str;
my @R1_fastq_files = split /,/, $R1_fastq_file_str;
my @R2_fastq_files = split /,/, $R2_fastq_file_str;

get_unique_aligned_reads(\@bam_files, \@R1_fastq_files, \@R2_fastq_files);


# Subroutines

sub get_unique_aligned_reads
{
	my ($bam_ref, $r1_ref, $r2_ref) = @_;
	my %recorder;
	
	# get the ids of mapped reads
	my @bam_files = @$bam_ref;
	foreach my $bam (@bam_files)
	{
		print_time_stamp("\tProcessing file ". $bam);
		open(IN, "$samtools_bin view -h $bam |") or die "$!: can't read file $bam\n";
		while(<IN>)
		{
			if(/^\@/){print $_; next}
			my @data = split /\s+/, $_; 
			my $id = $data[0];
			$recorder{$id} = 0 unless exists $recorder{$id};
			$recorder{$id}++;	
		}
		close IN;
	}

	# remove non-unique mapping reads
	foreach my $bam (@bam_files)
	{
		print_time_stamp("Get uniquely mapped reads from file ".$bam);
		open(IN, "$samtools_bin view -h $bam |") or die "$!: can't read file $bam\n";
		my $out = basename($bam) . "_rmMultiHits.sam";
		print_time_stamp("\tOutput to sam file " . $out);
		open(OUT, ">$out") or die "can't open file $out\n";
		while(<IN>)
		{
			if(/^\@/){print OUT $_; next}
			my @data = split /\s+/, $_; 
			#my $id = join("*",@data[0,1]);
			my $id = $data[0];
			print OUT $_ unless $recorder{$id} > 2;	
		}
		close IN;
		close OUT;
	}	
	

	# Combine sam files;
	my $picard_dir = "/homes/wangsc/Tools/picard/";
	my $picard_bin = "java -jar $picard_dir/MergeSamFiles.jar MSD=true";
	my @sam_files = map{basename($_)."_rmMultiHits.sam"}@bam_files;
	my @for_picard = map {"I=".$_;} @sam_files;
	my $out = "O=" . $acc_id . "_All_chromosome_MappedReads_rmMultiHits.bam";
	my $picard_cmd = join(" ", ($picard_bin, @for_picard, $out));
	print_time_stamp($picard_cmd);
	die $! if system($picard_cmd);

	# output fastq file for unmapped reads
	my $r1_unmap_fastq =  $acc_id . "_All_chromosome_unmapped_R1.fastq";
	open (R1, ">$r1_unmap_fastq") or die;
	my $r2_unmap_fastq =  $acc_id . "_All_chromosome_unmapped_R2.fastq";
	open (R2, ">$r2_unmap_fastq") or die;
	
	print_time_stamp("Generating R1 fastq file: " . $r1_unmap_fastq);
	foreach my $f (@$r1_ref)
	{
		open (IN, $f) or die;
		my $id;
		my $cnt_line = 0;
		my @temp_arr = ();
		while(<IN>)
		{
			chomp;
			$cnt_line++; 
			push @temp_arr, $_;
			if($cnt_line == 1)
			{
				$id = $1 if /\@(\S+)/;
											
			}
			elsif($cnt_line == 4)
			{
				print R1 join("\n", @temp_arr), "\n" unless exists $recorder{$id};
				@temp_arr = ();
				$cnt_line = 0;
			}
		}
	}
	close R1;

	print_time_stamp("Generating R2 fastq file: " . $r2_unmap_fastq);
	foreach my $f (@$r2_ref)
	{
		open (IN, $f) or die;
		my $id;
		my $cnt_line = 0;
		my @temp_arr = ();
		while(<IN>)
		{
			chomp;
			$cnt_line++; 
			push @temp_arr, $_;
			if($cnt_line == 1)
			{
				$id = $1 if /\@(\S+)/;
											
			}
			elsif($cnt_line == 4)
			{
				print R2 join("\n", @temp_arr), "\n" unless exists $recorder{$id};
				@temp_arr = ();
				$cnt_line = 0;
			}
		}	
	}
	close R2;			
	
}

sub print_time_stamp
{
        my $comments = shift;
        my $time = localtime();
        print STDERR $time, "\t", $comments, "\n";
}