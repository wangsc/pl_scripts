#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my $sUsage = <<USAGE;

Analysis steps:
# Step 1: quality filter and trim
# Step 2: barcode seperation
# Step 3: sequence cluster with cd-hit
# Step 4: read mapping
# Step 5: call variations in tags
# Step 6: call genotypes

GBS analysis pipeline options:

perl $0
*General options:
  --fastq_dir         directory where you put fastq files; all fastq files in the directory will be processed
  --gbs_pipeline_dir  directory where you installed gbs pipeline
  --output_dir        directory where all the output files will be put
  --barcode           plain text file that contained the barcodes for each accession: [accession_name  barcode]; 
                      if not provided, barcode seperation will not be processed.
  --call_snp          1 or 0, default is 1;  indicator for calling SNPs in tags
  --call_PA           1 or 0, default is 0;  indicator for calculating presence or absence of tags in each accession
  --start_step        the step you want to start, default is 1;

*Options for read quality filter and trimming  
  --min_qual          minimum quality score for bases, will be used for read quality filter and trim, default is 15
  --min_length        minimum length of reads after trimming to be retained, default 30

*Options for sequence cluster with cd-hit  
  --cd_hit_bin        must be provided unless it has beed included in PATH 
  --similarity        -c for cd-hit, default is 0.9 (90%);

*Options for mapping
  --mismatches        maximum number of mismatches allowed in mapped reads,default is 3
  --cpu               number of cpus used for mapping, default is 10;

*Options for snp calling
  --samtools_dir      directory that installed samtools
  --bcftools_dir      directory that installed bcftools
  
USAGE

my ($fastq_dir, $pipeline_dir, $output_dir, $barcode_file, $call_snp, $call_PA, $start_step, $min_qual, $min_length, $cd_hit_bin, $similarity);
my ($mismatches, $cpu) = (3, 10);
$min_qual = 15;
$min_length = 30;
$call_snp = 1;
$call_PA = 0;
$start_step = 1;
$similarity = 0.9;
my ($samtools_dir, $bcftools_dir);
GetOptions(
"fastq_dir=s"        => \$fastq_dir,
"gbs_pipeline_dir=s" => \$pipeline_dir,
"output_dir=s"       => \$output_dir,
"barcode=s"          => \$barcode_file,
"call_snp=i"         => \$call_snp,
"call_PA=i"          => \$call_PA,
"start_step=i"       => \$start_step,
"min_qual=i"         => \$min_qual,
"min_length=i"       => \$min_length,
"cd_hit_bin=s"       => \$cd_hit_bin,
"similarity=i"       => \$similarity,
"mismatches=i"       => \$mismatches,
"cpu=i"              => \$cpu,
"samtools_dir=s"     => \$samtools_dir,
"bcftools_fir=s"     => \$bcftools_dir
);
die "\n---Must define fastq_dir, gbs_pipeline_dir and output_dir\n", $sUsage 
		unless defined $fastq_dir and defined $pipeline_dir and defined $output_dir;
die "output_dir should NOT be the same as fastq_dir\n" if abs_path($fastq_dir) eq abs_path($output_dir);
if (-d $output_dir)
{
	warn "$output_dir already exist!\n";
}
else
{
	mkdir($output_dir);
}

my ($script, $cmd);

# Step 1: quality filter and trim
my $current_step = 1;
print_time_stamp('# Step 1: quality filter and trim');
$script = $pipeline_dir . "/1_quality_filter_and_trim.pl";
my $qc_out_dir = $output_dir ."/reads_after_quality_control";
mkdir($qc_out_dir) unless -d $qc_out_dir;
$cmd = "perl $script $fastq_dir $qc_out_dir $min_qual $min_length";
print STDERR "\t$cmd\n";
if($current_step>=$start_step){die "$cmd failed !!" if system($cmd);}

# Step 2: barcode seperation
$current_step++;
my $barcode_sep_dir;
if (defined $barcode_file)
{
	print_time_stamp('# Step 2: barcode seperation');
	$script = $pipeline_dir . '/' . "2_barcode_seperator.pl";
	$barcode_sep_dir = $output_dir . "/barcode_seperated";
	mkdir($barcode_sep_dir) unless -d $barcode_sep_dir;
	$cmd = "perl $script $qc_out_dir $barcode_file $barcode_sep_dir";
	print STDERR "\t$cmd\n";
	if($current_step>=$start_step){die "$cmd failed !!" if system($cmd);}
}
else
{
	print_time_stamp('# Step 2: barcode seperation skipped because of no barcode file provided')
}

# Step 3: sequence cluster with cd-hit
$current_step++;
print_time_stamp('Step 3: sequence cluster with cd-hit');
my $cd_hit = defined $cd_hit_bin?$cd_hit_bin:"";
$cd_hit = $cd_hit eq ""?'cd-hit-est':$cd_hit . '/' . 'cd-hit-est';
$script = $pipeline_dir . "/3_seq_cluster.pl";
my $cd_hit_out_dir = $output_dir . "/cd_hit_output";
mkdir($cd_hit_out_dir) unless -d $cd_hit_out_dir;
my $cd_hit_out_file = $cd_hit_out_dir . '/seq_cluster_' . $similarity . ".fasta";
$cmd = "perl $script $cd_hit $barcode_sep_dir $cd_hit_out_dir $cd_hit_out_file $similarity";
print STDERR "\t$cmd\n";
if($current_step>=$start_step){die "$cmd failed !!" if system($cmd);}

# Step 4: read mapping
$current_step++;
print_time_stamp('Step 4: read mapping with bowtie');
my $mapping_out_dir = $output_dir . "/mapping_output";
mkdir($mapping_out_dir) unless -d $mapping_out_dir;
my $bowtie_index_dir = $output_dir . "/bowtie_index";

unless (-d $bowtie_index_dir)
{
	print STDERR "\t Building bowtie index for $cd_hit_out_file\n";
	mkdir($bowtie_index_dir);
	my $base = basename($cd_hit_out_file, ".fasta");
	$base = $bowtie_index_dir . "/" . $base;
	my $build_index = "bowtie-build $cd_hit_out_file $base";
	die $! if system($build_index);
}
my @fastq_files = <$barcode_sep_dir/*.fastq>;
foreach my $file (@fastq_files)
{
	my $basename = basename($file, ".fastq");
	my $sam_out = $mapping_out_dir . "/". $basename . ".sam";
	$cmd = "bowtie --best -v $mismatches -p $cpu " . $bowtie_index_dir . "/" . basename($cd_hit_out_file, ".fasta") . " -q $file --sam $sam_out";
	last unless $current_step>=$start_step;
	print STDERR "\t$cmd\n";
	die "$cmd failed !!" if system($cmd);
}

# Step 5: call variations in tags
$current_step++;
print_time_stamp('Step 5: call variations in tags');
my @sam_files = <$mapping_out_dir/*.sam>;
my $samtools_bin = defined $samtools_dir?$samtools_dir."/samtools":"samtools";
my $bcftools_bin = defined $bcftools_dir?$bcftools_dir."/bcftools":"bcftools";
my $vcfutils = defined $bcftools_dir?$bcftools_dir."/vcfutils.pl":"vcfutils.pl";
foreach my $sam (@sam_files)
{
	my $u = $sam;
	$u =~ s/\.sam$//;
	$cmd="$samtools_bin view -bS $sam -o $u".".bam; ".  # sam to bam
			 "$samtools_bin sort $u".".bam $u"."_sorted; ".  # sort bam
			 "$samtools_bin mpileup -uf $cd_hit_out_file $u"."_sorted.bam |". 
			 "$bcftools_bin view -bvcgIN - >$u"."_raw_var.bcf;". 
			 "$bcftools_bin view $u"."_raw_var.bcf |". 
			 "$vcfutils varFilter -d 4  >$u"."_var.filtered.vcf";
	
	last unless $current_step>=$start_step;
	print STDERR "\t$cmd\n\n";
	die "$cmd failed !!" if system($cmd);
}
# combine all var.filter.vcf
my $all_vcf = $output_dir . "/all_snps.vcf"; 
my @vcf_files = <$mapping_out_dir/*var.filtered.vcf>;
combine_vcf_files($all_vcf, \@vcf_files);

# Step 6: call genotypes
$current_step++;
print_time_stamp('Step 6: call genotypes');
print_time_stamp('  Step 6-1 : count allele frequency');
$script  = $pipeline_dir . '/'. '6-1_count_allele_frequency_in_bam.pl';
my @bam_files = <$mapping_out_dir/*sorted.bam>;
foreach my $bam (@bam_files)
{
	my $tmp_out = $bam;
	$tmp_out =~ s/sorted.bam/allele_freq\.out/;
	$cmd = "perl $script $all_vcf $cd_hit_out_file $bam $tmp_out";
	last unless $current_step>=$start_step;
	print STDERR "\t$cmd\n";
	die "$cmd failed !!" if system($cmd);
}
print_time_stamp('  Step 6-2: call genotypes from sam files');
my $geno_out = $output_dir . "/" ."snp_genotype.csv";
$script = $pipeline_dir . '/6-2_call_genotype.pl';
$cmd = "perl $script $mapping_out_dir $geno_out";
print STDERR "\t$cmd\n";
if($current_step>=$start_step){die "$cmd failed !!" if system($cmd);}


# Subroutines
sub print_time_stamp
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}

sub combine_vcf_files
{
	my ($outfile, $in_files) = @_;
	my %recorder;
	foreach my $file (@$in_files)
	{
		open (IN, $file) or die $!;
		while(<IN>)
		{
			next if /^\#/; 
			my @t=split/\t/,$_;
			next if length $t[4] > 1;
			$recorder{join("\t", @t[0,1])} = $_;
		}
		close IN;
	}
	open (OUT, ">$outfile") or die $!;
	foreach (values %recorder)
	{
		print OUT $_;
	}
	close OUT;
	;
}




