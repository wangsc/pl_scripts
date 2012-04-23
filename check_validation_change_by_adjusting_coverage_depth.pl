#!/usr/bin/perl -w
use strict;

my @depth = 1..10;
my $run_snp_filter = 1;
foreach my $dp (@depth)
{
	# run snp_filter.pl
	my $cmd = "perl ~/snp_filter.pl ";
	my $out_file = "dp_".$dp.".out";
	my $params = "../all_var.filtered.vcf $dp acc_snp_freq ../*allele_frequency.out_formatted >". $out_file;	
	#if($run_snp_filter){die "$cmd failed\n" if system($cmd . $params)};
	
	# extract snps
	my $snp_file = "min".$dp."_2_2.snps";
	my $cmd_ex_snps = q"perl -ne '@t=split/\t/,$_; next unless @t>=5; $c=0; foreach $r(@t[2,3,4]){if($r=~/^\S{1}_(\d+)/){$c++ if $1>=2}} print $_ if $c==2' " . $out_file . ">". $snp_file;
	die "$cmd_ex_snps failed \n" if system($cmd_ex_snps);
	
	# extract flanking seq for snps
	my $flank_seq_file = "min".$dp."_2_2_50bp_flanking.fasta";
	my $oligo_file = "min".$dp."_2_2_50bp_oligo.fasta";
	my $cmd_flank = "perl ~/extract_flanking_seq_for_snps.pl ";
	$params = $snp_file . " ../all9lines_and_flcDNA.fasta2 50  $oligo_file >" . $flank_seq_file;
	die "$cmd_flank failed \n" if system($cmd_flank . $params);
	
	# run blat
	my $blat_out = '9k_blat_min'.$dp."_2_2.out";
	my $blat_cmd = "blat $flank_seq_file  PrivKSU_WheatCons_9k_11497518_A.fasta -fastMap -out=blast8 ". $blat_out;
	die "$blat_cmd failed \n" if system($blat_cmd);
	
	# get the info for SNPs
	my $shiaoman_file = "/home/swang/EduardWheat9Kproject_data.csv_convered2.csv";
	my $monomorphic_file = "/home/DNA/Data_analysis/SNP_discovery_in_Illuminar_data/compare_to_mira_snps/monomorphic_snps_in_9k_assay.list";
	my ($worked, $not_worked, $mono) = count($snp_file, $blat_out, $shiaoman_file, $monomorphic_file, $dp);
	print '*'x20,"\n";	
	print $dp,"\t", join("\t",($worked, $not_worked, $mono, sprintf("%.2f", $worked/($worked+$not_worked+$mono)), sprintf("%.2f", $mono/($worked+$not_worked+$mono)))),"\n";
	print '*'x20,"\n";
	# 
	
}

sub read_indication_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die "can't open file $file \n";
	while(<IN>)
	{
		chomp;
		my @t = split /\t/,$_; 
		$return{join("\t", @t[0,1])} = $t[-1];			
	}
	close IN;
	return %return;
}

sub count
{
	my ($snpfile, $blatfile, $shiaoman_file, $monomorphic_file, $dp) = @_;
	my %mono;
	open(IN, $monomorphic_file) or die;
	while(<IN>)
	{
		chomp; 
		my @t = split /\t/,$_;
		$mono{$t[0]} = 1;
	}
	close IN;
	
	my %shiaoman;
	open (IN, $shiaoman_file) or die "can't open $shiaoman_file\n";
	while(<IN>)
	{
		chomp; 
		s/\"//g;
		next unless /^\d/;
		my @t=split /\t/,$_;
		$shiaoman{$t[1]}=1;
	}
	close IN;
	
	my %blat;
	open (IN, $blatfile) or die "can't open $blatfile\n";
	while(<IN>)
	{
		chomp; 
		my @t=split /\t/,$_;
		my @m = split /-/,$t[0];
		$t[1] =~ s/\:(\d+)$//;
		push @{$blat{$m[0]}}, $t[1];
	}
	close IN;
	
	my ($worked, $not_worked, $monosnp) = (0, 0, 0);
	my %snp_info;
	foreach my $id (keys %blat)
	{
		if (exists $shiaoman{$id} and !(exists $mono{$id}))
		{
			$worked++;
			foreach (@{$blat{$id}}){$snp_info{$_}=1}
		}
		elsif (exists $shiaoman{$id} and (exists $mono{$id}))
		{
			$monosnp++;
			foreach (@{$blat{$id}}){$snp_info{$_}=0}
		}
		else
		{
			$not_worked++;
			foreach (@{$blat{$id}}){$snp_info{$_}=0}
		}	
	}
	
	my %snps;
	open (IN, $snpfile) or die "can't open $snpfile\n";
	my $temp_out = "min".$dp.'_2_2_incommon_info';
	open (OUT, ">$temp_out") or die "can't open file $temp_out\n";
	while(<IN>)
	{
		chomp; 
		my @t=split /\t/, $_;
		my $id = join(':', @t[0,1]);
		next unless exists $snp_info{$id};
		print OUT $_,"\t", $snp_info{$id},"\n";
		
	}
	close IN;	
	close OUT;
	return ($worked, $not_worked, $monosnp);
}






