#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<common_snps_worked_info>
<all_accessions_snp_frequency.out>
<all_var.filtered.vcf>
<output file prefix>
);
die $sUsage unless @ARGV >= 4;
my($snp_file, $frequency_file, $var_file, $out_prefix)  = @ARGV;

my %snp_info = read_snp_info_file($snp_file);
my %alleles = read_var_file($var_file);
my ($allele_count, $acc_count) = read_frequency_file($frequency_file, \%alleles);

my $out_ids = $out_prefix . "_ids";
my $out_features = $out_prefix . "_features";
open (I, ">$out_ids") or die "can't open file $out_ids";
open (F, ">$out_features") or die "can't open file $out_features";

foreach my $id (keys %$allele_count)
{
	print I $id,"\n";
	my @features = (@{$allele_count->{$id}}, @{$acc_count->{$id}});
	print F 0;
	if (exists $snp_info{$id}){print $snp_info{$id}?$snp_info{$id}:-1 }
	foreach (1..5)
	{
		print F "\t", $_,':',$features[$_-1];
		print "\t", $_,':',$features[$_-1] if exists $snp_info{$id};
	}
	print "\n" if exists $snp_info{$id};
	print F "\n";
	
}
close I;
close F;



# subroutines

sub read_snp_info_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die "can't open file $file\n";
	while(<IN>)
	{
		chomp;
		next if /^\s+$/;
    # Kukri_mira1_c10984:267:757_R    0
    my @t = split /\t/,$_;
    my @m = split /:/, $t[0];
    $return{join("\t", @m[0,1])} = $t[1];
	}
	close IN;
	return %return;
}

sub read_var_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die $!;
	my $debug = 1;
	while(<IN>)
	{
		# contig00066     679     .       A       C       55      .       DP=10;AF1=1;AC1=2;DP4=0,0,0,10;MQ=20;FQ=-57     GT:PL:GQ        1/1:88,30,0:57
		# contig00090     493     .       C       A       43      .       DP=9;AF1=1;AC1=2;DP4=0,0,9,0;MQ=20;FQ=-54       GT:PL:GQ        1/1:76,27,0:51
		# contig00101     170     .       T       C,A     37.3    .       DP=13;AF1=1;AC1=2;DP4=0,1,0,12;MQ=20;FQ=-41;PV4=1,0.12,1,1      GT:PL:GQ        1/1:71,15,1,56,0,61:25
    # contig00101     551     .       T       A       15.2    .       DP=6;AF1=0.5032;AC1=1;DP4=2,0,4,0;MQ=20;FQ=-8.64;PV4=1,6.4e-06,1,1      GT:PL:GQ        0/1:45,0,19:22		
		next if /^\#/;
		my @data = split /\t/,$_;
		next if length $data[4] > 1;
		$return{join("\t", @data[0,1])} = [@data[3,4]] if $data[5]>=20;
		print STDERR 'var ID: ', join("\t", @data[0,1]), "\n" if $debug; $debug=0;
	}
	close IN;
	return %return;
}

sub read_frequency_file
{
	my $file = shift;
	my $allele_ref = shift;
	my $min_depth = 4;
	my %allele_count;
	my %accession_count;
	my $debug = 1;
	open (IN, $file) or die "can't open file $file\n";
	while(<IN>)
	{
		#AC_Barrie       RAC875_mira1_rep_c73204 67      A_16    C_4
		chomp; 
		my @t = split /\t/,$_;
		my $id = join("\t", @t[1,2]);
		print STDERR 'freq ID: ', $id,"\n" if $debug; 
		next unless exists $allele_ref->{$id};
		$allele_count{$id} = [0, 0] unless exists $allele_count{$id};
		$accession_count{$id} = [0, 0, 0] unless exists $accession_count{$id};
		my %count;
		foreach (3..$#t)
		{
			if ($t[$_]=~/(\S+)_(\d+)/)
			{
				print STDERR $1,"\t",$2,"\n" if $debug; $debug=0;
				$count{$1} = $2
			}
		}
		
		my @alleles = @{$allele_ref->{$id}}; 
		map{ $count{$_} = 0 unless exists $count{$_} } @alleles;
		$allele_count{$id}[0] += $count{$alleles[0]};
		$allele_count{$id}[1] += $count{$alleles[1]};
		if($count{$alleles[0]} == 0)
		{
			if($count{$alleles[1]} > 0){ $accession_count{$id}[1]++ if $count{$alleles[1]}>=$min_depth }
		}
		else
		{			
			if($count{$alleles[1]} < $min_depth ){ $accession_count{$id}[0]++ if $count{$alleles[0]}>=$min_depth }
			else{$accession_count{$id}[2]++ if ($count{$alleles[1]}>=$min_depth and $count{$alleles[0]}>=$min_depth) }
		}

	}
	close IN;
	return(\%allele_count, \%accession_count)
}



