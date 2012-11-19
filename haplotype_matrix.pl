#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0 
<PHS input file>
<haplotype_matrx output>
<selection index file>
);
die $sUsage unless @ARGV >= 3;
my ($in_file, $hap_out, $idx_out) = @ARGV;
#my $in_file = shift or die $sUsage;
my %genotype = read_phs_input($in_file);
my @snp_files = <*.snps>;
die "No SNP files in current folder!\n" unless @snp_files;
my @snp_array;
my %haplotype;
my %select;
foreach my $snpfile (@snp_files)
{
	print STDERR "Processing file $snpfile ...\n";
	my %snp_dist = read_snp_file($snpfile);
	my $chr = $1 if $snpfile =~ /(\S+)\.snps/;
	my $phs_file = $chr . "_v3_phs.txt_addSNPs";
	open (IN, $phs_file) or die "can't open file $phs_file\n";
	while(<IN>)
	{
		chomp; 
		my @t=split /\s+/,$_;
		next unless $t[-1] =~ /wsnp/;
		my ($phs_a, $phs_b) = @t[4, 5];
		my $selected_allele = $phs_a > $phs_b?"A":"B";
		my @snps = ($t[2], (split /,/,$t[-1]));
		push @snp_array, $t[2];
		@snps = unique(@snps);
		map{print STDERR $_, "*\n" unless exists $snp_dist{$_}}@snps;
		@snps = sort{$snp_dist{$a}<=>$snp_dist{$b}} @snps;
		my %hap_tmp;
		#my %select_idx;
		foreach my $acc (keys %genotype)
		{
			
			my $geno_ref = $genotype{$acc};
			my @snps_data = grep {exists $geno_ref->{$_}} @snps;
			my $haplotype = join("", @{$geno_ref}{@snps_data});
			$hap_tmp{$acc} = $haplotype;
			#$select_idx{$haplotype} = $geno_ref->{$t[2]} eq $selected_allele?1:0;
			$select{$acc}{$t[2]} = $geno_ref->{$t[2]} eq $selected_allele?1:0;
		}
		my %hap_name = generate_hap_name(values %hap_tmp);
		
		foreach my $acc (keys %genotype)
		{
			$haplotype{$acc}{$t[2]} = $hap_name{$hap_tmp{$acc}};			
		}		
	}
	close IN;	
}

# output
open (H, ">$hap_out") or die "can't open file $hap_out\n";
open (D, ">$idx_out") or die "can't open file $idx_out\n";
print H join("\t", ("line", @snp_array)), "\n";
print D join("\t", ("line", @snp_array)), "\n";
foreach my $acc (keys %haplotype)
{
	my $hapref = $haplotype{$acc};
	my $idxref = $select{$acc};
	print H join("\t", ($acc, @{$hapref}{@snp_array})), "\n";
	print D join("\t", ($acc, @{$idxref}{@snp_array})), "\n";	
}
close H;
close D;

# 

sub read_phs_input
{
	my $file = shift;
	my %return;
	open (IN, $file) or die "can't open file \n";
	my $n = 1;
	my @array;
	while (<IN>)
	{
		chomp;
		my @t = split /\s+/, $_;
		if($n)
		{
			@array = @t;
			$n = 0;
			next;
		}
		
		foreach (1..$#t)
		{
			$return{$t[0]}{$array[$_]} = $t[$_] eq "NA"?"-":$t[$_];
		}
	}
	close IN;
	return %return;
}

sub read_snp_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die "can't open file \n";	
	while (<IN>)
	{
		chomp;
		my @t = split /\s+/, $_;
		$return{$t[1]} = $t[2];
	}
	close IN;
	
	return %return;
}

sub unique
{
	my %hash = map{$_, 1} @_;
	return keys %hash;
}

sub generate_hap_name
{
	my @arr = unique(@_);
	my $n = 0;
	my %name;
	foreach (1..(scalar @arr))
	{
		my $hap = $arr[$_-1];
		if ($hap =~ /\-/)
		{
			$name{$hap} = "NA";
		}
		else
		{
			$n++;
			$name{$hap} = "H".$n;
		}
	}
	return %name;
}


