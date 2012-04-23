#!/usr/bin/perl -w
use strict;


my $sUsage = qq(
perl $0
<all_var.filtered.vcf>
<coverage depth, 1..25 or 4..8 or ...>
<all_snps_matched_9k.out>
<allele frequency file>
);
die $sUsage unless @ARGV >= 4;

my ($all_var_file, $coverage_range, $matched_9k_file, @allele_freq_files) = @ARGV;
my %raw_snps = read_var_file($all_var_file);
my %matched_9k = read_matched_9k_file($matched_9k_file);
my %allele_freq = read_allele_freq_files(\%raw_snps, @allele_freq_files);
my @range = ($1, $2) if $coverage_range =~ /(\d+)\.\.(\d+)/;
print join("\t", qw(Min_coverage Not_work Monomorphic Polymorphic)),"\n";
foreach my $min_coverage ($range[0]..$range[1])
{
	my $total_snps_12 = 0;
	my $total_snps_22 = 0 ;
	my @snp_label_count_12 = (0,0,0);
	my @snp_label_count_22 = (0,0,0);
	my %genotypes = genotype(\%allele_freq, $min_coverage);
	foreach my $id (keys %genotypes)
	{
		my @types = @{$genotypes{$id}};
		my %count;
		map{$count{$_}++} @types;
		next if scalar (keys %count) <= 1;
		$total_snps_12++;
		my @unique_types = sort{$count{$b}<=>$count{$a}}keys %count;
		$snp_label_count_12[$matched_9k{$id}]++ if exists $matched_9k{$id};
		my $c = 0;
		foreach (@unique_types)
		{
			$c++ if $count{$_} >= 2;
		}
		if($c>=2){$total_snps_22++; $snp_label_count_22[$matched_9k{$id}]++ if exists $matched_9k{$id}}
	}
	print join("\t", ($min_coverage, @snp_label_count_12, $total_snps_12)),"\n";
	print join("\t", ($min_coverage, @snp_label_count_22, $total_snps_22)),"\n";
}

# subroutines

sub read_matched_9k_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die $!;
	while(<IN>)
	{
		chomp; 
		next if /^\s+$/;
		my @t = split/\t/,$_;
		$return{join(":", @t[0,1])} = $t[2];
	}	
	close IN;
	return %return;
}


sub read_var_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die $!;
	while(<IN>)
	{
		# contig00066     679     .       A       C       55      .       DP=10;AF1=1;AC1=2;DP4=0,0,0,10;MQ=20;FQ=-57     GT:PL:GQ        1/1:88,30,0:57
		# contig00090     493     .       C       A       43      .       DP=9;AF1=1;AC1=2;DP4=0,0,9,0;MQ=20;FQ=-54       GT:PL:GQ        1/1:76,27,0:51
		# contig00101     170     .       T       C,A     37.3    .       DP=13;AF1=1;AC1=2;DP4=0,1,0,12;MQ=20;FQ=-41;PV4=1,0.12,1,1      GT:PL:GQ        1/1:71,15,1,56,0,61:25
    # contig00101     551     .       T       A       15.2    .       DP=6;AF1=0.5032;AC1=1;DP4=2,0,4,0;MQ=20;FQ=-8.64;PV4=1,6.4e-06,1,1      GT:PL:GQ        0/1:45,0,19:22		
		next if /^\#/;
		my @data = split /\t/,$_;
		next if length $data[4] > 1;
		$return{join(":", @data[0,1])} = [@data[3,4]] if $data[5]>=20;
	}
	close IN;
	return %return;
}

sub genotype
{
	my ($allele_freq, $min_dep) = @_;
	my %return;
	foreach my $id (keys %$allele_freq)
	{
		foreach (@{$allele_freq->{$id}})
		{
			my ($aa, $ab, $ca, $cb) = @$_;
			my %freq;
			$freq{$aa} = $ca; $freq{$ab} = $cb;
			my @alleles = ($aa, $ab);
			my $type = determine_genotype(\@alleles, \%freq, $min_dep);
			push @{$return{$id}}, $type unless $type eq 'NA';
		}
	}
	return %return;
}

sub read_allele_freq_files
{
	my ($rawsnps, @files) = @_;
	my %return;

	foreach my $f (@files)
	{
		#print_time_comment("Processing file $f");
		open (IN, $f) or die $!;
		while(<IN>)
		{
			
			# Ok06114 contig1571304   83      C    T 0 1
			# Ok06114 contig1571304   61      A  T  1  14			
			chomp; 
			my @t = split /\t/,$_;
			my $snpid = join(":", @t[1,2]);
			next unless exists $rawsnps->{$snpid};
			push @{$return{$snpid}}, [@t[3..6]]
		}
		return %return;
		close IN;
	}

	return %return;
}

sub determine_genotype
{
	my ($allele, $freq, $min_coverage) = @_;
	map{$freq->{$_}=0 unless exists $freq->{$_}}@$allele;
	my $total_covered_depth = $freq->{$allele->[0]} + $freq->{$allele->[1]};
	return 'NA' if $total_covered_depth < $min_coverage;

	if((abs($freq->{$allele->[0]} - $freq->{$allele->[1]})/$total_covered_depth) <= 0.8) # minor allele should be at least 10% of total alleles for hetero
	{
		return join("", @$allele);
	}
	else
	{
		return $freq->{$allele->[0]}>$freq->{$allele->[1]}?$allele->[0]:$allele->[1]
	}		

}

sub print_time_comment
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}
