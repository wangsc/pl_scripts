#!/usr/bin/perl -w
use strict;


my $sUsage = qq(
perl $0
<Min num of accessions covered>
<all_unique_snps_matched_9k.out>
<count_called_genotype.out>
);
die $sUsage unless @ARGV >= 3;

my ($min_acc, $matched_9k_file, $genotype_count_file) = @ARGV;
my %matched_9k = read_matched_9k_file($matched_9k_file);

print join("\t", qw(Min_acc Not_work Monomorphic Polymorphic Total_snps)),"\n";
my @counts = (0, 0, 0);
my $total_snps; 
open (IN, $genotype_count_file) or die;
while(<IN>)
{
	chomp;
	my @t = split/\t/,$_;
	next unless @t == 4;
	my $id = join(":", @t[0,1]);
	my @c;
	foreach (2..3)
	{
		push @c, $1 if $t[$_]=~/_(\d+)/;
	}
	@c = sort{$a<=>$b} @c;
	if($c[0] >= $min_acc)
	{
		$total_snps++;
		$counts[$matched_9k{$id}]++ if exists $matched_9k{$id};
	}
}
print join("\t", ($min_acc, @counts, $total_snps)),"\n";

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
	my ($allele, $freq, $min_dep) = @_;
	my $min = $min_dep;
	$min_dep = 3;
	map{$freq->{$_}=0 unless exists $freq->{$_}}@$allele;
	return 'NA' if $freq->{$allele->[0]} < $min_dep and $freq->{$allele->[1]} < $min_dep;
	return join("", sort{$a cmp $b} @$allele) if $freq->{$allele->[0]} >= $min_dep and $freq->{$allele->[1]} >= $min_dep;
	
	my @alleles = sort{$freq->{$a}<=>$freq->{$b}} @$allele;
#	if ($freq->{$alleles[0]}/($freq->{$alleles[0]}+$freq->{$alleles[1]}) <= $min_percent/100)
#	{
#		return $alleles[1];
#	}
#	else
#	{
#		return join("", sort{$a cmp $b} @$allele)
#	}

	if($freq->{$alleles[0]} == 0 and $freq->{$alleles[1]} >= $min){return $alleles[1]} else{return 'NA'}
}

sub print_time_comment
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}
