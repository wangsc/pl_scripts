#!/usr/bin/perl -w
use strict;


my $sUsage = qq(
perl $0
<all_var.filtered.vcf>
<Min coverage depth range N..M >
<allele frequency file>
);
die $sUsage unless @ARGV >= 3;

my ($all_var_file, $percent_range,  @allele_freq_files) = @ARGV;
print_time_comment("processing $all_var_file");
my %raw_snps = read_var_file($all_var_file);
print_time_comment("processing allele frequency file");
my %allele_freq = read_allele_freq_files(\%raw_snps, @allele_freq_files);
my @range = ($1, $2) if $percent_range =~ /(\d+)\.\.(\d+)/;
#print join("\t", qw(Min_dep Not_work Monomorphic Polymorphic)),"\n";
foreach my $min_depth ($range[0]..$range[1])
{
	print_time_comment("Calculating for min coverage depth $min_depth ...");
	my $out = "dp_". $min_depth . ".out";
	open (OUT, ">$out") or die "can't open file $out \n";
	my %genotypes = genotype(\%allele_freq, $min_depth);
	foreach my $id (keys %genotypes)
	{
		my @types = @{$genotypes{$id}};
		my %count;
		map{$count{$_}++} @types;
		my @temp; map{push @temp, join("_", ($_, $count{$_}))} keys %count;
		print OUT join("\t", ($id, @temp)),"\n";
	}
	close OUT;
}

# subroutines

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
			
			# Ok06114 contig1571304:83      C    T 0 1
			# Ok06114 contig1571304:61      A  T  1  14			
			chomp; 
			my @t = split /\t/,$_;
			my $snpid = $t[1];
			next unless exists $rawsnps->{$snpid};
			my @alleles = @{$rawsnps->{$snpid}};
			my @count;
			if(@t==6) # Ok06114 contig1571304:83      C T 0 1
			{
				my %temp = @t[2,4,3,5];
				@count = map{$temp{$_}} @alleles;
				my %check_num_allele = map{$_,1} (@alleles, @t[2,3]);
				next if (scalar keys %check_num_allele) > 2;
			}
			elsif(@t==4) # Ok06114 contig1571304:80      C  1
			{
				my %temp = @t[2, 3];
				@count = map{exists $temp{$_}?$temp{$_}:0}@alleles;
				my %check_num_allele = map{$_,1} (@alleles, $t[2]);
				next if (scalar keys %check_num_allele) > 2;			
			}
			push @{$return{$snpid}}, [@alleles, @count]
		}
		return %return;
		close IN;
	}

	return %return;
}

sub determine_genotype
{
	my ($allele, $freq, $min_dep) = @_;
	#my $min = $min_dep;
	my $min_coverage_call_AB = 3; # call AB: A>=3 and B>=3; 
	map{$freq->{$_}=0 unless exists $freq->{$_}}@$allele;
	return 'NA' if $freq->{$allele->[0]} < $min_coverage_call_AB and $freq->{$allele->[1]} < $min_coverage_call_AB;
	return join("", sort{$a cmp $b} @$allele) if $freq->{$allele->[0]} >= $min_coverage_call_AB and $freq->{$allele->[1]} >= $min_coverage_call_AB;
	
	my @alleles = sort{$freq->{$a}<=>$freq->{$b}} @$allele;
	# call A: A>=$min_dep and B=0; call B: B >= $min_dep and A=0
	if($freq->{$alleles[0]} == 0 and $freq->{$alleles[1]} >= $min_dep){return $alleles[1]} else{return 'NA'}
}

sub print_time_comment
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}
