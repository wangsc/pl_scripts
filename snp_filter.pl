#!/usr/bin/perl -w
use strict;


my $sUsage = qq(
perl $0
<all_var.filtered.vcf>
<min_coverage depth>
<output file>
<allele frequency file>
);
die $sUsage unless @ARGV >= 4;

my ($all_var_file, $min_coverage, $out_file, @allele_freq_files) = @ARGV;
die "!!Error: Min coverage should be integer $sUsage" if $min_coverage =~ /\D/;
my %raw_snps = read_var_file($all_var_file);
my %haplotypes = read_allele_freq_files(\%raw_snps, $out_file, $min_coverage, @allele_freq_files);

foreach my $id (keys %haplotypes)
{
	my @types = @{$haplotypes{$id}};
	my %count;
	map{$count{$_}++} @types;
	my @unique_types = sort{$count{$a}<=>$count{$b}}keys %count;
	print $id;
	foreach (keys %count){next if /NA/; print "\t", $_,"_",$count{$_}}
	print "\t", sprintf("%.1f", $count{$unique_types[0]}/(scalar @types));
	print "\n";
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
		$return{join("_", @data[0,1])} = [@data[3,4]] if $data[5]>=20;
	}
	close IN;
	return %return;
}

sub read_allele_freq_files
{
	my ($rawsnps, $out_file, $min_coverage, @files) = @_;
	my %return;
	open (OUT, ">$out_file") or die;
	#$min_coverage = 2;
	foreach my $f (@files)
	{
		print_time_comment("Processing file $f");
		open (IN, $f) or die $!;
		while(<IN>)
		{
			
			# Ok06114 contig1571304   83      C_0     T_1
			# Ok06114 contig1571304   61      A_1			
			chomp; 
			my @t = split /\t/,$_;
			next unless exists $rawsnps->{join("_", @t[1,2])};
			my %frequency;
			while(/(\S)_(\d+)/g)
			{
				#print $1,"\t",$2,"\n";
				$frequency{$1} = $2;
			}
			my @alleles = @{$rawsnps->{join("_", @t[1,2])}};
			#print '@alleles ', join("\t", @alleles),"\n";
			die join("_", @t[1,2]),"\n" unless @alleles == 2;
			my $type = determine_genotype(\@alleles, \%frequency, $min_coverage);
			push @{$return{join("\t", @t[1,2])}}, $type;
			print OUT join("\t", (@t[0..2], @alleles)),"\t", 
							$frequency{$alleles[0]},"\t", exists $frequency{$alleles[1]}?$frequency{$alleles[1]}:0,"\n";	
		}
		close IN;
	}
	close OUT;
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
