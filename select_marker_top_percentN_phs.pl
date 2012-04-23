#!/usr/bin/perl -w
use strict;

# For the alleles which have top N percent PHS values, calculate the allele frequency in each sub population 
my $sUsage = "perl $0 <BScan_sp_pops.csv> <phs txt files>\n";
die $sUsage unless @ARGV;
my ($pop_file, @phs_files) = @ARGV;
my $percentile = 95;
my %phs = read_phs_files($percentile, @phs_files);
my ($allele_freq_hashref, $pops_ref) = read_pops_csv($pop_file);

# output
print join("\t", ("Chr\tSNP\tPHS", @$pops_ref)), "\n";
foreach my $chr (sort {$a cmp $b}keys %phs)
{
	foreach (@{$phs{$chr}})
	{
		my ($snp, $phs_value, $allele_index) = @$_;
		my @freq;
		foreach (@$pops_ref)
		{
			push @freq, $allele_freq_hashref->{$snp}->{$_}->[$allele_index -1];
		}
		print join("\t", ($chr, $snp, $phs_value, @freq)), "\n";
	}
}



# Subroutines
sub read_pops_csv
{
	my $file = shift;
	open (IN, $file) or die ;
	my %return;
	my %pops;
	while (<IN>)
	{
		chomp;
		s/\"//g;
		next if /^locus/i;
		# "locus","size","alNum","al1","al2","popNum","locusName","popName"
		# 1,322,2,138,184,1,"wsnp_AJ612027A_Ta_2_1","Australia_spring"
		my @t = split /,/, $_;
		my ($total_allele, $allele_a, $allele_b, $snp, $pop) = @t[1, 3,4, 6,7];
		$return{$snp}{$pop} = [0, 0];
		$return{$snp}{$pop} = [$allele_a/$total_allele, $allele_b/$total_allele] if $total_allele;
		$pops{$pop} = 1;
	}
	close IN;
	foreach my $snp (keys %return)
	{
		map{ $return{$snp}{$_}=["NA", "NA"] unless exists $return{$snp}{$_} } keys %pops;
	}
	
	return (\%return, [sort{$a cmp $b} keys %pops]);
}

sub read_phs_files
{
	my $percentile = shift;
	my @files = @_;
	my %return;
	my %percentile_values;
	foreach my $file (@files)
	{
		my $chr = $1 if $file =~/(\S{2})_\S+?_phs.txt/;
		my @phs;
		open (IN, $file) or die $!;
		while (<IN>)
		{
			# Chr	SNP_index	SNP_name	Genetic_dist	PHS_A	PHS_B	Freq_A	Freq_B	Block_A	Block_B
			next if /^chr/i;
			chomp;
			my @data = split /\s+/,$_;
			push @phs, (@data[4, 5]);
			push @{$return{$chr}}, [@data[2, 4, 5]];
		}
		close IN;
		my $percentile_value = get_percentile_value(@phs, $percentile);
		$percentile_values{$chr} = $percentile_value;
	}
	
	foreach my $chr (keys %return)
	{
		my $percentile_value = $percentile_values{$chr};
		my @data = @{$return{$chr}};
		my @filtered;
		foreach (@data)
		{
			my ($snp, $phsa, $phsb) = @$_;
			
			if($phsa > $percentile_value)
			{
				push @filtered, [$snp, $phsa, 1];
			}
			if($phsb > $percentile_value)
			{
				push @filtered, [$snp, $phsb, 2];
			}				
		}
		$return{$chr} = [@filtered];
	}
	return %return;
}

sub get_percentile_value
{
	my @values = @_;
	my $percent = pop @values;
	
	@values = sort {$b <=> $a} @values;
	my $index = int((scalar @values) * (100-$percent)/100);
	$index = 0 if $index < 0;
	return $values[$index]; 
}

