#!/usr/bin/perl -w
use strict;

#
my $file = shift or die "perl $0 [96selected.all.samples.out] [theat cutoff , default:0.4]\n";
my %sample_theta_each_snp = read_file($file);
my $theta_cutoff = shift;
$theta_cutoff = 0.4 unless defined $theta_cutoff;

my %pairwise_theta_diff;
foreach my $snpid (keys %sample_theta_each_snp)
{
	my @theta = @{$sample_theta_each_snp{$snpid}};
	my $total_sample = scalar @theta;
	foreach my $first_sam (1..$total_sample-1)
	{
		foreach my $second_sam ($first_sam+1.. $total_sample)
		{
			$pairwise_theta_diff{join('_', ($first_sam, $second_sam))} = [] unless exists $pairwise_theta_diff{join('_', ($first_sam, $second_sam))};
			if (abs($theta[$first_sam-1]-$theta[$second_sam-1])>=$theta_cutoff)
			{				
				push @{$pairwise_theta_diff{join('_', ($first_sam, $second_sam))}}, $snpid;
			}
		}
	}
}

foreach my $pair (keys %pairwise_theta_diff)
{
	#$pair =~ s/_/\t/;
	print join("\t", ($pair, @{$pairwise_theta_diff{$pair}}, scalar @{$pairwise_theta_diff{$pair}})),"\n";
}



sub read_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file \n";
	while (<IN>)
	{
		next if /^Index/;
		my @line_data = split /\t/, $_;
		my $snpid = $line_data[1];
		my $index = 12;
		while ($index < (scalar @line_data - 1))
		{
			push @{$return_hash{$snpid}}, $line_data[$index] if defined $line_data[$index];
			$index += 4;
		}
	}
	close IN;
	map { $return_hash{$_} = [sort {$a<=>$b} @{$return_hash{$_}}] } (keys %return_hash);
	return %return_hash;
}

sub min_max
{
	my $min = shift;
	my $max = $min;
	map{$min = $_ if $_<$min; $max = $_ if $_>$max}@_;
	return($min, $max);
}