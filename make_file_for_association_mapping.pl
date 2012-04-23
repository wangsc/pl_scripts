#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
Usage: (use >outfilename to redirect the output)
perl $0
<association_sample_list>
<ksu_sample_sheet.formatted>
<called genotype file>
<total number of SNPs used to genotype>
<show genotype as A:B (0) or 1:2 (1)>
);
die $sUsage unless @ARGV >= 5;

my($sample_list, $ksu_sample, $genotype_file, $total_snp, $format) = @ARGV;
my @ass_sample;
my %association_sample = extract_association_sample($sample_list, $ksu_sample, \@ass_sample);
my %genotype_data = read_genotype_file($genotype_file);
output(\%genotype_data, \%association_sample, \@ass_sample, $total_snp, $format);

# 
sub output
{
	my ($genotype_ref, $asso_hashref, $ass_arrayref, $total_snp, $format) = @_;
	foreach (1 .. $total_snp){print "\t", '9kIdnex_',$_}; print "\n";
	foreach my $name (@$ass_arrayref)
	{
		print $name;
		if (exists $asso_hashref->{$name})
		{
			my $sample_id = $asso_hashref->{$name};
			foreach my $index (1..$total_snp)
			{
				my $g = (defined $genotype_ref->{$sample_id}->[$index])?$genotype_ref->{$sample_id}->[$index]:'??';
				$g =~ s/(\w)(\w)/$1:$2/;
				$g = transform($g) if $format;
				print "\t", $g;
			}
			print "\n";
		}
		else
		{
			foreach my $index (1..$total_snp){print "\t", '?:?'}
			print "\n";
		}
	}
}

sub transform
{
	my $g = shift;
	return '?:?' if $g =~ /\?/;
	return '1:1' if $g =~ /A/;
	return '2:2';
}


sub read_genotype_file
{
	my $file = shift;
	my %return;
	open (IN, "$file") or die "can't open file $file\n";
	while (<IN>)
	{
		# wsnp_AJ612027A_Ta_2_1 5904341053_R01C01 A A 0.9289 1 A A A A A A - - 1.0000 AA AA
		chomp;
		next if /^\s+$/;
		my @t = split /\s+/, $_;
		my ($sample_id, $snp_index, $genotype) = @t[1, 5, 16];
		$genotype = '??' if is_hetero($genotype);
		$return{$sample_id}->[$snp_index] = $genotype;
	}
	close IN;
	return %return;
}

sub is_hetero
{
	my $g = shift;
	my @t = split //, $g;
	return ($t[0] eq $t[1])?0:1;
}


sub extract_association_sample
{
	my ($ass_file, $ksu_file, $ass_arrayref) = @_;
	my %return;
	my %ass_sample;
	open (IN, "$ass_file") or die "can't open file $ass_file\n";
	while (<IN>)
	{
		chomp; 
		next if /^\s+/;
		s/\s//g;
		$ass_sample{$_} = 1;
		push @{$ass_arrayref}, $_;
	}
	close IN;
	my %ksu_sample;
	open (IN, "$ksu_file") or die $!;
	while (<IN>)
	{
		chomp;
		my @t = split /\t/, $_;
		$t[0] =~ s/\s//g;
		$ksu_sample{$t[0]} = $t[-1];
	}
	close IN;
	foreach my $ass(keys %ass_sample)
	{
		#my $matched = 0;
		#print STDERR '$ass ', $ass,"\n";
		foreach my $ksu (keys %ksu_sample)
		{
			#print STDERR $ksu,"\n";
			$return{$ass} = $ksu_sample{$ksu} if $ksu =~ /^$ass/; 
		}
	}
	return %return;
}
