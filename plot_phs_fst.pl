#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <snp_fst.out> <phs txt files>\n";
die $sUsage unless @ARGV;

my ($fst_file, @phs_files) = @ARGV;
my %fst_values = read_fst_file($fst_file);
my $percentile = 95;
my @phs_values = read_phs_files(@phs_files);
my @phs = map{$_->[1], $_->[2]}@phs_values;
@phs = sort {$a<=>$b} @phs;
my $percent_index = int((scalar @phs)*$percentile/100) - 1;
my $percent_value = $phs[$percent_index];
my @fst_upper;
my @fst_lower;

foreach (@phs_values)
{
	my ($snp, $p1, $p2) = @$_;
	if($p1>$percent_value or $p2>$percent_value)
	{
		push @fst_upper, $fst_values{$snp};
	}
	else
	{
		push @fst_lower, $fst_values{$snp};
	}
}

# plot
my $r_script = "phs_fst.R";
open (R, ">$r_script") or die;
print R '#!/usr/bin/env Rscript', "\n";
my $pdffile = 'phs_fst' . ".pdf";
print R "pdf(\"$pdffile\")\n";
my ($min, $max) = min_max(@fst_upper, @fst_lower);
print R "Upper <- c(", join(",", @fst_upper), ")\n";
print R "Lower <- c(", join(",", @fst_lower), ")\n";
print R "boxplot(list(Upper=Upper, Lower=Lower), ylab=\"Fst\")\n";

print R "dev.off()\n";

# Run R script
my $cmd = "R CMD BATCH $r_script";
print "Running: $cmd", "\n";
print "Command: $cmd failed\n" if (system($cmd));

# Subroutiens
sub min_max
{
	my $min = shift;
	my $max = $min;
	foreach (@_)
	{
		next if $_ eq 'NA';
		$min = $_ if $min > $_;
		$max = $_ if $max < $_;
	}
	return ($min, $max);	
}

sub read_fst_file
{
	my $file = shift;
	open (IN, $file) or die $!;
	my %return;
	while(<IN>)
	{
		chomp; next if /^\s+$/;
		my @t = split /\s+/,$_;
		$return{$t[0]} = $t[1];
	}
	close IN;
	return %return;
}

sub read_phs_files
{
	my @files = @_;
	my @return;
	foreach my $file (@files)
	{
		open (IN, $file) or die $!;
		while (<IN>)
		{
			# Chr	SNP_index	SNP_name	Genetic_dist	PHS_A	PHS_B	Freq_A	Freq_B	Block_A	Block_B
			next if /^chr/i;
			chomp;
			my @data = split /\s+/,$_;
			push @return, ([@data[2, 4, 5]]);
		}
		close IN;
	}
	return @return;
}