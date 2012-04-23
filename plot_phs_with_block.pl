#!/usr/bin/perl -w
use strict;

die <<END unless @ARGV;

# This script will read one output file of calculate_PHS.pl, <prefix>_phs.txt, then invoke R to make a plot
# Make sure R is properly installed.
Usage: perl $0 <phs.txt> <1A_MAGIC_MAP.out> <chr, optional>
END

my $phs_file = shift;
my $map_file = shift;
my $chromosome = shift;
my @mapped_positions = read_map_file($map_file);
my $phs_cutoff = 0; # default cutoff for phs value is 0
open (IN, $phs_file) or die "can't open file $phs_file\n";

my @data_points;
my @phs_values;
while(<IN>)
{
	next if /^chr/i; # skip the header line
	chomp;
	my @data = split /\s+/,$_;
	# Chr     SNP_index       SNP_name   genetic_dist     PHS_A   PHS_B   Freq_A  Freq_B  Block_A   Block_B
	# 0				1                 2               3         4       5       6        7         8      9
	my ($chr, $gd, $phs_a, $phs_b, $freq_a, $freq_b, $block_a, $block_b) = @data[0, 3..9];
	push @phs_values, ($phs_a, $phs_b);
	$chromosome = $chr unless defined $chromosome;
	if($phs_a > $phs_cutoff)
	{
		push @data_points, [$gd, $freq_a, (split /_/, $block_a), $phs_a]
	}
	
	if($phs_b > $phs_cutoff)
	{
		push @data_points, [$gd, $freq_b, (split /_/, $block_b), $phs_b]
	}	
}
close IN;

my $percentile = 97;
my $phs_perc_value = get_percentile(@phs_values, $percentile);

# Generate R script
my $r_script = $chromosome . "_temp.R";
open (TR, ">$r_script") or die "can't open file $r_script\n";
print TR '#!/usr/bin/env Rscript', "\n";
my $pdffile = 'Chromosome_' . $chromosome . ".pdf";
print TR "pdf(\"$pdffile\")\n";
my @gd = map{$_->[2], $_->[3]} @data_points;
my ($min, $max) = min_max(@gd);
print TR "plot(1, type=\"n\", axes=T, xlim=c(0, $max), xlab=\"Position (cM)\", ylim=c(-0.02, 1), ylab=\"Frequency\")\n";
my @pointx;
my @pointy;
my @point_color;
my @line_seg;
foreach (@data_points)
{
	#print join(" ", @$_), "\n" unless $_->[0] >= $_->[2] and $_->[0] <= $_->[3];
	push @pointx, $_->[0];
	push @pointy, $_->[1];
	my $col = $_->[4] >= $phs_perc_value?"\"red\"":"\"black\"";
	push @point_color, $col;
	push @line_seg, [$_->[2], $_->[1], $_->[3], $_->[1], $_->[4]] unless $_->[2] eq 'NA' or $_->[3] eq 'NA';
}
print TR "x<-c(", join(",", @pointx), ")\n";
print TR "y<-c(", join(",", @pointy), ")\n";
print TR "p_col <- c(", join(",",@point_color), ")\n";
print TR "points(x, y, pch=20, cex=0.5, col=p_col)\n";
my $count = 0;
foreach (@line_seg)
{
	my @arr = @$_;
	my $phs = pop @arr;
	$count++ if $phs >= $phs_perc_value;
	my $col = $phs >= $phs_perc_value?"\"red\"":"\"black\"";
	print TR "segments(", join(",", @arr), ",col=", $col, ")\n";
}
print "Chromosome ", $chromosome, "\t", $count, "\n";
print TR "rugx <- c(", join(",", @mapped_positions), ")\n";
print TR "rug(rugx, col=\"blue\")\n";
print TR "dev.off()\n";

close TR;

# Run R script
my $cmd = "R CMD BATCH $r_script";
print "Running: $cmd", "\n";
print "Command: $cmd failed\n"if (system($cmd));

# subroutines
sub get_percentile
{
	my @values = @_;
	my $percent = pop @values;
	
	@values = sort {$b <=> $a} @values;
	my $index = int((scalar @values) * (100-$percent)/100);
	$index = 0 if $index < 0;
	return $values[$index]; 
}


sub read_map_file
{
	my $file = shift;
	open (IN, $file) or die;
	my @return;
	while (<IN>)
	{
		chomp;
		next if /^Chr/;
		my @t= split /\s+/,$_;
		push @return, $t[-1];
	}
	close IN;
	return @return;
}

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




