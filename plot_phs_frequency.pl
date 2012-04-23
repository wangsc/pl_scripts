#!/usr/bin/perl -w
use strict;
die <<END unless @ARGV;

# This script will read one output file of calculate_PHS.pl, <prefix>_phs.txt, then invoke R to make a plot
# Make sure R is properly installed.
Usage: perl $0 <*_phs.txt>
END

my @phs_files = @ARGV;

my @data_points = read_phs_files(@phs_files);
@data_points = sort {$a->[0] <=> $b->[0]} @data_points; # sort by PHS values in increasing order
my @percentiles = (95, 97.5, 99);

my @percent_data;
my %points;
map{push @{$points{sprintf("%.2f", $_->[0])}}, $_->[1] if $_->[1] > 0}@data_points;

foreach my $freq (sort {$a <=> $b} keys %points)
{
	my @phs = sort{$a <=> $b}@{$points{$freq}};
	my @percent_index = map{int($_*(scalar @phs)/100)-1}@percentiles;
	my @percent_phs = map{$phs[$_]}@percent_index;
	push @percent_data, [$freq, @percent_phs];
}

my $r_script = "phs_freq.R";
open (R, ">$r_script") or die;
print R '#!/usr/bin/env Rscript', "\n";
my $pdffile = 'phs_freq' . ".pdf";
print R "pdf(height=5, width=7, file=\"$pdffile\")\n";
my ($phs_min, $phs_max) = min_max(map{$_->[1]}@data_points);
print R "plot(1, type=\"n\", axes=T, xlim=c(0, 1), xlab=\"Frequency\", ylim=c($phs_min, $phs_max), ylab=\"PHS\")\n";
print R "freq<-c(", join(",", (map{$_->[0]}@data_points)), ")\n";
print R "phs<-c(", join(",", (map{$_->[1]}@data_points)), ")\n";
print R "points(freq, phs, pch=20, cex=0.3, col=\"grey\")\n";
my @avg_freq = map {$_->[0]} @percent_data;
print R "avgfreq <- c(", join(",", @avg_freq), ")\n";
my @cols = qw(blue darkgreen red);
foreach my $index (1..scalar @percentiles)
{
	#next unless $index ==2;
	my $name = "P".$index;
	print R "$name <- c(", join(",", (map{$_->[$index]}@percent_data)), ")\n";
	print R "lines(loess.smooth(avgfreq, $name), col=\"$cols[$index-1]\")\n";
}

print R "dev.off()\n";

# Run R script
my $cmd = "R CMD BATCH $r_script";
print "Running: $cmd", "\n";
print "Command: $cmd failed\n" if (system($cmd));

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
			push @return, ([@data[6, 4]], [@data[7, 5]]);
		}
		close IN;
	}
	return @return;
}

sub average {
	my $count = &count($_[0]);
	if ($count == 0) {
		return 0;
	} else {
		return (&sum($_[0]) / $count);
	}
}

sub sum {
	my $arr_ref = $_[0];

	my $sum = 0;
	foreach (@$arr_ref) {
		next unless defined($_);
		next if /^$/;
		$sum += $_;
	}
	return $sum;
}

sub count {
	my $arr_ref = $_[0];
	
	my $num = 0;
	foreach (@$arr_ref) {
		next unless defined($_);
		next if /^$/;
		$num++;
	}
	return $num;
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