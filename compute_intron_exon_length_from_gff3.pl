#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <gff3 file> <genewise, 0 or pasa, 1>\n";
die $sUsage unless @ARGV == 2;
my $gff_file = shift;
my $is_pasa = shift;
my ($intron_len, $exon_len) = $is_pasa?read_pasa_gff($gff_file):read_genewise_gff($gff_file);
print "** Intron length **\n";
output ($intron_len);
print "** Exon length **\n";
output ($exon_len);
print "***"x8,"\n";
my @plots = (50,100,200, 500, 1000, 1500, 2000, 2500, 3000);
my @intron_dis = calculate_distribution($intron_len, \@plots);
my @exon_dis = calculate_distribution($exon_len, \@plots);
output(\@intron_dis); print "***"x8,"\n";
output(\@exon_dis);


# Subroutines

sub calculate_distribution
{
	my $hashref = shift;
	my $plot_ref = shift; 
	my @plots = @$plot_ref;
	my $num = scalar @plots;
	my @record;
	foreach my $len (keys %$hashref)
	{
		foreach my $ind (0..$num)
		{
			$record[$ind] = 0 unless defined $record[$ind];
			if ($ind ==0)
			{
				 if ($len < $plots[0]){$record[$ind] += $hashref->{$len}; last}
			}
			elsif ($ind == $num)
			{
				if ($len >= $plots[$num-1]){$record[$ind] += $hashref->{$len}; last}
			}
			else
			{
				$record[$ind] += $hashref->{$len} if $len >= $plots[$ind-1] and $len < $plots[$ind];
			}
		}
	}
	return @record;
}

sub output
{
	my $inref = shift;
	if (ref($inref) eq 'HASH')
	{
		foreach (sort {$a<=>$b} keys %$inref)
		{
			print $_, "\t", $inref->{$_},"\n";
		}		
	}
	elsif (ref($inref) eq 'ARRAY')
	{
		print join("\n", @$inref),"\n";
	}

}

sub read_genewise_gff
{
	my $file = shift;
	open (IN,"$file") or die "can't open file $file\n";
	my (%intron, %exon);
	my $score = 0;
	my $score_cutoff = 35;
	while (<IN>)
	{
		#asmbl_3904	GeneWise	match	3336	8227	817.94	+	.	asmbl_3904-genewise-prediction-1
		#asmbl_3904	GeneWise	cds	3336	3454	0.00	+	0	asmbl_3904-genewise-prediction-1
		#asmbl_3904	GeneWise	intron	3455	3550	0.00	+	.	asmbl_3904-genewise-prediction-1
		next if /^\s+$/ or /^\//;
		chomp;
		my @data = split /\t/, $_;
		my ($type, $start, $end) = @data[2..4];
		if ($type =~ /match/){$score = $data[5]; next}
		my $length = abs($end-$start) + 1;
		if ($type =~ /cds/)
		{
			$exon{$length}++ if $score > $score_cutoff
		}
		elsif ($type =~ /intron/)
		{
			$intron{$length}++ if $score > $score_cutoff
		}
	}
	close IN;
	return (\%intron, \%exon)
}

sub read_pasa_gff
{
	my $file = shift;
	open (IN,"$file") or die "can't open file $file\n";
	my (%intron, %exon);
	my %target;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my $id = $1 if /Target=(\S+)/;
		my @data = split /\t/, $_;
		push @{$target{$id}}, [@data[3, 4]];
		my $length = abs($data[3] - $data[4]) + 1;
		$exon{$length}++;
	}
	close IN;
	foreach (keys %target)
	{
		my $num_exons = scalar @{$target{$_}};
		my @array = sort {$a->[0] <=> $b->[0]} @{$target{$_}};
		next unless ( $num_exons > 1);
		foreach my $index (0..($num_exons - 2))
		{
			my $len = abs($array[$index]->[1] - $array[$index+1]->[0]) + 1;
			$intron{$len}++;
		}
	}
	return (\%intron, \%exon);
}


