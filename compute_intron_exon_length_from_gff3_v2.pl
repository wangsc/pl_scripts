#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <gff3 file> < 0 for normal gff or 1 for pasa gff>\n";
die $sUsage unless @ARGV == 2;
my $gff_file = shift;
my $is_pasa = shift;
my ($intron_len, $exon_len) = $is_pasa?read_pasa_gff($gff_file):read_normal_gff($gff_file);
print "** Intron length **\n";
output ($intron_len);
print "** Exon length **\n";
output ($exon_len);
print "***"x8,"\n";
my @plots = (50,100,200, 500, 1000, 1500, 2000, 2500, 3000);
print join("\t", @plots),"\n","***"x8,"\n";
print 'Intron: ',"\n";
print 'Median and Mean lenght: ', join("\t", calcualte_median_and_mean($intron_len)),"\n";
my @intron_dis = calculate_distribution($intron_len, \@plots);
my @exon_dis = calculate_distribution($exon_len, \@plots);
output(\@intron_dis); 
print "***"x8,"\n";
print 'Exon: ',"\n";
print 'Median and Mean lenght: ', join("\t", calcualte_median_and_mean($exon_len)),"\n";
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

sub read_normal_gff
{
	my $file = shift;
	open (IN,"$file") or die "can't open file $file\n";
	my (%intron, %exon);
	my %exon_record;
	my $score = 0;
	while (<IN>)
	{
		next if /^\s+$/ or /^\#/;
		chomp;
		my @data = split /\t/, $_;
		my ($type, $start, $end) = @data[2..4];
		if ($type =~ /match/){$score = $data[5]; next}
		my $length = abs($end-$start) + 1;
		my $gene = $1 if /Parent=(\S+)/;
		
		if ($type =~ /cds/i)
		{
			$exon{$length}++;
			push @{$exon_record{$gene}}, [sort{$a<=>$b}($start, $end)];
		}
	}
	foreach (values %exon_record)
	{
		my @array = @$_;
		@array = sort{$a->[0]<=>$b->[0]} @array;
		foreach (0..($#array - 1))
		{
			my $length = $array[$_+1]->[0] - $array[$_]->[1] - 1;
			$intron{$length} ++;
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

sub calcualte_median_and_mean
{
	my $hashref = shift;
	my $total = 0;
	my $median; my $mean;
	map{$total += $_} values %$hashref;
	my $length_total = 0;
	my $count = 0;
	foreach (sort{$a<=>$b} keys %$hashref)
	{
		$length_total += $_*$hashref->{$_};
		$count += $hashref->{$_};
		$median = $_ if $count >= int($total/2) and (not defined $median);
	}
	$mean = int($length_total/$total);
	return($median, $mean);
}



