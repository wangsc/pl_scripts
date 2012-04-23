#!/usr/bin/perl -w

# detect local maximal and minimal value
my $file = shift or die "perl $0 <data file>\n";
my @array = read_file($file);
my @peaks = peak_detect(\@array, 0.1);
my $cutoff_percent = 0.9;
my @filtered_peaks = filter_peaks($cutoff_percent, @peaks);

output(@peaks);
output(@filtered_peaks);

sub output
{
	foreach (@_)
	{
		print join("\t", @$_),"\n";
	}
	print '*'x10,"\n";
}

sub read_file
{
	my $file = shift;
	my @array;
	open (IN, $file) or die;
	while(<IN>)
	{
		chomp;
		next if /^\s+$/;
		my @t= split /\t/,$_;
		$array[$t[0]] = $t[1];
		#print '!! ', $t[0], "\n";
	}
	close IN;
	return @array;
}

sub filter_peaks
{
	my ($cut_percent, @peaks) = @_;
	my @filtered;
	my $total = 0;
	my @peak_values;
	foreach (@peaks)
	{
		$total += $_->[1];
	}
	my $cutoff = ($total/(scalar @peaks))*$cut_percent;
	foreach (@peaks)
	{
		push @filtered, $_ if $_->[1] >= $cutoff;
	}
	return @filtered;
}



sub peak_detect
{
	my ($array_ref, $delta) = @_;
	my @peaks;
	map{$array_ref->[$_] = 0 unless defined $array_ref->[$_]}(0..(scalar @$array_ref));
	my ($min, $max) = (99999, -99999);
	my $maxpos;
	my $minpos;
	my $max_flag = 1;
	foreach my $i (0 .. (scalar @$array_ref - 1))
	{
		my $current = $array_ref->[$i];
		die "Not defined \n" unless defined $current;
		if ($current > $max){$max = $current; $maxpos = $i}
		if ($current < $min){$min = $current; $minpos = $i}		
		if ($max_flag)
		{
			if($current < ($max*(1-$delta)))
			{
				print $i,' ** ', $current," ** ", $max,"\n";
				push @peaks, [$maxpos, $max];
				$minpos = $i;
				$min = $current;
				$max_flag = 0;
			}		
		}
		else
		{
			if ($current > $min*(1+$delta))
			{
				$max = $current;
				$maxpos = $i;
				$max_flag = 1;
			}
		}

	}
	return @peaks;
}