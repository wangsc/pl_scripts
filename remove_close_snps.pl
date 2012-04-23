#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<snp file>
<min_distance between snps>
);
die $sUsage unless @ARGV >=2;
my($in_file, $min_dist) = @ARGV;

my($selected_snp_pos, $snp_record) = read_snp_file($in_file, $min_dist);

foreach my $key (keys %{$selected_snp_pos})
{
	die '$key: ', $key unless exists $selected_snp_pos->{$key};
	my @array = @{$selected_snp_pos->{$key}};
	foreach (@array)
	{
		my $id = join(":", ($key, $_));
		print $snp_record ->{$id};
	}	
}


sub read_snp_file
{
	my ($file, $min) = @_;
	open (IN, $file) or die $!;
	my (%snp_pos, %selected_pos);
	my %recorder;
	while(<IN>)
	{
		my @t= split /\t/,$_;
		my $id = join(':', @t[0,1]);
		push @{$snp_pos{$t[0]}}, $t[1];
		$recorder{$id} = $_;
	}
	close IN;
	foreach my $key (keys %snp_pos)
	{
		my @array = sort{$a<=>$b} @{$snp_pos{$key}};
		my $index = 0;
		$selected_pos{$key} = [@array] if @array == 1;
		while($index < $#array)
		{
			if (abs($array[$index] - $array[$index+1])< $min)
			{
				$index += 2;
				next
			}
			else
			{
				push @{$selected_pos{$key}}, $array[$index];
				$index++;
			}
			push @{$selected_pos{$key}}, $array[$index] if ($index == $#array and ($array[$index] - $array[$index-1])> $min);
		}
	}
	return(\%selected_pos, \%recorder)
}