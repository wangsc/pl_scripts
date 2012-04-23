#!/usr/bin/perl -w
use strict;

my $sUsage = "
perl $0
<combined_all_07.14.2011.csv>
<4k_unique_USDA_Result.csv>
";
die $sUsage unless @ARGV >= 2;
my($k9_file, $k4_file) = @ARGV;
my @names = read_k9_csv($k9_file);
read_k4_csv($k4_file, @names);


sub read_k4_csv
{
	my $file = shift;
	my @k9_names = @_;
	open (IN, $file) or die $!;
	my $lines = 0;
	my %names;
	my %name_index;
	my %return;
	my $snp_start = 9000;
	while(<IN>)
	{
		chomp;
		$lines++;
		s/\"//g;
		if ($lines == 2)
		{
			my @names = split /\t/,$_;
			foreach (0..$#names){$names[$_] =~ s/\s+//g;}
			%name_index = map{$names[$_],$_} (0..$#names);
			map{$names{$_}=1} @names;
			foreach (@k9_names)
			{
				print STDERR $_,"\n" unless exists $names{$_}
			}
			#exit;
		}
		next if $lines <= 4;
		$snp_start++;
		my @data = split /\t/, $_;
		print $snp_start,"\t", $data[1];
		foreach my $name (@k9_names)
		{
			if(exists $name_index{$name})
			{
				print "\t", $data[$name_index{$name}];
			}
			else
			{ print "\t--"}
		}
		print "\n";
	}
	close IN;
	return 1;	
}


sub read_k9_csv
{
	my $file = shift;
	open (IN, $file) or die $!;
	my(@heads, @names);
	while (<IN>)
	{
		s/\"//g;
		if(/^Index/)
		{
			chomp;
			my @data = split /\t/, $_;
			@names = @data[2..$#data];
			last;
		}
	}
	close IN;
	foreach (0..$#names){$names[$_] =~ s/\s//g;}
	return @names;
}