#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <345 combined csv file> <deletion_list file>";
die $sUsage unless @ARGV >= 2;
my $csv_file = shift;
my $del_list_file = shift;
my $num_missed_sample = 100;
my %del_lists = read_delete_list($del_list_file);
my ($max_sample, %records) = read_kian_csv_files($csv_file);
foreach (keys %records)
{
	next if exists $del_lists{$_};
	next if (scalar @{$records{$_}}) < ($max_sample-$num_missed_sample);
	print join("\n", @{$records{$_}})
}




sub read_kian_csv_files
{
	my @files = @_;
	my %records;
	my %index_name;
	my %name_index;
	my $max_sample = 0;
	foreach (@files)
	{
		open (IN, "$_") or die $!;
		while (<IN>)
		{
			next if /^\s+$/;
			next unless /^wsnp/;
			chomp;
			s/\s+$//;
			my @data = (/,/)?split /,/,$_:split /\s+/, $_;
			my ($snp_name, $sample_id, $snp_index, $genotype) = @data[0, 1, 2, 4];
			#print STDERR join("***", ($snp_name,$snp_index)),"\n" if /wsnp_Ex_rep_c67838_66536117/;
			$index_name{$snp_index} = $snp_name if $snp_index =~ /\d/;
			$name_index{$snp_name} = $snp_index if $snp_index =~ /\d/;
			die join('**', (@data)),"\n" unless defined $genotype;
			$snp_index = $name_index{$snp_name};
			$genotype = $data[3] if (length $genotype) > 2;
			push @{$records{$snp_index}}, $_;
		}
	}
	close IN;
	foreach (keys %records)
	{
		my $num_sample = scalar @{$records{$_}};
		$max_sample = $num_sample if $num_sample>$max_sample;
	}
	
	return ($max_sample, %records);
}

sub max
{
	
	my $m = shift;
	#map{$m = $_ if $_>$m} @_;
	foreach (@_){$m = $_ if $_>$m}
	return $m;
}


sub read_delete_list
{
	my $file = shift;
	my %return;
	open (IN, $file) or die $!;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		next if /^\D/;
		s/\s$//;
		print STDERR $_,"*\n";
		$return{$_} = 1;
	}
	close IN;
	return %return;
}