#!/usr/bin/perl -w
use strict;
my $sUsage = "perl $0 <ortholog list> <file_from_hanquan> \n";
die $sUsage unless @ARGV >= 2;
my ($ortho_list, $multi_hit_file) = @ARGV;

my ($record_array, $r_w_genes, $w_r_genes) = read_ortho_list($ortho_list);
my ($rice, $brachy) = read_multi_hit_file($multi_hit_file);

foreach (@$record_array)
{
	my $name;
	if (/\|/)
	{
		$name = $1 if /(LOC.*?)\.\d+\|/;
	}
	else
	{
		$name = $1 if /(Bradi.*?)\.\d/;	
	}
	next if $rice->{$name};
	next if scalar @{$r_w_genes->{$name}} >1;
	chomp;
	my @t= split /\t/, $_;
	if (scalar @{$w_r_genes->{$t[-1]}} >1){print STDERR $_,"\n";next }
	print $_,"\n";	
}


sub read_multi_hit_file
{
	my $file = shift;
	my %rice;
	my %brachy;
	open (IN, $file) or die $!;
	while(<IN>)
	{
		chomp;
		next if /^\s+$/;
		next unless /^wheat/;
		my @data = split /\t/, $_;
		if (scalar @data > 7)
		{
			$rice{$data[1]} = 1;
			$brachy{$data[2]} = 1;
		}
	}
	close IN;
	return (\%rice, \%brachy)
}


sub read_ortho_list
{
	my $file = shift;
	my @records;
	my %rice_wheat_genes;
	my %wheat_rice_genes;
	open (IN, $file) or die $!;
	while (<IN>)
	{
		next if /^\s+/;
		push @records, $_;
		chomp;		
		my @data = split /\t/, $_;
		my $name;
		if (/\|/)
		{
			$name = $1 if /(LOC.*?)\.\d+\|/;			
		}
		else
		{
			$name = $1 if /(Bradi.*?)\.\d/;	
		}
		push @{$rice_wheat_genes{$name}}, $data[-1];
		push @{$wheat_rice_genes{$data[-1]}}, $name;
	}
	close IN;
	map{$rice_wheat_genes{$_} = [unique(@{$rice_wheat_genes{$_}})]} keys %rice_wheat_genes;
	map{$wheat_rice_genes{$_} = [unique(@{$wheat_rice_genes{$_}})]} keys %wheat_rice_genes;
	return(\@records, \%rice_wheat_genes,\%wheat_rice_genes);
}

sub unique
{
	my %h=map{$_, 1} @_;
	return (keys %h);
}