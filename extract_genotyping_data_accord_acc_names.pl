#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <genotype_file> <accession names>\n";
die $sUsage unless @ARGV >= 2;
my ($genotype_file, $acc_file) = @ARGV;

my @acc_names = read_acc_file($acc_file);
print STDERR scalar @acc_names, "\n";
my %accname_hash = map{uc($_), 1} @acc_names;
open (IN, $genotype_file) or die;
my @acc_index;
my %values = ("AA", 1, "AB", 0, "BB", -1, "--", "NA");
while(<IN>)
{
	next unless /^\d+|Index/;
	chomp;
	s/\s+//g;	
	my @t = split /,/,$_;
	if(/^Index/)
	{	
		#print $_, "\n"; exit;
		my %h = map{uc($t[$_]), $_} (0..$#t);
		foreach my $name (@acc_names)
		{
			my $acc = uc($name);
			my $index = exists $h{$acc}?$h{$acc}:10000;
			push @acc_index, $index;
		}
		print STDERR scalar @acc_index, "\n";
		next;		
	}
	
	my @output;
	foreach (@acc_index)
	{
		
		if(defined $t[$_] and exists $values{$t[$_]})
		{
			#unless (exists $values{$t[$_]}){print STDERR $t[$_],"\n"; exit}
			push @output, $values{$t[$_]};
		}
		else
		{
			push @output, "NA";
		}		
	}
	print join("\t", ($t[1], @output) ),"\n";
}
close IN;


sub read_acc_file
{
	my $file = shift;
	my @return;
	open (IN, $file) or die;
	while(<IN>)
	{
		chomp;
		@return = split/\t/,$_;  
	}
	close IN;
	return @return[1..$#return];
}



