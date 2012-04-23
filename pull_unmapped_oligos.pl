#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <oligo seq in fasta> <BLAT output file>\n";
die $sUsage unless @ARGV >= 2;
my ($oligo_file, $blat_file) = @ARGV;

my %oligo = read_fasta_file($oligo_file);
my %mapped_oligo = read_blat_file($blat_file);
foreach (keys %oligo)
{
	next if exists $mapped_oligo{$_};
	print '>',$_,"\n", $oligo{$_}, "\n";
}


sub read_blat_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my @temp = split /\t/, $_;
		$return_hash{shift @temp} = 1;		
	}
	close IN;
	return %return_hash;
}



sub read_fasta_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file\n";
	my $id;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		if (/^>/)
		{
			$id = $_;
			$id =~ s/^>//;
			$return_hash{$id} = '';
			next;
		}
		else
		{
			$return_hash{$id} .= $_;
		}
	}
	close IN;
	return %return_hash; 
}