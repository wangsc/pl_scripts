#!/usr/bin/perl -w
use strict;
# concatenate contigs into one psudeo-chromosome
# 

my $sUsage = qq(
perl $0
<3AL fasta file>
<3AS fasta file>
<3BL fasta file>
<3BS fasta file>
<3DL fasta file>
<3DS fasta file>
);
die $sUsage unless @ARGV >= 2;
my @files = @ARGV;

my %chrs = map {($1, 1) if /^(\d)/} @files;
if(keys %chrs > 1)
{
	print STDERR 'Error: You provides more than one chromosomes: ', join("\t", keys %chrs), "\n";
	print $sUsage;
	exit;
}

my $output_file = "chr_" . (keys %chrs)[0] . "_concatenateX50.fasta";
print STDERR "!!\tThe output file is : ", $output_file, "\n";
open (OUT, ">$output_file") or die " can't open file $output_file\n";

my $linker_seq = "X" x 50;

foreach my $f (@files)
{
	my $id = $1 if $f=~/^(\S+?)\-/;
	print OUT ">", $id, "\n";
	open (IN, $f) or die ;
	my $n = 1;
	while (<IN>)
	{
		next if /^>/;
		chomp;
		s/\s//g;
		s/\d//g;
		if($n)
		{
			print OUT $_;
			$n = 0;
		}
		else
		{
			print OUT $linker_seq, $_
		}
	}
	print OUT "\n";
	close IN;
}
close OUT;
