#!/usr/bin/perl -w
use strict;

# calcualte coverage for target regions of sequence capture data
#

my $susage = qq(
perl $0
isoform GFF3 file for genes in chr1 and chr2 prediction by PASA
positional coverage file, generated by bedtools
);

die $susage unless @ARGV;

my ($gff_file, $pos_cov_file) = @ARGV;

my %isoform_vector = read_gff_file($gff_file);

my ($on_target_length, $off_target_length, $on_target_read_length, $off_target_read_length) = (0, 0, 0, 0);

open (IN, $pos_cov_file) or die;
while (<IN>)
{
	chomp;
	#1AL 23 4
	#1AS 24 5
	#1AS 25 1
	# ... 
	my @data = split /\s+/, $_;
	
	if (vec($isoform_vector{$data[0]}, $data[1], 1) == 1)
	{
		$on_target_length++;
		$on_target_read_length += $data[2];
		print $data[2], "\n";
	}
	else
	{
		$off_target_length++;
		$off_target_read_length += $data[2];
	}
	
}
close IN;

print STDERR  "Target sequences (bp): ", $on_target_length, "\n";
print STDERR "Target read sequences (bp): ", $on_target_read_length, "\n";
print STDERR "Off-target sequences (bp): ",  $off_target_length, "\n";
print STDERR "Off-target read sequences (bp): ", $off_target_read_length, "\n";


# Subroutines

sub read_gff_file
{
	my $file = shift or die;
	my %return_vec;
	open(IN, $file) or die;
	while(<IN>)
	{
		chomp; 
		next unless /\S/;
		# 1AL     PASA    cDNA_match      4385781 4386878 .       -       .       ID=S9-asmbl_9; Target=asmbl_9 1 1098 +
		my @t=split /\s+/, $_; 
		$return_vec{$t[0]} = '' unless exists $return_vec{$t[0]};
		map{vec($return_vec{$t[0]}, $_, 1) = 1}$t[3]..$t[4];
	}
	close IN;
	
	return %return_vec;	
}
