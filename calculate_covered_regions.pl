#!/usr/bin/perl -w
use strict;
# GC2 transcripts were queried against CS 5X data by BLAT,
# this script will check the covered length of each GC2 transcript.
# SW 08.11.2011


my $sUsage = "perl $0 <fasta file> <min similarity> <min length of fragment> <BLAT file 1> <BLAT file 2> ...\n";
die $sUsage unless @ARGV >= 4;
my ($fasta_file, $min_sim, $min_len, @blat_files) = @ARGV;
die $sUsage if $min_len =~ /\D/;

my %ctg_length = read_fasta_file($fasta_file);
my %covered_length = read_blat_files(\%ctg_length, $min_sim, $min_len, @blat_files);
foreach (keys %ctg_length)
{
	my $covered = exists $covered_length{$_}?$covered_length{$_}:0;
	print $_,"\t", $covered, "\t", $ctg_length{$_},"\t", $covered/$ctg_length{$_}, "\n";
}

# Subroutines

sub read_fasta_file
{
	my $file = shift;
	open (IN, $file) or die;
	my $id;
	my %return;
	my $debug=1;
	while(<IN>)
	{
		next if /^\s+$/;
		if(/^>(\S+)/)
		{
			$id=$1;
			#print STDERR '$id in fasta: ', $id,"\n" if $debug; $debug=0;
			next;
		}
		chomp;
		$return{$id} += length $_;
	}
	close IN;
	return %return;
}

sub read_blat_files
{
	my $ctg_length_ref = shift;
	my $min_sim = shift;
	my $min_len = shift;
	my @files = @_;
	my %covered_regions;
	my $debug=1;
	foreach my $file (@files)
	{
		open (IN, $file) or die;
		while (<IN>)
		{
	    # BobWhite_mira1_c1       3934886 74.84   465     1     0       230     694     4287    3823    1.9e-149        526.0
			next if /^\s+$/;
			my @t = split /\s+/,$_;
			next unless $t[2]>= $min_sim and $t[3]>= $min_len;
			print STDERR '$id in blat: ', $t[1],"\n" if $debug; $debug=0;
			#push @{$covered_regions{$t[0]}}, [@t[6, 7]];
			construct_vec(\%covered_regions, $t[0], sort{$a<=>$b}@t[6, 7]);
		}
		close IN;		
	}
	return calculate_covered_length(\%covered_regions, $ctg_length_ref);
}

sub calculate_covered_length
{
	my ($covered_ref, $length_ref) = @_;
	my %return;
	foreach my $id (keys %$covered_ref)
	{
		my $vec = $covered_ref -> {$id};
		print STDERR "No length for $id \n" unless exists $length_ref->{$id};
		next unless exists $length_ref->{$id};
		my $len = 0;
		#print STDERR "length of contig $id :", $length_ref->{$id}, "\n";
		foreach (1..$length_ref->{$id})
		{
			$len ++ if (vec($vec, $_,1) == 1);
		}
		$return{$id} = $len;
	}
	return %return;
}

sub construct_vec
{
	my ($hashref, $id, $start, $end) = @_;
	my $vec = exists $hashref->{$id}?$hashref->{$id}:"";
	my $debug = 1;
	#print STDERR join("\t", ($start, $end)), "\n" if $debug; $debug=0;
	foreach ($start .. $end)
	{
		vec($vec, $_, 1) = 0b1;
	}
	$hashref->{$id} = $vec;
}









