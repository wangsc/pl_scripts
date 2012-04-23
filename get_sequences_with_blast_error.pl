#!/usr/bin/perl -w
use strict;

# while running the script run_remote_blastx.pl, some sequences were returned with error which maybe caused by internet
# connection problems.
# This script will compare the fasta file used to do blast and the result file of blast to get those sequences and store 
# them in a fasta file.
# SW 05.04.2011

my $sUsage = qq(
# while running the script run_remote_blastx.pl, some sequences were returned with error which maybe caused by internet
# connection problems.
# This script will compare the fasta file used to do blast and the result file of blast to get those sequences and store 
# them in a fasta file.
# SW 05.04.2011

Usage: perl $0 
<fasta filename>
<blast result filename>
<output fasta filename>
);
die $sUsage unless @ARGV >= 3;

my ($fasta_file, $blast_outfile, $output_file) = @ARGV;

my %fasta_hash = read_fasta_file($fasta_file);
my %blast_hash = read_blast_output($blast_outfile);
output(\%fasta_hash, \%blast_hash, $output_file);

# subroutines;
sub read_fasta_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file\n";
	my $id;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		if (/^>/)
		{
			print 'fasta id ', $_,"\n" unless defined $id;
			$id = $_;			
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

sub read_blast_output
{
	my $file = shift;
	my %return_hash;
	my $flag = 1;
	open (IN, "<$file") or die "can't open file $file \n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		if (/Query\b?\=\b?(.*)/)
		{
			my $hashkey = '>'.$1;
			$hashkey =~ s/\s//g;
			print $hashkey,"\n" if $flag;
			$flag=0;
			$return_hash{$hashkey} = 1;
		}
	}
	close IN;
	return %return_hash;
}

sub output
{
	my ($fasta_hashref, $blast_hashref, $outfile) = @_;
	open (OUT, ">$outfile") or die "can't open file $outfile\n";
	foreach (keys %$fasta_hashref)
	{
		next if exists $blast_hashref->{$_};
		print OUT $_,"\n", $fasta_hashref->{$_},"\n";
	}
	close OUT;
}


