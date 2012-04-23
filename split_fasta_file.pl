#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $sUsage = qq(
perl $0
<original fasta file>
<the number of sub files>
);
die $sUsage unless @ARGV >= 2;
my ($orig_fasta, $num_files) = @ARGV;

my $total_seqs =  get_num_total_seq($orig_fasta);
my $seq_each = int($total_seqs/$num_files);
$seq_each += 1 if ($total_seqs%$num_files)>0;
print '$seq_each: ', $seq_each, "\n";
my @sub_files = map{join('_',($orig_fasta, $_))} 1..$num_files;
my %sub_file_handles = get_file_handles(@sub_files);
open (IN, "$orig_fasta") or die "can't open file $orig_fasta\n";
my $count = 0;
while(<IN>)
{
	next if /^\s+$/;
	$count++ if /^>/;
	my $file_index = int($count/$seq_each) + (($count%$seq_each == 0)?0:1) - 1;
	$file_index  = $num_files - 1 if $file_index >= $num_files;
	#print $count,"\t", $file_index,"\n";
	print {$sub_file_handles{$file_index}} $_;
}
foreach (values %sub_file_handles)
{
	close $_ or print STDERR  "Can't close File handle $_\n";
}

sub get_num_total_seq
{
	my $file = shift;
	my $total;
	open (IN, "$file") or die $!;
	while (<IN>)
	{
		$total++ if /^>/;
	}
	close IN;
	return $total;
}

sub get_file_handles
{
	my @files = @_;
	my %return_hash;
	foreach (0..$#files)
	{
		local *FH;
		open (FH, ">$files[$_]") or die "can't open file $files[$_]\n";
		$return_hash{$_} = *FH{IO};
	}
	return %return_hash;
}