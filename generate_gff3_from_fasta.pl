#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 <fasta_file> <gff3_file>\n";
die $sUsage unless @ARGV == 2;

my ($fasta_file, $gff_file) = @ARGV;

open (IN, "<$fasta_file") or die;
open (OUT, ">$gff_file") or die;

my %seq_length;
my $id;
while(<IN>)
{
	if(/>(\S+)/)
	{
		$id = $1;
		next;
	}
	chomp; 
	$seq_length{$id} += length $_;
}
close IN;

#output gff3
##gff-version 3
#ctg123  .  exon  1300  1500  .  +  .  ID=exon00001
#ctg123  .  exon  1050  1500  .  +  .  ID=exon00002
#ctg123  .  exon  3000  3902  .  +  .  ID=exon00003

foreach my $id (keys %seq_length)
{
	my $gene_id = $id . ".1";
	my $transcript_id = $id . ".1.1";
	my $len = $seq_length{$id};
	my $dot = ".";
	my $str_1 = join("\t", ($id, "rnaseq", "transcript", 1, $len, $dot, "+", $dot, "gene_id "."\"".$gene_id."\";"." transcript_id ". "\"".$transcript_id."\""));
	print OUT $str_1, "\n";
	my $str_2 = $str_1;
	$str_2 =~ s/transcript/exon/;
	print OUT $str_2, "; exon_number \"1\"\n";	 
	
}
close OUT;