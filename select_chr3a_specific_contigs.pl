#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $sUsage = qq(
The 454 RNA-seq contigs from all 9 accessions were mapped to chr3A genomic sequences using BLAT.
This script will select the contigs based on the similarity and mapped length in the BLAT output file.
Usage:
perl $0 
-in blat_output_file
-ctg RNA-seq contig fasta file
-ms 98 [minimum similarity (%), default is 98]
-mg 30 [minimum length of mapped fragments, smaller frag will not count, default 30]
-ml 300 [minimum total length of fragments of contig that mapped, default is 300]
-out output_file
);
my $ms = 98;
my $mg = 30;
my $ml = 300;
my ($input, $output, $contig);
GetOptions("in=s" => \$input, 
					 "ctg=s" => \$contig,
					 "ms=i" => \$ms,
					 "mg=i" => \$mg,
					 "ml=i" => \$ml,
					 "out=s" => \$output
);
die $sUsage unless defined $input and defined $output;

my %blat_results = read_blat_result($input, $ms, $mg);
select_contigs_and_output(\%blat_results, $ml, $output, $contig);

#
sub read_blat_result
{
	my ($file, $ms, $mg) = @_;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file\n";
	while (<IN>)
	{
		chomp;
		next if /^\s+$/;
		my @data = split /\s+/, $_;
		# BobWhite_mira1_c51	contig17660	100.00	117	0	0	943	1059	501	385	4.6e-60	229.0
		my ($id, $similarity, $length, $start, $end) = @data[0, 2, 3, 6, 7];
		$id = $1 if /(gi\|.*?)\|/;
		next unless $similarity>$ms and $length>$mg;
		($start, $end) = ($end, $start) if $start > $end;
		push @{$return_hash{$id}}, [$start, $end];
	}
	close IN;
	return %return_hash;
}

sub select_contigs_and_output
{
	my($blat_ref, $ml, $output, $contig) = @_;
	open (OUT, ">$output") or die "can't open file $output\n";
	my %selected_contigs = select_contigs($blat_ref, $ml);
	open (CTG, "$contig") or die "can't open file $contig\n";
	my $flag = 0;
	while (my $line = <CTG>)
	{
		next if $line =~ /^\s+$/;
		chomp $line;
		if ($line =~ /^>(.*)/)
		{
			
			#$flag = exists $selected_contigs{$1}?1:0;
			my $id = $1;
			$id = $1 if $line =~ /(gi\|.*?)\|/;
			if (exists $selected_contigs{$1})
			{
				$flag = 1;
			}

		}
	print OUT $line, "\n" if $flag == 1;
	}
	close CTG;
	close OUT;	
}

sub select_contigs
{
	my ($blat_ref, $min_len) = @_;
	my %return_hash;
	foreach my $id (keys %{$blat_ref})
	{
		my ($min, $max);
		my $vector = '';
		foreach (@{$blat_ref->{$id}})
		{
			my ($frag_start, $frag_end) = @$_;
			$min = $frag_start unless defined $min;
			$max = $frag_end unless defined $max;
			$min = $frag_start if $frag_start<$min;
			$max = $frag_end if $frag_end>$max;
			foreach ($frag_start .. $frag_end)
			{
				vec($vector, $_, 1) = 0b1;
			}
		}
		my $length = 0;
		foreach( $min.. $max)
		{
			$length++ if vec($vector, $_, 1) == 1;
		}
		$return_hash{$id} = 1 if $length >= $min_len;
	}
	return %return_hash;
}










