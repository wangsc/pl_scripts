#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<barcode file>
<fasta files, name must be accession_id.fasta>
);

die $sUsage unless @ARGV >= 2;
my ($barcode_file, @fasta_files) = @ARGV;
my %barcode = read_barcode_file($barcode_file);

foreach my $file (@fasta_files)
{
	print STDERR "Processing file $file ...\n";
	my $acc_id = $1 if $file =~ /(\S+)\.fasta/;
	my $output = $file . "_nobarcode";
	
	open (IN, $file) or die;
	open (OUT, ">$output") or die;
	while (<IN>)
	{
		if(/>/)
		{
			print OUT $_;
			next;
		}
		chomp;
		next if /^\s+$/;
		my $seq = $_;
		my @barcodes = @{$barcode{$acc_id}};
		foreach my $code (@barcodes)
		{
			my $count = $seq=~s/^$code//;
			last if $count;
		}
		print OUT $seq, "\n";
	}
	close IN;
	close OUT;
}

sub read_barcode_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while (<IN>)
	{
		chomp; 
		my @t = split /\s+/, $_;
		push @{$return{$t[0]}}, uc($t[1]);
	}
	close IN;
	return %return;
}


