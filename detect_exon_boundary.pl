#!/usr/bin/perl -w
use strict;
use Tie::File::AsHash;
use Search::Dict;

my $sUsage = qq(
# CS genomic sequences were mapped to transcript assemblies
# This script will detect potential exon boundaries in the transcript assemblies.

Usage: 
perl $0
<CS read length file>
<output file>
<minimum similarity, 95>
<minimum leghtn of fragment, 50>
<blat result files>
);
die $sUsage unless @ARGV >= 3;
my ($cs_read_len_file, $out_file, $min_sim, $min_len, @blat_files) = @ARGV;
$min_sim = 95 unless defined $min_sim;
$min_len = 50 unless defined $min_len;

print_time_comment("Reading CS read length file ...");
my %read_len = read_length_file($cs_read_len_file);
#my %read_len;
#tie %read_len, 'Tie::File::AsHash', $cs_read_len_file, split => qr(\t), join=>"\t" or die "Problem tying hash read_len: $!";

print_time_comment("Reading BLAT files ...");
my %exon_boundary = read_blat_files(\@blat_files, $min_sim, $min_len, \%read_len);

print_time_comment("Output results in file $out_file ...");
open (OUT, ">$out_file") or die $!;
foreach my $id(keys %exon_boundary)
{
	print OUT ">", $id,"\n";
	foreach (keys %{$exon_boundary{$id}})
	{
		print OUT join("\t", ($_, $exon_boundary{$id}->{$_})),"\n";
	}
	print OUT "\n";	
}
close OUT;

untie %read_len;

# Subroutines
sub read_length_file
{
	my $file = shift;
	my @return;
	open (IN, $file) or die "can't open file $file \n";
	while(<IN>)
	{
		chomp; 
		next if /^\s+$/;
		my @t=split /\t/,$_;
		push @return, @t;
	}
	close IN;
	return @return;
}



sub read_blat_files
{
	my($files_ref, $min_sim, $min_len, $read_len) = @_;
	my %gene_boundary;
	my $debug = 1;
	foreach my $file (@$files_ref)
	{
		print_time_comment("\tReading $file ...");
		open (IN, "$file") or die $file;
		while(<IN>)
		{
			#	asmbl_2 GJ20S3U01A02OO  97.90   238     2       1       90      324     1       238     6.9e-124        441.0
			# asmbl_2 GJ20S3U01A02OO  99.39   164     0       1       325     488     255     417     8.5e-87 317.0	
			chomp;
			next unless /^\S+/;
			my @data = split /\t/, $_;
			my ($similarity, $length, $num_gap) = @data[2, 3, 5];
			next unless $similarity >= $min_sim and $length >= $min_len and $num_gap == 0;
			my ($gene_id, $read_id, $gene_start, $gene_end, $read_start, $read_end) = @data[0, 1, 6..9];
			my $read_length = $read_len->{$read_id};
			print STDERR "\t$read_id\t", $read_length,"\n" if $debug; $debug=0;
			$gene_boundary{$gene_id}{$gene_start}++ unless $read_start == 1 or $read_start == $read_length;
			$gene_boundary{$gene_id}{$gene_end}++ unless $read_end == 1 or $read_end == $read_length;		
		}
		close IN;		
	}
	return %gene_boundary;
}

sub print_time_comment
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}
