#!/usr/bin/perl -w
use strict;
use Bio::DB::Sam;
use File::Basename;

my $sUsage = qq(
perl $0
<bam file>
<SNP file, combined output of "vcfutils.pl varFilter" from each accession>
<reference fasta file>
<output file>
);
die $sUsage unless @ARGV >= 4;

my($bamfile, $snpfile, $ref_fasta_file, $outfile) = @ARGV;
open (OUT, ">$outfile") or die "$!\n";
print_time_comment("Reading SNP file $snpfile ...");
my %snps = read_snp_file($snpfile);
print_time_comment("Reading file $bamfile ...");
my $bam_index = $bamfile . ".bai";
my $autoindex = -e $bam_index?1:0;
my $sam = Bio::DB::Sam->new(-bam  => $bamfile,
                            -fasta=> $ref_fasta_file,
                            -autoindex => $autoindex,
                            );
my $alignments = $sam->features(-iterator=>1);

my %allele_count;
my %ref_allele;
#my %covered_snps;
my $counter = 0;
while(my $aln = $alignments->next_seq())
{
	$counter++;
	print_time_comment("\t$counter alignments ...") unless $counter%1000000;
	my $id = $aln->seq_id();
	next unless exists $snps{$id};
	my $ref_start = $aln->start;
	next unless defined $ref_start;
	my $ref_end = $aln->end;
	my @snps = grep {$_>=$ref_start and $_<=$ref_end} @{$snps{$id}};
	#push @{$covered_snps{$id}}, @snps;
	my $ref_dna   = $aln->dna;
	my $query_dna = $aln->query->dna;
	my @qscores = $aln->qscore;
	foreach my $snppos (@snps)
	{
		next unless $qscores[$snppos-$ref_start] >= 20;
		$ref_allele{$id}{$snppos} = substr($ref_dna, $snppos-$ref_start, 1) unless exists $ref_allele{$id}{$snppos};
		$allele_count{join("\t",($id,$snppos))}{substr($query_dna, $snppos-$ref_start, 1)}++;
	}
}

# Output
print_time_comment("Output results to file $outfile...");
foreach my $id (keys %ref_allele)
{
	#my @covered_snps = @{$covered_snps{$id}};
	my @snps = @{$snps{$id}};
	foreach my $snppos (unique(@snps))
	{
		next unless exists $ref_allele{$id}{$snppos};
		print OUT $id, "\t", $snppos, "\t";
		print OUT $ref_allele{$id}{$snppos},"_", 
		     exists $allele_count{join("\t",($id,$snppos))}{$ref_allele{$id}{$snppos}}?$allele_count{join("\t",($id,$snppos))}{$ref_allele{$id}{$snppos}}:0;
		foreach (keys %{$allele_count{join("\t",($id,$snppos))}})
		{
			unless($_ eq  $ref_allele{$id}{$snppos})
			{
				print OUT "\t", $_,"_",exists $allele_count{join("\t",($id,$snppos))}{$_}?$allele_count{join("\t",($id,$snppos))}{$_}:0;
			}
		}
		print OUT "\n";
	}	
}

close OUT;

# subroutines

sub min_max
{
	my $min = shift;
	my $max = $min;
	foreach (@_)
	{
		$min = $_ if $_ < $min;
		$max = $_ if $_ > $max;
	}
	return ($min, $max);
}

sub unique
{
	my %h = map{$_, 1}@_;
	return keys %h;	
}

sub read_snp_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		next if /^\#/;
		next if /^\s+$/;
		chomp;
		my @t = split /\s+/,$_;
		push @{$return{$t[0]}}, $t[1];
	}
	close IN;
	map{ $return{$_} = [unique(@{$return{$_}})] }keys %return;
	return %return;
}

sub print_time_comment
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}
