#!/usr/bin/perl -w 
use strict;
use Getopt::Long;
use List::Util qw(min max sum);

my $usage = <<'USAGE';

################ CNV-seq ################

	usage: cnv.seq.pl [options]

		--test = test.hits.file
		--ref = ref.hits.file
    --gene = gene_position_in_ref.file
		--log2-threshold = number
		  (default=0.6)
		--p-value = number
		  (default=0.001)

		--bigger-window = number
		  (default=1.5)
		--window-size
		  (default: determined by log2 and p-value;
		   if set, then log2 and p-value are not used)

		--genome = (human, chicken, chrom1, autosome, sex, chromX, chromY)
		  (higher priority than --genome-size)
		--genome-size = number
		  (no use when --genome is avaliable)

		--annotate
		--no-annotate
		  (default: do annotation, ie --annotate)
		--minimum-windows-required = number 
		  (default=4; only for annotation of CNV)
		--global-normalization
		  (if used, normalization on whole genome,
		  instead of oneach chromosome)

		--Rexe = path to your R program
		  (default: R)

		--help
		--debug


#########################################

USAGE

my ($test_hit, $ref_hit, $genome_size);
my $log2 = 0.6;
my $pvalue = 0.001;
my $genome;
my $min_window=4;
my $chromosomal_normalization = 'TRUE';
my $annotation='TRUE';
my $window_size;
my $bigger = 1.5;
my $debug = 0;
my $Rexe = 'R';
my $gene_pos_file = "";

my $result = GetOptions(
	"test=s" => \$test_hit,
	"ref=s" => \$ref_hit,
	"gene=s" => \$gene_pos_file,
	"log2-threshold=f" => \$log2,
	"p-value=f" => \$pvalue,
	"bigger-window=f" => \$bigger,
	'window-size=i' => \$window_size,
	"genome=s" => \$genome,
	"genome-size=o" => \$genome_size,
	"global-normalization" => sub{$chromosomal_normalization='FALSE'},
	"no-annotate" => sub{$annotation='FALSE'},
	"annotate" => sub{$annotation='TRUE'},
	"Rexe=s" => \$Rexe,
	"debug" => sub{$debug='TRUE'},
	"minimum-windows-required=i" => \$min_window,
	"help|?" => sub{print $usage; exit}
);

die $usage unless($test_hit && $ref_hit);
open(TEST,$test_hit) or die "can not find $test_hit file\n";
open(REF, $ref_hit) or die "can not find $ref_hit file\n";
die "must set --genome or --genome-size or --window-size\n" unless($genome_size||$genome||$window_size);

my $out;
my $temp = $test_hit;
$temp =~ s/.+\///;
$out = $temp;
$temp = $ref_hit;
$temp =~ s/.+\///;
$out .= "-vs-$temp";
my ($total_test, $total_ref);
$out .= ".log2-$log2.pvalue-$pvalue";
my $cnvout = $out;
$cnvout .= ".minw-$min_window" if $annotation eq 'TRUE';

print "start counting test hits ... \n" if $debug;
my $size = 100000; # for speed up
my @gene_positions = read_gene_position_file($gene_pos_file,$size);
my (%count, %test_total, %ref_total);
my $pat = qr/^\s*(\S+)\s+([\d.eE\+-]+)\s*$/;
my %chrom;

my $n=0;
open(TEST,$test_hit) or die;
open(REF, $ref_hit) or die;
while(<TEST>)
{
	if(m/$pat/)
	{
		$n++;
		$chrom{$1}=1;
		my $aln_pos = $2+0;
		my $index = int($aln_pos/$size);
		my $window_pos;
		next unless defined $gene_positions[$index];
		my @positions = @{$gene_positions[$index]};
		foreach (0..$#positions)
		{
			if($aln_pos>=$positions[$_]->[0] and $aln_pos<$positions[$_]->[1])
			{
				$window_pos = $_;last
			}
		}
		next unless (defined $window_pos) ;#{print $_, $aln_pos, "\n"; exit;}
		$count{$1}->[0][$index*$size + $window_pos]++;
	}
}
print "read $n test reads, out of $total_test lines\n";
print "start counting ref hits ... \n" if $debug;
$n=0;
while(<REF>)
{
	if(m/$pat/)
	{
		$n++;
		$chrom{$1}=1;
		my $aln_pos = $2+0;
		my $index = int($aln_pos/$size);
		my $window_pos;
		next unless defined $gene_positions[$index];
		my @positions = @{$gene_positions[$index]};
		foreach (0..$#positions)
		{
			if($aln_pos>=$positions[$_]->[0] and $aln_pos<$positions[$_]->[1])
			{
				$window_pos = $_;last;
			}
		}
		next unless (defined $window_pos);
		$count{$1}->[1][$index*$size + $window_pos]++;
	}
}
print "read $n ref reads, out of $total_ref lines\n";

if($debug)
{
	print "done counting:\n";
	use Data::Dumper;
	my $chrom = Dumper(%chrom);
	my $count = Dumper(%count);
	print '%chrom: ', substr($chrom, 0, 100), " ...\n";
	print '%count: ', substr($count, 0, 100), " ...\n";
}

open(OUT, ">$cnvout.count") or die;
print "write read-counts into file: $out.count\n";
print OUT "chromosome\tstart\tend\ttest\tref\n";
foreach my $chr(sort{$a<=>$b} keys %chrom)
{
	my @test = @{$count{$chr}->[0]};
	my @ref = @{$count{$chr}->[1]};
	for my $id (0..min($#test,$#ref))
	{
		next unless defined $test[$id] or defined $ref[$id];
		my ($start, $end) = @{$gene_positions[int(($id)/$size)][($id)%$size]};
		$test[$id] = 0 unless defined $test[$id];
		my $test_raw = $test[$id]+0;
		$ref[$id] = 0 unless defined  $ref[$id];
		my $ref_raw = $ref[$id]+0;
		print OUT "$chr\t$start\t$end\t$test_raw\t$ref_raw\n";
	}
}
close OUT;

print "$Rexe package cnv output: $cnvout.cnv\n";
system(qq`echo 'library(cnv);\ncnv<-cnv.cal("$cnvout.count", log2=$log2, min=$min_window, chromosomal=$chromosomal_normalization, annota=$annotation);\nwrite.table(cnv,"$cnvout.cnv",row.names=FALSE,sep="\\t");\n' | $Rexe --vanilla --slave`);


sub read_gene_position_file
{
	my $file = shift;
	my $size = shift;
	open (IN, $file) or die;
	my @return;
	while (<IN>)
	{
		chomp; 
		my @t=split/\t/,$_;
		my $win = int($t[2]/$size);
		push @{$return[$win]}, [@t[1,2]];
	}
	close IN;
	return @return;
}



