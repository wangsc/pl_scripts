#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage = qq(
perl $0 
<est_transcript_combined.fasta>
<bac_specific_est_transcript>
<BAC sequences in fasta>
<output 1: cDNA for genes>
<output 2: gene structure in gff>
);
die $sUsage unless @ARGV >= 5;
my ($est_ctg_fasta, $bac_specific_file, $bac_fasta, $cdna_out, $gff_out) = @ARGV;
unlink($gff_out) if -e $gff_out;
my $est_obj = Bio::DB::Fasta->new($est_ctg_fasta);
my $bac_obj = Bio::DB::Fasta->new($bac_fasta);
my %bac_specific = read_bac_specific_file($bac_specific_file);

foreach my $bac_id (keys %bac_specific)
{
	my @ests = @{$bac_specific{$bac_id}};
	next unless @ests;
	generate_fasta("est.fasta", $est_obj, @ests);
	generate_fasta("bac.fasta", $bac_obj, $bac_id);
	run_geneseqer();
	parse_geneseqer_output($bac_id, $gff_out);	
}
output_cds($bac_obj, $gff_out, $cdna_out);

# Subroutines;
sub read_bac_specific_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		chomp;
		my @t = split /\s+/,$_;
		my $id = shift @t; 
		$return{$id} = [@t];
	}
	close IN;
	return %return;
}

sub generate_fasta
{
	my ($file, $gnobj, @ids) = @_;
	open (OUT, ">$file") or die;
	foreach (@ids)
	{
		my $seq = $gnobj->get_Seq_by_id($_)->seq;
		print OUT ">$_\n", $seq, "\n" if defined $seq;
	}
	close OUT;
	return 1;
}

sub run_geneseqer
{
	my $geneseqer_cmd = "/usr/local/bin/GeneSeqer";
	my $params = " -s rice -E est.fasta -L bac.fasta >geneseqer.out";
	warn "!! Geneseqer failed \n" if system($geneseqer_cmd . $params);
}

sub parse_geneseqer_output
{
	my $id = shift;
	my $gff = shift;
	open (OUT, ">>$gff") or die;
	open (IN, "geneseqer.out") or die;
	my $flag = 0;
	while(<IN>)
	{
		if(/Predicted gene structure/)
		{
			$flag++;
			next;
		}
	#	#  Exon  1  35628  35712 (  85 n);  cDNA      1     85 (  85 n); score: 1.000
		if(/^\s+Exon\s+\d+\s+(\d+)\s+(\d+).*cDNA\s+(\d+)\s+(\d+)/)
		{
			print STDERR $_;
			my $gene_id = $id . "_" . $flag;
			print OUT $id, "\tgeneseqer\texon\t$1\t$2\t.\t", $1<$2?'+':'-',"\t.\tID=$gene_id; Target=$gene_id $3 $4", "\n";
		}
	}
	close IN;
	close OUT;	
}

# gi|299109310|emb|FN564426.1|	PASA	cDNA_match	97191	98148	.	-	.	ID=S5-asmbl_6; Target=asmbl_6 2147 3104 +
# gi|299109310|emb|FN564426.1|	PASA	cDNA_match	98243	98596	.	-	.	ID=S5-asmbl_6; Target=asmbl_6 1793 2146 +


sub output_cds
{
	my ($gnobj, $gff, $output) = @_;
	my %cdna;
	open (G, "$gff") or die;
	while (<G>)
	{
		next if /^\s+$/;
		my $id = $1 if /ID=(\S+)\;/;
		my @t = split /\s+/, $_;
		push @{$cdna{$id}}, [@t[3, 4]];
	}
	close G;
	
	open (OUT, ">$output") or die;
	foreach my $gene (keys %cdna)
	{
		my $bac = $1 if $gene=~/^(\S+)_\d+$/;
		my $seq = "";
		foreach (@{$cdna{$gene}})
		{
			my ($start, $end) = @$_;
			$seq .= $gnobj->seq($bac, $start => $end);
		}
		print OUT ">$gene\n", $seq, "\n";
	}
	close OUT;
}















