#!/usr/bin/perl -w
# SW
use strict;
use Bio::DB::Fasta;

die qq(
perl $0 
<SNP position file; format:contig_id SNP_pos> 
<contig fasta> 
<flanking length: 20, 50, whatever you like> 
<oligo output file>) 
unless @ARGV >= 3;

my ($snp_file, $ctg_fasta_file, $flank_len, $oligo_file) = @ARGV;
open (OL, ">$oligo_file") or die $!;
my $ctg_fasta = Bio::DB::Fasta->new($ctg_fasta_file);
my %snp_pos = read_snp_file($snp_file);

foreach my $id (keys %snp_pos)
{
	my $obj = $ctg_fasta->get_Seq_by_id($id);
	my $length = $obj->length;
	foreach my $pos (@{$snp_pos{$id}})
	{
		my $start = $pos>$flank_len?$pos-$flank_len:0;		
		my $end = ($length - $pos)>=$flank_len?$pos+$flank_len:$length;		
		my $subseq = $obj->subseq($start => $end);
		print STDERR  join("\t",($start, $end)),"\n" if $id eq "tplb0059a11";
		print ">", join(":", ($id, $pos, $length)),"\n", $subseq,"\n";
		
		if($pos >= $start+$flank_len)
		{
			my $left_seq = $obj->subseq($start=>$pos-1);
			print OL ">", join(":", ($id, $pos, $length)),"_L\n", $left_seq,"\n";
		}
		if($end >= $pos + $flank_len)
		{
			my $right_seq = $obj->subseq($pos+1 => $end);
			print OL ">", join(":", ($id, $pos, $length)),"_R\n", $right_seq,"\n";
		}		
	}
}

sub read_snp_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die $!;
	while(<IN>)
	{
		chomp; 
		next if /^\s+$/;
		s/^>//;
		my @t=split /\s+/,$_;
		push @{$return{$t[0]}}, $t[1];
	}
	close IN;
	return %return;
}