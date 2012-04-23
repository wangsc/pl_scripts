#!/usr/bin/perl -w
use strict;
use Bio::DB::Sam;

# Just for test!

my $bam_file = shift;
my $ref_file = shift;

my $sam = Bio::DB::Sam->new(-bam  => $bam_file,
                            -fasta=> $ref_file,
                            );

my @targets = $sam->seq_ids;
#print join("\n", @targets),"\n";
my $id = "BobWhite_mira1_rep_c51791";
my $snp_pos = 1234;
#my @alignments = $sam->get_features_by_location(-seq_id => $id,
#                                                -start  => $snp_pos,
#                                                -end    => $snp_pos);
my $alignments = $sam->features(-iterator=>1);
while( my $a = $alignments->next_seq())
{
	# where does the alignment start in the query sequence
	print 'seq_id: ', $a->seq_id,"\n";
	my $query_start = $a->query->start;
	my $query_end   = $a->query->end;	
	my $ref_dna   = $a->dna;        # reference sequence bases
	my $ref_base = substr($ref_dna, $snp_pos-$a->start,1);
	my $query_dna = $a->query->dna; # query sequence bases
	my $alt_base = substr($query_dna, $snp_pos-$a->start,1);
	my @scores = $a->qscore;
	my @query_scores = $a->query->qscore;
	print $ref_dna,"\t", $a->start,"\t", $ref_base, "\t", $scores[$snp_pos-$a->start], "\n", 
	      $query_dna,"\t", $query_start, "\t",$alt_base,"\t", $query_scores[$snp_pos-$a->start], "\n", 
	      "*"x20,"\n";
}

my $segment = $sam->segment(-seq_id=>$id,-start=>$snp_pos,-end=>$snp_pos);
print $segment->dna,"\n", $segment->start,"\n", $segment->end,"\n";
my @all_alignments = $segment->features;
 # get an iterator across the alignments
my $iterator     = $segment->features(-iterator=>1);
while (my $align = $iterator->next_seq) 
{ 
	#print $align->dna,"\t", $align->query->dna,"\n";
}
