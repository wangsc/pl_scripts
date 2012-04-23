#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::DB::Sam;

my $sUsage = qq(
perl $0
<gbs tag id file>
<CS ref sequence in fasta>
<bam files>
<output file>
);
die $sUsage unless @ARGV >= 3;

my ($tagid_file, $refseq_file, @bam_files) = @ARGV;
my $outfile = pop @bam_files;
if (-e $outfile){print "Are you sure $outfile will be the output file?\n"; exit}
open (OUT, ">$outfile") or die "can't open file $outfile\n";
my $genome_obj = Bio::DB::Fasta->new($refseq_file);

my %tag_ids = read_tagid_file($tagid_file);

my %recorder;
foreach my $bamfile (@bam_files)
{
	print_time_comment("processing file $bamfile ...");
	my $bam_index = $bamfile . ".bai";
	my $autoindex = -e $bam_index?1:0;
	my $sam = Bio::DB::Sam->new(-bam  => $bamfile,
	                            -fasta=> $refseq_file,
	                            -autoindex => $autoindex,
	                            );
	my $alignments = $sam->features(-iterator=>1);
	
	while(my $aln = $alignments->next_seq())
	{
		my $ref_seq = $aln->dna;
		my $query_seq = $aln->query->dna;
		my $ori_len = length $ref_seq;
		my $strand = $aln->strand;
		my $seq_id = $aln->seq_id;
		my $start = $aln->start;
		my $end = $aln->end;
		next unless defined $start and defined $end;
		my @snp_pos = check_snp_pos($seq_id, $start, $end, \%tag_ids);
		next unless @snp_pos > 0;
		# concatenate 10bp before segment
		
		my $ten_bp_seq = $genome_obj->seq($seq_id, $start-10 => $start-1);
		$ten_bp_seq = $genome_obj->seq($seq_id, $end+1 => $end+10) if $strand < 0;
		
		if($strand < 0){$ref_seq .= $ten_bp_seq; $ref_seq = rev_comp($ref_seq)} 
		else{$ref_seq = $ten_bp_seq . $ref_seq;}
		
		my $rs_pos = 200;
		while($ref_seq =~ /TGCAG/g)
		{
			$rs_pos = pos($ref_seq);
			last;
		}
		next unless $rs_pos < 16;
		substr($ref_seq, 0, $rs_pos-5)="";
		next if exists $recorder{$ref_seq};
		my $increased_len = (length $ref_seq) - $ori_len;
		my $new_start = $strand<0?-1*$start:$start-$increased_len;
		print OUT join("\t",($ref_seq, @snp_pos, $new_start)),"\n";
		$recorder{$ref_seq} = 1;
	}
}

sub rev_comp
{
	my $seq = shift;
	$seq = reverse $seq;
	$seq =~ tr/[ATGC]/[TACG]/;
	return $seq;
}

sub check_snp_pos
{
	my ($id, $s, $e, $taghash) = @_;
	my @return;
	foreach ($s..$e)
	{
		my $new_id = $id . ":". $_;
		push @return, $new_id if exists $taghash->{$new_id}
	}
	return @return;
}

sub read_tagid_file
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while(<IN>)
	{
		chomp; 
		$return{$_}=1
	}
	close IN;
	return %return;
}


sub print_time_comment
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}
