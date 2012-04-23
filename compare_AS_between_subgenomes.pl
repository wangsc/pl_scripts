#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage = <<END;
To compare the AS between two subgenomes:
1. Extract genomic sequences of genes on 3B subgenome and cDNA assemblies for the correspondence ortholog genes on 3A subgenome;
2. Map the cDNA sequences to its genomic sequences by GMAP;
3. Compare the gene structure generated by GMAP to gene structures annotated for the genomic sequences

Usage:
perl $0
<ortholog gene pairs list file, 3B_S293	3A_S1531>
<3B contig/BAC sequences in fasta>
<3B annotation file generated by PASA, gff3 format>
<3A cDNA sequences in fasta>
<3A annotation file generated by PASA, gff3 format>
<Output file>
END

die $sUsage unless @ARGV >= 6;
my ($ortho_pair_file, $genome_file, $anno_gff3_3B, $cdna_file, $anno_gff3_3A, $output) = @ARGV;

my %ortho_pairs = read_ortho_pair_file($ortho_pair_file);
my $genome_obj = Bio::DB::Fasta->new($genome_file);
my @ids = $genome_obj->ids;
print $ids[0], "\n";
my ($genomic_pos_ref, $gene_structure_ref, $gene_ctg_ref) = read_gff($anno_gff3_3B);
my %cdna_seqs = read_cdna_fasta($cdna_file, $anno_gff3_3A);

############################################
# Get the predicted gene structure by GMAP #
############################################
my %predicted_structure; my $n = 0;
foreach my $gene_id_3A (keys %cdna_seqs)
{
	my $gene_id_3B = $ortho_pairs{$gene_id_3A};
	next unless defined $gene_id_3B;
	$n++; print STDERR '******', $n, " ****\n";
	open (T, ">transcripts") or die;
	my @seqs = @{$cdna_seqs{$gene_id_3A}};
	foreach (0..$#seqs)
	{
		print T ">$gene_id_3A", "_3A", $_, "\n";
		print T $seqs[$_], "\n";
	}
	close T;		
	
	my $ctg = $gene_ctg_ref->{$gene_id_3B};
	unless (defined $ctg){print $gene_id_3B, "\n"; die "** no ctg for $gene_id_3B !!\n"}
	my $genomic_seq = $genome_obj->seq("$ctg", $genomic_pos_ref->{$gene_id_3B}[0] => $genomic_pos_ref->{$gene_id_3B}[1]);
	print $ctg, "\n" unless defined $genomic_seq;
	open (G, '>genome.db') or die;
	print G ">$gene_id_3B", "_genomic", "\n", $genomic_seq, "\n";
	close G;
	
	# GMAP
	open (GMAP, "perl ./process_GMAP_alignments.pl genome.db transcripts|") or die $!;
	while (<GMAP>)
	{
		my @t = split /\t/, $_;
		next unless @t==9 and $t[2] eq 'exon';
		#my $name = $1 if /Name=(\S+?)\;/;
		push @{$predicted_structure{$gene_id_3B}}, join("_", sort{$a<=>$b} @t[3,4]);
	}
	close GMAP;
	map{$predicted_structure{$_} = [unique(@{$predicted_structure{$_}})] } keys %predicted_structure;
	sleep(1);
	system("rm -rf ./gmap_db_dir/*; rm ./Makefile.genome.db; rm ./coords.genome.db");	
}

############################################
# Compare gene structure                   #
############################################
my @extra_exons;
my %count;
foreach my $gene_id (keys %predicted_structure)
{
	my ($altA, $altD, $extra_exon) = (0, 0, 0);
	my @extra_exons_tmp;
	my $flag; 
	my @predic_struc = @{$predicted_structure{$gene_id}};
	foreach my $ind (0 .. $#predic_struc)
	{
		my $overlap_ind = check_overlapping($predic_struc[$ind], $gene_structure_ref->{$gene_id});
		if ($overlap_ind == -1) # no overlapping
		{
			push @extra_exons_tmp, join("-", ($gene_id, $predic_struc[$ind], $ind)) if $ind > 0;
			next;
		}
		$flag = $ind;
		my ($p_s, $p_e) = split /_/, $predic_struc[$ind];
		my ($o_s, $o_e) = split /_/, $gene_structure_ref->{$gene_id}->[$overlap_ind];
		$altD++ if( ($p_e != $o_e) and ($ind != $#predic_struc) and ($overlap_ind != ((scalar @{$gene_structure_ref->{$gene_id}})-1)) );
		$altA++ if ($p_s != $o_s) and ($ind != 0) and ($overlap_ind != 0);		
	}
	$count{$gene_id} = [$altA, $altD, $extra_exon];
	
	foreach (@extra_exons_tmp)
	{
		my @t = split /-/, $_;
		push @extra_exons, $_ if $t[-1] < $flag;
	}
	
}

# Output
open (OUT, ">$output") or die "can't write $output \n";
map{print OUT join("\t", ($_, @{$count{$_}})), "\n" } keys %count;
close OUT;

############################################
# Subroutines                              #
############################################
sub check_overlapping
{
	my ($exon, $array_ref) = @_;
	my ($start, $end) = split /_/, $exon;
	my @array = @$array_ref;
	foreach my $ind (0..$#array)
	{
		my ($s, $e)  = split /_/, $array[$ind];
		return $ind if( ($start >= $s and $start < $e) or ($s >= $start and $s < $end) ) 
	}
	return -1;
}



sub read_ortho_pair_file
{
	my $file = shift;
	open (IN, $file) or die "Error opening file $file\n";
	
	my %return;
	while(<IN>)
	{
		# 3B_S293	3A_S1531
		chomp;
		next if /^\s+$/;
		s/3B_//; s/3A_//;
		my @t=split /\s+/,$_;
		$return{$t[1]} = $t[0];
	}
	close IN;
	return %return;
}

sub read_gff
{
	my $file = shift;
	open (IN, $file) or die "Error opening file $file\n";
	my (%genomic_pos, %gene_structure, %gene_ctg);
	my %recorder;
	while(<IN>)
	{
		chomp; 
		next unless /\S/;
		# gi|300681572|emb|FN645450.1|	PASA	cDNA_match	16861	17032	.	+	.	ID=S319-asmbl_397; Target=asmbl_397 1 172 +
		# gi|300681572|emb|FN645450.1|	PASA	cDNA_match	17453	17519	.	+	.	ID=S319-asmbl_397; Target=asmbl_397 173 239 +
		# gi|300681572|emb|FN645450.1|	PASA	cDNA_match	17597	17670	.	+	.	ID=S319-asmbl_397; Target=asmbl_397 240 313 +
		my $gene_id = $1 if /ID=(\S+?)\-/;
		my @t = split /\s+/, $_;
		$gene_ctg{$gene_id} = $t[0];
		push @{$recorder{$gene_id}}, join("_", @t[3, 4]); 
	}
	close IN;
	
	foreach my $id (keys %recorder)
	{
		my ($start, $end);
		my @exons = @{$recorder{$id}};
		my @tmp;
		foreach (@exons){push @tmp, (split/_/,$_)}
		@tmp = sort {$a<=>$b} @tmp;
		($start, $end) = @tmp[0, -1];
		$genomic_pos{$id} = [$start, $end];
		
		my @new_exons;
		foreach (@exons)
		{
			my @t = split /_/,$_; 
			@t=map{$_-$start+1}@t;
			push @new_exons, join('_', sort{$a<=>$b}@t);
		}
		$gene_structure{$id} = [unique(@new_exons)]
	}
	
	return(\%genomic_pos, \%gene_structure, \%gene_ctg)
}

sub unique
{
	my %hash = map {$_, 1} @_;
	return keys %hash;
}

sub read_cdna_fasta
{
	my $file = shift; 
	open (IN, $file) or die;
	my %return;
	my ($tmp, $id);
	while (<IN>)
	{
		# >ID=S1;Name=S1%20asmbl_2%20(S1);
    # GGTTCGTCGCTCTGCCGTCTCCCTCCCTCCCTTCGCCCCCGC
    # GGTTCGTCGCTCTGCCGTCTCCCTCCCTCCCTTCGCCCCCGC
    # GGTTCGTCGCTCTGCCGTCTCCCTCCCTCCCTTCGCCCCCGC
    chomp;
    next if /^\s+/;
    if(/>ID=(\S+?)\;/)
    {
    	$id = $1;
    	push @{$return{$id}}, $tmp if defined $tmp;
    	$tmp='';
    }
    else
    {
    	$tmp .= $_;
    }
	}
	close IN;
	push @{$return{$id}}, $tmp;
	return %return;
}



