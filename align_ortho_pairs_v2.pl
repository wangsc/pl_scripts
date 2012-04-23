#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SearchIO;

my $sUsage = <<END;
Usage:
perl $0
<ortholog gene pairs list file, 3B_S293	3A_S1531>
<3B contig/BAC sequences in fasta>
<3B annotation file generated by PASA, gff3 format>
<3A contig/scafold sequence in fasta>
<3A annotation file generated by PASA, gff3 format>
<output file>
END

die $sUsage unless @ARGV >= 5;

my ($ortho_file, $B_ctg_file, $B_gff, $A_ctg_file, $A_gff, $outfile) = @ARGV;
open (OUT, ">$outfile") or die;

my @ortho_A_B = read_ortho_file($ortho_file);
my $gn_3A = Bio::DB::Fasta->new($A_ctg_file);
my $gn_3B = Bio::DB::Fasta->new($B_ctg_file);

my ($B_genomic_pos, $B_gene_struc, $B_gene_ctg) = read_gff($B_gff);
my ($A_genomic_pos, $A_gene_struc, $A_gene_ctg) = read_gff($A_gff);

foreach (@ortho_A_B)
{
	chomp;
	my @out=split /\s+/,$_;
	my $out_file = join("-", @out) . ".fasta";

	my ($b_gene, $a_gene) = $_=~/3B_(\S+)\s+3A_(\S+)/;
	my $b_gene_seq = get_seq($gn_3B, $b_gene, $B_genomic_pos, $B_gene_struc, $B_gene_ctg);
	my $a_gene_seq = get_seq($gn_3A, $a_gene, $A_genomic_pos, $A_gene_struc, $A_gene_ctg);
	
	open (Q, ">query") or die;
	print Q ">$out[0]\n", $b_gene_seq, "\n";
	close Q;
	
	open (S, ">subject") or die;
	print S ">$out[1]\n", $a_gene_seq, "\n";
	close S;
	
	#open (my $BLAST, "blastn -query query -subject subject |") or die "can't open BLAST\n";
	#sleep(1);
	my $blast_cmd = "bl2seq -i query -j subject -p blastn -o blastn2seq.out";
	system($blast_cmd);
	print STDERR "Finish bl2seq ...\n";
	
	my $searchio = Bio::SearchIO->new(-format => 'blast', -file => "blastn2seq.out", -report_type => 'blastn');
	while (my $result = $searchio->next_result())
	{
		#last unless defined $result;
		my $hit = $result->next_hit;
		#last unless defined $hit;
		while (my $hsp = $hit->next_hsp)
		{
			my $query_string = $hsp->query_string;
			my $query_start = $hsp->start('query');
			my $query_string_formatted = format_sequence($query_string, $query_start, $B_gene_struc->{$b_gene});
			
			my $hit_string = $hsp->hit_string;
			my $hit_start = $hsp->start('hit');
			my $hit_string_formatted = format_sequence($hit_string, $hit_start, $A_gene_struc->{$a_gene});
			
			print $out[0], " "x (25-(length $out[0])), $query_string_formatted, "\n";
			print $out[1], " "x (25-(length $out[1])), $hit_string_formatted, "\n";
			print "\n";
			my @count_structural_alter = compare_structure($query_string_formatted, $hit_string_formatted);
			print OUT join("-", @out), "\t", join("\t",(@count_structural_alter)), "\n";
			last;
		}		
		print '*'x 20, "\n";		
	}	
}
close OUT;


# Subroutines
sub compare_structure
{
	my ($seqa, $seqb) = @_;
	print STDERR length  $seqa, " *** ", length  $seqb, "\n";
	my @struc_a = get_structure($seqa);
	my @struc_b = get_structure($seqb);
	my @count_seg = count_coding_segments(\@struc_a, \@struc_b, length $seqa);
	return @count_seg;
}

sub get_structure
{
	my $seq = shift; 
	my @return;
	my @bases = split //, $seq;
	print STDERR 'length of seq: ', length $seq, "\t", 'length of array: ', scalar @bases, "\n";
	my $start;
	my $end;
	foreach my $index (0..$#bases)
	{
		if($index == $#bases and defined $start)
		{
			if($bases[$index] eq '-')
			{
				$end = $index;
				while($bases[$end] eq '-')
				{
					$end--;
				}
				push @return, [$start, $end];
			}
			next
		}
		
		if($bases[$index] =~ /[ATGCN]/)
		{
			$start = $index unless defined $start;
			if ($index == $#bases){push @return, [$start, $end]}
		}
		elsif($bases[$index] =~ /[atgcn]/)
		{
			if(defined $start)
			{
				$end = $index - 1;
				while($bases[$end] eq '-')
				{
					$end--;
				}
				push @return, [$start, $end];
				print STDERR join(",", ($start, $end)), "\n";
				undef($start);
			}
		}
	}
	return @return;
}

sub format_sequence
{
	my ($seq, $start, $structures) = @_;
	my %exon_pos;
	foreach (@$structures)
	{
		my ($start, $end) = split /_/, $_;
		map{ $exon_pos{$_}=1 } $start .. $end;
	}
	$seq =~ s/\s+$//;
	$seq = lc($seq);
	my $pos = $start - 1; 
	foreach my $index (0 .. (length $seq)-1)
	{
	 	my $base = substr $seq, $index, 1;
	 	next if $base eq '-';
	 	$pos++;
	 	substr ($seq, $index, 1) = uc($base) if exists $exon_pos{$pos};
	}
	
	return 'aaaaa' . $seq . 'aaaaa';	
}




sub get_seq
{
	my ($gn, $gene, $genomic_pos, $gene_struc, $gene_ctg) = @_;
	my $seq = $gn->seq($gene_ctg->{$gene}, $genomic_pos->{$gene}[0] => $genomic_pos->{$gene}[1]);
	$seq = lc ($seq);
	foreach (@{$gene_struc->{$gene}})
	{
		my ($start, $end) = split /_/, $_;
		my $subseq = substr ($seq, $start-1, ($end-$start+1));
		substr ($seq, $start-1, ($end-$start+1)) = uc($subseq);
	}
	return $seq;
}


sub read_ortho_file
{
	my $file = shift;
	open (IN, $file) or die;
	
	my @return;
	while (<IN>)
	{
		# 3B_S293	3A_S1531
		# 3B_S206	3A_S2560
		# 3B_S356	3A_S5665		
		chomp; next if /^\s+$/;

		push @return, $_;
	}
	close IN;
	return @return;
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

sub count_coding_segments
{
	my ($wise_ref, $pasa_ref, $seq_len) = @_;
	my @return_array;
	my $num_exons = 0;
	my $num_compare = 0;
	my ($ade_l, $ade_g, $aae_l, $aae_g, $ri, $si, $se, $re, $ce) = (0,0,0,0,0,0,0,0,0);
	my ($alt_acc_exon, $alt_donor_exon) = (0, 0);
	my ($wise_vec) = construct_vec($wise_ref);
	my ($pasa_vec) = construct_vec($pasa_ref);
#	my $coverage = calculate_coverage($wise_vec, $pasa_vec, $wise_max);
	my @coding_segments = calculate_coding_segment($wise_vec, $pasa_vec, $seq_len);
	my $total_segs = scalar @coding_segments;
	my $alternative_segs = 0;
	foreach  my $index (0..$#coding_segments)
	{
		my $segment = $coding_segments[$index];
		my $status = $segment->[2];
		next if $status == 0;
		if ($status == 11){$ce++ if $coding_segments[$index-1]->[2] == 0 and $coding_segments[$index+1]->[2] == 0; next}
		$alternative_segs++;
		$ade_g++ if $coding_segments[$index-1]->[2] == 11 and $coding_segments[$index+1]->[2] == 0 and $status == 1;
		$ade_l++ if $coding_segments[$index-1]->[2] == 11 and $coding_segments[$index+1]->[2] == 0 and $status == 10;
		$aae_g++ if $coding_segments[$index-1]->[2] == 00 and $coding_segments[$index+1]->[2] == 11 and $status == 1;
		$aae_l++ if $coding_segments[$index-1]->[2] == 00 and $coding_segments[$index+1]->[2] == 11 and $status == 10;
		$ri++ if $coding_segments[$index-1]->[2] == 11 and $coding_segments[$index+1]->[2] == 11 and $status == 10;
		$si++ if $coding_segments[$index-1]->[2] == 11 and $coding_segments[$index+1]->[2] == 11 and $status == 1;
		$se++ if $coding_segments[$index-1]->[2] == 0 and $coding_segments[$index+1]->[2] == 0 and $status == 1;
		$re++ if $coding_segments[$index-1]->[2] == 0 and $coding_segments[$index+1]->[2] == 0 and $status == 10;
	}
	@return_array = ($ade_g, $ade_l, $aae_g, $aae_l, $ri, $si, $se, $re, $ce);
	#@return_array = map{$_/$total_segs} @return_array;
	# return ($alternative_segs, $total_segs, \@coding_segments, @return_array);
	return @return_array;
}

sub calculate_coding_segment
{
	my ($wise_vec, $pasa_vec, $max)= @_;
	my @coding_segs;
	my ($seg_start, $seg_end) = (0, 0);
	my $pre_status = 0;
	my $current_status;
	foreach ( 0..$max+3) # 3 is used to make a small blank
	{
		my $w = vec($wise_vec, $_,1);
		$w = 0 unless defined $w;
		my $p =  vec($pasa_vec, $_,1);
		$p = 0 unless defined $p;
		$current_status = $w*10+$p;
		if ($current_status == $pre_status)
		{
			if ($_ == ($max+3))
			{
				push @coding_segs, [$seg_start, $max+3, $pre_status]
			}
			next;
		}
		$seg_end = $_ - 1;
		push @coding_segs, [$seg_start, $seg_end, $pre_status];
		$seg_start = $_;
		$pre_status = $current_status;		
	}
	return (@coding_segs);
}

sub construct_vec
{
	my $arrayref = shift;
	#my $flag= shift;
	my $vec = '';
	my $max;
	my $total;
	my $debug =1 ;
	my @array = @$arrayref;
	#if($flag){shift @array; pop @array} # remove terminal exons SW 10.06.2011
	foreach (@array)
	{
		my @d = sort{$a<=>$b}@$_;
#		print '@d: ', join("\t", @d),"\n" if $debug; $debug=0;
		foreach ($d[0]..$d[1])
		{
			vec($vec,$_,1) = 0b1;
		#	$total++;
			$max = $_ unless defined $max;
			$max = $_ if $_ > $max;
			
		}
	}
	return ($vec, $max);
}