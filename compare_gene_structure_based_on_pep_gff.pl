#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

die <<END unless @ARGV >= 5;
**************************************************************************
This scritp will compare the gene structure of orthologs of wheat 
and rice (or brachypodium).
1. read the output of blastx (blast wheat assemblies to rice protein db)
   to get the orthologs;
2. read the rice gff file and wheat gff file (from PASA)
3. compare gene structure

Usage:
perl $0
<blastx output>
<rice or brachypodium gff file>
<wheat gff file>
<flag, 1 for rice, 0 for brachypodium>
<output file>
**************************************************************************
END

my ($blast_result, $orth_gff, $wheat_gff, $flag, $output) = @ARGV;

my ($blast_ref, $best_frame_ref) = parse_blast_result($blast_result);

my %other_gff = read_other_gff($orth_gff, $flag);
my %wheat_gff = read_wheat_gff($wheat_gff);
my %exon_structure = transform_into_exon_structure($blast_ref, $best_frame_ref, \%other_gff, \%wheat_gff);
my ($total, %compare_results) = compare_exon_structure(\%exon_structure);
output($total, \%compare_results, $output);

# Subroutines
sub output
{
	# $return_hash{$id} = [$left_alt, $right_alt, $conserved, $hc_exon, $retained_intron, $join_exon, $splice_exon, $skip_exon];
	my($total, $hashref, $file) = @_;
	open (OUT, ">$file") or die "$!\n";
	my @record;
	my @names = qw(left_alt right_alt conserved hc_exon retained_intron join_exon splice_exon skip_exon);
	map{$record[$_] = 0}(0..$#names);
	print OUT 'ID',"\t", join("\t", @names),"\n";
	foreach my $id (keys %$hashref)
	{
		print OUT $id,"\t", join("\t", @{$hashref->{$id}}),"\n";
		my $num = scalar (@{$hashref->{$id}});
		foreach (0 .. $num-1)
		{
			$record[$_] ++ if $hashref->{$id}->[$_] > 0;
		}
	}
	close OUT;
	
	print 'Total',"\t", $total,"\n";
	foreach (0..$#names)
	{
		print $names[$_],"\t",$record[$_],"\n";
	}
}


sub compare_exon_structure
{
	my ($exon_ref) = @_;
	my %return_hash;
	my($total_genes);
	foreach my $id (keys %$exon_ref)
	{
		$total_genes++;
		my($left_alt, $right_alt, $conserved, $hc_exon, $retained_intron, $join_exon, $splice_exon, $skip_exon) = (0,0,0,0,0,0,0,0);
		my ($orth_vec, $orth_max, $orth_total) = construct_vec($exon_ref->{$id}->[0]);
		my ($wheat_vec, $wheat_max, $wheat_total) = construct_vec($exon_ref->{$id}->[1]);
		foreach (@{$exon_ref->{$id}->[0]})
		{
			my ($start, $end) = @$_;
		#	print STDERR '$start, $end ', $start, "\t", $end, "\n" if $id =~ /asmbl_1222/;
			my $covered_by_wheat = 0;
			foreach ($start..$end)
			{
				$covered_by_wheat++ if (vec($wheat_vec, $_,1) == 1);
				#print $_, "\t",  $covered_by_wheat,"\n" if (vec($wheat_vec, $_,1) == 1);
			}
		#	print '$covered_by_wheat ', $covered_by_wheat,"\n" if $id =~ /asmbl_1222/;
			#
			if ($covered_by_wheat >= ($end-$start+1)*0.9)
			{
				$hc_exon++ 
			}
			elsif ($covered_by_wheat > 5)
			{
				$retained_intron++ if (vec($wheat_vec, $start-1,1) == 1 or vec($wheat_vec, $end+1,1) == 1);
				$splice_exon++;
			}
			else
			{
				$skip_exon++
			}
			#
			foreach my $p_ref (@{$exon_ref->{$id}->[1]})
			{
				$conserved++ if $p_ref->[0] == $start and $p_ref->[1] == $end;
				$left_alt ++ if $p_ref->[0] != $start and $p_ref->[1] == $end;
				$right_alt++ if $p_ref->[0] == $start and $p_ref->[1] != $end;
			}
		}
		#
		foreach my $i (0..( (scalar @{$exon_ref->{$id}->[0]})-2 ))
		{
			my ($t, $c) = (0, 0);
			foreach ($i .. $i+1)
			{
				my ($start, $end) = @{$exon_ref->{$id}->[0]->[$_]};
				$t += ($end-$start+1);
				map{$c++ if (vec($wheat_vec, $_,1) == 1)} ($start, $end);
			}
			$join_exon++ if $c >= $t*0.9;
		}
		
		foreach (@{$exon_ref->{$id}->[1]})
		{
			my ($start, $end) = @$_;
	#		print STDERR '$start, $end ', $start, "\t", $end, "\n" if $id =~ /asmbl_1222/;
			my $covered_by_orth = 0;
			map{$covered_by_orth++ if (vec($orth_vec, $_,1) == 1)} ($start..$end);
			$retained_intron++ if $covered_by_orth <= ($end-$start+1)*0.9;
	#		print '$covered_by_orth ', $covered_by_orth,"\n" if $id =~ /asmbl_1222/;
		}
		$return_hash{$id} = [$left_alt, $right_alt, $conserved, $hc_exon, $retained_intron, $join_exon, $splice_exon, $skip_exon];
	}
	return ($total_genes, %return_hash);
}

sub construct_vec
{
	my $arrayref = shift;
	my $vec = '';
	my $max;
	my $total;
	my $debug =1 ;
	foreach (@$arrayref)
	{
		my @d = sort{$a<=>$b}@$_;
#		print '@d: ', join("\t", @d),"\n" if $debug; $debug=0;
		foreach ($d[0]..$d[1])
		{
			vec($vec,$_,1) = 0b1;
			$total++;
			$max = $_ unless defined $max;
			$max = $_ if $_ > $max;
			
		}
	}
	return ($vec, $max, $total);
}

sub transform_into_exon_structure
{
	my ($blast_ref, $frame_ref, $other_gffref, $wheat_gffref) = @_;
	my %return_hash;
	foreach my $ortho_pair (keys %{$blast_ref})
	{
		my ($wheat_id, $ortho_id) = split /__/, $ortho_pair;
		my @ortho_exons = transfrom_ortholog_exon($blast_ref, $other_gffref, $ortho_id, $ortho_pair);
		my @wheat_exons = tranform_wheat_exon($blast_ref, $wheat_gffref, $wheat_id, $ortho_pair, $frame_ref);
		$return_hash{$ortho_pair} = [[@ortho_exons],[@wheat_exons]];
	}
	return %return_hash;
}

sub transfrom_ortholog_exon
{
	my ($blast_ref, $other_gffref, $ortho_id, $ortho_pair) = @_;
	my @range = @{$blast_ref->{$ortho_pair}}[2,3];
	my @exons;
	my @return;
	my @cds = @{$other_gffref->{$ortho_id}{'CDS'}};
	@cds = sort{$a->[0] <=> $b->[0]} @cds;
	@cds = reverse @cds if $other_gffref->{$ortho_id}{'strand'} eq '-';
	my $count = 0;
	foreach (@cds)
	{
		my ($start, $end, $clip) = @$_;
		my $num_base = (abs($end - $start) + 1 -$clip);
		my $num_aa = int($num_base/3);
		push @exons, [$count+1, $count+$num_aa];
		$count += $num_aa;
	}
	foreach (@exons)
	{
		my ($start, $end) = @$_;
		last if $start > $range[1];
		if ($start <= $range[0] and $end >= $range[0])
		{
			push @return, [max($start,$range[0]), min($end,$range[1])];
		}
		elsif($start <= $range[1] and $end >= $range[1])
		{
			push @return, [$start, $range[1]];
		}
		
	}
	return @return;
}

sub max
{
	my $max = shift;
	foreach(@_){$max = $_ if $_ >= $max}
	return $max;
}

sub min
{
	my $min = shift;
	foreach(@_){$min = $_ if $_ < $min}
	return $min;
}

sub tranform_wheat_exon
{
	my ($blast_ref, $wheat_gffref, $wheat_id, $ortho_pair, $frame_ref) = @_;
	my @return;
	my @range = @{$blast_ref->{$ortho_pair}}[2,3];
	my $offset;
	my $frame = $frame_ref->{$wheat_id};
	my $s = int($blast_ref->{$ortho_pair}->[0]- abs($frame))/3 + 1;
	$offset = $blast_ref->{$ortho_pair}->[2] - $s;
	@range = map{$_ - $offset} @range;
	my @cds  = @{$wheat_gffref->{$wheat_id}};
	@cds = sort {$a->[0] <=> $b->[0]} @cds;
	@cds = reverse @cds if $frame < 0;
	my $count = 0;
	my @exons;
	foreach (@cds)
	{
		my ($start, $end) = @$_;
		my $num_base = (abs($end - $start) + 1 );
		$num_base -= abs($frame) unless $count;
		my $num_aa = int($num_base/3);
		push @exons, [$count+1, $count+$num_aa];
		$count += $num_aa;		
	}
	foreach (@exons)
	{
		my ($start, $end) = @$_;
		last if $start > $range[1];
		if ($start <= $range[0] and $end >= $range[0])
		{
			push @return, [max($start,$range[0]) + $offset, min($end,$range[1]) + $offset];
		}
		elsif($start <= $range[1] and $end >= $range[1])
		{
			push @return, [$start + $offset, $range[1] + $offset];
		}
		
	}
	return @return;
}



sub parse_blast_result
{
	my $file = shift;
	my $min_len = 500; #
	my %return_hash;
	my %frame;
	my $blast = Bio::SearchIO->new(-foramt=>'blast', -file => $file);
	my $debug = 1;
	while (my $result = $blast->next_result())
	{
		my $query_name = $result->query_name;
		my $hit = $result->next_hit(); next unless defined $hit;
		next if $hit->length < $min_len;
		my $hit_name = $hit->name;
		if ($hit_name =~ /(.*?)\|/){$hit_name = $1 } else{my @t=split /\s+/,$hit_name; $hit_name = $t[0]}
		print '$query_name ', $query_name,"\n", '$hit_name ', $hit_name,"\n" if $debug; $debug = 0;
		my $hsp = $hit->next_hsp();
		
		$return_hash{join('__', ($query_name, $hit_name))} = [$hsp->range('query'), $hsp->range('hit')];
		$frame{$query_name} = $hsp->frame('query');
	}
	return (\%return_hash, \%frame);
}

sub read_other_gff
{
	my $file = shift;
	my $rice_flag = shift;
	return read_rice_gff($file) if $rice_flag;
	return read_bd_gff($file);
}

sub read_rice_gff
{
	my $file = shift;
	my %return_hash;
	open (IN, $file) or die $!;
	
	#Chr1	MSU_osa1r6	mRNA	15321	19323	.	+	.	ID=13101.m00007;Parent=13101.t00004;Alias=LOC_Os01g01040.3
	#Chr1	MSU_osa1r6	five_prime_UTR	15321	15598	.	+	.	Parent=13101.m00007
	#Chr1	MSU_osa1r6	CDS	15599	15976	.	+	0	Parent=13101.m00007
	#Chr1	MSU_osa1r6	CDS	16383	16474	.	+	0	Parent=13101.m00007
	#Chr1	MSU_osa1r6	CDS	16558	17258	.	+	2	Parent=13101.m00007
	#Chr1	MSU_osa1r6	CDS	17501	17571	.	+	1	Parent=13101.m00007
	#Chr1	MSU_osa1r6	CDS	17968	18057	.	+	0	Parent=13101.m00007
	#Chr1	MSU_osa1r6	CDS	18142	18321	.	+	0	Parent=13101.m00007
	#Chr1	MSU_osa1r6	CDS	18531	18593	.	+	0	Parent=13101.m00007
	#Chr1	MSU_osa1r6	three_prime_UTR	18594	19323	.	+	.	Parent=13101.m00007
	my %id_alias;
	my $debug = 1;
	while (<IN>)
	{
		next if /^\s+$/ or /^\#/;
		chomp;
		my @data = split /\t/, $_;
		if ($data[2] eq 'mRNA')
		{
			my ($id, $alias) = $data[8] =~ /ID=(.*?)\;.*Alias=(.*)/;
			#die '$id, $alias ', join("  ", ($id, $alias)),"\n" if $debug;
			$id_alias{$id} = $alias if defined $alias;
			next;
		}
		if ($data[2] eq 'CDS')
		{
			my $id = $1 if $data[8] =~ /Parent=(.*)/;
			next unless exists $id_alias{$id};
			print 'Rice ID ', $id,"\n" if $debug; $debug = 0;
			my $strand = $data[6];
			$return_hash{$id_alias{$id}}{'strand'} = $strand unless exists $return_hash{$id_alias{$id}}{'strand'};
			push @{$return_hash{$id_alias{$id}}{'CDS'}}, [@data[3,4,7]];
		}
	}
	close IN;
	return %return_hash;
}

sub read_bd_gff
{
	my $file = shift;
	my %return_hash;
	open (IN, $file) or die $!;
	my $debug = 1;
# Bd1	MIPS_v1.2	predicted_gene	55668	56319	.	+	.	ID=Bradi1g00220;Name=BdGRAS1;Description=GRAS  transcription factor
# Bd1	MIPS_v1.2	mRNA	55668	56319	.	+	.	ID=Bradi1g00220.1;Parent=Bradi1g00220
# Bd1	MIPS_v1.2	exon	55668	56319	.	+	.	ID=Bradi1g00220.1_exon_1;Parent=Bradi1g00220.1
# Bd1	MIPS_v1.2	CDS_predicted	55668	56120	.	+	0	ID=Bradi1g00220.1_cds_1;Parent=Bradi1g00220.1

	while (<IN>)
	{
		next if /^\s+$/ or /^\#/;
		chomp;
		my @data = split /\t/, $_;
#		if ($data[2] eq 'mRNA')
#		{
#			$id = $1 if $data[8] =~ /ID=(.*?)\;/;
#			next;
#		}
		if ($data[2] =~ /CDS/)
		{
			my $id = $1 if $data[8] =~ /ID=(.*)_cds/;
			print 'BD ID ', $id,"\n" if $debug; $debug = 0;
			my $strand = $data[6];
			$return_hash{$id}{'strand'} = $strand unless exists $return_hash{'strand'};
			push @{$return_hash{$id}{'CDS'}}, [@data[3,4,7]];
		}
	}
	close IN;
	return %return_hash;
}

sub read_wheat_gff
{
	my $file = shift;
	my %return_hash;
	my $debug = 1;
	open (IN, "$file") or die $!;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my ($id, $start, $end) = $_=~ /Target=(\S+)\s(\S+)\s(\S+)/;
		print 'wheat ID ', $id, "\n" if $debug;$debug=0;
		push @{$return_hash{$id}}, [$start, $end];
	}
	close IN;
	return %return_hash;
}
