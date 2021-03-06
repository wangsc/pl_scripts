#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
# This script will read two specifed gff-like file (which is different with normal gff)
# and compare the exon structures in these two files.
# 07.06.2011 SW
perl $0
<gff3 file generated by script run_genewise_wrapper.pl>
<assemblies gff3 file generate by PASA pipeline>
<pasa_chr3A_clean_masked.pasa_assemblies.denovo_transcript_isoforms.gff3>
<exon groups output file>
);
die $sUsage unless @ARGV >= 4;
my($wise_gff, $pasa_gff, $iso_gff, $exon_outfile) = @ARGV;

my %wise_str = read_wise_gff($wise_gff);
my %pasa_str = read_pasa_gff($pasa_gff);
my %iso_gene = read_isoform_gff($iso_gff);
my %gene_iso;
foreach (keys %iso_gene)
{
	push @{ $gene_iso{$iso_gene{$_}} }, $_;
}
my %summary;
my %num_coding_each_gene;
my %acquir_alter_exon;
my %average_coverage;
foreach my $gene (keys %gene_iso)
{
	my $total_coding = 0;
	my $total_alt_coding = 0;
	my $total_cov = 0;
	my $num_comparison = 0;
	foreach my $iso_id (@{$gene_iso{$gene}})
	{
		my $pasa_exon_ref = $pasa_str{$iso_id};
		#print STDERR $iso_id,"\n" unless exists $wise_str{$iso_id};
		foreach my $exon_ref (@{$wise_str{$iso_id}})
		{
			my ($alt_seg, $total_seg, $coverage, $coding_seq_ref, @compare_result) = compare_structure($pasa_exon_ref, $exon_ref);
			$total_coding += $total_seg;
			$total_alt_coding += $alt_seg;
			$total_cov += $coverage;
			$num_comparison++;
			push @{$summary{$gene}}, [@compare_result];
		}		
		
	}
	next unless $num_comparison;
	$average_coverage{$gene} = $total_cov/$num_comparison;
	$num_coding_each_gene{$gene} = $total_coding;
}

output (\%summary, \%average_coverage, $exon_outfile);



# Subroutines
sub output
{
	# $return_hash{$id} = [$left_alt, $right_alt, $conserved, $hc_exon, $retained_intron, $join_exon, $splice_exon, $skip_exon];
	my($summary_ref, $cov_ref, $exon_out) = @_;
	open (E, ">$exon_out") or die "$!\n";
	my @record;
	my @names = qw(ALternative_donor_exon Alternative_acceptor_exon Retained_intron Skipped_exon Skipped_exon Retained_exon Conserved_exon Coverage);
	print E 'ID',"\t", join("\t", @names),"\n";
	# ($ade, $aae, $ri, $si, $se, $re, $ce);
	foreach my $gene (keys %{$summary_ref})
	{
		my $t = scalar (@{$summary_ref->{$gene}});
		my @count;
		foreach my $iso (@{$summary_ref->{$gene}})
		{
			foreach (0..(scalar @$iso -1))
			{
				$count[$_] += $iso->[$_];
			}
		}
		@count = map{$_/$t}@count;
		print E join("\t", ($gene, @count)),"\t", $cov_ref->{$gene}, "\n";
	}
	close E;
}

sub read_isoform_gff
{
	my $file = shift;
	my %iso_gene;
	my %iso_contig;
	open (IN, "$file") or die "$! $file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my ($gene, $iso_id) = $_=~/ID=(.*?)\-(.*?)\;/;
		$iso_gene{$iso_id} = $gene;
		my @data = split /\t/, $_;
		$iso_contig{$gene} = shift @data;
	}
	close IN;
	return %iso_gene;
}


sub read_wise_gff
{
	my $file = shift;
	my %return_hash;
	open (IN, $file) or die $!;
	my $score;
	my $genewise_cutoff = 0;
	my @array;
	my $iso_id;
	while (<IN>)
	{
		next if /^\s+$/;
		if (/^\/\//)
		{
			@array = sort{$a->[0] <=> $b->[0]} @array;
			push @{$return_hash{$iso_id}}, [@array] if @array;
			@array = ();
			next;
		}
		my @t = split /\t/, $_;
		$iso_id = $t[0];
		if ($t[2] =~ /match/){$score = $t[5]}
		next unless $t[2] =~ /cds/i;
		next unless $score > $genewise_cutoff;
		push @array,[sort{$a<=>$b}@t[3, 4]];		
	}
	close IN;
	return %return_hash;
}

sub read_pasa_gff
{
	my $file = shift;
	my %return_hash;
	open (IN, $file) or die $!;
	while (<IN>)
	{
		next if /^\s+$/;
		my $id = $1 if /Target=(\S+)\s/;
		my @t = split /\t/, $_;
	#	print $id, "\t", join("\t", @t[3, 4]),"\n" if $id =~ /asmbl_1222/;
		push @{$return_hash{$id}},[@t[3, 4]];		
	}
	close IN;
	map{ $return_hash{$_} = [ sort{$a->[0]<=>$b->[0]} @{$return_hash{$_}} ] }keys %return_hash;
	return %return_hash;
}

sub max
{
	
	my $m = shift;
	#map{$m = $_ if $_>$m} @_;
	foreach (@_){$m = $_ if $_>$m}
	return $m;
}

sub min
{
	my $m = shift;
	map{$m = $_ if $_<$m} @_;
	return $m;
}

sub calculate_coverage
{
	my ($wvec, $pvec, $max) = @_;
	my $n = 0;
	foreach (1..$max)
	{
		$n++ if (vec($wvec, $_, 1)==1 and vec($pvec, $_, 1)==1);
	}
	return $n/$max;
}

sub compare_structure
{
	my ($wise_ref, $pasa_ref) = @_;
	my @return_array;
	my $num_exons = 0;
	my $num_compare = 0;
	my ($ade, $aae, $ri, $si, $se, $re, $ce) = (0,0,0,0,0,0,0);
	my ($alt_acc_exon, $alt_donor_exon) = (0, 0);
	my ($wise_vec, $wise_max) = construct_vec($wise_ref);
	my ($pasa_vec, $pasa_max) = construct_vec($pasa_ref);
	my $coverage = calculate_coverage($wise_vec, $pasa_vec, $wise_max);
	my @coding_segments = calculate_coding_segment($wise_vec, $pasa_vec, $wise_max>$pasa_max?$wise_max:$pasa_max);
	my $total_segs = scalar @coding_segments;
	my $alternative_segs = 0;
	foreach  my $index (0..$#coding_segments)
	{
		my $segment = $coding_segments[$index];
		my $status = $segment->[2];
		next if $status == 0;
		if ($status == 11){$ce++ if $coding_segments[$index-1]->[2] == 0 and $coding_segments[$index+1]->[2] == 0; next}
		$alternative_segs++;
		$ade++ if $coding_segments[$index-1]->[2] == 11 and $coding_segments[$index+1]->[2] == 0;
		$aae++ if $coding_segments[$index-1]->[2] == 00 and $coding_segments[$index+1]->[2] == 11;
		$ri++ if $coding_segments[$index-1]->[2] == 11 and $coding_segments[$index+1]->[2] == 11 and $status == 10;
		$si++ if $coding_segments[$index-1]->[2] == 11 and $coding_segments[$index+1]->[2] == 11 and $status == 1;
		$se++ if $coding_segments[$index-1]->[2] == 0 and $coding_segments[$index+1]->[2] == 0 and $status == 1;
		$re++ if $coding_segments[$index-1]->[2] == 0 and $coding_segments[$index+1]->[2] == 0 and $status == 10;
	}
	@return_array = ($ade, $aae, $ri, $si, $se, $re, $ce);
	@return_array = map{$_/$total_segs} @return_array;
	return ($alternative_segs, $total_segs, $coverage, \@coding_segments, @return_array);
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
		#	$total++;
			$max = $_ unless defined $max;
			$max = $_ if $_ > $max;
			
		}
	}
	return ($vec, $max);
}