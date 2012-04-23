#!/usr/bin/perl -w
use strict;

# detect conserved exon(CE), acquired exon(AE) and missing Exon(ME) by 
# comparing rice or brachypodium gene structures (predicted by genewise) to wheat gene structures (predicted by PASA);

my $sUsage = qq(
perl $0
<rice protein genewise output>
<brachypodium gff3>
<output file>
);
die $sUsage unless @ARGV >= 3;

my ($genewise_file, $pasa_gff_file, $out_file) = @ARGV;
my %genewise_results = read_wise_gff($genewise_file);
my ($pasa_iso_gene, $pasa_gene_vec) = read_brachy_gff($pasa_gff_file);

open (OUT, ">$out_file") or die;
# Compare

foreach my $asmbl (keys %genewise_results)
{
	my ($genewise_vec, $genewise_min, $genewise_max) = @{$genewise_results{$asmbl}};
	my $gene = $pasa_iso_gene->{$asmbl};
	my ($pasa_vec, $pasa_min, $pasa_max) = @{$pasa_gene_vec->{$gene}};
	my @coding_segs = calculate_coding_segment($genewise_vec, $pasa_vec, $genewise_min, $genewise_max);
	my $effective_coding_seg = 0;
	foreach (@coding_segs)
	{
		$effective_coding_seg ++ if $_->[2] != 0;
	}
	$effective_coding_seg-- if $coding_segs[0]->[2]==1;
	$effective_coding_seg-- if $coding_segs[-1]->[2]==1;
	
	# Acquired exons [.. 01 ..]
	my @AE;
	foreach (1..$#coding_segs-1)
	{
		push @AE, $coding_segs[$_] if $coding_segs[$_]->[2] == 1;
	}
	
	# Conserved exons [11]
	my @CE;
	foreach (1..$#coding_segs-1)
	{
		if($_==1){push @CE, $coding_segs[$_] if $coding_segs[$_]->[2]==11 and $coding_segs[$_-1]->[2]==0 and $coding_segs[$_+1]->[2]==0}
		elsif( $_ == $#coding_segs-1){push @CE, $coding_segs[$_] if $coding_segs[$_]->[2]==11 and $coding_segs[$_-1]->[2]==0 and $coding_segs[$_+1]->[2]==0}
		else
		{
			push @CE, $coding_segs[$_] if $coding_segs[$_]->[2]==11;
		}
	}
	
	# Missing exons [.. 10 ..]
	my @ME;
	foreach (1..$#coding_segs-1)
	{
		push @ME, $coding_segs[$_] if $coding_segs[$_]->[2] == 10;
	}	
	print OUT join("\t", ($gene, (scalar @AE)/$effective_coding_seg, (scalar @CE)/$effective_coding_seg,  (scalar @ME)/$effective_coding_seg)),"\n";
	print ">", join("\t", ($gene, $genewise_max-$genewise_min)), "\n";
	my @arr = (); my $length = 0;
	foreach (@AE){push @arr, join(":", @{$_}[0,1]); $length += ($_->[1] - $_->[0] + 1)}print join("\t", ($length, @arr)), "\n";
	@arr = (); $length = 0;
	foreach (@CE){push @arr, join(":", @{$_}[0,1]); $length += ($_->[1] - $_->[0] + 1)}print join("\t", ($length, @arr)), "\n";
	@arr = (); $length = 0;
	foreach (@ME){push @arr, join(":", @{$_}[0,1]); $length += ($_->[1] - $_->[0] + 1)}print join("\t", ($length, @arr)), "\n";
}

close OUT;


# Subroutiens
sub calculate_coding_segment
{
	my ($wise_vec, $pasa_vec, $min, $max)= @_;
	my @coding_segs;
	my ($seg_start, $seg_end) = (0, 0);
	my $pre_status = 0;
	my $current_status;
	foreach ( $min-3..$max+3) # 3 is used to make a small blank
	{
		my $w = vec($wise_vec, $_,1);
		$w = 0 unless defined $w;
		my $p =  vec($pasa_vec, $_,1);
		$p = 0 unless defined $p;
		$current_status = $w*10+$p;
		if ($_ == $min-3) {$pre_status = $current_status; next}
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

sub read_brachy_gff
{
	my $file = shift;
	my %iso_gene;
	my %gene_exons;
	my ($start, $end);
	open (IN, "$file") or die "can't open file $file\n";
	my $debug = 1;
	
	while (<IN>)
	{
		next if /^\s+$/;
		my @data = split /\s+/, $_;
		if($data[2] eq "mRNA")
		{
			my $id = $1 if /ID=(\S{14})/;
			my $gene = $1 if /ID=(\S{12})/;
			$iso_gene{$id} = $gene;
			($start, $end) = @data[3,4];
			next;
		}
		next unless $data[2] eq 'exon';
		my $gene = $1 if /Parent=(\S{12})/;
		push @{$gene_exons{$gene}}, [ $data[3]-$start+1, $data[4]-$start+1 ];
	}
	close IN;
	map{$gene_exons{$_} = construct_vec($gene_exons{$_})} keys %gene_exons;
	return (\%iso_gene, \%gene_exons);
}

sub construct_vec
{
	my $arrayref = shift;
	my $vec = '';
	my $max;
	my $min;
	my $debug =1 ;
	foreach (@$arrayref)
	{
		my @d = sort{$a<=>$b}@$_;
#		print '@d: ', join("\t", @d),"\n" if $debug; $debug=0;
		foreach ($d[0]..$d[1])
		{
			$max = $_ unless defined $max;
			$max = $_ if $_ > $max;
			$min = $_ unless defined $min;
			$min = $_ if $_ < $min;
			vec($vec,$_,1) = 0b1;
		}
	}
	foreach ($min..$max)
	{
		my $value = vec($vec, $_, 1);
		vec($vec, $_, 1) = $value==1?$value:0;
	}
	return [$vec, $min, $max];
}

sub read_wise_gff
{
	my $file = shift;
	my %return_hash;
	open (IN, $file) or die $!;
	my $score;
	my $genewise_cutoff = 35;
	my @array;
	my $iso_id;
	my $flag;
	while (<IN>)
	{
		next if /^\s+$/;
		next if /^\/\//;
		my @t = split /\s+/,$_;
		if (/match/){$flag = 0; $flag = 1 if $t[5] >= $genewise_cutoff; next}
		push @{$return_hash{$t[0]}}, [sort {$a<=>$b} @t[3,4]] if $flag and $t[2] eq 'cds';
	}
	close IN;
	map{$return_hash{$_} = construct_vec($return_hash{$_})} keys %return_hash;
	return %return_hash;
}

sub check_overlap
{
	my ($region_a, $region_b) = @_;
	return 1 if $region_a->[0] >= $region_b->[0] and $region_a->[0] <= $region_b->[1];
	return 1 if $region_b->[0] >= $region_a->[0] and $region_b->[0] <= $region_a->[1];
	return 0;
}

sub max
{
	my $m = shift;
	map{$m = $_ if $_>$m} @_;
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

