#!/usr/bin/perl -w
use strict;

die <<END unless @ARGV >= 3;

 *******************************************************************************************************************************
  SW 12.16.2011
  This script will calcualte the haplotype share statistics (PHS) as reported in paper:																				
  A Nonparametric test reveals selection for rapid flowering in the arabidopsis genome, Christopher T et al, PLoS Biology 2006
  and generate plots of haplotype extent frequency for each chromosome.
 
  Input files needed:
  1. Mapped marker(SNP) distance for each chromosome
     chr_id SNP_name genetic_distance
     1A wsnp_Ex_c12345 5.8
     ...
  2. Genotype file for each accession
     accession_id snpname_1 snpname_2 ...
     acc1          A     B  B  NA...
     acc2          A     B  A  A...
  
  Output files:
  Two files will be generated: <prefix>_phs.txt and <prefix>_block_freq.txt
  
 
  Usage:
    perl $0
    <Mapped marker file>
    <Genotype file>
    <Output file prefix>
 
 *******************************************************************************************************************************
END

my ($mapped_marker_file, $genotype_file, $out_file_prefix) = @ARGV;

print_time("Reading marker file $mapped_marker_file");
my ($mapped_markers_ref, $marker_genetic_dist_ref) = read_marker_file($mapped_marker_file);

print_time("Reading genotype file $genotype_file");
my %genotype_marker = read_genotype_file($genotype_file, $mapped_markers_ref);

my %allele_freq_each_marker = count_allele_freq(\%genotype_marker);
my @accessions = keys %genotype_marker;
print STDERR "No. accessions: ", scalar @accessions, "\n";

my %haplotype_extent;
my %marker_blocks;
my %chr_block_dist;
my %acc_pairs;
print_time("Calculating haplotype extent");
foreach my $chr (keys %$mapped_markers_ref)
{
	my $count = 0;
	print_time("\tChromosome $chr");
	my @markers = @{$mapped_markers_ref->{$chr}};	
	foreach my $acc_a (0..($#accessions-1))
	{
		foreach my $acc_b ($acc_a+1 .. $#accessions)
		{
			$count++;
			print_time("Finished $count pairs ..") unless $count%1000000;
			$acc_pairs{join(":", (sort{$a cmp $b} ($accessions[$acc_a], $accessions[$acc_b])))} = 1;
			pairwise_haplotype_extent_calc
			($chr, \%genotype_marker, \%haplotype_extent, \%chr_block_dist, \%marker_blocks, \@markers, $marker_genetic_dist_ref, $accessions[$acc_a], $accessions[$acc_b]);
		}
	}
}

# calcuate PHS for each marker
print_time("Calculating PHS...");
my %PHS;
my %average_block;
my @pairs = keys %acc_pairs;
#print "No. Pairs: ", scalar @pairs,"\t", $pairs[0], "\n";
%acc_pairs = ();
foreach my $chr (keys %haplotype_extent)
{	
	my @markers = @{$mapped_markers_ref->{$chr}};
	foreach my $index (0..$#markers)
	{
		my $marker_name = $markers[$index];
		my %unique_acc;
		my @pair_with_blocks;
		foreach my $pair (@pairs)
		{
			#print $pair, "\n";
			#print join("*", (keys %{$marker_blocks{$chr}{$index}})), "\n";
			next unless exists $marker_blocks{$chr}{$index}{$pair};
			map{$unique_acc{$_}=1} (split/:/, $pair);
			push @pair_with_blocks, $pair;
		}
		next unless @pair_with_blocks > 0;
		# calculate N Z(xjx)
		my $N_sum_Z = sum_Z(\%haplotype_extent, \@pair_with_blocks, \%marker_blocks, $chr, $index, $marker_genetic_dist_ref, $mapped_markers_ref);
		my $N_sumZ_avg = $N_sum_Z/(scalar @pair_with_blocks);
		
		# calculate Z(xjx) for allele A and B;
		my ($pairs_AlleleA, $pairs_AlleleB) = group_acc_pair_according_to_allele([keys %unique_acc], $chr, $index, \%genotype_marker, $mapped_markers_ref);
		my $A_sumZ_avg;
		my $A_avg_block;
		if(@$pairs_AlleleA > 1)
		{
		my $A_sum_Z = sum_Z(\%haplotype_extent, $pairs_AlleleA, \%marker_blocks, $chr, $index, $marker_genetic_dist_ref, $mapped_markers_ref);
		$A_sumZ_avg = $A_sum_Z/(scalar @$pairs_AlleleA);
		
		$A_avg_block = calculate_average_block_left_right($pairs_AlleleA, \%marker_blocks, $chr, $index, $marker_genetic_dist_ref, $mapped_markers_ref);
		}
		else
		{
			$A_sumZ_avg = 0;
			$A_avg_block = "NA_NA";
		}
		
		my $B_sumZ_avg;
		my $B_avg_block;
		if(@$pairs_AlleleB > 1)
		{
		my $B_sum_Z = sum_Z(\%haplotype_extent, $pairs_AlleleB, \%marker_blocks, $chr, $index, $marker_genetic_dist_ref, $mapped_markers_ref);
		$B_sumZ_avg = $B_sum_Z/(scalar @$pairs_AlleleB);
		$B_avg_block = calculate_average_block_left_right($pairs_AlleleB, \%marker_blocks, $chr, $index, $marker_genetic_dist_ref, $mapped_markers_ref);	
		}
		else
		{
			$B_sumZ_avg = 0;
			$B_avg_block = "NA_NA";
		}

		$PHS{$chr}{$index} = [$A_sumZ_avg - $N_sumZ_avg, $B_sumZ_avg - $N_sumZ_avg];
		$average_block{$chr}{$index} = [$A_avg_block, $B_avg_block];
	}
}

# calcualte haplotype extent block frequency
print_time("Calculating haplotype block frequency...");
my %block_freq;
foreach my $chr (keys %haplotype_extent)
{
	my @markers = @{$mapped_markers_ref->{$chr}};		
	foreach my $pair (@pairs)
	{
		my @blocks_pairacc;
		foreach my $index (0..$#markers)
		{
			next unless exists $marker_blocks{$chr}{$index}{$pair};
			push @blocks_pairacc, $marker_blocks{$chr}{$index}{$pair};				
		}
		@blocks_pairacc = unique(@blocks_pairacc);
		map {$block_freq{$chr}{$_}++} @blocks_pairacc;
	}	
}

# Output
my $phs_out_file = $out_file_prefix . "_phs.txt";
open (PHS, ">$phs_out_file") or die $!;
print PHS join("\t", qw(Chr SNP_index SNP_name Genetic_dist PHS_A PHS_B Freq_A Freq_B Block_A Block_B)),"\n";
foreach my $chr (keys %PHS)
{
	my @markers = @{$mapped_markers_ref->{$chr}};
	foreach my $index (0..$#markers)
	{
		next unless exists $PHS{$chr}{$index};
		my $marker = $markers[$index];
		my $genetic_dist = $marker_genetic_dist_ref->{$marker};
		my @phs = @{$PHS{$chr}{$index}};
		my @avg_left_right_block = @{$average_block{$chr}{$index}};
		next if (not exists $allele_freq_each_marker{$marker});
		my @allele_freq = @{$allele_freq_each_marker{$marker}};
		print PHS join("\t",($chr, $index, $marker, $genetic_dist, @phs, @allele_freq, @avg_left_right_block)),"\n";
	}
}
close PHS;

my $block_freq_file = $out_file_prefix. "_block_freq.txt";
open (BLK, ">$block_freq_file") or die $!;
foreach my $chr (keys %block_freq)
{
	my @blocks = keys %{$block_freq{$chr}};
	@blocks = sort{(split/_/,$a)[0] <=> (split/_/,$b)[0]} @blocks;
	#map{print BLK join("\t",($chr, $_, $block_freq{$chr}{$_}/(scalar @pairs))),"\n" }@blocks;
	foreach (@blocks)
	{
		my ($start, $end) = split /_/, $_;
		#$mapped_markers_ref, $marker_genetic_dist_ref;
		my @genetic_dist = map{$marker_genetic_dist_ref->{$mapped_markers_ref->{$chr}[$_]}}($start, $end);
		unless (@genetic_dist){die $_,"\n"}
		print BLK join("\t",($chr, $_, @genetic_dist, $block_freq{$chr}{$_}/(scalar @pairs))),"\n" 
	}
}
close BLK;

print_time("Done...");

# Subroutines
sub calculate_average_block_left_right
{
	# $pairs_AlleleB, \%marker_blocks, $chr, $index, $marker_genetic_dist_ref, $mapped_markers_ref
	my ($pairs_ref, $markers_blocks_ref, $chr, $index, $marker_dist_ref, $mapped_markers_ref) = @_;
	my @avg_block = (0, 0);
	my $count = 0;
	foreach my $pair (@$pairs_ref)
	{		
		my $block = $markers_blocks_ref->{$chr}{$index}{$pair};
		next unless defined $block;
		$count++;
		my ($start, $end ) = split /_/, $block;
		unless ($index >= $start and $index <= $end){die join("\t", ($chr, $index, $block)), "\n"}
		my @genetic_dist = map{$marker_dist_ref->{$mapped_markers_ref->{$chr}[$_]}} ($start, $end);
		map{ $avg_block[$_] = ($avg_block[$_]*($count-1)+$genetic_dist[$_])/$count }(0..1);
	}
	my $index_dist = $marker_dist_ref->{$mapped_markers_ref->{$chr}[$index]};
	#unless ($index_dist >= $avg_block[0] and $index_dist <= $avg_block[1]){die join("\t", ($chr, $index, $index_dist, @avg_block)), "\n"}
	return join("_", @avg_block);
}

sub sum_Z
{
	my ($haplotype_extent, $pairs_ref, $marker_blocks_ref, $chr, $index, $marker_dist_ref, $mapped_markers_ref) = @_;
	my $no_pairs = scalar @$pairs_ref;
	my ($mean, $sd) = calculate_genome_wide_Mean_Sigma($haplotype_extent, $pairs_ref, $marker_dist_ref, $mapped_markers_ref);
	unless ($sd){print join("\t", @$pairs_ref), "\n"}
	my $sum_z = 0;
	foreach (@$pairs_ref)
	{
		my $block = $marker_blocks_ref->{$chr}->{$index}->{$_};
		next unless defined $block;
		my ($start, $end) = (split/_/, $block);
		my $dist = abs($marker_dist_ref->{$mapped_markers_ref->{$chr}->[$start]} - $marker_dist_ref->{$mapped_markers_ref->{$chr}->[$end]});
		$sum_z += ($dist - $mean)/$sd;
	}
	return $sum_z;
}


sub CN2
{
	my $N = shift;
	# calculate C(N,2);
	# C(4,2) = 6; C(5,2) = 10;
	return ($N-1)*$N/2;
}
sub group_acc_pair_according_to_allele
{
	my ($accs_ref, $chr, $index, $genotype_hashref, $mapped_marker) = @_;
	my %group_hash;
	foreach my $acc (@$accs_ref)
	{
		my $allele = $genotype_hashref->{$acc}{$mapped_marker->{$chr}->[$index]};
		next if $allele =~ /NA/;
		push @{$group_hash{$allele}}, $acc;
	}
	#if(keys %group_hash > 2){die "More than 2 alleles\n", join("\n",(keys %group_hash)),"\n"; }
	my @alleles = sort{$a cmp $b} keys %group_hash;
	$group_hash{"A"} = [] unless exists $group_hash{"A"};
	$group_hash{"B"} = [] unless exists $group_hash{"B"};
	return([pairs($group_hash{"A"})], [pairs($group_hash{"B"})]);
}

sub calculate_genome_wide_Mean_Sigma
# return (mean, stand_deviation);
{
	my ($hap_extend_freq, $accession_pairs_ref, $marker_genetic_dist_ref, $mapped_markers_ref) = @_;
	my @acc_pairs = @$accession_pairs_ref;
	my @block_length;
	foreach my $chr (keys %$hap_extend_freq)
	{
		my @markers = @{$mapped_markers_ref->{$chr}};
		foreach my $pair (@acc_pairs)
		{
			next unless defined $hap_extend_freq->{$chr}->{$pair};
			my @blocks = @{$hap_extend_freq->{$chr}->{$pair}};
			foreach (@blocks)
			{
				my @two_markers = @markers[(split/_/,$_)];
				my $dist  = abs($marker_genetic_dist_ref->{$two_markers[0]} - $marker_genetic_dist_ref->{$two_markers[1]});
				push @block_length, $dist;
			}
		}
	}
	#print join("*", @block_length), "\n"; exit;
	return(average(\@block_length), stdev(\@block_length));
}

sub pairs
{
	my $arr_ref = shift;
	return unless @$arr_ref > 1;
	my @return;
	foreach my $index (0..(scalar @$arr_ref)-2)
	{
		foreach ($index .. (scalar @$arr_ref)-1)
		{
			push @return, join(":", sort {$a cmp $b} ($arr_ref->[$index], $arr_ref->[$_]));
		}
	}
	return @return;
}


sub read_marker_file
{
	my $file = shift;
	my %markers_each_chr;
	my %distance_each_marker;
	open (IN, $file) or die;
	my $line = 0;
	while(<IN>)
	{
		chomp;
		$line++;
		next if $line==1; #/chr_id/i; # skip first line
		my ($chr, $marker, $distance) = split /\s+/, $_;
		push @{$markers_each_chr{$chr}}, $marker;
		$distance_each_marker{$marker} = $distance;
	}
	close IN;
	return(\%markers_each_chr, \%distance_each_marker);
}

sub read_genotype_file
{
	my ($genotype_file, $markers_each_chr) = @_;
	my %return;
	my @mapped_markers;
	map{push @mapped_markers, @$_} values %{$markers_each_chr};
	my %mapped_markers_hash = map{$_, 1} @mapped_markers;
	open (IN, $genotype_file) or die;
	my $line = 0;
	my @total_markers;
	my @index_for_mapped_markers;
	while(<IN>)
	{
		next if /^\s+$/;
		$line++;
		chomp;
		s/\s+$//;
		my @t = split /\s+/,$_;
		if ($line == 1) # get marker names from the first line
		{
			@total_markers = @t;
			foreach (1..$#t)
			{
				push @index_for_mapped_markers, $_ if exists $mapped_markers_hash{$t[$_]}
			}
			next;
		}
		my $acc = $t[0];
		foreach (@index_for_mapped_markers)
		{
			die $_, "\t", scalar @total_markers, "\n" unless defined $total_markers[$_];
			$return{$acc}{$total_markers[$_]} = $t[$_];
		}		
	}
	close IN;
	return %return;
}

sub count_allele_freq
{
	my $genotype_hashref = shift;
	my %return;
	my $num_accessions = scalar (keys %$genotype_hashref);
	foreach my $acc (keys %$genotype_hashref)
	{
		foreach my $mrk (keys %{$genotype_hashref->{$acc}})
		{
			my $allele = $genotype_hashref->{$acc}{$mrk};
			$return{$mrk}{$allele}++;
		}
	}
	
	foreach my $mrk (keys %return)
	{
		my @alleles = keys %{$return{$mrk}};
		@alleles = sort {$a cmp $b} @alleles;
		$return{$mrk}{"A"} = 0 unless exists $return{$mrk}{"A"};
		$return{$mrk}{"B"} = 0 unless exists $return{$mrk}{"B"};
		my $total = $return{$mrk}{"A"} + $return{$mrk}{"B"};
		if ($total > 0)
		{
			$return{$mrk} = [$return{$mrk}{"A"}/$total, $return{$mrk}{"B"}/$total];
		}
		else
		{
			$return{$mrk} = [0, 0];
		}
		
	}
	return %return;
}

sub pairwise_haplotype_extent_calc
{
	my ($chr, $gn_marker, $hap_ext_freq, $chr_block_dist, $marker_blocks, $markers, $marker_genetic_dist_ref, $acc1, $acc2) = @_;
	my $pair = join(":", sort{$a cmp $b}($acc1, $acc2));
	my @status_record;
	foreach my $index (0 ..  scalar @$markers - 1)
	{
		my $mk = $markers->[$index];
		$gn_marker->{$acc1}->{$mk} = "NA" unless exists $gn_marker->{$acc1}->{$mk};
		$gn_marker->{$acc2}->{$mk} = "NA" unless exists $gn_marker->{$acc2}->{$mk};
		if ($gn_marker->{$acc1}->{$mk} eq $gn_marker->{$acc2}->{$mk} and $gn_marker->{$acc2}->{$mk} ne "NA")
		{
			$status_record[$index] = 1;
		}
		elsif($gn_marker->{$acc1}->{$mk} eq "NA" or $gn_marker->{$acc2}->{$mk} eq "NA")
		{
			$status_record[$index] = -1;
		}
		else
		{
			$status_record[$index] = 0;
		}
		die "status_record $index not defiend\n" unless defined $status_record[$index];
	}
	my @blocks = calculate_blocks(@status_record);
	#if(@blocks){print join("\t", @blocks), "\n"; exit 1}
	#print join("\t", @blocks), "\n" if $pair eq "AUS108:AUS12";
	$hap_ext_freq->{$chr}->{$pair} = [@blocks] if @blocks > 0;
	
	foreach my $blk (@blocks)
	{
		my ($start, $end) = split /_/, $blk;
		#print $start, "\t", $end, "\n";
		foreach ($start .. $end)
		{
			$marker_blocks->{$chr}{$_}{$pair} = $blk;
			#$chr_block_dist{$chr}{$_}{$pair} = abs($marker_genetic_dist_ref->{$markers->[$start]} - $marker_genetic_dist_ref->{$markers->[$end]})
		}
	#	foreach ($start .. $end){$chr_block_dist{$chr}{$_}{$pair} = abs($marker_genetic_dist_ref->{$markers->[$start]} - $marker_genetic_dist_ref->{$markers->[$end]}) }
	}
}

sub calculate_blocks
{
	my @records = @_;
	my @return;
	my $pre; 
	foreach my $index (0..$#records)
	{
		if($records[$index] == 1)
		{
			$pre = $index unless defined $pre;
			if($index == $#records)
			{
				push @return, [$pre, $index] unless $pre == $index;
			}
		}
		elsif($records[$index] == 0)
		{
			if(defined $pre)
			{
				push @return, [$pre, $index-1] unless $pre == $index-1;
				undef $pre;
			}
		}
		elsif($records[$index] == -1)
		{
			if(defined $pre)
			{
				if(($index < $#records and $records[$index+1] != 1) or ($index == $#records))
				{
					push @return, [$pre, $index-1] unless $pre == $index-1;
					undef $pre;
				}
			}
		}
	}
	@return = map{join("_", @$_)}@return;
	#print join("", @records), "\n";
	return @return;
}

sub stdev {
	my $arr_ref = $_[0];

	my $n = 0;
	my $sumsq = 0;
	my $sum = 0;
	foreach (@$arr_ref) {
		next unless (defined($_) && ($_ !~ /^$/));
		$n++;
		$sumsq += $_ * $_;
		$sum += $_;
	}
	if ($n <= 1) { return 0; }
	my $stdev = (($n * $sumsq) - ($sum * $sum))/($n * ($n - 1));
	$stdev = ($stdev < 0) ? 0 : $stdev;

	return (sqrt($stdev));
}

sub average {
	my $count = &count($_[0]);
	if ($count == 0) {
		return 0;
	} else {
		return (&sum($_[0]) / $count);
	}
}

sub sum {
	my $arr_ref = $_[0];

	my $sum = 0;
	foreach (@$arr_ref) {
		next unless defined($_);
		next if /^$/;
		$sum += $_;
	}
	return $sum;
}

sub count {
	my $arr_ref = $_[0];
	
	my $num = 0;
	foreach (@$arr_ref) {
		next unless defined($_);
		next if /^$/;
		$num++;
	}
	return $num;
}

sub unique
{
	my %temp = map{$_, 1}@_;
	return keys %temp;
}

sub print_time
{
	my $str = shift;
	my $time = localtime(time);
	print STDERR join("  ", ($time, $str)), "\n";
}
