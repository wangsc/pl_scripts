#!/usr/bin/perl -w
use strict;

# By SWANG on May 9 2012
# Fast version

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
my ($mapped_markers_ref, $marker_genetic_dist_ref, $chr) = read_marker_file($mapped_marker_file);

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

my $count = 0;
print_time("\tChromosome $chr");
my @markers = @$mapped_markers_ref;	
foreach my $acc_a (0..($#accessions-1))
{
	foreach my $acc_b ($acc_a+1 .. $#accessions)
	{
		$count++;
		print_time("Finished $count pairs ..") unless $count%1000000;
		my $pair = join(":", (sort{$a cmp $b} ($accessions[$acc_a], $accessions[$acc_b])));
		$acc_pairs{$pair} = 1;
		pairwise_haplotype_extent_calc
		(\%genotype_marker, \%haplotype_extent, \%chr_block_dist, \%marker_blocks, \@markers, $marker_genetic_dist_ref, $accessions[$acc_a], $accessions[$acc_b], $pair);
	}
}


# calcuate PHS for each marker
print_time("Calculating PHS...");
my %PHS;
my %average_block;
my @pairs = keys %acc_pairs;
#print "No. Pairs: ", scalar @pairs,"\t", $pairs[0], "\n";
%acc_pairs = ();

my %mean_sd_hash = calculate_genome_wide_Mean_Sigma(\%haplotype_extent, \@pairs, $marker_genetic_dist_ref, $mapped_markers_ref);
foreach my $index (0..$#markers)
{
	my $marker_name = $markers[$index];
	
	#my ($mean, $sd) = calculate_genome_wide_Mean_Sigma(\%haplotype_extent, \@pair_with_blocks, $marker_genetic_dist_ref, $mapped_markers_ref);
  
	# calculate N Z(xjx)
	my ($N_sumZ_avg, $A_sumZ_avg, $B_sumZ_avg, $A_avg_block, $B_avg_block) = sum_Z($marker_name, \%genotype_marker, \%haplotype_extent, \@pairs, \%marker_blocks, $index, $marker_genetic_dist_ref, $mapped_markers_ref, \%mean_sd_hash);
	
	$PHS{$index} = [$A_sumZ_avg - $N_sumZ_avg, $B_sumZ_avg - $N_sumZ_avg];
	$average_block{$index} = [$A_avg_block, $B_avg_block];
}

=head
# calcualte haplotype extent block frequency
#print_time("Calculating haplotype block frequency...");
#my %block_freq;
#foreach my $chr (keys %haplotype_extent)
#{
#	my @markers = @{$mapped_markers_ref->{$chr}};		
#	foreach my $pair (@pairs)
#	{
#		my @blocks_pairacc;
#		foreach my $index (0..$#markers)
#		{
#			next unless exists $marker_blocks{$chr}{$index}{$pair};
#			push @blocks_pairacc, $marker_blocks{$chr}{$index}{$pair};				
#		}
#		@blocks_pairacc = unique(@blocks_pairacc);
#		map {$block_freq{$chr}{$_}++} @blocks_pairacc;
#	}	
#}
=cut
# Output
my $phs_out_file = $out_file_prefix . "_phs.txt";
open (PHS, ">$phs_out_file") or die $!;
print PHS join("\t", qw(Chr SNP_index SNP_name Genetic_dist PHS_A PHS_B Freq_A Freq_B Block_A Block_B)),"\n";

foreach my $index (0..$#markers)
{
	next unless exists $PHS{$index};
	my $marker = $markers[$index];
	my $genetic_dist = $marker_genetic_dist_ref->{$marker};
	my @phs = @{$PHS{$index}};
	next if $phs[0] == 0 and $phs[1] == 0;
	my @avg_left_right_block = @{$average_block{$index}};
	next if $avg_left_right_block[0] eq "0_0" or $avg_left_right_block[1] eq "0_0";
	next if (not exists $allele_freq_each_marker{$marker});
	my @allele_freq = @{$allele_freq_each_marker{$marker}};
	print PHS join("\t",($chr, $index, $marker, $genetic_dist, @phs, @allele_freq, @avg_left_right_block)),"\n";
}

close PHS;

=head
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
=cut
print_time("Done...");

# Subroutines

sub sum_Z
{
	my ($marker_name, $genotype_ref, $haplotype_extent, $pairs_ref, $marker_blocks_ref, $index, $marker_dist_ref, $mapped_markers_ref, $mean_sd_hashref) = @_;
	my $no_pairs = scalar @$pairs_ref;
	
	my @sum_z = (0,0,0);
	my @avg_left_right_block_A = (0,0);
	my @avg_left_right_block_B = (0,0);
	my @count = (0,0,0);
	foreach my $pair (@$pairs_ref)
	{
		next unless exists $marker_blocks{$index}{$pair};
		next unless exists $mean_sd_hashref->{$pair};
		my ($mean, $sd) = @{$mean_sd_hashref->{$pair}};
		unless (defined $sd){print STDERR $pair, "\n"}
		
		my ($acc1, $acc2) = split /:/, $pair;
		my $genotype1 = $genotype_ref->{$acc1}{$marker_name};
		my $genotype2 = $genotype_ref->{$acc2}{$marker_name};
		next if $genotype1 eq "NA" or $genotype2 eq "NA";
		if($genotype1 ne $genotype2){print STDERR "Err, strange ...\n"; next}
	
		my $block = $marker_blocks_ref->{$index}->{$pair};
		next unless defined $block;
		my ($start, $end) = (split/_/, $block);
		my ($left_dist, $right_dist) = sort{$a<=>$b} ($marker_dist_ref->{$mapped_markers_ref->[$start]}, $marker_dist_ref->{$mapped_markers_ref->[$end]});
		my $dist = $right_dist - $left_dist;
		$sum_z[0] += $sd>0?($dist - $mean)/$sd:0;
		$count[0]++;
		if($genotype1 eq "A"){$count[1]++; $sum_z[1] += $sd>0?($dist - $mean)/$sd:0; $avg_left_right_block_A[0] += $left_dist; $avg_left_right_block_A[1] += $right_dist}
		if($genotype1 eq "B"){$count[2]++; $sum_z[2] += $sd>0?($dist - $mean)/$sd:0; $avg_left_right_block_B[0] += $left_dist; $avg_left_right_block_B[1] += $right_dist}
		
	}
	
	map{$sum_z[$_] =  $sum_z[$_]/$count[$_] if $count[$_]>0}0..$#sum_z;
	map{$avg_left_right_block_A[$_] = $avg_left_right_block_A[$_]/$count[1] if $count[1]>0} 0 .. $#avg_left_right_block_A;
	map{$avg_left_right_block_B[$_] = $avg_left_right_block_B[$_]/$count[2] if $count[2]>0} 0 .. $#avg_left_right_block_B;
	return (@sum_z, join("_", @avg_left_right_block_A), join("_", @avg_left_right_block_B));
}


sub CN2
{
	my $N = shift;
	# calculate C(N,2);
	# C(4,2) = 6; C(5,2) = 10;
	return ($N-1)*$N/2;
}


sub calculate_genome_wide_Mean_Sigma
# return (mean, stand_deviation);
{
	my ($hap_extend_freq, $accession_pairs_ref, $marker_genetic_dist_ref, $mapped_markers_ref) = @_;
	my @acc_pairs = @$accession_pairs_ref;
	my %return;

	my @markers = @{$mapped_markers_ref};
	foreach my $pair (@acc_pairs)
	{
		my @block_length;
		next unless defined $hap_extend_freq->{$pair};
		my @blocks = @{$hap_extend_freq->{$pair}};
		foreach (@blocks)
		{
			my @two_markers = @markers[(split/_/,$_)];
			my $dist  = abs($marker_genetic_dist_ref->{$two_markers[0]} - $marker_genetic_dist_ref->{$two_markers[1]});
			push @block_length, $dist;
		}
		$return{$pair} = [average(\@block_length), stdev(\@block_length)];
	}
	
	#print join("*", @block_length), "\n"; exit;
	#return(average(\@block_length), stdev(\@block_length));
	return %return;
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
	my @markers_each_chr;
	my %distance_each_marker;
	my $chromosome;
	open (IN, $file) or die;
	my $line = 0;
	while(<IN>)
	{
		chomp;
		$line++;
		next if $line==1; #/chr_id/i; # skip first line
		my ($chr, $marker, $distance) = split /\s+/, $_;
		$chromosome = $chr;
		push @markers_each_chr, $marker;
		$distance_each_marker{$marker} = $distance;
	}
	close IN;
	return(\@markers_each_chr, \%distance_each_marker, $chromosome);
}

sub read_genotype_file
{
	my ($genotype_file, $markers_each_chr) = @_;
	my %return;
	my @mapped_markers = @$markers_each_chr;
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
	my ($gn_marker, $hap_ext_freq, $chr_block_dist, $marker_blocks, $markers, $marker_genetic_dist_ref, $acc1, $acc2, $pair) = @_;
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
	$hap_ext_freq->{$pair} = [@blocks] if @blocks > 0;
	
	foreach my $blk (@blocks)
	{
		my ($start, $end) = split /_/, $blk;
		#print $start, "\t", $end, "\n";
		foreach ($start .. $end)
		{
			$marker_blocks->{$_}{$pair} = $blk;
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
