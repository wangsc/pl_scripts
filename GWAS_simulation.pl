#!/usr/bin/perl -w
use strict;
use warnings FATAL => "all";
use Math::Random 'random_normal';
use Statistics::R;
use Statistics::Distrib::Normal;

# Simulation is based on the method presented in paper (James Cockram et al. PNAS(2010) 21611-21616)
# Steps:

# 1. Simulate phenotypes which is controled by N loci
#    i) phenotypes were allocated a genetic component, for which allele 1 was considered positive (contribute 9/N units), whereas allele 0 made no contribution;
#   ii) environmental variation(Ve): for certain heritability( Vg/(Vg+Ve), Vg was determined by step i and ii), Ve could be calculated; simulated Ve were
#       drawn from a distribution N(0, Ve);

# 2. Remove causative markers, then select marker in certain windowsize (cM, 1, 2, 4, 5, 10)
#
# 3. Run GWAS by tessal
# SW 09.21.2011


#
sub read_genotype_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die "can't open file $file\n";
	while(<IN>)
	{
		chomp;
		my @t = split /\s+/,$_;
		my $id = shift @t;
		$return{$id} = [@t];
	}
	close IN;
	return %return;
}

sub read_genetic_distance_file
{
	my $file = shift;
	my %distance;
	my @markers;
	open (IN, $file) or die "can't open file $file \n";
	while(<IN>)
	{
		chomp;
		#next unless /^SNP/;
		my @t = split/\t/,$_;
		push @markers, $t[0];
		$distance{$t[0]} = [@t[2,1]];		
	}
	close IN;
	return (\%distance, \@markers);
}


sub simulate_phenotypes
{
	my ($gc_total, $N, $herit, $genotype_hash, $markers) = @_;
	my %phenotypes;
	my $Vg = $gc_total/$N;
	my @causative_marker_index;
	foreach my $id (keys %$genotype_hash)
	{
		#print 'genotype $id ', $id,"\n";
		my $num_markers = scalar @{$genotype_hash->{$id}};
		@causative_marker_index = generate_random_index($num_markers, $N) unless @causative_marker_index > 0;
		my $sim_Vg;
		map{$sim_Vg += $genotype_hash->{$id}->[$_] eq '1:1'?$Vg:0} @causative_marker_index;
		my $Ve = (1 - $herit)*$sim_Vg/$herit;
		my $sim_Ve = random_normal(1, 0, sqrt($Ve));
		my $sim_pheno = $sim_Vg + $sim_Ve;
		$phenotypes{$id} = $sim_pheno;
	}
	return (\%phenotypes, [@{$markers}[@causative_marker_index]])
}

#sub random_normal
#{
#	my ($mu, $sigma) = @_;
#	my $dist = new Statistics::Distrib::Normal;
#	$dist->mu($mu);
#	$dist->sigma($sigma);
#	return $dist->rand(1);
#}

sub generate_random_index
{
	my ($max, $total_want) = @_;
	my %h; 
	while (scalar (keys %h) < $total_want)
	{
		my $r = int(rand($max));
		$h{$r} = 1;
	}
	return (keys %h);
}

sub generate_input_files
{
	my ($pheno_hashref, $causa_markers, $genotypes_ref, $markers, $herit) = @_;
	my %new_genotypes;
	my @ids = keys %$pheno_hashref;
	@ids = sort{$a cmp $b} @ids;
	my $total_acc = scalar @ids;
	my $total_marker = 0;
	my @index_remained = get_remained_indexes($causa_markers, $markers);
	foreach my $id (@ids)
	{
		$new_genotypes{$id} = [@{$genotypes_ref->{$id}}[@index_remained]];
	}
	my $geno_out = "genotype" . "_H" . $herit*100; 
	open (G, ">$geno_out") or die "can't open file $geno_out\n";
	my $pheno_out = "phenotype". "_H" . $herit*100; ; 
	open (P, ">$pheno_out") or die "can't open file $pheno_out\n";
	print G $total_acc, "\t",scalar @index_remained,":2","\n";	
	print G "\t", join("\t", @{$markers}[@index_remained]),"\n";
	print P join("\t", ($total_acc, 1, 1)),"\n";
	print P "\tsim","\n";
	foreach my $id (@ids)
	{
		print P $id,"\t", $pheno_hashref->{$id},"\n";
		print G $id, "\t", join("\t", @{$new_genotypes{$id}}),"\n";
	}
	close P;
	close G;
}

sub get_remained_indexes
{
	my ($cas_markers, $all_markers) = @_;
	my %cas = map{$_, 1} @$cas_markers;
	my @remained_indexes;
	foreach (0..scalar @$all_markers - 1)
	{
		push @remained_indexes, $_ unless exists $cas{$all_markers->[$_]}
	}
	return @remained_indexes;
}

sub run_gwas
{
	my ($tassel_dir, $herit) = @_;
	my $gwas_file_prefix = "tassel_result_H". $herit*100;
	my @files = <tassel_result*>;
	my $gwas_file;
	foreach (@files){$gwas_file = $_ if /^$gwas_file_prefix/}
	my $libdir = $tassel_dir . '/'. "lib";
	my @jars = <$libdir/*.jar>;
	push @jars , $tassel_dir . '/'. 'sTASSEL.jar';
	my $CP = join(":", @jars);
	my $geno_out = "genotype" . "_H" . $herit*100; 
	my $pheno_out = "phenotype". "_H" . $herit*100; ; 
	my $parameters = '-p "'. $geno_out .'" -t "'. $pheno_out. '" -k "K_matrix_237IND.txt" -q "Q_matrix_237IND.txt" -glm -o tassel_result'. "_H". $herit*100 .' -txt';
	my $cmd = "java -classpath '$CP' -Xmx2000m net.maizegenetics.pipeline.MLMGLMFileInputPipeline $parameters";

	#print $cmd,"\n";
	unlink ($gwas_file) if defined $gwas_file;
	die $! if system($cmd); 
}

sub parse_gwas_results
{
	my ($record_hashref, $genetic_dist_hashref, $window_sizes_arrary, $causa_markers, $sim_num, $Ni, $herit) = @_;
	my $fp = 0.1 ; # 10% false positive rate
	my $gwas_file_prefix = "tassel_result_H". $herit*100;
	my @files = <*.txt>;
	my $gwas_file;
	foreach (@files){$gwas_file = $_ if /^$gwas_file_prefix/}
	#unless ($gwas_file){system("ls -lh tassel_result*.txt"); exit}
	open (IN, "$gwas_file") or die "can't open file $gwas_file\n";
	my %pvalues;
	my $n = 0;
	while(<IN>)
	{
		$n++;
		next if $n == 1;
		next if /^\s+$/;
		my @t = split /\t/,$_;
		my $snp_id = $t[2];
		my $pvalue = $t[8];
		next unless $pvalue =~ /\d+/;
		eval{$pvalue = $t[8]+0;};
		if($@){print STDERR $@, "\n", $_; exit};
		
		$pvalues{$snp_id} = $pvalue;		
	}
	close IN;
	my @ids  = sort{$pvalues{$a} <=> $pvalues{$b}} keys %pvalues;
	my $total = scalar @ids;
	# find significant snps
	my @significant_snps;
	foreach my $index (0.. $#ids)
	{
		my $rank = $index +1;
		push @significant_snps, $ids[$index] if($pvalues{$ids[$index]} <= $rank*$fp/$total);
	}
	# calculate distance between causative markers and significant snps
	#my @distances_between;
	my %distance_between;
	foreach my $snp_id (@significant_snps)
	{
		my %dist_to_causa_marker;
		foreach my $causa_marker (@$causa_markers)
		{
			next unless $genetic_dist_hashref->{$snp_id}->[0] eq $genetic_dist_hashref->{$causa_marker}->[0];
			my $d = abs($genetic_dist_hashref->{$snp_id}->[1] - $genetic_dist_hashref->{$causa_marker}->[1]);
			$dist_to_causa_marker{$causa_marker} = $d;
			#push @distances_between, $d;
			#push @{$distance_between{$causa_marker}}, $d; 
		}
		my $detected_causa_marker = (sort{$dist_to_causa_marker{$a}<=>$dist_to_causa_marker{$b}} keys %dist_to_causa_marker)[0];
		next unless defined $detected_causa_marker;
		push @{$distance_between{$detected_causa_marker}}, $dist_to_causa_marker{$detected_causa_marker};
		#print STDERR '$detected_causa_marker ', $detected_causa_marker, "\t", $dist_to_causa_marker{$detected_causa_marker}, "\n"; exit;
	}
	my %check_each_window; map{$check_each_window{$_}=0}@$window_sizes_arrary;
	foreach my $win_size (@$window_sizes_arrary)
	{
		# map{$check_each_window{$win_size}++ if $_<= $win_size} @distances_between
		foreach (keys %distance_between)
		{
			my @d = @{$distance_between{$_}};
			foreach (@d){ if($_<= $win_size){$check_each_window{$win_size}++; last} }
		}
	}
	map{ $check_each_window{$_} /= $Ni } keys %check_each_window;
	map{ $record_hashref->{$_} = ($record_hashref->{$_}* ($sim_num - 1) + $check_each_window{$_})/$sim_num } @$window_sizes_arrary;
}

sub print_timestamp
{
	my $comment = shift;
	my $now = localtime time;
	print STDERR $now,"\t", $comment,"\n";
}

# Initial
my @N_units = (1, 2, 4, 6, 8, 10);
my $gc_total = 9;
my @window_sizes = (1..10);
#my $simulation_times = 100;
my $R = Statistics::R->new();
# Main
my $sUsage = qq(
perl $0 
<marker genetic distance file>
<genotype file>
<heritability, <=1 >
<tassel directory>
<number of simulations, 100 or ...>
);
die $sUsage unless @ARGV >= 5;

my ($mark_distance_file, $genotype_file, $herit, $tassel_directory, $simulation_times) = @ARGV;
my ($marker_genetic_distance, $markers) = read_genetic_distance_file($mark_distance_file);
my %genotypes = read_genotype_file($genotype_file);

foreach my $Ni (@N_units)
{
	print_timestamp("Simulation trait controled by $Ni loci\n");
	my %average_proportion_identified;
	map{$average_proportion_identified{$_}=0}@window_sizes;
	my $sim_num = 1;
	while( $sim_num <= $simulation_times)
	{
		my ($sim_pheno_hashref, $causa_marker_arrarref) = simulate_phenotypes($gc_total, $Ni, $herit, \%genotypes, $markers);
		generate_input_files($sim_pheno_hashref, $causa_marker_arrarref, \%genotypes, $markers, $herit);
		run_gwas($tassel_directory, $herit);
		parse_gwas_results(\%average_proportion_identified, $marker_genetic_distance, \@window_sizes, $causa_marker_arrarref, $sim_num, $Ni, $herit);
		
		print "Number simulation: $sim_num\n";
		my @proportions;
		print join ("\t", @window_sizes),"\n";
		map{push @proportions, $average_proportion_identified{$_}} @window_sizes;
		print join("\t", @proportions),"\n";
		
		$sim_num++;
		#$sim_num  = $simulation_times+1; # for debug
	}
	print_timestamp("  Finished $simulation_times simulations ...\n\n");
	print ">Simulation trait controled by $Ni loci\n";
	print join ("\t", @window_sizes),"\n";
	my @proportions;
	map{push @proportions, $average_proportion_identified{$_}} @window_sizes;
	print join("\t", @proportions),"\n";
	print "="x30, "\n";
}



