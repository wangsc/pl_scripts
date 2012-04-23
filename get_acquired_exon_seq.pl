#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<pasa assembly gff3>
<pasa transcripts isoform gff3>
<contig fasta>
<rice genewise outfile>
<result fasta file>
);
die $sUsage unless @ARGV >= 5;
my ($pasa_gff, $isoform_gff, $contig_fasta_file, $rice_wise, $output) = @ARGV;

my ($iso_contig, %pasa_structure) = read_pasa_gff($pasa_gff);
my %wise_structure = read_wise_gff($rice_wise);
my ($iso_gene) = read_isoform_gff($isoform_gff);
my %gene_iso;
foreach (keys %$iso_gene)
{
	push @{ $gene_iso{$iso_gene->{$_}} }, $_;
}
my %contig_fasta = read_fasta($contig_fasta_file);

my %acquired_exons = calculate_acquired_exons(\%pasa_structure, \%wise_structure, \%gene_iso);

output(\%acquired_exons, \%contig_fasta, $iso_contig, $output);



# Subroutines

sub output
{
	my ($exon_ref, $fasta_ref, $iso_contig, $outfile) = @_;
	open (OUT, ">$outfile") or die "$! $outfile\n";
	my %fasta_ids;
	foreach my $gene (keys %{$exon_ref})
	{
		foreach my $id (keys %{$exon_ref->{$gene}})
		{
			my $ctg = $iso_contig->{$id};
			foreach (@{$exon_ref->{$gene}{$id}})
			{
				my ($start, $end) = @$_;
				($start, $end) = ($end, $start) if $start > $end;
				die 'No value ' unless @$_;
				die '$id, $ctg, $start, $end ', join("\t", ($id, $ctg, $start, $end) ),"\n" if $end > (length $fasta_ref->{$ctg});
				my $seq = substr($fasta_ref->{$ctg}, $start-1, $end-$start+1);
				my $seqid = join('_', ($gene, $start, $end));
				print OUT '>', $seqid,"\n",$seq,"\n" unless exists $fasta_ids{$seqid};
				$fasta_ids{$seqid} = 1;
			}			
		}
	}
	close OUT;
}

sub read_fasta
{
	my $file = shift;
	open (IN, "$file") or die "$! $file\n";
	my %return_hash;
	my $id;
	my $debug = 1;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		if (/^>(\S+)/)
		{
			$id = $1;
			print 'fasta id: ', $id ,"\n" if $debug; $debug = 0;
			$return_hash{$id} = '';
			next;
		}
		$return_hash{$id} .= $_;
	}
	close IN;
	return %return_hash;
}

sub calculate_acquired_exons
{
	my ($pasa_ref, $wise_ref, $geneiso_ref) = @_;
	my %acquired_exons;
	my $debug = 1;
	foreach my $gene (keys %$geneiso_ref)
	{
		my @wise_vecs;
		foreach my $iso (@{$geneiso_ref->{$gene}})
		{
			#print 'ref $wise_ref->{$iso} ', ref $wise_ref->{$iso},"\n";
			next unless exists $wise_ref->{$iso};
			next unless  @{$wise_ref->{$iso}};
			print '$iso ', $iso,"\n" if $debug; $debug = 0;
			my ($wise_vec) = construct_vec($wise_ref->{$iso});
			push @wise_vecs, $wise_vec;
		}
		next unless @wise_vecs;
		#print '@wise_vecs ', scalar @wise_vecs,"\n";
		foreach my $iso (@{$geneiso_ref->{$gene}})
		{
			foreach (@{$pasa_ref->{$iso}})
			{
				my ($start, $end) = @$_;
				print '$iso ', $iso,"\n" if $start == 4191;
			#	die join("\t", ($start, $end)),"\n" if $iso =~ /asmbl_1081/;
				my $covered_flag = check_covered([$start, $end], @wise_vecs);
				push @{$acquired_exons{$gene}{$iso}}, [$start, $end] unless $covered_flag;
				unless ($covered_flag) {print join("\t",($gene, $start, $end)),"\n" if $gene eq 'S5719'}
			}
		}		
	}
	return %acquired_exons;
}


sub check_covered
{
	my ($arrayref, @vecs) = @_;
	my ($start, $end) = @$arrayref;
	foreach my $vec (@vecs)
	{
		foreach ($start .. $end)
		{
			return 1 if (vec($vec, $_, 1) ==1)
		}
	}
	return 0;
}


sub read_isoform_gff
{
	my $file = shift;
	my %iso_gene;
	my %iso_contig;
	my $debug = 1;
	open (IN, "$file") or die "$! $file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my ($gene, $iso_id) = $_=~/ID=(.*?)\-(.*?)\;/;
		$iso_gene{$iso_id} = $gene;
		my @data = split /\t/, $_;
		my $contig_id = $data[0];
		print 'contig id: ', $contig_id,"\n" if $debug; $debug = 0;
		$iso_contig{$iso_id} = $contig_id;
	}
	close IN;
	return (\%iso_gene, \%iso_contig, \%pasa_structure);
}

sub read_wise_gff
{
	my $file = shift;
	my %return_hash;
	open (IN, $file) or die $!;
	my $score;
	my $genewise_cutoff = 0; # 35;
	my @array;
	my $iso_id;
	my $debug = 1;
	while (<IN>)
	{
		next if /^\s+$/;
		if (/^\/\//)
		{
			next unless @array;
			@array = sort{$a->[0] <=> $b->[0]} @array;
			push @{$return_hash{$iso_id}}, [@array] if @array > 1;# at least two exons
			print 'wise iso_id ',$iso_id,"\n" if $debug; $debug =0;
			@array = ();
			next;
		}
		my @t = split /\t/, $_;
		if ($t[2] =~ /match/){$score = $t[5]}
		next unless $t[2] =~ /cds/i;
		next unless $score > $genewise_cutoff;
		$iso_id = $t[0];
		print join("\t", @t[3,4]),"\n" if $iso_id =~ /asmbl_1081/;
		push @array,[sort{$a<=>$b}@t[3, 4]];		
	}
	close IN;
	return %return_hash;
}

sub read_pasa_gff
{
	my $file = shift;
	my %return_hash;
	my %iso_contig;
	open (IN, $file) or die $!;
	while (<IN>)
	{
		next if /^\s+$/;
		my $id = $1 if /Target=(\S+)\s/;
		my @t = split /\t/, $_;
	#	print $id, "\t", join("\t", @t[3, 4]),"\n" if $id =~ /asmbl_1222/;
		push @{$return_hash{$id}},[@t[3, 4]];
		$iso_contig{$id} = $t[0];
	}
	close IN;
	map{ $return_hash{$_} = [ sort{$a->[0]<=>$b->[0]} @{$return_hash{$_}} ] }keys %return_hash;
	return (\%iso_contig, %return_hash);
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
	my ($wvec, $pvec, $max, $total) = @_;
	my $n = 0;
	foreach (1..$max)
	{
		$n++ if (vec($wvec, $_, 1)==1 and vec($pvec, $_, 1)==1);
	}
	return $n/$total;
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
		my @d = sort{$a->[0]<=>$b->[0]}@$_;
#		print '@d: ', join("\t", @d),"\n" if $debug; $debug=0;
		foreach (@d)
		{
			foreach ($_->[0]..$_->[1])
			{
				vec($vec,$_,1) = 0b1;
				$total++;
				$max = $_ unless defined $max;
				$max = $_ if $_ > $max;
				
			}			
		}

	}
	return ($vec, $max, $total);
}
