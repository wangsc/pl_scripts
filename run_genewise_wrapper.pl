#!usr/bin/perl -w
use strict;
my $sUsage = qq(
# run genewise wraper
perl $0 
<ortholog list file>
<pasa_chr3A.pasa_assemblies.gff3>
<rice or brachypodium protein seqences in fasta>
<3A_contigs_scaffold.fa>
<genewise gff output file name>
<output file>
);
die $sUsage unless @ARGV >= 6;
my ($orth_file, $ass_gff, $pep_fasta, $genome_fasta, $genewise_gff, $out_file) = @ARGV;

my %transcript_orth_pair = read_orthology_file($orth_file);
my ($transcript_genomic_pair, $exon_hash) = read_ass_gff($ass_gff);
my %pep = read_fasta_file($pep_fasta);
my %genome = read_fasta_file($genome_fasta);
my %predict_exons;
if (-e "$genewise_gff")
{
	%predict_exons = parse_genewise_all($genewise_gff);
}
else
{
	my %record;
	foreach my $id (keys %transcript_orth_pair)
	{
		my $pep_id = $transcript_orth_pair{$id};
		$pep_id =~ s/^(.*)\.\d+/$1/;
		my $geno_id = $transcript_genomic_pair->{$id};
		next unless exists $genome{$geno_id};
		foreach my $fasta_id (keys %pep)
		{
			next if exists $record{$fasta_id};
			next unless $fasta_id =~ /$pep_id/;
			$record{$fasta_id} = 1;
			#print $fasta_id ,"\n";
			write_pep_genomic($pep{$fasta_id}, $genome{$geno_id}, $id);
			run_genewise($id);
			#system("echo \"ID:$id\" >>genewise_all.out & cat genewise.out >>genewise_all.out");
			my $ind = $1 if $orth_file=~/_(\d+)$/;
			system("cat tmp$ind.gff >>$genewise_gff");
			parse_genewise_result(\%predict_exons);				
		}

	}
}

my %exon_overlapping = calc_exon_overlaping(\%predict_exons,$exon_hash);
output(\%exon_overlapping, $out_file);


########################
sub parse_genewise_all
{
	my $file = shift;
	open (IN, "$file") or die $!;
	my $genewise_score_cutoff = 35; # added 06.27.2011 SW
	my %return_hash;
	my $id;
	my $debug = 1;
	my $score;
	while(<IN>)
	{
		chomp;
		next if /^\s+$/;
		next if /^\//;
		my @t = split /\t/, $_;
		if ($t[2] =~ /match/){$score = $t[5]}
		next unless $t[2] =~ /cds/i;
		next unless $score > $genewise_score_cutoff;
		my ($id, $start, $end) = (@t[0,3,4]);
		push @{$return_hash{$id}}, [$start, $end];		
	}
	close IN;	
	return %return_hash;
}

sub read_orthology_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file\n";
	my $debug = 1;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my $line = $_;
		my @data = split /\s+/, $line;
		$data[1] =~ s/(.*?)\|.*/$1/;
		print 'Orth pair: ', join("\t", @data),"\n" if $debug==1; $debug =0;
		$return_hash{$data[0]} = $data[1];
	}
	close IN;
	return %return_hash;
}

sub read_ass_gff
{
	my $file = shift;
	my %return_hash;
	my %exon_hash;
	open (IN, "$file") or die "can't open file $file\n";
	my $debug = 1;
	while (<IN>)
	{
		next if /^\s+$/;
		my $id = $1 if /Target=(.*?) /;
		my @data = split /\s+/, $_;
		my $geno_id = $data[0];
		$return_hash{$id} = $geno_id;
		print join("\t",($id, $geno_id)),"\n" if $debug ==1; $debug =0;
		push @{$exon_hash{$id}}, [ @data[3,4] ];
	}
	close IN;
	return (\%return_hash, \%exon_hash);
}

sub read_fasta_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file\n";
	my $id;
	my $debug =1;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		if (/^>(\S+)/)
		{
			$id = $1;
			if (/mitochondrion/i){undef $id; next}
			$id =~ s/(.*?)\|.*/$1/ if $id =~ /LOC_/;
			my @tmp = split /\s+/, $id;
			$id = shift @tmp;
			$id =~ s/^>//;
			print 'fasta ID: ', $id,"\n" if $debug ==1; $debug=0;
			next;
		}
		$return_hash{$id} .= $_ if defined $id;
	}
	return %return_hash;
}

sub write_pep_genomic
{
	my ($pep, $genomic, $id) = @_;
	my $ind = $1 if $orth_file=~/_(\d+)$/;
	open (P, ">pep".$ind) or die $!;
	print P '>pep',"\n",$pep;
	close P;
	open (G, ">genomic".$ind) or die $!;
	print G ">$id","\n",$genomic,"\n";
	close G;
}

sub run_genewise
{
	my $ind = $1 if $orth_file=~/_(\d+)$/;
	my $command = "genewise -gff pep$ind genomic$ind -splice flat -both >tmp$ind.gff";
	#print $command,"\n";
	die "genewise failed \n" if system($command);
}

sub parse_genewise_result
{
	my ($predict_exons) = @_;
	my $ind = $1 if $orth_file=~/_(\d+)$/;
	open (IN, "tmp$ind.gff") or die $!;
	my $genewise_score_cutoff = 35; # added 06.27.2011 SW
	my $score;
	while(<IN>)
	{
		next if /^\s+$/;
		next if /^\//;
		my @t = split /\t/, $_;
		if ($t[2] =~ /match/){$score = $t[5]}
		next unless $t[2] =~ /cds/i;
		next unless $score > $genewise_score_cutoff;
		my $id = $t[0];
		my ($start, $end) = @t[3, 4];
		($start, $end) = ($end, $start) if $start > $end;
		push @{$predict_exons->{$id}}, [$start, $end];		
	}
	close IN;
}

sub calc_exon_overlaping
{
	my ($predic, $orig) = @_;
	my %return_hash;
	foreach my $id (keys %$predic)
	{
		next unless exists $orig->{$id};
		my ($p_vec, $p_max, $p_len) = construct_vec($predic->{$id});
		my ($o_vec, $o_max, $o_len) = construct_vec($orig->{$id});
		my $overlap_len = cal_vec_overlap($p_vec, $o_vec, $p_max>$o_max?$p_max:$o_max);
		$return_hash{$id} = [$o_len, $p_len, $overlap_len];
	}
	return %return_hash;
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
			$total++;
			$max = $_ unless defined $max;
			$max = $_ if $_ > $max;
			vec($vec,$_,1) = 0b1;
		}
	}
	return ($vec, $max, $total);
}

sub cal_vec_overlap
{
	my($pvec, $ovec, $len) = @_;
	my $count = 0;
	foreach my $pos ( 1..$len)
	{
		$count++ if vec($pvec, $pos,1)==1 and vec($ovec, $pos,1)==1;
	}
	return $count;
}

sub output
{
	my ($hashref, $file) = @_;
	open (OUT, ">$file") or die "can't open file $file\n";
	print OUT join("\t", qw(Asmbl_ID Wheat_exon Predicted_exon Overlapping) ),"\n";
	foreach my $id (keys %$hashref)
	{
		print OUT join("\t",($id, @{$hashref->{$id}})),"\n";
	}
	close OUT;
}

sub print_time
{
	my $str = shift;
	my $time = localtime(time);
	print STDERR join("  ", ($time, $str)), "\n";
}


