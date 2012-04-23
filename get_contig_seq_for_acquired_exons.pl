#!/usr/bin/perl -w
use strict;
my $sUsage = qq(
perl $0
<acquired exon fasta>
<pasa assembly gff3>
<pasa isoform gff3>
<3A contig fasta>
);
die $sUsage unless @ARGV >= 3;
my ($ac_fasta, $asmbl_gff, $iso_gff, $ctg_fasta) = @ARGV;

my %ac_genes = read_ac_fasta($ac_fasta);
my %iso_gene = rad_isoform_gff($iso_gff);
my ($iso_ctg, $pasa_exon, $region_ref) = read_pasa_gff($asmbl_gff);
my ($gene_ctg, $contig_regions, $exon_regions) = get_contig_name($asmbl_gff, \%ac_genes, \%iso_gene, $region_ref);
my %ctg_fasta = read_fasta($ctg_fasta);
foreach my $gene (keys %{$gene_ctg})
{
	my $ctg = $gene_ctg->{$gene};
	my @ctg_regs = sort{$a<=>$b}@{$contig_regions->{$gene}};
	my ($start, $end) = @ctg_regs[0, -1];
	my $seq = substr($ctg_fasta{$ctg}, $start-1, $end-$start+1);
	print '>', join("_", ($gene, @{$exon_regions->{$gene}})),"\n",$seq,"\n";
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
			print STDERR 'fasta id: ', $id ,"\n" if $debug; $debug = 0;
			$return_hash{$id} = '';
			next;
		}
		$return_hash{$id} .= $_;
	}
	close IN;
	return %return_hash;
}

sub read_ac_fasta
{
	my $file = shift;
	open(IN, $file) or die;
	my %return;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		next unless /^>/;
		s/>//;
		my @data = split /_/, $_;
		push @{$return{$data[0]}}, join('_', @data[1,2]);
	}
	close IN;
	return %return;
}

sub get_contig_name
{
	my $asmbl_gff = shift;
	my $ac_gene_ref = shift;
	my $iso_gene = shift;
	my $region = shift;
	my %return;
	open(IN, $asmbl_gff) or die;
	my %contig_regions;
	my %exon_regions;
	my %gene_ctg;
	while (<IN>)
	{
		next if /^\s+$/;
		my $line = $_;
		my @data = split /\t/, $line;
		my $asmbl = $1 if $line =~ /Target=(\S+)/;
		my $gene = $iso_gene->{$asmbl};
		next unless exists $ac_gene_ref->{$gene};
		$gene_ctg{$gene} = $data[0];
		push @{$contig_regions{$gene}}, (@data[3,4]);
		my $region = join('_', @data[3,4]);
		my $flag = 0;
		foreach ( @{$ac_gene_ref->{$gene}})
		{
			if($region eq $_)
			{
				$flag = 1;
				print STDERR $region,"\t",$_,"\n";
				my ($s, $e) = $line =~ /Target=\S+\s(\d+)\s+(\d+)/;
				print STDERR '$s, $e ', $s, "\t", $e, "\n";
				push @{$exon_regions{$gene}}, join('_', ($s, $e))
			}
		}
		print STDERR 'Gene: ', $gene,"\n" unless $flag;
	}	
	close IN;
	return(\%gene_ctg, \%contig_regions, \%exon_regions);
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
	return %iso_gene;
}

sub read_pasa_gff
{
	my $file = shift;
	my %return_hash;
	my %iso_contig;
	my %regions;
	open (IN, $file) or die $!;
	while (<IN>)
	{
		next if /^\s+$/;
		my $id = $1 if /Target=(\S+)\s/;
		my @t = split /\t/, $_;
	#	print $id, "\t", join("\t", @t[3, 4]),"\n" if $id =~ /asmbl_1222/;
		push @{$return_hash{$id}},[@t[3, 4]];
		$iso_contig{$id} = $t[0];
		my ($s, $e) = /Target=\S+\s(\d+)\s+(\d+)/;
		$regions{$id}{join('_',(@t[3,4]))} = join('_', ($s, $e));
	}
	close IN;
	map{ $return_hash{$_} = [ sort{$a->[0]<=>$b->[0]} @{$return_hash{$_}} ] }keys %return_hash;
	return (\%iso_contig, %return_hash,\%regions);
}



