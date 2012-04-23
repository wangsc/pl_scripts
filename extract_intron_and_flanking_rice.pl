#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<rice gff>
<rice gene fasta file>
<flanking lenght, 30>
<output file>
);
die $sUsage unless @ARGV >= 3;
my ($gff, $fasta, $fl, $outfile) = @ARGV;

my ($strand, $contigs, $gene_starts, $gene_pos)= read_gff($gff);
my %fasta_seq = read_fasta($fasta);

open (OUT ,">$outfile") or die $!;
foreach my $gene(keys %$gene_pos)
{
	next unless exists $fasta_seq{$gene};

	foreach my $asmbl (keys %{$gene_pos->{$gene}})
	{
		my @introns;
	  my @ctg_pos;
		next unless @{$gene_pos->{$gene}{$asmbl}} >1;
		my $arrref = $gene_pos->{$gene}{$asmbl};
		foreach my $index ( 0..(scalar @$arrref -2))
		{
			my @pos = ($arrref->[$index]->[0], $arrref->[$index]->[1], $arrref->[$index+1]->[0], $arrref->[$index+1]->[1]);
			@pos = sort{$a<=>$b} @pos;
			push @ctg_pos, [$pos[1]+1, $pos[2]-1];
			push @introns, 
			join("_",($arrref->[$index]->[0]+$gene_starts->{$gene}-1, $arrref->[$index]->[1]+$gene_starts->{$gene}-1, 
			     $arrref->[$index+1]->[0]+$gene_starts->{$gene}-1, $arrref->[$index+1]->[1]+$gene_starts->{$gene}-1))			
		}
		foreach my $ind (0..$#ctg_pos)
		{
			my ($start, $end) = @{$ctg_pos[$ind]};
			($start, $end) = ($end, $start) if $start > $end;
			my $seq = substr($fasta_seq{$gene}, $start-1-$fl, ($end-$start+1+$fl*2));
			print OUT '>', join("_",($gene, $asmbl, $introns[$ind])),"\n", $seq, "\n";
		}
	}
}
close OUT;



sub read_gff
{
	my $file = shift;
	open (IN, $gff) or die;
	my %strand;
	my %contigs;
	my %gene_starts;
	my %return;
	my ($gene_id, $iso_id, $gene_start);
	my $debug = 1;
	while (<IN>)
	{
		next unless /^Chr\d/;
		chomp;
		my @t = split /\t/,$_;
		if ($t[2] =~ /mRNA/)
		{
			($gene_id, $iso_id ) = ($1, $2) if /Alias=(.*?)\.(\d)/;
			print '$gene_id: ', $gene_id,"*\n" if $debug; $debug=0;
			$gene_start = $t[3];
			$gene_starts{$gene_id} = $gene_start;
			next;
		}
		next unless $t[2] =~ /cds/i;
		$strand{$gene_id} = $t[6] unless exists $strand{$gene_id};
		$contigs{$gene_id} = $gene_id;
		my @pos = map{$_-$gene_start+1}@t[3,4];
		push @{$return{$gene_id}{$iso_id}}, [@pos];
	}
	close IN;
	return (\%strand, \%contigs, \%gene_starts, \%return)
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
		if (/^>(\S+?)\|/)
		{
			$id = $1;
			print 'fasta id: ', $id ,"*\n" if $debug; $debug = 0;
			$return_hash{$id} = '';
			next;
		}
		$return_hash{$id} .= $_;
	}
	close IN;
	return %return_hash;
}