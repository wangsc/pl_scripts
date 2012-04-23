#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<pasa transcript isoform gff>
<ctg fasta file>
<flanking lenght, 30>
<output file>
);
die $sUsage unless @ARGV >= 3;
my ($gff, $fasta, $fl, $outfile) = @ARGV;

my ($strand, $contigs, $gene_pos)= read_gff($gff);
my %fasta_seq = read_fasta($fasta);

open (OUT ,">$outfile") or die $!;
foreach my $gene(keys %$gene_pos)
{
	foreach my $asmbl (keys %{$gene_pos->{$gene}})
	{
		my @introns;
		my @ctg_pos;
		next unless @{$gene_pos->{$gene}{$asmbl}} >1;
		my $arrref = $gene_pos->{$gene}{$asmbl};
		foreach my $index ( 0..(scalar @$arrref -2))
		{
			my @pos = ($arrref->[$index]->[2], $arrref->[$index]->[3], $arrref->[$index+1]->[2], $arrref->[$index+1]->[3]);
			@pos = sort{$a<=>$b} @pos;
			push @ctg_pos, [$pos[1]+1, $pos[2]-1];
			push @introns, join("_",($arrref->[$index]->[0], $arrref->[$index]->[1], $arrref->[$index+1]->[0], $arrref->[$index+1]->[1]));
			print STDERR $introns[$index],"\n" if $asmbl eq 'asmbl_7587';
			print STDERR join("\t",($arrref->[$index]->[0], $arrref->[$index]->[1])),"\n" if $asmbl eq 'asmbl_7587';
			print STDERR $gene_pos->{$gene}{$asmbl}->[$index]->[0],"\t", $gene_pos->{$gene}{$asmbl}->[$index]->[1],"\n" if $asmbl eq 'asmbl_7587';
		}
		foreach my $ind (0..$#ctg_pos)
		{
			my ($start, $end) = @{$ctg_pos[$ind]};
			($start, $end) = ($end, $start) if $start > $end;
			my $seq = substr($fasta_seq{$contigs->{$asmbl}}, $start-1-$fl, ($end-$start+1+$fl*2));
			print OUT '>', join("_",($gene, $asmbl, $introns[$ind])),"\n", $seq, "\n";
		#	print STDERR $introns[$ind],"\n" if $asmbl eq 'asmbl_7587';
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
	my %return;
	while (<IN>)
	{
		next if /^\s+$/;
		#	ID=S10-asmbl_10; Target=asmbl_10 1 399 +
		my ($gene, $asmbl, $start, $end) = $_ =~ /ID=(S\d+)\-(asmbl_\d+).*Target\S+\s(\d+)\s(\d+)/;
		my @t=split /\t/,$_;
		$strand{$asmbl} = $t[6];
		$contigs{$asmbl} = $t[0];
		print STDERR join("\t", ($start, $end)),"\n" if $asmbl eq 'asmbl_7587';
		push @{$return{$gene}{$asmbl}}, [$start, $end, @t[3, 4]];
	}
	close IN;
	return (\%strand, \%contigs, \%return)
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