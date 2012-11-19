#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage= qq(
perl $0
<stop_codon_New.out>
<asembl fasta file>
length cutoff for homopolymer, default 6
);
die $sUsage unless @ARGV;
my ($pst_file, $fasta_file, $cutoff) = @ARGV;
$cutoff = 6 unless defined $cutoff;

my $gn = Bio::DB::Fasta->new($fasta_file);

my %homopolymer_pos = get_homopolymer_pos($gn, $cutoff);

open (IN, $pst_file) or die;
print join("\t", qw(Gene Homopolymer PTC)), "\n";
while (<IN>)
{
	# ID=S1730;Name=S1730%20asmbl_1809%20(S1730);     No Stop_codon
	# ID=S4349;Name=S4349%20asmbl_5119%20(S4349);     160     164 
	my $id = $1 if /(asmbl_\d+)\%/;
	my @t = split /\s+/, $_;
	my $hp_pos = $homopolymer_pos{$id};
	if(/Stop_codon/)
	{
		print $id, "\t", $hp_pos>0?1:0, "\t", 0, "\n";
		next
	}
	my $flag = 0;
	foreach (@t[1..$#t])
	{
		my $pos = 3 * $_;
		$flag = 1 if $pos > $hp_pos;
	}
	print $id, "\t", $hp_pos>0?1:0, "\t", $flag, "\n";
}
close IN;

sub get_homopolymer_pos
{
	my $gn = shift;
	my $cutoff = shift;
	my @ids = $gn->ids;
	my %return;
	foreach my $id (@ids)
	{
		my $seq = $gn->get_Seq_by_id($id)->seq;
		$seq = uc($seq);
		my @hp_pos;
		while ($seq =~ /(A{$cutoff,})/g)
		{
			 my $p = pos($seq) - length($1)  + 1;
			 push @hp_pos , $p;
		}
		while ($seq =~ /(T{$cutoff,})/g)
		{
			 my $p = pos($seq) - length($1)  + 1;
			 push @hp_pos , $p;
		}
		while ($seq =~ /(G{$cutoff,})/g)
		{
			 my $p = pos($seq) - length($1)  + 1;
			 push @hp_pos , $p;
		}
		while ($seq =~ /(C{$cutoff,})/g)
		{
			 my $p = pos($seq) - length($1)  + 1;
			 push @hp_pos , $p;
		}
		
		if(@hp_pos == 0)
		{
			$return{$id} = -1;
		}
		else
		{
			@hp_pos = sort {$a <=> $b} @hp_pos;
			$return{$id} = $hp_pos[0];
		}
	}
	
	return %return;
}



