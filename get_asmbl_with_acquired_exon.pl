#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage = qq(
perl $0
<isoform gff3 file>
<acquired_exon list>
<asmbl_fasta file>
);
die $sUsage unless @ARGV;

my ($gff_file, $list_file, $fasta_file) = @ARGV;

my $gn_obj = Bio::DB::Fasta->new($fasta_file);

my ($asmbl_pos_hashref, $gene_asmbl_hashref, $strand_hashref) = read_gff($gff_file);

open (IN, $list_file) or die ;

while (<IN>)
{
	# >S6146_10522_10811
	next unless /^>/;
	chomp; 
	my $in = $_;
	s/>//;
	my $asmbl_id;
	my ($gene_id, $ae_start, $ae_end) = split /_/, $_;
	my @asmbls = @{$gene_asmbl_hashref->{$gene_id}};
	my @ids = map{$gene_id . "-" . $_}@asmbls;
	foreach my $id (@ids)
	{
		my @positions = @{$asmbl_pos_hashref->{$id}};
		foreach my $p (@positions)
		{
			my ($start, $end) = split /\t/, $p;
			if($start <= $ae_start and $end >= $ae_end)
			{
				$asmbl_id = $id;
				last;
			}
		}
		last if defined $asmbl_id;
	}
	
	my @pos_in_asmbl = calc_new_pos($asmbl_pos_hashref->{$asmbl_id}, $ae_start, $ae_end, $strand_hashref->{$asmbl_id});
	print join("\t", $in, $asmbl_id, @pos_in_asmbl), "\n";
	my $asm = $1 if $asmbl_id =~ /(asmbl_\d+)/;
	my $seq = $gn_obj->get_Seq_by_id($asm)->seq;
	print $seq, "\n";
}

# Subroutines
sub read_gff
{
	my $file  = shift;
	open (IN, $file) or die;
	my (%asmbl_pos, %gene_asmbl, %strand);
	while (<IN>)
	{
		# contig00869     PASA    cDNA_match      542     790     .       +       .       ID=S11-asmbl_11; Target=asmbl_11 1 249 +
		next unless /\S/;
		my @t = split /\s+/, $_;
		my $id = $1 if /ID=(\S+?)\;/;
		push @{$asmbl_pos{$id}}, join("\t", @t[3,4]);
		$strand{$id} = $t[6];
		my @temp = split /-/, $id;
		push @{$gene_asmbl{$temp[0]}}, $temp[1];		
	}
	close IN;
	
	return (\%asmbl_pos, \%gene_asmbl, \%strand)
}

sub calc_new_pos
{
	my ($pos_ref, $ae_start, $ae_end, $strand) = @_;
	my ($total_length, $left, $new_start, $new_end);
	foreach (@$pos_ref)
	{
		my ($s, $e) = split /\s+/, $_;
		$total_length += ($e - $s + 1);
		if($s <= $ae_start and $e >= $ae_end)
		{
			$left += ($ae_start - $s + 1);
			$new_start = $left;
			#$new_end = $new_start 
		}
		else
		{
			$left += ($e - $s + 1);
		}
	}
	
	$new_end = $new_start + abs($ae_start - $ae_end);
	if($strand eq '-')
	{
		$new_start = $total_length - $new_end  + 1;
		$new_end = $new_start + abs($ae_start - $ae_end);
	}
	
	return ($new_start, $new_end, $total_length);
	
}





