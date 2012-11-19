#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::DB::Fasta;
use Bio::Tools::Run::Alignment::Muscle;

my $sUsage = qq(
perl $0
<cluster file with chromosome assigned>
<all9line fasta file>
<output file>
);

die $sUsage unless @ARGV >= 3;

my ($cluster_file, $infasta, $output) = @ARGV;
my $gn_obj = Bio::DB::Fasta->new($infasta);

open (OUT, ">$output") or die;

my %ctg_homeologs = read_cluster_file($cluster_file);

foreach my $ctg (keys %ctg_homeologs)
{
	my @sub_chrs = keys %{$ctg_homeologs{$ctg}};
	next if @sub_chrs > 3;
	foreach my $sub_chr (@sub_chrs)
	{
		my @seq_ids = @{$ctg_homeologs{$ctg}{$sub_chr}};
		my $consensus = get_consensus($gn_obj, \@seq_ids);
		print OUT ">", $ctg,":", substr($sub_chr,0,2), "\n";
		print OUT $consensus, "\n";
	}
}
close OUT;

# Subroutines
sub read_cluster_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die "can't open file $file\n";
	while (<IN>)
	{
		chomp; 
		# Kukri_mira1_c11075  Excalibur_mira1_c28699:3AL+Kukri_mira1_c11075:3DL+Kukri_mira1_c39215:3AL
		my @t = split /\s+/, $_;
		my @p = split /\+/, $t[1]; 
		foreach (@p)
		{
			my @arr = split /:/, $_;
			push @{$return{$t[0]}{substr($arr[1],0,2)}}, $arr[0] unless $arr[1] eq "N";
		}
	}
	close IN;
	return %return;
}

sub get_consensus
{
	my ($gn, $lists) = @_;
	if(@$lists == 1)
	{
		my $seq = $gn->get_Seq_by_id($lists->[0])->seq;
		return $seq;
	}
	else
	{
		my $tmp_file = "tmp.fa";
		open (T, ">$tmp_file") or die "can't open file $tmp_file\n";
		foreach my $id (@$lists)
		{
			my $seq = $gn->get_Seq_by_id($id)->seq;
			print T ">$id\n$seq\n";
		}
		close T;
		
		my $aln_factory = Bio::Tools::Run::Alignment::Muscle->new();
		my $aln = $aln_factory->align($tmp_file);
		my $consensus = $aln->consensus_string();
		return $consensus;
	}
	
}












