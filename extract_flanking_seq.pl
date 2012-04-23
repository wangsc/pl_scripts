#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage = "perl $0 <rice.fasta> <rice_id.list> <flanking length>\n";
die $sUsage unless @ARGV >= 3;
my($fasta, $list, $flank) = @ARGV;

my $gn = Bio::DB::Fasta->new($fasta);
my %list = read_list($list);
foreach my $id (keys %list)
{
	my $seq = $gn->get_Seq_by_id($id)->seq;
	my ($start, $end) = @{$list{$id}};
	$start = ($start>$flank)?$start-$flank:0;
	$end = ($end+$flank)<=(length $seq)?$end+$flank:(length $seq);
	print '>',$id,'_', $start,"_",$end,"\n";
	print $gn->seq($id, $start=>$end),"\n";
}

sub read_list
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my @t=split /\s+/,$_;
		$return{$t[0]} = [sort{$a<=>$b}@t[1..2]];
	}
	close IN;
	return %return;
}