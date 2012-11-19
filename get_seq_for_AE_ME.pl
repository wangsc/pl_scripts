#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage = qq(
perl $0
<gtf3 file>
<fasta file>
<AE_CE_ME.out>
);
die $sUsage unless @ARGV >= 3;

my ($gtf_file, $fasta_file, $exon_file) = @ARGV;

my %ctg_id = read_gtf_file($gtf_file);
my $gn_obj = Bio::DB::Fasta->new($fasta_file);

my ($acquired_exon_neighbor, $missing_exons) = read_exon_file($exon_file);

# output
my $ae_out = 'AE_neighbor.fasta';
open (AE, ">$ae_out") or die "can't open file $ae_out\n";
foreach my $sid (keys %$acquired_exon_neighbor)
{
	my $ctg = $ctg_id{$sid};
	print STDERR '$ctg ', $ctg, "\t", '$sid ', $sid, "\n" unless defined $ctg;
	foreach (@{$acquired_exon_neighbor->{$sid}})
	{
		my ($start, $end) = split /:/, $_;
		next unless ($end-$start) >= 50;
		my $seq = $gn_obj->seq($ctg, $start => $end);
		print AE ">", join(':', ($sid, $start, $end)), "\n", $seq, "\n";		
	}

}
close AE;

my $me_out = 'ME.fasta';
open (ME, ">$me_out") or die;
foreach my $sid (keys %$missing_exons)
{
	my $ctg = $ctg_id{$sid};
	foreach (@{$missing_exons->{$sid}})
	{
		my ($start, $end) = split /:/, $_;
		next unless ($end-$start) >= 50;
		my $seq = $gn_obj->seq($ctg, $start => $end);
		print ME ">", join(':', ($sid, $start, $end)), "\n", $seq, "\n";		
	}
	

}
close ME;

# Subroutines
sub read_gtf_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while (<IN>)
	{
		next if /^\s+$/;
		my @t = split /\s+/, $_; 
		my $sid = $1 if /ID=(S\d+)/;
		$return{$sid} = $t[0];
	}
	close IN;
	return %return;
}

sub read_exon_file
{
	my $file = shift;
	open (IN, $file) or die;
	my (%ae_neighbor, %me);
	my $line_cnt = 0;
	my $id;
	my %record;
	while(<IN>)
	{
		next if /^\s+$/;
		$line_cnt++; 
		if(/^>(\S+)/)
		{
			$id = $1;
			$line_cnt = 0;
			next;
		}
		chomp;
		my @t = split /\s+/, $_;
		$record{$id}[$line_cnt-1] = [@t];		
	}
	close IN;
	
	foreach my $sid (keys %record)
	{
		my @ae = @{$record{$sid}[0]};
		my @ce = @{$record{$sid}[1]};
		my @me = @{$record{$sid}[2]};
		
		if (@ae > 1)
		{
			my @neighbors = check_neighbors(\@ae, \@ce);
			$ae_neighbor{$sid} = [@neighbors] if @neighbors;
		}
		
		if(@me > 1)
		{
			$me{$sid} = [@me[1..$#me]]
		}
	}
	return(\%ae_neighbor, \%me);
}

sub check_neighbors
{
	my ($arr1, $arr2) = @_;
	my @return;
	
	foreach (@$arr1)
	{
		next unless /:/;
		my ($s1, $e1) = split /:/, $_;
		foreach (@$arr2)
		{
			next unless /:/;
			my ($s2, $e2) = split /:/, $_;
			if( abs($s2-$e1) <= 3 )
			{
				push @return, $_;
			}
		}
	}
	return unique(@return);
}

sub unique
{
	my %h = map{$_, 1} @_;
	return (keys %h);	
}






 