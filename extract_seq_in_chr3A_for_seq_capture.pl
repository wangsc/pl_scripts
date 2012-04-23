#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<pasa transcript isoform gff>
<contig masked fasta>
<output>
<RBS_ortho_positions_in_chr3A, optional>
);
die $sUsage unless @ARGV >= 3;
my($iso_gff, $ctg_file, $outfile, $rbs_file) = @ARGV;

my %contig_gene_positions = read_iso_gff($iso_gff);
my %rbs_pos = defined $rbs_file?read_rbs_file($rbs_file):();
my %contig_fasta = read_fasta($ctg_file);
my %nonN_pos = detect_nonN_pos(\%contig_fasta);
my %ctg_pos_for_output = calculate_output_regions(\%nonN_pos, \%contig_gene_positions, \%rbs_pos);
output(\%ctg_pos_for_output, \%contig_fasta, $outfile);

# Subroutines

sub output
{
	my ($ctg_pos_out, $ctg_fasta, $outfile) = @_;
	open (OUT ,">$outfile") or die "can't open file $outfile\n";
	foreach my $ctg (keys %$ctg_pos_out)
	{
		foreach my $pos (@{$ctg_pos_out->{$ctg}})
		{
			my ($start, $end) = split /_/, $pos;
			die 'No contig: ', $ctg,"\n" unless exists $ctg_fasta->{$ctg};
			my $seq = substr($ctg_fasta->{$ctg}, $start-1, ($end-$start+1));
			print OUT '>',$ctg,'_',$pos,"\n", $seq, "\n";
		}
	}
	close OUT;
}

sub read_rbs_file
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while(<IN>)
	{
		chomp;
		next if /^\s+$/;
		my @t = split /\t/,$_;
		my $id = shift @t;
		$return{$id} = [@t];
	}	
	close IN;
	return %return;
}

sub calculate_output_regions
{
	my ($nonN_pos, $ctg_gen_pos, $rbs_pos) = @_;
	my %return;
	foreach my $ctg (keys %$ctg_gen_pos)
	{
		my @nonN_pos = @{$nonN_pos->{$ctg}};
		foreach my $gene (keys %{$ctg_gen_pos->{$ctg}})
		{
			print STDERR $ctg,"\t", $gene,"\n" unless exists $ctg_gen_pos->{$ctg}->{$gene};
			my @array = sort{$a<=>$b} (@{$ctg_gen_pos->{$ctg}->{$gene}});
			my ($start, $end) = @array[0, -1];
			foreach (@nonN_pos)
			{
				if($start >= $_->[0] and $end <= $_->[1])
				{
					push @{$return{$ctg}}, join('_', @$_);
				}
			}
		}
	}
	
	foreach my $ctg (keys %$rbs_pos)
	{
		my @array = sort{$a<=>$b}@{$rbs_pos->{$ctg}};
		my ($start, $end) = @array[0, -1];
		my @nonN_pos = @{$nonN_pos->{$ctg}};
		foreach (@nonN_pos)
		{
			if($start >= $_->[0] and $end <= $_->[1])
			{
				push @{$return{$ctg}}, join('_', @$_);
			}
		}
	}
	
	map{$return{$_} = [unique(@{$return{$_}})]} keys %return;
	
	return %return;	
}

sub unique
{
	my %h = map{$_, 1} @_;
	return keys %h;
}


sub detect_nonN_pos
{
	my $ctg_fasta = shift;
	my %return;
	foreach (keys %$ctg_fasta)
	{
		my $seq = $ctg_fasta->{$_};
		my @array;
		while($seq=~/([^N]+)/g)
		{
			my ($pos, $length) = (pos($seq), length $1);
			my ($start, $end) = ($pos-$length+1, $pos);
			push @array, [$start, $end];
		}
		$return{$_} = [@array];
	}
	return %return;
}

sub read_iso_gff
{
	my $file = shift;
	my $rbs_file = shift;
	open (IN, $file) or die;
	my %return;
	while(<IN>)
	{
		next unless /^\S+/;
		my $gene = $1 if /ID=(.*?)\-/;
		my @t = split /\t/,$_;
		my($ctg_id, $start, $end) = @t[0, 3, 4];
		push @{$return{$ctg_id}{$gene}}, ($start, $end);
	}
	close IN;	
	return %return;
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