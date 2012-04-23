#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<fasta file>
<simrank result>
<similarity cutoff  80,85...>
<output file>
);
die $sUsage unless @ARGV >= 4;

my ($fastafile, $simrank_file, $sim_cutoff, $output) = @ARGV;
my %read_length = get_length($fastafile);
my @remained = process_simrank_result($simrank_file, $sim_cutoff, \%read_length);
output($output, @remained);

sub output
{
	my $file = shift;
	my @arr= @_;
	open(OUT, ">$file") or die;
	foreach (@arr){print OUT join("\t", @$_),"\n"}
	close OUT;
}

sub process_simrank_result
{
	my $file = shift;
	my $cutoff = shift;
	my $read_length = shift;
	open (IN, $file) or die "can't open file $file \n";
	my (@recorder, %duplicates);
	my @return;
	my $id;
	while(<IN>)
	{
		chomp; next if /^\s+$/;
		if(/Processing (\S+)/)
		{
			$id = $1;
			push @recorder, $id;
			next;
		}
		next unless defined $id;
		next if exists $duplicates{$id};
		if(/^$id/)
		{
			my @dup;
			while(/(\S+):(\d+)/g)
			{
				my ($sim_id, $sim) = ($1, $2);
				last if $sim < $cutoff;
				push @dup, $1;
			}
			my @array = sort{$read_length->{$b}<=>$read_length->{$a}} (@dup);
			push @return, [@array, $read_length->{$array[0]}];
			map{$duplicates{$_}=1}@array;
		}
	}
	close IN;
#	my @remained = grep {not exists $duplicates{$_}} @recorder;
#	return @remained;
  return @return;
}

sub get_length
{
	my $file = shift;
	open (IN, $file) or die "can't open file $file \n";
	my $id;
	my %return;
	while(<IN>)
	{
		if(/^>(\S+)/)
		{
			$id = $1;
			next;
		}
		chomp;
		$return{$id} += length $_;
	}
	close IN;
	return %return;
}




