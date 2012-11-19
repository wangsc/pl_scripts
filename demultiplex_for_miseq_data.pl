#!/usr/bin/perl -w
use strict;
#
# Authored by Shichen Wang Apr.03.2012

my $sUsage = qq(
perl $0 
< Sample Index file, tab-delimited txt file generated from sample sheet; Foramt example: Sample_831_TU	GCCAAT >
< Index fastq files; if multiple files, separated them with comma: I1.fastq,I2.fastq >
< Reads fastq files; if multiple files, separated them with comma: R1.fastq,R2.fastq >
);
die $sUsage unless @ARGV >= 3;

my ($sample_index_file, $index_fastq, $read_fastq) = @ARGV;

my %sample_index = read_sample_index_file($sample_index_file);
my %out_fh = get_file_handles(values %sample_index);
my %read_id_to_sampe = read_index_fastq($index_fastq, \%sample_index);

separate_read_fastq_files($read_fastq, \%read_id_to_sampe, \%out_fh);


sub get_file_handles
{
	my @accs = unique(@_);
	my %return_hash;
	foreach my $acc_id (@accs)
	{
		my $file = $acc_id . '.fastq';
		local *FH;
		open (FH, ">$file") or die "can't open file $file\n";
		$return_hash{$acc_id} = *FH{IO};
	}
	return %return_hash;
}

sub unique
{
	my %return;
	map{$return{$_}=1} @_;
	return keys %return;
}

sub read_sample_index_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die $!;
	while(<IN>)
	{
		chomp; 
		my @t = split /\s+/,$_;
		die "Check the sample index file, it should be tab or space delimited\n" if @t>2;
		$return{$t[1]} = $t[0];
	}
	close IN;
	return %return;
}

sub read_index_fastq
{
	my @files = split /,/, shift;
	my $sample_index_ref = shift;
	my %return;
	foreach my $f (@files)
	{
		open (IN, $f) or die $!;
		my $count = 0;
		my $id;
		while (<IN>)
		{
			$count++;
			if($count == 1)
			{
				$id=$1 if /^(\S+)/;
			}
			elsif($count == 2)
			{
				chomp; 
				$return{$id} = $_;
			}
			elsif($count == 4)
			{
				$count=0;
			}
		}
		close IN;
	}
	
	foreach my $read_id (keys %return)
	{
		my $seq = $return{$read_id};
		$seq = uc($seq);
		my @acc_ids;
		my %acc_match;
		my $max_match = 0;
		foreach my $index_seq (keys %{$sample_index_ref})
		{
			my $acc = $sample_index_ref->{$index_seq};
			$index_seq = uc($index_seq);
			my $m = match($index_seq, $seq);
			$max_match = $m if $m >$max_match;
			$acc_match{$acc} = $m;
		}
		
		foreach (keys %acc_match)
		{
			push @acc_ids, $_ if $acc_match{$_} == $max_match and $max_match >= ((length $seq) - 2);
		}
		
		if(@acc_ids >1)
		{
			print STDERR "! Read $read_id has ambiguous accession ID: ", join(" ", @acc_ids), "\n";
			#$return{$read_id} = "NA";
		}
		else
		{
			$return{$read_id} = $acc_ids[0] if @acc_ids;
		}
	}
	
	return %return;
}

sub match
{
	my ($sa, $sb) = @_;
	my @sa=split //, $sa;
	my @sb=split //, $sb;
	
	my $return = 0;
	foreach (0..$#sa)
	{
		$return ++ if $sa[$_] eq $sb[$_]
	}
	return $return;
}

sub separate_read_fastq_files
{
	# ($read_fastq, \%read_id_to_sampe, \%out_fh);
	my @files = split /,/, shift;
	my $read_id_to_sample = shift;
	my $out_fh = shift;
	foreach my $f (@files)
	{
		print STDERR "Processing file $f ... \n";
		open (IN, $f) or die $!;
		my $count = 0;
		my @arr;
		my $fh;
		while (<IN>)
		{
			chomp;
			my $line = $_;
			$count++;
			if($count == 1)
			{
				my ($id, $p_index) = $line=~/^(\S+)\s(\d)/;
				push @arr, $id."/".$p_index;
				$fh = $out_fh->{$read_id_to_sample->{$id}} if exists $read_id_to_sample->{$id};
			}
			elsif($count == 4)
			{
				push @arr, $line;
				print {$fh} join("\n", @arr), "\n" if defined $fh;
				@arr = ();
				$count = 0;
				undef $fh;
			}
			else
			{
				push @arr, $line;
			}
		}
		close IN;
	}
	
}








