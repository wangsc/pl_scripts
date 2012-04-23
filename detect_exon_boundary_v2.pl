#!/usr/bin/perl -w
use strict;
use Search::Dict;

my $sUsage = qq(
# The transcript assemblies were blat against CS genomic sequences at first;
# This script will then detect potential exon boundaries in the transcript assemblies.

Usage: 
perl $0
<CS read length file>
<output file>
<minimum similarity, 95>
<minimum leghtn of fragment, 50>
<blat result files>
);
die $sUsage unless @ARGV >= 3;
my ($cs_read_len_file, $out_file, $min_sim, $min_len, @blat_files) = @ARGV;
$min_sim = 95 unless defined $min_sim;
$min_len = 50 unless defined $min_len;

print_time_comment("Reading CS read length file ...");
open(my $LEN, "$cs_read_len_file") or die $!;
print_time_comment("Reading BLAT files ...");
my %exon_boundary = read_blat_files(\@blat_files, $min_sim, $min_len, $LEN);
my %potential_boundaries = detect_boundaries(\%exon_boundary);
print_time_comment("Output results in file $out_file ...");
output(\%potential_boundaries, $out_file);


# Subroutines
sub output
{
	my($hashref, $outfile) = @_;
	open (OUT, ">$outfile") or die $!;
	foreach my $id (keys %{$hashref})
	{
		print OUT ">", $id,"\n";
		foreach (@{$hashref->{$id}})
		{
			print OUT join("\t", @$_),"\n";
		}
		print OUT "\n";	
	}
	close OUT;
}

sub detect_boundaries
{
	my $hashref = shift;
	my $delta = 0.1;
	my $cut_percent = 0.9;
	my %return;
	foreach my $id (keys %{$hashref})
	{
		my @array;
		foreach (keys %{$hashref->{$id}})
		{
			$array[$_] = $hashref->{$id}->{$_};
		}
		my @peaks = peak_detect(\@array, $delta);
		my @filtered_peaks = filter_peaks($cut_percent, @peaks);
		my $cut_percent_second = 0.6;
		my @filtered_peaks_second = filter_peaks($cut_percent_second, @filtered_peaks);
		$return{$id} = [@filtered_peaks_second];
	}
	return %return;
}

sub filter_peaks
{
	my ($cut_percent, @peaks) = @_;
	my @filtered;
	my $total = 0;
	my @peak_values;
	foreach (@peaks)
	{
		$total += $_->[1];
	}
	my $cutoff = ($total/(scalar @peaks))*$cut_percent;
	foreach (@peaks)
	{
		push @filtered, $_ if $_->[1] >= $cutoff;
	}
	return @filtered;
}

sub peak_detect
{
	my ($array_ref, $delta) = @_;
	my @peaks;
	map{$array_ref->[$_] = 0 unless defined $array_ref->[$_]}(0..(scalar @$array_ref));
	my ($min, $max) = (99999, -99999);
	my $maxpos;
	my $minpos;
	my $max_flag = 1;
	foreach my $i (0 .. (scalar @$array_ref - 1))
	{
		my $current = $array_ref->[$i];
		die "Not defined \n" unless defined $current;
		if ($current > $max){$max = $current; $maxpos = $i}
		if ($current < $min){$min = $current; $minpos = $i}		
		if ($max_flag)
		{
			if($current < ($max*(1-$delta)))
			{
				print $i,' ** ', $current," ** ", $max,"\n";
				push @peaks, [$maxpos, $max];
				$minpos = $i;
				$min = $current;
				$max_flag = 0;
			}		
		}
		else
		{
			if ($current > $min*(1+$delta))
			{
				$max = $current;
				$maxpos = $i;
				$max_flag = 1;
			}
		}

	}
	return @peaks;
}

sub read_length_file
{
	my $file = shift;
	my @return;
	open (IN, $file) or die "can't open file $file \n";
	while(<IN>)
	{
		chomp; 
		next if /^\s+$/;
		my @t=split /\t/,$_;
		push @return, @t;
	}
	close IN;
	return @return;
}



sub read_blat_files
{
	my($files_ref, $min_sim, $min_len, $LEN) = @_;
	my %gene_boundary;
	my $debug = 1;
	foreach my $file (@$files_ref)
	{
		print_time_comment("\tReading $file ...");
		open (IN, "$file") or die $file;
		while(<IN>)
		{
			#	asmbl_2 GJ20S3U01A02OO  97.90   238     2       1       90      324     1       238     6.9e-124        441.0
			# asmbl_2 GJ20S3U01A02OO  99.39   164     0       1       325     488     255     417     8.5e-87 317.0	
			chomp;
			next unless /^\S+/;
			my @data = split /\t/, $_;
			my ($similarity, $length, $num_gap) = @data[2, 3, 5];
			next unless $similarity >= $min_sim and $length >= $min_len and $num_gap == 0;
			my ($gene_id, $read_id, $gene_start, $gene_end, $read_start, $read_end) = @data[0, 1, 6..9];
			#print STDERR '$read_id ', $read_id,"\n" if $debug;
			look $LEN, $read_id;
			my $line = <$LEN>;
			#print STDERR $line,"\n";
			my $read_length = $1 if $line =~ /\s+(\d+)/;
			#print STDERR $read_id, " ** ", $read_length,"\n" if $debug; $debug=0;#exit 1;
			$gene_boundary{$gene_id}{$gene_start}++ unless $read_start == 1 or $read_start == $read_length;
			$gene_boundary{$gene_id}{$gene_end}++ unless $read_end == 1 or $read_end == $read_length;		
		}
		close IN;		
	}
	return %gene_boundary;
}

sub print_time_comment
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}
