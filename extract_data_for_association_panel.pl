#!/usr/bin/perl -w
use strict;

my $sUsage = "perl $0 
<combined file> 
<association samp list>
<output>
\n";
die $sUsage unless @ARGV >= 3;
my ($combined_file, $sample_list, $output) = @ARGV;
my ($sample_hash, @samples) = read_sample_file($sample_list);
read_combined_file($combined_file, $sample_hash, $output, @samples);

sub read_sample_file
{
	my $file = shift;
	my %return;
	my @samples;
	open (IN, "$file") or die;
	while (<IN>)
	{
		chomp;
		next if /^\s+$/;
		s/\s//g;
		push @samples, $_;
		$return{$_} = 1;
	}
	close IN;
	return (\%return, @samples);
}

sub read_combined_file
{
	my $file = shift;
	my $asso_sample = shift;
	my $outfile = shift;
	my @asso_sam = @_;
	open (IN, "$file") or die;
	my @record;
	my %all_samples;
	my $max_sample;
	my %associ_index;
	my %chr_pos;
	my %return;
	while (<IN>)
	{
		s/\"//g;
		s/,/\t/g;
		s/\s+$//;
		chomp;
		my @data = split /\t/,$_;
		if (/^Index/)
		{
			#%all_samples = map{ $_, $data[$_]} (2..$#data);
			foreach (2..$#data)
			{
				$data[$_] =~ s/ //g;
				$all_samples{$_} = $data[$_];
			}
			my $debug = 1; 
			#print scalar keys %all_samples,"***\n"; foreach (keys %all_samples){print $_,"**", $all_samples{$_},"**\n" if $debug; $debug=0}
			$max_sample = (scalar @data) - 2;		
			foreach (keys %all_samples)
			{
				$associ_index{$_} = 1 if exists $asso_sample->{$all_samples{$_}};
			}
		}
		next unless /^\d/;
		print STDERR join('**', @data),"\n" if $data[0] == 575;
		my %genotypes;
		foreach (2..$max_sample+1)
		{
			next unless exists $associ_index{$_};
			$genotypes{$data[$_]}++;
			#$genotypes{transform_genotype($data[$_])}++;
		}
		my $num_genotype = 0;
		foreach (keys %genotypes)
		{
			next unless /[ATCGB]/;
			$num_genotype++;
		}
		print STDERR join("\t", keys %genotypes),"\n" if $data[0] == 575;
		unless ($num_genotype > 1){print $data[0],"\n"; next}; # non-polymorphic snps
		
		push @record, [@data];
		if (@data > $max_sample+2)
		{
			$chr_pos{$data[0]} = [@data[$max_sample+2, $max_sample+3]]
		}
	}
	close IN;
	my %all_samples_reverse = reverse %all_samples; # name->index
	my $debug = 1;
	my @snps_array;
	my @snps_array_all;
	my %for_ld;
	my $flag = 1;
	foreach my $samp_name (@asso_sam)
	{
		
		my $index = $all_samples_reverse{$samp_name};
		next unless (defined $index);
		foreach (@record)
		{
			my $snp_index = $_->[0];
			push @snps_array, $snp_index if $flag;
			my $genotype = $_->[$index];
			push @{$return{$samp_name}}, transform_genotype($genotype);
			#print STDERR $genotype,"**", transform_genotype($genotype),"\n";
			if (exists $chr_pos{$snp_index})
			{
				push @{$for_ld{$samp_name}}, transform_genotype($genotype);
			}
		}
		$flag = 0;		
	}
#	
	open (OUT, ">$outfile") or die ;
	print OUT "\tSNP", join("\tSNP", @snps_array),"\n";
	#print STDERR 'scalar @record ', scalar @record, "\n";
	#print STDERR 'scalar @snps_array ', scalar @snps_array,"\n";
	foreach (@asso_sam)
	{
		unless (exists $return{$_}){print OUT $_, "\t--"x(scalar @snps_array), "\n"; next}
		#print STDERR $_,"\n";
		print OUT $_, "\t", join("\t", @{$return{$_}}),"\n";
		#print STDERR $_, join("\t", @{$return{$_}}),"\n";
	}
	close OUT;
#	
	open (OUT, ">chr_position.csv") or die;
	my @mapped_snps;
	foreach (@snps_array)
	{
		print OUT 'SNP', $_, "\t", join("\t", @{$chr_pos{$_}}),"\n" if exists $chr_pos{$_};
		push @mapped_snps, $_ if exists $chr_pos{$_};
	}
	close OUT;
#	
	open (OUT, ">LD_file.csv") or die $!;
	print OUT "\tSNP", join("\tSNP", @mapped_snps),"\n";
	foreach (@asso_sam)
	{
		unless (exists $for_ld{$_}){print OUT $_, "\t--"x(scalar @mapped_snps), "\n"; next}
		#print STDERR $_,"\n";
		print OUT $_, "\t", join("\t", @{$for_ld{$_}}),"\n";
	}
	close OUT;	
}

sub transform_genotype
{
	my $g = shift;
	return '--' unless $g =~ /[ATCGB]/;
	my @gs = split //, $g;
	if ($gs[0] eq $gs[1])
	{
		return ($g=~/A/)?'AA':'BB';
	}
	else{return 'AB'}
}

sub unique
{
	my %h;
	map{$h{$_}++} @_; 
	return (keys %h);
}

