#!/usr/bin/perl -w 
use strict;
use Bio::DB::Fasta;
use File::Basename;

my $sUsage = qq
(
perl $0 
<folder containing BLAT filtered files>
<folder containing survey sequence fasta files>
<snps_mapped_chr.out>
<dp_8_ref_alt_allele.out>
);

die $sUsage unless @ARGV >= 4;

my ($blatfile_folder, $fasta_folder, $mapped_chr_file, $ref_alt_file) = @ARGV;

my @blat_files = <$blatfile_folder/*filtered>;
my @fasta_files = <$fasta_folder/*.fa>;

my %snp_alleles = read_allele_file($ref_alt_file);
my %snp_chr = read_mapped_chr_file($mapped_chr_file);

my %recorder;

foreach my $blat_file (@blat_files)
{
	my $basefile = basename($blat_file);
	my $chr = $1 if $basefile =~ /blat_(\w+)/;
	print STDERR '$chr ', $chr, "\n";
	my $fasta;
	foreach (@fasta_files){if(/$chr/){$fasta=$_; last} }
	die "NO fasta file\n" unless defined $fasta;
	print STDERR '$fasta ', $fasta, "\n";
	my $fasta_db = Bio::DB::Fasta->new($fasta);
	
	open (BLAT, $blat_file) or die "can't open file $blat_file\n";
	while (<BLAT>)
	{
		# BobWhite_c1027_360      3976072 100.00  101     0       0       1       101     4327    4427    8.4e-51 198.0
		chomp;
		my @line_data = split /\s+/, $_;
		next unless exists $snp_chr{$line_data[0]};
		next if $line_data[5] > 0;
		my $snp_pos_relative = 51 - $line_data[6];
		my $snp_pos = $line_data[8] + ($line_data[8]>$line_data[9]?-1:1)*$snp_pos_relative;
		my $frag_seq = $fasta_db->seq($line_data[1], $line_data[8]=>$line_data[9]);
		my $survey_seq_allele = uc(substr($frag_seq, $snp_pos_relative, 1));
		push @{$recorder{$line_data[0]}}, join(":", ($chr, $line_data[1], $snp_pos, $survey_seq_allele));
	}
	close BLAT;
}

foreach my $snp (keys %recorder)
{
	my $mapped_chr = "ND";
	my @data = @{$recorder{$snp}};
	unless (exists $snp_alleles{$snp})
	{
		if (@data == 3)
		{
			my %tmp_hash;
			foreach (@data)
			{
				my ($chr, $contig_id, $snp_pos, $allele) = split /:/, $_;
				push @{$tmp_hash{$allele}}, $chr;
			}
			foreach (keys %tmp_hash)
			{
				if(@{$tmp_hash{$_}} == 1)
				{
					$mapped_chr = $tmp_hash{$_}[0];
					last;
				}
			}
		}

		print $snp, "\t", join("+", @data), "\t", "NA", "\t", $mapped_chr, "\n";
		next;
	}	
	
	foreach (@data)
	{
		my ($chr, $contig_id, $snp_pos, $allele) = split /:/, $_;
		if($allele eq $snp_alleles{$snp}[1])
		{
			$mapped_chr = $chr;
			last
		}
	}
		
	print $snp, "\t", join("+", @data), "\t", exists $snp_alleles{$snp}?join("/", @{$snp_alleles{$snp}}):"NA", "\t", $mapped_chr, "\n";	
}


# Subroutines
sub read_allele_file
{
	my $file = shift or die;
	open (IN, $file) or die;
	my %return;
	while (<IN>)
	{
		chomp;
		my @t = split /\s+/, $_;
		$return{$t[0]} = [split //, $t[1]];
	}
	close IN;
	return %return;
}

sub read_mapped_chr_file
{
	my $file = shift or die;
	open (IN, $file) or die;
	my %return;
	while (<IN>)
	{
		chomp;
		my @t = split /\s+/, $_;
		$return{$t[0]} = $t[1];
	}
	close IN;
	return %return;
}









