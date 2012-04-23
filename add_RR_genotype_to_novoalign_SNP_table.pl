#!/usr/bin/perl -w
use strict;
use CN::DB_Utilities;

# read the SAM files seperated for each chrommosome in each accession
# 
my $sUsage = qq(
	perl $0
	<database>
	<username>
	<password>
	<snp table>
	<directory containing SAM files>
	<output file>
);

die $sUsage unless @ARGV >= 6;
my ($db, $user, $pass, $snp_table, $sam_direc, $out_file) = @ARGV;
open (OUT,">$out_file") or die "can't open file $out_file\n";
my $accessions = 8;
my $chr_num = 10;

# DB connection
my $sUser = $user;
my $sPasswd = $pass;
my $dbh = DB_Utilities::connect_to_db($db, $sUser, $sPasswd);


foreach my $chr_id (1..$chr_num)
{
	print_time_comment("query chromosome ". $chr_id. " ...");
	my %snps_hash = query_snp_table($dbh, $snp_table, $chr_id);
	foreach my $acc_id (sort{$a<=>$b} keys %snps_hash)	
	{
		print_time_comment("\t accession: $acc_id");
		my $window_size = 200;
		my %hash_snppos_in_window = assign_snppos_to_window($snps_hash{$acc_id}, $window_size);
		my %snppos_coverage_qual_refbase;
		my $samfile = 'acc'. $acc_id. "_".$chr_id . '.sam';
		$samfile  = $sam_direc . '/'. $samfile;
		open (IN, "<$samfile") or die "can't open file $samfile\n";
		while (<IN>)
		{
			next if /^\s+$/;
			my @line_data = split /\t/, $_;
			my $flag = $line_data[1] #SW 09.20.2011
			my ($align_start, $seq, $quality_string) = @line_data[3, 9, 10];
			my $read_length = length $seq;
			my $num_win_read_start = int($align_start/$window_size);
			my $num_win_read_end = int(($align_start+$read_length-1)/$window_size);
			my @snps_array = @{$hash_snppos_in_window{$num_win_read_start}} if exists $hash_snppos_in_window{$num_win_read_start};
			@snps_array = (@snps_array, @{$hash_snppos_in_window{$num_win_read_end}}) if exists $hash_snppos_in_window{$num_win_read_end};
			next unless @snps_array > 0;
			foreach my $pos (@snps_array)
			{
				if ($pos >= $align_start and $pos <=($align_start + $read_length - 1))
				{
					$snppos_coverage_qual_refbase{$pos}->[0]++;
					my $qual_base = ' ';
					$qual_base = substr $quality_string, $pos-$align_start, 1;			
					my $qual_score = ord($qual_base) - 33;
					$snppos_coverage_qual_refbase{$pos}->[1] += $qual_score;
					my $ref_base = substr $seq, $pos-$align_start, 1;
					$snppos_coverage_qual_refbase{$pos}->[2] = $ref_base unless defined $snppos_coverage_qual_refbase{$pos}->[2];
				}			
			}
		}
		close IN;
		foreach my $snppos (keys %snppos_coverage_qual_refbase)
		{
			my ($coverage, $total_qual, $refbase) = @{$snppos_coverage_qual_refbase{$snppos}};
			my $avg_qual = int($total_qual/$coverage);
			print OUT join("\t", ($acc_id, $chr_id, $snppos, $refbase, $refbase,$avg_qual, $coverage) ),"\n";
		}
	}
}
close OUT;

# Subroutines
sub assign_snppos_to_window
{
	my ($array_ref, $win_size) = @_;
	my %return_hash;
	foreach (@$array_ref)
	{
		my $num_win = int($_/$win_size);
		push @{$return_hash{$num_win}}, $_;
	}
	return %return_hash;
}


sub print_time_comment
{
	my $comment = shift;
	my $now = localtime time;
	print STDERR $now, "\t",$comment,"\n";
}

# query snp table to get the snps which were not covered by all 8 accessions
sub query_snp_table
{
	my ($dbh, $table, $chr) = @_;
	my %snps_hash;
	my $query_txt = "select t1.accession_id, t1.gen_position from $table as t1, 
									 (SELECT chromosome, gen_position, count(*) as ACC FROM $table where chromosome = $chr
									 group by `chromosome`,`gen_position` having ACC<8 ) as t2 
								   where t1.chromosome = $chr and t1.gen_position = t2.gen_position";
	my $query = DB_Utilities::run_query($dbh, $query_txt);
	while (my($acc, $position) = $query->fetchrow_array)
	{
		push @{$snps_hash{$position}}, $acc;
	}
	
	my %return_hash;
	foreach my $snp_pos (keys %snps_hash)
	{
		my @accessions = @{$snps_hash{$snp_pos}};
		my @miss_accs = get_miss_acc(@accessions);
		#if ($snp_pos == 2346626){print join("\t", @accessions),"\n",join("\t", @miss_accs),"\n"; die}
		foreach my $missacc (@miss_accs)
		{
			push @{$return_hash{$missacc}}, $snp_pos;
		}
	}
	
	return %return_hash;
}

sub get_miss_acc
{
	my %array = map{$_, 0} @_;
	my @accs = 1..8;
	my @miss_accs;
	foreach( @accs)
	{
		push @miss_accs, $_ unless exists $array{$_};
	}	
	return @miss_accs;
}












