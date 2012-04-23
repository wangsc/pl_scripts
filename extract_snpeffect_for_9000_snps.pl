#!/usr/bin/perl -w
use strict;
# compare the SNP effect file genreated by blastx and summary to the 9000 SNPs which were used to design the SNP
# genotyping assay.
# Get the SNP for the 9000 SNPs, if SNPs were not reported in the SNP effect file,
# output the flanking sequencings of those SNPs to a fasta file.
# SW 05.10.2011

my $sUsage = qq(
perl $0
<snp effect file>
<9000 csv file>
<output file for snp effect>
<fasta file for snps not recorded in the snp effect file>
);
die $sUsage unless @ARGV >= 4;
my ($snpeffect_file, $csv_file, $output_snp, $output_fasta) = @ARGV;

# main
my %snpeffect = read_snpeffect_file($snpeffect_file);
my %snps_in_csv = read_csv_file($csv_file);
output(\%snpeffect, \%snps_in_csv, $output_snp, $output_fasta);

# subroutines

sub read_snpeffect_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file \n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my @line_data = split /\t/, $_;
		my $id  = shift @line_data;
		#print 'ID: ', $id,"\n";
		#print 'DATA: ', join("\t", @line_data);
		$return_hash{$id} = [@line_data];
	}
	return %return_hash;
}

sub read_csv_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file \n";
	while (<IN>)
	{
		next unless /^wsnp/;
		my($id, $seq) = split /,/, $_; # id; wsnp_Ku_c29287_39194579
		$return_hash{$id} = $seq;  
	}
	return %return_hash;
}

sub output
{
	my ($snpeffect_hashref, $csv_snps_hashref, $out_snp, $out_fasta) = @_;
	open (SNP, ">$out_snp") or die "can't open file $out_snp\n";
	open (FASTA, ">$out_fasta") or die "can't open file $out_fasta\n";
	my %exists_snps;
	foreach my $id (keys %{$snpeffect_hashref}) 
	{
		# id: CAP8_mira1_rep_c3917_1928346_1928869_1928366_1928366_[C/T]
    #     Ra_mira1_c39488_47184687_47185205_47185072_47185072_[A/G]
		my @data = split /_/, $id;
		my $prefix = 'wsnp';
		my ($acc_name, $contig, $snp_pos) = @data[0, -6, -2];
		unless ($data[2] =~ /c\d+/){$acc_name = $acc_name . '_' . $data[2]}
		my $csv_id = join('_', ($prefix, $acc_name, $contig, $snp_pos));
		if (exists $csv_snps_hashref->{$csv_id})
		{
			next if (exists $exists_snps{$csv_id});
			print SNP $csv_id,"\t", join("\t", @{$snpeffect_hashref->{$id}}),"\n" ;
			$exists_snps{$csv_id} = 1;
		}
	}
	foreach (keys %$csv_snps_hashref)
	{
		next if exists $exists_snps{$_};
		my $seq = $csv_snps_hashref->{$_};
		my $snp_pos;
		if ($seq =~ /\[/g){$snp_pos = pos($seq)};
		my $snp_info = $1 if $seq =~ /(\[.*\])/;
		$seq =~ s/(\[(\w)\/\w\])/$2/;
		my $fasta_id = join('_', ($_, $snp_pos, $snp_info));
		print FASTA '>',$fasta_id,"\n", $seq,"\n";
	}
	close SNP;
	close FASTA;
}





