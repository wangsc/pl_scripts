#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<min20_1_2.snps_nonredundant_mindistance100bp_noLowcomplex_removeNotworked9kSNPs_9k_bristol>
<min20_1_2_Nbp_flanking.fasta>
<length of flanking , 50,60 ...>
<output file>
);
die $sUsage unless @ARGV >= 3;

my($snpfile, $fasta_file, $length_flank, $outfile) = @ARGV;

my %snps = read_snp_file($snpfile);
my %fasta = read_fasta($fasta_file, $length_flank);

open (OUT, ">$outfile") or die "can't open file $outfile\n";
foreach my $id (keys %snps)
{
	my @alleles = (split/\t/,$snps{$id})[1,2];
	die "Not exist in fasta: ", $id,"\n" unless exists $fasta{$id};
	my $seq = $fasta{$id}->[1];
	my $len = length $seq;
	my $snppos = $len - $fasta{$id}->[0] -1;
	my $sub = '['.join('/', @alleles).']';
	substr($seq, $snppos, 1) = $sub;
	print OUT $snps{$id},"\t", $seq,"\n";
}
close OUT;

# subroutines;
sub read_snp_file
{
	my $file = shift;
	# RAC875_mira1_c18915     71      A_1     AC_7    C_1     0.1     BS00002550
	# RAC875_mira1_c18915     284     T_3     TA_7    0.1
	# Excalibur_mira1_s113903 79      A_2     CA_4    0.1     wsnp_Ex_rep_c68570_67412807     BS00004026
	open (IN, "$file") or die;
	my %return;
	while(<IN>)
	{
		chomp;
		my @t = split /\t/,$_;
		print join("**", @t) if /RAC875_mira1_c35013/;
		my $id = join(":", @t[0,1]);
		my $k9_id = $1 if /(wsnp\S+)/;
		my $bristol = $1 if /(BS\d+)/;
		my $k;
		if(defined $k9_id and defined $bristol){$k=join(",",($k9_id, $bristol))}
		elsif(defined $k9_id){$k=$k9_id}
		elsif(defined $bristol){$k=$bristol}
		else{$k=""}
		my @alleles  =get_allele(@t[2,3]);
		$return{$id} = join("\t",($id, @alleles, $k) );
	}
	close IN;
	return %return;	
}

sub get_allele
{
	my $t= join("",@_);
	$t=~s/[^ATGC]//g;
	my %h = map{$_, 1} (split//,$t);
	return sort{$a cmp $b}keys %h;
}



sub read_fasta
{
	my $file = shift;
	my $flank = shift;
	my %return;
	open (IN, $file ) or die;
	my $id;
	while(<IN>)
	{
		chomp;
		next if /^\s+$/;
		if(/^>(\S+)/)
		{
			my $line = $1;
			my @t=split /:/,$line;
			my $len = ($t[2]-$t[1])>=$flank?$flank:($t[2]-$t[1]);
			$id = join(":", @t[0,1]);
			$return{$id}->[0] = $len;
			next;
		}
		else
		{
			$return{$id}->[1] .= $_
		}
	}
	close IN;
	return %return;
}



