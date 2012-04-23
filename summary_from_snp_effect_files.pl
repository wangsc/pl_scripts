#!/usr/bin/perl -w
use strict;
#
# 
my $sUsage = qq(
perl $0
<fasta file used for blast>
<snp effect files>
);
die $sUsage unless @ARGV >= 2;
my $fasta_file = shift;
my @files = @ARGV;
#
my %fasta_hash = read_fasta_file($fasta_file);
my %snpeffect_hash = read_snpeffect_files(@files);
summarize_and_output(\%fasta_hash, \%snpeffect_hash);

# subroutiens

sub summarize_and_output
{
	my ($fasta_hashref, $snpeffect_hashref) = @_;
	my %summary;
	my $summary_out = 'snp_effect_summary.out';
	if (-e $summary_out) {$summary_out .= '.new'}
	open (SU, ">$summary_out") or die "can't open file $summary_out\n";
	my $nohit_fasta_out = 'sequences_no_blast_hit.fasta';
	if (-e $nohit_fasta_out) {$nohit_fasta_out .= '.new'};
	open (NH, ">$nohit_fasta_out") or die "can't open file $nohit_fasta_out\n";
	open (IG, ">integrated_snpeffect.out") or die "can't open file integrated_snpeffect.out\n";
	foreach my $id (keys %$snpeffect_hashref)
	{
		print IG $id, "\t", join("\t", @{$snpeffect_hashref->{$id}}),"\n";
		if (scalar @{$snpeffect_hashref->{$id}} == 1) # means no_hit
		{
			print NH '>'.$id,"\n", $fasta_hash{'>'.$id},"\n";
			$summary{'no_hit'}++;
			next;
		}
		my $aa_change = $snpeffect_hashref->{$id}->[1];
		print join("\t", @{$snpeffect_hashref->{$id}}),"\n" unless defined $aa_change;
		my ($ref_aa, $alt_aa) = $aa_change =~ /(\S)->(\S)/;
		if ($ref_aa eq '?'){$summary{'unknown'}++; next}
		if ($ref_aa eq $alt_aa)
		{
			$summary{'syn'} ++;
		}
		else
		{
			my $ref_prop = AA_prop($ref_aa); 
			my $alt_prop = AA_prop($alt_aa);
			$summary{join('<->', sort{$a cmp $b}($ref_prop, $alt_prop))}++;
		}
	}
	while (my @data = each %summary)
	{
		print SU join("\t", @data),"\n";
	}
	close SU;
	close NH;
}

sub AA_prop
{
	my $aa = shift;
	return 'stop' if $aa eq '_';
	my %basic = map{$_, 1} qw(K R H);
	my %acidic = map{$_, 1} qw(D E);
	my %polar = map{$_, 1} qw(G S T C Y N Q);
	my %nonpolar = map{$_, 1} qw(A V L I P F W M);
	return 'basic' if exists $basic{$aa};
	return 'acidic' if exists $acidic{$aa};
	return 'polar' if exists $polar{$aa};
	return 'nonpolar' if exists $nonpolar{$aa};
}


sub read_fasta_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file\n";
	my $id;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		if (/^>/)
		{
			$id = $_;
			$return_hash{$id} = '';
			next;
		}
		else
		{
			$return_hash{$id} .= $_;
		}
	}
	close IN;
	return %return_hash; 
}

sub read_snpeffect_files
{
	my @files = @_;
	my %snpeffect_hash;
	foreach my $file (@files)
	{
		open (IN,"$file") or die "can't open file $file \n";
		while (<IN>)
		{
			next if /^\s+$/;
			chomp;
			my @line_data = split /\s+/, $_;
			print $_,"\n" if (/no_hit/ and scalar @line_data >= 3);
			my $key = shift @line_data;
			if (not exists $snpeffect_hash{$key} or (scalar @{$snpeffect_hash{$key}}) == 1)
			{
				$snpeffect_hash{$key} = [@line_data];
			}
		}
		close IN;
	}
	return %snpeffect_hash;
}
