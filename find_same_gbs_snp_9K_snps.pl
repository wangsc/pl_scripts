#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<gbs_tag_blat_9K_contig output>
<gbs_snp_genotype>
<9K_genotype>
<9K_SNP_ID_NAME>
<9k_contig_start_end.out>
<output>
);
die $sUsage unless @ARGV >= 5;
my ($blat_out, $gbs_geno, $k9_geno, $id_name_file, $ctg_start_end_file, $out_file) = @ARGV;

my %k9_ctg_start = read_start_end_file($ctg_start_end_file);
my ($gbs_to_9k, $gbs_offset) = read_blat_file($blat_out);
my ($gbs_genotype, $ctg_snp) = read_gbs_genotype($gbs_geno);
my ($k9_genotype, $k9_ctg_snps, $id_name) = read_9k_genotype($k9_geno, $id_name_file);

# output
my @accs;
open (OUT, ">$out_file") or die;
foreach my $gbs_ctg (keys %$gbs_to_9k)
{
	next unless exists $ctg_snp->{$gbs_ctg};
	my @gbs_snps = @{ $ctg_snp->{$gbs_ctg}};
	
	foreach my $gbs_snp (@gbs_snps)
	{
		@accs = keys %{$gbs_genotype->{$gbs_snp}} unless @accs;
		#print STDERR $gbs_snp, "\n" unless exists $gbs_to_9k{$gbs_ctg};
		next unless exists $gbs_to_9k->{$gbs_ctg};
		my ($gbs_ctg, $gbs_snp_pos) = split /:/, $gbs_snp;
		my @k9_ctgs = @{$gbs_to_9k->{$gbs_ctg}};
		print OUT join("\t", ($gbs_snp, (map{$gbs_genotype->{$gbs_snp}{$_}}@accs))), "\n";
		foreach my $ctg (@k9_ctgs)
		{
			next unless exists $k9_ctg_snps->{$ctg};
			next unless exists $k9_ctg_start{$ctg};
			#print STDERR join("\t",($ctg, $k9_ctg_start{$ctg})), "\n";
			my @k9_snps = @{$k9_ctg_snps->{$ctg}};
			foreach my $k9_snp( @k9_snps)
			{
				my $new_gbs_pos = abs($gbs_offset->{$gbs_ctg}{$ctg}) + 1 + ($gbs_offset->{$gbs_ctg}{$ctg}<0?-1:1)*($gbs_snp_pos-1);
				my $orig_9k_pos = $1 if $id_name->{$k9_snp} =~ /_(\d+)$/;
				my $new_9k_pos = $orig_9k_pos - $k9_ctg_start{$ctg} + 1;
				print STDERR join("\t", ($gbs_snp, $new_gbs_pos, $k9_snp, $new_9k_pos)), "\n";
				next unless $new_gbs_pos == $new_9k_pos;
				print OUT join("\t",($k9_snp, (map{exists $k9_genotype->{$k9_snp}{$_}?$k9_genotype->{$k9_snp}{$_}:'?' }@accs))), "\n";
			}
		}
	}	
}
print OUT "\t", join("\t", @accs), "\n";
close OUT;

sub read_blat_file
{
	my $file = shift or die;
	open (IN, $file) or die;
	my %return;
	my %gbs_offset;
	while(<IN>)
	{
		chomp;
		my @t = split /\t/, $_;
		next unless $t[6] == 1;
		push @{$return{$t[0]}}, $t[1];
		my $offset = $t[8]-$t[6];
		$offset *= -1 if $t[8]>$t[9];
		$gbs_offset{$t[0]}{$t[1]} = $offset;
	}
	
	close IN;
	return (\%return, \%gbs_offset);
}

sub read_start_end_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		chomp;
		my @t = split /,/, $_;
		$return{$t[0]} = $t[1];
	}
	close IN;
	return %return;	
}



sub read_gbs_genotype
{	
	my $file = shift or die;
	open (IN, $file) or die;
	my %return;
	my %snp_to_tag;
	my @id_arr;
	my $n = 0;
	while(<IN>)
	{
		chomp;
		$n++;
		my @t = split /\s+/,$_; 
		if($n==1)
		{
			@id_arr = @t;
			next;
		}
		
		foreach (1..$#t)
		{
			my $acc = $id_arr[$_];
			$return{$t[0]}{$acc} = $t[$_];
		}
		
		my $ctg = $1 if $t[0] =~ /(\S+):/;
		push @{$snp_to_tag{$ctg}}, $t[0];
	}
	close IN;
	return (\%return, \%snp_to_tag);	
}



sub read_9k_genotype
{	
	my $file = shift or die;
	my $id_name_file = shift;
	my %id_name = read_id_name($id_name_file);
	my %name_id = reverse %id_name;
	open (IN, $file) or die;
	my %return;
	my %k9_ctg_snps;
	my @snp_arr;
	my $n = 0;
	while(<IN>)
	{
		chomp;
		$n++;
		my @t = split /\s+/,$_; 
		if($n==1)
		{
			@snp_arr = @t;
			foreach (1..$#snp_arr)
			{
				next unless (exists $id_name{$snp_arr[$_]});
				my $ctg = $1 if $id_name{$snp_arr[$_]} =~ /wsnp_(\S+)_(\d+)$/;
				if(not defined $ctg)
				{
					$snp_arr[$_] = 'NA';
				}
				else
				{
					push @{$k9_ctg_snps{$ctg}}, $snp_arr[$_];
				}
			}
			
			next;
		}
		
		foreach (1..$#t)
		{
			my $snp = $snp_arr[$_];
			next if $snp eq 'NA';
			next unless exists $id_name{$snp};
			$return{$snp}{$t[0]} = $t[$_];
		}		
	}
	close IN;
	return (\%return, \%k9_ctg_snps, \%id_name);
}

sub read_id_name
{
	my $file = shift;
	open (IN, $file) or die;
	my %return;
	while(<IN>)
	{
		chomp; 
		my @t=split/\s+/,$_;
		$return{$t[0]} = $t[1];
	}
	close IN;
	return %return;
}



