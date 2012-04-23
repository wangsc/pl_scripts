#!/usr/bin/perl -w
use strict;
my $sUsage = qq(
perl $0
<ortholist rice>
<ortholist bd>
<pasa_chr3A_cdna.pasa_assemblies.denovo_transcript_isoforms.gff3>
);

die $sUsage unless @ARGV >= 3;
my ($rice_list, $bd_list, $gff_file) = @ARGV;
my %rice_record = read_ortho_file($rice_list);
my %bd_record =  read_ortho_file($bd_list);
my %rice_com; my %bd_com;
foreach (keys %rice_record)
{
	#print 'keys %rice_record: ', $_,"\n";
	if (exists $bd_record{$_})
	{
		$rice_com{$_} = $rice_record{$_};
		$bd_com{$_} = $bd_record{$_};
	}
}
print 'Rice:',"\t", join("\t", (count(values %rice_com))),"\n";
print 'BD:',"\t", join("\t", (count(values %bd_com))),"\n";

#
my ($iso_gene_ref, $gene_ctg_ref)  = read_gff_file($gff_file);
my %rice_gene_ctgs = detect_multiple_match(\%rice_record, $iso_gene_ref, $gene_ctg_ref);
my %bd_gene_ctgs = detect_multiple_match(\%bd_record, $iso_gene_ref, $gene_ctg_ref);
output(\%rice_gene_ctgs);
output(\%bd_gene_ctgs);


######

sub output
{
	my $hashref = shift;
	print join("\t", qw(ID Contigs Num_ctg Same_ctg)  ),"\n";
	my ($num_same, $num_diff, $num_one) = (0, 0 , 0);
	foreach (keys %$hashref)
	{
		my @ctgs = @{$hashref->{$_}};
		my %h = map{$_, 1} @ctgs;
		my $num_ctg = scalar @ctgs;
		my $unique_ctg = scalar keys %h;
		my $same_ctg;
		if ($num_ctg == 1)
		{
			$same_ctg = 1;
		}
		else
		{
			$same_ctg = $unique_ctg>1?0:1; 
		}
		$num_one++ if $num_ctg==1;
		$num_same++ if $same_ctg==1 and $num_ctg>1;
		$num_diff++ if $same_ctg==0 and $num_ctg>1;
		print $_,"\t", join("**", (@ctgs)),"\t", $num_ctg, "\t", $same_ctg,"\n";
	}
	print '$num_one ',$num_one,"\t", '$num_same ', $num_same,"\n", '$num_diff ', $num_diff,"\n";
}

sub detect_multiple_match
{
	my ($rice, $iso_gene, $gene_ctg) = @_;
	my %record;
	foreach my $id (keys %$rice)
	{
		#print '$id ', $id,"\n", '$rice->{$id} ', $rice->{$id},"\n";
		$record{$rice->{$id}} = [] unless exists $record{$rice->{$id}};
		push @{$record{$rice->{$id}}}, $iso_gene->{$id};
	}
	#%record = map{$_, [get_unique($record{$_})]} keys %record;
	foreach (keys %record)
	{
		#print $_, "\t", $record{$_},"\n";
		$record{$_} = [get_unique($record{$_})]
	}
	my %return;
	foreach my $gene (keys %record)
	{
		my @ctgs;
		foreach (@{$record{$gene}})
		{
			push @ctgs, $gene_ctg->{$_};
		}
		$return{$gene} = [@ctgs];
		#print '$gene ', $gene,"\n";
	}
	return %return;
}

sub get_unique
{
	my $array_ref = shift;
	my %h = map {$_, 1}@$array_ref;
	return keys %h;
}

sub read_gff_file
{
	my $file = shift;
	my %iso_gene;
	my %gene_contig;
	open (IN, "$file") or die "$! $file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my ($gene, $iso_id) = $_=~/ID=(.*?)\-(.*?)\;/;
		$iso_gene{$iso_id} = $gene;
		my @data = split /\t/, $_;
		$gene_contig{$gene} = $data[0];
		#print '$iso_id ', $iso_id,"\n", '$gene ', $gene,"\n", 'contig_ID ', $data[0],"\n";
	}
	close IN;
	return (\%iso_gene, \%gene_contig);
}


sub count
{
	my @array = @_;
	my $total = scalar @array;
	my %hash = map{$_, 1} @array;
	my $unique = scalar keys %hash;
	return($total, $unique);
}

sub read_ortho_file
{
	my $file = shift;
	my %record;
	open (IN, "$file") or die "$! $file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		next if /Overlap/;
		chomp;
		my @data  = split /\s+/, $_;
		#print join(" **\n ", @data),"\n";
		my $id = $data[1];
		#print '**$id ', $id,"\n";
		if ($id =~ /(.*?)\|/)
		{
			$id = $1 if $data[1] =~ /^(.*?)\|/; 
		}
		$id = $1 if $id=~/(.*)\.\d/;
		#print '$id ', $id, "\n";
		$record{$data[0]} = $id;
	}
	close IN;
	return %record;
}