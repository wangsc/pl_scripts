#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<shiaoman file>
<ksu sample sheet>
<list to be deleted>
<kian csv files>
>
);
die $sUsage unless @ARGV >= 4;
my($shiaoman_file, $sample_sheet, $del_list, @csv_files) = @ARGV;

my ($sample_names, $snp_index_ref) = read_shiaoman_file($shiaoman_file);
my %sample_name_id = read_sample_sheet($sample_sheet);
my ($kian_results, $index_name) = read_kian_csv_files(@csv_files);
my %delete_list = read_delete_list($del_list);

open (ID, ">Names_ID.out") or die $!;
print ID "--\t--\t";
print ID join("\t", @$sample_names),"\n";
print ID "--\t--";
foreach (@$sample_names)
{
	my $id = (exists $sample_name_id{$_})?$sample_name_id{$_}:'--';
	print ID "\t", $id;
}
print ID "\n";
close ID; 

foreach my $index (keys %$kian_results)
{
	next if exists $snp_index_ref->{$index};
	next if exists $delete_list{$index};
	print $index,"\t", $index_name->{$index};
	foreach my $sample (@$sample_names)
	{
		$sample =~ s/\s+$//;
		unless (exists $sample_name_id{$sample}){print STDERR $sample,"**\n"; next;}
		next unless (exists $sample_name_id{$sample});
		my $genotype = (exists $kian_results->{$index}{$sample_name_id{$sample}})?$kian_results->{$index}{$sample_name_id{$sample}}:'--';
		print "\t", $genotype;
	}
	print "\n";
}

# Subroutines
sub read_delete_list
{
	my $file = shift;
	my %return;
	open (IN, $file) or die $!;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		next if /\D/;
		$return{$_} = 1;
	}
	close IN;
	return %return;
}


sub read_shiaoman_file
{
	my $file = shift;
	my (@names, %snp_index);
	open (IN, "$file") or die $!;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		s/\"//g;
		my @data = (/,/)?split /,/,$_:split /\t/, $_;
		if (/Index/)
		{
			@names = @data[2..$#data];
#			map{print STDERR $_,"*\n"} @names;
			next;
		}
		
		next unless $data[0] =~ /\d/;
		$snp_index{$data[0]} = 1;
	}
	close IN;
	return (\@names, \%snp_index);
}

sub read_sample_sheet
{
	my $file = shift;
	my %return;
	open (IN, "$file") or die $!;
	while (<IN>)
	{
		chomp;
		next if /^\s+$/;
		s/\"//g;
		my @data = split /\t/,$_;
		my $name = shift @data;
		$name =~ s/\s+$//;
#		print '$name ', $name,"\n";
		my $id = pop @data;
#		print '$id ', $id,"\n";
		#if (/Heyne/){print STDERR join("**", ($name, $id)),"\n"; die}
		$return{$name} = $id;
	}
	close IN;
	return %return;
}

sub read_kian_csv_files
{
	my @files = @_;
	my %records;
	my %index_name;
	my %name_index;
	foreach (@files)
	{
		open (IN, "$_") or die $!;
		while (<IN>)
		{
			next if /^\s+$/;
			next unless /^wsnp/;
			my @data = (/,/)?split /,/,$_:split /\s+/, $_;
			my ($snp_name, $sample_id, $snp_index, $genotype) = @data[0, 1, 2, 4];
			#print STDERR join("***", ($snp_name,$snp_index)),"\n" if /wsnp_Ex_rep_c67838_66536117/;
			$index_name{$snp_index} = $snp_name if $snp_index =~ /\d/;
			$name_index{$snp_name} = $snp_index if $snp_index =~ /\d/;
			die join('**', (@data)),"\n" unless defined $genotype;
			$snp_index = $name_index{$snp_name};
			$genotype = $data[3] if (length $genotype) > 2;
			$records{$snp_index}{$sample_id} = $genotype;
		}
	}
	close IN;
	return (\%records,\%index_name);
}





