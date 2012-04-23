#!/usr/bin/perl -w
use strict;
my $sUsage = "perl $0 <KSU_SNP_AllSNPAllSamples_PCRSelection.txt> <deletion line list>\n";
die $sUsage unless @ARGV >= 2;

my($bigfile, $del_list) = @ARGV;

my $cs_sample = '5904367093_R11C01';

my %del_sample = read_deletion_list($del_list);
$del_sample{$cs_sample} = 'CS';

read_big_file_and_output(\%del_sample, $cs_sample, $bigfile);


# 
sub read_deletion_list
{
	my $file = shift;
	my %return;
	open (IN, "$file") or die "can't open file $file \n";
	while (<IN>)
	{
		next if /^\s+$/;
		my @t = split /\s+/, $_;
		$return{$t[1]} = $t[0];
	}
	close IN;
	return %return;
}

sub read_big_file_and_output
{
	my ($delref, $cs_sample, $file) = @_;
	open (IN, "$file") or die "can't open file $file \n";
	my %index_hash = map{$_, []}(keys %$delref);
	my @samples = sort{$b cmp $a}keys %$delref;
	while (<IN>)
	{
		chomp;
		next if /^\s+$/;
		my @t = split /\t/, $_;
		#print 'scalar @t ', scalar @t,"\n";
		if (/^Index/)
		{
			foreach my $index (0 ..$#t)
			{
				my $id = $t[$index];
				if ( $id =~ /Theta/)
				{
					my ($sample) = split /\./, $id;
					#print '$id: ', $id,"\t", '$sample: ', $sample,"\n";
					$index_hash{$sample}->[0] = $index;
				}
				if ($id =~ /R$/)
				{
					my ($sample) = split /\./, $id;
					$index_hash{$sample}->[1] = $index;
					#print '$id: ', $id,"\t", '$sample: ', $sample,"\n";
				}
			}
			#exit;
			next;
		}
		else
		{
			print join("\t",(@t[0, 1])),"\n";
			print join("\t", qw(Acc Theta R) ),"\n";
			foreach (@samples)
			{
				my ($theta_index, $R_index) = @{$index_hash{$_}};
				print $delref->{$_}, "\t", join("\t", ($t[$theta_index], $t[$R_index])),"\n";
			}
			print "\n\n";
		}
	}
	close IN;	
}