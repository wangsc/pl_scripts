#!/usr/bin/perl -w
use strict;

my $sUsage = qq"
Usage:
perl $0 <input file, 9KSNP_AMPanel_FinalReport_nohead> <outfile>

Input file format:
wsnp_AJ612027A_Ta_2_5 5904341053_R08C01 A C 0.5376 2 T G T G A B - - 1.0000
wsnp_AJ612027A_Ta_2_5 5904341053_R09C01 A C 0.5376 2 T G T G A B - - 1.0000
wsnp_AJ612027A_Ta_2_5 5904341053_R09C02 C C 0.5376 2 G G G G B B - - 1.0000
wsnp_AJ612027A_Ta_2_5 5904341053_R10C01 A C 0.5376 2 T G T G A B - - 1.0000
wsnp_AJ612027A_Ta_2_5 5904341053_R12C01 C C 0.5376 2 G G G G B B - - 1.0000
...

Output file format:
wsnp_AJ612027A_Ta_2_5 5904341053_R08C01 A C 0.5376 2 T G T G A B - - 1.0000 AC AA
wsnp_AJ612027A_Ta_2_5 5904341053_R09C01 A C 0.5376 2 T G T G A B - - 1.0000 AC AA
wsnp_AJ612027A_Ta_2_5 5904341053_R09C02 C C 0.5376 2 G G G G B B - - 1.0000 CC CC
wsnp_AJ612027A_Ta_2_5 5904341053_R10C01 A C 0.5376 2 T G T G A B - - 1.0000 AC AA
wsnp_AJ612027A_Ta_2_5 5904341053_R12C01 C C 0.5376 2 G G G G B B - - 1.0000 CC CC
";
die $sUsage unless @ARGV;
my $out_file = pop;
my @files = @ARGV;
my (%file_content, %count_allele);
read_files(\%file_content, \%count_allele, @files);
my %genotype = call_genotype(\%count_allele);
output(\%file_content, \%genotype, \%count_allele, $out_file);

sub output
{
	my ($file_content, $geno_ref, $count_allele, $out) = @_;
	open (OUT, ">$out") or die $!;
	my %flag = map{$_, 1}keys %{$file_content};
	foreach my $snpid (keys %{$file_content})
	{

		foreach my $content (@{$file_content->{$snpid}})
		{
			my @array = @$content;
			my ($allea, $alleb) = @array[1,2];
			my $allel_id = join('', (sort{$a cmp $b} ($allea, $alleb)));
			my $change_allele = $geno_ref->{$snpid}{$allel_id};
			print OUT $snpid,"\t";
			print OUT join("\t", (@array)),"\t";
			print OUT $allel_id,"\t", $geno_ref->{$snpid}{$allel_id},"\t";
 #  	print $snpid,"\t";
#  		print join("\t", ($sid, $allea, $alleb, $score)),"\n";
 # 		print $allel_id,"\t", $change_allele,"\t";

			if ($flag{$snpid})
			{
				while (my @pair = each %{$count_allele->{$snpid}})
				{
					print OUT join("_", (@pair)),"\t";
#					print  join("_", (@pair)),"\t";
				}
				$flag{$snpid} = 0;		
			}
			print OUT "\n";
#			print "\n";
		}
	}
	close OUT;
}

sub read_files
{
	my $return_a = shift;
	my $return_b = shift;
	my @files = @_;
	foreach my $file (@files)
	{
		open (IN, "$file") or die "can't open file $file\n";
		while (<IN>)
		{
		#wsnp_AJ612027A_Ta_2_5 5904341053_R08C01 A C 0.5376 2 T G T G A B - - 1.0000
		#wsnp_AJ612027A_Ta_2_5 5904341053_R09C01 A C 0.5376 2 T G T G A B - - 1.0000
			
			chomp;
			next if /^\s+$/;
			next unless /^wsnp/;
			my @line_data = split /\s+/, $_;
			next if $line_data[2] eq '-';
			my $snpid = shift @line_data;			
			push @{$return_a->{$snpid}}, [@line_data];
			my ($allele_a, $allele_b) = @line_data[1,2];
			$return_b->{$snpid}{join('',(sort{$a cmp $b}($allele_a, $allele_b)))}++;
		}
		close IN;
	}
}

sub call_genotype
{
	my $count_allele = shift;
	my %return_hash;
	foreach my $snpid (keys %{$count_allele})
	{
		my $num_genotype = scalar(keys %{$count_allele->{$snpid}});
		if ($num_genotype == 1 or $num_genotype == 3)
		{
			foreach (keys %{$count_allele->{$snpid}}){$return_hash{$snpid}{$_} = $_}
		}
		else
		{
			 my %genotype_change = check_genotype(keys %{$count_allele->{$snpid}});
			 foreach (keys %genotype_change )
			 {
			 	$return_hash{$snpid}{$_} = $genotype_change{$_};
			 }
		}
	}
	return %return_hash;
}

sub check_genotype
{
	my($genoa, $genob) = @_;
	my %return;
	my @geno_a = split //, $genoa;
	my @geno_b = split //, $genob;
	if ($geno_a[0] eq $geno_a[1])
	{
		$return{$genoa} = $genoa;
		$return{$genob} = ($geno_b[0] eq $geno_a[0])?$geno_b[1] . $geno_b[1]:$geno_b[0] . $geno_b[0];
	}
	else
	{
		$return{$genob} = $genob;
		$return{$genoa} = ($geno_b[0] eq $geno_a[0])?$geno_a[1] . $geno_a[1]:$geno_a[0] . $geno_a[0];
	}
	return %return;
}



