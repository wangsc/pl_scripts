#!/usr/bin/perl -w
use strict;
die <<END if (@ARGV < 2 );
**********************************************************************************************
To use this script, excel files have to be saved as CSV file at first.
How to do this: open the file with excel -> select file -> save as,
then select CSV type to save.

Usage: 
Perl $0 [Genotype file] [guide file]

Note: 
Genotype file (CSV type); format is: snp_index, snp_name, genotype_acc1, genotype_acc2, ...  
convert genotype guide file (CSV type); format looks like: 34,wsnp_BE403597B_Ta_1_1,AB to BB, [A/G], TOP(or BOT)
* First step: 
  SNPs in the guide file with information like "AB to AA" will be converted;
  If no AB in file, script will just go to the second step.
* Second step:
  When [A/C] corresponds to TOP, that means AA genotype in the input data file should be converted to nucleotide AA, 
  and BB genotype should be converted to nucleotide CC.
  When [A/C] corresponds to BOT, that means AA genotype in the input file should be converted to nucleotide CC, 
  and BB genotype should be converted to nucleotide AA.

The new file would be named as [original_genotype_file]_converted.csv in the same folder

--- Updated on Jan.19.2012, modified by Chao Jan 20, 2012
convert all BOT to TOP in the Good SNP guide table, and also convert the SNP accordingly.  
For example, SNP index 2 has [T/G] on the SNP column, and BOT on the ILMN Strand column, 
we will change that to [A/C] and TOP. Append and write the changes to a new column. 
Then we convert all the A/B to nucleotides following only one rule,
When [A/C] corresponds to TOP, that means AA genotype in the input data file should 
be converted to nucleotide AA, and BB genotype should be converted to nucleotide CC.
--- end
***********************************************************************************************
END

my ($genotype_file, $guide_file) = @ARGV;

my %guide = read_guide_file($guide_file);
convert_genotype($genotype_file, \%guide);

# Subroutines
sub read_guide_file
{
	my $file = shift;
	my $out_file = $file . "_new.csv";
	my %return_hash;
	open (IN, $file) or die $!;
	open (OUT, ">$out_file") or die $!;
	while (<IN>)
	{
		next if /^\s+$/;
		s/\s+$//;
		if(/Index/){print OUT $_, "\n"; next} #skip Index line
		# 1	wsnp_AJ612027A_Ta_2_1	good	[A/G]	TOP
		# 2	wsnp_AJ612027A_Ta_2_5	poly	[T/G]	BOT
		# 3	wsnp_be352570A_Ta_1_1	AB to AA	[A/C]	TOP
		chomp;
		s/\"//g;
		s/\s+$//;
		my @data = split /,/,$_;
		my $index = $data[0];
		my ($old_geno, $new_geno);
		if ($data[2] =~ /to/)
		{
			($old_geno, $new_geno) = split /to/, $data[2];
			$old_geno =~ s/\s//g;
			$new_geno =~ s/\s//g;		
			$return_hash{$index}{$old_geno} = $new_geno;			
		}
		else # 08.31.2011 added
		{
			my $geno = $data[3];
			$geno =~ s/[^ATGC]//g;
			$return_hash{$index}{'AB'} = $geno;
		}
		my ($AA, $BB) = get_genotype(@data[-2,-1]);
		print join("\t", ($AA, $BB)),"\n" if $data[1] eq "wsnp_BE442776A_Ta_2_1";
		$return_hash{$index}{'AA'} = $AA if defined $AA;
		$return_hash{$index}{'BB'} = $BB if defined $BB;
		#print $return_hash{$index}{'AA'}, " * ", $return_hash{$index}{'BB'}, "\n" if $index ==109;
		# added Jan.19.2012
		if($data[-1] eq "BOT")
		{
			my $t = $data[-2];
			$t =~ s/[^ATGC]//g;
			$t=~ tr/[ATGC]/[TACG]/;
			$t = "[". join("\/", (split//, $t)). "]";
			print OUT join(",", (@data, $t, "TOP")),"\n";
		}
		else
		{
			print OUT join(",", @data), "\n"
		}	
	}
	close IN;
	close OUT;
	return %return_hash;
}

sub get_genotype
{
	my $seq = shift;
	my $flag = shift; # Top or Bot
	my $bot_flag = ($flag =~ /bot/i)?1:0;
	$seq =~ s/[^ATGC]//g;
	$seq =~ tr/[ATGC]/[TACG]/ if $bot_flag; # added Jan.19.2012
	my @seqs = split //, $seq;
	return ($seqs[0].$seqs[0], $seqs[1].$seqs[1])
}


sub convert_genotype
{
	my $file = shift;	
	my $guide_ref = shift;
	open (IN, $file) or die $!;
	my $out_file = $file . '_converted.csv';
	open (OUT, ">$out_file") or die "can't open file $out_file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		s/\"//g;
		if (/^Index/)
		{
			print OUT $_;
			next;
		}
		chomp;
		my @data = split /,/, $_;
		if (exists $guide_ref->{$data[0]})
		{
			foreach my $index (2..$#data)
			{
				#print $data[$index], " * ", $guide_ref->{$data[0]}->{$data[$index]}, "\n" if $data[0] == 109;
				if($data[$index] eq "AB")
				{
					$data[$index] =~ s/\s//g;
					$data[$index] = $guide_ref->{$data[0]}->{$data[$index]} if exists $guide_ref->{$data[0]}->{$data[$index]};
					$data[$index] = $guide_ref->{$data[0]}->{$data[$index]} if exists $guide_ref->{$data[0]}->{$data[$index]};
					next;		
				}
				$data[$index] =~ s/\s//g;
				$data[$index] = $guide_ref->{$data[0]}->{$data[$index]} if exists $guide_ref->{$data[0]}->{$data[$index]};
			}
		}
		print OUT join(',', @data),"\n";
	}
	close IN;
	close OUT;
}

# added Jan.19.2012
sub reverse_comp
{
	my $seq = shift;
	$seq =~ tr/[ATGC]/[TACG]/;
	return reverse $seq;
}

