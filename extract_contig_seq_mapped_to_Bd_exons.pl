#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

my $sUsage = qq(
perl $0
<cotig seq fasta file>
<blat result file>
<gff3 file>
<output file>
);
die $sUsage unless @ARGV >= 4;
my ($ctg_file, $blat_file, $gff_file, $out_file) = @ARGV;

print_time_comment("Read contig file ...");
my %contig = read_contig_file($ctg_file);

print_time_comment("Read blat file ...");
my %blat = read_blat_file($blat_file);

print_time_comment("Read gff file ...");
my %gff = read_gff_file($gff_file);

print_time_comment("Calculating_regions_mapped_to_exon ...");
my %regions_mapped_to_exons = calculate_regions_mapped_to_exon(\%contig, \%blat, \%gff);

output(\%contig, \%gff, \%blat, \%regions_mapped_to_exons, $out_file);

# Subroutines
sub print_time_comment
{
	my $comment = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $comment,"\n";
}

sub output
{
	my ($contig_ref, $gff_ref, $blat_ref, $mapped_region_ref, $out_file) = @_;
	open (OUT, ">$out_file") or die "can't open file $out_file\n";
	foreach my $ctg_id (keys %{$contig_ref})
	{
		my $count = 0;
#		if (not exists $blat_ref->{$ctg_id})
#		{
#			print OUT $ctg_id,"\t",'Not mapped to any genic regions',"\n";
#			next;
#		}		
		foreach my $chrid (keys %{$blat_ref->{$ctg_id}})
		{
			foreach (sort {$a->[0]<=>$b->[0]} @{$blat_ref->{$ctg_id}->{$chrid}})
			{
				my ($ctg_start, $ctg_end, $gene_start, $gene_end, $evalue) = @$_;
				#print "\t\t", join("\t", @$_),"\n";
				my ($frag_start, $frag_end);
				foreach ($ctg_start..$ctg_end)
				{
					if (vec($mapped_region_ref->{$ctg_id}, $_, 1) == 1)
					{
						$frag_start = $_ if vec($mapped_region_ref->{$ctg_id}, $_-1, 1) == 0 or $_ == $ctg_start;
						if ($_ == $ctg_end)
						{
							$frag_end = $_;
							$count++;
							if (($frag_end-$frag_start+1) >= 100)
							{
								print OUT '>', $ctg_id, '_', $count,"\n";
								print OUT substr($contig_ref->{$ctg_id}, $frag_start-1, ($frag_end-$frag_start+1)),"\n" ;								
							}

						}
					}
					else
					{
						next if $_ == $ctg_start;
						if (vec($mapped_region_ref->{$ctg_id}, $_-1, 1) == 1 )
						{
							$frag_end = $_-1;
							$count++;
							if (($frag_end-$frag_start+1) >= 100)
							{
								print OUT '>', $ctg_id, '_', $count,"\n";
								print OUT substr($contig_ref->{$ctg_id}, $frag_start-1, ($frag_end-$frag_start+1)),"\n" ;								
							}

						}
					}
				}
			}
		}
	}
	close OUT;
}

sub read_contig_file
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
			$id =~ s/^>//;
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

sub read_blat_file
{
	my $file = shift;
	open (IN,"$file") or die "can't open file $file\n";
	my %return_hash;
	my $evalue_cutoff = 1e-10;
	my $debug = 1;
	while (<IN>)
	{
		next if /^\s+$/;
		my @line_data = split /\s+/,$_;
		next if $line_data[-2] > $evalue_cutoff;
		my ($ctg_id, $chrid, $ctg_start, $ctg_end, $gene_start, $gene_end, $evalue) = @line_data[0..1, 6..10];
#		$chrid =~ s/(^LOC.*?)\|.*/$1/;
		print 'blat chr id: ', $chrid, "\n" if $debug;		$debug = 0;
		next unless $chrid =~ /bd\d/i;
		push @{$return_hash{$ctg_id}{$chrid} }, [$ctg_start, $ctg_end, $gene_start, $gene_end, $evalue];
	}
	close IN;
	return %return_hash;
}

sub read_gff_file
{
	my $file = shift;
	open (IN,"$file") or die "can't open file $file\n";
	my %exon_vector;
	my $debug = 1;
	while (<IN>)
	{
		next if /^\s+$/ or /^\#/;
		next unless /^Bd\d/;
		chomp;
		my @line_data = split /\s+/, $_;
		my ($chrid, $type, $start, $end, $description) = @line_data[0, 2..4, 8];
		print 'gff chr id: ', $chrid,"\n" if $debug; $debug = 0;
		next unless $type =~ /mRNA|exon/;
		if ($type =~ /mRNA/)
		{
			if ($description =~ /Name=(.*?)\;/)
			{
				
				push @{$exon_vector{$chrid}->[0]}, [$1, $start, $end];
				$exon_vector{$chrid}->[1] = '';
			}
		}
		else
		{
			foreach ($start .. $end)
			{
				vec($exon_vector{$chrid}->[1], $_, 1) = 0b1;
			}
		}
	}
	close IN;
	return %exon_vector;	
}

sub calculate_regions_mapped_to_exon
{
	my ($ctg_ref, $blat_hashref, $gff_hashref)= @_;
	my %return_hash;
	foreach my $ctg_id (keys %$ctg_ref)
	{
		my $mapped_region = '';
		my $overlap_region = '';
		foreach my $chrid (keys %{$blat_hashref->{$ctg_id}})
		{
			foreach (@{$blat_hashref->{$ctg_id}->{$chrid}})
			{
				my $flag = 1;
				my ($ctg_start, $ctg_end, $gen_start, $gen_end) = @$_;
				foreach my $pos ($ctg_start..$ctg_end)
				{
					if (vec($mapped_region, $pos,1) == 1)
					{
						vec($overlap_region, $pos,1) = 0b1;
						vec($mapped_region, $pos,1) = 0b0;
					}
					else
					{
						if (vec($overlap_region, $pos,1) == 1)
						{
							vec($mapped_region, $pos,1) = 0b0;
						}
						else
						{
							vec($mapped_region, $pos, 1) = vec($gff_hashref->{$chrid}->[1], $pos+$gen_start-$ctg_start, 1);
							vec($mapped_region, $pos-1, 1) = 0b0 if vec($mapped_region, $pos-1,1) == 1 and $flag == 1;
						}
					}
					$flag=0;
				}
			}			
		}
		$return_hash{$ctg_id} = $mapped_region;
	}
	return %return_hash;
}











