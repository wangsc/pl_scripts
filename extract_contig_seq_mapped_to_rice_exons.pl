#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<cotig seq csv file>
<blat result file>
<rice gff3 file>
<output file>
);
die $sUsage unless @ARGV >= 4;
my ($csv_file, $blat_file, $gff_file, $out_file) = @ARGV;
print_time_comment("Read csv file ...");
my %csv = read_csv_file($csv_file);
print_time_comment("Read blat file ...");
my %blat = read_blat_file($blat_file);
print_time_comment("Read gff file ...");
my %gff = read_gff_file($gff_file);
print_time_comment("Calculating_regions_mapped_to_exon ...");
my %regions_mapped_to_exons = calculate_regions_mapped_to_exon(\%csv, \%blat, \%gff);
output(\%csv, \%blat, \%regions_mapped_to_exons, $out_file);

# Subroutines
sub print_time_comment
{
	my $comment = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $comment,"\n";
}

sub output
{
	my ($csv_ref, $blat_ref, $mapped_region_ref, $out_file) = @_;
	open (OUT, ">$out_file") or die "can't open file $out_file\n";
	foreach my $snpid (keys %{$csv_ref})
	{
		if (not exists $blat_ref->{$snpid})
		{
			print OUT $snpid,"\t",'Not mapped to any rice genic regions',"\n";
			next;
		}
		my ($ctg_seq, $snppos, $flank_seq) = @{$csv_ref->{$snpid}};
		$ctg_seq = lc($ctg_seq);
		foreach (1..(length $ctg_seq))
		{
			my $base =substr($ctg_seq, $_-1, 1);
			substr($ctg_seq, $_-1, 1) = uc($base) if vec($mapped_region_ref->{$snpid}, $_, 1) == 1;
		}
		substr($ctg_seq, $snppos-1, 1) = $flank_seq;
		print OUT $snpid,"\t", $ctg_seq,"\n";
		print $snpid,"\n";
		foreach my $aliasname (keys %{$blat_ref->{$snpid}})
		{
			print "\t", $aliasname,"\n";
			foreach (sort {$a->[0]<=>$b->[0]} @{$blat_ref->{$snpid}->{$aliasname}})
			{
				print "\t\t", join("\t", @$_),"\n";
			}
		}
	}
	close OUT;
}

sub read_csv_file
{
	my $file = shift;
	open (IN,"$file") or die "can't open file $file\n";
	my %return_hash;
	while (<IN>)
	{
		next unless /wsnp/;
		next if /^\s+$/;
		s/\"//g;
		chomp;
		my @data = split /,/,$_;
		my ($id, $ctg_seq, $snp_pos, $flank_seq) = @data;
	#	print STDERR '$flank_seq _ori: ', $flank_seq,"\n";
		$flank_seq = $1 if $flank_seq =~ /.*(\[\w)/;
		$flank_seq = $flank_seq . ']';
	#	print STDERR '$flank_seq: ', $flank_seq,"\n";
		$return_hash{$id} = [$ctg_seq, $snp_pos, $flank_seq];
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
	while (<IN>)
	{
		next if /^\s+$/;
		my @line_data = split /\s+/,$_;
		next if $line_data[-2] > $evalue_cutoff;
		my ($snpid, $geneid, $ctg_start, $ctg_end, $gene_start, $gene_end, $evalue) = @line_data[0..1, 6..10];
		$geneid =~ s/(^LOC.*?)\|.*/$1/;
#		print '$geneid: ', $geneid, "\n";
		push @{$return_hash{$snpid}{$geneid} }, [$ctg_start, $ctg_end, $gene_start, $gene_end, $evalue];
	}
	close IN;
	return %return_hash;
}

sub read_gff_file
{
	my $file = shift;
	open (IN,"$file") or die "can't open file $file\n";
	my %gene_start;
	my %cds_vector;
	while (<IN>)
	{
		next if /^\s+$/ or /^\#/;
		next unless /^Chr\d/;
		chomp;
		my @line_data = split /\s+/, $_;
		my ($type, $start, $end, $description) = @line_data[2..4, 8];
		next unless $type =~ /mRNA|CDS/;
		if ($type =~ /mRNA/)
		{
#			print join("\t", @line_data[2..4, 8]),"\n";
#			print '$description: ', $description,"\n";
			if ($description =~ /ID=(.*?)\;.*Alias=(.{14})/)
			{
#				print '$1: ', $1, '$2: ',$2,"\n";
				$gene_start{$1} = [$start, $2];
			}
			next;
		}
		else
		{
			if ($description =~ /Parent\=(.*)/)
			{
				print '$1: ', $1, "*\n" if (not exists $gene_start{$1});
				my ($genstart, $aliasname) = @{$gene_start{$1}};
				$cds_vector{$aliasname} = '' unless exists $cds_vector{$aliasname};
				foreach ($start-$genstart+1..$end-$genstart+1)
				{
					vec($cds_vector{$aliasname}, $_, 1) = 1;
				}
			}
			else{print "What is this: $type $description","\n"}
		}
	}
	close IN;
	return %cds_vector;	
}

sub calculate_regions_mapped_to_exon
{
	my ($csv_hashref, $blat_hashref, $gff_hashref)= @_;
	my %return_hash;
	foreach my $snpid (keys %$csv_hashref)
	{
		my ($ctg_seq, $snppos, $flank_seq) = $csv_hashref->{$snpid};
		my $ctg_length = length $ctg_seq;
		
		my $mapped_region = '';
		my $overlap_region = '';
		foreach my $aliasname (keys %{$blat_hashref->{$snpid}})
		{
			foreach (@{$blat_hashref->{$snpid}->{$aliasname}})
			{
				my $flag = 1;
				my ($ctg_start, $ctg_end, $gen_start, $gen_end) = @$_;
				foreach my $pos ($ctg_start..$ctg_end)
				{
					if (vec($mapped_region, $pos,1) == 1)
					{
						vec($overlap_region, $pos,1) = 1;
						vec($mapped_region, $pos,1) = 0;
					}
					else
					{
						if (vec($overlap_region, $pos,1) == 1)
						{
							vec($mapped_region, $pos,1) = 0;
						}
						else
						{
							vec($mapped_region, $pos, 1) = vec($gff_hashref->{$aliasname}, $pos+$gen_start-$ctg_start, 1);
							vec($mapped_region, $pos-1, 1) = 0 if vec($mapped_region, $pos-1,1) == 1 and $flag == 1;
						}
					}
					$flag=0;
				}
				
			}			
		}
		$return_hash{$snpid} = $mapped_region;
	}
	return %return_hash;
}











