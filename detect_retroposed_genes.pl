#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $sUsage = qq(
perl $0 
<assemblies blast rice.cds file>
<rice.gff3>
<assemblies gff3>
<output>
);
die $sUsage unless @ARGV >= 4;
my ($blast_file, $rice_gff, $asmbl_gff, $output) = @ARGV;
open (OUT, ">$output") or die;

my ($rice_cds_structure, $alias_id) = read_rice_gff($rice_gff);
my ($asmbl_pos_to_rice, $query_start_end) = read_blast($blast_file);
my %asmbl_cds_structure = read_asmbl_gff($asmbl_gff, $asmbl_pos_to_rice, $query_start_end);

# check overlapping
foreach my $asmbl (keys %asmbl_cds_structure)
{
	my @retroposed_exons;
	my $rice_gene_alias = $asmbl_pos_to_rice->{$asmbl}->[0]; 
	my $rice_cds = $alias_id->{$rice_gene_alias};
	unless(defined $rice_cds){next ; print STDERR $rice_gene_alias, "\n";} 
	foreach my $exon (@{$asmbl_cds_structure{$asmbl}})
	{
		my @overlappings = check_overlapping($exon, $rice_cds_structure->{$rice_cds});
		push @retroposed_exons, join("_", @$exon) if @overlappings >= 2;
	}
	print OUT join("\t", ($asmbl, @retroposed_exons)), "\n";
	
	print join(":", ($asmbl, output_for_record($asmbl_cds_structure{$asmbl})));
	print "\t",join(":", ($rice_gene_alias,output_for_record($rice_cds_structure->{$rice_cds}))), "\n";	
}

# Subroutines
sub output_for_record
{
	my $arrayref = shift;
	my @return = map{join("_", @$_)}@$arrayref;
	return @return;
}

sub check_overlapping
{
	my ($exon, $rice_exons) = @_;
	my @return;
	foreach (@$rice_exons)
	{
		my ($start, $end) = @$_;
		if ( ($start >= $exon->[0] and $start<= $exon->[1]) or ($exon->[0]>=$start and $exon->[1]<=$end))
		{
			my @tmp = sort{$a<=>$b} (@$exon, $start, $end);
			my $overlap = $tmp[2] - $tmp[1];
			push @return, $overlap if $overlap >= 10;
		}	
	}
	return @return;
}


sub read_asmbl_gff
{
	my $file = shift;
	my $pos_transform = shift;
	my $query_start_end = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		#contig117559	alignAssembly-pasa_chr3A_clean_masked	cDNA_match	280	713	.	+	.	ID=chain_151;Target=asmbl_151 1 434 +
		if(/Target=(\S+)\s+(\d+)\s+(\d+)/)
		{
			my ($id, $start, $end) = ($1, $2, $3);
			next unless exists $pos_transform->{$id};
			
			my ($query_start, $query_end) = @{$query_start_end->{$id}};
			if( ($start >= $query_start and $start <= $query_end) or ($query_start>=$start and $query_start <= $end) )
			{
				my @tmp = sort{$a<=>$b}($query_start, $query_end, $start, $end);
				my $overlap = $tmp[2] - $tmp[1] + 1;
				next unless $overlap >= 0.8 * abs($end - $start);
			}
			else
			{
				next;
			}
						
			my $off_pos = $pos_transform->{$id}[1];			
			if($off_pos =~ /N(\S+)/)
			{
				$start = abs($1) - $start;
				$start = 1 if $start <=  0;
				$end = abs($1) - $end;
				$end = 1 if $end <= 0;
			}
			elsif($off_pos =~ /S(\S+)/)
			{
				$start -= $1;
				$start = 1 if $start <= 0;
				$end -= $1;
				$end = 1 if $end <= 0;
			}
			#if($start<0 or $end<0){print STDERR join("\t", ($id, $start, $end, $off_pos)), "\n"; exit}
			push @{$return{$id}}, [sort{$a<=>$b}($start, $end)] unless $start == $end;
		}
	}
	close IN;
	return %return;
}

sub read_blast
{
	my $file = shift;
	my %return;
	my %query_start_end;
	my $searchio = Bio::SearchIO->new(-file => $file, -type => 'blast');
	while (my $result = $searchio->next_result)
	{
		my $name = $result->query_name;
		my $hit = $result->next_hit;
		next unless defined $hit;
		my $hit_name = $hit->name;
		#print STDERR $hit_name, "\n"; exit;
		my $hsp = $hit->next_hsp;
		next unless defined $hsp;
		my $query_strand = $hsp->strand('query');
		my $hit_strand = $hsp->strand('hit');
		my $off_pos;
		if($name eq 'asmbl_71'){print STDERR $hsp->start('query'), "\t",$hsp->start('hit'), "\n"}
		if($query_strand eq $hit_strand)
		{
			$off_pos = $hsp->start('query') - $hsp->start('hit');
			$off_pos  = "S" . $off_pos;
		}
		else
		{
			$off_pos = -1 * ($hsp->start('query')>$hsp->end('query')?$hsp->start('query'):$hsp->end('query') +  $hsp->start('hit') < $hsp->end('hit')?$hsp->start('hit'):$hsp->end('hit'));
			$off_pos = "N" . $off_pos;
		}
		$return{$name} = [$hit_name, $off_pos];
		$query_start_end{$name} = [sort{$a<=>$b} ($hsp->start('query'), $hsp->end('query'))];
	}
	return (\%return, \%query_start_end);
}

sub read_rice_gff
{
	my $file = shift;
	my %return;
	my %alias_id;
	open (IN, $file) or die;
	while (<IN>)
	{
		chomp; next if /^\#/;
		my @t = split /\s+/,$_;
		if($t[2] eq 'mRNA')
		{
			if(/ID=(\S+?)\;.*Alias=(\S{16})/)
			{
				$alias_id{$2} = $1;
			}
			next
		}
		
		if($t[2] eq 'CDS')
		{
			# Chr1	MSU_osa1r6	CDS	2449	2616	.	+	0	Parent=13101.m00002
			my $id = $1 if /Parent=(\S{12})/;
			push @{$return{$id}}, [@t[3, 4]];
		}
	}
	close IN;
	foreach (keys %return)
	{
		$return{$_} = [ transform_pos($return{$_}) ];
	}
	return(\%return, \%alias_id)
}

sub transform_pos
{
	my $arrayref = shift;
	my @return;
	
	foreach my $index (0 .. (scalar @$arrayref)-1)
	{
		my @arr = sort{$a <=> $b} @{$arrayref->[$index]};
		my @new_pos = map{$_ - $arr[0] + 1}@arr;
		if($index)
		{
			@new_pos = map{$_ + $return[$index-1][1]}@new_pos;
		}
		push @return, [@new_pos];		
	}
	
	return @return;
}




