#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $sUsage = qq(
perl $0
<wheat pasa transcript isoform gff3>
<rice gff3>
<ortholog pairs for wheat and rice>
);
die $sUsage unless @ARGV >= 3;

my ($pasa_gff, $rice_gff, $ortho_pair_file, $blastx_file) = @ARGV;
#my ($asmbl_gene, %wheat_rice) = read_ortho_pair_file($ortho_pair_file);
#my %rice_wheat = reverse %wheat_rice;

my %wheat_as_events = compare_AS_within_wheat($pasa_gff);
my %rice_as_enents = compare_AS_within_rice($rice_gff);
output('wheat_introns_AS.out', %wheat_as_events);
output('rice_introns_AS.out', %rice_as_enents);


# Subroutines
sub output
{
	my ($file, %result) = @_;
	open (OUT, ">$file") or die "can't open file $file\n";
	foreach my $gene(keys %result)
	{
		print OUT $gene;
		foreach (keys %{$result{$gene}})
		{
			my @arr = unique(@{$result{$gene}{$_}});			
			print OUT "\t", join("_", ($_, @arr));
		}
		print OUT "\n";
	}
	close OUT;
}
sub unique
{
	my %h = map{$_, 1} @_;
	return keys %h;
}


sub read_ortho_pair_file
{
	my $file = shift;
	open (IN, $file) or die $!;
	my %return;
	my %asmbl_gene;
	while(<IN>)
	{
		chomp;
		next if /^\s+$/;
		my @t = split /\s+/, $_;
		my $rice_gene = $1 if $t[1] =~ /^(.*?)\.\d/;
		$return{$t[2]} = $rice_gene;
		$asmbl_gene{$t[0]} = $t[2];
	}
	close IN;
	return (\%asmbl_gene, %return);
}

sub compare_AS_within_wheat
{
	my ($gff) = @_;
	my (%ctg_pos, %gene_pos, %strand);
	open (IN, $gff) or die "$!\n";
	while (<IN>)
	{
		next if /^\s+$/;
		my $gene_id = $1 if /ID=(.*?)\-/;
	#	next unless exists $wheat_rice{$gene_id};
		my @t = split /\t/, $_;
		my @ctg_pos = @t[3,4];
		$strand{$gene_id} = $t[6] unless exists $strand{$gene_id};
		my @asmbl_pos = ($1, $2) if /Target=\S+\s+(\d+)\s+(\d+)/;
		my $asmbl_id = $1 if /Target=(\S+)/;
		push @{$ctg_pos{$gene_id}{$asmbl_id}}, [@ctg_pos];
		push @{$gene_pos{$gene_id}{$asmbl_id}}, [@asmbl_pos];
	}
	close IN;
	my %gene_as_type;
	foreach my $gene (keys %gene_pos)
	{
		$gene_as_type{$gene} = compare_iso(values %{$ctg_pos{$gene}}, $strand{$gene});
	}
	return (%gene_as_type);
}

sub compare_AS_within_rice
{
	my ($gff) = @_;
	my (%gene_pos, %strand);
	open (IN, $gff) or die "$!\n";
	my $gene_id; my $iso_id;
	while (<IN>)
	{
		next unless /^Chr\d/;
		chomp;
		my @t = split /\t/,$_;
		if ($t[2] =~ /mRNA/)
		{
			($gene_id, $iso_id ) = ($1, $2) if /Alias=(.*?)\.(\d)/;
			next;
		}
		next unless $t[2] =~ /cds/i;
		#next unless exists $rice_wheat{$gene_id};		
		$strand{$gene_id} = $t[6] unless exists $strand{$gene_id};
		push @{$gene_pos{$gene_id}{$iso_id}}, [@t[3,4]];
	}
	close IN;
	my %gene_as_type;
	foreach my $gene (keys %gene_pos)
	{
		$gene_as_type{$gene} = compare_iso(values %{$gene_pos{$gene}}, $strand{$gene});
	}
	return (%gene_as_type);
}

sub compare_iso
{
	my $strand = pop;
	my @array = @_;
	my @exon_vecs = map{construct_vec($_)} @array;
	my @introns;
	foreach (@array)
	{
		my $arrref = $_;
		my @iso_introns;
		next unless @$_ > 1;
		foreach my $index ( 0..(scalar @$arrref -2))
		{
			if ($strand eq '-')
			{
				push @iso_introns, [$arrref->[$index+1]->[1]+1, $arrref->[$index]->[0]-1]
			}
			else
			{
				push @iso_introns, [$arrref->[$index]->[1]+1,  $arrref->[$index+1]->[0]-1]
			}
		}
		push @introns, [@iso_introns];
	}
	my %as_types = compare_introns(@introns, \@exon_vecs);
	return \%as_types;
}

sub compare_introns
{
	my $exon_vecs = pop;
	my @introns = @_;
	my @vects = map{construct_vec($_)}@introns;
	my %return;
	foreach my $index (0..$#introns)
	{
		my $vec_out = $vects[$index];
		foreach my $ind (0..$#introns)
		{
			my $vec_in = $vects[$ind];
			my %comp_two = compare_two_iso_introns($vec_out, $vec_in, $introns[$index], $introns[$ind], $exon_vecs->[$ind]);
			foreach my $k (keys %comp_two)
			{
				if(exists $return{$k})
				{
					push @{$return{$k}}, @{$comp_two{$k}}
				}
				else
				{
					$return{$k} = $comp_two{$k};
				}
			}
		}
	}
	return %return;
}

sub compare_two_iso_introns
{
	my ($va, $vb, $ina, $inb, $exon_vec_b) = @_;
	my %types;
	foreach (@$ina)
	{ 
		my ($s_a, $e_a) = @$_;
		my $id = join('_', @$_);
		my @overlap;
		foreach (@$inb)
		{
			my ($s_b, $e_b) = @$_;
			next unless ($s_a >=$s_b and $s_a <=$e_b) or ($s_b >=$s_a and $s_b <=$e_a);
			push @overlap, $_;
			push @{$types{$id}}, 'AltD' if ($s_a != $s_b) and ($e_a == $e_b);
			push @{$types{$id}}, 'AltA' if ($s_a == $s_b) and ($e_a != $e_b);
			push @{$types{$id}}, 'AltP' if (($s_a != $s_b) and ($e_a != $e_b));
		}
		push @{$types{$id}}, 'ExonS' if @overlap >=2;
		my $covered_by_exon = 0;
		foreach ($s_a, $e_a)
		{
			$covered_by_exon ++ if vec($exon_vec_b, $_, 1) == 1;
		}
		push @{$types{$id}}, 'IntronR' if $covered_by_exon >= (abs($s_a-$e_a)+1)*90;
	}
	return %types;
}


sub construct_vec
{
	my $arrayref = shift;
	my $vec = '';
	my $max;
	my $total;
	my $debug =1 ;
	foreach (@$arrayref)
	{
		my @d = sort{$a<=>$b}@$_;
#		print '@d: ', join("\t", @d),"\n" if $debug; $debug=0;
		foreach ($d[0]..$d[1])
		{
			vec($vec,$_,1) = 0b1;
			$total++;
			$max = $_ unless defined $max;
			$max = $_ if $_ > $max;
			
		}
	}
	return ($vec, $max, $total);
}



