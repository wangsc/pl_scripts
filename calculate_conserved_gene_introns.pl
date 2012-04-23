#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $sUsage = qq(
perl $0
<ortholog_pair>
<tblastx result>
<pasa transcript isoform gff>
<wheat_intron_AS>
<rice_intron_AS>
);
die $sUsage unless @ARGV >= 5;
my ($ortho_pair_file, $tblastx_file, $pasa_gff, $wheat_AS_file, $rice_AS_file) = @ARGV;

my %ortho_pairs = read_ortho_pair_file($ortho_pair_file);
my %conserved_introns_pairs = read_tblastx_result($tblastx_file, \%ortho_pairs);
while(my ($k, $v) = each %conserved_introns_pairs){print join("\t", ($k, $v)),"\n"} exit 1;
my %wheat_intron_pos = transform_wheat_intron_pos($pasa_gff, keys %conserved_introns_pairs);
print 'Total number of conserved introns: ', scalar unique(values %wheat_intron_pos),"\n";
my %rice_intron_pos = transform_rice_intron_pos(\%conserved_introns_pairs);

my %wheat_AS = read_AS_file($wheat_AS_file);
my %rice_AS = read_AS_file($rice_AS_file);

# Conserved gene
my %cons_gene;
my $total_cons_genes = 0;
my %wheat_genes;
foreach (values %wheat_intron_pos)
{
	my @t = split /_/, $_;
	$wheat_genes{$t[0]}=1;
}
my $total_wheat_genes_with_cons_intron = scalar keys %wheat_genes;
print '$total_wheat_genes_with_cons_intron: ', $total_wheat_genes_with_cons_intron,"\n";
my %wheat_genes_with_cons_intron_AS;
my $num_rice_gene_with_as = 0;
my $num_wheat_gene_with_as = 0;
my %w_type_count;
foreach my $gene (keys %ortho_pairs)
{
	my $flag = 0;
	my $wheat_gene = $ortho_pairs{$gene};
	my %rice_types = types($rice_AS{$gene});
	my %wheat_types = types($wheat_AS{$wheat_gene});
	$num_rice_gene_with_as ++ if (scalar keys %rice_types) > 0;
	$num_wheat_gene_with_as++ if (scalar keys %wheat_types) > 0;
	$wheat_genes_with_cons_intron_AS{$gene}=1 if (scalar keys %wheat_types) > 0 and exists $wheat_genes{$wheat_gene};
	foreach (keys %wheat_types)
	{
		$w_type_count{$_}++;
		if (exists $rice_types{$_})
		{
			$cons_gene{$gene}{$_}=1 ;
			$flag=1;
		}		
	}
	#print STDERR join("\t", ($gene, $wheat_gene)),"\n" if $flag;
	$total_cons_genes++ if $flag;
}
print '$total_wheat_genes_with_cons_intron_AS: ', scalar keys %wheat_genes_with_cons_intron_AS,"\n";
print '$num_rice_gene_with_cons_intron_AS: ', $num_rice_gene_with_as,"\n";
# conserved intron
my %cons_intron_both;
my %cons_intron_wheat;
my $total_cons_introns = 0;
my $total_cons_events = 0;
my %conserved_introns;
foreach (keys %conserved_introns_pairs)
{
	my ($rice_gene, $r_start, $r_end) = split /_/, $rice_intron_pos{$_};
	my ($wheat_gene, $w_start, $w_end) = split /_/, $wheat_intron_pos{$_};
	#print STDERR $rice_gene,"\t", $wheat_gene,"\n";exit;
	next unless exists $cons_gene{$rice_gene};
	#$total_cons_introns++;
	my @rice_types = @{$rice_AS{$rice_gene}{join('_',($r_start, $r_end))}} if exists $rice_AS{$rice_gene}{join('_',($r_start, $r_end))};
	my @wheat_types = @{$wheat_AS{$wheat_gene}{join('_',($w_start, $w_end))}} if exists $wheat_AS{$wheat_gene}{join('_',($w_start, $w_end))};
	$conserved_introns{$wheat_intron_pos{$_}} = 1 if @wheat_types>0;
	$total_cons_introns++ if @wheat_types>0;
	next if @wheat_types==0;
	my %w_type_hash = map{$_, 1}@wheat_types;
	my $flag = 0;
	foreach (@rice_types)
	{
		$cons_intron_both{$_}++ if exists $w_type_hash{$_};
		$flag = 1 if exists $w_type_hash{$_};
	}
	$total_cons_events++ if $flag;
	map{$cons_intron_wheat{$_}++} @wheat_types;
}
#$total_cons_introns = scalar keys %conserved_introns;
# Output
my @types = qw(AltD AltA AltP ExonS IntronR);
foreach my $type (@types)
{
	my $num_gene = 0;
	foreach(keys %cons_gene){$num_gene++ if exists $cons_gene{$_}{$type}}
	print $type, "\t", 
	$num_gene,"\t",
	$w_type_count{$type},"\t",
	exists $cons_intron_wheat{$type}?$cons_intron_wheat{$type}:0,"\t",
	exists $cons_intron_both{$type}?$cons_intron_both{$type}:0,"\n";
}
print 'Total',"\t", $total_cons_genes, "\t", $num_wheat_gene_with_as,"\t", $total_cons_introns, "\t", $total_cons_events,"\n";

# Subroutines
sub unique
{
	my %h = map{$_, 1} @_;
	return keys %h;
}

sub types
{
	my $hashref = shift;
	my %return;
	foreach( values %{$hashref})
	{
		map{$return{$_}=1} @$_;
	}
	return %return;
}

sub transform_rice_intron_pos
{
	my $hashref = shift;
	my %hash = %{$hashref};
	my %return;
	foreach my $id (keys %hash)
	{
		my $flag = 0;
		$flag = 1 if $hash{$id}=~/LOC/;
		my @t = split /_/,$hash{$id};
		my @pos = $flag?@t[3..6]:@t[2..5];
		@pos = sort{$a<=>$b} @pos;
		#print STDERR join('_', ($t[1], $pos[1]+1, $pos[2]-1)),"\n";
		$return{$id} = $flag?join('_', ($t[1], $pos[1]+1, $pos[2]-1)):join('_', ($t[0], $pos[1]+1, $pos[2]-1));
	}
	return %return;
}


sub transform_wheat_intron_pos
{
	my $file = shift;
	my @ids = @_;
	my %asmbl_pos = read_isoform_gff($file);
	my %return;
	foreach my $id (@ids)
	{
		#S4824_asmbl_5699_1_42_43_175
		my @t = split /_/,$id;
		die $id unless defined $asmbl_pos{join('_', (@t[1..4]))};
		die $id unless defined $asmbl_pos{join('_', (@t[1,2,5,6]))};
		my @pos = (@{$asmbl_pos{join('_', (@t[1..4]))}}, @{$asmbl_pos{join('_', (@t[1,2,5,6]))}});
		@pos = sort{$a<=>$b}@pos;
		$return{$id} = join('_', ($t[0], $pos[1]+1, $pos[2]-1));
	}
	return %return;
}

sub read_isoform_gff
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while (<IN>)
	{
		chomp;
		next unless /^\S/;
		my ($asmbl_id, $start, $end) = ($1, $2, $3) if /Target=(\S+)\s+(\d+)\s+(\d+)/;
		die $_,"*\n" unless defined $asmbl_id;
		my @t= split /\t/, $_;
		$return{join('_', ($asmbl_id, $start, $end))} = [@t[3,4]];
	}
	close IN;
	return %return;
}



sub read_AS_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		chomp;
		my @t = split /\t/, $_;
		next unless @t > 1;
		my $id = shift @t;
		$id =~ s/LOC_//;
		foreach (@t)
		{
			my @pos_types = split /_/,$_;
			my @types = @pos_types[2..$#pos_types];
			my @pos = sort{$a<=>$b}@pos_types[0,1];
			push @{$return{$id}{join('_',@pos)}}, @types;
		}
	}
	close IN;
	return %return;
}


sub read_tblastx_result
{
	my $file = shift;
	my $ortho_pair = shift;
	my %return_hash;
	my $searchio = Bio::SearchIO->new(-format => 'blast', file => "$file" );
	my $total_matched = 0;
	while (my $result = $searchio->next_result())
	{
		last unless defined $result;
		my $query_name = $result->query_name;
		my $query_id = $1 if $query_name=~ /^(\S+?)\_/;
		my $query_length = $result->query_length;
		my $start_range = [1,30];
		my $end_range = [$query_length-29, $query_length];
	#	print $query_name, "\t";
		my $hit = $result->next_hit();
		next unless defined $hit;
		my @sp_names = split /_/, $hit->name;
		my $gene_id = $sp_names[0];
		$gene_id = $sp_names[1] if $hit->name =~ /LOC/;
		next unless $ortho_pair->{$gene_id} eq $query_id;
		my ($flag_a, $flag_b) = (0, 0);
		while(my $hsp = $hit->next_hsp)
		{
			my @range = $hsp->range('query');
			@range = sort{$a<=>$b} @range;
			$flag_a=1 if overlapping([@range], $start_range)>=10;
			$flag_b=1 if overlapping([@range], $end_range)>=10;
		}
		
		$return_hash{$query_name} = $hit->name if ($flag_a and $flag_b);
	}
	return (%return_hash);
}

sub overlapping
{
	my ($arraya, $arrayb) = @_;
	if(($arraya->[0] >= $arrayb->[0] and $arraya->[0] <= $arrayb->[1]) or ($arrayb->[0] >= $arraya->[0] and $arrayb->[0] <= $arraya->[1]))
	{
		my @t = sort{$a<=>$b}(@$arraya, @$arrayb);
		return $t[2]-$t[1]+1;
	}
	else {return 0;}
}

sub read_ortho_pair_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		chomp;
		my @t = split /\t/, $_;
		my $id;
		$id = $1 if /LOC_(.*?)\.\d/;
		$id = $1 if /(Brad.*?)\.\d/;
		$return{$id} = $t[-1];
	}
	close IN;
	return %return;
}