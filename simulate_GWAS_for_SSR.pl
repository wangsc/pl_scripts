#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<genotype input file>
<file contain SSR markers>
<tassel directory>
);
die $sUsage unless @ARGV >= 3;
my($genotype_input_file, $ssr_infor_file, $tassel_direc) = @ARGV;

my %ssr_info = read_ssr_info_file($ssr_infor_file);

foreach my $ssr_id (keys %ssr_info)
{
	my $phenotype_file = $ssr_id . "_phenotype";
	my @header;
	my $num_acc = scalar @{$ssr_info{$ssr_id}};
	open (PE, ">$phenotype_file") or die;
	print PE $num_acc,"\t1\t1\n\t$ssr_id\n";
	foreach (@{$ssr_info{$ssr_id}})
	{
		print PE join("\t", @$_),"\n";
	}
	close PE;
	run_gwas($ssr_id, $genotype_input_file, $phenotype_file, $tassel_direc);
	
}

sub read_ssr_info_file
{
	my $file = shift;
	open (IN ,$file) or die;
	my $f = 0;
	my %ids_index;
	my %record;
	while(<IN>)
	{
		chomp;
		$f++;
		if($f==1)
		{
			chomp; 
			my @t = split/\t/,$_;
			foreach (1..$#t)
			{
				$ids_index{$_} = $t[$_];
			}
			next;
		}
		chomp;
		my @t = split/\t/,$_;
		foreach (1..$#t)
		{
			my $pheno = $t[$_];
			$pheno =~ s/:\d+//;
			$pheno = -999 if $pheno =~ /\?/;
			push @{$record{$ids_index{$_}}}, [$t[0], $pheno]
		}
	}
	close IN;
	return %record;
}


sub run_gwas
{
	my ($ssr_id, $genotype_input, $phenotype_input, $tassel_dir) = @_;

	my $libdir = $tassel_dir . '/'. "lib";
	my @jars = <$libdir/*.jar>;
	push @jars , $tassel_dir . '/'. 'sTASSEL.jar';
	my $CP = join(":", @jars);

	my $parameters = '-p "'. $genotype_input .'" -t "'. $phenotype_input. '" -k "K_matrix_237IND.txt" -q "Q_matrix_237IND.txt" -glm -o tassel_result_'. $ssr_id . ' -txt';
	my $cmd = "java -classpath '$CP' -Xmx2000m net.maizegenetics.pipeline.MLMGLMFileInputPipeline $parameters";

	#print $cmd,"\n";
	die $! if system($cmd); 
}