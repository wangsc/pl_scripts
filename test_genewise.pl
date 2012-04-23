#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<rice pep file>
<rice gene seq file>
);
die $sUsage unless @ARGV >= 2;
my ($pep_file, $gene_seq_file) = @ARGV;
my %pep = read_fasta($pep_file);
my %gene_seq = read_fasta($gene_seq_file);
my $n = 5;
foreach my $id (keys %pep)
{
	open (P, ">pep") or die $!;
	print P '>',$id,"\n",$pep{$id};
	close P;
	open (G, ">genomic") or die $!;
	my $gene = $1 if $id=~/(.*?)\.\d/;
	print $id,"\n";
	print G ">",$gene,"\n",$gene_seq{$gene},"\n";
	close G;
	run_genewise();
	system("cat tmp.gff >>test_genewise_in_rice.gff");
	$n--;
	last unless $n;
}



sub read_fasta
{
	my $file = shift;
	open (IN, "$file") or die "$! $file\n";
	my %return_hash;
	my $id;
	my $debug = 1;
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		if (/^>(\S+?)\|/)
		{
			$id = $1;
			#print 'fasta id: ', $id ,"\n" if $debug; $debug = 0;
			$return_hash{$id} = '';
			next;
		}
		print STDERR $_,"\n" if not defined $id;
		$return_hash{$id} .= $_;
	}
	close IN;
	return %return_hash;
}

sub run_genewise
{
	my $command = "genewise -gff pep genomic -splice flat -alg 2193 -null flat -both >tmp.gff";
	#print $command,"\n";
	die "genewise failed \n" if system($command);
}