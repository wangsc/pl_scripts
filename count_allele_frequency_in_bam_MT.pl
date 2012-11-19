#!/usr/bin/perl -w
use strict;
use Bio::DB::Sam;
use File::Basename;
use threads;
use threads::shared;
use Thread::Semaphore;
use Getopt::Long;

my $sUsage = qq(
perl $0
--bam file1,file2,file3,...
--snp SNP file, combined output of "vcfutils.pl varFilter" from each accession
--ref reference fasta file
--out output file
);
die $sUsage unless @ARGV >= 4;

my($bamfile, $snpfile, $ref_fasta_file, $outfile);
GetOptions
(
	'bam=s' => \$bamfile,
	'snp=s' => \$snpfile,
	'ref=s' => \$ref_fasta_file,
	'out=s' => \$outfile
);
die $sUsage unless defined $bamfile and defined $snpfile and defined $ref_fasta_file and defined $outfile;

print_time_comment("Reading SNP file $snpfile ...");
my %snps = read_snp_file($snpfile);
my @bam_files = split /,/, $bamfile;

# shared variables
my $allele_count :shared;
my $ref_allele :shared;
my $covered_snps :shared;

# semaphore
my $semaphore = Thread::Semaphore->new();
my $num_thread = 0;
my @thr_names;
foreach my $file (@bam_files)
{
	$num_thread++;
	push @thr_names, "thr".$num_thread;
	$main::{"thr".$num_thread} = threads->create(\&count_allele_freq_in_bam, $file, $ref_fasta_file, $num_thread);
	print STDERR  "Startgin Thread ", $num_thread, "\n";
	#$main::{"thr".$num_thread}->detach();
}

foreach my $thr (threads->list())
{
	$thr->join();
}

# Output
open (OUT, ">$outfile") or die "$!\n";
print_time_comment("Output results to file $outfile \n");
foreach my $id (keys %{$allele_count})
{
	my @covered_snps = @{$covered_snps->{$id}};
	foreach my $snppos (unique(@covered_snps))
	{
		print OUT $id, "\t", $snppos, "\t";
		print OUT $ref_allele->{$id}{$snppos},"_", 
		     exists $allele_count->{$id}{$snppos}{$ref_allele->{$id}{$snppos}}?$allele_count->{$id}{$snppos}{$ref_allele->{$id}{$snppos}}:0;
		foreach (keys %{$allele_count->{$id}{$snppos}})
		{
			unless($_ eq  $ref_allele->{$id}{$snppos})
			{
				print OUT "\t", $_,"_",exists $allele_count->{$id}{$snppos}{$_}?$allele_count->{$id}{$snppos}{$_}:0;
			}
		}
		print OUT "\n";
	}	
}

close OUT;

# subroutines

sub make_shared_hash
{
	my $snps_ref = shift;
	my $allele_count  = &shared({});
	my $ref_allele = &shared({});
	my $covered_snps = &shared({});	
	
	foreach my $id(keys %$snps_ref)
	{
		my %p = &shared({});
		$allele_count->{$id} = \%p;
		my @c = &shared([]);
		$covered_snps->{$id} = \@c;
		
	}
}

sub count_allele_freq_in_bam
{
	my ($file, $ref_fasta, $num_thread) = @_;
	my $bam_index = $file . ".bai";
	my $autoindex = -e $bam_index?1:0;
	my $sam = Bio::DB::Sam->new(-bam  => $file,
	                            -fasta=> $ref_fasta,
	                            -autoindex => $autoindex,
	                            );
	my $alignments = $sam->features(-iterator=>1);	
	my $counter = 0;
	my %covered_snps_tmp;
	my %ref_allele_tmp;
	my %allele_count_tmp;
	while(my $aln = $alignments->next_seq())
	{
		$counter++;
		my $id = $aln->seq_id();
		next unless exists $snps{$id};
		my $ref_start = $aln->start;
		next unless defined $ref_start;
		my $ref_end = $aln->end;
		my @snps = grep {$_>=$ref_start and $_<=$ref_end} @{$snps{$id}};
		push @{$covered_snps_tmp{$id}}, @snps;
		$covered_snps_tmp{$id} = [unique(@{$covered_snps_tmp{$id}})];
		my $ref_dna   = $aln->dna;
		my $query_dna = $aln->query->dna;
		my @qscores = $aln->qscore;
		foreach my $snppos (@snps)
		{
			next unless $qscores[$snppos-$ref_start] >= 20;
			$ref_allele_tmp{$id}{$snppos} = substr($ref_dna, $snppos-$ref_start, 1) unless exists $ref_allele_tmp{$id}{$snppos};
			$allele_count_tmp{$id}{$snppos}{substr($query_dna, $snppos-$ref_start, 1)}++;
		}
		unless ($counter%1000)
		{
			#print STDERR "\tThread $num_thread finished $counter alignments ...\n";
			flush(\%covered_snps_tmp, %ref_allele_tmp, %allele_count_tmp);
		}		
	}
	flush(\%covered_snps_tmp, %ref_allele_tmp, %allele_count_tmp);
}

sub flush
{
	$semaphore->down();
	print STDERR "Thread which is flusing data: ", threads->tid(), "\n";
	my ($snps_ref, $ref_allele_ref, $allele_count_ref) = @_;	
	foreach my $id (keys %$snps_ref)
	{
		my @snps = @{$snps_ref->{$id}};
		my %p :shared;
		$covered_snps = \%p unless defined $covered_snps;
		push @{$covered_snps->{$id}}, @snps;
		$covered_snps->{$id} = [unique(@{$covered_snps->{$id}})];	
		foreach my $snppos (@snps)
		{
			$ref_allele->{$id}{$snppos} = $ref_allele_ref->{$id}{$snppos} unless exists $ref_allele->{$id}{$snppos};
			foreach my $allele (keys %{$allele_count_ref->{$id}{$snppos}})
			{
				$allele_count->{$id}{$snppos}{$allele} = $allele_count_ref->{$id}{$snppos}{$allele}
			}
		}
	}
	%$snps_ref = ();
	%$ref_allele_ref = ();
	%$allele_count_ref = ();
	$semaphore->up();
}


sub min_max
{
	my $min = shift;
	my $max = $min;
	foreach (@_)
	{
		$min = $_ if $_ < $min;
		$max = $_ if $_ > $max;
	}
	return ($min, $max);
}

sub unique
{
	my %h = map{$_, 1}@_;
	return keys %h;	
}

sub read_snp_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die;
	while(<IN>)
	{
		next if /^\#/;
		next if /^\s+$/;
		chomp;
		my @t = split /\s+/,$_;
		push @{$return{$t[0]}}, $t[1];
	}
	close IN;
	map{ $return{$_} = [unique(@{$return{$_}})] }keys %return;
	return %return;
}

sub print_time_comment
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}
