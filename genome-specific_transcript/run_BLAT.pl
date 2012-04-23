#!/usr/bin/perl -w
use strict;
use Parallel::ForkManager;
use Getopt::Long;


my $sUsage = <<USAGE;

perl $0
	--seqType  format of query file, fa or fq
	--query    query file
	--db       database file, must be fasta format
	--cpu      to speed up the blat, use multiple CPUs; default is 1;
	--blat     path to executable blat, default blat
  
USAGE

my ($seqtype, $query_file, $db_file);
my $cpu = 1; #default
my $blat_bin = "blat";

GetOptions(
"seqType=s" => \$seqtype,
"query=s"   => \$query_file,
"db=s"      => \$db_file,
"cpu=i"     => \$cpu,
"blat=s"    => \$blat_bin,
);
die $sUsage unless defined $seqtype and defined $db_file and defined $query_file;

my $output_type = "-out=blast8";
my $output_name = $query_file . ".blat.out";

my @blat_cmds;
if ($cpu <= 1)
{
	my $nf = $query_file;
	if($seqtype eq 'fq')
	{
		$nf = $query_file . ".fasta";
		fastq_to_fasta($query_file, $nf);
	}
	my $cmd = join(" ", ($blat_bin, $db_file, $nf, $output_type, $output_name));
	print_time_stamp("Running command: $cmd");
}
else
{
	my @sub_files = split_file($query_file, $cpu, $seqtype);
	my @sub_out_files;
	foreach my $sub (@sub_files)
	{
		my $sub_output = $sub . '.blat.out';
		push @sub_out_files, $sub_output;
		my $cmd =  join(" ", ($blat_bin, $db_file, $sub, $output_type, $sub_output));
		push @blat_cmds, $cmd;
	}
	my $pm = new Parallel::ForkManager(scalar @blat_cmds);
	foreach my $cmd (@blat_cmds)
	{
		my $pid = $pm->start and next;
		print_time_stamp("Running command: $cmd");
		die if system($cmd);
		$pm->finish;
	}
  $pm->wait_all_children;
	
	my $rec = system("cat @sub_out_files > $output_name");
	die "Something is wrong while running: \n cat @sub_out_files > $output_name" if $rec;
	system("rm -rf @sub_files @sub_out_files");
}

# Subroutines
sub split_file
{
	my ($file, $num, $type) = @_;
	my $nf = $file; 
	if($type eq "fq")
	{
		$nf = $file . ".fasta";
		fastq_to_fasta($file, $nf);
	}
	
	my @sub_files = map{join('_',($nf, $_))} 1..$num;
	my $out_dir = './tmp';
	mkdir($out_dir) unless -d $out_dir;
	map{$sub_files[$_]= $out_dir .'/'. $sub_files[$_] }0..$#sub_files;
	
	split_fasta_file($nf, $num, @sub_files);

	return @sub_files;	
}

sub split_fasta_file
{
	my ($orig_fasta, $num_files, @sub_files) = @_;
	my $total_seqs =  get_num_total_seq($orig_fasta);
	my $seq_each = int($total_seqs/$num_files);
	$seq_each += 1 if ($total_seqs%$num_files)>0;
	print '$seq_each: ', $seq_each, "\n";
	my %sub_file_handles = get_file_handles(@sub_files);
	open (IN, "$orig_fasta") or die "can't open file $orig_fasta\n";
	my $count = 0;
	while(<IN>)
	{
		next if /^\s+$/;
		$count++ if /^>/;
		my $file_index = int($count/$seq_each) + (($count%$seq_each == 0)?0:1) - 1;
		$file_index  = $num_files - 1 if $file_index >= $num_files;
		#print $count,"\t", $file_index,"\n";
		print {$sub_file_handles{$file_index}} $_;
	}
	foreach (values %sub_file_handles)
	{
		close $_ or print STDERR  "Can't close File handle $_\n";
	}

}

sub fastq_to_fasta
{
	my($fastq_file, $fasta_file) = @_;
	open (IN, "$fastq_file") or die;
	open (OUT, ">$fasta_file") or die;
	my $count=0;
	my @temp;
	while(<IN>)
	{
		next if /^\s+$/;
		$count++;
		if($count==1){s/^\@/\>/; print OUT $_}
		elsif($count==2){print OUT $_}
		elsif($count==4){$count=0;next}
	}
	close IN;
	close OUT;
}

sub get_num_total_seq
{
	my $file = shift;
	my $total;
	open (IN, "$file") or die $!;
	while (<IN>)
	{
		$total++ if /^>/;
	}
	close IN;
	return $total;
}

sub get_file_handles
{
	my @files = @_;
	my %return_hash;
	foreach (0..$#files)
	{
		local *FH;
		open (FH, ">$files[$_]") or die "can't open file $files[$_]\n";
		$return_hash{$_} = *FH{IO};
	}
	return %return_hash;
}

sub print_time_stamp
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $c,"\n";
}