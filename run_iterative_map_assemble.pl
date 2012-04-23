#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use threads;
use threads::shared;
use Getopt::Long qw(GetOptions);

my $sUsage  = qq(
************************************************************************
# This script will run mapping (BLAT) and assemble (Cap3) iteratively 
# until no new reads are added or reach the max_run.
# SW Jun 1 2011
Usage:
	perl $0
	-blat executive path of blat
	-cap3 executive path of cap3
	-read fasta file containing reads
	-db   fasta file which were used as database in blat
	-min_len  minimum length of mapped fragment in blat, optional
	-threads  the number of multiple process of blat, must be less than 10, default is 1;
	-max_run  maximum number of interations, optional
************************************************************************
);
die $sUsage unless @ARGV >= 4;
my ($blat, $cap3, $read_file, $db_file, $max_run, $min_len, $thread);
GetOptions('blat=s' => \$blat, 
					 'cap3=s' => \$cap3, 
					 'read=s' => \$read_file, 
					 'db=s' => \$db_file, 
					 'max=i' => \$max_run, 
					 'min=i' => \$min_len,
#				   'threads=i' => \$thread
					 );
$min_len = 20 unless defined $min_len;
$thread = 1 unless defined $thread;
print_time_comment("Generating Bio::DB::Fasta object ...");
my $reads_obj = Bio::DB::Fasta->new($read_file);
my $total_reads = scalar $reads_obj->ids;
my %mapped_reads_recorder;
my $count_run = 0;
my $flag = 1;
my @sub_read_files =  map{join('_', ($read_file, $_))}1..$thread;
print_time_comment("Spliting read file ...");
&split_file($total_reads, $read_file, $thread); #unless -e $sub_read_files[-1];

while ($flag) # keep running
{
	$count_run++;
	print_time_comment("Start runing $count_run round...");
	my $blat_out = $read_file . '_run_' . $count_run . '_blat.out';
	run_blat_threads($blat, $db_file, $read_file, $blat_out, \@sub_read_files, $thread);
	my @mapped_reads = parse_blat_result($blat_out, $min_len);
	my @new_reads;
	foreach (@mapped_reads)
	{
		if ( not exists $mapped_reads_recorder{$_})
		{
			$mapped_reads_recorder{$_} = 0;
			push @new_reads, $_;
		}
	}
	
	my $new_db_file;
	if (scalar @new_reads > 0)
	{
		my $mapped_read_file = extract_reads_seq($reads_obj, \@new_reads, $count_run);
		print_time_comment("Runing cap3 ...");
		$new_db_file = run_cap3($cap3, $mapped_read_file);		
	}
	else
	{
		print_time_comment("No new reads added. Run ". ($count_run-1) . ' is the last run!' );
		print_time_comment("File reads_for_assembly_run" . ($count_run-1) . '.cap.contigs is the last contig file!');
		$flag = 0;
		next;
	}

	if (defined $max_run and $count_run == $max_run)
	{
		$flag = 0;
	}
	else
	{
		$db_file = $new_db_file;
	}	
}

# subroutines

sub run_blat_threads
{
	my ($blat_bin, $db, $query, $out_file, $sub_files_ref, $thread) = @_;
	my $parameters = '-tileSize=12 -minScore=50 -minScore=50 -fastMap -out=blast8';
	if ($thread == 1)
	{
		print_time_comment("Start single process of blat ...");
		my $command = join(' ', ($blat_bin, $db, $query, $parameters, $out_file));
		print_time_comment("Runing command: $command ...");
		system($command) == 0 or die "Command: $command failed !\n";		
	}
	else
	{
		my @query_files = @$sub_files_ref;
		print_time_comment("\tStart multiple processes of blat (process \= $thread) ...");
		my @multi_blat;
		foreach my $count (1..$thread)
		{
			my $query_file = $query_files[$count-1];
			my $t_out_file = $query_file . '.blat';
			my $command = join(' ', ($blat_bin, $db, $query_file, $parameters, $t_out_file));
			my $thr = threads->new(\&run_blat, $command, $count);
			push @multi_blat, $thr;
		}
		foreach (@multi_blat)
		{
			my $count = $_->join();
			print_time_comment("\tprocess $count of blat finished ...");
		}
		print_time_comment("\tCombine the blat results ...");
		combine_blat_results_and_clean(\@query_files, $out_file);
	}
	return 1;
}

sub split_file
{
	my ($total_seqs, $orig_fasta, $num_files) = @_;
	my $seq_each = int($total_seqs/$num_files);
	$seq_each += 1 if $seq_each ==0;
#	print '$seq_each: ', $seq_each, "\n";
	my @sub_files = map{join('_',($orig_fasta, $_))} 1..$num_files;
	my %sub_file_handles = get_file_handles(@sub_files);
	open (IN, "$orig_fasta") or die "can't open file $orig_fasta\n";
	my $count = 0;
	while(<IN>)
	{
		next if /^\s+$/;
		chomp;
		$count++ if /^>/;
		my $file_index = int($count/$seq_each) + (($count%$seq_each == 0)?0:1) - 1;
		$file_index  = $num_files - 1 if $file_index >= $num_files;
	#	print $count,"\t", $file_index,"\n";
		print {$sub_file_handles{$file_index}} $_,"\n";
	}
	foreach (values %sub_file_handles)
	{
		close $_ or print STDERR  "Can't close File handle $_\n";
	}
	
	return @sub_files;
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

sub combine_blat_results_and_clean
{
	my ($sub_file, $out) = @_;
	foreach my $file (@$sub_file)
	{
		unlink($file) unless system("cat $file\.blat >>$out");
	}
}


sub run_blat
{
	my $cmd = shift;
	my $count_proc = shift;
	system($cmd) == 0 or die "Command: $cmd failed !\n";
	return $count_proc;
}


sub parse_blat_result
{
	my ($infile, $min_len) = @_;
	my @return_array;
	open (IN, "$infile") or die "can't open file $infile\n";
	while (<IN>)
	{
		next if /^\s+$/;
		my @line_data = split /\t/, $_;
		next if $line_data[3] < $min_len;
		push @return_array, $line_data[0]
	}
	close IN;
	return unique(@return_array);
}

sub extract_reads_seq
{
	my ($read_obj, $reads_ref, $count_run) = @_;
	my $read_file = 'reads_for_assembly_run'. $count_run;
	open (OUT, ">>$read_file") or die "can't open file $read_file\n";
	foreach (@{$reads_ref})
	{
		print OUT '>',$_,"\n",
					    $read_obj->get_Seq_by_id($_)->seq,"\n";
	}
	close OUT;
	return $read_file;
}

sub run_cap3
{
	my ($cap3_bin, $reads) = @_;
	my $command = join(" ", ($cap3_bin, $reads, ">>$reads.log"));
	system($command) == 0 or die "Command: $command failed !\n";
	my $cap_out = $reads . '.cap.contigs';
	return $cap_out;
}

sub unique
{
	my @array = @_;
	my %hash = map{$_, 1} @array;
	return keys %hash;
}

sub print_time_comment
{
	my $comment = shift;
	my $t = localtime(time);
	print STDERR $t,"\t", $comment,"\n";
}








