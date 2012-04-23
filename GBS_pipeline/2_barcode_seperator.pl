#!/usr/bin/perl -w
use strict;
use File::Basename;
use Parallel::ForkManager;

my $sUsage = qq(
**********************************************************************************
Seperate fastq files into accession-specific files based on barcode file.

The barcode file should look like this:
ID1  ARCAGT
ID2  ATCGATG
....
One ID could have several barcodes, but each barcode shoud be unique in this file.

Usage:
perl $0
<directory contained the fastq files>
<barcode file>
<direcotry for output files>

**********************************************************************************
);
die $sUsage unless @ARGV >= 3;
my ($qc_out_dir, $barcode_file, $barcode_sep_dir) = @ARGV;
my @fastq_files = <$qc_out_dir/*_trimmed_filtered>;
my $keep_barcode = 0;
my %barcode_ids = read_barcode_file($barcode_file);
my @accessions = unique(values %barcode_ids);
my %acc_fhs = get_file_handles($barcode_sep_dir, @accessions);
my ($read_code, $read_noncode) = seperate_barcodes(\%acc_fhs, \%barcode_ids, $keep_barcode, @fastq_files);
print "\t#Reads with valid barcodes: ", $read_code,"\n";
print "\t#Reads without valid barcodes: ", $read_noncode,"\n";

# fastq to fasta
my $MAX_PROCESSES = 5;
my $pm = new Parallel::ForkManager($MAX_PROCESSES);
my @sep_files = <$barcode_sep_dir/*.fastq>;
foreach my $file (@sep_files)
{

	my $pid = $pm->start and next;
	my $out = $file;
	$out =~ s/fastq$/fasta/;
	&fastq_to_fasta($file, $out);
	$pm->finish;
}

# Subroutines
sub fastq_to_fasta
{
	my($fastq_file, $fasta_file) = @_;
	my $base = basename($fastq_file);
	my $acc = $1 if $base =~ /(\S+)\.fastq/;
	open (IN, "$fastq_file") or die;
	open (OUT, ">$fasta_file") or die;
	my $count=0;
	my @temp;
	my $read_num = 0;
	while(<IN>)
	{
		next if /^\s+$/;
		$count++;
		if($count==1){$read_num++; print OUT ">$acc", "_", $read_num, "\n";}
		elsif($count==2){print OUT $_}
		elsif($count==4){$count=0;next}
	}
	close IN;
	close OUT;
}


sub read_barcode_file
{
	my $file = shift;
	my %return;
	open (IN, $file) or die ;
	while(<IN>)
	{
		chomp;
		next if /^\s+$/;
		s/^\s+//;
		my @t=split /\s+/,$_;
		$return{uc $t[1]} = $t[0];
	}
	close IN;
	return %return;
}

sub unique
{
	my %h = map{$_, 1} @_;
	return keys %h;
}


sub get_file_handles
{
	my $dir = shift;
	my @accs = @_;
	my %return_hash;
	foreach my $acc_id (@accs)
	{
		my $file = $dir . "/". $acc_id . '.fastq';
		local *FH;
		open (FH, ">$file") or die "can't open file $file\n";
		$return_hash{$acc_id} = *FH{IO};
	}
	return %return_hash;
}

sub seperate_barcodes
{
	my ($fhs_ref, $barcode_ref, $keep_barcode, @files) = @_;
	my ($read_code, $read_noncode) = (0, 0);
	foreach my $file (@files)
	{
		print_time_comment("Processing file $file");
		open (IN, $file) or die "can't open file $file\n";
		my @temp_array;
		my $index = -1;
		while(<IN>)
		{
			next if /^\s+$/;
			chomp;
			$index++;		
			$temp_array[$index] = $_;
			if($index == 3)
			{
				$index = -1;
				my $seq = $temp_array[1];
				my $barcode = detect_barcode($seq, $barcode_ref);
				unless (defined $barcode){$read_noncode++; next}
				my $code_length = length $barcode;
				$read_code++;
				unless($keep_barcode)
				{
					$temp_array[1] =~ s/^.{$code_length}//;
					$temp_array[3] =~ s/^.{$code_length}//;
				}
				# Only output reads with length >= 30 after barcode removed
				print {$fhs_ref->{$barcode_ref->{$barcode}}} join("\n", @temp_array),"\n" if (length $temp_array[1]) >= 30;
			}
		}
		close IN;
	}
	return($read_code, $read_noncode);
}

sub detect_barcode
{
	my ($seq, $hashref) = @_;
	my @barcodes = keys %$hashref;
	my @return;
	foreach my $code (@barcodes)
	{
		my $code_rc = rev_comp($code);
		my $code_reverse = reverse $code;
		my $code_rc_reverse = reverse $code_rc;
		if($seq =~ /^$code/)# or $seq =~ /^$code_rc/ or $seq =~/$code_reverse$/ or $seq =~/$code_rc_reverse$/)
		{
			push @return, $code;
		}
	}
	if(@return > 1)
	{
		print "Two barcode for one seq: \n";
		print $seq,"\n";
		print "Barcodes: ", join("\t", @return),"\n";
		#exit;
	}
	else
	{
		return $return[0];
	}
}

sub rev_comp
{
	my $seq = shift;
	$seq = reverse $seq;
	$seq =~ tr/ATGCatgc/TACGtacg/;
	return $seq;
}

sub print_time_comment
{
	my $c = shift;
	my $t = localtime(time);
	print STDERR "\t", $t,"\t", $c,"\n";
}
