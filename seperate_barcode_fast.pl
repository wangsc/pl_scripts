#!/usr/bin/perl -w
use strict;

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
<barcode file>
<keep barcode in file (1) or not keep (0)>
<fastq files>

**********************************************************************************
);
die $sUsage unless @ARGV >= 3;

my ($barcode_file, $keep_barcode, @fastq_files) = @ARGV;
my %barcode_ids = read_barcode_file($barcode_file);
my @accessions = unique(values %barcode_ids);
my %acc_fhs = get_file_handles(@accessions);
my ($read_code, $read_noncode) = seperate_barcodes(\%acc_fhs, \%barcode_ids, $keep_barcode, @fastq_files);
print '#Reads with barcodes: ', $read_code,"\n";
print '#Reads without barcodes: ', $read_noncode,"\n";

# Subroutines
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
	my @accs = @_;
	my %return_hash;
	foreach my $acc_id (@accs)
	{
		my $file = $acc_id . '.fastq';
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
		my @aux = undef;
		while(my ($name, $seq, $qual) = readfq(*IN, \@aux))
		{
			$name = '@'.$name;
			my $barcode = detect_barcode($seq, $barcode_ref);
			unless (defined $barcode){$read_noncode++; next}
			my $code_length = length $barcode;
			$read_code++;
			unless($keep_barcode)
			{
				#$seq =~ s/^.{$code_length}//;
				substr($seq, 0, $code_length) = '';
				#$qual =~ s/^.{$code_length}//;
				substr($qual, 0, $code_length) = '';
			}
			print {$fhs_ref->{$barcode_ref->{$barcode}}} join("\n", ($name, $seq, '+', $qual)),"\n" if (length $seq) >= 30;
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
		#	my $code_rc = rev_comp($code);
		#	my $code_reverse = reverse $code;
		#	my $code_rc_reverse = reverse $code_rc;
		my $prefix_seq = substr($seq, 0, length($code));	
		if($prefix_seq eq $code)
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
	print STDERR $t,"\t", $c,"\n";
}


# from lh3, https://github.com/lh3/readfq/blob/master/readfq.pl
sub readfq 
{
	my ($fh, $aux) = @_;
	@$aux = [undef, 0] if (!defined(@$aux));
	return if ($aux->[1]);
	if (!defined($aux->[0])) {
		while (<$fh>) {
			chomp;
			if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
				$aux->[0] = $_;
				last;
			}
		}
		if (!defined($aux->[0])) {
			$aux->[1] = 1;
			return;
		}
	}
	my $name = /^.(\S+)/? $1 : '';
	my $seq = '';
	my $c;
	$aux->[0] = undef;
	while (<$fh>) {
		chomp;
		$c = substr($_, 0, 1);
		last if ($c eq '>' || $c eq '@' || $c eq '+');
		$seq .= $_;
	}
	$aux->[0] = $_;
	$aux->[1] = 1 if (!defined($aux->[0]));
	return ($name, $seq) if ($c ne '+');
	my $qual = '';
	while (<$fh>) {
			chomp;
			$qual .= $_;
			if (length($qual) >= length($seq)) {
			$aux->[0] = undef;
			return ($name, $seq, $qual);
		}
	}
	$aux->[1] = 1;
	return ($name, $seq);
}



