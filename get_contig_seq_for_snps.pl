#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
my $sUsage = "perl $0 
<snp list file> 
<directory contained the contig fasta files> 
<all9SNPs_after_Cover_sel.txt>\n";
die $sUsage unless @ARGV >= 3;
my ($snplist, $ctg_direc, $all9snps_file) = @ARGV;

my @snps = read_snp_list($snplist);
my %all9snps = read_all9snps_file($all9snps_file);

my $ctg_file_suffix = 'cont_in.sanger.fasta';
foreach my $snpid (@snps)
{
	# wsnp_RFL_Contig4080_4587633; wsnp_CAP12_rep_c5274_2398727;
	chomp $snpid;
	next if $snpid =~ /^\s+$/;
	next if $snpid =~ /wsnp_b/i;
	next if $snpid =~ /wsnp_rfl/i;
	next if $snpid =~ /wsnp_Ta/;
	my @data = split /_/, $snpid;
	shift @data;
	my $acc_id = shift @data;
#	print STDERR 'acc_id: ', $acc_id,"\n";
	my $filename = join('_', ($acc_id, $ctg_file_suffix));
	$filename = $ctg_direc . '/' . $filename;
	next unless -e $filename;
	#print STDERR "not exists file $filename\n" unless -e $filename;
	pop @data;
	my $contig_id = (join('_', ('mira1',@data)));
#	print STDERR '$contig_id: ', $contig_id,"\n";
	my $ctg_obj = Bio::DB::Fasta->new($filename);
#	print STDERR 'ID number: ', scalar($ctg_obj->ids()),"\n";
	my $contig_seq = $ctg_obj->get_Seq_by_id($contig_id)->seq;
	my $id_nowsnp = $1 if $snpid =~ /wsnp_(.*)/;
#	print STDERR '$id_nowsnp: ', $id_nowsnp,"\n";
	my $flag =1;
	foreach (keys %all9snps)
	{
		
		if (/$id_nowsnp/ and /$acc_id/)
		{
	#		print STDERR 'all9snps key: ', $_,"\n" if $flag;
			my $flanking_seq = $all9snps{$_};
			my @data = split /_/, $_; #Ex_c35_77710_75789
			my $snp_pos = $data[-2] - $data[-1] + 1;
			print join("\t", ($snpid, $contig_seq, $snp_pos, $flanking_seq)),"\n";
			last;
		}
	}
}

sub read_snp_list
{
	my $file = shift;
	my @return;
	open (IN,"<$file") or die "can't open file $file\n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		push @return, $_;
	}
	close IN;
	return @return;
}

sub read_all9snps_file
{
	my $file = shift;
	open (IN,"<$file") or die "can't open file $file\n";
	my %return_hash;
	my $key;
	while (<IN>)
	{
		chomp;
		next if /^\s+$/;
		if (/mira1/)
		{
			my @data = split /\s+/, $_;
			my $snp_pos = $data[3] - $data[1] + 1;
			$key = join('_', (@data[0, 3, 1]));
			$key =~ s/mira1_//;
			next; 
		}
		else
		{
			$return_hash{$key} = $_;
		}
	}
	close IN;
	return %return_hash;
}














