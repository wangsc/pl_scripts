#!/usr/bin/perl -w
use strict;

my ($cd_hit, $barcode_sep_dir, $cd_hit_out_dir, $cd_hit_out_file, $similarity) = @ARGV;

mkdir ($cd_hit_out_dir) unless -d $cd_hit_out_dir;

my @fasta_files = <$barcode_sep_dir/*.fasta>;
my $all_in_one = $cd_hit_out_dir . "/all_in_one.fasta";
system("cat @fasta_files > $all_in_one");

# check cd-hit
my $t = `which $cd_hit`;
unless ($t)
{
	die <<END;
	Error: cannot call $cd_hit!
	Either cd-hit not installed or not in the PATH!
END
}

die if system("$cd_hit -i $all_in_one -c $similarity -o $cd_hit_out_file -M 0 -T 0 -d 0");