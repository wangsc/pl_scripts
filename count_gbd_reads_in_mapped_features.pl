#!/usr/bin/perl -w
use strict;

my $sUsage = qq(
perl $0
<gff3 file>
<read alignment file>
<output>
);
die $sUsage unless @ARGV >= 3;
my ($gff_file, $align_file, $out_file) = @ARGV;
