#!/usr/bin/perl -w
use strict;
use GD::Simple;

my $sUsage = qq(
perl $0 
<syntenic position file>
<total length of each chromosome>
);
die $sUsage unless @ARGV >=2;
my ($syntenic_file, $chr_length_file) = @ARGV;