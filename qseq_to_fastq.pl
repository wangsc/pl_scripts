#!/usr/bin/perl

use warnings;
use strict;

my $usage = "Usage: $0 <qseq input files> \n";
my @infiles = @ARGV;
die $usage unless @infiles ;
my $min_len = 30;
foreach my $infile (@infiles)
{
   open(IN,'<',$infile) || die "Could not open $infile: $!\n";
   while (<IN>) {
     chomp;
     #MachineID     run#     lane#     tile#     x-coord     y-coord     index     read#     sequence     q-scores    p/f flag
     my ($machine_id,$run,$lane,$tile,$x,$y,$index,$read,$seq,$qscore,$flag) = split(/\t/);
     $seq =~ s/\./N/g;
     next unless length $seq >= 30;
     print  "@","$machine_id:$lane:$tile:$x:$y#$index/$read\n";
     print  "$seq\n";
     print  "+","$machine_id:$lane:$tile:$x:$y#$index/$read\n";
     print  "$qscore\n";
   }

close(IN);
}
exit(0);
