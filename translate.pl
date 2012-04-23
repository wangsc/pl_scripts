#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;

my ($format) = 'fasta';
my $frame;
GetOptions(
'format:s' => \$format,
'frame:i' => $frame
);

my $oformat = 'fasta';
$frame = 0 unless defined $frame;
# this implicity uses the <> file stream
my $seqin = Bio::SeqIO->new( -format => $format, -file => shift);
my $seqout = Bio::SeqIO->new( -format => $oformat, -file => ">-" );


while( (my $seq = $seqin->next_seq()) ) {
my $pseq = $seq->translate(-frame=>$frame, -terminator=>'-');
$seqout->write_seq($pseq);
}