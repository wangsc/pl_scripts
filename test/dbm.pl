#!/usr/bin/perl -w
use strict;

my %test;
my $mode = '0666';
dbmopen(%test, "test_dbmfile", $mode);
foreach (1..100000)
{
	$test{$_} = [1, 2];
}
sleep(50);
dbmclose %test;