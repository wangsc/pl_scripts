#!/usr/bin/perl 
use strict;
use Parallel::ForkManager;

my @array;

my $pm = Parallel::ForkManager->new(10);
foreach my $ind (1..100)
{
	$pm->start and next;
	push @array, $ind;
	$pm->finish
}

print join("\t", @array), "\n"
