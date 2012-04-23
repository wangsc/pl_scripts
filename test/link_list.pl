#!/usr/bin/perl -w
use strict;
use constant NEXT => 0;
use constant VAL => 1;

my $list;
my $tail = \$list;
foreach (1..5)
{
	my $node = [undef, $_*$_];
	$$tail = $node;
	$tail = \$node->[NEXT];
}
print $list->[NEXT]->[NEXT]->[VAL], "\n";