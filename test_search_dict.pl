#!/usr/bin/perl -w
use strict;
use Search::Dict;

my $file = shift;
my $str = shift;
open (*IN, $file) or die;
look *IN, "$str";
my $line = <IN>;
print $line,"\n";
