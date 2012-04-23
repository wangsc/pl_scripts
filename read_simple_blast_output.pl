#!/usr/bin/perl -w
use strict;

my $file = shift;
open (IN, "$file") or die $!;
my %record;
my $name;
while (<IN>)
{
	next if /^\s+$/;
	if (/Query= (\S+)/)
	{
		$name = $1;
		$record{$name} = [];
		next;
	}
	my $spec = $1 if /\| (\S+ \S+) /;
	unless (defined $spec){$spec = $1 if /\|TAU76745 (\S+ \S+) /}
	print STDERR $_ unless defined $spec;
	
	push @{$record{$name}}, $spec;
}
close IN;

my %venn;
foreach my $id (keys %record)
{
	foreach my $spec (@{$record{$id}})
	{
		push @{$venn{$spec}}, $id;
	}
}
map{print $_,"\n", join("\n", unique(@{$venn{$_}})),"\n"}keys %venn;
print "\n";

my %count;
foreach my $id (keys %record)
{
	foreach (unique(@{$record{$id}}) )
	{
		$count{$_}++;
		#print STDERR $id,"\n" if $_ =~ /Oryza/;
	}
}

foreach (keys %count)
{
	print $_, "\t", $count{$_},"\n";
}

sub unique
{
	my %h = map{$_, 1} @_;
	return (keys %h);
}