#!/usr/bin/perl -w
use strict;


my @files = qw(A.test B.test C.test);
my %fh = get_file_handles(@files);


foreach my $index (0..$#files)
{
	print {$fh{$index}} $index;
}


sub get_file_handles
{
	my @files = @_;
	my %return_hash;
	foreach (0..$#files)
	{
		local *FH;
		open (FH, ">$files[$_]") or die "can't open file $files[$_]\n";
		$return_hash{$_} = *FH{IO};
	}
	return %return_hash;
}

foreach (values %fh)
{
	close $_ or print "Can't close File handle $_\n";
}