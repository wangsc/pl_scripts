#!/usr/bin/perl -w
use strict;
use Statistics::R;
use File::Basename;
use Cwd;

my $sUsage = "perl $0 [matrix files]\n";
die $sUsage unless @ARGV;
my @files = @ARGV;
my $cwd = getcwd;
my $plot_R_file = "plot.R";
open (IN, ">$plot_R_file") or die;
foreach my $file (@files)
{
	my $cmd = qq(perl -i -pe \047 s/^\\d+\$//g; s/^\\s+\$// \047 $file);
	system($cmd);
	#print $cmd,"\n";exit;
	my $basename = basename($file, ".matrix");
	my $full_path = $cwd . '/' . $file;
	my $read_table = "mt <- read.table($full_path)";
	my $r_cmds = ""
	
}

sub read_matrix_file
{
	
}