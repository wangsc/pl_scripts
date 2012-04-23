#!/usr/bin/perl -w
use strict;
use DBI;
# Update the SNP table by adding AA_change value
#
my $sUsage = qq(
perl $0
<database>
<table>
<user>
<passwd>
<aa_change file>
);
die $sUsage unless @ARGV >= 5;
my ($db, $table, $user, $pass, $aa_file) = @ARGV;
my $dbh = connect_to_db($db, $user, $pass);
my $update_query = "update $table set AA_change=? where name=?";
my $update = $dbh->prepare($update_query) or die "Couldn't prepare query: " . $dbh->errstr;

my %aa_change = read_file($aa_file);
while(my ($name, $aa_change) = each %aa_change)
{
	$update->execute($aa_change, $name);
}


# Sbroutines
sub read_file
{
	my $file = shift;
	my %return_hash;
	open (IN, "$file") or die "can't open file $file \n";
	while (<IN>)
	{
		next if /^\s+$/;
		chomp;
		my @line_data = split /\t/, $_;
		my $name = get_name($line_data[0]);
		my $aa_change = $line_data[-1];
		$return_hash{$name} = $aa_change;
	}
	close IN;
	return %return_hash;
}

sub get_name
{
	my $str = shift;#Ra_mira1_c4096_7494789_7496067_7495709_7495709_[A/G]
	my @array = split /_/, $str;
	return join('_', (@array[0, -6, -2]));	
}


sub connect_to_db
{
	my ($db, $login, $password) = @_;
	my $dbh = DBI->connect("dbi:mysql:". $db, $login, $password)
	  or die "Can't connect to MySQL database: $DBI::errstr\n";
}