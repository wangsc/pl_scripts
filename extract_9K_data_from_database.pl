#!/usr/bin/perl -w
use strict;
use DBI;

my $sUsage = qq(
perl $0
<9K design file, PrivKSU_WheatCons_9k_11497518_A.csv>
<Database Name>
<table Name>
<username for DB connection>
<password for DB connection>
);

die $sUsage unless @ARGV >= 5;

my ($design_file, $db, $table, $username, $passwd) = @ARGV;
my @k9_ids = read_design_file($design_file);
my $dbh = connect_to_db($db, $username, $passwd);

my $prep_query = prep_query ($dbh, "select * from $table where name = ?");

foreach my $id ( @k9_ids)
{
	$prep_query->execute($id);
	while(my @tmp = $prep_query->fetchrow_array)
	{
		print join("\t", @tmp), "\n";
	}
}



# 
sub read_design_file
{
	my $file = shift;
	open (IN, $file) or die "$file: ", $!;
	my @return;
	while (<IN>)
	{
		s/\"//g;
		next if /^ID/;
		next if /^\s+$/;
		my @t = split /\s+/, $_;
		next unless /wsnp/;
		my $rep = $t[2] =~ s/wsnp_//g;
		push @return, $t[2] if $rep;
	}
	close IN;
	return @return;
}



sub connect_to_db
{
	my ($db, $login, $password) = @_;
	my $dbh = DBI->connect("dbi:mysql:". $db, $login, $password)
	  or die "Can't connect to MySQL database: $DBI::errstr\n";
}

sub prep_query
{
	my ($dbh, $query_text) = @_;
	my $query = $dbh->prepare($query_text) or die "Couldn't prepare query: " . $dbh->errstr;
}

sub run_query
{
	my ($dbh, $query_text) = @_;
	my $query = prep_query(@_);
	$query->execute() or die "Couldn't execute query: $query_text\t" . $query->errstr;
	return $query;
}