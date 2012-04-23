#!/usr/bin/perl -w
use strict;
use DBI;

my $sUsage = "perl $0 <db> <table> <user> <pass>\n";
die $sUsage unless @ARGV >= 4;
my($db, $table, $user, $pass) = @ARGV;
my $dbh = connect_to_db($db, $user, $pass);

my @usda_accs = qw(CAP11 	CAP12 	CAP7 	CAP8 	CS 	Ex 	JD 	JG 	Ku 	Ra);
my @usda_index = 9..18;
my @australia_accs = qw(JUP 	ACT 	GCT 	ID0 	JAY 	TCT 	GAT 	GRA 	TAT 	TRU 	ATT 	STP 	CAT 	TGT 	STE 	AGT 	CGT 	CAL 	GTT 	CTT);
my @australia_index = 19..38;

my $query = run_query($dbh, "select * from $table");
my ($in_us, $in_au, $in_all) = (0, 0, 0);
my %record;
while (my @data = $query->fetchrow_array)
{
	my $snp_id = $data[0];
	my @us = @data[@usda_index];
	my @au = @data[@australia_index];
	foreach (@us)
	{
		next unless /x/;
		my ($allel_a, $allel_b) = split /x/, $_;
		if ($allel_a >0 or $allel_b>0)
		{
			$in_us++, $in_all++;
			$record{$snp_id}->[0]++;
			$record{$snp_id}->[2]++;
		}
	}
	foreach (@au)
	{
		next unless defined $_;
		next unless /x/;
		my ($allel_a, $allel_b) = split /x/, $_;
		if ($allel_a >0 or $allel_b>0)
		{
			$in_au++, $in_all++;
			$record{$snp_id}->[1]++;
			$record{$snp_id}->[2]++;			
		}
	}	
}
foreach my $snp (keys %record)
{
	foreach(0..2){$record{$snp}->[$_] = 0 unless defined $record{$snp}->[$_]}
	print $snp, "\t", join("\t",(@{$record{$snp}})),"\n";
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



