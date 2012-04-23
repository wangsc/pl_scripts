#!/usr/bin/perl -w
use DBI;

my ($db, $login, $password) = ('snp_database', 'browse_only', 'browse_only');
my $dbh = connect_to_db($db, $login, $password);
my $tables = run_query($dbh, "show tables");
my @tbls;
while (my @t = $tables->fetchrow_array)
{
	push @tbls, shift @t;
}
my @head = qw(Field Type Null Key);
foreach my $t (@tbls)
{
  print "\n", 'Table: ', $t,"\n";
  print join("\t", @head),"\n";
	my $query = "show columns from $t";
	my $run = run_query($dbh, $query);
	while (my @t = $run->fetchrow_array)
	{
		next unless @t;
		print join("\t", @t),"\n";
	}

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