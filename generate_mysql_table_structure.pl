#!/usr/bin/perl -w
use strict;
use DBI;
use CGI qw(:standard);

my $sUsage = "perl $0 <database> <user> <passwd>\n";
die $sUsage unless @ARGV >= 3;

my ($db, $user, $pass) = @ARGV;
my $dbh = connect_to_db($db, $user, $pass);

my $query = run_query($dbh, "show tables");
my @tables;
while (my ($table_name) = $query->fetchrow_array)
{
	push @tables, $table_name;
}
@tables = sort{$a cmp $b} @tables;
print start_html("Tables in database [$db]"),
			h4("Tables in database [$db]")			
			;
print "<ul>","\n";
foreach my $table (@tables)
{
	print "<li>","\n";
	
	print $table,"\t";
	print "<a href = /show_columns.php?table=$table target=\"showframe\">", 'Show columns', '</a><p>',"\n";
	
	print "</li>","\n";
}
print "<form action=\"show_columns.php\" method=\'get\'> </form>", "\n";
print "</ul>","\n";
print end_html;




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