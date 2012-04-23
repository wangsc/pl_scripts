#!/usr/bin/perl -w
use strict;
use DBI;
use File::Basename;
#
my $sUsage = qq(
perl $0
<db>
<table>
<mysql user>
<mysql passwd>
<SNPCounts.out files>
);
die $sUsage unless @ARGV > 4;
my $db = shift;
my $table = shift;
my $my_user = shift;
my $my_pass = shift;
my @files = @ARGV;
my $dbh = connect_to_db($db, $my_user, $my_pass);
# Main
grep {print $_ ,"\n"} @files;
my %table_snp_allele = get_snp_allele($dbh, $table);

my %acc_hash = read_snpcount_files(\%table_snp_allele, @files);

update_table($dbh, $table, \%acc_hash);



# subroutines

sub update_table
{
	my ($dbh, $table, $acc_hashref) = @_;
	my @accs = keys %$acc_hashref;
#	add_colums_to_table($dbh, $table, @accs);
	
	foreach my $acc (@accs)
	{
		
		foreach (@{$acc_hashref->{$acc}})
		{
			my ($snpname, $count) = @$_;
			my $query_txt = "UPDATE $table SET  \`$acc\` = \'$count\' WHERE \`snpName\` = \'$snpname\'";
			print $query_txt, "\n";
			run_query($dbh, $query_txt);
		}
	}
}

sub get_snp_allele
{
	my ($dbh, $table) = @_;
	my %return;
	my $query = run_query($dbh, "select snpName, SNP1, allele from $table");
	
	while (my @data = $query->fetchrow_array)
	{
		my ($snpname, $snppos, $allele) = @data;
		$allele =~ s/\[|\]//g;
		$allele =~ s/x/\//;
#		print '$snpname ', $snpname,"\n";
		my @t= split /_/, $snpname;
		my $c = pop @t;	
#		print 'New snpname ', join('_', (@t, $snppos)), "\n";
		$return{join('_', (@t, $snppos))} = [$allele, $c];
	}
	return %return;
}


sub add_colums_to_table
{
	my ($dbh, $table, @accs) = @_;
#	my $query = prep_query($dbh, "alter table $table add column ? varchar(20) NULL after 'Ra'") or die $!;
	foreach my $id (@accs)
	{
#		my $query_txt = " if not exists (SELECT * FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = $table AND COLUMN_NAME = $id)
#											BEGIN alter table \`$table\` add \`$id\` varchar(20) NUll after \`Ra\`;
#											END";
		my $query_txt = "alter table \`$table\` add \`$id\` varchar(20) NUll after \`Ra\`";
#		print $query_txt,"\n";
		run_query($dbh, $query_txt);
#		$query->execute($id);
	}
	return 1;
}


sub read_snpcount_files
{
	my $table_snp_allele = shift;
	my @files = @_;
	my %acc_hash;
	open (T, ">temp.out") or die $!;
	foreach my $file (@files)
	{
		my $basename = basename($file);
		my @t = split /_/, $basename;
		my $acc_id = shift @t;
		open (IN,"$file") or die "can't open file $file \n";
		while (<IN>)
		{
			next if /^\s+$/;
			print T $_;
			chomp;
			my @line_data = split /\s+/, $_; # Ra_mira1_c6065_10655997	T/C	0/0
		#	print join('***', @line_data),"\n";
			$acc_hash{$acc_id} = [] unless exists $acc_hash{$acc_id};
			my @t = split /_/, $line_data[0];
			my $acc = shift @t;
			my $snp_pos = pop @t;
			my $contig = join('_', @t);
			unless (exists $table_snp_allele->{$line_data[0]}) 
			{
		#		print 'Not exists ', $line_data[0],"\n"; 
				next
			}
	#		print 'Exists: ', $line_data[0],"\n" if exists $table_snp_allele->{$line_data[0]};
			$line_data[2] = reverse $line_data[2] unless $line_data[1] eq $table_snp_allele->{$line_data[0]}->[0];
			$line_data[2] =~ s:\/:x:;
			my $id = join('_', ($acc, $contig, $table_snp_allele->{$line_data[0]}->[1]));
	#		print 'snpName: ', $id,"\n";
			push @{$acc_hash{$acc_id}}, [$id, $line_data[2]];
		}
		close IN;
	}
	close T;
	return %acc_hash;
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