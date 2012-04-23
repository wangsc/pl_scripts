#!/usr/bin/perl -w
use strict;
use Cwd;
use CN::DB_Utilities;
use Bio::DB::Fasta;

my $sUsage = qq(
perl $0
<database name>
<talbe name>
<username for db>
<passwd for db>
<sorghum genome file, /bioinfo/sorghum/pipeline/genome/sbi1.fasta>
<output filename for the corrected db (could be loaded into table later)>
);
die $sUsage unless @ARGV;

my($db, $table, $user, $pass, $sorghum_fasta, $output) = @ARGV;
my $dbh = DB_Utilities::connect_to_db($db, $user, $pass);

# read genome file
my $sb_genome = Bio::DB::Fasta->new($sorghum_fasta);

# get records from table
# accession_id 	chromosome 	gen_position 	ref_base 	alt_base 	snp_score 	total_reads
my @records  =get_data_from_table($dbh, $table);

# correction
open (OUT, ">$output") or die "can't open file $output\n";
foreach my $line (@records)
{
	my @t = split/\t/,$line;
	my ($chr, $pos, $ref) = @t[1..3];
	my $correct_ref = $sb_genome->seq("chromosome_".$chr, $pos => $pos);
	@t[3,4] = ($correct_ref, $correct_ref) unless $ref eq $correct_ref;
	print OUT join("\t", @t),"\n";
}
close OUT;

# Make a new table and load data into database
my $newtable = $table . '_corrected_RR';
DB_Utilities::run_query($dbh, "create table $newtable like $table");
my $cwd = getcwd;
my $file = join('/', ($cwd, $output));
DB_Utilities::run_query($dbh, "load data local infile \'$file\' into table $newtable fields terminated by \'\t\' lines terminated by \'\n\'");

sub get_data_from_table
{
	my ($dbh, $table) = @_;
	my @records;
	my $query = DB_Utilities::run_query($dbh, "select * from $table");
	while(my @temp = $query->fetchrow_array)
	{
		push @records, join("\t", @temp);
	}
	return @records;
}

