#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Tools::Run::RemoteBlast;
use File::Basename;
#
my $sUsage = qq"perl $0 
<blast program>
<blast database, nt, nr, refseq_protein, swissprot, etc al>
<input fasta file>
<output type, default is 0, set 5 for xml output> 
< e value (optional), default is 1e-10>\n";
die $sUsage unless @ARGV >= 3;
my $prog = shift;
my $db = shift;
my $input_fastafile = shift;
my $outtype = shift;
$outtype = 0 unless defined $outtype;
my $e = shift;
my $e_val = (defined $e)?$e:'1e-10';
my @params = ('-prog' => $prog, '-data' => $db, '-expect' => $e_val, '-readmethod' => 'SeachIO');
@params = (@params, '-readmethod' => 'xml') if $outtype == 5;
my $factory = Bio::Tools::Run::RemoteBlast->new(@params);
$factory->retrieve_parameter('FORMAT_TYPE', 'XML') if $outtype == 5;
# Main
my $output_blastfile = basename($input_fastafile, ".fasta") . "\.$prog". "_". $db . '.out';
run_blast($factory, $prog, $input_fastafile, $output_blastfile);



# Subroutines
sub run_blast
{
	my ($factory, $prog, $in_file, $output_filename) = @_;
	my $str = Bio::SeqIO->new(-file=>$in_file , -format => 'fasta' );
	my $count_seq = 0;
	open (my $OUT, ">$output_filename") or die "can't open file $output_filename\n";
	while (my $input = $str->next_seq())
	{
		my $r = $factory->submit_blast($input);
		while (my @rid = $factory->each_rid)
		{
			foreach my $rid (@rid)
			{
				
				my $rc = $factory->retrieve_blast($rid);
				if (! ref($rc))
				{
					$factory->remove_rid($rid) if $rc < 0;
					sleep 5;
				}
				else
				{
					my $result = $rc->next_result();
					next unless defined $result;
					my $query_name = $result->query_name();
#					print 'Query_name: ', $query_name, "\n";
					my $filename = 'temp.out';
					$factory->save_output($filename);
					$factory->remove_rid($rid);
					save_tmp_out($OUT, $filename);
				}
			}
			$count_seq++;
			print_time("Finished blasting $count_seq sequences...") unless $count_seq%500;
		}		
	}

}

sub save_tmp_out
{
	my ($OUT, $tmpfile) = @_;
	open (IN, "$tmpfile") or die "can't open file $tmpfile\n";
	while (<IN>)
	{
		print $OUT $_;
	}
	close IN;
}

sub print_time
{
	my $text = shift;
	my $time = localtime;
	print $time,"\t", $text, "\n";	
}

