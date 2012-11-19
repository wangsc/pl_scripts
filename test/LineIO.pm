package LineIO;
use strict;
use Tie::File;


sub new
{
	my ($class, $file) = @_;
	#print STDERR join("\t", @_), "\n";
	#print STDERR $file, "\n";
	open (my $FH, $file) or die "can't open file $file\n";
	my $ref = {FH => $FH, cnt => 0};
	bless $ref, $class;
	return $ref;
}


sub next_line
{
	my $class = shift;
	my $FH = $class->{"FH"};
	while(<$FH>)
	{
		return $_;
	}
}

1;