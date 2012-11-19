# Test package
package Own::Test;

sub new
{
	my $class = shift;
	my $name = shift;
	my $r = {'name' => $name};
	bless $r, $class;
}

sub smile
{
	my $class = shift;
	print "I am happy :)\n";
}

sub cry
{
	print "I feel sad :(\n";
}

sub wheatever
{
	print "I am oook\n";
}

1;
