#!/usr/bin/perl -w
use strict;

my $num_pairs = shift or die "perl $0 num_pairs\n";
my  @str;
if($num_pairs > 0)
{
	printParenthesis(0, $num_pairs, 0, 0);
}

sub printParenthesis
{
	my ($pos, $num, $open, $close) = @_;	
	print join("\t", @_), "\n";
	if ($close == $num)
	{
		print "\t", join("", @str), "\n";
		return;
	}
	else
	{
		
		if($open > $close)
		{
			#print "\t", join("*", ($open, $close, $pos)), "\n";
			$str[$pos] = '}';
			print "*\t", join("", @str), "\n";
			printParenthesis($pos+1, $num, $open, $close+1);
		}
		if($open < $num)
		{
			#print "\t", join("+", ($open, $close, $pos)), "\n";
			$str[$pos] = '{';
			print "#\t", join("", @str), "\n";
			printParenthesis($pos+1, $num, $open+1, $close);
		}
	}
	
}



=head
# C code:
void _printParenthesis(int pos, int n, int open, int close)
{
  static char str[MAX_SIZE];

  if(close == n)
  {
    printf("%s \n", str);
    return;
  }
  else
  {
    if(open > close) {
        str[pos] = '}';
        _printParenthesis(pos+1, n, open, close+1);
    }
    if(open < n) {
       str[pos] = '{';
       _printParenthesis(pos+1, n, open+1, close);
    }
  }
}
=cut