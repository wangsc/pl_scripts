#!/usr/bin/perl -w
use strict;

# given a number, output the next higher number with the same digits
# For example, given 15432, output 21345;



while(1)
{
	print "Please input a number (>=10) or quit: \n";
	chomp (my $input = <STDIN>);
	last if $input =~ /^q/i;
	unless ($input=~/\d+/ and $input >= 10){print "Error: Number must >= 10\n"; next}
	$input += 0;
	my @digits = split //, $input;
	print "Next higher number is: ", next_highest_num(@digits), "\n\n";
}


sub next_highest_num
{
	my @digits = @_;
	#print @digits, "\n";
	my @new_number_digits;
	my $index = 0;
	if(@digits == 2)
	{
		return $digits[1]>$digits[0]?@digits[1,0]:@digits;
	}
	else
	{
		if(descendent(@digits[$index+1..$#digits]))
		{
			print @digits[$index+1, $#digits], "\tDescendent\n";
			my $swap = find_swap_index(@digits);
			@digits[0, $swap] = @digits[$swap, 0];
			return ($digits[0], sort{$a<=>$b} @digits[1..$#digits]);
		}
		else
		{
			my $d = shift @digits;			
			push @new_number_digits, ($d, next_highest_num(@digits));
		}		
	}
	return @new_number_digits;

}

sub find_swap_index
{
	my @arr = @_;
	my %hash;
	foreach my $index (1..$#arr)
	{
		$hash{$index} = $arr[$index] if $arr[$index] > $arr[0];
	}
	my @t = (sort {$hash{$a}<=>$hash{$b}} keys %hash);
	return $t[0];
}

sub descendent
{
	my @arr = @_;
	foreach (0..$#arr-1)
	{
		return 0 unless $arr[$_+1] <= $arr[$_];
	}
	return 1;
}