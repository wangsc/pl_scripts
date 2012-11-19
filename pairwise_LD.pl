#!/usr/bin/perl -w
# by SW
use strict;
my $sUsage = qq(
perl $0
<genotype file>
);
die $sUsage unless @ARGV;
my ($genotype_file, $out_file) = @ARGV;

my %genotype_vector = read_genotype_file($genotype_file);
my %r2 = calculate_LD(\%genotype_vector);

# Output
foreach my $ma (keys %r2)
{
	foreach my $mb (keys %{$r2{$ma}})
	{
		print join("\t", ($ma, $mb, $r2{$ma}{$mb})), "\n";
	}
}


sub read_genotype_file
{
	my $file = shift;
	open (IN, $file) or die "can't open file $file\n";
	my %return;
	#my @markers;
	my %marker_index;
	my %code = ("A", 0, "B", 2, "H", 1);
	my $n = 0;
	my $marker_start = 1;
	while(<IN>)
	{
		chomp;
		$n++;
		my @t = split /\s+/,$_; 
		if(/lines/ or $n == 1) #wsnp_AJ612027A_Ta_2_1   wsnp_AJ612027A_Ta_2_5   wsnp_be352570B_Ta_2_1   wsnp_BE398523A_Ta_2_1 
		{			
			my %markers = map{$t[$_],$_}0..$#t;
			%marker_index = reverse %markers;
			next
		}
		
		#0 1 2 0 2 0 2 2 0 NA 1 0 0 2
		foreach ($marker_start..$#t)
		{
			next unless exists $marker_index{$_};
			my $marker = $marker_index{$_};
			
			my $code;
			if($t[$_] eq 'NA')
			{
				$code = 4;
			}
			elsif ($t[$_] == 0 or $t[$_] eq 'A')
			{
				$code = 1 
			}
			elsif($t[$_] == 1 or $t[$_] eq 'H')
			{
				$code = 3
			}
			elsif($t[$_] == 2 or $t[$_] eq 'B')
			{
				$code = 2
			}
			else
			{
				$code = 4;
			}
			
			push @{$return{$marker}}, $code;
		}
	}
	close IN;
	return %return;
}

sub calculate_LD {
	my $hashref = shift;
	my %genotypeVec = %$hashref;
	my @sites = sort(keys(%genotypeVec));
	my @vector1 = ();
	my @vector2 = ();
	my %r2;
	my $genotype1 = 0;
	my $genotype2 = 0;
	my $i = 0;
	my $j = 0;
	my $k = 0;
	for ($i = 0; $i < @sites; $i++)
	{
		$r2{$sites[$i]}{$sites[$i]} = 1;
		@vector1 = @{$genotypeVec{$sites[$i]}};
		for($j = $i+1; $j < @sites; $j++) {
			@vector2 = @{$genotypeVec{$sites[$j]}};

			if (@vector1 != @vector2) {
				print STDERR join("\t", @vector1), "\n", join("\t", @vector2), "\n";
				die "Something is wrong with the prettybase input\n";
			}

			my $a1a1b1b1 = 0;
			my $a1a1b2b2 = 0;
			my $a1a1b1b2 = 0;
			my $a2a2b1b1 = 0;
			my $a2a2b2b2 = 0;
			my $a2a2b1b2 = 0;
			my $a1a2b1b1 = 0;
			my $a1a2b2b2 = 0;
			my $a1a2b1b2 = 0;
			for ($k = 0; $k < @vector1; $k++) {
				$genotype1 = $vector1[$k];
				$genotype2 = $vector2[$k];

				if ($genotype1 == 1) {
					if ($genotype2 == 1) {
						$a1a1b1b1++;
						# some3By3[0][0]++;
					}
					elsif ($genotype2 == 2) {
						$a1a1b2b2++;
						# some3By3[2][0]++;
					}
					elsif ($genotype2 == 3) {
						$a1a1b1b2++;
						# some3By3[1][0]++;
					}
				}
				elsif ($genotype1 == 2) {
					if ($genotype2 == 1) {
						$a2a2b1b1++;
						# some3By3[0][2]++;
					}
					elsif ($genotype2 == 2) {
						$a2a2b2b2++;
						# some3By3[2][2]++;
					}
					elsif ($genotype2 == 3) {
						$a2a2b1b2++;
						# some3By3[1][2]++;
					}
				}
				elsif ($genotype1 == 3) {
					if ($genotype2 == 1) {
						$a1a2b1b1++;
						# some3By3[0][1]++;
					}
					elsif ($genotype2 == 2) {
						$a1a2b2b2++;
						# some3By3[2][1]++;
					}
					elsif ($genotype2 == 3) {
						$a1a2b1b2++;
						# some3By3[1][1]++;
					}
				}
			}

			my $n = $a1a1b1b1 + $a1a1b1b2 + $a1a1b2b2 +
			        $a1a2b1b1 + $a1a2b1b2 + $a1a2b2b2 +
			        $a2a2b1b1 + $a2a2b1b2 + $a2a2b2b2; 
			next if $n == 0;
			my $x11 = 2*$a1a1b1b1 + $a1a1b1b2 + $a1a2b1b1;
			my $x12 = 2*$a1a1b2b2 + $a1a1b1b2 + $a1a2b2b2;
			my $x21 = 2*$a2a2b1b1 + $a1a2b1b1 + $a2a2b1b2;
			my $x22 = 2*$a2a2b2b2 + $a1a2b2b2 + $a2a2b1b2;

			my $p = ($x11 + $x12 + $a1a2b1b2) / (2*$n);
			my $q = ($x11 + $x21 + $a1a2b1b2) / (2*$n);

			my $p11 = $p * $q;

			my $convergentCounter = 0;
			my $oldP11 = $p11;
			my $range = 0.0;
			my $converged = "false";
			if ($p11 > 0) 
			{
				while ($converged eq "false" && $convergentCounter < 100) 
				{
					if ((1.0 - $p - $q + $p11) != 0.0 && $oldP11 != 0.0) 
					{
						$p11 = ($x11 + (($a1a2b1b2 * $p11 * (1.0 - $p - $q + $p11))/($p11 * (1.0 - $p - $q + $p11) + ($p - $p11)*($q - $p11))))/(2.0*$n);
						$range = $p11/$oldP11;
						if (($range >= 0.9999) && ($range <= 1.001)) 
						{
							$converged = "true";
						}
						$oldP11 = $p11;
						$convergentCounter++;
					}
					else 
					{
						#$p11 = $p * $q;
						$converged = "true";
					}
				}
			}

			# calculate D
			my $dValue = 0.0;
			if ($converged eq "true")
			{
				$dValue = $p11 - ($p * $q);
			}

			#calculate r2
			if ($dValue != 0.0)
			{
				$r2{$sites[$i]}{$sites[$j]} = ($dValue**2) / ($p * $q * (1 - $p) * (1 - $q));
			}
			else
			{
				$r2{$sites[$i]}{$sites[$j]} = 0;
			}
		}
	}
	return %r2;
}