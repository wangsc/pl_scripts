#!/usr/bin/perl -w
use strict;
use GD::Simple;

my $sUsage = qq(
perl $0
<Tag aligned postion file, each line should be: TagId ContigID AlignPosition >
<Contig length file, each line should be: ContigID ContigLength >
<mapped TagID file; One TagID per line>
<output pic name (*.png)>
);
die $sUsage unless @ARGV >= 4;
my ($alignPos_file, $ctg_len_file, $mapped_tags_file, $out) = @ARGV;
open (OUT, ">$out") or die;
binmode(OUT);

my ($align_positions, $pos_tags) = read_position_file($alignPos_file);
my %ctg_length = read_ctg_length($ctg_len_file);
my %mapped_tags = rad_mapped_tags_file($mapped_tags_file);
my $ctg_num = scalar keys %ctg_length;
my $max_length = max(values %ctg_length);

# plot size
my $block_height = 50;
my $tick_size = $block_height/4;
my ($left_margin, $right_margin, $up_margin, $down_margin) = (10, 40, 10, 10);
my @plot_size = (800, $block_height*$ctg_num + $up_margin + $down_margin);
my $bp_pix = int(($plot_size[0] - $left_margin - $right_margin)/$max_length);

my $img = GD::Simple->new(@plot_size);

my $num = 0;
foreach my $ctg (sort {$a cmp $b} keys %ctg_length)
{
	my $block_len = $ctg_length{$ctg} * $bp_pix;
	my $solid_line_y = $up_margin + $num*$block_height + 0.75*$block_height;
	$img->fgcolor('blue');
	$img->rectangle($left_margin, $solid_line_y, $left_margin+$block_len-1, $solid_line_y+1);
	
	foreach my $pos (@{$align_positions->{$ctg}})
	{
		$img->fgcolor('red');
		$img->rectangle($pos*$bp_pix, $solid_line_y-$tick_size, $pos*$bp_pix, $solid_line_y+$tick_size);
		my $tag = $pos_tags->{$ctg}{$pos};
		if(exists $mapped_tags{$tag})
		{
			$img->fgcolor('black');
			$img->moveTo($pos*$bp_pix, $solid_line_y-$tick_size-3);
			$img->lineTo($pos*$bp_pix-int($tick_size/4), $solid_line_y-int(1.5*$tick_size));
			$img->lineTo($pos*$bp_pix+int($tick_size/4), $solid_line_y-int(1.5*$tick_size));
			$img->lineTo($pos*$bp_pix, $solid_line_y-$tick_size-2);
			my $red = $img->colorAllocate(0,0,0);
			$img->fill($pos*$bp_pix, $solid_line_y-$tick_size-4, $red);
			#$img->moveTo($pos*$bp_pix-($tick_size/2), $solid_line_y-1.7*$tick_size);
			#$img->lineTo($pos*$bp_pix+($tick_size/2), $solid_line_y-1.7*$tick_size);
			
		}
		
	}
	
	$img->moveTo($plot_size[0]-$right_margin*1.3, $solid_line_y+$tick_size/4);
	$img->fgcolor('black');
	$img->font('Times:bold');
	$img->fontsize($tick_size);
	$img->string($ctg);
	$num++;
}

print OUT $img->png;

#
sub read_position_file
{
	my $file = shift;
	my %return;
	my %pos_tag;
	open (IN, "$file") or die;
	while(<IN>)
	{
		# TagId ContigID AlignPosition
		chomp; 
		my @t = split /\s+/,$_; 
		push @{$return{$t[1]}}, $t[2];
		$pos_tag{$t[1]}{$t[2]} = $t[0];
	}
	close IN;
	return (\%return, \%pos_tag);
}

sub read_ctg_length
{
	my $file = shift;
	my %return;
	open (IN, "$file") or die;
	while(<IN>)
	{
		# ContigID ContigLength
		chomp; 
		my @t = split /\s+/,$_; 
		$return{$t[0]} = $t[1];
	}
	close IN;
	return %return;	
}

sub rad_mapped_tags_file
{
	my $file = shift;
	my %return;
	open (IN, "$file") or die;
	while(<IN>)
	{
		# ContigID ContigLength
		chomp; 
		$return{$_} = 1;
	}
	close IN;
	return %return;		
}

sub max
{
	my $max = shift;
	foreach (@_)
	{
		$max = $_ if $_  >$max;
	}
	return $max;
}

