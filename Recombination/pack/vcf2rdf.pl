#!/usr/bin/perl -w

use strict;
use warnings;

my $line;
my @column;
my $number_of_sample;
my @sample_id_name;
my @position = ('       ' , '       ' , '       ' , '       ' , '       ' , '       ' );
my %haplotype;
my $weight = '';
my $offset = 9; #number of column before each sample's data.
my $line_break = "\r\n";

for($line = <>; $line =~ /^##/; $line = <>) {}

if ($line =~ /^#CHROM/) {
	chomp($line);
	@column = split(/\t/, $line);
	@sample_id_name = @column[$offset .. $#column];
	$number_of_sample = @sample_id_name;
}
else {
	die "Input file error.";
}

my $pos;
while ($line = <>) {
	chomp($line);
	@column = split(/\t/, $line);
	$pos = sprintf( "%9d", $column[1]);
	for (my $i = 0; $i < 6; $i ++ ) {
		$position[$i] .= substr( $pos, $i + 3, 1 );
	}
	for (my $i = 0; $i < $number_of_sample; $i ++) {
		my $gt;
		if ($column[$offset + $i] =~ /^([^:]*)/ ) {
			$gt = $1;
		}
		my @vcf_gt = split ('\|', $gt);
		if (@vcf_gt == 0) {
			@vcf_gt = ('.', '.');
		}
		elsif (@vcf_gt == 1) {
			if ($vcf_gt[0] eq '.') {
				@vcf_gt = (0, 0);
			}
			else {
				$vcf_gt[1] = '.';
			}
		}
		elsif (@vcf_gt >= 2) {
			@vcf_gt = map { $_ eq '.' ? 0 : $_ } @vcf_gt;
		}

		$haplotype{$sample_id_name[$i] . 'L'} .= $vcf_gt[0];
		$haplotype{$sample_id_name[$i] . 'R'} .= $vcf_gt[1];
	}
	$weight .= '10'
}

for (my $i = 0; $i < 6; $i ++ ) {
	print ( $position[$i], $line_break );
}
for (my $i = 0; $i < $number_of_sample; $i ++) {
	print ( substr ( $sample_id_name[$i], 2), 'L ', $haplotype{$sample_id_name[$i] . 'L'}, '  1', $line_break);
	print ( substr ( $sample_id_name[$i], 2), 'R ', $haplotype{$sample_id_name[$i] . 'R'}, '  1', $line_break);
}

# printing lines of "weight," those are one empty line and ROUNDUP({number of sites} / 125) line(s).
# Maybe each weight consists of two letters.
print ($line_break);
while ( length ( $weight ) > 250 ) {
	print ( substr ( $weight, 0, 250, ''), $line_break);
}
print ($weight, $line_break);

