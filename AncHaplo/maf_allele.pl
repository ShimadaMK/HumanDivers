#!/usr/bin/perl -w

use strict;
use warnings;

if (@ARGV != 5) {
	print ("Usage: maf_allele.pl [in_VCF_file] [in_MAF_file] [chr_name] [reference_name] [output_names]\n");
	print ("[output_names] should be separated by comma without space.\n");
	exit;
}

my $filename_vcf = $ARGV[0];
my $filename_maf = $ARGV[1];
my $chr_name = uc ($ARGV[2]);
my @out_names = ($ARGV[3], split (',', $ARGV[4]));

my $line;
my @columns;

# load MAF alignment
open (my $fh_in, '<', $filename_maf) or die "Cannot open $filename_maf: $!";

for($line = <$fh_in>; $line =~ /^#/; $line = <$fh_in>) {}
if (!$line) {
	print STDERR ("No alignment in $filename_maf.\n");
	close ($fh_in);
	exit;
}

my %alignment;
my @align_pos;
my @align_length;
my $align_number = -1;

my $status = 'separator';
while ($line) {
	chomp ($line);
	@columns = split (/ +/, $line, 7);
	if ($status eq 'separator') {
		if ($line =~ /^a /) {
			$align_number ++;
			$status = 'alignment';
		}
	}
	elsif ($status eq 'alignment') {
		if ($line eq '') {
			$status = 'separator';
		}
		elsif ($columns[0] eq 's') {
			foreach my $name (@out_names) {
				if ($columns[1] eq $name) {
					$alignment{$name}[$align_number] = uc ($columns[6]);
				}
			}
			if ($columns[1] eq $out_names[0]) {
				$align_pos[$align_number] = $columns[2];
				$align_length[$align_number] = $columns[3];
			}
		}
	}
	$line = <$fh_in>;
}
close ($fh_in);
if ($align_number < 0) {
	print STDERR ("No alignment in $filename_maf.\n");
	exit;
}

# load vcf header
open (my $fh_vcf, '<', $filename_vcf) or die "Cannot open $filename_vcf: $!";
for ($line = <$fh_vcf>; $line =~ /^##/; $line = <$fh_vcf>) {}
if ($line !~ /^#CHROM/) {
	print STDERR ("No title line in $filename_vcf.\n");
	close ($fh_vcf);
	exit;
}

my $chr_vcf;
my $pos_vcf;
my $ref_vcf;
my $alt_vcf;
my $pos_now;
my $align_number_now = 0;
my $align_number_end = 0;
my %out_maf;
my $letter_start;
my $letter_end;
# load vcf body before specified chromosome
while ($line = <$fh_vcf>) {
	chomp ($line);
	($chr_vcf, undef) = split ("\t", $line, 2);
	if ($chr_vcf eq $chr_name) {
		last;
	}
}

# load vcf body of specified chromosome
print ("CHROM\tPOS\tREF\tALT");
foreach my $name (@out_names) {
	print ("\t", $name);
}
print ("\n");
while ($line) {
	chomp ($line);
	($chr_vcf, $pos_vcf, undef, $ref_vcf, $alt_vcf, undef) = split ("\t", $line, 6);
	if ($chr_vcf ne $chr_name) {
		last;
	}
	($align_number_end, $letter_start, $letter_end) = &find_align_number ($pos_vcf, $ref_vcf);
	if ($align_number_now > $align_number) {
		print STDERR ($filename_maf, " is shorter than ", $filename_vcf, "\n");
		last;
	}
	if ($align_number_end == -1) {
		foreach my $name (@out_names) {
			$out_maf{$name} = ".";
		}
	}
	elsif ($align_number_now == $align_number_end) {
		foreach my $name (@out_names) {
			$out_maf{$name} = substr ($alignment{$name}[$align_number_now], $letter_start, $letter_end - $letter_start);
		}
	}
	else {
		foreach my $name (@out_names) {
			$out_maf{$name} = substr ($alignment{$name}[$align_number_now], $letter_start);
		}
		for (my $i = $align_number_now + 1; $i < $align_number_end; $i ++) {
			foreach my $name (@out_names) {
				$out_maf{$name} .= $alignment{$name}[$i];
			}
		}
		foreach my $name (@out_names) {
			$out_maf{$name} .= substr ($alignment{$name}[$align_number_end], 0, $letter_end);
		}
	}
	print ($chr_vcf, "\t", $pos_vcf, "\t", $ref_vcf, "\t", $alt_vcf);
	foreach my $name (@out_names) {
		$out_maf{$name} =~ s/-//g;
		print ("\t", $out_maf{$name});
	}
	print ("\n");
	$line = <$fh_vcf>;
}

close ($fh_vcf);
exit;


sub find_align_number {
	my $pos_start = $_[0];
	my $ref = $_[1];
	my $pos_end = $pos_start + length ($ref) - 1;

	if ($align_pos[$align_number_now] >= $pos_start) {
		return (-1, 0, 0);
	}
	while ($align_pos[$align_number_now] + $align_length[$align_number_now] < $pos_start) {
		$align_number_now ++;
		if ($align_number_now > $align_number) {
			return (-1, 0, 0);
		}
	}
	my $align_number_end = $align_number_now;
	while ($align_pos[$align_number_end] + $align_length[$align_number_end] < $pos_end) {
		$align_number_end ++;
		if ($align_number_end > $align_number) {
			return (-1, 0, 0);
		}
	}

	my $alignment_ref = $alignment{$out_names[0]}[$align_number_now];
	my $letter_start = 0;
	for (my $base = $pos_start - $align_pos[$align_number_now] - 1; $base > 0; $base --) {
		$alignment_ref = substr ($alignment_ref, 1);
		$letter_start ++;
		while ($alignment_ref =~ /^-/) {
			$alignment_ref = substr ($alignment_ref, 1);
			$letter_start ++;
		}
	}
	$alignment_ref = $alignment{$out_names[0]}[$align_number_end];
	my $letter_end = 0;
	for (my $base = $pos_end - $align_pos[$align_number_end] - 1; $base >= 0; $base --) {
		$alignment_ref = substr ($alignment_ref, 1);
		$letter_end ++;
		while ($alignment_ref =~ /^-/) {
			$alignment_ref = substr ($alignment_ref, 1);
			$letter_end ++;
		}
	}
	return ($align_number_end, $letter_start, $letter_end);
}
