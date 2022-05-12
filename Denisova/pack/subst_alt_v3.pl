#!/usr/bin/perl -w

# Version 3

# For more information
#  see /media/usb1/mak2/helper/nishida/human_genome/spec_subst_alt_v3.txt

use strict;
use warnings;

if (@ARGV < 4 || @ARGV > 5) {
	print ("Usage: subst_alt_v3.pl [in_file] [vcf_file] [out_file] [sample_name] ([range])\n");
	exit;
}

my $filename_in = $ARGV[0];
my $filename_vcf = $ARGV[1];
my $filename_out = $ARGV[2];
my @sample_ids = split (',', $ARGV[3]);
my $start_pos = '';
my $end_pos = '';
if (@ARGV == 5 && $ARGV[4] =~ /^([\d]*)-([\d]*)$/) {
	$start_pos = $1;
	$end_pos = $2;
}
my $line;
my $line_vcf;
my @column_vcf;
my $chr_fa;
my $pos_min;
my $pos_max;
my $offset = 9; #number of column before each sample's data in vcf.

# load DNA sequence 
open (my $fh_in, '<', $filename_in) or die "Cannot open $filename_in: $!";

while ($line = <$fh_in>) {
	if ($line =~ /^>/) {
		last;
	}
}
if (!$line) {
	print ("No title line in $filename_in.\n");
	close ($fh_in);
	exit;
}
if ($line =~ /range=chr(..?):([\d]*)-([\d]*)/ ) {
	$chr_fa  = $1;
	$pos_min = $2;
	$pos_max = $3;
}
else {
	print ("The title line of $filename_in is invalid.\nIt must contain range=chr_:____-____\n");
	close ($fh_in);
	exit;
}

# set start and end position
if ($start_pos !~ /^\d+$/ || $start_pos < $pos_min) {
	$start_pos = $pos_min;
}
if ($end_pos !~ /^\d+$/ || $end_pos > $pos_max) {
	$end_pos = $pos_max;
}
if ($start_pos > $end_pos) {
	close ($fh_in);
	exit;
}

my $padding = $start_pos - $pos_min;
my $interval_length = $end_pos - $start_pos + 1;

# load all sequences
my $haplotype = '';
while (length($haplotype) < $padding) {
	$line = <$fh_in>;
	if (!$line || $line eq "\n") {
		print ("$filename_in is too short.\n");
		close ($fh_in);
		exit;
	}
	else {
		chomp ($line);
		$haplotype .= $line;
	}
}
$haplotype = substr($haplotype, $padding);
while (length($haplotype) < $interval_length) {
	$line = <$fh_in>;
	if (!$line || $line eq "\n") {
		print ("$filename_in is too short.\n");
		close ($fh_in);
		exit;
	}
	else {
		chomp ($line);
		$haplotype .= $line;
	}
}
$haplotype = substr($haplotype, 0, $interval_length);
close ($fh_in);
if ($haplotype =~ /([^ACGTNacgt])/) {
	print ("Invalid letter \"$1\" in $filename_in.\n");
	exit;
}

# load vcf header
open (my $fh_vcf, '<', $filename_vcf) or die "Cannot open $filename_vcf: $!";
for($line_vcf = <$fh_vcf>; $line_vcf =~ /^##/; $line_vcf = <$fh_vcf>) {}
if ($line_vcf !~ /^#CHROM/) {
	print ("No title line in $filename_vcf.\n");
	close ($fh_vcf);
	exit;
}
chomp ($line_vcf);
my %sample_id_column_number;
@column_vcf = split ("\t", $line_vcf);
for (my $col = $offset; $col < @column_vcf; $col ++) {
	$sample_id_column_number{$column_vcf[ $col ]} = $col;
}

if ((@sample_ids == 1) && ($sample_ids[0] eq '--all')) {
	@sample_ids = @column_vcf[$offset..$#column_vcf];
}

my @valid_sample_ids;
for (my $id=0; $id < @sample_ids; $id ++) {
	if (exists $sample_id_column_number{$sample_ids[$id]}) {
		push (@valid_sample_ids, $sample_ids[$id]);
	}
	else {
		print ("$sample_ids[$id] does not exist in $filename_vcf.\n");
	}
}
if (@valid_sample_ids == 0) {
	print ("No specified sample in $filename_vcf.\n");
	close ($fh_vcf);
	exit;
}

# load vcf body before start position
my $pos = 0;
my $chr = 0;
while ($chr ne $chr_fa || $pos < $start_pos) {
	$line_vcf = <$fh_vcf>;
	if (!$line_vcf) {
		close ($fh_vcf);
		exit;
	}
	@column_vcf = split ("\t", $line_vcf);
	$chr = $column_vcf[0];
	$pos = $column_vcf[1];
}

# alternate allele
my $cur_pos;
my $ref;
my @alts;
my $alt;
my $org_base;
my $indel_length;
my $haplotypes;
my %plos;
foreach my $id (@valid_sample_ids) {
	$cur_pos->{$id}[0] = $start_pos - 1;
	$cur_pos->{$id}[1] = $start_pos - 1;
	$indel_length->{$id}[0] = 0;
	$indel_length->{$id}[1] = 0;
	$haplotypes->{$id}[0] = $haplotype;
	$haplotypes->{$id}[1] = $haplotype;
	$plos{$id} = 'N';
}
for ( ; $line_vcf ; $line_vcf = <$fh_vcf>) {
	chomp ($line_vcf);
	@column_vcf = split ("\t", $line_vcf);
	$chr = $column_vcf[0];
	$pos = $column_vcf[1];
	$ref = $column_vcf[3];
	$alt = $column_vcf[4];
	@alts = split (',', $alt);
	if ($pos > $end_pos || $chr ne $chr_fa) {
		last;
	}
	if ($ref eq '.' || $alt eq '.' || $alt =~ /^<.*>$/) {
		next;
	}
	if ($ref =~ /[^ACGTN]/) {
		print ("Invalid letter(s) in REF at pos=$pos.\nREF = $ref\n");
		close ($fh_vcf);
		exit;
	}
	foreach my $id (@valid_sample_ids) {
		if ($plos{$id} eq 'U') {
			next;
		}
		my $genotype = $column_vcf[$sample_id_column_number{$id}];
		($genotype, undef) = split (':', $genotype, 2);
		my ($plo, @is_alt) = &read_genotype ($genotype);
		if ($plo eq 'E') {
			print ("$id: Invalid genotype at pos=$pos.\nGT=$genotype\n");
			close ($fh_vcf);
			exit;
		}
		elsif ($plo eq '.') {
			next;
		}
		elsif ($plo eq 'U') {
			$plos{$id} = 'U';
			next;
		}
		elsif ($plo eq 'D') {
			if ($plos{$id} eq 'H') {
				print ("$id: Haploid diploid confusing at pos=$pos.\n");
				exit;
			}
			$plos{$id} = $plo;
		}
		elsif ($plo eq 'H') {
			if ($plos{$id} eq 'D') {
				print ("$id: Diploid haploid confusing at pos=$pos.\n");
				exit;
			}
			$plos{$id} = $plo;
		}
		for (my $phase = 0; $phase < 2; $phase ++) {
			if ($plo eq 'H' && $phase == 1) {
				next;
			}
			if ($is_alt[$phase] eq '.' || $is_alt[$phase] eq '0') {
				next;
			}
			if ($is_alt[$phase] <= @alts) {
				$alt = $alts[$is_alt[$phase] - 1];
			}
			else {
				$alt = $ref;
				$alt =~ s/./N/g;
			}
			if ($pos <= $cur_pos->{$id}[$phase]) {
				print ("$id-", $phase + 1, ": duplicate alternation at pos=$pos. It is ignored.\n");
				next;
			}
			$org_base = substr ($haplotypes->{$id}[$phase], $pos - $start_pos + 	$indel_length->{$id}[$phase], length($ref), $alt);
			$indel_length->{$id}[$phase] += length($alt) - length($ref);
			$cur_pos->{$id}[$phase] = $pos + length($ref) - 1;
			if ($pos == $cur_pos->{$id}[$phase]) {
				print ("$id-", $phase + 1, ": $org_base is altered $ref -> $alt at pos $pos.\n");
			}
			else {
				print ("$id-", $phase + 1, ": $org_base is altered $ref -> $alt at pos $pos - $cur_pos->{$id}[$phase].\n");
			}
		}
	}
}
close ($fh_vcf);

# output
open (my $fh_out, '>', $filename_out) or die "Cannot open $filename_out: $!";
foreach my $id (@valid_sample_ids) {
	if ($plos{$id} eq 'U') {
		print ("$id\tunphased.\n");
		next;
	}
	for (my $phase = 0; $phase < 2; $phase ++) {
		if ($plos{$id} eq 'H' && $phase == 1) {
			next;
		}
		print ("$id-", $phase + 1, "\tlength=", $interval_length + $indel_length->{$id}[$phase], "\n");
		print $fh_out (">$id-", $phase + 1, "\n");
		while (length($haplotypes->{$id}[$phase]) > 50) {
			print $fh_out (substr ($haplotypes->{$id}[$phase], 0, 50, ''), "\n");
		}
		print $fh_out ($haplotypes->{$id}[$phase], "\n");
	}
}

close ($fh_out);
exit (0);

sub read_genotype {
	my $gt = $_[0];
	if ($gt eq '.') {
		return ('.');
	}
	if ($gt =~ /^(\d+)$/) {
		return ('H', $gt);
	}
	if ($gt =~ /^(\d+|\.)([\|\/])(\d+|\.)$/) {
		if ($2 eq '/' && $1 ne $3) {
			return ('U');
		}
		return ('D', $1, $3);
	}
	return ('E');
}
