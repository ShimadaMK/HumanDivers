#!/usr/bin/perl -w

use strict;
use warnings;

my $len_seq = 6; # number of letters of a sequence name in RDF file
my $len_char = 6; # number of letters of a character name in RDF file

if (@ARGV < 1 || 3 < @ARGV) {
	print STDERR ("Usage: uniq_rdf.pl in_RDF_file [out_RDF_file [out_group_file]]\n");
	exit (1);
}
my $in_rdf_filename = $ARGV[0];
my $out_rdf_filename;
my $group_filename;
if (@ARGV == 1) {
	if ($in_rdf_filename =~ /^(.*)\.rdf$/) {
		$out_rdf_filename = $1 . '.uniq.rdf';
	}
	else {
		$out_rdf_filename = $in_rdf_filename . '.uniq.rdf';
	}
}
else {
	$out_rdf_filename = $ARGV[1];
}
if (@ARGV <= 2) {
	if ($out_rdf_filename =~ /^(.*)\.uniq\.rdf$/) {
		$group_filename = $1 . '.group';
	}
	elsif ($out_rdf_filename =~ /^(.*)\.rdf$/) {
		$group_filename = $1 . '.group';
	}
	else {
		$group_filename = $out_rdf_filename . '.group';
	}
}
else {
	$group_filename = $ARGV[2];
}

print ("Input RDF file: $in_rdf_filename\n");
print ("Output RDF file: $out_rdf_filename\n");
print ("Output group file: $group_filename\n\n");

# check filenames
my @op_check = &search_same_values ($in_rdf_filename, $out_rdf_filename, $group_filename);
if ($op_check[1] > 0) {
	my @err_filename = ('in_RDF_file', 'out_RDF_file', 'out_group_file');
	print STDERR ($err_filename[$op_check[0]], ' and ', $err_filename[$op_check[1]], " cannot be the same file.\n");
	exit (1);
}

# load input RDF file and make groups
open (my $fh_in, '<', $in_rdf_filename) or die "Cannot open $in_rdf_filename: $!";
my @header_lines;
my @footer_lines;
my $len_line = 0;
for (my $line_num = 0; $line_num < $len_char; $line_num ++) {
	if ($header_lines[$line_num] = <$fh_in>){
		$header_lines[$line_num] =~ s/\r?\n$//;
	}
	else {
		&rdf_error ($in_rdf_filename, $fh_in);
	}
	my $len_line_now = length ($header_lines[$line_num]);
	if ($len_line_now <= $len_seq) {
		&rdf_error ($in_rdf_filename, $fh_in);
	}
	elsif ($len_line == 0) {
		$len_line = $len_line_now;
	}
	elsif ($len_line_now != $len_line) {
		&rdf_error ($in_rdf_filename, $fh_in);
	}
}

my $name;
my @haplotypes;
my $haplotype;
my $group_otus;
my $freq;
my $line_break;
while (my $line = <$fh_in>) {
	$line =~ s/(\r?\n)$//;
	$line_break = $1;
	if ($line eq '') {
		last;
	}
	if (length ($line) < $len_line) {
		&rdf_error ($in_rdf_filename, $fh_in);
	}
	($name, $haplotype, $freq) = split (' ', $line, 3);
	if (!exists ($group_otus->{$haplotype})) {
		push (@haplotypes, $haplotype);
	}
	push (@{$group_otus->{$haplotype}}, $name);
}
while (my $line = <$fh_in>) {
	$line =~ s/\r?\n$//;
	push (@footer_lines, $line);
}

# output
open (my $fh_out_rdf, '>', $out_rdf_filename) or die "Cannot open $out_rdf_filename: $!";
for (my $line_num = 0; $line_num < $len_char; $line_num ++) {
	print $fh_out_rdf ($header_lines[$line_num], $line_break);
}
open (my $fh_out_group, '>', $group_filename) or die "Cannot open $group_filename: $!";
print $fh_out_group ("#Group name: members\n");

my $otu_name;
my $group_count = 1;
my $num_otu;
foreach $haplotype (@haplotypes) {
	if (($num_otu = @{$group_otus->{$haplotype}}) > 1) {
		$otu_name = sprintf ('H%d', $group_count ++);
		$otu_name = &check_special_member ($group_otus->{$haplotype}) . $otu_name;
		if ($haplotype !~ /[^0]/) {
			$otu_name = 'HR' . $otu_name;
		}
		print $fh_out_group ($otu_name, '-', $num_otu, ': ', join (',', @{$group_otus->{$haplotype}}), "\n");
	}
	else {
		$otu_name = $group_otus->{$haplotype}[0];
	}
	printf $fh_out_rdf ("%-7s%s  1%s", $otu_name, $haplotype, $line_break);
}
print $fh_out_rdf ($line_break);
foreach my $line (@footer_lines) {
	print $fh_out_rdf ($line, $line_break);
}

close ($fh_out_group);
close ($fh_out_rdf);

exit (0);

sub search_same_values
{
	my @in_v = @_;
	my $num_v = @in_v;
	if ($num_v <= 1) {
		return (0, 0);
	}
	my $i1;
	my $i2;
	for ($i1 = 0; $i1 < $num_v; $i1 ++) {
		for ($i2 = 0; $i2 < $i1; $i2 ++) {
			if ($in_v[$i1] eq $in_v[$i2]) {
				return ($i2, $i1);
			}
		}
	}
	return (0, 0);
}

sub rdf_error
{
	print STDERR ($_[0], " is invalid.\n");
	close ($_[1]);
	exit (2);
}

sub check_special_member
{
	my $members = $_[0];
	my @special_names  = ('DENISO');
	my @special_prefix = ('D'     );
	my $num_special = @special_names;
	my $prefix = '';
	for (my $num = 0; $num < $num_special; $num ++) {
		foreach my $member (@$members) {
			if ($special_names[$num] eq $member) {
				$prefix .= $special_prefix[$num];
				last;
			}
		}
	}
	return $prefix;
}
