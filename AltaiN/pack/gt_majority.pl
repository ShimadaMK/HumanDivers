#!/usr/bin/perl -w

use strict;
use warnings;

if (@ARGV > 1) {
	print ("Usage: gt_majority.pl [VCF_file]\n");
	exit (1);
}

my $vcf_file;
if (@ARGV == 0) {
	$vcf_file = &read_vcf_file (*STDIN);
}
else {
	my $in_vcf_filename = $ARGV[0];
	open (my $fh_vcf, '<', $in_vcf_filename) or die "Cannot open $in_vcf_filename:$!";
	$vcf_file = &read_vcf_file ($fh_vcf);
	close ($fh_vcf);
}
if ($vcf_file->{'error'}) {
	print STDERR ($vcf_file->{'error'}, " in the VCF file\n");
	exit (2);
}
my $number_of_sites = @{$vcf_file->{'ID'}};
my $number_of_samples = @{$vcf_file->{'sample_id'}};
my @site_valid = (1) x $number_of_sites;

# convert genotype (maybe unphased diploid) to haploid
for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
	if ($vcf_file->{'FILTER'}[$site_count] eq 'LowQual') {
		$site_valid[$site_count] = 0;
		next;
	}
	my $exist_alt = 0;
	for (my $sample_count = 0; $sample_count < $number_of_samples; $sample_count ++) {
		my %gt_field = &decode_gt ($vcf_file->{'FORMAT'}[$site_count], $vcf_file->{'GT'}[$site_count][$sample_count]);
		my $major_gt;
		if ($gt_field{'GT'} =~ /^(\d+)[\/\|](\d+)$/) {
			if ($1 eq $2) {
				$major_gt = $1;
			}
			else {
				$major_gt = &majority_gt ($vcf_file->{'REF_ALT'}[$site_count], \%gt_field, $1, $2);
			}
		}
		elsif ($gt_field{'GT'} =~ /^(\d+)$/) {
			$major_gt = $1;
		}
		else {
			$major_gt =  '.';
		}
		$vcf_file->{'GT'}[$site_count][$sample_count] = $major_gt . '|' . $major_gt;
		if ($major_gt !~ /^[0\.]$/) {
			$exist_alt = 1;
		}
	}
	$site_valid[$site_count] = $exist_alt;
}

print ("##fileformat=VCFv4.1\n");
print ("##source=gt_majority.pl\n");
print ("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Majority genotype based on number of reads, all described in homozygous\">\n");
print ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", join ("\t", @{$vcf_file->{'sample_id'}}), "\n");
for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
	if (!$site_valid[$site_count]) {
		next;
	}
	my $ref = $vcf_file->{'REF_ALT'}[$site_count][0];
	my $alt = (@{$vcf_file->{'REF_ALT'}[$site_count]} > 1) ? join (',', @{$vcf_file->{'REF_ALT'}[$site_count]}[1..$#{$vcf_file->{'REF_ALT'}[$site_count]}]) : '.';
	print ($vcf_file->{'CHR'}[$site_count], "\t", $vcf_file->{'POS'}[$site_count], "\t", $vcf_file->{'ID'}[$site_count], "\t");
	print ($ref, "\t", $alt, "\t.\t", $vcf_file->{'FILTER'}[$site_count], "\t.\tGT\t", join ("\t", @{$vcf_file->{'GT'}[$site_count]}));
	print ("\n");
}

exit (0);

sub read_vcf_file
{
	my $fh_vcf = $_[0];
	my $line_vcf;
	for ($line_vcf = <$fh_vcf>; $line_vcf =~ /^##/; $line_vcf = <$fh_vcf>) {}
	if ($line_vcf !~ /^#CHROM/) {
		my %ret_val = ('error' => 'No title line');
		return \%ret_val;
	}
	chomp ($line_vcf);
	my @columns = split ("\t", $line_vcf);
	if (@columns < 10) {
		my %ret_val = ('error' => 'No sample identifiers');
		return \%ret_val;
	}
	my @sample_ids = @columns[9 .. $#columns];
	my @chrs;
	my @poss;
	my @rsids;
	my @formats;
	my @filters;
	my @alleles;
	my @genotypes;
	while ($line_vcf = <$fh_vcf>) {
		chomp ($line_vcf);
		my ($chr, $pos, $rsid, $ref, $alt, undef, $filter, undef, $format, @gt) = split ("\t", $line_vcf);
		my @ref_alt = ($alt ne '.') ? ($ref, split (',', $alt)) : ($ref);
		push (@chrs, $chr);
		push (@poss, $pos);
		push (@rsids, $rsid);
		push (@alleles, [@ref_alt]);
		push (@filters, $filter);
		push (@formats, $format);
		push (@genotypes, [@gt]);
	}
	my %ret_val = ('sample_id' => \@sample_ids, 'CHR' => \@chrs, 'POS' => \@poss, 'ID' => \@rsids, 'REF_ALT' => \@alleles, 'FILTER' => \@filters, 'FORMAT' => \@formats, 'GT' => \@genotypes, 'error' => '');
	return \%ret_val;
}

sub decode_gt
{
	my ($format, $gt) = @_;
	my @formats = split (':', $format);
	my @gts = split (':', $gt);
	my %ret_val;
	while (@formats && @gts) {
		$ret_val{shift (@formats)} = shift (@gts);
	}
	return %ret_val;
}

sub majority_gt
{
	my ($ref_alt, $gt_field, @gts) = @_;
	if ($gts[0] >= @$ref_alt || $gts[1] >= @$ref_alt) {
		return '.';
	}
	elsif (length ($ref_alt->[$gts[0]]) > 1 || length ($ref_alt->[$gts[1]]) > 1) {
		return 0;
	}
	elsif (exists ($gt_field->{uc ($ref_alt->[$gts[0]])}) && exists ($gt_field->{uc ($ref_alt->[$gts[1]])})) {
		my @number_of_reads = (0, 0);
		foreach my $phase (0, 1) {
			if ($gt_field->{uc ($ref_alt->[$gts[$phase]])} =~ /^(\d+),(\d+)$/) {
				$number_of_reads[$phase] = $1 + $2;
			}
		}
		if ($number_of_reads[0] > $number_of_reads[1]) {
			return $gts[0];
		}
		elsif ($number_of_reads[0] < $number_of_reads[1]) {
			return $gts[1];
		}
		else {
			return ($gts[0] < $gts[1]) ? $gts[1] : $gts[0];
		}
	}
	else {
		return '.';
	}
}

