#!/usr/bin/perl -w

use strict;
use warnings;

if (@ARGV < 1 || 2 < @ARGV) {
	print STDERR ("Usage: vcf2beagle.pl VCF_file [Beagle_file_prefix]\n");
	exit (1);
}
my $in_vcf_filename = $ARGV[0];
my $out_file_prefix;
if (@ARGV == 1) {
	if ($in_vcf_filename =~ /^(.*)\.vcf$/) {
		$out_file_prefix = $1;
	}
	else {
		$out_file_prefix = $in_vcf_filename;
	}
}
else {
	$out_file_prefix = $ARGV[1];
}

print ("Input VCF file: $in_vcf_filename\n");

my $vcf_file = &read_vcf_file ($in_vcf_filename);

# sample grouping
my @haploid_samples_num = ();
my @phased_samples_num = ();
my @unphased_samples_num = ();
my @invalid_samples_num = ();
my $number_of_sites = @{$vcf_file->{'ID'}};
my $number_of_samples = @{$vcf_file->{'sample_id'}};
for (my $sample_count = 0; $sample_count < $number_of_samples; $sample_count ++) {
	my $phased = '0';# 0, Haploid, Phased, Unphased, X
	for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
		if ($vcf_file->{'GT'}[$site_count][$sample_count] =~ /^[\d+\.]$/) {
			if ($phased eq '0') {
				$phased = 'H';
			}
			elsif ($phased =~ /^[PU]$/) {
				$phased = 'X';
				last;
			}
		}
		elsif ($vcf_file->{'GT'}[$site_count][$sample_count] =~ /^[\d+\.]\|[\d+\.]$/) {
			if ($phased eq '0') {
				$phased = 'P';
			}
			elsif ($phased eq 'H') {
				$phased = 'X';
				last;
			}
		}
		elsif ($vcf_file->{'GT'}[$site_count][$sample_count] =~ /^[\d+\.]\/[\d+\.]$/) {
			if ($phased eq 'P' || $phased eq '0') {
				$phased = 'U';
			}
			elsif ($phased eq 'H') {
				$phased = 'X';
				last;
			}
		}
		else {
			$phased = 'X';
			last;
		}
	}
	if ($phased eq 'H') {
		push (@haploid_samples_num, $sample_count);
	}
	elsif ($phased eq 'P') {
		push (@phased_samples_num, $sample_count);
	}
	elsif ($phased eq 'U') {
		push (@unphased_samples_num, $sample_count);
	}
	elsif ($phased eq 'X') {
		push (@invalid_samples_num, $sample_count);
	}
}

# Output phased file
if (@haploid_samples_num + @phased_samples_num) {
	print ("Outputting $out_file_prefix.phased ...\n");
	open (my $fh_out, '>', $out_file_prefix . '.phased') or die "Cannot open $out_file_prefix.phased : $!";
	print $fh_out ("I\tid\t", join ("\t", map {($vcf_file->{'sample_id'}[$_], $vcf_file->{'sample_id'}[$_])} (@haploid_samples_num, @phased_samples_num)), "\n");
	for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
		my $rsid = $vcf_file->{'ID'}[$site_count];
		if ($rsid eq '.') {
			$rsid = $vcf_file->{'CHR'}[$site_count] . ':' . $vcf_file->{'POS'}[$site_count];
		}
		print $fh_out ("M\t$rsid");
		foreach my $sample_num (@haploid_samples_num) {
			my @allele = &gt_to_allele ($vcf_file->{'GT'}[$site_count][$sample_num], @{$vcf_file->{'REF_ALT'}[$site_count]});
			print $fh_out ("\t", $allele[0], "\t", $allele[0]);
		}
		foreach my $sample_num (@phased_samples_num) {
			my @allele = &gt_to_allele ($vcf_file->{'GT'}[$site_count][$sample_num], @{$vcf_file->{'REF_ALT'}[$site_count]});
			print $fh_out ("\t", $allele[0], "\t", $allele[1]);
		}
		print $fh_out ("\n");
	}
	close ($fh_out);
}

# Output unphased file
if (@unphased_samples_num > 0) {
	print ("Outputting $out_file_prefix.unphased ...\n");
	open (my $fh_out, '>', $out_file_prefix . '.unphased') or die "Cannot open $out_file_prefix.unphased : $!";
	print $fh_out ("I\tid\t", join ("\t", map {($vcf_file->{'sample_id'}[$_], $vcf_file->{'sample_id'}[$_])} (@unphased_samples_num)), "\n");
	for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
		my $rsid = $vcf_file->{'ID'}[$site_count];
		if ($rsid eq '.') {
			$rsid = $vcf_file->{'CHR'}[$site_count] . ':' . $vcf_file->{'POS'}[$site_count];
		}
		print $fh_out ("M\t$rsid");
		foreach my $sample_num (@unphased_samples_num) {
			my @allele = &gt_to_allele ($vcf_file->{'GT'}[$site_count][$sample_num], @{$vcf_file->{'REF_ALT'}[$site_count]});
			print $fh_out ("\t", $allele[0], "\t", $allele[1]);
		}
		print $fh_out ("\n");
	}
	close ($fh_out);
}

# Output marker file
if (@{$vcf_file->{'ID'}}) {
	print ("Outputting $out_file_prefix.marker ...\n");
	open (my $fh_out, '>', $out_file_prefix . '.marker') or die "Cannot open $out_file_prefix.marker : $!";
	for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
		my $rsid = $vcf_file->{'ID'}[$site_count];
			if ($rsid eq '.') {
			$rsid = $vcf_file->{'CHR'}[$site_count] . ':' . $vcf_file->{'POS'}[$site_count];
		}
		print $fh_out ($rsid, "\t", $vcf_file->{'POS'}[$site_count], "\t", join ("\t", @{$vcf_file->{'REF_ALT'}[$site_count]}), "\n");
	}
	close ($fh_out);
}

# print warning
if (@invalid_samples_num) {
	print ('Invalid genotype found in ', join (',', map {$vcf_file->{'sample_id'}[$_]} @invalid_samples_num), "\nThey have been ignored.\n");
}

exit (0);

sub read_vcf_file
{
	my $in_vcf_filename = $_[0];
	open (my $fh_vcf, '<', $in_vcf_filename) or die "Cannot open $in_vcf_filename: $!";
	my $line_vcf;
	for ($line_vcf = <$fh_vcf>; $line_vcf =~ /^##/; $line_vcf = <$fh_vcf>) {}
	if ($line_vcf !~ /^#CHROM/) {
		print STDERR ("No title line in $in_vcf_filename.\n");
		close ($fh_vcf);
		exit (5);
	}
	chomp ($line_vcf);
	my @columns = split ("\t", $line_vcf);
	if (@columns < 10) {
		print STDERR ("No sample identifiers in $in_vcf_filename.\n");
		close ($fh_vcf);
		exit (5);
	}
	my @sample_ids = @columns[9 .. $#columns];
	my @chrs;
	my @poss;
	my @rsids;
	my @alleles;
	my @genotypes;
	while ($line_vcf = <$fh_vcf>) {
		chomp ($line_vcf);
		my ($chr, $pos, $rsid, $ref, $alt, undef, undef, undef, undef, @gt) = split ("\t", $line_vcf);
		my @ref_alt = ($ref, split (',', $alt));
		push (@chrs, $chr);
		push (@poss, $pos);
		push (@rsids, $rsid);
		push (@alleles, [@ref_alt]);
		push (@genotypes, [@gt]);
	}
	my %ret_val = ('sample_id' => \@sample_ids, 'CHR' => \@chrs, 'POS' => \@poss, 'ID' => \@rsids, 'REF_ALT' => \@alleles, 'GT' => \@genotypes);
	return \%ret_val;
}

sub gt_to_allele
{
	my $gt_str = shift (@_);
	my @ref_alt = @_;
	my @genotype = split (/[\|\/]/, $gt_str);
	foreach my $gt (@genotype) {
		if ($gt =~ /^\d+$/ && $gt < @ref_alt) {
			$gt = $ref_alt[$gt];
		}
		else {
			$gt = 'N';
		}
	}
	return @genotype;
}
