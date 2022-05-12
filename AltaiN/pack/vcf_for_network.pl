#!/usr/bin/perl -w

# Removing (replacing genotype to missing value '.') whole haplotypes in the chromosomes including missing or unphased genotype.
# Then, (1) removing (replacing genotype to missing value '.') whole haplotypes in the chromosomes including the '3rd' genotype,
#   Or, (2) removing the sites that is triallelic or more.
# So the all sites in the output VCF is biallelic.

use strict;
use warnings;

if (@ARGV == 0 || 2 < @ARGV) {
	&usage;
	exit (1);
}
my $mode = shift (@ARGV);
if ($mode !~ /^(site|hap)$/) {
	&usage;
	exit (1);
}

# load VCF file
my $vcf_file;
if (@ARGV == 1) {
	my $in_vcf_filename = $ARGV[0];
	open (my $fh_vcf, '<', $in_vcf_filename) or die "Cannot open $in_vcf_filename: $!";
	$vcf_file = &read_vcf_file ($fh_vcf);
	close ($fh_vcf);
}
else {
	$vcf_file = &read_vcf_file (*STDIN);
}
if ($vcf_file->{'error'}) {
	print STDERR ($vcf_file->{'error'}, " in the VCF file\n");
	exit (2);
}
my $number_of_sites = @{$vcf_file->{'ID'}};
my @site_valid;
for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
	if (grep (!/^[ATGCN]+$/, @{$vcf_file->{'REF_ALT'}[$site_count]})) {
		$site_valid[$site_count] = 0;
		print STDERR ('chr', $vcf_file->{'CHR'}[$site_count], ':', $vcf_file->{'POS'}[$site_count], " is skipped.\n");
	}
	else {
		$site_valid[$site_count] = 1;
	}
}

# detect chromosomes including missing or unphased genotype
my $number_of_samples = @{$vcf_file->{'sample_id'}};
my $sample_remove = [];
for (my $sample_count = 0; $sample_count < $number_of_samples; $sample_count ++) {
	my @missing = (0, 0);
	my $remove = '';
	for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
		if (!$site_valid[$site_count]) {
			next;
		}
		my @alleles;
		(@alleles[0, 1], my $phase) = &genotype_allele ($vcf_file->{'GT'}[$site_count][$sample_count]);
		if ($alleles[0] eq '' || $alleles[0] gt $#{$vcf_file->{'REF_ALT'}[$site_count]}
		                      || $alleles[1] gt $#{$vcf_file->{'REF_ALT'}[$site_count]}) {
			$remove = 'X';
			last;
		}
		elsif ($phase eq '/') {
			$remove = 'U';
			last;
		}
		if ($alleles[0] eq '.') {
			$missing[0] = 1;
		}
		if ($alleles[1] eq '.') {
			$missing[1] = 1;
		}
	}
	$sample_remove->[$sample_count] = [0, 0];
	if ($remove) {
		$sample_remove->[$sample_count] = [1, 1];
		if ($remove eq 'X') {
			print STDERR ($vcf_file->{'sample_id'}[$sample_count], " contains an invalid genotype.\n");
		}
		elsif ($remove eq 'U') {
			print STDERR ($vcf_file->{'sample_id'}[$sample_count], " contains an unphased genotype.\n");
		}
		next;
	}
	if ($missing[0] == 1) {
		$sample_remove->[$sample_count][0] = 1;
		print STDERR ($vcf_file->{'sample_id'}[$sample_count], "-1 contains a missing genotype.\n");
	}
	if ($missing[1] == 1) {
		$sample_remove->[$sample_count][1] = 1;
		print STDERR ($vcf_file->{'sample_id'}[$sample_count], "-2 contains a missing genotype.\n");
	}
}

# replace whole genotypes to missing values, in chromosomes including missing or unphased genotype
&remove_gt;

# detect sites having more than 2 alleles
for (my $sample_count = 0; $sample_count < $number_of_samples; $sample_count ++) {
	@{$sample_remove->[$sample_count]}[0, 1] = (0, 0);
}
my @site_remove = (0) x $number_of_sites;
for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
	my $number_of_alleles = @{$vcf_file->{'REF_ALT'}[$site_count]};
	if (!$site_valid[$site_count] || $number_of_alleles <= 2) {
		next;
	}
	print STDERR ('chr', $vcf_file->{'CHR'}[$site_count], ':', $vcf_file->{'POS'}[$site_count], ' The "ALT" alleles are more than 1.', "\n");
	my @allele_counts = (0) x $number_of_alleles;
	for (my $sample_count = 0; $sample_count < $number_of_samples; $sample_count ++) {
		my @alleles;
		(@alleles[0, 1], undef) = &genotype_allele ($vcf_file->{'GT'}[$site_count][$sample_count]);
		if ($alleles[0] =~ /^\d+$/) {
			$allele_counts[$alleles[0]] ++;
		}
		if ($alleles[1] =~ /^\d+$/) {
			$allele_counts[$alleles[1]] ++;
		}
	}
	my @alleles_dec = sort {$allele_counts[$b] <=> $allele_counts[$a]} (0 .. $number_of_alleles - 1);
	if ($alleles_dec[0] > $alleles_dec[1]) {
		@alleles_dec[0, 1] = @alleles_dec[1, 0];
	}
	if ($alleles_dec[0] > 0) {
		print STDERR ('-- Warning: The "REF" allele has been changed: ', $vcf_file->{'REF_ALT'}[$site_count][0],
		                                                         ' -> ', $vcf_file->{'REF_ALT'}[$site_count][$alleles_dec[0]], ".\n");
	}
	print STDERR ('-- The "ALT" alleles have been changed: ', join (',', @{$vcf_file->{'REF_ALT'}[$site_count]}[1 .. $number_of_alleles - 1]),
	                                                               ' -> ', $vcf_file->{'REF_ALT'}[$site_count][$alleles_dec[1]], ".\n");
	@{$vcf_file->{'REF_ALT'}[$site_count]} = map {$vcf_file->{'REF_ALT'}[$site_count][$_]} @alleles_dec[0, 1];
	for (my $sample_count = 0; $sample_count < $number_of_samples; $sample_count ++) {
		my @alleles = 0;
		(@alleles[0, 1], my $phase) = &genotype_allele ($vcf_file->{'GT'}[$site_count][$sample_count]);
		if ($alleles[0] eq $alleles_dec[0]) {
			$alleles[0] = 0;
		}
		elsif ($alleles[0] eq $alleles_dec[1]) {
			$alleles[0] = 1;
		}
		elsif ($alleles[0] ne '.') {
			$alleles[0] = '.';
			$sample_remove->[$sample_count][0] = 1;
			print STDERR ('-- ', $vcf_file->{'sample_id'}[$sample_count], "-1 has the 3rd allele.\n");
			$site_remove[$site_count] = 1;
		}
		if ($alleles[1] eq $alleles_dec[0]) {
			$alleles[1] = 0;
		}
		elsif ($alleles[1] eq $alleles_dec[1]) {
			$alleles[1] = 1;
		}
		elsif ($alleles[1] ne '.' && $alleles[1] ne '') {
			$alleles[1] = '.';
			$sample_remove->[$sample_count][1] = 1;
			print STDERR ('-- ', $vcf_file->{'sample_id'}[$sample_count], "-2 has the 3rd allele.\n");
			$site_remove[$site_count] = 1;
		}
		$vcf_file->{'GT'}[$site_count][$sample_count] = $alleles[0] . $phase . $alleles[1];
	}
}

# replace whole genotypes to missing values, in chromosomes including the 3rd genotype
if ($mode eq 'hap') {
	&remove_gt;
	@site_remove = (0) x $number_of_sites;
}

# output VCF
print ("##fileformat=VCFv4.1\n");
print ("##source=vcf_for_network.pl\n");
print ("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Allele Count\">\n");
print ("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Alternate Allele Count\">\n");
print ("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
print ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", join ("\t", @{$vcf_file->{'sample_id'}}), "\n");
for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
	if (!$site_valid[$site_count]) {
		next;
	}
	if ($site_remove[$site_count]) {
		print STDERR ('More than 2 alleles appears at chr', $vcf_file->{'CHR'}[$site_count], ':', $vcf_file->{'POS'}[$site_count], ". It has been removed.\n");
		next;
	}
	my @allele_counts = (0, 0);
	for (my $sample_count = 0; $sample_count < $number_of_samples; $sample_count ++) {
		my @alleles;
		(@alleles[0, 1], undef) = &genotype_allele ($vcf_file->{'GT'}[$site_count][$sample_count]);
		if ($alleles[0] =~ /^(0|1)$/) {
			$allele_counts[$1] ++;
		}
		if ($alleles[1] =~ /^(0|1)$/) {
			$allele_counts[$1] ++;
		}
	}
	print STDOUT ($vcf_file->{'CHR'}[$site_count], "\t",
	              $vcf_file->{'POS'}[$site_count], "\t",
	              $vcf_file->{'ID'}[$site_count], "\t",
	              $vcf_file->{'REF_ALT'}[$site_count][0], "\t",
	              $vcf_file->{'REF_ALT'}[$site_count][1], "\t.\t.\tAN=", $allele_counts[0] + $allele_counts[1], ";AC=", $allele_counts[1], "\tGT\t",
	              join ("\t", @{$vcf_file->{'GT'}[$site_count]}), "\n");
}

exit 0;

sub usage
{
	print STDERR ("Usage: vcf_for_network.pl <site|hap> [VCF_file]\n");
	print STDERR ("<site> remove the triallelic (or more) sites.\n");
	print STDERR ("<hap> remove whole haplotypes in the chromosomes including neither REF nor ALT genotype.\n");
	return;
}

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
	my %ret_val = ('sample_id' => \@sample_ids, 'CHR' => \@chrs, 'POS' => \@poss, 'ID' => \@rsids, 'REF_ALT' => \@alleles, 'GT' => \@genotypes, 'error' => '');
	return \%ret_val;
}

sub genotype_allele
{
	my $genotype = $_[0];
	my @alleles;
	my $phase;
	if ($genotype =~ /^(\.|[\d]+)(\||\/)(\.|[\d]+)(:.*)?$/) {
		@alleles = ($1, $3);
		$phase = $2;
	}
	elsif ($genotype =~ /^(\.|[\d]+)(:.*)?$/) {
		@alleles = ($1, '');
		$phase = '';
	}
	else {
		@alleles = ('', '');
		$phase = '';
	}
	return (@alleles, $phase);
}

sub remove_gt
{
	for (my $sample_count = 0; $sample_count < $number_of_samples; $sample_count ++) {
		if ($sample_remove->[$sample_count][0] + $sample_remove->[$sample_count][1] == 2) {
			for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
				if (!$site_valid[$site_count]) {
					next;
				}
				my @alleles;
				(@alleles[0, 1], my $phase) = &genotype_allele ($vcf_file->{'GT'}[$site_count][$sample_count]);
				$vcf_file->{'GT'}[$site_count][$sample_count] = ($phase eq '') ? '.' : '.' . $phase . '.';
			}
		}
		elsif ($sample_remove->[$sample_count][0] == 1) {
			for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
				if (!$site_valid[$site_count]) {
					next;
				}
				my @alleles;
				(@alleles[0, 1], my $phase) = &genotype_allele ($vcf_file->{'GT'}[$site_count][$sample_count]);
				$vcf_file->{'GT'}[$site_count][$sample_count] = '.' . $phase . $alleles[1];
			}
		}
		elsif ($sample_remove->[$sample_count][1] == 1) {
			for (my $site_count = 0; $site_count < $number_of_sites; $site_count ++) {
				if (!$site_valid[$site_count]) {
					next;
				}
				my @alleles;
				(@alleles[0, 1], my $phase) = &genotype_allele ($vcf_file->{'GT'}[$site_count][$sample_count]);
				$vcf_file->{'GT'}[$site_count][$sample_count] = ($phase eq '') ? $alleles[0] : $alleles[0] . $phase . '.';
			}
		}
	}
}

__END__
error code
1: argument error
2: VCF file error
